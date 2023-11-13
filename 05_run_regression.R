## Run regression models
library(tidyverse)
## note that this is the version of multinma from 
# install_github("git@github.com:dmcalli2/multinma.git")
library(multinma)

## read in network characterisation ----
arm_regime <- readRDS("Scratch_data/arms_assign_drug_regime.Rds")
write_csv(arm_regime, "Outputs/arms_drugs_regimes.csv")

## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- read_csv("Data/agg.csv")
## note this is used for examining readon for exclusions only
drug <- readRDS("Scratch_data/drug_names_doses_regimen.Rds")
## correct the IPD arm label in one IPD trial uaa10735 instead of uaa10734
## Add arm label to IPD only one not in is where set to placebo
ipd <- ipd %>% 
  filter(!nct_id %in% c("NCT00306384", "NCT01368081"))
# One trial dropped as both arms are the same drug, other dropped as same drug (empa) against an open label extension
sum(!duplicated(ipd$nct_id))
ipd %>% 
  anti_join(arm_regime %>% select(nct_id)) %>% 
  count(nct_id, arm_id_unq)
## note, many to many expected as same trial can appear across multiple arms
ipd <- ipd %>% 
  inner_join(arm_regime, relationship = "many-to-many")
sum(!duplicated(ipd$nct_id))

## add arm label to aggregate - 29 trials dropped need to review again. Talk with ELB
sum(!duplicated(agg$nct_id))
agg %>% 
  mutate(in_arm_regime = arm_id_unq %in% arm_regime$arm_id_unq) %>%
  group_by(nct_id) %>% 
  mutate(in_arm_regime = any(in_arm_regime)) %>% 
  ungroup() %>% 
  filter(!in_arm_regime) %>% 
  distinct(nct_id, arm_id_unq) %>% 
  left_join(drug) %>% 
  group_by(nct_id) %>% 
  summarise(n_arms = length(arm_id_unq),
            drug_name  = paste(drug_name  %>% unique() %>% sort, collapse = ","),
            drug2_name = paste(drug2_name %>% unique() %>% sort, collapse = ",")) %>% 
  ungroup()
agg <- agg %>% 
  inner_join(arm_regime) 

## recode data into required format for IPD
ipd <- ipd %>% 
  rename(result = value_2, base = value_1) %>% 
  select(-age10)

# review NCT00838903 as has negative standard errors AND they are the same as the effect estimates
# no longer any need as these arms in the IPD
# correct directly from CTG
# drop one trial with no standard errors (was a trial registered on a Japanese register. presented mean and SD, 
# no participant count at all. No P-value found)
# drop another trial with no standard errors - UMIN000019022 - not available anywhere
agg <- agg %>% 
  filter(!nct_id %in% c("UMIN000019022") )

## impute N based on standard errors for 4 trials where do not have these
agg %>% 
  filter(is.na(n)) %>% 
  count(nct_id)
# 24 trials without age
agg %>% 
  filter(is.na(age_m)) %>% 
  count(nct_id)
# 38 trials without sex
agg %>% 
  filter(is.na(male)) %>% 
  count(nct_id)
aggsd <- agg %>% 
  mutate(sd = se*n^0.5) %>% 
  pull(sd)
hist(aggsd, breaks =50)
median(aggsd, na.rm = TRUE)
# 0.93
agg <- agg %>% 
  mutate(n_prov = if_else(is.na(n),
                     "imputed",
                     "observed"),
         n = if_else(!is.na(n), 
                     n,
                     (0.93/se)^2) %>% round(),
         n = pmax(1, n))

## split into networks and plot ----
RptNetwork <- function(ipd_choose, agg_choose){
  combine_network(
    set_ipd(ipd_choose,
            study = nct_id2,
            trt = drug_code,
            y = result,
            trt_class = trtcls5),
    set_agd_arm(agg_choose %>% 
                  filter(treat_or_ref == "arm_level_outcome"),
                study = nct_id,
                trt = drug_code,
                y = result, 
                se = se,
                trt_class = trtcls5,
                sample_size = n),
    set_agd_contrast(agg_choose %>% 
                       filter(!treat_or_ref == "arm_level_outcome"),
                     study = nct_id,
                     trt = drug_code,
                     y = result, 
                     se = se,
                     trt_class = trtcls5,
                     sample_size = n))  
}
ipd <- ipd %>% 
  group_by(drug_regime_smpl) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(ipd = data)
agg <- agg %>% 
  group_by(drug_regime_smpl) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(agg = data)
agg <- agg %>% 
  filter(!drug_regime_smpl == "unclear or missing")
tot <- agg %>% 
  inner_join(ipd)
rm(agg, ipd)

## identify where age missing and impute the median values
tot$agg <- map(tot$agg, ~ .x %>% 
                 mutate(across(c(age_m, age_sd), function(x) if_else(is.na(x), mean(x, na.rm = TRUE), x))))
saveRDS(tot, "Scratch_data/data_for_mars.Rds")


## draw network
tot$dm_nets <- map2(tot$ipd, tot$agg, ~ RptNetwork(.x, .y))
pdf("Outputs/network.pdf", height = 20, width = 20)
map2(tot$dm_nets, tot$drug_regime_smpl, ~ plot(.x, layout = "auto") +
       ggtitle(.y))
dev.off()
tot$ipd_trials <- map_int(tot$ipd, ~ sum(!duplicated(.x$nct_id)))
tot$agg_trials <- map_int(tot$agg, ~ sum(!duplicated(.x$nct_id)))

tot$dm_nets <- map(tot$dm_nets, ~ add_integration(.x,
                                        age = distr(qnorm, mean = age_m, sd = age_sd),
                                        n_int = 10))
tot$dm_net_fe <- map(tot$dm_nets, ~ nma(.x, 
                 regression = ~ age*.trt,
                 trt_effects = "fixed", link = "identity", likelihood = "normal", 
                 init_r = 0.1,
                 QR = TRUE, cores = 8))
tot$dm_net_fe <- map(tot$dm_net_fe, ~ {
  .x$code <-  rstan::get_stancode(.x$object)
  .x
})
saveRDS(tot, "Scratch_data/fe_model.Rds")
