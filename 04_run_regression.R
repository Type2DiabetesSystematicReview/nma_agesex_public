## Run regression models
library(tidyverse)
## note that this is the version of multinma from 
# install_github("git@github.com:dmcalli2/multinma.git")
library(multinma)

## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- read_csv("Data/agg.csv")
## read in arm drug names, doses, regiments
drug <- readRDS("Scratch_data/drug_names_doses_regimen.Rds")
## read in arm assignment
arm_assign <- readRDS("Scratch_data/arm_labels_final.Rds")

## correct the IPD arm label in one IPD trial uaa10735 instead of uaa10734

## Add arm label to IPD only one not in is where set to placebo
# One trial dropped as both amrs are the same drug
sum(!duplicated(ipd$nct_id))
ipd %>% 
  anti_join(arm_assign) %>% 
  count(nct_id, arm_id_unq)
ipd <- ipd %>% 
  inner_join(arm_assign)
sum(!duplicated(ipd$nct_id))

## add arm label to aggregate - 14 trials dropped are all trials with same drug spread across arms
sum(!duplicated(agg$nct_id))
agg %>% 
  mutate(in_arm_assign = arm_id_unq %in% arm_assign$arm_id_unq) %>%
  group_by(nct_id) %>% 
  mutate(in_arm_assign = any(in_arm_assign)) %>% 
  ungroup() %>% 
  filter(!in_arm_assign) %>% 
  distinct(nct_id, arm_id_unq) %>% 
  left_join(drug) %>% 
  group_by(nct_id) %>% 
  summarise(n_arms = length(arm_id_unq),
            drug_name  = paste(drug_name  %>% unique() %>% sort, collapse = ","),
            drug2_name = paste(drug2_name %>% unique() %>% sort, collapse = ",")) %>% 
  ungroup()
agg <- agg %>% 
  inner_join(arm_assign) 

## recode data into required format for IPD
ipd <- ipd %>% 
  rename(result = value_2, base = value_1) %>% 
  select(-age10)


## One IPD trial did not pull out "open label" comparison from IPD. Exct design of this trial unclear
rmv1 <- c("NCT01368081")
ipd <- ipd %>%
  filter(!nct_id %in% rmv1) 

## Remove ones which have only single arm in data. Is an issue in "raw" data. Ask ELB to review
rmv2 <- c("NCT03387683", "ChiCTR1800015296", "jRCTs031180159", "NCT02875821", "TCTR20170511001", "UMIN000021177", "UMIN000031451", "IRCT20200722048176N1", "NCT03313752")
agg <- agg %>%
  filter(!nct_id %in% rmv2)
# review NCT00838903 as has negative standard errors AND they are the same as the effect estimates
# correct directly from CTG
agg <- agg %>% 
  mutate(se = case_when(
    arm_id_unq == "uaa10383" ~ 0.113,
    arm_id_unq == "uaa10384" ~ 0.065,
    arm_id_unq == "uaa10382" ~ 0.064,
    arm_id_unq == "uaa10381" ~ 0.065,
    TRUE ~ se
  ))
# drop one trial with no standard errors (was a trial registered on a Japanese register. presented mean and SD, no participant count at all. No P-value found)
# drop another trial with no standard errors - NCT00101673 - not available anywhere
agg <- agg %>% 
  filter(!nct_id %in% c("NCT00101673", "UMIN000019022") )

## impute N based on standard errors for 4 trials where do not have these
agg %>% 
  filter(is.na(n)) %>% 
  count(nct_id)
# 31 trials without age
agg %>% 
  filter(is.na(age_m)) %>% 
  count(nct_id)
# 39 trials without sex
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

dm_net <- RptNetwork(ipd, agg)
pdf("Outputs/Initial_network.pdf", height = 20, width = 20)
plot(dm_net, layout = "auto")
dev.off()

# NCT02477865 and NCT02477969 have placebo arms but no data for these
agg <- agg %>% 
  filter(!nct_id %in% c("NCT04196231",
                        "NCT02477969",
                        "NCT02477865"))

dm_net <- RptNetwork(ipd, agg)
pdf("Outputs/Connected_network.pdf", height = 20, width = 20)
plot(dm_net, layout = "auto")
dev.off()

## drop missing age and sex
agg2 <- agg %>% 
  filter(!is.na(age_m),
         !is.na(age_sd))
dm_net <- RptNetwork(ipd, agg2)
dm_net <- add_integration(dm_net,
                          age = distr(qnorm, mean = age_m, sd = age_sd),
                          n_int = 10)
dm_net_fe <- nma(dm_net, 
                 regression = ~ age*.trt,
                 trt_effects = "fixed", link = "identity", likelihood = "normal", 
                 init_r = 0.1,
                 QR = TRUE, cores = 8)
dm_net_fe$code <- rstan::get_stancode(dm_net_fe$object)
saveRDS(dm_net_fe, "Scratch_data/fe_model.Rds")
