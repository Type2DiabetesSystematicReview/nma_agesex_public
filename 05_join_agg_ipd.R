library(tidyverse)
source("Scripts/00_functions.R")

## note that this is the version of multinma from 
# install_github("git@github.com:dmcalli2/multinma.git")
library(multinma)

source("../common_functions/Scripts/combine_sd.R")

## read in network characterisation ----
arm_regime <- readRDS("Scratch_data/arms_assign_drug_regime.Rds")

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

## identify where age missing and impute the median values
# 5 trials with missing sex data for a single arm, but not the whole trial
msng <- agg %>% 
  group_by(nct_id) %>% 
  mutate(some_missing = any(is.na(male)) & !all(is.na(male))) %>% 
  ungroup() %>% 
  filter(some_missing)
# these are also missing from the raw data so not a coding error
msng <- msng %>% left_join(readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
                             select(nct_id = trial_id, arm_id_unq = arm_id, first, second, variable) %>% filter( variable == "male") )
# 5 trials with missing age data for a singole arm, but not the whole trial
msng2 <- agg %>% 
  group_by(nct_id) %>% 
  mutate(some_missing = any(is.na(age_m)) & !all(is.na(age_m))) %>% 
  ungroup() %>% 
  filter(some_missing)
msng2 <- msng2 %>% left_join(readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>%
                               select(nct_id = trial_id, arm_id_unq = arm_id, first, second, variable) %>% filter( variable == "age") )
## 3 trials this is true for age adn sex, 2 trials only age and two trials only sex
union(msng2$nct_id, msng$nct_id)
intersect(msng2$nct_id, msng$nct_id)
setdiff(union(msng2$nct_id, msng$nct_id), intersect(msng2$nct_id, msng$nct_id))

## impute first within trial, then across trials for age
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(across(c(age_m, age_sd), function(x) if_else(is.na(x),
                                                      mean(x, na.rm = TRUE), 
                                                      x))) %>% 
  ungroup() %>% 
  mutate(across(c(age_m, age_sd), function(x) if_else(is.na(x),
                                                      mean(x, na.rm = TRUE), 
                                                      x)))
## repeat for sex
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(male = if_else(is.na(male), n*mean(male/n, na.rm = TRUE), male)) %>%
  ungroup() %>% 
  mutate(male = if_else(is.na(male), n*mean(male/n, na.rm = TRUE), male))
## pool arms which have same drug (eg metformin 500mg and metforming 1000g) into a single arm ----
## use cochrane formula. This gets rid of 161 rows
nrow(agg)
agg <- agg %>% 
  mutate(sd = se*n^0.5) %>% 
  group_by(nct_id, treat_or_ref, drug_code, trtcls5, trtcls4, drug_regime, drug_regime_smpl) %>% 
  summarise(pooled_arm = length(n),
            sd = CombSdVectorised(n = n, m = result, s = sd),
            result = weighted.mean(result, n),
            age_sd = CombSdVectorised(n = n, m = age_m, s = age_sd),
            age_m = weighted.mean(age_m, n),
            male = sum(male),
            n = sum(n),
            arm_id_unq = paste(arm_id_unq, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(se = sd/n^0.5)
nrow(agg)
## 115 trials/122 arms are arms that were pooled
agg %>% 
  filter(pooled_arm >=2) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = sum(!duplicated(arm_id_unq)))

## next need to drop resultant single arm trials ----
agg <- agg %>% 
  group_by(drug_regime_smpl, nct_id) %>% 
  mutate(n_arms = sum(!duplicated(arm_id_unq))) %>% 
  ungroup()
agg %>% 
  group_by(n_arms) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = length(nct_id))
single <- agg %>% 
  filter(n_arms ==1)
## all are trials now collapsed to a single arm
single_rv <- arm_regime %>% 
  filter(nct_id %in% single$nct_id)
agg <- agg %>% 
  filter(n_arms >=2)
## split into networks and plot ----

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
rm(ipd, agg)

tot$dm_nets <- map2(tot$ipd, tot$agg, ~ RptNetwork(.x, .y))
pdf("Outputs/disconnected_network.pdf", height = 20, width = 20)
map2(tot$dm_nets, tot$drug_regime_smpl, ~ plot(.x, layout = "auto") +
       ggtitle(.y))
dev.off()

tot$agg <- map(tot$agg, ~ .x %>% 
                 filter(!nct_id %in% c("NCT02752113", "NCT04196231", "NCT02787551")))
saveRDS(tot, "data_for_mars.Rds")
saveRDS(tot, "Scratch_data/agg_ipd_hba1c.Rds")

## generate stancode and standata
tot$mdl <- map(tot$dm_nets, ~ nma(.x, 
                                # regression = ~ base + age*.trt + sex*.trt,
                                trt_effects = "random", 
                                link = "identity", 
                                likelihood = "normal",
                                class_interactions = "common",
                                prior_intercept = normal(scale = 20),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                QR = TRUE, 
                                cores = 4))


## create network plot just showing classes
tot$cls_nets <- map2(tot$ipd, tot$agg, ~ {
  x <- .x %>% 
    distinct(nct_id, nct_id2, trtcls5, .keep_all = TRUE) %>% 
    group_by(nct_id, nct_id2) %>% 
    mutate(arm_n = length(age)) %>% 
    ungroup() %>% 
    filter(arm_n >=2)
  RptNetworkClass(x, .y)
})

tot$cls_nets_plt <- map(tot$cls_nets, plot)
saveRDS(tot %>% select(drug_regime_smpl, cls_nets_plt), "Scratch_data/network_class_plots.Rds")

                       


