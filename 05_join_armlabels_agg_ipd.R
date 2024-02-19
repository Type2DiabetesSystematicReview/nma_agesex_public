library(tidyverse)
library(multinma)

ipd_method <- "rsd"

source("Scripts/00_functions.R")
source("../common_functions/Scripts/combine_sd.R")

## read in arm labels ----
arm_meta <- read_csv("Data/arm_labels_hba1c.csv") 
## arm_id_subgroup has mismatched ordering between arm meta and agg so rename
arm_meta <- arm_meta %>% 
  rename(arm_id_subgroup_ordering = arm_id_subgroup)

## read in agg data - drop those not in arm label metadata ---
agg <- read_csv("Data/agg.csv") 
agg_orig <- agg 
sum(!duplicated(agg$nct_id))
agg <- agg %>% 
  semi_join(arm_meta %>% select(nct_id, arm_id_unq))
## 15 agg not in arm meta, expected due to combinations etc
agg_drop <- agg_orig %>% 
  filter(!nct_id %in% agg$nct_id) %>% 
  group_by(nct_id) %>% 
  summarise(arms = sum(!duplicated(arm_id_unq)),
            participants = sum(n)) %>% 
  ungroup()
write_csv(agg_drop, "Outputs/Trials_agg_not_in_meta.csv")
rm(agg_drop)
## no duplicates, nil dropped  on join
agg <- agg %>% 
  inner_join(arm_meta)
agg %>% filter(!is.na(arm_id_subgroup) | !is.na(arm_id_subgroup_ordering)) %>% select(starts_with("arm_id"))
agg <- agg %>% 
  select(-arm_id_subgroup, -arm_id_subgroup_ordering)
sum(!duplicated(agg$nct_id))

## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
sum(!duplicated(ipd$nct_id))
ipd <- ipd %>% 
  semi_join(arm_meta)
sum(!duplicated(ipd$nct_id))
## read in reg
reg <- readRDS("Scratch_data/ipd_coefs_frmttd.Rds")
sum(!duplicated(reg$nct_id))
reg <- reg %>% 
  rename(arm_id_unq = trt) %>% 
  semi_join(ipd %>% select(nct_id))
sum(!duplicated(reg$nct_id))

## check same in ipd and reg
ipd %>% count(nct_id, nct_id2, arm_id_unq, arm_f)
reg %>% filter(!is.na(arm_id_unq)) %>% count(nct_id, nct_id2, arm_id_unq, arm_f)

## all consistent in IPD and ipd reg
## note expect to see NAs in reg as these are for non-treatment covariates (ie not treatment effects or treamtnet covariate interactions)
ipd <- ipd %>% 
  inner_join(arm_meta)
sum(!duplicated(ipd$nct_id))
reg <- reg %>% 
  semi_join(arm_meta %>% select(nct_id)) 
reg <- reg %>% 
  left_join(arm_meta %>% select(nct_id, arm_id_unq, arm_lvl, trtcls5, trtcls4, drug_regime_smpl))
sum(!duplicated(ipd$nct_id))
sum(!duplicated(reg$nct_id))
setdiff(ipd$arm_id_unq, reg$arm_id_unq)
setdiff(reg$arm_id_unq, ipd$arm_id_unq)

## recode data into required format for IPD
if (ipd_method == "rds") {
  ipd <- ipd %>% 
    rename(result = value_2_rsd, base = value_1_rsd) %>% 
    select(-age10)
} else {
  ipd <- ipd %>% 
    rename(result = value_2_se, base = value_1_se) %>% 
    select(-age10)
}

## impute N based on standard errors for 3 trials where do not have these
agg %>% 
  filter(is.na(n)) %>% 
  count(nct_id)
# 24 trials without age
agg %>% 
  filter(is.na(age_m)) %>% 
  count(nct_id)
# 32 trials without sex
agg %>% 
  filter(is.na(male)) %>% 
  count(nct_id)
aggsd <- agg %>% 
  mutate(sd = se*n^0.5) %>% 
  pull(sd)
write_lines(paste0("The median standard deviation for HbA1c is", median(aggsd, na.rm = TRUE)), "Outputs/impute_missing_n_using_sd.txt")
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

## impute first within trial, then across trials for age. taking the mean
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(across(c(age_m, age_sd), function(x) if_else(is.na(x),
                                                      mean(x, na.rm = TRUE), 
                                                      x))) %>% 
  ungroup() %>% 
  mutate(across(c(age_m, age_sd), function(x) if_else(is.na(x),
                                                      mean(x, na.rm = TRUE), 
                                                      x))) 

## repeat this imputation for sex
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(male = if_else(is.na(male), n*mean(male/n, na.rm = TRUE), male)) %>%
  ungroup() %>% 
  mutate(male = if_else(is.na(male), n*mean(male/n, na.rm = TRUE), male)) 


## Aggregate arms where have multiple arms with same drug or drug and dose in case of metformin and novel antidiabetics) or class in case of insulin
## but not where results were supplied as subgroups as this was done in earlier script. Use Cochrane formula
# this applies to 25 trials.
nrow(agg)
arm_pre <- agg %>% 
  group_by(nct_id, arm_lvl) %>% 
  summarise(n = length(arm_lvl),
            across(c(arm_id_unq, drug_dose_orig), ~ paste(.x, collapse = ","))) %>% 
  ungroup() %>% 
  filter(n >=2)
write_csv(arm_pre, "Outputs/aggregate_arms_collapsed.csv")
rm(arm_pre)
agg <- agg %>% 
  mutate(sd = se*n^0.5) %>% 
  group_by(nct_id, arm_lvl, trtcls5, trtcls4, drug_regime_smpl) %>% 
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

## next need to drop any single arm trials ----
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(n_arms = sum(!duplicated(arm_id_unq))) %>% 
  ungroup()
agg %>% 
  group_by(n_arms) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = length(nct_id))
## Note are 3 single arm trials now. Drop
single <- agg %>% 
  filter(n_arms ==1)
write_csv(single, "Outputs/Trials_arms_dropped_single_after_aggregation.csv")
agg <- agg %>% 
  filter(n_arms >=2)

## split into networks and plot ----
ipd <- ipd %>% 
  group_by(drug_regime_smpl) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(ipd = data)

## following code needed to deal with main effects not having a drug regime drug code etc
reg2 <- reg %>% 
  group_by(nct_id) %>% 
  mutate(drug_regime_smpl = drug_regime_smpl[!is.na(drug_regime_smpl)][1]) %>% 
  ungroup

## drop main effects which are not for the drug regime smpl models
## None dropped
reg2 <- reg2 %>% 
                 group_by(nct_id, nct_id2, models) %>% 
                 mutate(anytreat = any(!is.na(arm_id_unq))) %>% 
                 ungroup() %>% 
                 filter(anytreat) %>% 
                 select(-anytreat)
## drop contrasts which are not in the main models. Both from correlation
## and from terms
reg2_arm_drp <- reg2 %>% 
  filter(models == "f1", !is.na(arm_id_unq)) %>% 
  distinct(nct_id, term, arm_id_unq, arm_f, arm_lvl)
reg2_arm_drp <- reg2_arm_drp %>% 
  filter(is.na(arm_lvl)) %>% 
  pull(arm_f) %>% 
  paste(collapse = "|")
reg2 <- reg2 %>% 
  filter(is.na(term) | !str_detect(term, reg2_arm_drp))
reg2$crl<- map(reg2$crl, ~ .x[!str_detect(colnames(.x), reg2_arm_drp),
                              !str_detect(rownames(.x), reg2_arm_drp)])
reg <- reg2
rm(reg2)
reg <- reg %>% 
  group_by(drug_regime_smpl) %>% 
  nest() %>% 
  ungroup()
agg <- agg %>% 
  group_by(drug_regime_smpl) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(agg = data)
agg <- agg %>% 
  filter(!drug_regime_smpl == "unclear or missing")
tot <- agg %>% 
  inner_join(ipd) %>% 
  inner_join(reg %>% rename(reg = data))
rm(ipd, agg)

saveRDS(tot, "Scratch_data/agg_ipd_hba1c.Rds")

                       


