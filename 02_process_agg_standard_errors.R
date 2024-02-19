# 01c_process_to_standard_errors
library(tidyverse)

## read data -----
droptrial <- read_csv("../cleaned_data/Data/excluded_trials_look_up.csv")
out <- readRDS("Scratch_data/agg_hba1c.Rds")
bas <- readRDS("Scratch_data/agg_hba1c_base.Rds")
bas <- map(bas, ~ .x %>% 
             rename(nct_id = trial_id, arm_id_unq  = arm_id))
comporig <- read_csv("../cleaned_data/Data/example_comparisons.csv")
comp2019 <- read_csv("../cleaned_data/Data/comparisons_2019.csv")
comp2022 <- read_csv("../cleaned_data/Data/comparisons_2022_update.csv")
comp <- bind_rows(comporig = comporig,
                      comp2019 = comp2019 %>% mutate(unique_id = as.character(unique_id)),
                      comp2022 = comp2022 %>% mutate(unique_id = as.character(unique_id)), .id = "type")
comp %>% 
  distinct(type, comp_id, arm_id)
comp <- comp %>% 
  distinct(comp_id, arm_id, .keep_all = TRUE)
arm <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv")

## NCT01541956  is actually a standard deviation not a standard error. checked paper
out$arm$meta <- out$arm$meta %>% 
  mutate(dispersion_type = if_else(nct_id == "NCT01541956", "sd", dispersion_type))

## calculate standard errors for arm data (change and end) ----
out$unav <- NULL
out <- transpose(out)
dt <- out$data %>% 
  bind_rows(.id = "stat_type")
meta <- out$meta %>% 
  bind_rows(.id = "stat_type")
rm(out)

## note all with SD are either arm or end
hba1c_sd <- dt %>% 
    semi_join(meta %>% 
                 filter(dispersion_type %in% c("sd")) %>% 
                 select(dispersion_type, result_id) %>% 
                 unnest(result_id))
hba1c_sd %>% count(stat_type)
hba1c_sd <- hba1c_sd %>% 
  left_join(bas$ns) %>% 
  mutate(participants = if_else(is.na(participants), n, participants)) %>% 
  select(-n) %>% 
  left_join(bas$ns_sg)
hba1c_sd %>% 
  filter(is.na(participants))
## 4 trials with missing n data. Pull manually. one where still NA could not locate
hba1c_sd_msng <- read_csv("nct_id,arm_id_unq,arm_label,result,dispersion,participants
UMIN000019022,ac0182,alogliptin,-0.5,0.7,NA
UMIN000019022,ac0182,alogliptin,-0.5,0.7,NA
UMIN000019022,ac0183,vildagliptin,-0.7,0.9,NA
UMIN000019022,ac0183,vildagliptin,-0.7,0.9,NA
UMIN000022953,ac0198,linagliptin,-0.5,0.4,21
UMIN000022953,ac0199,metformin,-0.5,0.9,22
NCT04183868,updab0089,glimepiride,-0.46,0.6,13
NCT04183868,updab0090,empagliflozin,-0.52,0.7,13
UMIN000024663,updac0190,vildagliptin,6.86,1.2,20
UMIN000024663,updac0191,metformin,7.12,0.8,20")
hba1c_sd <- hba1c_sd %>% 
  left_join(hba1c_sd_msng %>% 
              filter(!nct_id == "UMIN000019022") %>% 
              select(arm_id_unq, ppts_new = participants)) %>% 
  mutate(participants = if_else(is.na(participants), ppts_new, participants)) %>% 
  select(-ppts_new)

hba1c_sd <- hba1c_sd %>% 
  mutate(participants = if_else(is.na(participants), n, participants)) %>% 
  select(-n) %>% 
  mutate(se = as.double(dispersion)/participants^0.5) %>% 
  select(-dispersion)
hba1c_se <- dt %>% 
  semi_join(meta %>% 
              filter(dispersion_type %in% c("se", "sem")) %>% 
              select(dispersion_type, result_id) %>% 
              unnest(result_id)) %>% 
  mutate(se = as.double(dispersion)) %>% 
  select(-dispersion)
hba1c_se %>% count(stat_type)

#
hba1c_rng <- dt %>% 
  semi_join(meta %>% 
              filter(dispersion_type %in% c("95%ci", "95% ci", "ci")) %>% 
              select(dispersion_type, result_id) %>% 
              unnest(result_id))
# arm, comparison and end
hba1c_rng %>% 
  count(stat_type)
hba1c_rng <- hba1c_rng %>% 
  separate(dispersion, into = c("ll", "ul"), sep = ",|\\;")
hba1c_rng <- hba1c_rng %>% 
  mutate(across(c(ul, ll), as.double)) %>% 
  mutate(se = (ul-ll)/(2*1.96))  %>% 
  select(-ll, -ul)

nodisp <- dt %>% 
  inner_join(meta %>% 
  filter(!dispersion_type %in% c("se", "sd", "95%ci", "95% ci", "ci", "sem")) %>% 
    unnest(result_id) %>% 
    select(stat_type, nct_id, result_id, dispersion_type))
## 19 trials with no dispersion information
nodisp %>% 
  filter(is.na(dispersion_type)) %>% 
  count(nct_id)
## two trials with nonconvertible dispersion information
nodisp %>% 
  filter(!is.na(dispersion))
## so 21 trials in total dropped

## Join all data together as same structure baseline characteristics
hba1c <- bind_rows(rng = hba1c_rng,
                   sd = hba1c_sd,
                   se = hba1c_se,
                   .id = "dispersion_type_original")
rm(hba1c_rng, hba1c_sd, hba1c_se, dt)
## Pull arm-level data
hba1c_arm <- hba1c %>% 
  filter(stat_type %in% c("arm", "end"))
## aggregate over subgroups - 6 trials in CTG
hba1c_arm_sg <- hba1c_arm %>% 
  filter(!is.na(arm_id_subgroup)) %>% 
  group_by(dispersion_type_original, stat_type, nct_id, arm_id_unq, comp_id, ancova, units_label, timepoint) %>% 
  summarise(result_id = paste(result_id, collapse = ","),
              arm_id_subgroup = paste(arm_id_subgroup, collapse = ","),
              result = weighted.mean(result, participants),
              se = weighted.mean(se^2, participants)^0.5,
              participants = sum(participants, na.rm = TRUE)) %>% 
  ungroup()
hba1c_arm <- bind_rows(hba1c_arm %>% filter(is.na(arm_id_subgroup)),
                       hba1c_arm_sg)
rm(hba1c_arm_sg)

hba1c_arm <- hba1c_arm %>% 
  left_join(bas$ns) %>% 
  left_join(bas$sex %>% select(-n)) %>% 
  left_join(bas$age)
hba1c_arm <- hba1c_arm %>% 
  mutate(participants = if_else(!is.na(participants), participants, n),
         n = if_else(!is.na(n), n, participants))
## identify where has trial, but not arm with baseline data
## One trial has data for two arms but not two arms (these are called "ignore_diet" 
## remaining 5 trials missing for all trials in set)
xmn <- hba1c_arm %>% filter(is.na(n)) %>% pull(nct_id) %>% unique()
hba1c_arm %>% filter(nct_id %in% xmn) %>% inner_join(arm %>% select(arm_id_unq, arm_label))
hba1c_arm %>% filter(nct_id %in% xmn)  %>% count(nct_id)

## create structures with arm-level for baseline characteristics and arm-level style for comparison, where
## the result for one arm (the reference arm) is NA
# age and sex over treatment comparisons, to create a single row per comparison
hba1c_comp <- hba1c %>% 
  filter(stat_type == "comp")
comp <- comp %>% 
  filter(comp_id %in% hba1c_comp$comp_id)
comp2 <- comp %>% 
  rename(arm_id_unq = arm_id, nct_id = trial_id) %>% 
  left_join(bas$ns) %>% 
  left_join(bas$sex) %>% 
  left_join(bas$age)
## Two trial missing data. Add. age_sd missing for the latter
comp_msng <- read_csv(
"arm_id_unq,arm_label,n,male,age_m,age_sd
unq_updaa10055,dapagliflozin,12,7,63,4.7
unq_updaa10057,placebo,14,13,64.6,4.7
noarmid1,placebo,87,44,52,
noarmid2,gemigliptin,87,57,54,")
comp_msng_meta <- comp2 %>% 
  semi_join(comp_msng %>% select(arm_id_unq)) %>% 
  select(unique_id:treat_or_ref) %>% 
  inner_join(comp_msng)
comp3 <- bind_rows(comp2 %>% 
                     filter(!nct_id %in% comp_msng_meta$nct_id),
                   comp_msng_meta)

#differences in mean outcomes between arms)."
comp4 <- comp3 %>% 
  inner_join(hba1c_comp %>% 
               select(-arm_id_unq)) %>% 
  mutate(result = if_else(treat_or_ref == "reference", NA_real_, result),
         se = if_else(treat_or_ref == "reference", NA_real_, se))

## as per help menu on set_agd_contrast "Each study should have a single
#reference/baseline treatment, against which relative effects in the other
#arm(s) are given. For the reference arm, include a data row with continuous
#outcome y equal to NA. If a study has three or more arms (so two or more
#relative effects), set the standard error se for the reference arm data row
#equal to the standard error of the mean outcome on the reference arm (this
#determines the covariance of the relative effects, when expressed as # only one
#trial has three or more arms - NCT03159052 
# assume standard deviation arm same as median for other trials with a placebo arm
## note drop double reference for this trial
comp4 <- comp4 %>% 
  filter(!(nct_id == "NCT03159052" & unique_id == "0058"))

## convert all onto common units
## note error  in assigning unit label to one trial - Netherlands Trial Register NTR6709
## Need to add N for a trial as well updac0148 arm 
# 15 participants, updac0149  arm 18 participants (from https://doi.org/10.1111/dom.12936)
comp5 <- comp4 %>% 
  mutate(result = if_else(units_label == "mmol/mol", 
                          result * 0.09148 + 2.152,
                          result),
         ## note that conversion needs to happen at level of se
         se     = if_else(units_label == "mmol/mol", 
                          (0.09148*(se*n^0.5) + 2.152)/n^0.5,
                          se),
         units_label = if_else(units_label == "mmol/mol", "percentage_converted", units_label))
hba1c_arm2 <- hba1c_arm %>% 
  mutate(n = case_when(
    arm_id_unq == "updac0148" ~ 15L,
    arm_id_unq == "updac0149" ~ 18L,
    TRUE ~ n),
    result = if_else(units_label == "mmol/mol", 
                          result * 0.09148 + 2.152,
                          result),
         ## note that conversion needs to happen at level of se
         se     = if_else(units_label == "mmol/mol", 
                          (0.09148*(se*n^0.5) + 2.152)/n^0.5,
                          se),
         units_label = if_else(units_label == "mmol/mol", "percentage_converted", units_label))

placebo_arms <- arm %>% 
  filter(arm_label == "placebo") %>% 
  pull(arm_id_unq) %>% 
  unique()
placebo_arms <- hba1c_arm2 %>% 
  filter(arm_id_unq %in% placebo_arms) %>% 
  select(nct_id, se, n) %>% 
  distinct(nct_id, .keep_all = TRUE) %>% 
  mutate(sd = se*n^0.5) %>% 
  summarise(sd = median(sd, na.rm = TRUE)) %>% 
  pull(sd)
comp5 <- comp5 %>% 
  mutate(se = if_else(treat_or_ref == "reference", placebo_arms/n^0.5, se))

## set final datasets for saving as csv files
comp6 <- comp5 %>% 
  select(nct_id, arm_id_unq, treat_or_ref, result, se, n, male, age_m, age_sd, arm_id_subgroup)
hba1c_arm3 <- hba1c_arm2 %>% 
  select(nct_id, arm_id_unq, result, se, n, male, age_m, age_sd, arm_id_subgroup) %>% 
  mutate(treat_or_ref = "arm_level_outcome")

final <- bind_rows(comp5,
                   hba1c_arm3)
write_csv(final, "Data/agg.csv")