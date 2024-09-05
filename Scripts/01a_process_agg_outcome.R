#01a_process_hba1c_initial
library(tidyverse)
source("Scripts/common_functions/Scripts/misc.R")
source("Scripts/00_functions.R")
## read in all data ----
## read in drop trial list ----
# droplist <- read_csv("Data/cleaned_data/Data/excluded_trials_look_up.csv")
## Read in exclusions table so can add to within data
exclusions <- read_csv("Data/exclusions.csv")
droplist <- exclusions %>% 
  filter(exclude == 1L)

## Create folder structure if it does already not exist ----
foldermake <- c("Outputs", "FromVM", "Scratch_data")
walk(foldermake, ~ if (!file.exists(.x)) dir.create(.x))

## read in ipd so can drop from aggregate ----
ipd1 <- read_csv("Data/agesexhba1c_6115/hba1c_base_change_overall.csv") %>% 
  pull(nct_id)
ipd2 <- read.csv("Data/gsk/hba1c_base_change_overall.csv") %>% 
  pull(nct_id)
ipd3 <- read.csv("Data/agesexhba1c_8697/hba1c_base_change_overall.csv") %>% 
  pull(nct_id)
ipd4 <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv") %>% 
  pull(nct_id)
ipd <- c(ipd1, ipd2, ipd3, ipd4) %>% unique()
## read in aggregate level data for each hba1c trials and baseline data
hba1c_agg1 <- read_csv("Data/cleaned_data/Data/hba1c_outcome_data_2019.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg2 <- read_csv("Data/cleaned_data/Data/hba1c_outcome_data_2022.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg3 <- read_csv("Data/cleaned_data/Data/hba1c_outcome_data_2024.csv") %>% 
  rename(nct_id = trial_id)
arm <- read_csv("Data/cleaned_data/Data/arm_data_all_cleaned.csv")

## clean arm data 
arm <- arm %>% 
  mutate(nct_id = trial_id,
         drug_name = case_when(
    arm_label == "placebo" & is.na(drug_name) ~ "placebo",
    TRUE ~ drug_name
  ))

## read in outcome counts
oc_orig <- read_csv("Data/aact/oc_orig.csv")
oc_new <- read_csv("Data/aact/oc_new.csv")
## read in IDs for mapping counts onto outcomes
hba1c_ids_new <- read_csv("Data/extract_transform/hba1c_ids_new.csv")

# assign each type of outcome to a wider category. Eg mean change, baseline etc
result_determination <- read_csv("Created_metadata/resolve_result_type.csv")

## initial cleaning ----
hba1c_agg <- bind_rows(extract1 = hba1c_agg1,
                       extract2 = hba1c_agg2,
                       extract3 = hba1c_agg3,
                       .id = "extract")
rm(hba1c_agg1, hba1c_agg2, hba1c_agg3)

## pull duration data
duration <- hba1c_agg %>% 
  select(nct_id, timepoint, timepoint_units) 
duration <- duration %>% 
  mutate(timepoint_units = if_else(nct_id == "NCT02104804" & is.na(timepoint_units), "weeks", timepoint_units))
setdiff(ipd, duration$nct_id)
duration <- duration %>% 
  mutate(weeks = case_when(
    timepoint_units == "weeks" ~ timepoint,
    timepoint_units == "days" ~ timepoint/7,
    timepoint_units == "months" ~ 52*timepoint/12
  )) 
duration <- duration  %>% 
  arrange(desc(weeks)) %>% 
  distinct(nct_id, .keep_all = TRUE)
## 4 IPD mace trials without duration recorded here
addmis <- tibble(nct_id = c("NCT00968708", "NCT01989754", "NCT02065791", "NCT01131676"),
                 weeks = c(3.4, 3, 4.6, 4.6)*52) 
duration <- bind_rows(duration, addmis) %>% 
  mutate(timeraw = paste(timepoint, timepoint_units)) %>% 
  select(-timepoint, -timepoint_units)
shrt <- duration %>% 
  filter(weeks <12) %>% 
  pull(nct_id)
saveRDS(duration, "Scratch_data/trial_duration.Rds")
rm(duration, addmis)

# drop ipd trials from agg
hba1c_agg <- hba1c_agg %>% 
  filter(!nct_id %in% ipd)
# drop if on excluded list - one trial
hba1c_agg %>% 
  filter(nct_id %in% droplist$trial_id)
hba1c_agg <- hba1c_agg %>% 
  filter(!nct_id %in% droplist$trial_id)
# drop if not a selected arm
arm_orig <- arm
arm <- arm %>% 
  filter(is.na(arm_label) | arm_label != "total") %>% 
  select(trial_id, arm_id_unq, arm_id_subgroup, arm_id) %>% 
  gather("armids", "value", -trial_id) %>% 
  select(-armids) %>% 
  distinct()
hba1c_agg <- hba1c_agg %>% 
  semi_join(arm %>% select(nct_id = trial_id,
                           arm_id_unq = value))

## correct errors where result type missing
hba1c_agg <- hba1c_agg %>% 
  mutate(result_type = if_else(nct_id == "NCT00509262" & result_description == "change from baseline", "mean", result_type))
# cant use if_else here or drops dispersion where arm_id_unq is missing. so use case_when
hba1c_agg <- hba1c_agg %>% 
  mutate(dispersion = case_when(
    arm_id_unq == "unq_updaa10037" ~ "0,0.5",
    TRUE ~ dispersion))

## separate out result metadata ----
hba1c_meta <- hba1c_agg %>% 
  select(nct_id, extract, analysis_population, outcome_harm_label, result_description, result_type, 
         dispersion_type, timepoint_units, source, result_id) %>% 
  distinct() %>% 
  nest(result_id = result_id)
hba1c_agg <- hba1c_agg %>% 
  select(nct_id, result_id, arm_id_unq, comp_id, ancova, result, dispersion,
         units_label, timepoint, arm_id_subgroup) %>% 
  distinct()
hba1c_meta %>% 
  count(result_type)
hba1c_meta <- hba1c_meta %>% 
  inner_join(result_determination)
result_determination_rv <- hba1c_meta %>% 
  count(nct_id, result_type_smry) %>% 
  spread(result_type_smry, n, fill = 0L)
result_type_per_trial <- hba1c_meta %>% 
  group_by(nct_id) %>% 
  summarise(result_type_best = case_when(
    any(result_type_smry == "mean_base") & any(result_type_smry == "mean_end") ~ "mean_change_calculate",
    any(result_type_smry == "mean_change") ~ "mean_change",
    any(result_type_smry == "between_arm_mean") ~ "between_arm_mean",
    TRUE ~ "other"
  ))
result_type_per_trial %>% 
  count(result_type_best) %>% 
  write_csv("Outputs/Types of Hba1c result in each trial.csv")

## errors in result determination
## note is not error in processing but difference in way CTG data translated to aact. New extract performed to address this
hba1c_agg <- hba1c_agg %>% 
  mutate(dispersion = if_else(nct_id == "NCT02973477" & dispersion == "0.01-0.08",
                              "0.01;0.08", dispersion))
ci_error <- hba1c_meta %>% 
  filter(dispersion_type %in% c("95%ci", "95% ci")) %>% 
  unnest(result_id) %>% 
  inner_join(hba1c_agg) %>% 
  filter(!str_detect(dispersion, ",|\\;"))

ci_corr <- read_csv("Data/extract_transform/ci_corr.csv")
ci_error2 <- ci_error %>% 
  inner_join(ci_corr %>% 
               select(nct_id, result = param_value, dispersion_lower_limit, dispersion_upper_limit))
ci_error2 <- ci_error2 %>% 
  mutate(dispersion = paste(dispersion_lower_limit, dispersion_upper_limit, sep = ";")) %>% 
  select(-dispersion_lower_limit, -dispersion_upper_limit)
hba1c_agg <- bind_rows(hba1c_agg %>% 
                         anti_join(ci_error %>% select(result_id)),
                        ci_error2 %>% 
                          select(nct_id, result_id, arm_id_unq, comp_id, ancova, result, dispersion, units_label, timepoint, arm_id_subgroup))
rm(ci_error, ci_corr, ci_error2)

### Add in number of participants where available based on result_id ----
## pull other missing ns from results  note the IDs are AACT-extract specific
# hba1c_continuous_aact.csv
hba1c_ids_orig <-  read_csv("Data/extract_transform/hba1c_continuous_aact.csv")
oc_orig <- oc_orig %>% 
  semi_join(hba1c_ids_orig %>% select(nct_id))
oc_orig <- oc_orig %>% 
  inner_join(hba1c_ids_orig %>% rename(id_hba1c_ids = id))
oc_new <- oc_new %>% 
  semi_join(hba1c_ids_new %>% select(nct_id))
oc_new <- oc_new %>% 
  inner_join(hba1c_ids_new %>% rename(id_hba1c_ids = id))

hba1c_meta_orig <- hba1c_meta %>% 
  filter(nct_id %in% oc_orig$nct_id)
hba1c_meta_new <- hba1c_meta %>% 
  filter(nct_id %in% oc_new$nct_id, !nct_id %in% oc_orig$nct_id)
hba1c_meta_non <- hba1c_meta %>% 
  filter(!nct_id %in% c(oc_orig$nct_id, oc_new$nct_id))
hba1c_meta_orig <- hba1c_meta_orig %>% 
  unnest(result_id) %>% 
  mutate(id_hba1c_ids = as.integer(result_id)) %>% 
  left_join(oc_orig) 
## there are none with matching IDs for this set of result IDs
hba1c_meta_new <- hba1c_meta_new %>% 
  unnest(result_id) %>% 
  mutate(id_hba1c_ids = as.integer(result_id)) %>% 
  semi_join(oc_new) 
## no duplicates
hba1c_meta_orig %>% 
  select(nct_id, result_id, id_hba1c_ids, id, outcome_id, result_group_id) %>% 
  duplicated() %>% 
  any()

## check N's
hba1c_meta_ns_orig <- hba1c_meta_orig %>% 
  rename(participants = count) %>% 
  filter(!is.na(participants)) 
hba1c_meta_ns_orig <- hba1c_meta_ns_orig %>% 
  group_by(nct_id, result_id)  %>% 
  summarise(participants_multi = participants %>% unique() %>% sort() %>% paste(collapse = ","),
            participants_min = min(participants),
            participants_max = max(participants),
            participants = mean(participants)) %>% 
  ungroup()
# all 827 have a single set of participants so can add to hba1c results
# need to do the same for the new extract
hba1c_meta_ns_orig %>%
  filter(!participants_min == participants_max)
hba1c_meta_ns_orig <- hba1c_meta_ns_orig %>% 
  group_by(nct_id, result_id)  %>% 
  summarise(participants_multi = participants %>% unique() %>% sort() %>% paste(collapse = ","),
            participants_min = min(participants),
            participants_max = max(participants),
            participants = mean(participants)) %>% 
  ungroup()
hba1c_agg <- hba1c_agg %>% 
  left_join(hba1c_meta_ns_orig %>% select(nct_id, result_id, participants))

### Select outcomes from most to least informative ----
## Initially take arm-level mean change
hba1c_meta_mean <- hba1c_meta %>% 
  filter(result_type_smry %in% c("mean_change"))
hba1c_agg_mean <- hba1c_agg %>% 
  semi_join(hba1c_meta_mean %>% unnest(result_id))
hba1c_meta <- hba1c_meta %>% 
  filter(!nct_id %in% hba1c_meta_mean$nct_id)

## Next take endpoint-measure and if available start. Note that 
## Store the baseline as the value 1
hba1c_meta_end <- hba1c_meta %>% 
  filter(result_type_smry %in% c("mean_end"))
hba1c_meta_strt <- hba1c_meta %>% 
  filter(result_type_smry %in% c("mean_base"))
hba1c_agg_end <- hba1c_agg %>% 
  semi_join(hba1c_meta_end %>% unnest(result_id))
hba1c_agg_strt <- hba1c_agg %>% 
  semi_join(hba1c_meta_strt %>% unnest(result_id))
## Note will later need to impute and convert value_1
hba1c_agg_end <- hba1c_agg_end %>% 
  left_join(hba1c_agg_strt %>% select(arm_id_unq, value_1 = result, value_1_disp = dispersion, value_1_units = units_label))
hba1c_meta <- hba1c_meta %>% 
  filter(!nct_id %in% hba1c_meta_end$nct_id)
# Next take contrast in mean change
hba1c_meta_comp <- hba1c_meta %>% 
  filter(result_type_smry == "between_arm_mean") 
hba1c_agg_comp <- hba1c_agg %>% 
  semi_join(hba1c_meta_comp %>% unnest(result_id))
hba1c_meta <- hba1c_meta %>% 
  filter(!nct_id %in% hba1c_meta_comp$nct_id)

## leaves some trials with medians and percentage change. Will need to drop these
exclude <- tibble(trial_id =  hba1c_meta$nct_id %>% unique() ,
                  exclusion_reason2 = "only reported outcomes as medians or percentage change")
exclusions <- ExcludeRun(exclude = exclude)

saveRDS(list(arm  = list(data = hba1c_agg_mean,
                        meta = hba1c_meta_mean),
             comp = list(data = hba1c_agg_comp,
                         meta = hba1c_meta_comp),
             end  = list(data = hba1c_agg_end,
                         meta = hba1c_meta_end),
             unav = list(meta = hba1c_meta)),
        "Scratch_data/agg_hba1c.Rds")

