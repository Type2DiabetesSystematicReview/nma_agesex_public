#01a_process_hba1c_initial
library(tidyverse)

## read in all data ----
## read in ipd so can drop from aggregate
ipd <- read_csv("../from_vivli/Data/agesex/hba1c_base_change_overall.csv")
## read in aggregate level data for each hba1c trials and baseline data
hba1c_agg1 <- read_csv("../cleaned_data/Data/hba1c_outcome_data_2019.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg2 <- read_csv("../cleaned_data/Data/hba1c_outcome_data_2022.csv") %>% 
  rename(nct_id = trial_id)
arm <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv")

## read in outcome counts
oc_orig <- readRDS("../extract_transform/aact/data/aact_extract.Rds")$outcome_counts
oc_new <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$outcome_counts 
## read in IDs for mapping counts onto outcomes
## the first file is all hba1c
hba1c_ids_orig <- read_csv("../t2dm_nma_data_extract/Created_metadata/hba1c_continuous_aact.csv") %>% 
  distinct(nct_id, id, outcome_id, result_group_id, ctgov_group_code)
## note lacking IDs. Need to add from aact
hba1c_ids_new <- read_csv("../t2dm_nma_data_extract/aact/2022_update/Created_metadata/unique_outcomes_for_harmonisation_RVd_EB.csv") %>%
  filter(variable == "hba1c") 
hba1c_ids_new2 <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$outcome_measurements
hba1c_ids_new <- hba1c_ids_new %>% 
  inner_join(hba1c_ids_new2 %>% 
              select(nct_id, title, description, outcome_id, units, id, result_group_id, ctgov_group_code) %>% 
               distinct())
rm(hba1c_ids_new2)

# assign each type of outcome to a wider category. Eg mean change, baseline etc
result_determination <- read_csv("Created_metadata/resolve_result_type.csv")

## initial cleaning ----
hba1c_agg <- bind_rows(extract1 = hba1c_agg1,
                       extract2 = hba1c_agg2,
                       .id = "extract")
rm(hba1c_agg1, hba1c_agg2)
# drop ipd
hba1c_agg <- hba1c_agg %>% 
  filter(!nct_id %in% ipd$nct_id)
# drop if not a selected arm
arm <- arm %>% 
  filter(!arm_label == "total") %>% 
  select(trial_id, arm_id_unq, arm_id_subgroup, arm_id) %>% 
  gather("armids", "value", -trial_id) %>% 
  select(-armids) %>% 
  distinct()
hba1c_agg <- hba1c_agg %>% 
  semi_join(arm %>% select(nct_id = trial_id,
                           arm_id_unq = value))

## correct errors
hba1c_agg <- hba1c_agg %>% 
  mutate(result_type = if_else(nct_id == "NCT00509262" & result_description == "change from baseline", "mean", result_type))
# drop trials not in 
hba1c_agg %>% 
  count(nct_id)

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
## 633 have mean change or mean base and end. 38 do not. 
result_type_per_trial <- hba1c_meta %>% 
  group_by(nct_id) %>% 
  summarise(result_type_best = case_when(
    any(result_type_smry == "mean_change") ~ "mean_change",
    any(result_type_smry == "mean_base") & any(result_type_smry == "mean_end") ~ "mean_change_calculate",
    any(result_type_smry == "between_arm_mean") ~ "between_arm_mean",
    TRUE ~ "other"
  ))
## 524 mean change, 56 mean change calculate, 15 between arm means, 23 "others"
result_type_per_trial %>% 
  count(result_type_best)

### Add in number of participants where available based on result_id ----
## pull other missing ns from results  note the IDs are extract specific
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

# Next take contrast in mean change
hba1c_meta_comp <- hba1c_meta %>% 
  filter(result_type_smry == "between_arm_mean") 
hba1c_agg_comp <- hba1c_agg %>% 
  semi_join(hba1c_meta_comp %>% unnest(result_id))
hba1c_meta <- hba1c_meta %>% 
  filter(!nct_id %in% hba1c_meta_comp$nct_id)

## Next take endpoint-measure
hba1c_meta_end <- hba1c_meta %>% 
  filter(result_type_smry %in% c("mean_end"))
hba1c_agg_end <- hba1c_agg %>% 
  semi_join(hba1c_meta_end %>% unnest(result_id))
hba1c_meta <- hba1c_meta %>% 
  filter(!nct_id %in% hba1c_meta_end$nct_id)

## leaves 17 trials with medians and percentage change. Will need to drop these

saveRDS(list(arm  = list(data = hba1c_agg_mean,
                        meta = hba1c_meta_mean),
             comp = list(data = hba1c_agg_comp,
                         meta = hba1c_meta_comp),
             end  = list(data = hba1c_agg_end,
                         meta = hba1c_meta_end),
             unav = list(meta = hba1c_meta)),
        "Scratch_data/agg_hba1c.Rds")

