library(tidyverse)
## deal with combo trials
## sometimes these are dual or mono depending which arems we choose I think we should split the networks accordingly

arm_assign <- read_csv("Data/arm_labels_hba1c.csv")
whichnwork <- read_csv("../cleaned_data/Data/ancillary_drugs_data_all_cleaned.csv")
arm_meta_orig <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- read_csv("Data/agg.csv")

## re-run ancillary algorithm after dropping combination arms
## drops 2 trials only
combo <- arm_assign %>% 
  filter(str_detect(trtcls5, "_")) %>% 
  distinct(nct_id)
combo <- arm_assign %>% 
  semi_join(combo)
## assign format , placebo, drug class as A, B, A+B
combo$lng <- map(combo$drug_code, ~ str_split(.x, patt = "_") %>% unlist)
combo <- combo %>% 
  unnest(lng)
combo <- combo %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  ungroup()
combo$data <- map(combo$data, ~ .x %>% 
                    arrange(lng) %>% 
                    mutate(ltr = case_when(
                             lng == "placebo" ~ "P", 
                             lng == "A10BA02" ~ "M",
                             TRUE ~ LETTERS[cumsum(!duplicated(lng))])) %>% 
                    arrange(arm_id_unq, ltr) %>% 
                    group_by(arm_id_unq, drug_code) %>% 
                    summarise(ltr = paste(unique(ltr), collapse = "")) %>% 
                    ungroup())
combo <- combo %>% 
  unnest(data)
combo_smry <- combo %>% 
  arrange(nct_id, ltr) %>% 
  group_by(nct_id) %>% 
  summarise(cmpr = paste(ltr, collapse = ",")) %>% 
  ungroup()
combo <- combo %>% 
  left_join(combo_smry)
combo <- combo %>% 
  left_join(whichnwork %>% select(nct_id = trial_id, drug_regime, mandatory_ancillary1))
ns <- bind_rows(agg = agg %>% 
                  select(nct_id, arm_id_unq, n),
                ipd = ipd %>% select(nct_id, arm_id_unq) %>% 
                  count(nct_id, arm_id_unq), .id = "data_type")
combo <- combo %>% 
  left_join(ns)
combo <- combo %>% 
  # select(-ltr) %>% 
  left_join(arm_meta_orig %>% select(nct_id, arm_id_unq, drug_name:drug2_freq) %>% distinct())
combo <- combo %>% 
  count(cmpr, sort = TRUE) %>% 
  select(-n) %>% 
  left_join(combo)
# write_csv(combo, "Created_metadata/to_simplify_combination_trials.csv")
## note renamed to "simplify_combination_trials.csv" so do not accidentally overwrite

