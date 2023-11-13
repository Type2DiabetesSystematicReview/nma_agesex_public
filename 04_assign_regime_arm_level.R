library(tidyverse)
## deal with combo trials
## sometimes these are dual or mono depending which arems we choose I think we should split the networks accordingly
source("Scripts/00_functions.R")

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
rm(combo, combo_smry)
combo <- read_csv("Created_metadata/simplify_combination_trials.csv")
whichnwork_simple <- whichnwork %>% 
  rename(nct_id = trial_id) %>% 
  anti_join(combo %>% select(nct_id))
whichnwork_multi <- whichnwork %>% 
  rename(nct_id = trial_id) %>% 
  semi_join(combo %>% select(nct_id))
arm_assign_simple <- arm_assign %>% 
  anti_join(combo %>% select(nct_id))


## add which network to arms with one drug per arm
arm_assign_simple <- arm_assign_simple %>% 
  inner_join(whichnwork_simple %>% select(nct_id, drug_regime))

## add which network to arms with more than one drug per arm. need to select dual and triple separately
dual <- combo %>% 
  filter(!is.na(dual)) %>% 
  select(nct_id, arm_id_unq, drug_code) 
## get rid of double underscore as can cause problems
dual <- dual %>% 
  mutate(drug_code = str_replace_all(drug_code, "__", "$$"))
## simplify dual and combo
dual$drug_codes <- map(dual$drug_code, ~ str_split(.x, "_") %>% unlist())
dual <- dual %>% 
  select(-drug_code) %>% 
  rename(drug_code = drug_codes) %>% 
  unnest(drug_code)
dual2 <- SimplifyDrugs(dual)$keep
dual_cmpr <- dual %>% 
  anti_join(dual2) %>% 
  inner_join(dual2 %>% rename(drug_code2 = drug_code))
rm(dual, dual_cmpr)
dual <- dual2 %>% 
  mutate(drug_code = str_replace_all(drug_code, fixed("$$"), "__"),
         drug_code = if_else(drug_code == "implicit_control", "placebo", drug_code)) %>% 
  select(-sameacross)
rm(dual2)

triple <- combo %>% 
  filter(!is.na(triple)) %>% 
  select(nct_id, arm_id_unq, drug_code) 
triple <- triple %>% 
  mutate(drug_code = str_replace_all(drug_code, "__", "$$"))
## simplify triple and combo
triple$drug_codes <- map(triple$drug_code, ~ str_split(.x, "_") %>% unlist())
triple <- triple %>% 
  select(-drug_code) %>% 
  rename(drug_code = drug_codes) %>% 
  unnest(drug_code)
triple2 <- SimplifyDrugs(triple)$keep
triple_cmpr <- triple %>% 
  anti_join(triple2) %>% 
  inner_join(triple2 %>% rename(drug_code2 = drug_code))
rm(triple, triple_cmpr)
triple <- triple2 %>% 
  mutate(drug_code = str_replace_all(drug_code, fixed("$$"), "__"),
         drug_code = if_else(drug_code == "implicit_control", "placebo", drug_code)) %>% 
  select(-sameacross)
rm(triple2)

arm_assign_multi <- bind_rows(dual %>% 
                                mutate(drug_regime = "dual"),
                              triple %>% 
                                mutate(drug_regime = "triple")) %>% 
  mutate(trtcls5  = str_sub(drug_code, 1, 5),
         trtcls4 = str_sub(drug_code, 1, 4))
arm_assign_final <- bind_rows(arm_assign_simple,
                              arm_assign_multi)
## simplify regime names
arm_assign_final <- arm_assign_final %>% 
  mutate(drug_regime_smpl = case_when(
  drug_regime %in% c("triple", "triple+", "dual|triple+", "mono|dual|triple+") ~ "triple",
  drug_regime %in% c("dual", "mono|dual") ~ "dual",
  drug_regime %in% "mono" ~ "mono",
  TRUE ~ "unclear or missing"))

arm_assign_final %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = length(arm_id_unq),
            dups = sum(duplicated(arm_id_unq)))
arm_assign_final %>% 
  group_by(drug_regime_smpl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = length(arm_id_unq),
            dups = sum(duplicated(arm_id_unq)))
rv_trls <- arm_assign_final %>% 
  distinct(nct_id, drug_regime_smpl) %>% 
  arrange(drug_regime_smpl) %>% 
  group_by(nct_id) %>% 
  summarise(trl_type = paste(drug_regime_smpl, collapse = ",")) %>% 
  ungroup() 
rv_trls %>% 
  count(trl_type)
## note that trials (and arms) can appear in more than one network
saveRDS(arm_assign_final, "Scratch_data/arms_assign_drug_regime.Rds")
