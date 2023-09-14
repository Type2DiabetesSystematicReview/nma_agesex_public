## agesex
library(tidyverse)
library(multinma)

## read in aggregate level data for each hba1c trials and baselien data
hba1c_agg1 <- read_csv("../cleaned_data/Data/hba1c_outcome_data_2019.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg2 <- read_csv("../cleaned_data/Data/hba1c_outcome_data_2022.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg <- bind_rows(extract1 = hba1c_agg1,
                       extract2 = hba1c_agg2,
                       .id = "extract")
rm(hba1c_agg1, hba1c_agg2)
base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
  rename(nct_id = trial_id)
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds") %>% 
  rename(nct_id = trial_id)
## pull participants data from baseline
participants <- base_dsp %>% 
  filter(variable == "n")  %>% 
  select(unq_id, nct_id, arm_id, id, participants = first)
hba1c_agg <- hba1c_agg %>% 
  left_join(participants %>% rename(arm_id_unq = arm_id))
rm(participants)

## Only 10 trials with missing n where have sd instead of se. Will need to pull from results analyse in aact
hba1c_agg %>% 
  filter(dispersion_type == "sd" & is.na(participants)) %>% 
  count(nct_id)
warning("Still to add trials missing N's to this analysis")

## separate out result metadata 
hba1c_meta <- hba1c_agg %>% 
  select(nct_id, extract, analysis_population, outcome_harm_label, result_description, result_type, 
                                     dispersion_type, timepoint_units, source, result_id) %>% 
  nest(result_id = result_id)
hba1c_agg <- hba1c_agg %>% 
  select(nct_id, result_id, arm_id_unq, comp_id, ancova, result, dispersion,
         units_label, timepoint, arm_id_subgroup, unq_id, id, participants) %>% 
  distinct()

hba1c_meta %>% 
  count(result_type)
result_determination <- read_csv("Created_metadata/resolve_result_type.csv")
hba1c_meta <- hba1c_meta %>% 
  inner_join(result_determination)
result_determination_rv <- hba1c_meta %>% 
  count(nct_id, result_type_smry) %>% 
  spread(result_type_smry, n, fill = 0L)
## 633 have mean change or mean base and end. 38 do not. drop for now
## codify highest level of data
result_type_per_trial <- hba1c_meta %>% 
  group_by(nct_id) %>% 
  summarise(result_type_best = case_when(
    any(result_type_smry == "mean_change") ~ "mean_change",
    any(result_type_smry == "mean_base") & any(result_type_smry == "mean_end") ~ "mean_change_calculate",
    any(result_type_smry == "between_arm_mean") ~ "between_arm_mean",
    TRUE ~ "other"
  ))
## 577 mean change, 65 mean change calculate, 15 between arm means, 23 "others"
result_type_per_trial %>% 
  count(result_type_best)
## for now limit to mean change. Can calculate later
hba1c_meta <- hba1c_meta %>% 
  filter(result_type_smry %in% c("mean_change", "mean_base", "mean_end"))

hba1c_agg <- hba1c_agg %>% 
  semi_join(hba1c_meta %>% unnest(result_id))
warning("Dropped 38 arms without some kind of mean HbA1c change, or baseline and end. includes dropping difference between arms")
## calculate standard errors
hba1c_dsp <- hba1c_agg %>% 
  inner_join(hba1c_meta %>% 
               filter(dispersion_type %in% c("se", "sd")) %>% 
               select(dispersion_type, result_id) %>% 
               unnest(result_id))
## Note. 10 trials without standard error wihtou participant number. Need to find
hba1c_dsp <- hba1c_dsp %>% 
  mutate(se = case_when(
    dispersion_type == "se" ~ as.double(dispersion),
    dispersion_type == "sd" & !is.na(participants) ~ as.double(dispersion)/participants^0.5,
    TRUE ~ NA_real_
  )) %>% 
  select(-dispersion, -dispersion_type) 

hba1c_rng <- hba1c_agg %>% 
  inner_join(hba1c_meta %>% 
               filter(dispersion_type %in% c("95%ci", "95% ci")) %>% 
               select(dispersion_type, result_id) %>% 
               unnest(result_id)) %>% 
  separate(dispersion, into = c("ll", "ul"), sep = ",|\\;")

warning("Need to fix trials with missing ll/ul")
## where no missing UL and LL seems fine.
## where is missing seems to be a pasting (or similar) error I need to correct
hba1c_rng %>% filter(is.na(ll) | is.na(ul))
# NCT00813995 looking at CTG result is correct but dispersion is wrong
# NCT00837577 looking at CTG result is correct but dispersion is wrong
# NCT00885352 ditto and for other in this set
hba1c_rng <- hba1c_rng %>% 
  filter(!is.na(ul), !is.na(ll)) %>% 
  mutate(across(c(ul, ll), as.double)) %>% 
  mutate(se = (ul-ll)/(2*1.96))
hba1c_agg2 <- bind_rows(dsp = hba1c_dsp,
                       rng = hba1c_rng %>%  select(-ll, -ul), .id = "dsp_rng")
## it is quite complete WRT standard errors
hba1c_agg2 %>% 
  group_by(dsp_rng) %>% 
  summarise(across(c(result, se), ~ mean(!is.na(.x)))) %>% 
  ungroup()

## Calculate change in hba1c where have baseline and end ----
hba1c_agg %>% 
  count(result_description)
hba1c_chng <- hba1c_agg %>% 
  filter(result_description %in% c("adjusted mean change from baseline", 
                                   "change from baseline"))
hba1c_agg <- hba1c_agg %>% 
  filter(!nct_id %in% hba1c_chng$nct_id)
hba1c_base <- hba1c_agg %>% 
  filter(result_description %in% c("baseline", "baseline adjusted mean"))
hba1c_base %>% 
  count(result_type, dispersion_type)
hba1c_end <- hba1c_agg %>% 
  filter(result_description %in% c("endline", "endpoint"))
hba1c_end %>% 
  count(result_type, dispersion_type)
hba1c_diff <- bind_rows(hba1c_base,
                         hba1c_end) 
## drop the remaining base and end ones from agg
hba1c_agg <- hba1c_agg %>% 
  anti_join(hba1c_diff %>% 
              select(nct_id, arm_id_unq)) 
hba1c_diff <- hba1c_diff %>% 
  group_by(nct_id, arm_id_unq) %>% 
  summarise(result = mean(result),
            se = sum(se^2)^0.5,
            result_id  = paste(result_id %>% unique(), collapse = ","),
            across(c(dsp_rng, extract, analysis_population, comp_id:ancova, units_label:participants, dispersion_type), ~
                     .x %>% unique() %>% paste(collapse = ","))) %>% 
  mutate(result_description = "David calculated change") %>% 
  ungroup() %>% 
  mutate(participants = as.double(participants))
## 3 trials without participant data
hba1c_agg <- bind_rows(hba1c_agg %>% mutate(timepoint = as.character(timepoint)),
                       hba1c_diff)

## Pull age data for model ----
base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
  rename(nct_id = trial_id)
## almost all mean and sd
warning("Need to process other (non mean and sd)")
base_dsp_age <- base_dsp %>% 
  filter(variable == "age", first_format == "mean", second_format == "sd") %>% 
  mutate(age_m = first,
         age_s = second)

