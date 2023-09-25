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
# 671 trials with hba1c measures
hba1c_agg %>% 
  count(nct_id)
rm(hba1c_agg1, hba1c_agg2)
base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
  rename(nct_id = trial_id)
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds") %>% 
  rename(nct_id = trial_id)
## pull participants data from baseline
participants <- read_csv("Data/ns.csv") %>% 
  rename(nct_id = trial_id)
hba1c_agg <- hba1c_agg %>% 
  left_join(participants %>% rename(arm_id_unq = arm_id, participants = n))
rm(participants)
## 10 trials with missing n where have sd instead of se. 5 are in clinicaltrials.gov 3 of which are in sex database (obtained by multiplying n by %)
aact <- readRDS("../extract_transform/aact/data/aact_extract.Rds")
no_n <- hba1c_agg %>% 
  filter(dispersion_type == "sd" & is.na(participants)) %>% 
  distinct(nct_id, arm_id_unq) 
sex <- read_csv("Data/sex.csv") %>% 
  rename(nct_id = trial_id, arm_id_unq2 = arm_id)
no_n2 <- no_n %>% 
  left_join(sex)
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
  filter(result_type_smry %in% c("mean_change"))
hba1c_agg <- hba1c_agg %>% 
  semi_join(hba1c_meta %>% unnest(result_id))
warning("Dropped arms without some mean HbA1c change")

## calculate standard errors
hba1c_dsp <- hba1c_agg %>% 
  inner_join(hba1c_meta %>% 
               filter(dispersion_type %in% c("se", "sd")) %>% 
               select(dispersion_type, result_id) %>% 
               unnest(result_id))
## Note. 10 trials without standard error without participant number. Need to find
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

hba1c_agg <- bind_rows(dsp = hba1c_dsp,
                       rng = hba1c_rng %>%  select(-ll, -ul, -dispersion_type), .id = "dsp_rng")

## Pull age data for model ----
## almost all mean and sd
warning("Need to process other (non mean and sd)")
age <- base_dsp %>% 
  filter(variable == "age", first_format == "mean", second_format == "sd") %>% 
  mutate(age_m = first,
         age_s = second) %>% 
  select(nct_id, id_source, arm_id_unq = arm_id, age_m, age_s)
setdiff(union(age$arm_id_unq, hba1c_agg$arm_id_unq), age$arm_id_unq)
setdiff(union(age$arm_id_unq, hba1c_agg$arm_id_unq), hba1c_agg$arm_id_unq)
age <- age %>% 
  mutate(inhba1c = arm_id_unq %in% hba1c_agg$arm_id_unq)
smry1 <- age %>% 
  group_by(id_source, nct_id) %>% 
  summarise(inhba1c = if_else(any(inhba1c), "got", "not")) %>% 
  ungroup() %>% 
  count(id_source, inhba1c) %>% 
  spread(inhba1c, n, fill = 0L) 
bind_rows(smry1,
         smry1 %>% 
           mutate(id_source = "total") %>% 
           group_by(id_source) %>% 
           summarise(across(c(got, not), sum)) %>% 
           ungroup())
hba1c_agg %>% 
  count(nct_id)
# 543 trials with hba1c in current set; so 36 of these mismatching
# 40 trials where age is not merged in
# 18 are mean and SD. checked one has subgroup labels on arms
# some are age in other formats, eg mean and range or median and IQR 
# 9 are definitively without any age data 
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
ipd <- ipd %>% 
  distinct(nct_id)

hba1c_agg_noage <- hba1c_agg %>% 
  filter(!arm_id_unq %in% age$arm_id_unq)
hba1c_meta_noage <- hba1c_meta %>% 
  unnest(result_id) %>% 
  semi_join(hba1c_agg_noage %>% select(result_id)) %>% 
  nest(data = c(result_id)) %>% 
  left_join(base_dsp %>% 
              filter(variable == "age") %>% 
              select(nct_id, first_format, second_format) %>% 
              distinct(nct_id, .keep_all = TRUE)) %>% 
  left_join(base_rng %>% 
              filter(variable == "age") %>% 
              select(nct_id, first_format_rng = first_format, second_format_rng = second_format) %>% 
              distinct(nct_id, .keep_all = TRUE))
ipdage <- intersect(hba1c_meta_noage$nct_id, ipd$nct_id)
## 4 trials without age have IPD. Can ignore. 36 trials need to resolve age data
hba1c_meta_noage <- hba1c_agg_noage %>% 
  filter(!nct_id %in% ipdage)
ageother <- base_dsp %>%
  filter(nct_id %in% hba1c_meta_noage$nct_id) %>% 
  filter(variable == "age")


## join what is already matching for purpose of running model
hba1c_agg <- hba1c_agg %>% 
  inner_join(age)
# 513 trials
hba1c_agg <- hba1c_agg %>% 
  select(nct_id, arm_id_unq, participants, result, se, age_m, age_s) %>% 
  distinct()
dups <- hba1c_agg %>% 
  select(nct_id, arm_id_unq) %>% 
  duplicated()
dups <-  hba1c_agg %>% 
  select(nct_id, arm_id_unq) %>% 
  filter(dups) %>% 
  distinct()
## 14 rows with duplicates. Will need to resolve. Appears that issue is presence of subgroups. for now. drop
## need to deal with 
dups <- hba1c_agg %>% 
  semi_join(dups)
hba1c_agg <- hba1c_agg %>% 
  filter(!arm_id_unq %in% dups$arm_id_unq)
saveRDS(hba1c_agg, "Scratch_data/agg_hba1c.Rds")
