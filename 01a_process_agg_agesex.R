library(tidyverse)

source("../common_functions/Scripts/combine_sd.R")
rm(a, ab, b, mean_res, means, mymeans, myns, mysds, sd_res, sds)
source("../common_functions/Scripts/convert_iqr_to_sd.R")

base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds")
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds")

## read in arm data to allow to drop subgroups (to allow aggregating over these)
arm <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv")

## replace subgroup arm id with arm id - create look-up vector
arm_lkp <- arm %>% 
  filter(!is.na(arm_id_subgroup)) %>% 
  distinct(arm_id_unq, arm_id_subgroup, subgroup_name, subgroup_name2)
arm_lkp_vct1 <- arm_lkp$arm_id_unq
names(arm_lkp_vct1) <- arm_lkp$arm_id_subgroup
arm_lkp_vct2 <- setdiff(c(base_dsp$arm_id, base_rng$arm_id), names(arm_lkp_vct1))
names(arm_lkp_vct2) <- arm_lkp_vct2
arm_lkp_vct <- c(arm_lkp_vct1, arm_lkp_vct2)
sum(duplicated(names(arm_lkp_vct)))
rm(arm, arm_lkp_vct1, arm_lkp_vct2)
write_csv(arm_lkp, "Data/arm2_armsg.csv")

## pull number in each arm for linking
ns <- base_dsp %>% 
  filter(variable == "n") %>% 
  select(trial_id, arm_id, n = first)

## Pull age as continuous variable
age <- base_dsp %>% 
  filter(variable == "age") 
## 4 trials with age data but no number of participants
age %>% 
  anti_join(ns)
age <- age %>% 
  left_join(ns)
age %>% 
  count(first_format, second_format)
# Age mean and standard deviation
age_ms <- age %>% 
  filter(first_format == "mean", second_format == "sd") %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(age_sd = CombSdVectorised(n = n, m = first, s = second),
            age_m = weighted.mean(x = first, w = n)) %>% 
  ungroup()
age <- age %>% 
  anti_join(age_ms)
# age mean without SD
age_m <- age %>% 
  filter(first_format == "mean") %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(age_m = weighted.mean(x = first, w = n)) %>% 
  ungroup()
## one trial with median only (no range or anything)
age_med <- age %>% 
  anti_join(bind_rows(age_ms, age_m) %>% select(trial_id))
rm(age)

## age as a categorical variable, derive mean and sd 
age_ms_ctg <- base_rng %>% 
  filter(variable == "age") %>% 
  left_join(ns) %>% 
  mutate(
    age_sd = 
      case_when(
        second_format == "iqr" ~ EstSD2(low, upp, n),
        second_format == "range" ~ EstSD(low, upp, n)),
    age_m = case_when(
      first_format == "mean" ~ first,
      first_format == "median" & second_format == "range" ~ EstMean(first, low, upp),
      first_format == "median" & second_format == "iqr" ~ EstMean2(low, first, upp)
    )) %>% 
  select(trial_id, age_m, age_sd, arm_id, n) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(age_sd = CombSdVectorised(n = n, m = age_m, s = age_sd),
            age_m = weighted.mean(x = age_m, w = n)) %>% 
  ungroup()

age <- bind_rows(age_m,
                 age_ms,
                 age_ms_ctg)
## 14 trials with no age data
bind_rows(base_dsp, base_rng) %>% 
  select(trial_id) %>% 
  anti_join(age) %>% 
  distinct(trial_id)

## pull data on sex
sex <- base_dsp %>% 
  filter(variable == "male")
## all with missing number in each arm have sex as both a count and a percentage so calculate n
## then caluclate percentage if not already present and count if not already present
sex %>% 
  anti_join(ns) %>% 
  count(first_format, second_format)
sex <- sex %>% 
  left_join(ns)
sex <- sex %>% 
  mutate(ns = if_else(is.na(n) & second_format == "percentage", (second/100)/first, n))
sex <- sex %>% 
  mutate(male = if_else(is.na(first), (second/100)*n, first)) %>% 
  select(trial_id, arm_id, male, n) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(male = sum(male),
            n = sum(n)) %>% 
  ungroup()

## race_ethnic
race <- base_dsp %>% 
  filter(variable %in% c("race", "ethnicity"))
## as with sex, all ones without Ns have count and %
race %>% 
  anti_join(ns) %>% 
  count(first_format, second_format)
race <- race %>% 
  left_join(ns) 
race <- race %>% 
  mutate(first = if_else(is.na(first), (second/100)*n, first)) %>% 
  select(trial_id, arm_id, variable, category_level, n, first) %>% 
  group_by(trial_id, arm_id, variable, category_level) %>% 
  summarise(x = sum(first),
            n = sum(n)) %>% 
  ungroup()

## aggregate over subgroups ----
## note important that do this for ns last as used to weight
age <- age %>% 
  left_join(ns) %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  ## give equal weighting if n is missing
  mutate(n = if_else(is.na(n), 1L, n)) %>% 
  summarise(age_sd = CombSdVectorised(m = age_m, s = age_sd, n = n),
            age_m = weighted.mean(age_m, n)) %>% 
  ungroup()
write_csv(age, "Data/age.csv")

sex <- sex %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(n = sum(n),
            male = sum(male)) %>% 
  ungroup()
write_csv(sex, "Data/sex.csv")

race <- race %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id, variable, category_level) %>% 
  summarise(n = sum(n),
            x = sum(x)) %>% 
  ungroup()
write_csv(race, "Data/race_ethnicity.csv")

ns <- ns %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(n = sum(n)) %>% 
  ungroup()
write_csv(ns, "Data/ns.csv")


mydfs <- list.files("Data/")
res <- map(mydfs, ~ paste0("Data/", .x))
res <- map(res, read_csv)
names(res) <- list.files("Data/") %>% str_sub(1, -5)
