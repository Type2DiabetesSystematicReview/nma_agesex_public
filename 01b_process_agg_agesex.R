library(tidyverse)

source("../common_functions/Scripts/combine_sd.R")
# rm(a, ab, b, mean_res, means, mymeans, myns, mysds, sd_res, sds)
source("../common_functions/Scripts/convert_iqr_to_sd.R")

base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds")
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds")
hba1c_agg <- readRDS("Scratch_data/agg_hba1c.Rds")
comp2arm <- read_csv("../cleaned_data/Data/example_comparisons.csv")

mnl_n <- read_csv("trial_id,drug_name,n,arm_id
100-IRMI/PRI 16/6/2 (007/2017),DAPA,36,updac0091
100-IRMI/PRI 16/6/2 (007/2017),Placebo,36,updac0092
Hawler Medical University records of the clinical trials: No.276,Group A Glimepride,26,updac0231
Hawler Medical University records of the clinical trials: No.276,Group B Sitagliptin,28,updac0230
Hawler Medical University records of the clinical trials: No.276,Group C Canagliflozin,24,updac0229
NCT02730377,OAD,995,uaa11213
NCT02730377,liraglutide,996,uaa11212
Netherlands Trial Register NTR6709,Alirocumab,6,updac0150
Netherlands Trial Register NTR6709,Placebo,6,updac0151
UMIN000022953,Linagliptin group,21,ac0198
UMIN000022953,Metformin group,22,ac0199")

## create arm-level lookup based on comparison data and arm-level data
arm_in1 <- bind_rows(hba1c_agg$arm$data,
                  hba1c_agg$end$data) %>% 
  distinct(nct_id, arm_id_unq)
warning("6 trials without comp_id to arm_id lookup")
arm_in2 <- hba1c_agg$comp$data %>% 
  select(nct_id, comp_id) %>% 
  inner_join(comp2arm %>% select(comp_id, nct_id = trial_id, arm_id_unq = arm_id)) %>% 
  select(nct_id, arm_id_unq)
arm_in <- bind_rows(arm_in1, arm_in2) %>% 
  rename(trial_id = nct_id, arm_id = arm_id_unq)
rm(arm_in1, arm_in2)

## read in arm data to allow to drop subgroups (to allow aggregating over these)
arm <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv")
arm_ctg <- arm %>% 
  filter(!is.na(ctgov_group_code)) %>% 
  rename(nct_id = trial_id) %>% 
  distinct(nct_id, arm_id_unq, arm_id_subgroup, ctgov_group_code)
## replace subgroup arm id with arm id - create look-up vector
arm_lkp <- arm %>% 
  filter(!is.na(arm_id_subgroup)) %>% 
  distinct(arm_id_unq, arm_id_subgroup, subgroup_name, subgroup_name2)
arm_lkp_vct1 <- arm_lkp$arm_id_unq
names(arm_lkp_vct1) <- arm_lkp$arm_id_subgroup
arm_lkp_vct2 <- setdiff(c(base_dsp$arm_id, base_rng$arm_id, arm_in$arm_id), names(arm_lkp_vct1))
names(arm_lkp_vct2) <- arm_lkp_vct2
arm_lkp_vct <- c(arm_lkp_vct1, arm_lkp_vct2)
sum(duplicated(names(arm_lkp_vct)))
sum(is.na(arm_lkp_vct))
rm(arm, arm_lkp_vct1, arm_lkp_vct2)
write_csv(arm_lkp, "Data/arm2_armsg.csv")

## pull number in each arm for linking
## use number linked to outcome as first preference, note may not have subgroups
ns <- base_dsp %>% 
  filter(variable == "n") %>% 
  select(trial_id, arm_id, n = first)
ns2 <- hba1c_agg$arm$data %>% 
  filter(!is.na(participants)) %>% 
  select(trial_id = nct_id, arm_id = arm_id_unq, n = participants)
ns <- bind_rows(ns2,
                ns,
                mnl_n %>% select(-drug_name)) %>% 
  distinct(trial_id, arm_id, .keep_all = TRUE)
rm(ns2, mnl_n)

## Pull age as continuous variable
age <- base_dsp %>% 
  filter(variable == "age") 
age <- age %>% 
  left_join(ns)

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

## pull data on sex
sex <- base_dsp %>% 
  filter(variable == "male")
## all with missing number in each arm have sex as both a count and a percentage so calculate n
## then calculate percentage if not already present and count if not already present
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
basedata <- list(
  age = age %>% 
  left_join(ns) %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  ## give equal weighting if n is missing
  mutate(n = if_else(is.na(n), 1L, n)) %>% 
  summarise(age_sd = CombSdVectorised(m = age_m, s = age_sd, n = n),
            age_m = weighted.mean(age_m, n)) %>% 
  ungroup(),
 sex =  sex %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(n = sum(n),
            male = sum(male)) %>% 
  ungroup(),
 race = race %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id, variable, category_level) %>% 
  summarise(n = sum(n),
            x = sum(x)) %>% 
  ungroup(), 
 ns_sg = ns,
 ns = ns %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  group_by(trial_id, arm_id) %>% 
  summarise(n = sum(n)) %>% 
  ungroup())
rm(age, sex, race, ns)

## convert arm_in to arm_id_unq
arm_in <- arm_in %>% 
  mutate(arm_id = arm_lkp_vct[arm_id]) %>% 
  distinct(trial_id, arm_id)
  
basedata2 <- map(basedata, ~ .x %>% 
                   semi_join(arm_in))
map_int(basedata, nrow)
map_int(basedata2, nrow)


saveRDS(basedata2, "Scratch_data/agg_hba1c_base.Rds")
## Addressing missing ones n's where are using se's ----

## data complete wherever we have hba1c data
basedata2$age %>% 
  anti_join(arm_in)
basedata2$sex %>% 
  anti_join(arm_in)
basedata2$race %>% 
  anti_join(arm_in)
basedata2$ns %>% 
  anti_join(arm_in)


