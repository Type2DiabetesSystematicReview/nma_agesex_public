library(tidyverse)

source("../common_functions/Scripts/combine_sd.R")
rm(a, ab, b, mean_res, means, mymeans, myns, mysds, sd_res, sds)
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
arm_lkp_vct2 <- setdiff(c(base_dsp$arm_id, base_rng$arm_id), names(arm_lkp_vct1))
names(arm_lkp_vct2) <- arm_lkp_vct2
arm_lkp_vct <- c(arm_lkp_vct1, arm_lkp_vct2)
sum(duplicated(names(arm_lkp_vct)))
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

## pull participants data from baseline with subgroups for merging across ----
participants <- read_csv("Data/ns_by_sg.csv") %>% 
  rename(nct_id = trial_id,
         arm_id_unq = arm_id,
         participants = n)

## reviewed missing ones manually ----
function() {
  ## 5 trials without N data where have SD rather than SE. Following code helps extract and map these across
  se_present <- hba1c_agg %>% 
    filter(dispersion_type %in% c("se", "95% ci", "95%ci", "ci")) %>% 
    distinct(nct_id)
  mnl <- read_delim("nct_id;arm1;arm2;arm3
  100-IRMI/PRI 16/6/2 (007/2017);DAPA (n = 36);Placebo (n = 36);
  Hawler Medical University records of the clinical trials: No.276;Group A Glimepride (n = 26); Group B Sitagliptin (n = 28); Group C Canagliflozin (n = 24)
  NCT02730377;liraglutide, (n = 996); OAD, (n = 995);
  Netherlands Trial Register NTR6709;Alirocumab (n = 6); Placebo (n = 6);
  UMIN000022953;Linagliptin group (n = 21);Metformin group (n = 22);", delim = ";")
  mnl2 <- mnl %>%
    gather("arm", "value", -nct_id, na.rm = TRUE) %>%
    mutate(value = str_trim(value)) %>%
    separate(value, into = c("drug_name", "n"), sep = "n = ")  %>%
    mutate(drug_name = str_remove(drug_name, "\\(") %>%
             str_trim(),
           n = str_remove(n, "\\)") %>%
             str_trim() %>%
             as.integer()) %>%
    arrange(nct_id, drug_name) %>%
    select(nct_id, drug_name, n)
  ## none lost by merge
  addmnl <- hba1c_agg %>%
    semi_join(mnl2) %>%
    distinct(nct_id, arm_id_unq)
  armmnl <- arm %>%
    rename(nct_id = trial_id) %>%
    semi_join(addmnl) %>%
    select(nct_id, arm_id_unq, arm_id_subgroup, arm_label) %>%
    arrange(nct_id, arm_id_unq)
  write_csv(mnl2, "Scratch_data/manual_ns_missing.csv")
}


# mnln2 <- 

## review arm mean ----
## drop ones with se and ones with existing participant data
xmn <- hba1c_agg$arm$meta %>% 
  unnest(result_id) %>% 
  inner_join(hba1c_agg$arm$data) %>% 
  filter(is.na(participants)) 
se_present <- xmn %>%
  filter(dispersion_type %in% c("se", "95% ci", "95%ci", "ci")) %>%
  distinct(nct_id) %>% 
  pull()

# xmn <- xmn %>% 
#   filter(!dispersion_type %in% )
xmn <- xmn %>% 
  select(-participants) %>% 
  anti_join(participants)
# arm <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv")

hba1c_agg <- hba1c_agg %>% 
  mutate(arm_id_lnk = if_else(!is.na(arm_id_subgroup), arm_id_subgroup, arm_id_unq)) %>% 
  left_join(participants %>% rename(arm_id_lnk = arm_id_unq))
rm(participants)
## 10 trials with missing n where have sd instead of se. 5 are in clinicaltrials.gov 3 of which are in sex database (obtained by multiplying n by %)

no_n <- hba1c_agg %>% 
  filter(!nct_id %in% se_present, dispersion_type == "sd" & is.na(participants))


warning("Still to add trials missing N's to this analysis")


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

