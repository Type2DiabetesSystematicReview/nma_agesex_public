library(tidyverse)
library(multinma)

source("Scripts/00_functions.R")
source("Scripts/common_functions/Scripts/misc.R")
## read data ----
exclusions <- read_csv("Data/exclusions_update.csv")
## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in reg
reg <- readRDS("Scratch_data/ipd_coefs_frmttd.Rds")
## read in agg
agg <- read_csv("Scratch_data/agg.csv")
# read in arm metadata
whichnwork <- read_csv("Data/cleaned_data/Data/ancillary_drugs_data_all_cleaned.csv")
# as uploaded to vivli
arm_meta_orig <- read_csv("Data/cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)

arm_meta_orig <- arm_meta_orig %>% 
  semi_join(bind_rows(agg %>% select(nct_id, arm_id_unq),
                      ipd %>% select(nct_id, arm_id_unq)))
arm_meta_orig <- arm_meta_orig %>% 
  mutate(drug_name = case_when(
    arm_id_unq == "updac0239" ~ "sitagliptin|linagliptin|alogliptin|vildagliptin|saxagliptin |teneligliptin",
    TRUE ~ drug_name
  ))

# additions to trial metadata made within vivli ----
arm_meta_new <- read_csv("Data/agesexhba1c_6115/reference_arm_data_all_cleaned.csv")
arm_meta_new_ipd <- arm_meta_new %>% 
  semi_join(ipd %>% select(nct_id))
arm_meta_new_ipd <- arm_meta_new_ipd %>% 
  filter(!nct_id == "NCT01778049")
rm(arm_meta_new)
arm_meta_orig_ipd <- arm_meta_orig %>% 
  semi_join(ipd %>% select(nct_id))
currnames <- c("nct_id", "arm_id_unq", "drug_name", "drug_dose", "drug_unit", 
"drug_freq", "drug2_name", "drug2_dose", "drug2_unit", "drug2_freq", 
"new_note")
arm_meta_new_cmpr <- bind_rows(new = arm_meta_new_ipd[,currnames],
                               orig = arm_meta_orig_ipd %>% select(contains(currnames)), .id = "vivli_local")
arm_meta_new_cmpr <- arm_meta_new_cmpr %>% 
  gather("var", "value", -vivli_local, -nct_id, -arm_id_unq, na.rm = TRUE) %>% 
  pivot_wider(names_from = vivli_local, values_from = value)
arm_meta_new_cmpr$identical <- map2_lgl(arm_meta_new_cmpr$new, arm_meta_new_cmpr$orig, identical)
arm_meta_new_cmpr <- arm_meta_new_cmpr %>% 
  filter(!identical)
## all different ones are present in new but not old. No changes or subtractions. Additions only
arm_meta_new_cmpr$new <- map_chr(arm_meta_new_cmpr$new, ~ unique(.x))
arm_meta_new_cmpr <- arm_meta_new_cmpr %>% 
  select(-orig, -identical) %>% 
  spread(var, new)
# note that one drug has an unspecified dose but will be dropped anyway as has multiple drugs
arm_meta_new_cmpr <- arm_meta_new_cmpr %>% 
  mutate(drug_name = case_when(
                      drug_name == "linagliptinOL" ~ "linagliptin",
                      TRUE ~ drug_name))
arm_meta <- bind_rows(arm_meta_orig,
                      arm_meta_new_cmpr)
rm(arm_meta_new_cmpr, arm_meta_new_ipd, arm_meta_orig, arm_meta_orig_ipd)

## Drop combination arms ----
combo_drop <- arm_meta %>% 
  filter(!is.na(drug2_name)) %>% 
  distinct(nct_id, arm_id_unq) %>% 
  group_by(nct_id) %>% 
  summarise(dropped_arms = paste(arm_id_unq, collapse = ",")) %>% 
  ungroup()
## drop COMBINATION arms (this also directly drops 4 trials with ONLY combination arms)
arm_meta <- arm_meta %>% 
  filter(is.na(drug2_name)) 
combo_drop <- combo_drop %>% 
  left_join(arm_meta %>% 
              distinct(nct_id, arm_id_unq) %>% 
              count(nct_id) %>% 
              rename(other_arms = n)) %>% 
  mutate(other_arms = if_else(is.na(other_arms), 0L, other_arms))
## drops only 17 trials as 40 have two or more other arms
write_csv(combo_drop, "Outputs/Trials_arm_dropped_combination.csv")
exclude <- tibble(trial_id = combo_drop$nct_id[combo_drop$other_arms <2] %>% unique(),
                  exclusion_reason2 = "combination treatments only")
exclusions <- ExcludeRun(exclude = exclude, saveexclusion = TRUE)

## drop (17) TRIALS with one or fewer arms left after dropping combinations - already dropped arms
arm_meta <- arm_meta %>% 
  filter(!nct_id %in% combo_drop$nct_id[combo_drop$other_arms <2])
## Drop dietary arms from one trial. Trial not dropped
sum(!duplicated(arm_meta$nct_id))
write_csv(arm_meta %>% 
            filter(nct_id == "NCT02798744") %>% 
            select(nct_id, arm_id_unq, arm_label) %>% 
            filter(arm_label == "ignore_diet") %>% 
            group_by(nct_id) %>% 
            summarise(dropped_arms = paste(arm_id_unq, collapse = ",")) %>% 
            mutate(other_arms = 2), "Outputs/Trials_arm_dropped_dietary.csv")
arm_meta <- arm_meta %>% 
  filter(is.na(arm_label) | !arm_label == "ignore_diet")

## drop where an IPD trial is against an open label extension only
arm_meta <- arm_meta %>% 
  filter(!nct_id == "NCT00306384")
exclude <- tibble(trial_id = "NCT00306384",
                  exclusion_reason2 = "wrong phase/study design")
exclusions2 <- ExcludeRun(exclude, FALSE)
exclusions2 <- exclusions2 %>% 
  mutate(exclusion_note = if_else(trial_id == "NCT00306384",
                        "Open label extension study data only",
                        exclusion_note))
write_csv(exclusions2, "Data/exclusions_update.csv")
exclusions <- exclusions2
rm(exclusions2)

## Add placebo as a drug name when it is only listed as an arm label ----
arm_meta <- arm_meta %>% 
  mutate(drug_name = case_when(
    arm_label %in% c("control", "placebo") ~ "placebo",
    TRUE ~ drug_name))
## add arm label where it is missing
arm_meta <- arm_meta %>% 
  mutate(arm_label = if_else(is.na(arm_label), drug_name, arm_label))

## relabel drug names with ATC codes ----
rename_drug <- readxl::read_excel("Created_metadata/non_atc_drugs.xlsx", sheet = 1)
arm_meta <- arm_meta %>% 
  select(nct_id, arm_id_unq, arm_id_subgroup, arm_id, 
        drug_name, drug_dose, drug_unit, 
         drug_freq, note)
arm_meta <- arm_meta %>% 
  left_join(rename_drug) %>% 
  mutate(drug_name = if_else(!is.na(rename), rename, drug_name)) %>% 
  select(-rename)

## new arm labels in WHO database and own manual lookup
whoatc <- read_csv("Data/whoatcdiabetesnodose.csv")
whoatc <- whoatc %>% 
  select(nm = `ATC level name`,
         atc_code = `ATC code`) %>% 
  mutate(nm = str_to_lower(nm)) %>% 
  filter(str_sub(atc_code, 1, 3) == "A10") %>% 
  distinct()
nonwho <- readxl::read_excel("Created_metadata/non_atc_drugs.xlsx", sheet = 2)
nonwho <- nonwho %>% 
  mutate(nm = str_to_lower(drug_name)) %>% 
  select(nm, atc_code = code) %>% 
  distinct()
## where alternative (ie not multiple) drugs allowed but within the same class collapse these into the class name
alternatives <- read_csv(
  "nm,atc_code
pioglitazone|rosiglitazone,A10BG
liraglutide|exenatide|albiglutide|dulaglutide,A10BJ
aspart| lispro|,A10A
liraglutide| lixisenatide | glulisine,A10BJ
canagliflozin| dapagliflozin |empagliflozin,A10BK
mitiglinide|nateglinide,A10BX
voglibose|miglitol,A10BF
sitagliptin|tenegliptin|linagliptin|alogliptin,A10BH
glimperide| gliclazide,A10BB
metformin|glinide|alpha glucosidase inh,A10B
sitagliptin|linagliptin|alogliptin|vildagliptin|saxagliptin |teneligliptin,A10BH
sitagliptin|vildagliptin|alogliptin|linagliptin|teneligliptin|anagliptin|saxagliptin,A10BH")

whoatc <- bind_rows(whoatc,
                    nonwho,
                    alternatives) %>% 
  distinct()
whoatc_lkp <- whoatc$atc_code
names(whoatc_lkp) <- whoatc$nm
setdiff(str_to_lower(arm_meta$drug_name), names(whoatc_lkp))
saveRDS(whoatc_lkp, "Scratch_data/who_atc_lkp.Rds")
arm_meta <- arm_meta %>% 
  mutate(drug_code  = whoatc_lkp[drug_name %>% str_to_lower()]) 

## Change OAD to A10B which is all oral antidiabetic drugs
arm_meta <- arm_meta %>% 
  mutate(drug_code = case_when(
    drug_code == "OAD" ~ "A10B",
    TRUE ~ drug_code))

## Policy for network simplification ----
## Want to preserve drug dose comparisons as i) these are in the IPD regression model exports and ii) this preserves information.
## However, if done indiscrimnately this will lead to a very large complex network and possibility of disconnection.
## Therefore perform the following simplification (discussed with GP (PH) and diabetes consultant (KH))
## Simplify insulins to a single code A10A
## Keep dose for metformin, SGLT2, GLP1 and DPP4
## drop dose for drug class level comparisons and drugs not in key classes
## Reduce to drug (ie remove dose and regimen) for the trials in these remaining classes. This involves an aggregation step
## which will do when merge data onto arms information.
## Need dose data. Where a trial is (rarely) missing drug dose information take the median
Commonest <- function(x) {
  a <- table(x)
  sort(a, decreasing = TRUE)[1] %>% as.integer()
}
## preserve original dose information
arm_meta <- arm_meta %>% 
  mutate(drug_dose_orig = drug_dose)
arm_meta <- arm_meta %>% 
  group_by(drug_code) %>% 
  mutate(commonest_dose =  Commonest(drug_dose) %>% as.character()) %>% 
  ungroup()
arm_meta <- arm_meta %>% 
  mutate(drug_dose_simplify = case_when(
    str_length(drug_code) == 4 ~ "remove",
    str_length(drug_code) == 5 ~ "remove",
    !str_sub(drug_code, 1, 5) %in% c("A10BA", "A10BH", "A10BJ", "A10BK") ~ "remove",
    drug_dose %in% c("missing", "unclear", "") | is.na(drug_dose) ~ commonest_dose,
    TRUE ~ drug_dose
  ))

## where multiple doses - if three take the middle, if 2 take the commonest
## relevant drugs for this are sitagliptin, dapagliflozin, canagliflozin  and dulaglutide
dose_commonest <- arm_meta %>% 
  filter(drug_name %in% c("sitagliptin", "dapagliflozin", "canagliflozin",  "dulaglutide")) %>% 
  count(drug_name, drug_dose_simplify) %>% 
  arrange(drug_name, desc(n)) %>% 
  group_by(drug_name) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(-n) %>% 
  rename(drug_dose_commonest = drug_dose_simplify)
arm_meta <- arm_meta %>% 
  left_join(dose_commonest)
arm_meta %>% filter(str_detect(drug_dose_simplify, "\\|"))
arm_meta <- arm_meta %>% 
  mutate(drug_dose_simplify = case_when(
    arm_id_unq == "uaa10193" ~ "100",
    arm_id_unq == "uaa10272" ~ "10",
    arm_id_unq == "uaa10513" ~ "100",
    arm_id_unq %in% c("uaa10528", "uaa10534", "uaa10594", "uaa10667") ~ "1.5",
    TRUE ~ drug_dose_simplify
  )) 

## convert all insulins to a single code
arm_meta <- arm_meta %>% 
  mutate(drug_code_simplify = case_when(
    str_detect(drug_code, "^A10A") ~ "A10A",
    TRUE ~ drug_code
  ))
## all fine. merge into a single drug_code drug dose variable
arm_meta <- arm_meta %>% 
  mutate(arm_lvl = if_else(drug_dose_simplify == "remove",
                           drug_code_simplify,
                           paste0(drug_code_simplify, "_d", drug_dose_simplify)))

arm_meta <- arm_meta %>% 
  mutate(trtcls5= str_sub(arm_lvl, 1, 5),
         trtcls4 = str_sub(arm_lvl, 1, 4))

## 112 unique labels
map_int(arm_meta %>% select(-nct_id, -arm_id_unq, -arm_id_subgroup, -arm_id), ~ sum(!duplicated(.x)))
# nct_id arm_id_unq  drug_code  drug_dose  drug_unit  drug_freq    trtcls5    trtcls4 
# 632       1552         90         68          8          9         27          5 
arm_meta <- arm_meta %>% 
            select(-drug_dose_commonest, -commonest_dose)

## check for duplicate arms. There are some as expected. Can collapse these ONLY after merging
## back in baseline and outcome data
arm_meta %>% 
  distinct(nct_id, arm_id_unq, arm_lvl)

## simplify regime names
arm_meta <- arm_meta %>% 
  left_join(whichnwork %>% distinct(trial_id, drug_regime) %>% rename(nct_id = trial_id))
arm_meta <- arm_meta %>% 
  mutate(drug_regime_smpl = case_when(
    drug_regime %in% c("triple", "triple+", "dual|triple+", "mono|dual|triple+") ~ "triple",
    drug_regime %in% c("dual", "mono|dual") ~ "dual",
    drug_regime %in% "mono" ~ "mono",
    TRUE ~ "unclear or missing"))

## add in trial where obtained arm data manually separately due to mismatch in arm labels----
regime_add <- read_csv(
  "nct_id,arm_id_unq,drug_code,drug_dose,drug_regime,drug_regime_smpl,drug_dose_simplify,arm_lvl
NCT01601990,noarmid1,placebo,'50',mono,mono,'50',placebo
NCT01601990,noarmid2,A10BH06,'50',mono,mono,'50',A10BH06_d50")
regime_toadd <- arm_meta %>% 
  filter(drug_code %in% c("placebo", "A10BH06")) %>% 
  distinct(drug_code, .keep_all = TRUE) %>% 
  select(-nct_id, -arm_id_unq, -drug_dose, -drug_regime,-drug_regime_smpl, -drug_dose_simplify,-arm_lvl) %>% 
  distinct()
regime_add <- regime_add %>% 
  left_join(regime_toadd)
arm_meta <- bind_rows(arm_meta,
                        regime_add)
rm(regime_add, regime_toadd)

## Collapse arm subgroup IDs so only one per trial/arm id combination----
dups <- arm_meta %>% 
  select(nct_id, arm_lvl, arm_id_unq) %>% 
  duplicated()
dups <- arm_meta %>% 
  select(nct_id, arm_lvl, arm_id_unq) %>% 
  filter(dups)
dups <- arm_meta %>% 
  semi_join(dups)
dups <- dups %>% 
  mutate(drug_dose = if_else(drug_name == "placebo", "", drug_dose)) %>% 
  group_by(nct_id, arm_id_unq, arm_id, drug_name) %>% 
  summarise_all(~ .x %>% unique() %>% sort() %>% paste(collapse = ",")) %>% 
  ungroup()
arm_meta <- arm_meta %>% 
  anti_join(dups %>% select(nct_id, arm_id_unq)) %>% 
  bind_rows(dups)

## Note. Need to address tirzepatide currently under ATC code A10BX16 change to A10BJ
## as it is a type of GLP1 (sort of) "Tirzepatide is a GIP-analogue that activates both the GLP-1 and GIP receptors"
arm_meta <- arm_meta %>% 
  mutate(trtcls5 = if_else(drug_name == "tirzepatide", "A10BJ", trtcls5))
arm_meta <- arm_meta %>% 
  mutate(across(everything(), ~ if_else(is.na(.x) | .x == "", NA_character_, .x )))

## need to add arm labels for empareg hba1c
arm_meta_2 <- arm_meta %>% filter(nct_id == "NCT01131676", is.na(drug_dose) | drug_dose != "10|25")  
arm_meta_2 <- arm_meta_2 %>% 
  mutate(arm_id_unq = case_when(
    arm_lvl == "A10BK03_d10" ~ "BI 10773 10mg",
    arm_lvl == "A10BK03_d25" ~ "BI 10773 25mg",
    arm_lvl == "placebo" ~ "placebo"
  ))
arm_meta <- bind_rows(arm_meta,
                      arm_meta_2)
write_csv(arm_meta, "Scratch_data/arm_labels_hba1c.csv")
## no change to arm_meta on adding bocf and rcs
