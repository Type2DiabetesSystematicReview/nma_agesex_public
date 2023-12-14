library(tidyverse)
library(multinma)

source("Scripts/00_functions.R")
## read data ----
## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- read_csv("Data/agg.csv")
# read in arm metadata
# as uploaded to vivli
arm_meta_orig <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)
arm_meta_orig <- arm_meta_orig %>% 
  semi_join(bind_rows(agg %>% select(nct_id, arm_id_unq),
                      ipd %>% select(nct_id, arm_id_unq)))

# additions to trial metadata made within vivli ----
arm_meta_new <- read_csv("../from_vivli/Data/agesexhba1c_6115/reference_arm_data_all_cleaned.csv")
arm_meta_new %>% anti_join(arm_meta_orig)
arm_meta_orig %>% anti_join(arm_meta_new)
arm_meta <- bind_rows(arm_meta_orig,
                      arm_meta_new %>% 
                        anti_join(arm_meta_orig) %>% 
                        mutate(drug_name = case_when(
                          arm_id_unq %in% c("ipd00002", "ipd00004") ~ "placebo",
                          arm_id_unq %in% c("ipd00003", "ipd00001") ~ "linagliptin"),
                               drug_dose = case_when(
                                 arm_id_unq %in% c("ipd00001") ~ "25",
                                 arm_id_unq %in% c("ipd00003") ~ "10"))
                        )
rm(arm_meta_orig, arm_meta_new)

## Add placebo as a drug name when it is only listed as a label ----
arm_meta <- arm_meta %>% 
  mutate(drug_name = if_else(str_sub(arm_label) %>% str_to_lower() == "placebo" & 
                               (is.na(drug_name) | drug_name == ""), "placebo", drug_name))
# Correct error in coding for one drug (shifted along in excel by one cell)
arm_meta <- arm_meta %>% 
  mutate(drug2_dose = if_else(arm_id_unq == "updac0011" & drug2_name == 10,
                              "10",
                              drug2_dose),
         drug2_unit = if_else(arm_id_unq == "updac0011" & drug2_name == 10,
                              "mg",
                              drug2_unit),
         drug2_freq = if_else(arm_id_unq == "updac0011" & drug2_name == 10,
                              "od",
                              drug2_freq),
         drug2_name = if_else(arm_id_unq == "updac0011" & drug2_name == 10,
                              "dapagaliflozin",
                              drug2_name))
arm_meta <- arm_meta %>% 
  mutate(drug_name = if_else(arm_id_unq %in% c("updac0096", "updac0153", "updac0177", "updac0126"),
                             "placebo",
                             drug_name))
saveRDS(arm_meta, "Scratch_data/drug_names_doses_regimen.Rds")

## Drop drugs which are common to all groups and label arms uniquely ----
# note that where a single drug is present in multiple arms,
# and where one arm has only that drug, we need to add the name "control" or else the whole arm
# will be dropped from the data (When it is in fact a control arm)
# eg:-
# arm 1 metformin + SGLT2
# arm 2 metformin
# becomes:-
# arm 1 SGLT2
# arm 2 control
re_arm <- arm_meta %>% 
  select(nct_id, arm_id_unq, drug_name, drug2_name) %>% 
  gather("drug_order", "drug_name", drug_name, drug2_name, na.rm = TRUE)
## Split out where have pipe separator
re_arm$drug_name <- map(re_arm$drug_name, ~ str_split(.x, "\\|")[[1]])
re_arm <- re_arm %>% 
  unnest(drug_name) %>% 
  mutate(drug_name = drug_name %>% 
           str_remove("�") %>% 
           stringi::stri_enc_toutf8() %>% 
           str_trim() %>% 
           str_to_lower()) %>% 
  filter(!is.na(drug_name), !drug_name == "")
rename_drug <- readxl::read_excel("Created_metadata/non_atc_drugs.xlsx", sheet = 1)

## note this increases rows as some drugs were separated using underscores not pipes in data
re_arm <- re_arm %>% 
  left_join(rename_drug, relationship = "many-to-many") %>% 
  mutate(drug_name = if_else(!is.na(rename), rename, drug_name)) %>% 
  select(-rename)

## new arm labels in WHO database and own manual lookup
whoatc <- readxl::read_excel("../../../Medications_resources/WHO_ATC/2018 ATC index with DDDs.xlsx", sheet = 1)
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
whoatc <- bind_rows(whoatc,
                    nonwho) %>% 
  distinct()
whoatc_lkp <- whoatc$atc_code
names(whoatc_lkp) <- whoatc$nm
setdiff(str_to_lower(re_arm$drug_name), names(whoatc_lkp))
re_arm <- re_arm %>% 
  mutate(drug_code  = whoatc_lkp[drug_name %>% str_to_lower()]) %>% 
  distinct(nct_id, arm_id_unq, drug_code)

## Simplify insulins to a single code A10A; all are single drugs
re_arm <- re_arm %>% 
  mutate(drug_code = if_else(str_detect(drug_code, "^A10A"), "A10A", drug_code))
re_arm <- re_arm %>% 
  distinct()

## Where class is specified rather than drug (ie multiple drugs in the same class in an arm)
## simplify the drug code to the class. eg multiple glp1s
re_arm <- re_arm %>% 
  mutate(trtcls5 = str_sub(drug_code, 1, 5)) %>% 
  group_by(nct_id, arm_id_unq, trtcls5) %>% 
  mutate(appears = length(trtcls5),
         drug_code = if_else(appears >=2, trtcls5, drug_code)) %>% 
  ungroup()  %>% 
  distinct(nct_id, arm_id_unq, drug_code)

mylst <- SimplifyDrugs(re_arm)
## Reviewed ancillary drugs detected wiht comparison in code versus the manual assignment with ELB.
## 9th November 2023
## Consistent results. No need to change designation (eg mono, dual triple etc)
saveRDS(mylst, 
        "Scratch_data/arm_codes_sameacross_retain.Rds")
arm_assign <- re_arm %>% 
  select(nct_id, arm_id_unq, drug_code) %>% 
  mutate(drug_code = if_else(drug_code == "implicit_control",
                             "placebo",
                             drug_code))
## 113 unique labels
arm_assign <- arm_assign %>% 
  mutate(trtcls5= str_sub(drug_code, 1, 5),
         trtcls4 = str_sub(drug_code, 1, 4)) %>% 
  arrange(nct_id, arm_id_unq, drug_code) %>% 
  group_by(nct_id, arm_id_unq) %>% 
  summarise(across(c(drug_code, trtcls5, trtcls4), ~
              paste(.x %>% unique(), collapse = "_"))) %>% 
  ungroup()

map_int(arm_assign, ~ sum(!duplicated(.x)))
# nct_id arm_id_unq  drug_code    trtcls5    trtcls4 
# 602       1489         84         26          5 
write_csv(arm_assign, "Data/arm_labels_hba1c.csv")
