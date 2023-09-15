library(tidyverse)
library(multinma)

## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- readRDS("Scratch_data/agg_hba1c.Rds")

## drop agg where have ipd - 49 trials
agg %>% 
  semi_join(ipd %>% select(nct_id)) %>% 
  count(nct_id)
agg <- agg %>% 
  filter(!nct_id %in% ipd$nct_id)

# read in arm metadata
# as uploaded to vivli
arm_meta_orig <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)
# additions made within vivli
arm_meta_new<- read_csv("../from_vivli/Data/reference_arm_data_all_cleaned.csv")
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
## Add placebo as a drug name when it is only listed as a lable
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
forcode <- arm_meta %>% 
  select(nct_id, arm_id_unq, drug_name, drug2_name) %>% 
  gather("drug_order", "drug_name", drug_name, drug2_name, na.rm = TRUE)
## Split out where have pipe separator
forcode$drug_name <- map(forcode$drug_name, ~ str_split(.x, "\\|")[[1]])
forcode <- forcode %>% 
  unnest(drug_name) %>% 
  mutate(drug_name = drug_name %>% 
           str_remove("�") %>% 
           stringi::stri_enc_toutf8() %>% 
           str_trim() %>% 
           str_to_lower()) %>% 
  filter(!is.na(drug_name), !drug_name == "")
rename_drug <- readxl::read_excel("Created_metadata/non_atc_drugs.xlsx", sheet = 1)
## note this increases rows as some drugs were separated using underscores not pipes in data
forcode <- forcode %>% 
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
forcode <- forcode %>% 
  mutate(drug_code  = whoatc_lkp[drug_name]) %>% 
  distinct(nct_id, arm_id_unq, drug_code)

## drop any where same drug in every arm
forcode_drp <- forcode %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  ungroup()
IdentifyCrossArms <- function(arm_drug){
  x <- arm_drug %>% 
    select(arm_id_unq, drug_code) %>%
    mutate(v = 1L) %>% 
    spread(drug_code, v, fill = 0L) %>% 
    select(-arm_id_unq) %>% 
    summarise_all(all)
  x <- unlist(x)
  names(x[x])  
}
forcode_drp$sameacross <- map(forcode_drp$data, IdentifyCrossArms)
# forcode_drp$data <- NULL
forcode_drp$n <- map_int(forcode_drp$sameacross, length)
forcode_drp <- forcode_drp %>% 
  filter(n >0) %>% 
  unnest(sameacross)
## idetify implict controls and where present add a placevo
forcode_drp$impcntrl <- map(forcode_drp$data, function(a) {
  a %>% 
    mutate(v = 1L) %>% 
    spread(drug_code, v, fill = 0L) %>% 
    gather(key = "drug_code", cntrl, -arm_id_unq) %>% 
    mutate(drug_code = if_else(cntrl == 0L, 
                               "implicit_control",
                               drug_code))
})
forcode_drp$sameacross2 <- map(forcode_drp$impcntrl, IdentifyCrossArms)
forcode_drp$retain <- map2(forcode_drp$impcntrl, forcode_drp$sameacross2, ~ {(  
  .x %>% 
    filter(!drug_code %in% .y))
})
forcode_drp$retained <- map_int(forcode_drp$retain, nrow) 
forcode_keep <- forcode_drp %>% 
  filter(!retained ==0)
forcode_drp <- forcode_drp %>% 
  filter(retained ==0)
forcode_drp <- forcode_drp %>% 
  select(nct_id, sameacross, data) %>% 
  unnest(data) 
## 5 trials with implicit controls
forcode_keep <- forcode_keep %>% 
  select(nct_id, sameacross, retain) %>% 
  unnest(retain) %>% 
  select(-cntrl) 

forcode2 <- bind_rows(forcode %>% 
                        filter(!nct_id %in% forcode_keep$nct_id,
                               !nct_id %in% forcode_drp$nct_id) %>% 
                        mutate(sameacross = "nil here see ancillary"),
                      forcode_keep)

saveRDS(list(keep = forcode2, drop = forcode_drp), 
        "Scratch_data/arm_codes_sameacross_retain.Rds")

arm_assign <- forcode2 %>% 
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
              paste(.x, collapse = "_"))) %>% 
  ungroup()
map_int(arm_assign, ~ sum(!duplicated(.x)))
# nct_id arm_id_unq  drug_code    trtcls5    trtcls4 
# 782       1899        116         51         17 
saveRDS(arm_assign,  "Scratch_data/arm_labels_final.Rds")

