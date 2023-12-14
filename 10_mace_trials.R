#10_mace_trials

library(tidyverse)
library(multinma)

a <- read_csv("../cleaned_data/Data/example_trial_level_info.csv")
arm_regime <- readRDS("Scratch_data/arms_assign_drug_regime.Rds")
whichnwork <- read_csv("../cleaned_data/Data/ancillary_drugs_data_all_cleaned.csv")
arm_meta <- read_csv("../cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)
ipd <- read_csv("../ipd_overview/Outputs/MACE_trials.csv")

base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
  rename(nct_id = trial_id)
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds") %>% 
  rename(nct_id = trial_id)

# resgot1 <- read_csv("Data/availability_of_results_2022_update.csv")
# resgot2 <- read_csv("Data/availability_of_results_by_trial_id.csv") %>% 
#   mutate(trial_id = if_else(is.na(nct_id), other_id, nct_id))
# all(resgot2$trial_id %in% resgot1$trial_id)
# resgot <- resgot2
# rm(resgot1, resgot2)
# resgot <- resgot %>% 
#   mutate(nct_id = if_else(is.na(nct_id), other_id, nct_id)) %>% 
#   select(-other_id)
mace <- a %>% 
  filter(outcome_to_be_modelled == "MACE")  %>% 
  distinct(nct_id, other_id, outcome_to_be_modelled) %>% 
  mutate(nct_id = if_else(is.na(nct_id), other_id, nct_id)) %>% 
  select(nct_id)
mace <- mace %>% 
  filter(!nct_id %in% c("NCT02655757", "UMIN000025102", "NCT01018173"))

## have baseline data for all MACE
base_dsp <- base_dsp %>% 
  semi_join(mace)
base_rng <- base_rng %>% 
  semi_join(mace)
mace %>% 
  anti_join(base_dsp)


## all but one is triple
mace <- mace %>% 
  left_join(arm_regime %>% 
              distinct(nct_id, drug_regime_smpl) %>% 
              group_by(nct_id) %>% 
              summarise(drug_regime_smpl = 
                          paste(drug_regime_smpl %>% unique(), collapse = ",")) %>% 
              ungroup())
mace <- mace %>% 
  left_join(whichnwork %>% select(nct_id = trial_id, drug_regime))
mace <- mace %>% 
  mutate(drug_regime_smpl = case_when(
    !is.na(drug_regime_smpl) ~ drug_regime_smpl,
    str_detect(drug_regime, "triple") ~ "triple",
    str_detect(drug_regime, "dual") ~ "dual",
    str_detect(drug_regime, "mono") ~ "mono"))

## all but one is triple so I think we should just lump them
## There are no double drug arms
arm_meta <- arm_meta %>% 
  semi_join(mace)
mace <- mace %>% 
  inner_join(arm_meta %>% 
               filter(!is.na(drug_name)) %>% 
               select(nct_id, arm_id_unq, drug_name, drug_dose) %>% 
               distinct())

## how many are placebo controlled
mace <- mace %>% 
  group_by(nct_id) %>% 
  mutate(plac = any(drug_name %in% c("placebo", "control"))) %>% 
  ungroup() %>% 
  arrange(plac)
## 21 are placebo controlled, 1 is "control" and 1 is not. It is linked to the network via linagliptin
## count arms
mace_rv <- mace %>% 
  group_by(nct_id) %>% 
  summarise(arms = length(arm_id_unq),
         drug_unq = sum(!duplicated(drug_name))) %>% 
  ungroup()

## 18 two arm trials, 2 entries (and one is placbo)
## 4 four arm trials but with 2 unique drugs. Review thiese
## note ITCA_650 is an implantable GLP-1
mace4arm <- mace %>% 
  semi_join(mace_rv %>% filter(arms ==4)) 

# new arm labels in WHO database and own manual lookup
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
setdiff(str_to_lower(mace$drug_name), names(whoatc_lkp))
mace <- mace %>% 
  mutate(atc = whoatc_lkp[drug_name %>% str_to_lower()])
mace <- mace %>% 
  mutate(trtcl5 = str_sub(atc, 1, 5))
mace <- mace %>% 
  left_join(whoatc %>% 
               distinct(atc_code, .keep_all = TRUE) %>% 
               rename(trtcl5 = atc_code)) %>% 
  mutate(dc = if_else(trtcl5 == "place", 
                      "placebo",
                      paste(trtcl5, nm, sep = " ")))
mace <- mace %>% 
  mutate(ipd = if_else(nct_id %in% ipd$nct_id, 1L, 0L))
ns <- read_csv("../cleaned_data/Data/example_trial_level_info.csv")
mace <- mace %>% 
  group_by(nct_id) %>% 
  mutate(n_arms = sum(!duplicated(arm_id_unq))) %>% 
  ungroup() 
mace <- mace %>% 
  left_join(ns %>% select(nct_id, enrollment))
mace <- mace %>% 
  mutate(per_arm = as.integer(enrollment/n_arms))

## make table summary
plac_tbl <- mace %>% 
  filter(plac, !
           drug_name %in% c("placebo", "control")) %>% 
  arrange(dc, drug_name, drug_dose, nct_id, enrollment, ipd) %>% 
  select(nct_id, dc, drug_name, drug_dose, enrollment, ipd) %>% 
  group_by(nct_id, dc, drug_name, enrollment, ipd) %>% 
  summarise(drug_dose = paste(drug_dose, collapse = " and ")) %>% 
  ungroup() %>% 
  left_join(ipd %>% 
              select(nct_id, repo = data_request_id)) %>% 
  mutate(repo = if_else(is.na(repo), "", repo)) %>% 
  select(nct_id, dc, drug_name, drug_dose, enrollment, ipd, repo) 
noplac <- mace %>% 
  filter(!plac)
noplac_tbl <- noplac %>% 
  slice(2) %>% 
  inner_join(noplac %>% 
               slice(1) %>% 
               select(nct_id, drug_cmpr = drug_name)) %>%
  mutate(drug_name = paste0(drug_name, " vs ", drug_cmpr),
         repo = "") %>% 
  select(nct_id, dc, drug_name, drug_dose, enrollment, ipd, repo) 
mace_tbl <- bind_rows(plac_tbl,
                      noplac_tbl) %>% 
  arrange(dc, drug_name, drug_dose, nct_id, enrollment, ipd) 
sponsor <-  read_csv("../cleaned_data/Data/example_trial_level_info.csv")
sponsor <- sponsor %>% 
  mutate(nct_id = if_else(is.na(nct_id), other_id, nct_id)) %>% 
  semi_join(mace)
## nct_id, sponsor, primary_outcome, other_relevant_outcomes, primary_or_secondary_outcome, 
# phase, brief_title, starts_with("subgroup"), time_to_outcome_weeks, blinding, population_specifics, minimum_age, maximum_age
## mace_tbl - classes and drugs
mace_tbl1 <- mace_tbl
mace_tbl_outs <- mace_tbl %>% 
  select(nct_id, dc, ipd) %>% 
  left_join(sponsor %>% select(nct_id, outcome_to_be_modelled, primary_or_secondary_outcome, 
                               time_to_outcome_weeks, 
                               primary_outcome, other_relevant_outcomes))
# sponsor %>% 
#   filter(nct_id == "NCT01455896")

NCT01455896 <- read_csv("sg,present
age,1
bp,1
hba1c,1
duration_dm,1
race,1
ethnicity,1
gender,1
dm_therapy,1
bmi,1
cpeptide,0
egfr,1
geolocation,1
CVdisease,1
smoking,1
other,1")
NCT01455896_other <- "ldl, trig, hsCRP"

sg_smry <- mace_tbl %>% 
  select(nct_id, dc, ipd) %>% 
  left_join(sponsor %>% 
              select(nct_id, subgroups_yn) %>% 
              filter(!nct_id == "NCT01455896") %>% 
              bind_rows(tibble(nct_id = "NCT01455896", subgroups_yn = "y")))
sg_smry <- sg_smry %>% 
  count(ipd, subgroups_yn)
write_csv(sg_smry, "Outputs/mace_subgroup_summary.csv")
mace_tbl_sg <- mace_tbl %>% 
  select(nct_id, dc, ipd) %>% 
  left_join(sponsor %>% 
              select(nct_id, subgroups_yn, subgroup_names) %>% 
              filter(!nct_id == "NCT01455896") %>% 
              bind_rows(tibble(nct_id = "NCT01455896", subgroups_yn = "y", 
                               subgroup_names = paste(
                                 paste0(NCT01455896$sg[NCT01455896$present == 1L], collapse = ", "),
                                 ", ",
                                 NCT01455896_other))))
mace_plt_sg <- mace_tbl %>% 
  # filter(ipd == 0L) %>% 
  select(nct_id, dc, ipd) %>% 
  left_join(sponsor %>% 
              select(nct_id, subgroups_yn, contains("subgroup_")) %>% 
              select(-subgroup_names, -subgroup_cpeptide)) %>% 
  filter(subgroups_yn == "y") %>% 
  select(-subgroups_yn) %>% 
  gather("subgroup", "present", -nct_id, -dc, -ipd, na.rm = FALSE) %>% 
  mutate(present = if_else(is.na(present), "No", "Yes"),
         subgroup = str_sub(subgroup, 10))
NCT01455896b <- mace_tbl %>% 
  filter(nct_id == "NCT01455896") %>% 
  select(nct_id, dc, ipd) %>% 
  cross_join(NCT01455896) %>% 
  rename(subgroup = sg) %>% 
  mutate(present = if_else(present == 1, "Yes", "No")) %>% 
  filter(!subgroup == "cpeptide")
mace_plt_sg <- bind_rows(mace_plt_sg,
                         NCT01455896b)
mace_plt_sg <- mace_plt_sg %>% 
  mutate(dc2 = if_else(present == "Yes", dc, NA_character_))
lvls <- mace_plt_sg %>% 
  arrange(dc, nct_id) %>% 
  distinct(dc, ipd, nct_id) %>% 
  group_by(ipd) %>% 
  mutate(rdr = letters[seq_along(dc)]) %>% 
  ungroup()
mace_plt_sg <- mace_plt_sg %>% 
  inner_join(lvls)
plot1 <- ggplot(mace_plt_sg %>% 
                  mutate(ipd = if_else(ipd ==1L, "IPD", "Aggregate")), aes(x = rdr, y = subgroup, fill = dc2)) +
  geom_tile(colour = "white") +
  scale_x_discrete("", labels = NULL) +
  facet_wrap(~ipd) +
  scale_fill_discrete("") +
  scale_y_discrete("") +
  theme_minimal()
plot1

mace_tbl_other <- mace_tbl %>% 
  select(nct_id, dc, ipd) %>% 
  inner_join(sponsor %>% select(nct_id, phase, brief_title, population_specifics, minimum_age, maximum_age, other_relevant_outcomes))
# phase, brief_title, blinding, population_specifics, minimum_age, maximum_age
# starts_with("subgroup"),
## generate network

write_csv(mace_tbl1, "Outputs/mace_tbl_drugs.csv")
write_csv(mace_tbl_outs, "Outputs/mace_tbl_outcomes.csv")
write_csv(mace_tbl_sg, "Outputs/mace_tbl_subgroups.csv")
write_csv(mace_tbl_other, "Outputs/mace_tbl_other.csv")
saveRDS(plot1, "Scratch_data/heatmap_subgroups.Rds")


### draw network ----
mace <- mace %>% 
  mutate(drug_name = if_else(drug_name == "control", "placebo", drug_name),
         drug_dose = if_else(is.na(drug_dose), "", drug_dose))
mace_nst <- mace %>% 
  inner_join(mace_tbl %>% select(nct_id, dc_main = dc)) %>% 
  arrange(nct_id, drug_name) %>% 
  group_by(dc_main) %>% 
  nest() %>% 
  ungroup()
mace_nst$plt <- map(mace_nst$data, function(mace){
  mace_ipd <- mace %>% 
    filter(ipd == 1L)
  mace_agg <- mace %>%
    filter(!ipd == 1L)
  mace_agg <- mace_agg %>%
    mutate(per_arm = if_else(is.na(per_arm), 1, per_arm),
           r = rbinom(nrow(.), size = per_arm, p = 0.1))
  mace_ipd$ipd <- map(mace_ipd$per_arm, ~ tibble(r = rbinom(.x, size = 1, p = 0.1))) 
  mace_ipd <- mace_ipd %>%
    unnest(ipd)
  makenetdrugdose <- combine_network(
    set_ipd(
      mace_ipd,
      study = nct_id,
      trt = paste(drug_name, drug_dose),
      r = r
    ),
    set_agd_arm(
      mace_agg,
      study = nct_id,
      trt = paste(drug_name, drug_dose),
      r = r,
      sample_size = per_arm,
      n = per_arm
    )
  )
  nwork <- plot(makenetdrugdose)
  nwork
})
saveRDS(mace_nst %>% select(-data), "Scratch_data/network_by_class.Rds")

## got 16 of 20 results for MACE here
## read in data extract
agg <- read_tsv("Data/mace_example_extract")
main <- agg %>% 
  filter(subgroup == "main", !arm_label == "total") %>% 
  select(-p_interaction, -contains("level"), -subgroup)


## pull mace ipd aact
# mace_clean <- read_csv("../cleaned_data/Data/mace_extraction_example.csv")
# mace_clean_smry <- mace_clean %>% 
#   distinct(trial_id, result_type) %>% 
#   mutate(v = 1L) %>%  
#   spread(result_type, v, fill = 0L) %>% 
#   rename(nct_id = trial_id)
# mace_clean2 <- mace_clean %>% 
#   select(nct_id = trial_id, 
#          outcome, 
#          arm_id_unq,
#          arm_label,
#          result_type, result,
#          dispersion_type, dispersion)
# mace_clean2_d <- mace_clean2 %>% 
#   select(-result_type,-result) %>% 
#   spread(dispersion_type, dispersion)

## more outcome data ----
aact <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")
aact <- map(aact, ~{
  if("nct_id" %in% names(.x)) {
    .x %>% 
      semi_join(mace)
  } else .x
})

## have events at arm level for all 19 trials in AACT 
events <- aact$outcome_measurements %>% 
  filter( (title %>% str_to_lower() %>% str_detect("event") |
             category == "Number of patients with an event"))
percentage <- aact$outcome_measurements %>% 
                    anti_join(events %>% select(nct_id)) 
result_groups <- aact$result_groups
## all are treatment arms, not totals
events2 <- events %>% 
  inner_join(result_groups %>% 
               rename(result_group_id = id), by = "result_group_id")
percentage2 <- percentage %>% 
  inner_join(result_groups %>% 
             rename(result_group_id = id), by = "result_group_id")

## remaining trial has events even in abstract "We randomized 4,156 patients (2,075 assigned to
#receive ITCA 650 and 2,081 assigned to receive placebo) who were followed for a
#median of 16 months. The primary outcome occurred in 4.6% (95/2,075) of
#patients in the ITCA 650 group and 3.8% (79/2,081) of patients in the placebo
#group, meeting the pre-specified non-inferiority criterion (HR = 1.21, 95% CI,
#0.90-1.63, Pnon-inferiority = 0.004)."

mace %>% 
  filter(!nct_id %in% c(events$nct_id, percentage$nct_id))

## Note canvas and canvas-R both stratified according to prev CVD versus CV risk
## analysis stratified by trial (both trials together) and CVD/CVrisk
#NCT01032629 - CVD/CVrisk 
# NCT01131676 - stratified according to the glycated
#  hemoglobin level at screening (<8.5% or ≥8.5%), body-mass index at
#  randomization (<30 or ≥30), renal function at screening (eGFR, 30 to 59 ml, 60
#  to 89 ml, or ≥90 ml per minute per 1.73 m2), and geographic region (North
#  America [plus Australia and New Zealand], Latin America, Europe, Africa, or
#  Asia).
# NCT01989754 - CVD/CVrisk 
# NCT02065791 category of estimated GFR (30 to <45 ml, 45 to <60 ml, or 60 to <90 ml per minute per 1.73 m2) at screening
# NCT02465515 - no stratification

