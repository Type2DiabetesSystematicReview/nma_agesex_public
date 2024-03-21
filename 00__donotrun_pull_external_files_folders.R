# 000_pull_external_files_folders
## Pull external files and folders. Only re-run when want to incorporate updates in data 
## external to this folder.
## Will only be possible to run this if the files/folders are on the users machine.
## All subsequent scripts should run anywhere the github repository ahs been downloaded

library(tidyverse)

## read in hba1c vivli and gsk results and copy to the Data folder ----
allvivli <- list.files("../from_vivli/Data/agesexhba1c_6115/", patt = "csv$")
from <- paste0("../from_vivli/Data/agesexhba1c_6115/", allvivli)
dir.create("Data/agesexhba1c_6115")
to <- paste0("Data/agesexhba1c_6115/", allvivli)
file.copy(from, to)
allvivli <- list.files("../from_vivli/Data/agesexhba1c_8697/", patt = "csv$")
from <- paste0("../from_vivli/Data/agesexhba1c_8697/", allvivli)
dir.create("Data/agesexhba1c_8697")
to <- paste0("Data/agesexhba1c_8697/", allvivli)
file.copy(from, to)
allgsk <- list.files("../from_gsk//Data/agesex/", patt = "csv$")
from <- paste0("../from_gsk//Data/agesex/", allgsk)
dir.create("Data/gsk")
to <- paste0("Data/gsk/", allgsk)
file.copy(from, to)

## Read in MACE data, vivli plus a single file giving an overview of the IPD ---- 
## event times separate for single trial and other 4. Censoring times brought together early because sample with sd = 0 for categorical age
## trials on new vivli repository - 8697
fls <- c(
  # "event_time_distribution_linear.csv",
  "event_time_distribution_fp.csv",
  "event_time_distribution_single_trial.csv",
  "censoring_distribution.csv",
  "censoring_distribution_single_trial.csv")
res <- map(fls, ~ read_csv(paste0("../from_vivli/Data/agesexmace_8697/", .x)))
names(res) <- str_sub(fls, 1, -5)
list2env(res, envir = .GlobalEnv)
rm(res)
## trials on new previous repository - 6115
censoring_distribution_fp_6115 <- read_csv("../from_vivli/Data/agesexmace_6115/censoring_distribution_fp.csv")
event_time_distribution_fp_6115 <- read_csv("../from_vivli/Data/agesexmace_6115/event_time_distribution_fp.csv")
censoring_distribution_fp_6115 <- censoring_distribution_fp_6115 %>% 
  rename(arm = arm_label,
         sex = sex_decoded) %>% 
  mutate(sex = if_else(sex == "female", "F", "M"),
         `Number of quantiles` = str_count(`Censored at quantiles (%)`, "\\="))
censoring_distribution <- bind_rows(censoring_distribution,
                                    censoring_distribution_fp_6115)
rm(censoring_distribution_fp_6115)
event_time_distribution_fp_6115 <- event_time_distribution_fp_6115 %>% 
  rename(arm = arm_label,
         sex = sex_decoded) %>% 
  mutate(sex = if_else(sex == "female", "F", "M")) %>% 
  select(-fu_m, -fu_t) %>% 
  mutate(est_age2 = 0,
         se_age2 = 0)
event_time_distribution_fp <- bind_rows(event_time_distribution_fp,
                                        event_time_distribution_fp_6115 %>% mutate(r = as.character(r)))
rm(event_time_distribution_fp_6115)
dir.create("Data/vivli_mace/")
write_csv(censoring_distribution, "Data/vivli_mace/censoring_distribution.csv")
write_csv(censoring_distribution_single_trial, "Data/vivli_mace/censoring_distribution_single_trial.csv")
write_csv(event_time_distribution_fp, "Data/vivli_mace/event_time_distribution_fp.csv")
write_csv(event_time_distribution_single_trial, "Data/vivli_mace/event_time_distribution_single_trial.csv")

cfs <- bind_rows(`6115` = read_csv("../from_vivli/Data/agesexmacecentred_6115/age_sex_model_coefs.csv"),
                 `8697` = read_csv("../from_vivli/Data/agesexmacecentred_8697/age_sex_model_coefs.csv"),
                 .id = "repo")
write_csv(cfs, "Data/vivli_mace/model_coefficients.csv")
vcv <- bind_rows(`6115` = read_csv("../from_vivli/Data/agesexmacecentred_6115/age_sex_model_vcov.csv"),
                 `8697` = read_csv("../from_vivli/Data/agesexmacecentred_8697/age_sex_model_vcov.csv"),
                 .id = "repo")
write_csv(vcv, "Data/vivli_mace/model_vcv.csv")


ipd_meta <- read_csv("../ipd_overview/Outputs/MACE_trials.csv")
write_csv(ipd_meta, "Data/vivli_mace/Overview_MACE_trials.csv")

ipd_cnsr <- bind_rows(`6115` = read_csv("../from_vivli/Data/agesexmace_6115/censoring_distribution_fp.csv") %>% 
                        rename(arm = arm_label,
                               sex = sex_decoded),
                      `8697a` = read_csv("../from_vivli/Data/agesexmace_8697/censoring_distribution.csv"),
                      `8697b` = read_csv("../from_vivli/Data/agesexmace_8697/censoring_distribution_single_trial.csv"),
                      .id = "repo") %>% 
  mutate(male = if_else(sex %in% c("F", "Female"), 0L, 1L)) %>% 
  select(-sex)   %>% 
  mutate(arm = str_to_lower(arm))
write_csv(ipd_cnsr, "Data/vivli_mace/ipd_cnsr.csv")

## read in AACT results and copy to the data folder ----
oc_orig <- readRDS("../extract_transform/aact/data/aact_extract.Rds")$outcome_counts
dir.create("Data/aact")
write_csv(oc_orig, "Data/aact/oc_orig.csv")
oc_new <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$outcome_counts 
write_csv(oc_new, "Data/aact/oc_new.csv")
aact_age <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$baseline_measurements
aact_age <- aact_age %>% 
  filter(title == "Age, Continuous",
         param_type == "Mean",
         dispersion_type == "Standard Deviation",
         units %in% c("years", "Years")) %>% 
  select(nct_id, result_group_id, ctgov_group_code, 
         param_value_num, dispersion_value_num, number_analyzed)
write_csv(aact_age, "Data/aact/age.csv")
## Following function was to pull results and write to a file temp.csv
## This was then reviewed manually, renamed and read in, inline in a script using read_csv witha  string file
DoNotDeleteDoNotRun <- function () {aact_result_groups <- 
  readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$result_groups  %>% 
  semi_join(aact_age %>% select(nct_id, ctgov_group_code)) 
aact_result_groups %>% 
  distinct(nct_id, ctgov_group_code, title) %>% 
  mutate(arm_id = "") %>% 
  arrange(nct_id, arm_description) %>% 
  write_csv("temp.csv") }
NCT03496298 <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")
NCT03496298 <- map(NCT03496298, ~{
  if("nct_id" %in% names(.x)) {
    .x %>% 
      filter(nct_id == "NCT03496298")
  } else .x
})
saveRDS(NCT03496298, "Data/aact/NCT03496298.Rds")

## Read in from extract_transform data cleaning of raw data and data from AACT ----
hba1c_ids_orig <- read_csv("../extract_transform/Created_metadata/hba1c_continuous_aact.csv") %>% 
   distinct(nct_id, id, outcome_id, result_group_id, ctgov_group_code)
write_csv(hba1c_ids_orig, "Data/extract_transform/hba1c_continuous_aact.csv")
## note lacking IDs. Need to add from aact
hba1c_ids_new <- read_csv("../extract_transform/aact/2022_update/Created_metadata/unique_outcomes_for_harmonisation_RVd_EB.csv") %>%
  filter(variable == "hba1c") 
hba1c_ids_new2 <- readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$outcome_measurements
hba1c_ids_new <- hba1c_ids_new %>% 
  inner_join(hba1c_ids_new2 %>% 
               select(nct_id, title, description, outcome_id, units, id, result_group_id, ctgov_group_code) %>% 
               distinct())
rm(hba1c_ids_new2)
dir.create("Data/extract_transform/")
write_csv(hba1c_ids_new, "Data/extract_transform/hba1c_ids_new.csv")
ci_corr <- read_csv("../extract_transform/data/aact_hba1c_outcome_data_correct_aact_difference.csv")
write_csv(ci_corr, "Data/extract_transform/ci_corr.csv")

### Read in files from cleaned_data/ repository. This is downstream of extract_data/ ----
cleaned <- read_csv("Data/table_for_copying_files_into_directory.csv")
cleaned <- cleaned %>% 
  filter(top_level == "cleaned_data")
dir.create("Data/cleaned_data/")
## processed data
base_dsp <- readRDS("../cleaned_data/Processed_data/base_dsp.Rds") %>% 
  filter(variable %in% c("n", "age", "male", "race", "ethnicity"))
write_csv(base_dsp, "Data/cleaned_data/base_dsp.csv")
base_rng <- readRDS("../cleaned_data/Processed_data/base_rng.Rds") %>% 
  filter(variable %in% c("age"))
write_csv(base_rng, "Data/cleaned_data/base_rng.csv")
cleaned <- cleaned %>% 
  filter(!str_detect(subsequent, "^Processed"))
tocopy <- cleaned %>% 
  distinct(within, top_level, subsequent) %>% 
  separate(subsequent, into = c("datafolder", "filename"), sep = "/")
dir.create("Data/cleaned_data/Data/")
tocopy <- tocopy %>% 
  mutate(newname = paste0("Data/cleaned_data/Data/", filename))
file.copy(tocopy$within, tocopy$newname, overwrite = TRUE)

## read in scripts from "common_functions" that we use in analysis
common <- read_csv("Data/table_for_copying_files_into_directory.csv")
common <- common %>% 
  filter(top_level == "common_functions")
fromcopy <- common %>% distinct(within) %>% pull(within)
dir.create("Scripts/common_functions/")
dir.create("Scripts/common_functions/Scripts")
tocopy <- str_replace(fromcopy, "\\.{2,2}/", "Scripts/")
file.copy(fromcopy, tocopy)

## Read in WHO ATC table. Only pull atc codes and names. Not doses or routes ----
whoatc <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx", sheet = 1)
whoatc <- whoatc %>% 
  filter(str_detect(`ATC code`, "^A10")) %>% 
  select(1:2)
write_csv(whoatc, "Data/whoatcdiabetesnodose.csv")

## copy across readmes
file.copy("../from_vivli/Data/agesexhba1c_6115/00_readme.txt",
          "Data/agesexhba1c_6115/00_readme.txt")
file.copy("../from_vivli/Data/agesexhba1c_8697/readme.txt",
          "Data/agesexhba1c_8697/readme.txt")
file.copy("../from_vivli/Data/agesexmace_6115/00_readme.txt",
          "Data/vivli_mace/6115_00_readme.txt")
file.copy("../from_vivli/Data/agesexmace_8697/00_readme.txt",
          "Data/vivli_mace/8697_readme.txt")
file.copy("../from_vivli/Data/agesexmacecentred_6115/00_readme.txt",
          "Data/vivli_mace/6115_centred_readme.txt", overwrite = TRUE)
file.copy("../from_vivli/Data/agesexmacecentred_8697/00_readme.txt",
          "Data/vivli_mace/8697_centred_readme.txt", overwrite = TRUE)
