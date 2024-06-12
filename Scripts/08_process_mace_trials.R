#10_mace_trials

library(tidyverse)
library(multinma)
source("Scripts/common_functions/Scripts/convert_iqr_to_sd.R")
source("Scripts/common_functions/Scripts/combine_sd.R")
source("Scripts/common_functions/Scripts/truncated_normal.R")
## read in age criteria ----
age_elig <- read_csv("Data/cleaned_data/Data/age_max_min_elig.csv")
## assume no limits for unusual trial with minimal information (dropped later anyway as no outcomes)
age_elig <- bind_rows(age_elig,
                      tibble(nct_id ="UMIN000018395", min_age = 10L, max_age = 150L, age_unit = "Years"))
## read in association between median and mean fu (in days) ----
med_mean <- read_csv("Outputs/association_between_median_mean_fu.csv")

## read in arm data ----
arm_meta <- read_csv("Data/cleaned_data/Data/arm_data_all_cleaned.csv") %>% 
  rename(nct_id = trial_id)

## read in IPD regression outputs ----
ipd_meta <- read_csv("Data/vivli_mace/Overview_MACE_trials.csv")
ipd <- read_csv("Data/vivli_mace/model_coefficients.csv")
ipd %>% group_by(repo, nct_id) %>% summarise(n = length(term), arms = sum(!duplicated(term[str_detect(term, "^arm")])))

## read in IPD arms mapped across to ctg labelled in database ----
ipd_arms <- read_csv("nct_id,arm,arm_id
NCT00968708,alogliptin,aa0696|aa0697
NCT00968708,placebo,aa0698|aa0699
NCT01032629,jnj-28431754-100 mg,aa0791
NCT01032629,jnj-28431754-300 mg,aa0792
NCT01032629,placebo,aa0793
NCT01131676,bi 10773 10mg,aa0972
NCT01131676,bi 10773 25mg,aa0973
NCT01131676,placebo,aa0974
NCT01989754,canagliflozin,aa1604
NCT01989754,placebo,aa1605|aa1606
NCT02065791,cana 100 mg,aa1663
NCT02065791,placebo,aa1664
NCT02465515,albiglutide,aa1808|aa1809
NCT02465515,placebo,aa1810")
ipd_arms <- ipd_arms %>% 
  inner_join(arm_meta %>% select(-arm_id_subgroup, -(drug2_name:design_group_id), -(n_rand:source)))

## read in IPD trial-level information ----
ipd_cnsr <- read_csv("Data/vivli_mace/ipd_cnsr.csv")

ipd_cnsr %>% 
  distinct(nct_id, arm) %>% 
  count(nct_id, arm) %>% 
  spread(nct_id, n)

## recover trial follow-up

## read in aggregate-level outcome data and arrange into canonical format ----
maceout <- bind_rows(`19` = read_csv("Data/cleaned_data/Data/mace_2019_results.csv"),
                     `22` = read_csv("Data/cleaned_data/Data/mace_2022_results.csv"),
                     .id = "srch_round")
maceout <- maceout %>% 
  rename(nct_id = trial_id)

## arm ID flipped for NCT01897532. Label is correct but arm id is wrong
maceout <- maceout %>% 
  mutate(arm_id_unq = case_when(
    arm_id_unq == "uaa10962" & arm_description == "placebo" ~ "uaa10963",
    arm_id_unq == "uaa10963" & arm_description == "linagliptin" ~ "uaa10962",
    TRUE ~ arm_id_unq))

## Do not YET drop trials when have IPD - WILL BE down to 18 instead of 24
maceout_ipd <- maceout %>% 
  filter(nct_id %in% ipd$nct_id) %>% 
  distinct(srch_round, nct_id, trial_name, max_follow_up, timepoint_unit)
# maceout <- maceout %>% 
#   filter(!nct_id %in% ipd$nct_id)
maceout_trl <- map_int(maceout, ~ sum(!duplicated(paste0(maceout$nct_id, .x))))
maceout_trl <- maceout_trl[maceout_trl == sum(!duplicated(maceout$nct_id))]
maceout_trl <- maceout[ , names(maceout_trl)] %>% 
  distinct()
maceout <- maceout [ , c("nct_id", "timepoint_unit", setdiff(names(maceout), names(maceout_trl)))]
maceout_outcome <- maceout %>% 
  select(nct_id, outcome, outcome_definition) %>% 
  group_by(nct_id, outcome) %>% 
  summarise(outcome_definition = paste(outcome_definition[!is.na(outcome_definition)], collapse = ";")) %>% 
  ungroup()
maceout_cmpr <- maceout %>% 
  filter(!is.na(comparison_id)) %>% 
  select(-arm_id_unq, -arm_description, -timepoint_unit)
## all have these. No multiples
maceout_cmpr <- maceout_cmpr  %>% 
  filter(outcome %in% c("MACE occurrence", 
                        "Percentage of participants with MACE outcome", 
                        "Time to first occurrence of MACE"))
maceout <- maceout %>% 
  filter(!is.na(arm_id_unq))
maceout %>% count(result_type, nct_id)  %>% spread(result_type, n)
maceout %>% count(outcome, nct_id) %>% spread(outcome, n) 
## all mace have count, % or rate. None have multiples of these
maceout <- maceout %>% 
  filter(outcome %in% c("MACE occurrence", 
                        "Percentage of participants with MACE outcome", 
                        "Time to first occurrence of MACE"))
maceout <- maceout %>% 
  select(nct_id, result_id, arm_id_unq, arm_description, outcome, outcome_definition, result_description, result_type, result, dispersion_type,
         dispersion, median_follow_up, timepoint_unit)
maceout <- maceout %>% 
  mutate(median_fu_days = case_when(
    timepoint_unit == "years" ~ median_follow_up * 365.25,
    timepoint_unit == "months" ~ median_follow_up * 30,
    timepoint_unit == "weeks" ~ median_follow_up * 7)) %>% 
  select(-median_follow_up)
rm(maceout_outcome)

## read in baseline data for aggregate trials ----
base_dsp <- read_csv("Data/cleaned_data/base_dsp.csv") %>% 
  rename(nct_id = trial_id) %>% 
  semi_join(maceout) %>% 
  filter(variable %in% c("age", "male", "n"))
base_rng <- read_csv("Data/cleaned_data/base_rng.csv") %>% 
  rename(nct_id = trial_id)%>% 
  semi_join(maceout) %>% 
  filter(variable %in% c("age"))
base_rng <- base_rng %>% 
  inner_join(maceout %>% select(arm_id = arm_id_unq, arm_description) %>% distinct())
## pull mean age and STD from aact if present
aact_age <- read_csv("Data/aact/age.csv") %>% 
  semi_join(base_rng %>% 
              select(nct_id)) 
arm_lkp <- read_csv("nct_id,ctgov_group_code,title,arm_id
NCT03315143,BG000,Sotagliflozin,unq_updaa10053
NCT03315143,BG001,Placebo,unq_updaa10052
NCT03521934,BG000,Sotagliflozin,unq_updaa10073
NCT03521934,BG001,Placebo,unq_updaa10072
NCT01144338,BG000,Placebo,uaa10618
NCT01144338,BG001,Exenatide Once Weekly,uaa10617")

aact_age <- aact_age %>% 
  inner_join(arm_lkp)
base_rng <- base_rng %>% 
  filter(!nct_id %in% aact_age$nct_id)
base_rng <- base_rng %>% 
  inner_join(base_dsp %>% filter(variable == "n") %>% 
               select(nct_id, arm_id, n = first))
base_rng <- base_rng %>% 
  ## confusing naming here because "first" refers to first quartile
  # in function but first value in baseline table
  mutate(age_m = EstMean2(first = low, medn = first, third = upp),
         age_s = EstSD2(first = low, third = upp, smpl = n))

## Harmonise age data for aggregate trials ----
age_dsp <- base_dsp %>% 
  filter(variable == "age") %>% 
  select(unq_id:id, age_m = first, age_s = second) %>% 
  inner_join(base_dsp %>%
               filter(variable == "n") %>%
               select(nct_id, arm_id, participants = first))
age_aact <- aact_age %>% 
  select(nct_id, arm_id, id = result_group_id, age_m = param_value_num, age_s = dispersion_value_num, participants = number_analyzed)
age_cnvrt <- base_rng %>% 
  select(nct_id, arm_id, id, age_m, age_s, participants = n)
age_agg <- bind_rows(age_dsp = age_dsp %>% mutate(id = as.character(id)),
                     age_aact = age_aact %>% mutate(id = as.character(id)),
                     age_cnvrt = age_cnvrt, .id = "age_source")
rm(aact_age, age_aact, age_cnvrt, age_dsp, arm_lkp)

## Pull sex data for aggregate trials ----
male <- base_dsp %>% 
  filter(variable == "male") %>% 
  select(nct_id, arm_id, male_id = id, male_prcnt = second) 
base_agg <- age_agg %>% 
  inner_join(male)
rm(age_agg, male)

# Combine baseline and outcome data into a wide format ----
agg <- bind_rows(base_agg %>% distinct(nct_id, arm_id) %>% mutate(type = "base"),
                 maceout %>% select(nct_id, arm_id = arm_id_unq, arm_description) 
                 %>% distinct() %>% mutate(type = "outcome")) %>% 
  mutate(v = 1L) %>% 
  group_by(arm_id) %>% 
  mutate(arm_description = arm_description[!is.na(arm_description)][1]) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  left_join(arm_meta %>% select(-arm_id) %>% rename(arm_id = arm_id_unq) %>% 
              select(nct_id, arm_id, arm_label:drug_freq)) %>% 
  arrange(nct_id, arm_id)
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(msmtch = any(is.na(base) | any(is.na(outcome)))) %>% 
  ungroup()
## 3 trials, resolve mismatch appear to be different splits of data in outcome and base ----
agg_msmtch <- agg %>% 
  filter(msmtch)
## Review CTG directly
# NCT01720446 - subgroup uaa10891a and uaa10891b are both placebo arms with different dummy pills. Collapse base to
# NCT01720446 - subgroup uaa10892 and uaa10894 are both semaglutide arms with different doses. Collapse base to uaa10893
base_agg_simplify <- base_agg %>% 
  filter(nct_id == "NCT01720446") %>% 
  mutate(arm_id = case_when(
    arm_id %in% c("uaa10891a", "uaa10891b") ~ "uaa10891",
    arm_id %in% c("uaa10892", "uaa10894") ~ "uaa10893",
    TRUE ~ arm_id
  )) %>% 
  group_by(nct_id, arm_id) %>% 
  summarise(age_s = CombSdVectorised(participants, age_m, age_s),
            age_m = weighted.mean(age_m, participants),
            male_prcnt = weighted.mean(male_prcnt, participants),
            participants = sum(participants)) %>% 
  ungroup()
base_agg <- base_agg %>% 
  filter(!nct_id %in% "NCT01720446") %>% 
  bind_rows(base_agg_simplify)

# NCT01986881 - CTG has 3 arms in baseline Ertugliflozin 5, Ertugliflozin 15 and placebo and the same
# 3 arms in the outcome data, but with a collapsed arm for Ertugliflozin. Plan. Collapse baseline Ertugliflozin. 
# by reading in the following data. Note collapse rows 2 (5mg) and 3 (15mg)
NCT01986881_base <- read_csv("arm_description,arm_id,age_m,age_s,male_prcnt,participants
placebo,uaa10992,64.4,8.0,69.3,2747
ertugliflozin,uaa77777,64.3,8.2,70.9,2752
ertugliflozin,uaa77777,64.4,8.0,69.7,2747")
NCT01986881_base <- NCT01986881_base %>% 
  group_by(arm_description, arm_id) %>% 
  summarise(age_s = CombSdVectorised(participants, age_m, age_s),
            age_m = weighted.mean(age_m, participants),
            male_prcnt = weighted.mean(male_prcnt, participants),
            participants = sum(participants)) %>% 
  ungroup()
base_agg <- base_agg %>% 
  filter(!nct_id == "NCT01986881") %>% 
  bind_rows(NCT01986881_base %>% mutate(nct_id = "NCT01986881"))
## NCT03496298  can't figure out from data. Pull again from aact for BOTH outcome and base
NCT03496298 <- readRDS("Data/aact/NCT03496298.Rds")
male <- NCT03496298$baseline_measurements %>% 
  filter(title == "Sex: Female, Male") %>% 
  select(id, nct_id, result_group_id, ctgov_group_code, category, param_value_num) 
names(male) <- str_to_lower(names(male))
age <- NCT03496298$baseline_measurements %>% 
  filter(title == "Age, Continuous") %>% 
  select(id, nct_id, result_group_id, ctgov_group_code, age_m = param_value_num, age_s = dispersion_value_num, participants = number_analyzed)
out <- NCT03496298$outcome_measurements %>% 
  filter(title == "Time to First Occurrence of Major Adverse Cardiovascular Events: Event Rate Per 100 Participant-years for First Occurrence of Major Cardiovascular Event - Superiority Analysis")
out_counts <- NCT03496298$outcome_counts
# note none dropped in this join
out <- out %>% 
  inner_join(out_counts %>% select(-id, -units)) %>% 
  select(count, everything())
rg <- NCT03496298$result_groups %>% 
  select(nct_id, ctgov_group_code, title)
age <- age %>% 
  inner_join(rg)
male <- male %>% 
  inner_join(rg)
out <- out %>% 
  rename(title_outcome = title) %>% 
  inner_join(rg) %>% 
  select(nct_id, title, title_outcome, rate = param_value_num, participants = count)
male <- male %>% 
  filter(category == "Male") %>% 
  mutate(title = if_else(title %in% c("Efpeglenatide 4 mg", "Efpeglenatide 6 mg"), "Efpeglenatide 4 mg+6 mg", title)) %>% 
  group_by(nct_id, title) %>% 
  summarise(male = sum(param_value_num)) %>% 
  ungroup()
age <- age %>% 
  mutate(title = if_else(title %in% c("Efpeglenatide 4 mg", "Efpeglenatide 6 mg"), "Efpeglenatide 4 mg+6 mg", title)) %>% 
  group_by(nct_id, title) %>% 
  summarise(age_s = CombSdVectorised(participants, age_m, age_s),
            age_m = weighted.mean(age_m, participants),
            participants = sum(participants)) %>% 
  ungroup()
NCT03496298 <- age %>% 
  inner_join(male) %>% 
  inner_join(out)

## median FU is 19.9 and 19.8 months for tx and placebo respectively
NCT03496298 <- NCT03496298 %>% 
  mutate(median_fu_days = if_else(title == "Placebo", 19.8*30, 19.9*30),
         result_type = "rate",
         result_description = "events per 100 patient-years")
rm(age, agg, agg_msmtch, base_agg_simplify, NCT01986881_base, out, out_counts, rg)

base_out <- base_agg %>% 
  filter(!nct_id == "NCT03496298")
maceout <- maceout %>% 
  filter(!nct_id == "NCT03496298")
## all appear to match apart from NCT01131676. This is because the subgroup data has two arms
## and the baseline data is individual arms. Create an aggregate one for merging
base_out2 <- base_out %>% 
  filter(nct_id == "NCT01131676", arm_id %in% c("uaa10604", "uaa10606")) %>% 
  group_by(age_source, nct_id) %>% 
  summarise(age_s = CombSdVectorised(participants, age_m, age_s),
            age_m = weighted.mean(age_m, participants),
            male_prcnt = weighted.mean(male_prcnt, participants),
            participants = sum(participants),
            across(c(id, male_id, unq_id, arm_description), ~ paste(.x, collapse = ";"))) %>% 
  ungroup() %>% 
  mutate(arm_id = "uaa10605")
base_out <- bind_rows(base_out %>% 
                        filter(!arm_id %in% c("uaa10604", "uaa10606")),
                      base_out2)

rm(base_out2)
base_out %>% 
  anti_join(maceout %>% rename(arm_id = arm_id_unq) %>% select(-arm_description))

base_out <- base_out %>% 
  select(-arm_description) %>% 
  inner_join(maceout %>% rename(arm_id = arm_id_unq))
NCT03496298 <- NCT03496298 %>% 
  rename(outcome = title_outcome,
         arm_description = title,
         result = rate) %>% 
  mutate(male_prcnt = 100*male/participants,
         result_type = "rate",
         arm_id = paste0("DVD", 1:2)) %>% 
  select(-male)
base_out <- base_out[, names(NCT03496298)] %>% 
  bind_rows(NCT03496298)
rm(NCT03496298)
rm(base_agg, base_dsp, base_rng, maceout, male)
base_out <- base_out %>% 
  mutate(arm_description = if_else(arm_description == "Efpeglenatide 4 mg+6 mg", "efpeglenatide", arm_description),
         arm_description = str_to_lower(arm_description))
# base_out %>% count(arm_description, nct_id) %>% spread(nct_id, n)  %>% View()
## Note no arms with same drug, all 18 are placebo controlled, all are two arm trials (in this analysis)
base_out %>% count(arm_description, nct_id) %>% filter(n >2)

base_out <- base_out %>% 
  mutate(bythis = paste(nct_id, arm_id, sep = "__"))
base_out <- base_out %>% 
  inner_join(age_elig %>% select(-age_unit))
siml_estimates <- by(base_out %>% rename(age_sd = age_s), base_out$bythis, EstimateMuDispAge)
## Extract and convert to list
siml_estimates <- map(siml_estimates, identity)
siml_estimates <- map(siml_estimates, ~ .x %>% mutate(pick = seq_along(trial_mean)))
siml_estimates <- bind_rows(siml_estimates, .id = "bythis")
siml_estimates <- siml_estimates %>% 
  filter(pick ==1)
siml_estimates <- siml_estimates %>% 
  select(bythis, mu_x, sd_x) 

## update simulation
SafeMuSigmaCalc <- safely(MuSigmaCalc)

base_out$res <- pmap( list(base_out$min_age, base_out$max_age, base_out$age_m, base_out$age_s),
                          function(min_age, max_age, age_m, age_sd) {
                            # browser()
                            ## If no age limits return mean and SD
                            if(max_age >=100 & min_age <= 20) {
                              b <- tibble(mu = age_m, 
                                          sigma = age_sd, 
                                          convergence = 999, 
                                          value = 0,
                                          method = "ignore_limits")
                              return(b)  
                            } 
                            ## If age limits first try L-BFGS note this sometimes fails so need to run safely 
                            a <- SafeMuSigmaCalc(min_age, max_age, age_m, age_sd, mymethod = "L-BFGS-B", sigma_ranges = c(5, 20), mu_range = c(40,70))
                            if(is.null(a$error)) {
                              ## double if statement here as this object is only present if it runs
                              if(a$result$convergence[[1]] == 0) {
                                a <- a$result
                                b <- tibble(mu = a$par[1],
                                            sigma = a$par[2],
                                            convergence = a$convergence[[1]],
                                            value = a$value[[1]],
                                            method = "L-BFGS-B")
                                return(b) 
                              } ## end condition only if converges
                            } ## end condition only if not null
                            ## If L-BFGS-B fails or does not converge run Nelder-Mead Nelder-Mead
                            a <- SafeMuSigmaCalc(min_age, max_age, age_m, age_sd, mymethod = "CG")
                            if(is.null(a$error)) {
                              a <- a$result
                              b <- tibble(mu = a$par[1],
                                          sigma = a$par[2],
                                          convergence = a$convergence[[1]],
                                          value = a$value[[1]],
                                          method = "CG")
                              return(b) 
                            } 
                          })
if(any(map_lgl(base_out$res, is_null))) warning ("Cycled through 3 algorithms without running succesfully")
base_out <- base_out %>% 
  unnest(res)
if(any(!base_out$convergence %in% c(0, 999))) warning ("Optim did not converge for at least one row")
## very similar results with both approaches. Use optim preferentially
base_out <- base_out %>% 
  inner_join(siml_estimates) %>% 
  select(-bythis)  
## Compare grid approach with one of the optimisation algorithms
base_out %>%
  select(mu, mu_x) %>% 
  mutate(mu_diff = abs(mu - mu_x)) %>% 
  arrange(desc(mu_diff))
base_out %>% 
  select(sigma, sd_x) %>% 
  mutate(sigma_diff = abs(sigma - sd_x)) %>%
  arrange(desc(sigma_diff))

base_out <- base_out %>% 
  rename(age_mu = mu,
          age_sigma = sigma) %>% 
  select(-mu_x, -sd_x)
## count aggregate and IPD trials ----
mace <- base_out
rm(base_out)
mace <- mace %>% 
  rename(drug_name = arm_description)
mace %>% count(drug_name, nct_id) %>% spread(nct_id, n)
agg_arms <- mace %>% 
  select(nct_id, arm_id, drug_name_current = drug_name) %>% 
  left_join(arm_meta %>% 
              select(-arm_id) %>% 
              rename(arm_id = arm_id_unq) %>% 
              select(-arm_id_subgroup, -(drug2_name:design_group_id), -(n_rand:source))) 
## one trial where arm IDs appear flipped.
## check in both sources against CTG
## outcome and baseline data
mace %>% 
  filter(nct_id == "NCT01897532") %>% 
  select(arm_id, drug_name, everything())
## arm meta data
arm_meta %>% 
  filter(nct_id == "NCT01897532")
agg_arms %>% filter(!drug_name_current == arm_label)

agg_arms <- agg_arms %>% 
  mutate(drug_name = case_when(
    nct_id == "NCT03496298" & drug_name_current == "efpeglenatide" ~ "Efpeglenatide 4 mg+6 mg",
    nct_id == "NCT03496298" & drug_name_current == "placebo" ~ "placebo",
    drug_name_current == "placebo" ~ "placebo",
    TRUE ~ drug_name
  ))

bth_arm <- bind_rows(agg = agg_arms,
                     ipd = ipd_arms, .id = "data_lvl") %>% 
  select(-drug_name_current)
## For semaglutide note BNF statement "Oral semaglutide 14 mg once daily is comparable to subcutaneous semaglutide 0.5 mg once weekly"
bth_arm %>% 
  filter(drug_name %in% c("canagliflozin", "empagliflozin", "semaglutide")) %>% 
  arrange(nct_id, drug_name)

bth_arm <- bth_arm %>% 
  mutate(arm_lvl = 
           case_when(
             drug_name == "Efpeglenatide 4 mg+6 mg" ~ "efpeglenatide_4_6",
             drug_name == "ITCA_650" ~ "ITCA",
             drug_name == "semaglutide" & drug_dose == "0.5|1" ~ paste0(drug_name, "_", 14), 
             !drug_name == "placebo" ~ paste0(drug_name, 
                                              "_",
                                              drug_dose %>% 
                                                str_replace("\\|", "_") %>% 
                                                str_replace("\\.", "p") %>% 
                                                str_remove("pooled") %>% 
                                                str_trim()),
             drug_name == "placebo" ~ "placebo")) 
bth_arm <- bth_arm %>% 
  mutate(drug_dose = if_else(
           arm_lvl == "efpeglenatide_4_6", "4_6", drug_dose),
         drug_name = if_else(
           arm_lvl == "efpeglenatide_4_6", "efpeglenatide", drug_name))

# new arm labels in WHO database and own manual lookup ----
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
whoatc <- bind_rows(whoatc,
                    nonwho) %>% 
  distinct()
whoatc_lkp <- whoatc$atc_code
names(whoatc_lkp) <- whoatc$nm
setdiff(str_to_lower(bth_arm$drug_name), names(whoatc_lkp))
bth_arm <- bth_arm %>% 
  mutate(atc = whoatc_lkp[drug_name %>% str_to_lower()])
bth_arm <- bth_arm %>% 
  mutate(trtcls5 = str_sub(atc, 1, 5))
bth_arm <- bth_arm %>% 
  left_join(whoatc %>% 
               distinct(atc_code, .keep_all = TRUE) %>% 
               rename(trtcls5 = atc_code)) %>% 
  mutate(dc = if_else(trtcls5 == "place", 
                      "placebo",
                      paste(trtcls5, nm, sep = " "))) %>% 
  distinct()

## Estimate mean fu from median ----
## Prior to calculating reviewed NCT03521934 as an order of magnitude higher than every other trial and % more than 100% in placebo arm
## the following were the rates and follow-up times given on CTG but according to paper (https://www.nejm.org/doi/full/10.1056/NEJMoa2030183) r is 
## 247 in the sotagliflozin group and 330 in the placebo group as 
## row " per Deaths from cardiovascular causes, hospitalizations for heart failure, nonfatal myocardial infarctions, and nonfatal strokes — total no. of events (rate)†"
## in table 2 and "In the sotagliflozin group, the median duration of follow-up was 9.2 months, the median duration of treatment was 7.8 months ...
## the corresponding values in the placebo group were 8.9 months, 7.6 months" set to this in person years and use these events
## the reason the % comes out higher is that the follow-up in CTG is NOT the median follow-up it is the maximum follow-up 
## and there is a large difference between the two here as the trial stopped early
py_sota <- (608*9.2/12) 
py_plac <- (614*8.9/12)
rate_sota <- 100*247/py_sota
rate_plac <- 100*330/py_sota
## correct manually
NCT03521934 <- mace %>% 
  filter(nct_id == "NCT03521934")

NCT03521934 <- NCT03521934 %>% 
  mutate(r = case_when(
    drug_name == "sotagliflozin" ~ 245,
    drug_name == "placebo" ~ 355),
    median_fu_days = case_when(
      drug_name == "sotagliflozin" ~ 9.2*30,
      drug_name == "placebo" ~ 8.9*30))
mace <- mace %>% 
  filter(!nct_id == "NCT03521934") %>% 
  bind_rows(NCT03521934)
rm(NCT03521934)
mace <- mace %>% 
  mutate(mean_fu_days = med_mean$estimate[med_mean$term == "(Intercept)"] +
           med_mean$estimate[med_mean$term == "fu_median"] * median_fu_days)


## Read in subgroup data for age and sex ----
## Note. In IPD similar age between men and women (around 2 years different) so treat as independent in subgroup data
## ie assume mean age is same in both sexes. assume percent male is same in all age categories
mace_sg <- read_csv("Data/cleaned_data/Data/mace_subgroups_age_sex.csv")
## drop IPD ones
# mace_sg <- mace_sg %>% 
#   filter(!nct_id %in% ipd$nct_id)
# split into main, age and sex
mace_main <- mace_sg %>% 
  filter(subgroup == "main")
mace_age <- mace_sg %>% 
  filter(subgroup == "age")
mace_sex <- mace_sg %>% 
  filter(subgroup == "sex")

## merge on nct_id and treatment arm (not same arm_unq_id)
## resolve age. Drop total as no useful information
mace_age <- mace_age %>% 
  filter(!arm_label == "total") %>% 
  select(nct_id, arm_id_unq_ph = arm_id_unq, 
         arm_label, outcome, subgroup, level_min, level_max, n, hr, hr_lci, hr_uci) 
setdiff(mace_age$arm_label, mace$drug_name)
## one trial has two different cut-points. one at 65 and one at 75. Choose the latter
mace_age <- mace_age %>% 
  filter(!(nct_id == "NCT01144338" & (level_max == 64 | level_min == 65) ))
## One trial has additional rows without information. Drop
mace_age <- mace_age %>% 
  filter(!(nct_id == "NCT03496298" & is.na(n)))
## One trial has additional rows now collapsed Drop
mace_age <- mace_age %>%  
  filter(!(nct_id == "NCT01986881" & arm_id_unq_ph %in% c("uaa10990", "uaa10991")),
         !(nct_id == "NCT03496298" & arm_id_unq_ph == "unq_updaa10067" & is.na(n)))
mace_age %>% 
  filter(nct_id == "NCT01243424")
mace_age <- mace_age %>% 
  rename(drug_name = arm_label) %>% 
  inner_join(mace %>% 
               select(nct_id, drug_name, arm_id, male_prcnt, 
                      age_m, age_s, age_mu, age_sigma, max_age, min_age, n_overall = participants))
## of 10 trials with age subgroups, only 2 have a maximum age in the eligibility criteria
mace_age <- mace_age %>% 
  mutate(level_min = if_else(level_min == "min", min_age, as.double(level_min)),
         level_max = if_else(level_max == "max", max_age, as.double(level_max)))
## calculate number in each age group where this is missing
mace_age <- mace_age %>% 
  mutate(pfrom = truncnorm::ptruncnorm(level_min-1, a = min_age, b = max_age, age_mu, age_sigma),
         pto   = truncnorm::ptruncnorm(level_max, a = min_age, b = max_age, age_mu, age_sigma),
         pin = pto - pfrom,
         n_impute = round(n_overall * pin)) %>% 
  select(-pfrom, -pto)
## all sum to one
mace_age %>% 
  group_by(nct_id, arm_id_unq_ph) %>% 
  summarise(pin = sum(pin)) %>% 
  ungroup() %>% 
  filter(pin != 1)
## all imputed n's reasonably close to real N's (where these are not missing)
mace_age <- mace_age %>% 
  mutate(n = if_else(is.na(n), n_impute, n)) %>% 
  select(-n_impute)
mace_sex <- mace_sex %>% 
  filter(!arm_label == "total") %>% 
  select(nct_id, arm_id_unq_ph = arm_id_unq, arm_label, outcome, subgroup, 
         level_cat, n, hr, hr_lci, hr_uci) 
mace_sex <- mace_sex %>% 
  filter(!(nct_id == "NCT01986881" & arm_id_unq_ph %in% c("uaa10990", "uaa10991")),
         !(nct_id == "NCT03496298" & arm_id_unq_ph == "unq_updaa10067" & is.na(n)))
setdiff(mace_sex$arm_label, mace$drug_name)
mace_sex <- mace_sex %>% 
  rename(drug_name = arm_label) %>% 
  inner_join(mace %>% select(nct_id, drug_name, arm_id, male_prcnt, n_overall = participants,
                             age_m, age_s, age_mu, age_sigma, max_age, min_age)) 
mace_sex <- mace_sex %>% 
  mutate(male_prpn = male_prcnt/100,
         n = case_when(
    !is.na(n) ~ n,
    level_cat == "male" ~ n_overall*male_prpn,
    level_cat == "female" ~  n_overall*(1-male_prpn)),
    n = round(n)) %>% 
  select(-n_overall, -male_prpn) %>% 
  mutate(male_prcnt = if_else(level_cat == "male", 100, 0))
rm(mace_sg)

## estimate events from mean FU ----
## Note not appropriate to correct for events here as the mean is the mean time in the arm in those with and without events
mace <- mace %>% 
  mutate(pt = mean_fu_days * participants,
         r = case_when(
           result_type == "count" ~ result,
           result_type == "percentage" ~ result * participants/100,
           result_description == "events per 100 patient-years" ~ result * pt/(100*365),
           result_description == "events per 1000 patient-years" ~ result * pt/(1000*365)),
         r = round(r))
## note one trial has zero events in both arms, I think drop for all subsequent analyses
## review all trials
# select(nct_id, drug_name, outcome, participants, result, result_type, result_description, median_fu_days, mean_fu_days, pt, r, percentage_check) %>% 
mace <- mace %>% 
  mutate(percentage_check = round(100* r/participants,1),
         rate_check = round(100*365 *r/pt,2))
mace %>% 
  select(nct_id, drug_name, outcome, participants, result, result_type, result_description, median_fu_days, mean_fu_days, pt, r, 
         percentage_check, rate_check) 
## One trial stands out considerably, the one that was stopped early
mace_plot_check_data <- mace %>% 
  filter(!nct_id %in% c("UMIN000018395", "NCT01131676")) %>% 
  group_by(nct_id) %>% 
  mutate(arm = c("a", "b")) %>% 
  ungroup()
mace %>% 
  filter(nct_id %in% c("UMIN000018395", "NCT01131676"))

plot_rate_check <- ggplot(mace_plot_check_data %>% 
                      mutate(mydiff = nct_id == "NCT03521934"), aes(x = nct_id, y = rate_check, colour = arm, shape = arm)) +
  geom_point(position = position_dodge(width = 0.1)) +
  coord_flip()+
  facet_wrap(~mydiff, scales = "free")
plot_rate_check
plot_percentage_check <- ggplot(mace_plot_check_data %>% 
                            mutate(mydiff = nct_id == "NCT03521934"), aes(x = nct_id, y = percentage_check, colour = arm, shape = arm)) +
  geom_point(position = position_dodge(width = 0.1)) +
  coord_flip() +
  facet_wrap(~mydiff, scales = "free")
plot_percentage_check

## note all aggregate data trials are two arm analyses (some arms collapsed before we extracted data)
mace <- mace %>% 
  select(nct_id, drug_name, participants, age_m, age_s, male_prcnt, outcome, r, pt, mean_fu_days,
         min_age, max_age, age_mu, age_sigma)
## drop arm label as no additional information
bth_arm <- bth_arm %>% 
  mutate(drug_name = if_else(drug_name == "ITCA_650", "itca650", drug_name)) %>% 
  select(-arm_label)
bth_arm <- bth_arm %>% 
  rename(ipd_arm = arm)

## all arm_id_unq are NA for agg but present for IPD
mace2 <- mace %>% 
  inner_join(bth_arm %>% 
    select(nct_id, drug_name, arm_id, drug_dose, drug_unit, drug_freq, arm_lvl, atc, trtcls5, dc, note))


## add in treatment comparison data from aggregate ----
maceout_cmpr2 <- maceout_cmpr %>% 
  mutate(comparison_label = if_else(comparison_label == "saxaliptin_placebo", "saxagliptin_placebo", comparison_label)) %>% 
  separate(comparison_label, into = c("treat", "control"), sep = "_") %>% 
  select(nct_id, treat, control, hr = result, ci = dispersion) %>% 
  mutate(ci = if_else(ci == "0.57:1.11", "0.57;1.11", ci)) %>% 
  separate(ci, into = c("l", "u"), sep = ";", remove = FALSE) %>% 
  mutate(across(c(l, u), as.double),
         across(c(hr, l, u), log),
         se = (u-l)/(2*1.96)) %>% 
  select(nct_id, treat, control, loghr = hr, se = se) %>% 
  gather("treat_cmpr", "drug_name", treat, control) %>% 
  mutate(across(c(loghr, se), ~ if_else(treat_cmpr == "control", NA_real_, .x))) %>% 
  arrange(nct_id, drug_name) %>% 
  select(nct_id, drug_name, treat_cmpr, loghr, se)
maceout_cmpr2_add <- mace2 %>% 
  anti_join(maceout_cmpr2) %>%
  mutate(treat_cmpr = if_else(drug_name == "luseogliflozin", "treat", "control")) %>% 
  select(nct_id, drug_name, treat_cmpr)
maceout_cmpr2 <- bind_rows(maceout_cmpr2,
                           maceout_cmpr2_add)
## note none not joining
mace3 <- mace2 %>% 
  inner_join(maceout_cmpr2) %>% 
  select(nct_id:mean_fu_days, loghr, se, everything())

## conver hr to log-scale for subgroup data----
## when add IPD trials back in needs distinct rather than select statement. Not an issue for non-IPD aggregate subgroup trials (checked)
mace_age2 <- mace_age %>% 
  inner_join(mace3 %>% distinct(nct_id, drug_name, trtcls5))
mace_sex2 <- mace_sex %>% 
  inner_join(mace3 %>% distinct(nct_id, drug_name, trtcls5))

mace_age2 <- mace_age2  %>% 
  mutate(across(c(hr_lci, hr_uci ), as.double),
         across(c(hr, hr_lci, hr_uci), log),
         se = (hr_uci-hr_lci)/(2*1.96)) %>% 
  rename(loghr = hr, se = se) %>% 
  select(-hr_lci, -hr_uci)
mace_sex2 <- mace_sex2  %>% 
  mutate(across(c(hr_lci, hr_uci ), as.double),
         across(c(hr, hr_lci, hr_uci), log),
         se = (hr_uci-hr_lci)/(2*1.96)) %>% 
  rename(loghr = hr, se = se) %>% 
  select(-hr_lci, -hr_uci)

## review where have aggregate and IPD
# maceout_cmpr_ipd
ipd_cmpr_agg <- ipd %>% 
  filter(models == "f1") %>% 
  select(repo, nct_id, term, f1_est = estimate, f1_se = std.error) %>% 
  mutate(term = str_sub(term, 6)) %>% 
  inner_join(ipd_arms %>% select(nct_id, term = arm, drug_name, drug_dose))
mace3_cmpr_ipd_agg <- maceout_cmpr2 %>% 
  filter(nct_id %in% ipd_cmpr_agg$nct_id) %>% 
  select(nct_id, drug_name, treat_cmpr, loghr, se) %>% 
  distinct()
## allowing for collapsing of trial arms in aggregate data, results HIGLY CONSISTENT BETWEEN IPD and aggregate data
ipd_cmpr_agg <- ipd_cmpr_agg %>% arrange(nct_id) %>% 
  left_join( mace3_cmpr_ipd_agg %>% filter(!treat_cmpr == "control")) %>% 
  mutate(note = "Some arms collapsed for aggregate hence identical") %>% 
  rename(agg_est = loghr,
         agg_se = se)
write_csv(ipd_cmpr_agg, "Outputs/similar_results_ipd_agg_hrs.csv")

### save for subsequent analysis ----
## note there is some redundancy here as have merged some variables from mace_arms into the aggregate
## data for mace. But did so for clarity as there are multiple arm IDs
mace4 <- mace3 %>% 
  filter(! (nct_id == "NCT01131676" & arm_lvl %in% c("empagliflozin_10", "empagliflozin_25"))) %>% 
  distinct(nct_id, drug_name, arm_lvl)
lst <- list(mace_arms = bth_arm,
            mace_agg = mace3,
            mace_agg_age = mace_age2 %>% 
              rename(participants = n) %>% 
              inner_join(mace4),
            mace_agg_sex = mace_sex2 %>% 
              rename(participants = n) %>% 
              inner_join(mace4),
            mace_agg_trial_level = maceout_trl,
            mace_ipd_trial_level = maceout_ipd)
## drop trials with IPD from aggregate data
map(lst, ~ intersect(.x$nct_id, ipd$nct_id))
lst$mace_agg <- lst$mace_agg %>% 
  filter(!nct_id %in% ipd$nct_id)
lst$mace_agg_trial_level <- lst$mace_agg_trial_level %>% 
  filter(!nct_id %in% ipd$nct_id)
## The following drops aggregate arms from metadata where we already have IPD as these are duplicates
ipdandagg <- lst$mace_arms %>% 
  distinct(nct_id, data_lvl) %>% 
  arrange(desc(data_lvl)) %>% 
  distinct(nct_id, .keep_all = TRUE)
lst$mace_arms <- lst$mace_arms  %>% 
  semi_join(ipdandagg)
## checked no change in exiting subgroup data when update code to include two trials with SG data that also have IPD
## note only runs if existing file
ipdwsggot <- c("NCT00968708", "NCT02465515", "NCT01131676")
if (file.exists("Scratch_data/mace_arms_agg_data.Rds")){
 xmn <- readRDS("Scratch_data/mace_arms_agg_data.Rds")
 identical(lst$mace_agg_age %>% filter(!nct_id %in% ipdwsggot) ,xmn$mace_agg_age)
 identical(lst$mace_agg_sex %>% filter(!nct_id %in% ipdwsggot) ,xmn$mace_agg_sex)
 # separate subgroup data into separate table when there is IPD
}

lst$mace_agg_age_sens <- lst$mace_agg_age 
lst$mace_agg_age <- lst$mace_agg_age %>% 
  filter(!nct_id %in% ipdwsggot)
lst$mace_agg_sex_sens <- lst$mace_agg_sex 
lst$mace_agg_sex <- lst$mace_agg_sex %>% 
  filter(!nct_id %in% ipdwsggot)
if (file.exists("Scratch_data/mace_arms_agg_data.Rds")){
identical(lst$mace_agg_age, xmn$mace_agg_age)
identical(lst$mace_agg_sex, xmn$mace_agg_sex)
}

## set HR in placebo to NA rather than 1 (so consistent with other trials)
lst$mace_agg_sex_sens <- lst$mace_agg_sex_sens  %>% 
  mutate(loghr = if_else(nct_id == "NCT01131676" & drug_name == "placebo", NA_real_, loghr))
lst$mace_agg_age_sens <- lst$mace_agg_age_sens  %>% 
  mutate(loghr = if_else(nct_id == "NCT01131676" & drug_name == "placebo", NA_real_, loghr))

saveRDS(lst, "Scratch_data/mace_arms_agg_data.Rds")

