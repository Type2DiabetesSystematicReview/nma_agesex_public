#02_process_vivli
library(tidyverse)

## read in exclusions ----
exclusions <- read_csv("Data/exclusions_update.csv")

## read in IDs where have IPD ----
ipd1 <- read_csv("Data/agesexhba1c_6115/hba1c_base_change_overall.csv")
ipd2 <- read.csv("Data/gsk/hba1c_base_change_overall.csv")
ipd3 <- read.csv("Data/agesexhba1c_8697/hba1c_base_change_overall.csv")
warning("Update to local folder")
ipd4 <- read.csv("Data/agesexhba1c_9492/hba1c_base_change_overall.csv")
ipd_nct_id <- bind_rows(ipd1, ipd2, ipd3, ipd4) %>% 
  distinct(nct_id) %>% 
  pull()
newtrials <- setdiff(ipd4$nct_id, c(ipd1$nct_id, ipd2$nct_id, ipd3$nct_id))
saveRDS(newtrials, "Scratch_data/newtrials.Rds")
rm(ipd1, ipd2, ipd3, ipd4)

## read in new reference to rename terms, rows, cols and age_continuous.
newref <- read_csv("Data/cleaned_data/Data/link_new_arm_lbls_to_arm_ids.csv")
newref <- newref %>% 
  mutate(arm_id_unq = if_else(!is.na(arm_id_subgroup), arm_id_subgroup, arm_id_unq)) %>% 
  select(nct_id, trial_lbl = arm_f_cleaned, arm_id_unq, arm_label, dc, arm_type)

## read in vivli agesex results ----
allvivli <- list.files("Data/agesexhba1c_6115/", patt = "csv$")
res1 <- map(allvivli, ~ read_csv(paste0("Data/agesexhba1c_6115/", .x)))
names(res1) <- allvivli %>% str_sub(1, -5)
# modelspec <- c("age_sex_model_coefs",
#                "age_sex_model_diag",
#                "age_sex_model_vcov")
# res1 <- res1[setdiff(names(res1), modelspec)]

# bind categorical and continuous age data together. Only have continuous for gsk
res1$age_distribution_baseline_continuous <- bind_rows(res1$age_distribution_baseline_continuous,
                       res1$age_distribution_baseline_categorical) 

## read in gsk agesex results ----
allgsk <- list.files("Data/gsk/", patt = "csv$")
res2 <- map(allgsk, ~ read_csv(paste0("Data/gsk/", .x)))
names(res2) <- allgsk %>% str_sub(1, -5)
# res2 <- res2[setdiff(names(res2), modelspec)]
res1 <- res1[names(res2)]

## read in vivli agesex results from second vivli repository ----
allvivli <- list.files("Data/agesexhba1c_8697/", patt = "csv$")
res3 <- map(allvivli, ~ read_csv(paste0("Data/agesexhba1c_8697/", .x)))
names(res3) <- allvivli %>% str_sub(1, -5)
# res3 <- res3[setdiff(names(res3), modelspec)]

# drop csv only in second_vivli
res3$reference_trial_arm_to_arm_data_all_cleaned <- NULL
res <- pmap(list(res1, res2, res3), function(x, y, z) bind_rows(vivli1 = x, gsk = y, vivli2 = z))

## read in new gsk results with ae and rcs and bocf ----
allgsk <- list.files("Data/gsk_ae/", patt = "csv$")
allgsk <- list(
  age_sex_model_coefs = allgsk[str_detect(allgsk, "coefs") & !str_detect(allgsk, "_ae")],
  age_sex_model_vcov = allgsk[str_detect(allgsk, "vcov") & !str_detect(allgsk, "_ae")],
  age_sex_model_diag = allgsk[str_detect(allgsk, "diag") & !str_detect(allgsk, "_ae")],
  smry = allgsk[str_detect(allgsk, "smry") & !str_detect(allgsk, "_ae")])
res_gsk <- map(allgsk, ~ map(.x, ~ read_csv(paste0("Data/gsk_ae/", .x))) %>% 
                bind_rows())

## Add new age distribution data
age_distribution_baseline_continuous <- read_csv("Data/agesexhba1c_9492/age_distribution_baseline_continuous.csv") 
age_distribution_baseline_continuous <- age_distribution_baseline_continuous %>% 
  mutate(trial_lbl = str_sub(arm_id_unq, 4)) %>% 
  select(-arm_id_unq) %>% 
  inner_join(newref %>% select(nct_id, trial_lbl, arm_id_unq)) %>% 
  select(nct_id, arm_id_unq, sex, participants:`Age (years) at quantiles (%)`)

res$age_distribution_baseline_continuous <- bind_rows(res$age_distribution_baseline_continuous %>% 
                                                        filter(!nct_id %in% age_distribution_baseline_continuous$nct_id),
                                                      age_distribution_baseline_continuous)
rm(age_distribution_baseline_continuous)

## add in new data from 9492 including rcs with fixed knots and empareg ----
# first export, locf, bocf
age_sex_model_coefs_9492_new1 <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv") %>% 
  filter(outcome %in% c("hba1c"), !models == "rcs") %>% 
  select(-trials)
# read in rcs with fixed knots
age_sex_model_coefs_9492_new2 <- read_csv("Data/agesexhba1crcsfixedknots_9492/coef.csv") %>% 
  select(-trials)
## read in empareg
age_sex_model_coefs_9492_new3 <- read_csv("Data/agesexhba1cempareg_9492/coef.csv") %>% 
  select(-trials)
age_sex_model_coefs_9492_new <- bind_rows(age_sex_model_coefs_9492_new1 ,
                                          age_sex_model_coefs_9492_new2,
                                          age_sex_model_coefs_9492_new3)
newref_nst <- newref %>% 
  select(nct_id, trial_lbl, arm_id_unq) %>% 
  nest(.by = nct_id) %>% 
  rename(lkp = data)

maketerms <- age_sex_model_coefs_9492_new %>% 
  filter(nct_id %in% newtrials, str_detect(term, "arm_f")) %>% 
  distinct(nct_id, term) %>% 
  nest(.by = nct_id) %>% 
  rename(coefs = data)
maketerms <- maketerms %>% 
  inner_join(newref_nst)
maketerms$renamedf <- map2(maketerms$coefs,
                           maketerms$lkp, function(cfs, lkp) {
                             lkp$old_term <- map(lkp$trial_lbl, ~ cfs %>% filter(str_detect(term, .x)) %>% pull(term))
                             lkp %>% 
                               unnest(old_term) %>% 
                               distinct()
                           })
maketerms <- maketerms %>% 
  select(nct_id, renamedf) %>% 
  unnest(renamedf) %>% 
  distinct()
maketerms$new_term <- pmap_chr(list(maketerms$old_term, maketerms$trial_lbl, maketerms$arm_id_unq), function(old, lbl, id) {
  str_replace(old, lbl, id)
})
maketerms <- maketerms %>% 
  select(nct_id, old_term, new_term) %>% 
  distinct()
## drop where a new term matches both a longer and a shorter new term
## eg arm_f=_2_Dula 0.75mg  matches arm_f=_2_uaa11218 and arm_f=_2_uaa11218 0.75mg
maketerms <- maketerms %>% 
  mutate(rplc_shrt = str_length(new_term)) %>% 
  arrange(rplc_shrt) %>% 
  distinct(nct_id, old_term, .keep_all = TRUE) %>% 
  select(-rplc_shrt)
## one trial still needs replacing because of a "plus"
maketerms <- maketerms %>% 
  mutate(new_term = str_replace(new_term, 
                                "Empagliflozin 10mg \\+ uaa11133 5mg",
                                "uaa11132"))

#write_csv(maketerms, "Data/cleaned_data/Data/assign_armids_newtrials.csv")
age_sex_model_coefs_9492_new <- age_sex_model_coefs_9492_new %>% 
  left_join(maketerms %>% rename(term = old_term)) %>% 
  mutate(term = if_else(!is.na(new_term), new_term, term)) %>% 
  select(-new_term)
age_sex_model_vcov_9492_new1 <- read_csv("Data/agesexhba1cmaceupdate_9492/vcov.csv") %>% 
  select(-nct_id2) %>% 
  semi_join(age_sex_model_coefs_9492_new1)
# rcs
age_sex_model_vcov_9492_new2 <- read_csv("Data/agesexhba1crcsfixedknots_9492/vcov.csv") %>% 
  select(-nct_id2) %>% 
  semi_join(age_sex_model_coefs_9492_new2)
# empareg
age_sex_model_vcov_9492_new3 <- read_csv("Data/agesexhba1cempareg_9492/vcov.csv") %>% 
  select(-nct_id2) %>% 
  semi_join(age_sex_model_coefs_9492_new3)
age_sex_model_vcov_9492_new <- bind_rows(age_sex_model_vcov_9492_new1,
                                          age_sex_model_vcov_9492_new2,
                                          age_sex_model_vcov_9492_new3)
age_sex_model_vcov_9492_new <- age_sex_model_vcov_9492_new %>% 
  left_join(maketerms %>% rename(row = old_term)) %>% 
  mutate(row = if_else(!is.na(new_term), new_term, row)) %>% 
  select(-new_term) %>% 
  left_join(maketerms %>% rename(col = old_term)) %>% 
  mutate(col = if_else(!is.na(new_term), new_term, col)) %>% 
  select(-new_term)
age_sex_model_diag_9492_new1 <- read_csv("Data/agesexhba1cmaceupdate_9492/diag_lm.csv") %>% 
  select(-nct_id2, -trials) %>% 
  semi_join(age_sex_model_coefs_9492_new1)
age_sex_model_diag_9492_new3 <- read_csv("Data/agesexhba1cempareg_9492/diag_lm.csv") %>% 
  select(-trials) %>% 
  semi_join(age_sex_model_coefs_9492_new3)
age_sex_model_diag_9492_new <- bind_rows(age_sex_model_diag_9492_new1,
                                         age_sex_model_diag_9492_new3)

rm(age_sex_model_coefs_9492_new1, age_sex_model_coefs_9492_new2, age_sex_model_coefs_9492_new3,
   age_sex_model_vcov_9492_new1, age_sex_model_vcov_9492_new2, age_sex_model_vcov_9492_new3,
   age_sex_model_diag_9492_new1,
   age_sex_model_diag_9492_new3)
res_9492 <- list(age_sex_model_coefs = age_sex_model_coefs_9492_new,
                 age_sex_model_vcov = age_sex_model_vcov_9492_new,
                 age_sex_model_diag = age_sex_model_diag_9492_new)
rm(age_sex_model_coefs_9492_new,
   age_sex_model_diag_9492_new,
   age_sex_model_vcov_9492_new)

## add new coefs, vcov and diags For 6115; in this one nct_id2 is informative as it where we split the trial NCT01778049_a and NCT01778049_b  ----
## no missing nct_id2 for HbA1c outcomes
age_sex_model_coefs_6115_new1 <- read_csv("Data/agesexhba1cmaceupdate_6115/coef.csv") %>% 
  filter(outcome == "hba1c", !models == "rcs") %>% 
  select(-trials) 
age_sex_model_coefs_6115_new2 <- read_csv("Data/agesexhba1crcsfixedknots_6115/coef.csv") %>% 
  select(-trials) 
age_sex_model_coefs_6115_new <- bind_rows(age_sex_model_coefs_6115_new1,
                                          age_sex_model_coefs_6115_new2)
## Join all vcov, first join new 6115 which has bocf and rcs
age_sex_model_vcov_6115_new1 <- read_csv("Data/agesexhba1cmaceupdate_6115/vcov.csv") %>% 
  semi_join(age_sex_model_coefs_6115_new1)
age_sex_model_vcov_6115_new2 <- read_csv("Data/agesexhba1crcsfixedknots_6115/vcov.csv") %>% 
  semi_join(age_sex_model_coefs_6115_new2)
age_sex_model_vcov_6115_new <- bind_rows(age_sex_model_vcov_6115_new1,
                                         age_sex_model_vcov_6115_new2)

age_sex_model_diag_6115_new1 <- read_csv("Data/agesexhba1cmaceupdate_6115/diag_lm.csv") %>% 
  select(-trials) %>% 
  semi_join(age_sex_model_coefs_6115_new1)
## no second one as there are no diagnostics for the rcs
age_sex_model_diag_6115_new <- age_sex_model_diag_6115_new1
rm(age_sex_model_coefs_6115_new1, age_sex_model_coefs_6115_new2,
   age_sex_model_vcov_6115_new1, age_sex_model_vcov_6115_new2, 
   age_sex_model_diag_6115_new1)

res_6115 <- list(age_sex_model_coefs = age_sex_model_coefs_6115_new,
                 age_sex_model_vcov = age_sex_model_vcov_6115_new,
                 age_sex_model_diag = age_sex_model_diag_6115_new)
rm(age_sex_model_coefs_6115_new,
   age_sex_model_vcov_6115_new,
   age_sex_model_diag_6115_new)

## Change from base only required for 9492 ----
newbase <- read.csv("Data/agesexhba1c_9492/hba1c_base_change_overall.csv")  %>% 
  filter(!nct_id %in% res$hba1c_base_change_overall$nct_id) 
newbase <- newbase %>% 
  mutate(trial_lbl = str_sub(arm_id_unq, 4)) %>% 
  select(-arm_id_unq) %>% 
  inner_join(newref %>% select(nct_id, trial_lbl, arm_id_unq))
res$hba1c_base_change_overall <- bind_rows(res$hba1c_base_change_overall,
                                           newbase %>% 
                                             filter(!nct_id %in% res$hba1c_base_change_overall$nct_id))
## add  new reference arms ----
# Trial reference arms for new. as expected 28 trials. 
# res$reference_arms_slct_reference <- 
# uaa10604 empa 10
# uaa10606 empa 25
# uaa10607 placebo
res$reference_arms_slct_reference <-  bind_rows(res$reference_arms_slct_reference,
            newref %>% 
              filter(arm_type == "ref") %>% 
              select(-arm_type) %>% 
              mutate(nct_id2 = nct_id))
empaarm <-  tibble(nct_id = rep("NCT01131676", 3)) %>% 
  mutate(trial_lbl = c("placebo",
                       "BI 10773 10mg",
                       "BI 10773 25mg"),
         arm_id_unq = c("placebo", "BI 10773 10mg", "BI 10773 25mg"),
         arm_label = c("placebo", "empagliflozin10", "empagliflozin25"),
         dc = c("placebo","sglt2","sglt2"),
         nct_id2 = nct_id)
res$reference_arms_slct_reference <- bind_rows(
  res$reference_arms_slct_reference,
  empaarm %>% filter(arm_label == "placebo"))

## Check new res against old ----
res_6115$age_sex_model_coefs 
res$age_sex_model_coefs

#6115. as expected no rcs
res_6115$age_sex_model_coefs %>% 
  anti_join(res$age_sex_model_coefs, by = c("nct_id", "nct_id2", "models", "term")) %>% 
  count(models)
# expecting f5 to f8 to differ; f1 to f4 match on terms. Same number of measures
res_6115$age_sex_model_coefs %>% 
  filter(models %in% c("b1", "f1", "f2", "f3", "f4")) %>% 
  anti_join(res$age_sex_model_coefs, by = c("nct_id", "nct_id2", "models", "term")) 
cmpr <- res_6115$age_sex_model_coefs %>% 
  filter(models %in% c("f4"), str_detect(term, "age10\\:")) %>% 
  left_join(res$age_sex_model_coefs, by = c("nct_id", "nct_id2", "models", "term")) %>% 
  mutate(estimate_diff = (estimate.x - estimate.y)) %>% 
  select(starts_with("estimate"), everything())
## 100% match no change in estimates or SEs. Use new for 6115 for consistency with LOCF and RCS
### No change for 9492
### No change for GSK
res_9492$age_sex_model_coefs %>% 
  ## only join old trials to check
  semi_join(res$age_sex_model_coefs %>% select(nct_id)) %>% 
  filter(models %in% c("b1", "f1", "f2", "f3", "f4")) %>% 
  anti_join(res$age_sex_model_coefs) 

res_gsk$age_sex_model_coefs %>% 
  semi_join(res$age_sex_model_coefs %>% select(nct_id)) %>% 
  filter(models %in% c("b1", "f1", "f2", "f3", "f4")) %>% 
  anti_join(res$age_sex_model_coefs) 

## Combine all coefs and diags data ----
res_gsk$smry <- NULL
res_new <- list(v6115 = map(res_6115, ~ .x),
                v9492 = map(res_9492, ~ .x %>% mutate(nct_id2 = nct_id)),
                gsk =   map(res_gsk,  ~ .x)
                )
res_new <- transpose(res_new)
res_new <- map(res_new, ~ bind_rows(.x, .id = "datasource"))
res_new$age_sex_model_vcov <- res_new$age_sex_model_vcov %>% 
  select(-trials)
map(res_new, ~ .x %>% group_by(datasource) %>% slice(1:3))
# two missing values only both for rcs models; otherwise no issues. Data is same structure
anymissing <- map(res_new, ~ .x %>% 
                    summarise(across(everything(), ~ any(is.na(.x)))))
res$age_sex_model_coefs <- res_new$age_sex_model_coefs %>% select(-datasource)
res$age_sex_model_vcov <- res_new$age_sex_model_vcov %>% select(-datasource)
res$age_sex_model_diags <- res_new$age_sex_model_diag %>% select(-datasource)

## fix age continuous so has nct_id2 ----
res$age_distribution_baseline_continuous <- res$age_distribution_baseline_continuous %>% 
  mutate(nct_id2 = case_when(
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00001", "ipd00002") ~ paste0(nct_id, "_a"),
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00003", "ipd00004") ~ paste0(nct_id, "_b"),
    TRUE ~ nct_id)
  )

saveRDS(res, "Scratch_data/ipd_consolidated.Rds")
write_csv(empaarm, "Scratch_data/empaarms.csv")
