library(tidyverse)
library(multinma)
library(truncnorm)
## Functions ----
source("Scripts/common_functions/Scripts/combine_sd.R")
CnvrtCorrMatrix <- function(a){
  ## recovery whole matrix by duplication
  allnames <- union(a$row, a$col)
  a <- bind_rows(a, 
                 a %>% 
                   rename(row = col, col = row),
                 tibble(row = allnames, col = allnames, r = 1)) %>% 
    distinct()
  # convert into matrix format
  a <- a %>% 
    spread(col, r)
  a_row <- a$row
  a$row <- NULL
  a <- as.matrix(a)
  if (any(is.na(a))) warning("Missing values in matrix")
  rownames(a) <- a_row
  a
}

exclusions <- read_csv("Data/exclusions_update.csv")
exclusions <- exclusions %>% 
  mutate(exclusion_reason = if_else(trial_id == "UMIN000018395",
                                  paste0(exclusion_reason, ". MACE no events"),
                                  exclusion_reason))
write_csv(exclusions, "Data/exclusions_update.csv")

## read in aggregate data and arm metadata ----
readRDS("Scratch_data/mace_arms_agg_data.Rds") %>% 
  list2env(envir = .GlobalEnv)
mace_tbl_age <- read_csv("Outputs/age_summary_trials_mace.csv")
## check HRs in agg model against calculated (approx as using rate ratio) ----
## all are very similar
mace_agg <- mace_agg %>% 
  mutate(rate = r/pt,
         hr = exp(loghr)) %>% 
  arrange(treat_cmpr) %>% 
  group_by(nct_id) %>% 
  mutate(rr = rate/rate[1],
         lrr = log(rr)) %>% 
  ungroup() %>% 
  arrange(nct_id, treat_cmpr)

## read in pseudo ipd. Use for estimating correlations between covariates ----
pseudo <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
pseudo <- pseudo %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, arm_lvl, trtcls5))

## Create table 1a wiht a full listing of MACE trials ----
mace_tbl_name <- bind_rows(agg = mace_agg_trial_level,
                      ipd = mace_ipd_trial_level,
                      .id = "data_lvl") 
mace_tbl_name <- mace_tbl_name %>% 
  mutate(fu_years = case_when(
    timepoint_unit == "months" ~ max_follow_up/12,
    timepoint_unit == "weeks" ~ max_follow_up/52,
    timepoint_unit == "years" ~ max_follow_up
  )) %>% 
  select(data_lvl, nct_id, trial_name, fu_years)
mace_tbl_dc <- mace_arms %>% 
  filter(trtcls5 %in% c("A10BH","A10BJ","A10BK")) %>% 
  distinct(dc, nct_id)
mace_tbl_plac <- mace_arms %>% 
  group_by(nct_id) %>% 
  summarise(placebo = if_else(any(drug_name == "placebo"), 1L, 0L)) %>% 
  ungroup()

mace_tbl_arm <- mace_arms %>% 
  filter(!atc == "placebo") %>% 
  mutate(drug_tot = paste(drug_name, drug_dose, drug_unit)) %>% 
  group_by(nct_id) %>% 
  summarise(drug_tot = paste(drug_tot, collapse = " vs ")) %>% 
  ungroup() %>% 
  mutate(drug_tot = str_replace(drug_tot, "\\|", "/"))
mace_tbl_age
mace_tbl_sex_agg <- mace_agg %>% 
  group_by(nct_id) %>% 
  summarise(male_prcnt = weighted.mean(male_prcnt, participants),
            participants = sum(participants)) %>% 
  ungroup()
mace_tbl_sex_ipd <- pseudo %>% 
  group_by(nct_id) %>% 
  summarise(male_prcnt = 100* mean(sex == "M")) %>% 
  ungroup()
mace_tbl_sex <- bind_rows(mace_tbl_sex_agg,
                          mace_tbl_sex_ipd)
rm(mace_tbl_sex_agg,
   mace_tbl_sex_ipd)
names(mace_tbl_age)[-(1:2)] <- paste0("age_", names(mace_tbl_age)[-(1:2)])
mace_tbl <- mace_tbl_name %>% 
  inner_join(mace_tbl_dc) %>% 
  inner_join(mace_tbl_arm) %>% 
  inner_join(mace_tbl_plac) %>%
  inner_join(mace_tbl_age) %>% 
  inner_join(mace_tbl_sex) %>% 
  ## note participants is same as n
  select(-participants) %>% 
  rename(activ_tx = drug_tot,
         participants = n) %>% 
  select(dc, nct_id, trial_name, data_lvl, activ_tx, placebo, participants, fu_years, everything()) %>% 
  arrange(dc, data_lvl, nct_id) 
mace_tbl <- mace_tbl %>% 
  mutate(data_lvl = if_else(nct_id %in% mace_agg_age$nct_id, "sg", data_lvl))
write_csv(mace_tbl, "Outputs/manuscript_table1b_machine_readable.csv", na = "")

mace_tbl %>% 
  filter(!nct_id == "UMIN000018395") %>% 
  summarise(n = length(nct_id),
            age_s = CombSdVectorised(participants, age_m, age_s),
            age_m = weighted.mean(age_m, participants),
            male_prcnt1 = weighted.mean(male_prcnt, participants),
            male_prcnt2 = 0.01* sum(100*male_prcnt*participants)/sum(participants)) %>% 
  mutate(female_prcnt = 100 - male_prcnt1)

mace_tbl_neat <- mace_tbl %>% 
  mutate(across(starts_with("age"), ~ round(.x, 1) %>% formatC(digits = 1, format = "f")),
        age = paste0(age_m , " (", age_s, ") [", age_q05, "-", age_q95, "]"),
        data_lvl = str_to_upper(data_lvl),
        dc = dc %>% 
          str_replace("dpp-4", "DPP-4") %>% 
          str_replace("sglt-2", "SGLT2") %>% 
          str_replace("glp-1", "GLP-1")) %>% 
  select(-starts_with("age_"))
write_csv(mace_tbl_neat, "Outputs/manuscript_table1b.csv", na = "")

## Make plots of "raw" coefficeints ----
mace_agg_plt_df <- mace_agg %>% 
  mutate(
         active = !is.na(loghr),
         noveldc = trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  group_by(nct_id) %>% 
  mutate(drug_lbl = paste0(arm_lvl[active], " vs ", arm_lvl[!active]),
         trtcls_trt = intersect(dc, c(
           "A10BH dipeptidyl peptidase 4 (dpp-4) inhibitors", 
           "A10BJ glucagon-like peptide-1 (glp-1) analogues",  
           "A10BK sodium-glucose co-transporter 2 (sglt2) inhibitors"))) %>% 
  ungroup()
mace_agg_plt <- ggplot(mace_agg_plt_df,
                    aes(y = loghr,
                        ymin = loghr -2*se,
                        ymax = loghr + 2*se,
                        x = paste0(nct_id, ":", drug_lbl),
                        shape = noveldc))
mace_agg_dc <- mace_agg %>% 
  distinct(nct_id, arm_lvl, dc)
mace_agg_sex_plt_df <- mace_agg_sex %>% 
  inner_join(mace_agg_dc) %>% 
  mutate(
    active = !is.na(loghr),
    noveldc = trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  group_by(nct_id, subgroup, level_cat) %>% 
  mutate(drug_lbl = paste0(arm_lvl[active], " vs ", arm_lvl[!active]),
         trtcls_trt = intersect(dc, c(
           "A10BH dipeptidyl peptidase 4 (dpp-4) inhibitors", 
           "A10BJ glucagon-like peptide-1 (glp-1) analogues",  
           "A10BK sodium-glucose co-transporter 2 (sglt2) inhibitors"))) %>% 
  ungroup()  
mace_agg_sex_plt <- ggplot(mace_agg_sex_plt_df,
       aes(y = loghr,
           ymin = loghr -2*se,
           ymax = loghr + 2*se,
           x = paste0(nct_id, ":", drug_lbl),
           shape = noveldc,
           colour = level_cat))  
mace_agg_age_plt_df <- mace_agg_age %>% 
  inner_join(mace_agg_dc) %>% 
  mutate(active = !is.na(loghr),
         noveldc = trtcls5 %in% c("A10BH", "A10BJ", "A10BK"),
         level_max = if_else(level_max == 150, "-", as.character(level_max))) %>% 
  group_by(nct_id, arm_lvl) %>% 
  mutate(age_rdr = rev(order(level_min)) %>% as.character()) %>% 
  ungroup() %>% 
  group_by(nct_id, subgroup, level_min, level_max) %>% 
  mutate(drug_lbl = paste0(arm_lvl[active], " vs ", arm_lvl[!active]),
         trtcls_trt = intersect(dc, c(
           "A10BH dipeptidyl peptidase 4 (dpp-4) inhibitors", 
           "A10BJ glucagon-like peptide-1 (glp-1) analogues",  
           "A10BK sodium-glucose co-transporter 2 (sglt2) inhibitors"))) %>% 
  ungroup()
mace_agg_age_plt <- ggplot(mace_agg_age_plt_df,
                           aes(y = loghr,
                               ymin = loghr -2*se,
                               ymax = loghr + 2*se,
                               x = paste0(nct_id, ":", drug_lbl),
                               shape = noveldc,
                               colour = age_rdr)) 
plotlst <- list(main = mace_agg_plt,
                age = mace_agg_age_plt,
                sex = mace_agg_sex_plt)
plotlst <- map(plotlst, ~ .x + 
                 geom_hline(yintercept = 0, linetype = "dashed") +
                 geom_point(position = position_dodge(0.5)) +
                 geom_linerange(position = position_dodge(0.5)) + 
                 coord_flip() +
                 scale_x_discrete("")  +
                 facet_wrap(~trtcls_trt, scales = "free_y", ncol = 1) )
saveRDS(plotlst, "Scratch_data/raw_hrs_mace_plot.Rds")
pdf("Outputs/raw_mace_hrs_plots.pdf", width = 30, height = 10)
plotlst[[1]] 
plotlst[[1]] <- plotlst[[1]] %+% (mace_agg_plt_df %>% filter(nct_id %in% 
                                                               c(mace_agg_age_plt_df$nct_id,
                                                                 mace_agg_sex_plt_df$nct_id)))
cowplot::plot_grid(plotlist = plotlst, nrow = 1)
dev.off()
rm(mace_agg_plt, mace_agg_sex_plt, mace_agg_age_plt, plotlst,
   mace_agg_plt_df, mace_agg_sex_plt_df, mace_agg_age_plt_df)

## read in coefficients and variance/covariance matrix ----
cfs <- read_csv("Data/vivli_mace/model_coefficients.csv")
vcv <- read_csv("Data/vivli_mace/model_vcv.csv")
## note 4 more coefficients trials than VCV. Reason for this is that in one model (for 4 trials) 
## there is only a single parameter. For the other two trials there are 3 arms 
## create vcov with a single value of 1 for this
## so still have a vcv
cfs <- cfs %>% 
  group_by(repo, nct_id, models) %>% 
  nest() %>% 
  ungroup()
vcv <- vcv %>% 
  group_by(repo, nct_id, models) %>% 
  nest() %>% 
  ungroup()
cfs_no_r <- cfs %>% 
  anti_join(vcv %>% select(-data))
vcv_add <- cfs_no_r
vcv_add$data <- map(vcv_add$data, ~ .x %>% 
  mutate(row = term, 
         col = term,
         r = 1) %>% 
    select(row, col, r))
vcv <- bind_rows(vcv_add, vcv)  
rm(vcv_add)
cfs <- cfs %>% 
  semi_join(vcv %>% select(-data))
cfs <- cfs %>% 
  inner_join(vcv %>% rename(vcv = data))
cfs$r <- map(cfs$vcv, CnvrtCorrMatrix)
cors <- cfs %>% 
  select(repo, nct_id, models, r)
cfs$r <- NULL

## Convert regression data in multinma format ----
MakeMaceData <- function(modeln) {
  
  cfs <- cfs %>% 
    filter(models == modeln)
  cfs <- cfs %>% 
    unnest(data)
  ## set up names as per set_agd_regression
  cfs <- cfs %>% 
    mutate(age10c = if_else(str_detect(term, "age10c"), 10L, NA_integer_),
           male = if_else(str_detect(term, "sexM"), TRUE, NA),
           trt = term %>% 
             str_remove("\\:") %>% 
             str_remove("age10c") %>% 
             str_remove("sexMALE") %>% 
             str_remove("sexM") %>% 
             str_remove("arm_f"))
  ## check treatment arms
  cfs %>% 
    filter(!trt == "") %>% 
    count(trt) %>% 
    left_join(mace_arms %>% 
                rename(trt = ipd_arm) %>% 
                select(nct_id, trt, arm_lvl, trtcls5))
  ## Note all reference arms for IPD trials are placebo so can simplify joining by assigning
  ## reference treatment to "placebo"
  cfs <- cfs %>% 
    left_join(mace_arms %>% 
                rename(trt = ipd_arm) %>% 
                select(nct_id, trt, arm_lvl, trtcls5)) %>% 
    group_by(nct_id) %>% 
    nest() %>% 
    ungroup()
  cfs$reference <- map(cfs$data, ~ .x %>% 
                         slice(1) %>% 
                         mutate(estimate = NA_real_, std.error = NA_real_,
                                term = NA_character_, age10c = 0L, male = FALSE,
                                trt = "placebo", arm_lvl = "placebo", trtcls5 = "place"))
  cfs$data <- map2(cfs$data, cfs$reference, ~ bind_rows(ref = .y, notref = .x, .id = "reference"))
  cfs$reference <- NULL
  cfs <- cfs %>% 
    unnest(data)
  cors <- cors %>% 
    filter(models == modeln)
  cors_lst <- cors$r
  names(cors_lst) <- cors$nct_id
  cors_lst[[1]]
  cfs %>% 
    filter(nct_id == names(cors_lst)[1])
  cfs <- cfs %>%
    mutate(ltime = 0)

  ## transform aggregate data and pseudo IPD into correct format for multinma ----
  pseudo <- pseudo %>% 
    mutate(male = if_else(sex == "M", 1L, 0L),
           time = time/365,
           ltime = log(time + 1))
  mace_agg <- mace_agg %>% 
    mutate(male_p = male_prcnt/100,
           time = mean_fu_days*participants,
           time = time/365,
           ltime = log(time))

  ## differs from mace_agg and mace_agg_sex because age subgroup cut-points
  # browser()
  mace_agg_age <- mace_agg_age %>% 
    select(-min_age, -max_age) %>% 
    rename(min_age = level_min, max_age = level_max) %>% 
    mutate(male_p = male_prcnt/100) %>% 
    inner_join(mace_agg %>% 
                 distinct(arm_id, ltime))

  mace_agg_age_sens <- mace_agg_age_sens %>%
    select(-min_age, -max_age) %>%
    rename(min_age = level_min, max_age = level_max) %>%
    mutate(male_p = male_prcnt/100) %>%
    left_join(mace_agg %>%
                 distinct(arm_id, ltime))
  
  mace_agg_sex <- mace_agg_sex %>% 
    mutate(male_p = male_prcnt/100) %>% 
    inner_join(mace_agg %>% 
                 distinct(arm_id, ltime))
  mace_agg_sex_sens <- mace_agg_sex_sens %>% 
    mutate(male_p = male_prcnt/100) %>% 
    left_join(mace_agg %>% 
                 distinct(arm_id, ltime))
  
  mace_agg <- mace_agg %>% 
    filter(!nct_id == "UMIN000018395")
  list(mace_agg = mace_agg,
               mace_agg_age = mace_agg_age,
               mace_agg_sex = mace_agg_sex,
       mace_agg_age_sens = mace_agg_age_sens,
       mace_agg_sex_sens = mace_agg_sex_sens,
               pseudo = pseudo,
               cors_lst = cors_lst,
               cfs = cfs)
}
f5 <- MakeMaceData("f5")
f1 <- MakeMaceData("f1")
f2 <- MakeMaceData("f2")

## make plots of "raw" coefficients for main effects and each interaction from IPD ----
cfsf5 <- f5$cfs %>% 
  filter(!is.na(term)) %>% 
  mutate(term_lbl = case_when(
    str_detect(term, "arm") & str_detect(term, "\\:") & male == TRUE ~ "arm_male_inter",
    str_detect(term, "arm") & str_detect(term, "\\:") & age10c == 10L ~ "arm_age_inter",
    str_detect(term, "arm") ~ "arm_main_age60_female",
    male == TRUE ~ "male",
    age10c == 10L ~ "age")) %>%
  group_by(nct_id) %>% 
  mutate(trtcls5 = intersect(trtcls5, c("A10BH", "A10BJ", "A10BK"))) %>% 
  ungroup() %>% 
  inner_join(mace_agg %>% distinct(trtcls5, dc)) %>% 
  mutate(arm_lvl = if_else(is.na(arm_lvl), "None", arm_lvl))
plotcfsf5 <- ggplot(cfsf5,
                  aes(x = paste0(trtcls5, ":",nct_id), y = estimate,
                      ymin = estimate - 2*std.error,
                      ymax = estimate + 2*std.error,
                      colour = arm_lvl,
                      shape = dc)) +
  geom_point(position = position_dodge(0.25)) +
  geom_linerange(position = position_dodge(0.25)) +
  facet_wrap(~term_lbl, scales = "free_y") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("Coefficients, full models")
plotcfsf5

cfsf1 <- f1$cfs %>% 
  filter(!is.na(term)) %>% 
  inner_join(mace_agg %>% distinct(trtcls5, dc)) 
plotcfsf1 <- ggplot(cfsf1,
       aes(x = paste0(trtcls5, ":",nct_id), y = estimate,
           ymin = estimate - 2*std.error,
           ymax = estimate + 2*std.error,
           shape = dc,
           colour =term)) +
  geom_point(position = position_dodge(0.25)) +
  geom_linerange(position = position_dodge(0.25)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggtitle("Coefficients, Treatment effect only models")
plotcfsf1 

saveRDS( list(simple = plotcfsf1, full = plotcfsf5), "Scratch_data/raw_coefs_mace_plot.Rds")
pdf("Outputs/raw_coefs_mace_plot.pdf", width = 20, height = 10)
plotcfsf1
plotcfsf5
dev.off()
saveRDS(f5, "Scratch_data/for_mace_regression_inter.Rds")
saveRDS(f1, "Scratch_data/for_mace_regression_nointer.Rds")
saveRDS(f2, "Scratch_data/for_mace_regression_nointercovs.Rds")

## get contrasts from model wiht no age, sex or interactions
f1 <- readRDS("Scratch_data/for_mace_regression_nointer.Rds")
aggedipd <- f1$cfs
aggedipd <- aggedipd %>% 
  select(nct_id, arm_lvl, loghr = estimate, se = std.error, trt, trtcls5, ltime)
censoring_distribution <- read_csv("Data/vivli_mace/censoring_distribution.csv")
male <- censoring_distribution %>% 
  select(nct_id:age_s) %>% 
  select(nct_id:participants) %>% 
  distinct() %>% 
  spread(sex, participants) %>% 
  mutate(male_p = M/(F+M),
         participants = F + M)  %>% 
  select(nct_id, trt = arm, male_p, participants) %>% 
  mutate(trt = str_to_lower(trt))
age <- censoring_distribution %>% 
  select(nct_id:age_s) %>% 
  group_by(nct_id, arm) %>% 
  summarise(age_s = CombSdVectorised(participants, age_m, age_s),
         age_m = weighted.mean(age_m, participants)) %>% 
  ungroup() %>% 
  select(nct_id, trt = arm, age_m, age_s) %>% 
  mutate(min_age = 0, max_age = 100,
         trt = str_to_lower(trt),
         age_mu = age_m,
         age_sigma = age_s)
aggedipd2 <- aggedipd %>% 
  inner_join(age) %>% 
  inner_join(male)
setdiff(names(f1$mace_agg), names(aggedipd2))
## need to drop 3 arm trials from no IPD analysis as cannot accomodate without SE in control arm
aggedipd2 <- aggedipd2 %>% 
  filter(! (arm_lvl %in% c("empagliflozin_10","canagliflozin_100") &
              nct_id %in% c("NCT01131676", "NCT01032629")))

f1$mace_agg <- bind_rows(f1$mace_agg,
                         aggedipd2)
saveRDS(f1, "Scratch_data/for_mace_regression_noipd.Rds")
