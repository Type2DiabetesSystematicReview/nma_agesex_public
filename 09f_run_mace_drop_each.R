# 09f_run_mace_dropping_each
library(dplyr)
library(tidyr)
library(multinma)
library(truncnorm)

alltrials <- c("NCT00968708", "NCT02465515", "NCT01131676", "NCT01032629", 
               "NCT01989754", "NCT02065791", "NCT00790205", "NCT01107886", 
               "NCT01144338", "NCT01147250", "NCT01179048", "NCT01243424",
               "NCT01394952", "NCT01455896", "NCT01703208", "NCT01720446", 
               "NCT01730534", "NCT01897532", "NCT01986881", "NCT02692716",
               "NCT03315143", "NCT03496298", "NCT03521934")

args <- commandArgs(trailingOnly=TRUE)
sg <- args[[1]]
whichtrial <- as.integer(args[[2]])
excld <- alltrials[whichtrial]
filename <- paste0("fixed_mace_", "agesex", "_", sg, "_", excld,".Rds")
print(filename)

list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
myreg <- ~ (male + age)*.trt 

cfs <- cfs %>% 
  rename(age = age10c) %>% 
  filter(!nct_id %in% excld)
cors_lst <- cors_lst[!names(cors_lst) == excld]
mace_agg <- mace_agg %>% 
  filter(!nct_id %in% excld)
mace_agg_age <- mace_agg_age %>% 
  filter(!nct_id %in% excld)
mace_agg_sex <- mace_agg_sex %>% 
  filter(!nct_id %in% excld)
pseudo <- pseudo %>% 
  filter(!nct_id %in% excld)


## Create networks with different combinations (aggregate = event/hr, main/sg, IPD = pseudo/coefs) ----
pseudo <- combine_network(set_agd_contrast(data = mace_agg,
                                           study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                           trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                          set_ipd(data = pseudo, study = nct_id, trt = arm_lvl, r = event,
                                  trt_ref = "placebo", trt_class = trtcls5))
# Note - approximately centre but do not scale age by subtracting 60
pseudo <- pseudo %>%  add_integration(
  age = distr(qtruncnorm, a = min_age-60, b = max_age-60, mean = age_mu-60, sd = age_sigma),
  male = distr(qbern, prob = male_p))
mycor <- pseudo$int_cor
rm(pseudo)

if(sg == "main") {
  nwork <- combine_network(set_agd_contrast(data = mace_agg,
                                            study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_regression(cfs,
                                              study = nct_id,
                                              trt = arm_lvl,
                                              estimate = estimate,
                                              se = std.error,
                                              cor = cors_lst,
                                              trt_ref = "placebo",
                                              trt_class = trtcls5,
                                              regression = myreg))
}
if(sg == "age") {
  nwork <- combine_network(set_agd_contrast(data = mace_agg_age,
                                            study = paste0(nct_id, "_", min_age, "_", max_age), trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_contrast(data = mace_agg %>% filter(!nct_id %in% mace_agg_age$nct_id),
                                            study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_regression(cfs,
                                              study = nct_id,
                                              trt = arm_lvl,
                                              estimate = estimate,
                                              se = std.error,
                                              cor = cors_lst,
                                              trt_ref = "placebo",
                                              trt_class = trtcls5,
                                              regression = myreg))
}
if(sg == "sex") {
  nwork <- combine_network(set_agd_contrast(data = mace_agg_sex,
                                            study = paste0(nct_id, "_", level_cat), trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_contrast(data = mace_agg %>% filter(!nct_id %in% mace_agg_sex$nct_id),
                                            study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_regression(cfs,
                                              study = nct_id,
                                              trt = arm_lvl,
                                              estimate = estimate,
                                              se = std.error,
                                              cor = cors_lst,
                                              trt_ref = "placebo",
                                              trt_class = trtcls5,
                                              regression = myreg))
}
## Add integration points. Note taking correlations from pseudo IPD networks
# Note - approximately centre but do not scale age by subtracting 60
nwork <-  add_integration(nwork,
                          age = distr(qtruncnorm,  a = min_age-60, b = max_age-60, mean = age_mu-60, sd = age_sigma),
                          male = distr(qbern, prob = male_p),
                          cor = mycor)
mdl <- nma(nwork,
           trt_effects = "fixed",
           link = "identity",
           regression = myreg,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10),
           chains = 4, cores = 4,
           control = list(max_treedepth = 15))
saveRDS(mdl, filename)
print(filename)
