# 09b_run_mace
library(dplyr)
library(tidyr)
library(multinma)
library(truncnorm)

args <- commandArgs(trailingOnly=TRUE)
sg <- args[[1]]
refe <- args[[2]]
filename <- paste0(refe, "_mace_", "agesex", "_", sg, ".Rds")
print(filename)

list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
myreg <- ~ (male + age15)*.trt 


## Create networks with different combinations (aggregate = event/hr, main/sg, IPD = pseudo/coefs) ----
pseudo <- combine_network(set_agd_contrast(data = mace_agg,
                                          study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                          trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                         set_ipd(data = pseudo, study = nct_id, trt = arm_lvl, r = event,
                                 trt_ref = "placebo", trt_class = trtcls5))
pseudo <- pseudo %>%  add_integration(
                                 age15 = distr(qtruncnorm, a = min_age, b = max_age, mean = age_mu, sd = age_sigma),
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
if(sg == "sens") {
  nwork <- combine_network(set_agd_contrast(data = mace_agg,
                                            study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                            trt_ref = "placebo", trt_class = trtcls5, sample_size = participants), 
                           set_agd_regression(cfs %>% filter(!nct_id == "NCT01131676"),
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
                                           study = paste0(nct_id, "_", level_min, "_", level_max), trt = arm_lvl, y = loghr, se = se, 
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
nwork <-  add_integration(nwork,
                                age15 = distr(qtruncnorm,  a = min_age, b = max_age, mean = age_mu, sd = age_sigma),
                                male = distr(qbern, prob = male_p),
                               cor = mycor)
mdl <- nma(nwork,
           trt_effects = refe,
           link = "identity",
           regression = myreg,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), 
           chains = 4, cores = 4,
           control = list(max_treedepth = 15))
saveRDS(mdl, filename)

