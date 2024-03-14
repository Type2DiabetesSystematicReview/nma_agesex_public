# 09e_runclass
library(dplyr)
library(tidyr)
library(multinma)
library(truncnorm)

args <- commandArgs(trailingOnly=TRUE)
refe <- args[[1]]
filename <- paste0(refe, "_mace_", "classoverall", ".Rds")
print(filename)

list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
myreg <- ~ (male + age)*.trt 

cfs <- cfs %>% 
  rename(age = age10c)

mace_agg <- mace_agg %>% 
  mutate(arm_lvl = trtcls5)
mace_agg_sex <- mace_agg_sex %>% 
  mutate(arm_lvl = trtcls5)
pseudo <- pseudo %>% 
  mutate(arm_lvl = trtcls5)
cfs <- cfs %>% 
  mutate(arm_lvl = trtcls5)

## Create networks with different combinations (aggregate = event/hr, main/sg, IPD = pseudo/coefs) ----
pseudo <- combine_network(set_agd_contrast(data = mace_agg,
                                           study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                           trt_ref = "place", trt_class = trtcls5, sample_size = participants), 
                          set_ipd(data = pseudo,
                                  study = nct_id, trt = arm_lvl, r = event,
                                  trt_ref = "place", trt_class = trtcls5))
# Note - approximately centre but do not scale age by subtracting 60
pseudo <- pseudo %>%  add_integration(
  age = distr(qtruncnorm, a = min_age-60, b = max_age-60, mean = age_mu-60, sd = age_sigma),
  male = distr(qbern, prob = male_p))
mycor <- pseudo$int_cor
rm(pseudo)

nwork <- combine_network(set_agd_contrast(data = mace_agg_sex,
                                          study = paste0(nct_id, "_", level_cat), 
                                          trt = arm_lvl, y = loghr, se = se, 
                                          trt_ref = "place", 
                                          trt_class = trtcls5, sample_size = participants), 
                         set_agd_contrast(data = mace_agg %>% filter(!nct_id %in% mace_agg_sex$nct_id),
                                          study = nct_id, 
                                          trt = arm_lvl, y = loghr, se = se, 
                                          trt_ref = "place", 
                                          trt_class = trtcls5, sample_size = participants), 
                           set_agd_regression(cfs,
                                              study = nct_id,
                                              trt = arm_lvl,
                                              estimate = estimate,
                                              se = std.error,
                                              cor = cors_lst,
                                              trt_ref = "place",
                                              trt_class = trtcls5,
                                              regression = myreg))


## Add integration points. Note taking correlations from pseudo IPD networks
# Note - approximately centre but do not scale age by subtracting 60
nwork <-  add_integration(nwork,
                          age = distr(qtruncnorm,  
                                      a = min_age-60, 
                                      b = max_age-60,
                                      mean = age_mu-60, 
                                      sd = age_sigma),
                          male = distr(qbern, prob = male_p),
                          n_int = 128,
                          cor = mycor)
mdl <- nma(nwork,
           trt_effects = refe,
           link = "identity",
           regression = myreg,
           class_interactions = "independent",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), 
           chains = 4, cores = 4,
           control = list(max_treedepth = 15))
saveRDS(mdl, filename)

