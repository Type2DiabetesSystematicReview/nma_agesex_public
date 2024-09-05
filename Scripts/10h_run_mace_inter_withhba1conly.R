# 09b_run_mace
library(dplyr)
library(tidyr)
library(multinma)
library(truncnorm)

list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
myreg <- ~ (male + age)*.trt 

cfs <- cfs %>% 
  rename(age = age10c)

hba1c <-   c("NCT00790205", "NCT01243424", "NCT01455896", "NCT01703208", 
            "NCT01720446", "NCT01986881", "NCT02692716", "NCT03315143", "NCT00968708", 
            "NCT02465515", "NCT01032629", "NCT01989754", "NCT02065791", "NCT01131676")
## drop one disconnected trial
hba1c <- setdiff(hba1c, "NCT01243424")
cfs <- cfs %>% 
  filter(nct_id %in% hba1c)
mace_agg_sex <- mace_agg_sex %>% 
  filter(nct_id %in% hba1c)
mace_agg <- mace_agg %>% 
  filter(nct_id %in% hba1c)
psuedo <- pseudo %>% 
  filter(nct_id %in% hba1c)

## Create networks  ----
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

## Add integration points. Note taking correlations from pseudo IPD networks
# Note - approximately centre but do not scale age by subtracting 60
nwork <-  add_integration(nwork,
                          age = distr(qtruncnorm,  a = min_age-60, b = max_age-60, mean = age_mu-60, sd = age_sigma),
                          male = distr(qbern, prob = male_p),
                          cor = mycor)
mdl_fe <- nma(nwork,
           trt_effects = "fixed",
           link = "identity",
           regression = myreg,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), 
           chains = 4, cores = 4, 
           control = list(max_treedepth = 15))
saveRDS(mdl_fe, "mace_sensitivity_fe_hba1conly.Rds")
rm(mdl_fe)
mdl_re <- nma(nwork,
              trt_effects = "random",
              link = "identity",
              regression = myreg,
              class_interactions = "common",
              prior_intercept = normal(scale = 10),
              prior_trt = normal(scale = 10),
              prior_reg = normal(scale = 10), 
              chains = 4, cores = 4, 
              control = list(max_treedepth = 15))
saveRDS(mdl_re, "mace_sensitivity_re_hba1conly.Rds")
rm(mdl_re)