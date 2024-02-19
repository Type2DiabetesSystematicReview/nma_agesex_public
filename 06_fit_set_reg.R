# 06_fit_set_reg
library(tidyverse)
library(multinma)

tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")

dual_reg <- tot$reg[tot$drug_regime_smpl == "dual"][[1]]
dual_reg <- dual_reg %>% 
  filter(models == "f4")
## to check code runs drop the three trials in multiple settings. Will need to update
dual_agg <- tot$agg[tot$drug_regime_smpl == "dual"][[1]]
dual_agg <- dual_agg %>% 
  mutate(nct_id2 = nct_id)


crl <- dual_reg %>% 
  filter(is.na(term))
crl_final <- crl$crl
names(crl_final) <- crl$nct_id2
dual_net <- combine_network(
  set_agd_regression(dual_reg,
                     study = nct_id2,
                     trt = arm_lvl,
                     estimate = estimate,
                     se = std.error,
                     cor = crl_final,
                     # cov = pso_reg_cov,
                     trt_class = trtcls5,
                     regression = ~ (sex + age10)*.trt))


dual_FE <- nma(dual_net,
                  trt_effects = "fixed",
                  regression = ~ (sex + age10)*.trt,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10))
