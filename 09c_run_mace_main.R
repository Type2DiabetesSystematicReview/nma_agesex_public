# 09c_run_mace
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(truncnorm)
library(purrr)

list2env(readRDS("Scratch_data/for_mace_regression_nointer.Rds"), envir = .GlobalEnv)
myregl <- ~ .trt + offset(ltime)
myreg <- ~ .trt 

## Set-up network ----
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

## Model running function
fe <-   nma(nwork,
      trt_effects = "fixed",
      link = "identity",
      regression = myreg,
      class_interactions = "common",
      prior_intercept = normal(scale = 10),
      prior_trt = normal(scale = 10),
      prior_reg = normal(scale = 10), chains = 4, cores = 4)
re <-   nma(nwork,
            trt_effects = "fixed",
            link = "identity",
            regression = myreg,
            class_interactions = "common",
            prior_intercept = normal(scale = 10),
            prior_trt = normal(scale = 10),
            prior_reg = normal(scale = 10), chains = 4, cores = 4)
saveRDS(list(fe = fe,
               re = re),
          paste0("Scratch_data/mace_", "nointer", ".Rds"))

