# 09e_run_mace_ipd
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(truncnorm)
library(purrr)

list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
myregl <- ~ (male + age15)*.trt + offset(ltime)
myreg <- ~ (male + age15)*.trt 


## Set-up aggregate and IPD data in different formatsfor combining into networks ----
## event and person-time data
## HR data
agg_hrs    <- set_agd_contrast(data = mace_agg,
                               study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                               trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
## Regression coefficients for use with event aggregate data (added persontime)
## Regression coefficients for use with HR aggregate data (no persontime)
nwork <- set_agd_regression(cfs,
                                        study = nct_id,
                                        trt = arm_lvl,
                                        estimate = estimate,
                                        se = std.error,
                                        cor = cors_lst,
                                        trt_ref = "placebo",
                                        trt_class = trtcls5,
                                        regression = myreg)
## Note that RE will not run for this
mdl_ipd <- nma(nwork,
           trt_effects = "fixed",
           link = "identity",
           regression = myreg,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), 
           chains = 2, cores = 4)
