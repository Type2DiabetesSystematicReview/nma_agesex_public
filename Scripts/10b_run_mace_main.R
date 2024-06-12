# 09c_run_mace
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(truncnorm)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop("Need to pass arguments (eg in terminal via Rscript Scripts/10b_run_mace_main.R ARGS HERE) to indicate which model wish to run")

refe <- args[[1]]
filename <- paste0(refe, "_mace_", "nointer", ".Rds")
print(filename)

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
mdl <-   nma(nwork,
      trt_effects = refe,
      link = "identity",
      regression = myreg,
      class_interactions = "common",
      prior_intercept = normal(scale = 10),
      prior_trt = normal(scale = 10),
      prior_reg = normal(scale = 10), chains = 4, cores = 4)
saveRDS(mdl, filename)

