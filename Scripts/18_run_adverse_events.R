# 18_run_ae_ipd
library(dplyr)
library(tidyr)
library(multinma)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop("Need to pass argument to indicate which model wish to run")
i <- as.integer(args[1])

mydata <- readRDS("Scratch_data/for_ae_regression.Rds")
# myreg <- ~ age10:.trt + sex:.trt

cfs <- mydata$cfs[[i]]
vcv <- mydata$crl[[i]]
mdl_chs <- paste0(mydata$outcome[[i]], "_", mydata$models[[i]])

## Set-up IPD data ----
nwork <- set_agd_regression(cfs,
                            study = nct_id2,
                            trt = arm_lvl,
                            estimate = estimate, 
                            cov = vcv,
                            trt_ref = "placebo",
                            trt_class = trtcls5,
                            regression = ~ age10:.trt + sex:.trt)

mdl <- nma(nwork,
               trt_effects = "fixed",
               link = "identity",
               regression = ~ age10:.trt + sex:.trt,
               class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10), 
               chains = 4, cores = 4)
saveRDS(mdl, paste0(mdl_chs, "_re.Rds"))
rm(mdl)
