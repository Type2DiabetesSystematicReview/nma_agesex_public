# 09b_run_mace
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)
library(truncnorm)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) ==0) model <- "f1" else {
  model <- args[[1]]
}

print(model)

if(model == "f5"){
  list2env(readRDS("Scratch_data/for_mace_regression_inter.Rds"), envir = .GlobalEnv)
  myregl <- ~ (male + age15)*.trt + offset(ltime)
  myreg <- ~ (male + age15)*.trt 
}

if(model == "f1"){
  list2env(readRDS("Scratch_data/for_mace_regression_nointer.Rds"), envir = .GlobalEnv)
  myregl <- ~ .trt + offset(ltime)
  myreg <- ~ .trt 
}


## Set-up aggregate and IPD data in different formatsfor combining into networks ----
## event and person-time data
agg_events <- set_agd_arm(data = mace_agg, 
                          study = nct_id, trt = arm_lvl, r = r, n = participants, 
                          trt_ref = "placebo", trt_class = trtcls5)
## HR data
agg_hrs    <- set_agd_contrast(data = mace_agg,
                               study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                               trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
## Age and sex subgroups (HRs only)
agg_hrs_age   <- set_agd_contrast(data = mace_agg_age,
                                  study = paste0(nct_id, "_", level_min, "_", level_max), trt = arm_lvl, y = loghr, se = se, 
                                  trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
agg_hrs_noage <- set_agd_contrast(data = mace_agg %>% filter(!nct_id %in% mace_agg_age$nct_id),
                                  study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                  trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
agg_hrs_sex   <- set_agd_contrast(data = mace_agg_sex,
                                  study = paste0(nct_id, "_", level_cat), trt = arm_lvl, y = loghr, se = se, 
                                  trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
agg_hrs_nosex <- set_agd_contrast(data = mace_agg %>% filter(!nct_id %in% mace_agg_sex$nct_id),
                                  study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                                  trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
## Pseudo IPD
ipd_pseudo <- set_ipd(data = pseudo, study = nct_id, trt = arm_lvl, r = event,
                      trt_ref = "placebo", trt_class = trtcls5)
## Regression coefficients for use with event aggregate data (added persontime)
ipd_regres_w_events <- set_agd_regression(cfs,
                                          study = nct_id,
                                          trt = arm_lvl,
                                          estimate = estimate,
                                          se = std.error,
                                          cor = cors_lst,
                                          trt_ref = "placebo",
                                          trt_class = trtcls5,
                                          regression = myregl)
## Regression coefficients for use with HR aggregate data (no persontime)
ipd_regres_w_coef <- set_agd_regression(cfs,
                                        study = nct_id,
                                        trt = arm_lvl,
                                        estimate = estimate,
                                        se = std.error,
                                        cor = cors_lst,
                                        trt_ref = "placebo",
                                        trt_class = trtcls5,
                                        regression = myreg)

## Create networks with different combinations (aggregate = event/hr, main/sg, IPD = pseudo/coefs) ----
net_lst <- list(
  e_p = combine_network(agg_events, ipd_pseudo),
  e_c = combine_network(agg_events, ipd_regres_w_events),
  h_p = combine_network(agg_hrs, ipd_pseudo),
  h_c = combine_network(agg_hrs, ipd_regres_w_coef),
  h_c_age = combine_network(agg_hrs_age, agg_hrs_noage, ipd_regres_w_coef),
  h_c_sex = combine_network(agg_hrs_sex, agg_hrs_nosex, ipd_regres_w_coef)
)
## plot networks ----
# map(net_lst, plot, layout = "kk")

## Add integration points. Note taking correlations from pseudo IPD networks
pseudos <- c("e_p", "h_p")
regs <- c("e_c", "h_c", "h_c_age", "h_c_sex")
net_lst[pseudos] <- map(net_lst[pseudos], ~ add_integration(.x,
                                                            age15 = distr(qtruncnorm, a = min_age, b = max_age, mean = age_mu, sd = age_sigma),
                                                            male = distr(qbern, prob = male_p),
                                                            ltime = distr(qnorm, ltime, 0)))
net_lst[regs] <- map2(net_lst[regs], net_lst[pseudos[c(1, 2, 2, 2)]], ~ add_integration(.x,
                                                                                        age15 = distr(qtruncnorm,  a = min_age, b = max_age, mean = age_mu, sd = age_sigma),
                                                                                        male = distr(qbern, prob = male_p),
                                                                                        ltime = distr(qnorm, ltime, 0), 
                                                                                        cor = .y$int_cor))
## Model running function
MyNMA <- function(mynet, mylink, myreg, fe_re = "fixed", ...) {
  nma(mynet,
      trt_effects = fe_re,
      link = mylink,
      regression = myreg,
      class_interactions = "common",
      prior_intercept = normal(scale = 10),
      prior_trt = normal(scale = 10),
      prior_reg = normal(scale = 10), chains = 4, cores = 4)
}

net_lst <- net_lst[c("h_c", "h_c_age", "h_c_sex")]

if(model == "f1") {
  fe <- MyNMA(net_lst[[1]], 
              mylink = "identity", 
              myreg = myreg, 
              fe_re = "fixed")
  re <- MyNMA(net_lst[[1]],
              mylink = "identity",
              myreg = myreg,
              fe_re = "random")
  saveRDS(list(fe = fe,
               re = re),
          "Scratch_data/mace_nointer.Rds")
} else {
  for(chs in 1:3){
    print(chs)
    nwork <- net_lst[[chs]]
    fe <- MyNMA(nwork, 
                mylink = "identity", 
                myreg = myreg, 
                fe_re = "fixed")
    re <- MyNMA(nwork,
                mylink = "identity",
                myreg = myreg,
                fe_re = "random")
    saveRDS(list(fe = fe,
                 re = re),
            paste0("Scratch_data/mace_", model, "_", chs, ".Rds"))
}}
