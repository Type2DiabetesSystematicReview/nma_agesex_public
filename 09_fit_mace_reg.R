library(tidyverse)
library(multinma)
library(truncnorm)
## Functions ----
CnvrtCorrMatrix <- function(a){
  ## recovery whole matrix by duplication
  allnames <- union(a$row, a$col)
  a <- bind_rows(a, 
                 a %>% 
                   rename(row = col, col = row),
                 tibble(row = allnames, col = allnames, r = 1)) %>% 
    distinct()
  # convert into matrix format
  a <- a %>% 
    spread(col, r)
  a_row <- a$row
  a$row <- NULL
  a <- as.matrix(a)
  if (any(is.na(a))) warning("Missing values in matrix")
  rownames(a) <- a_row
  a
}

## read in aggregate data and arm metadata ----
readRDS("Scratch_data/mace_arms_agg_data.Rds") %>% 
  list2env(envir = .GlobalEnv)

## check hrs in agg model against calculated (approx as using rate ratio) ----
mace_agg <- mace_agg %>% 
  mutate(rate = r/pt,
         hr = exp(loghr)) %>% 
  arrange(treat_cmpr) %>% 
  group_by(nct_id) %>% 
  mutate(rr = rate/rate[1],
         lrr = log(rr)) %>% 
  ungroup() %>% 
  arrange(nct_id, treat_cmpr)
## all are very similar

## read in pseudo ipd ----
pseudo <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
pseudo <- pseudo %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, arm_lvl, trtcls5))

## read in coefficients and variance/covariance matrix ----
cfs <- bind_rows(`6115` = read_csv("../from_vivli/Data/agesexmace_6115/age_sex_model_coefs.csv"),
                 `8697` = read_csv("../from_vivli/Data/agesexmace_8697/age_sex_model_coefs.csv"),
                 .id = "repo")
vcv <- bind_rows(`6115` = read_csv("../from_vivli/Data/agesexmace_6115/age_sex_model_vcov.csv"),
                 `8697` = read_csv("../from_vivli/Data/agesexmace_8697/age_sex_model_vcov.csv"),
                 .id = "repo")
## note 4 more coefficients tables than VCV as one model (4 trials) with only the treatment effect.
## The two trials where we have 3 arms still have a vcv
cfs <- cfs %>% 
  group_by(repo, nct_id, models) %>% 
  nest() %>% 
  ungroup()
vcv <- vcv %>% 
  group_by(repo, nct_id, models) %>% 
  nest() %>% 
  ungroup()
cfs_no_r <- cfs %>% 
  anti_join(vcv %>% select(-data))
cfs <- cfs %>% 
  semi_join(vcv %>% select(-data))
cfs$r <- map(vcv$data, CnvrtCorrMatrix)
cors <- cfs %>% 
  select(repo, nct_id, models, r)
cfs$r <- NULL

## Convert regression data ----
cfs <- cfs %>% 
  filter(models == "f5")
cfs <- cfs %>% 
  unnest(data)
## set up names as per set_agd_regression
cfs <- cfs %>% 
  mutate(age15 = if_else(str_detect(term, "age15"), 1, NA_real_),
         male = if_else(str_detect(term, "sexM"), TRUE, NA),
         trt = term %>% 
           str_remove("\\:") %>% 
           str_remove("age15") %>% 
           str_remove("sexMALE") %>% 
           str_remove("sexM") %>% 
           str_remove("arm_f"))
## check treatment arms
cfs %>% 
  filter(!trt == "") %>% 
  count(trt) %>% 
  left_join(mace_arms %>% 
              rename(trt = ipd_arm) %>% 
              select(nct_id, trt, arm_lvl, trtcls5))
## Note all reference arms are placebo so can simplify joining by assigning reference treatment to "placebo"
cfs <- cfs %>% 
  left_join(mace_arms %>% 
              rename(trt = ipd_arm) %>% 
              select(nct_id, trt, arm_lvl, trtcls5)) %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  ungroup()
cfs$reference <- map(cfs$data, ~ .x %>% 
                       slice(1) %>% 
                       mutate(estimate = NA_real_, std.error = NA_real_,
                              term = NA_character_, age15 = 0, male = FALSE,
                              trt = "placebo", arm_lvl = "placebo", trtcls5 = "place"))
cfs$data <- map2(cfs$data, cfs$reference, ~ bind_rows(ref = .y, notref = .x, .id = "reference"))
cfs$reference <- NULL
cfs <- cfs %>% 
  unnest(data)
cors <- cors %>% 
  filter(models == "f5")
cors_lst <- cors$r
names(cors_lst) <- cors$nct_id
cors_lst[[1]]
cfs %>% 
  filter(nct_id == names(cors_lst)[1])
cfs <- cfs %>%
  mutate(ltime = 0)
rm(cfs_no_r, cors, vcv)

## transform aggregate data and pseudo IPD ----
pseudo <- pseudo %>% 
  mutate(age15 = age/15,
         male = if_else(sex == "M", 1L, 0L),
         time = time/365,
         ltime = log(time + 1))
mace_agg <- mace_agg %>% 
  mutate(age_mu = age_mu/15,
         age_sigma = age_sigma/15,
         male_p = male_prcnt/100,
         time = mean_fu_days*participants,
         time = time/365,
         ltime = log(time))
mace_agg_age <- mace_agg_age %>% 
  mutate(age_mu = age_mu/15,
         age_sigma = age_sigma/15,
         male_p = male_prcnt/100) %>% 
  inner_join(mace_agg %>% 
              distinct(arm_id, ltime))
mace_agg_sex <- mace_agg_sex %>% 
  mutate(age_mu = age_mu/15,
         age_sigma = age_sigma/15,
         male_p = male_prcnt/100) %>% 
  inner_join(mace_agg %>% 
               distinct(arm_id, ltime))
mace_agg <- mace_agg %>% 
  filter(!nct_id == "UMIN000018395")

## Set-up aggregate and IPD data in different formats ----
agg_events <- set_agd_arm(data = mace_agg, 
                          study = nct_id, trt = arm_lvl, r = r, n = participants, 
                          trt_ref = "placebo", trt_class = trtcls5)
agg_hrs    <- set_agd_contrast(data = mace_agg,
                          study = nct_id, trt = arm_lvl, y = loghr, se = se, 
                          trt_ref = "placebo", trt_class = trtcls5, sample_size = participants)
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

ipd_pseudo <- set_ipd(data = pseudo, study = nct_id, trt = arm_lvl, r = event,
                      trt_ref = "placebo", trt_class = trtcls5)
ipd_regres_w_events <- set_agd_regression(cfs,
                                 study = nct_id,
                                 trt = arm_lvl,
                                 estimate = estimate,
                                 se = std.error,
                                 cor = cors_lst,
                                 trt_ref = "placebo",
                                 trt_class = trtcls5,
                                 regression = ~ (male + age15)*.trt + offset(ltime))
ipd_regres_w_coef <- set_agd_regression(cfs,
                                          study = nct_id,
                                          trt = arm_lvl,
                                          estimate = estimate,
                                          se = std.error,
                                          cor = cors_lst,
                                          trt_ref = "placebo",
                                          trt_class = trtcls5,
                                          regression = ~ (male + age15)*.trt)

## Create networks with 2x2 different combinations (aggregate = event/hr, IPD = pseudo/coefs) ----
net_lst <- list(
  e_p = combine_network(agg_events, ipd_pseudo),
  e_c = combine_network(agg_events, ipd_regres_w_events),
  h_p = combine_network(agg_hrs, ipd_pseudo),
  h_c = combine_network(agg_hrs, ipd_regres_w_coef),
  h_c_age = combine_network(agg_hrs_age, agg_hrs_noage, ipd_regres_w_coef),
  h_c_sex = combine_network(agg_hrs_sex, agg_hrs_nosex, ipd_regres_w_coef)
)
map(net_lst, plot, layout = "kk")
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

## Run models ----
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
## Commented out ones won't run
# h_p <- MyNMA(net_lst$h_p, 
#              mylink = "cloglog", 
#              myreg = ~ (male + age15)*.trt + offset(ltime))
h_c <- MyNMA(net_lst$h_c, 
             mylink = "identity", 
             myreg = ~ (male + age15)*.trt)
h_c_age <- MyNMA(net_lst$h_c_age, 
             mylink = "identity", 
             myreg = ~ (male + age15)*.trt)
h_c_sex <- MyNMA(net_lst$h_c_sex, 
                 mylink = "identity", 
                 myreg = ~ (male + age15)*.trt)
# h_c_re <- MyNMA(net_lst$h_c, 
#              mylink = "identity", 
#              myreg = ~ (male + age15)*.trt, fe_re = "random")
# e_c <- MyNMA(net_lst$e_c, 
#              mylink = "cloglog", 
#              myreg = ~ (male + age15)*.trt + offset(ltime))
# e_p <- MyNMA(net_lst$e_p, 
#              mylink = "cloglog", 
#              myreg = ~ (male + age15)*.trt + offset(ltime))

## Model runs fast. No divergent transitions. Rhat low. ESS high
saveRDS(h_c_sex, "Scratch_data/macesex_h_c_model.Rds")
saveRDS(h_c, "Scratch_data/mace_h_c_model.Rds")
saveRDS(h_c_age, "Scratch_data/maceage_h_c_model.Rds")

## For random effects model
saveRDS(h_c_re, "Scratch_data/macerand_h_c_model.Rds")
# The following variables have undefined values:  delta[18],The following
# variables have undefined values:  delta[19],The following variables have
# undefined values:  delta[20],The following variables have undefined values:
# delta[21],The following variables have undefined values:  delta[22],The
# following variables have undefined values:  delta[23],The following variables
# have undefined values:  delta[24],The following variables have undefined
# values:  delta[25].
