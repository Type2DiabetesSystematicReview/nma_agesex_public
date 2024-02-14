library(tidyverse)
library(multinma)

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

## read in pseudo ipd ----
pseudo <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
pseudo <- pseudo %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, drug_mdl, trtcl5))

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

## transform data ----
pseudo <- pseudo %>% 
  mutate(age15 = age/15,
         male = if_else(sex == "M", 1L, 0L),
         ltime = log(time + 1))
mace_agg <- mace_agg %>% 
  mutate(age_m = age_m/15,
         age_s = age_s/15,
         male_p = male_prcnt/100,
         time = mean_fu_days,
         ltime = log(time))

mace_agg <- mace_agg %>% 
  filter(!nct_id == "UMIN000018395")

## first attempt to fit pseudo-ipd and aggregate data ----
mynet <- combine_network(
  set_agd_arm(data = mace_agg, 
              study = nct_id, trt = drug_mdl, r = r, n = participants, 
              trt_ref = "placebo", trt_class = trtcl5),
  set_ipd(data = pseudo, study = nct_id, trt = drug_mdl, r = event,
          trt_ref = "placebo", trt_class = trtcl5))
mynet <- add_integration(mynet,
                         age15 = distr(qnorm, mean = age_m, sd = age_s),
                         male = distr(qbern, prob = male_p),
                         ltime = distr(qnorm, ltime, 0))
pseudocor <- mynet$int_cor
plot(mynet, layout = "kk")

# res <- nma(mynet,
#     trt_effects = "fixed",
#     link = "cloglog",
#     regression = ~ (age15 + male)*.trt + offset(ltime),
#     class_interactions = "common",
#     prior_intercept = normal(scale = 10),
#     prior_trt = normal(scale = 10),
#     prior_reg = normal(scale = 10))

## Next attempt to fit using regression data and aggregate data ----
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
              select(nct_id, trt, drug_mdl, trtcl5))
## Note all reference arms are placebo so can simplify joining by assigning reference treatment to "placebo"
cfs <- cfs %>% 
  left_join(mace_arms %>% 
              rename(trt = ipd_arm) %>% 
              select(nct_id, trt, drug_mdl, trtcl5)) %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  ungroup()
cfs$reference <- map(cfs$data, ~ .x %>% 
                        slice(1) %>% 
                        mutate(estimate = NA_real_, std.error = NA_real_,
                               term = NA_character_, age15 = 0, male = FALSE,
                               trt = "placebo", drug_mdl = "placebo", trtcl5 = "place"))
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

## take inverse cloglog of estimate
cloglog <- function(p) {
  log(-log(1-p))
}
icloglog <- function(y) {
  1 - exp(-exp(y))
}
cfs <- cfs %>%
  mutate(estimate_tform = icloglog(estimate),
         std.error_tform = icloglog(std.error))

reg_net <- combine_network(
  set_agd_arm(data = mace_agg, 
              study = nct_id, trt = drug_mdl, r = r, n = participants, 
              trt_ref = "placebo", trt_class = trtcl5),
  set_agd_regression(cfs,
                     study = nct_id,
                     trt = drug_mdl,
                     estimate = estimate_tform,
                     se = std.error_tform,
                     cor = cors_lst,
                     trt_ref = "placebo",
                     trt_class = trtcl5,
                     regression = ~ (male + age15)*.trt + offset(ltime)))
plot(reg_net)

reg_net <- add_integration(reg_net,
                age15 = distr(qnorm, mean = age_m, sd = age_s),
                male = distr(qbern, prob = male_p),
                ltime = distr(qnorm, ltime, 0), cor = pseudocor)

mace_FE <- nma(reg_net,
               trt_effects = "fixed",
               link = "cloglog",
               regression = ~ (male + age15)*.trt + offset(ltime),
               class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10), iter = 100, chains = 2)

