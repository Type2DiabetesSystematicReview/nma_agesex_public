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

## transform data ----
pseudo <- pseudo %>% 
  mutate(age = age/15,
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
                         age = distr(qnorm, mean = age_m, sd = age_s),
                         male = distr(qbern, prob = male_p),
                         ltime = distr(qnorm, ltime, 0))

plot(mynet, layout = "kk")

res <- nma(mynet,
    trt_effects = "fixed",
    link = "cloglog",
    regression = ~ (age + male)*.trt + offset(ltime),
    class_interactions = "common",
    prior_intercept = normal(scale = 10),
    prior_trt = normal(scale = 10),
    prior_reg = normal(scale = 10))

## Next attempt to fit using regression data and aggregate data ----
cfs <- cfs %>% 
  filter(models == "f5")
## map names to names in 