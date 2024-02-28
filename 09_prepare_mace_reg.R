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

## check HRs in agg model against calculated (approx as using rate ratio) ----
## all are very similar
mace_agg <- mace_agg %>% 
  mutate(rate = r/pt,
         hr = exp(loghr)) %>% 
  arrange(treat_cmpr) %>% 
  group_by(nct_id) %>% 
  mutate(rr = rate/rate[1],
         lrr = log(rr)) %>% 
  ungroup() %>% 
  arrange(nct_id, treat_cmpr)

## read in pseudo ipd. Use for estimating correlations between covariates ----
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
## note 4 more coefficients trials than VCV. Reason for this is that in one model (for 4 trials) 
## there is only a single parameter. For the other two trials there are 3 arms 
## create vcov with a single value of 1 for this
## so still have a vcv
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
vcv_add <- cfs_no_r
vcv_add$data <- map(vcv_add$data, ~ .x %>% 
  mutate(row = term, 
         col = term,
         r = 1) %>% 
    select(row, col, r))
vcv <- bind_rows(vcv_add, vcv)  
rm(vcv_add)
cfs <- cfs %>% 
  semi_join(vcv %>% select(-data))
cfs <- cfs %>% 
  inner_join(vcv %>% rename(vcv = data))
cfs$r <- map(cfs$vcv, CnvrtCorrMatrix)
cors <- cfs %>% 
  select(repo, nct_id, models, r)
cfs$r <- NULL

## Convert regression data in multinma format ----
MakeMaceData <- function(modeln) {
  
  cfs <- cfs %>% 
    filter(models == modeln)
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
  ## Note all reference arms for IPD trials are placebo so can simplify joining by assigning
  ## reference treatment to "placebo"
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
    filter(models == modeln)
  cors_lst <- cors$r
  names(cors_lst) <- cors$nct_id
  cors_lst[[1]]
  cfs %>% 
    filter(nct_id == names(cors_lst)[1])
  cfs <- cfs %>%
    mutate(ltime = 0)

  ## transform aggregate data and pseudo IPD into correct format for multinma ----
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
  
  list(mace_agg = mace_agg,
               mace_agg_age = mace_agg_age,
               mace_agg_sex = mace_agg_sex,
               pseudo = pseudo,
               cors_lst = cors_lst,
               cfs = cfs)
}
f5 <- MakeMaceData("f5")
f1 <- MakeMaceData("f1")

saveRDS(f5, "Scratch_data/for_mace_regression_inter.Rds")
saveRDS(f1, "Scratch_data/for_mace_regression_nointer.Rds")
