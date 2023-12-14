#02_process_vivli
library(tidyverse)

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


cor2cov <- function(term, vr, r) {
  r_nm <- rownames(r)
  c_nm <- colnames(r)
  if (! all(r_nm == c_nm)) stop("rownames and column names of matrix dont match")
  if (! all(term %in% r_nm)) stop ("Have standard errors without matching terms in matrix")
  if (! all(r_nm %in% term)) stop ("Terms in matrix for which dont have standard errors/variances")
  # arrange matrix to match term
  if (! all(term == r_nm)) 
    { warning("Rearranging matrix so same order as standard error terms")
    r <- r[term, term]
  }
  vr %*% r %*% t(vr)
}

## read in IDs where have IPD ----
ipd1 <- read_csv("../from_vivli/Data/agesexhba1c_6115/hba1c_base_change_overall.csv")
ipd2 <- read.csv("../from_gsk/Data/agesex/hba1c_base_change_overall.csv")
ipd3 <- read.csv("../from_vivli/Data/agesexhba1c_8697/hba1c_base_change_overall.csv")
ipd_nct_id <- bind_rows(ipd1, ipd2, ipd3) %>% 
  distinct(nct_id) %>% 
  pull()

## read in vivli agesex results ----
allvivli <- list.files("../from_vivli/Data/agesexhba1c_6115/", patt = "csv$")
read_lines("../from_vivli/Data/agesexhba1c_6115/00_readme.txt") 
res1 <- map(allvivli, ~ read_csv(paste0("../from_vivli/Data/agesexhba1c_6115/", .x)))
names(res1) <- allvivli %>% str_sub(1, -5)
# bind categorical and continuous age data together. Only have continuous for gsk
res1$age_distribution_baseline_continuous <- bind_rows(res1$age_distribution_baseline_continuous,
                       res1$age_distribution_baseline_categorical) 
## read in gsk agesex results ----
allgsk <- list.files("../from_gsk//Data/agesex/", patt = "csv$")
read_lines("../from_gsk/README.md") 
res2 <- map(allgsk, ~ read_csv(paste0("../from_gsk/Data/agesex/", .x)))
names(res2) <- allgsk %>% str_sub(1, -5)
res1 <- res1[names(res2)]

## read in vivli agesex results from second vivli repository ----
allvivli <- list.files("../from_vivli/Data/agesexhba1c_8697/", patt = "csv$")
read_lines("../from_vivli/Data/agesexhba1c_8697/readme.txt") 
res3 <- map(allvivli, ~ read_csv(paste0("../from_vivli/Data/agesexhba1c_8697/", .x)))
names(res3) <- allvivli %>% str_sub(1, -5)
# drop csv only in second_vivli
res3$reference_trial_arm_to_arm_data_all_cleaned <- NULL
res <- pmap(list(res1, res2, res3), function(x, y, z) bind_rows(vivli1 = x, gsk = y, vivli2 = z))
list2env(res, envir = .GlobalEnv)
rm(res, res1, res2, allvivli, allgsk)

## simulate age and sex data ----
## crude simulation, can improve using ecdf
saveRDS(age_distribution_baseline_continuous, "Scratch_data/ipd_age_sex.Rds")
age_distr <- age_distribution_baseline_continuous %>% 
  select(nct_id:sd)
rm(age_distribution_baseline_continuous)
age_distr$sim <- pmap(list(age_distr$participants, age_distr$mean, age_distr$sd), function(p, m, s) rnorm(p, m, s))
age_distr <- age_distr %>% 
  select(nct_id, arm_id_unq, sex, sim) %>% 
  unnest(sim) %>% 
  rename(age = sim) 
# deal with NCT01778049 population a and population b
age_distr <- age_distr %>% 
  mutate(nct_id2 = case_when(
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00001", "ipd00002") ~ paste0(nct_id, "_a"),
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00003", "ipd00004") ~ paste0(nct_id, "_b"),
    TRUE ~ nct_id)
  )
reference_arms_slct_reference <- reference_arms_slct_reference %>% 
  mutate(nct_id2 = case_when(
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00001", "ipd00002") ~ paste0(nct_id, "_a"),
    nct_id == "NCT01778049" & arm_id_unq %in% c("ipd00003", "ipd00004") ~ paste0(nct_id, "_b"),
    TRUE ~ nct_id)
  )
age_distr <- age_distr %>% 
  left_join(reference_arms_slct_reference %>% 
              select(nct_id, nct_id2, arm_id_unq) %>% 
              distinct() %>% 
              mutate(reference_arm = 1)) %>% 
  mutate(reference_arm = if_else(is.na(reference_arm), 2L, reference_arm),
         arm_f = paste0("_", reference_arm, "_", arm_id_unq)) 

## Recover variance-covariance matrix ----
age_sex_model_coefs <- age_sex_model_coefs %>% 
  group_by(nct_id, nct_id2, models) %>% 
  nest() %>% 
  ungroup()
age_sex_model_vcov <- age_sex_model_vcov %>% 
  group_by(nct_id, nct_id2, models) %>% 
  nest() %>% 
  ungroup()
age_sex_model_vcov$r <- map(age_sex_model_vcov$data, CnvrtCorrMatrix)
age_sex_model_vcov$data <- NULL
age_sex_model_vcov <- age_sex_model_vcov %>% 
  inner_join(age_sex_model_coefs)

## convert to variance-covariance matrix
age_sex_model_vcov$vcv <- map2(age_sex_model_vcov$data, age_sex_model_vcov$r, ~ {
  term <- .x$term
  vr <- diag(.x$std.error^2)
  r <- .y
  res <- cor2cov(term = term, vr = vr, r = r)
  print(round(res, 3))
  
 res <- Matrix::nearPD(res, ensureSymmetry = TRUE, base.matrix = TRUE)$mat
 colnames(res) <- term
 rownames(res) <- term
res
})
age_sex_model_vcov$r <- NULL
saveRDS(age_sex_model_vcov, "Scratch_data/combined_cfs_vcov.Rds" )

## plot outputs of simpler age-sex treatment models ----
ast <- age_sex_model_vcov %>% 
  filter(models == "f2")
ast <- ast %>% 
  unnest(data) 
ast <- ast %>% 
  mutate(term = case_when(
    str_detect(term, "age\\:") ~ "age_interaction",
    str_detect(term, "\\:sex") ~ "sex_interaction",
    str_detect(term, "arm") ~ "arm",
    TRUE ~ term))
ast <- ast %>% 
  select(nct_id, nct_id2, term, estimate, std.error) %>% 
  gather("var", "value", estimate, std.error)
saveRDS(ast, "Scratch_data/ipd_raw_coefs.Rds")

##sample from distribution to get set of sampled coefficients. Sample same number as observations for convenience
age_sex_model_vcov <- age_sex_model_vcov %>% 
  left_join(age_sex_model_diags %>% select(nct_id, nct_id2, models, nobs)) 

age_sex_model_vcov$sim <- pmap(list(age_sex_model_vcov$data, age_sex_model_vcov$vcv, age_sex_model_vcov$nobs), function(.x, .y, .z) {
  # browser()
                                 vcv <- .y
                                 cf <- .x$estimate
                                 term <- .x$term
                                 if(!all(term == rownames(vcv))) stop("Matrix and coefficients dont match")
                                 res <- mvtnorm::rmvnorm(.z, cf, vcv)
                                 colnames(res) <- term
                                 res
})
age_sex_smpl_cfs <- age_sex_model_vcov %>% 
  select(nct_id, nct_id2, models, sim) 
rm(age_sex_model_coefs, age_sex_model_diags, age_sex_model_vcov)

age_distr <- age_distr %>% 
  mutate(sex = if_else(sex == "male", 1L, 0L),
         age10 = age/10) %>% 
  group_by(nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup()
age_distr <- age_distr %>% 
  inner_join(age_sex_smpl_cfs %>% 
               filter(models == "b1"))
age_distr$mm <- map(age_distr$data, ~ model.matrix(~ sex*age10*arm_f, data = .x))
a <- age_distr$mm[[1]] %>% head(2)
b <- age_distr$sim[[1]] %>% head(2)
all(colnames(a) == colnames(b))

age_distr$data <- pmap(list(age_distr$data, age_distr$mm, age_distr$sim), function(df, mm, cf) {
  ## deal with lower number of observation in model than in data
  chsrows <- sample(1:nrow(mm), size = nrow(cf))
  df <- df[chsrows,]
  mm <- mm[chsrows,]
  # calculate base
  df$value_1 <- rowSums(mm * cf)
  df
})

age_distr <- age_distr %>% 
  select(-models, -sim) %>% 
  inner_join(age_sex_smpl_cfs %>% 
               filter(models == "f8"))
age_distr$mm <- map(age_distr$data, ~ model.matrix(~ value_1 + age10 + arm_f + sex + age10:arm_f + arm_f:sex + value_1:arm_f, data = .x))
a <- age_distr$mm[[1]] %>% head(2)
b <- age_distr$sim[[1]] %>% head(2)
all(colnames(a) == colnames(b))
age_distr$data <- pmap(list(age_distr$data, age_distr$mm, age_distr$sim), function(df, mm, cf) {
  ## deal with lower number of observation in model than in data
  chsrows <- sample(1:nrow(mm), size = nrow(cf))
  df <- df[chsrows,]
  mm <- mm[chsrows,]
  # calculate base
  df$value_2 <- rowSums(mm * cf)
  df
})

age_distr_lng <- age_distr %>% 
  select(nct_id, nct_id2, data) %>% 
  unnest(data)
saveRDS(age_distr_lng, "Scratch_data/simulated_ipd.Rds")
