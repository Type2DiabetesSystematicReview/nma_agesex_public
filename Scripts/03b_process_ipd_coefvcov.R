#03b_process__ipd_create_
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

cor2cov <- function(V, sigma) {
  p <- (d <- dim(V))[1L]
  if (!is.numeric(V) || length(d) != 2L || p != d[2L])
    stop("'V' is not a square numeric matrix")
  if (length(sigma) != p)
    stop("'sigma' is not a vector comformable as the standard deviations of 'V'")
  if (any(diag(V) != 1))
    warning("diag(.) contained non 1 entries.  Did you pass a correlation matrix?")
  sigma * V * rep(sigma, each = p)
}

res <- readRDS("Scratch_data/ipd_consolidated.Rds")
list2env(res, envir = .GlobalEnv)
empaarm <- read_csv("Scratch_data/empaarms.csv")

## simulate age and sex data ----
## Used ECDF. 6 trials still required sampling from normal distribution using mean and sd as these were too small to do so without compromising privacy
saveRDS(age_distribution_baseline_continuous, "Scratch_data/ipd_age_sex.Rds")

ColnamesPipe <- function(x, vct){
  colnames(x) <- vct
  x
}
# Add zero to two trials where not age zero. Set as same age as next value as same as earlier age NCT01381900 and NCT02065791 
age_distribution_baseline_continuous$splt <-  str_split(age_distribution_baseline_continuous$`Age (years) at quantiles (%)`, pattern = "\\,")
age_distribution_baseline_continuous$splt <- map(age_distribution_baseline_continuous$splt, ~ str_split_fixed(.x, "\\=", n = 2) %>% ColnamesPipe(c("x", "y")) %>% 
                                     as_tibble() %>% 
                                     mutate(across(everything(), ~ str_trim(.x) %>% as.integer()),
                                            x = x/100))
age_distribution_baseline_continuous$nozero <- map_lgl(age_distribution_baseline_continuous$splt, ~ ! any(.x$x ==0))
age_distribution_baseline_continuous$splt <- map2(age_distribution_baseline_continuous$nozero,
                                                  age_distribution_baseline_continuous$splt, function(condition, mydf){
                                                    if(is.na(condition) | !condition) mydf else {
                                                      a <- mydf %>% 
                                                        slice(1) %>% 
                                                        mutate(x = 0)
                                                      bind_rows(a, mydf)
                                                    }
                                                  })
age_distribution_baseline_continuous$ecdf_version <- map_lgl(age_distribution_baseline_continuous$splt, ~ !nrow(.x) ==1)
age_distribution_baseline_continuous$sim <- pmap(list(
  age_distribution_baseline_continuous$ecdf_version,
  age_distribution_baseline_continuous$splt, 
  age_distribution_baseline_continuous$participants,
  age_distribution_baseline_continuous$mean, 
  age_distribution_baseline_continuous$sd),
  function(ecdf_version, a, p, m, s) {
    if(ecdf_version) {
      fromdec <- approx(x = a$x, y =a$y, xout = seq(0, 1, 1/(p-1))) %>% 
        as_tibble() 
      result <- fromdec$y
    } else {
      result <- rnorm(p, m, s)
    }
    result
})
age_distribution_baseline_continuous$splt <- NULL
age_distr <- age_distribution_baseline_continuous %>% 
  select(nct_id, nct_id2, arm_id_unq, sex, participants, sim) %>% 
  unnest(sim) %>% 
  rename(age = sim) 

## add in empareg which we know is normally distributed
empa_base <- read_csv("Data/vivli_mace/censoring_distribution.csv") %>% 
  filter(nct_id == "NCT01131676") %>% 
  select(nct_id, arm, sex, participants, age_m, age_s) %>% 
  mutate(trial_lbl = if_else(arm == "Placebo", "placebo", arm))
empa_base$age <- pmap(list(empa_base$participants, empa_base$age_m, empa_base$age_s), function(n, m, s) {
  rnorm(n, m, s)
})
empa_base <- empa_base %>% 
  select(-age_m, -age_s) %>% 
  mutate(sex = if_else(sex == "F", "female", "male"),
         nct_id2 = nct_id) %>% 
  unnest(age)
empa_base <- empa_base %>% 
  inner_join(empaarm %>% 
               select(nct_id, trial_lbl, arm_id_unq))
empa_base <- empa_base %>% 
  select(nct_id, nct_id2, arm_id_unq, sex, participants, age)
age_distr <- bind_rows(age_distr,
                        empa_base)
# arm_f_2_uaa10253

age_distr <- age_distr %>% 
  # left_join(reference_arms_slct_reference %>% 
  #             select(nct_id, nct_id2, arm_id_unq) %>% 
  #             distinct() %>% 
  #             mutate(reference_arm = 1)) %>% 
  mutate(reference_arm = case_when(
    arm_id_unq %in% reference_arms_slct_reference$arm_id_unq ~ 1L,
    TRUE ~ 2L),
    arm_f = paste0("_", reference_arm, "_", arm_id_unq)) 

age_distr_smry <- age_distr %>% 
  group_by(nct_id, sex) %>% 
  summarise(m = mean(age),
            q5 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() 
age_distr_smry2 <- age_distr_smry %>% 
  arrange(sex) %>% 
  group_by(nct_id) %>% 
  summarise(across(c(m, q5, q95), ~ .x[1] - .x[2])) %>% 
  ungroup() %>% 
  gather("measure", "value", -nct_id) %>% 
  arrange(nct_id)

age_sex_plot <- ggplot(age_distr %>% 
                         filter(nct_id %in% sample(age_distr$nct_id %>% unique(), size = 20)),
                       aes(x = nct_id, 
                           y = age, fill = sex)) +
  geom_violin() +
  scale_fill_discrete("") +
  scale_x_discrete("Trial",guide = "none") +
  scale_y_continuous("Age (years)") +
  theme_minimal()
age_sex_plot2 <- ggplot(age_distr_smry2,
                    aes(x = measure, y = value)) +
  geom_violin() +
  scale_x_discrete("Mean, upper 5th and 95th centiles for age") +
  scale_y_continuous("Difference in age summary statistic women - men, in years")

saveRDS(age_sex_plot, "Scratch_data/age_comparison_by_sex_hba1c.Rds")
saveRDS(age_sex_plot2, "Scratch_data/age_comparison_by_sex_hba1c_smry.Rds")

## Obtain coefficients -----
age_sex_model_coefs <- age_sex_model_coefs %>% 
  group_by(nct_id, nct_id2, models, outcome, outcome_method) %>% 
  nest() %>% 
  ungroup()
## Recover variance-covariance matrix ----
age_sex_model_vcov <- age_sex_model_vcov %>% 
  group_by(nct_id, nct_id2, models, outcome, outcome_method) %>% 
  nest() %>% 
  ungroup()
# mysf <- safely(CnvrtCorrMatrix)
age_sex_model_vcov$r <- map(age_sex_model_vcov$data, CnvrtCorrMatrix)
# age_sex_model_vcov$error <- map_lgl(age_sex_model_vcov$r, ~ !is.null(.x$error))

age_sex_model_vcov$data <- NULL
age_sex_model_vcov <- age_sex_model_vcov %>% 
  inner_join(age_sex_model_coefs)

## rearrange terms so matches r in ordering
age_sex_model_vcov$data <- map2(age_sex_model_vcov$data, age_sex_model_vcov$r, ~ {
  term_order <- colnames(.y)
  .x %>% 
    arrange(match(term, term_order))
})

## convert to variance-covariance matrix
age_sex_model_vcov$vcv <- map2(age_sex_model_vcov$data, age_sex_model_vcov$r, ~ {
  r <- .y
  cfs <- .x
  if(! all(colnames(r) == cfs$term)) stop("Coefficient terms do not match correlation matrix")
  se <- cfs$std.error
  a <- cov.mat <- sweep(sweep(r, 1L, se, "*"), 2L, se, "*")
  a
  })

## Convert both R and vcv to positive definite (not currently due to rounding when exported)
age_sex_model_vcov$vcv_pd <- map(age_sex_model_vcov$vcv, ~ {
  as.matrix(Matrix::nearPD(.x)$mat)
})
age_sex_model_vcov$r_pd <- map(age_sex_model_vcov$r, ~ {
  as.matrix(Matrix::nearPD(.x, corr = TRUE)$mat)
})

## compare original and converted pd, very small differences
age_sex_model_vcov$vcv_diff <- map2_dbl(age_sex_model_vcov$vcv_pd, age_sex_model_vcov$vcv, ~ sum(abs(.x-.y)))
age_sex_model_vcov$r_diff <- map2_dbl(age_sex_model_vcov$r_pd, age_sex_model_vcov$r, ~ sum(abs(.x-.y)))
age_sex_model_vcov$vcv_diff %>% hist(breaks = 100)
age_sex_model_vcov$r_diff %>% hist(breaks = 100)

age_sex_model_vcov <- age_sex_model_vcov %>% 
  select(-c(r, vcv, r_diff, vcv_diff)) %>% 
  rename(vcv = vcv_pd,
         r = r_pd)

## want R as well as can use to run set_agd_regression
# age_sex_model_vcov$r <- NULL
## note that this is the same as the count at the end
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

##sample from distribution to get set of sampled coefficients based on vcov. 
## Sample from nobs *1.5 for convenience
age_sex_model_diags <- age_sex_model_diags %>% 
  mutate(rsd = (deviance/(nobs-1))^0.5)
## Not got diagnostics for all models (rcs dropped), so assume is same as for other nct_id/nct_id2 with f1
age_sex_model_diags_imp <- age_sex_model_diags %>% 
  filter(models == "f1") %>% 
  semi_join(age_sex_model_vcov %>% 
              filter(models == "rcs") %>% 
              select(nct_id, nct_id2)) %>% 
  mutate(models = "rcs", outcome_method = "locf", outcome = "hba1c") %>% 
  select(nct_id, nct_id2, models, outcome_method, outcome, nobs, rsd)
age_sex_model_diags_imp <- bind_rows(age_sex_model_diags,
                                     age_sex_model_diags_imp)

age_sex_model_vcov <- age_sex_model_vcov %>% 
  inner_join(age_sex_model_diags_imp %>% select(nct_id, nct_id2, models, outcome_method, outcome, nobs, rsd))
age_sex_model_vcov$sim <- pmap(list(age_sex_model_vcov$data, 
                                    age_sex_model_vcov$vcv, 
                                    age_sex_model_vcov$nobs,
                                    seq_along(age_sex_model_vcov$nct_id)), 
                               function(.x, .y, .z, i) {
                                 # if(i == 48) browser()
                                 smpls <- round(.z*1.5)
                                 vcv <- .y
                                 cf <- .x$estimate
                                 term <- .x$term
                                 if(!all(term == rownames(vcv))) stop("Matrix and coefficients dont match")
                                 res <- mvtnorm::rmvnorm(smpls, cf, vcv)
                                 colnames(res) <- term
                                 res
                               })
## Repeat mean coefficient across the range of the data. 
## Not efficient in terms of computing but re-uses code from simulation so easier to implement
age_sex_model_vcov$sngl_cf <- pmap(list(age_sex_model_vcov$data, 
                                        age_sex_model_vcov$vcv, 
                                        age_sex_model_vcov$nobs), 
                                   function(.x, .y, .z) {
                                     # browser()
                                     smpls <- round(.z*1.5)
                                     
                                     vcv <- .y
                                     vcv[] <- 0
                                     cf <- .x$estimate
                                     term <- .x$term
                                     if(!all(term == rownames(vcv))) stop("Matrix and coefficients dont match")
                                     res <- mvtnorm::rmvnorm(smpls, cf, vcv)
                                     colnames(res) <- term
                                     res
                                   })

age_sex_smpl_cfs <- age_sex_model_vcov %>% 
  select(nct_id, nct_id2, models, sim, sngl_cf, rsd, outcome_method, outcome) 
rm(age_sex_model_coefs, age_sex_model_vcov, age_sex_model_diags)

age_distr <- age_distr %>% 
  select(-participants) %>% 
  mutate(sex = if_else(sex == "male", 1L, 0L),
         age10 = age/10) %>% 
  group_by(nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup()


## adds in residual standard deviation
# In R's lm output "deviance" (as obtained from broom::glance) is the sum of squared residuals
# Therefore you can calculate the standard deviation of the residuals by (dev/(n-1))^0.5. This is identical to sd(residuals)
age_distr <- age_distr %>% 
  inner_join(age_sex_smpl_cfs %>% 
               filter(models == "b1"))
age_distr$mm <- map(age_distr$data, ~ model.matrix(~ sex*age10*arm_f, data = .x, ))

a <- age_distr$mm[[1]] %>% head(2)
b <- age_distr$sim[[1]] %>% head(2)
b <- b[ , c(intersect(colnames(a), colnames(b)), setdiff(colnames(b), colnames(a)))]
a <- a[ , c(intersect(colnames(a), colnames(b)), setdiff(colnames(a), colnames(b)))]
all(colnames(a) == colnames(b))

## Note that value_1 and value_1_sngl are the estimated baseline hba1c where the variation is
## based on the standard errors and residual standard deviation respectively. Both should have the same mean
## allowing for sampling variation
age_distr$data <- pmap(list(age_distr$data, 
                            age_distr$mm, 
                            age_distr$sim, 
                            age_distr$sngl_cf,
                            age_distr$rsd,
                            seq_along(age_distr$nct_id)), function(df, mm, cf, cf_sngl, rsd, i) {
                              # browser()
                              # if(i ==108) browser()
                              # i <<- i + 1L
                              if(!all(colnames(mm) %in% colnames(cf))) stop("Column names don't match")
                              if(!all(colnames(cf) %in% colnames(mm))) stop("Column names don't match")
                              cf <- cf[, colnames(mm)]
                              # if(nrow(cf) < nrow(mm)) stop("Not enough rows in data to match model matrix")
                              cf <- cf[sample(1:nrow(cf), nrow(mm), replace = TRUE), colnames(mm)]
                              if(!all(colnames(mm) %in% colnames(cf_sngl))) stop("Column names don't match")
                              if(!all(colnames(cf_sngl) %in% colnames(mm))) stop("Column names don't match")
                              cf_sngl <- cf_sngl[, colnames(mm)]
                              # if(nrow(cf_sngl) < nrow(mm)) stop("Not enough rows in data to match model matrix")
                              cf_sngl <- cf_sngl[sample(1:nrow(cf_sngl), nrow(mm), replace = TRUE), colnames(mm)]
                              
                              # calculate base
                              df$value_1 <- rowSums(mm * cf)
                              df$value_1_sngl <- rowSums(mm * cf_sngl)
                              df$value_1_sngl <- rnorm(length(df$value_1_sngl),
                                                       mean = df$value_1_sngl,
                                                       sd = rsd)
                              df
                            })
## Note needed to change prediction of value_2 to f4 model as didn't export f8 from Vivli as not needed in later modelling
age_distr <- age_distr %>% 
  select(-models, -sim, -sngl_cf, -rsd, -outcome_method) %>% 
  inner_join(age_sex_smpl_cfs %>% 
               filter(models == "f4"))
# age_distr$mm <- map(age_distr$data, ~ model.matrix(~ value_1 + age10 + arm_f + sex + age10:arm_f + arm_f:sex + value_1:arm_f, data = .x))
age_distr$mm <- map(age_distr$data, ~ model.matrix(~ value_1 + age10 + arm_f + sex + age10:arm_f + arm_f:sex, data = .x))

a <- age_distr$mm[[1]] %>% head(2)
b <- age_distr$sim[[1]] %>% head(2)
all(colnames(a) %in% colnames(b))

age_distr$data <- pmap(list(age_distr$data, 
                            age_distr$mm, age_distr$sim,
                            age_distr$sngl_cf,
                            age_distr$rsd), function(df, mm, cf, cf_sngl, rsd) {
                              
                              ## align each column
                              cf <- cf[ , colnames(mm)]
                              ## deal with lower number of observation in model than in data
                              cf <- cf[sample(1:nrow(cf), nrow(mm), replace = TRUE), colnames(mm)]
                              cf_sngl <- cf_sngl[sample(1:nrow(cf_sngl), nrow(mm), replace = TRUE), colnames(mm)]
                              
                              # calculate base
                              df$value_2 <- rowSums(mm * cf)
                              df$value_2_sngl <- rowSums(mm * cf_sngl)
                              df$value_2_sngl <- rnorm(length(df$value_2_sngl),
                                                       mean = df$value_2_sngl,
                                                       sd = rsd)
                              df
                            })

age_distr_lng <- age_distr %>% 
  select(nct_id, nct_id2, data) %>% 
  unnest(data)
age_distr_lng <- age_distr_lng %>% 
  rename(value_1_se = value_1,
         value_1_rsd = value_1_sngl,
         value_2_se = value_2,
         value_2_rsd = value_2_sngl)
saveRDS(age_distr_lng, "Scratch_data/simulated_ipd.Rds")
age_distr_lng <- readRDS("Scratch_data/simulated_ipd.Rds")
hba1c_base_smry <- age_distr_lng %>% 
  group_by(nct_id, arm_f) %>% 
  summarise(hba1cbase_m = mean(value_1_rsd),
            hba1cbase_s = sd(value_1_rsd)) %>% 
  ungroup()
write_csv(hba1c_base_smry, "Outputs/hba1c_base_summary.csv")
## check same with correlation/covariance data and append arms
cfs_cors <- readRDS("Scratch_data/combined_cfs_vcov.Rds")
all(age_distr_lng$nct_id %in% cfs_cors$nct_id) && all(cfs_cors$nct_id %in% age_distr_lng$nct_id)
age_distr_lng <- age_distr_lng %>% 
  distinct(nct_id, nct_id2, arm_id_unq, arm_f, reference_arm)
age_distr_lng <- age_distr_lng %>% 
  group_by(nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup()
cfs_cors <- cfs_cors %>% 
  inner_join(age_distr_lng %>% rename(arm_lkp = data))

## reformat for correct shape for set_agd_regression ----
## a bit slow but code easy to follow
reg_frmt <- pmap(
  list(cfs_cors$nct_id,
       cfs_cors$nct_id2,
       cfs_cors$models,
       cfs_cors$data,
       cfs_cors$arm_lkp,
       cfs_cors$vcv,
       cfs_cors$r),
  function(nct_id, nct_id2, models, cfs, lkp, vcv, crl) {
    # cfs <- cfs_cors$data[[1]]
    # lkp <- cfs_cors$arm_lkp[[1]]
    # vcv <- cfs_cors$vcv[[1]]
    # crl <- cfs_cors$r[[1]]
    ## add treatment variable
    
    ## arrange covariance and correlation matrix to match terms
    if(!all(colnames(crl) == cfs$term )) stop("coefficient terms and correlation matrix do not match")
    if(!all(colnames(vcv) == cfs$term )) stop("coefficient terms and correlation matrix do not match")
    
    # crl <- crl[cfs$term, cfs$term]
    # vcv <- vcv[cfs$term, cfs$term]
    
    ## add NA row for no treatments (note the correlation for this is set to NULL)
    cfs <- cfs %>% 
      mutate(trt = NA_character_)
    for (i in lkp$arm_id_unq) {
      cfs <- cfs %>% 
        mutate(trt = case_when(
          !is.na(trt) ~ trt,
          str_detect(term, i) ~ i,
          TRUE ~ trt))
    }
    addref <- setdiff(lkp$arm_id_unq, cfs$trt)
    cfs <- bind_rows(tibble(trt = addref),
                     cfs)
    cfs <- cfs %>% 
      left_join(lkp %>% rename(trt = arm_id_unq))
    
    ## add covariate columns
    ## Note than manual suggests listing as NA but example data has zeros and ones.
    cfs <- cfs %>% 
      mutate(sex = case_when(
        str_detect(term, "sex") ~ TRUE,
        reference_arm ==1 ~ FALSE,
        TRUE ~ NA))
    
    ## add covariate columns
    cfs <- cfs %>% 
      mutate(age10 = case_when(
        str_detect(term, "age10") ~ 10,
        reference_arm ==1 ~ 0,
        TRUE ~ NA_real_))
    
    tibble(nct_id = nct_id,
           nct_id2 = nct_id2,
           models = models,
           cfs = list(cfs),
           crl = list(crl),
           vcv = list(vcv))
  })
reg_frmt <- bind_rows(reg_frmt)
reg_frmt <- reg_frmt %>% 
  unnest(cfs)
saveRDS(reg_frmt, "Scratch_data/ipd_coefs_frmttd.Rds")

