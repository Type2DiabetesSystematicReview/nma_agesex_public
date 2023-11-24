# 05_ipdtrial_regression

#fit models multivariate
library(brms)
library(tidyverse)
count <- 1

## Select priors based on https://mc-stan.org/docs/stan-users-guide/multivariate-hierarchical-priors.html
oth <- readRDS("data_for_mars.Rds")
cfs_vcv <- readRDS("Scratch_data/combined_cfs_vcov.Rds" )
cfs_vcv <- cfs_vcv %>% 
  filter(models == "f4")

## reduce to pairwise with just the relevant interaction estimates
cfs_vcv$justarms <- map(cfs_vcv$data, ~ .x %>% 
                          filter(!str_detect(term, "age|sex|Intercept|value_1")) %>% 
                          pull(term) %>% 
                          unique())
justarms <- cfs_vcv %>% 
  select(nct_id, nct_id2, justarms) %>% 
  unnest(justarms) 
cfs_vcv <- cfs_vcv %>% 
  select(-justarms)
cfs_vcv <- justarms %>% 
  inner_join(cfs_vcv)
cfs_vcv$data <- map2(cfs_vcv$data, cfs_vcv$justarms, ~ .x %>% 
                      filter(str_detect(term, .y)))
cfs_vcv$vcv <- map2(cfs_vcv$vcv, cfs_vcv$justarms, ~ {
  .x[str_detect(rownames(.x), .y),
     str_detect(colnames(.x), .y)]
})

oth <- oth %>% 
  select(-agg, -dm_nets)
ipd1 <- oth$ipd[[1]]
mod1 <- lme4::lmer(result ~ base + age + sex + (drug_code|trtcls5), data = ipd1)


oth$ipd <- map(oth$ipd, ~ .x %>% 
                 distinct(nct_id, nct_id2, arm_f, drug_code, trtcls5, trtcls4, reference_arm))
oth <- oth %>% 
  unnest(ipd)
oth <- oth %>% 
  group_by(drug_regime_smpl, reference_arm, nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(lkp = data)
actv <- oth %>% 
  filter(reference_arm != 1)
refr <- oth %>% 
  filter(reference_arm == 1)
cfs_vcv <- cfs_vcv %>% 
  inner_join(actv)
cfs_vcv$data <- map2(cfs_vcv$data, cfs_vcv$lkp, ~ {
  a <- .x
  for (i in seq_along(.y$arm_f)) {
    a <- a %>% 
      mutate(term = str_replace_all(term, .y$arm_f[i], .y$drug_code[i]))
  }
  a
})
cfs_vcv$vcv <- map2(cfs_vcv$vcv, cfs_vcv$lkp, ~ {
  a <- colnames(.x)
  b <- rownames(.x)
  for (i in seq_along(.y$arm_f)) {
    a <- str_replace_all(a, .y$arm_f[i], .y$drug_code[i])
    b <- str_replace_all(b, .y$arm_f[i], .y$drug_code[i])
  }
  colnames(.x) <- a
  rownames(.x) <- b
  .x
})
## add in reference
cfs_vcv <- cfs_vcv %>% 
  select(-lkp, -reference_arm) %>% 
  inner_join(refr)
cfs_vcv$data <- map2(cfs_vcv$data, cfs_vcv$lkp, ~ {
  a <- .x
  for (i in seq_along(.y$arm_f)) {
    a <- a %>% 
      mutate(term = str_replace_all(term, "arm_f", paste0(.y$drug_code[i], "_")))
  }
  a
})
cfs_vcv$vcv <- map2(cfs_vcv$vcv, cfs_vcv$lkp, ~ {
  a <- colnames(.x)
  b <- rownames(.x)
  for (i in seq_along(.y$arm_f)) {
    a <- str_replace_all(a, "arm_f", paste0(.y$drug_code[i], "_"))
    b <- str_replace_all(b, "arm_f", paste0(.y$drug_code[i], "_"))
  }
  colnames(.x) <- a
  rownames(.x) <- b
  .x
})

## pull out what the pairwise comparison is
cfs_vcv$pairwise <- map_chr(cfs_vcv$data, ~ .x %>% 
                 filter(!str_detect(term, "age|sex")) %>% 
                 pull(term))
cfs_vcv <- cfs_vcv %>% 
  filter(!str_detect(pairwise, "uaa"))
pairwise <- tibble(orig = cfs_vcv$pairwise %>% 
  unique()) %>% 
  separate(col = orig, sep = "_|#", into = c(paste0("v", 1:3)), remove = FALSE) %>% 
  gather("order", "treat", -orig, na.rm = TRUE) %>% 
  mutate(trtcls5 = str_sub(treat, 1, 5),
         trtcls4 = str_sub(treat, 1, 4)) %>% 
  arrange(order) %>% 
  group_by(orig) %>% 
  summarise(across(c(trtcls5, trtcls4), ~ paste(.x, collapse = "_"))) %>% 
  ungroup() %>% 
  rename(pairwise = orig)
cfs_vcv <- cfs_vcv %>% 
  inner_join(pairwise)
## can add after running models
oth$mypriors <- vector(mode = "list", length = nrow(oth))

prior_orig <- TRUE

if (prior_orig) {
  oth$mypriors[oth$formtype == "simple"] <- map(oth$mypriors[oth$formtype == "simple"] , function(x) c(
    prior(constant(1), class = "sigma"),
    set_prior("normal(0, 1)", class = "b")))
  oth$mypriors[!oth$formtype == "simple"] <- map(oth$mypriors[!oth$formtype == "simple"] , function(x) c(
    prior(constant(1), class = "sigma"),
    set_prior("normal(0, 1)", class = "b"),
    set_prior("normal(0, 1)", class = "sd"),
    set_prior("lkj_corr_cholesky(1)", class = "L")))
} else {
  oth$mypriors[oth$formtype == "simple"] <- map(oth$mypriors[oth$formtype == "simple"] , function(x) c(
    prior(constant(1), class = "sigma"),
    set_prior("normal(0, 2)", class = "b")))
  oth$mypriors[!oth$formtype == "simple"] <- map(oth$mypriors[!oth$formtype == "simple"] , function(x) c(
    prior(constant(1), class = "sigma"),
    set_prior("normal(0, 2)", class = "b"),
    set_prior("normal(0, 2)", class = "sd"),
    set_prior("lkj_corr_cholesky(2)", class = "L")))
}


# oth_orig <- oth
# oth <- oth %>% 
#   filter(condition_cb == "Hypertension")
# oth <- oth %>% 
#   filter(model_type == "cnt_sq")

brms_make <- function(...){
  list(code = make_stancode(...),
       data = make_standata(...))
}


oth$res <- pmap(list(
  oth$myform,
  oth$cfs,
  oth$bd,
  oth$mypriors), function(myform, cfs, bd, mypriors) {
    print(count)
    count <<- count + 1
    brms_make(myform, 
              data = cfs, 
              data2 = list(covs = bd),
              prior = mypriors,
              control = list(adapt_delta = 0.999),
              cores = 4)  
  })

## save single trials for post-processing (no meta-analysis)
if (!prior_orig) {
  saveRDS(oth, "Scratch_data/mdls_mvn_vmbox.Rds")
} else {
  saveRDS(oth, "Scratch_data/mdls_mvn_vmbox_prior_orig.Rds")
}