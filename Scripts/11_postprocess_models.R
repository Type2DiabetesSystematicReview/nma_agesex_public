# 06b_read_models_vm

library(tidyverse)
library(multinma)
library(rstan)

## read in data
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
tot <- tot %>% 
  slice(3, 1:2)

mdl_names <- list.files("FromVM/hba1c_agesex/", patt = "Rds$")
res1 <- map(mdl_names, ~ readRDS(paste0("FromVM/hba1c_agesex/", .x)))
names(res1) <- mdl_names %>% str_sub(1, -5)

mdl_names <- list.files("FromVM/mace_agesex/", patt = "Rds$")
res3 <-  map(mdl_names, ~ readRDS(paste0("FromVM/mace_agesex/", .x)))
names(res3) <- mdl_names %>% str_sub(1, -5)

mdl_names <- list.files("FromVM/mace_nointer/", patt = "Rds$")
res4 <- map(mdl_names, ~ readRDS(paste0("FromVM/mace_nointer/", .x)))
names(res4) <- mdl_names %>% str_sub(1, -5)

mdl_names <- list.files("FromVM/mace_sensitivity/", patt = "Rds$")
res5 <- map(mdl_names, ~ readRDS(paste0("FromVM/mace_sensitivity/", .x)))
names(res5) <- mdl_names %>% str_sub(1, -5)

mdl_names <- list.files("FromVM/hba1c_noipd/", patt = "Rds$")
res6 <- map(mdl_names, ~ readRDS(paste0("FromVM/hba1c_noipd/", .x)))
names(res6) <- mdl_names %>% str_sub(1, -5)
names(res6) <- str_replace(names(res6), "^hba1c_", "hba1c_agesex_")

mdl_names <- list.files("FromVM/mace_noipd/", patt = "Rds$")
res7 <- map(mdl_names, ~ readRDS(paste0("FromVM/mace_noipd/", .x)))
names(res7) <- mdl_names %>% str_sub(1, -5)
names(res7) <- str_replace(names(res7), "main", "noipd")

res <- c(res1, res3, res4, res5, res6, res7)
rm(res1, res3, res4, res5, res6, res7)

smpls <- map(res, ~ as.data.frame(.x$stanfit))
betas <- map(smpls, ~ .x[str_detect(names(.x), "beta")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "beta")])
ds  <- map(smpls, ~ .x[!str_detect(names(.x), "^delta") & str_detect(names(.x), "^d")])
rm(smpls)
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^d")])
## the following drops all columns
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^mu")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^res")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^fitted")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^log_lik")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^lp__")])
# smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^tau")])
# map(smpls, colnames)
# rm(smpls)
saveRDS(ds, "Scratch_data/tx_samples.Rds")
saveRDS(betas, "Scratch_data/cov_nter_samples.Rds")


beta <- map(res, ~ summary(.x$stanfit)$summary %>% 
  as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "tosep")
names(beta) <- str_to_lower(names(beta))

divergent <- map(res, ~ get_sampler_params(.x$stanfit, inc_warmup = FALSE))
divergent <- map(divergent, ~ map_int(.x, ~ sum(.x[, "divergent__"])))
divergent <- bind_rows(divergent, .id = "tosep") %>% 
   write_csv("Outputs/divergent_transition_count.csv")


## Check with DP RE missing deltas
missings <- beta %>% 
  filter(is.na(rhat)) %>% 
  count(tosep, params)
write_csv(missings, "Outputs/nma_results_missing_delta_values.csv")
beta <- beta %>% 
  filter(!is.na(rhat))
high_rhat <- beta %>% 
  filter(rhat > 1.05)

write_csv(high_rhat, "Outputs/nma_results_high_rhats.csv")
write_csv(beta, "Outputs/betas_meta_analysis.csv")
tau <- beta %>% 
  filter(params == "tau")
write_csv(tau, "Outputs/taus_meta_analysis.csv")

## Aggregated IPD models did not converge, so also did the BOCF model for monotherapy 
failed_model <- high_rhat %>% 
  distinct(tosep)

## produce table of priors
priors <- map(res, ~ tibble(nm = names(.x$priors), value = .x$priors))
priors <- bind_rows(priors, .id = "tosep")
priors$value <- map_chr(priors$value, capture.output)
priors <- priors %>% 
  filter(!value == "NULL")
priors <- priors %>% 
  mutate(value = str_replace(value, '\\[1\\] \\\"sd\\\"', "standard deviation."))
priors <- priors %>% 
  mutate(value = if_else(nm %in% c("prior_het"), str_sub(value, 1, -2), value),
         nm = if_else(nm %in% c("prior_het", "prior_het_type"), "prior_het", nm)) %>% 
  group_by(tosep, nm) %>% 
  summarise(value = paste(value, collapse = ", on the ")) %>% 
  ungroup() %>% 
  arrange(tosep, nm)
rnm <- c("prior_het" = "D: Heterogeneity", 
         "prior_intercept" = "A: Intercept", 
         "prior_reg" = "C: Covariates main effects and treatment interactions", 
         "prior_trt" = "B: Main treatment effects")
priors <- priors %>% 
  mutate(nm = rnm[nm])
priors <- priors %>% 
  spread(nm, value, fill = "")
write_csv(priors, "Scratch_data/priors_meta_analysis.csv")

