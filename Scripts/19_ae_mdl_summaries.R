# 19_plot_model_summaries_ae
library(tidyverse)
library(multinma)
library(rstan)

## read in data
tot <- readRDS("Scratch_data/for_ae_regression.Rds")

mdl_names <- list.files("FromVM/ae_agesex/", patt = "Rds$")
res1 <- map(mdl_names, ~ readRDS(paste0("FromVM/ae_agesex/", .x)))
names(res1) <- mdl_names %>% str_sub(1, -5)
names(res1) <- str_replace(names(res1), "^hba1c_", "hba1c_agesex_")

smpls <- map(res1, ~ as.data.frame(.x$stanfit))
betas <- map(smpls, ~ .x[str_detect(names(.x), "beta")])
smpls <- map(smpls, ~ .x[!str_detect(names(.x), "beta")])
ds <- map(smpls, ~ .x[!str_detect(names(.x), "^delta") & str_detect(names(.x), "^d")])
smpls <- map(smpls, ~ .x[!str_detect(names(.x), "^d")])

rm(smpls)
saveRDS(ds, "Scratch_data/tx_samples_ae.Rds")
saveRDS(betas, "Scratch_data/cov_nter_samples_ae.Rds")

beta <- map(res1, ~ summary(.x$stanfit)$summary %>% 
              as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "tosep")
names(beta) <- str_to_lower(names(beta))

divergent <- map(res1, ~ get_sampler_params(.x$stanfit, inc_warmup = FALSE))
divergent <- map(divergent, ~ map_int(.x, ~ sum(.x[, "divergent__"])))
divergent <- bind_rows(divergent, .id = "tosep") %>% 
  write_csv("Outputs/divergent_transition_count_ae.csv")

high_rhat <- beta %>% 
  filter(rhat > 1.05)
write_csv(high_rhat, "Outputs/nma_results_high_rhats_ae.csv")
write_csv(beta, "Outputs/betas_meta_analysis_ae.csv")

## produce table of priors
priors <- map(res1, ~ tibble(nm = names(.x$priors), value = .x$priors))
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
write_csv(priors, "Scratch_data/priors_meta_analysis_ae.csv")

## plots and tables 
