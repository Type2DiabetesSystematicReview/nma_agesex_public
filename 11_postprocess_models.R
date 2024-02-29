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
names(res1) <- str_replace(names(res1), "^hba1c_", "hba1c_agesex_")

mdl_names <- list.files("FromVM/hba1c_nocovariates/", patt = "Rds$")
res2 <- map(mdl_names, ~ readRDS(paste0("FromVM/hba1c_nocovariates/", .x)))
names(res2) <- mdl_names %>% str_sub(1, -5)
names(res2) <- str_replace(names(res2), "no_inter", "nointer")

mdl_names <- list.files("FromVM/mace_agesex/", patt = "Rds$")
res3 <-  map(mdl_names, ~ readRDS(paste0("FromVM/mace_agesex/", .x)))
names(res3) <- mdl_names %>% str_sub(1, -5)
res4 <- readRDS("Scratch_data/mace_nointer.Rds")
names(res4) <- c("fixed", "random")
names(res4) <- paste0(names(res4), "_mace_nointer_main")
res <- c(res1, res2, res3, res4)
rm(res1, res2, res3, res4)
beta <- map(res, ~ summary(.x$stanfit)$summary %>% 
  as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "tosep")
names(beta) <- str_to_lower(names(beta))

divergent <- map(res, ~ get_sampler_params(.x$stanfit, inc_warmup = FALSE))
divergent <- map(divergent, ~ map_int(.x, ~ sum(.x[, "divergent__"])))
divergent <- bind_rows(divergent, .id = "tosep") %>% 
   write_csv("Outputs/divergent_transition_count.csv")
## no divergent transitions
rm(res)

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
write_lines("", "Outputs/trials_with_high_rhats.txt")
write_csv(beta, "Outputs/betas_meta_analysis.csv")







