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

mace_nms <- list(mace_main = "mace_1.Rds", 
     mace_age = "mace_2.Rds", 
     mace_sex = "mace_3.Rds")
mace_nms <- list(mace_main_f1 = "mace_f1_1.Rds", 
                 mace_age_f1 = "mace_f1_2.Rds", 
                 mace_sex_f1 = "mace_f1_3.Rds", 
                 mace_main_f5 = "mace_f5_1.Rds", 
                 mace_age_f5 = "mace_f5_2.Rds",
                 mace_main_f5 = "mace_f5_3.Rds")

mace <- map(mace_nms, ~ readRDS(paste0("Scratch_data/", .x)))
mace <- transpose(mace)
mace <- c(fe = mace$fe,
          re = mace$re)
res <- c(res1, res2, mace)
rm(res1, res2, mace)
beta <- map(res, ~ summary(.x$stanfit)$summary %>% 
  as_tibble(rownames = "params")) %>% 
  bind_rows(.id = "tosep")
names(beta) <- str_to_lower(names(beta))
divergent <- map(res, ~ get_sampler_params(.x$stanfit, inc_warmup = FALSE))
divergent <- map(divergent, ~ map_int(.x, ~ sum(.x[, "divergent__"])))
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

## relabel outputs for subsequent plotting ----
hba1c <- beta %>% 
  filter(str_detect(tosep, "^hba1c")) %>% 
  separate(tosep, into = c("outcome",
                           "mainorinter",
                           "modelnum",
                           "datalevel",
                           "fixedrand",
                           "network"),
           sep = "_")  %>% 
  mutate(sg = "main")
mace <- beta %>% 
  filter(str_detect(tosep, "mace")) %>% 
  mutate(modelnum = paste0("m", cumsum(!duplicated(tosep) + 1)),
         network = "triple",
         # mainorinter = "agesex",
         datalevel = "aggipd") %>% 
  separate(tosep, into = c("fixedrand",
                           "outcome",
                           "sg",
                           "mainorinter"),
           sep = "_|\\.")  %>% 
  mutate(mainorinter = if_else(mainorinter == "f5", "agesex", "nointer"),
         fixedrand = if_else(fixedrand == "fe", "fixed", "random"))
beta <- bind_rows(hba1c, mace)
rm(hba1c, mace)
write_csv(beta, "Outputs/betas_meta_analysis.csv")





