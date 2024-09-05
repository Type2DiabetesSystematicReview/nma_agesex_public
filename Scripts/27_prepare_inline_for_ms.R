# 27_prepare_inline_forms
library(tidyverse)
source("Scripts/00_functions.R")
hb <- readRDS("Scratch_data/pre_hba1c_results_ms.Rds")
mc <- readRDS("Scratch_data/pre_mace_results_ms.Rds")

## Create summary statistics ----
hba1c_lst <- hb %>% 
  mutate(network = as.character(network)) %>% 
  mutate(across(c(mean, se_mean, x2_5_percent:x97_5_percent), ~ round(.x, 2) %>% formatC(digits = 2, format = "f")),
         est = mean,
         rngto = paste0(x2_5_percent, " to ", x97_5_percent),
         rngdash = paste0(x2_5_percent, "-", x97_5_percent),
         flto = paste0(est, " (95% CI ", rngto, ")"),
         fldash = paste0(est, " (95% CI ", rngdash, ")"),
         fltonopar = paste0(est, "; 95% CI ", rngto, ""),
         fldashnopar = paste0(est, "; 95% CI ", rngdash, "")) %>% 
  select(fixedrand, network, covariate, trtclass, est, rngto, rngdash, flto, fldash, fltonopar, fldashnopar) 
mace_lst <- mc %>% 
  filter(nct_id == "none") %>% 
  mutate(covariate = if_else(covariate == "Male sex", "male", "age30")) %>% 
  mutate(across(c(mean, se_mean, x2_5_percent:x97_5_percent), ~ round(.x %>% exp(), 2) %>% formatC(digits = 2, format = "f")),
         est = mean,
         rngto = paste0(x2_5_percent, " to ", x97_5_percent),
         rngdash = paste0(x2_5_percent, "-", x97_5_percent),
         flto = paste0(est, " (95% CI ", rngto, ")"),
         fldash = paste0(est, " (95% CI ", rngdash, ")"),
         fltonopar = paste0(est, "; 95% CI ", rngto, ""),
         fldashnopar = paste0(est, "; 95% CI ", rngdash, "")) %>% 
  select(fixedrand, sg, covariate, trtclass, est, rngto, rngdash, flto, fldash, fltonopar, fldashnopar)

## Convert data to a list format ----
hba1c_lst <- hba1c_lst %>% 
  nest(data = c("est", "rngto", "rngdash", "flto", "fldash", "fltonopar", "fldashnopar"))
hba1c_lst$data <- map(hba1c_lst$data, as.list)
mace_lst <- mace_lst %>% 
  nest(data = c("est", "rngto", "rngdash", "flto", "fldash", "fltonopar", "fldashnopar"))
mace_lst$data <- map(mace_lst$data, as.list)

## Convert whole structure into nested list format ----
hba1c_lst <- nest_dataframe(df = hba1c_lst, columns = c("fixedrand", "network", "covariate", "trtclass"), data_col = "data")
hba1c_res <- map(#fixedrand
  hba1c_lst, ~ 
                  #network
             map(.x, ~ 
                   #covariate
                   map(.x, ~ 
                         #trtclass
                         map(.x, ~ 
                               #estimates
                               .x[[1]]))))
saveRDS(hba1c_res, "Scratch_data/hba1c_results_ms.Rds")


mace_lst <- nest_dataframe(df = mace_lst, columns = c("fixedrand", "sg","covariate", "trtclass"), data_col = "data")
mace_res <- map(#fixedrand
  mace_lst, ~ 
    # subgroup
    map(.x, ~ 
          #covariate
          map(.x, ~ 
                #trtclass
                map(.x, ~ 
                      #estimates
                      .x[[1]]))))

saveRDS(mace_res, "Scratch_data/mace_results_ms.Rds")

## Mace relative results (figure 3)
macef3 <- readRDS("Scratch_data/pre_mace_relative_results_ms.Rds")
macef3 <- macef3 %>% 
  mutate(sex = if_else(male ==1, "male", "female"),
         age = paste0("age", age)) %>% 
  mutate(across(c(m, s, q2_5, q97_5), ~ round(.x %>% exp(), 2) %>% formatC(digits = 2, format = "f")),
         est = m,
         rngto = paste0(q2_5, " to ", q97_5),
         rngdash = paste0(q2_5, "-", q97_5),
         flto = paste0(est, " (95% CI ", rngto, ")"),
         fldash = paste0(est, " (95% CI ", rngdash, ")"),
         fltonopar = paste0(est, "; 95% CI ", rngto, ""),
         fldashnopar = paste0(est, "; 95% CI ", rngdash, ""))
macef3 <- macef3 %>% 
  select(trtclass, sex, age, est:fldashnopar) %>% 
  nest(data = est:fldashnopar)
macef3_lst <- nest_dataframe(macef3, c("trtclass", "sex", "age"), data_col = "data")
macef3_res <- map(#trtclass
  macef3_lst, ~ 
    map(.x, ~ 
          map(.x, ~ 
                .x[[1]])))
saveRDS(macef3_res, "Scratch_data/mace_relative_results_ms.Rds")
