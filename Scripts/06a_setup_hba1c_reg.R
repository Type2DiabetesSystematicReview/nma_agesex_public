library(tidyverse)
source("Scripts//common_functions/Scripts/misc.R")
source("Scripts/00_functions.R")

exclusions <- read_csv("Data/exclusions_update.csv")

a <- expand_grid(data_lvl = c("aggipd", "ipd"),
                 fe_re = c("fixed", "random"),
                 nwork = c("mono", "dual", "triple")) 

a2 <- a %>% 
  filter(data_lvl == "aggipd", fe_re == "random") %>% 
  mutate(f_mdl = "f8", , durn = "ge12")
a3 <- a %>% 
  filter(data_lvl == "aggipd", fe_re == "random") %>% 
  mutate(f_mdl = "f4", durn = "ge26")
a <- a %>% 
  mutate(f_mdl = "f4", durn = "ge12")
a <- bind_rows(a,
               a2,
               a3) %>% 
  arrange(data_lvl, nwork)

a4 <- a %>% 
  filter(data_lvl == "aggipd", fe_re == "random", durn == "ge12", f_mdl == "f4") %>% 
  mutate(f_mdl = "f1")

a5 <- a %>% 
  filter(data_lvl == "aggipd", fe_re == "fixed", durn == "ge12", f_mdl == "f4") %>% 
  mutate(f_mdl = "f1")

a <- bind_rows(a,
               a4,
               a5)


a <- a %>% 
  mutate(model_n = seq_along(a$data_lvl)) %>% 
  select(model_n, everything()) %>% 
  mutate(model_n = paste0("m", str_pad(model_n, width = 2, pad = "0")))

write_csv(a, file = "Outputs/model_order.csv")

dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
exclude <- dropdisconnect
exclude <- tibble(exclusion_reason2 = "Disconnected from network.",
                  trial_id = exclude)
exclusions <- ExcludeRun(exclude)
