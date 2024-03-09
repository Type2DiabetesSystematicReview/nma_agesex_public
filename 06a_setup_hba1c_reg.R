library(tidyverse)
source("../common_functions/Scripts/misc.R")
a <- expand_grid(data_lvl = c("aggipd", "ipd"),
                 fe_re = c("fixed", "random"),
                 nwork = c("mono", "dual", "triple")) 
a <- a %>% 
  mutate(model_n = seq_along(a$data_lvl)) %>% 
  select(model_n, everything())
write_csv(a, file = "Outputs/model_order.csv")


dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
exclude <- dropdisconnect

ipd_nct <- bind_rows(
  read_csv("../from_vivli/Data/agesexhba1c_6115/hba1c_base_change_overall.csv"),
  read.csv("../from_gsk/Data/agesex/hba1c_base_change_overall.csv"),
  read.csv("../from_vivli/Data/agesexhba1c_8697/hba1c_base_change_overall.csv")) %>% 
  pull(nct_id) %>% 
  unique()

exclude <- tibble(reason = "Disconnected from network.",
                  trials = length(exclude),
                  trials_ipd = length(intersect(exclude, ipd_nct)),
                  nct_ids = exclude %>% paste(collapse = ";"),
                  nct_ids_ipd = intersect(exclude, ipd_nct) %>% paste(collapse = ";"))
write_tsv(exclude, "Outputs/Trial_exclusion_during_cleaning.txt", append = TRUE)