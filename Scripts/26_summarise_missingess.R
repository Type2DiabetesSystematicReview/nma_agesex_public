## Summarise missingness
library(tidyverse)

## Missing Hba1c ----
res <- readRDS("Scratch_data/ipd_consolidated.Rds")
dgn <- res$age_sex_model_diags
cfs <- res$age_sex_model_coefs
rm(res)

exclude <- read_csv("Data/exclusions_update.csv")
dgn <- dgn %>% 
  semi_join(exclude %>% filter(!exclude == 1L) %>% select(nct_id = trial_id))

cfs <- cfs %>% 
  semi_join(dgn) %>% 
  mutate(term2 = term %>% 
           str_remove("arm_f_[0-9]_uaa[0-9]*[a-z]*") %>% 
           str_remove("arm_f_[0-9]_ipd[0-9]*") %>% 
           str_remove("arm_f_[0-9]_ab[0-9]*") %>% 
           str_remove("arm_f_2_BI 10773 10mg|arm_f_2_BI 10773 25mg"))
cfs <- cfs %>%
  distinct(nct_id, models, term2) %>% 
  arrange(term2) %>% 
  group_by(nct_id, models) %>% 
  summarise(terms = paste(term2, collapse = ";")) %>% 
  ungroup()
cfs %>% count(models, terms)
models <- list(b1 = c("base", "sex", "age"),
               f1 = c("base", "final"),
               f2 = c("base", "final", "age"),
               f3 = c("base", "final", "sex"),
               f4 = c("base", "final", "age", "sex"))

dgn <- dgn %>% 
  filter(models %in% c("b1", "f1", "f2", "f3", "f4")) %>% 
  select(nct_id, nct_id2, models, nobs) %>% 
  distinct()
dgn <- dgn %>% 
  spread(models, nobs)
dgn <- dgn %>% 
  group_by(nct_id) %>% 
  summarise((across(b1:f4, sum))) %>% 
  ungroup()

dgn$sameall <- dgn %>% 
  select(b1:f4) %>% 
  apply(MARGIN = 1, function(x) sum(!duplicated(x)))
dgn$samef <- dgn %>% 
  select(f1:f4) %>% 
  apply(MARGIN = 1, function(x) sum(!duplicated(x)))
sameall <- dgn %>% 
  filter(sameall ==1)
dgn <- dgn %>% 
  filter(!sameall ==1)
samef <- dgn %>% 
  filter(samef ==1)
dgn <- dgn %>% 
  filter(!samef ==1)
## none with missing sex data
msng_sex <- dgn %>% 
  filter(f1 == f2,
         f1 != f3)
## all remaining have missing age, data calculate missing final by subtracting b1 from f4
msng_age <- dgn %>% 
  filter(!f1 == f2,
         f1 == f3)
## f2 has same informaion as f4. Drop
msng_age %>% 
  filter(!f2 == f4)
msng_age <- msng_age %>% 
  select(-f4) 
msng_age <- msng_age %>% 
  mutate(age_mis = f1 - f2,
         sex_mis = f1 - f3,
         last_hba1c_mis = b1-f2) %>% 
  mutate(age_mis_p = 100*age_mis/b1,
         last_hba1c_mis_p = 100*last_hba1c_mis/b1)
table(msng_age$age_mis_p >=1)
table(msng_age$last_hba1c_mis_p <1)
msng_age$last_hba1c_mis_p [msng_age$last_hba1c_mis_p >=1 ] %>% 
  round(2)
rm(list = ls())
## Missing MACE ----
diag <- read_csv("Data/agesexhba1cmaceupdate_9492/diag_cox.csv")
diag <- diag %>% 
  filter(models %in% paste0("f", 2:4)) %>% 
  mutate(models = case_match(models,
                             "f2" ~ "agesex",
                             "f3" ~ "age",
                             "f4" ~ "sex")) %>% 
  select(nct_id, models, nobs) %>% 
  spread(models, nobs)
diag <- diag %>% 
  mutate(misage = sex - age,
         misage_p = round(100*misage/sex, 2))
## Summary
# Of 103 HbA1c trials, none had missing sex data. Twelve had missing age data; of these nine had missing age data for <1% participants, with the remaining three having missing age data for
# 1.29%, 1.34% and 1.95% respectively. Nine had no HbA1c at follow-up; of these six had missing HbA1c data for less than 1% of participants, with the  remaining having missing HbA1c data for
# 1.11%, 1.17% and 4.81% respectively. The trials/participants with missing HbA1c at follow-up were a subset of the trials/participants with missing age data.
# Of the six IPD MACE trials three had missing age data; there were three (0.06%), two (0.02%) and one (0.01%) of participants with missing age data respectively.
# Missing age data was due to the trial sponsors redaction algorithm.