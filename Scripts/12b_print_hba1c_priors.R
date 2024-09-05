#12b_print_hba1c_priors
library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
priors <- read_csv("Scratch_data/priors_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
source("Scripts/00_functions.R")

priors <- priors %>% 
  filter(str_detect(tosep, "^m")) %>% 
  separate(tosep, into = c("modelnum",
                           "datalevel",
                           "fixedrand",
                           "network",
                           "f_model",
                           "duration"),
           sep = "_", remove = FALSE)  %>% 
  mutate(outcome = "hba1c")

priors <- priors %>% 
  mutate(datalevel = factor(datalevel,
                            levels = c("aggipd", "ipd", "cnvrtagged"),
                            labels = c("All data", "IPD only", "IPD collapsed")),
         mygrp = paste(network, datalevel, fixedrand, f_model, duration, sep = ", "),
         mygrp = mygrp %>% 
           # str_remove("All data ")  %>% 
           # str_to_sentence() %>% 
           str_replace(" ge12", ">= 12 weeks") %>% 
           str_replace("ge26", ">= 26 weeks") %>% 
           str_replace("f4", "LOCF") %>% 
           str_replace("f8", "BOCF") %>% 
           str_replace("Ipd", "IPD") %>% 
           str_trim()) %>% 
  select(mygrp, `A: Intercept`:`D: Heterogeneity`)
write_csv(priors, "Outputs/priors_meta_analysis_hba1c.csv", na = "")
