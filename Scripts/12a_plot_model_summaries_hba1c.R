library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
source("Scripts/00_functions.R")

## read in information on which models failed to converge ----
high_rhat <- read_csv("Outputs/nma_results_high_rhats.csv") %>% 
  distinct(tosep)

## read in data to label trials and count number in each group ----
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
tot$ipd <- map(tot$ipd, ~ .x  %>% 
                 count(nct_id, arm_lvl, trtcls5, trtcls4))
tot$agg <- map(tot$agg, ~ .x %>% 
                 select(nct_id, arm_lvl, trtcls5, trtcls4, n))
tot <- tot %>% 
  select(drug_regime_smpl, agg, ipd) %>% 
  gather("datatype", "data", -drug_regime_smpl) %>% 
  unnest(data)
tot <- tot %>% 
  filter(!nct_id %in% dropdisconnect)
nbyclass <- tot %>% 
  group_by(drug_regime_smpl, trtcls5) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            participants = sum(n)) %>% 
  ungroup()
nbyregime <- tot %>% 
  group_by(drug_regime_smpl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            participants = sum(n)) %>% 
  ungroup()

## relabel outputs for subsequent plotting ----
beta <- beta %>% 
  anti_join(high_rhat)
hba1c <- beta %>% 
  filter(str_detect(tosep, "^m")) %>% 
  separate(tosep, into = c("modelnum",
                           "datalevel",
                           "fixedrand",
                           "network",
                           "f_model",
                           "duration"),
           sep = "_", remove = FALSE)  %>% 
  mutate(outcome = "hba1c")
rm(beta)
hba1c <- hba1c %>% 
  mutate(params = str_replace(params, "age10", "age"))

## Relabel parameters with drug class and harmonise names ----
hba1c <- hba1c %>% 
  filter(params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:")
hba1c <- hba1c %>% 
  inner_join(who_atc %>% 
               select(trtclass = `ATC code`,
                      cls = `ATC level name`) %>% 
               distinct()) 
hba1c <- hba1c %>% 
  mutate(cls_bare = cls,
         cls = paste0(trtclass, ":", cls))
cls_bare_levels <- hba1c %>% 
  distinct(cls, cls_bare) 
hba1c <- hba1c %>% 
  mutate(cls_bare = factor(cls_bare, levels = cls_bare_levels$cls_bare))
rm(cls_bare_levels)
hba1c <- hba1c %>% 
  janitor::clean_names() %>% 
  mutate(network = factor(network, levels = c("mono", "dual", "triple")))
hba1c <- hba1c %>% 
  mutate(myalpha = case_when(
    (trtclass == "A10BX") ~ 0.2,
    (trtclass == "A10A" & network == "mono") ~ 0.2,
    (trtclass %in% c("A10BX", "A10BA") & network == "triple") ~ 0.2,
    TRUE ~ 1),
    across(mean:x97_5_percent, ~ case_when(
      covariate %in% c("age") ~ .x*30,
      TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age"), "age30", covariate))

## Interactions for sensitivity analysis ----
hba1c_sens <- hba1c %>% 
  filter(trtclass %in% paste0("A10B", c("H", "J", "K")),
         !f_model == "f1") %>% 
  mutate(datalevel = factor(datalevel,
                            levels = c("aggipd", "ipd", "cnvrtagged"),
                            labels = c("All data", "IPD only", "IPD collapsed")),
         mygrp = paste(datalevel, fixedrand, f_model, duration),
         mygrp = mygrp %>% 
           str_remove("All data ")  %>% 
           str_to_sentence() %>% 
           str_remove(" ge12") %>% 
           str_replace("ge26", "longer duration") %>% 
           str_replace("f4", "LOCF") %>% 
           str_replace("f8", "BOCF") %>% 
           str_replace("Ipd", "IPD") %>% 
           str_trim())
hba1c_sens_plt <- ggplot(hba1c_sens %>% 
                           mutate(cls_shrt = str_extract(cls_bare, "\\(.*\\)") %>% str_sub(2, -2)), 
                              aes(x = cls_shrt, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                                  colour = mygrp, 
                                  alpha = myalpha)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid( network ~ covariate) + 
  scale_x_discrete(limits = rev) +
  scale_alpha_identity() +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal2() +
  scale_y_continuous("HbA1c (%)") +
  scale_x_discrete("") +
  scale_color_discrete("")

hba1c_main <- hba1c %>% 
  filter(datalevel == "aggipd",
         f_model == "f4",
         duration == "ge12",
         trtclass %in% paste0("A10B", c("A", "B", "H", "J", "K")),
         !(network == "triple" & trtclass == "A10BA")) 
saveRDS(hba1c_main, "Scratch_data/pre_hba1c_results_ms.Rds")

hba1c_main <- hba1c_main %>% 
  inner_join(nbyregime %>% select(network = drug_regime_smpl,
                                 nwork_trials = trials,
                                 nwork_n = participants))
## add class level labels
hba1c_main <- hba1c_main %>% 
  inner_join(nbyclass %>% select(network = drug_regime_smpl,
                                 trtclass = trtcls5,
                                 cls_trials = trials,
                                 cls_n = participants))
hba1c_main <- hba1c_main %>% 
  mutate(nwork_p = cls_n/nwork_n)

hba1c_main <- hba1c_main%>% 
  mutate(datalevel = factor(datalevel,
                            levels = c("aggipd", "ipd"),
                            labels = c("All data", "IPD only")),
         covariate = factor(covariate,
                            levels = c("age30", "male"),
                            labels = c("Age per 30 years",
                                       "Male sex")),
         network2 = case_match(network,
                          "mono" ~ "Monotherapy",
                          "dual" ~ "Dual therapy",
                          "triple" ~ "Triple therapy"))
## Quite complex approach to creating factor labelling so will update if numbers change
hba1c_main <- hba1c_main %>% 
         mutate(network_fac = paste0(
           network2,
           "\n",
           nwork_trials,
           " trials. ",
           formatC(nwork_n, format = "d", big.mark = ","),
           " participants"))
nworkforfac <- hba1c_main %>% 
  distinct(network, network2, network_fac) %>% 
  arrange(factor(network, levels = c("mono", "dual", "triple")))
hba1c_main <- hba1c_main %>% 
  mutate(network_fac = factor(network,
                              levels = nworkforfac$network,
                              labels = nworkforfac$network_fac))

hba1c_main_plt <- ggplot(hba1c_main, 
                         aes(x = cls_bare, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                             colour = factor(fixedrand, levels = c("random", "fixed"), labels = c("Random", "Fixed")), 
                             alpha = myalpha)) +
  geom_point(position = position_dodge(0.5), mapping = aes(size = nwork_p), shape = "square") +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ network_fac) + 
  scale_x_discrete("", limits = rev) +
  scale_alpha_identity() +
  scale_color_discrete("Parameterisation of trial estimates within drug classes", limits = rev) +
  coord_flip(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("HbA1c (%)") +
  scale_size_area(guide = "none", limits = c(0, 1)) +
  theme_minimal5()
hba1c_main_plt
saveRDS(list(main = hba1c_main_plt, sens = hba1c_sens_plt), "Scratch_data/hba1c_plots.Rds")
