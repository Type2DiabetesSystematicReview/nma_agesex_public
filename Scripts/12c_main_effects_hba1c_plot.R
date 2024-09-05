library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
source("Scripts/00_functions.R")

## read in information on which models failed to converge ----
high_rhat <- read_csv("Outputs/nma_results_high_rhats.csv") %>% 
  distinct(tosep)

## main effects hba1c ----
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
hba1c <- hba1c %>% 
  filter(f_model == "f1")
# pull main effects
hba1c <- hba1c %>% 
  filter(!str_detect(params, "^delta"),
         str_detect(params, "^d")) 
hba1c <- hba1c %>% 
  mutate(hba1c = str_sub(params, 3, -2)) %>% 
  separate_wider_delim(params, names = c("drug", "dose"), delim = "_d",too_few = "align_start" ) 
hba1c <- hba1c %>% 
  mutate(drg_lkp = drug %>% 
           str_remove("\\[") %>% 
           str_remove("\\]") %>% 
           str_remove("^d") %>% 
           str_trim())
hba1c <- hba1c %>% 
  mutate(atc = drug,
         drug = who_lkp_rev[drg_lkp])
hba1c <- hba1c %>% 
  separate_wider_delim(atc, names = c("atc2", "todrop"), delim = "_", too_many = "merge", too_few = "align_start") %>% 
  mutate(lbl= if_else(!is.na(dose),
                      paste0(atc2, "-", str_to_sentence(drug), " ", dose),
                      paste0(atc2, "-", str_to_sentence(drug))))
hba1c <- hba1c %>% 
  janitor::clean_names()

mainhba1cplot <- ggplot(hba1c %>% 
                          filter(datalevel == "aggipd") %>% 
                          mutate(network = factor(network,
                                                  levels = c("mono", "dual", "triple"),
                                                  labels = c("Monotherapy", "Dual therapy", "Triple therapy"))),
                        aes(x = lbl, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                            colour = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete("", limits = rev) +
  scale_color_discrete("") +
  theme_bw() +
  scale_y_continuous("HbA1c (%)") +
  coord_flip() +
  facet_wrap( ~ network, scales = "free_y", ncol = 3)
hba1c <- hba1c %>% 
  mutate(grp = network) %>% 
  group_by(grp) %>% 
  nest() %>% 
  ungroup()
mainhba1cplot_lst <- map(hba1c$data, ~ {ggplot(.x %>% 
                                                          filter(datalevel == "aggipd") %>% 
                                                          mutate(network = factor(network,
                                                                                  levels = c("mono", "dual", "triple"),
                                                                                  labels = c("Monotherapy", "Dual therapy", "Triple therapy"))),
                                                        aes(x = lbl, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                                                            colour = fixedrand)) +
    geom_point(position = position_dodge(0.5)) +
    geom_linerange(position = position_dodge(0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_discrete("", limits = rev) +
    scale_color_discrete("") +
    theme_bw() +
    scale_y_continuous("HbA1c (%)") +
    coord_flip() +
    facet_wrap( ~ network, scales = "free_y", ncol = 1)})
saveRDS(mainhba1cplot_lst, "Scratch_data/main_eff_lst.Rds")
