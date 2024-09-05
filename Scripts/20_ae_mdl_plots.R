# 20_ae_mdl_plts
library(tidyverse)

## run this to obtain whoatc
source("Scripts/04b_arm_label_reverse.R")

## read in data to label trials ----
tot <- readRDS("Scratch_data/for_ae_regression.Rds")
cfs <- tot %>%  
  unnest(cfs) %>% 
  distinct(nct_id, trtcls5)
nct_id_lkp <- cfs %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  distinct(nct_id, trtcls5, .keep_all = TRUE)

## relabel outputs for subsequent plotting ----
beta <- read_csv("Outputs/betas_meta_analysis_ae.csv")
beta <- beta %>% 
  separate(tosep, into = c("ae",
                           "outcome",
                           "modeltype"),
           sep = "_", remove = FALSE)  
beta <- beta %>% 
  mutate(params = str_replace(params, "age10", "age") %>% 
           str_replace("sexTRUE", "male"))
## additional relabelling
beta_age_sex <- beta %>% 
  filter(params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:") %>% 
  mutate(covariate2 = if_else(covariate %in% c("age", "male"), covariate, trtclass),
         trtclass2 = if_else(covariate %in% c("age", "male"), trtclass, covariate)) %>% 
  select(-covariate, -trtclass) %>% 
  rename(covariate = covariate2, trtclass = trtclass2)

beta_age_sex <- beta_age_sex %>% 
  inner_join(who_atc %>% 
               select(trtclass = `ATC code`,
                      cls = `ATC level name`) %>% 
               distinct()) %>% 
  mutate(cls = paste0(trtclass, ":", cls))

beta_age_sex <- beta_age_sex %>% 
  janitor::clean_names()

beta_age_sex <- beta_age_sex %>% 
  mutate(
    across(mean:x97_5_percent, ~ case_when(
      covariate %in% c("age") ~ .x*3,
      TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age"), "age30", covariate))

interaeplot <- ggplot(beta_age_sex %>% 
                           filter(trtclass %in% paste0("A10B", c("A", "B", "H", "J", "K"))) %>% 
                            mutate(covariate = factor(covariate,
                                                     levels = c("age30", "male"),
                                                     labels = c("Age per 30 years",
                                                                "Male sex")),
                                   outcome_f = factor(outcome,
                                                      levels = c("sae", "gi", "hypog", "uti"),
                                                      labels = c("SAEs", "Gastrointenstinal",
                                                                 "Hypoglycaemia",
                                                                 "UTIs")),
                                   modeltype_f = factor(modeltype,
                                                      levels = c("pois", "nb"),
                                                      labels = c("Quasi-Poisson", "Negative binomial"))) %>% 
                        mutate(cls_shrt = case_when(
                          str_detect(cls, "\\(") ~ str_extract(cls, "\\(.*\\)") %>% str_sub(2, -2),
                          TRUE ~ cls)), 
                         aes(x = cls_shrt, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                             colour = outcome_f, 
                             linetype = modeltype_f)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_wrap(~covariate, ncol = 1) + 
  scale_x_discrete("",limits = rev) +
  coord_flip(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Log rate ratio") +
  scale_linetype("") +
  scale_color_discrete("")
interaeplot
saveRDS(interaeplot, "Scratch_data/interaeplot.Rds")
pdf("Outputs/ae_nter_with_pois.pdf", height = 10, width = 20)
interaeplot
dev.off()

tiff("Outputs/ae_nter_with_pois.tiff", height = 10, width = 20, res = 300, unit = "in", compression = "lzw")
interaeplot
dev.off()
