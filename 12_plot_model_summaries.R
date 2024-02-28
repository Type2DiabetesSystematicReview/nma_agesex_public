library(tidyverse)

beta <- read_csv("Outputs/betas_meta_analysis.csv")
if(sessionInfo()$platform == "x86_64-pc-linux-gnu (64-bit)") {
  whoatc <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx", sheet = 1) 
} else {
  whoatc <- readxl::read_excel("../../../Medications_resources/WHO_ATC/2018 ATC index with DDDs.xlsx", sheet = 1)
}
whoatc <- whoatc %>% 
  select(trtclass = `ATC code`,
         cls = `ATC level name`) %>% 
  distinct()

## interaction plots hba1c ----
beta_age_sex <- beta %>% 
  filter(params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:")
beta_age_sex <- beta_age_sex %>% 
  inner_join(whoatc) %>% 
  mutate(cls = paste0(trtclass, ":", cls))
beta_age_sex <- beta_age_sex %>% 
  janitor::clean_names() %>% 
  mutate(network = factor(network, levels = c("mono", "dual", "triple")))
beta_age_sex <- beta_age_sex %>% 
  mutate(myalpha = case_when(
    (trtclass == "A10BX") ~ 0.2,
    (trtclass == "A10A" & network == "mono") ~ 0.2,
    (trtclass %in% c("A10BX", "A10BA") & network == "triple") ~ 0.2,
    TRUE ~ 1),
    across(mean:x97_5_percent, ~ case_when(
      covariate == "age10" ~ .x*3,
      covariate == "age15" ~ .x*2,
      TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age10", "age15"), "age30", covariate))

plotbeta <- ggplot(beta_age_sex %>% 
                     filter(!is.na(network), outcome == "hba1c") %>% 
                     mutate(datalevel = factor(datalevel,
                                               levels = c("aggipd", "ipd"),
                                               labels = c("All data", "IPD only"))), 
                   aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                       colour = fixedrand, 
                       shape = datalevel, 
                       linetype = datalevel,
                       alpha = myalpha)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ network) + 
  # scale_y_continuous(limits = c(-1, 1)) +
  scale_x_discrete(limits = rev) +
  scale_alpha_identity() +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("HbA1c (%)")
plotbeta

plotbetapaper <- ggplot(beta_age_sex %>% 
                          filter(!is.na(network), outcome == "hba1c",
                                 datalevel == "aggipd",
                                 trtclass %in% paste0("A10B", c("A", "B", "H", "J", "K")),
                                 !(network == "triple" & trtclass == "A10BA")) %>% 
                          mutate(datalevel = factor(datalevel,
                                                    levels = c("aggipd", "ipd"),
                                                    labels = c("All data", "IPD only")),
                                 covariate = factor(covariate,
                                                    levels = c("age30", "male"),
                                                    labels = c("Age per 30 years",
                                                               "Male sex")),
                                 network = factor(network,
                                                  levels = c("mono", "dual", "triple"),
                                                  labels = c("Monotherapy", "Dual therapy", "Triple therapy"))), 
                        aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                            colour = fixedrand, 
                            alpha = myalpha)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ network) + 
  # scale_y_continuous(limits = c(-1, 1)) +
  scale_x_discrete("",limits = rev) +
  scale_alpha_identity() +
  scale_color_discrete("") +
  coord_flip(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("HbA1c (%)")
plotbetapaper

## interaction plots mace ----
plotmace <- ggplot(beta_age_sex %>% 
                     filter(outcome == "mace") %>% 
                     mutate(datalevel = factor(datalevel,
                                               levels = c("aggipd", "ipd"),
                                               labels = c("All data", "IPD only")),
                            sg = factor(sg,
                                        levels = c("main", "age", "sex"),
                                        labels = c("None", "Age", "Sex")),
                            covariate = factor(covariate,
                                               levels = c("age30", "male"),
                                               labels = c("Age per 30 years",
                                                          "Male sex"))), 
                   aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                       colour = sg,
                       shape = fixedrand,
                       linetype = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(~ covariate) + 
  scale_x_discrete(limits = rev) +
  scale_color_discrete("Subgroup data") +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Hazard ratio (log-scale)")
plotmace

plotmacepaper <- ggplot(beta_age_sex %>% 
                          filter(outcome == "mace", sg == "main",
                                 trtclass %in% c("A10BH", "A10BJ", "A10BK")) %>% 
                          mutate(covariate = factor(covariate,
                                                    levels = c("age30", "male"),
                                                    labels = c("Age per 30 years",
                                                               "Male sex"))), 
                        aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                            colour = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(~ covariate) + 
  scale_x_discrete("", limits = rev) +
  scale_color_discrete("") +
  coord_flip(ylim = c(-0.7, 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Hazard ratio", 
                     breaks = seq(-0.5, 0.5, 0.25), 
                     labels = round(exp(seq(-0.5, 0.5, 0.25)),2))
plotmacepaper

## main effects hba1c ----
main_hba1c <- beta %>% 
  filter(mainorinter == "nointer", 
         outcome == "hba1c",
         !str_detect(params, "^delta"),
         str_detect(params, "^d")) 
main_hba1c <- main_hba1c %>% 
  mutate(params = str_sub(params, 3, -2)) %>% 
  separate_wider_delim(params, names = c("drug", "dose"), delim = "_d",too_few = "align_start" ) %>% 
  mutate(cls = str_sub(drug, 1, 5))
who_atc <- readxl::read_excel("../2018 ATC index with DDDs.xlsx")
who_lkp_rev1 <- who_atc$`ATC level name`
names(who_lkp_rev1) <- who_atc$`ATC code`
who_lkp <- readRDS("Scratch_data/who_atc_lkp.Rds")
who_lkp_rev2 <- names(who_lkp)
names(who_lkp_rev2) <- who_lkp
rm(who_atc, who_lkp)
who_lkp_rev <- c(who_lkp_rev1, who_lkp_rev2)
who_lkp_rev <- who_lkp_rev[!duplicated(names(who_lkp_rev))]
who_lkp_rev <- who_lkp_rev[!str_detect(who_lkp_rev, "\\|")]
who_lkp_rev <- who_lkp_rev[!duplicated(who_lkp_rev)]
setdiff(main_hba1c$drug, names(who_lkp_rev))
rm(who_lkp_rev1, who_lkp_rev2)
main_hba1c <- main_hba1c %>% 
  mutate(atc = drug,
         drug = who_lkp_rev[drug])
main_hba1c <- main_hba1c %>% 
  separate_wider_delim(atc, names = c("atc2", "todrop"), delim = "_", too_many = "merge", too_few = "align_start") %>% 
  mutate(lbl= if_else(!is.na(dose),
                      paste0(atc2, "-", str_to_sentence(drug), " ", dose),
                      paste0(atc2, "-", str_to_sentence(drug))))
main_hba1c <- main_hba1c %>% 
  janitor::clean_names()
main_hba1c <- main_hba1c %>% 
  mutate(cls_facet = case_when(
    cls %in% c("A10BA", "A10BB", "A10BH", "A10BJ", "A10BK") ~ cls,
    TRUE ~ "Other"
  ))
mainhba1cplot <- ggplot(main_hba1c %>% 
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
mainhba1cplot


pdf("Outputs/hba1c_mace.pdf", height = 10, width = 20)
mainhba1cplot + ggtitle("HbA1c meta-analysis main effects for paper")
plotbetapaper  + ggtitle("HbA1c meta-analysis for paper")
plotmacepaper  + ggtitle("MACE meta-analysis for paper")
plotbeta + ggtitle("HbA1c meta-analysis additional information")
plotmace + ggtitle("MACE meta-analysisadditional information")
dev.off()