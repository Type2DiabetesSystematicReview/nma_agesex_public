library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
## relabel outputs for subsequent plotting ----
hba1c <- beta %>% 
  filter(str_detect(tosep, "^hba1c")) %>% 
  separate(tosep, into = c("outcome",
                           "mainorinter",
                           "modelnum",
                           "datalevel",
                           "fixedrand",
                           "network"),
           sep = "_", remove = FALSE)  %>% 
  mutate(sg = "main")
mace <- beta %>% 
  filter(str_detect(tosep, "mace")) %>% 
  mutate(modelnum = paste0("m", cumsum(!duplicated(tosep))+1),
         network = "triple",
         # mainorinter = "agesex",
         datalevel = "aggipd") %>% 
  separate(tosep, into = c("fixedrand",
                           "outcome",
                           "mainorinter",
                           "sg"),
           sep = "_|\\.", remove = FALSE)  

x <- names(mace)
y <- names(hba1c)
setdiff(union(x, y), intersect(x, y))
cmpr <- map2(hba1c %>% select(outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg),
     mace  %>% select(outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg), ~ {
       list(hba1conly = setdiff(.x, .y),
            maceonly  = setdiff(.y, .x))
     })
beta <- bind_rows(hba1c, mace)
rm(hba1c, mace, x, y, cmpr)
write_csv(beta %>% select(tosep:network, sg) %>% distinct(), "Scratch_data/modelname_content_lkp.csv")

## interaction plots hba1c ----
beta_age_sex <- beta %>% 
  filter(mainorinter == "agesex",
         datalevel == "aggipd",
         # sg == "main",
         params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:")
beta_age_sex <- beta_age_sex %>% 
  inner_join(who_atc %>% 
               select(trtclass = `ATC code`,
                      cls = `ATC level name`) %>% 
               distinct()) %>% 
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

interhba1cplotappen <- ggplot(beta_age_sex %>% 
                     filter(outcome == "hba1c") %>% 
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
# interhba1cplotappen

interhba1cplot <- ggplot(beta_age_sex %>% 
                          filter(outcome == "hba1c",
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
# interhba1cplot

## interaction plots mace ----
intermaceplotappen <- ggplot(beta_age_sex %>% 
                     filter(outcome == "mace") %>% 
                     mutate(datalevel = factor(datalevel,
                                               levels = c("aggipd", "ipd"),
                                               labels = c("All data", "IPD only")),
                            sg = factor(sg,
                                        levels = c("main", "age", "sex", "sens", "sens2"),
                                        labels = c("None (full set)", "Age", "Sex", 
                                                   "None (one SGLT2 IPD trial excluded)",
                                                   "None (one DPP-4 IPD trial excluded)")),
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
  scale_color_discrete("Subgroup data/ Sensitivity analysis") +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Hazard ratio (log-scale)")
# intermaceplotappen

intermaceplot <- ggplot(beta_age_sex %>% 
                          filter(outcome == "mace", sg == "main",
                                 trtclass %in% c("A10BH", "A10BJ", "A10BK")) %>% 
                          mutate(covariate = factor(covariate,
                                                    levels = c("age30", "male"),
                                                    labels = c("Age per 30 years",
                                                               "Male sex")),
                         top = "Triple therapy"), 
                        aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                            colour = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ top) + 
  scale_x_discrete("", limits = rev) +
  scale_color_discrete("") +
  coord_flip(ylim = c(-0.7, 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Hazard ratio", 
                     breaks = seq(-0.5, 0.5, 0.25), 
                     labels = round(exp(seq(-0.5, 0.5, 0.25)),2))


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
# mainhba1cplot

## main effects hba1c ----
main_mace <- beta %>% 
  mutate(sg = if_else(is.na(sg), "main", sg)) %>% 
  filter(mainorinter == "nointer", 
         outcome == "mace",
         sg == "main",
         !str_detect(params, "^delta"),
         str_detect(params, "^d")) 
main_mace <- main_mace %>% 
  mutate(params = str_sub(params, 3, -2)) %>% 
  separate_wider_delim(params, names = c("drug", "dose"), delim = "_", too_few = "align_start", too_many = "merge") %>% 
  mutate(cls = str_sub(drug, 1, 5))

main_mace <- main_mace %>% 
  mutate(atc = who_lkp_mace[drug])
main_mace <- main_mace %>% 
  separate_wider_delim(atc, names = c("atc2", "todrop"), delim = "_", too_many = "merge", too_few = "align_start") %>% 
  mutate(lbl= if_else(!is.na(dose),
                      paste0(atc2, "-", str_to_sentence(drug), " ", dose),
                      paste0(atc2, "-", str_to_sentence(drug))))
main_mace <- main_mace %>% 
  janitor::clean_names() 
main_mace <- main_mace %>% 
  mutate(cls = str_sub(atc2, 1, 5))
main_mace <- main_mace %>% 
  mutate(cls_lbl = who_lkp_rev[cls],
         cls_lbl = paste0(cls, "-", cls_lbl))
mainmacecplot <- ggplot(main_mace, 
                        aes(x = lbl, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, 
                            colour = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_discrete("", limits = rev) +
  scale_color_discrete("") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_y_continuous("Hazard ratio", 
                     breaks = seq(-0.5, 0.5, 0.25), 
                     labels = round(exp(seq(-0.5, 0.5, 0.25)),2)) +
  coord_flip(ylim = c(-0.61, 0.61)) +
  facet_wrap(~cls_lbl, scales = "free_y", ncol = 1)

pdf("Outputs/hba1c_mace.pdf", height = 10, width = 20)
cowplot::plot_grid(interhba1cplot + 
                     scale_color_discrete(guide = "none") + 
                     ggtitle("Hba1c"), 
                   intermaceplot + 
                     ggtitle("MACE") + 
                     scale_color_discrete(guide = "none"), 
                   nrow = 1, rel_widths = c(3, 1.9))
mainhba1cplot + ggtitle("HbA1c main effects meta-analysis for paper")
interhba1cplot  + ggtitle("HbA1c interactions meta-analysis for paper")
interhba1cplotappen + ggtitle("HbA1c interactions meta-analysis for appendix")
mainmacecplot + ggtitle("MACE main effects meta-analysis for paper")
intermaceplot  + ggtitle("MACE interactions meta-analysis for paper")
intermaceplotappen + ggtitle("MACE interactions meta-analysis for appendix")
dev.off()

regplots <- list(interhba1cmaceplot = cowplot::plot_grid(interhba1cplot + 
                                      scale_color_discrete(guide = "none") + 
                                      ggtitle("Hba1c"), 
                                    intermaceplot + 
                                      ggtitle("MACE") + 
                                      scale_color_discrete(guide = "none"), 
                                    nrow = 1, rel_widths = c(3, 1.9)),
                 mainhba1cplot = mainhba1cplot,
                 interhba1cplot = interhba1cplot,
                 interhba1cplotappen = interhba1cplotappen,
                 mainmaceplot = mainmacecplot,
                 intermaceplot = intermaceplot,
                 intermaceplotappen = intermaceplotappen)
names(regplots) <- names(regplots) %>% str_remove("plot$") %>% 
  str_replace("^main", "m") %>% 
  str_replace("^inter", "nt") %>% 
  str_replace("^hba1c", "hb")%>% 
  str_replace("^mace", "mc")
saveRDS(regplots, "Scratch_data/regplots.Rds")
