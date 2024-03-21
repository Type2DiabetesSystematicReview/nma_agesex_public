library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
source("Scripts/00_functions.R")
## read in data to label trials ----
tot <- readRDS("Scratch_data/for_mace_regression_inter.Rds")
cfs <- tot$cfs %>% 
  distinct(nct_id, trtcls5)
mace_agg_sex <- tot$mace_agg_sex %>% 
  distinct(nct_id, trtcls5)
mace_agg_age <- tot$mace_agg_age %>% 
  distinct(nct_id, trtcls5)
mace_agg <- tot$mace_agg %>% 
  distinct(nct_id, trtcls5)
nct_id_lkp <- bind_rows(ipd = cfs,
                        sg = mace_agg_age,
                        sg = mace_agg_sex,
                        agg = mace_agg, 
                        .id = "data_lvl") %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  distinct(nct_id, trtcls5, .keep_all = TRUE)
rm(tot, cfs, mace_agg_age, mace_agg_sex, mace_agg)

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
  mutate(sg = "main",
         nct_id = "none")
hba1c <- hba1c %>% 
  mutate(params = str_replace(params, "age10", "age"))

mace <- beta %>% 
  filter(str_detect(tosep, "mace")) %>% 
  mutate(modelnum = paste0("m", cumsum(!duplicated(tosep))+1),
         network = "triple",
         # mainorinter = "agesex",
         datalevel = "aggipd") %>% 
  separate(tosep, into = c("fixedrand",
                           "outcome",
                           "mainorinter",
                           "sg",
                           "nct_id"),
           sep = "_|\\.", remove = FALSE)  %>% 
  mutate(nct_id = if_else(is.na(nct_id), "none", nct_id))

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
      covariate %in% c("age") ~ .x*30,
      TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age"), "age30", covariate))

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
interhba1cplot

## interaction plots mace ----
# "sens", "sens2"
# "None (one SGLT2 IPD trial excluded)",
# "None (one DPP-4 IPD trial excluded)"
intermaceplot <- ggplot(beta_age_sex %>% 
                          filter(outcome == "mace", sg == "sex",
                                 nct_id == "none",
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

## sensitivity analyses for MACE ----
nct_id_lkp <- nct_id_lkp %>% 
  mutate(trtcls5 = case_when(
    trtcls5 == "A10BH" ~ "DPP-4",
    trtcls5 == "A10BJ" ~ "GLP-1",
    trtcls5 == "A10BK" ~ "SGLT-2",
  ))

macesens <- beta_age_sex %>% 
  mutate(issg = str_sub(nct_id, 12),
         nct_id = str_sub(nct_id, 1, 11)) %>% 
  filter(!(sg == "main" & issg == "issg"),
         fixedrand == "fixed",
         outcome == "mace",
         trtclass %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  mutate(trtclass = case_when(
    trtclass == "A10BH" ~ "DPP-4",
    trtclass == "A10BJ" ~ "GLP-1",
    trtclass == "A10BK" ~ "SGLT-2",
  )) %>% 
  left_join(nct_id_lkp) %>% 
  mutate(sg_lbl = case_when(
    sg %in% c("age", "sex") & issg == "issg" ~ "Yes, ipdreplace",
    sg %in% c("age", "sex") ~ "Yes, noipdreplace",
    TRUE ~ "No"),
         sg_lbl = factor(sg_lbl, levels = c("Yes, ipdreplace", "Yes, noipdreplace", "No")),
         data_lvl = if_else(is.na(data_lvl), 
                            "na", data_lvl),
         data_lvl = factor(data_lvl,
                           levels = c("na",
                                      "agg",
                                      "sg",
                                      "ipd"),
                           labels = c("Non-applicable",
                                      "Aggregate only",
                                      "Subgroup",
                                      "IPD")),
         xvar = case_when(
          nct_id == "none" ~ "All trials included",
          TRUE ~ paste0(trtcls5, ": ", nct_id)
         ))
macesensage <- macesens %>% 
  filter(sg %in% c("main", "age"),
         covariate == "age30")
macesenssex <- macesens %>% 
  filter(sg %in% c("main", "sex"),
         covariate == "male")
## examine impact of using sex-subgroup on age analysis
macesensage_contra <- macesens %>% 
  filter(sg %in% c("main", "sex"),
         covariate == "age30")
## examine impact of using age-subgroup on sex analysis
macesenssex_contra  <- macesens %>% 
  filter(sg %in% c("main", "age"),
         covariate == "male")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

macesensageplt <- ggplot(macesensage,
                         aes(x = xvar,
                             y = mean, 
                             ymin = x2_5_percent, 
                             ymax = x97_5_percent, 
                             shape = sg_lbl,
                             linetype = sg_lbl,
                             colour = data_lvl)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_wrap(~trtclass) + 
  scale_x_discrete("", limits = rev) +
  scale_color_manual("Trial data type", 
                     values = cbbPalette[1:4]) +
  coord_flip(ylim = c(-0.7, 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal3() +
  scale_y_continuous("Hazard ratio", 
                     breaks = seq(-0.5, 0.5, 0.25), 
                     labels = round(exp(seq(-0.5, 0.5, 0.25)),2)) +
  scale_shape("Subgroup data included") +
  scale_linetype("Subgroup data included")
macesensageplt
macesenssexplt <- macesensageplt %+% macesenssex
macesensageplt_contra <- macesensageplt %+% macesensage_contra
macesenssexplt_contra <- macesenssexplt %+% macesenssex_contra
ipdwsggot <- c("NCT00968708", "NCT02465515", "NCT01131676")

macesensageplt_paper <- macesensageplt %+% 
  (macesensage %>% 
     filter( (nct_id %in% ipdwsggot & trtclass == trtcls5) |
     nct_id == "none") %>% 
     mutate(sg_lbl = case_when(
       sg_lbl == "Yes, ipdreplace" ~ "IPD, AGG and SG data, including SG for this trial",
       sg_lbl == "Yes, noipdreplace" ~ "IPD, AGG and SG data, SG only for dropped/downgraded trials",
       sg_lbl == "No" ~ "IPD and AGG data only"),
       xvar = if_else(nct_id == "none", 
                      xvar,
                      paste0(xvar, " dropped/downgraded")))) +
  aes(colour = sg_lbl, linetype = NULL, shape = NULL) +
  theme_minimal4() 
macesensageplt_paper
macesenssexplt_paper <- macesenssexplt %+% 
  (macesenssex %>% 
     filter( (nct_id %in% ipdwsggot & trtclass == trtcls5) |
               nct_id == "none") %>% 
     mutate(sg_lbl = case_when(
       sg_lbl == "Yes, ipdreplace" ~ "IPD, AGG and SG data, no data (IPD or SG) for this trial",
       sg_lbl == "Yes, noipdreplace" ~ "IPD, AGG and SG data, SG only for dropped/downgraded trials",
       sg_lbl == "No" ~ "IPD and AGG data only"),
       xvar = if_else(nct_id == "none", 
                      xvar,
                      paste0(xvar, " dropped/downgraded")))) +
  aes(colour = sg_lbl, linetype = NULL, shape = NULL) +
  theme_minimal4() 
macesenssexplt_paper
pdf("Outputs/sens_onetrial.pdf", width = 15, height = 8)
macesenssexplt_paper + ggtitle("Sex-treatment interaction with/without sex subgroup data")
macesensageplt_paper + ggtitle("Age-treatment interaction with/without age subgroup data")
dev.off()

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

regplots <- list(interhba1cplot = interhba1cplot,
                 interhba1cplotappen = interhba1cplotappen,
                 intermaceplot = intermaceplot,
                 macesensageplt = macesensageplt,
                 macesenssexplt = macesenssexplt,
                 macesensageplt_contra = macesensageplt_contra,
                 macesenssexplt_contra = macesenssexplt_contra,
                 mainhba1cplot = mainhba1cplot,
                 mainmacecplot = mainmacecplot)
names(regplots) <- names(regplots) %>% str_remove("plot$") %>% 
  str_replace("^main", "m") %>% 
  str_replace("^inter", "nt") %>% 
  str_replace("^hba1c", "hb")%>% 
  str_replace("^mace", "mc") %>% 
  str_replace("plot", "plt")
pdf("Outputs/hba1c_mace.pdf", height = 10, width = 20)
regplots
dev.off()
saveRDS(regplots, "Scratch_data/regplots.Rds")
