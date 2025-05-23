library(tidyverse)
beta <- read_csv("Outputs/betas_meta_analysis.csv")
priors <- read_csv("Scratch_data/priors_meta_analysis.csv")
source("Scripts/04b_arm_label_reverse.R")
source("Scripts/00_functions.R")

## read in data to label trials and get number in each group ----
tot <- readRDS("Scratch_data/for_mace_regression_inter.Rds")
cfs <- tot$cfs %>% 
  distinct(nct_id, trtcls5)
mace_lng <- read_csv("Outputs/manuscript_table1b_machine_readable.csv")
mace_lng <- mace_lng %>% 
  filter(!nct_id == "UMIN000018395")
mace_lng <- mace_lng %>% 
  mutate(trtcls5 = str_sub(dc, 1, 5)) %>% 
  select(nct_id, trtcls5, participants)
mace_tot <- mace_lng %>% 
  group_by(trtcls5) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            participants = sum(participants)) %>% 
  ungroup()
nbyclass <- mace_tot
nbyregime <- nbyclass %>% 
  summarise(trials = sum(trials),
            participants = sum(participants))
nbyclass <- bind_cols(nbyregime %>% select(nwork_trials = trials,
                                           nwork_n = participants),
                      nbyclass %>% select(trtclass = trtcls5,
                                          cls_trials = trials,
                                          cls_n = participants)) %>% 
  mutate(nwork_p = cls_n/nwork_n)

mace_agg_sex <- tot$mace_agg_sex %>% 
  distinct(nct_id, trtcls5)
mace_agg_age <- tot$mace_agg_age %>% 
  distinct(nct_id, trtcls5)
mace_agg <- tot$mace_agg %>% 
  distinct(nct_id, trtcls5)

write_csv(nbyclass, "Scratch_data/mace_n_by_class.csv")
nct_id_lkp <- bind_rows(ipd = cfs,
                        sg = mace_agg_age,
                        sg = mace_agg_sex,
                        agg = mace_agg, 
                        .id = "data_lvl") %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  distinct(nct_id, trtcls5, .keep_all = TRUE)
rm(tot, cfs, mace_agg_age, mace_agg_sex, mace_agg, agg, psd)

## count number of trials and participants in each class ----

## relabel outputs for subsequent plotting ----
mace <- beta %>% 
  filter(str_detect(tosep, "mace")) 
mace <- mace %>% 
  mutate(tosep = case_when(
    tosep == "fixed_mace_agesex_noipd" ~ "fixed_mace_agesex_main_noipd",
    tosep == "random_mace_agesex_noipd" ~ "random_mace_agesex_main_noipd",
    TRUE ~ tosep
  ))
mace <- mace %>% 
  mutate(modelnum = paste0("m", cumsum(!duplicated(tosep))+1),
         network = "triple") %>% 
  separate(tosep, into = c("fixedrand",
                           "outcome",
                           "mainorinter",
                           "sg",
                           "nct_id"),
           sep = "_|\\.", remove = FALSE)  %>% 
  mutate(nct_id = if_else(is.na(nct_id), "none", nct_id),
         sg = if_else(mainorinter == "nointer", "main", sg),
         datalevel = if_else(nct_id == "noipd", "noipd", "aggipd"))  

## Mace priors tabulation ----
mace_priors <- priors %>% 
  filter(str_detect(tosep, "mace")) %>% 
  mutate(modelnum = paste0("m", cumsum(!duplicated(tosep))+1),
         network = "triple") %>% 
  separate(tosep, into = c("fixedrand",
                           "outcome",
                           "mainorinter",
                           "sg",
                           "nct_id"),
           sep = "_|\\.", remove = TRUE)  %>% 
  mutate(nct_id = if_else(is.na(nct_id), "none", nct_id),
         sg = if_else(mainorinter == "nointer", "main", sg))   %>% 
  select(-outcome) %>% 
  mutate(modelnum = if_else(str_length(modelnum) ==2, 
                            str_replace(modelnum, "m", "m0"),
                            modelnum),
         datalevel = if_else(nct_id == "noipd", "noipd","aggipd"))

mace_priors <- mace_priors %>% 
  mutate(`C: Covariates main effects and treatment interactions` = if_else(mainorinter == "nointer", "", 
                                                                           `C: Covariates main effects and treatment interactions`),
         mainorinter = if_else(mainorinter == "nointer", "No covariates", "Age, sex and treatment interactions"),
         sg = case_when(
           sg == "age" ~ "IPD, age-subgroup and aggregate",
           sg == "sex" ~ "IPD, sex-subgroup and aggregate",
           TRUE ~ "IPD and aggregate")) %>% 
  rename(covariates = mainorinter,
         fe_or_re = fixedrand,
         data_used = sg, 
         trials_dropped_or_downgraded = nct_id) %>% 
  select(modelnum, everything()) %>% 
  arrange(modelnum, fe_or_re, data_used)
write_csv(mace_priors, "Outputs/priors_meta_analysis_mace.csv", na = "")

## process data
mace_age_sex <- mace %>% 
  filter(mainorinter == "agesex",
         # datalevel == "aggipd",
         # sg == "main",
         params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:")
mace_age_sex <- mace_age_sex %>% 
  inner_join(who_atc %>% 
               select(trtclass = `ATC code`,
                      cls_raw = `ATC level name`,
                      cls = `ATC level name`) %>% 
               distinct()) %>% 
  mutate(cls = paste0(trtclass, ":", cls))
mace_age_sex <- mace_age_sex %>% 
  janitor::clean_names() 
mace_age_sex <- mace_age_sex %>% 
  mutate(across(mean:x97_5_percent, ~ case_when(
      covariate %in% c("age") ~ .x*30,
      TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age"), "age30", covariate))  %>% 
  mutate(covariate = factor(covariate,
                            levels = c("age30", "male"),
                            labels = c("Age per 30 years",
                                       "Male sex")))  %>% 
  inner_join(nbyclass) 

## create nested list for embedding results
saveRDS(mace_age_sex, "Scratch_data/pre_mace_results_ms.Rds")

## interaction plots mace ----
# "sens", "sens2"
# "None (one SGLT2 IPD trial excluded)",
# "None (one DPP-4 IPD trial excluded)"
a <- c(0, log(0.6), log(0.7), log(0.8))
b <- c(0, log(1.25), log(1.6), log(2.1))
ab <- union(a,b) %>% sort()
ab_lbl <- exp(ab) %>% round(2) %>% formatC(format = "f", digits = 2)
ab_lbl[ab_lbl == "1.00"] <- "1"
mace_age_sex <- mace_age_sex %>% 
  mutate(top = paste0("Triple therapy",
                      "\n",
                      nwork_trials,
                      " trials. ",
                      formatC(nwork_n, format = "d", big.mark = ","),
                      " participants")) %>% 
  mutate(cls_raw_n = paste0(cls_raw,  "\n(", cls_trials, " trials)"))
intermaceplot <- ggplot(mace_age_sex %>% 
                          filter(outcome == "mace", sg == "main",
                                 nct_id == "none",
                                 trtclass %in% c("A10BH", "A10BJ", "A10BK"),
                                 datalevel == "aggipd"), 
                        aes(x = cls_raw_n, y = mean, ymin = x2_5_percent, ymax = x97_5_percent,
                            colour = factor(fixedrand, levels = c("random", "fixed"), labels = c("Random", "Fixed")))) +
  geom_point(position = position_dodge(0.5),  shape = "square", mapping = aes(size = nwork_p)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ top) + 
  scale_x_discrete("", limits = rev) +
  scale_color_discrete("Parameterisation of trial estimates within drug classes", limits = rev) +
  coord_flip(ylim = c(-0.75, 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_minimal5() +
  scale_y_continuous("Hazard ratio", 
                     breaks = ab,
                     labels = ab_lbl) +
  scale_size_area(guide = "none", limits = c(0, 1))
intermaceplot

## sensitivity analyses for MACE ----
intermaceplotnoipd <- intermaceplot %+% (mace_age_sex %>% 
  filter(outcome == "mace", 
         nct_id == "noipd",
         trtclass %in% c("A10BH", "A10BJ", "A10BK"),
         datalevel == "noipd")) +
  coord_flip(ylim = c(-5, 5))  +
  scale_y_continuous("Hazard ratio", 
                     breaks = seq(-5, 5, 2), 
                     labels = round(exp(seq(-5, 5, 2)),2))

mace_age_sex <- mace_age_sex %>% 
  filter(!datalevel == "noipd")

nct_id_lkp <- nct_id_lkp %>% 
  mutate(trtcls5 = case_when(
    trtcls5 == "A10BH" ~ "DPP-4",
    trtcls5 == "A10BJ" ~ "GLP-1",
    trtcls5 == "A10BK" ~ "SGLT-2",
  ))

macesens <- mace_age_sex %>% 
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
         covariate == "Age per 30 years")
macesenssex <- macesens %>% 
  filter(sg %in% c("main", "sex"),
         covariate == "Male sex")
## examine impact of using sex-subgroup on age analysis
macesensage_contra <- macesens %>% 
  filter(sg %in% c("main", "sex"),
         covariate == "Age per 30 years")
## examine impact of using age-subgroup on sex analysis
macesenssex_contra  <- macesens %>% 
  filter(sg %in% c("main", "age"),
         covariate == "Male sex")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

macesensageplt <- ggplot(macesensage,
                         aes(x = xvar,
                             y = mean, 
                             ymin = x2_5_percent, 
                             ymax = x97_5_percent, 
                             colour = sg_lbl)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_wrap(~trtclass, scales = "free_y", ncol = 3) + 
  scale_x_discrete("", limits = rev) +
  scale_color_manual("", values = cbbPalette[1:4]) +
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

macesensage_forplot <- macesensage %>% 
  filter( (nct_id %in% ipdwsggot & trtclass == trtcls5) |
            nct_id == "none") %>% 
  mutate(sg_lbl = case_when(
    sg_lbl == "Yes, ipdreplace" ~ "SG data rather than IPD for dropped/downgraded trial
SG data where available for other trials",
    sg_lbl == "Yes, noipdreplace" ~ "No data for dropped/downgraded trial
SG data where available for other trials",
    sg_lbl == "No" ~ "No data for dropped/downgraded trial
No SG data for any trial"),
    xvar = if_else(nct_id == "none", 
                   xvar,
                   paste0(xvar, "\ndropped/downgraded")))
macesenssex_forplot  <- macesenssex %>% 
  filter( (nct_id %in% ipdwsggot & trtclass == trtcls5) |
            nct_id == "none") %>% 
  mutate(sg_lbl = case_when(
    sg_lbl == "Yes, ipdreplace" ~ "SG data rather than IPD for dropped/downgraded trial
SG data where available for other trials",
    sg_lbl == "Yes, noipdreplace" ~ "No data for dropped/downgraded trial
SG data where available for other trials",
    sg_lbl == "No" ~ "No data for dropped/downgraded trial
No SG data for any trial"),
    xvar = if_else(nct_id == "none", 
                   xvar,
                   paste0(xvar, "\ndropped/downgraded")))

macesensageplt_paper <- macesensageplt %+% 
  macesensage_forplot +
  aes(colour = sg_lbl, linetype = NULL, shape = NULL) +
  theme_minimal4()  +
  theme(legend.position="bottom")
macesenssexplt_paper <- macesenssexplt %+%  
  macesenssex_forplot +
  aes(colour = sg_lbl, linetype = NULL, shape = NULL) +
  theme_minimal4() +
  theme(legend.position="bottom")
pdf("Outputs/sens_onetrial.pdf", width = 15, height = 8)
macesensageplt_paper + ggtitle("Age-treatment interaction with/without age subgroup data")
macesenssexplt_paper + ggtitle("Sex-treatment interaction with/without sex subgroup data")
dev.off()


## main effects MACE ----
main_mace <- mace %>% 
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
write_csv(main_mace, "Outputs/mace_no_inter_ma_res.csv")
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

regplots <- list(intermaceplot = intermaceplot,
                 intermaceplotnoipd = intermaceplotnoipd,
                 macesensageplt = macesensageplt,
                 macesenssexplt = macesenssexplt,
                 macesensageplt_contra = macesensageplt_contra,
                 macesenssexplt_contra = macesenssexplt_contra,
                 mainmacecplot = mainmacecplot,
                 macesenssexplt_paper = macesenssexplt_paper,
                 macesensageplt_paper = macesensageplt_paper
                 )
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
