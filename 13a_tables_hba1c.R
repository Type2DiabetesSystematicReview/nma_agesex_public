# 13_tables
library(tidyverse)
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
source("Scripts/04b_arm_label_reverse.R")
source("../common_functions/Scripts/misc.R")

## arrange into same format as summaries
ages <- read_csv("Outputs/age_summary_hba1c.csv")
trials <- ages %>% 
  select(trl_lbl, data_lvl, res = trials) %>% 
  mutate(measure = "count",
         var = "trials")
ages <- ages %>% 
  select(-trials) %>% 
  gather("measure", "res", -trl_lbl, -data_lvl, -cls) %>% 
  mutate(var = "age")

## Hba1c tables
agg <- tot %>% 
  select(drug_regime_smpl, agg) %>% 
  unnest(agg) %>% 
  select(drug_regime_smpl, nct_id, arm_lvl, trtcls5, participants = n, n_arms, male)
ipd <- tot %>% 
  select(drug_regime_smpl, ipd) %>% 
  unnest(ipd) %>% 
  select(drug_regime_smpl, nct_id, arm_lvl, trtcls5, sex) %>% 
  group_by(drug_regime_smpl, nct_id, trtcls5, arm_lvl) %>% 
  summarise(participants = length(sex),
            male = sum(sex)) %>% 
  group_by(drug_regime_smpl, nct_id) %>% 
  mutate(n_arms = sum(!duplicated(arm_lvl))) %>% 
  ungroup()
hba1c <- bind_rows(agg = agg,
                   ipd = ipd, .id = "data_lvl")
rm(agg, ipd)

## Separate each into classes
hba1c_cls <- map(c("A10BK", "A10BH", "A10BJ"), ~ {
              hba1c %>% 
                 group_by(nct_id) %>% 
                 mutate(cls = if_else(any(trtcls5 == .x),
                                      .x,
                                      ""))  %>% 
    ungroup()
              }) 
hba1c_cls <- bind_rows(hba1c_cls) %>% 
  filter(!cls == "")
hba1c_cls <- bind_rows(hba1c_cls, 
                       hba1c %>% 
                         mutate(cls = "Any"))

## Add labels for treatment and comparison arms
hba1c_cls <- hba1c_cls %>% 
  left_join(who_atc %>% select(cls = `ATC code`, trl_lbl = `ATC level name`)) %>% 
  mutate(trl_lbl = if_else(is.na(trl_lbl), "Total trials", trl_lbl)) %>% 
  left_join(who_atc %>% select(trtcls5 = `ATC code`, nm = `ATC level name`)) %>% 
  mutate(nm = if_else(trtcls5 == "place", "Placebo", nm)) 
comparisons <- hba1c_cls %>%
  filter(!cls == trtcls5) %>% 
  group_by(trl_lbl, data_lvl, nm) %>% 
  summarise(n = sum(!duplicated(nct_id))) %>% 
  ungroup() %>% 
  mutate(res = paste0(nm, " (", n, ")")) %>%
  arrange(trl_lbl, desc(n)) %>% 
  group_by(trl_lbl, data_lvl) %>%
  summarise(res_chr = PasteAnd(res)) %>% 
  ungroup() %>% 
  mutate(var = "comparisons",
         measure = "description_n")
male <- hba1c_cls  %>% 
  group_by(trl_lbl, data_lvl) %>% 
  summarise(male = round(sum(male)),
            participants = sum(participants),
            male_prcnt = round(100*male/participants, 1)) %>% 
  ungroup() %>% 
  gather("cmplx", "res", male, participants, male_prcnt) %>% 
  mutate(measure = if_else(cmplx %in% c("male", "participants"), "count", "percent"),
         var = if_else(cmplx %in% c("male", "male_prcnt"), "male", "participants")) %>%
  select(-cmplx)
n_arms <- hba1c_cls %>% 
  mutate(n_arms = if_else(n_arms %in% 2:3, as.character(n_arms), "4 or 5")) %>% 
  group_by(trl_lbl, data_lvl, n_arms) %>% 
  summarise(res = sum(!duplicated(nct_id))) %>% 
  ungroup() %>% 
  rename(lvls =n_arms) %>% 
  mutate(measure = "count",
         var = "arms")

tbl_lng <- bind_rows(trials,
                     n_arms,
                     comparisons,
                     male,
                     ages ,
                     .id = "orig_tbl") %>% 
  select(var, trl_lbl, data_lvl, lvls, measure, res, res_chr)
write_csv(tbl_lng, "Outputs/manuscript_table1a_machine_readable.csv", na = "", col_names = FALSE)

## Note participant count from ages summary and male summary are, as expected, identical
tbl_lng %>% filter(var == "age" & measure == "n" | var == "participants")

tbl_wide <- tbl_lng %>% 
  filter(!measure == "n") %>% 
  mutate(res_chr = case_when(
    measure == "description_n" ~ res_chr,
    measure %in% "count" ~ round(res) %>% formatC(digits = 0, format = "f"),
    measure %in% c("m", "s", "q05", "q95", "percent") ~ round(res, 1)  %>% formatC(digits = 1, format = "f")
  )) %>% 
  select(-res) %>% 
  pivot_wider(names_from = c(trl_lbl, data_lvl), values_from = res_chr) %>% 
  arrange(match(var, c("trials",  "arms","comparisons",
                       "participants","male", "age")))
tbl_lbl1 <- tbl_wide %>% slice(1)
tbl_lbl1[] <- map(names(tbl_wide), ~ str_remove(.x, "_ipd|_agg") )
tbl_lbl2 <- tbl_wide %>%slice(1)
tbl_lbl2[] <- map(names(tbl_wide), ~ str_extract(.x, "ipd|agg") )
tbl_wide_char <- bind_rows(tbl_lbl1,
                      tbl_lbl2,
                      tbl_wide)

write_csv(tbl_wide_char, "Outputs/manuscript_table1a.csv", na = "", col_names = FALSE)

