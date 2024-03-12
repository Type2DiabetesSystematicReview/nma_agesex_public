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
  count(trl_lbl, nm, data_lvl) %>%
  mutate(n = if_else(trl_lbl == nm, "-", as.character(n))) %>% 
  pivot_wider(names_from = c(trl_lbl, data_lvl), values_from = n, values_fill = "0")   %>% 
  arrange(desc(as.integer(`Total trials_agg`)))
comparisons2 <- hba1c_cls %>% 
  distinct(nct_id, trl_lbl, data_lvl) %>% 
  count(trl_lbl, data_lvl) %>%
  mutate(n = as.character(n)) %>% 
  pivot_wider(names_from = c(trl_lbl, data_lvl), values_from = n, values_fill = "0")  %>% 
  mutate(nm = "Total")
comparisons <- bind_rows(comparisons2[, names(comparisons)],
                         comparisons)
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
                     male,
                     ages ,
                     .id = "orig_tbl") %>% 
  select(var, trl_lbl, data_lvl, lvls, measure, res)
write_csv(tbl_lng, "Outputs/manuscript_table1a_machine_readable.csv", na = "")

## Note participant count from ages summary and male summary are, as expected, identical
tbl_lng %>% filter(var == "age" & measure == "n" | var == "participants")

## produce "Nice" format for report/paper
tbl_lng2 <- tbl_lng %>% 
  filter(!measure == "n") %>% 
  mutate(res = case_when(
    measure %in% "count" ~ round(res) %>% formatC(digits = 0, format = "f"),
    measure %in% c("m", "s", "q05", "q95", "percent") ~ round(res, 1)  %>% formatC(digits = 1, format = "f")
  ))
age <- tbl_lng2 %>% 
  filter(var == "age") %>%
  spread(measure, res) %>% 
  mutate(res = paste0(m, " (",s, ")", " [",q05, "-", q95, "]")) %>% 
  select(var, trl_lbl, data_lvl, res)
male <- tbl_lng2 %>% 
  filter(var == "male") %>%
  spread(measure, res) %>% 
  mutate(res = paste0(count, " (", percent, "%)")) %>% 
  select(var, trl_lbl, data_lvl, res)
tbl_lng3 <- bind_rows(tbl_lng2 %>% 
                        filter(var %in% c("trials", "arms", "participants")),
                      age,
                      male) 
tbl_wide <- tbl_lng3 %>% 
  pivot_wider(names_from = c(trl_lbl, data_lvl), values_from = res) %>% 
  arrange(match(var, c("trials",  "arms",
                       "participants","male", "age")))

tbl_wide2 <- tbl_wide %>% 
  mutate(collbl = if_else(is.na(lvls),
                          str_to_sentence(var),
                          paste0(lvls, " ", var)),
         collbl = if_else(!duplicated(collbl), collbl, "")) %>% 
  select(-var, -lvls, -measure) %>% 
  select(collbl, everything())

write_csv(tbl_wide2, "Outputs/manuscript_table1a.csv", na = "")
write_csv(comparisons, "Outputs/manuscript_hba1c_comparisons.csv", na = "")

