# 13_tables
library(tidyverse)
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
source("Scripts/04b_arm_label_reverse.R")
source("../common_functions/Scripts/misc.R")
## Hba1c tables
agg <- tot %>% 
  select(drug_regime_smpl, agg) %>% 
  unnest(agg)

## Count trials with each novel class (no double counting for simplicity and consistency of reporting)
agg <- agg %>% 
  group_by(nct_id) %>% 
  mutate(cls_slct = case_when(
    any(trtcls5 == "A10BK") ~ "A10BK",
    any(trtcls5 == "A10BH") ~ "A10BH",
    any(trtcls5 == "A10BJ") ~ "A10BJ",
  )) %>% 
  ungroup()
## drop treatment arm (other than in cls slct)
agg <- agg %>% 
  filter(!trtcls5 == cls_slct)
## count trials per arm
trials <- agg %>% 
  distinct(cls_slct, nct_id) %>% 
  count(cls_slct) 

## Count number of comparison arms
cmpr_arms <- agg %>% 
  group_by(cls_slct, nct_id) %>% 
  summarise(cmpr_arms = length(nct_id)) %>% 
  ungroup() %>% 
  count(cls_slct,
        cmpr_arms) %>% 
  spread(cls_slct, n, fill = 0L)

## organise comparison arms
agg %>% 
  count(trtcls5)
comparisons <- agg %>% 
  group_by(cls_slct, trtcls5) %>% 
  summarise(n = sum(!duplicated(nct_id))) %>% 
  ungroup() %>% 
  left_join(who_atc %>% select(trtcls5 = `ATC code`, nm = `ATC level name`)) %>% 
  arrange(cls_slct, desc(n)) %>% 
  mutate(nm = if_else(trtcls5 == "place", "Placebo", nm)) %>% 
  mutate(res = paste0(nm, " (", n, ")")) %>% 
  select(cls_slct, res) %>% 
  group_by(cls_slct) %>% 
  summarise(res = PasteAnd(res)) %>% 
  mutate(Comparisons = "") %>% 
  spread(cls_slct, res)

NewBindRows <- function(x) {
  map(x, ~ {
      .x[1,1] <- names(.x)[1]
      names(.x)[1] <- "col1"
  }) %>% 
    bind_rows()
}

