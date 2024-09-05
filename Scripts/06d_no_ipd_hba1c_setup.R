#26_noIPD_model

library(tidyverse)
ColnamesPipe <- function(x, vct){
  colnames(x) <- vct
  x
}

source("Scripts/common_functions/Scripts/combine_sd.R")
## reformat age and sex data ----
age <- readRDS("Scratch_data/ipd_consolidated.Rds")$age_distribution_baseline_continuous
age$splt <-  str_split(age$`Age (years) at quantiles (%)`, pattern = "\\,")
age$splt <- map(age$splt, ~ str_split_fixed(.x, "\\=", n = 2) %>% ColnamesPipe(c("x", "y")) %>% 
                  as_tibble() %>% 
                  mutate(across(everything(), ~ str_trim(.x) %>% as.integer()),
                         x = x/100))
# Add zero to two trials where not age zero. Set as same age as next value as same as earlier age NCT01381900 and NCT02065791 
age$nozero <- map_lgl(age$splt, ~ ! any(.x$x ==0))
age$splt <- map2(age$nozero,
                 age$splt, function(condition, mydf){
                   if(is.na(condition) | !condition) mydf else {
                     a <- mydf %>% 
                       slice(1) %>% 
                       mutate(x = 0)
                     bind_rows(a, mydf)
                   }
                 })
age$min <- map_dbl(age$splt, ~ .x$y[1])
age$max <- map_dbl(age$splt, ~ .x$y[length((.x$y))])
age <- age %>% 
  select(-splt, -nozero)
age <- age %>% 
  group_by(nct_id) %>% 
  mutate(min  = min(min),
         max = max(max)) %>% 
  ungroup()
age <- age %>% 
  group_by(nct_id, nct_id2, arm_id_unq) %>% 
  mutate(age_m = weighted.mean(mean, participants),
         age_sd = CombSdVectorised(participants, mean, sd)) %>% 
  ungroup()
age <- age %>%
  group_by(nct_id, nct_id2, arm_id_unq) %>% 
  mutate(sex_prcnt = participants/sum(participants),
         n = sum(participants)) %>% 
  ungroup()
age <- age %>% 
  filter(sex == "male") %>% 
  select(nct_id, nct_id2, arm_id_unq, age_sd = age_sd, age_m = age_m, male = sex_prcnt, n = n) %>% 
  distinct()

## add values at baseline and end ----
base_end <- readRDS("Scratch_data/ipd_consolidated.Rds")$hba1c_base_change_overall
base_end <- base_end %>% 
  mutate(se = sd/n^0.5) %>% 
  select(-result) %>% 
  select(nct_id, nct_id2, arm_id_unq, value_1 = base, result = end, se) %>% 
  distinct()

## add age and sex and value2 to treatment information ----
## note took value1_sd from regression (b1) but not other values
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
tot$aggedipd <- map(tot$ipd, ~ .x %>% 
                      group_by(nct_id, nct_id2, arm_lvl, trtcls5, trtcls4, arm_id_unq) %>% 
                      summarise(value_1_sd = sd(value_1_rsd)) %>% 
                      ungroup())
tot$aggedipd <- map(tot$aggedipd, ~ .x %>% 
                      inner_join(age) %>% 
                      inner_join(base_end) %>% 
                    select(-nct_id) %>% 
                    rename(nct_id = nct_id2))
tot <- tot %>% 
  select(drug_regime_smpl, agg, aggedipd) %>% 
  gather("aggtype", "agg", -drug_regime_smpl)
tot <- tot %>% 
  unnest(agg) %>% 
  select(-aggtype) %>% 
  nest(.by = drug_regime_smpl)
saveRDS(tot, "Scratch_data/agg_agged_ipd_hba1c.Rds")


  