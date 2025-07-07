# 28_confirm_new_data
library(tidyverse)
library(metafor)

tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")

tot <- tot %>% 
  select(drug_regime_smpl, reg) %>% 
  unnest(reg)
plac <- tot %>% 
  filter(reference_arm ==1, trtcls5 == "place") %>% 
  select(nct_id) %>% 
  pull()

agef4 <- tot %>% 
  filter(nct_id %in% plac,
         models == "f4",
         age10 == 10, 
         is.na(sex)) %>% 
  select(nct_id, estimate, std.error, trtcls5, arm_lvl, term)
a10bh_age <- agef4 %>% 
  filter(trtcls5 == "A10BH")
a10bj_age <- agef4 %>% 
  filter(trtcls5 == "A10BJ")
a10bk_age <- agef4 %>% 
  filter(trtcls5 == "A10BK")
# mod1_bh <- rma(data = a10bh_age, yi = estimate, sei = std.error)
mod1_bj_age <- rma(data = a10bj_age, yi = estimate, sei = std.error)
mod1_bk_age <- rma(data = a10bk_age, yi = estimate, sei = std.error)

## same result in updated data for lm model
# estimate      se    zval    pval   ci.lb   ci.ub      
# 0.0973  0.0153  6.3503  <.0001  0.0673  0.1274  *** 
## updated
# 0.1011 0.0161 6.28 p < 0.001 0.0696 to 0.1327


sexf4 <- tot %>% 
  filter(nct_id %in% plac,
         models == "f4",
         is.na(age10), 
         sex) %>% 
  select(nct_id, estimate, std.error, trtcls5, arm_lvl, term)
a10bh_sex <- sexf4 %>% 
  filter(trtcls5 == "A10BH")
a10bj_sex <- sexf4 %>% 
  filter(trtcls5 == "A10BJ")
a10bk_sex <- sexf4 %>% 
  filter(trtcls5 == "A10BK")
mod1_bh_sex <- rma(data = a10bh_sex, yi = estimate, sei = std.error)
mod1_bj_sex <- rma(data = a10bj_sex, yi = estimate, sei = std.error)
mod1_bk_sex <- rma(data = a10bk_sex, yi = estimate, sei = std.error)
