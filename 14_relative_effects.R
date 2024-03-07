library(tidyverse)

ds <- readRDS("Scratch_data/tx_samples.Rds")
betas <- readRDS("Scratch_data/cov_nter_samples.Rds")
lkp <- read_csv("Scratch_data/modelname_content_lkp.csv")
mace_arms <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_arms
mace_arms <- mace_arms %>% 
  distinct(arm_lvl, trtcls5, atc, drug_name, drug_dose)
hba1c_arms <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
hba1c_arms <- hba1c_arms %>% 
  gather("data_type", "data", -drug_regime_smpl)
hba1c_arms$data <- map(hba1c_arms$data, ~ .x %>% 
                         select(matches("drug_name|drug_dose|arm_lvl|trtcls5")) %>% 
                  distinct())
hba1c_arms <- hba1c_arms %>% 
  select(data) %>% 
  unnest(data) %>% 
  distinct()
hba1c_arms_nomis <- hba1c_arms %>% 
  filter(!is.na(drug_name))
hba1c_arms_nomis %>% 
  anti_join(hba1c_arms %>% select(arm_lvl))
rm(hba1c_arms)
hba1c_arms <- hba1c_arms_nomis
hba1c_arms <- hba1c_arms %>% 
  separate(arm_lvl, into = c("atc"), sep = "_", extra = "drop", remove = FALSE) %>% 
  select(-drug_dose_orig, -drug_dose_simplify)
rm(hba1c_arms_nomis)
arms <- bind_rows(mace = mace_arms,
                  hba1c = hba1c_arms, .id = "outcome")
arms <- arms %>% 
  select(-outcome) %>% 
  distinct()
rm(mace_arms, hba1c_arms)


ds <- map(ds, function(a) {
  a %>% 
    mutate(smpls = 1:nrow(a)) %>% 
    gather("param", "value", -smpls)
})
ds <- bind_rows(ds, .id = "tosep") %>% 
  as_tibble()
betas <- map(betas, function(a) {
  a %>% 
    mutate(smpls = 1:nrow(a)) %>% 
    gather("param", "value", -smpls)
})
betas <- bind_rows(betas, .id = "tosep") %>% 
  as_tibble()
## drop models without betas
betas <- betas %>% 
  filter(!tosep %in% c("fixed_mace_nointer",
                       "random_mace_nointer"))

## reduce number of iterations
warning("Need to drop this code, just doing while developing")
betas <- betas %>% 
  group_by(tosep, param) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()
ds <- ds %>% 
  group_by(tosep, param) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()
## separate betas into covariate treatment etc
betas_dist <- betas %>% 
  distinct(param)
betas_dist <- betas_dist %>% 
  mutate(param_new = param %>% 
          str_remove("^beta\\[") %>% 
          str_remove("\\]") %>%
          str_remove("\\.trtclass")) %>% 
  separate(param_new, into = c("covariate", "trtclass"), sep = "\\:")

## repeat for ds
ds_dist <- ds %>% 
  distinct(param)
ds_dist <- ds_dist %>% 
  separate(param, into = c("d", "trt"), sep = "\\[", remove = FALSE) %>% 
  mutate(trt = str_sub(trt, 1, -2))
## chekd dups in arm, just where simplified semaglutide dose etc
ds_dist <- ds_dist %>% 
  left_join(arms %>% 
              distinct(arm_lvl, trtcls5, atc, .keep_all = TRUE) %>% 
              mutate(trt = arm_lvl))
## reviewed where trtcls5 missing, all ATC codes
ds_dist <- ds_dist %>% 
  mutate(trtcls5 = if_else(is.na(trtcls5), str_sub(trt, 1, 5), trtcls5))

bth_dist <- ds_dist %>% 
  left_join(betas_dist %>% rename(trtcls5 = trtclass) %>% 
              spread(covariate, param)) %>% 
  select(param_ds = param, trt, trtcls5, age10, age15, male, value_1) %>% 
  distinct()

ds <- ds %>% 
  inner_join(bth_dist %>% rename(param = param_ds))

## limit to mace trials
ds2 <- ds %>% 
  semi_join(lkp %>% 
              filter(outcome == "mace",
                     mainorinter== "agesex") %>%
              select(tosep))

## Checked and have correctly joined
ds3 <- ds2 %>% 
  # left_join(betas %>% rename(age10 = param, age10_value = value)) %>% 
  # left_join(betas %>% rename(value_1 = param,  value_1_value = value)) %>% 
  left_join(betas %>% rename(age15 = param, age15_value = value)) %>% 
  left_join(betas %>% rename(male = param,  male_value = value)) %>% 
  select(-value_1, -age10)

ds3 <- ds3 %>% 
  mutate(loghr = value + age15_value*(60/15) + 0*male_value)



