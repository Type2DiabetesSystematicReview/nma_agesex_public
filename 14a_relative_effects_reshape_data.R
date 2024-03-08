library(tidyverse)
## Code designed to link estimates for classes to drugs within those classes

ds <- readRDS("Scratch_data/tx_samples.Rds")
betas <- readRDS("Scratch_data/cov_nter_samples.Rds")
lkp <- read_csv("Scratch_data/modelname_content_lkp.csv")


## select only models of interest - 4 models. HbA1c for mono, dual and triple ----
lkp2 <- lkp %>% 
  filter(mainorinter == "agesex",
         datalevel == "aggipd",
         fixedrand == "random",
         sg == "main")
ds <- ds[lkp2$tosep]
betas <- betas[lkp2$tosep]

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

## pull in arm names to link ----
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

## Pull whole model level coefficients unrelated to drug or drug class (value_1, age, and sex) ----
## and drop from other betas
wholemodel <- c("beta[male]", "beta[value_1]", "beta[age10]", "beta[age15]")
beta_agesexvalue <- betas %>% 
  filter(param %in% wholemodel)
betas <-  betas %>% 
  filter(!param %in% wholemodel)
beta_agesexvalue <- beta_agesexvalue %>% 
  separate(param, into = c("beta", "main_term"), sep = "\\[") %>% 
  mutate(main_term = str_sub(main_term, 1, -2)) %>% 
  spread(main_term, value) %>% 
  select(-beta)

## separate betas and ds params into covariate treatment etc and link ds to betas using a distinct table ----
## Then map bacross to raw data
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
## reviewed where trtcls5 missing, all ATC codes for hba1c. paste across. Deal with arm labelling later
ds_dist <- ds_dist %>% 
  mutate(trtcls5 = if_else(is.na(trtcls5), str_sub(trt, 1, 5), trtcls5)) %>% 
  select(param, d, trt, trtcls5)
## rearrange betas so wide format, this is mapping between class and parameterss
betas_dist_wide <- betas_dist %>% 
  rename(trtcls5 = trtclass) %>% 
  spread(covariate, param)
bth_dist <- ds_dist %>% 
  left_join(betas_dist_wide) 
rm(ds_dist)

## separate ds and betas into mace and hba1c and join  onto drug effects ----
ds_mace <- ds %>% 
  semi_join(lkp %>% 
              filter(outcome == "mace") %>%
              select(tosep))
ds_mace <- ds_mace %>% 
  left_join(bth_dist %>% select(param, trtcls5, age15, male))
beta_mace <- betas %>% 
  semi_join(lkp %>% 
              filter(outcome == "mace") %>%
              select(tosep)) 
ds_mace <- ds_mace %>% 
  rename(d_value = value) %>% 
  left_join(beta_mace %>% rename(age15 = param, age15_inter = value)) %>% 
  left_join(beta_mace %>% rename(male = param, male_inter = value)) %>% 
  left_join(beta_agesexvalue %>% select(tosep, smpls, age15_main = age15, male_main = male)) %>% 
  select(tosep, trtcls5, d = param, age15, male, smpls, d_value, age15_inter, male_inter, age15_main, male_main)
rm(beta_mace)

## drop mace models and columns from overall effects
beta_agesexvalue <- beta_agesexvalue %>% 
  filter(tosep != "random_mace_agesex_main") %>% 
  select(-age15)
ds_hba1c <- ds %>% 
  semi_join(lkp %>% 
              filter(outcome == "hba1c") %>%
              select(tosep)) 
ds_hba1c <- ds_hba1c %>% 
  left_join(bth_dist %>% select(param, trtcls5, age10, male))
beta_hba1c <- betas %>% 
  semi_join(lkp %>% 
              filter(outcome == "hba1c") %>%
              select(tosep))
ds_hba1c <- ds_hba1c %>% 
  rename(d_value = value) %>% 
  left_join(beta_hba1c %>% rename(age10 = param, age10_inter = value)) %>% 
  left_join(beta_hba1c %>% rename(male = param, male_inter = value)) %>% 
  left_join(beta_agesexvalue %>% select(tosep, smpls, age10_main = age10, male_main = male)) %>% 
  select(tosep, trtcls5, d = param, age10, male, smpls, d_value, age10_inter, male_inter, age10_main, male_main)
rm(betas, ds)

## rescale age so same for both ----
ds_mace <- ds_mace %>% 
  rename(modelname = tosep) %>% 
  mutate(age10_inter = 10 * age15_inter/15,
         age10_main  = 10 * age15_main/15,
         age10 = age15) %>% 
  select(-age15, -age15_inter, -age15_main)
ds_hba1c <- ds_hba1c %>% 
  rename(modelname = tosep)
ds_tot <- bind_rows(hba1c = ds_hba1c,
                    mace = ds_mace,
                    .id = "outcome_rel") 
rm(ds_mace, ds_hba1c, beta_hba1c, beta_agesexvalue)
ds_tot <- inner_join(lkp %>% rename(modelname = tosep),
                      ds_tot)
## save objects ----
saveRDS(ds_tot, "Scratch_data/random_effects_smpls.Rds")

