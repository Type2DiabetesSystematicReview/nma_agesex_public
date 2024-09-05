library(tidyverse)
## Code designed to link estimates for classes to drugs within those classes
## for simplicity run separately for MACE and hba1c

ds <- readRDS("Scratch_data/tx_samples.Rds")
betas <- readRDS("Scratch_data/cov_nter_samples.Rds")

## select only models of interest - 4 models. HbA1c for mono, dual and triple ----
lkp <- read_csv("tosep,outcome,mainorinter,modelnum,datalevel,fixedrand,network,sg
m06_aggipd_random_mono_f4_ge12,hba1c,agesex,m06,aggipd,random,mono,main
m02_aggipd_random_dual_f4_ge12,hba1c,agesex,m02,aggipd,random,dual,main
m10_aggipd_random_triple_f4_ge12,hba1c,agesex,m10,aggipd,random,triple,main
random_mace_agesex_main,mace,agesex,m6,aggipd,random,triple,main
random_mace_agesex_sex,mace,agesex,m7,aggipd,random,triple,sex")

ds <- ds[lkp$tosep]
betas <- betas[lkp$tosep]

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

## split into hba1c and mace as differnt names ----
## ## three trials
lkp_hba1c <- lkp %>% 
  filter(outcome == "hba1c")
## Single trial
lkp_mace <- lkp %>% 
  filter(outcome == "mace", sg == "sex")
betas_mace <- betas %>% 
  filter(tosep %in% lkp_mace$tosep)
betas_hba1c <- betas %>% 
  filter(tosep %in% lkp_hba1c$tosep)
rm(betas)
ds_mace <- ds %>% 
  filter(tosep %in% lkp_mace$tosep)
ds_hba1c <- ds %>% 
  filter(tosep %in% lkp_hba1c$tosep)
rm(ds)

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

## Pull whole model level coefficients unrelated to drug or drug class (value_1, age, and sex) ----
## and drop from other betas
wholemodel <- c("beta[male]", "beta[value_1]", "beta[age]", "beta[age10]")

beta_agesexvalue_mace <- betas_mace %>% 
  filter(param %in% wholemodel)
betas_mace <-  betas_mace %>% 
  filter(!param %in% wholemodel)
beta_agesexvalue_hba1c <- betas_hba1c %>% 
  filter(param %in% wholemodel)
betas_hba1c <-  betas_hba1c %>% 
  filter(!param %in% wholemodel)

beta_agesexvalue_mace <- beta_agesexvalue_mace %>% 
  separate(param, into = c("beta", "main_term"), sep = "\\[") %>% 
  mutate(main_term = str_sub(main_term, 1, -2)) %>% 
  spread(main_term, value) %>% 
  select(-beta)
beta_agesexvalue_hba1c <- beta_agesexvalue_hba1c %>% 
  separate(param, into = c("beta", "main_term"), sep = "\\[") %>% 
  mutate(main_term = str_sub(main_term, 1, -2)) %>% 
  spread(main_term, value) %>% 
  select(-beta)

## separate betas and ds params into covariate treatment etc and link ds to betas using a distinct table ----
## Then map bacross to raw data
betas_dist_hba1c <- betas_hba1c %>% 
  distinct(param)
betas_dist_mace <- betas_mace %>% 
  distinct(param)
betas_dist_hba1c <- betas_dist_hba1c %>% 
  mutate(param_new = param %>% 
          str_remove("^beta\\[") %>% 
          str_remove("\\]") %>%
          str_remove("\\.trtclass")) %>% 
  separate(param_new, into = c("covariate", "trtclass"), sep = "\\:")
betas_dist_mace <- betas_dist_mace %>% 
  mutate(param_new = param %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(param_new, into = c("covariate", "trtclass"), sep = "\\:")

## repeat for ds
ds_dist_hba1c <- ds_hba1c %>% 
  distinct(param)
ds_dist_mace <- ds_mace %>% 
  distinct(param)
ds_dist_hba1c <- ds_dist_hba1c %>% 
  separate(param, into = c("d", "trt"), sep = "\\[", remove = FALSE) %>% 
  mutate(trt = str_sub(trt, 1, -2))
ds_dist_mace <- ds_dist_mace %>% 
  separate(param, into = c("d", "trt"), sep = "\\[", remove = FALSE) %>% 
  mutate(trt = str_sub(trt, 1, -2))
## chekd dups in arm, just where simplified semaglutide dose etc
ds_dist_hba1c <- ds_dist_hba1c %>% 
  left_join(mace_arms %>% 
              distinct(arm_lvl, trtcls5, atc, .keep_all = TRUE) %>% 
              mutate(trt = arm_lvl))
ds_dist_mace <- ds_dist_mace %>% 
  left_join(mace_arms %>% 
              distinct(arm_lvl, trtcls5, atc, .keep_all = TRUE) %>% 
              mutate(trt = arm_lvl))

## reviewed where trtcls5 missing, all ATC codes for hba1c. paste across. Deal with arm labelling later
ds_dist_hba1c <- ds_dist_hba1c %>% 
  mutate(trtcls5 = if_else(is.na(trtcls5), str_sub(trt, 1, 5), trtcls5)) %>% 
  select(param, d, trt, trtcls5)
ds_dist_mace <- ds_dist_mace %>% 
  select(param, d, trt, trtcls5)

## rearrange betas so wide format, this is mapping between class and parameterss
betas_dist_wide_hba1c <- betas_dist_hba1c %>% 
  rename(trtcls5 = trtclass) %>% 
  spread(covariate, param)
betas_dist_wide_mace <- betas_dist_mace %>% 
  rename(trtcls5 = trtclass) %>% 
  spread(covariate, param)
bth_dist_hba1c <- ds_dist_hba1c %>% 
  left_join(betas_dist_wide_hba1c) 
bth_dist_mace <- ds_dist_mace %>% 
  left_join(betas_dist_wide_mace) 
rm(ds_dist_hba1c, ds_dist_mace)

## separate ds and betas into mace and hba1c and join  onto drug effects ----
ds_mace <- ds_mace %>% 
  left_join(bth_dist_mace %>% select(param, trtcls5, age, male))
ds_hba1c <- ds_hba1c %>% 
  left_join(bth_dist_hba1c %>% select(param, trtcls5, age10, male))

ds_mace <- ds_mace %>% 
  rename(d_value = value) %>% 
  left_join(betas_mace %>% rename(age = param, age_inter = value)) %>% 
  left_join(betas_mace %>% rename(male = param, male_inter = value)) %>% 
  left_join(beta_agesexvalue_mace %>% select(tosep, smpls, age_main = age, male_main = male)) %>% 
  select(tosep, trtcls5, d = param, age, male, smpls, d_value, age_inter, male_inter, age_main, male_main)
rm(betas_mace)

ds_hba1c <- ds_hba1c %>% 
  rename(d_value = value) %>% 
  left_join(betas_hba1c %>% rename(age10 = param, age_inter = value)) %>% 
  left_join(betas_hba1c %>% rename(male = param, male_inter = value)) %>% 
  left_join(beta_agesexvalue_hba1c %>% select(tosep, smpls, age_main = age10, male_main = male, bline_hba1c = value_1)) %>% 
  select(tosep, trtcls5, d = param, age = age10, male, smpls, d_value, age_inter, male_inter, age_main, male_main, bline_hba1c)
rm(betas_hba1c)

## drop a10bx from triple therapy as does not fit
ds_hba1c <- ds_hba1c %>% 
  filter(! (tosep == "hba1c_agesex_m6_aggipd_random_triple" & trtcls5 == "A10BX"))

## save objects ----
ds_mace <- ds_mace %>% 
  rename(modelname = tosep) 
ds_hba1c <- ds_hba1c %>% 
  rename(modelname = tosep)
ds_tot <- bind_rows(hba1c = ds_hba1c,
                    mace = ds_mace,
                    .id = "outcome_rel") 
rm(ds_mace, ds_hba1c)
ds_tot <- inner_join(lkp %>% rename(modelname = tosep),
                      ds_tot)
## save objects ----
saveRDS(ds_tot, "Scratch_data/random_effects_smpls.Rds")

