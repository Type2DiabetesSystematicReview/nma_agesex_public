library(tidyverse)
source("../common_functions/Scripts/combine_sd.R")
tot <- readRDS("data_for_mars.Rds")
agg <- tot %>% 
  select(drug_regime_smpl, agg)
ipd <- tot %>% 
  select(drug_regime_smpl, ipd) %>% 
  unnest(ipd) %>% 
  distinct(drug_regime_smpl, nct_id, nct_id2, arm_id_unq, drug_code, trtcls5, trtcls4)
ipd_age_sex <- readRDS("Scratch_data/ipd_age_sex.Rds")
ipd_sex <- ipd_age_sex %>%
  select(nct_id:participants) %>% 
  spread(sex, participants, fill = 0L)
ipd_age <- ipd_age_sex %>% 
  group_by(nct_id, arm_id_unq) %>% 
  summarise(age_m = weighted.mean(mean, participants),
            age_s = CombSdVectorised(participants, mean, sd)) %>% 
  ungroup()
ipd <- ipd %>% 
  inner_join(ipd_sex) %>% 
  inner_join(ipd_age)
rm(ipd_age, ipd_sex, ipd_age_sex)

ipd <- bind_rows(ipd,
                 ipd %>% 
                   mutate(drug_regime_smpl = "total") %>% 
                   distinct())
agg <- tot %>% 
  select(drug_regime_smpl, agg) %>% 
  unnest(agg) %>% 
  distinct()
agg <- bind_rows(agg,
                 agg %>% 
                   mutate(drug_regime_smpl = "total", 
                          drug_regime = "tot", treat_or_ref = "", pooled_arm = 0L) %>% 
                   distinct())
agg_smry <- agg %>% 
  mutate(drug_regime_smpl = factor(drug_regime_smpl,
                                   levels = c("mono", "dual", "triple", "total"))) %>% 
  group_by(drug_regime_smpl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = sum(!duplicated(arm_id_unq)),
            participants = sum((n)),
            male = sum(male),
            age_m_min  = min(age_m),
            age_m_max = max(age_m),
            age_q1  = quantile(age_m, probs = 0.25),
            age_q3 = quantile(age_m, probs = 0.75),
            age_sd = CombSdVectorised(n, age_m, age_sd),
            age_m = weighted.mean(age_m, n)) %>% 
  ungroup() 

ipd_smry <- ipd %>% 
  mutate(drug_regime_smpl = factor(drug_regime_smpl,
                                   levels = c("mono", "dual", "triple", "total")),
         n = male + female) %>% 
  group_by(drug_regime_smpl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            arms = sum(!duplicated(arm_id_unq)),
            participants = sum((n)),
            male = sum(male),
            age_m_min  = min(age_m),
            age_m_max = max(age_m),
            age_q1  = quantile(age_m, probs = 0.25),
            age_q3 = quantile(age_m, probs = 0.75),
            age_sd = CombSdVectorised(n, age_m, age_s),
            age_m = weighted.mean(age_m, n)) %>% 
  ungroup() 
tot_smry <- bind_rows(agg = agg_smry,
                      ipd = ipd_smry, .id = "agg_lvl")
write_csv(tot_smry, "Outputs/demo_included_trials.csv")


## save R print 
tot <- readRDS("Scratch_data/fe_model.Rds")
tot$txt <- map(tot$dm_nets, ~ {
  sink("my_data.txt")
  #write this string to file
  print(.x, n = 1000)
  #close the external connection
  sink() 
  
  a <- read_lines("my_data.txt")
  a <- tibble(raw = a)
  a <- a %>% 
    mutate(grp = cumsum(as.integer(str_detect(raw, "^-"))))
  a <- a %>% 
    group_by(grp) %>% 
    nest()
  a$fl <- map_chr(a$data, ~ .x$raw[1])
  a <- a %>% 
    mutate(fl = str_remove_all(fl, "\\-") %>% str_trim())
  a$data[a$fl == "AgD studies (armbased)"][[1]] <- a$data[a$fl == "AgD studies (armbased)"][[1]] %>% 
    head(21)
  a$data <- map(a$data, pull)
  
  do.call(c, a$data)
})
saveRDS(tot %>% select(drug_regime_smpl, txt), "Scratch_data/net_summary")
