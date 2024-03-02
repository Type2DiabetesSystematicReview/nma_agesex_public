# summarise mace separately

mace_agg <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_agg
pseudo_mace <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
mace_arms <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_arms
pseudo_mace <- pseudo_mace %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, arm_lvl, trtcls5))
maceinhba1c_ipd <- intersect(pseudo_mace$nct_id, pseudo_hba1c$nct_id)

pseudo_hba1c <- pseudo_hba1c %>% 
  anti_join(pseudo_mace %>% select(nct_id))

mace = pseudo_mace %>% select(nct_id, sex, age, arm_lvl, trtcls5, arm_lvl)

mace_agg <- mace_agg %>% 
  select(nct_id, age_m, age_sd = age_s, n = participants,
         trtcls5,
         arm_lvl = arm_lvl)
maceinhba1c_agg <- intersect(mace_agg$nct_id, hba1c_agg$nct_id)
maceinhba1c <- union(maceinhba1c_agg, maceinhba1c_ipd)
hba1c_agg <- hba1c_agg %>% 
  anti_join(mace_agg %>% select(nct_id))

mace = mace_agg

## note that MACE novel antidiabetic trials are all against non novel antidiabetic comparators
## spo easier to code (all mutually exclusive)
mace <- ages %>% 
  filter(outcome == "mace")
mace <- mace %>% 
  group_by(nct_id) %>% 
  mutate(trtcls5 = case_when(
    any(trtcls5 %in% novel[1]) ~ novel[1],
    any(trtcls5 %in% novel[2]) ~ novel[2],
    any(trtcls5 %in% novel[3]) ~ novel[3],
  )) %>% 
  ungroup()
mace_smry1 <- mace %>% 
  group_by(trtcls5, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup()
mace_smry2 <- mace %>% 
  group_by(data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() %>% 
  mutate(trtcls5 = "Any")
mace_smry <- bind_rows(mace_smry1,
                       mace_smry2) %>% 
  rename(cls = trtcls5) %>% 
  left_join(ages_smry %>% distinct(cls, trl_lbl))
mace_smry <- mace_smry[, names(ages_smry)]
write_csv(mace_smry, "Outputs/age_summary_mace.csv")
