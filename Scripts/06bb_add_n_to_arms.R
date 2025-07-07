# 06bb_add_n_to_arms
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
getarms <- tot %>% 
  select(reg) 
getarms$reg <- map(getarms$reg, ~ .x %>% 
                     filter(models == "f1") %>% 
                     select(nct_id, nct_id2, arm_id_unq, arm_lvl, reference_arm, 
                            trtcls5, trtcls4, models) %>%
                     distinct())
getarms <- getarms %>% 
  unnest(reg) %>% 
  filter(!is.na(arm_lvl))
## identical arms so join and get new n
totipd <- tot %>% 
  select(ipd) %>% 
  unnest(ipd)
totipd <- totipd %>% 
  count(nct_id, arm_id_unq) %>% 
  rename(narm = n)
getarms <- getarms %>% 
  inner_join(totipd)

## note three arms not included in NMA as complex doses
drpd <- c("NCT00688701", "NCT00712673", "NCT00763451")
## examine which reg have a placebo arm only (as other combination?)
# NCT00688701 plac only

## add n to arms
dgn <- bind_rows(
  dgn1 = read_csv("Data/agesexhba1c_9492/age_sex_model_diags_fixed.csv"),
  dgn2 = read_csv("Data/agesexhba1c_6115/age_sex_model_diags.csv"), 
  dgn3 = read_csv("Data/gsk/age_sex_model_diags.csv"), 
  dgn4 = read_csv("Data/agesexhba1c_8697/age_sex_model_diags.csv"),
  dgn5 = read_csv("Data/agesexhba1cempareg_9492/diag_lm.csv"),
  .id = "src") %>% 
  filter(models == "f1") %>% 
  distinct(src, nct_id, nct_id2, nobs)
dgn %>% 
  filter(nct_id %in% drpd)
getarms %>% 
  anti_join(dgn)
getarms <- getarms %>% 
  inner_join(dgn)
getarms_smry <- getarms %>% 
  group_by(nct_id) %>% 
  summarise(nobs = unique(nobs),
            nobs_slctd = sum(narm)) %>% 
  ungroup()
## note that nobs is for trial not arm
write_csv(getarms, "Outputs/arm_lvl_info_as_analysed.csv")
getn <- getarms %>% 
  distinct(nct_id, nobs) %>% 
  arrange(nct_id)

d <- getn %>% 
  arrange(nct_id) %>% 
  mutate(nct_id_as_n = str_sub(nct_id, 4) %>% as.integer(),
         tst = as.integer(str_sub(nct_id, -2)) * round(nobs, -1))
sum(d$nct_id_as_n)
sum(d$tst)
sum(d$tst[1:50])
sum(d$tst[-(1:50)])

d2 <- d %>% 
  filter(!nct_id %in% drpd)
sum(d2$nct_id_as_n)
sum(d2$nobs)
