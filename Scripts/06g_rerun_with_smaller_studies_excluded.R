library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(multinma)
# 
# xmn <- readRDS("Scratch_data/exclude_by_size.Rds")
# slct <- xmn$armbycut %>% filter(nwork == "triple", grp == "arms_connected_to_placebo") %>% distinct(ncut, arm_lvl) 
# slct <- slct %>% 
#   mutate(ncut = str_remove(ncut, "ncut") %>% as.integer)
# slct <- slct %>% 
#   filter(ncut %in% c(50, 100, 150, 200, 300))
# write_csv(slct, "Outputs/select_cutpoints.csv")
slct <- read_csv("Outputs/select_cutpoints.csv")

## allow passing arguments
if(interactive()) {
  cutchoose <- 50
} else {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args) == 0) stop("Need to pass argument")
  cutchoose <- args[1]
}

ncut_arms <- slct %>% 
  filter(ncut == cutchoose) %>% 
  distinct(arm_lvl) %>% 
  pull(1)

mymod <- ~ value_1 + (male + age10)*.trt
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
tot <- tot %>% 
  filter(drug_regime_smpl == "triple") %>% 
  select(-drug_regime_smpl)
chsn_agg <- tot$agg[[1]] %>% 
    group_by(nct_id) %>% 
    filter(any(arm_lvl %in% ncut_arms)) %>% 
    ungroup()
chsn_reg <- tot$reg[[1]] %>% 
  group_by(nct_id) %>% 
  filter(any(arm_lvl %in% ncut_arms)) %>% 
  ungroup()

chsn_reg <- chsn_reg %>% 
  mutate(value_1 = case_when(
    is.na(term) ~ 0,
    str_detect(term, "value_1") ~ 1,
    TRUE ~ NA_real_))
chsn_reg <- chsn_reg %>% 
  filter(models == "f4")
chsn_agg <- chsn_agg %>% 
  mutate(nct_id2 = nct_id) %>% 
  mutate(male_p = male/n) %>% 
  select(-male)
crl <- chsn_reg %>% 
  filter(is.na(term))
crl_final <- crl$crl
names(crl_final) <- crl$nct_id2
rm(crl)
vcv <- chsn_reg %>% 
  filter(is.na(term))
vcv_final<- vcv$vcv
names(vcv_final) <- vcv$nct_id2
rm(vcv)

chsn_reg <- chsn_reg %>% 
  select(nct_id2, term, estimate, std.error, arm_lvl, trtcls5, age10, sex, value_1)  %>% 
  mutate(male = case_when(
    sex == TRUE ~ 1L,
    sex == FALSE ~ 0L,
    TRUE ~ NA_integer_))

chsn_agg <- chsn_agg %>% 
  group_by(nct_id) %>% 
  mutate(cntrst = any(is.na(result))) %>% 
  ungroup()
chsn_cntrst <- chsn_agg %>% 
  filter(cntrst)
chsn_agg <- chsn_agg %>% 
  filter(!cntrst)
sd_plac <- chsn_agg %>% 
  filter(arm_lvl == "placebo") %>% 
  mutate(sd = se*n^0.5) %>% 
  pull(sd) %>% 
  median()
chsn_cntrst <- chsn_cntrst %>% 
  mutate(se = if_else(is.na(se) & nct_id == "NCT01744236" & arm_lvl == "placebo", sd_plac/n^0.5, se))

regpart <- set_agd_regression(chsn_reg,
                              study = nct_id2,
                              trt = arm_lvl,
                              # trt_ref = "placebo",
                              estimate = estimate,
                              cov = vcv_final,
                              trt_class = trtcls5,
                              regression = ~ value_1 + (male + age10)*.trt)
armpart <- set_agd_arm(chsn_agg,
                       study = nct_id,
                       trt = arm_lvl,
                       trt_ref = "placebo",
                       y = result,
                       se = se,
                       sample_size = n,
                       trt_class = trtcls5)
cntrstpart <- set_agd_contrast(chsn_cntrst,
                               study = nct_id,
                               trt = arm_lvl,
                               # trt_ref = "placebo",
                               y = result,
                               se = se,
                               sample_size = n,
                               trt_class = trtcls5)
chsn_net <- combine_network(
  regpart,
  armpart,
  cntrstpart,
  trt_ref = "placebo")

if (nrow(chsn_cntrst) >=  1) {
  chsn_net <- combine_network(
    regpart,
    armpart,
    cntrstpart,
    trt_ref = "placebo")
} else {
  chsn_net <- combine_network(
    regpart,
    armpart,
    trt_ref = "placebo")
}

mycor <- matrix(rep(0, 9), nrow = 3)
diag(mycor) <- 1
chsn_net <- add_integration(chsn_net,
                            age10 = distr(qnorm, mean = age_m, sd = age_sd),
                            male = distr(qbern, prob = male_p),
                            value_1 = distr(qnorm, value_1, value_1_sd), n_int = 1L, cor = mycor)
mdl <- nma(chsn_net,
           trt_effects = "fixed",
           regression = ~ value_1 + (male + age10)*.trt,
           class_interactions = "common",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), chains = 4, cores = 4)
