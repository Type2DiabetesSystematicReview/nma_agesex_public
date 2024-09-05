# 06_fit_hba1c_reg
## Designed to run inside virtual machine for speed using Rscript with arguments
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)

## read in metadata 
mdl_order <- read.csv("Outputs/model_order.csv", stringsAsFactors = FALSE)

## allow passing arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop("Need to pass argument  (eg in terminal via Rscript Scripts/06b_fit_hba1c_reg.R 1) to indicate which model wish to run")
if(length(args) == 2) testrun <- TRUE else testrun <- FALSE
rowchoose <- args[1]

dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")

mdl_chs <- mdl_order[rowchoose,] %>% 
  as.list()
print(paste(as.vector(mdl_chs), collapse = ":"))
filename <- paste(mdl_chs, collapse = "_")

if (mdl_chs$f_mdl == "f1") {
  mymod <- ~ value_1 + .trt
} else if (mdl_chs$f_mdl %in% c("f4", "f8")) {
  mymod <- ~ value_1 + (male + age10)*.trt
}
print(as.character(mymod)[2])

## For running on VM use conditional statements to read/write from top-level folder
## note want to read in data at lowest point in each loop
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")

## Limit to longer duration
if(mdl_chs$durn == "ge26") {
  tot$agg[] <- lapply(tot$agg, function(x) x %>% filter(ge26 ==1L))
  tot$reg[] <- lapply(tot$reg, function(x) x %>% filter(ge26 ==1L))
}

## choose network
chsn_reg <- tot$reg[tot$drug_regime_smpl == mdl_chs$nwork][[1]]
chsn_reg <- chsn_reg %>% 
  mutate(value_1 = case_when(
    is.na(term) ~ 0,
    str_detect(term, "value_1") ~ 1,
    TRUE ~ NA_real_))
chsn_reg <- chsn_reg %>% 
  filter(models == mdl_chs$f_mdl)
chsn_agg <- tot$agg[tot$drug_regime_smpl == mdl_chs$nwork][[1]]
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

## For dual drop two disconnected trials for aggregate NCT02477969 and JapicCTI-101352
## and one for contrast NCT03508323
## For triple one for aggregate - UMIN000007051
## for mono two for aggregate - NCT02477865 and JapicCTI-101351
chsn_agg <- chsn_agg %>% 
  filter(!nct_id %in% dropdisconnect)
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

## Set-up data for network
if (mdl_chs$f_mdl == "f1") {
  regpart <- set_agd_regression(chsn_reg,
                                study = nct_id2,
                                trt = arm_lvl,
                                # trt_ref = "placebo",
                                estimate = estimate,
                                cov = vcv_final,
                                trt_class = trtcls5,
                                regression = ~ value_1 + .trt)
} else if (mdl_chs$f_mdl %in% c("f4", "f8")) {
  regpart <- set_agd_regression(chsn_reg,
                                study = nct_id2,
                                trt = arm_lvl,
                                # trt_ref = "placebo",
                                estimate = estimate,
                                cov = vcv_final,
                                trt_class = trtcls5,
                                regression = ~ value_1 + (male + age10)*.trt)
}


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
if(mdl_chs$data_lvl == "aggipd") {
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
  
} else {
  chsn_net <- combine_network(
    regpart,
    trt_ref = "placebo")
}

if(testrun) {
  mdl <- "somestring" 
  if(!dir.exists("Outputs/nwork_plots")) dir.create("Outputs/nwork_plots")
  ggplot2::ggsave(plot = plot(chsn_net, layout = "auto"),
         filename = paste0("Outputs/nwork_plots/", filename, ".pdf"), height = 20, width = 20)
  } else {
    if (mdl_chs$f_mdl == "f1") {
    mdl <- nma(chsn_net,
               trt_effects = mdl_chs$fe_re,
               regression = ~ value_1 + .trt,
               class_interactions = "common",
               prior_intercept = normal(scale = 10),
               prior_trt = normal(scale = 10),
               prior_reg = normal(scale = 10), chains = 4, cores = 4)
    } else if (mdl_chs$f_mdl %in% c("f4", "f8")) {
      mdl <- nma(chsn_net,
                 trt_effects = mdl_chs$fe_re,
                 regression = ~ value_1 + (male + age10)*.trt,
                 class_interactions = "common",
                 prior_intercept = normal(scale = 10),
                 prior_trt = normal(scale = 10),
                 prior_reg = normal(scale = 10), chains = 4, cores = 4)
    }
  }
saveRDS(mdl, paste0(filename, ".Rds"))
rm(mdl)
gc()
