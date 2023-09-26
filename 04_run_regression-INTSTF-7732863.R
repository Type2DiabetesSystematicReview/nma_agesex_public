## Run regression models
library(tidyverse)
library(multinma)

## read in simulated IPD
ipd <- readRDS("Scratch_data/simulated_ipd.Rds")
## read in agg
agg <- readRDS("Scratch_data/agg_hba1c.Rds")
## read in arm drug names, doses, regiments
drug <- readRDS("Scratch_data/drug_names_doses_regimen.Rds")
## read in arm assignment
arm_assign <- readRDS("Scratch_data/arm_labels_final.Rds")



## Drop agg if already in ipd, is so for 50 trials ----
intersect(agg$nct_id, ipd$nct_id) %>% length()
agg <- agg %>% 
  filter(!nct_id %in% ipd$nct_id)

## Add arm label to IPD only one not in is where set to placebo
# as expected drop some arms from two trials but keep trials
sum(!duplicated(ipd$nct_id))
ipd <- ipd %>% 
  inner_join(arm_assign)
sum(!duplicated(ipd$nct_id))

## add arm label to aggregate - 4 trials dropped
sum(!duplicated(agg$nct_id))
agg <- agg %>% 
  inner_join(arm_assign) 

## recode data as required
ipd <- ipd %>% 
  rename(result = value_2, base = value_1) %>% 
  select(-age10)


## Drop apparently single arm trials. Need to review
## Seems to be some problem wiht metformin and a missing. Review
rmv1 <- c("NCT01368081")
# rmv1 <- ""
## Two clearly are single arm drugs.
## one appears to be linagliptin versus placebo, with different doses of empa (10 and 25). Appears to be coding error in
## sameclass algorithm. Need to review all "sameclass" algorithm again
ipd %>% 
  filter(nct_id %in% rmv1) %>% 
  count(nct_id, arm_id_unq) %>% 
  left_join(drug) %>% 
  select(nct_id, drug_name:drug2_freq) %>% 
  distinct()
ipd <- ipd %>%
  filter(!nct_id %in% rmv1) 

## need to review some of these again same class algorithm is not quite working.
rmv2 <- c("NCT00316082", "NCT00688701", "NCT00712673", "NCT00763451", "NCT01032629", "ChiCTR1800015296", "jRCTs031180159", "NCT02875821", "TCTR20170511001", "UMIN000021177",
"UMIN000030514", "UMIN000031451", "NCT01042977", "NCT03387683")
agg <- agg %>%
  filter(!nct_id %in% rmv2)
# review NCT00838903 as has negative standard errors AND they are the same as the effect estimates
agg <- agg %>% 
  filter(!nct_id %in% "NCT00838903")
dm_net <- combine_network(
  set_ipd(ipd,
          study = nct_id2,
          trt = drug_code,
          y = result,
          trt_class = trtcls5),
  set_agd_arm(agg,
              study = nct_id,
              trt = drug_code,
              y = result, 
              se = se,
              trt_class = trtcls5,
              sample_size = participants))
pdf("Outputs/Initial_network.pdf", height = 20, width = 20)
plot(dm_net, layout = "auto")
dev.off()
# NCT02477865 and NCT02477969 have placebo arms but no data for these
agg <- agg %>% 
  filter(!nct_id %in% c("NCT04196231",
                        "NCT02477969",
                        "NCT02477865"))
dm_net <- combine_network(
  set_ipd(ipd,
          study = nct_id2,
          trt = drug_code,
          y = result,
          trt_class = trtcls5),
  set_agd_arm(agg,
              study = nct_id,
              trt = drug_code,
              y = result, 
              se = se,
              trt_class = trtcls5,
              sample_size = participants))
plot(dm_net, layout = "auto")
pdf("Outputs/Connected_network.pdf", height = 20, width = 20)
plot(dm_net, layout = "auto")
dev.off()
dm_net <- add_integration(dm_net,
                          age = distr(qnorm, mean = age_m, sd = age_s),
                          n_int = 100)
dm_net_fe <- nma(dm_net,
                 trt_effects = "fixed", link = "identity", likelihood = "normal", 
                 init_r = 0.1,
                 QR = TRUE, cores = 4)
dm_net_fe$code <- rstan::get_stancode(dm_net_fe$object)
saveRDS(dm_net_fe, "Scratch_data/fe_model.Rds")
