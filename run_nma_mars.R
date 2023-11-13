# rn_nma_mars.R

## functions
RptNetwork <- function(ipd_choose, agg_choose){
  combine_network(
    set_ipd(ipd_choose,
            study = nct_id2,
            trt = drug_code,
            y = result,
            trt_class = trtcls5),
    set_agd_arm(agg_choose %>% 
                  filter(treat_or_ref == "arm_level_outcome"),
                study = nct_id,
                trt = drug_code,
                y = result, 
                se = se,
                trt_class = trtcls5,
                sample_size = n),
    set_agd_contrast(agg_choose %>% 
                       filter(!treat_or_ref == "arm_level_outcome"),
                     study = nct_id,
                     trt = drug_code,
                     y = result, 
                     se = se,
                     trt_class = trtcls5,
                     sample_size = n))  
}


install.packages("multinma", repos = c("https://dmphillippo.r-universe.dev", 
                                       "https://cloud.r-project.org/"))
library(multinma)

tot <- readRDS("data_for_mars.Rds")

ipd <- tot$ipd[[1]]
agg <- tot$agg[[1]]
rm(tot)

print(nrow(ipd))
print(nrow(agg))

dm_nets <- RptNetwork(ipd, agg)
dm_nets <- add_integration(dm_nets, age = distr(qnorm, mean = age_m, sd = age_sd))
dm_net_fe <- nma(dm_nets, 
                                    regression = ~ age*.trt,
                                    trt_effects = "fixed", link = "identity", likelihood = "normal", 
                                    init_r = 0.1,
                                    QR = TRUE, cores = 4)

saveRDS(dm_net_fe, "nma_res.Rds")
