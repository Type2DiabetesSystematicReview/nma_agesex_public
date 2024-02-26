# 06_fit_hba1c_reg
library(dplyr)
library(tidyr)
library(stringr)
library(multinma)

testrun <- FALSE
dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
exclude <- dropdisconnect
chsn <- "global"
sens <- "full"
for(mdl_type in c("fixed", "random")){
       ## For running on VM use conditional statements to read/write from top-level folder
      ## note want to read in data at lowest point in each loop
      if(dir.exists("Scratch_data/")){
        tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
      } else {
        tot <- readRDS("agg_ipd_hba1c.Rds")
      }
      
      print(paste(sens, mdl_type, chsn, sep = ":"))
      ## Run once after excluding trials from dual with a high Rhat
      if(sens == "exclude") {
        if(mdl_type == "random" & chsn == "dual") {
          problem_trials <- c(
            "NCT00622284", "NCT00856284", "NCT00968812", "NCT01167881", 
            "NCT01289119","NCT01624259", "NCT01644500")
          tot$agg[] <- lapply(tot$agg, function(x) x %>% filter(!nct_id %in% problem_trials))
          tot$ipd[] <- lapply(tot$ipd, function(x) x %>% filter(!nct_id %in% problem_trials))
          tot$reg[] <- lapply(tot$reg, function(x) x %>% filter(!nct_id %in% problem_trials)) } else {
            next()
          }
      }
      print(tot)
      chsn_reg <- tot %>% 
        select(drug_regime_smpl, reg) %>% 
        unnest(cols = c(reg))
      ## need to set first row to the reference for baseline value and ones with baseline value to 1L
      chsn_reg <- chsn_reg %>% 
        mutate(value_1 = case_when(
          is.na(term) ~ 0,
          str_detect(term, "value_1") ~ 1,
          TRUE ~ NA_real_))
      chsn_reg <- chsn_reg %>% 
        filter(models == "f4")
      chsn_agg <- tot %>% 
        select(drug_regime_smpl, agg) %>% 
        unnest(cols = c(agg))
      chsn_agg <- chsn_agg %>% 
        mutate(nct_id2 = nct_id) %>% 
        mutate(male_p = male/n) %>% 
        select(-male)
      crl <- chsn_reg %>% 
        filter(is.na(term))
      crl_final <- crl$crl
      names(crl_final) <- crl$nct_id2
      rm(crl)
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
      chsn_net <- combine_network(
        set_agd_regression(chsn_reg,
                           study = nct_id2,
                           trt = arm_lvl,
                           trt_ref = "placebo",
                           estimate = estimate,
                           se = std.error,
                           cor = crl_final,
                           # cov = pso_reg_cov,
                           trt_class = trtcls5,
                           regression = ~ value_1 + (male + age10)*.trt),
        set_agd_arm(chsn_agg,
                    study = nct_id,
                    trt = arm_lvl,
                    trt_ref = "placebo",
                    y = result,
                    se = se,
                    sample_size = n,
                    trt_class = trtcls5),
        set_agd_contrast(chsn_cntrst,
                         study = nct_id,
                         trt = arm_lvl,
                         trt_ref = "placebo",
                         y = result,
                         se = se,
                         sample_size = n,
                         trt_class = trtcls5))
      
      mycor <- matrix(rep(0, 9), nrow = 3)
      diag(mycor) <- 1
      chsn_net <- add_integration(chsn_net,
                                  age10 = distr(qnorm, mean = age_m/10, sd = age_sd/10),
                                  male = distr(qbern, prob = male_p),
                                  value_1 = distr(qnorm, 0, 0), n_int = 1L, cor = mycor)
      # plot(chsn_net)
      
      if(testrun) {
        mdl <- "somestring" } else {
          mdl <- nma(chsn_net,
                     trt_effects = mdl_type,
                     regression = ~ value_1 + (male + age10)*.trt,
                     class_interactions = "common",
                     prior_intercept = normal(scale = 10),
                     prior_trt = normal(scale = 10),
                     prior_reg = normal(scale = 10), chains = 4, cores = 4)
        }
      
      
      filename <- paste(mdl_type, "hba1c", chsn, sens, sep = "_")
      if(dir.exists("Scratch_data/")) {
        saveRDS(mdl, paste0("Scratch_data/", filename, ".Rds"))
      } else {
        saveRDS(mdl, paste0(filename, ".Rds"))
      }
      rm(mdl)
      gc()
    }
  

