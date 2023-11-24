# 06_ipd_trial_psuedoipd
library(tidyverse)
library(multinma)
source("Scripts/00_functions.R")
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
ipd <- tot %>% 
  select(drug_regime_smpl, ipd) %>% 
  unnest(ipd)
rm(tot)


ipd_nst <- ipd %>% 
  group_by(drug_regime_smpl, nct_id, nct_id2, arm_id_unq, arm_f, trtcls4, reference_arm, trtcls5) %>% 
  nest() %>% 
  ungroup()
ipd_nst$n <- map_int(ipd_nst$data, nrow)
ipd_nst$reg <- map(ipd_nst$data, ~ lm(result ~ base + age + sex, data = .x))
ipd_nst$res <- map(ipd_nst$reg, ~ broom::tidy(.x) %>% select(term, estimate, std.error))
ipd_nst <- ipd_nst %>% 
  select(-data, -reg) %>% 
  unnest(res)
## 76 trials. One row per regime trial arm combo
ipd_nst <- ipd_nst %>% 
  filter(term %in% c("age", "sex"))
ipd_nst <- ipd_nst %>% 
  group_by(drug_regime_smpl, term) %>% 
  nest() %>% 
  ungroup()

ipd_nst$nwork <- map(ipd_nst$data, ~ {
  set_agd_arm(.x, 
              study = nct_id2,
              trt = trtcls5, 
              y = estimate,
              se = std.error, 
              sample_size = n, 
              trt_ref = "place", 
              trt_class = trtcls4)
})
rm(ipd)
# ipd_nst$stan_nma <- map(ipd_nst$nwork, ~ nma(.x, runmodel = TRUE, trt_effects = "random", cores = 2, chains = 2))

# res2 <- nma(ipd_nst$nwork[[1]],  runmodel = TRUE, trt_effects = "random", cores = 3)

# saveRDS(ipd_nst, "Scratch_data/ipd_only_res.Rds")
ipd_nst <- readRDS("Scratch_data/ipd_only_res.Rds")
ipd_nst$res <- map(ipd_nst$stan_nma, ~ summary(.x, pars = "d")$summary)
ipd_nst$plt <- map2(ipd_nst$nwork, ipd_nst$drug_regime_smpl, ~ plot(.x) + ggtitle(.y))
ipd_nst$plt 
nwork3 <- cowplot::plot_grid(plotlist = ipd_nst$plt[ipd_nst$term == "age"][c(3,1,2)])
saveRDS(nwork3, "Scratch_data/IPD_nwork.Rds")
res <- ipd_nst %>% 
  select(drug_regime_smpl, term, res) %>% 
  unnest(res)
res <- res %>% 
  mutate(across(c(`mean`, `2.5%`, `97.5%`), ~ if_else(term == "age", 30*.x, .x)))
dc_lkp <- ipd_nst %>% 
  select(drug_regime_smpl, data) %>% 
  unnest(data) %>% 
  distinct(drug_regime_smpl, drug_code, trtcls5)
res <- res %>% 
  mutate(drug_code = str_sub(parameter, 3, -2)) %>% 
  inner_join(dc_lkp) 

plot1 <- ggplot(res, aes(x = paste0(drug_regime_smpl, "_", drug_code),
                         colour = drug_regime_smpl,
                         y = mean, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_point() +
  geom_linerange() +
  facet_grid(trtcls5 ~ term, scales = "free_y") +
  geom_hline(yintercept = 0) +
  coord_flip(ylim = c(-1, 1)) +
  scale_x_discrete("", guide = "none")
saveRDS(plot1, "Scratch_data/pseudo_ipd_tx_compare.Rds")
