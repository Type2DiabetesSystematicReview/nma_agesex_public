# 06_fit_set_reg
library(tidyverse)
library(multinma)

tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")

dual_reg <- tot$reg_f4[tot$drug_regime_smpl == "dual"][[1]]
## to check code runs drop the three trials in multiple settings. Will need to update
dual_reg <- dual_reg %>% 
  filter(!nct_id %in% c("NCT00798161", "NCT01023581", "NCT01890122")) %>% 
  slice(1:24)

dual_agg <- tot$agg[tot$drug_regime_smpl == "dual"][[1]]
dual_agg <- dual_agg %>% 
  mutate(nct_id2 = nct_id)
terms <- dual_reg %>% 
  filter(!is.na(arm_id_unq)) %>% 
  distinct(nct_id2, arm_id_unq)
terms <- terms %>% 
  nest(data = arm_id_unq) %>% 
  rename(terms = data)
crl <- dual_reg %>% 
  distinct(nct_id2, crl)
crl <- crl %>% 
  inner_join(terms)
crl$terms <- map(crl$terms, ~ .x$arm_id_unq)
crl$crl2 <- map2(crl$crl, crl$terms, function(mycor, terms) {
  mycols <- colnames(mycor)
  termsrch <- paste(terms, collapse = "|")
  chs <- !str_detect(mycols, "^arm_f") | str_detect(mycols, termsrch)
  mycor[chs, chs]
})

crl_final <- crl$crl2
names(crl_final) <- crl$nct_id2
dual_net <- combine_network(
  set_agd_regression(dual_reg,
                     study = nct_id2,
                     trt = drug_code,
                     estimate = estimate,
                     se = std.error,
                     cor = crl_final,
                     # cov = pso_reg_cov,
                     trt_class = trtcls5,
                     regression = ~ (sex + age10)*.trt))


dual_FE <- nma(dual_net,
                  trt_effects = "fixed",
                  regression = ~ (sex + age10)*.trt,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10))
