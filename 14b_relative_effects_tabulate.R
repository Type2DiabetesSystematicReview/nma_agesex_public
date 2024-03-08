# 14b_relative_effects_tabulate
library(tidyverse)
source("Scripts/04b_arm_label_reverse.R")

mydf <- readRDS("Scratch_data/random_effects_smpls.Rds")
## drop A10BX from hba1c_m6_aggipd_random_triple
mydf <- mydf %>% 
  filter(!(modelname == "hba1c_agesex_m6_aggipd_random_triple" & trtcls5 == "A10BX"))

mtrx <- mydf %>% 
  select(d_value, male_inter, age10_inter) %>% 
  as.matrix()

lpfx <- function(trt, male, age10, mymtrx = mtrx) {
  as.vector(mymtrx %*% c(trt, male, age10))
}
lpfx(1, 1, 5) %>% head()

predictfor <- expand_grid(trt = 1L, male = 0:1, age10 = c(4, 6, 8))
predictforlst <- as.list(predictfor)
predictfor$res <- pmap(predictforlst, function(trt, male, age10) {
  a <- lpfx(trt, male, age10)
  name_a <- paste0(if_else(male ==0, "f", "m"),
                     "_",
                     age10*10)
  a <- tibble(a)
  names(a) <- name_a
  a
  })
mydf <- bind_cols(mydf,
                  bind_cols(predictfor$res))

mydf <- mydf %>% 
  select(modelname, outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg, 
         outcome_rel, trtcls5, d, f_40:m_80) %>% 
  group_by(modelname, outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg, 
           outcome_rel, trtcls5, d) %>% 
  nest() %>% 
  ungroup()

mydf$smry <- map(mydf$data, ~ {
  map(.x, ~ {
  tibble(m = mean(.x),
         s = sd(.x),
         q2_5 = quantile(.x, probs = 0.025),
         q97_5 = quantile(.x, probs = 0.975)) }) %>% 
    bind_rows(.id = "grp")
})

smry <- mydf %>% 
  select(-data) %>% 
  unnest(smry)
smry_hba1c <- smry %>% 
  filter(outcome == "hba1c") %>% 
  mutate(across(c(m, s, q2_5, q97_5), ~ round(.x, 1) %>% formatC(digits = 1, format = "f")),
         res = paste0(m, " (", q2_5, " to ", q97_5, ")")) %>% 
  select(-m, -s, -q2_5, -q97_5) %>% 
  spread(grp, res)

smry_mace <- smry %>% 
  filter(outcome == "mace") %>% 
  mutate(across(c(m, s, q2_5, q97_5), ~ round(exp(.x), 2) %>% formatC(digits = 2, format = "f") %>% 
                  str_replace("0.00", "0.004")),
         res = paste0(m, " (", q2_5, " to ", q97_5, ")")) %>% 
  select(-m, -s, -q2_5, -q97_5) %>% 
  spread(grp, res)

mace_novel <- smry_mace %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  select(network, trtcls5, d, f_40:m_80)  %>% 
  mutate(d = str_sub(d, 3, -2)) %>% 
  separate(d, into = c("drug", "dose"), sep = "_", extra = "merge", fill = "right")

## relabel hba1c WHO ATCs to the names
hba1c_novel <- smry_hba1c %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  select(network, trtcls5, d, f_40:m_80) %>% 
  mutate(d = str_sub(d, 3, -2))  %>% 
  separate(d, into = c("drug", "dose"), sep = "_d", extra = "merge", fill = "right", remove = FALSE) %>% 
  mutate(drug = who_lkp_rev[drug])
rm(mtrx)
gc()

