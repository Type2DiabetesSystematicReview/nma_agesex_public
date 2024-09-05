# 14b_relative_effects_tabulate
library(tidyverse)
source("Scripts/04b_arm_label_reverse.R")

mydf <- readRDS("Scratch_data/random_effects_smpls.Rds")


hba1c <- mydf %>% 
  filter(outcome == "hba1c")
mtrx_hba1c <- hba1c %>% 
  select(d_value, male_inter, age_inter) %>% 
  as.matrix()
mace <- mydf %>% 
  filter(outcome == "mace")
mtrx_mace <- mace %>% 
  select(d_value, male_inter, age_inter) %>% 
  as.matrix()
rm(mydf)

lpfx <- function(trt, male, age10, mymtrx) {
  as.vector(mymtrx %*% c(trt, male, age10))
}
lpfx(1, 1, 50-60, mtrx_mace) %>% head()

predictfor <- expand_grid(trt = 1L, male = 0:1, age = c(40, 60, 80))
predictforlst <- as.list(predictfor)
predictfor$hba1c <- pmap(predictforlst, function(trt, male, age) {
  lpfx(trt, male, age, mtrx_hba1c)
  })
## Note for mace centring age first
predictfor$mace <- pmap(predictforlst, function(trt, male, age) {
  lpfx(trt, male, age-60, mtrx_mace)
})
rm(mtrx_hba1c, mtrx_mace)
gc()

predictfor$mace <- map(predictfor$mace, ~ mace %>% 
                         mutate(res = .x))
predictfor$hba1c <- map(predictfor$hba1c, ~ hba1c %>% 
                         mutate(res = .x))

## get summary stats for hba1c
hba1c <- predictfor %>% 
  rename(male_predict = male, age_predict = age) %>% 
  select(-mace) %>% 
  unnest(hba1c) %>% 
  filter(!trtcls5 == "A10BX") %>% 
  group_by(modelname, outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg, 
           outcome_rel, trtcls5, d,
           trt, male_predict, age_predict) %>% 
  nest() %>% 
  ungroup()
hba1c$smry <- map(hba1c$data, ~ {
  .x %>% 
    summarise(
      m = mean(res),
      s = sd(res),
      q2_5 = quantile(res, probs = 0.025),
      q97_5 = quantile(res, probs = 0.975))
    })
hba1c <- hba1c %>% 
  select(-data) %>% 
  unnest(smry)

## get summary stats for MACE
mace <- predictfor %>% 
  rename(male_predict = male, age_predict = age) %>% 
  select(-hba1c) %>% 
  unnest(mace) %>% 
  group_by(modelname, outcome, mainorinter, modelnum, datalevel, fixedrand, network, sg, 
           outcome_rel, trtcls5, d,
           trt, male_predict, age_predict) %>% 
  nest() %>% 
  ungroup()
mace$smry <- map(mace$data, ~ {
  .x %>% 
    summarise(
      m = mean(res),
      s = sd(res),
      q2_5 = quantile(res, probs = 0.025),
      q97_5 = quantile(res, probs = 0.975))
})
mace <- mace %>% 
  select(-data) %>% 
  unnest(smry)

mace_plot <- ggplot(mace %>% 
                      filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
                             mutate(sex = if_else(male_predict ==0, "Female", "Male"),
                                    dc = case_when(
                                      trtcls5 == "A10BH" ~ "A10BH - DPP4",
                                      trtcls5 == "A10BJ" ~ "A10BJ - GLP1",
                                      trtcls5 == "A10BK" ~ "A10BK - SGLT2")) %>% 
                      group_by(dc) %>% 
                      mutate(d_recode = cumsum(!duplicated(d)) %>% as.character()) %>% 
                      ungroup(),
                    aes(x = age_predict, y = m, ymin = q2_5, ymax = q97_5, 
                              group = interaction(sex, d), linetype = sex, colour = d)) +
  geom_point(position = position_dodge(10)) +
  geom_linerange(position = position_dodge(10)) +
  facet_wrap(~dc) +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") 
mace_plot
pdf("Outputs/relative_mace.pdf")
mace_plot
dev.off()

hba1c <- hba1c %>% 
  mutate(across(c(m, s, q2_5, q97_5), ~ round(.x, 1) %>% formatC(digits = 1, format = "f")),
         res = paste0(m, " (", q2_5, " to ", q97_5, ")")) %>% 
  select(-m, -s, -q2_5, -q97_5) 
mace <- mace %>% 
  mutate(across(c(m, s, q2_5, q97_5), ~ round(exp(.x), 2) %>% formatC(digits = 2, format = "f") %>% 
                  str_replace("0.00", "0.004")),
         res = paste0(m, " (", q2_5, " to ", q97_5, ")")) %>% 
  select(-m, -s, -q2_5, -q97_5) 

## arrange into formatting for paper
mace_novel <- mace %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  select(network, trtcls5, male_predict, age_predict, d, res)  %>% 
  mutate(d = str_sub(d, 3, -2)) %>% 
  separate(d, into = c("drug", "dose"), sep = "_", extra = "merge", fill = "right")
mace_novel <- mace_novel %>% 
  mutate(male_predict = if_else(male_predict ==0L, "Female", "Male")) %>% 
  pivot_wider(names_from = c(male_predict, age_predict), values_from = res) %>% 
  arrange(trtcls5, drug, dose)

## relabel hba1c WHO ATCs to the names
hba1c_novel <- hba1c %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK")) %>% 
  select(network, trtcls5, male_predict, age_predict, d, res)  %>% 
  mutate(d = str_sub(d, 3, -2))  %>% 
  separate(d, into = c("drug", "dose"), sep = "_d", extra = "merge", fill = "right", remove = FALSE) %>% 
  mutate(drug = who_lkp_rev[drug])
hba1c_novel <- hba1c_novel %>% 
  mutate(male_predict = if_else(male_predict ==0L, "Female", "Male")) %>% 
  pivot_wider(names_from = c(male_predict, age_predict), values_from = res) %>% 
  arrange(trtcls5, drug, dose)

write_csv(hba1c_novel, "Outputs/hba1c_relative_effects.csv")
write_csv(mace_novel, "Outputs/mace_relative_effects.csv")
