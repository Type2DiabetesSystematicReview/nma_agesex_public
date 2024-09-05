# ssensitivity analysis confined to mace trials (14 have mace and ipd)

library(tidyverse)
library(multinma)
library(rstan)
hba1c_m6_main <- readRDS("FromVM/hba1c_agesex/m12_aggipd_random_triple_f4_ge26.Rds")
hba1c_m6_sens <- readRDS("FromVM/hba1c_in_mace_trials.Rds")

hba1c_m6_main <- as.data.frame(hba1c_m6_main$stanfit)
hba1c_m6_sens <- as.data.frame(hba1c_m6_sens$stanfit)

hba1c_m6_main <- hba1c_m6_main[str_detect(names(hba1c_m6_main), "beta")]
hba1c_m6_sens <- hba1c_m6_sens[str_detect(names(hba1c_m6_sens), "beta")]

hba1c_m6_main <- hba1c_m6_main[ , "beta[age10:.trtclassA10BK]"]
hba1c_m6_sens <- hba1c_m6_sens[ , "beta[age10:.trtclassA10BK]"]

res <- tibble(main = hba1c_m6_main,
              sens = hba1c_m6_sens) %>% 
  mutate(iter = 1:4000)
res <- res %>% 
  gather("model", "value", - iter)
rm(hba1c_m6_main, hba1c_m6_sens)

res <- res %>% 
  mutate(value = value * 30)
plot1 <- ggplot(res %>% 
                  mutate(model = if_else(model == "main", "Main", "Sensitivity")), aes(x = value, fill = model)) +
  geom_density(alpha = 0.1) +
  geom_vline(xintercept = 0L, linetype = "dashed") +
  scale_x_continuous("Hba1c effect (% units) per 30-year increment in age\n(Higher number indicates worse treatment effect)") +
  scale_y_continuous("Probability density") +
  scale_fill_discrete("") 
plot1
pdf("Outputs/hba1c_sensitivity_mace.pdf", page = "a4r")
plot1
dev.off()
tiff("Outputs/hba1c_sensitivity_mace.tiff", height = 7, width = 7, unit = "in", compression = "lzw", res = 300)
plot1
dev.off()
## Also look for MACE
cmpr <- list(
  sens_fe = readRDS("FromVM/mace_withhba1c/mace_sensitivity_fe_hba1conly.Rds"),
  sens_re = readRDS("FromVM/mace_withhba1c/mace_sensitivity_re_hba1conly.Rds"),
  main_fe = readRDS("FromVM/mace_agesex/fixed_mace_agesex_sex.Rds"),
  main_re = readRDS("FromVM/mace_agesex/random_mace_agesex_sex.Rds"))
cmpr <- map(cmpr, ~ as.data.frame(.x$stanfit))
cmpr <- map(cmpr, ~ .x[str_detect(names(.x), "beta")])
cmpr <- map(cmpr, ~ .x[ , "beta[age:.trtclassA10BK]"])
cmpr <- bind_rows(cmpr, .id = "tosep")
cmpr <- cmpr %>% 
  mutate(iter = 1:4000) %>% 
  gather("model", "value", - iter)
res <- cmpr
rm(cmpr)
res <- res %>% 
  mutate(value = value * 30)
res <- res %>% 
  separate(model, into = c("model", "fere"), sep = "_")
plot2 <- ggplot(res %>% 
                  mutate(model = if_else(model == "main", "Main", "Sensitivity"),
                         fere = if_else(fere == "fe", "Fixed effects", "Random effects")),
                aes(x = value, fill = model)) +
  geom_density(alpha = 0.1) +
  geom_vline(xintercept = 0L, linetype = "dashed") +
  scale_x_continuous("MACE effect (hazard ratio) per 30-year increment in age\n(Higher number indicates worse treatment effect)") +
  scale_y_continuous("Probability density") +
  scale_fill_discrete("") +
  facet_wrap(~fere, ncol = 1)
plot2
pdf("Outputs/mace_sensitivity_hba1c.pdf", page = "a4")
plot2
dev.off()

tiff("Outputs/mace_sensitivity_hba1c.tiff", height = 7, width = 7, unit = "in", compression = "lzw", res = 300)
plot2
dev.off()
res %>% 
  group_by(model, fere) %>% 
  summarise(p = mean(value < 0)) %>% 
  ungroup() %>% 
  spread(model, p)

saveRDS(list(plot1, plot2), "Scratch_data/mace_hba1c_confined_both.Rds")
