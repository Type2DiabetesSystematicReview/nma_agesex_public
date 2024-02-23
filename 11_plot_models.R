# 06b_read_models_vm

library(tidyverse)
library(multinma)

if(sessionInfo()$platform == "x86_64-pc-linux-gnu (64-bit)") {
  whoatc <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx", sheet = 1) 
} else {
  whoatc <- readxl::read_excel("../../../Medications_resources/WHO_ATC/2018 ATC index with DDDs.xlsx", sheet = 1)
}
whoatc <- whoatc %>% 
  select(trtclass = `ATC code`,
         cls = `ATC level name`) %>% 
  distinct()

mdl_names <- list.files("FromVM/", patt = "Rds$")
res <- map(mdl_names, ~ readRDS(paste0("FromVM/", .x)))
names(res) <- mdl_names %>% str_sub(1, -5)
res$mace <- readRDS("Scratch_data/mace_h_c_model.Rds")

beta <- map(res, ~ summary(.x$stanfit, pars = "beta"))
beta <- map(beta, ~ .x$summary)
beta <- map(beta, ~ as_tibble(.x, rownames = "params"))
beta <- bind_rows(beta, .id = "tosep")
beta <- beta %>% 
  mutate(tosep = if_else(tosep == "mace", "fe_model_mace", tosep))
beta <- beta %>% 
  separate(tosep, into = c("fixedrand", "model", "network"), sep = "_")
write_csv(beta, "Outputs/hba1c_meta_analysis.csv")
divergent <- map(res, ~ get_sampler_params(.x$stanfit, inc_warmup = FALSE))
divergent <- map(divergent, ~ map_int(.x, ~ sum(.x[, "divergent__"])))
## no divergent transitions
rm(res)
beta_age_sex <- beta %>% 
  filter(params %>% str_detect("age|male"),
         params %>% str_detect("\\:")) %>% 
  mutate(params = params %>% 
           str_remove("^beta\\[") %>% 
           str_remove("\\]") %>%
           str_remove("\\.trtclass")) %>% 
  separate(params, into = c("covariate", "trtclass"), sep = "\\:")
beta_age_sex <- beta_age_sex %>% 
  inner_join(whoatc) %>% 
  mutate(cls = paste0(trtclass, ":", cls))
beta_age_sex <- beta_age_sex %>% 
  janitor::clean_names() %>% 
  mutate(network = factor(network, levels = c("mono", "dual", "triple")))
beta_age_sex <- beta_age_sex %>% 
   mutate(myalpha = case_when(
     (trtclass == "A10BX") ~ 0.2,
     (trtclass == "A10A" & network == "mono") ~ 0.2,
     (trtclass %in% c("A10BX", "A10BA") & network == "triple") ~ 0.2,
     TRUE ~ 1),
     across(mean:x97_5_percent, ~ case_when(
       covariate == "age10" ~ .x*3,
       covariate == "age15" ~ .x*2,
       TRUE ~ .x))) %>% 
  mutate(covariate = if_else(covariate %in% c("age10", "age15"), "age30", covariate))

plotbeta <- ggplot(beta_age_sex %>% filter(!is.na(network)), 
                   aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, colour = fixedrand, alpha = myalpha)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(covariate ~ network) + 
  # scale_y_continuous(limits = c(-1, 1)) +
  scale_x_discrete(limits = rev) +
  scale_alpha_identity() +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("HbA1c (mmol/mol")
plotbeta

plotmace <- ggplot(beta_age_sex %>% filter(is.na(network)), 
                   aes(x = cls, y = mean, ymin = x2_5_percent, ymax = x97_5_percent, colour = fixedrand)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_grid(~ covariate) + 
  scale_x_discrete(limits = rev) +
  coord_flip(ylim = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  scale_y_continuous("Hazard ratio (log-scale)")
plotmace

pdf("Outputs/hba1c_mace.pdf", height = 10, width = 20)
plotbeta + ggtitle("HbA1c meta-analysis")
plotmace + ggtitle("MACE meta-analysis")
dev.off()




