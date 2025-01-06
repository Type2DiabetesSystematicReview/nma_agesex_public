# 09b_run_mace
library(dplyr)
library(tidyr)
library(multinma)
library(truncnorm)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0) stop("Need to pass arguments (eg in terminal via Rscript Scripts/10c_run_mace_inter.R ARGS HERE) to indicate which model wish to run")

refe <- args[[1]]
filename <- paste0("Scratch_data/mace_class_nocovs_", refe, ".Rds")
print(filename)

list2env(readRDS("Scratch_data/for_mace_regression_noipd.Rds"), envir = .GlobalEnv)

cfs <- cfs %>% 
  rename(age = age10c)
# mace_agg_from_cfs <- cfs %>% 
#   filter(!is.na(term)) %>% 
#   group_by(nct_id, trtcls5) %>% 
#   summarise(estimate = weighted.mean(estimate, 1/std.error^2),
#             std.error = sum(std.error^2)^0.5) %>% 
#   ungroup()
mace_agg <- mace_agg %>% 
  # filter(!is.na(loghr)) %>% 
  group_by(nct_id, trtcls5) %>% 
  summarise(estimate = weighted.mean(loghr, 1/se^2),
            std.error = sum(se^2)^0.5,
            participants = sum(participants)) %>% 
  ungroup()

nwork <- combine_network(set_agd_contrast(data = mace_agg,
                                            study = nct_id, trt = trtcls5, y = estimate, se = std.error, 
                                            trt_ref = "place", sample_size = participants))

## Add integration points. Note taking correlations from pseudo IPD networks

mdl <- nma(nwork,
           trt_effects = refe,
           link = "identity",
           prior_intercept = normal(scale = 10),
           prior_trt = normal(scale = 10),
           prior_reg = normal(scale = 10), 
           chains = 4, cores = 4,
           control = list(max_treedepth = 15))
saveRDS(mdl, filename)

## create plots
mdl_fe <- readRDS("Scratch_data/mace_class_nocovs_fixed.Rds")
mdl_re <- readRDS("Scratch_data/mace_class_nocovs_random.Rds")

res_fe <- summary(mdl_fe)
res_re <- summary(mdl_re)
library(tidyverse)
res <- bind_rows(fixed = res_fe$summary,
                 random = res_re$summary %>% filter(parameter %>% str_detect("^d\\[")),
                 .id = "mdl")
res <- res %>% 
  gather("var", "value", -mdl, -parameter) 
res <- res %>% 
  filter(var %in% c("mean", "2.5%", "97.5%")) %>% 
  mutate(var = str_replace(var, "\\.", "_") %>% 
           str_remove("\\%"))
res <- res %>% 
  mutate(drgcls = case_match(parameter,
                             "d[A10BB]" ~ "Metformin",
                             "d[A10BH]" ~ "DPP4",
                             "d[A10BJ]" ~ "GLP1",
                             "d[A10BK]" ~ "SGLT2"))
res <- res %>% 
  mutate(var = if_else(var == "mean", "estimate", paste0("q", var)))
res <- res %>% 
  spread(var, value)
source("Scripts/00_functions.R")
mylbl <- c(0.67, 0.8, 1, 1.25, 1.5)
plt_new <- ggplot(res %>% filter(!drgcls == "Metformin"), aes(x = drgcls, y = estimate, ymin = q2_5, ymax = q97_5, colour = mdl)) +
  geom_point(position = position_dodge(0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  geom_linerange(position = position_dodge(0.1)) +
  scale_x_discrete("", limits = rev) +
  scale_y_continuous("", limits = log(c(0.6, 1/0.6)),breaks = log(mylbl), labels = mylbl) +
  scale_color_discrete("") +
  coord_flip() +
  theme_minimal5()
plt_new
tiff("Outputs/fig_mace_nocov.tiff", res = 300, compression = "lzw", unit = "in", height = 8, width = 6)
plt_new
dev.off()
