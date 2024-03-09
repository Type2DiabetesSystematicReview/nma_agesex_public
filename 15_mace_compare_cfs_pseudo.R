# 14_compare_cfs_pseudo
library(broom)
library(tidyverse)
library(survival)

## test cfs against pseudo IPD as an independent check on export

## Read in coefficients and pseudodata and arrange into similar format with similar names ----
f1 <- readRDS("Scratch_data/for_mace_regression_nointer.Rds")
pseudo <- f1$pseudo
f1 <- f1$cfs
## pseudo same for all 3
# identical(f1$pseudo, f5$pseudo)
# identical(f1$pseudo, f2$pseudo)
f2 <- readRDS("Scratch_data/for_mace_regression_nointercovs.Rds")$cfs
f5 <- readRDS("Scratch_data/for_mace_regression_inter.Rds")$cfs
cfs <- bind_rows(f1 = f1,
                 f2 = f2,
                 f5 = f5)
cfs <- cfs %>% 
  mutate(term = str_replace(term, "sexMALE", "sexM"))

rm(f1, f2, f5)

## note that IPD arm corresponds to cfs "trt" which was used to create arm_f
pseudo <- pseudo %>% 
  mutate(arm_f = if_else(str_sub(ipd_arm) == "placebo", paste0("a_", ipd_arm), ipd_arm),
         age10c = (age-60)/10,
         sex = if_else(male ==1, "M", "F")) %>% 
  group_by(nct_id) %>% 
  nest() %>% 
  ungroup()
pseudo$data <- map(pseudo$data, ~ .x %>% 
                     mutate(arm_f = factor(arm_f)))

## create formulas ----
fs <- list(f1 =  Surv(time, event) ~ arm_f,
     f2 = Surv(time, event) ~ arm_f + age10c + sex,
     f5 = Surv(time, event) ~ age10c*arm_f + sex*arm_f)

## Run models on pseudodata and obtain coefficeints ----
res <- map(fs, function(frm) {
  map(pseudo$data, ~ coxph(formula = frm, data = .x))
})
res <- as_tibble(res)
pseudo <- bind_cols(pseudo, res)
rm(res)
pseudo <- pseudo %>% 
  gather("frm", "mdl", -nct_id, -data)
pseudo$cfs <- map(pseudo$mdl, tidy)
pseudo <- pseudo %>% 
  select(-data, -mdl) %>% 
  unnest(cfs)


## Check terms match and then join ----
intersect(cfs$term, pseudo$term)
setdiff(cfs$term, pseudo$term)
setdiff(pseudo$term, cfs$term)
pseudo <- pseudo %>% 
  rename(estimate_psd = estimate,
         std.error_psd = std.error) %>% 
  inner_join(cfs %>% select(nct_id, reference, repo, frm = models, term, estimate_cfs = estimate, std.error_cfs = std.error, trtcls5, trt, arm_lvl))
## relabel so using drug name
cleanup <- pseudo %>% 
  filter(!is.na(arm_lvl)) %>% 
  distinct(trt, arm_lvl)
for (i in seq_along(cleanup$trt)) {
  pseudo$term <- str_replace(pseudo$term, cleanup$trt[[i]], cleanup$arm_lvl[[i]])
}
pseudo <- pseudo %>% 
  mutate(term = str_remove(term, "arm_f"))
## flip so that age10c in interaction comes after arm for plotting
pseudo <- pseudo %>% 
  mutate(term = 
           if_else(str_detect(term, "^age10c\\:"),
                   paste0(str_remove(term, "^age10c\\:"), ":age10c"),
                   term))
## plot differences ----
pseudo_lng <- pseudo %>% 
  gather("measure_src", "value", estimate_psd, estimate_cfs, std.error_psd, std.error_cfs) %>% 
  separate(measure_src, into = c("measure", "src"), sep = "_") %>% 
  spread(measure, value) %>% 
  mutate(lci = estimate - 1.96*std.error,
         uci = estimate + 1.96*std.error)

psd_cfs_plt <- ggplot(pseudo_lng, aes(x = interaction(nct_id, term), y = estimate, ymin = lci, ymax = uci, colour = src)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap( ~ frm, scales = "free", ncol = 3) 
psd_cfs_plt
pdf("Outputs/compare_hr_pseudo_cfs.pdf", width = 20, height = 15)
psd_cfs_plt
dev.off()

## review correlations
