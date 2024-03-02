library(tidyverse)
library(truncnorm)
source("../common_functions/Scripts/truncated_normal.R")

## read and transform data ----
who <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx") 
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
mace_agg <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_agg
elig <- read_csv("../cleaned_data/Data/age_max_min_elig.csv")
source("../common_functions/Scripts/combine_sd.R")
agg <- tot %>% 
  select(drug_regime_smpl, agg) %>% 
  unnest(agg)
pseudo_hba1c <- tot %>% 
  select(drug_regime_smpl, ipd) %>% 
  unnest(ipd) %>% 
  mutate(sex = if_else(sex == 1L, "M", "F"))
pseudo_mace <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
mace_arms <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_arms
pseudo_mace <- pseudo_mace %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, arm_lvl, trtcls5))
maceinhba1c_ipd <- intersect(pseudo_mace$nct_id, pseudo_hba1c$nct_id)
pseudo_hba1c <- pseudo_hba1c %>% 
  anti_join(pseudo_mace %>% select(nct_id))
ipd <- bind_rows(hba1c = pseudo_hba1c %>% select(nct_id, sex, age, arm_lvl, trtcls5),
                 mace = pseudo_mace %>% select(nct_id, sex, age, arm_lvl, trtcls5),
                 .id = "outcome")
rm(tot)

### Pull in aggregate age Hba1c and MACE data
mace_agg <- mace_agg %>% 
  select(nct_id, age_m, age_sd = age_s, n = participants,
         trtcls5,
         arm_lvl = arm_lvl)
maceinhba1c_agg <- intersect(mace_agg$nct_id, agg$nct_id)
maceinhba1c <- union(maceinhba1c_agg, maceinhba1c_ipd)
agg <- agg %>% 
  anti_join(mace_agg %>% select(nct_id))
bth <- bind_rows(hba1c = agg[,names(mace_agg)],
                 mace = mace_agg,
                 .id = "outcome") 
bth <- bth %>% 
  filter(!nct_id %in% ipd$nct_id)
ipd_cls <- ipd %>% 
  select(nct_id, trtcls5) %>% 
  distinct()
ipd <- ipd %>% 
  select(-trtcls5)
ipd_agg <- ipd %>% 
  group_by(nct_id, outcome, arm_lvl) %>% 
  summarise(age_m = mean(age),
            age_sd = sd(age),
            n = length(age)) %>% 
  ungroup()
bth_cls <- bth %>% 
  distinct(nct_id, trtcls5)
bth <- bth %>% 
  select(-trtcls5)
bth <- bind_rows(yes = ipd_agg, 
                no = bth, .id = "for_ipd_chk")
ipd_agg$age_rnorm <- pmap(list(ipd_agg$n, ipd_agg$age_m, ipd_agg$age_sd), function(n, m, s) rnorm(n, m, s))
ipd_agg <- ipd_agg %>% 
  unnest(age_rnorm)
elig <- elig %>% 
  semi_join(bth)
## 56 trials without age eligibility data. All are non NCT
bth %>% 
  anti_join(elig) %>% 
  distinct(nct_id)
bth <- bth %>% 
  group_by(for_ipd_chk, outcome, nct_id) %>% 
  summarise(age_sd2 = CombSdVectorised(n = n, m = age_m, s = age_sd),
            age_m2 = weighted.mean(age_m, n),
            n = sum(n)) %>% 
  ungroup() %>% 
  rename(age_m = age_m2, age_sd = age_sd2)

bth <- bth %>% 
  left_join(elig) %>% 
  mutate(max_imp = if_else(is.na(max_age), "imp", "known"),
         max_age = if_else(is.na(as.integer(max_age)), 150L, as.integer(max_age)),
         min_age = if_else(is.na(as.integer(min_age)), 10L, as.integer(min_age)))

## Obtain mu and dispersion parameter ----
# Functions for truncted normal distributions
# https://www.r-bloggers.com/truncated-normal-distribution/

## Find mu and s from truncated normal distributions
bth <- bth %>% 
  mutate(bythis = paste0(for_ipd_chk, "__", nct_id))
siml_estimates <- by(bth, bth$bythis, EstimateMuDispAge)

## Extract and convert to list
siml_estimates <- map(siml_estimates, identity)
siml_estimates <- map(siml_estimates, ~ .x %>% mutate(pick = seq_along(trial_mean)))
siml_estimates <- bind_rows(siml_estimates, .id = "bythis")
siml_estimates <- siml_estimates %>% 
  filter(pick ==1)
siml_estimates <- siml_estimates %>% 
  separate(bythis, into = c("for_ipd_chk", "nct_id"), sep = "__")
bth <- bth %>% 
  inner_join(siml_estimates %>% 
               select(for_ipd_chk, nct_id, mu_x, sd_x))
## Sample from truncated normal distributions
bth$sim <- pmap(list(bth$n,
                     bth$min_age,
                     bth$max_age,
                     bth$mu_x,
                     bth$sd_x),
                function(n, a, b, m, s) {
                  x <- truncnorm::rtruncnorm(n = n,
                                        a = a, b = b, mean = m, sd = s)
                  x
                  })
bth <- bth %>% 
  unnest(sim)
names(ipd)
names(bth)
 
ages <- bind_rows(agg = bth %>% 
                          select(outcome, nct_id, 
                                 age = sim, 
                                 for_ipd_chk, max_imp),
                  ipd = ipd %>% 
                          select(outcome, nct_id, age) %>% 
                    mutate(for_ipd_chk = "both"),
                  .id = "data_lvl")
rm(ipd, bth, elig, mace_agg, mace_arms, pseudo_hba1c, pseudo_mace, agg, siml_estimates)
ipd_chk <- ages %>% 
  filter(for_ipd_chk %in% c("both", "yes"))
ages <- ages %>% 
  filter(for_ipd_chk %in% c("both", "no"))

## plot age distributions where have mean, sd and (for Hba1c trials) truncations ----
theme_minimal2 <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                            base_rect_size = base_size/22) 
{
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = 'bottom',
          complete = TRUE)
}

## Plot age distribution for ipd
ages <- ages %>% 
  mutate(category = 
           case_when(data_lvl == "agg" & outcome == "mace" ~ "MACE, aggregate",
                     data_lvl == "agg" & outcome == "hba1c" ~ "HbA1c, aggregate (random sample)",
                     data_lvl == "ipd" & outcome == "mace" ~ "MACE, IPD",
                     data_lvl == "ipd" & outcome == "hba1c" ~ "HbA1c, IPD"))

cls <- bind_rows(ipd_cls,
                 bth_cls) %>% 
  distinct(nct_id, trtcls5) %>% 
  filter(trtcls5 %in% c("A10BH", "A10BJ", "A10BK"))
ages2 <- ages %>% 
  left_join(cls, relationship =
              "many-to-many")
who <- who %>% 
  filter(`ATC code` %in% ages2$trtcls5)
who_vct <- who$`ATC level name`
names(who_vct) <- who$`ATC code`
ages2 <- ages2 %>% 
  mutate(trl_lbl = who_vct[trtcls5])
agg_trials <- ages2 %>% 
  filter(data_lvl == "agg", outcome == "hba1c", max_imp == "known") %>%
  distinct(trl_lbl, trtcls5, nct_id) %>% 
  group_by(trl_lbl, trtcls5) %>% 
  sample_frac(0.3) %>% 
  ungroup() %>% 
  pull(nct_id)

age_plot <- ggplot(ages2 %>% 
                     filter(nct_id %in% agg_trials |
                              data_lvl == "ipd" |
                              outcome == "mace"),
                       aes(x = interaction(nct_id, category), y = age, 
                           colour = category)) +
  geom_violin(draw_quantiles = c(0.025,0.975), scale = "width") +
  scale_x_discrete("Trials") +
  scale_y_continuous("Age (years) density plots with 2.5th and 97.5th centiles", 
                     breaks = seq(20, 100, 20)) +
  theme_minimal2() +
  facet_wrap(~trl_lbl, ncol = 1, scales = "free_x") +
  scale_color_discrete("")  + 
  ggtitle("Age distribution by trial") +
  geom_hline(yintercept = c(40, 80), linetype = "dashed", colour = "grey")
age_plot

ipd_chk <- bind_rows(ipd_chk %>% 
                        select(data_lvl, outcome, nct_id, age),
                      ipd_agg %>% mutate(data_lvl = "normal") %>% 
                        select(data_lvl, outcome, nct_id, age = age_rnorm)) 

## Drop trial with categorical age
# NCT02065791
ipd_chk <- ipd_chk %>% 
  mutate(distrb = case_when(
    data_lvl == "normal" ~ "normal",
    data_lvl == "agg" ~ "truncated normal",
    data_lvl == "ipd" ~ "empirical"),
    distrb = factor(distrb,
                    levels = c("normal", "truncated normal", "empirical")))
chk_sim <- ggplot(ipd_chk %>% 
                    filter(!nct_id == "NCT02065791"), aes(x = distrb, y = age, fill = distrb)) +
  geom_violin() +
  facet_wrap(~nct_id) +
  scale_x_discrete(guide = "none") +
  ggtitle("For IPD compare normal, truncated normal and empirical distributions")
chk_sim
pdf("Outputs/age_plots.pdf", height = 10, width = 20)
age_plot
chk_sim
dev.off()

## Summarise ages by trial class and ipd or agg
ages3 <- ages2 %>% 
  filter(nct_id %in% maceinhba1c) %>% 
  mutate(outcome = "hba1c",
         trl_lbl = str_replace(trl_lbl, "MACE", "HbA1c"))
ages3 <- ages2 %>% 
  bind_rows(ages3)
## Check on coding as quite complex against CTG website 
## n's all match, m and sd all match
ages_smry_chk <- ages2 %>% 
  bind_rows(ages3) %>% 
  group_by(nct_id) %>% 
  summarise(n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() %>%
  filter(str_detect(nct_id, "^N")) %>% sample_n(3)

ages_smry <- ages3 %>% 
  group_by(outcome, trl_lbl, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup()
ages3_trl_grp <- ages3 %>% 
  distinct(nct_id, trl_lbl) %>% 
  distinct(nct_id, .keep_all = TRUE)
ages4 <- ages3 %>% 
  semi_join(ages3_trl_grp) %>% 
  select(-trl_lbl)
ages_smry_tot <- ages4 %>% 
  group_by(outcome, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() %>% 
  mutate(trl_lbl = "Total trials")
ages_smry_final <- bind_rows(ages_smry,
                       ages_smry_tot) %>% 
  arrange(outcome, trl_lbl, data_lvl)
write_csv(ages_smry_final, "Outputs/age_summary.csv")
