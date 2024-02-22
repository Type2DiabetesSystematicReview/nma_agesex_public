library(tidyverse)
library(truncnorm)

## read and transform data ----
who <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx") 
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
mace_agg <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_agg
eliga <- 
  readRDS("../extract_transform/aact/data/April2023_extract/aact_extract_April_2023.Rds")$eligibilities %>% 
  mutate(across(c(minimum_age, maximum_age), ~ if_else(.x == "N/A", "NA Years", .x))) %>% 
  separate(minimum_age, into = c("min_age", "min_age_unit"), sep = "\\s") %>% 
  separate(maximum_age, into = c("max_age", "max_age_unit"), sep = "\\s") %>% 
  select(nct_id, min_age, max_age, min_age_unit, max_age_unit)
eligb <- read_csv("Data/elig_cri_all_tr.csv") %>% 
  select(nct_id = trial_id, 
         min_age = less_than_age_excluded, 
         min_age_unit = minimum_age_unit, 
         max_age = more_than_age_excluded, 
         max_age_unit = maximum_age_unit) %>% 
  mutate(across(c(min_age, max_age), as.character))
elig <- bind_rows(eliga,
                  eligb) %>% 
  distinct(nct_id, .keep_all = TRUE)
rm(eliga, eligb)
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
pseudo_hba1c <- pseudo_hba1c %>% 
  anti_join(pseudo_mace %>% select(nct_id))
ipd <- bind_rows(hba1c = pseudo_hba1c %>% select(nct_id, sex, age, arm_lvl, trtcls5),
                 mace = pseudo_mace %>% select(nct_id, sex, age, arm_lvl, trtcls5),
                 .id = "outcome")
ipd <- ipd %>% 
  group_by(nct_id) %>% 
  mutate(trtcls5_trial = 
           case_when(
             any(trtcls5 == "A10BK") ~ "A10BK",
             any(trtcls5 == "A10BJ") ~ "A10BJ",
             any(trtcls5 == "A10BH") ~ "A10BH",
             TRUE ~ "Other"
               )) %>% 
  ungroup()
rm(tot)

### Pull in aggregate age Hba1c and MACE data
mace_agg <- mace_agg %>% 
  select(nct_id, age_m, age_sd = age_s, n = participants,
         trtcls5,
         arm_lvl = arm_lvl)
agg <- agg %>% 
  anti_join(mace_agg %>% select(nct_id))

bth <- bind_rows(hba1c = agg[,names(mace_agg)],
                 mace = mace_agg,
                 .id = "outcome") 
bth <- bth %>% 
  filter(!nct_id %in% ipd$nct_id)
ipd_agg <- ipd %>% 
  group_by(nct_id, outcome, arm_lvl, trtcls5) %>% 
  summarise(age_m = mean(age),
            age_sd = sd(age),
            n = length(age)) %>% 
  ungroup()
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
            n = sum(n),
            trtcls5_trial = 
               case_when(
                   any(trtcls5 == "A10BK") ~ "A10BK",
                   any(trtcls5 == "A10BJ") ~ "A10BJ",
                   any(trtcls5 == "A10BH") ~ "A10BH",
                   TRUE ~ "Other")) %>% 
  ungroup() %>% 
  rename(age_m = age_m2, age_sd = age_sd2)
# Note includes 5 more trials with nonoe of the three novel anti-diabetics. Leave in in case helps connect network
elig <- elig %>% 
  mutate(max_age = if_else(is.na(as.integer(max_age)), 150L, as.integer(max_age)),
         min_age = if_else(is.na(as.integer(min_age)), 10L, as.integer(min_age)))
bth <- bth %>% 
  left_join(elig) %>% 
  mutate(max_age = if_else(is.na(as.integer(max_age)), 150L, as.integer(max_age)),
         min_age = if_else(is.na(as.integer(min_age)), 10L, as.integer(min_age)))

## Obtain mu and dispersion parameter ----
# Functions for truncted normal distributions
# https://www.r-bloggers.com/truncated-normal-distribution/
mean.tnorm<-function(mu,sd,lower,upper){
  ##return the expectation of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  mean=mu+sd*(dnorm(lower.std)-dnorm(upper.std))/
    (pnorm(upper.std)-pnorm(lower.std))
  return(mean)
}
var.tnorm<-function(mu,sd,lower,upper){
  ##return the variance of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  variance=sd^2*(1+(lower.std*dnorm(lower.std)-upper.std*dnorm(upper.std))/
                   (pnorm(upper.std)-pnorm(lower.std))-((dnorm(lower.std)-dnorm(upper.std))/
                                                          (pnorm(upper.std)-pnorm(lower.std)))^2)
  return(variance)
}

## Find mu and s from truncated normal distributions
bth <- bth %>% 
  mutate(bythis = paste0(for_ipd_chk, "__", nct_id))
siml_estimates <- by(bth, bth$bythis, function (x) {
  # Loop through trials
  if(str_sub(x$nct_id, -1) == 1) print(x$nct_id)
  
  lower <- x$min_age
  upper <- x$max_age
  trial_mean <- x$age_m
  trial_sd <- x$age_sd
  trial_var <- x$age_sd^2
  
  ## Create grid
  mu_x <- seq(lower, upper, 0.5)
  sd_x <- seq(1, upper-lower, 0.5)
  full_grid <- expand.grid(mu_x = mu_x, sd_x = sd_x)
  
  ## Calculate for all values of grid, is vectorised so is fast, is faster than one in truncnorm package
  full_grid$mean_x <- mean.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
  full_grid$var_x <- var.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
  
  # print(nrow(full_grid))
  # browser()  
  ## Identify closest values
  full_grid <- full_grid %>%
    as_tibble() %>%
    mutate(mean_diff = abs(trial_mean - mean_x),
           var_diff = abs(trial_var - var_x),
           total_diff = mean_diff + var_diff) %>%
    arrange(total_diff, mean_diff, var_diff)
  ## Append original parameters
  estimate <- full_grid %>%
    slice(1:10) %>%
    mutate(trial_mean = trial_mean,
           trial_var = trial_var,
           trial_lower = lower,
           trial_upper = upper,
           trial_sd = trial_sd) %>%
    select(trial_mean, mean_x, trial_var, var_x, mu_x, sd_x, trial_sd, everything())
  estimate 
})

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
               select(for_ipd_chk, nct_id, mu_x, var_x))
## Sample from truncated normal distributions
# bth$sim <- pmap(list(bth$min_age, 
#                      bth$max_age,
#                      bth$mu_x,
#                      bth$var_x), 
#                 function(a, b, m, v) {
#                   x <- truncnorm::qtruncnorm(p = c(0.025, 0.975),
#                                         a = a, b = b, mean = m, sd = v^0.5)
#                   tibble(lci = x[1], uci = x[2])
#                   })

bth$sim <- pmap(list(bth$n,
                     bth$min_age,
                     bth$max_age,
                     bth$mu_x,
                     bth$var_x),
                function(n, a, b, m, v) {
                  x <- truncnorm::rtruncnorm(n = n,
                                        a = a, b = b, mean = m, sd = v^0.5)
                  x
                  })

bth <- bth %>% 
  unnest(sim)

names(ipd)
names(bth)
 
ages <- bind_rows(agg = bth %>% 
                          select(outcome, nct_id, 
                                 age = sim, 
                                 trtcls5_trial, for_ipd_chk),
                  ipd = ipd %>% 
                          select(outcome, nct_id, age,
                                 trtcls5_trial) %>% 
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
agg_trials <- ages %>% 
  filter(data_lvl == "agg", outcome == "hba1c") %>%
  pull(nct_id) %>% 
  unique()
agg_trials_smpl <- sample(agg_trials %>% unique(), size = 0.3*length(agg_trials))
ages <- ages %>% 
  mutate(category = 
           case_when(data_lvl == "agg" & outcome == "mace" ~ "MACE, aggregate",
                     data_lvl == "agg" & outcome == "hba1c" ~ "HbA1c, aggregate (random sample)",
                     data_lvl == "ipd" & outcome == "mace" ~ "MACE, IPD",
                     data_lvl == "ipd" & outcome == "hba1c" ~ "HbA1c, IPD"))
who <- who %>% 
  filter(`ATC code` %in% ages$trtcls5_trial)
who_vct <- who$`ATC level name`
names(who_vct) <- who$`ATC code`

ages <- ages %>% 
  mutate(trl_lbl = who_vct[trtcls5_trial])

age_plot <- ggplot(ages %>% 
                     filter(!trtcls5_trial == "Other",
                            nct_id %in% agg_trials_smpl |
                              data_lvl == "ipd" |
                              outcome == "mace"),
                       aes(x = interaction(nct_id, category), y = age, 
                           colour = category)) +
  geom_violin(draw_quantiles = c(0.025,0.975), scale = "width") +
  scale_x_discrete("Trials") +
  scale_y_continuous("Age (years) density plots with 2.5th and 97.5th centiles", breaks = seq(20, 100, 10)) +
  theme_minimal2() +
  facet_wrap(~trl_lbl, ncol = 1, scales = "free_x") +
  scale_color_discrete("")  + 
  ggtitle("Age distribution by trial") +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey")
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
