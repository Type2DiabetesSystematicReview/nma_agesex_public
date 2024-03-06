# summarise age for mace separately
library(tidyverse)
library(truncnorm)
source("../common_functions/Scripts/truncated_normal.R")

## read and transform data ----
who <- readxl::read_excel("~/2018 ATC index with DDDs.xlsx") 
elig <- read_csv("../cleaned_data/Data/age_max_min_elig.csv")
source("../common_functions/Scripts/combine_sd.R")

mace_agg <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_agg
pseudo_mace <- readRDS("Scratch_data/ipd_age_sex_mace.Rds")
mace_arms <- readRDS("Scratch_data/mace_arms_agg_data.Rds")$mace_arms
pseudo_mace <- pseudo_mace %>% 
  mutate(ipd_arm = str_to_lower(arm)) %>% 
  inner_join(mace_arms %>% select(nct_id, ipd_arm, arm_lvl, trtcls5))
mace <- pseudo_mace %>% select(nct_id, sex, age, arm_lvl, trtcls5, arm_lvl)
rm(pseudo_mace)
mace_agg <- mace_agg %>% 
  select(nct_id, age_m, age_sd = age_s, n = participants,
         trtcls5,
         arm_lvl = arm_lvl)

## one trial without mace eligibility
mace_agg <- mace_agg %>% 
  left_join(elig) %>% 
  mutate(max_imp = if_else(is.na(max_age), "imp", "known"),
         max_age = if_else(is.na(as.integer(max_age)), 150L, as.integer(max_age)),
         min_age = if_else(is.na(as.integer(min_age)), 10L, as.integer(min_age)))
  
## sample from truncated normal distribution for mace
SafeMuSigmaCalc <- safely(MuSigmaCalc)
mace_agg$res <- pmap( list(mace_agg$min_age, mace_agg$max_age, mace_agg$age_m, mace_agg$age_sd),
                 function(min_age, max_age, age_m, age_sd) {
                   ## If no age limits return mean and SD
                   if(max_age >= 100 & min_age <= 20) {
                     b <- tibble(mu = age_m, 
                                 sigma = age_sd, 
                                 convergence = 999, 
                                 value = 0,
                                 method = "ignore_limits")
                     return(b)  
                   } 
                   ## If age limits first try L-BFGS note this sometimes fails so need to run safely 
                   a <- SafeMuSigmaCalc(min_age, max_age, age_m, age_sd, mymethod = "L-BFGS-B")
                   if(is.null(a$error)) {
                     ## double if statement here as this object is only present if it runs
                     if(a$result$convergence[[1]] == 0) {
                       a <- a$result
                       b <- tibble(mu = a$par[1],
                                   sigma = a$par[2],
                                   convergence = a$convergence[[1]],
                                   value = a$value[[1]],
                                   method = "L-BFGS-B")
                       return(b) 
                     } ## end condition only if converges
                   } ## end condition only if not null
                   ## If L-BFGS-B fails or does not converge run Nelder-Mead Nelder-Mead
                   a <- SafeMuSigmaCalc(min_age, max_age, age_m, age_sd, mymethod = "CG")
                   if(is.null(a$error)) {
                     a <- a$result
                     b <- tibble(mu = a$par[1],
                                 sigma = a$par[2],
                                 convergence = a$convergence[[1]],
                                 value = a$value[[1]],
                                 method = "CG")
                     return(b) 
                   } 
                 })
## all converge (as with HbA1c)
if(any(map_lgl(mace_agg$res, is_null))) warning ("Cycled through 3 algorithms without running succesfully")
mace_agg <- mace_agg %>% 
  unnest(res)
if(any(!mace_agg$convergence %in% c(0, 999))) warning ("Optim did not converge for at least one row")

## sample from distribution
## Sample from truncated normal distributions
mace_agg$sim <- pmap(list(mace_agg$n,
                          mace_agg$min_age,
                          mace_agg$max_age,
                          mace_agg$mu,
                          mace_agg$sigma),
                function(n, a, b, mu, sigma) {
                  x <- truncnorm::rtruncnorm(n = n,
                                             a = a, b = b, mean = mu, sd = sigma)
                  x
                })
mace_agg <- mace_agg %>% 
  unnest(sim)
mace <- bind_rows(agg = mace_agg %>% 
                    select(nct_id, 
                           age = sim, 
                           max_imp,
                           arm_lvl,
                           trtcls5),
                  ipd = mace %>% 
                    select(nct_id, 
                           age = age, 
                           # max_imp,
                           arm_lvl,
                           trtcls5),
                  .id = "data_lvl")
rm(mace_agg, mace_arms, elig)
## note that MACE novel antidiabetic trials are all against non novel antidiabetic comparators
## so easier to code (all mutually exclusive)
novel <- c("A10BK", "A10BH", "A10BJ")
mace <- mace %>% 
  group_by(nct_id) %>% 
  mutate(trtcls5 = case_when(
    any(trtcls5 %in% novel[1]) ~ novel[1],
    any(trtcls5 %in% novel[2]) ~ novel[2],
    any(trtcls5 %in% novel[3]) ~ novel[3],
  )) %>% 
  ungroup()
mace_smry1 <- mace %>% 
  group_by(trtcls5, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup()
mace_smry2 <- mace %>% 
  group_by(data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() %>% 
  mutate(trtcls5 = "Any")
mace_smry_trial <- mace %>% 
  group_by(nct_id) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) 
write_csv(mace_smry_trial %>% select(-trials), "Outputs/age_summary_trials_mace.csv")
hba1c_smry <- read_csv("Outputs/age_summary_hba1c.csv")
mace_smry <- bind_rows(mace_smry1,
                       mace_smry2) %>% 
  rename(cls = trtcls5) %>% 
  left_join(hba1c_smry %>% distinct(cls, trl_lbl))
mace_smry <- mace_smry[, names(hba1c_smry)]
write_csv(mace_smry, "Outputs/age_summary_mace.csv")
rm(hba1c_smry, mace_smry, mace_smry1, mace_smry2)

## read in hba1c ages for plotting
hba1c <- readRDS("Scratch_data/for_age_violin_hba1c_add_mace.Rds")

theme_minimal2 <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                            base_rect_size = base_size/22) {
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
mace <- mace %>% 
  mutate(outcome = "mace") %>% 
  rename(cls = trtcls5)
cls_lbl <- hba1c %>% 
  distinct(cls, trl_lbl)
mace <- mace %>% 
  inner_join(cls_lbl)
hba1c <- hba1c[, names(mace)]
hba1c <- hba1c %>% 
  filter(!nct_id %in% mace$nct_id)
## Drop mace trials from hba1c set to prevent double counting in plot
forplot <- bind_rows(hba1c, 
                     mace)
forplot <- forplot  %>% 
  mutate(category = 
           case_when(data_lvl == "agg" & outcome == "mace" ~ "MACE, aggregate",
                     data_lvl == "agg" & outcome == "hba1c" ~ "HbA1c, aggregate (random sample)",
                     data_lvl == "ipd" & outcome == "mace" ~ "MACE, IPD",
                     data_lvl == "ipd" & outcome == "hba1c" ~ "HbA1c, IPD"))
age_plot <- ggplot(forplot,
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
saveRDS(age_plot, "Scratch_data/age_plots.Rds")
pdf("Outputs/age_plots.pdf", height = 10, width = 20)
age_plot
dev.off()