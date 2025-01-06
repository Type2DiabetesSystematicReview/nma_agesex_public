library(tidyverse)
library(truncnorm)
source("Scripts/common_functions/Scripts/truncated_normal.R")
source("Scripts/00_functions.R")
## read and transform data ----
who <- read_csv("Data/whoatcdiabetesnodose.csv")
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
elig <- read_csv("Data/cleaned_data/Data/age_max_min_elig.csv")
source("Scripts/common_functions/Scripts/combine_sd.R")

## Set exclusion to hba1c for ones with dispersion /percentage issues
exclusions <- read_csv("Data/exclusions_update.csv")
exclusions <- exclusions %>% 
  mutate(across(c(agg_hba1c, any_hba1c), ~ if_else(exclusion_reason %in% c(
    "only reported outcomes as medians or percentage change",
    "no dispersion statistics reported"), 1L, .x)))
write_csv(exclusions, "Data/exclusions_update.csv")

exclusions <- read_csv("Data/exclusions_update.csv")
keephba1c <- exclusions %>% 
  filter(any_hba1c ==1L | any_mace ==1L, exclude ==0L)

hba1c_agg <- tot %>% 
  select(drug_regime_smpl, agg) %>% 
  unnest(agg)
pseudo_hba1c <- tot %>% 
  select(drug_regime_smpl, ipd) %>% 
  unnest(ipd) %>% 
  mutate(sex = if_else(sex == 1L, "M", "F"))
ipd_ecdf <- bind_rows(hba1c = pseudo_hba1c %>% select(nct_id, sex, age, arm_lvl, trtcls5, arm_lvl),
                 .id = "outcome")
rm(tot)

### Pull in aggregate age Hba1c and MACE data
bth <- bind_rows(hba1c = hba1c_agg,
                 .id = "outcome") 
ipd_agg <- ipd_ecdf %>% 
  group_by(nct_id, outcome, arm_lvl, trtcls5) %>% 
  summarise(age_m = mean(age),
            age_sd = sd(age),
            n = length(age)) %>% 
  ungroup()
bth <- bind_rows(ipd = ipd_agg, 
                 agg = bth[, names(ipd_agg)], 
                 .id = "orig_data_lvl")
elig <- elig %>% 
  semi_join(bth)
## 55 trials without age eligibility data. All are non NCT
bth %>% 
  anti_join(elig) %>% 
  distinct(nct_id)
bth <- bth %>% 
  left_join(elig) %>% 
  mutate(max_imp = if_else(is.na(max_age), "imp", "known"),
         max_age = if_else(is.na(as.integer(max_age)), 150L, as.integer(max_age)),
         min_age = if_else(is.na(as.integer(min_age)), 10L, as.integer(min_age)))
bth <- bth %>% 
  semi_join(keephba1c %>% rename(nct_id = trial_id))

## Obtain mu and dispersion parameter ----
# Functions for truncted normal distributions
# https://www.r-bloggers.com/truncated-normal-distribution/

## Find mu and s from truncated normal distributions
# siml_estimates <- by(bth, bth$nct_id, EstimateMuDispAge)
SafeMuSigmaCalc <- safely(MuSigmaCalc)
bth$res <- pmap( list(bth$min_age, bth$max_age, bth$age_m, bth$age_sd),
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
if(any(map_lgl(bth$res, is_null))) warning ("Cycled through 3 algorithms without running succesfully")
bth <- bth %>% 
  unnest(res)
if(any(!bth$convergence %in% c(0, 999))) warning ("Optim did not converge for at least one row")
## Sample from truncated normal distributions
bth$sim <- pmap(list(bth$n,
                     bth$min_age,
                     bth$max_age,
                     bth$mu,
                     bth$sigma),
                function(n, a, b, mu, sigma) {
                  x <- truncnorm::rtruncnorm(n = n,
                                        a = a, b = b, mean = mu, sd = sigma)
                  x
                  })
bth <- bth %>% 
  unnest(sim)

## For plotting ages, select ecdf from ipd (npte fpr MACE IPD is actually normal distribution, but checked this with plots inside vivli)
ages <- bind_rows(agg = bth %>% 
                    filter(orig_data_lvl == "agg") %>% 
                          select(outcome,
                                 nct_id, 
                                 age = sim, 
                                 max_imp,
                                 arm_lvl,
                                 trtcls5),
                  ipd = ipd_ecdf %>% 
                    select(outcome,
                           nct_id, 
                           age = age, 
                           # max_imp,
                           arm_lvl,
                           trtcls5),
                  .id = "data_lvl")

## Using IPD compare normal, truncnormal and ECDF
## sample from normal distribution (wihtout truncation) for later comparison
ipd_agg$age_rnorm <- pmap(list(ipd_agg$n, ipd_agg$age_m, ipd_agg$age_sd), function(n, m, s) rnorm(n, m, s))
ipd_agg <- ipd_agg %>% 
  unnest(age_rnorm)

ipd_chk <- bind_rows(truncnorm = bth %>% 
                       filter(orig_data_lvl == "ipd") %>% 
                       select(outcome,
                              nct_id, 
                              age = sim, 
                              max_imp,
                              arm_lvl,
                              trtcls5),
                     ecdf = ipd_ecdf %>% 
                       select(outcome,
                              nct_id, 
                              age = age, 
                              # max_imp,
                              arm_lvl,
                              trtcls5),
                     normal = ipd_agg %>% 
                       select(outcome,
                              nct_id, 
                              age = age_rnorm, 
                              # max_imp,
                              arm_lvl,
                              trtcls5),
                     .id = "age_method")
rm(bth, elig, hba1c_agg, ipd_agg, ipd_ecdf, pseudo_hba1c)

## plot age distributions where have mean, sd and (for Hba1c trials) truncations ----


## Replicate trials to allow for fact that some trials are in multiple classes
ages <- ages %>% 
  mutate(category = 
           case_when(data_lvl == "agg" & outcome == "mace" ~ "MACE, aggregate",
                     data_lvl == "agg" & outcome == "hba1c" ~ "HbA1c, aggregate (random sample)",
                     data_lvl == "ipd" & outcome == "mace" ~ "MACE, IPD",
                     data_lvl == "ipd" & outcome == "hba1c" ~ "HbA1c, IPD"))
write_csv(ages %>% count(data_lvl, nct_id), "Outputs/nt_trl_lvl_hba1c.csv")
novel <- c("A10BK", "A10BH", "A10BJ")
ages_cls <- map(novel, ~ {
  ages %>% 
    group_by(nct_id) %>% 
    mutate(cls = if_else(any(trtcls5 == .x),
                         .x,
                         ""))  %>% 
    ungroup()
}) 
ages_cls <- bind_rows(ages_cls) %>% 
  filter(!cls == "")
who <- who %>% 
  filter(`ATC code` %in% ages_cls$cls)
who_vct <- who$`ATC level name`
names(who_vct) <- who$`ATC code`
ages_cls <- ages_cls %>% 
  mutate(trl_lbl = who_vct[cls])
## Take random sample of aggregte level trials for plot
agg_trials <- ages_cls %>% 
  filter(data_lvl == "agg", outcome == "hba1c", max_imp == "known") %>%
  distinct(trl_lbl, trtcls5, nct_id) %>% 
  group_by(trl_lbl, trtcls5) %>% 
  sample_frac(0.3) %>% 
  ungroup() %>% 
  pull(nct_id)
saveRDS(ages_cls %>% 
          filter(nct_id %in% agg_trials |
                   data_lvl == "ipd" |
                   outcome == "mace"),
        "Scratch_data/for_age_violin_hba1c_add_mace.Rds")

## Drop trial with categorical age
# NCT02065791
ipd_chk <- ipd_chk %>% 
  mutate(distrb = case_when(
    age_method == "normal" ~ "normal",
    age_method == "truncnorm" ~ "truncated normal",
    age_method == "ecdf" ~ "empirical"),
    distrb = factor(distrb,
                    levels = c("normal", "truncated normal", "empirical")))
chk_sim <- ggplot(ipd_chk %>% 
                    filter(!nct_id == "NCT02065791"), aes(x = distrb, y = age, fill = distrb)) +
  geom_violin() +
  facet_wrap(~nct_id) +
  scale_x_discrete(guide = "none") 
saveRDS(chk_sim, "Scratch_data/chk_sim.Rds")
pdf("Outputs/age_plots_check_simulation.pdf", height = 10, width = 20)
chk_sim+
  ggtitle("For IPD compare normal, truncated normal and empirical distributions")
dev.off()

chk_sim2 <- ggplot(ipd_chk %>% 
                    filter(!nct_id == "NCT02065791", !distrb == "normal"), aes(x = age, fill = distrb)) +
  geom_density(alpha = 0.2) +
  facet_wrap(~nct_id, scales = "free_y") 
ggsave(plot = chk_sim2, filename = "Outputs/age_plots_check_simulation_density_plots.pdf", height = 20, width = 20)

## Summarise ages by trial class and ipd or agg
ages_cls2 <- bind_rows(ages_cls,
                       ages %>% mutate(cls = "Any",
                                       trl_lbl = "Total trials")) 
## Generate age summaries
## n's all match, m and sd all match
ages_smry_chk <- ages_cls2 %>% 
  group_by(nct_id) %>% 
  summarise(n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup() %>%
  filter(str_detect(nct_id, "^N")) %>% sample_n(3)

ages_smry <- ages_cls2 %>% 
  group_by(cls, trl_lbl, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup()
fi_trials <- read_lines("Data/fi_trials.txt")
ages_smry_fi <- ages_cls2 %>% 
  mutate(data_lvl = if_else(nct_id %in% fi_trials, "ipd_fi", data_lvl)) %>% 
  group_by(cls, trl_lbl, data_lvl) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            n = length(nct_id),
            m = mean(age),
            s = sd(age),
            q05 = quantile(age, probs = 0.05),
            q95 = quantile(age, probs = 0.95)) %>% 
  ungroup()
write_csv(ages_smry, "Outputs/age_summary_hba1c.csv")
write_csv(ages_smry_fi, "Outputs/age_summary_hba1c_fi.csv")

## examine with missing ages - 10 (9 IPD and one aggregate) of the 12 are mace trials so adding later
rv_msng <- exclusions %>% 
  filter(exclude ==0) %>% 
  filter(!trial_id %in% ages_cls2$nct_id)
rv_msng %>% 
  count(any_mace, ipd_mace)

## review remainder, two are our "add not sure why left outs
basedata <- readRDS("Scratch_data/agg_hba1c_base.Rds")$age
rv_msng_nomace <- rv_msng %>% 
  filter(!any_mace ==1L)
rv_msng_nomace %>% 
  left_join(basedata %>% distinct(trial_id, .keep_all = TRUE))

## examine with got ages but in exclusions; none
exclusions %>% 
  filter(exclude == 1L) %>% 
  filter(trial_id %in% ages_cls2$nct_id)
