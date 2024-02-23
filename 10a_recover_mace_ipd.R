#13
library(tidyverse)
source("Scripts/00_functions.R")
source("../common_functions/Scripts/fractional_polynomial.R")
## event times separate for single trial and other 4. Censoring times brought together early because sample with sd = 0 for categorical age
## trials on new vivli repository - 8697
fls <- c(
  # "event_time_distribution_linear.csv",
  "event_time_distribution_fp.csv",
  "event_time_distribution_single_trial.csv",
  "censoring_distribution.csv",
  "censoring_distribution_single_trial.csv")
res <- map(fls, ~ read_csv(paste0("../from_vivli/Data/agesexmace_8697/", .x)))
names(res) <- str_sub(fls, 1, -5)
list2env(res, envir = .GlobalEnv)
rm(res)
## trials on new previous repository - 6115
censoring_distribution_fp_6115 <- read_csv("../from_vivli/Data/agesexmace_6115/censoring_distribution_fp.csv")
event_time_distribution_fp_6115 <- read_csv("../from_vivli/Data/agesexmace_6115/event_time_distribution_fp.csv")
censoring_distribution_fp_6115 <- censoring_distribution_fp_6115 %>% 
  rename(arm = arm_label,
         sex = sex_decoded) %>% 
  mutate(sex = if_else(sex == "female", "F", "M"),
         `Number of quantiles` = str_count(`Censored at quantiles (%)`, "\\="))
censoring_distribution <- bind_rows(censoring_distribution,
                                    censoring_distribution_fp_6115)
rm(censoring_distribution_fp_6115)
event_time_distribution_fp_6115 <- event_time_distribution_fp_6115 %>% 
  rename(arm = arm_label,
         sex = sex_decoded) %>% 
  mutate(sex = if_else(sex == "female", "F", "M")) %>% 
  select(-fu_m, -fu_t) %>% 
  mutate(est_age2 = 0,
         se_age2 = 0)
event_time_distribution_fp <- bind_rows(event_time_distribution_fp,
                                     event_time_distribution_fp_6115 %>% mutate(r = as.character(r)))
rm(event_time_distribution_fp_6115)
## note there was no need to export both linear and FP as gave same result

## sample from censoring distributions (not dependent on age - 5 trials) ----
ColnamesPipe <- function(x, vct){
  colnames(x) <- vct
  x
}
censoring_distribution <- censoring_distribution %>% 
  anti_join(censoring_distribution_single_trial %>% select(nct_id))
censoring_distribution$splt <-  str_split(censoring_distribution$`Censored at quantiles (%)`, pattern = "\\,")
censoring_distribution$splt <- map(censoring_distribution$splt, ~ str_split_fixed(.x, "\\=", n = 2) %>% ColnamesPipe(c("x", "y")) %>% 
                                     as_tibble() %>% 
                                     mutate(across(everything(), ~ str_trim(.x) %>% as.integer()),
                                            x = x/100))
censoring_distribution$cnsr_time <- map2(censoring_distribution$splt, censoring_distribution$participants, function(a, b) {
  fromdec <- approx(x = a$x, y =a$y, xout = seq(0, 1, 1/(b-1))) %>% 
    as_tibble() 
  fromdec$y
})
censoring_distribution$splt <- NULL

## single trial
censoring_distribution_single_trial$splt <-  str_split(censoring_distribution_single_trial$`Censored at quantiles (%)`, pattern = "\\,")
censoring_distribution_single_trial$splt <- map(censoring_distribution_single_trial$splt, ~ str_split_fixed(.x, "\\=", n = 2) %>% ColnamesPipe(c("x", "y")) %>% 
                                                  as_tibble() %>% 
                                                  mutate(across(everything(), ~ str_trim(.x) %>% as.integer()),
                                                         x = x/100))
censoring_distribution_single_trial$cnsr_time <- map2(censoring_distribution_single_trial$splt, censoring_distribution_single_trial$participants, function(a, b) {
  fromdec <- approx(x = a$x, y =a$y, xout = seq(0, 1, 1/(b-1))) %>% 
    as_tibble() 
  fromdec$y
})
censoring_distribution_single_trial$splt <- NULL

## note that can run from here as set age_sd to zero for the trial where age is solely categorical
censoring_distribution <- bind_rows(censoring_distribution,
                                    censoring_distribution_single_trial)
cnsr_max <- censoring_distribution %>% 
  select(nct_id, arm, sex, cnsr_time) %>% 
  unnest(cnsr_time) %>% 
  arrange(desc(cnsr_time)) %>% 
  distinct(nct_id, arm, sex, .keep_all = TRUE)

## Loop here at 3 seeds
seed <- c(1111)
  print(seed)
  set.seed(seed)
  ## note can do MVN here if want to use correlation
  censoring_distribution$age <- pmap(list(censoring_distribution$participants,
                                          censoring_distribution$age_m,
                                          censoring_distribution$age_s), function(n, m, s) rnorm(n, m, s))
  censoring_distribution <- censoring_distribution %>% 
    select(nct_id, arm, sex, participants, age, cnsr_time) %>% 
    unnest(c(age, cnsr_time))
  censoring_distribution %>% 
    count(nct_id, arm, sex, participants) %>% 
    filter(!participants == n)
  plot_ct <- ggplot(censoring_distribution, aes(x = interaction(arm, sex), y = cnsr_time/365)) +
    geom_violin() +
    facet_wrap(~nct_id, scales = "free")
  plot_ct
  rm(plot_ct)
  
  ## sample from event time distributions (4 trials) ----
  # event time distribution model with fractional polynomials 
  set.seed(seed)
  event_time_distribution_fp$age <-  pmap(list(event_time_distribution_fp$participants,
                                               event_time_distribution_fp$age_m,
                                               event_time_distribution_fp$age_s), function(n, m, s) rnorm(n, m, s))
  event_time_distribution_fp <- event_time_distribution_fp %>% 
    inner_join(cnsr_max)
  
  ## Only 2-years difference in mean age between women and men
  event_time_distribution_fp %>%
    select(nct_id, arm, sex, age_m) %>% 
    spread(sex, age_m) %>% 
    mutate(sex_diff = F-M) %>% 
    write_csv("Outputs/mean_age_difference_men_women.csv")
  ## simulate 
  # mydf <- mydf %>%
  #   mutate(age1 = MP1(age/mod_fp$scale, mod_fp$power1),
  #          age2 = MP2(age/mod_fp$scale, mod_fp$power1, mod_fp$power2))
  # time_sim <- rnorm(nrow(mydf), mod_fp$est_cept +
  #                     mod_fp$est_age1*mydf$age1 +
  #                     mod_fp$est_age2*mydf$age2,
  # mod_fp$residsd)
  
  event_time_distribution_fp$time <- pmap(list(
    age = event_time_distribution_fp$age, 
    scale = event_time_distribution_fp$scale,
    power1 = event_time_distribution_fp$power1,
    power2 = event_time_distribution_fp$power2,
    residsd = event_time_distribution_fp$residsd,
    est_cept = event_time_distribution_fp$est_cept,
    est_age1 = event_time_distribution_fp$est_age1,
    est_age2 = event_time_distribution_fp$est_age2,
    cnsr_time = event_time_distribution_fp$cnsr_time), 
    function(age, scale, power1, power2, residsd, est_cept, est_age1, est_age2, cnsr_time){
      age1 = MP1(age/scale, power1)
      age2 = MP2(age/scale, power1, power2)
      time_mt = est_cept +
        est_age1*age1 +
        est_age2*age2
      time_t = rnorm(length(age1), time_mt, residsd)
      time_p = plogis(time_t)
      time = cnsr_time*time_p
      time})
  event_time_distribution_fp <- event_time_distribution_fp %>% 
    unnest(c(age, time))
  
  plot_et <- ggplot(event_time_distribution_fp, aes(x = interaction(arm, sex), y = time/365)) +
    geom_violin() +
    facet_wrap(~nct_id, scales = "free")
  plot_et
  rm(plot_et)
  
  ## join together censored and uncensored pseudo-data for 4 trials ----
  pseudo <- bind_rows(censoring_distribution %>% 
                        select(nct_id, arm, sex, age, time = cnsr_time) %>% 
                        mutate(event = 0L,
                               event_method = "nonapplicable"),
                      event_time_distribution_fp %>% 
                        select(nct_id, arm, sex, age, time = time) %>% 
                        mutate(event = 1L,
                               event_method = "age_fp"))
  rm(censoring_distribution, event_time_distribution_fp)
  
  ## Add in event data for event data fir the single trial ----
  et_single <- event_time_distribution_single_trial
  
  et_single <- et_single %>% 
    inner_join(cnsr_max)
  et_single$time <- pmap(list(et_single$events,
                              et_single$time_t_m,
                              et_single$time_t_s), function(n, m, s) rnorm(n, m, s))
  et_single <- et_single %>% 
    mutate(event = 1L) %>% 
    select(nct_id, arm, sex, age, time, cnsr_time) %>% 
    unnest(cols = c(time))
  et_single <- et_single %>% 
    mutate(time = plogis(time)*cnsr_time) %>% 
    select(-cnsr_time)
  et_single <- et_single %>% 
    mutate(event = 1L,
           event_method = "single_trial")
  
  pseudo <- bind_rows(pseudo,
                      et_single)

  saveRDS(pseudo, "Scratch_data/ipd_age_sex_mace.Rds")
  pseudo_smry <- pseudo %>% 
    group_by(nct_id, arm) %>% 
    summarise(fu_mean = mean(time),
              fu_max = max(time), 
              fu_median = median(time)) %>% 
    ungroup()

  pseudo_smry <- pseudo_smry %>% 
    mutate(arm_simple = if_else(arm %in% c("placebo", "Placebo", "PLACEBO"), "placebo", "active")) %>% 
    arrange(nct_id, arm_simple) %>% 
    group_by(nct_id) %>%
    mutate(arm_simple = if_else(arm_simple != "placebo", 
                                paste0(arm_simple, 1+cumsum(duplicated(arm_simple))),
                                "placebo")) %>% 
    ungroup()
    
  pseudo_smry <- pseudo_smry %>% 
    group_by(arm_simple) %>% 
    nest() %>% 
    ungroup()
  pseudo_smry$r_u_med <- map_dbl(pseudo_smry$data, ~ cor(.x$fu_median, .x$fu_mean))
  pseudo_smry$r_u_max <- map_dbl(pseudo_smry$data, ~ cor(.x$fu_max, .x$fu_mean))
  pseudo_smry$r_med_max <- map_dbl(pseudo_smry$data, ~ cor(.x$fu_median, .x$fu_max))
  pseudo_smry$lm <- map(pseudo_smry$data, ~ lm(fu_mean ~ fu_median + fu_max, data = .x))
  pseudo_smry$pred <- map2(pseudo_smry$lm, pseudo_smry$data, ~ broom::augment(.x, newdata = .y, se_fit = TRUE, 
                                                                              interval = "prediction"))
  pseudo_smry <- pseudo_smry %>% 
    unnest(c(pred))
  ## correlation between median and 
  lm_global <- lm(fu_mean ~ fu_median, data = pseudo_smry)
  pseudo_smry$overall <- broom::augment(lm_global) %>% 
    pull(.fitted)
  r_med_max_global <- cor(pseudo_smry$fu_median, pseudo_smry$fu_mean)
  pseudo_plot1 <- ggplot(pseudo_smry, aes(x = fu_median, colour = fu_max, 
                                          y = fu_mean,
                                          group = arm_simple)) +
    geom_point(shape = 1) +
    geom_line(mapping = aes(y = .fitted), col = "red",
              linetype = "dotted") +
    geom_line(mapping = aes(y = overall), col = "blue", linetype = "dotted") +
    geom_abline(intercept = 0, slope = 1, col = "black") +
    facet_grid(~arm_simple) +
    theme_bw()
  pseudo_plot1
  coef(lm_global)
write_csv(broom::tidy(lm_global) %>% 
            mutate(overall_r = r_med_max_global), "Outputs/association_between_median_mean_fu.csv")

