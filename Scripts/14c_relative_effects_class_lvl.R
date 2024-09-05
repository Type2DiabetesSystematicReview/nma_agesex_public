library(tidyverse)
library(rstan)
library(multinma)

source("Scripts/00_functions.R")

## Read in age and arm-level information for plots ----
tot <- readRDS("Scratch_data/for_mace_regression_inter.Rds")
pseud <- tot$pseudo %>% 
  select(nct_id, arm_lvl, sex, age, trtcls5)
aggsex <- tot$mace_agg_sex %>% 
  select(nct_id, arm_lvl, sex = level_cat, participants, age_mu, age_sigma, max_age, min_age, trtcls5)
agg <- tot$mace_agg %>% 
  anti_join(aggsex %>% select(nct_id)) %>% 
  mutate(male = round(male_p * participants),
         female = participants - male) %>% 
  select(nct_id, arm_lvl, male, female, age_mu, age_sigma, max_age, min_age, trtcls5) 
agg <- agg %>% 
  gather("sex", "participants", male, female) %>% 
  arrange(nct_id, arm_lvl)
agg <- bind_rows(agg,
                 aggsex)
rm(aggsex)
agg$age_sim <- pmap(list(
  agg$participants*10, agg$age_mu, agg$age_sigma, agg$min_age, agg$max_age), function(n, m, s, l, u){
    truncnorm::rtruncnorm(n, a = l, b = u, mean = m, sd = s)  
  })
agg <- agg %>% 
  select(nct_id, arm_lvl, sex, age_sim, min_age, max_age, trtcls5) %>% 
  unnest(age_sim) %>% 
  rename(age = age_sim)
age_dist <- bind_rows(agg,
                      pseud %>% 
                        mutate(sex = if_else(sex == "F", "female", "male")))
rm(agg, pseud, tot)

## Make density curves ----
age_dist <- age_dist %>% 
  mutate(arm_type = if_else(arm_lvl == "placebo", "placebo", "active")) %>% 
  nest(data = c(age, nct_id, arm_lvl, min_age, max_age))
age_dist$n <- map_int(age_dist$data, nrow)
age_dist$m <- map_dbl(age_dist$data, ~ mean(.x$age))
age_dist$dns <- map(age_dist$data, ~ density(.x$age, adjust = 2)[c("x", "y")] %>% as_tibble())
age_dist$dns_area <- map_dbl(age_dist$dns, ~ sum(.x$y))
# map so all on the same scale, so total area sums to one
age_dist$dns <- map2(age_dist$dns, age_dist$dns_area, ~ .x %>% 
                       mutate(y = y/.y))
age_dist$dns_area <- map_dbl(age_dist$dns, ~ sum(.x$y))
## re-weight areas based on number male and female
age_dist$dns <- map2(age_dist$dns, age_dist$n, ~ .x %>% 
                       mutate(y = y*.y))
age_dist <- age_dist %>% 
  mutate(male = if_else(sex == "male", 1L, 0L))

## get nbybclass ----
mace_lng <- read_csv("Outputs/manuscript_table1b_machine_readable.csv")
mace_lng <- mace_lng %>% 
  filter(!nct_id == "UMIN000018395")
mace_lng <- mace_lng %>% 
  mutate(trtcls5 = str_sub(dc, 1, 5)) %>% 
  select(nct_id, trtcls5, participants)
nbyclass <- mace_lng %>% 
  group_by(trtcls5) %>% 
  summarise(trials = sum(!duplicated(nct_id)),
            participants = sum(participants)) %>% 
  ungroup()

## read in model ouptuts and process ----
mdl <- readRDS("FromVM/mace_classlvl/fixed_mace_classoverall.Rds")
smpl <- as.data.frame(mdl$stanfit)
rm(mdl)
betas <- smpl[, str_detect(names(smpl), "^beta")]
ds <- smpl[,!str_detect(names(smpl), "^delta") & str_detect(names(smpl), "^d")]
smpl <-smpl[, setdiff(names(smpl), c(names(betas), names(ds)))]
rm(smpl)

betas <- betas %>% 
  as_tibble() %>% 
  mutate(i = 1:4000)
ds <- ds %>% 
  as_tibble() %>% 
  mutate(i = 1:4000)
betas <- betas %>% 
  gather("param", "value", -i)
ds <- ds %>% 
  gather("param", "value", -i)
tot <- bind_rows(ds, betas)
rm(ds, betas)

novels <- c("A10BH", "A10BJ", "A10BK")
res <- map(novels, ~ tot %>% 
             filter(str_detect(param, .x)))
names(res) <- novels
res <- map(res, ~ .x %>%  spread(param, value))
map(res, names)

res <- map(res, ~ .x %>% 
             select(starts_with("d"), contains("male"), contains("age")) %>% 
             as.matrix())
map(res, colnames)
map(res, colMeans)
lpfx <- function(trt, male, age10, mymtrx) {
  as.vector(mymtrx %*% c(trt, male, age10))
}
lpfx(1, 0, 20, res$A10BJ) %>% mean

predictfor <- expand_grid(trt = 1L, male = 0:1, age = c(55, 65, 75))
predictforlst <- as.list(predictfor)
predictfor$mace <- pmap(predictforlst, function(trt, male, age) {
  b <- map2(names(res), res, ~ {
    a <- lpfx(trt, male, age-60, .y)
    tibble(cls = rep(.x, length(a)), value = a )
  }) %>% 
    bind_rows()
  b
})
predictfor <- predictfor %>% 
  unnest(mace)
smry <- predictfor %>% 
  group_by(cls, trt, male, age) %>% 
  summarise(
    m = mean(value),
    s = sd(value),
    q2_5 = quantile(value, probs = 0.025),
    q97_5 = quantile(value, probs = 0.975)) %>% 
  ungroup()

smry <- smry %>% 
  inner_join(nbyclass %>% rename(cls = trtcls5, cls_trials = trials, nwork_n = participants))

whoatc <- read_csv("Data/whoatcdiabetesnodose.csv")
whoatclkp <- whoatc$`ATC level name`
names(whoatclkp) <- whoatc$`ATC code`
smry <- smry %>% 
  mutate(trtclass = cls, 
         cls = whoatclkp[cls])
smry <- smry %>% 
  mutate(forfacet = paste0(cls, " ", cls_trials, " trials. ", 
                           formatC(nwork_n, format = "d", big.mark = ","),
                           " participants"))

smry <- smry %>% 
  mutate(sex = factor(male, levels = 1:0, 
                      labels = c("Male", "Female")))
age_dist <- age_dist %>% 
  mutate(sex = factor(male, levels = 1:0, 
                      labels = c("Male", "Female")))
saveRDS(smry, "Scratch_data/pre_mace_relative_results_ms.Rds")
  
plot_dc <- ggplot(smry,
                        aes(x = age, 
                            y = m, 
                            ymin = q2_5, 
                            ymax =  q97_5, 
                            colour = sex)) +
  geom_point(position = position_dodge(5)) +
  geom_linerange(position = position_dodge(5)) +
  facet_wrap(~forfacet, ncol = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # coord_cartesian(ylim = c(-1, 1)) +
  scale_y_continuous(name = "Treatment efficacy by age and sex (hazard ratio)",
                     breaks = log(c(1/1.5, 0.8, 1, 1.25, 1.5)) ,
                     labels = c(0.67, 0.8, 1, 1.25, 1.5)) +
  scale_x_continuous("Age (years)", breaks = c(55, 65, 75)) +
  coord_cartesian(xlim = c(50, 80), ylim = c(log(0.8)*2, log(1.5))) +
  scale_color_discrete(guide = "none") +
  theme_minimal() +
  scale_color_manual("", values = rev(c("#E69F00", "#56B4E9")), )  +
  scale_fill_manual("", values = rev(c("#E69F00", "#56B4E9")))
plot_dc
age_dist <- age_dist %>% 
  select(sex, trtcls5, dns) %>% 
  unnest(dns) 
age_dist <- age_dist %>% 
  filter(!trtcls5 %in% c("place", "A10BB"))
age_dist <- age_dist %>% 
  inner_join(smry %>% distinct(trtclass, forfacet) %>% rename(trtcls5 = trtclass))
age_dist <- age_dist %>% 
  mutate(upr = log(0.78),
         lwr = -0.4,
         upr_lwr = upr-lwr) %>% 
  ## do not group by sex here as want different distributions for men and women
  group_by(trtcls5, forfacet) %>% 
  mutate(y01 = y/(max(y)),
         upr = y01*upr_lwr + lwr) %>% 
  ungroup()
plot_dc <- plot_dc +
  geom_ribbon(data = age_dist, mapping = aes(x = x, ymin = lwr, ymax = upr, fill = sex, y = NULL), 
              alpha = 0.2,
              colour = NA) 
plot_dc
saveRDS(plot_dc, "Scratch_data/relative_mace_class_level.Rds")
ggsave(plot_dc, filename = "Outputs/rel_effects_ms.svg", height = 12, width = 7, bg = "transparent")
ggsave(plot_dc, filename = "Outputs/rel_effects_ms.pdf", height = 12, width = 7, bg = "transparent")

pdf("Outputs/relative_mace_class_level.pdf")
plot_dc
dev.off()

## make components
smry_cls <- smry %>% 
  mutate(cls2 = cls) %>% 
  group_by(cls2) %>% 
  nest() %>% 
  ungroup()
plot_dc_lst <- map(smry_cls$data, ~ ggplot(.x %>% 
                    mutate(sex = factor(male, levels = 0:1, 
                                        labels = c("Female", "Male"))),
                  aes(x = age, 
                      y = m, 
                      ymin = q2_5, 
                      ymax =  q97_5, 
                      colour = sex)) +
  geom_point(position = position_dodge(5)) +
  geom_linerange(position = position_dodge(5)) +
  facet_wrap(~cls, ncol = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # coord_cartesian(ylim = c(-1, 1)) +
  scale_y_continuous(name = "Treatment efficacy by age and sex (hazard ratio)",
                     breaks = log(c(0.6, 0.8, 1, 1.25, 1.67)) ,
                     labels = c(0.6, 0.8, 1, 1.25, 1.67),
                     limits = log(c(0.50, 1.77))) +
  scale_x_continuous("Age (years)") +
  scale_color_manual("Sex", values = c("#E69F00", "#56B4E9")) +
  theme_minimal())



saveRDS(plot_dc_lst, "Scratch_data/relative_mace_class_level_components.Rds")

