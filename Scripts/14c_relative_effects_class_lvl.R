library(tidyverse)
library(rstan)
library(multinma)

source("Scripts/00_functions.R")

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

predictfor <- expand_grid(trt = 1L, male = 0:1, age = seq(40, 80, 10))
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
whoatc <- read_csv("Data/whoatcdiabetesnodose.csv")
whoatclkp <- whoatc$`ATC level name`
names(whoatclkp) <- whoatc$`ATC code`
smry <- smry %>% 
  mutate(cls = whoatclkp[cls])
plot_dc <- ggplot(smry %>% 
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
                     labels = c(0.6, 0.8, 1, 1.25, 1.67)) +
  scale_x_continuous("Age (years)") +
  # scale_color_discrete(guide = "none") +
  theme_minimal()
plot_dc
saveRDS(plot_dc, "Scratch_data/relative_mace_class_level.Rds")

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

