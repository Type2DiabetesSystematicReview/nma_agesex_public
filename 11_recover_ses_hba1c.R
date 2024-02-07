# 11_recover_ses
## very close, but I still think we need to sample from residuals
library(tidyverse)

ipd <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
ipd <- ipd %>% 
  select(-agg, -dm_nets) %>% 
  filter(drug_regime_smpl == "triple")
ipd <- ipd %>% 
  unnest(ipd) 
ipd <- ipd %>% 
  mutate(age10 = age/10) %>% 
  rename(arm = arm_f,
         value_1 = base) %>% 
  select(nct_id, nct_id2, result, age10, arm, value_1)
ipd1 <- ipd %>% 
  semi_join(ipd %>% slice(1) %>% select(nct_id))
ipd <- ipd %>% 
  group_by(nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup()
ipd$mdl <- map(ipd$data, ~ {
  lm(result ~ value_1 + age10 + arm, data = .x)
})
ipd$res <- map(ipd$mdl, broom::tidy)
ipd <- ipd %>% 
  select(nct_id, nct_id2, res) %>% 
  unnest(res) 
ipd <- ipd %>% 
  select(nct_id, nct_id2, term, estimate, std.error) %>% 
  gather("var", "pseudo", -nct_id2, -nct_id, -term)
mdls <- readRDS("Scratch_data/ipd_raw_coefs.Rds")
ipd2 <- ipd %>% 
  inner_join(mdls)

ipd2 <- ipd2 %>% 
  mutate(ps_rl = pseudo - value)

plt1 <- ggplot(ipd2 %>% filter(var == "std.error"),
               aes(x = ps_rl)) +
  geom_density() +
  facet_wrap(~ term, scales = "free")
plt1

ipd2 %>% 
  group_by(term, var) %>% 
  summarise(m = mean(ps_rl),
            s = sd(ps_rl)) %>% 
  ungroup() %>% 
  arrange(var)

## explore recovering with residuals
mod1 <- lm(result ~ age10 + value_1 + arm, data = ipd1)
library(broom)
cfs <- tidy(mod1)
dgn <- glance(mod1)
a <- augment(mod1, se_fit = TRUE)
b <- rnorm(length(ipd1$result), a$.fitted, dgn$sigma[1])
d <- rnorm(length(ipd1$result), a$.fitted, a$.se.fit)
par(mfrow = c(2,2))
hist(ipd1$result)
hist(b)
hist(d)
