library(tidyverse)
tot <- readRDS("Scratch_data/rcs_mdl_data.Rds")
ref <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")

## Identify age ranges
ages <- ref %>% 
  select(ipd)
ages$ipd <- map(ages$ipd, ~ .x %>% 
                 group_by(nct_id) %>% 
                 summarise(min_age = min(age),
                           max_age = max(age)) %>% 
                             ungroup())
ages <- ages %>% 
  unnest(ipd)

## Identify placebo controlled trials
ref <- ref %>% 
  select(reg) %>% 
  unnest(reg) %>% 
  filter(reference_arm ==1) %>% 
  distinct(nct_id, trtcls5) %>% 
  mutate(reference = if_else(trtcls5 == "place", "placebo", trtcls5)) %>% 
  select(-trtcls5)

## sample from regression models 
reg <- tot %>% 
  select(drug_regime_smpl, reg2) %>% 
  unnest(reg2) %>% 
  filter(!is.na(term))

## note all nct_id2 = nct_id as NCT01778049 not in this set. 
reg_nst <- reg %>% 
  select(drug_regime_smpl, nct_id, term, estimate, arm_lvl, trtcls4, trtcls5) %>% 
  nest(cf = c(arm_lvl, trtcls4, trtcls5, term, estimate))
vcv_nst <- tot %>% 
  select(vcv2) %>% 
  unnest(vcv2) %>% 
  rename(vcv = vcv2)
reg_nst <- bind_cols(reg_nst, vcv_nst)
reg_nst$chk <- map2_lgl(reg_nst$cf, reg_nst$vcv, ~ all(.x$term == colnames(.y)))
reg_nst %>% 
  count(chk)

reg_nst$smpl <- map2(reg_nst$cf, reg_nst$vcv, ~ MASS::mvrnorm(100, .x$estimate, .y) %>% 
                       as_tibble() %>% 
                       mutate(iter = 1:nrow(.)) %>% 
                       gather("term", "value", -iter))
reg_nst <- reg_nst %>% 
  select(nct_id, smpl) %>% 
  unnest(smpl)
rm(tot, vcv_nst)  
reg <- reg %>% 
  distinct(nct_id, term, arm_lvl, trtcls5, trtcls4)
reg_nst %>% 
  anti_join(reg)
reg_nst <- reg_nst %>% 
  inner_join(reg)
rm(reg)
smpls <- reg_nst
rm(reg_nst)

## separate age to determine zs
smpls <- smpls %>% 
  mutate(z = case_when(
    str_detect(term, "^age10 ") ~ "z1",
    str_detect(term, "^age10\\' ") ~ "z2",
    str_detect(term, "^age10\\'\\'") ~ "z3"))

## estimate z_Xs ----
vals <- Hmisc::rcspline.eval(seq(-2, 1.5, 0.1), knots = c(-2, -1, 0, 1.5), inclx = TRUE)
colnames(vals) <- paste0("z", 1:3)
vals <- as_tibble(vals)
vals <- vals %>% 
  mutate(age = z1*10 + 60) %>% 
  gather("z", "zx", -age)

## add to smpls
smpls <- smpls %>% 
  inner_join(vals, by = "z", relationship = "many-to-many")


## perform calculation to obtain total effect
smpls <- smpls %>% 
  mutate(res = value*zx) %>% 
  group_by(nct_id, iter, arm_lvl, trtcls5, trtcls4, age) %>% 
  summarise(res =sum(res)) %>% 
  ungroup()

## limit to classes of interest, add reference then nest
smpls <- smpls %>% 
  filter(trtcls5 %in% paste0("A10B", c("H", "J", "K"))) %>% 
  inner_join(ref) %>% 
  nest(.by = trtcls5)

smpls$data <-map(smpls$data, ~ .x %>% 
                   arrange(nct_id, arm_lvl) %>% 
                   group_by(nct_id, reference) %>% 
                   mutate(arm_cnt = cumsum(!duplicated(arm_lvl)),
                          arm_cnt = letters[arm_cnt]) %>% 
                   ungroup())
who <- read_csv("Data/whoatcdiabetesnodose.csv")
smpls <- smpls %>% 
  inner_join(who %>% select(trtcls5= `ATC code`, dc = `ATC level name`))
smpls$data <- map(smpls$data, ~ .x %>% 
                    inner_join(ages) %>% 
                    filter(age >= min_age,
                           age <= max_age))
## drop NCT00707993 as didn't converge
smpls$plt <- map2(smpls$data, smpls$dc, ~ ggplot(.x %>% 
                                                   filter(!nct_id == "NCT00707993") %>% 
                                                   mutate(forfac = paste(nct_id, "\n vs", reference)), 
                                                 aes(x = age, y = res, colour = arm_cnt, 
                                                     group = interaction(arm_cnt, iter))) +
                   geom_line(alpha = 0.05) +
                    geom_smooth(mapping = aes(x = age, y = res, colour = arm_cnt, group = arm_cnt), method = "lm", 
                                colour = "black", linetype = "dashed", linewidth = 0.2) +
                   facet_wrap(~ forfac) +
                   theme_bw() +
                    scale_color_discrete(guide = "none") +
                    ggtitle(.y) +
                    coord_cartesian(ylim = c(-2.5, 2.5)))
saveRDS(smpls$plt, "Scratch_data/rcs_hba1c_plots.Rds")
walk2(smpls$trtcls5, smpls$plt, ~ {
  tiff(paste0("Outputs/rcs_plot_", .x, ".tiff"), res = 300, height = 12, width = 12, unit = "in", compression = "lzw")
  print(.y)
  dev.off()
})
