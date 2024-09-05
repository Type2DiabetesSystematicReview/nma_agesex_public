# 23_rcs_hba1c
library(tidyverse)
tot <- read_rds("Scratch_data/agg_ipd_hba1c.Rds") %>% 
  select(drug_regime_smpl, reg)

## review number of arms in each trial
rv <- tot %>% 
  select(reg) %>% 
  unnest(reg) %>%
  mutate(rcs = if_else(models == "rcs", 1L, 0L)) %>% 
  select(nct_id, nct_id2, arm_id_unq, term, arm_lvl, rcs)
rv_smry <- rv %>% 
  group_by(nct_id, nct_id2, rcs) %>% 
  summarise(arm_n = sum(!duplicated(arm_lvl[!is.na(arm_lvl)])),
            arm_lvls = paste(arm_lvl %>% sort() %>% unique(), collapse = ";")) %>% 
  ungroup() %>% 
  arrange(arm_n)
rv_smry <- rv_smry %>% 
  pivot_wider(names_from = rcs, values_from = c(arm_n, arm_lvls))
rv %>% 
  filter(nct_id == "NCT00688701")

## drop trials with only single arms left
drp <- c("NCT00688701", "NCT00712673", "NCT00763451")
tot$reg <- map(tot$reg, ~ .x %>% 
                  filter(models %in% "rcs", !nct_id %in% drp))

# 
## 92 trials; fewer than main analysis due to exclusions
tot$reg <- map(tot$reg, ~ .x %>% 
        separate(term, into = c("tform", "others"), sep = " ", remove = FALSE, extra = "merge") %>% 
          mutate(tform = 
                   case_match(tform,
                              "age10" ~ "a",
                              "age10'" ~ "b",
                              "age10''" ~ "c",
                              "age10'''" ~ "d")))


tot$reg <- map(tot$reg, ~ .x %>% 
                 mutate(a = if_else(tform == "a", 1L, NA_integer_),
                        b = if_else(tform ==  "b", 1L, NA_integer_),
                        c = if_else(tform == "c", 1L, NA_integer_),
                        d = if_else(tform == "d", 1L, NA_integer_),
                        across(c(a, b, c, d), ~ if_else(is.na(term), 0L, .x)),
                        value_1 = case_when(
                          is.na(term) ~ 0,
                          str_detect(term, "value_1") ~ 1,
                          TRUE ~ NA_real_),
                        male = case_when(
                          sex == TRUE ~ 1L,
                          sex == FALSE ~ 0L,
                          TRUE ~ NA_integer_)))

tot$reg <- map(tot$reg, ~ .x %>% 
                 select(nct_id, nct_id2, models, arm_id_unq, term, estimate, std.error, arm_lvl, trtcls5, trtcls4, male, a, b, c, d, value_1, crl, vcv))
tot$vcv <- map(tot$reg, ~ .x %>% 
                 filter(is.na(term)) %>% 
                 pull(vcv))
tot$vcv <- map2(tot$vcv, tot$reg, ~ {
  names(.x) <- unique(.y$nct_id[!is.na(.y$nct_id)])
  .x
})
tot$reg <- map(tot$reg, ~ .x %>% 
                 select(-vcv, -crl))
tot$reg2 <- map(tot$reg, ~ .x %>% 
                  filter(is.na(term) | 
                           (str_detect(term, "age10") & str_detect(term, "arm"))))
tot$vcv2 <- map(tot$vcv, function(vcv1) {
vcv2 <- map(vcv1, ~ {
    term <- rownames(.x)
    slct <- str_detect(term, "age10") & str_detect(term, "arm")
    .x[slct,slct]
    })
vcv2
})

## limit to just the age:arm interactions
saveRDS(tot, "Scratch_data/rcs_mdl_data.Rds")
