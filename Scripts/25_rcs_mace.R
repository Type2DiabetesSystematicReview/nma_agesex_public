# 27_rcs_mace
library(tidyverse)

CnvrtCorrMatrix <- function(a){
  ## recovery whole matrix by duplication
  # browser()
  allnames <- union(a$row, a$col)
  a <- bind_rows(a,
                 a %>%
                   rename(row = col, col = row),
                 tibble(row = allnames, col = allnames, r = 1)) %>%
    distinct()
  # convert into matrix format
  a <- a %>%
    spread(col, r)
  a_row <- a$row
  a$row <- NULL
  a <- as.matrix(a)
  if (any(is.na(a))) warning("Missing values in matrix")
  rownames(a) <- a_row
  a
}
cor2cov <- function(V, sigma) {
  p <- (d <- dim(V))[1L]
  if (!is.numeric(V) || length(d) != 2L || p != d[2L])
    stop("'V' is not a square numeric matrix")
  if (length(sigma) != p)
    stop("'sigma' is not a vector comformable as the standard deviations of 'V'")
  if (any(diag(V) != 1))
    warning("diag(.) contained non 1 entries.  Did you pass a correlation matrix?")
  sigma * V * rep(sigma, each = p)
}

## armlkp 
armlkp <- (readRDS("Scratch_data/for_mace_regression_inter.Rds")$cfs) %>%
  filter(!is.na(arm_lvl), !arm_lvl == "placebo")  %>% distinct(nct_id, trt, arm_lvl, trtcls5)
who <- read_csv("Data/whoatcdiabetesnodose.csv")
armlkp <- armlkp %>% 
  inner_join(who %>% select(trtcls5= `ATC code`, dc = `ATC level name`))
rm(who)

cfs <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv")
cfs <- cfs %>% 
  filter(outcome == "mace", models == "rcs")
cfs <- cfs %>% 
  separate(term, into = c("cvrt", "trt"), sep = " \\*", remove = FALSE, fill = "right") %>% 
  mutate(tform = 
  case_match(cvrt,
             "age10c" ~ "a",
             "age10c'" ~ "b",
             "age10c''" ~ "c",
             "age10c'''" ~ "d")) %>% 
  filter(tform %in% letters[1:4], !is.na(trt)) 
cfs <- cfs %>% 
  mutate(trt = str_remove(trt, " arm_f\\="))
vcv <- read_csv("Data/agesexhba1cmaceupdate_9492/vcov.csv") %>% 
  semi_join( cfs %>% rename(row = term)) %>% 
  semi_join( cfs %>% rename(col = term))
vcv <- vcv %>% 
  select(nct_id, row, col, r) %>% 
  nest(.by = nct_id)
cfs <- cfs %>% 
  select(nct_id, term, tform, estimate, std.error) %>% 
  nest(.by = nct_id) %>% 
  rename(cfs = data)
cfs <- cfs %>% 
  inner_join(vcv) %>% 
  rename(vcv = data)
cfs$r <- map(cfs$vcv, CnvrtCorrMatrix)

## Convert to variance-covariance matrix
cfs$vcv <- NULL
cfs$chk <- map2_lgl(cfs$cfs, cfs$r, ~ all(rownames(.y) %in% .x$term))
cfs$r <- map2(cfs$r, cfs$cfs, ~ .x[.y$term, .y$term])
cfs$chk <- map2_lgl(cfs$cfs, cfs$r, ~ all(rownames(.y) == .x$term))
cfs$vcv <- map2(cfs$r, cfs$cfs, ~ cor2cov(.x, .y$std.error))
cfs$vcv <- map(cfs$vcv, ~ Matrix::nearPD(.x)$mat)

## sample from distribution ----
cfs$smpl <- map2(cfs$cfs, cfs$vcv, ~ MASS::mvrnorm(100, .x$estimate, .y) %>% 
                       as_tibble() %>% 
                       mutate(iter = 1:nrow(.)) %>% 
                       gather("term", "value", -iter))
smpls <- cfs %>% 
  select(nct_id, smpl) %>% 
  unnest(smpl)
smpls <- smpls %>% 
  separate(term, into = c("term", "trt"), sep = " \\* arm_f\\=")
rm(cfs, tot, vcv, vcv_nst)

smpls <- smpls %>% 
  mutate(z = case_match(term,
    "age10c" ~ "z1",
    "age10c'" ~ "z2",
    "age10c''" ~ "z3",
    "age10c'''" ~ "z4"))
## obtain age distribution
age_dist <- readRDS("Scratch_data/for_mace_regression_inter.Rds")$pseudo
age_dist <- age_dist %>% 
  select(nct_id, age) %>% 
  nest(.by = nct_id)
age_dist <- age_dist %>% 
  semi_join(smpls %>% select(nct_id))
age_dist$ageq <- map(age_dist$data, ~ tibble(q = c(0.05, 0.27, 0.5, 0.72, 0.95),
                                             age = quantile(.x$age, c(0.05, 0.27, 0.5, 0.72, 0.95))))
age_dist <- age_dist %>% 
  select(nct_id, ageq) %>% 
  unnest(ageq)
age_dist <- age_dist %>% 
  mutate(age10c = round((age-60)/10, 1))
age_dist %>% 
  select(-age) %>% 
  spread(nct_id, age10c)
age_dist <- age_dist %>% 
  select(nct_id, q, age10c) %>% 
  distinct()
age_dist <- age_dist %>% 
  select(-q) %>% 
  nest(.by = nct_id)
age_dist$vals <- map(age_dist$data, ~ {
  vals <- Hmisc::rcspline.eval(seq(-1.5, 1.5, 0.1), knots = .x$age10c, inclx = TRUE)
  colnames(vals) <- paste0("z", 1:4)
  vals <- as_tibble(vals)
  vals <- vals %>% 
    mutate(age = z1*10 + 60) %>% 
    gather("z", "zx", -age)
  vals
})
vals <- age_dist %>% 
  select(nct_id, vals) %>% 
  unnest(vals)

## add to smpls
smpls <- smpls %>% 
  inner_join(vals, by = c("z", "nct_id"), relationship = "many-to-many")

## perform calculation to obtain total effect
smpls <- smpls %>% 
  mutate(res = value*zx) %>% 
  group_by(nct_id, trt, iter, age) %>% 
  summarise(res =sum(res)) %>% 
  ungroup()
smpls <- smpls %>% 
  inner_join(armlkp)
smpls <- smpls %>% 
  mutate(forfacet = paste0(dc, "\n", nct_id),
         arm_lbl = arm_lvl %>% 
           str_replace("_", " ") %>% 
           paste0(" mg") %>% 
           str_to_sentence())

plt <- ggplot(smpls, aes(x = age, y = res, group = interaction(arm_lbl, iter), colour = arm_lbl)) +
  geom_line(alpha = 0.05) +
  geom_smooth(mapping = aes(x = age, y = res, colour = arm_lbl, group = arm_lbl), method = "lm", 
              colour = "black", linetype = "dashed", linewidth = 0.2) +
  facet_wrap(~dc + nct_id, ncol = 2, scales = "free_y") +
  scale_colour_discrete("", guide = guide_legend(override.aes = list(alpha = 1))) +
  theme_minimal()
plt
saveRDS(plt, "Scratch_data/mace_rcs_plots.Rds")
