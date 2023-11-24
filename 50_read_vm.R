# 50_read_vm
library(multinma)
library(tidyverse)
library(cowplot)

ReadSummary <- function(choosedir = "mono_common_random_20231116_045721_417484") {
  res <- readRDS(paste0("FromVM/", choosedir, "/summary.Rds"))
  res2 <- res
  res2$summary <- res2$summary  %>% 
    filter(str_detect(parameter, "^beta|^d"), !str_detect(parameter, "^delta"))
  res2$sims <- res2$sims[, , res2$summary$parameter]
  # pdf(paste0("Outputs/effect_estimates_", choosedir, ".pdf"))
  # a <- plot(res2) +
  #   scale_x_continuous(limits = c(-2, 2))
  # dev.off()
  res2
}

## Get timings ----
fldrs <- c("dual_common_fixed_20231121_022309_581732/dualcomf-12528.out", 
             "dual_common_random_20231117_143718_465782/dualcommon-12044.out", 
             "dual_independent_random_20231117_213609_237076/dualindep-12047.out", 
             "mono_common_fixed_20231120_163959_242145/monocomf-12530.out", 
             "mono_common_random_20231116_045721_417484/monocommon-12043.out", 
             "mono_independent_random_20231116_145814_243569/monoindep-12046.out", 
             "triple_common_fixed_20231119_034234_118584/tripcomf-12157.out", 
             "triple_common_random_20231120_112004_635204/tripcomr-12158.out")
timings <- tibble(fldrs = fldrs)
timings$res <- map(timings$fldrs, ~ read_lines(paste0("FromVM/", .x)))
timings$res <- map(timings$res, ~ .x[str_detect(.x, "seconds") &
                                       str_detect(.x, "Warm|Sampling|Total")])
timings <- timings %>% 
  unnest(res)
timings <- timings %>% 
  mutate(res = str_remove_all(res, "Elapsed Time\\:"))
timings <- timings %>% 
  separate(res, into = c("chain", "times", "type"), sep = "\\:|\\(") %>% 
  mutate(across(c(chain, times), str_trim)) %>% 
  separate(fldrs, into = c("nwork", "comind", "fixran", "date", "time", "time2"),
           sep = "_", remove = FALSE)
timings <- timings %>% 
  select(nwork, comind, fixran, chain, type, times)
timings <- timings %>% 
  mutate(times = as.double(times %>% 
                             str_remove_all("[a-z]") %>% 
                             str_trim()),
         hours = times/60^2)
timings <- timings %>% 
  group_by(nwork, comind, fixran, type) %>% 
  summarise(max_hours = max(hours) %>% round(2),
            min_hours = min(hours) %>% round(2),
            hours = paste0(min_hours, " to ", max_hours)) %>% 
  ungroup() %>% 
  select(-min_hours, -max_hours)
timings_wide <- timings %>% 
  mutate(type = str_sub(type, 1, -2) %>% 
           str_to_lower() %>% 
           str_replace_all("\\-", "_"),
         type = case_when(
           type == "warm_up" ~ "a_warm_up",
           type == "sampling" ~ "b_sampling",
           type == "total" ~ "c_total"
         )) %>% 
  arrange(type) %>% 
  pivot_wider(names_from = type, values_from = hours)
write_csv(timings_wide, "Outputs/Model_timings.csv")
rm(timings, timings_wide)

### results. Common models
common <- list.files("FromVM/", patt = "common")
res <- map(common, ReadSummary)
res <- tibble(fldrs = common, res = res)
res <- res %>% 
  separate(fldrs, into = c("nwork", "comind", "fixran", "date", "time", "time2"),
           sep = "_", remove = FALSE)
res$summary <- map(res$res, ~.x$summary %>% 
  separate(parameter, into = c("stem", "branch"), sep = "\\[") %>% 
  mutate(branch = str_sub(branch, 1, -2) %>% 
           str_remove("\\.trtclass"),
         branch = case_when(
           str_detect(branch, "\\:") ~ branch,
           branch %in% c("base", "age", "sex") ~ paste0(branch, ":none"),
           TRUE ~ paste0("main:", branch))) %>% 
  separate(branch, into = c("covariate", "arm"), sep = "\\:"))
res$summary <- map(res$summary, function(mydf) mydf %>% 
  mutate(across(c(`mean`, `2.5%`, `97.5%`), ~ if_else(covariate == "age", 30*.x, .x))))
res <- res %>% 
  select(-res) %>% 
  unnest(summary)
DCExtract <- function(x) {
  x2 <- str_replace_all(x, pattern = "_{1,2}[0-1]{1,1}", "9")
  x2 <- str_replace(x2, "A10BK_05", "A10BK")
  x2 <- str_replace_all(x2, pattern = "#", "")
  x2 <- str_split(x2, pattern = "_")
  map_chr(x2, ~ str_sub(.x, 1, 5) %>% paste(collapse = "_"))
}
res <- res  %>% 
  mutate(covariate = factor(covariate, levels = c("none", "base","main", "age", "sex"))) %>% 
  mutate(dc = DCExtract(arm),
         dc = if_else(dc == "A10AA", "A10A", dc),
         dc = if_else(dc == "A10BK_06", "A10BK", dc))
all_classes <- res %>% 
  count(dc, nwork) %>% 
  spread(nwork, n)
res <- res %>% 
                 mutate(dc_lbl = factor(dc,
                                    levels = c('none',
                                               'A10BA',
                                               'A10BB',
                                               'A10BF',
                                               'A10BG',
                                               'A10BH',
                                               'A10BJ',
                                               'A10BK',
                                               'A10BX',
                                               'combi',
                                               'OAD',
                                               'A10A',
                                               'A10A_A10BB',
                                               'A10A_A10BJ',
                                               'A10BH_A10BK',
                                               "A10A_A10BH"),
                                    labels = c('covariate effects',
                                               'A10BA - Biguanides',
                                               'A10BB - Sulphonylureas',
                                               'A10BF - Alpha glucosidase inhibitors',
                                               'A10BG - Thiazolidinediones (glitzones)',
                                               'A10BH - DPP-4 inhibitors',
                                               'A10BJ - GLP-1 analogues',
                                               'A10BK - SGLT2 inhibitors',
                                               'A10BX - Other blood lowering drugs',
                                               'combi - combinations',
                                               'OAD - any oral antidiabetic',
                                               'A10A - insulins',
                                               'A10A_A10BB - insulins and SUs',
                                               'A10A_A10BJ - insulin and GLP1',
                                               'A10BH_A10BK - insulin and SGLT2',
                                               'A10BH_A10BK - DPP4 SGLT2')))
res <- res %>% 
  mutate(nwork = factor(nwork, levels = c("mono", "dual", "triple")))
res_nst <- res %>% 
  group_by(fixran, dc, dc_lbl) %>% 
  nest() %>% 
  ungroup()
res_nst$plot <- map2(res_nst$data, res_nst$dc_lbl, ~ ggplot(.x, aes(x = interaction(covariate, arm), 
                             y = mean, ymin = `2.5%`, ymax = `97.5%`, 
                             colour = covariate)) +
  geom_point() + 
  geom_linerange() +
  # scale_y_continuous(limits = c(-1, 1)) +
  facet_wrap(~nwork, ncol = 3) +
  coord_flip(ylim = c(-2, 2)) +
    scale_x_discrete("") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggtitle(.y))
lgnd <- res_nst$plot[[1]] %>% get_legend()
res_nst$plot <- map(res_nst$plot, ~ .x + theme(legend.position="none"))

res_nst <- res_nst %>% 
  mutate(plotset = 
           case_when(dc  == "none" ~ "age_sex_effects",
                     dc %in% c("A10",
                               "A10A_A10BJ",
                               "A10A_A10BB",
                               "A10A_A10BH") ~ "insulin_and_combos",
                     dc %in% paste0("A10B", LETTERS[1:6]) ~ "older_oads_class_single",
                     dc %in% paste0("A10B", LETTERS[7:11]) ~ "newer_oads_class_single",
                     TRUE ~ "other"))
res_nst2 <- res_nst %>% 
  select(fixran, plotset, plot) %>% 
  nest(data = plot)
## fixed
fixed <- res_nst2 %>% 
  filter(fixran == "fixed")
random <- res_nst2 %>% 
  filter(fixran == "random")
random$data <- map(random$data, pull)
random$cowplot <- map(random$data, ~ cowplot::plot_grid(plotlist = .x, ncol = 1, axis = "b", align = "v"))
fixed$data <- map(fixed$data, pull)
fixed$cowplot <- map(fixed$data, ~ cowplot::plot_grid(plotlist = .x, ncol = 1, axis = "b", align = "v"))

saveRDS(random, "Scratch_data/re_common_plots.Rds")
saveRDS(fixed, "Scratch_data/fe_common_plots.Rds")

## examine DIC for independent and commpn
foldic <- list.files("FromVM/")
dic <- map(foldic, ~ readRDS(paste0("FromVM/", .x, "/dic.Rds"))[c(1, 2, 4)] %>% 
             as_tibble())
names(dic) <- foldic
dic <- bind_rows(dic, .id = "model") %>% 
  separate(model, c("nwork", "sharing", "fixed_rand", "daterun", "a", "b")) %>% 
  select(-a, -b)
dic <- dic %>% 
  arrange(nwork, sharing, fixed_rand)
write_csv(dic, "Outputs/dics.csv")

## Independent plots ----
folci<- c("dual_common_random_20231117_143718_465782", 
            "dual_independent_random_20231117_213609_237076", 
             "mono_common_random_20231116_045721_417484", 
             "mono_independent_random_20231116_145814_243569")
ci <- map(folci, ~ readRDS(paste0("FromVM/", .x, "/summary.Rds"))$summary)
names(ci) <- folci
onlycommon <- setdiff(ci$dual_common_random_20231117_143718_465782$parameter, ci$dual_independent_random_20231117_213609_237076$parameter)
onlyindep <- setdiff(ci$dual_independent_random_20231117_213609_237076$parameter, ci$dual_common_random_20231117_143718_465782$parameter)
ci <- ci %>% 
  bind_rows(.id = "model")
ci <- ci %>% 
  separate(model, c("nwork", "sharing", "fixed_rand", "daterun", "a", "b")) %>% 
  select(-a, -b)
ci <- ci %>% 
  filter( (sharing == "common" & parameter %in% onlycommon) |
            (sharing == "independent"  & parameter %in% onlyindep))
ci <- ci %>% 
  mutate(varname = case_when(
    str_detect(parameter, "age") & str_detect(parameter, "\\:") ~ "age_inter",
    str_detect(parameter, "sex") & str_detect(parameter, "\\:") ~ "sex_inter",
    str_detect(parameter, "age") & !str_detect(parameter, "\\:") ~ "age",
    str_detect(parameter, "sex") & !str_detect(parameter, "\\:") ~ "sex"))
ci <- ci %>% 
  mutate(drug_or_dc = parameter %>% 
           str_remove("beta\\[\\.trt") %>% 
           str_remove("\\:(age|sex)\\]") %>% 
           str_remove("beta\\[(age|sex)\\:\\.trtclass") %>% 
           str_remove("\\]"),
         dc = str_sub(drug_or_dc, 1, 5),
         dc = if_else(str_detect(dc, "^A10A"), "A10A", dc),
         across(c(`mean`, `2.5%`, `97.5%`), ~ if_else(varname == "age_inter", 30*.x, .x)))
plot_ci <- ggplot(ci %>% filter(str_detect(dc, "A10B[H-K]"),
                                !str_detect(drug_or_dc, "_A")), aes(x = drug_or_dc, y = mean, ymin = `2.5%`, ymax = `97.5%`, 
                          colour = sharing)) +
  geom_point() + 
  geom_linerange(alpha = 0.5) +
  facet_grid(dc + nwork~varname, scales = "free_y") +
  coord_flip(ylim = c(-2, 2)) +
  scale_x_discrete("") +
  scale_color_discrete(guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)
plot_ci
saveRDS(plot_ci, "Scratch_data/plot_independent_vs_common.Rds")

## Model fit statistics that are available from summary ----
fldrs <- c("dual_common_fixed_20231121_022309_581732/summary.Rds", "dual_common_random_20231117_143718_465782/summary.Rds", 
           "dual_independent_random_20231117_213609_237076/summary.Rds", 
           "mono_common_fixed_20231120_163959_242145/summary.Rds", "mono_common_random_20231116_045721_417484/summary.Rds", 
           "mono_independent_random_20231116_145814_243569/summary.Rds", 
           "triple_common_fixed_20231119_034234_118584/summary.Rds", "triple_common_random_20231120_112004_635204/summary.Rds")
ss <- map(fldrs, ~ readRDS(paste0("FromVM/", .x))$summary)
names(ss) <- fldrs
ss <- ss %>% 
  bind_rows(.id = "model")
ss <- ss %>% 
  separate(model, c("nwork", "sharing", "fixed_rand"), sep = "_", extra = "drop") 
ss <- ss %>% 
  separate(parameter, into = c("trunc", "branc"), sep = "\\[")
ss_table <- ss %>% 
  filter(is.na(branc)) %>% 
  select(nwork, sharing, fixed_rand, parameter = trunc, Bulk_ESS, Tail_ESS, Rhat)
ss <- ss %>% 
  filter(!is.na(branc))
bulk_ess <- ggplot(ss, aes(x = trunc, 
                          y = Bulk_ESS, fill = sharing, colour = sharing)) +
  geom_violin() +
  facet_wrap(~nwork + fixed_rand, scales = "free_y") +
  coord_flip()
tail_ess <- bulk_ess + geom_violin(mapping = aes(y = Rhat))
rhat <-bulk_ess + geom_violin(mapping = aes(y = Rhat))
saveRDS(list(bulk_ess = bulk_ess, tail_ess = tail_ess, rhat= rhat, ss_table = ss_table), "Scratch_data/diagnostic_summaries.Rds")


