#22 death competing
library(tidyverse)

dth <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv") %>% 
  filter(outcome == "death")
arm_map <- read_csv("nct_id,drug_name,drug_dose,drug_unit,term
NCT00968708,alogliptin,25,milligram,arm_falogliptin
NCT00968708,placebo,,,placebo
NCT00968708,placebo,,,placebo
NCT01032629,canagliflozin,100,milligram,arm_fjnj-28431754-100 mg
NCT01032629,canagliflozin,300,milligram,arm_fjnj-28431754-300 mg
NCT01032629,placebo,,placebo,placebo
NCT01989754,canagliflozin,100,milligram,arm_fcanagliflozin
NCT01989754,placebo,,,placebo
NCT02065791,canagliflozin,100,milligram,arm_fcana 100 mg
NCT02065791,placebo,,,placebo
NCT02465515,albiglutide,30,milligram,arm_falbiglutide
NCT02465515,placebo,,,placebo
NCT01131676,Empagliflozin,10,milligram,arm_fbi 10773 10mg
NCT01131676,Empagliflozin,25,milligram,arm_fbi 10773 25mg")

arm_map <- arm_map %>%
  mutate(arm_lbl = case_when(
  term == "placebo" ~ "Placebo",
  TRUE ~ paste(drug_name, drug_dose)))

dth <- dth %>% 
  mutate(term_orig = term)
for (i in seq_along(arm_map$nct_id)) {
  dth <- dth %>% 
    mutate(term = str_replace(term, arm_map$term[[i]], arm_map$arm_lbl[[i]]))
}
dth <- dth %>% 
  mutate(term = term %>% 
           str_replace("sexM", "Male") %>% 
           str_replace("age10c", "Age (decades)"))
dth <- dth %>% 
  mutate(models_f = factor(models, levels = paste0("f", 6:10),
                           labels = c("arm",
                                      "arm + age + sex",
                                      "arm*age + sex",
                                      "arm + age*sex",
                                      "arm*age + arm*sex")))
dth <- dth %>% 
  mutate(lci = estimate - 1.96*std.error,
         uci = estimate + 1.96*std.error)
plots <- ggplot(dth, aes(x = term, y = estimate, ymin = lci, ymax = uci, colour = models_f)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(position = position_dodge(0.5)) +
  facet_wrap(~nct_id, scales = "free") +
  coord_flip() +
  scale_color_discrete("") 
plots
pdf("Outputs/compete_death.pdf", height = 12, width = 20)
plots
dev.off()

tiff("Outputs/compete_death.tiff", height = 12, width = 20, units = "in", res = 300, compression = "lzw")
plots
dev.off()
saveRDS(plots, "Scratch_data/competing_risk_model_deaths.Rds")

deaths1 <- read_csv("Data/agesexhba1cmaceupdate_6115/deaths_summary.csv") %>% 
  rename(pop = n)
deaths2 <- read_csv("Data/agesexhba1cmaceupdate_9492/death_summary.csv") %>% 
  rename(death = deaths)
deaths3 <- read_csv("Data/agesexhba1cmaceupdate_9492/diag_cox.csv") %>% 
  filter(outcome == "death", models == "f6") %>% 
  select(nct_id, pop = nobs, death = nevent) %>% 
  distinct() %>% 
  mutate(prcnt = 100*death/pop,
         prcnt = round(prcnt, 1),
         death = as.character(death))

deaths <- bind_rows(deaths3 %>% select(nct_id, pop, death, prcnt), 
                    deaths1 %>% select(nct_id, pop, death, prcnt), 
                    deaths2 %>% select(nct_id, pop, death, prcnt)) %>% 
  distinct(nct_id, .keep_all = TRUE)
deaths <- deaths %>% 
  mutate(mace = if_else(nct_id %in% dth$nct_id, "MACE", "Non-MACE"))
pltdeath <- ggplot(deaths, aes(x = prcnt, fill = mace)) +
  geom_histogram() +
  scale_y_continuous("Count of trials") +
  scale_x_continuous("Percentage of participants who died") +
  scale_fill_discrete("")
saveRDS(pltdeath, "Scratch_data/death_histogram.Rds")
tiff("Outputs/death_histogram.tiff", res = 300, height = 6, width = 6, units = "in")
pltdeath
dev.off()
