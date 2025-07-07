# 06_fit_hba1c_reg
## Designed to run inside virtual machine for speed using Rscript with arguments
library(tidyverse)
library(igraph)

dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
tot <- readRDS("Scratch_data/agg_ipd_hba1c.Rds")
tot <- tot %>% 
  select(drug_regime_smpl, agg, ipd)
tot$ipd <- map(tot$ipd, ~ .x %>%
                 group_by(nct_id) %>%
                 mutate(ntot = length(trtcls4)) %>%
                 ungroup() %>%
                 distinct(nct_id, arm_lvl, trtcls5, trtcls5, ntot))
tot$agg <- map(tot$agg, ~ .x %>%
                 group_by(nct_id) %>% 
                 mutate(ntot = sum(n)) %>% 
                 ungroup() %>% 
                 distinct(nct_id, arm_lvl, trtcls5, trtcls5, ntot))
tot <- tot %>% 
  gather("datalvl", "data", -drug_regime_smpl) %>% 
  unnest(data) 
tot <- tot %>% 
  filter(!nct_id %in% dropdisconnect)
neach <- tot$ntot %>% unique() %>% sort()
tot <- tot %>% 
  group_by(nct_id) %>% 
  filter(n() > 1) %>% 
  ungroup()
tot_unnest <- tot
tot <- tot %>% 
  nest(.by = drug_regime_smpl)
## create edge dataframe
tot$data <- map(tot$data, ~ .x %>% 
  group_by(nct_id, ntot) %>%
  summarise(pairs = list(as.data.frame(t(combn(arm_lvl, 2, simplify = TRUE)),
                                       stringsAsFactors = FALSE)),
            .groups = "drop") %>%
  unnest(pairs) %>% 
    select(from  = V1, to = V2, ntot))


res <- map(neach, function(ncut) {
  print(ncut)
connections <- map2(tot$drug_regime_smpl, tot$data, ~ {
## drop trials which are below cut off then
  edge_df <- .y %>% 
    filter(ntot >= ncut) %>% 
    select(-ntot)
  # For each trial, generate all unique treatment pairs
## create a graph
g <- graph_from_data_frame(edge_df, directed = FALSE)
## obtain components
comp <- components(g)
## obtain placebo component ID
placebo_component_id <- comp$membership[which(V(g)$name == "placebo")]
# Extract the names of all treatments in the placebo component.
arms_connected_to_placebo <- names(comp$membership)[comp$membership == placebo_component_id]
arms_not_connected_to_placebo = setdiff(names(comp$membership), arms_connected_to_placebo)

arms_connected_to_placebo <- tibble(nwork = .x, 
                                    grp = "arms_connected_to_placebo", 
                                    arm_lvl = arms_connected_to_placebo)
arms_not_connected_to_placebo <- tibble(nwork = .x, 
                                        grp = "arms_not_connected_to_placebo", 
                                        arm_lvl = arms_not_connected_to_placebo)
tot <- bind_rows(arms_connected_to_placebo,
             arms_not_connected_to_placebo)
# list(nwork = .x,
#      connected = arms_connected_to_placebo,
#      disconnected = arms_not_connected_to_placebo)
})
  bind_rows(connections)
})
names(res) <- paste0("ncut", neach)
res <- bind_rows(res, .id = "ncut")
saveRDS(list(tot = tot_unnest, armbycut = res), "Scratch_data/exclude_by_size.Rds")

res <- readRDS("Scratch_data/exclude_by_size.Rds")
tot <- res$tot
armbycut <- res$armbycut
armbycut <- armbycut %>% 
  inner_join(tot %>% distinct(arm_lvl, trtcls5))
NFunction <- function(arm_lvl_slct, nwork_slct, ncut) {
  tot %>% 
    filter(ntot >= ncut,
           arm_lvl %in% arm_lvl_slct,
           drug_regime_smpl %in% nwork_slct) %>% 
    distinct(nct_id, ntot) %>% 
    summarise(ntrials = sum(!duplicated(nct_id)),
              ntot = sum(ntot)) %>% 
    ungroup()
}
armbycut_smry <- armbycut %>% 
  group_by(ncut, nwork, grp, trtcls5) %>% 
  nest() %>% 
  ungroup()
armbycut_smry <- armbycut_smry %>% 
  mutate(ncut = str_remove(ncut, "ncut") %>% as.integer())
armbycut_smry$smry <- pmap(list(armbycut_smry$nwork,
                                armbycut_smry$data, 
                                armbycut_smry$ncut), 
                           function(nwork, data, ncut) {
                                  print(ncut)
                                  NFunction(arm_lvl_slct = data$arm_lvl, 
                                            nwork_slct = nwork, ncut = ncut)
                                  })
armbycut_smry <- armbycut_smry %>% 
  select(-data) %>% 
  unnest(smry)
armbycut_smry2 <- armbycut_smry %>% 
  rename(participants = ntot) %>% 
  mutate(grp = if_else(grp == "arms_connected_to_placebo", "connected", "unnconnected")) %>% 
  pivot_wider(names_from = c(trtcls5, grp), values_from = c(ntrials, participants), values_fill = 0L)
## take only cut-points which change the network. There are only 143 of these
armbycut_smry3 <- armbycut_smry2 %>% 
  distinct(nwork, ntrials_A10BK_connected, ntrials_A10BB_connected, ntrials_A10BH_connected, 
           ntrials_A10BF_connected, ntrials_A10A_connected, ntrials_A10BJ_connected, ntrials_A10BG_connected, 
           ntrials_A10BA_connected, ntrials_place_connected, ntrials_A10B_connected, ntrials_A10BX_connected, 
           ntrials_A10BB_unnconnected, ntrials_A10BH_unnconnected, ntrials_A10BG_unnconnected, 
           ntrials_A10BF_unnconnected, ntrials_A10BJ_unnconnected, ntrials_A10BA_unnconnected, 
           ntrials_A10A_unnconnected, ntrials_A10BK_unnconnected, ntrials_A10B_unnconnected, 
           participants_A10BK_connected, participants_A10BB_connected, participants_A10BH_connected, 
           participants_A10BF_connected, participants_A10A_connected, participants_A10BJ_connected, 
           participants_A10BG_connected, participants_A10BA_connected, participants_place_connected, 
           participants_A10B_connected, participants_A10BX_connected, participants_A10BB_unnconnected, 
           participants_A10BH_unnconnected, participants_A10BG_unnconnected, participants_A10BF_unnconnected, 
           participants_A10BJ_unnconnected, participants_A10BA_unnconnected,
           participants_A10A_unnconnected, participants_A10BK_unnconnected, participants_A10B_unnconnected,
           .keep_all = TRUE) 
write_csv(armbycut_smry3, "Outputs/trials_participants_by_class_cut_size.csv")
totals <- armbycut %>% 
  group_by(ncut, nwork, grp) %>% 
  nest() %>% 
  ungroup()
totals <- totals %>% 
  mutate(ncut = str_remove(ncut, "ncut") %>% as.integer()) 
totals$smry <- pmap(list(totals$nwork,
                         totals$data, 
                         totals$ncut), 
                    function(nwork, data, ncut) {
                      print(ncut)
                      NFunction(arm_lvl_slct = data$arm_lvl, 
                                nwork_slct = nwork, ncut = ncut)
                    })
totals_smry <- totals %>% 
  select(-data) %>% 
  unnest(smry)
totals_smry2 <- totals_smry %>% 
  rename(participants = ntot) %>% 
  mutate(grp = if_else(grp == "arms_connected_to_placebo", "connected", "unnconnected")) %>% 
  pivot_wider(names_from = c(grp), values_from = c(ntrials, participants), values_fill = 0L)
totals_smry2 <- totals_smry2 %>% 
  distinct(ntrials_connected, ntrials_unnconnected, participants_connected, participants_unnconnected, .keep_all = TRUE)

totals_smry2 <- totals_smry2 %>% 
  arrange(ncut)
write_csv(totals_smry2, "Outputs/trials_participants_cut_size.csv")
