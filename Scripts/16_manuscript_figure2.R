# 16_manuscript_figures
library(tidyverse)
library(patchwork)
ntmace <- readRDS("Scratch_data/regplots.Rds")$ntmace
ggsave(ntmace, filename = "Outputs/mace_for_ms.svg", height = 4, width = 8, bg = "transparent")

nthba1c <- readRDS("Scratch_data/hba1c_plots.Rds")$main
ggsave(nthba1c, filename = "Outputs/hba1c_for_ms.svg", height = 4, width = 12, bg = "transparent")
