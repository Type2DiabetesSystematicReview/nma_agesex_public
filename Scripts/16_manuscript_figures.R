# 16_manuscript_figures
library(tidyverse)
library(cowplot)
library(patchwork)
theme_minimal2 <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                            base_rect_size = base_size/22) 
{
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(
      legend.background = element_blank(), 
      legend.key = element_blank(), 
      panel.background = element_blank(), 
      panel.border = element_blank(), 
      strip.background = element_blank(), 
      plot.background = element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = 'bottom',
      complete = TRUE)
}

regplots <- readRDS("Scratch_data/regplots.Rds")
clsplts <- regplots[c("nthba1c", "ntmace")]
clsplts <- map(clsplts, ~ .x + scale_color_discrete("Parameterisation of trial estimates within drug classes"))

# + theme_minimal2()
leg <- map(clsplts, cowplot::get_legend)[[1]]
clsplts <- map(clsplts, ~ .x + scale_color_discrete(guide = "none"))

res <- clsplts[[1]] / (clsplts[[2]] + leg)
ggsave(plot = res, filename = "Outputs/hba1c_mace.tiff", dpi = 300, height = 7, width = 14, compression = "lzw")


releffects <- readRDS("Scratch_data/relative_mace_class_level_components.Rds")
# releffects[1] <- map(releffects[1], ~ .x + scale_x_continuous("", guide = "none"))
res2 <- releffects[[1]] + releffects[[2]] + releffects[[3]] + guide_area() +
  plot_layout(guides = "collect",
              axes = "collect",
              axis_titles = "collect")
res2
ggsave(plot = res2, filename = "Outputs/rel_effects.tiff", dpi = 300, height = 8, width = 14, compression = "lzw")
