#04b_arm_label_reverse
## ciode to reverse who_atc lookups back to drug names
library(tidyverse)
who_atc <- read_csv("../nma_agesex_public/Data/whoatcdiabetesnodose.csv")
who_lkp <- readRDS("../nma_agesex_public/Scratch_data/who_atc_lkp.Rds")

who_lkp_rev1 <- who_atc$`ATC level name`
names(who_lkp_rev1) <- who_atc$`ATC code`
who_lkp_rev2 <- names(who_lkp)
names(who_lkp_rev2) <- who_lkp

who_lkp_rev <- c(who_lkp_rev1, who_lkp_rev2)
who_lkp_rev <- who_lkp_rev[!duplicated(names(who_lkp_rev))]
who_lkp_rev <- who_lkp_rev[!str_detect(who_lkp_rev, "\\|")]
who_lkp_rev <- who_lkp_rev[!duplicated(who_lkp_rev)]
rm(who_lkp_rev1, who_lkp_rev2)

who_lkp_rev <- c(who_lkp_rev,
                 c("A10BJ01" = "ITCA",
                   "A10BH" = "omarigliptin"))
who_lkp_mace <- names(who_lkp_rev)
names(who_lkp_mace) <- who_lkp_rev

