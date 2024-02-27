library(tidyverse)
a <- expand_grid(data_lvl = c("aggipd", "ipd"),
                 fe_re = c("fixed", "random"),
                 nwork = c("mono", "dual", "triple")) 
a <- a %>% 
  mutate(model_n = seq_along(a$data_lvl)) %>% 
  select(model_n, everything())
write_csv(a, file = "Outputs/model_order.csv")


dropdisconnect <- c("NCT02477865", "JapicCTI-101351",
                    "NCT02477969", "JapicCTI-101352", "NCT03508323",
                    "UMIN000007051")
exclude <- dropdisconnect
PasteAnd <- function(x) {
  if(length(x) == 1) return(x)
  x1 <- x[-length(x)]
  x2 <- x[length(x)]
  x1 <- paste(x1, sep = ", ")
  paste0(x1, " and ", x2)
}
exclude <- tibble(reason = "Disconnected from network.",
                  trials = length(exclude),
                  nct_ids = exclude %>% PasteAnd(),
                  level = "Aggregate")
write_tsv(exclude, "Outputs/Trial_exclusion_during_cleaning.txt", append = TRUE)