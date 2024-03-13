## Find external dependencies - data and scripts
## Ran towards end of analysis so that could identify all data and scripts used by this project and make it self contained.

library(tidyverse)
names_a <- list.files("Scripts/", full.names = FALSE, pattern = "^[0-8]")
a <- list.files("Scripts/", full.names = TRUE, pattern = "^[0-8]")
names(a) <- names_a
rm(names_a)
## drop the dependency renamining file
a <- a[setdiff(names(a), c("00__pull_external_files_folders.R", "Scripts/00__identify_external_data_scripts.R"))]

rv <- tibble(filename = names(a), content = a)
rv$content <- map(rv$content, ~ {
  x <- read_lines(.x)
  x[str_detect(x, "\\.{2,2}/")]
})
rv <- rv %>% 
  unnest(cols = content) %>% 
  distinct()
rv <- rv %>% 
  mutate(list.files = str_detect(content, "list\\.files"))
## Now no folders. do not delete this as needs to be able to run again if create new external dependencies
fldrs <- rv %>% 
  filter(list.files)
rv <- rv %>% 
  filter(!list.files) %>% 
  select(-list.files)
rv <- rv %>% 
  separate(content, into = c("pre", "within", "post"), sep = '\\"', remove = FALSE)
rv %>% 
  count(within)

## Drop WHO ones
who <- rv %>% 
  filter(str_detect(within, "2018 ATC index with DDDs.xlsx"))
rv <- rv %>% 
  filter(!str_detect(within, "2018 ATC index with DDDs.xlsx"))

## create rename files for cleaned data, create folders for cleaned_data
rv <- rv %>% 
  separate(within, c("moveup","top_level", "subsequent"), sep = "/", extra = "merge", remove = FALSE)
# write_csv(rv, "Data/table_for_copying_files_into_directory.csv")

         