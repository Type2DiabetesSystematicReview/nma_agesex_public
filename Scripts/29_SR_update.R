#29_SR_update
library(tidyverse)

a <- read_csv("Data/update2025/fulltext_screen.csv")
a %>% 
  count(decision)
a <- a %>% 
  filter(decision == "INCLUDED")

## separate out trial IDs
a$trial_id <- str_split(a$trial_id, ";")
a <- a %>% 
  unnest(trial_id)
a <- a %>% 
  mutate(trial_id = str_trim(trial_id))

## get enrolled numbers
a <- a %>%
  # first, pull off the prefix (case‐insensitive)
  mutate(
    prefix = str_to_lower(trial_id) %>%
      str_extract("^(nct|chictr|ctri|eu|irct|jrct|umin)"),
    registry = case_match(
      prefix,
      "nct"    ~ "nct_id",
      "chictr" ~ "chictr",
      "ctri"   ~ "ctri",
      "eu"     ~ "eu",
      "irct"   ~ "icrt",
      "jrct"   ~ "jrct",
      "umin"   ~ "umin",
      .default = NA_character_
    )
  ) %>%
  select(-prefix)

## Add AACT CTTI database password. Will create a popout window
mypassword <- rstudioapi::askForPassword()

## Connect to database
library(RPostgreSQL)
drv <- dbDriver('PostgreSQL')
con <- dbConnect(drv, dbname="aact",host="aact-db.ctti-clinicaltrials.org",
                 port=5432, user="dmcalli",
                 password=mypassword)
# 127 unique strings
nct_id_srch <- a$trial_id[a$registry == "nct_id"] %>% 
  unique() %>% 
  paste(collapse = "','")
myq <- paste0("SELECT * FROM STUDIES where nct_id IN ('",
              nct_id_srch,
              "')")
nct_id_q <- dbGetQuery(con, myq)
saveRDS(nct_id_q, "Scratch_data/new_trial_aact.Rds")

umin_srch <- a$trial_id[a$registry == "umin"] %>% 
  unique() %>% 
  paste(collapse = " OR ")
# umin <- read_csv("https://upload.umin.ac.jp/ctr_csv/ctr_data_j.csv.gz")
# saveRDS(umin, "Scratch_data/umin.Rds")
umin <- readRDS("Scratch_data/umin.Rds")
umin2 <- umin %>% 
  select(2, 105) %>% 
  slice(-1)
names(umin2) <- c("trial_id", "enrollment")
umin_got <- umin2 %>% 
  semi_join(a)
uminb <- umin
names(uminb)[2] <- "trial_id"
uminb <- uminb %>% 
  semi_join(a)
for_aa <- bind_rows(nct_id_q %>% 
                       select(trial_id = nct_id, enrollment),
                    umin_got %>% mutate(enrollment = as.integer(enrollment)))
msng <- a %>% 
  anti_join(for_aa %>% filter(!is.na(enrollment)))
ictrp <- read_csv("Data/update2025/ictrp_srch_lkp.csv")
msng2 <- msng %>% 
  inner_join(ictrp)
for_aa2 <- bind_rows(for_aa,
                    msng2 %>% select(trial_id, enrollment = enrolled))
write_csv(for_aa2, "Scratch_data/enrolled_ad_az.csv")

msng_new <- a %>% 
  anti_join(for_aa2 %>% filter(!is.na(enrollment)))
msng_new %>% 
  count(trial_id, sort = TRUE)
## Compare to previous trials
xcld <- read_csv("Data/exclusions_update.csv")
a_dist <- a %>% 
  distinct(trial_id)
a_in <-  %>% 
  inner_join(xcld)
## Of 203, 58 already got. 
a_got <- a_in %>% 
  filter(exclude ==0) 
## remaining 21 in data 
a_in <- a_in %>% 
  filter(exclude ==1)
## of these 4 should be excluded. Other 17 are "no results reported"
a_in %>% 
  count(exclusion_reason)
# 124 plus 17 gives 141 new trials/results from update
a_not <- a_dist %>% 
  anti_join(xcld)
# Of these 28 have an unreported ID and can be excluded
a %>% 
  semi_join(a_not) %>% 
  count(trial_id == "Unreported ID")


