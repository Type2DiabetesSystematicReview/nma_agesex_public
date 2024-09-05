# 17_adverse_events
library(tidyverse)

## functions
CnvrtCorrMatrix <- function(a){
  ## recovery whole matrix by duplication
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

ae6115_nct_id <- read_csv("Data/agesexhba1cmaceupdate_6115/coef.csv")$nct_id %>% unique()
ae9492_cf_nct_id <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv")$nct_id %>% unique()
gsk_cf_nct_id <- read_csv("Data/gsk_ae/age_sex_model_coefs_ae.csv")$nct_id %>% unique()

ae6115_cf <- read_csv("Data/agesexhba1cmaceupdate_6115/coef.csv")
ae9492_cf <- read_csv("Data/agesexhba1cmaceupdate_9492/coef.csv") %>% 
  filter(!nct_id %in% ae6115_nct_id)
gsk_cf <- read_csv("Data/gsk_ae/age_sex_model_coefs_ae.csv")%>% 
  filter(!nct_id %in% c(ae6115_nct_id, ae9492_cf_nct_id))

ae6115_vcv <- read_csv("Data/agesexhba1cmaceupdate_6115/vcov.csv") %>% 
  semi_join(ae6115_cf)
ae9492_vcv <- read_csv("Data/agesexhba1cmaceupdate_9492/vcov.csv")%>% 
  semi_join(ae9492_cf)
gsk_vcv <- read_csv("Data/gsk_ae/age_sex_model_vcov_ae.csv") %>% 
  semi_join(gsk_cf)

## bind together ae ones
aecf <- bind_rows(ae6115_cf %>% 
                    filter(outcome_method == "incident_event"),
                  ae9492_cf %>% 
                    filter(outcome_method == "incident_event"),
                  gsk_cf)
aevcv <- bind_rows(ae6115_vcv <- ae6115_vcv %>% 
                     filter(outcome_method == "incident_event"),
                   ae9492_vcv %>% 
                     filter(outcome_method == "incident_event"),
                   gsk_vcv)
rm(ae6115_cf, ae9492_cf, gsk_cf, ae6115_vcv, ae9492_vcv, gsk_vcv)

## drop excluded, 5 trials dropped----
exclusions <- read_csv("Data/exclusions_update.csv")
aecf %>% 
  inner_join(exclusions %>% filter(exclude ==1L) %>% rename(nct_id = trial_id)) %>% 
  count(nct_id, exclusion_reason)
aecf <- aecf %>% 
  anti_join(exclusions %>% filter(exclude ==1L) %>% select(nct_id = trial_id))

## add arm_meta and drop arm_meta where not in aecf
arm_meta <- read_csv("Scratch_data/arm_labels_hba1c.csv") 
## arm_id_subgroup has mismatched ordering between arm meta and agg so rename
arm_meta <- arm_meta %>% 
  rename(arm_id_subgroup_ordering = arm_id_subgroup)
arm_meta <- arm_meta %>% 
  semi_join(aecf)
aecf %>% 
  anti_join(arm_meta)
arm_meta %>% 
  filter(nct_id %in% c("NCT00763451", "NCT00712673"))
aecf %>% 
  filter(nct_id %in% c("NCT00763451", "NCT00712673")) %>% 
  distinct(nct_id, term) %>% 
  filter(str_detect(term, "^arm"), !str_detect(term, "sex"))

## drop where model failed, Apparent from VERY high standard errors
aecf <- aecf %>% 
  mutate(fail = std.error > 100) %>% 
  group_by(nct_id, models, outcome) %>% 
  mutate(fail = any(fail)) %>% 
  ungroup()
aecf_fail <- aecf %>% 
  filter(fail) %>% 
  distinct(nct_id, models, outcome)
aecf_fail %>% 
  mutate(v = 1L) %>% 
  pivot_wider(names_from = c(models, outcome), values_from = v)
aecf <- aecf %>% 
  filter(!fail)

## map terms across to arm_id_unqs for the 28 new trials -----
maketerms <- read_csv("Data/cleaned_data/Data/assign_armids_newtrials.csv")
aecf <- aecf %>% 
  left_join(maketerms %>% rename(term = old_term)) %>% 
  mutate(term = if_else(!is.na(new_term), new_term, term)) %>% 
  select(-new_term)
aevcv <- aevcv %>% 
  left_join(maketerms %>% rename(row = old_term)) %>% 
  mutate(row = if_else(!is.na(new_term), new_term, row)) %>% 
  select(-new_term) %>% 
  left_join(maketerms %>% rename(col = old_term)) %>% 
  mutate(col = if_else(!is.na(new_term), new_term, col)) %>% 
  select(-new_term)
rm(maketerms)

# deal with NCT01778049 population a and population b
aecf <- aecf %>%
  mutate(nct_id2 = case_when(
    nct_id == "NCT01778049" & str_detect(term, "ipd00001|ipd00002") ~ paste0(nct_id, "_a"),
    nct_id == "NCT01778049" & str_detect(term, "ipd00003|ipd00004") ~ paste0(nct_id, "_b"),
    TRUE ~ nct_id)
  )
aevcv <- aevcv %>%
  mutate(nct_id2 = case_when(
    nct_id == "NCT01778049" & 
      (str_detect(row, "ipd00001|ipd00002") | str_detect(col, "ipd00001|ipd00002")) ~ paste0(nct_id, "_a"),
    nct_id == "NCT01778049" & 
      (str_detect(row, "ipd00003|ipd00004") | str_detect(col, "ipd00003|ipd00004")) ~ paste0(nct_id, "_b"),
    TRUE ~ nct_id)
  )

## select only variables of interest
aecf <- aecf %>% 
  filter(str_detect(term, "arm_f"), str_detect(term, "sex") | str_detect(term, "age")) %>% 
  nest(.by = c(outcome, nct_id, nct_id2, models)) %>% 
  rename(cf = data)
aevcv <- aevcv %>% 
  nest(.by = c(outcome, nct_id, nct_id2, models)) %>% 
  rename(vcv = data)
ae <- aecf %>% 
  inner_join(aevcv)
rm(aecf, aevcv)
ae$vcv <- map2(ae$cf, ae$vcv, ~ .y %>% 
                 filter(row %in% .x$term,
                        col %in% .x$term))
ae$vcv  <- map(ae$vcv, ~ CnvrtCorrMatrix(.x %>% select(row, col, r)))

## convert correlation to variance-covariance mtrix
ae$crl <- ae$vcv
ae$vcv <- NULL
ae$vcv <- map2(ae$cf, ae$crl, function(cf, crl){
  crl * outer(cf$std.error, cf$std.error)
})
## all are positive definite
ae$crl_pd <- map_lgl(ae$crl, matrixcalc::is.positive.definite)
ae$vcv_pd <- map_lgl(ae$vcv, matrixcalc::is.positive.definite)

## set-up for fitting regression models ----
fromhba1c <- readRDS("Scratch_data/simulated_ipd.Rds")

## check same with correlation/covariance data and append arms
all(ae$nct_id2 %in% fromhba1c$nct_id2)
fromhba1c <- fromhba1c %>% 
  semi_join(ae %>% select(nct_id2))
fromhba1c <- fromhba1c %>% 
  distinct(nct_id, nct_id2, arm_id_unq, arm_f, reference_arm)
fromhba1c <- fromhba1c %>% 
  group_by(nct_id, nct_id2) %>% 
  nest() %>% 
  ungroup()
ae <- ae %>% 
  inner_join(fromhba1c %>% rename(arm_lkp = data))

reg_ae <- pmap(
  list(ae$nct_id,
       ae$nct_id2,
       ae$outcome,
       ae$models,
       ae$cf,
       ae$crl,
       ae$vcv,
       ae$arm_lkp,
       seq_along(ae$nct_id)),
  function(nct_id, nct_id2, outcome, models, cfs, crl, vcv, lkp, i) {
    ## arrange covariance and correlation matrix to match terms
    # if(i == 181) browser()
    # browser()
    if(!all(colnames(crl) %in% cfs$term )) stop("coefficient terms and correlation matrix do not match")
    if(!all(colnames(vcv) %in% cfs$term )) stop("coefficient terms and correlation matrix do not match")
    crl <- crl[cfs$term,cfs$term]
    vcv <- vcv[cfs$term,cfs$term]
    
    ## add NA row for no treatments (note the correlation for this is set to NULL)
    cfs <- cfs %>% 
      mutate(trt = NA_character_)
    for (i in lkp$arm_id_unq) {
      cfs <- cfs %>% 
        mutate(trt = case_when(
          !is.na(trt) ~ trt,
          str_detect(term, i) ~ i,
          TRUE ~ trt))
    }
    addref <- setdiff(lkp$arm_id_unq, cfs$trt)
    cfs <- bind_rows(tibble(trt = addref),
                     cfs)
    cfs <- cfs %>% 
      left_join(lkp %>% rename(trt = arm_id_unq), by = join_by(trt))
    
    ## add covariate columns
    ## Note than manual suggests listing as NA but example data has zeros and ones.
    cfs <- cfs %>% 
      mutate(sex = case_when(
        str_detect(term, "sex") ~ TRUE,
        reference_arm ==1 ~ FALSE,
        TRUE ~ NA))
    
    ## add covariate columns
    cfs <- cfs %>% 
      mutate(age10 = case_when(
        str_detect(term, "age10") ~ 10,
        reference_arm ==1 ~ 0,
        TRUE ~ NA_real_))
    
    tibble(nct_id = nct_id,
           nct_id2 = nct_id2,
           outcome = outcome,
           models = models,
           cfs = list(cfs),
           crl = list(crl),
           vcv = list(vcv))
  })


reg_ae <- bind_rows(reg_ae)
reg_ae %>% 
  filter(nct_id %in% c("NCT00763451", "NCT00712673"), outcome == "ae_sae", models == "pois") %>% 
  unnest(cfs)
reg_ae <- reg_ae %>% 
  mutate(models = case_when(
    models %in% c("pois", "glm") ~ "pois",
    models %in% c("nb", "glmnb")  ~ "nb"))
## add in arm labels
## address two trials with different regimes. pick more aggressive regime
arm_meta <- arm_meta %>% 
  mutate(arm_id_unq = case_when(
    arm_id_unq %in% c("uaa10323", "uaa10349") ~ paste0(arm_id_unq, "b"),
    TRUE ~ arm_id_unq)
  )
arm_meta <- arm_meta %>% 
  nest(.by = nct_id) %>% rename(meta = data)
reg_ae <- reg_ae %>% 
  inner_join(arm_meta)
reg_ae$cfs <- map2(reg_ae$cfs, reg_ae$meta, function(cfs, meta) {
 cfs %>% 
    inner_join(meta %>% 
                 select(trt = arm_id_unq, arm_lvl, trtcls5, trtcls4) %>% 
                 distinct(),
               by = join_by(trt))})
reg_ae$crl <- map2(reg_ae$crl, reg_ae$cfs, ~ .x[.y$term[-1],
                                                .y$term[-1]])
reg_ae$vcv <- map2(reg_ae$vcv, reg_ae$cfs, ~ .x[.y$term[-1],
                                                .y$term[-1]])
reg_ae_nst <- reg_ae %>% 
  nest(.by = c(outcome, models))

reg_ae_nst$crl <- map(reg_ae_nst$data, function(mydf) {
  crl <- mydf$crl
  names(crl) <- mydf$nct_id
  crl
})
reg_ae_nst$vcv <- map(reg_ae_nst$data, function(mydf) {
  vcv <- mydf$vcv
  names(vcv) <- mydf$nct_id
  vcv
})
reg_ae_nst$data <- map(reg_ae_nst$data, ~ .x %>% 
                         select(nct_id, nct_id2, cfs) %>% 
                         unnest(cfs))
reg_ae_nst <- reg_ae_nst %>% 
  rename(cfs = data)

## trials to drop due to single arms or disconnected networks ----
drp <- read_csv("outcome,models,arm_drp,trl_drp
ae_sae,pois,A10BH03_d5,
ae_gi,pois,A10BH03_d5,NCT00688701
ae_hypog,pois,,
ae_uti,pois,A10BJ05_d1.5,
ae_sae,nb,A10BH03_d5,
ae_gi,nb,A10BH03_d5,NCT00688701
ae_hypog,nb,,
ae_uti,nb,A10BK01_d10;A10BK03_d25,")

drp$arm_drp <- str_split(drp$arm_drp, pattern = ";")
reg_ae_nst <- reg_ae_nst %>% 
  inner_join(drp)
reg_ae_nst$trl_arm_drp <- map2(reg_ae_nst$cfs, reg_ae_nst$arm_drp, ~ .x %>% 
                             filter(arm_lvl %in% .y) %>% 
                             distinct(nct_id) %>% 
                             pull())
reg_ae_nst$trl_drp <- map2(reg_ae_nst$trl_arm_drp, reg_ae_nst$trl_drp, union)
reg_ae_nst$trl_arm_drp <- NULL
reg_ae_nst$cfs <- map2(reg_ae_nst$cfs, reg_ae_nst$trl_drp, ~ .x %>% 
                     filter(!nct_id %in% .y))
reg_ae_nst$crl <- map2(reg_ae_nst$crl, reg_ae_nst$trl_drp, ~ .x[setdiff(names(.x), .y)])
reg_ae_nst$vcv <- map2(reg_ae_nst$vcv, reg_ae_nst$trl_drp, ~ .x[setdiff(names(.x), .y)])
reg_ae_nst <- reg_ae_nst %>% 
  select(-arm_drp, -trl_drp)
saveRDS(reg_ae_nst, "Scratch_data/for_ae_regression.Rds")
