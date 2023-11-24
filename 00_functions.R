# 00_functions


SimplifyDrugs <- function(re_arm) {
  ## drop any where same drug in every arm
  # expects nct_id, arm_id_unq and drug_code
  ## function is not very modular. Consider reorganising into smaller parts
  re_arm_drp <- re_arm %>% 
    group_by(nct_id) %>% 
    nest() %>% 
    ungroup()
  IdentifyCrossArms <- function(arm_drug){
    x <- arm_drug %>% 
      select(arm_id_unq, drug_code) %>%
      mutate(v = 1L) %>% 
      spread(drug_code, v, fill = 0L) %>% 
      select(-arm_id_unq) %>% 
      summarise_all(all)
    x <- unlist(x)
    names(x[x])  
  }
  re_arm_drp$sameacross <- map(re_arm_drp$data, IdentifyCrossArms)
  re_arm_drp$n <- map_int(re_arm_drp$sameacross, length)
  re_arm_drp <- re_arm_drp %>% 
    filter(n >0) %>% 
    unnest(sameacross)
  ## identify implicit controls and where present add a placebo
  re_arm_drp$impcntrl <- map(re_arm_drp$data, function(a) {
    a %>% 
      mutate(v = 1L) %>% 
      spread(drug_code, v, fill = 0L) %>% 
      gather(key = "drug_code", cntrl, -arm_id_unq) %>% 
      mutate(drug_code = if_else(cntrl == 0L, 
                                 "implicit_control",
                                 drug_code))
  })
  re_arm_drp$sameacross2 <- map(re_arm_drp$impcntrl, IdentifyCrossArms)
  re_arm_drp$retain <- map2(re_arm_drp$impcntrl, re_arm_drp$sameacross2, ~ {(  
    .x %>% 
      filter(!drug_code %in% .y))
  })
  re_arm_drp$retained <- map_int(re_arm_drp$retain, nrow) 
  re_arm_keep <- re_arm_drp %>% 
    filter(!retained ==0)
  re_arm_drp <- re_arm_drp %>% 
    filter(retained ==0)
  re_arm_drp <- re_arm_drp %>% 
    select(nct_id, sameacross, data) %>% 
    unnest(data) 
  ## 5 trials with implicit controls
  re_arm_keep <- re_arm_keep %>% 
    select(nct_id, sameacross, retain) %>% 
    unnest(retain) %>% 
    select(-cntrl) 
  
  re_arm2 <- bind_rows(re_arm %>% 
                         filter(!nct_id %in% re_arm_keep$nct_id,
                                !nct_id %in% re_arm_drp$nct_id) %>% 
                         mutate(sameacross = "nil here see ancillary"),
                       re_arm_keep)
  list(keep = re_arm2, drop = re_arm_drp)
}

RptNetwork <- function(ipd_choose, agg_choose){
  combine_network(
    set_ipd(ipd_choose,
            study = nct_id2,
            trt = drug_code,
            y = result,
            trt_class = trtcls5),
    set_agd_arm(agg_choose %>% 
                  filter(treat_or_ref == "arm_level_outcome"),
                study = nct_id,
                trt = drug_code,
                y = result, 
                se = se,
                trt_class = trtcls5,
                sample_size = n),
    set_agd_contrast(agg_choose %>% 
                       filter(!treat_or_ref == "arm_level_outcome"),
                     study = nct_id,
                     trt = drug_code,
                     y = result, 
                     se = se,
                     trt_class = trtcls5,
                     sample_size = n))  
}

RptNetworkClass <- function(ipd_choose, agg_choose){
  combine_network(
    set_ipd(ipd_choose,
            study = nct_id2,
            trt = trtcls5,
            y = result,
            trt_class = trtcls4),
    set_agd_arm(agg_choose %>% 
                  filter(treat_or_ref == "arm_level_outcome"),
                study = nct_id,
                trt = trtcls5,
                y = result, 
                se = se,
                trt_class = trtcls4,
                sample_size = n),
    set_agd_contrast(agg_choose %>% 
                       filter(!treat_or_ref == "arm_level_outcome"),
                     study = nct_id,
                     trt = trtcls5,
                     y = result, 
                     se = se,
                     trt_class = trtcls4,
                     sample_size = n))  
}