# rn_nma_mars.R
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Set up plaque psoriasis network combining IPD and AgD
library(dplyr)
library(stringr)
library(tibble)
library(multinma)

## functions
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

## data
tot <- readRDS("data_for_mars.Rds")
tot <- tot %>% 
  filter(drug_regime_smpl == args[1])

ipd <- tot$ipd[[1]]
agg <- tot$agg[[1]]
rm(tot)

## print trial numbers
print(nrow(ipd))
print(nrow(agg))

# set-up network
dm_net <- RptNetwork(ipd, agg)

dm_net <- add_integration(dm_net,
                          age = distr(qnorm, mean = age_m, sd = age_sd),
                          sex = distr(qbinom, size = 1, 
                                      prob = male/n),
                          base = distr(qnorm, mean = 0, sd = 0),
                          n_int = 64)

x <- nma(dm_net, 
        regression = ~ base + age*.trt + sex*.trt,
        trt_effects = args[[3]], 
        link = "identity", 
        likelihood = "normal",
        class_interactions = args[[2]],
        prior_intercept = normal(scale = 20),
        prior_trt = normal(scale = 10),
        prior_reg = normal(scale = 10),
        QR = TRUE, 
        cores = 4)


# save outputs
datetime <- Sys.time() %>% 
   as.character() %>% 
   str_replace_all("\\-", "") %>% 
   str_replace_all("\\:", "") %>% 
   str_replace_all("\\.", "_") %>% 
   str_replace_all("\\s{1,}", "_") 
 
mypath <- paste0(args[[1]], "_", args[[2]], "_", args[[3]], "_", datetime)
print(mypath)

if (!dir.exists(mypath)) dir.create(mypath)

MySave <- function(objname, filename){
   saveRDS(object = objname,
           file = paste0(mypath, "/", filename))
 }
 
 nms <- names(x$stanfit)
 pars <- nms[str_detect(nms, "^beta|^d|^delta|^fitted_agd|^lp|^mu|^tau")]
 
 a <- dic(x)
 MySave(a[c("dic", "pd", "pv", "resdev")], "dic.Rds")
 rm(a)
 ## obtain samples and other summaries for parameters of interest
 a <- summary(x,
              pars = pars,
              probs = c(0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975))
 MySave(a, "summary.Rds")
 a <- relative_effects(x,
                       probs = c(0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975),
                       predictive_distribution = TRUE)
 MySave(a, "relative_effects.Rds")
 a <- posterior_ranks(x)
 MySave(a, "posterior_ranks.Rds")
 a <- posterior_rank_probs(x)
 MySave(a, "posterior_rank_probs.Rds")
 a <- plot(dic(x))
 pdf(paste0(mypath, "/dicplots.pdf"))
 a
 dev.off()
 