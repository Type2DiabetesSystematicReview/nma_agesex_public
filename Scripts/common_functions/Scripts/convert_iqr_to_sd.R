
## Functions to get mean and sd from different data summaries ----
## Use published method to estimate the mean and sd from each set of statistics, note that none of these have other information (eg also categories)
## doi = {10.1186/1471-2288-14-135}
# see excel spreadhseet in supplement, using summary 1.
# =($D5-$B5)/(2*NORM.INV(($E5-0.375)/($E5+0.25),0,1))
# calculate median from mean re-writing equation from scenario 1
EstMed <- function(mu, minim, maxim){
  (4*mu - minim - maxim)/2
}
EstMean <- function(med, minim, maxim){
  (minim + 2*med + maxim)/4
}
EstSD <- function(minim, maxim, smpl){
  (maxim - minim)/ (2*qnorm((smpl-0.375)/ (smpl+0.25)))
}
EstSD2 <- function(first, third, smpl){
  numer <- (third - first)
  denom <- 2*qnorm((0.75*smpl-0.125)/ (smpl+0.25))
  numer/denom
}
EstMean2 <- function(first, medn, third){
  (first + medn + third)/3
}

