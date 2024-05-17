## 00_functions

## Can calcualte combined standard deviation of groups as per
## https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formulae_for_combining_groups.htm
## However, and while written code to allow this, it is probably unnecessary
## as these are randomised groups
## instead just assume same SD across groups and pool
## Note comb mean is the same for both

## Originally taken from above link. Latest link and text are here https://training.cochrane.org/handbook/current/chapter-06#section-6-5-2-10

CombMean <- function(n1, m1, n2, m2){
  top <- n1*m1 + n2*m2
  top/(n1+n2)
}

CombSd <- function(n1, m1, n2, m2, s1, s2){
  top = (n1-1)*s1^2 + (n2-1)*s2^2 + (n1*n2/(n1+n2)) * (m1^2 + m2^2 - 2*m1*m2)
  bottom = n1 + n2 -1
  res = (top/bottom)^0.5
  res
}


## test this, it is correct
function() {
a <- rnorm(10, sample(1:100, 1), sample(5:15, 1))
b <- rnorm(10, sample(1:100, 1), sample(5:15, 1))
ab <- c(a, b)
means <- map_dbl(list(a, b, ab), mean)
sds <- map_dbl(list(a, b, ab), sd)
mean_res <- CombMean(n1 = 10, m1 = means[1], n2 = 10, m2 = means[2])
sd_res <- CombSd(n1 = 10, m1 = means[1], n2 = 10, m2 = means[2], s1 = sds[1], s2 = sds[2])
}
# means
# sds
# mean_res
# sd_res
# 
# print(round(means[3],5) == round(mean_res,5))
# print(round(sds[3],5) == round(sd_res,5))

## Expand so that it is vectorised

CombSdVectorised <- function(n, m, s){
  myrows <- length(n)
  myrows_now <- myrows
  while(myrows_now >= 2){
    ## select first two values
    n1 <- n[1]
    n2 <- n[2]
    s1 <- s[1]
    s2 <- s[2]
    m1 <- m[1]
    m2 <- m[2]
    ## replace first value with combination of first two values and drop second value
    new_s1 <- CombSd(n1, m1, n2, m2, s1, s2)
    new_m1 <- CombMean(n1, m1, n2, m2)
    new_n1 <- n1 + n2
    s[1] <- new_s1
    m[1] <- new_m1
    n[1] <- new_n1
    s <- s[-2]
    m <- m[-2]
    n <- n[-2]
    # print(s)
    # print(m)
    ## recalculate the length
    myrows_now <- length(s)
  }
  s
}


## test this, it is correct
function() {
myvects <- map(rep(1, 10), ~ rnorm(sample(2:10, 1), sample(1:100, .x), sample(5:15, .x)))
mymeans <- map_dbl(myvects, mean)
mysds <- map_dbl(myvects, sd)
myns <- map_dbl(myvects, length)

a <- CombSdVectorised(n = myns, m = mymeans, s = mysds)
sd(unlist(myvects))
a}

##Vectorised and simple pooled means and sds, note that in the mean case it is 
## just a weighted means
PoolSD <- function(sds, ns){
  df <- ns-1
  sqrt( sum(sds^2 * df) / sum(df) )  
}

