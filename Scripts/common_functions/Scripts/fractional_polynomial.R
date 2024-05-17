MakePowers <- function(value, power1, power2) {
  value1 = if_else(power1 ==0, log(value), value^power1)
  value2 = case_when(
    is.na(power2) ~ 0,
    power2 ==0 ~ log(value),
    TRUE ~ value^power2)
  value2 = if_else(power1 == power2, value2 * log(value), value2)
  cbind(value1, value2)
}
MP1 <- function(value, power1) {
  if(power1 ==0) log(value) else value^power1
}
MP2 <- function(value, power1, power2) {
  if(is.na(power2)) value2 <- rep(0, length(value))
  if(!is.na(power2)){
    if(power2 ==0) value2 <- log(value)
    if(power2 !=0) value2 <- value^power2
    if(power1 == power2) value2 <- value2 * log(value)
  }
  value2
}