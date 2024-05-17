MakeDist <- function(m, s, l, u){
  freqs <- msm::ptnorm(seq(0, 100, 1),
                       mean = m, sd = s,
                       lower = l, upper = u)
  freqs <- freqs - dplyr::lag(freqs, default = 0)
  freqs <- freqs/sum(freqs)
  freqs
}

# Functions for truncted normal distributions written by asking copilot to re-write wikipedia equations as R code
# Checked and give correct results
truncated_normal_mean <- function(mu, sigma, a, b) {
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  Z <- pnorm(beta) - pnorm(alpha)
  return(mu + (dnorm(alpha) - dnorm(beta)) * sigma / Z)
}

truncated_normal_variance <- function(mu, sigma, a, b) {
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  Z <- pnorm(beta) - pnorm(alpha)
  variance_value <- sigma^2 * (1 - (beta * dnorm(beta) - alpha * dnorm(alpha)) / Z - ((dnorm(alpha) - dnorm(beta)) / Z)^2)
  return(variance_value)
}


EstimateMuDispAge <- function (x) {
  # Loop through trials
  if(str_sub(x$nct_id, -1) == 1) print(x$nct_id)
  
  lower <- x$min_age
  upper <- x$max_age
  trial_mean <- x$age_m
  trial_sd <- x$age_sd
  trial_var <- x$age_sd^2
  
  ## Create grid
  mu_x <- seq(lower, upper, 0.5)
  sd_x <- seq(1, upper-lower, 0.5)
  full_grid <- expand.grid(mu_x = mu_x, sd_x = sd_x)
  
  ## Calculate for all values of grid, is vectorised so is fast, is faster than one in truncnorm package
  full_grid$mean_x <- truncated_normal_mean(full_grid$mu_x, full_grid$sd_x, lower, upper)
  full_grid$var_x <- truncated_normal_variance(full_grid$mu_x, full_grid$sd_x, lower, upper)
  
  # print(nrow(full_grid))
  # browser()  
  ## Identify closest values
  full_grid <- full_grid %>%
    as_tibble() %>%
    mutate(mean_diff = abs(trial_mean - mean_x),
           var_diff = abs(trial_var - var_x),
           total_diff = mean_diff + var_diff) %>%
    arrange(total_diff, mean_diff, var_diff)
  ## Append original parameters
  estimate <- full_grid %>%
    slice(1:10) %>%
    mutate(trial_mean = trial_mean,
           trial_var = trial_var,
           trial_lower = lower,
           trial_upper = upper,
           trial_sd = trial_sd) %>%
    select(trial_mean, mean_x, trial_var, var_x, mu_x, sd_x, trial_sd, everything())
  estimate 
}

MuSigmaCalc <- function(x_min, x_max, mean, sd, mymethod) {
  ## Supply lower, upper, mean and sd
  ## function obtains values of mu and sigma which give
  ## values of mean and sd closest to the observed values
  ## mu and sigma are constrained (0-100) and 1-50 respectively
  # browser()
  ToMin <- function(mu_sigma, ...) {
    abs(mean - truncated_normal_mean(mu = mu_sigma[1], sigma = mu_sigma[2], a = x_min, b = x_max)) +
      abs(sd - truncated_normal_variance(mu = mu_sigma[1], sigma = mu_sigma[2], a = x_min, b = x_max)^0.5)
  }
  if (mymethod == "L-BFGS-B") {
    optim(par = c(mean, sd), fn = ToMin, method = mymethod, lower = c(0, 0.1*sd), upper = c(mean + 3*sd, 5*sd))
  } else {
    optim(par = c(mean, sd), fn = ToMin, method = mymethod)
    
  }
}

