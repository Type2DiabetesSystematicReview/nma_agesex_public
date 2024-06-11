README

The following export request contains summary statistics and regression models for the effect of each treatment on MACE adjusting for age and sex.

# Figures

## The Tiff files are all diagnostic plots checking the distribution/fit of summary statsitics to the data.

age_normally_distributed.tiff - shows density plots of age by trial, arm, sex and event status. Compares actual density to noraml distribution with same mean and standard deviation
check_distr_simulated_real_event_times_single_trial.tiff - shows density plots of event times by age category, sex and treatment arm for a single trial. Compares actual event time density (having dividied each value by the follow-up time and performed a logit transformation) to a normal distribution with the same mean and standard deviation.
check_ecdf_fit_censor_time.tiff - For trial adn treatment arm compares actual cumulative distribution of censoring times to simulation based on summary statistics. Removed axes, rescaled all and randomised order to prevent possibility of re-identifying true censoring times
lm_eventtime_fit.tiff -  shows density plots of event times by sex and treatment arm for a single trial. Compares actual event time density (having dividied each value by the follow-up time and performed a logit transformation) to a normal distribution with the same mean and standard deviation. Mean estimated using linear regression models.

# PDF file

test_for_non_phy_cox.pdf - is  graphical test for non proportional hazards for each trial and model. Horizontal lines show now evidence for non-proportional hazards. The residuals were not plotted to prrevent possibility of identifying true event times.

# R files

19_mace_recover_simulation.R, 20_mace_cox_regression_pseudo_real.R, 21_mace_cox_regression_real.R - are some analysis scripts we wrote and wish to re-use some of the code from

# CSV files
One trial had collaped age into four categories. As few individuals were in the lowest age category I further collapsed the lower two categories into one.

age_sex_model_coefs.csv - regression model outputs. One row per model and term. Contains term, estimate, standard error, test statistic, p-value and lower and upper confidence intervals
age_sex_model_diags.csv - regression model diagnostics. One row per model contains loglikelihood, AIC, BIC, number of observations and number of event.
age_sex_model_diags_fixed.csv - regression model diagnostics for linear HbA1c models. Contains various model diagnostics produced by linear model in R.
age_sex_model_vcov.csv - regression model outputs. Shows correlation between coefficients (not variables) for various model outputs. Obtained from vcov function in R.
censoring_distribution.csv - Summary of censoring times. One row per trial, arm and sex. Gives mean age, standard deviation, number of participants, rank correlation between age and censoring time and selected quantiles for the distribution of censoring times.
censoring_distribution_single_trial.csv - as above, but for a trial where age was given as three categories rather than as a continuous variable
event_time_distribution_fp.csv - Summary of event times. One row per trial, arm and sex. Among those with an event gives number of participants, age (as mean and standard deviation), and the ouputs of a linear regression model predicting censoring times on the logit scale. THe power1 and power2 columns indicate transformations used for age, and the columns starting with "est_" and "r se_" indicate the model coefficients and standard errors respectively. r indicates the correlation between the coefficients (not variables)
event_time_distribution_single_trial.csv - For a single trial, mean and standard deviation of the event time among those who had an event (having transformed the event time by dividing it by the censoring time and perfoming a logit transformation). One row per arm, sex and age cagegory.
