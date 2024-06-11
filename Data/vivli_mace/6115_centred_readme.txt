00_readme

This export request contains the results of regression models of MACE events for a single trial - NCT00968708 - sponsor Takeda. We previously exported similar models (in September 2023) but this time have centred age (subtracting 60 and dividing by 10) as including it in the model uncentred was causing some problems with the variance-covariance matrix and standard errors during our subsequent analysis.

CSV files
 age_sex_model_coefs.csv - regression model output. One row per model and term. Contains term, estimate, standard error, test statistic, p-value and lower and upper confidence intervals.
 age_sex_model_diags.csv - regression model output. One row per model. Contains log likelihood, AIC, BIC, number of observations and number of events.
 age_sex_model_vcov.csv  - regression model output. Shows correlation between coefficients (not variables) for coefficients in model. Obtained from vcov function in R.

 