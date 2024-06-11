README

The following export request contains regression models for the effect of each treatment on MACE adjusting for age and sex. We previously exported these summaries of regression models before but had not centred age. This led to some problems wiht subsequent analysis. This time we have scaled age before fitting the models (age-60)/10

# PDF file

test_for_non_phy_cox.pdf - is  graphical test for non proportional hazards for each trial and model. Horizontal lines show now evidence for non-proportional hazards. The residuals were not plotted to prrevent possibility of identifying true event times.

# CSV files

age_sex_model_coefs.csv - regression model outputs. One row per model and term. Contains term, estimate, standard error, test statistic, p-value and lower and upper confidence intervals
age_sex_model_diags.csv - regression model diagnostics. One row per model contains loglikelihood, AIC, BIC, number of observations and number of event.
age_sex_model_vcov.csv - regression model outputs. Shows correlation between coefficients (not variables) for various model outputs. Obtained from vcov function in R.

# We also wish to provide very approximate summaries of the number of people aged over 80 in each trial and do so below.
In NCT02065791 age is in bands so this calculation is not included

- redacted gives a count of the number of participants in whom age is redacted.
- percent gives a very approximate percentage of the percentage aged over 80 years old

nct_id	redacted	percentgt80
NCT01032629	0	<0.5%
NCT01131676	1	1.7%
NCT01989754	0	<1%
NCT02465515	2	2.3%