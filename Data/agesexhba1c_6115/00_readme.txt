# Overview

This is an export request comprising summary statistics, trial-level regression model outputs and some reference files (ie metadata) based on an analysis of treatment heterogeneity by age and sex.

# Specific files

For the purpose of this description, these files are split into summary statistics/plots, model outputs and reference files.

# Summary statistics/plots

## age_distribution_baseline_categorical.csv

This file gives the mean and standard deviation for age according to trial, sex and treatment arm for six trials. Since for these trials age was provided in categories rather than as an integer, the mean and standard deviation for each trial was estimated using a function in R which fits a univariate distribution (in this case normal) to censored data using maximum likelihood.

### Rows - one row per trial/arm/sex combination
### columns

- nct_id - trial ID
- arm_id_unq - arm ID assigned by our team to uniquely identify each arm
- sex - male or female
- participants - number of participants in stratum
- mean - mean age in years (rounded to 1 DP)
- sd - standard deviation for age in years (rounded to 1 DP)

## age_distribution_baseline_continuous.csv

As "age_distribution_baseline_categorical.csv" for the 48 trials with age as continuous data, EXCEPT i) the mean and standard deviation were calculated directly from data and ii) selected quantiles from the empirical distribution were also obtained in order to describe the age distribution. The quantiles were selected using the Ramer Douglas Peucker algorithm for simplifying curves. This allowed us to produce a non-parametric summary of the age distribution.

### rows - as above
### columns - as above AND

- Number of quantiles - number of quantiles included from empirical cumulative distribution
- Age (years) at quantiles (%) - the age in years at selected quantiles (eg in the first row at 26th quantile the age in years was 53)
- epsilon - the threshold for filtering outliers from the simplified curve

## age_ecdf_simplification_diagnostics.tiff

This is a tiff file comparing the empirical cumulative distribution function (ECDF) to a simplified representation of the ECDF using fewer points. This was done to check whether the simplified curve (produced using the Ramer Douglas Peucker algorithm) is farily representative of original curve. Each pair of lines (red and black) correspond to single trial/arm/sex combination, So that the original data CANNOT be recovered from this plot (even with digitisation software) we randomised the plotting order for the facets (trials, arms sexes are all mixed-up), allowed the scales to to be free (so it impossible to tell the actual values and unlabelled) for any plot and saved as the image a low resolution (150 dpi) tiff.

## hba1c_base_change_overall.csv

This file summarises the unadjusted mean HBA1c at baseline and the end-point (mean) and the change from baseline (mean and SD). There is one row per trial and arm. This was produced to allow a check of the quality of our analysis against similar statistics in clinicaltrials.gov.

### Rows - One row per trial and arm. 
### Cols
-  nct_id - trial ID
-  nct_id2 - trial ID with a letter appended for one trial which had two populations within the same NCT_ID
- arm_id_unq - arm ID assigned by our team to uniquely identify each arm
- base - mean HbA1c at baseline - rounded to 2 DP
- end - mean HbA1c at end - rounded to 2 DP
- result - mean change in HbA1c rounded to 2 DP
- n - number at baseline
- n_analysed - number with two or more HbA1c measurements
- sd - standard deviation for change in HbA1c

# Model outputs

The three model output files - all starting "age_sex_model... .csv" - have results (coefficients, diagnostics and variance-covariance matrix for the same set of models).

We first fit a model of hba1c at baseline on age, sex and treatment arm. We then fit the following ANCOVA-type models.
[1] final_hba1c ~ baseline_hba1c + treatment_arm
[2] final_hba1c ~ baseline_hba1c + age * treatment_arm                                                
[3] final_hba1c ~ baseline_hba1c + sex * treatment_arm                                                  
[4] final_hba1c ~ baseline_hba1c + age * treatment_arm + sex * treatment_arm                                  
[5] final_hba1c ~ baseline_hba1c + treatment_arm + baseline_hba1c:treatment_arm                                        
[6] final_hba1c ~ baseline_hba1c + age + treatment_arm + age:treatment_arm + baseline_hba1c:treatment_arm                
[7] final_hba1c ~ baseline_hba1c + sex + treatment_arm + sex:treatment_arm + baseline_hba1c:treatment_arm                 
[8] final_hba1c ~ baseline_hba1c + age + treatment_arm + sex + age:treatment_arm + treatment_arm:sex + baseline_hba1c:treatment_arm

The first is a standard ancova. The fifth is an anvova with a baseline:treatment arm interaction. The other six models aadd in interactions for age and sex. The following summaries were produced for all nine models (1 baseline HbA1c 8 ANCOVA).

## age_sex_model_coefs.csv

Coefficeints, standard errors, t.statistics and p-values for each model.

### Rows - one row per trial, model and coefficient
### Cols

-  nct_id - trial ID
-  nct_id2 - trial ID with a letter appended for one trial which had two populations within the same NCT_ID
-  models - model name (b1 for baseline HbA1c models and f1-8 for ANCOVA models)
-  term - name of parameter arm_f = treatment arm, sex = sex, age10 = age (divided by ten), arm_f = treatment arm
-  estimate - point estimate - rounded to 3 DP
-  std.error - standard error of point estimate - rounded to 3 DP
-  statistic - t-statistic - rounded to 3 DP
-  p.value - p_value - Scientific notation, rounded to 3 significant figures

## age_sex_model_diags.csv

Model diagnostics for fitted model.

### Rows - one row per trial/model combination.
### Cols

-  nct_id - as above
-  nct_id2 - as above
-  models - as above
The following are written in scientific notation and rounded to 3 DP (except df, df.residual and n.obs which are integers)
-  r.squared
-  adj.r.squared
-  sigma - residual standard deviation
-  statistic - F-statistic
-  p.value
-  df - degrees of freedom 
-  logLik - log likelihood
-  AIC - Akaike information criterion 
-  BIC - Bayesian information criterion 
-  deviance - model deviance 
-  df.residual - residual degrees of freedom
-  nobs - number of observations used to fit model

##  age_sex_model_vcov.csv

Variance covariance matrix of each model converted into a correlation matrix with each element rounded to 3 DP. Please note that this is the correltion of the parameter estimates NOT of the original variables. To avoid duplication only the upper triangle is included. The data is arranged into a long format to make review easier.

### Rows - one row per trial/model/parameter1/parameter2 combination.
### Cols
-  nct_id - as above
-  nct_id2 - as above
-  models - as above
-  row - parameter name of first parameter (eg interecept)
-  col - parameter name of second parameter (eg sex)
-  r - correlation between standard errors of each parameter.

# Reference files

These are reference tables indicating the meaning of the trial arm IDs and the selected reference arm. These are pure metadata with no data. When refered to in any analysis scripts the prefix "reference_" was not used. That is just used in this export for clarity.

## reference_arms_slct_reference.csv

Indicates which was the reference arm for each model.

### rows - one row per trial. Extra rows as some names are duplicated (eg placebo and Placebo)
### Cols
 -  nct_id - as above
 - trial_lbl - arm label used in IPD
 - arm_id_unq - arm ID assigned by our team to uniquely identify each arm
 - arm_label - short arm label 
 - dc - drug class

##  reference_trial_arm_to_arm_data_all_cleaned.csv

Allows mapping across between arm name in trial IPD and arm name we used.

### Rows - one row per trial/arm combination
### Cols
nct_id - as above
trial_lbl - arm name in trial IPD
arm_id_unq - unique arm ID assingned by my team
arm_label - label for unique arm ID assigned by my team (note not unique)

## reference_ancillary_drugs_all_data_cleaned.csv

For one trial only indicates names of ancillary drugs.

### Rows - one row per trial
### cols
ancillary_id - ID indicating the ancillary drug information
nct_id - as above
drug1 - first ancillary drug
drug2- second ancillary drug
ancillary_notes - notes on ancillary drugs


## reference_arm_data_all_cleaned.csv

Indicates the drugs, doses etc of each arm.

#### rows - one row per trial and arm
#### Cols
 - nct_id - as above
 - arm_id_unq - as above
 - drug_name - name of first drug in arm
 - drug_dose - dose of first drug in arm
 - drug_unit - unit of dose
 - drug_freq - frequency
 - drug2_name - as above but second drug in arm
 - drug2_dose - as above but second drug in arm
 - drug2_unit - as above but second drug in arm
 - drug2_freq - as above but second drug in arm
 - new note - note about trial with two populations (population A and population B)


