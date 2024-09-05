README

Please see contemporaneous request made on 6115 and note at the foot of this readme.

The contents of this results export request are needed to respond to the reviewer and editor comments from our recent journal submission (name of journal will be posted in Vivli chat).

The results are:-
1. Summaries and diagnostics for regression models fit to individual trials
2. Selected trial-level summary statistics.

These cover a total of 44 trials from multiple sponsors:-

astrazeneca	1
boehringer	10
eli lilly	3
gsk	1
janssen	13
sanofi	14
takeda	2


The models are all trial-level models and are as follows:-
- hba1c - linear models of hba1c on treatment age and sex. Modelled with last observation carried forward and baseline observation carried forward. Models include interactions and in one model a restricted cubic spline for age with 4 knots.
- mace - Cox models of MACE on treatment age and sex. Models include interactions and a restricted cubic spline for age with 4 knots.
- non-cardiovascular death - Cox models of non cardiovascular death before MACE on treatment age and sex. Models include interactions.
- adverse events - Outcomes are incident adverse events:- this was split into serious adverse events, hypoglycaemic events, urinary tract infection and gatrointestinal events. Outcomes were defined using combinations of MedDRA codes. Models are Poisson and negative binomial models of the outcome on age, sex and treatment with age-treatment and sex-treatment interactions.

The model outputs are contained in the following files:-
- coef.csv - estimates (ie betas) and standard errors for each term within each model. One row per trial/outcome/model type (eg poisson, locf, bocf) and term
- vcov.csv - correlations for the standard errors in "coefs.csv". Note that these are NOT correlations of participant charcteristics, they are correlations of standard errors.
- diag_cox.csv, diag_glm.csv and diag_lm.csv - Model diagnostic summaries. One row per trial/outcome/model type (eg poisson, locf, bocf)
- pred.csv - for the models with restricted cubic splines only, these are model-derived estimates and 95% confidence intervals for the outcome level for different combinations of the covariates. 

The summary statistics are as follows:-
- adverse_event_summary.csv - One row per trial. Summary statistics for adverse events in each column. Statistic shows events/person time (rate per 100 person years) [%]. If fewer than 10 events (but not zero) the number is suppressed (replaced with "<10"). Percentages are rounded to two decimal places.
- death_summary.csv - One row per trial. Statistic show population, deaths and percentage. If fewer than 10 deaths (but not zero) the number is suppressed (replaced with "<10""). Percentages are rounded to one decimal place.

Note on 6115
Please note that we are simultaneously making a similar request in repository 6115. That is for 54 other trials. These analyses are covered under both applications (6115 and 9492). We are making the request in both repositories as we have not had time to complete the data processing of the 6115 trials in the 9492 repository to our satisfaction in time for the journal revision submission.
