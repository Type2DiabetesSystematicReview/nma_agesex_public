README

THis export includes results of regression models of HbA1c on age.

The results are :-
1. Summaries for regression models fit to individual trials

These cover a total of 53 trials from multiple sponsors:-

astrazeneca	4
boehringer	28
eli_lilly	9
takeda	12

The models are all trial-level models and are of one outcome, HbA1c
- hba1c - linear models of hba1c on treatment age and sex. Modelled with last observation carried forward. There is one model for each trial and it includes age-treatment and sex-treatment interactions. Age is modelled using a restricted cubic spline for age with 4 knots. Unlike the previous analysis (exported last Monday), the knots were set at ages 40, 50, 60 and 75 years.

The model outputs are contained in the following files:-
- coefs.csv - estimates (ie betas) and standard errors for each term within each model. One row per trial and term
- vcov.csv - correlations for the standard errors in "coefs.csv". Note that these are NOT correlations of participant characteristics, they are correlations of standard errors.

Note on 9492
Please note that we are simultaneously making a similar request in repository 9492.