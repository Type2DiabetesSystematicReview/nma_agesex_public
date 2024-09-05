README

This export includes results of regression models of HbA1c on age.

The results are:-
1. Summaries for regression models fit to individual trials

These cover a total of 34 trials from multiple sponsors:-

sponsor         n
astrazeneca     1
boehringer      8
eli lilly       3
janssen        10
sanofi         11
takeda          1


The models are all trial-level models and are of one outcome. HbA1c:-
- hba1c - linear models of hba1c on treatment age and sex. Modelled with last observation carried forward. Models include interactions and a restricted cubic spline for age with 4 knots. Unlike the previous analysis (exported last Monday) the knots are set at ages 40, 50, 60 and 75 years.

The model outputs are contained in the following files:-
- coef.csv - estimates (ie betas) and standard errors for each term within each model. One row per trial and term
- vcov.csv - correlations for the standard errors in "coefs.csv". Note that these are NOT correlations of participant charcteristics, they are correlations of standard errors.

Note on 6115
Please note that we are simultaneously making a similar request in repository 6115. 