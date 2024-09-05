README

The results are summaries and diagnostics for regression models fit to a single trial - NCT01131676 (aka EMPAREG). The sponsor is Boehringer Ingelheim.

The models are all trial-level models and are as follows:-
- hba1c - linear models of hba1c on treatment age and sex. Modelled with last observation carried forward and baseline observation carried forward. Models include interactions and in one model a restricted cubic spline for age with 4 knots.

The model outputs are contained in the following files:-
- coef.csv - estimates (ie betas) and standard errors for each term within each model. One row per trial and model type (eg locf, bocf) and term
- vcov.csv - correlations for the standard errors in "coefs.csv". Note that these are NOT correlations of participant charcteristics, they are correlations of standard errors.
- diag_lm.csv - Model diagnostic summaries. One row per trial and model type 
