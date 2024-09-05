This request for results contains trial summary statistics and outputs of regression models. 
For a detailed description of the request, please see detailed readme from request 6115 made on the 13th of September 2023 as the results export is in ALMOST the same format (as is a request from 8697 relating to 11 YODA/JJ trials). The only two differences are listed below.

Difference 1.The current request is for 29 trials more recently provided. The new trials are for the following sponsors:-
Sponsor         Trials
astrazeneca     1
boehringer      9
eli lilly       3
sanofi         15
takeda          1

Difference 2. There is one additional file which in a different format to the 6115 request. It is called "hba1c_categorical_trial.csv" and is described below.

hba1c_categorical_trial.csv
This file include summary data for one trial - NCT01117350, sponsor Sanofi. There is one row per arm/age category/sex stratum (16 rows in total not including the header). The columns are as follows:-

nct_id  - clinicicaltrials.gov trial ID
arm  - treatment arm 
arm_f  - treatment arm with indicator showing which arm is the reference arm
sex  - Male/Female
age_cat  - age category ([35, 51], ]51, 58], ]58, 64], or ]64, 76])
n  - total number of participants in stratum
f_mean  - first hba1c mean rounded to 2 dp
f_sd  - first hba1c standard deviation rounded to 2 dp
f_msng  - first hba1c, number of participants with missing values
l_mean  - last hba1c mean rounded to 2 dp
l_sd  - last hba1c standard deviation rounded to 2 dp
l_msng  - last hba1c, number of participants with missing values
result_mean  - last-first hba1c mean rounded to 2 dp
result_sd  - last-first hba1c standard deviation rounded to 2 dp
result_msng  - last-first hba1c, number of participants with missing values
