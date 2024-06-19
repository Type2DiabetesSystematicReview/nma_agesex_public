# Age and sex differences in controlling glycaemia and preventing cardiovascular events in type 2 diabetes: a clinical trials meta-analysis using aggregate-level and individual-participant data

This repository contains all the code and non-individual-level data for a paper analysing heterogeneity of treatment effects.
The repository contains the entire git history of the analysis unchanged, except the first few commits which were squashed to remove an absolute path which had been included in an early iteration.

## Folders

### [Created_metadata](Created_metadata/) 

CSV files which encode information about data processing. These are called by [Scripts](Scripts/)

### [Data](Data/) 

Contains data used in the analysis. These are called by [Scripts](Scripts/). The types of data include:- 

- Aggregate-level trial data
- Trial-level model outputs (coefficients, variance-covariance matrices, model diagnostics)
- Trial summary statistics (marginal summary statistics and models describing relations between covariates)
- Lists of included/excluded trials
- Risk of bias results
- Metadata

### [FromVM](FromVM/)

Model objects transferred from ML-NMR models run on high-performance machine. Too large to be uploaded to github. Model objects will be created using code in [Scripts](Scripts/)

### [Outputs](Outputs/)

Tables and figures showing analyses used in the analyses. Includes both final tables and results of intermediate analyses. Includes summaries and plots and (in the case of the multi-level network meta regression (ML-NMR) models) samples from the posterior.

### [Scratch_data](Scratch_data/)

Empty folder. Stores objects created by one script then used by another.

### [Scripts](Scripts/)

Analysis code. Summary of scripts below.

- 00__'s pull data in from original cleaning/other repositories
- 00_functions.R functions re-used across analysis
- 01 to 06a cleaning, transforming and set-up of data for HbA1c ML-NMR models
- 06b to 06c fit regression models for HbA1c. Need to run this from the terminal using Rscripts. Quite time consuming so run on HPCC
- 07-10a prepare MACE data for ML-NMR models
- 10b-10f fit regression models for MACE. Need to run this from the terminal using Rscripts. A bit faster than the HbA1c models, but still best to run on HPCC
- 11-15 post-processing of ML-NMR models to produce summaries
- 16 additional figures for manuscript

### [Supporting](Supporting/)

Contains some additional exports from the trial-level analysis and summaries, illustrating, for example, how well the summaries of age and censoring time fitted true distributions.



