# Age and sex differences in controlling glycaemia and preventing cardiovascular events in type 2 diabetes: a clinical trials meta-analysis using aggregate-level and individual-participant data

This repository contains all the code and non-individual-level data for the paper "Age and Sex Differences in Efficacy of Treatments for Type 2 Diabetes A Network Meta-Analysis paper analysing heterogeneity of treatment effects" available at  [JAMA](https://jamanetwork.com/journals/jama/fullarticle/2829848), [Europe PMC](https://europepmc.org/article/MED/39899304) and [the University of Glasgow Enlighten repository](https://eprints.gla.ac.uk/343359/). The [pubmed link is here](https://pubmed.ncbi.nlm.nih.gov/39899304/).

This repository contains the entire git history of the analysis unchanged, except the first few commits which were squashed to remove an absolute path which had been included in an early iteration.

## Folders

### [Created_metadata](Created_metadata/)

CSV files which encode information about data processing. These are called by [Scripts](Scripts/)

### [Data](Data/)

Contains data used in the analysis. These are called by [Scripts](Scripts/). The types of data include:-

-   Aggregate-level trial data
-   Trial-level model outputs (coefficients, variance-covariance matrices, model diagnostics)
-   Trial summary statistics (marginal summary statistics and models describing relations between covariates)
-   Lists of included/excluded trials
-   Risk of bias results
-   Metadata

### [FromVM](FromVM/)

Model objects transferred from ML-NMR models run on high-performance machine. Too large to be uploaded to github. Model objects will be created using code in [Scripts](Scripts/)

### [Outputs](Outputs/)

Tables and figures showing analyses used in the analyses. Includes both final tables and results of intermediate analyses. Includes summaries and plots and (in the case of the multi-level network meta regression (ML-NMR) models) samples from the posterior.

### [Scratch_data](Scratch_data/)

Empty folder. Stores objects created by one script then used by another.

### [Scripts](Scripts/)

Analysis code. Summary of scripts below.

-   00\_\_'s pull data in from original cleaning/other repositories
-   00_functions.R functions re-used across analysis
-   01 to 06a cleaning, transforming and set-up of data for HbA1c ML-NMR models
-   06b to 06c fit regression models for HbA1c. Need to run this from the terminal using Rscripts. Quite time consuming so run on HPCC
-   07-10a prepare MACE data for ML-NMR models
-   10b-10f fit regression models for MACE. Need to run this from the terminal using Rscripts. A bit faster than the HbA1c models, but still best to run on HPCC
-   11-15 post-processing of ML-NMR models to produce summaries
-   16 additional figures for manuscript

### [Supporting](Supporting/)

Contains some additional exports from the trial-level analysis and summaries, illustrating, for example, how well the summaries of age and censoring time fitted true distributions.

# Running the code

The analysis code was written on Ubuntu 22.04.03, on a machine running Windows 11 Pro. The ML-NMR models were run on the [MVLS Advanced Research System (MARS) high performance computing resource](https://mars.ice.gla.ac.uk/).

The R version used was 4.4.0.

The following versions of R packages were used:-

-   broom - 1.0.5
-   cowplot - 1.1.2
-   dplyr - 1.1.4
-   multinma - development version install from devtools::install_github("dmphillippo/multinma", ref = "feature/set_agd_regression")
-   patchwork - 1.2.0
-   purrr - 1.0.2
-   rstan - 2.32.6
-   stringr - 1.5.1
-   survival - 3.2-13
-   tidyr - 1.3.1
-   tidyverse - 2.0.0
-   truncnorm - 1.0.9
