#! /usr/bin/bash

echo "Runs all scripts that are run locally (ie not ones running the nma in the HPCC)"
cd ~/nma_agesex_public
#ls

# Rscript Scripts/01a_process_agg_outcome.R
# Rscript Scripts/01b_process_agg_agesex.R
# Rscript Scripts/02_process_agg_standard_errors.R
# Rscript Scripts/03a_process_ipd_consolidate.R
# Rscript Scripts/03b_process_ipd_coefvcov.R
# Rscript Scripts/04_process_arm_labels.R
# Rscript Scripts/04b_arm_label_reverse.R
# Rscript Scripts/05_join_armlabels_agg_ipd.R
# Rscript Scripts/06a_setup_hba1c_reg.R
# 
# Rscript Scripts/07_create_mace_pseudoipd.R
# Rscript Scripts/08_process_mace_trials.R
# 
# Rscript Scripts/09a_summarise_plot_ages_hba1c.R
# Rscript Scripts/09b_summarise_plot_ages_mace.R
# Rscript Scripts/10a_plot_tabulate_prepare_mace.R

echo 11_postprocess_models.R
Rscript Scripts/11_postprocess_models.R
echo 12a_plot_model_summaries_hba1c.R
Rscript Scripts/12a_plot_model_summaries_hba1c.R
echo 12b_print_hba1c_priors.R
Rscript Scripts/12b_print_hba1c_priors.R
echo 12c_main_effects_hba1c_plot.R
Rscript Scripts/12c_main_effects_hba1c_plot.R
echo 12d_mace_plots.R
Rscript Scripts/12d_mace_plots.R
echo 13a_tables_hba1c.R
Rscript Scripts/13a_tables_hba1c.R
echo 14a_relative_effects_reshape_data.R
Rscript Scripts/14a_relative_effects_reshape_data.R
echo 14b_relative_effects_tabulate.R
Rscript Scripts/14b_relative_effects_tabulate.R
echo 14c_relative_effects_class_lvl.R
Rscript Scripts/14c_relative_effects_class_lvl.R
echo 15_mace_compare_cfs_pseudo.R
Rscript Scripts/15_mace_compare_cfs_pseudo.R
echo 16_manuscript_figure2.R
Rscript Scripts/16_manuscript_figure2.R
echo 17_adverse_events.R
Rscript Scripts/17_adverse_events.R

echo 19_ae_mdl_summaries.R
Rscript Scripts/19_ae_mdl_summaries.R
echo 20_ae_mdl_plots.R
Rscript Scripts/20_ae_mdl_plots.R
echo 21_hba1c_and_mace_trials_sens_analysis.R
Rscript Scripts/21_hba1c_and_mace_trials_sens_analysis.R
echo 22_death_summaries.R
Rscript Scripts/22_death_summaries.R
echo 23_rcs_hba1c_setup.R
Rscript Scripts/23_rcs_hba1c_setup.R
echo 24_rcs_plot_model_outputs.R
Rscript Scripts/24_rcs_plot_model_outputs.R
echo 25_rcs_mace.R
Rscript Scripts/25_rcs_mace.R
echo 26_summarise_missingess.R
Rscript Scripts/26_summarise_missingess.R
echo 27_prepare_inline_for_ms.R
Rscript Scripts/27_prepare_inline_for_ms.R
