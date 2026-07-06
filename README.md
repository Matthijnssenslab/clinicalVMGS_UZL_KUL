# Publication Code

Code and data for reproducing the CSF/BAL hybrid-capture viral metagenomics analyses.

## Structure

```text
data_input/   input data
scripts/      R scripts
tables/       generated tables
figures/      PDF figures
```

## Run

```bash
Rscript scripts/00_run_all.R
```

Required R packages: `ggplot2`, `lhs`, `readxl`, `patchwork`, `ggnewscale`.

## Scripts

```text
00_run_all.R
01_build_records_and_lhs.R
02_figure2_qcmd_viral_load_metrics.R
03_figure3bc_lhs_cutoff_distributions.R
04_figure3d_qcmd_top100_validation.R
05_methods_qcmd_selected_cutoff_validation.R
06_figure4a_clinical_sensitivity_specificity.R
07_figure4b_bal_cq_sensitivity.R
08_figure5_pcr_metagenomics_concordance.R
09_supplementary_figure1_top100_cutoff_robustness.R
```
