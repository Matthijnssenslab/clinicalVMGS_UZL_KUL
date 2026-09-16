# Clinical viral metagenomics in cerebrospinal fluid and bronchoalveolar lavage

Thank you for your interest in our study!

This repository contains input datasets and R scripts to reproduce the CSF/BAL hybrid-capture viral metagenomics analyses comparing AVITI and Oxford Nanopore Technologies (ONT).

## Analyses

The workflow uses Latin hypercube sampling (LHS) to evaluate detection cutoffs, validates them using QCMD reference panels, and compares metagenomic detections with clinical PCR results. It generates the tables and figures for viral-load comparisons, cutoff selection, clinical sensitivity and specificity, and PCR–metagenomics concordance.

The included inputs are processed metagenomics outputs and reference/comparison data. Raw-read processing is not part of this R workflow.

## Repository structure

```text
data_input/   Processed AVITI and ONT inputs, QCMD data, and PCR comparisons
scripts/      Numbered R analysis and figure scripts
tables/      Generated record-level tables, LHS results, and selected cutoffs
figures/      PDF figures created when the workflow runs
```

Generated tables are included in the repository; figure PDFs are generated locally.

## Requirements

Install the required R packages:

```r
install.packages(c("ggplot2", "lhs", "readxl", "patchwork", "ggnewscale"))
```

## Usage

Clone the repository and run the complete workflow from its root directory:

```bash
git clone https://github.com/Matthijnssenslab/clinicalVMGS_UZL_KUL.git
cd clinicalVMGS_UZL_KUL
Rscript scripts/00_run_all.R
```

Alternatively, open `publication_code_tbu.Rproj` in RStudio and run:

```r
source("scripts/00_run_all.R")
```

The workflow regenerates tables in `tables/` and saves PDF figures in `figures/`.

## Figure scripts

| Script | Analysis / output |
| --- | --- |
| `01_build_records_and_lhs.R` | Record-level tables and LHS cutoff selection |
| `02_figure2_qcmd_viral_load_metrics.R` | Figure 2: QCMD viral load, genome coverage, and reads per million |
| `03_figure3bc_lhs_cutoff_distributions.R` | Figure 3b–c: cutoff distributions |
| `04_figure3d_qcmd_top100_validation.R` | Figure 3d: QCMD validation of the top 100 cutoff combinations |
| `05_methods_qcmd_selected_cutoff_validation.R` | QCMD validation of the selected cutoffs |
| `06_figure4a_clinical_sensitivity_specificity.R` | Figure 4a: clinical sensitivity and specificity |
| `07_figure4b_bal_cq_sensitivity.R` | Figure 4b: BAL sensitivity by PCR Cq |
| `08_figure5_pcr_metagenomics_concordance.R` | Figure 5: PCR–metagenomics concordance |
| `09_supplementary_figure1_top100_cutoff_robustness.R` | Supplementary Figure 1: cutoff robustness |
