#!/usr/bin/env Rscript

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))

required_packages <- c("ggplot2", "lhs", "readxl", "patchwork", "ggnewscale")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(sprintf(
    "Missing required R package(s): %s. Install them before running this script.",
    paste(missing_packages, collapse = ", ")
  ))
}

run_script <- function(file) {
  message(sprintf("--- Running %s ---", file))
  env <- new.env(parent = globalenv())
  scratch_pdf <- tempfile(fileext = ".pdf")
  grDevices::pdf(scratch_pdf)
  on.exit({
    if (grDevices::dev.cur() > 1L) {
      grDevices::dev.off()
    }
    unlink(scratch_pdf)
  }, add = TRUE)
  source(file.path(scripts_dir, file), local = env)
  env
}

save_pdf <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

run_script("01_build_records_and_lhs.R")

qcmd <- run_script("02_figure2_qcmd_viral_load_metrics.R")
save_pdf(qcmd$p_cov, "Figure2a_QCMD_viral_load_vs_genome_coverage.pdf", 6.4, 4.8)
save_pdf(qcmd$p_rpm, "Figure2b_QCMD_viral_load_vs_RPM.pdf", 6.4, 4.8)

cutoffs <- run_script("03_figure3bc_lhs_cutoff_distributions.R")
save_pdf(cutoffs$fig, "Figure3bc_LHS_cutoff_distributions.pdf", 11.0, 5.2)

qcmd_top100 <- run_script("04_figure3d_qcmd_top100_validation.R")
save_pdf(qcmd_top100$p, "Figure3d_QCMD_top100_cutoff_validation.pdf", 10.0, 4.8)

qcmd_selected <- run_script("05_methods_qcmd_selected_cutoff_validation.R")
save_pdf(qcmd_selected$p, "Methods_QCMD_selected_cutoff_validation.pdf", 5.8, 4.2)

clinical <- run_script("06_figure4a_clinical_sensitivity_specificity.R")
save_pdf(clinical$p, "Figure4a_clinical_sensitivity_specificity.pdf", 8.5, 5.0)

cq <- run_script("07_figure4b_bal_cq_sensitivity.R")
save_pdf(cq$p, "Figure4b_BAL_Cq_sensitivity.pdf", 6.2, 4.8)

run_script("08_figure5_pcr_metagenomics_concordance.R")

supp <- run_script("09_supplementary_figure1_top100_cutoff_robustness.R")
save_pdf(supp$p, "Supplementary_Figure1_top100_cutoff_robustness.pdf", 10.0, 4.8)

message("Done. Outputs are in tables/ and figures/.")
