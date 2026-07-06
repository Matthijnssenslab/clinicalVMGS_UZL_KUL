#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplot2))

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))
source_dir <- data_dir

N_TOP <- as.integer(Sys.getenv("OPERATIONAL_TOP_N", Sys.getenv("FIG2A_N_TOP", "100")))
if (!is.finite(N_TOP) || N_TOP < 1L) N_TOP <- 100L
MIN_SENS_FLOOR <- 0
MIN_SPEC_FLOOR <- 0
FIXED_AVITI_COVERED_BASES <- 150
FIXED_AVITI_READS_ALIGNED <- 1
FIXED_AVITI_RPKMF <- 0
FIXED_ONT_ASSIGNMENT_RANK <- 1
DONOR_MIN_READS <- 1
DONOR_MIN_RPM <- 1000

OUTSET_DONOR_AVITI <- c(
  "Human mastadenovirus E" = 5000000
)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
trim_chr <- function(x) trimws(as.character(x))
safe_div <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)

lhs_param_burden <- function(lhs_df, method) {
  params <- if (method == "AVITI") {
    c("aviti_bleed_ratio_min_pct", "aviti_min_genome_coverage_pct", "aviti_min_reads_per_million_filtered")
  } else {
    c("ont_bleed_ratio_min_pct", "ont_min_nb_reads", "ont_min_genome_coverage_pct")
  }
  M <- sapply(params, function(nm) {
    v <- if (nm %in% names(lhs_df)) as_num(lhs_df[[nm]]) else rep(NA_real_, nrow(lhs_df))
    v[!is.finite(v)] <- 0
    v
  })
  rowSums(M, na.rm = TRUE)
}

pick_top <- function(df, n = 100L, method = c("AVITI", "ONT")) {
  method <- match.arg(method)
  sens <- as_num(df$sensitivity_pct)
  spec <- as_num(df$specificity_like_pct)
  ok <- is.finite(sens) & is.finite(spec) & sens >= MIN_SENS_FLOOR & spec >= MIN_SPEC_FLOOR
  pool <- if (any(ok)) df[ok, , drop = FALSE] else df
  spec100 <- pool[is.finite(as_num(pool$specificity_like_pct)) &
                    abs(as_num(pool$specificity_like_pct) - 100) < 1e-9, , drop = FALSE]
  if (nrow(spec100) > 0L) pool <- spec100
  br <- lhs_param_burden(pool, method)
  ord <- order(
    -as_num(pool$sensitivity_pct),
    -as_num(pool$min_sens_spec_like_pct),
    -as_num(pool$avg_sens_spec_like_pct),
    -as_num(pool$specificity_like_pct),
    br,
    as_num(pool$lhs_id)
  )
  pool[ord[seq_len(min(n, nrow(pool)))], , drop = FALSE]
}

any_match <- function(records, col, keyword) {
  any(grepl(keyword, records[[col]], ignore.case = TRUE))
}
sample_level_sensitivity <- function(det, detected_col) {
  keys <- paste(det$year, det$sample, sep = "||")
  pred_ok <- tapply(det[[detected_col]], keys, all)
  sum(pred_ok, na.rm = TRUE) / max(length(pred_ok), 1L)
}
EMPTY_SAMPLE_KEYS <- c("2024||Sample5", "2025||Sample2")
empty_sample_specificity <- function(spec_df, fp_col) {
  row_keys  <- paste(spec_df$year, spec_df$sample, sep = "||")
  df        <- spec_df[row_keys %in% EMPTY_SAMPLE_KEYS, , drop = FALSE]
  if (nrow(df) == 0L) return(list(pct = NA_real_, num = 0L, den = 0L))
  keys      <- paste(df$year, df$sample, sep = "||")
  fp_per    <- tapply(df[[fp_col]], keys, any)
  tn        <- as.integer(sum(!fp_per, na.rm = TRUE))
  total     <- as.integer(length(fp_per))
  list(pct = 100 * tn / max(total, 1L), num = tn, den = total)
}

score_detection <- function(truth, av_df, on_df, av_map, on_map, vkeys) {
  rows <- vector("list", nrow(truth))
  for (i in seq_len(nrow(truth))) {
    yr <- truth$year[i]; smp <- truth$sample[i]; vir <- truth$virus[i]
    av_sid <- av_map$sid[av_map$year == yr & av_map$sample == smp]
    on_sid <- on_map$sid[on_map$year == yr & on_map$sample == smp]
    kw <- vkeys[vkeys$virus == vir, , drop = FALSE]
    av_pass <- av_df[av_df$sample_ID == av_sid & av_df$pass, , drop = FALSE]
    on_pass <- on_df[on_df$SampleName == on_sid & on_df$pass, , drop = FALSE]
    rows[[i]] <- data.frame(
      year = yr, sample = smp, virus = vir,
      aviti_detected = if (nrow(kw) > 0 && nrow(av_pass) > 0) any_match(av_pass, "species", kw$aviti) else FALSE,
      ont_detected = if (nrow(kw) > 0 && nrow(on_pass) > 0) any_match(on_pass, "StrainName", kw$ont) else FALSE,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

score_specificity <- function(truth_pos, universe, av_df, on_df, av_map, on_map, vkeys) {
  all_samples <- unique(rbind(
    av_map[av_map$sample != "NC", c("year", "sample")],
    on_map[, c("year", "sample")]
  ))
  all_samples <- all_samples[all_samples$sample != "NC", , drop = FALSE]

  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(all_samples))) {
    yr <- all_samples$year[i]; smp <- all_samples$sample[i]
    expected <- truth_pos$virus[truth_pos$year == yr & truth_pos$sample == smp]
    non_expected <- setdiff(universe, expected)
    av_sid <- av_map$sid[av_map$year == yr & av_map$sample == smp]
    on_sid <- on_map$sid[on_map$year == yr & on_map$sample == smp]
    av_pass <- av_df[av_df$sample_ID == av_sid & av_df$pass, , drop = FALSE]
    on_pass <- on_df[on_df$SampleName == on_sid & on_df$pass, , drop = FALSE]
    for (vir in non_expected) {
      kw <- vkeys[vkeys$virus == vir, , drop = FALSE]
      k <- k + 1L
      rows[[k]] <- data.frame(
        year   = yr, sample = smp,
        aviti_fp = if (nrow(kw) > 0 && nrow(av_pass) > 0) any_match(av_pass, "species", kw$aviti) else FALSE,
        ont_fp   = if (nrow(kw) > 0 && nrow(on_pass) > 0) any_match(on_pass, "StrainName", kw$ont) else FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
ov_raw <- read.csv(file.path(source_dir, "QCMD_overview.csv"), skip = 1, stringsAsFactors = FALSE)
names(ov_raw)[1:4] <- c("year", "sample", "virus", "vl_log10")
ov_raw$year <- as.integer(as_num(ov_raw$year))
ov_raw$sample <- trim_chr(ov_raw$sample)
ov_raw$virus <- trim_chr(ov_raw$virus)
truth_pos <- ov_raw[!is.na(ov_raw$virus) & ov_raw$virus != "-", c("year", "sample", "virus")]
truth_pos$virus[grepl("Influenza A", truth_pos$virus)] <- "Influenza A"
truth_pos <- unique(truth_pos)
QCMD_UNIVERSE <- unique(truth_pos$virus)

aviti_sample_map <- data.frame(
  year = c(rep(2024, 7), rep(2025, 7)),
  sample = c(paste0("Sample", 1:6), "NC", paste0("Sample", 1:6), "NC"),
  sid = c("NGS_1T", "NGS_2T", "NGS_3T", "NGS_4T", "NGS_5T", "NGS_6T", "NGS_NC_T",
          "QCMD1_Q1", "QCMD2_Q2", "QCMD3_Q3", "QCMD4_Q4", "QCMD5_Q5", "QCMD6_Q6", "NC_QCMD"),
  stringsAsFactors = FALSE
)
aviti_sample_map$sample[aviti_sample_map$sample == "Sample6" & aviti_sample_map$year == 2025] <- "Sample6*"

ont_sample_map <- data.frame(
  year = c(rep(2024, 7), rep(2025, 6)),
  sample = c(paste0("Sample", 1:6), "NC", paste0("Sample", 1:6)),
  sid = c("QCMD2024_sample1", "QCMD2024_sample2", "QCMD2024_sample3",
          "QCMD2024-sample4", "QCMD2024-sample5", "QCMD2024_sample6", "QCMD2024-neg",
          "NGS_meta_25S_01", "NGS_meta_25S_02", "NGS_meta_25S_03",
          "NGS_meta_25S_04", "NGS_meta_25S_05", "NGS_meta_25S_06"),
  stringsAsFactors = FALSE
)
ont_sample_map$sample[ont_sample_map$sample == "Sample6" & ont_sample_map$year == 2025] <- "Sample6*"

virus_keys <- data.frame(
  virus = c("MPXV", "HIV-1", "EVD68", "HSV-1", "HSV-2", "CMV", "Influenza A", "HadV-4"),
  aviti = c("Monkeypox virus", "immunodeficiency virus 1", "Enterovirus D",
            "Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "betaherpesvirus 5",
            "Influenza A virus", "Human mastadenovirus E"),
  ont = c("Monkeypox virus", "immunodeficiency virus 1", "deconjuncti",
          "humanalpha1", "humanalpha2", "humanbeta5", "Alphainfluenzavirus", "exoticum"),
  stringsAsFactors = FALSE
)
av <- rbind(
  read.csv(file.path(source_dir, "2024_QCMD.csv"), stringsAsFactors = FALSE),
  read.csv(file.path(source_dir, "2025_QCMD.csv"), stringsAsFactors = FALSE)
)
av$year <- ifelse(av$sample_ID %in% aviti_sample_map$sid[aviti_sample_map$year == 2024], 2024, 2025)
av$genome_cov_pct <- safe_div(av$covered_bases, av$reference_length) * 100
av$rpm_filtered <- safe_div(av$reads_aligned, av$total_filtered_reads_in_sample) * 1e6
av_nonNC <- av[!grepl("NC|neg", av$sample_ID, ignore.case = TRUE), , drop = FALSE]
av_sum <- aggregate(cbind(reads_aligned, rpm_filtered) ~ sample_ID + year + species, data = av_nonNC, FUN = sum)
av_max <- aggregate(cbind(reads_aligned, rpm_filtered) ~ year + species, data = av_sum, FUN = max)
names(av_max)[3:4] <- c("taxon_max_reads", "taxon_max_rpm")
av_sum <- merge(av_sum, av_max, by = c("year", "species"), all.x = TRUE)
av_sum$bleed_pct <- 100 * safe_div(av_sum$reads_aligned, av_sum$taxon_max_reads)
av <- merge(av, av_sum[, c("sample_ID", "species", "taxon_max_reads", "taxon_max_rpm", "bleed_pct")],
            by = c("sample_ID", "species"), all.x = TRUE)
outset_reads <- OUTSET_DONOR_AVITI[av$species]
has_outset <- !is.na(outset_reads)
if (any(has_outset)) {
  av$bleed_pct[has_outset] <- av$bleed_pct[has_outset] * av$taxon_max_reads[has_outset] / outset_reads[has_outset]
}
av$donor_taxon <- is.finite(av$taxon_max_reads) & av$taxon_max_reads >= DONOR_MIN_READS
av_nc <- av[grepl("NC|neg", av$sample_ID, ignore.case = TRUE) & is.finite(av$genome_cov_pct), , drop = FALSE]
av_nc_max <- aggregate(genome_cov_pct ~ species, data = av_nc, FUN = max)
names(av_nc_max)[2] <- "nc_gc_max"
av <- merge(av, av_nc_max, by = "species", all.x = TRUE)
on <- as.data.frame(read_excel(file.path(source_dir, "Query.xls"), sheet = "Sheet0"))
qcmd_sids <- c(ont_sample_map$sid, "QCMD2024-neg")
on <- on[on$SampleName %in% qcmd_sids, , drop = FALSE]
on$year <- ifelse(on$SampleName %in% ont_sample_map$sid[ont_sample_map$year == 2024], 2024, 2025)
mapped_tsv <- read.delim(file.path(source_dir, "mapped_reads_per_sample.tsv"), stringsAsFactors = FALSE)
on$total_mapped <- mapped_tsv$total_mapped_reads[match(on$SampleName, mapped_tsv$sample_name)]
on$rpm_mapped <- safe_div(on$NbReads, on$total_mapped) * 1e6
on$gc_pct <- as_num(on$GenomeCoveragePct)
al <- tolower(trim_chr(on$AssignmentLevel))
on$assignment_rank <- 0L
on$assignment_rank[al == "discovery"] <- 1L
on$assignment_rank[al == "assignment"] <- 2L
on_nonNC <- on[!grepl("neg|NC", on$SampleName, ignore.case = TRUE), , drop = FALSE]
on_sum <- aggregate(cbind(NbReads, rpm_mapped) ~ SampleName + year + StrainName, data = on_nonNC, FUN = sum)
on_max <- aggregate(cbind(NbReads, rpm_mapped) ~ year + StrainName, data = on_sum, FUN = max)
names(on_max)[3:4] <- c("taxon_max_reads", "taxon_max_rpm")
on_sum <- merge(on_sum, on_max, by = c("year", "StrainName"), all.x = TRUE)
on_sum$bleed_pct <- 100 * safe_div(on_sum$NbReads, on_sum$taxon_max_reads)
on <- merge(on, on_sum[, c("SampleName", "StrainName", "taxon_max_reads", "taxon_max_rpm", "bleed_pct")],
            by = c("SampleName", "StrainName"), all.x = TRUE)
on$donor_taxon <- is.finite(on$taxon_max_reads) & is.finite(on$taxon_max_rpm) &
  on$taxon_max_reads >= DONOR_MIN_READS & on$taxon_max_rpm >= DONOR_MIN_RPM
on_nc <- on[on$SampleName == "QCMD2024-neg" & is.finite(on$gc_pct), , drop = FALSE]
on_nc_max <- aggregate(gc_pct ~ StrainName, data = on_nc, FUN = max)
names(on_nc_max)[2] <- "nc_gc_max"
on <- merge(on, on_nc_max, by = "StrainName", all.x = TRUE)
lhs_aviti_all <- read.csv(file.path(tables_dir, "lhs_aviti.csv"), stringsAsFactors = FALSE, check.names = FALSE)
lhs_ont_all <- read.csv(file.path(tables_dir, "lhs_ont.csv"), stringsAsFactors = FALSE, check.names = FALSE)
lhs_av <- pick_top(
  lhs_aviti_all,
  n = N_TOP, method = "AVITI"
)
lhs_on <- pick_top(
  lhs_ont_all,
  n = N_TOP, method = "ONT"
)

eval_rows <- list()
k <- 0L
for (i in seq_len(nrow(lhs_av))) {
  ac <- lhs_av[i, , drop = FALSE]
  av$pass <- with(av,
    is.finite(covered_bases) & covered_bases >= FIXED_AVITI_COVERED_BASES &
      is.finite(reads_aligned) & reads_aligned >= FIXED_AVITI_READS_ALIGNED &
      is.finite(RPKMF) & RPKMF >= FIXED_AVITI_RPKMF &
      is.finite(genome_cov_pct) & genome_cov_pct >= as_num(ac$aviti_min_genome_coverage_pct) &
      is.finite(rpm_filtered) & rpm_filtered >= as_num(ac$aviti_min_reads_per_million_filtered) &
      (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(ac$aviti_bleed_ratio_min_pct))) &
      (is.na(nc_gc_max) | genome_cov_pct >= nc_gc_max)
  )
  on$pass <- FALSE
  det <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  spec <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  sens_det  <- mean(det$aviti_detected)
  sens_samp <- sample_level_sensitivity(det, "aviti_detected")
  sp        <- 1 - mean(spec$aviti_fp)
  sp_nc     <- empty_sample_specificity(spec, "aviti_fp")$pct
  k <- k + 1L
  eval_rows[[k]] <- data.frame(method = "AVITI", lhs_id = as.integer(as_num(ac$lhs_id)),
                               sensitivity_detection_pct = 100 * sens_det,
                               sensitivity_sample_pct    = 100 * sens_samp,
                               specificity_like_pct      = 100 * sp,
                               specificity_sample_pct    = sp_nc,
                               stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(lhs_on))) {
  oc <- lhs_on[i, , drop = FALSE]
  on$pass <- with(on,
    is.finite(assignment_rank) & assignment_rank >= FIXED_ONT_ASSIGNMENT_RANK &
      is.finite(NbReads) & NbReads >= as_num(oc$ont_min_nb_reads) &
      is.finite(gc_pct) & gc_pct >= as_num(oc$ont_min_genome_coverage_pct) &
      (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(oc$ont_bleed_ratio_min_pct))) &
      (is.na(nc_gc_max) | gc_pct >= nc_gc_max)
  )
  av$pass <- FALSE
  det <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  spec <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  sens_det  <- mean(det$ont_detected)
  sens_samp <- sample_level_sensitivity(det, "ont_detected")
  sp        <- 1 - mean(spec$ont_fp)
  sp_nc     <- empty_sample_specificity(spec, "ont_fp")$pct
  k <- k + 1L
  eval_rows[[k]] <- data.frame(method = "ONT", lhs_id = as.integer(as_num(oc$lhs_id)),
                               sensitivity_detection_pct = 100 * sens_det,
                               sensitivity_sample_pct    = 100 * sens_samp,
                               specificity_like_pct      = 100 * sp,
                               specificity_sample_pct    = sp_nc,
                               stringsAsFactors = FALSE)
}

eval_df <- do.call(rbind, eval_rows)
eval_df <- eval_df[order(eval_df$method, eval_df$lhs_id), , drop = FALSE]
write.csv(eval_df, file.path(tables_dir, "qcmd_external_validation_top100.csv"), row.names = FALSE)

selected_co <- read.csv(file.path(tables_dir, "selected_cutoffs.csv"), stringsAsFactors = FALSE, check.names = FALSE)
aviti_sel_id <- as.integer(as_num(selected_co$lhs_id[selected_co$method == "AVITI"]))
ont_sel_id <- as.integer(as_num(selected_co$lhs_id[selected_co$method == "ONT"]))
best_ac <- lhs_aviti_all[as.integer(as_num(lhs_aviti_all$lhs_id)) == aviti_sel_id, , drop = FALSE]
best_oc <- lhs_ont_all[as.integer(as_num(lhs_ont_all$lhs_id)) == ont_sel_id, , drop = FALSE]
if (nrow(best_ac) != 1L || nrow(best_oc) != 1L) {
  stop("ExternalValidation_top100: selected_cutoffs lhs_id not found in lhs_aviti.csv / lhs_ont.csv")
}
ac <- best_ac
av$pass <- with(av,
  is.finite(covered_bases) & covered_bases >= FIXED_AVITI_COVERED_BASES &
    is.finite(reads_aligned) & reads_aligned >= FIXED_AVITI_READS_ALIGNED &
    is.finite(RPKMF) & RPKMF >= FIXED_AVITI_RPKMF &
    is.finite(genome_cov_pct) & genome_cov_pct >= as_num(ac$aviti_min_genome_coverage_pct) &
    is.finite(rpm_filtered) & rpm_filtered >= as_num(ac$aviti_min_reads_per_million_filtered) &
    (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(ac$aviti_bleed_ratio_min_pct))) &
    (is.na(nc_gc_max) | genome_cov_pct >= nc_gc_max)
)
on$pass <- FALSE
det_best_av  <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
spec_best_av <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
av_det_tp   <- sum(det_best_av$aviti_detected)
av_det_n    <- nrow(det_best_av)
av_samp_tmp <- tapply(det_best_av$aviti_detected,
                      paste(det_best_av$year, det_best_av$sample, sep = "||"), all)
av_samp_tp  <- as.integer(sum(av_samp_tmp, na.rm = TRUE))
av_samp_n   <- as.integer(length(av_samp_tmp))
av_splike_tn <- sum(!spec_best_av$aviti_fp)
av_splike_n  <- nrow(spec_best_av)
av_spsamp   <- empty_sample_specificity(spec_best_av, "aviti_fp")
best_eval_av <- data.frame(
  method = "AVITI", lhs_id = aviti_sel_id,
  sensitivity_detection_pct = 100 * av_det_tp  / av_det_n,
  sensitivity_sample_pct    = 100 * av_samp_tp  / av_samp_n,
  specificity_like_pct      = 100 * av_splike_tn / av_splike_n,
  specificity_sample_pct    = av_spsamp$pct,
  det_num = av_det_tp,    det_den = av_det_n,
  samp_num = av_samp_tp,  samp_den = av_samp_n,
  splike_num = av_splike_tn, splike_den = av_splike_n,
  spsamp_num = av_spsamp$num, spsamp_den = av_spsamp$den,
  stringsAsFactors = FALSE
)
oc <- best_oc
on$pass <- with(on,
  is.finite(assignment_rank) & assignment_rank >= FIXED_ONT_ASSIGNMENT_RANK &
    is.finite(NbReads) & NbReads >= as_num(oc$ont_min_nb_reads) &
    is.finite(gc_pct) & gc_pct >= as_num(oc$ont_min_genome_coverage_pct) &
    (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(oc$ont_bleed_ratio_min_pct))) &
    (is.na(nc_gc_max) | gc_pct >= nc_gc_max)
)
av$pass <- FALSE
det_best_on  <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
spec_best_on <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
on_det_tp   <- sum(det_best_on$ont_detected)
on_det_n    <- nrow(det_best_on)
on_samp_tmp <- tapply(det_best_on$ont_detected,
                      paste(det_best_on$year, det_best_on$sample, sep = "||"), all)
on_samp_tp  <- as.integer(sum(on_samp_tmp, na.rm = TRUE))
on_samp_n   <- as.integer(length(on_samp_tmp))
on_splike_tn <- sum(!spec_best_on$ont_fp)
on_splike_n  <- nrow(spec_best_on)
on_spsamp   <- empty_sample_specificity(spec_best_on, "ont_fp")
best_eval_on <- data.frame(
  method = "ONT", lhs_id = ont_sel_id,
  sensitivity_detection_pct = 100 * on_det_tp  / on_det_n,
  sensitivity_sample_pct    = 100 * on_samp_tp  / on_samp_n,
  specificity_like_pct      = 100 * on_splike_tn / on_splike_n,
  specificity_sample_pct    = on_spsamp$pct,
  det_num = on_det_tp,    det_den = on_det_n,
  samp_num = on_samp_tp,  samp_den = on_samp_n,
  splike_num = on_splike_tn, splike_den = on_splike_n,
  spsamp_num = on_spsamp$num, spsamp_den = on_spsamp$den,
  stringsAsFactors = FALSE
)
best_eval <- rbind(best_eval_av, best_eval_on)

median_cutoff_row <- function(feasible_df, method_label) {
  if (method_label == "AVITI") {
    params <- c("aviti_bleed_ratio_min_pct",
                "aviti_min_genome_coverage_pct",
                "aviti_min_reads_per_million_filtered")
  } else {
    params <- c("ont_bleed_ratio_min_pct",
                "ont_min_nb_reads",
                "ont_min_genome_coverage_pct")
  }
  out <- as.list(setNames(rep(NA_real_, length(params)), params))
  for (nm in params) out[[nm]] <- median(as_num(feasible_df[[nm]]), na.rm = TRUE)
  if ("ont_min_nb_reads" %in% names(out)) out[["ont_min_nb_reads"]] <- round(out[["ont_min_nb_reads"]])
  as.data.frame(out, stringsAsFactors = FALSE)
}

eval_av_cutoff <- function(ac) {
  av$pass <<- with(av,
    is.finite(covered_bases) & covered_bases >= FIXED_AVITI_COVERED_BASES &
      is.finite(reads_aligned) & reads_aligned >= FIXED_AVITI_READS_ALIGNED &
      is.finite(RPKMF) & RPKMF >= FIXED_AVITI_RPKMF &
      is.finite(genome_cov_pct) & genome_cov_pct >= as_num(ac$aviti_min_genome_coverage_pct) &
      is.finite(rpm_filtered) & rpm_filtered >= as_num(ac$aviti_min_reads_per_million_filtered) &
      (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(ac$aviti_bleed_ratio_min_pct))) &
      (is.na(nc_gc_max) | genome_cov_pct >= nc_gc_max)
  )
  on$pass <<- FALSE
  det_best_av  <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  spec_best_av <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  av_det_tp   <- sum(det_best_av$aviti_detected)
  av_det_n    <- nrow(det_best_av)
  av_samp_tmp <- tapply(det_best_av$aviti_detected,
                        paste(det_best_av$year, det_best_av$sample, sep = "||"), all)
  av_samp_tp  <- as.integer(sum(av_samp_tmp, na.rm = TRUE))
  av_samp_n   <- as.integer(length(av_samp_tmp))
  av_splike_tn <- sum(!spec_best_av$aviti_fp)
  av_splike_n  <- nrow(spec_best_av)
  av_spsamp   <- empty_sample_specificity(spec_best_av, "aviti_fp")
  data.frame(
    method = "AVITI",
    sensitivity_detection_pct = 100 * av_det_tp / av_det_n,
    sensitivity_sample_pct    = 100 * av_samp_tp / av_samp_n,
    specificity_like_pct      = 100 * av_splike_tn / av_splike_n,
    specificity_sample_pct    = av_spsamp$pct,
    det_num = av_det_tp, det_den = av_det_n,
    samp_num = av_samp_tp, samp_den = av_samp_n,
    splike_num = av_splike_tn, splike_den = av_splike_n,
    spsamp_num = av_spsamp$num, spsamp_den = av_spsamp$den,
    stringsAsFactors = FALSE
  )
}

eval_on_cutoff <- function(oc) {
  on$pass <<- with(on,
    is.finite(assignment_rank) & assignment_rank >= FIXED_ONT_ASSIGNMENT_RANK &
      is.finite(NbReads) & NbReads >= as_num(oc$ont_min_nb_reads) &
      is.finite(gc_pct) & gc_pct >= as_num(oc$ont_min_genome_coverage_pct) &
      (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= as_num(oc$ont_bleed_ratio_min_pct))) &
      (is.na(nc_gc_max) | gc_pct >= nc_gc_max)
  )
  av$pass <<- FALSE
  det_best_on  <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  spec_best_on <- score_specificity(truth_pos, QCMD_UNIVERSE, av, on, aviti_sample_map, ont_sample_map, virus_keys)
  on_det_tp   <- sum(det_best_on$ont_detected)
  on_det_n    <- nrow(det_best_on)
  on_samp_tmp <- tapply(det_best_on$ont_detected,
                        paste(det_best_on$year, det_best_on$sample, sep = "||"), all)
  on_samp_tp  <- as.integer(sum(on_samp_tmp, na.rm = TRUE))
  on_samp_n   <- as.integer(length(on_samp_tmp))
  on_splike_tn <- sum(!spec_best_on$ont_fp)
  on_splike_n  <- nrow(spec_best_on)
  on_spsamp   <- empty_sample_specificity(spec_best_on, "ont_fp")
  data.frame(
    method = "ONT",
    sensitivity_detection_pct = 100 * on_det_tp / on_det_n,
    sensitivity_sample_pct    = 100 * on_samp_tp / on_samp_n,
    specificity_like_pct      = 100 * on_splike_tn / on_splike_n,
    specificity_sample_pct    = on_spsamp$pct,
    det_num = on_det_tp, det_den = on_det_n,
    samp_num = on_samp_tp, samp_den = on_samp_n,
    splike_num = on_splike_tn, splike_den = on_splike_n,
    spsamp_num = on_spsamp$num, spsamp_den = on_spsamp$den,
    stringsAsFactors = FALSE
  )
}

white_ac <- median_cutoff_row(lhs_av, "AVITI")
white_oc <- median_cutoff_row(lhs_on, "ONT")
white_eval <- rbind(eval_av_cutoff(white_ac), eval_on_cutoff(white_oc))
metric_levels <- c("detection_sens", "sample_sens", "spec_like", "spec_sample")
metric_labels <- c(
  detection_sens = "Detection\nsensitivity",
  sample_sens    = "Sample\nsensitivity",
  spec_like      = "Specificity-like\n(detection-level)",
  spec_sample    = "Specificity\n(empty samples)"
)
long <- rbind(
  data.frame(method = eval_df$method, metric = "detection_sens", value = eval_df$sensitivity_detection_pct),
  data.frame(method = eval_df$method, metric = "sample_sens",    value = eval_df$sensitivity_sample_pct),
  data.frame(method = eval_df$method, metric = "spec_like",      value = eval_df$specificity_like_pct),
  data.frame(method = eval_df$method, metric = "spec_sample",    value = eval_df$specificity_sample_pct)
)
long$method <- factor(long$method, levels = c("AVITI", "ONT"))
long$metric <- factor(long$metric, levels = metric_levels)
best_long <- rbind(
  data.frame(method = best_eval$method, metric = "detection_sens",
             value  = best_eval$sensitivity_detection_pct,
             num    = best_eval$det_num,   den = best_eval$det_den),
  data.frame(method = best_eval$method, metric = "sample_sens",
             value  = best_eval$sensitivity_sample_pct,
             num    = best_eval$samp_num,  den = best_eval$samp_den),
  data.frame(method = best_eval$method, metric = "spec_like",
             value  = best_eval$specificity_like_pct,
             num    = best_eval$splike_num, den = best_eval$splike_den),
  data.frame(method = best_eval$method, metric = "spec_sample",
             value  = best_eval$specificity_sample_pct,
             num    = best_eval$spsamp_num, den = best_eval$spsamp_den)
)
best_long$method <- factor(best_long$method, levels = c("AVITI", "ONT"))
best_long$metric <- factor(best_long$metric, levels = metric_levels)
best_long$label  <- sprintf("%.1f%%\n(%d/%d)", best_long$value,
                             as.integer(best_long$num), as.integer(best_long$den))
best_long$label_y <- pmin(best_long$value + 2, 99)

white_long <- rbind(
  data.frame(method = white_eval$method, metric = "detection_sens",
             value  = white_eval$sensitivity_detection_pct,
             num    = white_eval$det_num,   den = white_eval$det_den),
  data.frame(method = white_eval$method, metric = "sample_sens",
             value  = white_eval$sensitivity_sample_pct,
             num    = white_eval$samp_num,  den = white_eval$samp_den),
  data.frame(method = white_eval$method, metric = "spec_like",
             value  = white_eval$specificity_like_pct,
             num    = white_eval$splike_num, den = white_eval$splike_den),
  data.frame(method = white_eval$method, metric = "spec_sample",
             value  = white_eval$specificity_sample_pct,
             num    = white_eval$spsamp_num, den = white_eval$spsamp_den)
)
white_long$method <- factor(white_long$method, levels = c("AVITI", "ONT"))
white_long$metric <- factor(white_long$metric, levels = metric_levels)
white_long$label  <- sprintf("%.1f%%\n(%d/%d)", white_long$value,
                             as.integer(white_long$num), as.integer(white_long$den))
white_long$label_y <- pmax(white_long$value - 7, 31)

cols <- c(AVITI = "#0072B2", ONT = "#D55E00")

p <- ggplot(long, aes(x = method, y = value, fill = method, color = method)) +
  geom_violin(alpha = 0.35, trim = TRUE, linewidth = 0.35) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.25, linewidth = 0.3) +
  geom_jitter(width = 0.08, alpha = 0.3, size = 0.9, color = "grey35") +
  geom_point(data = best_long, aes(x = method, y = value), inherit.aes = FALSE,
             shape = 21, fill = "red", color = "darkred", size = 3, stroke = 0.85) +
  geom_point(data = white_long, aes(x = method, y = value), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2.5, stroke = 0.8) +
  geom_text(
    data = white_long,
    aes(x = method, y = label_y, label = label),
    inherit.aes = FALSE,
    color      = "black",
    size       = 2.5,
    fontface   = "bold",
    vjust      = 1,
    lineheight = 0.9
  ) +
  geom_text(
    data = best_long,
    aes(x = method, y = label_y, label = label),
    inherit.aes = FALSE,
    color      = "darkred",
    size       = 2.8,
    fontface   = "bold",
    vjust      = 0,
    lineheight = 0.9
  ) +
  facet_wrap(~metric, nrow = 1, labeller = as_labeller(metric_labels)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(30, 108)) +
  labs(
    title = sprintf("QCMD external validation across top %d LHS candidates", N_TOP),
    subtitle = sprintf(
      "white = median cutoff-set result | red = operational cutoff result (AVITI lhs_id %d, ONT lhs_id %d) | labels show n/N",
      aviti_sel_id, ont_sel_id
    ),
    x = NULL, y = "Performance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line        = element_blank(),
    legend.position  = "none",
    strip.text       = element_text(face = "bold")
  )

print(p)
cat(sprintf("Saved: %s\n", file.path(tables_dir, "qcmd_external_validation_top100.csv")))
