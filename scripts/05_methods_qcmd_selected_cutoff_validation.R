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

cutoffs <- read.csv(file.path(tables_dir, "selected_cutoffs.csv"), stringsAsFactors = FALSE)
ac <- cutoffs[cutoffs$method == "AVITI", ]
oc <- cutoffs[cutoffs$method == "ONT",   ]
FIXED_AVITI_COVERED_BASES <- 150
FIXED_AVITI_READS_ALIGNED <- 1
FIXED_AVITI_RPKMF         <- 0
DONOR_MIN_READS           <- 1
DONOR_MIN_RPM             <- 1000
OUTSET_DONOR_AVITI <- c(
  "Human mastadenovirus E" = 5000000
)

ov_raw       <- read.csv(file.path(source_dir, "QCMD_overview.csv"), skip = 1, header = TRUE,
                         stringsAsFactors = FALSE)
names(ov_raw)[1:4] <- c("year", "sample", "virus", "vl_log10")
ov_raw$year   <- as.integer(suppressWarnings(as.numeric(ov_raw$year)))
ov_raw$sample <- trimws(ov_raw$sample)
ov_raw$virus  <- trimws(ov_raw$virus)
truth_pos <- ov_raw[!is.na(ov_raw$virus) & ov_raw$virus != "-", c("year", "sample", "virus")]
truth_pos$virus[grepl("Influenza A", truth_pos$virus)] <- "Influenza A"
truth_pos <- unique(truth_pos)
QCMD_VIRUS_UNIVERSE <- unique(truth_pos$virus)

aviti_sample_map <- data.frame(
  year   = c(rep(2024, 7),              rep(2025, 7)),
  sample = c(paste0("Sample", 1:6), "NC",   paste0("Sample", 1:6), "NC"),
  sid    = c("NGS_1T", "NGS_2T", "NGS_3T", "NGS_4T", "NGS_5T", "NGS_6T", "NGS_NC_T",
             "QCMD1_Q1", "QCMD2_Q2", "QCMD3_Q3", "QCMD4_Q4", "QCMD5_Q5", "QCMD6_Q6", "NC_QCMD"),
  stringsAsFactors = FALSE
)
aviti_sample_map$sample[aviti_sample_map$sample == "Sample6" & aviti_sample_map$year == 2025] <- "Sample6*"

ont_sample_map <- data.frame(
  year   = c(rep(2024, 7),              rep(2025, 6)),
  sample = c(paste0("Sample", 1:6), "NC",   paste0("Sample", 1:6)),
  sid    = c("QCMD2024_sample1", "QCMD2024_sample2", "QCMD2024_sample3",
             "QCMD2024-sample4", "QCMD2024-sample5", "QCMD2024_sample6", "QCMD2024-neg",
             "NGS_meta_25S_01", "NGS_meta_25S_02", "NGS_meta_25S_03",
             "NGS_meta_25S_04", "NGS_meta_25S_05", "NGS_meta_25S_06"),
  stringsAsFactors = FALSE
)
ont_sample_map$sample[ont_sample_map$sample == "Sample6" & ont_sample_map$year == 2025] <- "Sample6*"

virus_keys <- data.frame(
  virus = c("MPXV",            "HIV-1",                   "EVD68",
            "HSV-1",           "HSV-2",                   "CMV",
            "Influenza A",     "HadV-4"),
  aviti = c("Monkeypox virus", "immunodeficiency virus 1", "Enterovirus D",
            "Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "betaherpesvirus 5",
            "Influenza A virus",
            "Human mastadenovirus E"),
  ont   = c("Monkeypox virus", "immunodeficiency virus 1", "deconjuncti",
            "humanalpha1",     "humanalpha2",              "humanbeta5",
            "Alphainfluenzavirus",
            "exoticum"),
  stringsAsFactors = FALSE
)
any_match <- function(records, col, keyword) {
  any(grepl(keyword, records[[col]], ignore.case = TRUE))
}

av24 <- read.csv(file.path(source_dir, "2024_QCMD.csv"), stringsAsFactors = FALSE)
av25 <- read.csv(file.path(source_dir, "2025_QCMD.csv"), stringsAsFactors = FALSE)
av   <- rbind(av24, av25)

av$year           <- ifelse(av$sample_ID %in% unique(aviti_sample_map$sid[aviti_sample_map$year == 2024]),
                            2024, 2025)
av$genome_cov_pct <- av$covered_bases / av$reference_length * 100
av$rpm_filtered   <- av$reads_aligned / av$total_filtered_reads_in_sample * 1e6
av_nonNC <- av[!grepl("NC|neg", av$sample_ID, ignore.case = TRUE), ]
av_samp_total <- aggregate(reads_aligned ~ sample_ID + year + species, data = av_nonNC, FUN = sum)
av_rpm_total  <- aggregate(rpm_filtered  ~ sample_ID + year + species, data = av_nonNC, FUN = sum)
av_samp_total$rpm_sum <- av_rpm_total$rpm_filtered[
  match(paste(av_samp_total$sample_ID, av_samp_total$species),
        paste(av_rpm_total$sample_ID,  av_rpm_total$species))
]
av_max <- aggregate(cbind(reads_aligned, rpm_sum) ~ year + species, data = av_samp_total, FUN = max)
names(av_max)[3:4] <- c("taxon_max_reads", "taxon_max_rpm")
av_samp_total <- merge(av_samp_total, av_max, by = c("year", "species"), all.x = TRUE)
av_samp_total$bleed_pct <- 100 * av_samp_total$reads_aligned / av_samp_total$taxon_max_reads
av <- merge(av, av_samp_total[, c("sample_ID", "species", "taxon_max_reads", "taxon_max_rpm", "bleed_pct")],
            by = c("sample_ID", "species"), all.x = TRUE)
av$donor_taxon <- is.finite(av$taxon_max_reads) & av$taxon_max_reads >= DONOR_MIN_READS
outset_reads <- OUTSET_DONOR_AVITI[av$species]
has_outset   <- !is.na(outset_reads)
if (any(has_outset)) {
  av$bleed_pct[has_outset] <- av$bleed_pct[has_outset] *
                               av$taxon_max_reads[has_outset] /
                               outset_reads[has_outset]
  av$donor_taxon[has_outset] <- TRUE
}
av_nc <- av[grepl("NC|neg", av$sample_ID, ignore.case = TRUE) & is.finite(av$genome_cov_pct), ]
av_nc_max <- aggregate(genome_cov_pct ~ species, data = av_nc, FUN = max)
names(av_nc_max)[2] <- "nc_gc_max"
av <- merge(av, av_nc_max, by = "species", all.x = TRUE)
av$pass <- with(av,
  is.finite(covered_bases)    & covered_bases    >= FIXED_AVITI_COVERED_BASES &
  is.finite(reads_aligned)    & reads_aligned    >= FIXED_AVITI_READS_ALIGNED &
  is.finite(RPKMF)            & RPKMF            >= FIXED_AVITI_RPKMF         &
  is.finite(genome_cov_pct)   & genome_cov_pct   >= ac$aviti_min_genome_coverage_pct &
  is.finite(rpm_filtered)     & rpm_filtered      >= ac$aviti_min_reads_per_million_filtered &
  (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= ac$aviti_bleed_ratio_min_pct)) &
  (is.na(nc_gc_max) | genome_cov_pct >= nc_gc_max)
)

ont_raw <- as.data.frame(read_excel(file.path(source_dir, "Query.xls"), sheet = "Sheet0"))
qcmd_sids <- c(ont_sample_map$sid, "QCMD2024-neg")
on <- ont_raw[ont_raw$SampleName %in% qcmd_sids, ]
on$year <- ifelse(on$SampleName %in% ont_sample_map$sid[ont_sample_map$year == 2024], 2024, 2025)
mapped_tsv <- read.delim(file.path(source_dir, "mapped_reads_per_sample.tsv"), stringsAsFactors = FALSE)
on$total_mapped <- mapped_tsv$total_mapped_reads[match(on$SampleName, mapped_tsv$sample_name)]
on$rpm_mapped   <- on$NbReads / on$total_mapped * 1e6
on$gc_pct       <- as.numeric(on$GenomeCoveragePct)

al <- tolower(trimws(as.character(on$AssignmentLevel)))
on$assignment_rank <- 0L
on$assignment_rank[al == "discovery"]  <- 1L
on$assignment_rank[al == "assignment"] <- 2L
on_nonNC <- on[!grepl("neg|NC", on$SampleName, ignore.case = TRUE), ]
on_samp_total <- aggregate(NbReads ~ SampleName + year + StrainName, data = on_nonNC, FUN = sum)
on_rpm_total  <- aggregate(rpm_mapped ~ SampleName + year + StrainName, data = on_nonNC, FUN = sum)
on_samp_total$rpm_sum <- on_rpm_total$rpm_mapped[
  match(paste(on_samp_total$SampleName, on_samp_total$StrainName),
        paste(on_rpm_total$SampleName,  on_rpm_total$StrainName))
]
on_max <- aggregate(cbind(NbReads, rpm_sum) ~ year + StrainName, data = on_samp_total, FUN = max)
names(on_max)[3:4] <- c("taxon_max_reads", "taxon_max_rpm")
on_samp_total <- merge(on_samp_total, on_max, by = c("year", "StrainName"), all.x = TRUE)
on_samp_total$bleed_pct <- 100 * on_samp_total$NbReads / on_samp_total$taxon_max_reads
on <- merge(on, on_samp_total[, c("SampleName", "StrainName", "taxon_max_reads", "taxon_max_rpm", "bleed_pct")],
            by = c("SampleName", "StrainName"), all.x = TRUE)
on$donor_taxon <- is.finite(on$taxon_max_reads) & is.finite(on$taxon_max_rpm) &
                  on$taxon_max_reads >= DONOR_MIN_READS & on$taxon_max_rpm >= DONOR_MIN_RPM
on_nc <- on[on$SampleName == "QCMD2024-neg" & is.finite(on$gc_pct), ]
on_nc_max <- aggregate(gc_pct ~ StrainName, data = on_nc, FUN = max)
names(on_nc_max)[2] <- "nc_gc_max"
on <- merge(on, on_nc_max, by = "StrainName", all.x = TRUE)
on$pass <- with(on,
  is.finite(assignment_rank) & assignment_rank >= 1L                    &
  is.finite(NbReads)         & NbReads         >= oc$ont_min_nb_reads   &
  is.finite(gc_pct)          & gc_pct          >= oc$ont_min_genome_coverage_pct &
  (!donor_taxon | (is.finite(bleed_pct) & bleed_pct >= oc$ont_bleed_ratio_min_pct)) &
  (is.na(nc_gc_max) | gc_pct >= nc_gc_max)
)

score_detection <- function(truth, av_df, on_df, av_map, on_map, vkeys) {
  rows <- vector("list", nrow(truth))
  for (i in seq_len(nrow(truth))) {
    yr  <- truth$year[i]
    smp <- truth$sample[i]
    vir <- truth$virus[i]

    av_sid <- av_map$sid[av_map$year == yr & av_map$sample == smp]
    on_sid <- on_map$sid[on_map$year == yr & on_map$sample == smp]
    kw     <- vkeys[vkeys$virus == vir, , drop = FALSE]

    av_pass <- av_df[av_df$sample_ID == av_sid & av_df$pass == TRUE, ]
    on_pass <- on_df[on_df$SampleName == on_sid & on_df$pass == TRUE, ]

    aviti_hit <- if (nrow(kw) > 0 && nrow(av_pass) > 0) any_match(av_pass, "species",    kw$aviti) else FALSE
    ont_hit   <- if (nrow(kw) > 0 && nrow(on_pass) > 0) any_match(on_pass, "StrainName", kw$ont)   else FALSE

    rows[[i]] <- data.frame(year = yr, sample = smp, virus = vir,
                            aviti_detected = aviti_hit, ont_detected = ont_hit,
                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

detections <- score_detection(truth_pos, av, on, aviti_sample_map, ont_sample_map, virus_keys)

all_samples <- rbind(
  aviti_sample_map[aviti_sample_map$sample != "NC",  c("year", "sample", "sid")],
  ont_sample_map[,                                    c("year", "sample", "sid")]
)
all_samples <- unique(all_samples[, c("year", "sample")])
all_samples <- all_samples[all_samples$sample != "NC", ]

spec_rows <- list()
spec_k <- 0L
for (i in seq_len(nrow(all_samples))) {
  yr  <- all_samples$year[i]
  smp <- all_samples$sample[i]
  expected <- truth_pos$virus[truth_pos$year == yr & truth_pos$sample == smp]
  non_expected <- setdiff(QCMD_VIRUS_UNIVERSE, expected)
  if (length(non_expected) == 0) next

  av_sid <- aviti_sample_map$sid[aviti_sample_map$year == yr & aviti_sample_map$sample == smp]
  on_sid <- ont_sample_map$sid[ont_sample_map$year == yr & ont_sample_map$sample == smp]
  av_pass <- av[av$sample_ID == av_sid & av$pass == TRUE, ]
  on_pass <- on[on$SampleName == on_sid & on$pass == TRUE, ]

  for (vir in non_expected) {
    kw <- virus_keys[virus_keys$virus == vir, , drop = FALSE]
    aviti_fp <- if (nrow(kw) > 0 && nrow(av_pass) > 0) any_match(av_pass, "species",    kw$aviti) else FALSE
    ont_fp   <- if (nrow(kw) > 0 && nrow(on_pass) > 0) any_match(on_pass, "StrainName", kw$ont)   else FALSE
    spec_k <- spec_k + 1L
    spec_rows[[spec_k]] <- data.frame(
      year = yr, sample = smp, virus = vir,
      aviti_fp = aviti_fp, ont_fp = ont_fp,
      stringsAsFactors = FALSE
    )
  }
}
specificity_df <- do.call(rbind, spec_rows)

aviti_sens  <- mean(detections$aviti_detected)
ont_sens    <- mean(detections$ont_detected)
aviti_spec  <- 1 - mean(specificity_df$aviti_fp)
ont_spec    <- 1 - mean(specificity_df$ont_fp)

aviti_n_det  <- sum(detections$aviti_detected)
ont_n_det    <- sum(detections$ont_detected)
n_expected   <- nrow(detections)
aviti_n_tn   <- sum(!specificity_df$aviti_fp)
ont_n_tn     <- sum(!specificity_df$ont_fp)
n_neg        <- nrow(specificity_df)

summary_tbl <- data.frame(
  method           = c("AVITI", "ONT"),
  sensitivity      = round(c(aviti_sens, ont_sens) * 100, 1),
  specificity_like = round(c(aviti_spec, ont_spec) * 100, 1),
  n_expected       = c(n_expected, n_expected),
  n_detected       = c(aviti_n_det, ont_n_det),
  n_tn             = c(aviti_n_tn,  ont_n_tn),
  n_fp             = c(n_neg - aviti_n_tn, n_neg - ont_n_tn),
  stringsAsFactors = FALSE
)

cat("\n=== QCMD External Validation: cutoff accuracy ===\n\n")
print(summary_tbl, row.names = FALSE)

cat("\n--- Detection detail (sensitivity) ---\n")
print(detections, row.names = FALSE)

cat("\n--- Off-target hits (false positives within QCMD virus universe) ---\n")
fp_rows <- specificity_df[specificity_df$aviti_fp | specificity_df$ont_fp, ]
if (nrow(fp_rows) == 0) cat("None\n") else print(fp_rows, row.names = FALSE)

write.csv(summary_tbl, file.path(tables_dir, "qcmd_validation_results.csv"), row.names = FALSE)
write.csv(detections,  file.path(tables_dir, "qcmd_validation_detections.csv"), row.names = FALSE)
cat("\nSaved to tables/qcmd_validation_results.csv and qcmd_validation_detections.csv\n")

plot_df <- data.frame(
  method  = rep(c("AVITI", "ONT"), 2),
  metric  = c("Sensitivity", "Sensitivity", "Specificity-like", "Specificity-like"),
  value   = c(aviti_sens * 100, ont_sens * 100, aviti_spec * 100, ont_spec * 100),
  numer   = c(aviti_n_det, ont_n_det, aviti_n_tn, ont_n_tn),
  denom   = c(n_expected,  n_expected, n_neg,      n_neg),
  stringsAsFactors = FALSE
)
plot_df$metric <- factor(plot_df$metric, levels = c("Sensitivity", "Specificity-like"))
plot_df$label  <- paste0(round(plot_df$value, 1), "%\n(",
                         plot_df$numer, "/", plot_df$denom, ")")

method_colors <- c(AVITI = "#0072B2", ONT = "#D55E00")

p <- ggplot(plot_df, aes(x = metric, y = value, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "QCMD external validation — LHS-derived cutoffs",
       x = NULL, y = "Accuracy (%)", fill = "Method") +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

print(p)
