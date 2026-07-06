#!/usr/bin/env Rscript

library(ggplot2)

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))
source_dir <- data_dir
csv_path   <- file.path(source_dir, "QCMD_overview.csv")
if (!file.exists(csv_path)) {
  stop(sprintf(
    "Missing input: %s",
    csv_path
  ))
}
raw <- read.csv(csv_path, skip = 1, stringsAsFactors = FALSE,
                check.names = FALSE, na.strings = c("", "NA"))
names(raw) <- c(
  "year", "sample", "virus", "vl_log10",
  "aviti_reads", "aviti_depth", "aviti_genome_cov_pct", "aviti_total_filtered_reads",
  "ont_reads", "ont_depth", "ont_genome_cov_pct", "ont_total_mapped_reads",
  "sensitivity"
)
to_num <- function(x) {
  x <- gsub("%", "", trimws(as.character(x)), fixed = TRUE)
  x[is.na(x) | x == "-"] <- "0"
  suppressWarnings(as.numeric(x))
}

raw$vl_log10                 <- to_num(raw$vl_log10)
raw$aviti_reads              <- to_num(raw$aviti_reads)
raw$aviti_genome_cov_pct     <- to_num(raw$aviti_genome_cov_pct)
raw$aviti_total_filtered_reads <- to_num(raw$aviti_total_filtered_reads)
raw$ont_reads                <- to_num(raw$ont_reads)
raw$ont_genome_cov_pct       <- to_num(raw$ont_genome_cov_pct)
raw$ont_total_mapped_reads   <- to_num(raw$ont_total_mapped_reads)
raw$aviti_rpm <- ifelse(raw$aviti_total_filtered_reads > 0,
                        raw$aviti_reads / raw$aviti_total_filtered_reads * 1e6, 0)
raw$ont_rpm   <- ifelse(raw$ont_total_mapped_reads > 0,
                        raw$ont_reads   / raw$ont_total_mapped_reads   * 1e6, 0)
qcmd <- raw

message(sprintf("QCMD rows: %d total (across %d years, %d unique samples)",
                nrow(qcmd), length(unique(qcmd$year)), length(unique(qcmd$sample))))
long_df <- rbind(
  data.frame(year = qcmd$year, sample = qcmd$sample, virus = qcmd$virus,
             method = "AVITI",
             vl_log10 = qcmd$vl_log10,
             genome_cov_pct = qcmd$aviti_genome_cov_pct,
             rpm = qcmd$aviti_rpm,
             stringsAsFactors = FALSE),
  data.frame(year = qcmd$year, sample = qcmd$sample, virus = qcmd$virus,
             method = "ONT",
             vl_log10 = qcmd$vl_log10,
             genome_cov_pct = qcmd$ont_genome_cov_pct,
             rpm = qcmd$ont_rpm,
             stringsAsFactors = FALSE)
)
long_df <- long_df[!(long_df$vl_log10 > 0 & long_df$genome_cov_pct == 0 & long_df$rpm == 0), , drop = FALSE]
long_df$method   <- factor(long_df$method, levels = c("AVITI", "ONT"))
long_df$rpm_log1 <- log10(long_df$rpm + 1)

method_colors <- c(AVITI = "#0072B2", ONT = "#D55E00")
JITTER_X <- 0.05
JITTER_Y <- 0.05
spearman_label <- function(d, x_col, y_col, method_name) {
  keep <- is.finite(d[[x_col]]) & is.finite(d[[y_col]])
  n    <- sum(keep)
  if (n < 3 || length(unique(d[[x_col]][keep])) < 2) return(NA_character_)
  r <- suppressWarnings(cor.test(d[[x_col]][keep], d[[y_col]][keep], method = "spearman"))
  sprintf("%s: rho = %.2f (p = %.3f, n = %d)", method_name,
          unname(r$estimate), r$p.value, n)
}
make_pairs <- function(y_col) {
  merge(
    long_df[long_df$method == "AVITI", c("year", "sample", "virus", "vl_log10", y_col)],
    long_df[long_df$method == "ONT",   c("year", "sample", "virus", "vl_log10", y_col)],
    by = c("year", "sample", "virus", "vl_log10"),
    suffixes = c("_aviti", "_ont")
  )
}
pairs_cov <- make_pairs("genome_cov_pct")
pairs_rpm <- make_pairs("rpm_log1")
cor_cov_aviti <- spearman_label(long_df[long_df$method == "AVITI", ], "vl_log10", "genome_cov_pct", "AVITI")
cor_cov_ont   <- spearman_label(long_df[long_df$method == "ONT",   ], "vl_log10", "genome_cov_pct", "ONT")
cor_cov_label <- paste(c(cor_cov_aviti, cor_cov_ont), collapse = "\n")

p_cov <- ggplot(long_df, aes(x = vl_log10, y = genome_cov_pct, color = method)) +
  geom_segment(data = pairs_cov,
               aes(x = vl_log10, xend = vl_log10,
                   y = genome_cov_pct_aviti, yend = genome_cov_pct_ont),
               inherit.aes = FALSE, color = "grey60", linewidth = 0.4, alpha = 0.6) +
  geom_jitter(size = 2.5, alpha = 0.85, width = JITTER_X, height = JITTER_Y) +
  geom_smooth(method = "loess", span = 1.2, se = FALSE, linewidth = 0.9, na.rm = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = cor_cov_label,
           hjust = -0.05, vjust = 1.4, size = 3.0, color = "black") +
  scale_color_manual(values = method_colors) +
  coord_cartesian(xlim = c(NA, NA), ylim = c(0, 100)) +
  labs(title = "A  Viral load vs genome coverage",
       x = "Viral load (log10 copies/ml)", y = "Genome coverage (%)") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold"))
cor_rpm_aviti <- spearman_label(long_df[long_df$method == "AVITI", ], "vl_log10", "rpm_log1", "AVITI")
cor_rpm_ont   <- spearman_label(long_df[long_df$method == "ONT",   ], "vl_log10", "rpm_log1", "ONT")
cor_rpm_label <- paste(c(cor_rpm_aviti, cor_rpm_ont), collapse = "\n")
rpm_max    <- max(long_df$rpm, na.rm = TRUE)
rpm_ticks  <- c(0, 1, 10, 100, 1000, 10000, 100000, 1e6)
rpm_ticks  <- rpm_ticks[rpm_ticks <= rpm_max * 1.2]

p_rpm <- ggplot(long_df, aes(x = vl_log10, y = rpm_log1, color = method)) +
  geom_segment(data = pairs_rpm,
               aes(x = vl_log10, xend = vl_log10,
                   y = rpm_log1_aviti, yend = rpm_log1_ont),
               inherit.aes = FALSE, color = "grey60", linewidth = 0.4, alpha = 0.6) +
  geom_jitter(size = 2.5, alpha = 0.85, width = JITTER_X, height = JITTER_Y * 0.05) +
  geom_smooth(method = "loess", span = 1.2, se = FALSE, linewidth = 0.9, na.rm = TRUE) +
  annotate("text", x = -Inf, y = Inf, label = cor_rpm_label,
           hjust = -0.05, vjust = 1.4, size = 3.0, color = "black") +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(breaks = log10(rpm_ticks + 1),
                     labels = format(rpm_ticks, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(title = "B  Viral load vs reads per million mapped reads",
       subtitle = "AVITI denominator: total filtered reads | ONT denominator: total mapped reads",
       x = "Viral load (log10 copies/ml)",
       y = "Reads per million (log10 scale)") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))
print(p_cov)
print(p_rpm)
