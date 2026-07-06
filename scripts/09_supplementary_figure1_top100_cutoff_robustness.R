#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))

aviti <- read.csv(file.path(tables_dir, "aviti_records.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ont   <- read.csv(file.path(tables_dir, "ont_records.csv"),   stringsAsFactors = FALSE, check.names = FALSE)
lhs_aviti <- read.csv(file.path(tables_dir, "lhs_aviti.csv"), stringsAsFactors = FALSE, check.names = FALSE)
lhs_ont   <- read.csv(file.path(tables_dir, "lhs_ont.csv"),   stringsAsFactors = FALSE, check.names = FALSE)
selected_cutoffs <- read.csv(file.path(tables_dir, "selected_cutoffs.csv"), stringsAsFactors = FALSE, check.names = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

N_TOP <- as.integer(Sys.getenv("OPERATIONAL_TOP_N", Sys.getenv("FIG2A_N_TOP", "100")))
if (!is.finite(N_TOP) || N_TOP < 1L) N_TOP <- 100L
MIN_SENS_FLOOR <- 0
MIN_SPEC_FLOOR <- 0

lhs_param_burden <- function(lhs_df, method_label) {
  params <- if (method_label == "AVITI") {
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

pick_top <- function(lhs_df, method_label) {
  sens <- as_num(lhs_df$sensitivity_pct)
  spec <- as_num(lhs_df$specificity_like_pct)
  ok <- is.finite(sens) & is.finite(spec) & sens >= MIN_SENS_FLOOR & spec >= MIN_SPEC_FLOOR
  pool <- if (any(ok)) lhs_df[ok, , drop = FALSE] else lhs_df
  spec_pool <- pool[is.finite(as_num(pool$specificity_like_pct)) &
                      abs(as_num(pool$specificity_like_pct) - 100) < 1e-9, , drop = FALSE]
  if (nrow(spec_pool) > 0L) pool <- spec_pool
  br <- lhs_param_burden(pool, method_label)
  ord <- order(
    -as_num(pool$sensitivity_pct),
    -as_num(pool$min_sens_spec_like_pct),
    -as_num(pool$avg_sens_spec_like_pct),
    -as_num(pool$specificity_like_pct),
    br,
    as_num(pool$lhs_id)
  )
  top <- pool[ord[seq_len(min(N_TOP, nrow(pool)))], , drop = FALSE]
  top$method <- method_label
  top
}

prep_scored <- function(df) {
  df <- df[df$in_panel_grid == TRUE, , drop = FALSE]
  df$exclude_degraded <-
    !is.na(df$reflex_tested) & !is.na(df$reflex_positive) &
    tolower(df$reflex_tested) == "yes" & tolower(df$reflex_positive) == "no"
  df <- df[!df$exclude_degraded, , drop = FALSE]
  df
}

calc_metrics <- function(df, hit_vec, method_label, lhs_id_val) {
  groups <- c("ALL", "BAL", "CSF")
  out <- list()
  k <- 1L

  routine_yes <- tolower(trimws(as.character(df$routinely_tested))) == "yes"
  sample_id <- as.character(df$sample_id)
  sample_type <- as.character(df$sample_type)
  sample_is_pos <- tapply(routine_yes, sample_id, any)
  neg_samples <- names(sample_is_pos)[!sample_is_pos]
  is_neg_sample <- sample_id %in% neg_samples

  for (g in groups) {
    idx_g <- rep(TRUE, nrow(df))
    if (g != "ALL") idx_g <- sample_type == g

    d_idx <- idx_g & routine_yes
    tp <- sum(hit_vec[d_idx], na.rm = TRUE)
    fn <- sum(!hit_vec[d_idx], na.rm = TRUE)
    det_den <- tp + fn
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val, level = "detection",
      metric = "Sensitivity", sample_type = g,
      value = if (det_den > 0) tp / det_den else NA_real_,
      numerator = tp, denominator = det_den,
      stringsAsFactors = FALSE
    ); k <- k + 1L

    pos_rows <- df[idx_g & routine_yes, , drop = FALSE]
    if (nrow(pos_rows) > 0) {
      pos_hit <- hit_vec[idx_g & routine_yes]
      pred_pos <- tapply(pos_hit, pos_rows$sample_id, all)
      sens_num <- sum(pred_pos)
      sens_den <- length(pred_pos)
      sens_s <- sens_num / sens_den
    } else {
      sens_s <- NA_real_
      sens_num <- NA_integer_
      sens_den <- NA_integer_
    }
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val, level = "sample",
      metric = "Sensitivity", sample_type = g, value = sens_s,
      numerator = sens_num, denominator = sens_den,
      stringsAsFactors = FALSE
    ); k <- k + 1L

    neg_rows <- df[idx_g & is_neg_sample, , drop = FALSE]
    if (nrow(neg_rows) > 0) {
      neg_hit <- hit_vec[idx_g & is_neg_sample]
      any_hit <- tapply(neg_hit, neg_rows$sample_id, any)
      spec_num <- sum(!any_hit)
      spec_den <- length(any_hit)
      spec_s <- spec_num / spec_den
    } else {
      spec_s <- NA_real_
      spec_num <- NA_integer_
      spec_den <- NA_integer_
    }
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val, level = "sample",
      metric = "Specificity", sample_type = g, value = spec_s,
      numerator = spec_num, denominator = spec_den,
      stringsAsFactors = FALSE
    ); k <- k + 1L
  }
  do.call(rbind, out)
}

aviti_hit <- function(df, p) {
  rec <- as_bool(df$record_present)
  donor <- as_bool(df$aviti_donor_taxon)
  bleed <- as_num(df$aviti_taxon_bleed_ratio_pct)
  gc <- as_num(df$aviti_genome_cov_pct)
  nc <- as_num(df$aviti_nc_genome_cov_pct_max)
  cb <- as_num(df$aviti_covered_bases)
  ra <- as_num(df$aviti_reads_aligned)
  rpkmf <- as_num(df$aviti_rpkmf)
  rpm <- as_num(df$aviti_rpm_filtered)

  bleed_pass <- rep(TRUE, nrow(df))
  idx <- rec & donor
  bleed_pass[idx] <- is.finite(bleed[idx]) & bleed[idx] >= as_num(p$aviti_bleed_ratio_min_pct)
  nc_pass <- ifelse(is.na(nc), TRUE, is.finite(gc) & gc >= nc)

  rec & bleed_pass & nc_pass &
    is.finite(cb) & cb >= 150 &
    is.finite(ra) & ra >= 1 &
    is.finite(rpkmf) & rpkmf >= 0 &
    is.finite(gc) & gc >= as_num(p$aviti_min_genome_coverage_pct) &
    is.finite(rpm) & rpm >= as_num(p$aviti_min_reads_per_million_filtered)
}

ont_hit <- function(df, p) {
  rec <- as_bool(df$record_present)
  donor <- as_bool(df$ont_donor_taxon)
  bleed <- as_num(df$ont_taxon_bleed_ratio_pct)
  gc <- as_num(df$ont_genome_cov_pct)
  nc <- as_num(df$ont_nc_genome_cov_pct_max)
  nb <- as_num(df$ont_nb_reads)
  rank <- as_num(df$ont_assignment_rank)

  bleed_pass <- rep(TRUE, nrow(df))
  idx <- rec & donor
  bleed_pass[idx] <- is.finite(bleed[idx]) & bleed[idx] >= as_num(p$ont_bleed_ratio_min_pct)
  nc_pass <- ifelse(is.na(nc), TRUE, is.finite(gc) & gc >= nc)

  rec & bleed_pass & nc_pass &
    is.finite(nb) & nb >= as_num(p$ont_min_nb_reads) &
    is.finite(gc) & gc >= as_num(p$ont_min_genome_coverage_pct) &
    is.finite(rank) & rank >= 1
}

build_dist_metrics <- function(df_scored, feasible_df, method_label, hit_fun) {
  rows <- vector("list", nrow(feasible_df))
  for (i in seq_len(nrow(feasible_df))) {
    p <- feasible_df[i, , drop = FALSE]
    hit <- hit_fun(df_scored, p)
    rows[[i]] <- calc_metrics(df_scored, hit, method_label, as.integer(p$lhs_id))
  }
  do.call(rbind, rows)
}

feas_av <- pick_top(lhs_aviti, "AVITI")
feas_on <- pick_top(lhs_ont, "ONT")
av_sc <- prep_scored(aviti)
on_sc <- prep_scored(ont)

dist <- rbind(
  build_dist_metrics(av_sc, feas_av, "AVITI", aviti_hit),
  build_dist_metrics(on_sc, feas_on, "ONT", ont_hit)
)
dist <- dist[is.finite(dist$value), , drop = FALSE]

selected_id_aviti <- as.integer(as_num(selected_cutoffs$lhs_id[selected_cutoffs$method == "AVITI"][1]))
selected_id_ont <- as.integer(as_num(selected_cutoffs$lhs_id[selected_cutoffs$method == "ONT"][1]))
selected_row_aviti <- lhs_aviti[as.integer(as_num(lhs_aviti$lhs_id)) == selected_id_aviti, , drop = FALSE]
selected_row_ont <- lhs_ont[as.integer(as_num(lhs_ont$lhs_id)) == selected_id_ont, , drop = FALSE]
if (nrow(selected_row_aviti) != 1L || nrow(selected_row_ont) != 1L) {
  stop("selected_cutoffs lhs_id not found in lhs_aviti/lhs_ont.")
}

selected_metrics <- rbind(
  calc_metrics(av_sc, aviti_hit(av_sc, selected_row_aviti), "AVITI", selected_id_aviti),
  calc_metrics(on_sc, ont_hit(on_sc, selected_row_ont), "ONT", selected_id_ont)
)
selected_metrics <- selected_metrics[is.finite(selected_metrics$value), , drop = FALSE]

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

white_row_aviti <- median_cutoff_row(feas_av, "AVITI")
white_row_ont   <- median_cutoff_row(feas_on, "ONT")
white_metrics <- rbind(
  calc_metrics(av_sc, aviti_hit(av_sc, white_row_aviti), "AVITI", NA_integer_),
  calc_metrics(on_sc, ont_hit(on_sc, white_row_ont), "ONT", NA_integer_)
)
white_metrics <- white_metrics[is.finite(white_metrics$value), , drop = FALSE]

dist$sample_type <- factor(dist$sample_type, levels = c("ALL", "BAL", "CSF"))
dist$method <- factor(dist$method, levels = c("AVITI", "ONT"))
dist$panel_id <- factor(
  paste0(as.character(dist$level), "__", as.character(dist$metric)),
  levels = c("detection__Sensitivity", "sample__Sensitivity", "sample__Specificity")
)
white_metrics$sample_type <- factor(white_metrics$sample_type, levels = levels(dist$sample_type))
white_metrics$method <- factor(white_metrics$method, levels = levels(dist$method))
white_metrics$panel_id <- factor(paste0(white_metrics$level, "__", white_metrics$metric), levels = levels(dist$panel_id))
selected_metrics$sample_type <- factor(selected_metrics$sample_type, levels = levels(dist$sample_type))
selected_metrics$method <- factor(selected_metrics$method, levels = levels(dist$method))
selected_metrics$panel_id <- factor(paste0(selected_metrics$level, "__", selected_metrics$metric), levels = levels(dist$panel_id))

white_metrics$label <- ifelse(
  is.na(white_metrics$denominator) | white_metrics$denominator <= 0,
  sprintf("%d%%", as.integer(round(100 * white_metrics$value))),
  sprintf(
    "%d%% (%d/%d)",
    as.integer(round(100 * white_metrics$value)),
    as.integer(white_metrics$numerator),
    as.integer(white_metrics$denominator)
  )
)
white_metrics$label_y <- pmin(white_metrics$value + 0.03, 1.0)

selected_metrics$label <- ifelse(
  is.na(selected_metrics$denominator) | selected_metrics$denominator <= 0,
  sprintf("%d%%", as.integer(round(100 * selected_metrics$value))),
  sprintf(
    "%d%% (%d/%d)",
    as.integer(round(100 * selected_metrics$value)),
    as.integer(selected_metrics$numerator),
    as.integer(selected_metrics$denominator)
  )
)
selected_metrics$label_y <- ifelse(
  as.character(selected_metrics$method) == "AVITI",
  pmin(selected_metrics$value + 0.055, 1.0),
  pmax(selected_metrics$value - 0.06, 0.03)
)

method_colors <- c("AVITI" = "#0072B2", "ONT" = "#D55E00")
panel_labels <- c(
  "detection__Sensitivity" = "Detection\nSensitivity",
  "sample__Sensitivity" = "Sample\nSensitivity",
  "sample__Specificity" = "Sample\nSpecificity"
)

dodge_main <- position_dodge(width = 0.72)
jitter_dodge_main <- position_jitterdodge(jitter.width = 0.07, jitter.height = 0, dodge.width = 0.72)

p <- ggplot(dist, aes(x = sample_type, y = value, fill = method, color = method)) +
  geom_boxplot(position = dodge_main, width = 0.58, outlier.shape = NA, alpha = 0.35, linewidth = 0.35) +
  geom_jitter(alpha = 0.26, size = 0.75, position = jitter_dodge_main) +
  geom_point(data = white_metrics, aes(x = sample_type, y = value, group = method),
             position = dodge_main, inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2.5, stroke = 0.8) +
  geom_text(
    data = white_metrics,
    aes(x = sample_type, y = label_y, label = label, group = method),
    position = dodge_main, inherit.aes = FALSE,
    vjust = 0, size = 2.4, lineheight = 0.9, color = "black", fontface = "bold"
  ) +
  geom_point(data = selected_metrics, aes(x = sample_type, y = value, group = method),
             position = dodge_main, inherit.aes = FALSE,
             shape = 21, fill = "darkred", color = "darkred", size = 1.7, stroke = 0.9) +
  geom_text(
    data = selected_metrics,
    aes(x = sample_type, y = label_y, label = label, group = method),
    position = dodge_main, inherit.aes = FALSE,
    vjust = 0.5, size = 2.5, lineheight = 0.9, color = "darkred", fontface = "bold"
  ) +
  facet_wrap(~panel_id, nrow = 1, drop = TRUE, labeller = as_labeller(panel_labels)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = function(x) sprintf("%d%%", as.integer(round(100 * x))),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = c("ALL" = "BAL+CSF", "BAL" = "BAL", "CSF" = "CSF")) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    title = "Supplementary cutoff check (top-N by defined LHS ranking)",
    subtitle = sprintf(
      "N=%d | white=median cutoff-set result | red=selected cutoff-set result (AVITI lhs_id=%d, ONT lhs_id=%d)",
      N_TOP, selected_id_aviti, selected_id_ont
    ),
    x = "Sample type", y = "Performance", fill = "Method", color = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5, color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(face = "bold", color = "black", size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", color = "black", size = 11, margin = margin(r = 8)),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    legend.position = "top"
  )

print(p)
