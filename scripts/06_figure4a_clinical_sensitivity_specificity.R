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

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

MIN_SENS_FLOOR <- 0
MIN_SPEC_FLOOR <- 95

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

pick_best_real <- function(lhs_df, method_label) {
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
  pool[ord[1L], , drop = FALSE]
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
    det_denom <- tp + fn
    det_val <- if (det_denom > 0) tp / det_denom else NA_real_
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val,
      level = "detection", metric = "Sensitivity", sample_type = g,
      value = det_val, numerator = tp, denominator = det_denom,
      stringsAsFactors = FALSE
    ); k <- k + 1L

    pos_rows <- df[idx_g & routine_yes, , drop = FALSE]
    if (nrow(pos_rows) > 0) {
      pos_hit <- hit_vec[idx_g & routine_yes]
      pred_pos <- tapply(pos_hit, pos_rows$sample_id, all)
      tp_s <- sum(pred_pos)
      fn_s <- sum(!pred_pos)
      sens_s <- tp_s / (tp_s + fn_s)
      sens_num <- tp_s
      sens_den <- tp_s + fn_s
    } else {
      sens_s <- NA_real_
      sens_num <- NA_integer_
      sens_den <- NA_integer_
    }
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val,
      level = "sample", metric = "Sensitivity", sample_type = g,
      value = sens_s, numerator = sens_num, denominator = sens_den,
      stringsAsFactors = FALSE
    ); k <- k + 1L

    neg_rows <- df[idx_g & is_neg_sample, , drop = FALSE]
    if (nrow(neg_rows) > 0) {
      neg_hit <- hit_vec[idx_g & is_neg_sample]
      any_hit <- tapply(neg_hit, neg_rows$sample_id, any)
      tn <- sum(!any_hit)
      fp <- sum(any_hit)
      spec_s <- tn / (tn + fp)
      spec_num <- tn
      spec_den <- tn + fp
    } else {
      spec_s <- NA_real_
      spec_num <- NA_integer_
      spec_den <- NA_integer_
    }
    out[[k]] <- data.frame(
      method = method_label, lhs_id = lhs_id_val,
      level = "sample", metric = "Specificity", sample_type = g,
      value = spec_s, numerator = spec_num, denominator = spec_den,
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

av_sc <- prep_scored(aviti)
on_sc <- prep_scored(ont)

best_av <- pick_best_real(lhs_aviti, "AVITI")
best_on <- pick_best_real(lhs_ont, "ONT")

hit_av  <- aviti_hit(av_sc, best_av)
hit_on  <- ont_hit(on_sc, best_on)
metrics <- rbind(
  calc_metrics(av_sc, hit_av, "AVITI", as.integer(best_av$lhs_id[1])),
  calc_metrics(on_sc, hit_on, "ONT",   as.integer(best_on$lhs_id[1]))
)
metrics <- metrics[is.finite(metrics$value), , drop = FALSE]

metrics$sample_type <- factor(metrics$sample_type, levels = c("ALL", "BAL", "CSF"))
metrics$method <- factor(metrics$method, levels = c("AVITI", "ONT"))
metrics$panel_id <- factor(
  paste0(as.character(metrics$level), "__", as.character(metrics$metric)),
  levels = c("detection__Sensitivity", "sample__Sensitivity", "sample__Specificity")
)
panel_labels <- c(
  "detection__Sensitivity" = "Detection\nSensitivity",
  "sample__Sensitivity" = "Sample\nSensitivity",
  "sample__Specificity" = "Sample\nSpecificity"
)

metrics$percent_label <- sprintf("%.1f%%", 100 * metrics$value)
metrics$count_label <- ifelse(
  is.na(metrics$numerator) | is.na(metrics$denominator) | metrics$denominator <= 0,
  "",
  sprintf("%d/%d", as.integer(metrics$numerator), as.integer(metrics$denominator))
)
metrics$percent_y <- pmin(metrics$value + 0.015, 1.045)
metrics$count_y <- pmax(metrics$value * 0.50, 0.06)
lbl_key <- paste(metrics$panel_id, metrics$sample_type, sep = "||")
lbl_split <- split(metrics, lbl_key)
percent_label_rows <- lapply(lbl_split, function(d) {
  d <- d[order(as.character(d$method)), , drop = FALSE]
  uniq <- unique(d$percent_label)
  if (length(uniq) == 1L) {
    data.frame(
      panel_id = d$panel_id[1],
      sample_type = d$sample_type[1],
      method = NA_character_,
      percent_label = uniq[1],
      percent_y = max(d$percent_y, na.rm = TRUE),
      use_dodge = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      panel_id = d$panel_id,
      sample_type = d$sample_type,
      method = as.character(d$method),
      percent_label = d$percent_label,
      percent_y = d$percent_y,
      use_dodge = TRUE,
      stringsAsFactors = FALSE
    )
  }
})
percent_labels <- do.call(rbind, percent_label_rows)
percent_labels$sample_type <- factor(percent_labels$sample_type, levels = levels(metrics$sample_type))
percent_labels$panel_id <- factor(percent_labels$panel_id, levels = levels(metrics$panel_id))
percent_labels$method <- factor(percent_labels$method, levels = levels(metrics$method))

percent_labels_single <- percent_labels[!percent_labels$use_dodge, , drop = FALSE]
percent_labels_dodge <- percent_labels[percent_labels$use_dodge, , drop = FALSE]

method_colors <- c("AVITI" = "#0072B2", "ONT" = "#D55E00")

p <- ggplot(metrics, aes(x = sample_type, y = value, fill = method, color = method)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, alpha = 0.95) +
  geom_text(
    data = percent_labels_dodge,
    aes(x = sample_type, y = percent_y, label = percent_label, group = method, color = method),
    position = position_dodge(width = 0.72),
    inherit.aes = FALSE,
    vjust = 0, lineheight = 0.92, size = 3.4, fontface = "bold"
  ) +
  geom_text(
    data = percent_labels_single,
    aes(x = sample_type, y = percent_y, label = percent_label),
    inherit.aes = FALSE,
    vjust = 0, lineheight = 0.92, size = 3.4, fontface = "bold", color = "black"
  ) +
  geom_text(
    aes(label = count_label, y = count_y, group = method),
    position = position_dodge(width = 0.72),
    vjust = 0.5, angle = 90, lineheight = 0.9, size = 4.2, color = "black", fontface = "bold"
  ) +
  facet_wrap(~panel_id, nrow = 1, labeller = as_labeller(panel_labels)) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_x_discrete(labels = c("ALL" = "BAL+CSF", "BAL" = "BAL", "CSF" = "CSF")) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = function(x) sprintf("%d", as.integer(round(100 * x))),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +
  labs(
    x = NULL,
    y = "Performance (%)",
    tag = "a"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black", family = "sans"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12, color = "black", lineheight = 0.95),
    axis.text = element_text(color = "black", size = 12),
    axis.ticks = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = grid::unit(0.15, "cm"),
    axis.title.x = element_text(face = "bold", color = "black", size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", size = 16, margin = margin(r = 8)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 10, b = 8, l = 8),
    plot.tag = element_text(face = "bold", size = 42, color = "black"),
    plot.tag.position = c(0.005, 0.995)
  )

diagnose_aviti <- function(df, p) {
  rec    <- as_bool(df$record_present)
  donor  <- as_bool(df$aviti_donor_taxon)
  bleed  <- as_num(df$aviti_taxon_bleed_ratio_pct)
  gc     <- as_num(df$aviti_genome_cov_pct)
  nc     <- as_num(df$aviti_nc_genome_cov_pct_max)
  cb     <- as_num(df$aviti_covered_bases)
  ra     <- as_num(df$aviti_reads_aligned)
  rpkmf  <- as_num(df$aviti_rpkmf)
  rpm    <- as_num(df$aviti_rpm_filtered)
  gc_cut  <- as_num(p$aviti_min_genome_coverage_pct)
  rpm_cut <- as_num(p$aviti_min_reads_per_million_filtered)
  br_cut  <- as_num(p$aviti_bleed_ratio_min_pct)

  vapply(seq_len(nrow(df)), function(i) {
    r <- character(0)
    if (!rec[i]) {
      r <- c(r, "No sequencing record")
    } else {
      if (donor[i] && (is.na(bleed[i]) || bleed[i] < br_cut))
        r <- c(r, sprintf("Bleed ratio %.1f%% < cutoff %.1f%%", bleed[i], br_cut))
      if (!is.na(nc[i]) && !(is.finite(gc[i]) && gc[i] >= nc[i]))
        r <- c(r, sprintf("Genome cov %.1f%% < NC baseline %.1f%%", gc[i], nc[i]))
      if (!(is.finite(cb[i]) && cb[i] >= 150))
        r <- c(r, sprintf("Covered bases %d < 150", as.integer(cb[i])))
      if (!(is.finite(ra[i]) && ra[i] >= 1))
        r <- c(r, sprintf("Reads aligned %d < 1", as.integer(ra[i])))
      if (!(is.finite(rpkmf[i]) && rpkmf[i] >= 0))
        r <- c(r, "RPKMF < 0")
      if (!(is.finite(gc[i]) && gc[i] >= gc_cut))
        r <- c(r, sprintf("Genome cov %.1f%% < cutoff %.1f%%", gc[i], gc_cut))
      if (!(is.finite(rpm[i]) && rpm[i] >= rpm_cut))
        r <- c(r, sprintf("RPM %.2f < cutoff %.2f", rpm[i], rpm_cut))
    }
    if (length(r) == 0) "passed (unexpected)" else paste(r, collapse = "; ")
  }, character(1))
}

diagnose_ont <- function(df, p) {
  rec    <- as_bool(df$record_present)
  donor  <- as_bool(df$ont_donor_taxon)
  bleed  <- as_num(df$ont_taxon_bleed_ratio_pct)
  gc     <- as_num(df$ont_genome_cov_pct)
  nc     <- as_num(df$ont_nc_genome_cov_pct_max)
  nb     <- as_num(df$ont_nb_reads)
  rank   <- as_num(df$ont_assignment_rank)
  gc_cut  <- as_num(p$ont_min_genome_coverage_pct)
  nb_cut  <- as_num(p$ont_min_nb_reads)
  br_cut  <- as_num(p$ont_bleed_ratio_min_pct)

  vapply(seq_len(nrow(df)), function(i) {
    r <- character(0)
    if (!rec[i]) {
      r <- c(r, "No sequencing record")
    } else {
      if (donor[i] && (is.na(bleed[i]) || bleed[i] < br_cut))
        r <- c(r, sprintf("Bleed ratio %.1f%% < cutoff %.1f%%", bleed[i], br_cut))
      if (!is.na(nc[i]) && !(is.finite(gc[i]) && gc[i] >= nc[i]))
        r <- c(r, sprintf("Genome cov %.1f%% < NC baseline %.1f%%", gc[i], nc[i]))
      if (!(is.finite(nb[i]) && nb[i] >= nb_cut))
        r <- c(r, sprintf("Read count %d < cutoff %d", as.integer(nb[i]), as.integer(nb_cut)))
      if (!(is.finite(gc[i]) && gc[i] >= gc_cut))
        r <- c(r, sprintf("Genome cov %.1f%% < cutoff %.1f%%", gc[i], gc_cut))
      if (!(is.finite(rank[i]) && rank[i] >= 1))
        r <- c(r, "Assignment rank < 1")
    }
    if (length(r) == 0) "passed (unexpected)" else paste(r, collapse = "; ")
  }, character(1))
}

build_miss_table <- function(sc_df, hit_vec, diagnose_fn, best_row, method_label) {
  routine_yes <- tolower(trimws(as.character(sc_df$routinely_tested))) == "yes"
  fn_idx      <- routine_yes & !hit_vec
  fn_df       <- sc_df[fn_idx, , drop = FALSE]
  if (nrow(fn_df) == 0L) {
    return(data.frame(
      method = character(), sample_id = character(), sample_type = character(),
      virus = character(), genome_cov_pct = numeric(), nc_baseline_pct = numeric(),
      reads_aligned = numeric(), rpm_filtered = numeric(),
      gc_cutoff = numeric(), rpm_cutoff = numeric(), reads_cutoff = numeric(),
      reason_missed = character(), stringsAsFactors = FALSE
    ))
  }
  is_aviti <- method_label == "AVITI"
  data.frame(
    method          = method_label,
    sample_id       = as.character(fn_df$sample_id),
    sample_type     = as.character(fn_df$sample_type),
    virus           = as.character(fn_df$virus),
    genome_cov_pct  = round(as_num(if (is_aviti) fn_df$aviti_genome_cov_pct  else fn_df$ont_genome_cov_pct), 2),
    nc_baseline_pct = round(as_num(if (is_aviti) fn_df$aviti_nc_genome_cov_pct_max else fn_df$ont_nc_genome_cov_pct_max), 2),
    reads_aligned   = as.integer(round(as_num(if (is_aviti) fn_df$aviti_reads_aligned else fn_df$ont_nb_reads))),
    rpm_filtered    = ifelse(is_aviti, round(as_num(fn_df$aviti_rpm_filtered), 2), NA_real_),
    gc_cutoff       = round(as_num(if (is_aviti) best_row$aviti_min_genome_coverage_pct  else best_row$ont_min_genome_coverage_pct), 2),
    rpm_cutoff      = ifelse(is_aviti, round(as_num(best_row$aviti_min_reads_per_million_filtered), 2), NA_real_),
    reads_cutoff    = ifelse(!is_aviti, round(as_num(best_row$ont_min_nb_reads), 0), NA_real_),
    reason_missed   = diagnose_fn(fn_df, best_row),
    stringsAsFactors = FALSE
  )
}

miss_av  <- build_miss_table(av_sc, hit_av, diagnose_aviti, best_av, "AVITI")
miss_on  <- build_miss_table(on_sc, hit_on, diagnose_ont,   best_on, "ONT")
miss_all <- rbind(miss_av, miss_on)
miss_all <- miss_all[order(miss_all$method, miss_all$sample_type,
                           miss_all$sample_id, miss_all$virus), ]

write.csv(miss_all, file.path(tables_dir, "Figure4a_missed_detections.csv"),
          row.names = FALSE, na = "")
message(sprintf(
  "Wrote Figure4a_missed_detections.csv  (%d missed: %d AVITI, %d ONT)",
  nrow(miss_all), nrow(miss_av), nrow(miss_on)
))

print(p)
