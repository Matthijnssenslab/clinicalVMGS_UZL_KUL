#!/usr/bin/env Rscript

library(ggplot2)

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))

aviti <- read.csv(file.path(tables_dir, "aviti_records.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ont   <- read.csv(file.path(tables_dir, "ont_records.csv"),   stringsAsFactors = FALSE, check.names = FALSE)
shared <- intersect(names(aviti), names(ont))
both <- rbind(aviti[, shared], ont[, shared])
both$exclude_degraded <-
  !is.na(both$reflex_tested) & !is.na(both$reflex_positive) &
  tolower(both$reflex_tested) == "yes" & tolower(both$reflex_positive) == "no"

routine <- both[
  both$sample_type == "BAL" &
    both$routinely_tested == "yes" &
    !both$exclude_degraded &
    is.finite(both$cq_value),
  ,
  drop = FALSE
]
cq_grid <- seq(floor(min(routine$cq_value)), ceiling(max(routine$cq_value)), by = 1)
curve_for <- function(method) {
  d <- routine[routine$method == method, , drop = FALSE]
  out <- do.call(rbind, lapply(cq_grid, function(ct) {
    sub <- d[d$cq_value <= ct, , drop = FALSE]
    det_tp <- sum( sub$hit_at_best); det_fn <- sum(!sub$hit_at_best)
    if (nrow(sub) > 0) {
      pred_pos <- tapply(sub$hit_at_best, sub$sample_id, all)
      smp_tp <- sum(pred_pos); smp_fn <- sum(!pred_pos)
    } else {
      smp_tp <- 0; smp_fn <- 0
    }
    rbind(
      data.frame(method = method, level = "detection", cq_cutoff = ct,
                 sensitivity_num = det_tp, sensitivity_den = det_tp + det_fn,
                 sensitivity = if (det_tp + det_fn > 0) det_tp / (det_tp + det_fn) else NA_real_,
                 stringsAsFactors = FALSE),
      data.frame(method = method, level = "sample", cq_cutoff = ct,
                 sensitivity_num = smp_tp, sensitivity_den = smp_tp + smp_fn,
                 sensitivity = if (smp_tp + smp_fn > 0) smp_tp / (smp_tp + smp_fn) else NA_real_,
                 stringsAsFactors = FALSE)
    )
  }))
  out
}
curves <- rbind(curve_for("AVITI"), curve_for("ONT"))
curves$level <- factor(curves$level, levels = c("detection", "sample"))
x_low  <- 30
x_high <- ceiling(max(curves$cq_cutoff))
if (x_high < x_low) x_low <- x_high
x_step <- if ((x_high - x_low) > 12) 2 else 1
x_breaks <- seq(x_high, x_low, by = -x_step)

p <- ggplot(curves, aes(x = cq_cutoff, y = sensitivity, color = level,
                        linetype = method, shape = method,
                        group = interaction(method, level))) +
  geom_line(linewidth = 1.35, na.rm = TRUE) +
  geom_point(size = 3.0, stroke = 0.7, na.rm = TRUE) +
  scale_x_reverse(limits = c(x_high, x_low), breaks = x_breaks) +
  scale_y_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, 0.1),
    labels = function(x) sprintf("%d", as.integer(round(100 * x)))
  ) +
  scale_color_manual(
    values = c("detection" = "#009E73", "sample" = "#CC79A7"),
    labels = c("detection" = "Detection sensitivity", "sample" = "Sample sensitivity")
  ) +
  scale_linetype_manual(values = c("AVITI" = "solid", "ONT" = "longdash")) +
  scale_shape_manual(values = c("AVITI" = 16, "ONT" = 17)) +
  labs(
    x = "Cq cutoff",
    y = "Sensitivity (%)",
    color = "Sensitivity type",
    linetype = "Method",
    shape = "Method",
    tag = "b"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black", family = "sans"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.text = element_text(color = "black", size = 16),
    axis.title.x = element_text(face = "bold", color = "black", size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", color = "black", size = 20, margin = margin(r = 8)),
    axis.ticks = element_line(color = "black", linewidth = 1.0),
    axis.ticks.length = grid::unit(0.22, "cm"),
    legend.position = "top",
    legend.title = element_text(face = "bold", color = "black", size = 12),
    legend.text = element_text(color = "black", size = 11),
    legend.box = "vertical",
    legend.key.width = grid::unit(1.1, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 8, l = 8),
    plot.tag = element_text(face = "bold", size = 42, color = "black"),
    plot.tag.position = c(0.005, 0.995)
  )

print(p)
invisible(p)
