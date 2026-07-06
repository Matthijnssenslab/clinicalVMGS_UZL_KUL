#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))
N_TOP          <- as.integer(Sys.getenv("PLATEAU_N_TOP",
                                        Sys.getenv("OPERATIONAL_TOP_N", "100")))
if (!is.finite(N_TOP) || N_TOP < 1L) N_TOP <- 100L
MIN_SENS_FLOOR <- 0
MIN_SPEC_FLOOR <- 0
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
lhs_aviti <- read.csv(file.path(tables_dir, "lhs_aviti.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
lhs_ont   <- read.csv(file.path(tables_dir, "lhs_ont.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
selected  <- read.csv(file.path(tables_dir, "selected_cutoffs.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)

best_aviti <- selected[selected$method == "AVITI", , drop = FALSE]
best_ont   <- selected[selected$method == "ONT",   , drop = FALSE]

if (nrow(best_aviti) == 0L) stop("No AVITI row found in selected_cutoffs.csv")
if (nrow(best_ont)   == 0L) stop("No ONT row found in selected_cutoffs.csv")
aviti_params <- c(
  "aviti_bleed_ratio_min_pct",
  "aviti_min_genome_coverage_pct",
  "aviti_min_reads_per_million_filtered"
)
ont_params <- c(
  "ont_bleed_ratio_min_pct",
  "ont_min_genome_coverage_pct",
  "ont_min_nb_reads"
)

pretty_lbl <- c(
  aviti_bleed_ratio_min_pct            = "Bleed ratio min (%)",
  aviti_min_genome_coverage_pct        = "Min genome coverage (%)",
  aviti_min_reads_per_million_filtered = "Min RPM (filtered)",
  ont_bleed_ratio_min_pct              = "Bleed ratio min (%)",
  ont_min_genome_coverage_pct          = "Min genome coverage (%)",
  ont_min_nb_reads                     = "Min read count"
)
av_levels <- c("Bleed ratio min (%)", "Min RPM (filtered)", "Min genome coverage (%)")
on_levels <- c("Bleed ratio min (%)", "Min read count",     "Min genome coverage (%)")
pick_top_pool <- function(lhs_full, param_names, n_top = N_TOP) {
  sens <- as_num(lhs_full$sensitivity_pct)
  spec <- as_num(lhs_full$specificity_like_pct)
  ok   <- is.finite(sens) & is.finite(spec) &
          sens >= MIN_SENS_FLOOR & spec >= MIN_SPEC_FLOOR
  pool <- if (any(ok)) lhs_full[ok, , drop = FALSE] else {
    message("pick_top_pool: no rows pass floors; using full sweep.")
    lhs_full
  }
  spec100 <- pool[is.finite(as_num(pool$specificity_like_pct)) &
                    abs(as_num(pool$specificity_like_pct) - 100) < 1e-9, , drop = FALSE]
  if (nrow(spec100) > 0L) pool <- spec100
  burden <- rowSums(as.data.frame(lapply(param_names, function(nm) {
    x <- as_num(pool[[nm]]); x[!is.finite(x)] <- 0; x
  })), na.rm = TRUE)
  ord <- order(
    -as_num(pool$sensitivity_pct),
    -as_num(pool$min_sens_spec_like_pct),
    -as_num(pool$avg_sens_spec_like_pct),
    -as_num(pool$specificity_like_pct),
    burden,
    as_num(pool$lhs_id)
  )
  pool[ord[seq_len(min(n_top, nrow(pool)))], , drop = FALSE]
}

feasible_av <- pick_top_pool(lhs_aviti, aviti_params)
feasible_on <- pick_top_pool(lhs_ont,   ont_params)
message(sprintf("Feasible pool: AVITI n = %d, ONT n = %d",
                nrow(feasible_av), nrow(feasible_on)))
to_long <- function(feasible, params, lvls) {
  df <- do.call(rbind, lapply(params, function(p) {
    data.frame(parameter = pretty_lbl[[p]],
               value     = as_num(feasible[[p]]),
               stringsAsFactors = FALSE)
  }))
  df$parameter <- factor(df$parameter, levels = lvls)
  df
}

sel_long <- function(best_row, params, lvls) {
  df <- do.call(rbind, lapply(params, function(p) {
    data.frame(parameter = pretty_lbl[[p]],
               value     = as_num(best_row[[p]][1L]),
               stringsAsFactors = FALSE)
  }))
  df$parameter <- factor(df$parameter, levels = lvls)
  df
}
med_long <- function(feasible, params, lvls) {
  df <- do.call(rbind, lapply(params, function(p) {
    data.frame(parameter = pretty_lbl[[p]],
               value     = median(as_num(feasible[[p]]), na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  df$parameter <- factor(df$parameter, levels = lvls)
  df
}

av_long <- to_long(feasible_av, aviti_params, av_levels)
av_med  <- med_long(feasible_av, aviti_params, av_levels)
av_sel  <- sel_long(best_aviti,  aviti_params, av_levels)

on_long <- to_long(feasible_on, ont_params, on_levels)
on_med  <- med_long(feasible_on, ont_params, on_levels)
on_sel  <- sel_long(best_ont,    ont_params, on_levels)
av_ann <- sprintf(
  paste0("Selected cutoff (red dot)\n",
         "Bleed ratio min  :  %.2f%%\n",
         "Min genome cov.  :  %.2f%%\n",
         "Min RPM          :  %.2f"),
  as_num(best_aviti$aviti_bleed_ratio_min_pct[1L]),
  as_num(best_aviti$aviti_min_genome_coverage_pct[1L]),
  as_num(best_aviti$aviti_min_reads_per_million_filtered[1L])
)

on_ann <- sprintf(
  paste0("Selected cutoff (red dot)\n",
         "Bleed ratio min  :  %.2f%%\n",
         "Min genome cov.  :  %.2f%%\n",
         "Min read count   :  %d"),
  as_num(best_ont$ont_bleed_ratio_min_pct[1L]),
  as_num(best_ont$ont_min_genome_coverage_pct[1L]),
  as.integer(round(as_num(best_ont$ont_min_nb_reads[1L])))
)
theme_pub <- function() {
  theme_minimal(base_size = 12) %+replace%
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text         = element_text(color = "black", size = 10),
      axis.title.x      = element_text(face = "bold", color = "black", size = 11,
                                       margin = margin(t = 6)),
      axis.title.y      = element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(2, "mm"),
      plot.tag          = element_text(face = "bold", size = 14),
      plot.title        = element_text(face = "bold", color = "black", size = 12,
                                       margin = margin(b = 3)),
      plot.subtitle     = element_text(color = "grey35", size = 9),
      plot.margin       = margin(t = 5, r = 10, b = 80, l = 5)
    )
}
make_panel <- function(long_df, med_df, sel_df, ann_text,
                       fill_col, title_str, tag_str, n_feas) {
  ggplot(long_df, aes(x = .data$value, y = .data$parameter)) +
    geom_violin(fill       = fill_col,
                color      = NA,
                alpha      = 0.45,
                trim       = FALSE,
                scale      = "width") +
    geom_jitter(color      = "grey40",
                size       = 0.9,
                width      = 0,
                height     = 0.08,
                alpha      = 0.28) +
    geom_point(data        = med_df,
               aes(x = .data$value, y = .data$parameter),
               inherit.aes = FALSE,
               shape       = 21,
               fill        = "white",
               color       = "black",
               size        = 3.0,
               stroke      = 0.8) +
    geom_point(data        = sel_df,
               aes(x = .data$value, y = .data$parameter),
               inherit.aes = FALSE,
               color       = "#D62728",
               size        = 3.2) +
    annotate("label",
             x             = -Inf,
             y             = -Inf,
             label         = ann_text,
             hjust         = -0.03,
             vjust         = 1.18,
             size          = 3.0,
             lineheight    = 1.35,
             label.padding = unit(0.45, "lines"),
             label.r       = unit(0.18, "lines"),
             fill          = "white",
             color         = "black") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      tag      = tag_str,
      title    = title_str,
      subtitle = sprintf(
        "Top-%d feasible combinations | white = top-N median | red = selected cutoff",
        n_feas
      ),
      x = "Cutoff value"
    ) +
    coord_cartesian(clip = "off") +
    theme_pub()
}

pa <- make_panel(av_long, av_med, av_sel, av_ann,
                 fill_col  = "#0072B2",
                 title_str = "AVITI — LHS feasible plateau",
                 tag_str   = "A",
                 n_feas    = nrow(feasible_av))

pb <- make_panel(on_long, on_med, on_sel, on_ann,
                 fill_col  = "#D55E00",
                 title_str = "ONT — LHS feasible plateau",
                 tag_str   = "B",
                 n_feas    = nrow(feasible_on))
fig <- pa + pb + plot_layout(ncol = 2)


print(fig)
invisible(fig)
