#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggnewscale))

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))
fig_dir <- figures_dir

aviti <- read.csv(file.path(tables_dir, "aviti_records.csv"), stringsAsFactors = FALSE, check.names = FALSE)
ont   <- read.csv(file.path(tables_dir, "ont_records.csv"),   stringsAsFactors = FALSE, check.names = FALSE)

norm_virus <- function(x) {
  x <- as.character(x)
  x[x == "HHV7"] <- "HHV-7"
  x
}
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
yn <- function(x) tolower(trimws(as.character(x))) == "yes"
sample_num <- function(x) {
  out <- suppressWarnings(as.integer(gsub("^[A-Za-z]+", "", as.character(x))))
  out[!is.finite(out)] <- 9999L
  out
}
mg_text_colour <- function(x) {
  if (is.finite(x) && x >= 50) "white" else "black"
}

aviti$virus <- norm_virus(aviti$virus)
ont$virus   <- norm_virus(ont$virus)
record_common_cols <- intersect(names(aviti), names(ont))

key_cols <- c("sample_id", "sample_id_norm", "sample_type", "virus")
keep_event <- function(df) {
  cls <- as.character(df$comparison_class)
  routine_target <- yn(df$routine_targeted) | yn(df$routinely_tested)
  reflex_done <- yn(df$reflex_tested)
  pcr_pos <- yn(df$pcr_positive_for_class) | yn(df$pcr_detected_by_cq)
  nonroutine_pcr_pos <- grepl("Non-routine PCR[+]", cls)
  routine_target | reflex_done | pcr_pos | nonroutine_pcr_pos
}

events <- unique(rbind(
  aviti[keep_event(aviti), key_cols, drop = FALSE],
  ont[keep_event(ont),     key_cols, drop = FALSE]
))
events <- events[events$sample_type %in% c("BAL", "CSF"), , drop = FALSE]
events <- events[order(events$sample_type, sample_num(events$sample_id), events$virus), , drop = FALSE]
rownames(events) <- NULL

row_for <- function(df, ev) {
  z <- df[df$sample_id == ev$sample_id &
            df$sample_type == ev$sample_type &
            df$virus == ev$virus, , drop = FALSE]
  if (nrow(z) == 0L) return(NULL)
  z[1, , drop = FALSE]
}

pcr_tile <- function(ev) {
  candidates <- list(row_for(aviti, ev), row_for(ont, ev))
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  if (length(candidates) == 0L) {
    return(data.frame(
      assay = "PCR", fill_class = "Not tested", label = "NT",
      stringsAsFactors = FALSE
    ))
  }
  z <- do.call(rbind, lapply(candidates, function(x) x[, record_common_cols, drop = FALSE]))
  z_pcr <- z[keep_event(z), , drop = FALSE]
  r <- if (nrow(z_pcr) > 0L) z_pcr[1, , drop = FALSE] else z[1, , drop = FALSE]

  cq <- as_num(r$cq_value)
  cq_lab <- if (is.finite(cq)) sprintf(" (Cq %.1f)", cq) else ""
  routine_target <- yn(r$routine_targeted) | yn(r$routinely_tested)
  reflex_done <- yn(r$reflex_tested)
  pcr_pos <- yn(r$pcr_positive_for_class) | yn(r$pcr_detected_by_cq)
  reflex_pos <- yn(r$reflex_positive)
  nonroutine_pcr_pos <- grepl("Non-routine PCR[+]", as.character(r$comparison_class))

  if (pcr_pos && routine_target && reflex_done) {
    fill_class <- if (reflex_pos) "Routine PCR+ / reflex+" else "Routine PCR+ / reflex-"
    label <- paste0("Routine+", cq_lab, "; Reflex", if (reflex_pos) "+" else "-")
  } else if (pcr_pos && routine_target) {
    fill_class <- "Routine PCR+"
    label <- paste0("Routine+", cq_lab)
  } else if (nonroutine_pcr_pos) {
    fill_class <- "Non-routine PCR+"
    label <- paste0("Non-routine+", cq_lab)
  } else if (pcr_pos && reflex_done && reflex_pos) {
    fill_class <- "Reflex PCR+"
    label <- paste0("Reflex+", cq_lab)
  } else if (pcr_pos && reflex_done) {
    fill_class <- "Reflex PCR+"
    label <- paste0("Reflex+", cq_lab)
  } else if (routine_target) {
    fill_class <- "PCR-"
    label <- "Routine-"
  } else if (reflex_done) {
    fill_class <- "PCR-"
    label <- "Reflex-"
  } else {
    fill_class <- "Not tested"
    label <- "NT"
  }
  data.frame(assay = "PCR", fill_class = fill_class, label = label, stringsAsFactors = FALSE)
}

mg_tile <- function(ev, df, assay_label) {
  r <- row_for(df, ev)
  if (is.null(r)) {
    return(data.frame(
      assay = assay_label, fill_class = "MG-", label = "MG-", high_complete = FALSE,
      mg_coverage = NA_real_, text_colour = "black",
      stringsAsFactors = FALSE
    ))
  }
  detected <- yn(r$mg_detected_best_cutoff)
  if (detected) {
    gc <- as_num(r$mg_genome_coverage_pct)
    gc_lab <- if (is.finite(gc)) sprintf(" (%.1f%%)", gc) else ""
    complete <- is.finite(gc) && gc >= 90
    fill_class <- paste0(assay_label, " MG+")
    label <- paste0("MG+", gc_lab)
  } else {
    gc <- NA_real_
    complete <- FALSE
    fill_class <- "MG-"
    label <- "MG-"
  }
  data.frame(
    assay = assay_label,
    fill_class = fill_class,
    label = label,
    high_complete = complete,
    mg_coverage = gc,
    text_colour = mg_text_colour(gc),
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", nrow(events) * 3L)
k <- 0L
for (i in seq_len(nrow(events))) {
  ev <- events[i, , drop = FALSE]
  event_label <- paste(ev$sample_id, ev$virus, sep = " | ")
  pieces <- list(
    pcr_tile(ev),
    mg_tile(ev, aviti, "AVITI"),
    mg_tile(ev, ont,   "ONT")
  )
  for (piece in pieces) {
    k <- k + 1L
    if (!"high_complete" %in% names(piece)) piece$high_complete <- FALSE
    if (!"mg_coverage" %in% names(piece)) piece$mg_coverage <- NA_real_
    if (!"text_colour" %in% names(piece)) piece$text_colour <- "black"
    rows[[k]] <- cbind(
      sample_type = ev$sample_type,
      sample_id = ev$sample_id,
      virus = ev$virus,
      event_label = event_label,
      piece,
      stringsAsFactors = FALSE
    )
  }
}
plot_df <- do.call(rbind, rows)
nt_events <- unique(plot_df$event_label[plot_df$assay == "PCR" & plot_df$label == "NT"])
if (length(nt_events) > 0L) {
  plot_df <- plot_df[!plot_df$event_label %in% nt_events, , drop = FALSE]
}

pcr_sort_group <- function(fill_class) {
  ifelse(
    fill_class %in% c("Routine PCR+", "Routine PCR+ / reflex+"), "Routine PCR+",
    ifelse(fill_class == "Routine PCR+ / reflex-", "Routine PCR+ / reflex-",
           ifelse(fill_class == "Non-routine PCR+", "Non-routine PCR+", "Other"))
  )
}

event_order <- unique(plot_df[plot_df$assay == "PCR",
                              c("sample_type", "sample_id", "virus", "event_label", "fill_class"),
                              drop = FALSE])
event_order$pcr_sort_group <- pcr_sort_group(event_order$fill_class)
event_order$pcr_sort_rank <- match(
  event_order$pcr_sort_group,
  c("Routine PCR+", "Routine PCR+ / reflex-", "Non-routine PCR+", "Other")
)
event_order$sample_sort <- sample_num(event_order$sample_id)
event_order <- event_order[order(
  match(event_order$sample_type, c("BAL", "CSF")),
  ifelse(event_order$sample_type == "BAL", event_order$pcr_sort_rank, 0L),
  event_order$sample_sort,
  event_order$virus
), , drop = FALSE]

event_levels <- rev(event_order$event_label)
plot_df$event_label <- factor(plot_df$event_label, levels = event_levels)
plot_df$assay <- factor(plot_df$assay, levels = c("PCR", "AVITI", "ONT"))
plot_df$sample_type <- factor(plot_df$sample_type, levels = c("BAL", "CSF"))
plot_df$display_order <- match(as.character(plot_df$event_label), event_order$event_label)
plot_df$pcr_sort_group <- event_order$pcr_sort_group[
  match(as.character(plot_df$event_label), event_order$event_label)
]
plot_df <- plot_df[order(plot_df$display_order, plot_df$assay), , drop = FALSE]

tile_cols <- c(
  "Routine PCR+" = scales::alpha("#12441B",0.40),
  "Routine PCR+ / reflex+" = scales::alpha("#11A354",0.90),
  "Routine PCR+ / reflex-" = scales::alpha("darkred",0.80),
  "Reflex PCR+"  = scales::alpha("#6A3D9A",0.60),
  "Non-routine PCR+" = scales::alpha("#E7258A",0.40),
  "PCR-"         = scales::alpha("#D9D9D9",0.60),
  "Not tested"   = "white",
  "MG-"          = scales::alpha("gray95",0.60)
)

pcr_plot <- plot_df[plot_df$assay == "PCR", , drop = FALSE]
mg_neg_plot <- plot_df[plot_df$assay %in% c("AVITI", "ONT") & plot_df$label == "MG-", , drop = FALSE]
aviti_plot <- plot_df[plot_df$fill_class == "AVITI MG+", , drop = FALSE]
ont_plot <- plot_df[plot_df$fill_class == "ONT MG+", , drop = FALSE]

p <- ggplot(plot_df, aes(x = assay, y = event_label)) +
  geom_tile(
    data = pcr_plot,
    aes(fill = fill_class),
    color = "white",
    linewidth = 0.65
  ) +
  geom_tile(
    data = mg_neg_plot,
    aes(fill = fill_class),
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_manual(
    name = "Result",
    values = tile_cols,
    guide = guide_legend(order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = aviti_plot,
    aes(fill = mg_coverage),
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_gradientn(
    name = "AVITI coverage (%)",
    colours = c("#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    oob = scales::squish,
    na.value = "#C6DBEF",
    guide = guide_colorbar(order = 2)
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = ont_plot,
    aes(fill = mg_coverage),
    color = "white",
    linewidth = 0.65
  ) +
  scale_fill_gradientn(
    name = "ONT coverage (%)",
    colours = c("#FEE6CE", "#FDBE85", "#FD8D3C", "#E6550D", "#A63603"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    oob = scales::squish,
    na.value = "#FDD0A2",
    guide = guide_colorbar(order = 3)
  ) +
  geom_text(
    data = pcr_plot,
    aes(label = label),
    color = "black",
    size = 2.0,
    lineheight = 0.85,
    fontface = "bold"
  ) +
  geom_text(
    data = rbind(mg_neg_plot, aviti_plot, ont_plot),
    aes(label = label, color = text_colour),
    size = 2.0,
    lineheight = 0.85,
    fontface = "bold"
  ) +
  geom_point(
    data = subset(plot_df, high_complete),
    aes(x = assay, y = event_label, shape = ">=90% genome completeness"),
    size = 2,
    stroke = 1.0,
    color = "black",
    position = position_nudge(x = 0.31),
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
  scale_color_identity(guide = "none") +
  scale_shape_manual(
    name = "Genome completeness",
    values = c(">=90% genome completeness" = 8)
  ) +
  labs(
    x = NULL,
    y = "Sample | virus",
    fill = "Result",
    title = "Potential Figure 3: clinical PCR and cutoff-filtered metagenomic concordance",
    subtitle = "PCR labels indicate routine, reflex, or non-routine testing status; AVITI/ONT tile intensity indicates genome coverage; black asterisk indicates >=90% genome completeness"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(3, "pt"),
    strip.text.y = element_text(face = "bold", angle = 0),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_text(size = 7.5, color = "black"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

out_pdf <- file.path(fig_dir, "Figure5a_pcr_metagenomics_concordance.pdf")
height_in <- max(6, min(18, 0.24 * length(unique(plot_df$event_label)) + 2.2))
ggsave(out_pdf, p, width = 10.5, height = height_in, units = "in")
write.csv(plot_df, file.path(tables_dir, "Figure5a_concordance_long.csv"), row.names = FALSE)

print(p)
message("Wrote: ", out_pdf)
message("Wrote: ", file.path(tables_dir, "Figure5a_concordance_long.csv"))
pcr_df <- plot_df[plot_df$assay == "PCR",
                  c("sample_type", "sample_id", "virus", "event_label", "fill_class", "label"),
                  drop = FALSE]
aviti_df <- plot_df[plot_df$assay == "AVITI", c("event_label", "label", "high_complete"), drop = FALSE]
names(aviti_df)[2:3] <- c("AVITI_label", "AVITI_high_complete")
ont_df <- plot_df[plot_df$assay == "ONT", c("event_label", "label", "high_complete"), drop = FALSE]
names(ont_df)[2:3] <- c("ONT_label", "ONT_high_complete")
event_df <- merge(merge(pcr_df, aviti_df, by = "event_label"), ont_df, by = "event_label")
event_df$AVITI_detected <- event_df$AVITI_label != "MG-"
event_df$ONT_detected <- event_df$ONT_label != "MG-"

nonroutine_df <- event_df[event_df$fill_class == "Non-routine PCR+", , drop = FALSE]
nonroutine_df$overlap <- ifelse(
  nonroutine_df$AVITI_detected & nonroutine_df$ONT_detected, "Both",
  ifelse(nonroutine_df$AVITI_detected, "AVITI only",
         ifelse(nonroutine_df$ONT_detected, "ONT only", "Neither"))
)
overlap_levels <- c("Both", "AVITI only", "ONT only", "Neither")
overlap_counts <- as.data.frame(table(factor(nonroutine_df$overlap, levels = overlap_levels)),
                                stringsAsFactors = FALSE)
names(overlap_counts) <- c("category", "n")
overlap_counts$category <- factor(overlap_counts$category, levels = overlap_levels)

overlap_cols <- c(
  "Both" = "#4D4D4D",
  "AVITI only" = "#0072B2",
  "ONT only" = "#D55E00",
  "Neither" = "gray80"
)

p_overlap <- ggplot(overlap_counts, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.35, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.35, fontface = "bold", size = 3.4) +
  scale_fill_manual(values = overlap_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Non-routine PCR+ recovery",
    x = NULL,
    y = "Sample-virus rows"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(3, "pt"),
    axis.text.x = element_text(face = "bold", color = "black", angle = 25, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
event_df$pcr_group <- ifelse(
  event_df$fill_class == "Non-routine PCR+", "Non-routine PCR+",
  ifelse(grepl("^Routine PCR", event_df$fill_class) &
           event_df$fill_class != "Routine PCR+ / reflex-",
         "Routine PCR+", NA_character_)
)
complete_rows <- rbind(
  data.frame(
    pcr_group = event_df$pcr_group,
    method = "AVITI",
    near_complete = event_df$AVITI_high_complete,
    stringsAsFactors = FALSE
  ),
  data.frame(
    pcr_group = event_df$pcr_group,
    method = "ONT",
    near_complete = event_df$ONT_high_complete,
    stringsAsFactors = FALSE
  )
)
complete_rows <- complete_rows[!is.na(complete_rows$pcr_group), , drop = FALSE]
complete_counts <- aggregate(near_complete ~ pcr_group + method, complete_rows, sum)
names(complete_counts)[3] <- "n"
complete_counts$pcr_group <- factor(complete_counts$pcr_group,
                                    levels = c("Routine PCR+", "Non-routine PCR+"))
complete_counts$method <- factor(complete_counts$method, levels = c("AVITI", "ONT"))

p_complete <- ggplot(complete_counts, aes(x = pcr_group, y = n, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           color = "black", linewidth = 0.35) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7),
            vjust = -0.35, fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c("AVITI" = "#0072B2", "ONT" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Near-complete genomes",
    x = NULL,
    y = "Rows with >=90% genome completeness",
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(3, "pt"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

draw_summary_panels <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
  print(p_overlap, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_complete, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::popViewport()
}

out_panels_pdf <- file.path(fig_dir, "Figure5bc_summary_panels.pdf")
pdf(out_panels_pdf, width = 4.8, height = 6.4)
draw_summary_panels()
dev.off()
draw_summary_panels()

write.csv(overlap_counts, file.path(tables_dir, "Figure5b_nonroutine_overlap.csv"),
          row.names = FALSE)
write.csv(complete_counts, file.path(tables_dir, "Figure5c_near_complete_counts.csv"),
          row.names = FALSE)

message("Wrote: ", out_panels_pdf)
message("Wrote: ", file.path(tables_dir, "Figure5b_nonroutine_overlap.csv"))
message("Wrote: ", file.path(tables_dir, "Figure5c_near_complete_counts.csv"))
