#!/usr/bin/env Rscript

this_file <- local({
  ofiles <- vapply(sys.frames(), function(x) if (!is.null(x$ofile)) x$ofile else NA_character_, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) ofiles[length(ofiles)] else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(this_file, mustWork = TRUE)), "_paths.R"))
source_dir <- data_dir

if (!requireNamespace("lhs", quietly = TRUE)) {
  stop("Missing required R package: lhs. Install it before running this script.")
}
library(lhs)
N_AVITI <- 100000
N_ONT   <- 100000
SEED_AVITI <- 42
SEED_ONT   <- 43
MIN_SENS_FLOOR <- 0
MIN_SPEC_FLOOR <- 0
OPERATIONAL_TOP_N <- as.integer(Sys.getenv("OPERATIONAL_TOP_N", Sys.getenv("FIG2A_N_TOP", "100")))
if (!is.finite(OPERATIONAL_TOP_N) || OPERATIONAL_TOP_N < 1L) OPERATIONAL_TOP_N <- 100L

AVITI_DONOR_MIN_RPM   <- 1000
AVITI_DONOR_MIN_READS <- 100000
ONT_DONOR_MIN_RPM     <- 1000
ONT_DONOR_MIN_READS   <- 10000

FIXED_AVITI_COVERED_BASES <- 150
FIXED_AVITI_READS_ALIGNED <- 1
FIXED_AVITI_RPKMF         <- 0
FIXED_ONT_ASSIGNMENT_RANK <- 1
NC_FLOOR_OVERRIDES <- data.frame(
  species_norm = c(
    "human mastadenovirus c",
    "mastadenovirus caesari"
  ),
  sample_type           = "BAL",
  nc_genome_cov_pct_max = 10.0,
  stringsAsFactors = FALSE
)

PROGRESS_EVERY <- 1000L
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

norm_species <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x == ""] <- NA_character_
  gsub("[[:space:]]+", " ", x)
}
normalize_sample_id <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out == ""] <- NA_character_
  out <- gsub("[[:space:]-]+", "_", out)
  out <- gsub("_+", "_", out)
  fix <- function(pattern, extract, prefix) {
    idx <- grepl(pattern, out)
    if (any(idx, na.rm = TRUE)) {
      n <- as.integer(sub(extract, "\\1", out[idx], perl = TRUE))
      out[idx] <<- paste0(prefix, n)
    }
  }
  fix("^LDC_?BAL_?0*[0-9]+$", "^LDC_?BAL_?0*([0-9]+)$", "B")
  fix("^B_?0*[0-9]+$",        "^B_?0*([0-9]+)$",         "B")
  fix("^CSF_?0*[0-9]+$",      "^CSF_?0*([0-9]+)$",       "C")
  fix("^C_?0*[0-9]+$",        "^C_?0*([0-9]+)$",         "C")
  out
}
parse_ct <- function(x) {
  x <- gsub(",", ".", as.character(x), fixed = TRUE)
  n <- sub(".*?([0-9]+(?:\\.[0-9]+)?).*", "\\1", x, perl = TRUE)
  n[!grepl("[0-9]", x)] <- NA_character_
  suppressWarnings(as.numeric(n))
}
split_species_set <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0) return(character())
  toks <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  toks <- norm_species(toks)
  toks <- toks[!is.na(toks) & toks != ""]
  unique(toks)
}
normalize_ont_join_id <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out == ""] <- NA_character_
  out <- gsub("[[:space:]-]+", "_", out)
  out <- gsub("_+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  bal <- grepl("^LDC_BAL_0*[0-9]+(?:_POOL[0-9]+)?$", out, perl = TRUE)
  if (any(bal, na.rm = TRUE)) {
    n <- as.integer(sub("^LDC_BAL_0*([0-9]+)(?:_POOL[0-9]+)?$", "\\1", out[bal], perl = TRUE))
    out[bal] <- sprintf("LDC_BAL_%02d", n)
  }
  csf <- grepl("^CSF_?0*[0-9]+.*$", out, perl = TRUE)
  if (any(csf, na.rm = TRUE)) {
    n <- as.integer(sub("^CSF_?0*([0-9]+).*$", "\\1", out[csf], perl = TRUE))
    out[csf] <- paste0("CSF", n)
  }
  out
}
message("Loading raw inputs from data_input/ ...")

compare <- read.csv(
  file.path(source_dir, "CSF_BAL_compare.csv"),
  stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM"
)
names(compare) <- sub("^\xef\xbb\xbf", "", names(compare), perl = TRUE, useBytes = TRUE)

aviti_raw <- read.delim(
  file.path(source_dir, "aggregated_AVITI.tsv"),
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = c("", "NA"), quote = ""
)

ont_raw <- read.csv(
  file.path(source_dir, "combined_ONT.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
)

ont_total <- read.csv(
  file.path(source_dir, "ONT_total_mapped.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
)

ont_neg_lines <- readLines(
  file.path(source_dir, "neg_samples_CSF-BAL_GD.csv"),
  warn = FALSE, encoding = "UTF-8"
)
compare$sample_id_norm <- normalize_sample_id(compare$sample_ID)
compare$sample_type <- toupper(trimws(as.character(compare$sample_type)))
compare$sample_type[compare$sample_type == ""] <- NA_character_
need_st <- is.na(compare$sample_type)
compare$sample_type[need_st & grepl("^B[0-9]+$", compare$sample_id_norm)] <- "BAL"
compare$sample_type[need_st & grepl("^C[0-9]+$", compare$sample_id_norm)] <- "CSF"

compare$virus_ct_text       <- tolower(trimws(as.character(compare$virus_ct)))
compare$virus_ct_text[is.na(compare$virus_ct_text)] <- ""
compare$virus_ct_numeric    <- parse_ct(compare$virus_ct)
compare$cq_value            <- compare$virus_ct_numeric
compare$cq_value[grepl("log", compare$virus_ct_text, fixed = TRUE)] <- NA_real_
compare$reflex_tested_norm  <- tolower(trimws(as.character(compare$reflex_tested)))
compare$reflex_positive_norm <- tolower(trimws(as.character(compare$reflex_positive)))
compare$reflex_positive_norm[is.na(compare$reflex_positive_norm)] <- ""
compare$routinely_tested_norm <- tolower(trimws(as.character(compare$routinely_tested)))

compare$pcr_detected_by_cq <- ifelse(is.finite(compare$cq_value), "yes", "no")
compare$pcr_positive_for_class <- ifelse(
  is.finite(compare$cq_value) |
    grepl("positive", compare$virus_ct_text, fixed = TRUE) |
    grepl("log", compare$virus_ct_text, fixed = TRUE) |
    compare$reflex_positive_norm == "yes",
  "yes", "no"
)
cmp_yes <- compare[
  compare$routinely_tested_norm == "yes" &
    !is.na(compare$sample_id_norm) &
    !is.na(compare$sample_type) &
    !is.na(compare$virus),
  ,
  drop = FALSE
]
panel <- unique(cmp_yes[, c("sample_type", "virus")])
panel_species <- function(species_col) {
  out <- vector("list", nrow(panel))
  joined <- character(nrow(panel))
  for (i in seq_len(nrow(panel))) {
    vals <- cmp_yes[
      cmp_yes$sample_type == panel$sample_type[i] & cmp_yes$virus == panel$virus[i],
      species_col
    ]
    vals <- vals[!is.na(vals) & trimws(vals) != ""]
    out[[i]] <- split_species_set(vals)
    joined[i] <- if (length(vals) == 0) NA_character_ else paste(unique(vals), collapse = " | ")
  }
  list(species = out, joined = joined)
}
panel_aviti <- panel_species("equiv_AVITI_species")
panel_ont   <- panel_species("equiv_ONT_species")
compare$compare_sample_type_known <- !need_st
all_samples <- unique(compare[, c("sample_ID", "sample_id_norm",
                                  "sample_type", "compare_sample_type_known")])
all_samples <- all_samples[
  !is.na(all_samples$sample_id_norm) & !is.na(all_samples$sample_type),
  ,
  drop = FALSE
]
all_samples <- all_samples[order(-as.integer(all_samples$compare_sample_type_known)), , drop = FALSE]
all_samples <- all_samples[!duplicated(all_samples$sample_id_norm), , drop = FALSE]

grid_parts <- lapply(seq_len(nrow(all_samples)), function(i) {
  st <- all_samples$sample_type[i]
  pan_idx <- which(panel$sample_type == st)
  if (length(pan_idx) == 0) return(NULL)
  data.frame(
    sample_id        = all_samples$sample_ID[i],
    sample_id_norm   = all_samples$sample_id_norm[i],
    sample_type      = st,
    lhs_eligible     = TRUE,
    virus            = panel$virus[pan_idx],
    panel_idx        = pan_idx,
    in_panel_grid    = TRUE,
    stringsAsFactors = FALSE
  )
})
grid <- do.call(rbind, grid_parts)
grid_keys <- paste(grid$sample_id_norm, grid$virus, sep = "||")
extra <- compare[
  !is.na(compare$sample_id_norm) & !is.na(compare$sample_type) &
    !is.na(compare$virus) & tolower(compare$virus) != "negative" &
    !(paste(compare$sample_id_norm, compare$virus, sep = "||") %in% grid_keys),
  ,
  drop = FALSE
]
if (nrow(extra) > 0) {
  grid <- rbind(
    grid,
    data.frame(
      sample_id      = extra$sample_ID,
      sample_id_norm = extra$sample_id_norm,
      sample_type    = extra$sample_type,
      lhs_eligible   = FALSE,
      virus          = extra$virus,
      panel_idx      = NA_integer_,
      in_panel_grid  = FALSE,
      stringsAsFactors = FALSE
    )
  )
}
grid$row_id <- seq_len(nrow(grid))
truth_key <- paste(compare$sample_id_norm, compare$virus, sep = "||")
grid_key  <- paste(grid$sample_id_norm,    grid$virus,    sep = "||")
m <- match(grid_key, truth_key)
grid$compare_listed <- !is.na(m)
grid$routinely_tested      <- ifelse(is.na(compare$routinely_tested_norm[m]) |
                                       compare$routinely_tested_norm[m] != "yes",
                                     "no", "yes")
grid$reflex_tested         <- compare$reflex_tested_norm[m]
grid$reflex_positive       <- compare$reflex_positive_norm[m]
grid$reflex_positive[is.na(grid$reflex_positive)] <- ""
grid$virus_ct              <- as.character(compare$virus_ct[m])
grid$virus_ct_numeric      <- compare$virus_ct_numeric[m]
grid$cq_value              <- compare$cq_value[m]
grid$pcr_detected_by_cq    <- compare$pcr_detected_by_cq[m]
grid$pcr_positive_for_class <- compare$pcr_positive_for_class[m]
grid$pcr_detected_by_cq[is.na(grid$pcr_detected_by_cq)] <- "no"
grid$pcr_positive_for_class[is.na(grid$pcr_positive_for_class)] <- "no"
grid$routine_targeted <- grid$routinely_tested
grid$reflex_state <- ifelse(
  grid$routine_targeted == "yes", "routine_target",
  ifelse(grid$pcr_positive_for_class == "yes", "reflex_done", "reflex_not_done")
)
grid$ground_truth_status <- ifelse(
  grid$routine_targeted == "yes" & grid$pcr_positive_for_class == "yes", "routine_pcr_positive",
  ifelse(grid$routine_targeted == "yes" & grid$pcr_positive_for_class == "no", "routine_pcr_negative",
  ifelse(grid$routine_targeted == "no"  & grid$pcr_positive_for_class == "yes", "non_routine_reflex_positive",
                                                                                "non_routine_not_tested"))
)
grid$target_species_aviti <- ifelse(grid$in_panel_grid,
                                    panel_aviti$joined[grid$panel_idx],
                                    compare$equiv_AVITI_species[m])
grid$target_species_ont   <- ifelse(grid$in_panel_grid,
                                    panel_ont$joined[grid$panel_idx],
                                    compare$equiv_ONT_species[m])
grid$species_set_aviti <- vector("list", nrow(grid))
grid$species_set_ont   <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  if (isTRUE(grid$in_panel_grid[i])) {
    grid$species_set_aviti[[i]] <- panel_aviti$species[[grid$panel_idx[i]]]
    grid$species_set_ont  [[i]] <- panel_ont$species  [[grid$panel_idx[i]]]
  } else {
    grid$species_set_aviti[[i]] <- split_species_set(compare$equiv_AVITI_species[m[i]])
    grid$species_set_ont  [[i]] <- split_species_set(compare$equiv_ONT_species  [m[i]])
  }
}

message(sprintf("  Panel grid: %d rows (%d samples x panel viruses)",
                nrow(grid), nrow(all_samples)))
av <- data.frame(
  sample_id_norm = normalize_sample_id(aviti_raw$sample_ID),
  species        = aviti_raw$species,
  species_norm   = norm_species(aviti_raw$species),
  genus          = aviti_raw$genus,
  strain         = aviti_raw$strain,
  accession      = aviti_raw$accession,
  ref_length     = as_num(aviti_raw$reference_length),
  covered_bases  = as_num(aviti_raw$covered_bases),
  reads_aligned  = as_num(aviti_raw$reads_aligned),
  mean_coverage  = as_num(aviti_raw$mean_coverage),
  rpkmf          = as_num(aviti_raw$RPKMF),
  total_reads    = as_num(aviti_raw$total_filtered_reads_in_sample),
  stringsAsFactors = FALSE
)
sample_to_type <- unique(compare[, c("sample_id_norm", "sample_type")])
sample_to_type <- sample_to_type[!is.na(sample_to_type$sample_id_norm), , drop = FALSE]
av$sample_type <- sample_to_type$sample_type[match(av$sample_id_norm, sample_to_type$sample_id_norm)]

av$genome_cov_pct <- ifelse(is.finite(av$ref_length) & av$ref_length > 0,
                            av$covered_bases / av$ref_length * 100, NA_real_)
av$rpm_filtered   <- ifelse(is.finite(av$total_reads) & av$total_reads > 0,
                            av$reads_aligned / av$total_reads * 1e6, NA_real_)
av_bleed <- av[
  !is.na(av$sample_id_norm) & !is.na(av$sample_type) &
    !is.na(av$species_norm) & is.finite(av$reads_aligned),
  c("sample_id_norm", "sample_type", "species_norm", "reads_aligned", "rpm_filtered")
]
av_bleed_sum <- aggregate(
  cbind(reads_aligned, rpm_filtered) ~ sample_id_norm + sample_type + species_norm,
  data = av_bleed, FUN = sum
)
names(av_bleed_sum)[4:5] <- c("taxon_reads_in_sample", "taxon_rpm_in_sample")
av_bleed_max <- aggregate(
  cbind(taxon_reads_in_sample, taxon_rpm_in_sample) ~ sample_type + species_norm,
  data = av_bleed_sum, FUN = max
)
names(av_bleed_max)[3:4] <- c("taxon_max_reads_in_sample_type", "taxon_max_rpm_in_sample_type")
av_bleed_join <- merge(av_bleed_sum, av_bleed_max,
                       by = c("sample_type", "species_norm"), all.x = TRUE, sort = FALSE)
av_bleed_join$taxon_bleed_ratio_pct <- 100 *
  av_bleed_join$taxon_reads_in_sample / av_bleed_join$taxon_max_reads_in_sample_type
av <- merge(
  av,
  av_bleed_join[, c("sample_id_norm", "sample_type", "species_norm",
                    "taxon_max_reads_in_sample_type", "taxon_max_rpm_in_sample_type",
                    "taxon_bleed_ratio_pct")],
  by = c("sample_id_norm", "sample_type", "species_norm"),
  all.x = TRUE, sort = FALSE
)
av$sample_type_for_nc <- av$sample_type
av$sample_type_for_nc[grepl("^B_NC",   av$sample_id_norm)] <- "BAL"
av$sample_type_for_nc[grepl("^CSF_NC", av$sample_id_norm)] <- "CSF"

av_nc <- av[
  !is.na(av$sample_id_norm) & grepl("NC", av$sample_id_norm, fixed = TRUE) &
    !is.na(av$species_norm) & is.finite(av$genome_cov_pct) &
    !is.na(av$sample_type_for_nc),
  c("species_norm", "sample_type_for_nc", "genome_cov_pct")
]
if (nrow(av_nc) > 0) {
  av_nc_max <- aggregate(genome_cov_pct ~ species_norm + sample_type_for_nc,
                         data = av_nc, FUN = max)
  names(av_nc_max)[3] <- "nc_genome_cov_pct_max"
} else {
  av_nc_max <- data.frame(species_norm = character(),
                          sample_type_for_nc = character(),
                          nc_genome_cov_pct_max = numeric(),
                          stringsAsFactors = FALSE)
}
av$nc_genome_cov_pct_max <- av_nc_max$nc_genome_cov_pct_max[
  match(paste(av$sample_type_for_nc, av$species_norm),
        paste(av_nc_max$sample_type_for_nc, av_nc_max$species_norm))
]
av$sample_type_for_nc <- NULL
for (i in seq_len(nrow(NC_FLOOR_OVERRIDES))) {
  ov_idx <- av$sample_type == NC_FLOOR_OVERRIDES$sample_type[i] &
            av$species_norm == NC_FLOOR_OVERRIDES$species_norm[i] &
            !is.na(av$sample_type) & !is.na(av$species_norm)
  av$nc_genome_cov_pct_max[ov_idx] <- pmax(
    NC_FLOOR_OVERRIDES$nc_genome_cov_pct_max[i],
    av$nc_genome_cov_pct_max[ov_idx],
    na.rm = TRUE
  )
}

av$donor_taxon <- is.finite(av$taxon_max_reads_in_sample_type) &
  is.finite(av$taxon_max_rpm_in_sample_type) &
  av$taxon_max_reads_in_sample_type >= AVITI_DONOR_MIN_READS &
  av$taxon_max_rpm_in_sample_type   >= AVITI_DONOR_MIN_RPM
on <- data.frame(
  sample_id_norm    = normalize_sample_id(ont_raw$SAMPLE),
  ont_sample_raw    = ont_raw$SAMPLE,
  species           = ont_raw$species,
  species_norm      = norm_species(ont_raw$species),
  genus             = ont_raw$genus,
  nb_reads          = as_num(ont_raw$NbReads),
  depth_cov         = as_num(ont_raw$DepthOfCoverage),
  genome_cov_pct    = as_num(ont_raw$GenomeCoveragePct),
  assignment_lvl    = ont_raw$AssignmentLevel,
  stringsAsFactors  = FALSE
)
on$sample_type <- sample_to_type$sample_type[match(on$sample_id_norm, sample_to_type$sample_id_norm)]

al <- tolower(trimws(as.character(on$assignment_lvl)))
on$assignment_rank <- 0
on$assignment_rank[al == "discovery"]  <- 1
on$assignment_rank[al == "assignment"] <- 2
ont_total$study <- toupper(trimws(as.character(ont_total$study)))
ont_total <- ont_total[ont_total$study %in% c("BAL", "CSF"), , drop = FALSE]
ont_total$sample_join_id <- normalize_ont_join_id(ont_total$sample_name)
ont_total$total_mapped_reads <- as_num(ont_total$total_mapped_reads)
ont_total <- ont_total[
  !is.na(ont_total$sample_join_id) &
    !(ont_total$sample_join_id %in% c("BAL_NC", "CSFPANEL_NEG")) &
    is.finite(ont_total$total_mapped_reads) &
    ont_total$total_mapped_reads > 0,
  c("sample_join_id", "total_mapped_reads")
]

on$sample_join_id <- normalize_ont_join_id(on$ont_sample_raw)
on$total_mapped_reads <- ont_total$total_mapped_reads[
  match(on$sample_join_id, ont_total$sample_join_id)
]
on$rpm_mapped <- ifelse(
  is.finite(on$nb_reads) & is.finite(on$total_mapped_reads) & on$total_mapped_reads > 0,
  on$nb_reads / on$total_mapped_reads * 1e6, NA_real_
)
on_bleed <- on[
  !is.na(on$sample_id_norm) & !is.na(on$sample_type) &
    !is.na(on$species_norm) & is.finite(on$nb_reads),
  c("sample_id_norm", "sample_type", "species_norm", "nb_reads", "rpm_mapped")
]
on_bleed_sum <- aggregate(
  cbind(nb_reads, rpm_mapped) ~ sample_id_norm + sample_type + species_norm,
  data = on_bleed, FUN = sum
)
names(on_bleed_sum)[4:5] <- c("taxon_reads_in_sample", "taxon_rpm_in_sample")
on_bleed_max <- aggregate(
  cbind(taxon_reads_in_sample, taxon_rpm_in_sample) ~ sample_type + species_norm,
  data = on_bleed_sum, FUN = max
)
names(on_bleed_max)[3:4] <- c("taxon_max_reads_in_sample_type", "taxon_max_rpm_in_sample_type")
on_bleed_join <- merge(on_bleed_sum, on_bleed_max,
                       by = c("sample_type", "species_norm"), all.x = TRUE, sort = FALSE)
on_bleed_join$taxon_bleed_ratio_pct <- 100 *
  on_bleed_join$taxon_reads_in_sample / on_bleed_join$taxon_max_reads_in_sample_type
on <- merge(
  on,
  on_bleed_join[, c("sample_id_norm", "sample_type", "species_norm",
                    "taxon_max_reads_in_sample_type", "taxon_max_rpm_in_sample_type",
                    "taxon_bleed_ratio_pct")],
  by = c("sample_id_norm", "sample_type", "species_norm"),
  all.x = TRUE, sort = FALSE
)
ont_neg_parts   <- strsplit(ont_neg_lines, ",", fixed = TRUE)
strain_raw      <- vapply(ont_neg_parts, function(v) if (length(v) >= 1) v[1] else "", character(1))
gc_raw          <- vapply(ont_neg_parts, function(v) if (length(v) >= 7) v[7] else NA_character_, character(1))
is_header       <- grepl("^\\s*StrainName", strain_raw, ignore.case = TRUE)

section_type  <- rep(NA_character_, length(strain_raw))
current_type  <- NA_character_
for (i in seq_along(strain_raw)) {
  if (is_header[i]) {
    current_type <- if (grepl("BAL", strain_raw[i], ignore.case = TRUE)) "BAL" else "CSF"
  }
  section_type[i] <- current_type
}

gc_num          <- suppressWarnings(as.numeric(trimws(gc_raw)))
species_nc_norm <- norm_species(sub("\\s*/.*$", "", strain_raw, perl = TRUE))
keep            <- !is_header & !is.na(species_nc_norm) & species_nc_norm != "" &
                   is.finite(gc_num) & !is.na(section_type)
if (any(keep)) {
  on_nc_max <- aggregate(gc_num[keep] ~ species_nc_norm[keep] + section_type[keep],
                         FUN = max)
  names(on_nc_max) <- c("species_norm", "sample_type", "nc_genome_cov_pct_max")
} else {
  on_nc_max <- data.frame(species_norm = character(),
                          sample_type = character(),
                          nc_genome_cov_pct_max = numeric(),
                          stringsAsFactors = FALSE)
}
on$nc_genome_cov_pct_max <- on_nc_max$nc_genome_cov_pct_max[
  match(paste(on$sample_type, on$species_norm),
        paste(on_nc_max$sample_type, on_nc_max$species_norm))
]
for (i in seq_len(nrow(NC_FLOOR_OVERRIDES))) {
  ov_idx <- on$sample_type == NC_FLOOR_OVERRIDES$sample_type[i] &
            on$species_norm == NC_FLOOR_OVERRIDES$species_norm[i] &
            !is.na(on$sample_type) & !is.na(on$species_norm)
  on$nc_genome_cov_pct_max[ov_idx] <- pmax(
    NC_FLOOR_OVERRIDES$nc_genome_cov_pct_max[i],
    on$nc_genome_cov_pct_max[ov_idx],
    na.rm = TRUE
  )
}

on$donor_taxon <- is.finite(on$taxon_max_reads_in_sample_type) &
  is.finite(on$taxon_max_rpm_in_sample_type) &
  on$taxon_max_reads_in_sample_type >= ONT_DONOR_MIN_READS &
  on$taxon_max_rpm_in_sample_type   >= ONT_DONOR_MIN_RPM

message(sprintf("  AVITI raw records: %d | NC baseline species: %d", nrow(av), nrow(av_nc_max)))
message(sprintf("  ONT   raw records: %d | NC baseline species: %d", nrow(on), nrow(on_nc_max)))
pick_best <- function(sample_norm_vec, species_set_list, raw_df, read_col) {
  best <- vector("list", length(sample_norm_vec))
  for (i in seq_along(sample_norm_vec)) {
    sp <- species_set_list[[i]]
    if (is.na(sample_norm_vec[i]) || length(sp) == 0) next
    cand <- raw_df[
      raw_df$sample_id_norm == sample_norm_vec[i] & raw_df$species_norm %in% sp,
      ,
      drop = FALSE
    ]
    if (nrow(cand) == 0) next
    reads <- as_num(cand[[read_col]])
    ord <- order(reads, decreasing = TRUE, na.last = TRUE)
    best[[i]] <- cand[ord[1], , drop = FALSE]
  }
  best
}
take_col <- function(best_list, col) {
  out <- rep(NA, length(best_list))
  for (i in seq_along(best_list)) {
    if (!is.null(best_list[[i]]) && nrow(best_list[[i]]) > 0) {
      out[i] <- best_list[[i]][[col]][1]
    }
  }
  out
}

best_av <- pick_best(grid$sample_id_norm, grid$species_set_aviti, av, "reads_aligned")
av_panel <- data.frame(
  row_id           = grid$row_id,
  sample_id_norm   = grid$sample_id_norm,
  matched_species  = take_col(best_av, "species"),
  ref_length       = as_num(take_col(best_av, "ref_length")),
  covered_bases    = as_num(take_col(best_av, "covered_bases")),
  reads_aligned    = as_num(take_col(best_av, "reads_aligned")),
  rpkmf            = as_num(take_col(best_av, "rpkmf")),
  total_reads      = as_num(take_col(best_av, "total_reads")),
  genome_cov_pct   = as_num(take_col(best_av, "genome_cov_pct")),
  rpm_filtered     = as_num(take_col(best_av, "rpm_filtered")),
  taxon_max_reads_in_sample_type = as_num(take_col(best_av, "taxon_max_reads_in_sample_type")),
  taxon_max_rpm_in_sample_type   = as_num(take_col(best_av, "taxon_max_rpm_in_sample_type")),
  taxon_bleed_ratio_pct          = as_num(take_col(best_av, "taxon_bleed_ratio_pct")),
  donor_taxon                    = as.logical(take_col(best_av, "donor_taxon")),
  nc_genome_cov_pct_max          = as_num(take_col(best_av, "nc_genome_cov_pct_max")),
  stringsAsFactors = FALSE
)
av_panel$record_present <- !is.na(av_panel$matched_species)
av_panel$donor_taxon[is.na(av_panel$donor_taxon)] <- FALSE

best_on <- pick_best(grid$sample_id_norm, grid$species_set_ont, on, "nb_reads")
on_panel <- data.frame(
  row_id          = grid$row_id,
  sample_id_norm  = grid$sample_id_norm,
  matched_species = take_col(best_on, "species"),
  nb_reads        = as_num(take_col(best_on, "nb_reads")),
  total_mapped_reads = as_num(take_col(best_on, "total_mapped_reads")),
  rpm_mapped      = as_num(take_col(best_on, "rpm_mapped")),
  genome_cov_pct  = as_num(take_col(best_on, "genome_cov_pct")),
  depth_cov       = as_num(take_col(best_on, "depth_cov")),
  assignment_lvl  = take_col(best_on, "assignment_lvl"),
  assignment_rank = as_num(take_col(best_on, "assignment_rank")),
  taxon_max_reads_in_sample_type = as_num(take_col(best_on, "taxon_max_reads_in_sample_type")),
  taxon_max_rpm_in_sample_type   = as_num(take_col(best_on, "taxon_max_rpm_in_sample_type")),
  taxon_bleed_ratio_pct          = as_num(take_col(best_on, "taxon_bleed_ratio_pct")),
  donor_taxon                    = as.logical(take_col(best_on, "donor_taxon")),
  nc_genome_cov_pct_max          = as_num(take_col(best_on, "nc_genome_cov_pct_max")),
  stringsAsFactors = FALSE
)
on_panel$record_present <- !is.na(on_panel$matched_species)
on_panel$donor_taxon[is.na(on_panel$donor_taxon)] <- FALSE
included_for_scoring <- function(reflex_pos, routine) {
  rp <- tolower(trimws(as.character(reflex_pos)))
  base <- is.na(rp) | rp != "no"
  base | (tolower(trimws(as.character(routine))) == "yes")
}
truth_positive <- function(routine) {
  tolower(trimws(as.character(routine))) == "yes"
}

include_idx <- included_for_scoring(grid$reflex_positive, grid$routinely_tested) & grid$lhs_eligible
truth_idx   <- truth_positive(grid$routinely_tested)
truth_idx[is.na(truth_idx)] <- FALSE
scored_sids  <- unique(grid$sample_id_norm[include_idx])
pos_sids     <- unique(grid$sample_id_norm[include_idx & truth_idx])
neg_sids     <- scored_sids[!scored_sids %in% pos_sids]
neg_inc_idx  <- include_idx & (grid$sample_id_norm %in% neg_sids)
neg_inc_sids <- grid$sample_id_norm[neg_inc_idx]
message(sprintf("  Negative samples (0 truth-pos rows): %d  |  Positive samples: %d",
                length(neg_sids), length(pos_sids)))
draw_lhs <- function(n, ranges, seed) {
  set.seed(seed)
  m <- randomLHS(n, length(ranges))
  out <- as.data.frame(m)
  names(out) <- names(ranges)
  for (j in seq_along(ranges)) {
    pd <- ranges[[j]]
    val <- pd$min + m[, j] * (pd$max - pd$min)
    if (isTRUE(pd$integer)) val <- round(val)
    out[[j]] <- val
  }
  out
}

aviti_ranges <- list(
  aviti_bleed_ratio_min_pct            = list(min = 0, max = 10, integer = FALSE),
  aviti_min_genome_coverage_pct        = list(min = 0, max = 50, integer = FALSE),
  aviti_min_reads_per_million_filtered = list(min = 0, max = 20, integer = FALSE)
)
ont_ranges <- list(
  ont_bleed_ratio_min_pct      = list(min = 0, max = 10,  integer = FALSE),
  ont_min_nb_reads             = list(min = 0, max = 100, integer = TRUE),
  ont_min_genome_coverage_pct  = list(min = 0, max = 50,  integer = FALSE)
)
aviti_lhs <- draw_lhs(N_AVITI, aviti_ranges, SEED_AVITI)
ont_lhs   <- draw_lhs(N_ONT,   ont_ranges,   SEED_ONT)
aviti_hit <- function(p) {
  donor_idx <- av_panel$record_present & av_panel$donor_taxon
  bleed_pass <- rep(TRUE, nrow(av_panel))
  bleed_pass[donor_idx] <- !is.na(av_panel$taxon_bleed_ratio_pct[donor_idx]) &
    av_panel$taxon_bleed_ratio_pct[donor_idx] >= p$aviti_bleed_ratio_min_pct
  nc_pass <- ifelse(
    is.na(av_panel$nc_genome_cov_pct_max), TRUE,
    is.finite(av_panel$genome_cov_pct) &
      av_panel$genome_cov_pct >= av_panel$nc_genome_cov_pct_max
  )
  cb_pass    <- is.finite(av_panel$covered_bases)  & av_panel$covered_bases  >= FIXED_AVITI_COVERED_BASES
  ra_pass    <- is.finite(av_panel$reads_aligned)  & av_panel$reads_aligned  >= FIXED_AVITI_READS_ALIGNED
  rpkmf_pass <- is.finite(av_panel$rpkmf)          & av_panel$rpkmf          >= FIXED_AVITI_RPKMF
  cov_pass   <- is.finite(av_panel$genome_cov_pct) & av_panel$genome_cov_pct >= p$aviti_min_genome_coverage_pct
  rpm_pass   <- is.finite(av_panel$rpm_filtered)   & av_panel$rpm_filtered   >= p$aviti_min_reads_per_million_filtered
  av_panel$record_present & bleed_pass & nc_pass & cb_pass & ra_pass & rpkmf_pass & cov_pass & rpm_pass
}
ont_hit <- function(p) {
  donor_idx <- on_panel$record_present & on_panel$donor_taxon
  bleed_pass <- rep(TRUE, nrow(on_panel))
  bleed_pass[donor_idx] <- !is.na(on_panel$taxon_bleed_ratio_pct[donor_idx]) &
    on_panel$taxon_bleed_ratio_pct[donor_idx] >= p$ont_bleed_ratio_min_pct
  nc_pass <- ifelse(
    is.na(on_panel$nc_genome_cov_pct_max), TRUE,
    is.finite(on_panel$genome_cov_pct) &
      on_panel$genome_cov_pct >= on_panel$nc_genome_cov_pct_max
  )
  reads_pass <- is.finite(on_panel$nb_reads)        & on_panel$nb_reads        >= p$ont_min_nb_reads
  cov_pass   <- is.finite(on_panel$genome_cov_pct)  & on_panel$genome_cov_pct  >= p$ont_min_genome_coverage_pct
  rank_pass  <- is.finite(on_panel$assignment_rank) & on_panel$assignment_rank >= FIXED_ONT_ASSIGNMENT_RANK
  on_panel$record_present & bleed_pass & nc_pass & reads_pass & cov_pass & rank_pass
}

score <- function(hit) {
  hit[is.na(hit)] <- FALSE
  pos <- include_idx & truth_idx
  tp <- sum(hit[pos]); fn <- sum(!hit[pos])
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  n_neg <- length(neg_sids)
  if (n_neg > 0) {
    fp_s <- sum(tapply(hit[neg_inc_idx], neg_inc_sids, any), na.rm = TRUE)
    spec <- (n_neg - fp_s) / n_neg
  } else {
    fp_s <- 0
    spec <- NA_real_
  }
  neg_row_idx <- include_idx & !truth_idx
  ot <- sum(hit[neg_row_idx]); ca <- sum(!hit[neg_row_idx])

  c(TP = tp, FN = fn, OFF_TARGET = ot, CORRECT_ABSENCE = ca,
    NEG_SAMPLES_N = n_neg, NEG_SAMPLES_FP = fp_s,
    sensitivity = sens, specificity_like = spec,
    avg = if (!is.na(sens) && !is.na(spec)) (sens + spec) / 2 else NA_real_,
    minss = if (!is.na(sens) && !is.na(spec)) min(sens, spec) else NA_real_)
}

run_lhs <- function(lhs_df, hit_fun, label) {
  n <- nrow(lhs_df)
  results <- vector("list", n)
  for (i in seq_len(n)) {
    if (i == 1L || i == n || (i %% PROGRESS_EVERY) == 0L) {
      message(sprintf("  %s LHS: %d/%d", label, i, n))
    }
    s <- score(hit_fun(as.list(lhs_df[i, , drop = FALSE])))
    results[[i]] <- data.frame(
      lhs_id              = i,
      TP                  = unname(s["TP"]),
      FN                  = unname(s["FN"]),
      OFF_TARGET          = unname(s["OFF_TARGET"]),
      CORRECT_ABSENCE     = unname(s["CORRECT_ABSENCE"]),
      neg_samples_n       = unname(s["NEG_SAMPLES_N"]),
      neg_samples_fp      = unname(s["NEG_SAMPLES_FP"]),
      sensitivity_pct     = round(unname(s["sensitivity"])      * 100, 2),
      specificity_like_pct = round(unname(s["specificity_like"]) * 100, 2),
      avg_sens_spec_like_pct = round(unname(s["avg"])           * 100, 2),
      min_sens_spec_like_pct = round(unname(s["minss"])         * 100, 2),
      stringsAsFactors    = FALSE
    )
    for (nm in names(lhs_df)) results[[i]][[nm]] <- lhs_df[[nm]][i]
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

message(sprintf("Running AVITI LHS sweep (%d combinations)...", N_AVITI))
lhs_aviti <- run_lhs(aviti_lhs, aviti_hit, "AVITI")
message(sprintf("Running ONT LHS sweep (%d combinations)...", N_ONT))
lhs_ont   <- run_lhs(ont_lhs,   ont_hit,   "ONT")
pick_top_pool <- function(lhs_full, param_names, n_top) {
  sweep <- lhs_full
  sens <- as_num(sweep$sensitivity_pct)
  spec <- as_num(sweep$specificity_like_pct)
  ok <- is.finite(sens) & is.finite(spec) &
    sens >= MIN_SENS_FLOOR & spec >= MIN_SPEC_FLOOR
  pool <- if (any(ok)) sweep[ok, , drop = FALSE] else {
    message("pick_top_pool: no LHS rows pass sensitivity/spec floors; using full sweep table.")
    sweep
  }
  spec100 <- pool[is.finite(as_num(pool$specificity_like_pct)) &
                    abs(as_num(pool$specificity_like_pct) - 100) < 1e-9, , drop = FALSE]
  if (nrow(spec100) > 0L) pool <- spec100
  burden <- rowSums(as.data.frame(lapply(param_names, function(nm) {
    x <- as_num(pool[[nm]])
    x[!is.finite(x)] <- 0
    x
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

pool_aviti <- pick_top_pool(lhs_aviti, names(aviti_ranges), OPERATIONAL_TOP_N)
pool_ont   <- pick_top_pool(lhs_ont,   names(ont_ranges),   OPERATIONAL_TOP_N)

write.csv(lhs_aviti, file.path(tables_dir, "lhs_aviti.csv"), row.names = FALSE, na = "")
write.csv(lhs_ont,   file.path(tables_dir, "lhs_ont.csv"),   row.names = FALSE, na = "")
message(sprintf(
  "Wrote lhs_aviti.csv (%d rows), lhs_ont.csv (%d rows); top-N pool N=%d",
  nrow(lhs_aviti), nrow(lhs_ont),
  OPERATIONAL_TOP_N
))
pick_best_real <- function(lhs_df, param_names) {
  sens <- as_num(lhs_df$sensitivity_pct)
  spec <- as_num(lhs_df$specificity_like_pct)
  ok   <- is.finite(sens) & is.finite(spec) & sens >= MIN_SENS_FLOOR & spec >= MIN_SPEC_FLOOR
  pool <- if (any(ok)) lhs_df[ok, , drop = FALSE] else lhs_df
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
  pool[ord[1L], , drop = FALSE]
}

best_aviti <- pick_best_real(lhs_aviti, names(aviti_ranges))
best_ont   <- pick_best_real(lhs_ont,   names(ont_ranges))
feasible_aviti <- pool_aviti
feasible_ont <- pool_ont

plateau_summary <- function(method_label, feasible_df, best_row, params) {
  score_cols <- c("sensitivity_pct", "specificity_like_pct",
                  "avg_sens_spec_like_pct", "min_sens_spec_like_pct")
  out <- data.frame(method = method_label, n_feasible = nrow(feasible_df), stringsAsFactors = FALSE)
  for (nm in score_cols) {
    x <- as_num(feasible_df[[nm]])
    out[[paste0(nm, "_min")]] <- round(min(x, na.rm = TRUE), 2)
    out[[paste0(nm, "_median")]] <- round(median(x, na.rm = TRUE), 2)
    out[[paste0(nm, "_max")]] <- round(max(x, na.rm = TRUE), 2)
    out[[paste0("selected_", nm)]] <- round(as_num(best_row[[nm]]), 2)
  }
  for (nm in params) {
    x <- as_num(feasible_df[[nm]])
    out[[paste0(nm, "_min")]] <- round(min(x, na.rm = TRUE), 6)
    out[[paste0(nm, "_median")]] <- round(median(x, na.rm = TRUE), 6)
    out[[paste0(nm, "_max")]] <- round(max(x, na.rm = TRUE), 6)
    out[[paste0("selected_", nm)]] <- round(as_num(best_row[[nm]]), 6)
  }
  out
}

summary_aviti <- plateau_summary("AVITI", feasible_aviti, best_aviti, names(aviti_ranges))
summary_ont   <- plateau_summary("ONT",   feasible_ont,   best_ont,   names(ont_ranges))
write.csv(summary_aviti, file.path(tables_dir, "lhs_aviti_plateau_summary.csv"), row.names = FALSE, na = "")
write.csv(summary_ont,   file.path(tables_dir, "lhs_ont_plateau_summary.csv"),   row.names = FALSE, na = "")
message(sprintf("Wrote lhs_aviti_plateau_summary.csv (n=%d), lhs_ont_plateau_summary.csv (n=%d)",
                nrow(feasible_aviti), nrow(feasible_ont)))

selected_cutoffs <- data.frame(
  method = c("AVITI", "ONT"),
  lhs_id = c(best_aviti$lhs_id, best_ont$lhs_id),
  avg_sens_spec_like_pct  = c(best_aviti$avg_sens_spec_like_pct,  best_ont$avg_sens_spec_like_pct),
  specificity_like_pct    = c(best_aviti$specificity_like_pct,    best_ont$specificity_like_pct),
  sensitivity_pct         = c(best_aviti$sensitivity_pct,         best_ont$sensitivity_pct),
  aviti_bleed_ratio_min_pct            = c(best_aviti$aviti_bleed_ratio_min_pct, NA_real_),
  aviti_min_genome_coverage_pct        = c(best_aviti$aviti_min_genome_coverage_pct, NA_real_),
  aviti_min_reads_per_million_filtered = c(best_aviti$aviti_min_reads_per_million_filtered, NA_real_),
  ont_bleed_ratio_min_pct              = c(NA_real_, best_ont$ont_bleed_ratio_min_pct),
  ont_min_nb_reads                     = c(NA_real_, best_ont$ont_min_nb_reads),
  ont_min_genome_coverage_pct          = c(NA_real_, best_ont$ont_min_genome_coverage_pct),
  stringsAsFactors = FALSE
)
write.csv(selected_cutoffs,
          file.path(tables_dir, "selected_cutoffs.csv"),
          row.names = FALSE, na = "")
message(sprintf(
  "AVITI best: lhs_id=%d  avg=%.2f  sens=%.2f  spec=%.2f",
  best_aviti$lhs_id, best_aviti$avg_sens_spec_like_pct,
  best_aviti$sensitivity_pct, best_aviti$specificity_like_pct
))
message(sprintf(
  "ONT   best: lhs_id=%d  avg=%.2f  sens=%.2f  spec=%.2f",
  best_ont$lhs_id, best_ont$avg_sens_spec_like_pct,
  best_ont$sensitivity_pct, best_ont$specificity_like_pct
))
hit_av <- aviti_hit(as.list(best_aviti[1, , drop = FALSE]))
hit_on <- ont_hit  (as.list(best_ont  [1, , drop = FALSE]))

cov_av <- ifelse(hit_av, av_panel$genome_cov_pct, NA_real_)
cov_on <- ifelse(hit_on, on_panel$genome_cov_pct, NA_real_)
classify <- function(routine, pcr_pos, hit) {
  ifelse(routine == "yes" & pcr_pos == "yes" & hit, "Routine PCR+ & MG+",
  ifelse(routine == "yes" & pcr_pos == "yes" & !hit, "Routine PCR+ & MG-",
  ifelse(routine == "yes" & pcr_pos == "no"  & hit, "Routine PCR- & MG+",
  ifelse(routine == "yes" & pcr_pos == "no"  & !hit, "Routine PCR- & MG-",
  ifelse(routine == "no"  & pcr_pos == "yes" & hit, "Non-routine PCR+ & MG+",
  ifelse(routine == "no"  & pcr_pos == "yes" & !hit, "Non-routine PCR+ & MG-",
  ifelse(routine == "no"  & pcr_pos == "no"  & hit, "Non-routine MG+ only",
                                                    "Non-routine not tested")))))))
}

virus_family_for <- function(virus) {
  key <- gsub("[^a-z0-9]+", "", tolower(trimws(virus)))
  out <- rep("Other", length(key))
  out[key %in% c("adenovirus")]                                                  <- "Adenoviridae"
  out[key %in% c("cmv", "ebv", "hhv6", "hhv7", "hsv1", "hsv2", "vzv")]           <- "Herpesviridae"
  out[key %in% c("echovirus30", "enterorhino", "enterovirusd68")]                <- "Picornaviridae"
  out[key %in% c("hpyv6", "polyomavirus")]                                       <- "Polyomaviridae"
  out[key %in% c("hmpv", "parainfluenzatype1", "parainfluenzatype2",
                 "parainfluenzatype3", "rsv", "rsvb")]                           <- "Paramyxoviridae"
  out[key %in% c("influenzaa")]                                                  <- "Orthomyxoviridae"
  out[key %in% c("coronanel63", "coronaoc43", "sarscov2")]                       <- "Coronaviridae"
  out
}
build_records <- function(method_label, lhs_id_best, target_species_col,
                          panel_df, hit_vec, cov_vec) {
  d <- data.frame(
    method                  = method_label,
    lhs_id_best             = lhs_id_best,
    sample_id               = grid$sample_id,
    sample_id_norm          = grid$sample_id_norm,
    sample_type             = grid$sample_type,
    virus                   = grid$virus,
    in_panel_grid           = grid$in_panel_grid,
    compare_listed          = grid$compare_listed,
    target_species          = grid[[target_species_col]],
    matched_species         = panel_df$matched_species,
    routinely_tested        = grid$routinely_tested,
    routine_targeted        = grid$routine_targeted,
    reflex_tested           = grid$reflex_tested,
    reflex_positive         = grid$reflex_positive,
    virus_ct                = grid$virus_ct,
    virus_ct_numeric        = grid$virus_ct_numeric,
    cq_value                = grid$cq_value,
    pcr_detected_by_cq      = grid$pcr_detected_by_cq,
    pcr_positive_for_class  = grid$pcr_positive_for_class,
    reflex_state            = grid$reflex_state,
    ground_truth_status     = grid$ground_truth_status,
    record_present          = panel_df$record_present,
    hit_at_best             = hit_vec,
    mg_detected_best_cutoff = ifelse(hit_vec, "yes", "no"),
    mg_genome_coverage_pct  = cov_vec,
    stringsAsFactors        = FALSE
  )
  d$comparison_class <- classify(d$routine_targeted, d$pcr_positive_for_class, d$hit_at_best)
  d$virus_family    <- virus_family_for(d$virus)
  d$plot_include    <- d$compare_listed &
                       d$comparison_class != "Non-routine not tested" &
                       d$comparison_class != "Routine PCR- & MG-"
  d
}

aviti_records <- build_records(
  "AVITI", best_aviti$lhs_id, "target_species_aviti",
  av_panel, hit_av, cov_av
)
aviti_records$aviti_reads_aligned                   <- av_panel$reads_aligned
aviti_records$aviti_covered_bases                   <- av_panel$covered_bases
aviti_records$aviti_ref_length                      <- av_panel$ref_length
aviti_records$aviti_genome_cov_pct                  <- av_panel$genome_cov_pct
aviti_records$aviti_rpm_filtered                    <- av_panel$rpm_filtered
aviti_records$aviti_rpkmf                           <- av_panel$rpkmf
aviti_records$aviti_taxon_bleed_ratio_pct           <- av_panel$taxon_bleed_ratio_pct
aviti_records$aviti_taxon_max_reads_in_sample_type  <- av_panel$taxon_max_reads_in_sample_type
aviti_records$aviti_taxon_max_rpm_in_sample_type    <- av_panel$taxon_max_rpm_in_sample_type
aviti_records$aviti_donor_taxon                     <- av_panel$donor_taxon
aviti_records$aviti_nc_genome_cov_pct_max           <- av_panel$nc_genome_cov_pct_max

ont_records <- build_records(
  "ONT", best_ont$lhs_id, "target_species_ont",
  on_panel, hit_on, cov_on
)
ont_records$ont_nb_reads                          <- on_panel$nb_reads
ont_records$ont_total_mapped_reads                <- on_panel$total_mapped_reads
ont_records$ont_rpm_mapped                        <- on_panel$rpm_mapped
ont_records$ont_genome_cov_pct                    <- on_panel$genome_cov_pct
ont_records$ont_depth_cov                         <- on_panel$depth_cov
ont_records$ont_assignment_lvl                    <- on_panel$assignment_lvl
ont_records$ont_assignment_rank                   <- on_panel$assignment_rank
ont_records$ont_taxon_bleed_ratio_pct             <- on_panel$taxon_bleed_ratio_pct
ont_records$ont_taxon_max_reads_in_sample_type    <- on_panel$taxon_max_reads_in_sample_type
ont_records$ont_taxon_max_rpm_in_sample_type      <- on_panel$taxon_max_rpm_in_sample_type
ont_records$ont_donor_taxon                       <- on_panel$donor_taxon
ont_records$ont_nc_genome_cov_pct_max             <- on_panel$nc_genome_cov_pct_max

aviti_records <- aviti_records[
  order(aviti_records$sample_type, aviti_records$sample_id, aviti_records$virus),
  ,
  drop = FALSE
]
ont_records <- ont_records[
  order(ont_records$sample_type, ont_records$sample_id, ont_records$virus),
  ,
  drop = FALSE
]
rownames(aviti_records) <- NULL
rownames(ont_records)   <- NULL

write.csv(aviti_records, file.path(tables_dir, "aviti_records.csv"), row.names = FALSE, na = "")
write.csv(ont_records,   file.path(tables_dir, "ont_records.csv"),   row.names = FALSE, na = "")
message(sprintf("Wrote aviti_records.csv (%d rows), ont_records.csv (%d rows)",
                nrow(aviti_records), nrow(ont_records)))
if (interactive()) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    make_plateau_plot <- function(feasible, best_row, method_label, params) {
      long <- do.call(rbind, lapply(params, function(p) data.frame(
        method = method_label,
        parameter = p,
        value = as_num(feasible[[p]]),
        stringsAsFactors = FALSE
      )))
      sel <- do.call(rbind, lapply(params, function(p) data.frame(
        method = method_label,
        parameter = p,
        value = as_num(best_row[[p]]),
        stringsAsFactors = FALSE
      )))
      long$parameter <- factor(long$parameter, levels = params)
      sel$parameter  <- factor(sel$parameter,  levels = params)
      ggplot(long, aes(x = value, y = parameter)) +
        geom_violin(fill = "#8DA0CB", alpha = 0.6, trim = FALSE, color = NA) +
        geom_jitter(height = 0.08, width = 0, alpha = 0.25, size = 1.0, color = "grey35") +
        geom_point(data = sel, aes(x = value, y = parameter), color = "#D62728", size = 2.8) +
        labs(
          title = sprintf("%s feasible plateau (n = %d)", method_label, nrow(feasible)),
          subtitle = "Violin + points show near-best region; red point is selected cutoff",
          x = "Cutoff value", y = NULL
        ) +
        theme_minimal(base_size = 12)
    }
    print(make_plateau_plot(feasible_aviti, best_aviti, "AVITI", names(aviti_ranges)))
    print(make_plateau_plot(feasible_ont,   best_ont,   "ONT",   names(ont_ranges)))
  } else {
    message("interactive(): ggplot2 not installed; skipping display-only plateau plots.")
  }
}

invisible(list(aviti = aviti_records, ont = ont_records, cutoffs = selected_cutoffs))
