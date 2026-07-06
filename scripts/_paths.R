script_path <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(x) {
    if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = TRUE) else NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]
  if (length(ofiles) > 0L) return(ofiles[length(ofiles)])

  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  stop("Could not determine script path.")
}

script_dir <- dirname(script_path())
scripts_dir <- script_dir
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir <- file.path(project_dir, "data_input")
tables_dir <- file.path(project_dir, "tables")
figures_dir <- file.path(project_dir, "figures")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
