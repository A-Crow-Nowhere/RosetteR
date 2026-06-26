#!/usr/bin/env Rscript

# =============================================================================
# 04_correct_and_summarize_rosettes.R
#
# Step 4 correction + summary script for rosette/cell analysis.
#
# Supports:
#   1. Per-replicate / normal mode:
#        --input points at a Step 3 output folder or a 04_corrected folder.
#
#   2. Project-level mode:
#        --input points at a project output folder organized like:
#        project/sample/replicate/image-or-step-folder/<cell table>
#
# Default behavior:
#   - Final summaries use accepted cells only.
#   - Raw combined table preserves rejected cells.
#   - Rosettes with fewer than --min_cells_per_rosette accepted cells are demoted
#     into single/unassigned counts.
#   - Writes both modern summary names and wrapper-compatible *.corrected.tsv names.
#
# Important robustness patch:
#   - Drops transient valid_rosette / valid_rosette.x / valid_rosette.y columns before
#     recomputing correction status. This prevents merge collisions when the project-
#     level summary reads already-corrected per-image/per-replicate tables.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

SCRIPT_VERSION <- "04_summary_v04_valid_rosette_collision_safe"

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input directory. Either Step 3 output, 04_corrected output, or project-level folder."),

  make_option(c("--out"), type = "character", default = NULL,
              help = "Output directory for Step 4 correction/summary tables."),

  make_option(c("--project_level"), type = "logical", default = FALSE,
              help = "If TRUE, interpret input paths using sample/replicate/image path depths. Default: FALSE"),

  make_option(c("--sample_depth"), type = "integer", default = 1,
              help = "In project-level mode, path part under --input used as sample_id. Default: 1"),

  make_option(c("--replicate_depth"), type = "integer", default = 2,
              help = "In project-level mode, path part under --input used as replicate_id. Default: 2"),

  make_option(c("--image_depth"), type = "integer", default = 4,
              help = paste(
                "In project-level mode, path part under --input used as image_id.",
                "Default: 4, suitable for OUT/sample/replicate/04_corrected/image.cells.corrected.tsv.",
                "Use 3 for OUT/sample/replicate/image/all_cells.tsv."
              )),

  make_option(c("--sample_id"), type = "character", default = "",
              help = "Optional wrapper-provided sample ID. Overrides inferred sample_id when provided."),

  make_option(c("--replicate_id"), type = "character", default = "",
              help = "Optional wrapper-provided replicate ID. Overrides inferred replicate_id when provided."),

  make_option(c("--cell_table_name"), type = "character", default = "",
              help = "Optional exact cell table filename to search for. If blank, common names/patterns are tried."),

  make_option(c("--recursive"), type = "logical", default = TRUE,
              help = "Search recursively under --input for cell tables. Default: TRUE"),

  make_option(c("--accepted_only"), type = "logical", default = TRUE,
              help = "Use only accepted cells for final summaries. Default: TRUE"),

  make_option(c("--single_cell_rosette_correction"), type = "logical", default = TRUE,
              help = "Demote rosettes with fewer than --min_cells_per_rosette accepted cells. Default: TRUE"),

  make_option(c("--min_cells_per_rosette"), type = "integer", default = 2,
              help = "Minimum accepted cells required to keep an object as a rosette. Default: 2"),

  make_option(c("--corrected_suffix_outputs"), type = "logical", default = TRUE,
              help = "Also write wrapper-compatible *.corrected.tsv summary aliases. Default: TRUE"),

  make_option(c("--write_per_image_corrected_cells"), type = "logical", default = TRUE,
              help = "Write one *.cells.corrected.tsv table per image. Default: TRUE"),

  make_option(c("--sample_col"), type = "character", default = "",
              help = "Optional sample column name for non-project mode."),

  make_option(c("--rosette_col"), type = "character", default = "",
              help = "Optional rosette/center assignment column name. If blank, inferred."),

  make_option(c("--cell_col"), type = "character", default = "",
              help = "Optional cell ID column name. If blank, inferred."),

  make_option(c("--area_col"), type = "character", default = "",
              help = "Optional cell area column name. If blank, inferred.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$out)) {
  stop("Required arguments: --input and --out", call. = FALSE)
}

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

message("[04_summary] SCRIPT_VERSION: ", SCRIPT_VERSION)
message("[04_summary] input: ", opt$input)
message("[04_summary] out:   ", opt$out)
message("[04_summary] project_level: ", opt$project_level)
message("[04_summary] accepted_only: ", opt$accepted_only)
message("[04_summary] single_cell_rosette_correction: ", opt$single_cell_rosette_correction)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

escape_regex <- function(x) {
  gsub("([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1", x)
}

as_clean_logical <- function(x) {
  if (is.logical(x)) return(x)

  if (is.numeric(x) || is.integer(x)) {
    return(!is.na(x) & x != 0)
  }

  x_chr <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(x_chr))

  out[x_chr %in% c("true", "t", "yes", "y", "1", "accepted", "pass", "passed")] <- TRUE
  out[x_chr %in% c("false", "f", "no", "n", "0", "rejected", "fail", "failed")] <- FALSE

  out
}

safe_uniqueN <- function(x) {
  data.table::uniqueN(x[!is.na(x)])
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_sum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(0)
  sum(x, na.rm = TRUE)
}

pick_first_existing_col <- function(dt, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% names(dt)][1]

  if (is.na(hit) || length(hit) == 0) {
    if (required) {
      stop(
        "Could not infer ", label, ". Tried: ",
        paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
    return(NULL)
  }

  hit
}

drop_transient_summary_cols <- function(dt) {
  dt <- as.data.table(dt)

  transient_cols <- c(
    "valid_rosette",
    "valid_rosette.x",
    "valid_rosette.y"
  )

  drop_cols <- intersect(names(dt), transient_cols)

  if (length(drop_cols) > 0) {
    dt[, (drop_cols) := NULL]
  }

  dt[]
}

relative_path_parts <- function(path, input_root) {
  abs_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  abs_root <- normalizePath(input_root, winslash = "/", mustWork = FALSE)

  rel <- sub(paste0("^", escape_regex(abs_root), "/?"), "", abs_path)
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]

  # Remove filename.
  if (length(parts) > 1) {
    parts <- parts[-length(parts)]
  } else {
    parts <- character()
  }

  parts
}

part_or_na <- function(parts, idx) {
  if (length(parts) >= idx && idx >= 1) return(parts[[idx]])
  NA_character_
}

infer_image_id_from_file <- function(source_file, parts, image_depth) {
  # If the requested image path component is a processing folder like 04_corrected,
  # use the filename stem instead. This is common for:
  #   sample/replicate/04_corrected/Snap_001.cells.corrected.tsv
  image_id <- part_or_na(parts, image_depth)

  if (is.na(image_id) || image_id == "" || image_id %in% c("03_cells", "04_corrected", "04_summary")) {
    stem <- basename(source_file)
    stem <- sub("\\.cells\\.corrected\\.tsv$", "", stem)
    stem <- sub("\\.corrected\\.tsv$", "", stem)
    stem <- sub("\\.all_cells\\.tsv$", "", stem)
    stem <- sub("\\.tsv$", "", stem)
    image_id <- stem
  }

  image_id
}

infer_ids_from_path <- function(source_file, input_root, project_level,
                                sample_depth, replicate_depth, image_depth) {
  parts <- relative_path_parts(source_file, input_root)

  if (isTRUE(project_level)) {
    sample_id <- part_or_na(parts, sample_depth)
    replicate_id <- part_or_na(parts, replicate_depth)
    image_id <- infer_image_id_from_file(source_file, parts, image_depth)

    if (is.na(sample_id) || sample_id == "") sample_id <- "unknown_sample"
    if (is.na(replicate_id) || replicate_id == "") replicate_id <- "unknown_replicate"
    if (is.na(image_id) || image_id == "") image_id <- tools::file_path_sans_ext(basename(source_file))

    return(list(sample_id = sample_id, replicate_id = replicate_id, image_id = image_id))
  }

  # Non-project mode: Step 3 often outputs image tables directly or under image folders.
  fallback <- if (length(parts) >= 1) parts[[1]] else tools::file_path_sans_ext(basename(source_file))

  list(sample_id = fallback, replicate_id = "replicate_01", image_id = fallback)
}

find_cell_tables <- function(input_dir, exact_name = "", recursive = TRUE) {
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir, call. = FALSE)
  }

  common_names <- c(
    "all_cells.tsv",
    "cells.tsv",
    "cell_table.tsv",
    "assigned_cells.tsv",
    "cell_assignments.tsv",
    "rosette_cell_assignments.tsv",
    "all_assigned_cells.tsv",
    "per_cell.tsv"
  )

  if (!is.null(exact_name) && nzchar(exact_name)) {
    files <- list.files(
      input_dir,
      pattern = paste0("^", escape_regex(exact_name), "$"),
      recursive = recursive,
      full.names = TRUE
    )
  } else {
    all_tsv <- list.files(input_dir, pattern = "\\.tsv$", recursive = recursive, full.names = TRUE)

    # Exclude summary tables, but keep per-image/per-cell corrected tables.
    all_tsv <- all_tsv[
      !grepl("summary", basename(all_tsv), ignore.case = TRUE) &
        !grepl("project_summary", basename(all_tsv), ignore.case = TRUE) &
        !grepl("global_summary", basename(all_tsv), ignore.case = TRUE) &
        !grepl("sample_summary", basename(all_tsv), ignore.case = TRUE) &
        !grepl("replicate_summary", basename(all_tsv), ignore.case = TRUE) &
        !grepl("image_summary", basename(all_tsv), ignore.case = TRUE)
    ]

    # Prefer explicit corrected cell tables if present. This helps project-level summary
    # read the per-image corrected cell files in 04_corrected folders.
    files <- all_tsv[grepl("\\.cells\\.corrected\\.tsv$", basename(all_tsv))]

    if (length(files) == 0) {
      files <- all_tsv[basename(all_tsv) %in% common_names]
    }

    if (length(files) == 0) {
      files <- all_tsv[
        grepl("cell", basename(all_tsv), ignore.case = TRUE) &
          !grepl("summary", basename(all_tsv), ignore.case = TRUE)
      ]
    }
  }

  unique(files)
}

standardize_cell_table <- function(dt, source_file, input_root, opt) {
  dt <- as.data.table(dt)
  dt <- drop_transient_summary_cols(dt)

  ids <- infer_ids_from_path(
    source_file = source_file,
    input_root = input_root,
    project_level = opt$project_level,
    sample_depth = opt$sample_depth,
    replicate_depth = opt$replicate_depth,
    image_depth = opt$image_depth
  )

  dt[, source_file := source_file]
  dt[, sample_id := ids$sample_id]
  dt[, replicate_id := ids$replicate_id]
  dt[, image_id := ids$image_id]

  # Wrapper compatibility: explicit sample/replicate IDs override inference.
  if (!is.null(opt$sample_id) && nzchar(opt$sample_id)) {
    dt[, sample_id := opt$sample_id]
  }

  if (!is.null(opt$replicate_id) && nzchar(opt$replicate_id)) {
    dt[, replicate_id := opt$replicate_id]
  }

  # In non-project mode, allow user/table sample column to override folder sample
  # only if sample_id was not explicitly supplied by wrapper.
  if (!isTRUE(opt$project_level) && (!nzchar(opt$sample_id))) {
    sample_col <- NULL

    if (nzchar(opt$sample_col)) {
      if (!opt$sample_col %in% names(dt)) {
        stop("Requested --sample_col not found: ", opt$sample_col, call. = FALSE)
      }
      sample_col <- opt$sample_col
    } else {
      sample_col <- pick_first_existing_col(
        dt,
        c("sample_id", "sample", "sample_name", "image_id", "image", "folder", "snap", "file_id"),
        required = FALSE,
        label = "sample column"
      )
    }

    if (!is.null(sample_col)) {
      dt[, sample_id := as.character(get(sample_col))]
      dt[, image_id := as.character(get(sample_col))]
    }
  }

  # Cell ID.
  cell_col <- NULL
  if (nzchar(opt$cell_col)) {
    if (!opt$cell_col %in% names(dt)) stop("Requested --cell_col not found: ", opt$cell_col, call. = FALSE)
    cell_col <- opt$cell_col
  } else {
    cell_col <- pick_first_existing_col(
      dt,
      c("cell_id", "cell", "label", "cell_label", "object_id", "segment_id", "cell_index"),
      required = FALSE,
      label = "cell ID column"
    )
  }

  if (is.null(cell_col)) {
    dt[, cell_id := as.character(seq_len(.N))]
  } else if (cell_col != "cell_id") {
    dt[, cell_id := as.character(get(cell_col))]
  } else {
    dt[, cell_id := as.character(cell_id)]
  }

  dt[, cell_uid := paste(source_file, cell_id, sep = "__")]

  # Rosette / assignment ID.
  rosette_col <- NULL
  if (nzchar(opt$rosette_col)) {
    if (!opt$rosette_col %in% names(dt)) stop("Requested --rosette_col not found: ", opt$rosette_col, call. = FALSE)
    rosette_col <- opt$rosette_col
  } else {
    rosette_col <- pick_first_existing_col(
      dt,
      c(
        "rosette_id",
        "assigned_rosette_id",
        "center_id",
        "assigned_center_id",
        "candidate_id",
        "merged_center_id",
        "cluster_id",
        "blob_id",
        "object_id"
      ),
      required = FALSE,
      label = "rosette/assignment column"
    )
  }

  if (is.null(rosette_col)) {
    dt[, rosette_id := NA_character_]
  } else if (rosette_col != "rosette_id") {
    dt[, rosette_id := as.character(get(rosette_col))]
  } else {
    dt[, rosette_id := as.character(rosette_id)]
  }

  dt[
    is.na(rosette_id) |
      rosette_id == "" |
      tolower(rosette_id) %in% c("na", "nan", "none", "null", "unassigned"),
    rosette_id := NA_character_
  ]

  dt[, rosette_uid := ifelse(is.na(rosette_id), NA_character_, paste(source_file, rosette_id, sep = "__"))]

  # Area.
  area_col <- NULL
  if (nzchar(opt$area_col)) {
    if (!opt$area_col %in% names(dt)) stop("Requested --area_col not found: ", opt$area_col, call. = FALSE)
    area_col <- opt$area_col
  } else {
    area_col <- pick_first_existing_col(
      dt,
      c("area_px", "cell_area_px", "area", "cell_area", "n_px", "pixels"),
      required = FALSE,
      label = "area column"
    )
  }

  if (is.null(area_col)) {
    dt[, area_px := NA_real_]
  } else if (area_col != "area_px") {
    dt[, area_px := suppressWarnings(as.numeric(get(area_col)))]
  } else {
    dt[, area_px := suppressWarnings(as.numeric(area_px))]
  }

  # Accepted-cell status.
  accepted_col <- pick_first_existing_col(
    dt,
    c("accepted_cell", "accepted_cells", "cell_accepted", "accepted", "pass_cell", "cell_pass", "accepted_cell_for_summary"),
    required = FALSE,
    label = "accepted-cell column"
  )

  if (is.null(accepted_col)) {
    dt[, accepted_cell_for_summary := TRUE]
    dt[, accepted_cell_source_col := NA_character_]
  } else {
    dt[, accepted_cell_for_summary := as_clean_logical(get(accepted_col))]
    dt[is.na(accepted_cell_for_summary), accepted_cell_for_summary := FALSE]
    dt[, accepted_cell_source_col := accepted_col]
  }

  dt[]
}

empty_dt <- function() data.table()

write_empty_outputs <- function(out_dir) {
  fwrite(empty_dt(), file.path(out_dir, "all_cells_combined.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "all_cells_summary_input.accepted_only.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "cells_used_for_rosette_summary.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "cells_counted_as_single_or_unassigned.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "rosette_summary.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "image_summary.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "replicate_summary.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "sample_summary.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "global_summary.tsv"), sep = "\t")

  fwrite(empty_dt(), file.path(out_dir, "image_summary.corrected.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "replicate_summary.corrected.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "sample_summary.corrected.tsv"), sep = "\t")
  fwrite(empty_dt(), file.path(out_dir, "global_summary.corrected.tsv"), sep = "\t")
}

summarize_level <- function(level_name,
                            grouping_cols,
                            cells_all,
                            cells_summary_input,
                            rosette_summary,
                            cells_single_summary_input) {
  all_groups <- unique(cells_all[, ..grouping_cols])

  raw_counts <- cells_all[
    ,
    .(
      raw_cell_rows_all = as.integer(.N),
      raw_unique_cells_all = as.integer(safe_uniqueN(cell_uid)),
      raw_accepted_cell_rows = as.integer(sum(accepted_cell_for_summary == TRUE, na.rm = TRUE)),
      raw_rejected_cell_rows = as.integer(sum(accepted_cell_for_summary == FALSE | is.na(accepted_cell_for_summary), na.rm = TRUE)),
      n_cells_total = as.integer(safe_uniqueN(cell_uid))
    ),
    by = grouping_cols
  ]

  accepted_counts <- cells_summary_input[
    ,
    .(
      accepted_cell_rows_used_for_summary = as.integer(.N),
      accepted_unique_cells_used_for_summary = as.integer(safe_uniqueN(cell_uid)),
      accepted_mean_cell_area_px = safe_mean(area_px),
      accepted_median_cell_area_px = safe_median(area_px)
    ),
    by = grouping_cols
  ]

  if (nrow(rosette_summary) > 0) {
    rosette_counts <- rosette_summary[
      ,
      .(
        n_rosettes = as.integer(.N),
        n_rosettes_corrected = as.integer(.N),
        total_cells_in_rosettes = as.integer(sum(n_cells, na.rm = TRUE)),
        n_cells_in_rosettes_corrected = as.integer(sum(n_cells, na.rm = TRUE)),
        mean_cells_per_rosette = safe_mean(n_cells),
        median_cells_per_rosette = safe_median(n_cells),
        mean_cells_per_rosette_corrected = safe_mean(n_cells)
      ),
      by = grouping_cols
    ]
  } else {
    rosette_counts <- data.table()
  }

  if (nrow(cells_single_summary_input) > 0) {
    single_counts <- cells_single_summary_input[
      ,
      .(
        n_single_or_unassigned_cells = as.integer(safe_uniqueN(cell_uid)),
        n_single_or_unassigned_rows = as.integer(.N),
        n_single_cells_corrected = as.integer(safe_uniqueN(cell_uid))
      ),
      by = grouping_cols
    ]
  } else {
    single_counts <- data.table()
  }

  # Raw rosette count before single-cell correction.
  raw_rosette_counts <- cells_summary_input[!is.na(rosette_id)][
    ,
    .(n_rosettes_raw = as.integer(safe_uniqueN(rosette_uid))),
    by = grouping_cols
  ]

  # Demoted rosette count by grouping level.
  demoted_counts <- cells_single_summary_input[!is.na(rosette_id)][
    ,
    .(n_demoted_single_cell_rosettes = as.integer(safe_uniqueN(rosette_uid))),
    by = grouping_cols
  ]

  out <- Reduce(
    function(x, y) merge(x, y, by = grouping_cols, all.x = TRUE),
    list(all_groups, raw_counts, accepted_counts, raw_rosette_counts, rosette_counts, single_counts, demoted_counts)
  )

  count_cols <- c(
    "raw_cell_rows_all",
    "raw_unique_cells_all",
    "raw_accepted_cell_rows",
    "raw_rejected_cell_rows",
    "n_cells_total",
    "accepted_cell_rows_used_for_summary",
    "accepted_unique_cells_used_for_summary",
    "n_rosettes_raw",
    "n_rosettes",
    "n_rosettes_corrected",
    "total_cells_in_rosettes",
    "n_cells_in_rosettes_corrected",
    "n_single_or_unassigned_cells",
    "n_single_or_unassigned_rows",
    "n_single_cells_corrected",
    "n_demoted_single_cell_rosettes"
  )

  for (cc in count_cols) {
    if (cc %in% names(out)) {
      set(out, which(is.na(out[[cc]])), cc, 0L)
      out[, (cc) := as.integer(get(cc))]
    }
  }

  # Ensure expected columns exist even if absent from all joins.
  for (cc in count_cols) {
    if (!cc %in% names(out)) out[, (cc) := 0L]
  }

  out[, total_accepted_cells_after_correction := as.integer(n_cells_in_rosettes_corrected + n_single_cells_corrected)]
  out[, summary_level := level_name]

  setcolorder(out, c("summary_level", grouping_cols, setdiff(names(out), c("summary_level", grouping_cols))))
  setorderv(out, grouping_cols)

  out[]
}

write_per_image_corrected_cells <- function(cells_all, out_dir) {
  if (nrow(cells_all) == 0) return(invisible(NULL))

  # Write one file per image_id. This is useful for project-level re-summary and debugging.
  image_keys <- unique(cells_all[, .(sample_id, replicate_id, image_id)])

  for (ii in seq_len(nrow(image_keys))) {
    sid <- image_keys$sample_id[[ii]]
    rid <- image_keys$replicate_id[[ii]]
    iid <- image_keys$image_id[[ii]]

    sub <- cells_all[sample_id == sid & replicate_id == rid & image_id == iid]
    if (nrow(sub) == 0) next

    safe_iid <- gsub("[^A-Za-z0-9_.-]", "_", iid)
    out_file <- file.path(out_dir, paste0(safe_iid, ".cells.corrected.tsv"))
    fwrite(sub, out_file, sep = "\t")
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Discover and read input tables
# -----------------------------------------------------------------------------

cell_files <- find_cell_tables(
  input_dir = opt$input,
  exact_name = opt$cell_table_name,
  recursive = opt$recursive
)

message("[04_summary] Found ", length(cell_files), " candidate cell table(s).")

if (length(cell_files) == 0) {
  warning("[04_summary] No cell tables found. Writing empty outputs.")
  write_empty_outputs(opt$out)
  quit(status = 0)
}

cell_tables <- vector("list", length(cell_files))

for (i in seq_along(cell_files)) {
  f <- cell_files[[i]]
  message("[04_summary] Reading: ", f)

  dt <- tryCatch(
    fread(f),
    error = function(e) {
      warning("[04_summary] Failed to read ", f, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(dt) || nrow(dt) == 0) {
    warning("[04_summary] Skipping empty/unreadable table: ", f)
    next
  }

  cell_tables[[i]] <- standardize_cell_table(
    dt = dt,
    source_file = f,
    input_root = opt$input,
    opt = opt
  )
}

cell_tables <- Filter(Negate(is.null), cell_tables)

if (length(cell_tables) == 0) {
  warning("[04_summary] No readable non-empty cell tables found. Writing empty outputs.")
  write_empty_outputs(opt$out)
  quit(status = 0)
}

cells_all <- rbindlist(cell_tables, fill = TRUE)
cells_all <- drop_transient_summary_cols(cells_all)

message("[04_summary] Combined raw cell rows: ", nrow(cells_all))

fwrite(cells_all, file.path(opt$out, "all_cells_combined.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Accepted-cell filter for summary input
# -----------------------------------------------------------------------------

if (isTRUE(opt$accepted_only)) {
  before_n <- nrow(cells_all)

  cells_summary_input <- cells_all[
    !is.na(accepted_cell_for_summary) & accepted_cell_for_summary == TRUE
  ]

  after_n <- nrow(cells_summary_input)

  message(
    "[04_summary] Using accepted cells only for final summary: ",
    after_n, " / ", before_n, " cells retained."
  )
} else {
  cells_summary_input <- copy(cells_all)

  message(
    "[04_summary] accepted_only is FALSE; using all cells for final summary: ",
    nrow(cells_summary_input), " cells."
  )
}

# Critical collision-safety patch.
cells_summary_input <- drop_transient_summary_cols(cells_summary_input)

fwrite(cells_summary_input, file.path(opt$out, "all_cells_summary_input.accepted_only.tsv"), sep = "\t")

if (nrow(cells_summary_input) == 0) {
  warning("[04_summary] No cells remain after accepted-cell filtering. Writing empty summary outputs.")
  write_empty_outputs(opt$out)
  quit(status = 0)
}

# -----------------------------------------------------------------------------
# Single-cell rosette correction
# -----------------------------------------------------------------------------

assigned_cells <- cells_summary_input[!is.na(rosette_id)]
unassigned_cells <- cells_summary_input[is.na(rosette_id)]
assigned_cells <- drop_transient_summary_cols(assigned_cells)
unassigned_cells <- drop_transient_summary_cols(unassigned_cells)

rosette_group_cols <- c("sample_id", "replicate_id", "image_id", "source_file", "rosette_id", "rosette_uid")

if (nrow(assigned_cells) > 0) {
  rosette_precheck <- assigned_cells[
    ,
    .(
      n_accepted_cells_pre_correction = as.integer(safe_uniqueN(cell_uid)),
      n_rows_pre_correction = as.integer(.N)
    ),
    by = rosette_group_cols
  ]
} else {
  rosette_precheck <- data.table()
}

if (isTRUE(opt$single_cell_rosette_correction)) {
  if (nrow(rosette_precheck) > 0) {
    valid_rosettes <- rosette_precheck[
      n_accepted_cells_pre_correction >= opt$min_cells_per_rosette,
      ..rosette_group_cols
    ]

    if (nrow(valid_rosettes) > 0) {
      assigned_cells <- drop_transient_summary_cols(assigned_cells)

      valid_rosettes <- unique(valid_rosettes)
      valid_rosettes[, valid_rosette := TRUE]

      assigned_cells <- merge(
        assigned_cells,
        valid_rosettes,
        by = rosette_group_cols,
        all.x = TRUE,
        sort = FALSE
      )

      assigned_cells[is.na(valid_rosette), valid_rosette := FALSE]
    } else {
      assigned_cells <- drop_transient_summary_cols(assigned_cells)
      assigned_cells[, valid_rosette := FALSE]
    }
  } else {
    assigned_cells <- drop_transient_summary_cols(assigned_cells)
    assigned_cells[, valid_rosette := logical()]
  }

  cells_rosette_summary_input <- assigned_cells[valid_rosette == TRUE]

  cells_single_summary_input <- rbindlist(
    list(unassigned_cells, assigned_cells[valid_rosette == FALSE]),
    fill = TRUE
  )

  n_demoted <- nrow(assigned_cells[valid_rosette == FALSE])

  message(
    "[04_summary] Single-cell correction enabled. Demoted ",
    n_demoted,
    " accepted cell row(s) from undersized rosettes to single-cell/unassigned counts."
  )
} else {
  cells_rosette_summary_input <- assigned_cells
  cells_single_summary_input <- unassigned_cells

  if (nrow(cells_rosette_summary_input) > 0) {
    cells_rosette_summary_input[, valid_rosette := TRUE]
  }

  message("[04_summary] Single-cell correction disabled.")
}

fwrite(cells_rosette_summary_input, file.path(opt$out, "cells_used_for_rosette_summary.tsv"), sep = "\t")
fwrite(cells_single_summary_input, file.path(opt$out, "cells_counted_as_single_or_unassigned.tsv"), sep = "\t")

# Preserve per-image corrected cell tables after recomputed correction status.
if (isTRUE(opt$write_per_image_corrected_cells)) {
  corrected_all <- rbindlist(list(cells_rosette_summary_input, cells_single_summary_input), fill = TRUE)
  write_per_image_corrected_cells(corrected_all, opt$out)
}

# -----------------------------------------------------------------------------
# Rosette-level summary
# -----------------------------------------------------------------------------

if (nrow(cells_rosette_summary_input) > 0) {
  rosette_summary <- cells_rosette_summary_input[
    ,
    .(
      n_cells = as.integer(safe_uniqueN(cell_uid)),
      n_cell_rows = as.integer(.N),
      mean_cell_area_px = safe_mean(area_px),
      median_cell_area_px = safe_median(area_px),
      total_cell_area_px = safe_sum(area_px)
    ),
    by = rosette_group_cols
  ]

  setorder(rosette_summary, sample_id, replicate_id, image_id, rosette_id)
} else {
  rosette_summary <- data.table(
    sample_id = character(),
    replicate_id = character(),
    image_id = character(),
    source_file = character(),
    rosette_id = character(),
    rosette_uid = character(),
    n_cells = integer(),
    n_cell_rows = integer(),
    mean_cell_area_px = numeric(),
    median_cell_area_px = numeric(),
    total_cell_area_px = numeric()
  )
}

fwrite(rosette_summary, file.path(opt$out, "rosette_summary.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# Image-, replicate-, sample-level summaries
# -----------------------------------------------------------------------------

image_summary <- summarize_level(
  level_name = "image",
  grouping_cols = c("sample_id", "replicate_id", "image_id"),
  cells_all = cells_all,
  cells_summary_input = cells_summary_input,
  rosette_summary = rosette_summary,
  cells_single_summary_input = cells_single_summary_input
)

replicate_summary <- summarize_level(
  level_name = "replicate",
  grouping_cols = c("sample_id", "replicate_id"),
  cells_all = cells_all,
  cells_summary_input = cells_summary_input,
  rosette_summary = rosette_summary,
  cells_single_summary_input = cells_single_summary_input
)

sample_summary <- summarize_level(
  level_name = "sample",
  grouping_cols = c("sample_id"),
  cells_all = cells_all,
  cells_summary_input = cells_summary_input,
  rosette_summary = rosette_summary,
  cells_single_summary_input = cells_single_summary_input
)

fwrite(image_summary, file.path(opt$out, "image_summary.tsv"), sep = "\t")
fwrite(replicate_summary, file.path(opt$out, "replicate_summary.tsv"), sep = "\t")
fwrite(sample_summary, file.path(opt$out, "sample_summary.tsv"), sep = "\t")

if (isTRUE(opt$corrected_suffix_outputs)) {
  fwrite(image_summary, file.path(opt$out, "image_summary.corrected.tsv"), sep = "\t")
  fwrite(replicate_summary, file.path(opt$out, "replicate_summary.corrected.tsv"), sep = "\t")
  fwrite(sample_summary, file.path(opt$out, "sample_summary.corrected.tsv"), sep = "\t")
}

# -----------------------------------------------------------------------------
# Global / project summary
# -----------------------------------------------------------------------------

global_summary <- data.table(
  script_version = SCRIPT_VERSION,
  input = opt$input,
  project_level = opt$project_level,
  accepted_only = opt$accepted_only,
  single_cell_rosette_correction = opt$single_cell_rosette_correction,
  min_cells_per_rosette = opt$min_cells_per_rosette,

  n_samples = as.integer(safe_uniqueN(cells_all$sample_id)),
  n_replicates = as.integer(safe_uniqueN(paste(cells_all$sample_id, cells_all$replicate_id, sep = "__"))),
  n_images = as.integer(safe_uniqueN(paste(cells_all$sample_id, cells_all$replicate_id, cells_all$image_id, sep = "__"))),
  n_input_cell_tables = as.integer(length(cell_files)),

  raw_cell_rows_all = as.integer(nrow(cells_all)),
  raw_unique_cells_all = as.integer(safe_uniqueN(cells_all$cell_uid)),

  raw_accepted_cell_rows = as.integer(sum(cells_all$accepted_cell_for_summary == TRUE, na.rm = TRUE)),
  raw_rejected_cell_rows = as.integer(sum(cells_all$accepted_cell_for_summary == FALSE | is.na(cells_all$accepted_cell_for_summary), na.rm = TRUE)),

  accepted_cell_rows_used_for_summary = as.integer(nrow(cells_summary_input)),
  accepted_unique_cells_used_for_summary = as.integer(safe_uniqueN(cells_summary_input$cell_uid)),

  n_rosettes_raw = as.integer(safe_uniqueN(cells_summary_input[!is.na(rosette_id)]$rosette_uid)),
  n_rosettes = as.integer(nrow(rosette_summary)),
  n_rosettes_corrected = as.integer(nrow(rosette_summary)),
  total_cells_in_rosettes = as.integer(sum(rosette_summary$n_cells, na.rm = TRUE)),
  n_cells_in_rosettes_corrected = as.integer(sum(rosette_summary$n_cells, na.rm = TRUE)),

  n_single_or_unassigned_cells = as.integer(safe_uniqueN(cells_single_summary_input$cell_uid)),
  n_single_or_unassigned_rows = as.integer(nrow(cells_single_summary_input)),
  n_single_cells_corrected = as.integer(safe_uniqueN(cells_single_summary_input$cell_uid)),
  n_demoted_single_cell_rosettes = as.integer(safe_uniqueN(cells_single_summary_input[!is.na(rosette_id)]$rosette_uid)),

  n_cells_total = as.integer(safe_uniqueN(cells_all$cell_uid)),
  total_accepted_cells_after_correction = as.integer(
    sum(rosette_summary$n_cells, na.rm = TRUE) + safe_uniqueN(cells_single_summary_input$cell_uid)
  ),

  mean_rosettes_per_image_corrected = safe_mean(image_summary$n_rosettes_corrected),
  median_rosettes_per_image_corrected = safe_median(image_summary$n_rosettes_corrected),
  mean_cells_per_image = safe_mean(image_summary$n_cells_total),
  median_cells_per_image = safe_median(image_summary$n_cells_total),
  mean_cells_per_rosette_corrected = safe_mean(rosette_summary$n_cells)
)

fwrite(global_summary, file.path(opt$out, "global_summary.tsv"), sep = "\t")
fwrite(global_summary, file.path(opt$out, "project_summary.tsv"), sep = "\t")

if (isTRUE(opt$corrected_suffix_outputs)) {
  fwrite(global_summary, file.path(opt$out, "global_summary.corrected.tsv"), sep = "\t")
  fwrite(global_summary, file.path(opt$out, "project_summary.corrected.tsv"), sep = "\t")
}

# -----------------------------------------------------------------------------
# Text report
# -----------------------------------------------------------------------------

report_file <- file.path(opt$out, "summary_report.txt")

report_lines <- c(
  paste0("SCRIPT_VERSION: ", SCRIPT_VERSION),
  paste0("input: ", opt$input),
  paste0("out: ", opt$out),
  paste0("project_level: ", opt$project_level),
  paste0("accepted_only: ", opt$accepted_only),
  paste0("single_cell_rosette_correction: ", opt$single_cell_rosette_correction),
  paste0("min_cells_per_rosette: ", opt$min_cells_per_rosette),
  "",
  paste0("n_samples: ", global_summary$n_samples),
  paste0("n_replicates: ", global_summary$n_replicates),
  paste0("n_images: ", global_summary$n_images),
  paste0("n_input_cell_tables: ", global_summary$n_input_cell_tables),
  "",
  paste0("raw_cell_rows_all: ", global_summary$raw_cell_rows_all),
  paste0("raw_accepted_cell_rows: ", global_summary$raw_accepted_cell_rows),
  paste0("raw_rejected_cell_rows: ", global_summary$raw_rejected_cell_rows),
  paste0("accepted_cell_rows_used_for_summary: ", global_summary$accepted_cell_rows_used_for_summary),
  "",
  paste0("n_rosettes_raw: ", global_summary$n_rosettes_raw),
  paste0("n_rosettes_corrected: ", global_summary$n_rosettes_corrected),
  paste0("n_demoted_single_cell_rosettes: ", global_summary$n_demoted_single_cell_rosettes),
  paste0("n_cells_in_rosettes_corrected: ", global_summary$n_cells_in_rosettes_corrected),
  paste0("n_single_cells_corrected: ", global_summary$n_single_cells_corrected),
  paste0("total_accepted_cells_after_correction: ", global_summary$total_accepted_cells_after_correction)
)

writeLines(report_lines, report_file)

message("[04_summary] Wrote:")
message("[04_summary]   ", file.path(opt$out, "all_cells_combined.tsv"))
message("[04_summary]   ", file.path(opt$out, "all_cells_summary_input.accepted_only.tsv"))
message("[04_summary]   ", file.path(opt$out, "cells_used_for_rosette_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "cells_counted_as_single_or_unassigned.tsv"))
message("[04_summary]   ", file.path(opt$out, "rosette_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "image_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "replicate_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "sample_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "global_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "project_summary.tsv"))
message("[04_summary]   ", file.path(opt$out, "image_summary.corrected.tsv"))
message("[04_summary]   ", file.path(opt$out, "replicate_summary.corrected.tsv"))
message("[04_summary]   ", file.path(opt$out, "sample_summary.corrected.tsv"))
message("[04_summary]   ", file.path(opt$out, "global_summary.corrected.tsv"))
message("[04_summary]   ", report_file)
message("[04_summary] Done.")

