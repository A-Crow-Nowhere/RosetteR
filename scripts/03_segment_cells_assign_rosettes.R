#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(EBImage)
  library(data.table)
})

# ============================================================
# 03_segment_cells_assign_rosettes.R
#
# Purpose:
#   Segment cells/nuclei inside object masks from Step 1/2,
#   then uniquely assign detected cells to nearby rosette centers/radii.
#
# Expected project layout:
#
# RosetteCounter/
#   raw_images/
#   results/
#     outline_v01/ or outlines_v01/
#       <sample_name>/
#         mask.binary.png
#         gray.png
#     centers_v01/
#       <sample_name>/
#         candidate_centers.tsv
#     cells_debug/
#
# Step 2 compatibility:
#   Step 2 writes:
#     centers_v01/<sample>/candidate_centers.tsv
#
#   with columns:
#     rosette_candidate_id
#     parent_blob_id
#     center_x
#     center_y
#     fitted_radius_px
#     accepted
#     confidence
#
#   In Step 2:
#     center_x = EBImage matrix row
#     center_y = EBImage matrix column
#
#   In Step 3 internal coordinates:
#     x = image column
#     y = image row
#
#   Therefore this script maps:
#     x         <- center_y
#     y         <- center_x
#     radius    <- fitted_radius_px
#     object_id <- parent_blob_id
#
# Cell detection polarity:
#   --cell_signal_polarity bright
#     detects bright cell/nucleus signal.
#
#   --cell_signal_polarity dark
#     detects dark cell bodies/nuclei against brighter background.
#
# Speed/debug features:
#   - pairs outline and center folders by sample subfolder name
#   - crop-first segmentation
#   - max image/object/rosette caps
#   - skips irrelevant objects before expensive work
#   - writes progress continuously
#   - avoids global all-cell x all-rosette distance matrix
# ============================================================


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Project root, results directory, outline_v01/outlines_v01 directory, or one sample outline directory."),
  
  make_option(c("--out"), type = "character", default = "results/cells_v01",
              help = "Output directory."),
  
  make_option(c("--outlines_dir"), type = "character", default = NULL,
              help = "Directory containing outline sample folders. Example: results/outline_v01 or results/outlines_v01"),
  make_option(c("--centers_dir"), type = "character", default = NULL,
              help = "Directory containing center sample folders. Example: results/centers_v01"),
  make_option(c("--raw_dir"), type = "character", default = NULL,
              help = "Optional raw image directory. Example: raw_images"),
  
  make_option(c("--sample"), type = "character", default = NULL,
              help = "Optional sample name to process, e.g. Snap_246."),
  
  make_option(c("--object_mask_file"), type = "character", default = "mask.binary.png",
              help = "Exact object-mask filename inside each outline_v01/<sample> folder."),
  make_option(c("--nuclei_file"), type = "character", default = "gray.png",
              help = "Exact gray/nuclei/cell image filename inside each outline_v01/<sample> folder."),
  make_option(c("--rosette_file"), type = "character", default = "candidate_centers.tsv",
              help = "Exact rosette/center TSV filename inside each centers_v01/<sample> folder."),
  
  make_option(c("--accepted_only"), type = "logical", default = TRUE,
              help = "Use only accepted centers from candidate_centers.tsv when an accepted column exists."),
  
  make_option(c("--mask_channel"), type = "integer", default = 1,
              help = "Image channel to use for object mask if multichannel."),
  make_option(c("--nuclei_channel"), type = "integer", default = 1,
              help = "Image channel to use for nuclei/cell image if multichannel."),
  
  make_option(c("--max_images"), type = "integer", default = NA,
              help = "Debug: maximum number of image/sample folders to process."),
  make_option(c("--max_objects"), type = "integer", default = NA,
              help = "Debug: maximum number of objects/clumps per image."),
  make_option(c("--max_rosettes"), type = "integer", default = NA,
              help = "Debug: maximum number of rosettes per image."),
  
  make_option(c("--min_object_area_px"), type = "integer", default = 500,
              help = "Skip object masks smaller than this area."),
  make_option(c("--max_object_area_px"), type = "integer", default = NA,
              help = "Skip object masks larger than this area. NA disables."),
  
  make_option(c("--min_cell_area_px"), type = "integer", default = 10,
              help = "Minimum area for a segmented cell/nucleus object."),
  make_option(c("--max_cell_area_px"), type = "integer", default = 5000,
              help = "Maximum area for a segmented cell/nucleus object."),
  
  make_option(c("--threshold_quantile"), type = "double", default = 0.80,
              help = "Quantile threshold used inside object crops."),
  make_option(c("--cell_signal_polarity"), type = "character", default = "bright",
              help = "Whether cells/nuclei are bright or dark relative to local background. Use bright or dark."),
  make_option(c("--smooth_sigma"), type = "double", default = 1.0,
              help = "Gaussian blur sigma before thresholding. Use 0 to disable."),
  make_option(c("--opening_radius"), type = "integer", default = 1,
              help = "Morphological opening brush radius. Use 0 to disable."),
  
  make_option(c("--radius_multiplier"), type = "double", default = 1.15,
              help = "Cells are candidate members if within radius * this multiplier plus buffer."),
  make_option(c("--radius_buffer_px"), type = "double", default = 5,
              help = "Extra pixels added to rosette radius during assignment."),
  
  make_option(c("--write_debug_overlays"), type = "logical", default = FALSE,
              help = "Write debug overlay PNGs."),
  make_option(c("--debug_draw_n"), type = "integer", default = 10,
              help = "Only draw overlays for first N processed objects."),
  make_option(c("--draw_cell_ids"), type = "logical", default = FALSE,
              help = "Draw cell IDs on debug overlays. Slower."),
  
  make_option(c("--dry_run"), type = "logical", default = FALSE,
              help = "Only discover input files and print paired jobs; do not process images."),
  
  make_option(c("--verbose"), type = "logical", default = TRUE,
              help = "Print progress messages.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) && is.null(opt$outlines_dir)) {
  stop("Provide either --input or --outlines_dir.", call. = FALSE)
}

opt$cell_signal_polarity <- tolower(opt$cell_signal_polarity)

if (!opt$cell_signal_polarity %in% c("bright", "dark")) {
  stop("--cell_signal_polarity must be either 'bright' or 'dark'.", call. = FALSE)
}

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

msg <- function(...) {
  if (isTRUE(opt$verbose)) {
    message(...)
  }
}

tic <- function(label) {
  msg("[TIMER start] ", label)
  list(label = label, time = Sys.time())
}

toc <- function(t) {
  elapsed <- round(as.numeric(difftime(Sys.time(), t$time, units = "secs")), 3)
  msg("[TIMER done] ", t$label, " | ", elapsed, " sec")
  invisible(elapsed)
}

first_existing_dir <- function(paths) {
  paths <- paths[!is.na(paths)]
  hits <- paths[dir.exists(paths)]
  if (length(hits) == 0) {
    return(NULL)
  }
  hits[1]
}

truthy <- function(x) {
  if (is.logical(x)) {
    return(x %in% TRUE)
  }
  
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}


# ------------------------------------------------------------
# Filename patterns
# ------------------------------------------------------------

exclude_image_patterns <- c(
  "overlay",
  "debug",
  "qc",
  "preview",
  "annotated",
  "labels_on",
  "ids",
  "summary"
)

object_mask_patterns <- c(
  "^mask\\.binary\\.png$",
  "^object_mask\\.(tif|tiff|png)$",
  "^objects_mask\\.(tif|tiff|png)$",
  "^object_labels\\.(tif|tiff|png)$",
  "^objects_labels\\.(tif|tiff|png)$",
  "^labels\\.(tif|tiff|png)$",
  "^labeled_objects\\.(tif|tiff|png)$",
  "^clump_mask\\.(tif|tiff|png)$",
  "^clump_labels\\.(tif|tiff|png)$",
  "^outline_mask\\.(tif|tiff|png)$",
  "^mask_objects\\.(tif|tiff|png)$",
  "^binary_mask\\.(tif|tiff|png)$",
  "^segmentation_mask\\.(tif|tiff|png)$",
  
  "object.*mask.*\\.(tif|tiff|png)$",
  "objects.*mask.*\\.(tif|tiff|png)$",
  "object.*label.*\\.(tif|tiff|png)$",
  "objects.*label.*\\.(tif|tiff|png)$",
  "clump.*mask.*\\.(tif|tiff|png)$",
  "clump.*label.*\\.(tif|tiff|png)$",
  "outline.*mask.*\\.(tif|tiff|png)$",
  "mask.*object.*\\.(tif|tiff|png)$",
  "mask.*clump.*\\.(tif|tiff|png)$",
  "label.*object.*\\.(tif|tiff|png)$",
  "labeled.*object.*\\.(tif|tiff|png)$",
  "binary.*mask.*\\.(tif|tiff|png)$",
  "segmented.*mask.*\\.(tif|tiff|png)$",
  
  "objects.*\\.(tif|tiff|png)$",
  "object.*\\.(tif|tiff|png)$",
  "clumps.*\\.(tif|tiff|png)$",
  "clump.*\\.(tif|tiff|png)$",
  "mask.*\\.(tif|tiff|png)$"
)

nuclei_patterns <- c(
  "^gray\\.png$",
  "^gray\\.(tif|tiff|png|jpg|jpeg)$",
  "^nuclei\\.(tif|tiff|png|jpg|jpeg)$",
  "^nuclear\\.(tif|tiff|png|jpg|jpeg)$",
  "^cells\\.(tif|tiff|png|jpg|jpeg)$",
  "^cell_channel\\.(tif|tiff|png|jpg|jpeg)$",
  "^channel_B\\.(tif|tiff|png|jpg|jpeg)$",
  "^B\\.(tif|tiff|png|jpg|jpeg)$",
  
  "gray.*\\.(tif|tiff|png|jpg|jpeg)$",
  "nuc.*\\.(tif|tiff|png|jpg|jpeg)$",
  "nuclear.*\\.(tif|tiff|png|jpg|jpeg)$",
  "channel.*B.*\\.(tif|tiff|png|jpg|jpeg)$",
  "cell.*\\.(tif|tiff|png|jpg|jpeg)$"
)

rosette_patterns <- c(
  "^candidate_centers\\.tsv$",
  "^accepted_centers\\.tsv$",
  "^rosettes\\.tsv$",
  "^rosette_centers\\.tsv$",
  "^centers\\.tsv$",
  "^center_points\\.tsv$",
  "^step2_rosettes\\.tsv$",
  "^detected_rosettes\\.tsv$",
  
  "candidate_centers\\.tsv$",
  "accepted.*center.*\\.tsv$",
  "rosette.*\\.tsv$",
  "center.*\\.tsv$"
)

find_first_existing <- function(dir,
                                patterns,
                                exact_filename = NULL,
                                exclude_patterns = character()) {
  if (is.null(dir) || is.na(dir) || !dir.exists(dir)) {
    return(NA_character_)
  }
  
  if (!is.null(exact_filename) && !is.na(exact_filename)) {
    exact_path <- file.path(dir, exact_filename)
    if (file.exists(exact_path)) {
      return(exact_path)
    }
  }
  
  files <- list.files(dir, full.names = TRUE, recursive = FALSE)
  
  if (length(files) == 0) {
    return(NA_character_)
  }
  
  bn <- basename(files)
  
  if (length(exclude_patterns) > 0) {
    exclude_hit <- rep(FALSE, length(files))
    for (ep in exclude_patterns) {
      exclude_hit <- exclude_hit | grepl(ep, bn, ignore.case = TRUE)
    }
    files <- files[!exclude_hit]
    bn <- basename(files)
  }
  
  if (length(files) == 0) {
    return(NA_character_)
  }
  
  for (pat in patterns) {
    hit <- files[grepl(pat, bn, ignore.case = TRUE)]
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  NA_character_
}


# ------------------------------------------------------------
# Image reading and normalization
# ------------------------------------------------------------

read_gray_image <- function(path, channel = 1) {
  if (is.na(path) || !file.exists(path)) {
    stop("Image file not found: ", path, call. = FALSE)
  }
  
  img <- EBImage::readImage(path)
  arr <- EBImage::imageData(img)
  
  if (length(dim(arr)) == 2) {
    mat <- arr
  } else if (length(dim(arr)) == 3) {
    if (channel < 1 || channel > dim(arr)[3]) {
      stop(
        "Requested channel ", channel,
        " but image has ", dim(arr)[3], " channel(s): ",
        path,
        call. = FALSE
      )
    }
    mat <- arr[, , channel]
  } else if (length(dim(arr)) == 4) {
    if (channel < 1 || channel > dim(arr)[3]) {
      stop(
        "Requested channel ", channel,
        " but image has ", dim(arr)[3], " channel(s): ",
        path,
        call. = FALSE
      )
    }
    mat <- arr[, , channel, 1]
  } else {
    stop("Unsupported image dimensions for: ", path, call. = FALSE)
  }
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}

normalize01 <- function(x) {
  rng <- range(x, finite = TRUE, na.rm = TRUE)
  
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  
  (x - rng[1]) / diff(rng)
}


# ------------------------------------------------------------
# Mask/object helpers
# ------------------------------------------------------------

robust_label_mask <- function(mask_mat) {
  vals <- unique(as.integer(round(mask_mat[is.finite(mask_mat) & mask_mat > 0])))
  vals <- vals[vals > 0]
  
  if (length(vals) > 2) {
    lab <- as.integer(round(mask_mat))
    lab[!is.finite(lab)] <- 0
    dim(lab) <- dim(mask_mat)
    return(lab)
  }
  
  bw <- mask_mat > 0
  lab <- EBImage::bwlabel(bw)
  lab <- as.integer(lab)
  dim(lab) <- dim(mask_mat)
  lab
}

get_label_bbox <- function(label_mat, object_id, pad = 2) {
  coords <- which(label_mat == object_id, arr.ind = TRUE)
  coords <- as.matrix(coords)
  
  if (nrow(coords) == 0 || ncol(coords) < 2) {
    return(NULL)
  }
  
  ymin <- max(1, min(coords[, 1]) - pad)
  ymax <- min(nrow(label_mat), max(coords[, 1]) + pad)
  xmin <- max(1, min(coords[, 2]) - pad)
  xmax <- min(ncol(label_mat), max(coords[, 2]) + pad)
  
  list(
    ymin = ymin,
    ymax = ymax,
    xmin = xmin,
    xmax = xmax,
    area_px = nrow(coords)
  )
}


# ------------------------------------------------------------
# Rosette table handling
# ------------------------------------------------------------

standardize_rosettes <- function(path, accepted_only = TRUE) {
  if (is.na(path) || !file.exists(path)) {
    stop("Rosette table not found: ", path, call. = FALSE)
  }
  
  dt <- fread(path)
  setnames(dt, names(dt), tolower(names(dt)))
  
  if (isTRUE(accepted_only) && "accepted" %in% names(dt)) {
    before_n <- nrow(dt)
    dt <- dt[truthy(accepted)]
    
    message(
      "[03_cells] accepted_only=TRUE: kept ",
      nrow(dt),
      " / ",
      before_n,
      " centers from ",
      basename(path)
    )
  }
  
  if (nrow(dt) == 0) {
    warning("[03_cells] Rosette/center table has zero usable rows after filtering: ", path)
    
    return(data.table(
      rosette_id = character(),
      x = numeric(),
      y = numeric(),
      radius = numeric(),
      object_id = integer()
    ))
  }
  
  pick_col <- function(candidates) {
    hit <- candidates[candidates %in% names(dt)]
    if (length(hit) == 0) {
      return(NA_character_)
    }
    hit[1]
  }
  
  is_step2_center_table <- all(c("center_x", "center_y", "fitted_radius_px") %in% names(dt)) ||
    "parent_blob_id" %in% names(dt)
  
  if (is_step2_center_table) {
    id_col  <- pick_col(c("rosette_candidate_id", "rosette_id", "id", "center_id"))
    row_col <- pick_col(c("center_x", "weighted_center_x", "best_center_x"))
    col_col <- pick_col(c("center_y", "weighted_center_y", "best_center_y"))
    r_col   <- pick_col(c("fitted_radius_px", "weighted_fitted_radius_px", "best_fitted_radius_px", "radius", "r"))
    obj_col <- pick_col(c("parent_blob_id", "object_id", "cluster_id", "clump_id", "mask_id"))
    
    if (is.na(row_col) || is.na(col_col) || is.na(r_col)) {
      stop(
        "Step 2-style center table must contain center_x, center_y, and fitted_radius_px. Found columns: ",
        paste(names(dt), collapse = ", "),
        call. = FALSE
      )
    }
    
    out <- data.table(
      rosette_id = if (!is.na(id_col)) {
        as.character(dt[[id_col]])
      } else {
        paste0("R", seq_len(nrow(dt)))
      },
      x = as.numeric(dt[[col_col]]),
      y = as.numeric(dt[[row_col]]),
      radius = as.numeric(dt[[r_col]])
    )
    
    if (!is.na(obj_col)) {
      out[, object_id := as.integer(dt[[obj_col]])]
    } else {
      out[, object_id := NA_integer_]
    }
    
    if ("confidence" %in% names(dt)) {
      out[, confidence := as.numeric(dt[["confidence"]])]
    }
    
    if ("possible_overlap_break" %in% names(dt)) {
      out[, possible_overlap_break := truthy(dt[["possible_overlap_break"]])]
    }
    
    if ("accepted" %in% names(dt)) {
      out[, accepted := truthy(dt[["accepted"]])]
    }
    
  } else {
    id_col  <- pick_col(c("rosette_id", "rosette", "id", "center_id"))
    x_col   <- pick_col(c("x", "center_x", "cx", "x_center", "centroid_x"))
    y_col   <- pick_col(c("y", "center_y", "cy", "y_center", "centroid_y"))
    r_col   <- pick_col(c("radius", "r", "rad", "estimated_radius", "rosette_radius"))
    obj_col <- pick_col(c("object_id", "cluster_id", "clump_id", "mask_id", "parent_blob_id"))
    
    if (is.na(x_col) || is.na(y_col) || is.na(r_col)) {
      stop(
        "Rosette table must contain x/y/radius columns. Found columns: ",
        paste(names(dt), collapse = ", "),
        call. = FALSE
      )
    }
    
    out <- data.table(
      rosette_id = if (!is.na(id_col)) {
        as.character(dt[[id_col]])
      } else {
        paste0("R", seq_len(nrow(dt)))
      },
      x = as.numeric(dt[[x_col]]),
      y = as.numeric(dt[[y_col]]),
      radius = as.numeric(dt[[r_col]])
    )
    
    if (!is.na(obj_col)) {
      out[, object_id := as.integer(dt[[obj_col]])]
    } else {
      out[, object_id := NA_integer_]
    }
  }
  
  out <- out[is.finite(x) & is.finite(y) & is.finite(radius) & radius > 0]
  
  message(
    "[03_cells] Loaded ",
    nrow(out),
    " rosette/center entries from ",
    basename(path)
  )
  
  out
}


# ------------------------------------------------------------
# Progress and numeric-safe summaries
# ------------------------------------------------------------

write_progress <- function(progress_path, row) {
  row <- as.data.table(row)
  
  fwrite(
    row,
    file = progress_path,
    sep = "\t",
    append = file.exists(progress_path),
    col.names = !file.exists(progress_path)
  )
}

safe_mean_num <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(mean(x, na.rm = TRUE))
}

safe_median_num <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  as.numeric(stats::median(x, na.rm = TRUE))
}

safe_uniqueN_int <- function(x) {
  as.integer(data.table::uniqueN(x))
}


# ------------------------------------------------------------
# Cell segmentation inside cropped object
# ------------------------------------------------------------

segment_cells_in_object_crop <- function(nuc_crop,
                                         object_mask_crop,
                                         min_cell_area_px,
                                         max_cell_area_px,
                                         threshold_quantile,
                                         cell_signal_polarity,
                                         smooth_sigma,
                                         opening_radius) {
  if (!any(object_mask_crop)) {
    return(data.table())
  }
  
  x <- nuc_crop
  x[!object_mask_crop] <- NA_real_
  
  vals <- x[is.finite(x)]
  
  if (length(vals) < min_cell_area_px) {
    return(data.table())
  }
  
  x_norm <- normalize01(x)
  x_norm[!object_mask_crop] <- 0
  
  if (!is.na(smooth_sigma) && smooth_sigma > 0) {
    x_norm <- EBImage::gblur(x_norm, sigma = smooth_sigma)
  }
  
  vals2 <- x_norm[object_mask_crop & is.finite(x_norm)]
  
  if (length(vals2) < min_cell_area_px) {
    return(data.table())
  }
  
  if (cell_signal_polarity == "bright") {
    thr <- as.numeric(stats::quantile(vals2, probs = threshold_quantile, na.rm = TRUE))
    cell_bw <- x_norm >= thr
  } else if (cell_signal_polarity == "dark") {
    thr <- as.numeric(stats::quantile(vals2, probs = threshold_quantile, na.rm = TRUE))
    cell_bw <- x_norm <= thr
  } else {
    stop("Unknown cell_signal_polarity: ", cell_signal_polarity, call. = FALSE)
  }
  
  cell_bw[!object_mask_crop] <- FALSE
  cell_bw[is.na(cell_bw)] <- FALSE
  
  if (!is.na(opening_radius) && opening_radius > 0) {
    brush <- EBImage::makeBrush(size = opening_radius * 2 + 1, shape = "disc")
    cell_bw <- EBImage::opening(cell_bw, brush)
    cell_bw <- cell_bw > 0
  }
  
  lab <- EBImage::bwlabel(cell_bw)
  lab <- as.integer(lab)
  dim(lab) <- dim(cell_bw)
  
  idx <- which(lab > 0, arr.ind = TRUE)
  idx <- as.matrix(idx)
  
  if (nrow(idx) == 0 || ncol(idx) < 2) {
    return(data.table())
  }
  
  dt <- data.table(
    local_label = lab[idx],
    local_y = idx[, 1],
    local_x = idx[, 2]
  )
  
  cells <- dt[
    ,
    .(
      area_px = as.integer(.N),
      local_x = as.numeric(mean(local_x)),
      local_y = as.numeric(mean(local_y)),
      local_xmin = as.numeric(min(local_x)),
      local_xmax = as.numeric(max(local_x)),
      local_ymin = as.numeric(min(local_y)),
      local_ymax = as.numeric(max(local_y))
    ),
    by = local_label
  ]
  
  cells <- cells[
    area_px >= min_cell_area_px &
      area_px <= max_cell_area_px
  ]
  
  cells[, local_label := as.integer(local_label)]
  cells[, area_px := as.numeric(area_px)]
  
  cells
}


# ------------------------------------------------------------
# Rosette filtering and cell assignment
# ------------------------------------------------------------

filter_rosettes_for_object <- function(rosettes, object_id, bbox) {
  if (nrow(rosettes) == 0) {
    return(rosettes)
  }
  
  in_bbox <- (
    rosettes$x >= bbox$xmin &
      rosettes$x <= bbox$xmax &
      rosettes$y >= bbox$ymin &
      rosettes$y <= bbox$ymax
  )
  
  has_object_ids <- "object_id" %in% names(rosettes) && any(!is.na(rosettes$object_id))
  
  if (has_object_ids) {
    by_object <- (!is.na(rosettes$object_id)) & rosettes$object_id == object_id
    rosettes[by_object | in_bbox]
  } else {
    rosettes[in_bbox]
  }
}

assign_cells_to_rosettes <- function(cells,
                                     rosettes,
                                     radius_multiplier,
                                     radius_buffer_px) {
  if (nrow(cells) == 0 || nrow(rosettes) == 0) {
    return(data.table())
  }
  
  assignments <- vector("list", nrow(rosettes))
  
  for (i in seq_len(nrow(rosettes))) {
    r <- rosettes[i]
    
    search_radius <- as.numeric(r$radius) * radius_multiplier + radius_buffer_px
    
    cand <- cells[
      global_x >= r$x - search_radius &
        global_x <= r$x + search_radius &
        global_y >= r$y - search_radius &
        global_y <= r$y + search_radius
    ]
    
    if (nrow(cand) == 0) {
      next
    }
    
    cand <- copy(cand)
    
    cand[
      ,
      dist_to_center := sqrt((global_x - r$x)^2 + (global_y - r$y)^2)
    ]
    
    cand <- cand[dist_to_center <= search_radius]
    
    if (nrow(cand) == 0) {
      next
    }
    
    cand[, rosette_id := as.character(r$rosette_id)]
    cand[, rosette_x := as.numeric(r$x)]
    cand[, rosette_y := as.numeric(r$y)]
    cand[, rosette_radius := as.numeric(r$radius)]
    cand[, assignment_score := as.numeric(abs(dist_to_center - rosette_radius))]
    
    assignments[[i]] <- cand
  }
  
  all <- rbindlist(assignments, fill = TRUE)
  
  if (nrow(all) == 0) {
    return(data.table())
  }
  
  all[, cell_id := as.character(cell_id)]
  all[, rosette_id := as.character(rosette_id)]
  all[, area_px := as.numeric(area_px)]
  all[, global_x := as.numeric(global_x)]
  all[, global_y := as.numeric(global_y)]
  all[, dist_to_center := as.numeric(dist_to_center)]
  all[, rosette_x := as.numeric(rosette_x)]
  all[, rosette_y := as.numeric(rosette_y)]
  all[, rosette_radius := as.numeric(rosette_radius)]
  all[, assignment_score := as.numeric(assignment_score)]
  
  setorder(all, cell_id, assignment_score, dist_to_center)
  unique_assigned <- all[, .SD[1], by = cell_id]
  
  unique_assigned
}


# ------------------------------------------------------------
# Debug overlay
# ------------------------------------------------------------

make_debug_overlay <- function(nuc_crop,
                               object_mask_crop,
                               cells,
                               rosettes_local,
                               out_png,
                               draw_cell_ids = FALSE) {
  png(out_png, width = 1000, height = 1000)
  par(mar = c(1, 1, 2, 1))
  
  img <- normalize01(nuc_crop)
  
  image(
    t(img[nrow(img):1, ]),
    col = gray.colors(256),
    axes = FALSE,
    asp = 1,
    main = basename(out_png)
  )
  
  bidx <- which(object_mask_crop, arr.ind = TRUE)
  bidx <- as.matrix(bidx)
  
  if (nrow(bidx) > 0 && ncol(bidx) >= 2) {
    points(
      x = bidx[, 2] / ncol(img),
      y = 1 - bidx[, 1] / nrow(img),
      pch = ".",
      cex = 0.2
    )
  }
  
  if (nrow(cells) > 0) {
    points(
      x = cells$local_x / ncol(img),
      y = 1 - cells$local_y / nrow(img),
      pch = 16,
      cex = 0.8
    )
    
    if (isTRUE(draw_cell_ids)) {
      text(
        x = cells$local_x / ncol(img),
        y = 1 - cells$local_y / nrow(img),
        labels = cells$cell_id,
        cex = 0.5,
        pos = 3
      )
    }
  }
  
  if (nrow(rosettes_local) > 0) {
    symbols(
      x = rosettes_local$local_x / ncol(img),
      y = 1 - rosettes_local$local_y / nrow(img),
      circles = rosettes_local$radius / max(nrow(img), ncol(img)),
      inches = FALSE,
      add = TRUE,
      lwd = 2
    )
    
    points(
      x = rosettes_local$local_x / ncol(img),
      y = 1 - rosettes_local$local_y / nrow(img),
      pch = 3,
      cex = 1.2,
      lwd = 2
    )
  }
  
  dev.off()
}


# ------------------------------------------------------------
# Project path inference
# ------------------------------------------------------------

infer_project_paths <- function(opt) {
  input <- opt$input
  
  outlines_dir <- opt$outlines_dir
  centers_dir  <- opt$centers_dir
  raw_dir      <- opt$raw_dir
  
  if (!is.null(outlines_dir)) {
    if (is.null(centers_dir)) {
      parent <- dirname(outlines_dir)
      
      centers_dir <- first_existing_dir(c(
        file.path(parent, "centers_v01"),
        file.path(parent, "center_v01")
      ))
    }
    
    if (is.null(raw_dir)) {
      project_guess <- dirname(dirname(outlines_dir))
      
      raw_dir <- first_existing_dir(c(
        file.path(project_guess, "raw_images"),
        file.path(project_guess, "raw")
      ))
    }
    
    return(list(
      outlines_dir = outlines_dir,
      centers_dir = centers_dir,
      raw_dir = raw_dir
    ))
  }
  
  if (is.null(input)) {
    stop("No --input or --outlines_dir provided.", call. = FALSE)
  }
  
  if (!dir.exists(input)) {
    stop("Input directory does not exist: ", input, call. = FALSE)
  }
  
  input_norm <- normalizePath(input, mustWork = FALSE)
  input_base <- basename(input_norm)
  
  outlines_guess <- first_existing_dir(c(
    file.path(input, "results", "outline_v01"),
    file.path(input, "results", "outlines_v01")
  ))
  
  if (!is.null(outlines_guess)) {
    centers_guess <- first_existing_dir(c(
      file.path(input, "results", "centers_v01"),
      file.path(input, "results", "center_v01")
    ))
    
    raw_guess <- first_existing_dir(c(
      file.path(input, "raw_images"),
      file.path(input, "raw")
    ))
    
    return(list(
      outlines_dir = outlines_guess,
      centers_dir = centers_guess,
      raw_dir = raw_guess
    ))
  }
  
  outlines_guess <- first_existing_dir(c(
    file.path(input, "outline_v01"),
    file.path(input, "outlines_v01")
  ))
  
  if (!is.null(outlines_guess)) {
    centers_guess <- first_existing_dir(c(
      file.path(input, "centers_v01"),
      file.path(input, "center_v01")
    ))
    
    project_guess <- dirname(input)
    
    raw_guess <- first_existing_dir(c(
      file.path(project_guess, "raw_images"),
      file.path(project_guess, "raw")
    ))
    
    return(list(
      outlines_dir = outlines_guess,
      centers_dir = centers_guess,
      raw_dir = raw_guess
    ))
  }
  
  if (input_base %in% c("outline_v01", "outlines_v01")) {
    outlines_dir <- input
    
    centers_guess <- first_existing_dir(c(
      file.path(dirname(input), "centers_v01"),
      file.path(dirname(input), "center_v01")
    ))
    
    project_guess <- dirname(dirname(input))
    
    raw_guess <- first_existing_dir(c(
      file.path(project_guess, "raw_images"),
      file.path(project_guess, "raw")
    ))
    
    return(list(
      outlines_dir = outlines_dir,
      centers_dir = centers_guess,
      raw_dir = raw_guess
    ))
  }
  
  parent <- dirname(input)
  parent_base <- basename(normalizePath(parent, mustWork = FALSE))
  
  if (parent_base %in% c("outline_v01", "outlines_v01")) {
    outlines_dir <- parent
    
    centers_guess <- first_existing_dir(c(
      file.path(dirname(parent), "centers_v01"),
      file.path(dirname(parent), "center_v01")
    ))
    
    project_guess <- dirname(dirname(parent))
    
    raw_guess <- first_existing_dir(c(
      file.path(project_guess, "raw_images"),
      file.path(project_guess, "raw")
    ))
    
    return(list(
      outlines_dir = outlines_dir,
      centers_dir = centers_guess,
      raw_dir = raw_guess,
      single_sample_dir = input,
      single_sample = basename(input)
    ))
  }
  
  stop(
    "Could not infer project layout from --input: ", input, "\n",
    "Expected one of:\n",
    "  RosetteCounter/\n",
    "  RosetteCounter/results/\n",
    "  RosetteCounter/results/outline_v01/\n",
    "  RosetteCounter/results/outlines_v01/\n",
    "  RosetteCounter/results/outline_v01/<sample>/\n",
    "  RosetteCounter/results/outlines_v01/<sample>/\n",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Raw image discovery
# ------------------------------------------------------------

find_raw_image_for_sample <- function(sample_id, raw_dir) {
  if (is.null(raw_dir) || is.na(raw_dir) || !dir.exists(raw_dir)) {
    return(NA_character_)
  }
  
  files <- list.files(
    raw_dir,
    full.names = TRUE,
    recursive = FALSE,
    pattern = "\\.(tif|tiff|png|jpg|jpeg)$",
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    return(NA_character_)
  }
  
  bn <- basename(files)
  
  hit <- files[grepl(sample_id, bn, ignore.case = TRUE)]
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  if (length(files) == 1) {
    return(files[1])
  }
  
  NA_character_
}

find_nuclei_or_raw_image <- function(sample_id, outline_dir, raw_dir, exact_nuclei_file = NULL) {
  nuclei_path <- find_first_existing(
    dir = outline_dir,
    patterns = nuclei_patterns,
    exact_filename = exact_nuclei_file,
    exclude_patterns = exclude_image_patterns
  )
  
  if (!is.na(nuclei_path)) {
    return(nuclei_path)
  }
  
  raw_hit <- find_raw_image_for_sample(sample_id, raw_dir)
  
  if (!is.na(raw_hit)) {
    return(raw_hit)
  }
  
  NA_character_
}


# ------------------------------------------------------------
# Sample job discovery
# ------------------------------------------------------------

discover_sample_jobs <- function(opt) {
  paths <- infer_project_paths(opt)
  
  outlines_dir <- paths$outlines_dir
  centers_dir  <- paths$centers_dir
  raw_dir      <- paths$raw_dir
  
  msg("[03_cells] inferred outlines_dir: ", outlines_dir)
  msg("[03_cells] inferred centers_dir:  ", centers_dir)
  msg("[03_cells] inferred raw_dir:      ", raw_dir)
  
  if (is.null(outlines_dir) || !dir.exists(outlines_dir)) {
    stop("Outlines directory does not exist: ", outlines_dir, call. = FALSE)
  }
  
  if (is.null(centers_dir) || !dir.exists(centers_dir)) {
    stop("Centers directory does not exist: ", centers_dir, call. = FALSE)
  }
  
  outline_sample_dirs <- list.dirs(outlines_dir, recursive = FALSE, full.names = TRUE)
  center_sample_dirs  <- list.dirs(centers_dir, recursive = FALSE, full.names = TRUE)
  
  if (!is.null(paths$single_sample)) {
    outline_sample_dirs <- paths$single_sample_dir
  }
  
  outline_samples <- basename(outline_sample_dirs)
  center_samples  <- basename(center_sample_dirs)
  
  common_samples <- intersect(outline_samples, center_samples)
  
  if (!is.null(opt$sample)) {
    common_samples <- common_samples[common_samples == opt$sample]
  }
  
  if (length(common_samples) == 0) {
    message("[03_cells] No paired sample folders found.")
    message("[03_cells] Outline sample folders:")
    
    if (length(outline_samples) == 0) {
      message("  <none>")
    } else {
      message(paste0("  - ", outline_samples, collapse = "\n"))
    }
    
    message("[03_cells] Center sample folders:")
    
    if (length(center_samples) == 0) {
      message("  <none>")
    } else {
      message(paste0("  - ", center_samples, collapse = "\n"))
    }
    
    stop(
      "No paired sample folders between outline_v01/outlines_v01 and centers_v01.",
      call. = FALSE
    )
  }
  
  jobs <- vector("list", length(common_samples))
  
  for (i in seq_along(common_samples)) {
    sample_id <- common_samples[i]
    
    outline_dir <- file.path(outlines_dir, sample_id)
    center_dir  <- file.path(centers_dir, sample_id)
    
    object_mask_path <- find_first_existing(
      dir = outline_dir,
      patterns = object_mask_patterns,
      exact_filename = opt$object_mask_file,
      exclude_patterns = exclude_image_patterns
    )
    
    nuclei_path <- find_nuclei_or_raw_image(
      sample_id = sample_id,
      outline_dir = outline_dir,
      raw_dir = raw_dir,
      exact_nuclei_file = opt$nuclei_file
    )
    
    rosette_path <- find_first_existing(
      dir = center_dir,
      patterns = rosette_patterns,
      exact_filename = opt$rosette_file,
      exclude_patterns = c("\\.raw\\.", "overlay", "summary")
    )
    
    jobs[[i]] <- data.table(
      image_id = sample_id,
      outline_dir = outline_dir,
      center_dir = center_dir,
      raw_dir = ifelse(is.null(raw_dir), NA_character_, raw_dir),
      object_mask_path = object_mask_path,
      nuclei_path = nuclei_path,
      rosette_path = rosette_path
    )
  }
  
  jobs <- rbindlist(jobs, fill = TRUE)
  
  msg("[03_cells] Discovered sample jobs:")
  for (i in seq_len(nrow(jobs))) {
    msg("  sample:      ", jobs$image_id[i])
    msg("    outline:  ", jobs$outline_dir[i])
    msg("    center:   ", jobs$center_dir[i])
    msg("    mask:     ", jobs$object_mask_path[i])
    msg("    nuclei:   ", jobs$nuclei_path[i])
    msg("    rosettes: ", jobs$rosette_path[i])
  }
  
  valid <- !is.na(jobs$object_mask_path) &
    !is.na(jobs$nuclei_path) &
    !is.na(jobs$rosette_path)
  
  if (!any(valid)) {
    message("[03_cells] Found paired sample folders, but no complete valid jobs.")
    message("[03_cells] This usually means filename patterns are not matching.")
    
    message("\n[03_cells] Example files in outline folders:")
    for (d in head(jobs$outline_dir, 3)) {
      message("  ", d)
      ff <- list.files(d, recursive = FALSE)
      if (length(ff) == 0) {
        message("    <none>")
      } else {
        message(paste0("    - ", ff, collapse = "\n"))
      }
    }
    
    message("\n[03_cells] Example files in center folders:")
    for (d in head(jobs$center_dir, 3)) {
      message("  ", d)
      ff <- list.files(d, recursive = FALSE)
      if (length(ff) == 0) {
        message("    <none>")
      } else {
        message(paste0("    - ", ff, collapse = "\n"))
      }
    }
    
    stop(
      "No valid jobs. Need object mask from outline_v01/<sample>, gray/nuclei image, and candidate_centers.tsv from centers_v01/<sample>.",
      call. = FALSE
    )
  }
  
  jobs <- jobs[valid]
  
  if (!is.na(opt$max_images) && nrow(jobs) > opt$max_images) {
    jobs <- jobs[seq_len(opt$max_images)]
  }
  
  jobs
}

image_jobs <- discover_sample_jobs(opt)

msg("[03_cells] Found ", nrow(image_jobs), " image job(s).")

if (isTRUE(opt$dry_run)) {
  msg("[03_cells] dry run requested; exiting before image processing.")
  quit(save = "no", status = 0)
}


# ------------------------------------------------------------
# Main processing loop
# ------------------------------------------------------------

all_image_summaries <- list()

for (job_i in seq_len(nrow(image_jobs))) {
  job <- image_jobs[job_i]
  
  image_id <- job$image_id
  image_dir <- job$outline_dir
  center_dir <- job$center_dir
  
  object_mask_path <- job$object_mask_path
  nuclei_path <- job$nuclei_path
  rosette_path <- job$rosette_path
  
  msg("[03_cells] Processing: ", image_id)
  msg("[03_cells] outline dir: ", image_dir)
  msg("[03_cells] center dir:  ", center_dir)
  msg("[03_cells] object mask: ", object_mask_path)
  msg("[03_cells] nuclei/raw:  ", nuclei_path)
  msg("[03_cells] rosette TSV: ", rosette_path)
  
  image_out <- file.path(opt$out, image_id)
  debug_out <- file.path(image_out, "debug_overlays")
  
  dir.create(image_out, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(opt$write_debug_overlays)) {
    dir.create(debug_out, recursive = TRUE, showWarnings = FALSE)
  }
  
  progress_path <- file.path(image_out, "step3_progress.tsv")
  
  if (file.exists(progress_path)) {
    file.remove(progress_path)
  }
  
  if (is.na(object_mask_path) || is.na(nuclei_path) || is.na(rosette_path)) {
    warning("[03_cells] Skipping ", image_id, ": missing required files.")
    next
  }
  
  tt <- tic("read inputs")
  
  object_mask_raw <- read_gray_image(
    object_mask_path,
    channel = opt$mask_channel
  )
  
  nuclei_img <- read_gray_image(
    nuclei_path,
    channel = opt$nuclei_channel
  )
  
  rosettes <- standardize_rosettes(
    rosette_path,
    accepted_only = opt$accepted_only
  )
  
  toc(tt)
  
  if (!identical(dim(object_mask_raw), dim(nuclei_img))) {
    stop(
      "Object mask and nuclei image dimensions do not match for ",
      image_id,
      ". object mask dim=",
      paste(dim(object_mask_raw), collapse = "x"),
      "; nuclei dim=",
      paste(dim(nuclei_img), collapse = "x"),
      "\nobject mask: ",
      object_mask_path,
      "\nnuclei/raw: ",
      nuclei_path,
      call. = FALSE
    )
  }
  
  if (!is.na(opt$max_rosettes) && nrow(rosettes) > opt$max_rosettes) {
    rosettes <- rosettes[seq_len(opt$max_rosettes)]
    msg("[03_cells] Debug mode: keeping first ", opt$max_rosettes, " rosettes.")
  }
  
  tt <- tic("label object mask")
  
  object_labels <- robust_label_mask(object_mask_raw)
  object_ids <- sort(unique(as.integer(object_labels[object_labels > 0])))
  
  toc(tt)
  
  msg("[03_cells] Found ", length(object_ids), " object(s) before filters.")
  
  if (!is.na(opt$max_objects) && length(object_ids) > opt$max_objects) {
    object_ids <- object_ids[seq_len(opt$max_objects)]
    msg("[03_cells] Debug mode: keeping first ", opt$max_objects, " objects.")
  }
  
  cells_all <- list()
  assignments_all <- list()
  rosette_counts_all <- list()
  
  cell_counter <- 0L
  processed_counter <- 0L
  skipped_counter <- 0L
  overlay_counter <- 0L
  
  for (object_i in seq_along(object_ids)) {
    object_id <- object_ids[object_i]
    object_start_time <- Sys.time()
    
    bbox <- get_label_bbox(object_labels, object_id, pad = 3)
    
    if (is.null(bbox)) {
      skipped_counter <- skipped_counter + 1L
      next
    }
    
    area_px <- bbox$area_px
    skip_reason <- NA_character_
    
    if (area_px < opt$min_object_area_px) {
      skip_reason <- paste0("area_below_min_", opt$min_object_area_px)
    }
    
    if (!is.na(opt$max_object_area_px) && area_px > opt$max_object_area_px) {
      skip_reason <- paste0("area_above_max_", opt$max_object_area_px)
    }
    
    if (!is.na(skip_reason)) {
      skipped_counter <- skipped_counter + 1L
      
      write_progress(
        progress_path,
        data.table(
          image_id = as.character(image_id),
          object_i = as.integer(object_i),
          object_id = as.integer(object_id),
          status = "skipped",
          reason = as.character(skip_reason),
          area_px = as.numeric(area_px),
          n_rosettes_in_crop = as.integer(0L),
          n_cells = as.integer(0L),
          n_assigned = as.integer(0L),
          elapsed_sec = as.numeric(0),
          time = as.character(Sys.time())
        )
      )
      
      next
    }
    
    msg(
      "[03_cells] object ", object_i, "/", length(object_ids),
      " | object_id=", object_id,
      " | area=", area_px
    )
    
    object_crop_labels <- object_labels[
      bbox$ymin:bbox$ymax,
      bbox$xmin:bbox$xmax,
      drop = FALSE
    ]
    
    object_mask_crop <- object_crop_labels == object_id
    
    nuclei_crop <- nuclei_img[
      bbox$ymin:bbox$ymax,
      bbox$xmin:bbox$xmax,
      drop = FALSE
    ]
    
    ros_crop <- filter_rosettes_for_object(
      rosettes = rosettes,
      object_id = object_id,
      bbox = bbox
    )
    
    if (nrow(ros_crop) == 0) {
      skipped_counter <- skipped_counter + 1L
      
      elapsed <- round(as.numeric(difftime(Sys.time(), object_start_time, units = "secs")), 3)
      
      write_progress(
        progress_path,
        data.table(
          image_id = as.character(image_id),
          object_i = as.integer(object_i),
          object_id = as.integer(object_id),
          status = "skipped",
          reason = "no_rosettes_in_object_bbox",
          area_px = as.numeric(area_px),
          n_rosettes_in_crop = as.integer(0L),
          n_cells = as.integer(0L),
          n_assigned = as.integer(0L),
          elapsed_sec = as.numeric(elapsed),
          time = as.character(Sys.time())
        )
      )
      
      next
    }
    
    ros_crop <- copy(ros_crop)
    
    ros_crop[, local_x := as.numeric(x - bbox$xmin + 1)]
    ros_crop[, local_y := as.numeric(y - bbox$ymin + 1)]
    
    tt <- tic(paste0("segment cells in object ", object_id))
    
    cells <- segment_cells_in_object_crop(
      nuc_crop = nuclei_crop,
      object_mask_crop = object_mask_crop,
      min_cell_area_px = opt$min_cell_area_px,
      max_cell_area_px = opt$max_cell_area_px,
      threshold_quantile = opt$threshold_quantile,
      cell_signal_polarity = opt$cell_signal_polarity,
      smooth_sigma = opt$smooth_sigma,
      opening_radius = opt$opening_radius
    )
    
    toc(tt)
    
    if (nrow(cells) > 0) {
      cells[, object_id := as.integer(object_id)]
      cells[, image_id := as.character(image_id)]
      
      cells[, global_x := as.numeric(local_x + bbox$xmin - 1)]
      cells[, global_y := as.numeric(local_y + bbox$ymin - 1)]
      
      cells[, cell_id := as.character(paste0(image_id, "_obj", object_id, "_cell", seq_len(.N) + cell_counter))]
      
      cells[, local_label := as.integer(local_label)]
      cells[, area_px := as.numeric(area_px)]
      cells[, local_x := as.numeric(local_x)]
      cells[, local_y := as.numeric(local_y)]
      cells[, local_xmin := as.numeric(local_xmin)]
      cells[, local_xmax := as.numeric(local_xmax)]
      cells[, local_ymin := as.numeric(local_ymin)]
      cells[, local_ymax := as.numeric(local_ymax)]
      
      cell_counter <- cell_counter + nrow(cells)
      
      setcolorder(
        cells,
        c(
          "image_id",
          "object_id",
          "cell_id",
          "local_label",
          "area_px",
          "global_x",
          "global_y",
          "local_x",
          "local_y",
          "local_xmin",
          "local_xmax",
          "local_ymin",
          "local_ymax"
        )
      )
    }
    
    tt <- tic(paste0("assign cells in object ", object_id))
    
    assigned <- assign_cells_to_rosettes(
      cells = cells,
      rosettes = ros_crop,
      radius_multiplier = opt$radius_multiplier,
      radius_buffer_px = opt$radius_buffer_px
    )
    
    toc(tt)
    
    if (nrow(assigned) > 0) {
      assigned[, image_id := as.character(image_id)]
      assigned[, object_id := as.integer(object_id)]
      assigned[, rosette_id := as.character(rosette_id)]
      assigned[, cell_id := as.character(cell_id)]
      
      assigned[, global_x := as.numeric(global_x)]
      assigned[, global_y := as.numeric(global_y)]
      assigned[, area_px := as.numeric(area_px)]
      assigned[, dist_to_center := as.numeric(dist_to_center)]
      assigned[, rosette_x := as.numeric(rosette_x)]
      assigned[, rosette_y := as.numeric(rosette_y)]
      assigned[, rosette_radius := as.numeric(rosette_radius)]
      assigned[, assignment_score := as.numeric(assignment_score)]
      
      setcolorder(
        assigned,
        c(
          "image_id",
          "object_id",
          "rosette_id",
          "cell_id",
          "global_x",
          "global_y",
          "area_px",
          "dist_to_center",
          "rosette_x",
          "rosette_y",
          "rosette_radius",
          "assignment_score"
        )
      )
    }
    
    if (nrow(assigned) > 0) {
      assigned <- copy(assigned)
      
      assigned[, image_id := as.character(image_id)]
      assigned[, object_id := as.integer(object_id)]
      assigned[, rosette_id := as.character(rosette_id)]
      assigned[, rosette_x := as.numeric(rosette_x)]
      assigned[, rosette_y := as.numeric(rosette_y)]
      assigned[, rosette_radius := as.numeric(rosette_radius)]
      assigned[, area_px := as.numeric(area_px)]
      assigned[, dist_to_center := as.numeric(dist_to_center)]
      
      counts <- assigned[
        ,
        .(
          n_cells = safe_uniqueN_int(cell_id),
          mean_cell_area_px = safe_mean_num(area_px),
          median_cell_area_px = safe_median_num(area_px),
          mean_dist_to_center = safe_mean_num(dist_to_center),
          median_dist_to_center = safe_median_num(dist_to_center)
        ),
        by = .(
          image_id,
          object_id,
          rosette_id,
          rosette_x,
          rosette_y,
          rosette_radius
        )
      ]
      
      counts[, image_id := as.character(image_id)]
      counts[, object_id := as.integer(object_id)]
      counts[, rosette_id := as.character(rosette_id)]
      counts[, rosette_x := as.numeric(rosette_x)]
      counts[, rosette_y := as.numeric(rosette_y)]
      counts[, rosette_radius := as.numeric(rosette_radius)]
      counts[, n_cells := as.integer(n_cells)]
      counts[, mean_cell_area_px := as.numeric(mean_cell_area_px)]
      counts[, median_cell_area_px := as.numeric(median_cell_area_px)]
      counts[, mean_dist_to_center := as.numeric(mean_dist_to_center)]
      counts[, median_dist_to_center := as.numeric(median_dist_to_center)]
      
    } else {
      counts <- data.table(
        image_id = as.character(image_id),
        object_id = as.integer(object_id),
        rosette_id = as.character(ros_crop$rosette_id),
        rosette_x = as.numeric(ros_crop$x),
        rosette_y = as.numeric(ros_crop$y),
        rosette_radius = as.numeric(ros_crop$radius),
        n_cells = as.integer(0L),
        mean_cell_area_px = as.numeric(NA_real_),
        median_cell_area_px = as.numeric(NA_real_),
        mean_dist_to_center = as.numeric(NA_real_),
        median_dist_to_center = as.numeric(NA_real_)
      )
    }
    
    if (isTRUE(opt$write_debug_overlays) && overlay_counter < opt$debug_draw_n) {
      overlay_counter <- overlay_counter + 1L
      
      overlay_path <- file.path(
        debug_out,
        paste0(
          image_id,
          "_object_",
          sprintf("%04d", object_id),
          "_cells.png"
        )
      )
      
      make_debug_overlay(
        nuc_crop = nuclei_crop,
        object_mask_crop = object_mask_crop,
        cells = cells,
        rosettes_local = ros_crop,
        out_png = overlay_path,
        draw_cell_ids = opt$draw_cell_ids
      )
    }
    
    processed_counter <- processed_counter + 1L
    
    cells_all[[length(cells_all) + 1L]] <- cells
    assignments_all[[length(assignments_all) + 1L]] <- assigned
    rosette_counts_all[[length(rosette_counts_all) + 1L]] <- counts
    
    elapsed <- round(as.numeric(difftime(Sys.time(), object_start_time, units = "secs")), 3)
    
    write_progress(
      progress_path,
      data.table(
        image_id = as.character(image_id),
        object_i = as.integer(object_i),
        object_id = as.integer(object_id),
        status = "processed",
        reason = NA_character_,
        area_px = as.numeric(area_px),
        n_rosettes_in_crop = as.integer(nrow(ros_crop)),
        n_cells = as.integer(nrow(cells)),
        n_assigned = as.integer(nrow(assigned)),
        elapsed_sec = as.numeric(elapsed),
        time = as.character(Sys.time())
      )
    )
  }
  
  cells_dt <- rbindlist(cells_all, fill = TRUE)
  assigned_dt <- rbindlist(assignments_all, fill = TRUE)
  rosette_counts_dt <- rbindlist(rosette_counts_all, fill = TRUE)
  
  cells_path <- file.path(image_out, "cells.tsv")
  assigned_path <- file.path(image_out, "cell_rosette_assignments.tsv")
  counts_path <- file.path(image_out, "rosette_cell_counts.tsv")
  
  fwrite(cells_dt, cells_path, sep = "\t")
  fwrite(assigned_dt, assigned_path, sep = "\t")
  fwrite(rosette_counts_dt, counts_path, sep = "\t")
  
  image_summary <- data.table(
    image_id = as.character(image_id),
    n_objects_total = as.integer(length(unique(as.integer(object_labels[object_labels > 0])))),
    n_objects_considered = as.integer(length(object_ids)),
    n_objects_processed = as.integer(processed_counter),
    n_objects_skipped = as.integer(skipped_counter),
    n_rosettes = as.integer(nrow(rosettes)),
    n_cells_detected = as.integer(nrow(cells_dt)),
    n_cells_assigned = as.integer(if (nrow(assigned_dt) > 0) uniqueN(assigned_dt$cell_id) else 0L),
    cells_path = as.character(cells_path),
    assignments_path = as.character(assigned_path),
    counts_path = as.character(counts_path)
  )
  
  all_image_summaries[[length(all_image_summaries) + 1L]] <- image_summary
  
  fwrite(
    image_summary,
    file.path(image_out, "image_summary.tsv"),
    sep = "\t"
  )
  
  msg("[03_cells] Finished: ", image_id)
  msg("[03_cells]   cells: ", cells_path)
  msg("[03_cells]   assignments: ", assigned_path)
  msg("[03_cells]   counts: ", counts_path)
}

summary_dt <- rbindlist(all_image_summaries, fill = TRUE)
summary_path <- file.path(opt$out, "step3_summary.tsv")
fwrite(summary_dt, summary_path, sep = "\t")

msg("[03_cells] Complete.")
msg("[03_cells] Summary: ", summary_path)