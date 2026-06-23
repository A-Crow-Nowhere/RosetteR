#!/usr/bin/env Rscript

# 03b_segment_cells_by_cluster_geometry.R
#
# SCRIPT_VERSION: 03b_geometry_v12_sample_folder_diagnostics
#
# Alternate Step 3 for rosette/cell analysis.
#
# This version segments cells from the geometry of each Step 1 cluster mask.
# It is membrane-aware: dark membrane lines in gray.png are detected as local
# dark ridges and can be used as internal cut lines before watershed subdivision.
#
# Required R packages:
#   optparse, data.table, png, EBImage
#
# Bioconductor install example:
#   install.packages(c("optparse", "data.table", "png"))
#   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#   BiocManager::install("EBImage")

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Missing R package: optparse")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Missing R package: data.table")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Missing R package: png")
  }
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop("Missing Bioconductor package: EBImage")
  }
  
  library(optparse)
  library(data.table)
  library(png)
  library(EBImage)
})

# ----------------------------- arguments -------------------------------------

option_list <- list(
  make_option("--input", type = "character", default = ".",
              help = "Input folder. Can be a single image result folder or a parent folder containing result folders."),
  make_option("--out", type = "character", default = "results/03b_geometry_cells",
              help = "Output folder for alternate Step 3 results."),
  make_option("--image_id_override", type = "character", default = NA,
              help = "Optional image ID to use for output naming and filtering combined candidate-center tables."),
  
  make_option("--gray_name", type = "character", default = "gray.png",
              help = "Name of the grayscale image inside each result folder."),
  make_option("--cluster_mask_name", type = "character", default = NA,
              help = "Name of the Step 1 cluster/object mask. If omitted, the script searches common mask names."),
  make_option("--rosette_table_name", type = "character", default = NA,
              help = "Name of the Step 2 rosette center/radius TSV. If omitted, the script searches common rosette TSV names."),
  
  make_option("--gray_path", type = "character", default = NA,
              help = "Explicit path to gray image for single-folder runs. Overrides --gray_dir/--gray_name."),
  make_option("--cluster_mask_path", type = "character", default = NA,
              help = "Explicit path to cluster/object mask for single-folder runs. Overrides --cluster_mask_dir/--cluster_mask_name."),
  make_option("--rosette_table_path", type = "character", default = NA,
              help = "Explicit path to rosette/candidate-center TSV for single-folder runs. Overrides --rosette_table_dir/--rosette_table_name."),
  
  make_option("--gray_dir", type = "character", default = NA,
              help = "Optional base directory containing gray images. Batch mode tries <gray_dir>/<image_id>/<gray_name>, then <gray_dir>/<gray_name>."),
  make_option("--cluster_mask_dir", type = "character", default = NA,
              help = "Optional base directory containing cluster masks. Batch mode tries <cluster_mask_dir>/<image_id>/<cluster_mask_name>, then <cluster_mask_dir>/<cluster_mask_name>."),
  make_option("--rosette_table_dir", type = "character", default = NA,
              help = "Optional base directory containing rosette/candidate-center tables. Batch mode tries <rosette_table_dir>/<image_id>/<rosette_table_name>, then <rosette_table_dir>/<rosette_table_name>."),
  
  make_option("--candidate_center_source", type = "character", default = "weighted",
              help = "Which candidate-center columns to use when present: weighted, best, raw, or auto. For all_candidate_centers.tsv, weighted uses weighted_center_x/y and weighted_fitted_radius_px."),
  make_option("--accepted_only", type = "character", default = "TRUE",
              help = "TRUE/FALSE. If the rosette table has an accepted column, only use accepted rows."),
  make_option("--require_center_inside_blob", type = "character", default = "FALSE",
              help = "TRUE/FALSE. If matching columns exist, require selected center to be inside blob."),
  make_option("--min_candidate_confidence", type = "double", default = NA,
              help = "Optional minimum confidence threshold for candidate-center rows."),
  
  make_option("--mask_foreground", type = "character", default = "auto",
              help = "Mask foreground polarity: auto, bright, or dark."),
  make_option("--fill_cluster_holes", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Fill holes in the cluster mask before subdivision."),
  make_option("--min_cluster_area_px", type = "double", default = 100,
              help = "Minimum connected-component cluster area to process."),
  
  make_option("--use_membrane_cuts", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Use dark membrane lines from gray.png as cut lines before watershed."),
  make_option("--membrane_bg_radius_px", type = "double", default = 9,
              help = "Local median background radius for dark membrane enhancement."),
  make_option("--membrane_quantile", type = "double", default = 0.82,
              help = "Quantile of local dark-ridge score inside each cluster used as membrane threshold."),
  make_option("--membrane_min_score", type = "double", default = 0.025,
              help = "Minimum dark-ridge score required for a membrane candidate pixel."),
  make_option("--membrane_open_px", type = "double", default = 0,
              help = "Opening radius for membrane mask cleanup. Use 0 to disable."),
  make_option("--membrane_dilate_px", type = "double", default = 1,
              help = "Dilation radius for membrane cut lines."),
  
  make_option("--min_cell_area_px", type = "double", default = 50,
              help = "Minimum segmented cell-object area."),
  make_option("--max_cell_area_px", type = "double", default = 4000,
              help = "Maximum segmented cell-object area."),
  make_option("--min_cell_radius_px", type = "double", default = 2,
              help = "Minimum object max distance-transform radius."),
  make_option("--max_cell_radius_px", type = "double", default = 80,
              help = "Maximum object max distance-transform radius."),
  
  make_option("--seed_smooth_sigma", type = "double", default = 1.0,
              help = "Gaussian smoothing sigma applied to distance map before watershed."),
  make_option("--seed_min_distance_px", type = "double", default = 5,
              help = "Watershed extent/minimum spacing control. Larger values reduce over-splitting."),
  make_option("--watershed_tolerance", type = "double", default = 1.0,
              help = "Watershed tolerance. Larger values merge weak local maxima and reduce over-splitting."),
  
  make_option("--min_solidity", type = "double", default = 0.35,
              help = "Minimum solidity for accepted complete-ish cells."),
  make_option("--min_circularity", type = "double", default = 0.10,
              help = "Minimum circularity/compactness for accepted complete-ish cells."),
  make_option("--max_edge_contact_fraction", type = "double", default = 0.75,
              help = "Maximum fraction of cell boundary allowed to contact the outer cluster boundary."),
  
  make_option("--assignment_max_norm_distance", type = "double", default = 1.60,
              help = "Maximum centroid-to-rosette distance normalized by rosette radius for assignment."),
  make_option("--assignment_radius_weight", type = "double", default = 1.0,
              help = "Weight for normalized center distance in assignment score."),
  make_option("--assignment_overlap_weight", type = "double", default = 0.5,
              help = "Weight for cell overlap with rosette disk in assignment score."),
  make_option("--assignment_min_overlap", type = "double", default = 0.0,
              help = "Minimum overlap fraction that can rescue a cell assignment even if normalized distance is high."),
  
  make_option("--draw_cell_ids", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Draw cell IDs in overlay."),
  make_option("--draw_rejected", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Draw rejected fragments in overlay."),
  make_option("--label_mode", type = "character", default = "both",
              help = "Overlay label mode: none, cell, rosette, or both."),
  make_option("--overlay_scale", type = "double", default = 4,
              help = "Scale factor for PNG debug overlays."),
  make_option("--boundary_line_width_px", type = "double", default = 1,
              help = "Pixel dilation radius for drawn cell boundaries."),
  make_option("--debug", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Write intermediate membrane/distance/debug images."),
  make_option("--debug_parser", type = "character", default = "TRUE",
              help = "TRUE/FALSE. Print detailed candidate-center table parsing diagnostics."),
  make_option("--keep_going", type = "character", default = "FALSE",
              help = "TRUE/FALSE. Continue after an image fails. Default FALSE stops at the first real error for debugging."),
  
  make_option("--max_folders", type = "integer", default = NA,
              help = "Optional limit on number of image folders to process for fast debugging.")
)

args <- parse_args(OptionParser(option_list = option_list))

# ----------------------------- small helpers ---------------------------------

is_na_arg <- function(x) {
  length(x) == 0 || is.na(x) || identical(x, "NA") || identical(x, "") || identical(tolower(as.character(x)), "null")
}

parse_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}

args$fill_cluster_holes <- parse_bool(args$fill_cluster_holes)
args$use_membrane_cuts <- parse_bool(args$use_membrane_cuts)
args$draw_cell_ids <- parse_bool(args$draw_cell_ids)
args$draw_rejected <- parse_bool(args$draw_rejected)
args$debug <- parse_bool(args$debug)
args$debug_parser <- parse_bool(args$debug_parser)
args$keep_going <- parse_bool(args$keep_going)
args$accepted_only <- parse_bool(args$accepted_only)
args$require_center_inside_blob <- parse_bool(args$require_center_inside_blob)
args$candidate_center_source <- tolower(trimws(as.character(args$candidate_center_source)))

dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
message("[03b_geometry] SCRIPT_VERSION: 03b_geometry_v12_sample_folder_diagnostics")
message("[03b_geometry] INPUT_MODE: sample-folder-only; expects results/<experiment>/<samplefolder>/{gray,mask,centers}")
message("[03b_geometry] DISCOVERY: no recursive search; direct sample folders only")

messagef <- function(...) {
  message(sprintf(...))
}

safe_file_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

odd_size_from_radius <- function(radius_px) {
  s <- max(3L, as.integer(2L * ceiling(radius_px) + 1L))
  if (s %% 2L == 0L) s <- s + 1L
  s
}

normalize01 <- function(x) {
  x <- as.matrix(x)
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!all(is.finite(rng)) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

plain_matrix <- function(x, nr = NULL, nc = NULL, label = "matrix") {
  # EBImage sometimes returns Image/S4 arrays whose orientation can dispatch
  # oddly under arithmetic/logical operators. This helper converts to a plain
  # base R matrix and, when requested, forces it to match a reference dimension.
  if (inherits(x, "Image")) {
    x <- as.array(x)
  }
  
  if (length(dim(x)) == 3L) {
    x <- x[, , 1]
  }
  
  m <- as.matrix(x)
  
  if (!is.null(nr) && !is.null(nc)) {
    if (identical(dim(m), c(as.integer(nr), as.integer(nc)))) {
      return(matrix(as.vector(m), nrow = nr, ncol = nc))
    }
    
    if (identical(dim(m), c(as.integer(nc), as.integer(nr)))) {
      mt <- t(m)
      return(matrix(as.vector(mt), nrow = nr, ncol = nc))
    }
    
    if (length(m) == nr * nc) {
      return(matrix(as.vector(m), nrow = nr, ncol = nc))
    }
    
    stop(sprintf(
      "%s has dimensions %s but expected %d x %d",
      label,
      paste(dim(m), collapse = " x "),
      nr,
      nc
    ))
  }
  
  matrix(as.vector(m), nrow = nrow(m), ncol = ncol(m))
}

plain_logical_matrix <- function(x, nr = NULL, nc = NULL, label = "logical matrix") {
  m <- plain_matrix(x, nr = nr, nc = nc, label = label)
  matrix(as.logical(m), nrow = nrow(m), ncol = ncol(m))
}

make_disc <- function(radius_px) {
  radius_px <- as.integer(ceiling(radius_px))
  if (radius_px <= 0L) {
    return(matrix(1, nrow = 1, ncol = 1))
  }
  EBImage::makeBrush(size = 2L * radius_px + 1L, shape = "disc")
}

m_dilate <- function(mask, radius_px) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "dilate input")
  if (radius_px <= 0) return(mask)
  plain_logical_matrix(EBImage::dilate(EBImage::Image(mask), make_disc(radius_px)), nr = nr, nc = nc, label = "dilate output")
}

m_erode <- function(mask, radius_px) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "erode input")
  if (radius_px <= 0) return(mask)
  plain_logical_matrix(EBImage::erode(EBImage::Image(mask), make_disc(radius_px)), nr = nr, nc = nc, label = "erode output")
}

m_opening <- function(mask, radius_px) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "opening input")
  if (radius_px <= 0) return(mask)
  plain_logical_matrix(EBImage::opening(EBImage::Image(mask), make_disc(radius_px)), nr = nr, nc = nc, label = "opening output")
}

m_closing <- function(mask, radius_px) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "closing input")
  if (radius_px <= 0) return(mask)
  plain_logical_matrix(EBImage::closing(EBImage::Image(mask), make_disc(radius_px)), nr = nr, nc = nc, label = "closing output")
}

fill_hull <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "fillHull input")
  plain_logical_matrix(EBImage::fillHull(EBImage::Image(mask)), nr = nr, nc = nc, label = "fillHull output")
}

bwlabel_matrix <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "bwlabel input")
  lab <- plain_matrix(EBImage::bwlabel(EBImage::Image(mask)), nr = nr, nc = nc, label = "bwlabel output")
  storage.mode(lab) <- "integer"
  lab
}

relabel_sequential <- function(lab) {
  lab <- as.matrix(lab)
  ids <- sort(unique(as.integer(lab[lab > 0])))
  out <- matrix(0L, nrow = nrow(lab), ncol = ncol(lab))
  if (length(ids) == 0L) return(out)
  for (i in seq_along(ids)) {
    out[lab == ids[i]] <- as.integer(i)
  }
  out
}

remove_small_components <- function(mask, min_area_px) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "remove_small_components input")
  lab <- bwlabel_matrix(mask)
  ids <- sort(unique(as.integer(lab[lab > 0])))
  if (length(ids) == 0L) return(matrix(FALSE, nrow = nr, ncol = nc))
  tab <- tabulate(lab[lab > 0], nbins = max(lab))
  keep <- which(tab >= min_area_px)
  matrix(lab %in% keep, nrow = nr, ncol = nc)
}

boundary_from_mask <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  mask <- plain_logical_matrix(mask, nr = nr, nc = nc, label = "boundary_from_mask input")
  if (!any(mask)) return(matrix(FALSE, nrow = nr, ncol = nc))
  er <- m_erode(mask, 1)
  er <- plain_logical_matrix(er, nr = nr, nc = nc, label = "boundary_from_mask eroded")
  matrix(mask & !er, nrow = nr, ncol = nc)
}

boundary_from_label <- function(lab) {
  lab <- as.matrix(lab)
  nr <- nrow(lab)
  nc <- ncol(lab)
  b <- matrix(FALSE, nr, nc)
  
  if (nc > 1L) {
    diff_h <- lab[, -1L, drop = FALSE] != lab[, -nc, drop = FALSE]
    b[, -1L] <- b[, -1L] | (diff_h & lab[, -1L, drop = FALSE] > 0)
    b[, -nc] <- b[, -nc] | (diff_h & lab[, -nc, drop = FALSE] > 0)
  }
  if (nr > 1L) {
    diff_v <- lab[-1L, , drop = FALSE] != lab[-nr, , drop = FALSE]
    b[-1L, ] <- b[-1L, ] | (diff_v & lab[-1L, , drop = FALSE] > 0)
    b[-nr, ] <- b[-nr, ] | (diff_v & lab[-nr, , drop = FALSE] > 0)
  }
  b & lab > 0
}

polygon_area <- function(x, y) {
  if (length(x) < 3L) return(0)
  idx <- grDevices::chull(x, y)
  hx <- x[idx]
  hy <- y[idx]
  hx2 <- c(hx, hx[1])
  hy2 <- c(hy, hy[1])
  abs(sum(hx2[-1] * hy2[-length(hy2)] - hx2[-length(hx2)] * hy2[-1]) / 2)
}

read_gray_png <- function(path) {
  img <- png::readPNG(path)
  if (length(dim(img)) == 3L) {
    if (dim(img)[3] >= 3L) {
      img <- (img[, , 1] + img[, , 2] + img[, , 3]) / 3
    } else {
      img <- img[, , 1]
    }
  }
  img <- as.matrix(img)
  img <- normalize01(img)
  img
}

read_mask_png <- function(path, foreground = "auto") {
  x <- png::readPNG(path)
  if (length(dim(x)) == 3L) {
    if (dim(x)[3] >= 3L) {
      x <- (x[, , 1] + x[, , 2] + x[, , 3]) / 3
    } else {
      x <- x[, , 1]
    }
  }
  x <- as.matrix(x)
  x <- normalize01(x)
  
  foreground <- tolower(foreground)
  if (foreground == "bright") {
    return(x > 0.5)
  }
  if (foreground == "dark") {
    return(x <= 0.5)
  }
  
  bright <- x > 0.5
  dark <- x <= 0.5
  bright_frac <- mean(bright)
  dark_frac <- mean(dark)
  
  # Most masks are bright foreground on dark background. If that is nearly all
  # the image, flip to dark foreground. Otherwise keep bright foreground.
  if (bright_frac > 0.90 && dark_frac > 0.001) {
    dark
  } else {
    bright
  }
}

pick_col <- function(nms, exact = character(), regex = character()) {
  lower <- tolower(nms)
  for (x in exact) {
    hit <- which(lower == tolower(x))
    if (length(hit) > 0L) return(nms[hit[1]])
  }
  for (rx in regex) {
    hit <- grep(rx, lower, perl = TRUE)
    if (length(hit) > 0L) return(nms[hit[1]])
  }
  NA_character_
}

find_first_file <- function(folder, explicit_name = NA, patterns = character(), exclude_regex = NULL) {
  if (!dir.exists(folder)) return(NA_character_)
  
  if (!is_na_arg(explicit_name)) {
    p <- file.path(folder, explicit_name)
    if (file.exists(p)) return(normalizePath(p, mustWork = TRUE))
    p2 <- list.files(folder, pattern = paste0("^", gsub("([.])", "\\\\\\1", basename(explicit_name)), "$"),
                     full.names = TRUE, ignore.case = TRUE)
    if (length(p2) > 0L) return(normalizePath(p2[1], mustWork = TRUE))
    return(NA_character_)
  }
  
  files <- list.files(folder, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0L) return(NA_character_)
  bn <- basename(files)
  
  if (!is.null(exclude_regex)) {
    keep <- !grepl(exclude_regex, bn, ignore.case = TRUE, perl = TRUE)
    files <- files[keep]
    bn <- bn[keep]
  }
  
  for (pat in patterns) {
    hit <- files[grepl(pat, bn, ignore.case = TRUE, perl = TRUE)]
    if (length(hit) > 0L) return(normalizePath(hit[1], mustWork = TRUE))
  }
  NA_character_
}

looks_like_path <- function(x) {
  !is_na_arg(x) && grepl("[/\\\\]", as.character(x))
}

resolve_external_file <- function(folder,
                                  image_id,
                                  explicit_path = NA,
                                  base_dir = NA,
                                  file_name = NA,
                                  patterns = character(),
                                  exclude_regex = NULL) {
  # Priority 1: explicit one-off path, intended for single-image debugging.
  if (!is_na_arg(explicit_path)) {
    if (file.exists(explicit_path)) return(normalizePath(explicit_path, mustWork = TRUE))
    return(NA_character_)
  }
  
  # Priority 2: if the "name" argument itself contains a directory, treat it as a path.
  # This lets --rosette_table_name results/step2/Snap_246/all_candidate_centers.tsv work.
  if (looks_like_path(file_name) || (!is_na_arg(file_name) && file.exists(file_name))) {
    if (file.exists(file_name)) return(normalizePath(file_name, mustWork = TRUE))
    return(NA_character_)
  }
  
  # Priority 3: external base directory, useful when Step 1 and Step 2 outputs
  # are in different roots. In batch mode, the per-image subfolder is tried first.
  if (!is_na_arg(base_dir)) {
    per_image_dir <- file.path(base_dir, image_id)
    
    if (!is_na_arg(file_name)) {
      p1 <- file.path(per_image_dir, file_name)
      if (file.exists(p1)) return(normalizePath(p1, mustWork = TRUE))
      
      p2 <- file.path(base_dir, file_name)
      if (file.exists(p2)) return(normalizePath(p2, mustWork = TRUE))
    }
    
    p3 <- find_first_file(
      per_image_dir,
      explicit_name = file_name,
      patterns = patterns,
      exclude_regex = exclude_regex
    )
    if (!is.na(p3) && file.exists(p3)) return(normalizePath(p3, mustWork = TRUE))
    
    p4 <- find_first_file(
      base_dir,
      explicit_name = file_name,
      patterns = patterns,
      exclude_regex = exclude_regex
    )
    if (!is.na(p4) && file.exists(p4)) return(normalizePath(p4, mustWork = TRUE))
  }
  
  # Priority 4: original behavior: file lives inside the image result folder.
  if (!is_na_arg(file_name)) {
    p5 <- file.path(folder, file_name)
    if (file.exists(p5)) return(normalizePath(p5, mustWork = TRUE))
  }
  
  p6 <- find_first_file(
    folder,
    explicit_name = file_name,
    patterns = patterns,
    exclude_regex = exclude_regex
  )
  if (!is.na(p6) && file.exists(p6)) return(normalizePath(p6, mustWork = TRUE))
  
  NA_character_
}

find_result_folders <- function(input_dir, args) {
  # Simplified, strict behavior:
  #   --input may be either:
  #       1) one sample folder:
  #          results/<experiment>/<samplefolder>/
  #       2) one experiment folder whose direct children are sample folders:
  #          results/<experiment>/<samplefolder>/
  #
  # Each sample folder must directly contain the Step 1 and Step 2 outputs:
  #       gray.png
  #       mask.binary.png
  #       all_candidate_centers.tsv
  #
  # No recursive discovery, no external summary-table locations, and no alternate
  # output roots are interpreted here. The *_path and *_dir arguments are still
  # accepted for CLI compatibility but are intentionally ignored by this lookup.
  
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  
  mask_patterns <- c(
    "^cluster.*mask.*\\.png$",
    "^object.*mask.*\\.png$",
    "^outline.*mask.*\\.png$",
    "^filled.*mask.*\\.png$",
    "^mask.*\\.png$",
    "mask.*\\.png$",
    "binary.*\\.png$"
  )
  
  rosette_patterns <- c(
    "rosette.*center.*\\.tsv$",
    "rosette.*radius.*\\.tsv$",
    "rosette.*\\.tsv$",
    "candidate.*center.*\\.tsv$",
    "all_candidate_centers.*\\.tsv$",
    "center.*\\.tsv$",
    "radii.*\\.tsv$"
  )
  
  inspect_sample_dir <- function(sample_dir) {
    if (!dir.exists(sample_dir)) {
      return(list(
        sample_dir = sample_dir,
        is_dir = FALSE,
        gray_path = NA_character_,
        gray_ok = FALSE,
        mask_path = NA_character_,
        mask_ok = FALSE,
        rosette_path = NA_character_,
        rosette_ok = FALSE,
        valid = FALSE,
        local_files = character()
      ))
    }
    
    local_files <- list.files(sample_dir, full.names = FALSE, recursive = FALSE)
    local_files <- local_files[file.info(file.path(sample_dir, local_files))$isdir %in% FALSE]
    
    gray_path <- file.path(sample_dir, args$gray_name)
    gray_ok <- file.exists(gray_path)
    
    mask_path <- find_first_file(
      sample_dir,
      explicit_name = args$cluster_mask_name,
      patterns = mask_patterns,
      exclude_regex = "overlay|debug|cell|assignment|boundary|membrane|distance|thumb"
    )
    mask_ok <- !is.na(mask_path) && file.exists(mask_path)
    
    rosette_path <- find_first_file(
      sample_dir,
      explicit_name = args$rosette_table_name,
      patterns = rosette_patterns,
      exclude_regex = "cell_objects|cell_counts|geometry"
    )
    rosette_ok <- !is.na(rosette_path) && file.exists(rosette_path)
    
    list(
      sample_dir = sample_dir,
      is_dir = TRUE,
      gray_path = if (gray_ok) normalizePath(gray_path, mustWork = TRUE) else gray_path,
      gray_ok = gray_ok,
      mask_path = mask_path,
      mask_ok = mask_ok,
      rosette_path = rosette_path,
      rosette_ok = rosette_ok,
      valid = gray_ok && mask_ok && rosette_ok,
      local_files = local_files
    )
  }
  
  format_found <- function(ok, path) {
    if (isTRUE(ok)) {
      paste0("FOUND: ", path)
    } else {
      paste0("MISSING: ", path)
    }
  }
  
  print_inspection <- function(rec, max_files = 20L) {
    messagef("[03b_geometry] inspect: %s", rec$sample_dir)
    messagef("[03b_geometry]   gray    %s", format_found(rec$gray_ok, rec$gray_path))
    messagef("[03b_geometry]   mask    %s", format_found(rec$mask_ok, rec$mask_path))
    messagef("[03b_geometry]   centers %s", format_found(rec$rosette_ok, rec$rosette_path))
    
    if (length(rec$local_files) == 0L) {
      messagef("[03b_geometry]   local files: <none>")
    } else {
      shown <- utils::head(rec$local_files, max_files)
      suffix <- if (length(rec$local_files) > max_files) {
        sprintf(" ... (+%d more)", length(rec$local_files) - max_files)
      } else {
        ""
      }
      messagef("[03b_geometry]   local files: %s%s", paste(shown, collapse = ", "), suffix)
    }
  }
  
  ignored_args <- c(
    "gray_path", "gray_dir",
    "cluster_mask_path", "cluster_mask_dir",
    "rosette_table_path", "rosette_table_dir"
  )
  used_ignored <- ignored_args[vapply(ignored_args, function(x) !is_na_arg(args[[x]]), logical(1))]
  if (length(used_ignored) > 0L) {
    messagef(
      "[03b_geometry] NOTE: ignoring external path/root argument(s) in sample-folder-only mode: %s",
      paste(paste0("--", used_ignored), collapse = ", ")
    )
  }
  
  # Case 1: input is itself a sample folder.
  input_rec <- inspect_sample_dir(input_dir)
  if (isTRUE(input_rec$valid)) {
    messagef("[03b_geometry] Treating --input as one sample folder.")
    return(input_dir)
  }
  
  # Case 2: input is an experiment folder; process direct child sample folders.
  direct_children <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  direct_children <- direct_children[dir.exists(direct_children)]
  direct_children <- sort(unique(normalizePath(direct_children, mustWork = TRUE)))
  
  if (length(direct_children) == 0L) {
    messagef("[03b_geometry] No direct child folders found under: %s", input_dir)
    messagef("[03b_geometry] Since --input was not itself a valid sample folder, discovery failed.")
    print_inspection(input_rec)
    return(character())
  }
  
  inspections <- lapply(direct_children, inspect_sample_dir)
  valid <- vapply(inspections, function(x) isTRUE(x$valid), logical(1))
  valid_dirs <- direct_children[valid]
  
  if (length(valid_dirs) == 0L) {
    messagef("[03b_geometry] No valid sample folders found in: %s", input_dir)
    messagef("[03b_geometry] Checked %d direct child folder(s). Showing first %d:", length(inspections), min(10L, length(inspections)))
    
    for (rec in utils::head(inspections, 10L)) {
      print_inspection(rec)
    }
    
    messagef("[03b_geometry] Expected either:")
    messagef("[03b_geometry]   --input results/<experiment>/<samplefolder>")
    messagef("[03b_geometry] or:")
    messagef("[03b_geometry]   --input results/<experiment>")
    messagef("[03b_geometry] with direct child folders containing these exact/named files:")
    messagef("[03b_geometry]   gray image: %s", args$gray_name)
    messagef("[03b_geometry]   mask:       %s", ifelse(is_na_arg(args$cluster_mask_name), "<mask file matching mask/binary patterns>", args$cluster_mask_name))
    messagef("[03b_geometry]   centers:    %s", ifelse(is_na_arg(args$rosette_table_name), "<center table matching candidate/center patterns>", args$rosette_table_name))
    messagef("[03b_geometry] IMPORTANT: files must be directly inside each sample folder, not nested one level deeper.")
  } else {
    messagef("[03b_geometry] Valid sample folders found: %d/%d", length(valid_dirs), length(inspections))
    invalid_n <- length(inspections) - length(valid_dirs)
    if (invalid_n > 0L) {
      messagef("[03b_geometry] Skipping %d direct child folder(s) that did not contain the required files.", invalid_n)
    }
  }
  
  valid_dirs
}

ensure_same_dims <- function(a, b, label_a = "gray", label_b = "mask") {
  if (!identical(dim(a), dim(b))) {
    stop(sprintf("Dimension mismatch between %s %s and %s %s",
                 label_a, paste(dim(a), collapse = "x"),
                 label_b, paste(dim(b), collapse = "x")))
  }
}

# ----------------------------- rosette table ---------------------------------

standardize_rosettes <- function(path, h, w, cluster_lab = NULL, args = NULL, image_id = NULL) {
  parser_debug <- !is.null(args) && isTRUE(args$debug_parser)
  
  dbg <- function(...) {
    if (parser_debug) messagef(...)
  }
  
  dbg("[03b_geometry]   parser: reading candidate table with base read.delim")
  
  # Read with base R, not fread/data.table. This is intentionally conservative:
  # it avoids data.table's row-subsetting methods inside the parser, which is
  # where the previous "primitive next method: subscript out of bounds" was
  # being triggered.
  df <- tryCatch(
    utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      messagef("[03b_geometry]   parser: base read.delim failed: %s", conditionMessage(e))
      messagef("[03b_geometry]   parser: falling back to fread(data.table=FALSE)")
      data.table::fread(path, data.table = FALSE)
    }
  )
  
  # Force plain list-backed data.frame with unique object identity.
  df <- as.data.frame(as.list(df), stringsAsFactors = FALSE, check.names = FALSE)
  
  dbg("[03b_geometry]   parser: class=%s", paste(class(df), collapse = ","))
  dbg("[03b_geometry]   parser: dimensions=%d rows x %d columns", nrow(df), ncol(df))
  dbg("[03b_geometry]   parser: columns=%s", paste(names(df), collapse = ", "))
  
  if (nrow(df) == 0L) {
    return(data.table())
  }
  
  get_df_col <- function(d, col) {
    j <- match(col, names(d))
    if (is.na(j)) {
      stop(sprintf("Internal parser error: column '%s' not found. Available columns: %s",
                   col, paste(names(d), collapse = ", ")))
    }
    unclass(d)[[j]]
  }
  
  filter_df_rows <- function(d, keep, label) {
    keep <- as.logical(keep)
    keep[is.na(keep)] <- FALSE
    idx <- which(keep)
    dbg("[03b_geometry]   parser: %s row indices length=%d; first few=%s",
        label, length(idx), paste(utils::head(idx, 10), collapse = ","))
    
    if (length(idx) == 0L) {
      # Rebuild as empty data.frame without using d[idx, ].
      z <- lapply(unclass(d), function(x) x[FALSE])
      z <- as.data.frame(z, stringsAsFactors = FALSE, check.names = FALSE)
      names(z) <- names(d)
      return(z)
    }
    
    # Rebuild row-filtered data.frame column-by-column. This avoids calling
    # [.data.frame or [.data.table at all.
    z <- lapply(unclass(d), function(x) x[idx])
    z <- as.data.frame(z, stringsAsFactors = FALSE, check.names = FALSE)
    names(z) <- names(d)
    z
  }
  
  nms <- names(df)
  lower_nms <- tolower(nms)
  
  get_existing <- function(candidates) {
    for (cc in candidates) {
      hit <- which(lower_nms == tolower(cc))
      if (length(hit) > 0L) return(nms[hit[1]])
    }
    NA_character_
  }
  
  safe_bool <- function(x) {
    out <- tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
    out[is.na(out)] <- FALSE
    out
  }
  
  source <- "weighted"
  if (!is.null(args) && !is.null(args$candidate_center_source)) {
    source <- tolower(args$candidate_center_source)
  }
  if (!source %in% c("weighted", "best", "raw", "auto")) {
    stop("--candidate_center_source must be one of: weighted, best, raw, auto")
  }
  
  col_sets <- list(
    weighted = list(
      x = c("weighted_center_x", "weighted_x", "rosette_weighted_center_x"),
      y = c("weighted_center_y", "weighted_y", "rosette_weighted_center_y"),
      r = c("weighted_fitted_radius_px", "weighted_radius_px", "weighted_radius")
    ),
    best = list(
      x = c("best_center_x", "best_x", "rosette_best_center_x"),
      y = c("best_center_y", "best_y", "rosette_best_center_y"),
      r = c("best_fitted_radius_px", "best_radius_px", "best_radius")
    ),
    raw = list(
      x = c("center_x", "rosette_center_x", "cx", "x_center", "x", "center_col", "col", "column", "x_px", "center_x_px", "candidate_x"),
      y = c("center_y", "rosette_center_y", "cy", "y_center", "y", "center_row", "row", "y_px", "center_y_px", "candidate_y"),
      r = c("fitted_radius_px", "rosette_radius_px", "radius_px", "radius", "r", "r_px", "candidate_radius", "center_radius")
    )
  )
  
  try_order <- if (source == "auto") c("weighted", "best", "raw") else c(source, setdiff(c("weighted", "best", "raw"), source))
  
  xcol <- ycol <- rcol <- NA_character_
  chosen_source <- NA_character_
  for (src in try_order) {
    tx <- get_existing(col_sets[[src]]$x)
    ty <- get_existing(col_sets[[src]]$y)
    tr <- get_existing(col_sets[[src]]$r)
    if (!is.na(tx) && !is.na(ty)) {
      xcol <- tx
      ycol <- ty
      rcol <- tr
      chosen_source <- src
      break
    }
  }
  
  if (is.na(xcol) || is.na(ycol)) {
    xcol <- pick_col(
      nms,
      exact = c("rosette_center_x", "center_x", "cx", "x_center", "x",
                "center_col", "col", "column", "x_px", "center_x_px", "candidate_x"),
      regex = c("(^|_)center.*x$", "(^|_)center.*col$", "rosette.*x", "centroid.*x", "candidate.*x")
    )
    ycol <- pick_col(
      nms,
      exact = c("rosette_center_y", "center_y", "cy", "y_center", "y",
                "center_row", "row", "y_px", "center_y_px", "candidate_y"),
      regex = c("(^|_)center.*y$", "(^|_)center.*row$", "rosette.*y", "centroid.*y", "candidate.*y")
    )
    rcol <- pick_col(
      nms,
      exact = c("rosette_radius_px", "radius_px", "radius", "r", "r_px", "candidate_radius", "center_radius", "fitted_radius_px"),
      regex = c("radius", "(^|_)r_px$")
    )
    chosen_source <- "generic"
  }
  
  idcol <- pick_col(
    nms,
    exact = c("rosette_candidate_id", "rosette_id", "center_id", "id"),
    regex = c("rosette.*id", "center.*id", "candidate.*id")
  )
  clcol <- pick_col(
    nms,
    exact = c("parent_blob_id", "cluster_id", "object_id", "clump_id", "blob_id"),
    regex = c("parent.*blob.*id", "cluster.*id", "object.*id", "clump.*id", "blob.*id")
  )
  accepted_col <- pick_col(nms, exact = c("accepted"), regex = c("^accepted$"))
  confidence_col <- pick_col(nms, exact = c("confidence"), regex = c("^confidence$"))
  
  inside_col <- NA_character_
  if (!is.null(args) && isTRUE(args$require_center_inside_blob)) {
    if (chosen_source == "weighted") {
      inside_col <- get_existing(c("weighted_center_inside_blob", "center_inside_blob"))
    } else {
      inside_col <- get_existing(c("center_inside_blob", "weighted_center_inside_blob"))
    }
  }
  
  if (is.na(xcol) || is.na(ycol)) {
    stop(sprintf("Could not identify rosette x/y columns in %s. Columns: %s",
                 path, paste(nms, collapse = ", ")))
  }
  
  messagef("[03b_geometry]   centers table source: %s", chosen_source)
  messagef("[03b_geometry]   centers x/y/r columns: %s / %s / %s",
           xcol, ycol, ifelse(is.na(rcol), "fallback_radius", rcol))
  dbg("[03b_geometry]   parser: id column=%s; cluster/blob column=%s; accepted column=%s; confidence column=%s",
      ifelse(is.na(idcol), "none", idcol),
      ifelse(is.na(clcol), "none", clcol),
      ifelse(is.na(accepted_col), "none", accepted_col),
      ifelse(is.na(confidence_col), "none", confidence_col))
  
  if (!is.null(image_id) && "image_id" %in% names(df) && nrow(df) > 0L) {
    img_col <- as.character(get_df_col(df, "image_id"))
    m <- img_col == as.character(image_id)
    m[is.na(m)] <- FALSE
    if (any(m)) {
      messagef("[03b_geometry]   image_id filter: keeping %d/%d center rows for %s", sum(m), nrow(df), image_id)
      df <- filter_df_rows(df, m, "image_id filter")
    } else {
      dbg("[03b_geometry]   parser: no image_id rows exactly matched '%s'; not filtering by image_id", image_id)
      dbg("[03b_geometry]   parser: first image_id values: %s", paste(utils::head(unique(img_col), 10), collapse = ", "))
    }
  }
  
  if (!is.null(args) && isTRUE(args$accepted_only) && !is.na(accepted_col) && nrow(df) > 0L) {
    keep <- safe_bool(get_df_col(df, accepted_col))
    messagef("[03b_geometry]   accepted filter: keeping %d/%d center rows", sum(keep), nrow(df))
    df <- filter_df_rows(df, keep, "accepted filter")
    dbg("[03b_geometry]   parser: after accepted filter: %d rows x %d columns", nrow(df), ncol(df))
  }
  
  if (!is.null(args) && isTRUE(args$require_center_inside_blob) && !is.na(inside_col) && nrow(df) > 0L) {
    keep <- safe_bool(get_df_col(df, inside_col))
    messagef("[03b_geometry]   inside-blob filter: keeping %d/%d center rows", sum(keep), nrow(df))
    df <- filter_df_rows(df, keep, "inside-blob filter")
  }
  
  if (!is.null(args) && !is.na(args$min_candidate_confidence) && !is.na(confidence_col) && nrow(df) > 0L) {
    conf <- suppressWarnings(as.numeric(get_df_col(df, confidence_col)))
    keep <- is.finite(conf) & conf >= args$min_candidate_confidence
    keep[is.na(keep)] <- FALSE
    messagef("[03b_geometry]   confidence filter: keeping %d/%d center rows", sum(keep), nrow(df))
    df <- filter_df_rows(df, keep, "confidence filter")
  }
  
  if (nrow(df) == 0L) {
    stop(sprintf("No candidate center rows remain after filtering in %s", path))
  }
  
  dbg("[03b_geometry]   parser: extracting numeric x/y/r columns")
  x_raw <- suppressWarnings(as.numeric(get_df_col(df, xcol)))
  y_raw <- suppressWarnings(as.numeric(get_df_col(df, ycol)))
  r_raw <- if (!is.na(rcol) && rcol %in% names(df)) {
    suppressWarnings(as.numeric(get_df_col(df, rcol)))
  } else {
    rep(NA_real_, nrow(df))
  }
  
  dbg("[03b_geometry]   parser: x range raw=%s", paste(range(x_raw, na.rm = TRUE), collapse = " to "))
  dbg("[03b_geometry]   parser: y range raw=%s", paste(range(y_raw, na.rm = TRUE), collapse = " to "))
  dbg("[03b_geometry]   parser: image bounds h=%d, w=%d", h, w)
  
  id_raw <- if (!is.na(idcol) && idcol %in% names(df)) {
    as.character(get_df_col(df, idcol))
  } else {
    as.character(seq_len(nrow(df)))
  }
  
  cl_raw <- if (!is.na(clcol) && clcol %in% names(df)) {
    as.character(get_df_col(df, clcol))
  } else {
    rep(NA_character_, nrow(df))
  }
  
  bad_xy <- !is.finite(x_raw) | !is.finite(y_raw)
  if (any(bad_xy)) {
    warning(sprintf("Dropping %d candidate center row(s) with non-finite x/y coordinates from %s", sum(bad_xy), path))
  }
  
  keep_xy <- which(!bad_xy)
  dbg("[03b_geometry]   parser: finite coordinate rows=%d/%d", length(keep_xy), length(x_raw))
  if (length(keep_xy) == 0L) {
    stop(sprintf("No usable candidate center rows found in %s", path))
  }
  
  x_raw <- x_raw[keep_xy]
  y_raw <- y_raw[keep_xy]
  r_raw <- r_raw[keep_xy]
  id_raw <- id_raw[keep_xy]
  cl_raw <- cl_raw[keep_xy]
  
  x_clamped <- pmin(pmax(x_raw, 1), w)
  y_clamped <- pmin(pmax(y_raw, 1), h)
  
  if (all(is.na(r_raw))) {
    r_raw <- rep(25, length(x_raw))
  }
  
  good_r <- is.finite(r_raw) & r_raw > 0
  med_r <- if (any(good_r)) stats::median(r_raw[good_r], na.rm = TRUE) else 25
  r_raw[!is.finite(r_raw) | r_raw <= 0] <- med_r
  r_raw[!is.finite(r_raw) | r_raw <= 0] <- 25
  
  vals <- integer(length(x_clamped))
  if (!is.null(cluster_lab) && length(x_clamped) > 0L) {
    # Force the label image into a plain base integer matrix and use linear
    # indexing. This avoids S3/S4 [ dispatch problems seen on Windows with
    # cbind-style matrix indexing.
    cluster_lab_plain <- matrix(
      as.integer(as.vector(cluster_lab)),
      nrow = as.integer(h),
      ncol = as.integer(w)
    )
    cluster_lab_vec <- as.integer(as.vector(cluster_lab_plain))
    
    yy <- as.integer(round(y_clamped))
    xx <- as.integer(round(x_clamped))
    yy <- pmin(pmax(yy, 1L), as.integer(h))
    xx <- pmin(pmax(xx, 1L), as.integer(w))
    ok <- is.finite(yy) & is.finite(xx) & yy >= 1L & yy <= h & xx >= 1L & xx <= w
    dbg("[03b_geometry]   parser: cluster index ok rows=%d/%d", sum(ok), length(ok))
    if (any(ok)) {
      lin <- yy + (xx - 1L) * as.integer(h)
      vals[ok] <- cluster_lab_vec[lin[ok]]
    }
    vals[is.na(vals)] <- 0L
  }
  
  dbg("[03b_geometry]   parser: building standardized rosette table")
  out <- data.table::data.table(
    rosette_row = seq_along(x_clamped),
    rosette_id = id_raw,
    rosette_center_x_raw = as.numeric(x_raw),
    rosette_center_y_raw = as.numeric(y_raw),
    rosette_center_x = as.numeric(x_clamped),
    rosette_center_y = as.numeric(y_clamped),
    rosette_radius_px = as.numeric(r_raw),
    source_cluster_id = cl_raw,
    candidate_center_source = chosen_source,
    mask_cluster_id = as.integer(vals)
  )
  
  dbg("[03b_geometry]   parser: standardized center rows=%d", nrow(out))
  out
}

# ----------------------------- segmentation ----------------------------------

compute_membrane_mask <- function(gray, cluster_mask, args) {
  nr <- nrow(gray)
  nc <- ncol(gray)
  
  gray <- plain_matrix(gray, nr = nr, nc = nc, label = "gray in compute_membrane_mask")
  cluster_mask <- plain_logical_matrix(cluster_mask, nr = nr, nc = nc, label = "cluster mask in compute_membrane_mask")
  
  size <- odd_size_from_radius(args$membrane_bg_radius_px)
  
  bg <- EBImage::medianFilter(EBImage::Image(gray), size = size)
  bg <- plain_matrix(bg, nr = nr, nc = nc, label = "medianFilter output")
  
  # Use base matrices from here forward. This avoids EBImage/S4 logical dispatch
  # errors such as: e1@.Data & e2@.Data : non-conformable arrays.
  score <- matrix(pmax(0, as.numeric(bg) - as.numeric(gray)), nrow = nr, ncol = nc)
  score[!cluster_mask] <- 0
  score <- normalize01(score)
  score <- plain_matrix(score, nr = nr, nc = nc, label = "membrane score")
  
  inside <- score[cluster_mask]
  if (length(inside) == 0L || !any(is.finite(inside))) {
    return(list(
      mask = matrix(FALSE, nrow = nr, ncol = nc),
      score = score,
      threshold = NA_real_
    ))
  }
  
  qthr <- as.numeric(stats::quantile(inside, probs = args$membrane_quantile, na.rm = TRUE, names = FALSE))
  thr <- max(qthr, args$membrane_min_score, na.rm = TRUE)
  
  mem <- matrix((score >= thr) & cluster_mask, nrow = nr, ncol = nc)
  
  if (args$membrane_open_px > 0) {
    mem <- m_opening(mem, args$membrane_open_px)
  }
  if (args$membrane_dilate_px > 0) {
    mem <- m_dilate(mem, args$membrane_dilate_px)
  }
  
  mem <- matrix(mem & cluster_mask, nrow = nr, ncol = nc)
  
  list(mask = mem, score = score, threshold = thr)
}

segment_one_cluster <- function(gray, cluster_mask, args) {
  nr <- nrow(gray)
  nc <- ncol(gray)
  
  gray <- plain_matrix(gray, nr = nr, nc = nc, label = "gray in segment_one_cluster")
  cluster_mask <- plain_logical_matrix(cluster_mask, nr = nr, nc = nc, label = "cluster mask in segment_one_cluster")
  
  mem <- compute_membrane_mask(gray, cluster_mask, args)
  
  work_mask <- cluster_mask
  if (args$use_membrane_cuts) {
    work_mask <- matrix(cluster_mask & !plain_logical_matrix(mem$mask, nr = nr, nc = nc, label = "membrane mask"), nrow = nr, ncol = nc)
  }
  
  # Remove dust introduced by membrane cuts, but keep subcell fragments around
  # for rejection labels/debugging.
  work_mask <- remove_small_components(work_mask, min_area_px = max(3, floor(args$min_cell_area_px / 4)))
  work_mask <- plain_logical_matrix(work_mask, nr = nr, nc = nc, label = "work mask")
  
  if (!any(work_mask)) {
    return(list(
      label = matrix(0L, nrow = nr, ncol = nc),
      membrane_mask = plain_logical_matrix(mem$mask, nr = nr, nc = nc, label = "empty membrane mask"),
      membrane_score = plain_matrix(mem$score, nr = nr, nc = nc, label = "empty membrane score"),
      membrane_threshold = mem$threshold,
      distance = matrix(0, nrow = nr, ncol = nc)
    ))
  }
  
  dist <- plain_matrix(EBImage::distmap(EBImage::Image(work_mask)), nr = nr, nc = nc, label = "distmap output")
  dist[!work_mask] <- 0
  dist_s <- dist
  if (args$seed_smooth_sigma > 0) {
    dist_s <- plain_matrix(
      EBImage::gblur(EBImage::Image(dist), sigma = args$seed_smooth_sigma),
      nr = nr,
      nc = nc,
      label = "gblur distance output"
    )
    dist_s[!work_mask] <- 0
  }
  
  ws <- tryCatch(
    plain_matrix(
      EBImage::watershed(
        EBImage::Image(dist_s),
        tolerance = args$watershed_tolerance,
        ext = args$seed_min_distance_px
      ),
      nr = nr,
      nc = nc,
      label = "watershed output"
    ),
    error = function(e) {
      warning(sprintf("Watershed failed; falling back to connected components: %s", conditionMessage(e)))
      bwlabel_matrix(work_mask)
    }
  )
  
  ws[!work_mask] <- 0
  ws <- relabel_sequential(ws)
  
  if (max(ws, na.rm = TRUE) == 0L) {
    ws <- bwlabel_matrix(work_mask)
  }
  
  list(
    label = ws,
    membrane_mask = plain_logical_matrix(mem$mask, nr = nr, nc = nc, label = "membrane mask return"),
    membrane_score = plain_matrix(mem$score, nr = nr, nc = nc, label = "membrane score return"),
    membrane_threshold = mem$threshold,
    distance = plain_matrix(dist_s, nr = nr, nc = nc, label = "distance return")
  )
}

measure_cell <- function(cell_mask, cluster_mask, args) {
  nr <- nrow(cluster_mask)
  nc <- ncol(cluster_mask)
  
  cluster_mask <- plain_logical_matrix(cluster_mask, nr = nr, nc = nc, label = "measure_cell cluster_mask")
  cell_mask <- plain_logical_matrix(cell_mask, nr = nr, nc = nc, label = "measure_cell cell_mask")
  
  if (!any(cell_mask)) {
    return(list(
      area = 0,
      centroid_x = NA_real_,
      centroid_y = NA_real_,
      perimeter = 0,
      circularity = NA_real_,
      solidity = NA_real_,
      max_radius = NA_real_,
      mean_radius = NA_real_,
      edge_contact_fraction = NA_real_,
      bbox_width_px = NA_real_,
      bbox_height_px = NA_real_,
      aspect_ratio = NA_real_,
      completeness_score = 0,
      accepted = FALSE,
      reject_reason = "empty_cell_mask"
    ))
  }
  
  idx <- which(cell_mask, arr.ind = TRUE)
  area <- length(idx[, 1])
  centroid_x <- mean(idx[, 2])
  centroid_y <- mean(idx[, 1])
  
  cell_boundary <- boundary_from_mask(cell_mask)
  cell_boundary <- plain_logical_matrix(cell_boundary, nr = nr, nc = nc, label = "measure_cell cell_boundary")
  perimeter <- sum(cell_boundary)
  
  circularity <- if (perimeter > 0) {
    4 * pi * area / (perimeter^2)
  } else {
    NA_real_
  }
  if (is.finite(circularity)) circularity <- min(circularity, 1)
  
  hull_area <- polygon_area(idx[, 2], idx[, 1])
  hull_area <- max(hull_area, area)
  solidity <- if (hull_area > 0) area / hull_area else NA_real_
  if (is.finite(solidity)) solidity <- min(solidity, 1)
  
  cell_dt <- plain_matrix(
    EBImage::distmap(EBImage::Image(cell_mask)),
    nr = nr,
    nc = nc,
    label = "measure_cell cell distmap"
  )
  max_radius <- if (any(cell_mask)) max(cell_dt[cell_mask], na.rm = TRUE) else NA_real_
  mean_radius <- if (any(cell_mask)) mean(cell_dt[cell_mask], na.rm = TRUE) else NA_real_
  
  cl_boundary <- boundary_from_mask(cluster_mask)
  cl_boundary <- plain_logical_matrix(cl_boundary, nr = nr, nc = nc, label = "measure_cell cluster_boundary")
  cl_boundary_d <- m_dilate(cl_boundary, 1)
  cl_boundary_d <- plain_logical_matrix(cl_boundary_d, nr = nr, nc = nc, label = "measure_cell dilated_cluster_boundary")
  
  edge_contact <- if (perimeter > 0) {
    sum(cell_boundary & cl_boundary_d) / perimeter
  } else {
    NA_real_
  }
  
  bbox_w <- max(idx[, 2]) - min(idx[, 2]) + 1
  bbox_h <- max(idx[, 1]) - min(idx[, 1]) + 1
  aspect_ratio <- max(bbox_w, bbox_h) / max(1, min(bbox_w, bbox_h))
  
  reject <- character()
  
  if (!is.finite(area) || area < args$min_cell_area_px) reject <- c(reject, "area_low")
  if (!is.finite(area) || area > args$max_cell_area_px) reject <- c(reject, "area_high")
  if (!is.finite(max_radius) || max_radius < args$min_cell_radius_px) reject <- c(reject, "radius_low")
  if (!is.finite(max_radius) || max_radius > args$max_cell_radius_px) reject <- c(reject, "radius_high")
  if (!is.finite(solidity) || solidity < args$min_solidity) reject <- c(reject, "solidity_low")
  if (!is.finite(circularity) || circularity < args$min_circularity) reject <- c(reject, "circularity_low")
  if (!is.finite(edge_contact) || edge_contact > args$max_edge_contact_fraction) reject <- c(reject, "edge_contact_high")
  
  area_score <- if (area < args$min_cell_area_px) {
    area / max(1, args$min_cell_area_px)
  } else if (area > args$max_cell_area_px) {
    args$max_cell_area_px / max(1, area)
  } else {
    1
  }
  
  radius_score <- if (!is.finite(max_radius)) {
    0
  } else if (max_radius < args$min_cell_radius_px) {
    max_radius / max(1, args$min_cell_radius_px)
  } else if (max_radius > args$max_cell_radius_px) {
    args$max_cell_radius_px / max(1, max_radius)
  } else {
    1
  }
  
  solidity_score <- if (is.finite(solidity)) min(1, solidity / max(1e-6, args$min_solidity)) else 0
  circ_score <- if (is.finite(circularity)) min(1, circularity / max(1e-6, args$min_circularity)) else 0
  edge_score <- if (is.finite(edge_contact)) max(0, 1 - edge_contact / max(1e-6, args$max_edge_contact_fraction)) else 0
  
  completeness_score <- mean(c(area_score, radius_score, solidity_score, circ_score, edge_score), na.rm = TRUE)
  
  list(
    area = as.numeric(area),
    centroid_x = as.numeric(centroid_x),
    centroid_y = as.numeric(centroid_y),
    perimeter = as.numeric(perimeter),
    circularity = as.numeric(circularity),
    solidity = as.numeric(solidity),
    max_radius = as.numeric(max_radius),
    mean_radius = as.numeric(mean_radius),
    edge_contact_fraction = as.numeric(edge_contact),
    bbox_width_px = as.numeric(bbox_w),
    bbox_height_px = as.numeric(bbox_h),
    aspect_ratio = as.numeric(aspect_ratio),
    completeness_score = as.numeric(completeness_score),
    accepted = length(reject) == 0L,
    reject_reason = if (length(reject) == 0L) "" else paste(unique(reject), collapse = ";")
  )
}

assign_cells_to_rosettes <- function(cells, label_img, rosettes, args) {
  if (nrow(cells) == 0L) return(cells)
  
  label_img <- plain_matrix(label_img, label = "assign label_img")
  storage.mode(label_img) <- "integer"
  
  cells[, `:=`(
    assigned_rosette_id = NA_character_,
    assigned_rosette_row = NA_integer_,
    assignment_score = NA_real_,
    distance_to_rosette_center_px = NA_real_,
    distance_to_rosette_center_norm = NA_real_,
    rosette_overlap_fraction = NA_real_
  )]
  
  if (nrow(rosettes) == 0L) return(cells)
  
  h <- nrow(label_img)
  w <- ncol(label_img)
  
  for (i in seq_len(nrow(cells))) {
    cid <- cells$cell_num[i]
    cell_mask <- matrix(label_img == cid, nrow = h, ncol = w)
    
    candidates <- rosettes
    clid <- cells$mask_cluster_id[i]
    if ("mask_cluster_id" %in% names(rosettes) && is.finite(clid) && clid > 0) {
      tmp <- rosettes[mask_cluster_id == clid]
      if (nrow(tmp) > 0L) candidates <- tmp
    }
    
    if (nrow(candidates) == 0L) next
    
    cx <- cells$cell_centroid_x[i]
    cy <- cells$cell_centroid_y[i]
    if (!is.finite(cx) || !is.finite(cy)) next
    
    best <- NULL
    
    for (j in seq_len(nrow(candidates))) {
      rx <- candidates$rosette_center_x[j]
      ry <- candidates$rosette_center_y[j]
      rr <- max(1, candidates$rosette_radius_px[j])
      
      d <- sqrt((cx - rx)^2 + (cy - ry)^2)
      nd <- d / rr
      
      overlap_frac <- 0
      if (args$assignment_overlap_weight != 0 || args$assignment_min_overlap > 0) {
        idx <- which(cell_mask, arr.ind = TRUE)
        if (nrow(idx) > 0L) {
          dx <- idx[, 2] - rx
          dy <- idx[, 1] - ry
          overlap_frac <- mean((dx * dx + dy * dy) <= rr * rr)
        }
      }
      
      allowed <- is.finite(nd) &&
        (nd <= args$assignment_max_norm_distance || overlap_frac >= args$assignment_min_overlap)
      
      if (!allowed) next
      
      score <- args$assignment_radius_weight * nd -
        args$assignment_overlap_weight * overlap_frac
      
      rec <- list(
        row = candidates$rosette_row[j],
        id = candidates$rosette_id[j],
        score = score,
        d = d,
        nd = nd,
        overlap = overlap_frac
      )
      
      if (is.null(best) || rec$score < best$score) best <- rec
    }
    
    if (!is.null(best)) {
      cells$assigned_rosette_row[i] <- as.integer(best$row)
      cells$assigned_rosette_id[i] <- as.character(best$id)
      cells$assignment_score[i] <- as.numeric(best$score)
      cells$distance_to_rosette_center_px[i] <- as.numeric(best$d)
      cells$distance_to_rosette_center_norm[i] <- as.numeric(best$nd)
      cells$rosette_overlap_fraction[i] <- as.numeric(best$overlap)
    }
  }
  
  cells
}

# ----------------------------- overlay output --------------------------------

blend_mask <- function(rgb, mask, color, alpha = 0.45) {
  mask <- mask > 0
  if (!any(mask)) return(rgb)
  for (k in 1:3) {
    plane <- rgb[, , k]
    plane[mask] <- (1 - alpha) * plane[mask] + alpha * color[k]
    rgb[, , k] <- plane
  }
  rgb
}

set_mask_color <- function(rgb, mask, color) {
  mask <- mask > 0
  if (!any(mask)) return(rgb)
  for (k in 1:3) {
    plane <- rgb[, , k]
    plane[mask] <- color[k]
    rgb[, , k] <- plane
  }
  rgb
}

palette_for_ids <- function(ids) {
  ids <- sort(unique(as.character(ids[!is.na(ids)])))
  if (length(ids) == 0L) return(setNames(character(), character()))
  cols <- grDevices::rainbow(length(ids), s = 0.9, v = 1.0)
  setNames(cols, ids)
}

col_to_rgb01 <- function(col) {
  grDevices::col2rgb(col)[, 1] / 255
}

make_rgb_base <- function(gray) {
  gray <- normalize01(gray)
  array(rep(gray, 3), dim = c(nrow(gray), ncol(gray), 3))
}

write_png01 <- function(mat, path) {
  png::writePNG(normalize01(mat), target = path)
}

write_overlay <- function(gray, cluster_lab, cell_lab, membrane_mask, rosettes, cells, path, args, mode = "full") {
  h <- nrow(gray)
  w <- ncol(gray)
  
  rgb <- make_rgb_base(gray)
  
  cluster_boundary <- boundary_from_label(cluster_lab)
  cell_boundary <- boundary_from_label(cell_lab)
  if (args$boundary_line_width_px > 0) {
    cell_boundary <- m_dilate(cell_boundary, args$boundary_line_width_px)
  }
  
  # rejected fragments first, accepted cells second
  if (nrow(cells) > 0L) {
    if (args$draw_rejected) {
      rej_ids <- cells[accepted_cell == FALSE, cell_num]
      if (length(rej_ids) > 0L) {
        rgb <- blend_mask(rgb, cell_lab %in% rej_ids, c(1.0, 0.0, 0.0), alpha = 0.35)
      }
    }
    
    ros_cols <- palette_for_ids(cells$assigned_rosette_id)
    acc <- cells[accepted_cell == TRUE]
    if (nrow(acc) > 0L) {
      for (i in seq_len(nrow(acc))) {
        cid <- acc$cell_num[i]
        rid <- acc$assigned_rosette_id[i]
        col <- if (!is.na(rid) && rid %in% names(ros_cols)) {
          col_to_rgb01(ros_cols[[rid]])
        } else {
          c(0.65, 0.65, 0.65)
        }
        rgb <- blend_mask(rgb, cell_lab == cid, col, alpha = 0.28)
      }
    }
  }
  
  if (!is.null(membrane_mask) && any(membrane_mask)) {
    # cyan membrane candidates
    rgb <- blend_mask(rgb, membrane_mask, c(0.0, 1.0, 1.0), alpha = 0.45)
  }
  
  # white cluster outline, yellow cell boundaries
  rgb <- set_mask_color(rgb, m_dilate(cluster_boundary, 1), c(1, 1, 1))
  rgb <- set_mask_color(rgb, cell_boundary, c(1, 1, 0))
  
  png(
    filename = path,
    width = max(480, as.integer(w * args$overlay_scale)),
    height = max(480, as.integer(h * args$overlay_scale))
  )
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(1, w), ylim = c(h, 1), asp = 1, axes = FALSE, xlab = "", ylab = "")
  rasterImage(as.raster(rgb), 1, h, w, 1)
  
  if (nrow(rosettes) > 0L) {
    # Draw rosette radii and centers.
    symbols(
      rosettes$rosette_center_x,
      rosettes$rosette_center_y,
      circles = rosettes$rosette_radius_px,
      inches = FALSE,
      add = TRUE,
      fg = "deepskyblue",
      lwd = 2
    )
    points(rosettes$rosette_center_x, rosettes$rosette_center_y, pch = 3, col = "deepskyblue", lwd = 2)
    text(
      rosettes$rosette_center_x,
      rosettes$rosette_center_y,
      labels = paste0("R", rosettes$rosette_id),
      col = "deepskyblue",
      cex = 0.75,
      pos = 3
    )
  }
  
  if (nrow(cells) > 0L) {
    # Assignment lines for accepted cells.
    assigned <- cells[accepted_cell == TRUE & !is.na(assigned_rosette_row)]
    if (nrow(assigned) > 0L && nrow(rosettes) > 0L) {
      tmp <- merge(
        assigned,
        rosettes[, .(rosette_row, rosette_center_x, rosette_center_y)],
        by.x = "assigned_rosette_row",
        by.y = "rosette_row",
        all.x = TRUE,
        sort = FALSE
      )
      if (nrow(tmp) > 0L) {
        segments(
          tmp$cell_centroid_x,
          tmp$cell_centroid_y,
          tmp$rosette_center_x,
          tmp$rosette_center_y,
          col = grDevices::adjustcolor("white", alpha.f = 0.35),
          lwd = 1
        )
      }
    }
    
    # Cell centroids.
    acc <- cells[accepted_cell == TRUE]
    rej <- cells[accepted_cell == FALSE]
    if (nrow(acc) > 0L) {
      points(acc$cell_centroid_x, acc$cell_centroid_y, pch = 16, col = "limegreen", cex = 0.55)
    }
    if (args$draw_rejected && nrow(rej) > 0L) {
      points(rej$cell_centroid_x, rej$cell_centroid_y, pch = 4, col = "red", cex = 0.55)
    }
    
    if (args$draw_cell_ids && tolower(args$label_mode) != "none") {
      lab_dt <- if (args$draw_rejected) cells else cells[accepted_cell == TRUE]
      if (nrow(lab_dt) > 0L) {
        labels <- switch(
          tolower(args$label_mode),
          cell = paste0("C", lab_dt$cell_num),
          rosette = ifelse(is.na(lab_dt$assigned_rosette_id), "R?", paste0("R", lab_dt$assigned_rosette_id)),
          both = paste0("C", lab_dt$cell_num, ifelse(is.na(lab_dt$assigned_rosette_id), "/R?", paste0("/R", lab_dt$assigned_rosette_id))),
          paste0("C", lab_dt$cell_num, ifelse(is.na(lab_dt$assigned_rosette_id), "/R?", paste0("/R", lab_dt$assigned_rosette_id)))
        )
        text(
          lab_dt$cell_centroid_x,
          lab_dt$cell_centroid_y,
          labels = labels,
          col = ifelse(lab_dt$accepted_cell, "white", "red"),
          cex = 0.55
        )
      }
    }
  }
  
  legend(
    "bottomleft",
    legend = c("cluster outline", "cell boundary", "membrane candidate", "accepted centroid", "rejected fragment", "rosette radius"),
    col = c("white", "yellow", "cyan", "limegreen", "red", "deepskyblue"),
    pch = c(NA, NA, 15, 16, 4, NA),
    lty = c(1, 1, NA, NA, NA, 1),
    lwd = c(2, 2, NA, NA, NA, 2),
    bty = "n",
    text.col = "white",
    cex = 0.70
  )
  
  invisible(TRUE)
}

# ----------------------------- summaries -------------------------------------

summarize_rosettes <- function(rosettes, cells, image_id) {
  if (nrow(rosettes) == 0L) {
    return(data.table())
  }
  
  base <- copy(rosettes)
  base[, image_id := image_id]
  
  accepted_counts <- cells[
    accepted_cell == TRUE & !is.na(assigned_rosette_row),
    .(
      n_cells_geometry = as.integer(.N),
      mean_cell_area_px = as.numeric(mean(cell_area_px, na.rm = TRUE)),
      median_cell_area_px = as.numeric(stats::median(cell_area_px, na.rm = TRUE)),
      mean_cell_circularity = as.numeric(mean(cell_circularity, na.rm = TRUE)),
      mean_cell_solidity = as.numeric(mean(cell_solidity, na.rm = TRUE))
    ),
    by = .(assigned_rosette_row)
  ]
  
  rejected_counts <- cells[
    accepted_cell == FALSE & !is.na(assigned_rosette_row),
    .(n_rejected_fragments_nearby = as.integer(.N)),
    by = .(assigned_rosette_row)
  ]
  
  out <- merge(
    base,
    accepted_counts,
    by.x = "rosette_row",
    by.y = "assigned_rosette_row",
    all.x = TRUE,
    sort = FALSE
  )
  out <- merge(
    out,
    rejected_counts,
    by.x = "rosette_row",
    by.y = "assigned_rosette_row",
    all.x = TRUE,
    sort = FALSE
  )
  
  int_cols <- c("n_cells_geometry", "n_rejected_fragments_nearby")
  for (cc in int_cols) {
    if (!cc %in% names(out)) out[, (cc) := 0L]
    idx <- which(is.na(out[[cc]]))
    if (length(idx) > 0L) data.table::set(out, i = idx, j = cc, value = 0L)
    out[, (cc) := as.integer(get(cc))]
  }
  
  num_cols <- c("mean_cell_area_px", "median_cell_area_px", "mean_cell_circularity", "mean_cell_solidity")
  for (cc in num_cols) {
    if (!cc %in% names(out)) out[, (cc) := NA_real_]
    out[, (cc) := as.numeric(get(cc))]
  }
  
  preferred_cols <- c(
    "image_id",
    "rosette_row",
    "rosette_id",
    "mask_cluster_id",
    "source_cluster_id",
    "candidate_center_source",
    "rosette_center_x",
    "rosette_center_y",
    "rosette_radius_px",
    "n_cells_geometry",
    "n_rejected_fragments_nearby",
    "mean_cell_area_px",
    "median_cell_area_px",
    "mean_cell_circularity",
    "mean_cell_solidity"
  )
  data.table::setcolorder(out, c(intersect(preferred_cols, names(out)), setdiff(names(out), preferred_cols)))
  
  out
}

# ----------------------------- main processing -------------------------------

process_folder <- function(folder, args) {
  folder_id <- basename(normalizePath(folder, mustWork = TRUE))
  image_id <- folder_id
  if (!is_na_arg(args$image_id_override)) {
    image_id <- as.character(args$image_id_override)
  }
  image_safe <- safe_file_stem(image_id)
  
  messagef("[03b_geometry] Processing: %s", image_id)
  if (!identical(image_id, folder_id)) {
    messagef("[03b_geometry]   folder id: %s", folder_id)
  }
  
  mask_patterns <- c(
    "^cluster.*mask.*\\.png$",
    "^object.*mask.*\\.png$",
    "^outline.*mask.*\\.png$",
    "^filled.*mask.*\\.png$",
    "^mask.*\\.png$",
    "mask.*\\.png$",
    "binary.*\\.png$"
  )
  
  rosette_patterns <- c(
    "rosette.*center.*\\.tsv$",
    "rosette.*radius.*\\.tsv$",
    "rosette.*\\.tsv$",
    "candidate.*center.*\\.tsv$",
    "all_candidate_centers.*\\.tsv$",
    "center.*\\.tsv$",
    "radii.*\\.tsv$"
  )
  
  # Sample-folder-only behavior:
  # all required Step 1 and Step 2 files are expected directly in this folder.
  gray_path <- file.path(folder, args$gray_name)
  
  mask_path <- find_first_file(
    folder,
    explicit_name = args$cluster_mask_name,
    patterns = mask_patterns,
    exclude_regex = "overlay|debug|cell|assignment|boundary|membrane|distance|thumb"
  )
  
  rosette_path <- find_first_file(
    folder,
    explicit_name = args$rosette_table_name,
    patterns = rosette_patterns,
    exclude_regex = "cell_objects|cell_counts|geometry"
  )
  
  if (!file.exists(gray_path)) {
    stop(sprintf(
      "Could not find gray image for %s. Expected local file: %s",
      folder,
      gray_path
    ))
  }
  if (is.na(mask_path) || !file.exists(mask_path)) {
    stop(sprintf(
      "Could not find cluster mask for %s. Expected local file name/pattern in the sample folder.",
      folder
    ))
  }
  if (is.na(rosette_path) || !file.exists(rosette_path)) {
    stop(sprintf(
      "Could not find rosette/candidate-center table for %s. Expected local file name/pattern in the sample folder.",
      folder
    ))
  }
  
  messagef("[03b_geometry]   sample folder: %s", folder)
  messagef("[03b_geometry]   gray:          %s", gray_path)
  messagef("[03b_geometry]   mask:          %s", mask_path)
  messagef("[03b_geometry]   centers:       %s", rosette_path)
  
  out_dir <- file.path(args$out, image_safe)
  dbg_dir <- file.path(out_dir, "debug_geometry_cells")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dbg_dir, recursive = TRUE, showWarnings = FALSE)
  
  gray <- read_gray_png(gray_path)
  cluster_mask <- read_mask_png(mask_path, foreground = args$mask_foreground)
  ensure_same_dims(gray, cluster_mask, "gray", "cluster mask")
  
  if (args$fill_cluster_holes) {
    cluster_mask <- fill_hull(cluster_mask)
  }
  cluster_mask <- remove_small_components(cluster_mask, args$min_cluster_area_px)
  cluster_lab <- bwlabel_matrix(cluster_mask)
  
  h <- nrow(gray)
  w <- ncol(gray)
  
  rosettes <- standardize_rosettes(rosette_path, h = h, w = w, cluster_lab = cluster_lab, args = args, image_id = image_id)
  rosettes[, image_id := image_id]
  
  cluster_ids <- sort(unique(as.integer(cluster_lab[cluster_lab > 0])))
  
  if (length(cluster_ids) == 0L) {
    warning(sprintf("No clusters found in %s", folder))
  }
  
  all_cell_lab <- matrix(0L, nrow = h, ncol = w)
  all_membrane <- matrix(FALSE, nrow = h, ncol = w)
  all_membrane_score <- matrix(0, nrow = h, ncol = w)
  all_distance <- matrix(0, nrow = h, ncol = w)
  
  cell_rows <- list()
  cell_num <- 0L
  
  for (clid in cluster_ids) {
    cl_mask <- cluster_lab == clid
    if (sum(cl_mask) < args$min_cluster_area_px) next
    
    seg <- segment_one_cluster(gray, cl_mask, args)
    
    all_membrane <- all_membrane | seg$membrane_mask
    all_membrane_score <- pmax(all_membrane_score, seg$membrane_score)
    all_distance <- pmax(all_distance, seg$distance)
    
    lab <- seg$label
    local_ids <- sort(unique(as.integer(lab[lab > 0])))
    
    if (length(local_ids) == 0L) next
    
    for (lid in local_ids) {
      cmask <- lab == lid
      if (!any(cmask)) next
      
      cell_num <- cell_num + 1L
      all_cell_lab[cmask] <- cell_num
      
      m <- measure_cell(cmask, cl_mask, args)
      
      cell_rows[[length(cell_rows) + 1L]] <- data.table(
        image_id = image_id,
        mask_cluster_id = as.integer(clid),
        cell_num = as.integer(cell_num),
        cell_id = sprintf("%s_cl%03d_cell%04d", image_safe, as.integer(clid), as.integer(cell_num)),
        local_cell_label = as.integer(lid),
        cell_area_px = as.numeric(m$area),
        cell_centroid_x = as.numeric(m$centroid_x),
        cell_centroid_y = as.numeric(m$centroid_y),
        cell_perimeter_px = as.numeric(m$perimeter),
        cell_circularity = as.numeric(m$circularity),
        cell_solidity = as.numeric(m$solidity),
        cell_max_radius_px = as.numeric(m$max_radius),
        cell_mean_radius_px = as.numeric(m$mean_radius),
        cell_edge_contact_fraction = as.numeric(m$edge_contact_fraction),
        cell_bbox_width_px = as.numeric(m$bbox_width_px),
        cell_bbox_height_px = as.numeric(m$bbox_height_px),
        cell_aspect_ratio = as.numeric(m$aspect_ratio),
        cell_completeness_score = as.numeric(m$completeness_score),
        accepted_cell = as.logical(m$accepted),
        reject_reason = as.character(m$reject_reason)
      )
    }
  }
  
  cells <- if (length(cell_rows) > 0L) {
    data.table::rbindlist(cell_rows, fill = TRUE)
  } else {
    data.table(
      image_id = character(),
      mask_cluster_id = integer(),
      cell_num = integer(),
      cell_id = character(),
      local_cell_label = integer(),
      cell_area_px = numeric(),
      cell_centroid_x = numeric(),
      cell_centroid_y = numeric(),
      cell_perimeter_px = numeric(),
      cell_circularity = numeric(),
      cell_solidity = numeric(),
      cell_max_radius_px = numeric(),
      cell_mean_radius_px = numeric(),
      cell_edge_contact_fraction = numeric(),
      cell_bbox_width_px = numeric(),
      cell_bbox_height_px = numeric(),
      cell_aspect_ratio = numeric(),
      cell_completeness_score = numeric(),
      accepted_cell = logical(),
      reject_reason = character()
    )
  }
  
  cells <- assign_cells_to_rosettes(cells, all_cell_lab, rosettes, args)
  
  counts <- summarize_rosettes(rosettes, cells, image_id)
  
  cells_path <- file.path(out_dir, "cell_objects_geometry.tsv")
  counts_path <- file.path(out_dir, "rosette_cell_counts_geometry.tsv")
  rosettes_path <- file.path(out_dir, "rosettes_standardized_geometry.tsv")
  
  data.table::fwrite(cells, cells_path, sep = "\t")
  data.table::fwrite(counts, counts_path, sep = "\t")
  data.table::fwrite(rosettes, rosettes_path, sep = "\t")
  
  if (args$debug) {
    write_png01(all_membrane_score, file.path(dbg_dir, "membrane_dark_ridge_score.png"))
    png::writePNG(all_membrane * 1, target = file.path(dbg_dir, "membrane_candidate_mask.png"))
    write_png01(all_distance, file.path(dbg_dir, "distance_map_used_for_watershed.png"))
    png::writePNG((all_cell_lab > 0) * 1, target = file.path(dbg_dir, "segmented_cell_objects_binary.png"))
  }
  
  write_overlay(
    gray = gray,
    cluster_lab = cluster_lab,
    cell_lab = all_cell_lab,
    membrane_mask = all_membrane,
    rosettes = rosettes,
    cells = cells,
    path = file.path(dbg_dir, "geometry_cells_full_overlay.png"),
    args = args,
    mode = "full"
  )
  
  # A cleaner boundary-only overlay can be useful when labels get visually busy.
  old_draw <- args$draw_cell_ids
  old_rej <- args$draw_rejected
  args$draw_cell_ids <- FALSE
  args$draw_rejected <- TRUE
  write_overlay(
    gray = gray,
    cluster_lab = cluster_lab,
    cell_lab = all_cell_lab,
    membrane_mask = all_membrane,
    rosettes = rosettes,
    cells = cells,
    path = file.path(dbg_dir, "geometry_cells_boundaries_overlay.png"),
    args = args,
    mode = "boundaries"
  )
  args$draw_cell_ids <- old_draw
  args$draw_rejected <- old_rej
  
  messagef(
    "[03b_geometry] %s: clusters=%d, segmented_objects=%d, accepted_cells=%d, assigned_cells=%d",
    image_id,
    length(cluster_ids),
    nrow(cells),
    sum(cells$accepted_cell, na.rm = TRUE),
    sum(cells$accepted_cell & !is.na(cells$assigned_rosette_id), na.rm = TRUE)
  )
  
  list(cells = cells, counts = counts)
}

# ----------------------------- run -------------------------------------------

folders <- find_result_folders(
  input_dir = args$input,
  args = args
)

if (length(folders) == 0L) {
  stop(sprintf(
    paste0(
      "No image result folders found under %s.\n",
      "Expected --input to be results/<experiment> or results/<experiment>/<samplefolder>.\n",
      "Each sample folder must directly contain the gray image, cluster mask, and rosette/candidate-center TSV."
    ),
    args$input
  ))
}

if (!is.na(args$max_folders)) {
  folders <- folders[seq_len(min(length(folders), args$max_folders))]
}

messagef("[03b_geometry] Found %d image result folder(s).", length(folders))

all_cells <- list()
all_counts <- list()

for (folder in folders) {
  if (args$keep_going) {
    res <- tryCatch(
      process_folder(folder, args),
      error = function(e) {
        warning(sprintf("[03b_geometry] Failed on %s: %s", folder, conditionMessage(e)))
        NULL
      }
    )
  } else {
    res <- process_folder(folder, args)
  }
  
  if (!is.null(res)) {
    all_cells[[length(all_cells) + 1L]] <- res$cells
    all_counts[[length(all_counts) + 1L]] <- res$counts
  }
}

if (length(all_cells) > 0L) {
  data.table::fwrite(
    data.table::rbindlist(all_cells, fill = TRUE),
    file.path(args$out, "all_cell_objects_geometry.tsv"),
    sep = "\t"
  )
}

if (length(all_counts) > 0L) {
  data.table::fwrite(
    data.table::rbindlist(all_counts, fill = TRUE),
    file.path(args$out, "all_rosette_cell_counts_geometry.tsv"),
    sep = "\t"
  )
}

messagef("[03b_geometry] Done. Output: %s", args$out)
