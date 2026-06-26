#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(png)
})

SCRIPT_VERSION <- "05_mask_confluence_summary_v02_outline_only"

option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "Input directory to recursively search for binary mask images."
  ),
  make_option(
    c("--out"),
    type = "character",
    default = "results/confluence_summary",
    help = "Output directory for confluence TSV files."
  ),
  make_option(
    c("--mask_name"),
    type = "character",
    default = "mask.binary.png",
    help = "Exact filename of binary cluster mask images to use."
  ),
  make_option(
    c("--outline_folder"),
    type = "character",
    default = "01_outline",
    help = "Only use masks located inside this folder name, e.g. 01_outline."
  ),
  make_option(
    c("--foreground"),
    type = "character",
    default = "auto",
    help = "Foreground polarity: auto, bright, or dark. Use bright if masks are white shapes on black background."
  ),
  make_option(
    c("--threshold"),
    type = "double",
    default = 0.5,
    help = "Pixel threshold for binary mask calling after scaling image to 0-1."
  ),
  make_option(
    c("--default_replicate"),
    type = "character",
    default = "replicate_01",
    help = "Replicate ID to use when no replicate-like folder is found."
  ),
  make_option(
    c("--sample_regex"),
    type = "character",
    default = "",
    help = "Optional regex for sample folders. If empty, sample is inferred from path."
  ),
  make_option(
    c("--replicate_regex"),
    type = "character",
    default = "^replicate[_-]?[0-9]+$",
    help = "Regex used to identify replicate folders."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("[", SCRIPT_VERSION, "] ERROR: --input is required.", call. = FALSE)
}

if (!dir.exists(opt$input)) {
  stop("[", SCRIPT_VERSION, "] ERROR: input directory does not exist: ", opt$input, call. = FALSE)
}

if (!opt$foreground %in% c("auto", "bright", "dark")) {
  stop("[", SCRIPT_VERSION, "] ERROR: --foreground must be one of: auto, bright, dark", call. = FALSE)
}

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

message("[", SCRIPT_VERSION, "] Input: ", opt$input)
message("[", SCRIPT_VERSION, "] Output: ", opt$out)
message("[", SCRIPT_VERSION, "] Searching for masks named: ", opt$mask_name)
message("[", SCRIPT_VERSION, "] Restricting search to folder named: ", opt$outline_folder)

escape_regex <- function(x) {
  gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", x)
}

mask_files <- list.files(
  opt$input,
  pattern = paste0("^", escape_regex(opt$mask_name), "$"),
  recursive = TRUE,
  full.names = TRUE
)

if (length(mask_files) > 0) {
  mask_files <- mask_files[
    grepl(
      paste0("(^|/)", escape_regex(opt$outline_folder), "(/|$)"),
      normalizePath(mask_files, winslash = "/", mustWork = FALSE)
    )
  ]
}

if (length(mask_files) == 0) {
  stop(
    "[", SCRIPT_VERSION, "] ERROR: No mask files found named '", opt$mask_name,
    "' inside outline folder '", opt$outline_folder,
    "' under: ", opt$input,
    call. = FALSE
  )
}

message("[", SCRIPT_VERSION, "] Found ", length(mask_files), " mask file(s) in outline folder(s).")

normalize_path_parts <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parts <- unlist(strsplit(path, "/", fixed = TRUE))
  parts[nzchar(parts)]
}

infer_ids_from_path <- function(mask_path, input_root, replicate_regex, default_replicate, sample_regex = "") {
  mask_path_norm <- normalizePath(mask_path, winslash = "/", mustWork = FALSE)
  input_root_norm <- normalizePath(input_root, winslash = "/", mustWork = FALSE)

  rel <- sub(
    paste0("^", escape_regex(input_root_norm), "/?"),
    "",
    mask_path_norm
  )

  parts <- unlist(strsplit(rel, "/", fixed = TRUE))
  parts <- parts[nzchar(parts)]

  filename <- basename(mask_path_norm)
  image_id <- basename(dirname(mask_path_norm))

  replicate_id <- default_replicate
  sample_id <- NA_character_

  rep_idx <- which(grepl(replicate_regex, parts, ignore.case = TRUE))

  if (length(rep_idx) > 0) {
    rep_idx <- rep_idx[1]
    replicate_id <- parts[rep_idx]

    if (rep_idx > 1) {
      sample_id <- parts[rep_idx - 1]
    }
  }

  if (is.na(sample_id) || !nzchar(sample_id)) {
    if (nzchar(sample_regex)) {
      sample_hits <- parts[grepl(sample_regex, parts)]
      if (length(sample_hits) > 0) {
        sample_id <- sample_hits[1]
      }
    }
  }

  if (is.na(sample_id) || !nzchar(sample_id)) {
    # Common expected structures:
    #
    # input/sample/replicate_01/01_outline/image_id/mask.binary.png
    # input/sample/01_outline/image_id/mask.binary.png
    # input/01_outline/sample_or_image/mask.binary.png
    #
    # If no replicate folder exists, choose a sensible upstream folder.
    outline_idx <- which(parts == opt$outline_folder)

    if (length(outline_idx) > 0) {
      outline_idx <- outline_idx[1]

      if (outline_idx > 1) {
        sample_id <- parts[outline_idx - 1]
      } else if (length(parts) >= outline_idx + 1) {
        sample_id <- parts[outline_idx + 1]
      } else {
        sample_id <- "sample_unknown"
      }
    } else if (length(parts) >= 3) {
      sample_id <- parts[length(parts) - 2]
    } else if (length(parts) >= 2) {
      sample_id <- parts[length(parts) - 1]
    } else {
      sample_id <- "sample_unknown"
    }
  }

  list(
    sample_id = sample_id,
    replicate_id = replicate_id,
    image_id = image_id,
    relative_path = rel,
    filename = filename
  )
}

read_mask01 <- function(mask_file) {
  img <- png::readPNG(mask_file)

  # Convert RGB/RGBA to grayscale-like matrix if needed.
  if (length(dim(img)) == 3) {
    img <- img[, , 1]
  }

  img <- as.matrix(img)

  # png::readPNG usually returns 0-1, but keep this robust.
  maxv <- suppressWarnings(max(img, na.rm = TRUE))
  if (is.finite(maxv) && maxv > 1) {
    img <- img / maxv
  }

  img[is.na(img)] <- 0
  img
}

call_foreground <- function(img, foreground = "auto", threshold = 0.5) {
  if (foreground == "bright") {
    fg <- img >= threshold
  } else if (foreground == "dark") {
    fg <- img <= threshold
  } else {
    # Auto assumes the minority pixel class is foreground.
    bright_fg <- img >= threshold
    dark_fg <- img <= threshold

    bright_frac <- mean(bright_fg)
    dark_frac <- mean(dark_fg)

    # Filled cluster masks usually occupy less than the whole image.
    # Pick the smaller class unless nearly tied.
    if (bright_frac <= dark_frac) {
      fg <- bright_fg
    } else {
      fg <- dark_fg
    }
  }

  fg
}

rows <- vector("list", length(mask_files))

for (i in seq_along(mask_files)) {
  f <- mask_files[i]

  ids <- infer_ids_from_path(
    mask_path = f,
    input_root = opt$input,
    replicate_regex = opt$replicate_regex,
    default_replicate = opt$default_replicate,
    sample_regex = opt$sample_regex
  )

  img <- read_mask01(f)
  fg <- call_foreground(
    img,
    foreground = opt$foreground,
    threshold = opt$threshold
  )

  n_total_px <- length(fg)
  n_foreground_px <- sum(fg, na.rm = TRUE)
  confluence_fraction <- n_foreground_px / n_total_px
  confluence_percent <- 100 * confluence_fraction

  rows[[i]] <- data.table(
    sample_id = ids$sample_id,
    replicate_id = ids$replicate_id,
    image_id = ids$image_id,
    mask_file = f,
    relative_path = ids$relative_path,
    image_width_px = ncol(img),
    image_height_px = nrow(img),
    total_pixels = n_total_px,
    foreground_pixels = n_foreground_px,
    confluence_fraction = confluence_fraction,
    confluence_percent = confluence_percent
  )
}

per_image <- rbindlist(rows, fill = TRUE)

setorder(per_image, sample_id, replicate_id, image_id, relative_path)

summary_by_sample_replicate <- per_image[
  ,
  .(
    n_images = .N,
    mean_confluence_percent = mean(confluence_percent, na.rm = TRUE),
    median_confluence_percent = median(confluence_percent, na.rm = TRUE),
    sd_confluence_percent = sd(confluence_percent, na.rm = TRUE),
    min_confluence_percent = min(confluence_percent, na.rm = TRUE),
    max_confluence_percent = max(confluence_percent, na.rm = TRUE),
    total_foreground_pixels = sum(foreground_pixels, na.rm = TRUE),
    total_pixels = sum(total_pixels, na.rm = TRUE),
    pooled_confluence_percent = 100 * sum(foreground_pixels, na.rm = TRUE) / sum(total_pixels, na.rm = TRUE)
  ),
  by = .(sample_id, replicate_id)
]

setorder(summary_by_sample_replicate, sample_id, replicate_id)

per_image_out <- file.path(opt$out, "mask_confluence.per_image.tsv")
summary_out <- file.path(opt$out, "mask_confluence.by_sample_replicate.tsv")

fwrite(per_image, per_image_out, sep = "\t")
fwrite(summary_by_sample_replicate, summary_out, sep = "\t")

message("[", SCRIPT_VERSION, "] Wrote per-image table: ", per_image_out)
message("[", SCRIPT_VERSION, "] Wrote sample/replicate summary: ", summary_out)

message("[", SCRIPT_VERSION, "] Done.")
