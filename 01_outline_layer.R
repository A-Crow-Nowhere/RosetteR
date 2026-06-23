#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EBImage)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  opts <- list(
    input = NULL,
    out = "results/outline_v01",
    pattern = "\\.(png|jpg|jpeg|tif|tiff)$",
    recursive = FALSE,
    
    # foreground = "dark" is usually right for black phase-contrast outlines.
    # use "bright" if cells are bright on dark background.
    # use "auto" to try both and choose the one with foreground fraction
    # closest to target_fraction.
    foreground = "dark",
    target_fraction = 0.15,
    
    # threshold method: adaptive or otsu
    method = "adaptive",
    window = 41,
    offset = 0.03,
    
    # preprocessing
    bg_sigma = 0,
    smooth_sigma = 1,
    
    # morphology
    open_radius = 1,
    close_radius = 3,
    fill_holes = TRUE,
    
    # object filtering/classification
    min_area = 50,
    single_area_max = NA
  )
  
  if (length(args) == 0) {
    return(opts)
  }
  
  if (length(args) %% 2 != 0) {
    stop("Arguments must be supplied as --key value pairs.")
  }
  
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    val <- args[i + 1]
    
    if (!key %in% names(opts)) {
      stop("Unknown option: --", key)
    }
    
    old <- opts[[key]]
    
    if (is.logical(old)) {
      opts[[key]] <- tolower(val) %in% c("true", "t", "1", "yes", "y")
    } else if (is.numeric(old)) {
      opts[[key]] <- if (tolower(val) %in% c("na", "nan", "null")) NA else as.numeric(val)
    } else {
      opts[[key]] <- val
    }
  }
  
  if (is.null(opts$input)) {
    stop("Please provide --input as an image file or image directory.")
  }
  
  opts
}

safe_stem <- function(path) {
  x <- basename(path)
  x <- sub("\\.[^.]+$", "", x)
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

as_gray_image <- function(img) {
  if (length(dim(img)) == 3) {
    img <- channel(img, "gray")
  }
  normalize(img)
}

make_brush_safe <- function(radius) {
  radius <- as.integer(radius)
  if (is.na(radius) || radius <= 0) {
    return(NULL)
  }
  makeBrush(size = radius * 2 + 1, shape = "disc")
}

segment_once <- function(gray, opts, foreground = "dark") {
  work <- gray
  
  if (foreground == "dark") {
    work <- 1 - work
  }
  
  if (!is.na(opts$bg_sigma) && opts$bg_sigma > 0) {
    bg <- gblur(work, sigma = opts$bg_sigma)
    work <- normalize(work - bg)
  }
  
  if (!is.na(opts$smooth_sigma) && opts$smooth_sigma > 0) {
    work <- gblur(work, sigma = opts$smooth_sigma)
  }
  
  work <- normalize(work)
  
  if (opts$method == "adaptive") {
    mask <- thresh(
      work,
      w = as.integer(opts$window),
      h = as.integer(opts$window),
      offset = opts$offset
    )
  } else if (opts$method == "otsu") {
    mask <- work > otsu(work)
  } else {
    stop("Unknown --method. Use adaptive or otsu.")
  }
  
  open_brush <- make_brush_safe(opts$open_radius)
  close_brush <- make_brush_safe(opts$close_radius)
  
  if (!is.null(open_brush)) {
    mask <- opening(mask, open_brush)
  }
  
  if (!is.null(close_brush)) {
    mask <- closing(mask, close_brush)
  }
  
  if (isTRUE(opts$fill_holes)) {
    mask <- fillHull(mask)
  }
  
  lab <- bwlabel(mask)
  lab_mat <- imageData(lab)
  ids <- sort(unique(as.integer(lab_mat)))
  ids <- ids[ids > 0]
  
  if (length(ids) > 0) {
    areas <- tabulate(as.integer(lab_mat[lab_mat > 0]))
    keep <- which(areas >= opts$min_area)
    
    lab_mat[!(lab_mat %in% keep)] <- 0
    mask <- Image(lab_mat > 0, colormode = Grayscale)
  } else {
    mask <- Image(lab_mat > 0, colormode = Grayscale)
  }
  
  mask <- normalize(mask > 0)
  lab <- bwlabel(mask)
  
  list(
    processed = work,
    mask = mask,
    labels = lab,
    foreground = foreground,
    foreground_fraction = mean(imageData(mask) > 0)
  )
}

choose_segmentation <- function(gray, opts) {
  if (opts$foreground %in% c("dark", "bright")) {
    return(segment_once(gray, opts, opts$foreground))
  }
  
  if (opts$foreground != "auto") {
    stop("--foreground must be dark, bright, or auto.")
  }
  
  dark <- segment_once(gray, opts, "dark")
  bright <- segment_once(gray, opts, "bright")
  
  dark_score <- abs(dark$foreground_fraction - opts$target_fraction)
  bright_score <- abs(bright$foreground_fraction - opts$target_fraction)
  
  if (dark_score <= bright_score) dark else bright
}

make_outline <- function(mask) {
  brush <- makeBrush(3, shape = "disc")
  eroded <- erode(mask, brush)
  outline <- (mask > 0) & !(eroded > 0)
  normalize(outline)
}

make_overlay <- function(gray, outline) {
  g <- imageData(normalize(gray))
  o <- imageData(outline) > 0
  
  rgb <- array(0, dim = c(dim(g)[1], dim(g)[2], 3))
  rgb[, , 1] <- g
  rgb[, , 2] <- g
  rgb[, , 3] <- g
  
  # red outline
  rgb[, , 1][o] <- 1
  rgb[, , 2][o] <- 0
  rgb[, , 3][o] <- 0
  
  Image(rgb, colormode = Color)
}

object_table <- function(labels, gray, image_id, image_file, opts) {
  lab_mat <- imageData(labels)
  ids <- sort(unique(as.integer(lab_mat)))
  ids <- ids[ids > 0]
  
  if (length(ids) == 0) {
    return(data.frame())
  }
  
  shape <- as.data.frame(computeFeatures.shape(labels))
  moment <- as.data.frame(computeFeatures.moment(labels))
  basic <- as.data.frame(computeFeatures.basic(labels, gray))
  
  df <- cbind(shape, moment, basic)
  df$object_id <- as.integer(rownames(df))
  
  bbox_list <- lapply(df$object_id, function(id) {
    ij <- which(lab_mat == id, arr.ind = TRUE)
    
    data.frame(
      object_id = id,
      bbox_xmin = min(ij[, 1]),
      bbox_xmax = max(ij[, 1]),
      bbox_ymin = min(ij[, 2]),
      bbox_ymax = max(ij[, 2])
    )
  })
  
  bbox <- do.call(rbind, bbox_list)
  df <- merge(df, bbox, by = "object_id", all.x = TRUE)
  
  df$image_id <- image_id
  df$image_file <- image_file
  
  # Crude first-pass type assignment. This will become smarter later.
  if (!is.na(opts$single_area_max)) {
    df$object_type <- ifelse(
      df$s.area <= opts$single_area_max,
      "candidate_single_swimmer",
      "candidate_rosette_or_cluster"
    )
  } else {
    df$object_type <- "unclassified_connected_object"
  }
  
  # Cleaner column order
  front <- c(
    "image_id", "image_file", "object_id", "object_type",
    "s.area", "s.perimeter",
    "m.cx", "m.cy",
    "m.majoraxis", "m.eccentricity", "m.theta",
    "bbox_xmin", "bbox_xmax", "bbox_ymin", "bbox_ymax"
  )
  
  front <- front[front %in% names(df)]
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}

process_image <- function(image_file, opts) {
  image_id <- safe_stem(image_file)
  image_out <- file.path(opts$out, image_id)
  dir.create(image_out, recursive = TRUE, showWarnings = FALSE)
  
  message("[01_outline] Reading: ", image_file)
  
  img <- readImage(image_file)
  gray <- as_gray_image(img)
  
  seg <- choose_segmentation(gray, opts)
  outline <- make_outline(seg$mask)
  overlay <- make_overlay(gray, outline)
  
  writeImage(normalize(gray), file.path(image_out, "gray.png"))
  writeImage(normalize(seg$processed), file.path(image_out, "processed.png"))
  writeImage(normalize(seg$mask), file.path(image_out, "mask.binary.png"))
  writeImage(normalize(outline), file.path(image_out, "outline.png"))
  writeImage(overlay, file.path(image_out, "outline_overlay.png"))
  
  obj <- object_table(
    labels = seg$labels,
    gray = gray,
    image_id = image_id,
    image_file = image_file,
    opts = opts
  )
  
  write.table(
    obj,
    file = file.path(image_out, "objects.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  summary <- data.frame(
    image_id = image_id,
    image_file = image_file,
    width_px = dim(imageData(gray))[1],
    height_px = dim(imageData(gray))[2],
    foreground = seg$foreground,
    foreground_fraction = seg$foreground_fraction,
    n_objects = nrow(obj),
    n_candidate_single_swimmers = sum(obj$object_type == "candidate_single_swimmer", na.rm = TRUE),
    n_candidate_rosette_or_clusters = sum(obj$object_type == "candidate_rosette_or_cluster", na.rm = TRUE),
    min_area = opts$min_area,
    single_area_max = opts$single_area_max,
    threshold_method = opts$method,
    threshold_window = opts$window,
    threshold_offset = opts$offset
  )
  
  list(objects = obj, summary = summary)
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
  
  if (file.info(opts$input)$isdir) {
    files <- list.files(
      opts$input,
      pattern = opts$pattern,
      recursive = opts$recursive,
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    files <- opts$input
  }
  
  if (length(files) == 0) {
    stop("No image files found.")
  }
  
  message("[01_outline] Found ", length(files), " image(s).")
  
  res <- lapply(files, process_image, opts = opts)
  
  all_objects <- do.call(rbind, lapply(res, `[[`, "objects"))
  all_summary <- do.call(rbind, lapply(res, `[[`, "summary"))
  
  write.table(
    all_objects,
    file = file.path(opts$out, "all_objects.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    all_summary,
    file = file.path(opts$out, "image_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("[01_outline] Done.")
  message("[01_outline] Results: ", opts$out)
}

main()