#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EBImage)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  opts <- list(
    input = "results/outline_v01",
    out = "results/centers_v01",
    
    mask_name = "mask.binary.png",
    gray_name = "gray.png",
    
    min_blob_area = 500,
    
    max_boundary_points = 2500,
    max_seed_points = 600,
    
    neighborhood_radius = 60,
    min_arc_points = 25,
    min_radius = 15,
    max_radius = 300,
    min_arc_angle = 35,
    max_fit_error = 0.12,
    
    require_center_inside = TRUE,
    center_inside_erode_radius = 0,
    
    ray_angles = 32,
    radial_bin_width_deg = 12,
    
    center_merge_dist = 20,
    min_total_score = 0.45,
    min_merged_support = 2,
    
    merge_method = "weighted",
    merge_weight_power = 2,
    merge_min_score_fraction = 0.50,
    require_weighted_center_inside = TRUE,
    
    overlap_imbalance_threshold = 0.45,
    
    draw_rejected_top_n = 100
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
  
  if (!file.exists(opts$input)) {
    stop("Input does not exist: ", opts$input)
  }
  
  opts$ray_angles <- as.integer(opts$ray_angles)
  
  if (opts$ray_angles %% 2 != 0) {
    stop("--ray_angles must be an even number.")
  }
  
  if (!opts$merge_method %in% c("best", "weighted")) {
    stop("--merge_method must be either 'best' or 'weighted'.")
  }
  
  opts
}

clamp01 <- function(x) {
  pmax(0, pmin(1, x))
}

safe_stem <- function(x) {
  x <- basename(x)
  x <- sub("\\.[^.]+$", "", x)
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

as_gray_image <- function(img) {
  if (length(dim(img)) == 3) {
    img <- channel(img, "gray")
  }
  
  normalize(img)
}

make_outline <- function(mask) {
  brush <- makeBrush(3, shape = "disc")
  eroded <- erode(mask, brush)
  outline <- (mask > 0) & !(eroded > 0)
  normalize(outline)
}

make_center_valid_mask <- function(label_mat, blob_id, erode_radius = 0) {
  blob_mask <- label_mat == blob_id
  
  if (!is.na(erode_radius) && erode_radius > 0) {
    brush <- makeBrush(
      size = as.integer(erode_radius) * 2 + 1,
      shape = "disc"
    )
    
    blob_mask <- imageData(
      erode(Image(blob_mask, colormode = Grayscale), brush)
    ) > 0
  }
  
  blob_mask
}

point_inside_mask <- function(mask_mat, x, y) {
  xi <- round(x)
  yi <- round(y)
  
  if (
    xi < 1 || xi > nrow(mask_mat) ||
    yi < 1 || yi > ncol(mask_mat)
  ) {
    return(FALSE)
  }
  
  isTRUE(mask_mat[xi, yi])
}

center_is_inside_blob <- function(
    label_mat,
    blob_id,
    center_x,
    center_y,
    erode_radius = 0
) {
  valid_mask <- make_center_valid_mask(
    label_mat = label_mat,
    blob_id = blob_id,
    erode_radius = erode_radius
  )
  
  point_inside_mask(valid_mask, center_x, center_y)
}

subsample_df <- function(df, max_n) {
  if (nrow(df) <= max_n) {
    return(df)
  }
  
  idx <- unique(round(seq(1, nrow(df), length.out = max_n)))
  df[idx, , drop = FALSE]
}

angle_span_deg <- function(theta) {
  theta <- sort((theta + 2 * pi) %% (2 * pi))
  
  if (length(theta) < 2) {
    return(0)
  }
  
  gaps <- diff(c(theta, theta[1] + 2 * pi))
  span <- 2 * pi - max(gaps)
  
  span * 180 / pi
}

fit_circle_kasa <- function(x, y) {
  if (length(x) < 3) {
    return(NULL)
  }
  
  A <- cbind(x, y, 1)
  b <- -(x^2 + y^2)
  
  fit <- tryCatch(
    solve(t(A) %*% A, t(A) %*% b),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  D <- fit[1]
  E <- fit[2]
  F <- fit[3]
  
  cx <- -D / 2
  cy <- -E / 2
  
  r2 <- (D^2 + E^2) / 4 - F
  
  if (!is.finite(r2) || r2 <= 0) {
    return(NULL)
  }
  
  r <- sqrt(r2)
  d <- sqrt((x - cx)^2 + (y - cy)^2)
  
  rms <- sqrt(mean((d - r)^2))
  rms_norm <- rms / r
  
  list(
    center_x = as.numeric(cx),
    center_y = as.numeric(cy),
    radius = as.numeric(r),
    fit_error_px = as.numeric(rms),
    fit_error_norm = as.numeric(rms_norm)
  )
}

radial_metrics <- function(center_x, center_y, boundary_pts, opts) {
  n_angles <- opts$ray_angles
  
  angles <- seq(0, 2 * pi, length.out = n_angles + 1)
  angles <- angles[-length(angles)]
  
  dx <- boundary_pts$x - center_x
  dy <- boundary_pts$y - center_y
  
  pt_angles <- atan2(dy, dx)
  pt_dist <- sqrt(dx^2 + dy^2)
  
  half_width <- (opts$radial_bin_width_deg * pi / 180) / 2
  
  ray_dist <- rep(NA_real_, n_angles)
  
  for (i in seq_along(angles)) {
    a <- angles[i]
    adiff <- atan2(sin(pt_angles - a), cos(pt_angles - a))
    hit <- which(abs(adiff) <= half_width)
    
    if (length(hit) > 0) {
      ray_dist[i] <- median(pt_dist[hit], na.rm = TRUE)
    }
  }
  
  coverage <- mean(!is.na(ray_dist))
  
  n_pairs <- n_angles / 2
  pair_imbalance <- rep(NA_real_, n_pairs)
  
  for (i in seq_len(n_pairs)) {
    d1 <- ray_dist[i]
    d2 <- ray_dist[i + n_pairs]
    
    if (is.finite(d1) && is.finite(d2) && (d1 + d2) > 0) {
      pair_imbalance[i] <- abs(d1 - d2) / (d1 + d2)
    }
  }
  
  if (all(is.na(pair_imbalance))) {
    radial_balance_score <- NA_real_
    max_pair_imbalance <- NA_real_
    median_pair_imbalance <- NA_real_
  } else {
    median_pair_imbalance <- median(pair_imbalance, na.rm = TRUE)
    max_pair_imbalance <- max(pair_imbalance, na.rm = TRUE)
    radial_balance_score <- 1 - median_pair_imbalance
  }
  
  if (sum(is.finite(ray_dist)) >= 3) {
    radial_distance_cv <- sd(ray_dist, na.rm = TRUE) / mean(ray_dist, na.rm = TRUE)
  } else {
    radial_distance_cv <- NA_real_
  }
  
  list(
    radial_coverage = coverage,
    radial_balance_score = radial_balance_score,
    median_pair_imbalance = median_pair_imbalance,
    max_pair_imbalance = max_pair_imbalance,
    radial_distance_cv = radial_distance_cv
  )
}

empty_raw_table <- function() {
  data.frame(
    image_id = character(),
    parent_blob_id = integer(),
    seed_id = integer(),
    seed_x = numeric(),
    seed_y = numeric(),
    center_x = numeric(),
    center_y = numeric(),
    fitted_radius_px = numeric(),
    fit_error_px = numeric(),
    fit_error_norm = numeric(),
    arc_points = integer(),
    arc_span_deg = numeric(),
    center_inside_blob = logical(),
    radial_coverage = numeric(),
    radial_balance_score = numeric(),
    median_pair_imbalance = numeric(),
    max_pair_imbalance = numeric(),
    radial_distance_cv = numeric(),
    total_score = numeric(),
    stringsAsFactors = FALSE
  )
}

fit_candidates_for_blob <- function(
    image_id,
    blob_id,
    label_mat,
    boundary_pts,
    opts
) {
  if (nrow(boundary_pts) < opts$min_arc_points) {
    return(empty_raw_table())
  }
  
  boundary_fit <- subsample_df(boundary_pts, opts$max_boundary_points)
  seed_pts <- subsample_df(boundary_fit, opts$max_seed_points)
  
  valid_center_mask <- make_center_valid_mask(
    label_mat = label_mat,
    blob_id = blob_id,
    erode_radius = opts$center_inside_erode_radius
  )
  
  raw_list <- list()
  raw_i <- 0
  
  for (s in seq_len(nrow(seed_pts))) {
    sx <- seed_pts$x[s]
    sy <- seed_pts$y[s]
    
    d2 <- (boundary_fit$x - sx)^2 + (boundary_fit$y - sy)^2
    local <- boundary_fit[d2 <= opts$neighborhood_radius^2, , drop = FALSE]
    
    if (nrow(local) < opts$min_arc_points) {
      next
    }
    
    fit <- fit_circle_kasa(local$x, local$y)
    
    if (is.null(fit)) {
      next
    }
    
    if (!is.finite(fit$radius)) {
      next
    }
    
    if (fit$radius < opts$min_radius || fit$radius > opts$max_radius) {
      next
    }
    
    if (fit$fit_error_norm > opts$max_fit_error) {
      next
    }
    
    theta <- atan2(local$y - fit$center_y, local$x - fit$center_x)
    arc_span <- angle_span_deg(theta)
    
    if (arc_span < opts$min_arc_angle) {
      next
    }
    
    center_inside <- point_inside_mask(
      valid_center_mask,
      fit$center_x,
      fit$center_y
    )
    
    if (isTRUE(opts$require_center_inside) && !isTRUE(center_inside)) {
      next
    }
    
    rm <- radial_metrics(
      center_x = fit$center_x,
      center_y = fit$center_y,
      boundary_pts = boundary_fit,
      opts = opts
    )
    
    fit_score <- clamp01(1 - fit$fit_error_norm / opts$max_fit_error)
    arc_score <- clamp01(arc_span / 120)
    coverage_score <- clamp01(rm$radial_coverage)
    
    balance_score <- ifelse(
      is.finite(rm$radial_balance_score),
      clamp01(rm$radial_balance_score),
      0
    )
    
    inside_score <- ifelse(center_inside, 1, 0)
    
    total_score <- (
      0.40 * fit_score +
        0.25 * arc_score +
        0.15 * coverage_score +
        0.10 * inside_score +
        0.10 * balance_score
    )
    
    raw_i <- raw_i + 1
    
    raw_list[[raw_i]] <- data.frame(
      image_id = image_id,
      parent_blob_id = blob_id,
      seed_id = s,
      seed_x = sx,
      seed_y = sy,
      center_x = fit$center_x,
      center_y = fit$center_y,
      fitted_radius_px = fit$radius,
      fit_error_px = fit$fit_error_px,
      fit_error_norm = fit$fit_error_norm,
      arc_points = nrow(local),
      arc_span_deg = arc_span,
      center_inside_blob = center_inside,
      radial_coverage = rm$radial_coverage,
      radial_balance_score = rm$radial_balance_score,
      median_pair_imbalance = rm$median_pair_imbalance,
      max_pair_imbalance = rm$max_pair_imbalance,
      radial_distance_cv = rm$radial_distance_cv,
      total_score = total_score,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(raw_list) == 0) {
    return(empty_raw_table())
  }
  
  do.call(rbind, raw_list)
}

empty_merged_table <- function() {
  data.frame(
    image_id = character(),
    rosette_candidate_id = integer(),
    parent_blob_id = integer(),
    
    center_x = numeric(),
    center_y = numeric(),
    fitted_radius_px = numeric(),
    
    best_center_x = numeric(),
    best_center_y = numeric(),
    best_fitted_radius_px = numeric(),
    
    weighted_center_x = numeric(),
    weighted_center_y = numeric(),
    weighted_fitted_radius_px = numeric(),
    
    center_shift_from_best_px = numeric(),
    
    confidence = numeric(),
    accepted = logical(),
    possible_overlap_break = logical(),
    
    n_raw_support = integer(),
    n_weighted_support = integer(),
    representative_method = character(),
    
    best_fit_error_norm = numeric(),
    best_arc_span_deg = numeric(),
    best_arc_points = integer(),
    
    center_inside_blob = logical(),
    weighted_center_inside_blob = logical(),
    
    radial_coverage = numeric(),
    radial_balance_score = numeric(),
    median_pair_imbalance = numeric(),
    max_pair_imbalance = numeric(),
    radial_distance_cv = numeric(),
    
    weighted_mean_score = numeric(),
    best_total_score = numeric(),
    
    stringsAsFactors = FALSE
  )
}

merge_candidates <- function(raw, opts, label_mat = NULL) {
  if (nrow(raw) == 0) {
    return(empty_merged_table())
  }
  
  if (!opts$merge_method %in% c("best", "weighted")) {
    stop("--merge_method must be either 'best' or 'weighted'.")
  }
  
  raw <- raw[order(raw$parent_blob_id, -raw$total_score), , drop = FALSE]
  
  merged_list <- list()
  out_i <- 0
  
  for (blob_id in sort(unique(raw$parent_blob_id))) {
    sub <- raw[raw$parent_blob_id == blob_id, , drop = FALSE]
    
    if (nrow(sub) == 0) {
      next
    }
    
    order_idx <- order(-sub$total_score)
    cluster_id <- rep(0L, nrow(sub))
    cid <- 0L
    
    for (idx in order_idx) {
      if (cluster_id[idx] != 0L) {
        next
      }
      
      cid <- cid + 1L
      
      d <- sqrt(
        (sub$center_x - sub$center_x[idx])^2 +
          (sub$center_y - sub$center_y[idx])^2
      )
      
      members <- which(cluster_id == 0L & d <= opts$center_merge_dist)
      cluster_id[members] <- cid
    }
    
    sub$cluster_id <- cluster_id
    
    for (cluster in sort(unique(sub$cluster_id))) {
      cl <- sub[sub$cluster_id == cluster, , drop = FALSE]
      
      if (nrow(cl) == 0) {
        next
      }
      
      best <- cl[which.max(cl$total_score), , drop = FALSE]
      best_score <- best$total_score[1]
      
      score_cutoff <- best_score * opts$merge_min_score_fraction
      
      cl_weighted <- cl[
        cl$total_score >= score_cutoff,
        ,
        drop = FALSE
      ]
      
      if (nrow(cl_weighted) == 0) {
        cl_weighted <- best
      }
      
      w <- pmax(cl_weighted$total_score, 1e-6)^opts$merge_weight_power
      
      weighted_center_x <- weighted.mean(cl_weighted$center_x, w)
      weighted_center_y <- weighted.mean(cl_weighted$center_y, w)
      weighted_radius <- weighted.mean(cl_weighted$fitted_radius_px, w)
      
      weighted_mean_score <- weighted.mean(cl_weighted$total_score, w)
      
      center_shift_from_best_px <- sqrt(
        (weighted_center_x - best$center_x)^2 +
          (weighted_center_y - best$center_y)^2
      )
      
      weighted_inside <- NA
      
      if (!is.null(label_mat)) {
        weighted_inside <- center_is_inside_blob(
          label_mat = label_mat,
          blob_id = blob_id,
          center_x = weighted_center_x,
          center_y = weighted_center_y,
          erode_radius = opts$center_inside_erode_radius
        )
      }
      
      if (opts$merge_method == "best") {
        final_center_x <- best$center_x
        final_center_y <- best$center_y
        final_radius <- best$fitted_radius_px
        final_center_inside <- best$center_inside_blob
        representative_method <- "best"
      } else {
        final_center_x <- weighted_center_x
        final_center_y <- weighted_center_y
        final_radius <- weighted_radius
        final_center_inside <- weighted_inside
        representative_method <- "weighted"
        
        if (
          isTRUE(opts$require_weighted_center_inside) &&
          isFALSE(weighted_inside)
        ) {
          final_center_x <- best$center_x
          final_center_y <- best$center_y
          final_radius <- best$fitted_radius_px
          final_center_inside <- best$center_inside_blob
          representative_method <- "best_fallback_weighted_center_outside"
        }
      }
      
      support_score <- clamp01(log1p(nrow(cl)) / log1p(10))
      
      confidence <- clamp01(
        0.55 * weighted_mean_score +
          0.30 * best_score +
          0.15 * support_score
      )
      
      possible_overlap <- is.finite(best$max_pair_imbalance) &&
        best$max_pair_imbalance >= opts$overlap_imbalance_threshold
      
      accepted <- confidence >= opts$min_total_score &&
        nrow(cl) >= opts$min_merged_support
      
      out_i <- out_i + 1
      
      merged_list[[out_i]] <- data.frame(
        image_id = best$image_id,
        rosette_candidate_id = out_i,
        parent_blob_id = best$parent_blob_id,
        
        center_x = final_center_x,
        center_y = final_center_y,
        fitted_radius_px = final_radius,
        
        best_center_x = best$center_x,
        best_center_y = best$center_y,
        best_fitted_radius_px = best$fitted_radius_px,
        
        weighted_center_x = weighted_center_x,
        weighted_center_y = weighted_center_y,
        weighted_fitted_radius_px = weighted_radius,
        
        center_shift_from_best_px = center_shift_from_best_px,
        
        confidence = confidence,
        accepted = accepted,
        possible_overlap_break = possible_overlap,
        
        n_raw_support = nrow(cl),
        n_weighted_support = nrow(cl_weighted),
        representative_method = representative_method,
        
        best_fit_error_norm = best$fit_error_norm,
        best_arc_span_deg = best$arc_span_deg,
        best_arc_points = best$arc_points,
        
        center_inside_blob = final_center_inside,
        weighted_center_inside_blob = weighted_inside,
        
        radial_coverage = best$radial_coverage,
        radial_balance_score = best$radial_balance_score,
        median_pair_imbalance = best$median_pair_imbalance,
        max_pair_imbalance = best$max_pair_imbalance,
        radial_distance_cv = best$radial_distance_cv,
        
        weighted_mean_score = weighted_mean_score,
        best_total_score = best_score,
        
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(merged_list) == 0) {
    return(empty_merged_table())
  }
  
  merged <- do.call(rbind, merged_list)
  merged$rosette_candidate_id <- seq_len(nrow(merged))
  
  merged
}

set_pixel <- function(rgb, x, y, col) {
  nx <- dim(rgb)[1]
  ny <- dim(rgb)[2]
  
  ok <- x >= 1 & x <= nx & y >= 1 & y <= ny
  
  if (!any(ok)) {
    return(rgb)
  }
  
  x <- x[ok]
  y <- y[ok]
  
  rgb[cbind(x, y, 1)] <- col[1]
  rgb[cbind(x, y, 2)] <- col[2]
  rgb[cbind(x, y, 3)] <- col[3]
  
  rgb
}

draw_cross <- function(rgb, x, y, col = c(0, 1, 0), size = 7) {
  x <- round(x)
  y <- round(y)
  
  xs <- c((x - size):(x + size), rep(x, length((y - size):(y + size))))
  ys <- c(rep(y, length((x - size):(x + size))), (y - size):(y + size))
  
  set_pixel(rgb, xs, ys, col)
}

draw_circle <- function(rgb, x, y, r, col = c(0, 1, 0), n = 720, thickness = 2) {
  if (!is.finite(r) || r <= 0) {
    return(rgb)
  }
  
  theta <- seq(0, 2 * pi, length.out = n)
  xs <- round(x + r * cos(theta))
  ys <- round(y + r * sin(theta))
  
  for (dx in seq(-thickness, thickness)) {
    for (dy in seq(-thickness, thickness)) {
      rgb <- set_pixel(rgb, xs + dx, ys + dy, col)
    }
  }
  
  rgb
}

make_center_overlay <- function(gray, outline, merged, opts, accepted_only = FALSE) {
  g <- imageData(normalize(gray))
  o <- imageData(outline) > 0
  
  rgb <- array(0, dim = c(dim(g)[1], dim(g)[2], 3))
  
  rgb[, , 1] <- g
  rgb[, , 2] <- g
  rgb[, , 3] <- g
  
  rgb[, , 1][o] <- 1
  rgb[, , 2][o] <- 0
  rgb[, , 3][o] <- 0
  
  if (nrow(merged) == 0) {
    return(Image(rgb, colormode = Color))
  }
  
  accepted <- merged[merged$accepted, , drop = FALSE]
  rejected <- merged[!merged$accepted, , drop = FALSE]
  
  if (!accepted_only && nrow(rejected) > 0) {
    rejected <- rejected[order(-rejected$confidence), , drop = FALSE]
    rejected <- head(rejected, opts$draw_rejected_top_n)
    
    for (i in seq_len(nrow(rejected))) {
      rgb <- draw_cross(
        rgb,
        rejected$center_x[i],
        rejected$center_y[i],
        col = c(0, 0.8, 1),
        size = 5
      )
    }
  }
  
  if (nrow(accepted) > 0) {
    for (i in seq_len(nrow(accepted))) {
      rgb <- draw_circle(
        rgb,
        accepted$center_x[i],
        accepted$center_y[i],
        accepted$fitted_radius_px[i],
        col = c(0, 1, 0),
        thickness = 1
      )
      
      rgb <- draw_cross(
        rgb,
        accepted$center_x[i],
        accepted$center_y[i],
        col = c(0, 1, 0),
        size = 8
      )
      
      if (isTRUE(accepted$possible_overlap_break[i])) {
        rgb <- draw_cross(
          rgb,
          accepted$center_x[i],
          accepted$center_y[i],
          col = c(1, 0, 1),
          size = 12
        )
      }
    }
  }
  
  Image(rgb, colormode = Color)
}

find_image_dirs <- function(opts) {
  if (file.exists(file.path(opts$input, opts$mask_name))) {
    return(normalizePath(opts$input))
  }
  
  all_files <- list.files(
    opts$input,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
  )
  
  mask_files <- all_files[basename(all_files) == opts$mask_name]
  
  if (length(mask_files) == 0) {
    stop("No ", opts$mask_name, " files found under: ", opts$input)
  }
  
  unique(dirname(mask_files))
}

process_image_dir <- function(image_dir, opts) {
  image_id <- safe_stem(basename(image_dir))
  out_dir <- file.path(opts$out, image_id)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  mask_file <- file.path(image_dir, opts$mask_name)
  gray_file <- file.path(image_dir, opts$gray_name)
  
  message("[02_centers] Processing: ", image_id)
  
  if (!file.exists(mask_file)) {
    stop("Missing mask file: ", mask_file)
  }
  
  mask <- readImage(mask_file)
  mask <- Image(imageData(mask) > 0.5, colormode = Grayscale)
  
  if (file.exists(gray_file)) {
    gray <- as_gray_image(readImage(gray_file))
  } else {
    gray <- normalize(mask)
  }
  
  labels <- bwlabel(mask)
  label_mat <- imageData(labels)
  
  outline_all <- make_outline(mask)
  
  blob_ids <- sort(unique(as.integer(label_mat)))
  blob_ids <- blob_ids[blob_ids > 0]
  
  raw_list <- list()
  raw_i <- 0
  
  if (length(blob_ids) > 0) {
    for (blob_id in blob_ids) {
      blob_area <- sum(label_mat == blob_id)
      
      if (blob_area < opts$min_blob_area) {
        next
      }
      
      blob_mask <- Image(label_mat == blob_id, colormode = Grayscale)
      blob_outline <- make_outline(blob_mask)
      
      idx <- which(imageData(blob_outline) > 0, arr.ind = TRUE)
      
      if (nrow(idx) < opts$min_arc_points) {
        next
      }
      
      boundary_pts <- data.frame(
        x = idx[, 1],
        y = idx[, 2]
      )
      
      raw_blob <- fit_candidates_for_blob(
        image_id = image_id,
        blob_id = blob_id,
        label_mat = label_mat,
        boundary_pts = boundary_pts,
        opts = opts
      )
      
      if (nrow(raw_blob) > 0) {
        raw_i <- raw_i + 1
        raw_list[[raw_i]] <- raw_blob
      }
    }
  }
  
  raw <- if (length(raw_list) > 0) {
    do.call(rbind, raw_list)
  } else {
    empty_raw_table()
  }
  
  merged <- merge_candidates(
    raw = raw,
    opts = opts,
    label_mat = label_mat
  )
  
  write.table(
    raw,
    file.path(out_dir, "candidate_centers.raw.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    merged,
    file.path(out_dir, "candidate_centers.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  overlay_all <- make_center_overlay(
    gray = gray,
    outline = outline_all,
    merged = merged,
    opts = opts,
    accepted_only = FALSE
  )
  
  overlay_accepted <- make_center_overlay(
    gray = gray,
    outline = outline_all,
    merged = merged,
    opts = opts,
    accepted_only = TRUE
  )
  
  writeImage(
    overlay_all,
    file.path(out_dir, "candidate_centers_overlay.png")
  )
  
  writeImage(
    overlay_accepted,
    file.path(out_dir, "accepted_centers_overlay.png")
  )
  
  summary <- data.frame(
    image_id = image_id,
    image_dir = image_dir,
    n_blobs = length(blob_ids),
    n_raw_candidates = nrow(raw),
    n_merged_candidates = nrow(merged),
    n_accepted_candidates = sum(merged$accepted, na.rm = TRUE),
    n_possible_overlap_breaks = sum(
      merged$possible_overlap_break & merged$accepted,
      na.rm = TRUE
    ),
    require_center_inside = opts$require_center_inside,
    center_inside_erode_radius = opts$center_inside_erode_radius,
    merge_method = opts$merge_method,
    merge_weight_power = opts$merge_weight_power,
    merge_min_score_fraction = opts$merge_min_score_fraction,
    require_weighted_center_inside = opts$require_weighted_center_inside,
    stringsAsFactors = FALSE
  )
  
  write.table(
    summary,
    file.path(out_dir, "center_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  list(
    raw = raw,
    merged = merged,
    summary = summary
  )
}

main <- function() {
  opts <- parse_args()
  
  dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
  
  image_dirs <- find_image_dirs(opts)
  
  message("[02_centers] Found ", length(image_dirs), " image result folder(s).")
  
  res <- lapply(image_dirs, process_image_dir, opts = opts)
  
  all_raw <- do.call(rbind, lapply(res, `[[`, "raw"))
  all_merged <- do.call(rbind, lapply(res, `[[`, "merged"))
  all_summary <- do.call(rbind, lapply(res, `[[`, "summary"))
  
  write.table(
    all_raw,
    file.path(opts$out, "all_candidate_centers.raw.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    all_merged,
    file.path(opts$out, "all_candidate_centers.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    all_summary,
    file.path(opts$out, "center_image_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("[02_centers] Done.")
  message("[02_centers] Results: ", opts$out)
}

main()