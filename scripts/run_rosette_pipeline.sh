#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_rosette_pipeline.sh
#
# Stable Bash wrapper for rosette image-analysis pipeline.
#
# Expected input structure:
#
# Option A:
# raw_images/
#   Sample_A/
#     image1.png
#     image2.png
#   Sample_B/
#     image3.png
#
# Option B:
# raw_images/
#   Sample_A/
#     rep1/
#       image1.png
#       image2.png
#     rep2/
#       image3.png
#   Sample_B/
#     rep1/
#       image4.png
#
# If a sample folder directly contains images, it is treated as
# one replicate named "replicate_01".
#
# Pipeline:
#   01_outline_layer.R
#   02_candidate_rosette_centers.R
#   03b_segment_cells_by_cluster_geometry.R
#   04_correct_and_summarize_rosettes.R
#   05_mask_confluence_summary.R
#
# Final project-level output now also includes enriched post-correction
# rosette summaries with:
#   - cells-per-rosette mean / median / sd / sem / min / max
#   - cell_counts_per_rosette, e.g. r1=11;r2=8;r3=14
#   - mask confluence merged from Step 05
# ============================================================


# ----------------------------
# Defaults
# ----------------------------

INPUT=""
OUT="results/rosette_pipeline"

ENV_NAME="rosette_pipeline"
ENV_YML="envs/rosette_pipeline.yml"

STEP1_SCRIPT="scripts/01_outline_layer.R"
STEP2_SCRIPT="scripts/02_candidate_rosette_centers.R"
STEP3_SCRIPT="scripts/03b_segment_cells_by_cluster_geometry.R"
CORRECT_SCRIPT="scripts/04_correct_and_summarize_rosettes.R"
CONFLUENCE_SCRIPT="scripts/05_mask_confluence_summary.R"
ENV_CHECK_SCRIPT="scripts/check_rosette_R_env.R"

RUN_STEP1="TRUE"
RUN_STEP2="TRUE"
RUN_STEP3="TRUE"
RUN_CORRECTION="TRUE"
RUN_PROJECT_SUMMARY="TRUE"
RUN_CONFLUENCE="TRUE"
RUN_ENRICHED_SUMMARY="TRUE"

OVERWRITE="FALSE"
DRY_RUN="FALSE"

IMAGE_EXT_REGEX="png|jpg|jpeg|tif|tiff"

MIN_CELLS_PER_ROSETTE=2
CELL_TABLE_NAME="all_cells.tsv"
PROJECT_SUMMARY_DIR_NAME="04_project_summary"

CONFLUENCE_DIR_NAME="05_confluence"
CONFLUENCE_MASK_NAME="mask.binary.png"
CONFLUENCE_OUTLINE_FOLDER="01_outline"
CONFLUENCE_FOREGROUND="bright"
CONFLUENCE_THRESHOLD="0.5"

ENRICHED_SUMMARY_DIR_NAME="06_enriched_summary"
CORRECTED_CELL_PATTERN="\\.cells\\.corrected\\.tsv$"

# ------------------------------------------------------------
# Step 1 defaults: outline layer
# ------------------------------------------------------------

STEP1_EXTRA="\
--pattern \"\\\\.(png|jpg|jpeg|tif|tiff)$\" \
--recursive FALSE \
--foreground bright \
--target_fraction 0.30 \
--method adaptive \
--window 41 \
--offset 0.03 \
--bg_sigma 90 \
--smooth_sigma 1.2 \
--open_radius 2 \
--close_radius 5 \
--fill_holes TRUE \
--min_area 100 \
--single_area_max 1000"

# ------------------------------------------------------------
# Step 2 defaults: candidate rosette centers
# ------------------------------------------------------------

STEP2_EXTRA="\
--mask_name mask.binary.png \
--gray_name gray.png \
--min_blob_area 675 \
--max_boundary_points 1000 \
--max_seed_points 500 \
--neighborhood_radius 20 \
--min_arc_points 20 \
--min_radius 15 \
--max_radius 400 \
--min_arc_angle 45 \
--max_fit_error 0.12 \
--ray_angles 26 \
--radial_bin_width_deg 110 \
--center_merge_dist 45 \
--min_total_score 0.60 \
--require_center_inside TRUE \
--center_inside_erode_radius 0 \
--merge_method weighted \
--merge_weight_power 2 \
--merge_min_score_fraction 0.50 \
--min_merged_support 2 \
--overlap_imbalance_threshold 0.50 \
--draw_rejected_top_n 0"

# ------------------------------------------------------------
# Step 3 defaults: cell segmentation by cluster geometry
# ------------------------------------------------------------

STEP3_EXTRA="\
--gray_name gray.png \
--cluster_mask_name mask.binary.png \
--rosette_table_name all_candidate_centers.tsv \
--candidate_center_source weighted \
--accepted_only TRUE \
--require_center_inside_blob FALSE \
--min_candidate_confidence 0.1 \
--mask_foreground bright \
--fill_cluster_holes FALSE \
--min_cluster_area_px 40 \
--use_membrane_cuts TRUE \
--membrane_bg_radius_px 7 \
--membrane_quantile 0.80 \
--membrane_min_score 0.020 \
--membrane_open_px 0 \
--membrane_dilate_px 0 \
--min_cell_area_px 35 \
--max_cell_area_px 1200 \
--min_cell_radius_px 1.5 \
--max_cell_radius_px 80 \
--seed_smooth_sigma 0.6 \
--seed_min_distance_px 3 \
--watershed_tolerance 0.6 \
--min_solidity 0.25 \
--min_circularity 0.15 \
--max_edge_contact_fraction 0.5 \
--assignment_max_norm_distance 1.30 \
--assignment_radius_weight 1.0 \
--assignment_overlap_weight 0.5 \
--assignment_min_overlap 0.25 \
--draw_cell_ids TRUE \
--draw_rejected TRUE \
--label_mode both \
--overlay_scale 5 \
--boundary_line_width_px 0 \
--debug TRUE \
--debug_parser TRUE \
--keep_going FALSE"


# ----------------------------
# Usage
# ----------------------------

usage() {
cat <<EOF

Usage:

  bash scripts/run_rosette_pipeline.sh \\
    --input raw_images \\
    --out results/rosette_pipeline_v01

Required:
  --input DIR
      Folder containing sample folders.

Optional:
  --out DIR
      Output folder.
      Default: results/rosette_pipeline

  --env-name NAME
      Conda environment name.
      Default: rosette_pipeline

  --env-yml FILE
      Conda environment YAML.
      Default: envs/rosette_pipeline.yml

  --step1-script FILE
      Step 1 R script.
      Default: scripts/01_outline_layer.R

  --step2-script FILE
      Step 2 R script.
      Default: scripts/02_candidate_rosette_centers.R

  --step3-script FILE
      Step 3 R script.
      Default: scripts/03b_segment_cells_by_cluster_geometry.R

  --correct-script FILE
      Correction / summary R script.
      Default: scripts/04_correct_and_summarize_rosettes.R

  --confluence-script FILE
      Step 05 confluence R script.
      Default: scripts/05_mask_confluence_summary.R

  --run-step1 TRUE/FALSE
      Default: TRUE

  --run-step2 TRUE/FALSE
      Default: TRUE

  --run-step3 TRUE/FALSE
      Default: TRUE

  --run-correction TRUE/FALSE
      Run Step 4 per replicate after Step 3.
      Default: TRUE

  --run-project-summary TRUE/FALSE
      Run Step 4 once across the whole output project tree.
      Default: TRUE

  --run-confluence TRUE/FALSE
      Run Step 05 mask confluence after Step 4.
      Default: TRUE

  --run-enriched-summary TRUE/FALSE
      Run final enriched rosette summary after Step 05.
      Default: TRUE

  --cell-table-name FILE
      Cell table filename to use for final project-level Step 4 search.
      Default: all_cells.tsv

  --project-summary-dir-name NAME
      Name of final project-level summary folder inside --out.
      Default: 04_project_summary

  --confluence-dir-name NAME
      Name of Step 05 confluence folder inside --out.
      Default: 05_confluence

  --confluence-mask-name NAME
      Binary mask filename for confluence.
      Default: mask.binary.png

  --confluence-outline-folder NAME
      Only masks inside this folder are used.
      Default: 01_outline

  --confluence-foreground bright/dark/auto
      Foreground polarity for mask confluence.
      Default: bright

  --confluence-threshold FLOAT
      Binary mask threshold.
      Default: 0.5

  --enriched-summary-dir-name NAME
      Name of final enriched summary folder inside --out.
      Default: 06_enriched_summary

  --corrected-cell-pattern REGEX
      Regex for corrected per-cell TSVs used by enriched summary.
      Default: \\.cells\\.corrected\\.tsv$

  --min-cells-per-rosette INT
      Minimum cells required for a corrected rosette.
      Rosettes with fewer cells are demoted to single cells.
      Default: 2

  --overwrite TRUE/FALSE
      Re-run steps even if output folders already exist.
      Default: FALSE

  --dry-run TRUE/FALSE
      Print commands without executing.
      Default: FALSE

Examples:

  bash scripts/run_rosette_pipeline.sh \\
    --input raw_images \\
    --out results/rosette_pipeline_v01

  bash scripts/run_rosette_pipeline.sh \\
    --input raw_images \\
    --out results/rosette_pipeline_v01 \\
    --overwrite TRUE

EOF
}


# ----------------------------
# Argument parsing
# ----------------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT="$2"
      shift 2
      ;;

    --out)
      OUT="$2"
      shift 2
      ;;

    --env-name)
      ENV_NAME="$2"
      shift 2
      ;;

    --env-yml)
      ENV_YML="$2"
      shift 2
      ;;

    --step1-script)
      STEP1_SCRIPT="$2"
      shift 2
      ;;

    --step2-script)
      STEP2_SCRIPT="$2"
      shift 2
      ;;

    --step3-script)
      STEP3_SCRIPT="$2"
      shift 2
      ;;

    --correct-script)
      CORRECT_SCRIPT="$2"
      shift 2
      ;;

    --confluence-script)
      CONFLUENCE_SCRIPT="$2"
      shift 2
      ;;

    --run-step1)
      RUN_STEP1="$2"
      shift 2
      ;;

    --run-step2)
      RUN_STEP2="$2"
      shift 2
      ;;

    --run-step3)
      RUN_STEP3="$2"
      shift 2
      ;;

    --run-correction)
      RUN_CORRECTION="$2"
      shift 2
      ;;

    --run-project-summary)
      RUN_PROJECT_SUMMARY="$2"
      shift 2
      ;;

    --run-confluence)
      RUN_CONFLUENCE="$2"
      shift 2
      ;;

    --run-enriched-summary)
      RUN_ENRICHED_SUMMARY="$2"
      shift 2
      ;;

    --cell-table-name)
      CELL_TABLE_NAME="$2"
      shift 2
      ;;

    --project-summary-dir-name)
      PROJECT_SUMMARY_DIR_NAME="$2"
      shift 2
      ;;

    --confluence-dir-name)
      CONFLUENCE_DIR_NAME="$2"
      shift 2
      ;;

    --confluence-mask-name)
      CONFLUENCE_MASK_NAME="$2"
      shift 2
      ;;

    --confluence-outline-folder)
      CONFLUENCE_OUTLINE_FOLDER="$2"
      shift 2
      ;;

    --confluence-foreground)
      CONFLUENCE_FOREGROUND="$2"
      shift 2
      ;;

    --confluence-threshold)
      CONFLUENCE_THRESHOLD="$2"
      shift 2
      ;;

    --enriched-summary-dir-name)
      ENRICHED_SUMMARY_DIR_NAME="$2"
      shift 2
      ;;

    --corrected-cell-pattern)
      CORRECTED_CELL_PATTERN="$2"
      shift 2
      ;;

    --min-cells-per-rosette)
      MIN_CELLS_PER_ROSETTE="$2"
      shift 2
      ;;

    --overwrite)
      OVERWRITE="$2"
      shift 2
      ;;

    --dry-run)
      DRY_RUN="$2"
      shift 2
      ;;

    -h|--help)
      usage
      exit 0
      ;;

    *)
      echo "[wrapper] ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT" ]]; then
  echo "[wrapper] ERROR: --input is required." >&2
  usage
  exit 1
fi


# ----------------------------
# Helper functions
# ----------------------------

timestamp() {
  date +"%Y-%m-%d %H:%M:%S"
}

log() {
  echo "[$(timestamp)] $*"
}

run_cmd() {
  echo
  echo "[cmd] $*"

  if [[ "$DRY_RUN" != "TRUE" ]]; then
    eval "$@"
  fi
}

script_supports_option() {
  local script="$1"
  local opt_name="$2"

  if [[ "$DRY_RUN" == "TRUE" ]]; then
    return 0
  fi

  Rscript "$script" --help 2>&1 | grep -q -- "$opt_name"
}

append_arg_if_supported() {
  local script="$1"
  local opt_name="$2"
  local opt_value="$3"
  local -n arg_ref="$4"

  if script_supports_option "$script" "$opt_name"; then
    arg_ref+=" $opt_name '$opt_value'"
  else
    log "[wrapper] NOTE: $script does not advertise $opt_name; not passing it."
  fi
}

normalize_step4_outputs() {
  local correct_out="$1"

  local names=(image_summary replicate_summary sample_summary global_summary project_summary rosette_summary)
  local n

  for n in "${names[@]}"; do
    if [[ -f "$correct_out/${n}.tsv" && ! -f "$correct_out/${n}.corrected.tsv" ]]; then
      cp "$correct_out/${n}.tsv" "$correct_out/${n}.corrected.tsv"
    fi
  done
}

run_step4_summary() {
  local input_dir="$1"
  local out_dir="$2"
  local sample_id="$3"
  local replicate_id="$4"
  local project_level="$5"
  local sample_depth="$6"
  local replicate_depth="$7"
  local image_depth="$8"

  mkdir -p "$out_dir"

  local args="--input '$input_dir' --out '$out_dir' --min_cells_per_rosette '$MIN_CELLS_PER_ROSETTE'"

  append_arg_if_supported "$CORRECT_SCRIPT" "--accepted_only" "TRUE" args

  if [[ "$project_level" == "TRUE" ]]; then
    append_arg_if_supported "$CORRECT_SCRIPT" "--project_level" "TRUE" args
    append_arg_if_supported "$CORRECT_SCRIPT" "--sample_depth" "$sample_depth" args
    append_arg_if_supported "$CORRECT_SCRIPT" "--replicate_depth" "$replicate_depth" args
    append_arg_if_supported "$CORRECT_SCRIPT" "--image_depth" "$image_depth" args

    if [[ -n "$CELL_TABLE_NAME" ]]; then
      append_arg_if_supported "$CORRECT_SCRIPT" "--cell_table_name" "$CELL_TABLE_NAME" args
    fi
  else
    append_arg_if_supported "$CORRECT_SCRIPT" "--sample_id" "$sample_id" args
    append_arg_if_supported "$CORRECT_SCRIPT" "--replicate_id" "$replicate_id" args
  fi

  run_cmd "Rscript '$CORRECT_SCRIPT' $args"

  if [[ "$DRY_RUN" != "TRUE" ]]; then
    normalize_step4_outputs "$out_dir"
  fi
}

copy_project_summary_aliases() {
  local project_summary_dir="$1"

  local pairs=(
    "image_summary.tsv:image_summary.corrected.tsv"
    "replicate_summary.tsv:replicate_summary.corrected.tsv"
    "sample_summary.tsv:sample_summary.corrected.tsv"
    "global_summary.tsv:global_summary.corrected.tsv"
    "project_summary.tsv:project_summary.corrected.tsv"
    "rosette_summary.tsv:rosette_summary.corrected.tsv"
  )

  local pair src dst
  for pair in "${pairs[@]}"; do
    src="${pair%%:*}"
    dst="${pair##*:}"

    if [[ -f "$project_summary_dir/$src" ]]; then
      cp "$project_summary_dir/$src" "$OUT/$dst"
    elif [[ -f "$project_summary_dir/$dst" ]]; then
      cp "$project_summary_dir/$dst" "$OUT/$dst"
    fi
  done
}

copy_enriched_summary_aliases() {
  local enriched_dir="$1"

  if [[ -f "$enriched_dir/rosette_summary.enriched.per_image.tsv" ]]; then
    cp "$enriched_dir/rosette_summary.enriched.per_image.tsv" "$OUT/rosette_summary.enriched.per_image.tsv"
  fi

  if [[ -f "$enriched_dir/rosette_summary.enriched.by_sample_replicate.tsv" ]]; then
    cp "$enriched_dir/rosette_summary.enriched.by_sample_replicate.tsv" "$OUT/rosette_summary.enriched.by_sample_replicate.tsv"
  fi

  if [[ -f "$enriched_dir/rosette_counts.enriched.per_rosette.tsv" ]]; then
    cp "$enriched_dir/rosette_counts.enriched.per_rosette.tsv" "$OUT/rosette_counts.enriched.per_rosette.tsv"
  fi
}

has_images_directly() {
  local d="$1"

  find "$d" -maxdepth 1 -type f | grep -Eiq "\.(${IMAGE_EXT_REGEX})$"
}

safe_name() {
  basename "$1" | sed 's/[^A-Za-z0-9_.-]/_/g'
}

require_file() {
  local f="$1"

  if [[ ! -f "$f" ]]; then
    echo "[wrapper] ERROR: required file not found: $f" >&2
    exit 1
  fi
}

require_dir() {
  local d="$1"

  if [[ ! -d "$d" ]]; then
    echo "[wrapper] ERROR: required directory not found: $d" >&2
    exit 1
  fi
}

find_centers_table() {
  local step2_out="$1"
  local f=""

  f="$(find "$step2_out" -name "all_candidate_centers.tsv" -type f | head -n 1 || true)"

  if [[ -z "$f" ]]; then
    f="$(find "$step2_out" -name "*candidate*center*.tsv" -type f | head -n 1 || true)"
  fi

  if [[ -z "$f" ]]; then
    f="$(find "$step2_out" -name "*.tsv" -type f | grep -Ei "center|candidate|rosette" | head -n 1 || true)"
  fi

  echo "$f"
}

stage_step3_input() {
  local step1_out="$1"
  local step2_out="$2"
  local staged_out="$3"

  mkdir -p "$staged_out"

  echo "[$(timestamp)] Staging Step 3 input:"
  echo "[$(timestamp)]   Step 1 source: $step1_out"
  echo "[$(timestamp)]   Step 2 source: $step2_out"
  echo "[$(timestamp)]   Staged dir:   $staged_out"

  if [[ ! -d "$step1_out" ]]; then
    echo "[wrapper] ERROR: Step 1 output directory not found for staging: $step1_out" >&2
    exit 1
  fi

  if [[ ! -d "$step2_out" ]]; then
    echo "[wrapper] ERROR: Step 2 output directory not found for staging: $step2_out" >&2
    exit 1
  fi

  if command -v rsync >/dev/null 2>&1; then
    rsync -a "$step1_out"/ "$staged_out"/
  else
    cp -a "$step1_out"/. "$staged_out"/
  fi

  local centers_table=""
  centers_table="$(find_centers_table "$step2_out")"

  if [[ -z "$centers_table" ]]; then
    echo "[wrapper] ERROR: could not find candidate center table in: $step2_out" >&2
    exit 1
  fi

  echo "[$(timestamp)]   Centers table: $centers_table"

  cp "$centers_table" "$staged_out/all_candidate_centers.tsv"

  local n_staged=0

  for image_dir in "$staged_out"/*; do
    [[ -d "$image_dir" ]] || continue

    if [[ -f "$image_dir/gray.png" && -f "$image_dir/mask.binary.png" ]]; then
      cp "$centers_table" "$image_dir/all_candidate_centers.tsv"
      n_staged=$((n_staged + 1))
      echo "[$(timestamp)]   staged centers into: $image_dir/all_candidate_centers.tsv"
    fi
  done

  if [[ "$n_staged" -eq 0 ]]; then
    echo "[wrapper] ERROR: staged input contains no direct child folders with gray.png and mask.binary.png" >&2
    echo "[wrapper]        staged_out: $staged_out" >&2
    echo "[wrapper]        direct child folders/files:" >&2
    find "$staged_out" -maxdepth 2 -type f | head -50 >&2 || true
    exit 1
  fi

  echo "[$(timestamp)]   Staged centers table into $n_staged image folder(s)."
}

run_step5_confluence() {
  local input_root="$1"
  local out_dir="$2"

  mkdir -p "$out_dir"

  run_cmd "Rscript '$CONFLUENCE_SCRIPT' \
    --input '$input_root' \
    --out '$out_dir' \
    --mask_name '$CONFLUENCE_MASK_NAME' \
    --outline_folder '$CONFLUENCE_OUTLINE_FOLDER' \
    --foreground '$CONFLUENCE_FOREGROUND' \
    --threshold '$CONFLUENCE_THRESHOLD'"
}

run_enriched_summary() {
  local input_root="$1"
  local out_dir="$2"
  local confluence_per_image="$3"
  local confluence_by_rep="$4"

  mkdir -p "$out_dir"

  if [[ "$DRY_RUN" == "TRUE" ]]; then
    echo "[dry-run] run enriched rosette summary into '$out_dir'"
    return 0
  fi

  Rscript --vanilla - "$input_root" "$out_dir" "$CORRECTED_CELL_PATTERN" "$MIN_CELLS_PER_ROSETTE" "$confluence_per_image" "$confluence_by_rep" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

input_root <- args[[1]]
out_dir <- args[[2]]
corrected_pattern <- args[[3]]
min_cells_per_rosette <- as.integer(args[[4]])
confluence_per_image_path <- args[[5]]
confluence_by_rep_path <- args[[6]]

SCRIPT_VERSION <- "embedded_enriched_rosette_summary_v01"

message("[", SCRIPT_VERSION, "] Input root: ", input_root)
message("[", SCRIPT_VERSION, "] Output dir: ", out_dir)
message("[", SCRIPT_VERSION, "] Corrected pattern: ", corrected_pattern)
message("[", SCRIPT_VERSION, "] Min cells per rosette: ", min_cells_per_rosette)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

escape_regex <- function(x) {
  gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", x)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

bool_like_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  xl <- tolower(as.character(x))
  xl %in% c("true", "t", "1", "yes", "y", "accepted", "pass", "passed")
}

first_existing_col <- function(dt, candidates) {
  hits <- candidates[candidates %in% names(dt)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

infer_ids_from_path <- function(f, input_root) {
  f_norm <- normalizePath(f, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(input_root, winslash = "/", mustWork = FALSE)

  rel <- sub(paste0("^", escape_regex(root_norm), "/?"), "", f_norm)
  parts <- unlist(strsplit(rel, "/", fixed = TRUE))
  parts <- parts[nzchar(parts)]

  sample_id <- NA_character_
  replicate_id <- NA_character_

  rep_idx <- which(grepl("^replicate[_-]?[0-9]+$", parts, ignore.case = TRUE))
  if (length(rep_idx) > 0) {
    rep_idx <- rep_idx[[1]]
    replicate_id <- parts[[rep_idx]]
    if (rep_idx > 1) sample_id <- parts[[rep_idx - 1]]
  }

  if (is.na(sample_id) || !nzchar(sample_id)) {
    corr_idx <- which(parts == "04_corrected")
    if (length(corr_idx) > 0) {
      corr_idx <- corr_idx[[1]]
      if (corr_idx > 2) {
        sample_id <- parts[[corr_idx - 2]]
      } else if (corr_idx > 1) {
        sample_id <- parts[[corr_idx - 1]]
      }
    }
  }

  if (is.na(replicate_id) || !nzchar(replicate_id)) {
    replicate_id <- "replicate_01"
  }

  if (is.na(sample_id) || !nzchar(sample_id)) {
    if (length(parts) >= 3) {
      sample_id <- parts[[length(parts) - 2]]
    } else if (length(parts) >= 2) {
      sample_id <- parts[[length(parts) - 1]]
    } else {
      sample_id <- "sample_unknown"
    }
  }

  image_id <- basename(f_norm)
  image_id <- sub("\\.cells\\.corrected\\.tsv$", "", image_id)
  image_id <- sub("\\.corrected\\.tsv$", "", image_id)
  image_id <- sub("\\.tsv$", "", image_id)

  list(
    sample_id = sample_id,
    replicate_id = replicate_id,
    image_id = image_id,
    relative_path = rel,
    source_file = f_norm
  )
}

standardize_file <- function(f) {
  message("[", SCRIPT_VERSION, "] Reading: ", f)
  dt <- fread(f, sep = "\t", showProgress = FALSE)

  ids <- infer_ids_from_path(f, input_root)

  if (!"sample_id" %in% names(dt)) {
    dt[, sample_id := ids$sample_id]
  } else {
    dt[is.na(sample_id) | sample_id == "", sample_id := ids$sample_id]
  }

  if (!"replicate_id" %in% names(dt)) {
    dt[, replicate_id := ids$replicate_id]
  } else {
    dt[is.na(replicate_id) | replicate_id == "", replicate_id := ids$replicate_id]
  }

  if (!"image_id" %in% names(dt)) {
    dt[, image_id := ids$image_id]
  } else {
    dt[is.na(image_id) | image_id == "", image_id := ids$image_id]
  }

  dt[, source_file := ids$source_file]
  dt[, relative_path := ids$relative_path]

  dt
}

corrected_files <- list.files(
  input_root,
  pattern = corrected_pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (length(corrected_files) > 0) {
  corrected_files <- corrected_files[
    grepl("(^|/)04_corrected(/|$)", normalizePath(corrected_files, winslash = "/", mustWork = FALSE))
  ]
}

if (length(corrected_files) == 0) {
  stop(
    "[", SCRIPT_VERSION, "] ERROR: no corrected cell files found under ",
    input_root, " with pattern ", corrected_pattern,
    call. = FALSE
  )
}

message("[", SCRIPT_VERSION, "] Found corrected files: ", length(corrected_files))

cells <- rbindlist(lapply(corrected_files, standardize_file), fill = TRUE)
message("[", SCRIPT_VERSION, "] Combined corrected rows: ", nrow(cells))

accepted_col <- first_existing_col(
  cells,
  c("accepted_cell_for_summary", "accepted_cell", "cell_accepted", "accepted", "valid_cell", "keep_cell")
)

if (!is.na(accepted_col)) {
  before <- nrow(cells)
  cells <- cells[bool_like_true(get(accepted_col))]
  message("[", SCRIPT_VERSION, "] Accepted-cell filter via ", accepted_col, ": ", nrow(cells), " / ", before)
}

valid_rosette_col <- first_existing_col(
  cells,
  c("valid_rosette", "valid_rosette.x", "valid_rosette.y", "accepted_rosette", "rosette_accepted", "keep_rosette")
)

if (!is.na(valid_rosette_col)) {
  before <- nrow(cells)
  cells <- cells[bool_like_true(get(valid_rosette_col))]
  message("[", SCRIPT_VERSION, "] Valid-rosette filter via ", valid_rosette_col, ": ", nrow(cells), " / ", before)
}

rosette_col <- first_existing_col(
  cells,
  c("rosette_id", "corrected_rosette_id", "assigned_rosette_id", "rosette", "cluster_id")
)

if (is.na(rosette_col)) {
  stop("[", SCRIPT_VERSION, "] ERROR: could not identify rosette ID column.", call. = FALSE)
}

cell_id_col <- first_existing_col(
  cells,
  c("cell_id", "corrected_cell_id", "object_id", "label", "cell_label", "seg_cell_id")
)

n_cells_col <- first_existing_col(
  cells,
  c("n_cells", "n_cell", "cell_count", "rosette_n_cells")
)

timepoint_col <- first_existing_col(cells, c("timepoint", "time_point", "tp"))

cells[, rosette_id_internal := as.character(get(rosette_col))]

bad_rosette_values <- c(
  "", "NA", "NaN", "nan", "NULL", "null", "none", "None",
  "single", "singleton", "unassigned", "background", "0", "-1"
)

before <- nrow(cells)
cells <- cells[
  !is.na(rosette_id_internal) &
    !(rosette_id_internal %in% bad_rosette_values)
]
message("[", SCRIPT_VERSION, "] Rosette-assigned rows retained: ", nrow(cells), " / ", before)

id_cols <- c("sample_id", "replicate_id", "image_id")
if (!is.na(timepoint_col)) {
  if (timepoint_col != "timepoint") cells[, timepoint := get(timepoint_col)]
  id_cols <- c("sample_id", "replicate_id", "timepoint", "image_id")
}

if (!is.na(cell_id_col)) {
  rosette_counts <- cells[
    ,
    .(
      n_cells = uniqueN(get(cell_id_col)),
      n_rows_used = .N
    ),
    by = c(id_cols, "rosette_id_internal")
  ]
} else if (!is.na(n_cells_col)) {
  rosette_counts <- cells[
    ,
    .(
      n_cells = as.integer(round(max(safe_numeric(get(n_cells_col)), na.rm = TRUE))),
      n_rows_used = .N
    ),
    by = c(id_cols, "rosette_id_internal")
  ]
} else {
  rosette_counts <- cells[
    ,
    .(
      n_cells = .N,
      n_rows_used = .N
    ),
    by = c(id_cols, "rosette_id_internal")
  ]
}

setnames(rosette_counts, "rosette_id_internal", "rosette_id")
rosette_counts[, n_cells := as.integer(n_cells)]

before <- nrow(rosette_counts)
rosette_counts <- rosette_counts[n_cells >= min_cells_per_rosette]
message(
  "[", SCRIPT_VERSION, "] Rosettes passing min cell count: ",
  nrow(rosette_counts), " / ", before
)

rosette_counts[, rosette_sort := safe_numeric(rosette_id)]
rosette_counts[is.na(rosette_sort), rosette_sort := Inf]
setorderv(rosette_counts, c(id_cols, "rosette_sort", "rosette_id"))

rosette_counts[, rosette_label := paste0("r", rosette_id, "=", n_cells)]

per_image <- rosette_counts[
  ,
  .(
    n_rosettes = .N,
    total_cells_in_rosettes = sum(n_cells, na.rm = TRUE),
    mean_cells_per_rosette = mean(n_cells, na.rm = TRUE),
    median_cells_per_rosette = median(n_cells, na.rm = TRUE),
    sd_cells_per_rosette = if (.N > 1) stats::sd(n_cells, na.rm = TRUE) else NA_real_,
    sem_cells_per_rosette = sem(n_cells),
    min_cells_per_rosette = min(n_cells, na.rm = TRUE),
    max_cells_per_rosette = max(n_cells, na.rm = TRUE),
    cell_counts_per_rosette = paste(rosette_label, collapse = ";")
  ),
  by = id_cols
]

setorderv(per_image, id_cols)

sample_group_cols <- c("sample_id", "replicate_id")
if ("timepoint" %in% names(per_image)) {
  sample_group_cols <- c("sample_id", "replicate_id", "timepoint")
}

by_rep_image_stats <- per_image[
  ,
  .(
    n_images = .N,
    total_rosettes = sum(n_rosettes, na.rm = TRUE),
    mean_rosettes_per_image = mean(n_rosettes, na.rm = TRUE),
    median_rosettes_per_image = median(n_rosettes, na.rm = TRUE),
    sd_rosettes_per_image = if (.N > 1) stats::sd(n_rosettes, na.rm = TRUE) else NA_real_,
    sem_rosettes_per_image = sem(n_rosettes),
    min_rosettes_per_image = min(n_rosettes, na.rm = TRUE),
    max_rosettes_per_image = max(n_rosettes, na.rm = TRUE),
    total_cells_in_rosettes = sum(total_cells_in_rosettes, na.rm = TRUE),
    mean_image_mean_cells_per_rosette = mean(mean_cells_per_rosette, na.rm = TRUE),
    sd_image_mean_cells_per_rosette = if (.N > 1) stats::sd(mean_cells_per_rosette, na.rm = TRUE) else NA_real_,
    sem_image_mean_cells_per_rosette = sem(mean_cells_per_rosette)
  ),
  by = sample_group_cols
]

by_rep_rosette_stats <- rosette_counts[
  ,
  .(
    n_rosettes_pooled = .N,
    mean_cells_per_rosette_pooled = mean(n_cells, na.rm = TRUE),
    median_cells_per_rosette_pooled = median(n_cells, na.rm = TRUE),
    sd_cells_per_rosette_pooled = if (.N > 1) stats::sd(n_cells, na.rm = TRUE) else NA_real_,
    sem_cells_per_rosette_pooled = sem(n_cells),
    min_cells_per_rosette_pooled = min(n_cells, na.rm = TRUE),
    max_cells_per_rosette_pooled = max(n_cells, na.rm = TRUE)
  ),
  by = sample_group_cols
]

by_rep <- merge(
  by_rep_image_stats,
  by_rep_rosette_stats,
  by = sample_group_cols,
  all = TRUE
)

if (file.exists(confluence_per_image_path)) {
  conf_img <- fread(confluence_per_image_path, sep = "\t", showProgress = FALSE)

  keep <- intersect(
    c(
      "sample_id",
      "replicate_id",
      "timepoint",
      "image_id",
      "image_width_px",
      "image_height_px",
      "total_pixels",
      "foreground_pixels",
      "confluence_fraction",
      "confluence_percent"
    ),
    names(conf_img)
  )

  conf_img <- unique(conf_img[, ..keep])

  merge_keys <- intersect(
    c("sample_id", "replicate_id", "timepoint", "image_id"),
    intersect(names(per_image), names(conf_img))
  )

  message("[", SCRIPT_VERSION, "] Merging per-image confluence on: ", paste(merge_keys, collapse = ", "))
  per_image <- merge(per_image, conf_img, by = merge_keys, all.x = TRUE)
} else {
  warning("[", SCRIPT_VERSION, "] Confluence per-image file not found: ", confluence_per_image_path)
}

if (file.exists(confluence_by_rep_path)) {
  conf_rep <- fread(confluence_by_rep_path, sep = "\t", showProgress = FALSE)

  merge_keys <- intersect(
    c("sample_id", "replicate_id", "timepoint"),
    intersect(names(by_rep), names(conf_rep))
  )

  non_keys <- setdiff(names(conf_rep), merge_keys)
  setnames(conf_rep, non_keys, paste0("mask_", non_keys))

  message("[", SCRIPT_VERSION, "] Merging replicate confluence on: ", paste(merge_keys, collapse = ", "))
  by_rep <- merge(by_rep, conf_rep, by = merge_keys, all.x = TRUE)
} else {
  warning("[", SCRIPT_VERSION, "] Confluence replicate file not found: ", confluence_by_rep_path)
}

setorderv(per_image, intersect(c("sample_id", "replicate_id", "timepoint", "image_id"), names(per_image)))
setorderv(by_rep, intersect(c("sample_id", "replicate_id", "timepoint"), names(by_rep)))

per_rosette_out <- file.path(out_dir, "rosette_counts.enriched.per_rosette.tsv")
per_image_out <- file.path(out_dir, "rosette_summary.enriched.per_image.tsv")
by_rep_out <- file.path(out_dir, "rosette_summary.enriched.by_sample_replicate.tsv")

fwrite(rosette_counts, per_rosette_out, sep = "\t")
fwrite(per_image, per_image_out, sep = "\t")
fwrite(by_rep, by_rep_out, sep = "\t")

message("[", SCRIPT_VERSION, "] Wrote: ", per_rosette_out)
message("[", SCRIPT_VERSION, "] Wrote: ", per_image_out)
message("[", SCRIPT_VERSION, "] Wrote: ", by_rep_out)
message("[", SCRIPT_VERSION, "] Done.")
RSCRIPT
}


# ----------------------------
# Initial checks
# ----------------------------

require_dir "$INPUT"
require_file "$ENV_YML"
require_file "$STEP1_SCRIPT"
require_file "$STEP2_SCRIPT"
require_file "$STEP3_SCRIPT"
require_file "$ENV_CHECK_SCRIPT"

if [[ "$RUN_CORRECTION" == "TRUE" ]]; then
  require_file "$CORRECT_SCRIPT"
fi

if [[ "$RUN_CONFLUENCE" == "TRUE" ]]; then
  require_file "$CONFLUENCE_SCRIPT"
fi

mkdir -p "$OUT"
mkdir -p "$OUT/logs"

MASTER_LOG="$OUT/logs/wrapper.master.log"

exec > >(tee -a "$MASTER_LOG") 2>&1

log "[wrapper] Starting rosette pipeline."
log "[wrapper] Input: $INPUT"
log "[wrapper] Output: $OUT"
log "[wrapper] Environment name: $ENV_NAME"
log "[wrapper] Environment YAML: $ENV_YML"
log "[wrapper] Min cells per corrected rosette: $MIN_CELLS_PER_ROSETTE"
log "[wrapper] Run project-level summary: $RUN_PROJECT_SUMMARY"
log "[wrapper] Run Step 05 confluence: $RUN_CONFLUENCE"
log "[wrapper] Run enriched summary: $RUN_ENRICHED_SUMMARY"
log "[wrapper] Project summary dir name: $PROJECT_SUMMARY_DIR_NAME"
log "[wrapper] Confluence dir name: $CONFLUENCE_DIR_NAME"
log "[wrapper] Enriched summary dir name: $ENRICHED_SUMMARY_DIR_NAME"
log "[wrapper] Project summary cell table name: $CELL_TABLE_NAME"
log "[wrapper] Dry run: $DRY_RUN"
log "[wrapper] Overwrite: $OVERWRITE"


# ----------------------------
# Conda / mamba setup
# ----------------------------

if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_EXE="conda"
else
  echo "[wrapper] ERROR: neither conda nor mamba found on PATH." >&2
  exit 1
fi

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "$CONDA_BASE/etc/profile.d/conda.sh"
else
  echo "[wrapper] ERROR: conda is required for environment activation." >&2
  exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  log "[wrapper] Conda env exists: $ENV_NAME"
  log "[wrapper] Updating env from: $ENV_YML"
  run_cmd "$CONDA_EXE env update -n '$ENV_NAME' -f '$ENV_YML' --prune"
else
  log "[wrapper] Creating conda env: $ENV_NAME"
  run_cmd "$CONDA_EXE env create -n '$ENV_NAME' -f '$ENV_YML'"
fi

log "[wrapper] Activating env: $ENV_NAME"

if [[ "$DRY_RUN" != "TRUE" ]]; then
  conda activate "$ENV_NAME"
else
  echo "[dry-run] conda activate '$ENV_NAME'"
fi

log "[wrapper] Checking R environment."
run_cmd "Rscript '$ENV_CHECK_SCRIPT'"


# ----------------------------
# Write run manifest header
# ----------------------------

RUN_TABLE="$OUT/pipeline_run_table.tsv"

printf "sample_id\treplicate_id\tinput_dir\tstep1_out\tstep2_out\tstep3_input\tstep3_out\tstatus\n" > "$RUN_TABLE"


# ----------------------------
# Main loop
# ----------------------------

shopt -s nullglob

log "[wrapper] Sample folders discovered under input:"
find "$INPUT" -mindepth 1 -maxdepth 1 -type d -print | sort || true

for SAMPLE_DIR in "$INPUT"/*; do
  [[ -d "$SAMPLE_DIR" ]] || continue

  SAMPLE_ID="$(safe_name "$SAMPLE_DIR")"

  log "============================================================"
  log "[wrapper] Sample: $SAMPLE_ID"
  log "============================================================"

  REPLICATE_DIRS=()

  if has_images_directly "$SAMPLE_DIR"; then
    REPLICATE_DIRS+=("$SAMPLE_DIR")
  fi

  for REP_DIR in "$SAMPLE_DIR"/*; do
    [[ -d "$REP_DIR" ]] || continue

    if has_images_directly "$REP_DIR"; then
      REPLICATE_DIRS+=("$REP_DIR")
    fi
  done

  if [[ "${#REPLICATE_DIRS[@]}" -eq 0 ]]; then
    log "[wrapper] WARNING: no image-containing replicate folders found for sample: $SAMPLE_ID"
    continue
  fi

  REP_INDEX=0

  for REP_DIR in "${REPLICATE_DIRS[@]}"; do
    REP_INDEX=$((REP_INDEX + 1))

    if [[ "$REP_DIR" == "$SAMPLE_DIR" ]]; then
      REPLICATE_ID="replicate_$(printf "%02d" "$REP_INDEX")"
    else
      REPLICATE_ID="$(safe_name "$REP_DIR")"
    fi

    RUN_ROOT="$OUT/$SAMPLE_ID/$REPLICATE_ID"
    STEP1_OUT="$RUN_ROOT/01_outline"
    STEP2_OUT="$RUN_ROOT/02_centers"
    STEP3_INPUT="$RUN_ROOT/03_input_staged"
    STEP3_OUT="$RUN_ROOT/03_cells"
    CORRECT_OUT="$RUN_ROOT/04_corrected"

    REP_LOG="$OUT/logs/${SAMPLE_ID}__${REPLICATE_ID}.log"

    mkdir -p "$RUN_ROOT"

    log "[wrapper] Replicate: $REPLICATE_ID"
    log "[wrapper] Input: $REP_DIR"
    log "[wrapper] Log: $REP_LOG"

    {
      echo "[$(timestamp)] Sample: $SAMPLE_ID"
      echo "[$(timestamp)] Replicate: $REPLICATE_ID"
      echo "[$(timestamp)] Input: $REP_DIR"
      echo "[$(timestamp)] Run root: $RUN_ROOT"

      # ----------------------------
      # Step 1: outline layer
      # ----------------------------

      if [[ "$RUN_STEP1" == "TRUE" ]]; then
        if [[ -d "$STEP1_OUT" && "$OVERWRITE" != "TRUE" ]]; then
          echo "[$(timestamp)] Step 1 exists; skipping: $STEP1_OUT"
        else
          mkdir -p "$STEP1_OUT"
          run_cmd "Rscript '$STEP1_SCRIPT' --input '$REP_DIR' --out '$STEP1_OUT' $STEP1_EXTRA"
        fi
      fi

      # ----------------------------
      # Step 2: candidate centers
      # ----------------------------

      if [[ "$RUN_STEP2" == "TRUE" ]]; then
        if [[ -d "$STEP2_OUT" && "$OVERWRITE" != "TRUE" ]]; then
          echo "[$(timestamp)] Step 2 exists; skipping: $STEP2_OUT"
        else
          mkdir -p "$STEP2_OUT"
          run_cmd "Rscript '$STEP2_SCRIPT' --input '$STEP1_OUT' --out '$STEP2_OUT' $STEP2_EXTRA"
        fi
      fi

      # ----------------------------
      # Step 3 input staging
      # ----------------------------

      if [[ "$RUN_STEP3" == "TRUE" ]]; then
        if [[ -d "$STEP3_INPUT" && "$OVERWRITE" != "TRUE" ]]; then
          echo "[$(timestamp)] Step 3 staged input exists; skipping staging: $STEP3_INPUT"
        else
          rm -rf "$STEP3_INPUT"
          mkdir -p "$STEP3_INPUT"

          if [[ "$DRY_RUN" != "TRUE" ]]; then
            stage_step3_input "$STEP1_OUT" "$STEP2_OUT" "$STEP3_INPUT"
          else
            echo "[dry-run] stage_step3_input '$STEP1_OUT' '$STEP2_OUT' '$STEP3_INPUT'"
          fi
        fi
      fi

      # ----------------------------
      # Step 3: segment cells by cluster geometry
      # ----------------------------

      if [[ "$RUN_STEP3" == "TRUE" ]]; then
        if [[ -d "$STEP3_OUT" && "$OVERWRITE" != "TRUE" ]]; then
          echo "[$(timestamp)] Step 3 exists; skipping: $STEP3_OUT"
        else
          mkdir -p "$STEP3_OUT"

          run_cmd "Rscript '$STEP3_SCRIPT' --input '$STEP3_INPUT' --out '$STEP3_OUT' $STEP3_EXTRA"
        fi
      fi

      # ----------------------------
      # Step 4: correction and summary
      # ----------------------------

      if [[ "$RUN_CORRECTION" == "TRUE" ]]; then
        if [[ -d "$CORRECT_OUT" && "$OVERWRITE" != "TRUE" ]]; then
          echo "[$(timestamp)] Step 4 exists; skipping per-replicate correction/summary: $CORRECT_OUT"
        else
          rm -rf "$CORRECT_OUT"
          mkdir -p "$CORRECT_OUT"

          run_step4_summary \
            "$STEP3_OUT" \
            "$CORRECT_OUT" \
            "$SAMPLE_ID" \
            "$REPLICATE_ID" \
            "FALSE" \
            "1" \
            "2" \
            "3"
        fi
      fi

      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE_ID" \
        "$REPLICATE_ID" \
        "$REP_DIR" \
        "$STEP1_OUT" \
        "$STEP2_OUT" \
        "$STEP3_INPUT" \
        "$STEP3_OUT" \
        "complete" >> "$RUN_TABLE"

    } 2>&1 | tee -a "$REP_LOG"

  done
done


# ----------------------------
# Final project-level summary
# ----------------------------

PROJECT_SUMMARY_OUT="$OUT/$PROJECT_SUMMARY_DIR_NAME"

if [[ "$RUN_CORRECTION" == "TRUE" && "$RUN_PROJECT_SUMMARY" == "TRUE" ]]; then
  log "[wrapper] Running final project-level summary."

  if [[ -d "$PROJECT_SUMMARY_OUT" && "$OVERWRITE" != "TRUE" ]]; then
    log "[wrapper] Project-level summary exists; skipping: $PROJECT_SUMMARY_OUT"
  else
    rm -rf "$PROJECT_SUMMARY_OUT"
    mkdir -p "$PROJECT_SUMMARY_OUT"

    run_step4_summary \
      "$OUT" \
      "$PROJECT_SUMMARY_OUT" \
      "" \
      "" \
      "TRUE" \
      "1" \
      "2" \
      "4"
  fi

  if [[ "$DRY_RUN" != "TRUE" ]]; then
    normalize_step4_outputs "$PROJECT_SUMMARY_OUT"
    copy_project_summary_aliases "$PROJECT_SUMMARY_OUT"

    log "[wrapper] Final project-level summaries:"
    log "  $PROJECT_SUMMARY_OUT/image_summary.tsv"
    log "  $PROJECT_SUMMARY_OUT/replicate_summary.tsv"
    log "  $PROJECT_SUMMARY_OUT/sample_summary.tsv"
    log "  $PROJECT_SUMMARY_OUT/global_summary.tsv"
    log "[wrapper] Top-level compatibility aliases:"
    log "  $OUT/image_summary.corrected.tsv"
    log "  $OUT/replicate_summary.corrected.tsv"
    log "  $OUT/sample_summary.corrected.tsv"
    log "  $OUT/global_summary.corrected.tsv"
  fi
fi


# ----------------------------
# Step 05: mask confluence summary
# ----------------------------

CONFLUENCE_OUT="$OUT/$CONFLUENCE_DIR_NAME"

if [[ "$RUN_CONFLUENCE" == "TRUE" ]]; then
  log "[wrapper] Running Step 05 mask confluence summary."

  if [[ -d "$CONFLUENCE_OUT" && "$OVERWRITE" != "TRUE" ]]; then
    log "[wrapper] Step 05 confluence output exists; skipping: $CONFLUENCE_OUT"
  else
    rm -rf "$CONFLUENCE_OUT"
    mkdir -p "$CONFLUENCE_OUT"

    run_step5_confluence "$OUT" "$CONFLUENCE_OUT"
  fi

  if [[ "$DRY_RUN" != "TRUE" ]]; then
    log "[wrapper] Step 05 outputs:"
    log "  $CONFLUENCE_OUT/mask_confluence.per_image.tsv"
    log "  $CONFLUENCE_OUT/mask_confluence.by_sample_replicate.tsv"
  fi
else
  log "[wrapper] Skipping Step 05 mask confluence summary."
fi


# ----------------------------
# Final enriched rosette summary
# ----------------------------

ENRICHED_SUMMARY_OUT="$OUT/$ENRICHED_SUMMARY_DIR_NAME"

if [[ "$RUN_ENRICHED_SUMMARY" == "TRUE" ]]; then
  log "[wrapper] Running enriched post-correction rosette summary."

  if [[ -d "$ENRICHED_SUMMARY_OUT" && "$OVERWRITE" != "TRUE" ]]; then
    log "[wrapper] Enriched summary exists; skipping: $ENRICHED_SUMMARY_OUT"
  else
    rm -rf "$ENRICHED_SUMMARY_OUT"
    mkdir -p "$ENRICHED_SUMMARY_OUT"

    run_enriched_summary \
      "$OUT" \
      "$ENRICHED_SUMMARY_OUT" \
      "$CONFLUENCE_OUT/mask_confluence.per_image.tsv" \
      "$CONFLUENCE_OUT/mask_confluence.by_sample_replicate.tsv"
  fi

  if [[ "$DRY_RUN" != "TRUE" ]]; then
    copy_enriched_summary_aliases "$ENRICHED_SUMMARY_OUT"

    log "[wrapper] Enriched summary outputs:"
    log "  $ENRICHED_SUMMARY_OUT/rosette_counts.enriched.per_rosette.tsv"
    log "  $ENRICHED_SUMMARY_OUT/rosette_summary.enriched.per_image.tsv"
    log "  $ENRICHED_SUMMARY_OUT/rosette_summary.enriched.by_sample_replicate.tsv"
    log "[wrapper] Top-level enriched aliases:"
    log "  $OUT/rosette_counts.enriched.per_rosette.tsv"
    log "  $OUT/rosette_summary.enriched.per_image.tsv"
    log "  $OUT/rosette_summary.enriched.by_sample_replicate.tsv"
  fi
else
  log "[wrapper] Skipping enriched post-correction rosette summary."
fi

log "[wrapper] Done."
