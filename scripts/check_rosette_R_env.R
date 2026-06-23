#!/usr/bin/env Rscript

# ============================================================
# check_rosette_R_env.R
#
# Verifies that the active R environment has all packages needed
# for the rosette image-analysis pipeline.
#
# Important:
#   This script intentionally does NOT install packages.
#   EBImage/magick are best installed through conda/mamba using:
#     bioconductor-ebimage
#     r-magick
#
# Runtime package installation tends to fail because EBImage/magick
# have compiled system/image-library dependencies.
# ============================================================

required_packages <- c(
  "data.table",
  "optparse",
  "jsonlite",
  "stringr",
  "EBImage"
)

message("[env_check] R executable: ", file.path(R.home("bin"), "R"))
message("[env_check] Rscript executable: ", file.path(R.home("bin"), "Rscript"))
message("[env_check] R version: ", R.version.string)

message("[env_check] Library paths:")
message(paste("  -", .libPaths()), sep = "\n")

message("[env_check] Checking required packages...")

missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  message("[env_check] Missing packages:")
  message(paste("  -", missing), sep = "\n")

  message("")
  message("[env_check] Recommended fix:")
  message("  Install missing packages through conda/mamba, not install.packages() or BiocManager at runtime.")
  message("")
  message("  Your envs/rosette_pipeline.yml should include at least:")
  message("    - r-data.table")
  message("    - r-optparse")
  message("    - r-jsonlite")
  message("    - r-stringr")
  message("    - r-magick")
  message("    - imagemagick")
  message("    - bioconductor-ebimage")
  message("")
  message("  Then recreate/update the environment:")
  message("    conda env remove -n rosette_pipeline -y")
  message("    mamba env create -f envs/rosette_pipeline.yml")
  message("  or:")
  message("    conda env create -f envs/rosette_pipeline.yml")
  message("")

  stop(
    "[env_check] Missing required R packages: ",
    paste(missing, collapse = ", ")
  )
}

message("[env_check] Loading packages to verify shared libraries...")

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
  library(stringr)
  library(EBImage)
})

message("[env_check] Package versions:")

for (pkg in required_packages) {
  ver <- as.character(utils::packageVersion(pkg))
  message("  - ", pkg, ": ", ver)
}

message("[env_check] Session info:")
si <- sessionInfo()
message("  R version: ", si$R.version$version.string)
message("  Platform: ", si$platform)
message("  Running under: ", si$running)

message("[env_check] R environment looks good.")
