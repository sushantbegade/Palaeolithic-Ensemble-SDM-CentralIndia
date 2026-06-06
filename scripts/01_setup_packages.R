# ============================================================
# SCRIPT 01: SETUP, FILE PATHS, AND PACKAGE INSTALLATION
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 01 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1. Defines ALL file paths used across all 25 scripts
#   2. Installs and loads all required packages
#   3. Sets global parameters (CRS, resolution, random seed)
#   4. Initialises renv for reproducibility
# ============================================================
# HOW TO USE:
#   Run this script at the start of every R session.
#   All other scripts (02-25) source this file automatically.
#   Never hardcode paths in other scripts — use variables
#   defined here.
# ============================================================

# ── 0. REPRODUCIBILITY ──────────────────────────────────────

set.seed(42)  # Fixed random seed — state in Methods Section 5.15

# Redirect terra temp files to E drive — prevents C drive overflow
# during heavy raster processing (Scripts 11-22)
E_TEMP <- "E:/R_temp"
dir.create(E_TEMP, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = E_TEMP)

# ── 1. BASE DIRECTORY ───────────────────────────────────────
# Change ONLY this one path if the project moves to a
# different drive or computer. All other paths update
# automatically.

BASE <- "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset"

# ── 2. INPUT DATA PATHS (raw data — never modified) ─────────

# Site data
PATH_SITES_XLSX <- file.path(BASE,
                             "Nagpur Chandrapur Palaeolithic Site Raw Data.xlsx")
PATH_SITES_CSV  <- file.path(BASE,
                             "Nagpur Chandrapur Palaeolithic Site Raw Data.csv")

# Administrative boundary
PATH_DISTRICT   <- file.path(BASE, "District Boundary")

# Terrain (DEM)
PATH_DEM        <- file.path(BASE, "DEM")

# Geological predictors
PATH_GEOLOGY    <- file.path(BASE, "Geology")
PATH_GEOMORPH   <- file.path(BASE, "Geomorphology")
PATH_LITHOLOGY  <- file.path(BASE, "Lithology")

# Hydrological predictors
PATH_RIVERS     <- file.path(BASE, "Rivers")
PATH_WATERBODY  <- file.path(BASE, "Waterbody")

# Sentinel-2A (processed bands and indices)
PATH_S2         <- file.path(BASE, "Sentinel 2A Multispectral Data")

# Climate data (CHELSA-TraCE21k + modern + WorldClim)
PATH_CLIMATE    <- file.path(BASE, "Climate Data")

# Supplementary geological data
PATH_SOIL       <- file.path(BASE, "Soil Data")
PATH_FAULT      <- file.path(BASE, "Fault")
PATH_DYKE       <- file.path(BASE, "Dyke")
PATH_LINEAMENT  <- file.path(BASE, "Linement")

# ── 3. OUTPUT PATHS (R_Analysis — all generated outputs) ────

R_ANALYSIS      <- file.path(BASE, "R_Analysis")

SCRIPTS_DIR <- file.path(R_ANALYSIS,
                         "Palaeolithic-Ensemble-SDM-CentralIndia/scripts")

# Processed data
OUT_PREDICTORS  <- file.path(R_ANALYSIS, "data_processed/predictors")
OUT_SITES       <- file.path(R_ANALYSIS, "data_processed/sites")
OUT_BACKGROUND  <- file.path(R_ANALYSIS, "data_processed/background")
OUT_BIAS        <- file.path(R_ANALYSIS, "data_processed/bias_surface")
OUT_PALAEO      <- file.path(R_ANALYSIS, "data_processed/palaeochannel")
OUT_CLIMATE_PRO <- file.path(R_ANALYSIS, "data_processed/climate")
OUT_CV          <- file.path(R_ANALYSIS, "data_processed/cv_blocks")

# Models
OUT_MOD_IND     <- file.path(R_ANALYSIS, "models/individual")
OUT_MOD_ENS     <- file.path(R_ANALYSIS, "models/ensemble")
OUT_MOD_SUB     <- file.path(R_ANALYSIS, "models/submodels")

# Outputs
OUT_FIG_MAIN    <- file.path(R_ANALYSIS, "outputs/figures/main")
OUT_FIG_SUPP    <- file.path(R_ANALYSIS, "outputs/figures/supplementary")
OUT_TABLES      <- file.path(R_ANALYSIS, "outputs/tables")
OUT_SHAP        <- file.path(R_ANALYSIS, "outputs/shap")
OUT_EVAL        <- file.path(R_ANALYSIS, "outputs/evaluation")
OUT_SUPP_AN     <- file.path(R_ANALYSIS, "outputs/supplementary_analyses")

# ── 4. CREATE ALL OUTPUT FOLDERS (safe — skips if exists) ───

all_out_dirs <- c(
  OUT_PREDICTORS, OUT_SITES, OUT_BACKGROUND, OUT_BIAS,
  OUT_PALAEO, OUT_CLIMATE_PRO, OUT_CV,
  OUT_MOD_IND, OUT_MOD_ENS, OUT_MOD_SUB,
  OUT_FIG_MAIN, OUT_FIG_SUPP, OUT_TABLES,
  OUT_SHAP, OUT_EVAL, OUT_SUPP_AN
)

for (d in all_out_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
cat("✓ All output folders ready\n")

# ── 5. GLOBAL PARAMETERS ────────────────────────────────────

# Coordinate Reference System — WGS84 UTM Zone 44N
# EPSG:32644 — standard for central India
CRS_PROJECT <- "EPSG:32644"

# Target raster resolution (metres)
RESOLUTION <- 30

# Study area extent (WGS84 geographic — for initial loading)
# Nagpur + Chandrapur: 19.3N-21.8N, 78.1E-80.2E
STUDY_EXTENT_GEO <- c(
  xmin = 78.10, xmax = 80.20,
  ymin = 19.30, ymax = 21.80
)

# Number of background points for bias-corrected sampling
N_BACKGROUND <- 10000

# Spatial thinning distance (metres)
THIN_DIST_KM <- 1

# VIF threshold for multicollinearity screening
VIF_THRESHOLD <- 5

# SHAP grid resolution for dominant driver map (metres)
SHAP_RESOLUTION <- 250

# MaxEnt logistic output type — CRITICAL — never change this
MAXENT_OUTPUT_TYPE <- "logistic"

# ── 6. PACKAGE INSTALLATION ─────────────────────────────────
# This block checks if each package is installed.
# If not, it installs from CRAN automatically.
# Run this only once — takes 10-20 minutes first time.

packages_needed <- c(
  # ── Core spatial ──────────────────────────────────────────
  "terra",          # Raster and vector processing
  "sf",             # Vector data (shapefiles)
  "tidyterra",      # terra + ggplot2 integration
  
  # ── Data handling ─────────────────────────────────────────
  "readxl",         # Read Excel files
  "writexl",        # Write Excel files
  "dplyr",          # Data manipulation
  "tidyr",          # Data reshaping
  "tibble",         # Modern data frames
  
  # ── Predictive modelling ──────────────────────────────────
  "maxnet",         # MaxEnt (logistic output — no Java)
  "ENMeval",        # MaxEnt parameter tuning by AICc
  "randomForest",   # Random Forest
  "xgboost",        # XGBoost
  "gbm",            # Boosted Regression Trees
  "mgcv",           # Generalised Additive Models
  "kernlab",        # Support Vector Machine
  "glmnet",         # LASSO meta-learner (Tier 2 ensemble)
  "dismo",          # randomPoints() for background sampling
  "PresenceAbsence", # TSS and threshold optimisation
  
  # ── Cross-validation and thinning ─────────────────────────
  "blockCV",        # Spatial block cross-validation
  "spThin",         # Spatial thinning at 1km
  
  # ── Bias correction ───────────────────────────────────────
  "ks",             # KDE for survey effort surface
  
  # ── Explainability ────────────────────────────────────────
  "shapviz",        # SHAP values and visualisation
  
  # ── Model evaluation ──────────────────────────────────────
  "pROC",           # AUC, ROC, DeLong's test
  "ecospat",        # Boyce Index
  "usdm",           # VIF analysis
  
  # ── Visualisation ─────────────────────────────────────────
  "ggplot2",        # Publication-ready figures
  "patchwork",      # Multi-panel figure layout
  "cowplot",        # Figure alignment and themes
  "viridis",        # Colour-blind-safe palettes
  "RColorBrewer",   # Additional colour palettes
  "ggpubr",         # Figure annotation and formatting
  "ggspatial",      # Scale bars and north arrows
  "scales",         # Axis scaling
  
  # ── Reproducibility ───────────────────────────────────────
  "renv"            # Package version locking
)

# Install missing packages
installed <- rownames(installed.packages())
to_install <- packages_needed[!packages_needed %in% installed]

if (length(to_install) > 0) {
  cat("Installing", length(to_install), "packages...\n")
  install.packages(to_install, dependencies = TRUE)
  cat("✓ Installation complete\n")
} else {
  cat("✓ All packages already installed\n")
}

# ── 7. LOAD ALL PACKAGES ────────────────────────────────────

suppressPackageStartupMessages({
  for (pkg in packages_needed) {
    library(pkg, character.only = TRUE)
  }
})
cat("✓ All packages loaded\n")

# ── 8. VERIFY CRITICAL VERSIONS ─────────────────────────────
# These versions are the minimum required for this analysis

check_version <- function(pkg, min_ver) {
  v <- packageVersion(pkg)
  status <- if (v >= min_ver) "✓" else "✗ UPDATE NEEDED"
  cat(sprintf("  %s %-20s version %s (min: %s)\n",
              status, pkg, v, min_ver))
}

cat("\nCritical package versions:\n")
check_version("terra",         "1.7.0")
check_version("maxnet",        "0.1.4")
check_version("ENMeval",       "2.0.0")
check_version("randomForest",  "4.7.0")
check_version("xgboost",       "1.7.0")
check_version("gbm",           "2.1.8")
check_version("mgcv",          "1.9.0")
check_version("kernlab",       "0.9.0")
check_version("blockCV",       "3.1.0")
check_version("shapviz",       "0.9.0")
check_version("pROC",          "1.18.0")

# ── 9. INITIALISE RENV ──────────────────────────────────────
# Run renv::snapshot() after installing all packages
# to lock package versions for reproducibility.
# Only run once — after this, renv.lock is on GitHub.

# renv::init()       # Run once to initialise
# renv::snapshot()   # Run once after all packages confirmed working

cat("\n========================================\n")
cat("Script 01 complete. Environment ready.\n")
cat("CRS:        ", CRS_PROJECT, "\n")
cat("Resolution: ", RESOLUTION, "m\n")
cat("Background: ", N_BACKGROUND, "points\n")
cat("Seed:       42\n")
cat("========================================\n")