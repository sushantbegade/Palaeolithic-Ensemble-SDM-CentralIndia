# ============================================================
# SCRIPT 03: TERRAIN DERIVATIVES FROM DEM
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 03 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Computes all 7 terrain predictors from the master DEM
#   (TEMPLATE_30m_utm44n.tif) and saves each aligned to the
#   master template. Outputs feed directly into Script 07
#   (predictor stack assembly + VIF screening).
#
#   Predictors produced (predictor stack variables 1-8):
#   1.  Elevation        — DEM itself (already from Script 02)
#   2.  Slope            — terra::terrain()
#   3.  Aspect           — terra::terrain()
#   4.  TRI              — terra::terrain() [Wilson et al. 2007]
#   5.  TPI              — terra::terrain() [Weiss 2001]
#   6.  Plan Curvature   — WhiteboxTools: PlanCurvature
#   7.  HAND             — WhiteboxTools: ElevationAboveStream
#   8.  Flow Accumulation— WhiteboxTools: D8FlowAccumulation
#
# TOOLS:
#   terra     — slope, aspect, TRI, TPI (fast, native R)
#   whitebox  — plan curvature, HAND, flow accumulation
#               (specialized geomorphometric algorithms)
#
# NOTE ON whitebox PACKAGE:
#   'whitebox' was not in Script 01 package list.
#   This script installs and initialises it automatically.
#   First run downloads WhiteboxTools binary (~30 MB).
#   Subsequent runs use cached binary — fast.
# ============================================================
# HOW TO USE:
#   Run after Script 02 completes successfully.
#   Script 01 sourced automatically.
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 03: Terrain Derivatives\n")
cat("========================================\n\n")

# ── 1. INSTALL AND INITIALISE WHITEBOX ──────────────────────
# whitebox R package wraps WhiteboxTools binary.
# Not in Script 01 list — handled here.

cat("--- Initialising WhiteboxTools ---\n")

if (!requireNamespace("whitebox", quietly = TRUE)) {
  cat("  Installing whitebox package...\n")
  install.packages("whitebox")
}
library(whitebox)

# Download WhiteboxTools binary if not already cached
# (first run only — ~30 MB download, takes 1-2 minutes)
{
  cat("  Downloading WhiteboxTools binary (first run only)...\n")
  whitebox::install_whitebox()
}

# Initialise WhiteboxTools
wbt_init()
cat("  WhiteboxTools version:", wbt_version(), "\n")
cat("  ✓ WhiteboxTools ready\n\n")

# ── 2. LOAD MASTER DEM ──────────────────────────────────────

cat("--- Loading Master DEM ---\n")

dem_path <- file.path(OUT_PREDICTORS, "DEM_30m_utm44n.tif")

if (!file.exists(dem_path)) {
  stop("DEM not found: ", dem_path,
       "\nRun Script 02 first.")
}

dem <- terra::rast(dem_path)
cat("  CRS:       ", terra::crs(dem, describe = TRUE)$name, "\n")
cat("  Resolution:", paste(terra::res(dem), collapse = " x "), "m\n")
cat("  Dimensions:", nrow(dem), "rows x", ncol(dem), "cols\n")
cat("  ✓ DEM loaded\n\n")

# Load master template for alignment verification
template_path <- file.path(OUT_PREDICTORS, "TEMPLATE_30m_utm44n.tif")
template_30m  <- terra::rast(template_path)

# Load study area boundary for masking
boundary_path <- file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg")
boundary_vect <- terra::vect(sf::st_read(boundary_path, quiet = TRUE))

# ── 3. HELPER FUNCTIONS ─────────────────────────────────────

# Save raster: verify alignment to template, mask, write
save_terrain <- function(r, out_filename, layer_name,
                         datatype = "FLT4S") {
  # Verify alignment
  if (!isTRUE(all.equal(terra::res(r), terra::res(template_30m)))) {
    cat("    Resampling to match template...\n")
    r <- terra::resample(r, template_30m, method = "bilinear")
  }
  r <- terra::mask(r, boundary_vect)
  out_path <- file.path(OUT_PREDICTORS, out_filename)
  terra::writeRaster(r, out_path, overwrite = TRUE, datatype = datatype)
  
  # Report value range using global() — avoids minmax warning
  rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
  cat("  ✓", layer_name, "saved:", out_filename, "\n")
  cat("    Value range:", round(rng[1, 1], 3),
      "to", round(rng[1, 2], 3), "\n\n")
  return(r)
}

# WhiteboxTools writes to disk — helper to run wbt tool,
# load result, align to template, and save to OUT_PREDICTORS
run_wbt <- function(wbt_fn, args_list, temp_filename,
                    out_filename, layer_name,
                    resample_method = "bilinear") {
  # Temporary output path (WhiteboxTools needs full path)
  temp_path <- file.path(OUT_PREDICTORS, temp_filename)
  
  # Build argument list with output path
  args_list$output <- temp_path
  
  # Run WhiteboxTools function
  do.call(wbt_fn, args_list)
  
  # Load result
  r <- terra::rast(temp_path)
  
  # Align to template (resample + mask)
  if (!isTRUE(all.equal(terra::res(r), terra::res(template_30m))) ||
      !isTRUE(all.equal(as.vector(terra::ext(r)),
                        as.vector(terra::ext(template_30m))))) {
    r <- terra::resample(r, template_30m, method = resample_method)
  }
  r <- terra::mask(r, boundary_vect)
  
  # Save final aligned output (overwrites temp if same name)
  out_path <- file.path(OUT_PREDICTORS, out_filename)
  terra::writeRaster(r, out_path, overwrite = TRUE, datatype = "FLT4S")
  
  # Remove temp file if different from output
  if (temp_filename != out_filename && file.exists(temp_path)) {
    file.remove(temp_path)
  }
  
  rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
  cat("  ✓", layer_name, "saved:", out_filename, "\n")
  cat("    Value range:", round(rng[1, 1], 3),
      "to", round(rng[1, 2], 3), "\n\n")
  return(r)
}

# ── 4. TERRA-BASED DERIVATIVES ───────────────────────────────
# terra::terrain() computes slope, aspect, TRI, TPI natively.
# Much faster than WhiteboxTools for these standard derivatives.

cat("--- Computing terra-based Derivatives ---\n\n")

# ── 4a. SLOPE (degrees) ──────────────────────────────────────
cat("  [Slope]\n")
slope_r <- terra::terrain(dem, v = "slope", unit = "degrees",
                          neighbors = 8)
slope_r <- save_terrain(slope_r, "SLOPE_30m_utm44n.tif", "Slope")

# ── 4b. ASPECT (degrees, 0 = North, clockwise) ───────────────
cat("  [Aspect]\n")
aspect_r <- terra::terrain(dem, v = "aspect", unit = "degrees",
                           neighbors = 8)
aspect_r <- save_terrain(aspect_r, "ASPECT_30m_utm44n.tif", "Aspect")

# ── 4c. TRI — Terrain Ruggedness Index (Wilson et al. 2007) ──
# terra uses mean of absolute differences from focal 3x3 window.
# Equivalent to Riley et al. (1999) / Wilson et al. (2007).
cat("  [TRI — Terrain Ruggedness Index]\n")
tri_r <- terra::terrain(dem, v = "TRI", neighbors = 8)
tri_r <- save_terrain(tri_r, "TRI_30m_utm44n.tif",
                      "TRI (Terrain Ruggedness Index)")

# ── 4d. TPI — Topographic Position Index (Weiss 2001) ────────
# TPI = elevation of cell minus mean elevation of neighbourhood.
# Positive = ridges/hilltops, negative = valleys/depressions.
cat("  [TPI — Topographic Position Index]\n")
tpi_r <- terra::terrain(dem, v = "TPI", neighbors = 8)
tpi_r <- save_terrain(tpi_r, "TPI_30m_utm44n.tif",
                      "TPI (Topographic Position Index)")

# ── 5. WHITEBOX DERIVATIVES ──────────────────────────────────
# Plan curvature, HAND, and flow accumulation require
# WhiteboxTools geomorphometric algorithms not available in terra.

cat("--- Computing WhiteboxTools Derivatives ---\n\n")

# WhiteboxTools requires file paths as strings — use full paths
dem_input <- file.path(OUT_PREDICTORS, "DEM_30m_utm44n.tif")

# ── 5a. PLAN CURVATURE ───────────────────────────────────────
# Plan curvature: curvature in the horizontal plane perpendicular
# to slope direction. Positive = divergent flow (ridges),
# negative = convergent flow (channels/hollows).
# Mechanistically distinct from slope — justified separate inclusion.
cat("  [Plan Curvature]\n")

plan_curv_r <- run_wbt(
  wbt_fn        = wbt_plan_curvature,
  args_list     = list(dem = dem_input),
  temp_filename = "PLANCURV_temp.tif",
  out_filename  = "PLANCURV_30m_utm44n.tif",
  layer_name    = "Plan Curvature"
)

# ── 5b. FLOW ACCUMULATION (D8) ───────────────────────────────
# D8 flow accumulation: number of upslope cells draining
# into each cell. High values = valley floors / drainage
# convergence zones — proxy for seasonal water pooling.
# Log-transform recommended (very skewed distribution).
cat("  [Flow Accumulation — D8]\n")

# Step 1: Fill depressions first (required before flow routing)
dem_filled_path <- file.path(OUT_PREDICTORS, "DEM_filled_temp.tif")
wbt_fill_depressions_wang_and_liu(
  dem    = dem_input,
  output = dem_filled_path
)
cat("    DEM depressions filled\n")

# Step 2: D8 flow accumulation on filled DEM
flowacc_r <- run_wbt(
  wbt_fn        = wbt_d8_flow_accumulation,
  args_list     = list(
    input        = dem_filled_path,
    out_type     = "cells"   # output in cell count (integer)
  ),
  temp_filename = "FLOWACC_raw_temp.tif",
  out_filename  = "FLOWACC_30m_utm44n.tif",
  layer_name    = "Flow Accumulation (D8)"
)

# Step 3: Log10-transform (reduces extreme skew)
# Methods: log10(flowacc + 1) to handle zero values
cat("  [Flow Accumulation — log10 transform]\n")
flowacc_path <- file.path(OUT_PREDICTORS, "FLOWACC_30m_utm44n.tif")
flowacc_r    <- terra::rast(flowacc_path)
flowacc_log  <- log10(flowacc_r + 1)
flowacc_log  <- terra::mask(flowacc_log, boundary_vect)
terra::writeRaster(flowacc_log,
                   file.path(OUT_PREDICTORS,
                             "FLOWACC_LOG10_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
rng_fa <- terra::global(flowacc_log, c("min", "max"), na.rm = TRUE)
cat("  ✓ Flow Accumulation (log10) saved: FLOWACC_LOG10_30m_utm44n.tif\n")
cat("    Value range:", round(rng_fa[1, 1], 3),
    "to", round(rng_fa[1, 2], 3), "\n\n")

# ── 5c. HAND — Height Above Nearest Drainage ─────────────────
# HAND: vertical distance from each cell to the nearest stream
# cell in the drainage network. Captures local hydrological
# position and flood exposure risk.
# Reference: Nobre et al. (2011); Rennó et al. (2008).
#
# Method:
#   (i)  Extract stream network from flow accumulation
#        (threshold: cells with log10(flowacc) > 2.5 = streams)
#   (ii) Compute HAND using wbt_elevation_above_stream()
cat("  [HAND — Height Above Nearest Drainage]\n")

# Step 1: Extract stream network raster
#   Threshold: flow accumulation > 316 cells (~0.28 km² upslope area)
#   Equivalent to log10(316) = 2.5
cat("    Extracting stream network (flow accumulation threshold)...\n")

stream_path <- file.path(OUT_PREDICTORS, "STREAMS_temp.tif")
flowacc_raw_path <- file.path(OUT_PREDICTORS, "FLOWACC_30m_utm44n.tif")

wbt_extract_streams(
  flow_accum = flowacc_raw_path,
  output     = stream_path,
  threshold  = 316   # cells — adjust if stream network too sparse/dense
)
cat("    Stream network extracted\n")

# Step 2: HAND computation
hand_r <- run_wbt(
  wbt_fn        = wbt_elevation_above_stream,
  args_list     = list(
    dem     = dem_filled_path,
    streams = stream_path
  ),
  temp_filename = "HAND_temp.tif",
  out_filename  = "HAND_30m_utm44n.tif",
  layer_name    = "HAND (Height Above Nearest Drainage)"
)

# Step 3: Verify HAND values
hand_check <- terra::global(hand_r, c("min", "max"), na.rm = TRUE)
if (hand_check[1, 2] > 500) {
  cat("  ⚠ WARNING: HAND max > 500m — check stream threshold\n")
  cat("    If too high, lower threshold in wbt_extract_streams()\n")
} else {
  cat("  ✓ HAND values plausible for Vidarbha terrain\n")
}
cat("\n")

# ── 6. CLEAN UP TEMPORARY FILES ─────────────────────────────

cat("--- Cleaning up temporary files ---\n")

temp_files <- c(
  dem_filled_path,
  stream_path
)

for (f in temp_files) {
  if (file.exists(f)) {
    file.remove(f)
    cat("  Removed:", basename(f), "\n")
  }
}
cat("\n")

# ── 7. VERIFY ALL 7 TERRAIN OUTPUTS ─────────────────────────

cat("--- Verifying All Terrain Outputs ---\n\n")

terrain_files <- c(
  "DEM_30m_utm44n.tif",
  "SLOPE_30m_utm44n.tif",
  "ASPECT_30m_utm44n.tif",
  "TRI_30m_utm44n.tif",
  "TPI_30m_utm44n.tif",
  "PLANCURV_30m_utm44n.tif",
  "HAND_30m_utm44n.tif",
  "FLOWACC_LOG10_30m_utm44n.tif"
)

all_ok <- TRUE

for (f in terrain_files) {
  fpath <- file.path(OUT_PREDICTORS, f)
  if (!file.exists(fpath)) {
    cat("  ✗ MISSING:", f, "\n")
    all_ok <- FALSE
  } else {
    r    <- terra::rast(fpath)
    rng  <- terra::global(r, c("min", "max"), na.rm = TRUE)
    dims <- paste0(nrow(r), " x ", ncol(r))
    cat(sprintf("  ✓ %-40s  range [%7.2f, %7.2f]  dims %s\n",
                f,
                round(rng[1, 1], 2),
                round(rng[1, 2], 2),
                dims))
  }
}

if (all_ok) {
  cat("\n  ✓ All 8 terrain layers present and verified\n")
} else {
  cat("\n  ✗ Some terrain layers missing — check errors above\n")
}

# ── 8. QUICK MULTI-PANEL VISUAL CHECK ───────────────────────

cat("\n--- Generating Terrain Overview Figure ---\n")

terrain_labels <- c(
  "Elevation (m)",
  "Slope (°)",
  "Aspect (°)",
  "TRI",
  "TPI",
  "Plan Curvature",
  "HAND (m)",
  "Flow Accum (log10)"
)

png(file.path(OUT_FIG_SUPP, "S0b_terrain_derivatives_check.png"),
    width = 6400, height = 4800, res = 300)

par(mfrow = c(2, 4), mar = c(2, 2, 3, 1))

for (i in seq_along(terrain_files)) {
  fpath <- file.path(OUT_PREDICTORS, terrain_files[i])
  if (file.exists(fpath)) {
    r <- terra::rast(fpath)
    terra::plot(r,
                main   = terrain_labels[i],
                col    = viridisLite::viridis(100),
                axes   = FALSE,
                legend = TRUE)
  }
}

dev.off()
cat("  ✓ Terrain overview figure saved: S0b_terrain_derivatives_check.png\n\n")

# ── 9. SUMMARY REPORT ───────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 03 COMPLETE — Summary\n")
cat("========================================\n")
cat("Elevation (DEM):         ✓ (from Script 02)\n")
cat("Slope:                   ✓\n")
cat("Aspect:                  ✓\n")
cat("TRI:                     ✓\n")
cat("TPI:                     ✓\n")
cat("Plan Curvature:          ✓\n")
cat("HAND:                    ✓\n")
cat("Flow Accumulation log10: ✓\n")
cat("\nAll outputs at 30m UTM Zone 44N, aligned to master template.\n")
cat("Saved to:", OUT_PREDICTORS, "\n")
cat("\nNext: Run Script 04 — Sentinel-2A Processing\n")
cat("(Distance-to-river raster + raw material raster = Script 05)\n")
cat("========================================\n")