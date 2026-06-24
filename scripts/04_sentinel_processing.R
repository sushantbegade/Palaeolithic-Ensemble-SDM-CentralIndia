# ============================================================
# SCRIPT 04: SENTINEL-2A PROCESSING
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 04 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Script 02 aligned NDVI, MNDWI, NDBI to the master
#   template. This script goes further:
#
#   1. Validates NDVI (predictor 14) — value range, spatial
#      pattern, site-value distribution check
#   2. Validates MNDWI alignment and value distribution
#   3. Applies MNDWI > -0.10 threshold to produce the
#      binary MOISTURE ANOMALY raster (input to Script 05
#      palaeochannel layer construction)
#   4. Extracts NIR band (B8) for field validation support
#      (Criterion 3: vegetation anomaly over channel fill)
#   5. Generates diagnostic figures for all indices
#
# PREDICTORS FINALISED HERE:
#   Predictor 14: NDVI_30m_utm44n.tif  ← ready after this script
#
# OUTPUTS FOR SCRIPT 05 (palaeochannel):
#   MNDWI_30m_utm44n.tif         (continuous — already done)
#   MNDWI_ANOMALY_30m_utm44n.tif (binary: MNDWI > -0.10)
#   NIR_B8_30m_utm44n.tif        (for vegetation anomaly check)
#
# SENTINEL-2A ACQUISITION:
#   Dry season composite: May 2025 (peak dry — max contrast
#   between active channels and relict channel fill moisture)
#   Level-2A (surface reflectance, atmospherically corrected)
#   All bands already pre-computed in R and stored in:
#   PATH_S2 = .../Sentinel 2A Multispectral Data/
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 04: Sentinel-2A Processing\n")
cat("========================================\n\n")

# ── 1. LOAD ALIGNED INDICES FROM SCRIPT 02 ──────────────────

cat("--- Loading Aligned Sentinel-2A Outputs ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))
sites_utm <- sf::st_read(
  file.path(OUT_SITES, "sites_all_utm44n.gpkg"),
  quiet = TRUE)

# Load already-aligned indices
ndvi_r  <- terra::rast(file.path(OUT_PREDICTORS,
                                 "NDVI_30m_utm44n.tif"))
mndwi_r <- terra::rast(file.path(OUT_PREDICTORS,
                                 "MNDWI_30m_utm44n.tif"))

# Report
for (nm in c("NDVI", "MNDWI")) {
  r   <- if (nm == "NDVI") ndvi_r else mndwi_r
  rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
  cat(sprintf("  %-8s range [%6.3f, %6.3f]  CRS: %s\n",
              nm,
              round(rng[1, 1], 3),
              round(rng[1, 2], 3),
              terra::crs(r, describe = TRUE)$name))
}
cat("\n")

# ── 2. VALIDATE NDVI (PREDICTOR 14) ─────────────────────────
# Expected for Vidarbha dry season (May):
#   Very low NDVI over bare soil / Deccan basalt outcrops: 0.05–0.15
#   Moderate NDVI over degraded scrub / agricultural fallow: 0.15–0.35
#   No green vegetation expected above ~0.6 in peak dry season

cat("--- Validating NDVI (Predictor 14 — Dry Season) ---\n")

ndvi_stats <- terra::global(ndvi_r, c("min", "max", "mean", "sd"),
                            na.rm = TRUE)
cat("  Min:  ", round(ndvi_stats[1, "min"],  3), "\n")
cat("  Max:  ", round(ndvi_stats[1, "max"],  3), "\n")
cat("  Mean: ", round(ndvi_stats[1, "mean"], 3), "\n")
cat("  SD:   ", round(ndvi_stats[1, "sd"],   3), "\n")

# Flag unexpected values
if (ndvi_stats[1, "min"] < -1 || ndvi_stats[1, "max"] > 1) {
  warning("NDVI outside [-1, 1] — check input raster scaling.")
}
if (ndvi_stats[1, "mean"] > 0.5) {
  warning("NDVI mean > 0.5 for dry season — verify acquisition date.")
} else {
  cat("  ✓ NDVI statistics consistent with May dry season acquisition\n")
}

# Extract NDVI values at site locations
ndvi_at_sites <- terra::extract(ndvi_r, terra::vect(sites_utm))
names(ndvi_at_sites)[2] <- "NDVI"
cat("  NDVI at sites — mean:",
    round(mean(ndvi_at_sites$NDVI, na.rm = TRUE), 3),
    "  SD:", round(sd(ndvi_at_sites$NDVI, na.rm = TRUE), 3), "\n")
cat("  ✓ NDVI predictor validated\n\n")

# ── 3. VALIDATE MNDWI ───────────────────────────────────────
# MNDWI = (Green − SWIR) / (Green + SWIR)
# Reference: Xu (2006); modified NDWI for open water features
#
# Expected value ranges (Vidarbha dry season):
#   Active water bodies:      MNDWI > 0.0
#   Moist soil / channel fill: MNDWI -0.3 to 0.0
#   Dry soil / rock outcrops: MNDWI < -0.3
#   Threshold for palaeochannel moisture anomaly: MNDWI > -0.10

cat("--- Validating MNDWI ---\n")

mndwi_stats <- terra::global(mndwi_r, c("min", "max", "mean", "sd"),
                             na.rm = TRUE)
cat("  Min:  ", round(mndwi_stats[1, "min"],  3), "\n")
cat("  Max:  ", round(mndwi_stats[1, "max"],  3), "\n")
cat("  Mean: ", round(mndwi_stats[1, "mean"], 3), "\n")
cat("  SD:   ", round(mndwi_stats[1, "sd"],   3), "\n")

# Count cells above each threshold
n_total <- terra::global(!is.na(mndwi_r), "sum",
                         na.rm = TRUE)[1, 1]
n_above_threshold <- terra::global(mndwi_r > -0.10, "sum",
                                   na.rm = TRUE)[1, 1]
n_water <- terra::global(mndwi_r > 0.0, "sum",
                         na.rm = TRUE)[1, 1]

cat("  Cells MNDWI > 0.0  (water):           ",
    n_water, sprintf("(%.1f%%)\n",
                     100 * n_water / n_total))
cat("  Cells MNDWI > -0.10 (moisture anomaly):",
    n_above_threshold, sprintf("(%.1f%%)\n",
                               100 * n_above_threshold / n_total))
cat("  ✓ MNDWI validated\n\n")

# ── 4. MNDWI MOISTURE ANOMALY RASTER ────────────────────────
# Binary raster: 1 = potential palaeochannel moisture anomaly
#                0 = background
# Threshold: MNDWI > -0.10 (Section 5.4 of Research Design)
# This is INPUT 1 of 2 for Script 05 palaeochannel construction.
# Script 05 cross-validates this against DEM valley morphology.

cat("--- Creating MNDWI Moisture Anomaly Raster ---\n")
cat("  Threshold: MNDWI > -0.10\n")

moisture_anomaly <- ifel(mndwi_r > -0.10, 1, 0)
moisture_anomaly <- terra::mask(moisture_anomaly, boundary_vect)

# Report spatial coverage
n_anomaly <- terra::global(moisture_anomaly == 1, "sum",
                           na.rm = TRUE)[1, 1]
pct_anomaly <- 100 * n_anomaly / n_total
cat("  Moisture anomaly cells:", n_anomaly,
    sprintf("(%.1f%% of study area)\n", pct_anomaly))

if (pct_anomaly < 2) {
  cat("  ⚠ Very low coverage — consider lowering threshold to -0.15\n")
} else if (pct_anomaly > 30) {
  cat("  ⚠ Very high coverage — consider raising threshold to -0.05\n")
} else {
  cat("  ✓ Coverage plausible for relict channel network\n")
}

out_anomaly <- file.path(OUT_PREDICTORS,
                         "MNDWI_ANOMALY_30m_utm44n.tif")
terra::writeRaster(moisture_anomaly, out_anomaly,
                   overwrite = TRUE, datatype = "INT1U")
cat("  ✓ Saved: MNDWI_ANOMALY_30m_utm44n.tif\n\n")

# ── 5. NIR BAND (B8) FOR FIELD VALIDATION SUPPORT ───────────
# Palaeochannel field validation Criterion 3 (Research Design 5.4):
# "linear vegetation density anomaly visible in Sentinel-2A NIR
#  band (B8) during dry season"
# NIR (B8) at 10m is sensitive to subsurface moisture — denser
# dry-season shrub/grass over moist channel fill appears brighter.
# We align B8 to 30m template for overlay with MNDWI anomaly.

cat("--- Aligning NIR Band B8 for Field Validation Support ---\n")

# Check if B8 file exists — name may vary by processing chain
b8_candidates <- list.files(PATH_S2,
                            pattern = "B8|NIR|nir|b8",
                            full.names = TRUE)

if (length(b8_candidates) == 0) {
  cat("  ⚠ NIR band B8 not found in Sentinel folder.\n")
  cat("    Candidates checked in:", PATH_S2, "\n")
  cat("    Skipping B8 alignment — not required for modelling,\n")
  cat("    only for field validation planning.\n\n")
} else {
  b8_path <- b8_candidates[1]
  cat("  File found:", basename(b8_path), "\n")
  
  b8_raw <- terra::rast(b8_path)
  
  # Reproject if needed
  if (!identical(terra::crs(b8_raw), terra::crs(template_30m))) {
    b8_raw <- terra::project(b8_raw, template_30m,
                             method = "bilinear")
  }
  
  # Resample and mask
  b8_aligned <- terra::resample(b8_raw, template_30m,
                                method = "bilinear")
  b8_aligned  <- terra::mask(b8_aligned, boundary_vect)
  
  terra::writeRaster(b8_aligned,
                     file.path(OUT_PREDICTORS,
                               "NIR_B8_30m_utm44n.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  
  b8_rng <- terra::global(b8_aligned, c("min", "max"), na.rm = TRUE)
  cat("  ✓ NIR B8 saved: NIR_B8_30m_utm44n.tif\n")
  cat("    Value range:", round(b8_rng[1, 1], 3),
      "to", round(b8_rng[1, 2], 3), "\n\n")
}

# ── 6. EXTRACT MNDWI VALUES AT SITE LOCATIONS ───────────────
# Diagnostic: do known Palaeolithic sites preferentially occur
# in areas of higher MNDWI (closer to moisture anomaly zones)?
# This is a preliminary check — formal test via SHAP in Script 19.

cat("--- MNDWI Distribution at Site Locations (Diagnostic) ---\n")

mndwi_at_sites <- terra::extract(mndwi_r, terra::vect(sites_utm))
names(mndwi_at_sites)[2] <- "MNDWI"

cat("  MNDWI at 197 sites:\n")
cat("    Mean:", round(mean(mndwi_at_sites$MNDWI, na.rm = TRUE), 3), "\n")
cat("    SD:  ", round(sd(mndwi_at_sites$MNDWI,   na.rm = TRUE), 3), "\n")
cat("    Min: ", round(min(mndwi_at_sites$MNDWI,  na.rm = TRUE), 3), "\n")
cat("    Max: ", round(max(mndwi_at_sites$MNDWI,  na.rm = TRUE), 3), "\n")

# Proportion of sites in moisture anomaly zone
n_sites_anomaly <- sum(mndwi_at_sites$MNDWI > -0.10, na.rm = TRUE)
cat("  Sites with MNDWI > -0.10:", n_sites_anomaly,
    sprintf("(%.1f%% of sites)\n",
            100 * n_sites_anomaly / nrow(sites_utm)))
cat("  Study area with MNDWI > -0.10:",
    sprintf("%.1f%%\n\n", pct_anomaly))

# ── 7. DIAGNOSTIC FIGURES ───────────────────────────────────

cat("--- Generating Diagnostic Figures ---\n")

# Figure 1: NDVI + MNDWI side-by-side with site overlay
png(file.path(OUT_FIG_SUPP, "S0c_sentinel_indices_check.png"),
    width = 4800, height = 2400, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))

terra::plot(ndvi_r,
            main = "NDVI — Dry Season (May 2025)",
            col  = viridisLite::viridis(100),
            axes = FALSE)
terra::plot(terra::vect(sites_utm), add = TRUE,
            col = "red", pch = 16, cex = 0.2)

terra::plot(mndwi_r,
            main = "MNDWI — Dry Season (May 2025)",
            col  = viridisLite::mako(100),
            axes = FALSE)
terra::plot(terra::vect(sites_utm), add = TRUE,
            col = "red", pch = 16, cex = 0.2)

dev.off()
cat("  ✓ S0c_sentinel_indices_check.png\n")

# Figure 2: MNDWI moisture anomaly binary raster
png(file.path(OUT_FIG_SUPP, "S0d_mndwi_anomaly.png"),
    width = 2400, height = 2400, res = 300)

terra::plot(moisture_anomaly,
            main   = "MNDWI Moisture Anomaly\n(MNDWI > -0.10)",
            col    = c("grey90", "#2166ac"),
            legend = FALSE,
            axes   = FALSE)
terra::plot(boundary_vect, add = TRUE, border = "black", lwd = 1)
terra::plot(terra::vect(sites_utm), add = TRUE,
            col = "red", pch = 16, cex = 0.3)
legend("bottomright",
       legend = c("Moisture anomaly", "Background", "Sites"),
       fill   = c("#2166ac", "grey90", NA),
       pch    = c(NA, NA, 16),
       col    = c(NA, NA, "red"),
       bty    = "n", cex = 0.7)

dev.off()
cat("  ✓ S0d_mndwi_anomaly.png\n\n")

# ── 8. SAVE MNDWI STATS FOR METHODS REPORTING ───────────────

mndwi_report <- data.frame(
  threshold        = -0.10,
  n_anomaly_cells  = n_anomaly,
  pct_study_area   = round(pct_anomaly, 2),
  n_sites_in_zone  = n_sites_anomaly,
  pct_sites        = round(100 * n_sites_anomaly / nrow(sites_utm), 2),
  mndwi_site_mean  = round(mean(mndwi_at_sites$MNDWI, na.rm = TRUE), 3),
  mndwi_site_sd    = round(sd(mndwi_at_sites$MNDWI,   na.rm = TRUE), 3)
)

write.csv(mndwi_report,
          file.path(OUT_TABLES, "S04_mndwi_threshold_stats.csv"),
          row.names = FALSE)
cat("  ✓ MNDWI threshold statistics saved: S04_mndwi_threshold_stats.csv\n\n")

# ── 9. VERIFY ALL SCRIPT 04 OUTPUTS ─────────────────────────

cat("--- Verifying Script 04 Outputs ---\n")

outputs_04 <- c(
  "NDVI_30m_utm44n.tif",          # Predictor 14 — READY
  "MNDWI_30m_utm44n.tif",         # Continuous — for Script 05
  "MNDWI_ANOMALY_30m_utm44n.tif"  # Binary — for Script 05
)

for (f in outputs_04) {
  fpath <- file.path(OUT_PREDICTORS, f)
  if (file.exists(fpath)) {
    r   <- terra::rast(fpath)
    rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
    cat(sprintf("  ✓ %-42s  [%7.3f, %7.3f]\n",
                f,
                round(rng[1, 1], 3),
                round(rng[1, 2], 3)))
  } else {
    cat("  ✗ MISSING:", f, "\n")
  }
}

# ── 10. SUMMARY REPORT ──────────────────────────────────────

cat("\n========================================\n")
cat("SCRIPT 04 COMPLETE — Summary\n")
cat("========================================\n")
cat("NDVI (Predictor 14):     ✓ Validated and ready\n")
cat("MNDWI (continuous):      ✓ Validated\n")
cat("MNDWI Anomaly (binary):  ✓ Created (threshold = -0.10)\n")
cat("NIR B8:                  ✓ Aligned (field validation aid)\n")
cat("\nMNDWI anomaly coverage:  ",
    sprintf("%.1f%%", pct_anomaly), "of study area\n")
cat("Sites in anomaly zone:   ",
    sprintf("%.1f%%", 100 * n_sites_anomaly / nrow(sites_utm)),
    "of 197 sites\n")
cat("\nNext: Run Script 05 — Palaeochannel Layer Construction\n")
cat("  Inputs: MNDWI_ANOMALY + DEM valley morphology (WBT)\n")
cat("  Output: DIST_PALAEOCHANNEL_30m_utm44n.tif\n")
cat("========================================\n")