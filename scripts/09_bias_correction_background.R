# ============================================================
# SCRIPT 09: SURVEY-EFFORT BIAS CORRECTION + BACKGROUND SAMPLING
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 09 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Implements survey-effort bias correction following the
#   target-group approach (Phillips et al. 2009). Site records
#   concentrate where archaeologists have surveyed — without
#   correction, models partially predict survey history rather
#   than genuine hominin environmental preferences.
#
#   METHOD:
#   A 2D kernel density estimate (KDE) of site distribution
#   is used as a proxy for survey effort intensity. This is
#   the standard approach when explicit survey-effort polygon
#   shapefiles are unavailable (Fourcade et al. 2014). The KDE
#   surface is rasterized and used as sampling weights for
#   background point generation — background points are drawn
#   disproportionately from areas of higher survey effort,
#   ensuring models contrast presences against a background
#   with similar observation bias.
#
#   KDE implementation:
#   - Package: ks (Duong 2007)
#   - Bandwidth: ks::Hscv() — cross-validated smoothed estimator
#   - Evaluated on 1km grid, resampled to 30m template
#   - Input: ALL 197 site UTM coordinates (pre-thinning)
#     (thinned sites used for modelling; all sites used for
#      bias surface so KDE reflects full survey coverage)
#
#   BACKGROUND POINTS:
#   - N = 10,000 (Warton & Shepherd 2010 justification)
#   - Sampled proportional to KDE bias surface
#   - Spatial extent: Nagpur-Chandrapur district ONLY
#   - Package: terra::spatSample() with prob weighting
#
#   SENSITIVITY ANALYSIS (Supplementary Table S5):
#   N = 1,000 / 5,000 / 10,000 / 20,000 background points
#   Run full ensemble with each — report spatial CV AUC.
#   If AUC stable (variation < 0.02), N=10,000 confirmed.
#   Background sets for all four counts saved here.
#
# OUTPUTS (to data_processed/background/):
#   bias_surface_kde_30m_utm44n.tif   — KDE bias raster
#   background_N10000.gpkg            — primary (10,000 pts)
#   background_N1000.gpkg             — sensitivity S5
#   background_N5000.gpkg             — sensitivity S5
#   background_N20000.gpkg            — sensitivity S5
#   S09_bias_surface_check.png        — diagnostic figure
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

# Ensure terra temp on E drive
terra::terraOptions(tempdir = "E:/R_temp")

cat("\n========================================\n")
cat("SCRIPT 09: Bias Correction + Background\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# Load ALL 197 sites (pre-thinning) for bias surface
# Thinned sites used in modelling; all sites used here
# so KDE reflects full archaeological survey coverage
sites_sf <- sf::st_read(file.path(OUT_SITES,
                                  "sites_all_utm44n.gpkg"),
                        quiet = TRUE)
cat("  Sites for KDE (all, pre-thinning):", nrow(sites_sf), "\n")

coords_utm <- sf::st_coordinates(sites_sf)
cat("  UTM coordinate range:\n")
cat("    X:", round(range(coords_utm[,1])), "\n")
cat("    Y:", round(range(coords_utm[,2])), "\n\n")

# ── 2. COMPUTE KDE BIAS SURFACE ──────────────────────────────

cat("--- Computing KDE Bias Surface ---\n")
cat("  Package: ks\n")
cat("  Bandwidth: Hscv() cross-validated smoothed estimator\n\n")

# Coordinate matrix for KDE
X <- cbind(coords_utm[, 1], coords_utm[, 2])

# Bandwidth matrix via cross-validated smoothed estimator
# Hscv is recommended for archaeological spatial data
# (unconstrained — allows asymmetric bandwidth)
cat("  Estimating bandwidth (Hscv)...\n")
set.seed(42)
H_cv <- ks::Hscv(X)
cat("  Bandwidth matrix H:\n")
cat(sprintf("    [%10.1f  %10.1f]\n", H_cv[1,1], H_cv[1,2]))
cat(sprintf("    [%10.1f  %10.1f]\n", H_cv[2,1], H_cv[2,2]))
cat("\n")

# Define evaluation grid at 1km resolution
# (KDE at 30m would require ~23.7M eval points — impractical)
# Resampled to 30m template after KDE evaluation
template_1km <- terra::aggregate(template_30m, fact = 33,
                                 fun = "mean", na.rm = TRUE)
grid_cells   <- which(!is.na(terra::values(template_1km)))
grid_xy      <- terra::xyFromCell(template_1km, grid_cells)

cat("  Evaluating KDE on 1km grid...\n")
cat("  Grid points:", nrow(grid_xy), "\n")

set.seed(42)
kde_fit <- ks::kde(X, H = H_cv, eval.points = grid_xy)

# Extract KDE values and clip negatives to zero
kde_vals <- kde_fit$estimate
kde_vals[kde_vals < 0] <- 0

cat("  KDE value range:", round(min(kde_vals), 8),
    "to", round(max(kde_vals), 8), "\n")
cat("  Non-zero cells:", sum(kde_vals > 0), "\n\n")

# ── 3. RASTERIZE KDE TO 30M TEMPLATE ────────────────────────

cat("--- Rasterizing KDE to 30m Template ---\n")

# Place KDE values into 1km raster
kde_1km <- terra::rast(template_1km)
terra::values(kde_1km) <- NA_real_
kde_1km[grid_cells]    <- kde_vals

# Resample to 30m template (bilinear — smooth gradient)
kde_30m <- terra::resample(kde_1km, template_30m,
                           method = "bilinear")
kde_30m <- terra::mask(kde_30m, boundary_vect)

# Clip any negative interpolation artefacts
kde_30m <- terra::ifel(kde_30m < 0, 0, kde_30m)

# Normalise to [0, 1] — used as sampling probability weights
kde_max <- terra::global(kde_30m, "max", na.rm = TRUE)[1, 1]
kde_norm <- kde_30m / kde_max

kde_stats <- terra::global(kde_norm, c("min","max","mean"),
                           na.rm = TRUE)
cat(sprintf("  KDE normalised range: %.4f to %.4f  (mean %.4f)\n",
            kde_stats[1,"min"], kde_stats[1,"max"],
            kde_stats[1,"mean"]))

# Save bias surface
out_bias <- file.path(OUT_BIAS, "bias_surface_kde_30m_utm44n.tif")
terra::writeRaster(kde_norm, out_bias,
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Saved: bias_surface_kde_30m_utm44n.tif\n\n")

# ── 4. SAMPLE BACKGROUND POINTS ─────────────────────────────

cat("--- Sampling Background Points ---\n\n")
cat("  Extent: Nagpur-Chandrapur district boundary ONLY\n")
cat("  Method: probability proportional to KDE bias surface\n")
cat("  Justification: Warton & Shepherd (2010)\n\n")

# Helper: sample N background points from kde_norm
sample_background <- function(n, label, seed = 42) {
  set.seed(seed)
  cat(sprintf("  Sampling N = %6d [%s]...", n, label))
  
  # terra::spatSample with prob=TRUE samples cells
  # proportional to raster values (bias-weighted sampling)
  bg_pts <- terra::spatSample(
    kde_norm,
    size   = n,
    method = "weights",   # weighted by raster values
    na.rm  = TRUE,
    as.points = TRUE,
    exhaustive = FALSE
  )
  
  # Convert to sf with UTM coordinates
  bg_sf <- sf::st_as_sf(bg_pts)
  bg_coords <- sf::st_coordinates(bg_sf)
  bg_df <- data.frame(
    bg_id     = seq_len(nrow(bg_coords)),
    presence  = 0L,
    UTM_X     = bg_coords[, 1],
    UTM_Y     = bg_coords[, 2]
  )
  bg_sf2 <- sf::st_as_sf(bg_df,
                         coords = c("UTM_X","UTM_Y"),
                         crs    = 32644,
                         remove = FALSE)
  
  cat(sprintf(" → %d points\n", nrow(bg_sf2)))
  return(bg_sf2)
}

# Primary background set (N = 10,000)
bg_10000 <- sample_background(N_BACKGROUND, "primary")

# Sensitivity analysis sets (Supplementary Table S5)
bg_1000  <- sample_background(1000,  "sensitivity")
bg_5000  <- sample_background(5000,  "sensitivity")
bg_20000 <- sample_background(20000, "sensitivity")

cat("\n")

# ── 5. VERIFY BACKGROUND SPATIAL DISTRIBUTION ────────────────

cat("--- Verifying Background Distribution ---\n\n")

# Check: background points should broadly track site density
# Extract KDE values at background and site locations
kde_at_bg <- terra::extract(kde_norm,
                            terra::vect(bg_10000))[, 2]
kde_at_sites <- terra::extract(kde_norm,
                               terra::vect(sites_sf))[, 2]

cat(sprintf("  KDE at background pts: mean=%.4f  SD=%.4f\n",
            mean(kde_at_bg,    na.rm=TRUE),
            sd(kde_at_bg,      na.rm=TRUE)))
cat(sprintf("  KDE at site locs:      mean=%.4f  SD=%.4f\n",
            mean(kde_at_sites, na.rm=TRUE),
            sd(kde_at_sites,   na.rm=TRUE)))

# Both should have similar mean KDE — confirms bias correction
# is working (background tracks survey effort like sites do)
ratio <- mean(kde_at_bg, na.rm=TRUE) /
  mean(kde_at_sites, na.rm=TRUE)
cat(sprintf("  KDE ratio (bg/sites):  %.3f\n", ratio))
if (ratio > 0.5 && ratio < 2.0) {
  cat("  ✓ Background KDE distribution plausible\n\n")
} else {
  cat("  ⚠ Large KDE ratio — check bias surface\n\n")
}

# ── 6. SAVE ALL BACKGROUND SETS ─────────────────────────────

cat("--- Saving Background Points ---\n\n")

bg_list <- list(
  list(bg_10000, "background_N10000.gpkg", "Primary"),
  list(bg_1000,  "background_N1000.gpkg",  "Sensitivity S5"),
  list(bg_5000,  "background_N5000.gpkg",  "Sensitivity S5"),
  list(bg_20000, "background_N20000.gpkg", "Sensitivity S5")
)

for (item in bg_list) {
  bg_sf    <- item[[1]]
  filename <- item[[2]]
  label    <- item[[3]]
  out_path <- file.path(OUT_BACKGROUND, filename)
  sf::st_write(bg_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("  ✓ %-30s %6d pts  [%s]\n",
              filename, nrow(bg_sf), label))
}

# ── 7. DIAGNOSTIC FIGURE ────────────────────────────────────

cat("\n--- Generating Diagnostic Figure ---\n")

png(file.path(OUT_FIG_SUPP, "S09_bias_surface_check.png"),
    width = 7200, height = 2400, res = 300)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))

# Panel 1: KDE bias surface
terra::plot(kde_norm,
            main = "KDE Bias Surface\n(normalised survey effort)",
            col  = viridisLite::plasma(100),
            axes = FALSE)
terra::plot(boundary_vect, add = TRUE,
            border = "white", lwd = 0.8)
terra::plot(terra::vect(sites_sf), add = TRUE,
            col = "white", pch = 16, cex = 0.3)

# Panel 2: Background points (N=10,000)
terra::plot(kde_norm,
            main = "Background Points (N=10,000)\nBias-weighted",
            col  = viridisLite::plasma(100),
            axes = FALSE)
terra::plot(boundary_vect, add = TRUE,
            border = "white", lwd = 0.8)
terra::plot(terra::vect(bg_10000), add = TRUE,
            col = "#00FF0044", pch = 16, cex = 0.05)

# Panel 3: Sites vs background overlay
terra::plot(template_30m,
            main = "Sites (red) vs Background (green)\nN=10,000",
            col  = "grey92", legend = FALSE, axes = FALSE)
terra::plot(boundary_vect, add = TRUE,
            border = "grey50", lwd = 0.8)
terra::plot(terra::vect(bg_10000), add = TRUE,
            col = "#00AA0066", pch = 16, cex = 0.05)
terra::plot(terra::vect(sites_sf), add = TRUE,
            col = "red", pch = 16, cex = 0.4)

dev.off()
cat("  ✓ S09_bias_surface_check.png\n\n")

# ── 8. SAVE SUMMARY FOR METHODS ─────────────────────────────

summary_df <- data.frame(
  set           = c("Primary","Sensitivity","Sensitivity","Sensitivity"),
  N             = c(10000, 1000, 5000, 20000),
  filename      = c("background_N10000.gpkg","background_N1000.gpkg",
                    "background_N5000.gpkg","background_N20000.gpkg"),
  justification = c("Warton & Shepherd (2010)",
                    "Supplementary Table S5",
                    "Supplementary Table S5",
                    "Supplementary Table S5"),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(OUT_BACKGROUND,
                    "S09_background_summary.csv"),
          row.names = FALSE)
cat("  ✓ S09_background_summary.csv\n\n")

# ── 9. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 09 COMPLETE — Summary\n")
cat("========================================\n")
cat("KDE method:      ks::kde() + Hscv() bandwidth\n")
cat("Extent:          District boundary only\n")
cat("Primary N:       10,000 (Warton & Shepherd 2010)\n")
cat("Sensitivity:     N = 1,000 / 5,000 / 10,000 / 20,000\n")
cat("  (for Supplementary Table S5)\n")
cat("\nFiles saved to:", OUT_BACKGROUND, "\n")
cat("\nNext: Run Script 10 — Spatial Block CV Design\n")
cat("========================================\n")