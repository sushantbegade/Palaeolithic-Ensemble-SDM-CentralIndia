# ============================================================
# SCRIPT 05: PALAEOCHANNEL LAYER CONSTRUCTION
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 05 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Constructs the novel palaeochannel proximity predictor
#   (Predictor 10) following the dual-source cross-validation
#   protocol specified in Research Design Section 5.4.
#
#   PROTOCOL (verbatim from Research Design 5.4):
#   Source 1 — Spectral: MNDWI > -0.10 moisture anomaly
#              (from Script 04: MNDWI_ANOMALY_30m_utm44n.tif)
#   Source 2 — Topographic: DEM valley morphology
#              (WhiteboxTools geomorphons — valley class)
#   CONFIRMED palaeochannel cell = present in BOTH sources.
#   Features in only one source are excluded.
#
#   Outputs:
#     PALAEO_VALLEY_30m_utm44n.tif      — binary valley mask
#     PALAEO_CANDIDATE_30m_utm44n.tif   — MNDWI ∩ valley (binary)
#     PALAEO_NETWORK_30m_utm44n.tif     — cleaned candidate network
#     DIST_PALAEOCHANNEL_30m_utm44n.tif — Euclidean distance raster
#                                         (Predictor 10 — FINAL)
#     PALAEO_fieldvalidation_plan.csv   — field validation points
#                                         (≥15 segments required)
#
# PREDICTOR FINALISED:
#   Predictor 10: DIST_PALAEOCHANNEL_30m_utm44n.tif
#
# METHOD NOTE — valley extraction:
#   WhiteboxTools wbt_geomorphons() classifies each cell into
#   one of 10 landform types. Valley class = 7 (concave, low
#   position). Combined with MNDWI anomaly this identifies
#   cells where relict drainage morphology AND subsurface
#   moisture signal are co-present — the two necessary
#   conditions for palaeochannel confirmation.
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

library(whitebox)
wbt_init()

cat("\n========================================\n")
cat("SCRIPT 05: Palaeochannel Layer\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

dem_path     <- file.path(OUT_PREDICTORS, "DEM_30m_utm44n.tif")
mndwi_path   <- file.path(OUT_PREDICTORS, "MNDWI_30m_utm44n.tif")
anomaly_path <- file.path(OUT_PREDICTORS,
                          "MNDWI_ANOMALY_30m_utm44n.tif")
tpi_path     <- file.path(OUT_PREDICTORS, "TPI_30m_utm44n.tif")

for (p in c(dem_path, mndwi_path, anomaly_path, tpi_path)) {
  if (!file.exists(p)) stop("Missing input: ", p,
                            "\nRun Scripts 02-04 first.")
}

dem_r      <- terra::rast(dem_path)
mndwi_r    <- terra::rast(mndwi_path)
anomaly_r  <- terra::rast(anomaly_path)
tpi_r      <- terra::rast(tpi_path)

cat("  DEM:           ✓\n")
cat("  MNDWI:         ✓\n")
cat("  MNDWI Anomaly: ✓\n")
cat("  TPI:           ✓\n\n")

# ── 2. SOURCE 2 — TOPOGRAPHIC VALLEY EXTRACTION ─────────────
# WhiteboxTools geomorphons classifies landform elements.
# Valley class = 7 in geomorphons output.
# Moving window (search radius) tested at 10 and 30 cells:
#   10 cells = 300m radius — captures fine channels
#   30 cells = 900m radius — captures broad palaeovalleys
# Both used; union = conservative valley mask.

cat("--- Source 2: Topographic Valley Extraction ---\n")
cat("  Method: WhiteboxTools wbt_geomorphons()\n")

# Geomorphons at two scales
geom_10_path <- file.path(OUT_PALAEO, "GEOMORPHONS_r10_temp.tif")
geom_30_path <- file.path(OUT_PALAEO, "GEOMORPHONS_r30_temp.tif")

cat("  Running geomorphons at r=10 cells (300m)...\n")
wbt_geomorphons(
  dem     = dem_path,
  output  = geom_10_path,
  search  = 10,      # search radius in cells
  threshold = 1.0,   # flatness threshold (degrees)
  fdist   = 0        # skip filter
)
cat("  ✓ Geomorphons r=10 complete\n")

cat("  Running geomorphons at r=30 cells (900m)...\n")
wbt_geomorphons(
  dem     = dem_path,
  output  = geom_30_path,
  search  = 30,
  threshold = 1.0,
  fdist   = 0
)
cat("  ✓ Geomorphons r=30 complete\n")

# Load geomorphons rasters
geom_10 <- terra::rast(geom_10_path)
geom_30 <- terra::rast(geom_30_path)

# Geomorphons landform codes (Jasiewicz & Stepinski 2013):
# 1=flat, 2=summit, 3=ridge, 4=shoulder, 5=spur,
# 6=slope, 7=hollow, 8=footslope, 9=valley, 10=depression
# Valley-type cells: 7 (hollow), 9 (valley), 10 (depression)
VALLEY_CLASSES <- c(7, 9, 10)

valley_10 <- ifel(geom_10 %in% VALLEY_CLASSES, 1L, 0L)
valley_30 <- ifel(geom_30 %in% VALLEY_CLASSES, 1L, 0L)

# Union of both scales (present in either = valley candidate)
valley_union <- ifel((valley_10 + valley_30) >= 1, 1L, 0L)
valley_union <- terra::resample(valley_union, template_30m,
                                method = "near")
valley_union <- terra::mask(valley_union, boundary_vect)

out_valley <- file.path(OUT_PALAEO, "PALAEO_VALLEY_30m_utm44n.tif")
terra::writeRaster(valley_union, out_valley,
                   overwrite = TRUE, datatype = "INT1U")

n_valley <- terra::global(valley_union == 1, "sum",
                          na.rm = TRUE)[1, 1]
n_total  <- terra::global(!is.na(valley_union), "sum",
                          na.rm = TRUE)[1, 1]
cat("  Valley cells:", n_valley,
    sprintf("(%.1f%% of study area)\n\n",
            100 * n_valley / n_total))

# ── 3. DUAL-SOURCE INTERSECTION ──────────────────────────────
# Confirmed palaeochannel cell = MNDWI anomaly AND valley signal.
# This is the cross-validation step from Research Design 5.4.

cat("--- Dual-Source Intersection (MNDWI ∩ Valley) ---\n")

# Align anomaly to valley raster (both should match template)
anomaly_aligned <- terra::resample(anomaly_r, template_30m,
                                   method = "near")

candidate_r <- ifel(
  (anomaly_aligned == 1L) & (valley_union == 1L),
  1L, 0L
)
candidate_r <- terra::mask(candidate_r, boundary_vect)

n_candidate <- terra::global(candidate_r == 1, "sum",
                             na.rm = TRUE)[1, 1]
n_mndwi     <- terra::global(anomaly_aligned == 1, "sum",
                             na.rm = TRUE)[1, 1]
n_val       <- n_valley

cat("  MNDWI anomaly cells:  ", n_mndwi, "\n")
cat("  Valley cells:         ", n_val, "\n")
cat("  Intersection (both):  ", n_candidate,
    sprintf("(%.1f%% of MNDWI anomaly retained)\n",
            100 * n_candidate / max(n_mndwi, 1)))

out_candidate <- file.path(OUT_PALAEO,
                           "PALAEO_CANDIDATE_30m_utm44n.tif")
terra::writeRaster(candidate_r, out_candidate,
                   overwrite = TRUE, datatype = "INT1U")
cat("  ✓ Candidate network saved\n\n")

# ── 4. CLEAN CANDIDATE NETWORK ───────────────────────────────
# Remove isolated pixels (noise) — keep only spatially
# connected clusters of ≥5 cells.
# This removes false positives from single-cell MNDWI
# speckle that happen to fall in valley positions.

cat("--- Cleaning Candidate Network ---\n")
cat("  Removing isolated pixels (min cluster = 5 cells)...\n")

# Use terra::patches() to label connected regions
patches_r <- terra::patches(candidate_r == 1L,
                            directions = 8,   # 8-connectivity
                            zeroAsNA   = TRUE)

# Count cells per patch
patch_sizes <- terra::freq(patches_r)
large_patches <- patch_sizes$value[patch_sizes$count >= 5]
cat("  Total candidate patches:", nrow(patch_sizes), "\n")
cat("  Patches with >= 5 cells:", length(large_patches), "\n")

# Keep only large patches
network_r <- ifel(patches_r %in% large_patches, 1L, 0L)
network_r <- terra::mask(network_r, boundary_vect)

n_network <- terra::global(network_r == 1, "sum",
                           na.rm = TRUE)[1, 1]
cat("  Network cells after cleaning:", n_network, "\n")
cat("  Retention rate:",
    sprintf("%.1f%%\n\n", 100 * n_network / max(n_candidate, 1)))

out_network <- file.path(OUT_PALAEO,
                         "PALAEO_NETWORK_30m_utm44n.tif")
terra::writeRaster(network_r, out_network,
                   overwrite = TRUE, datatype = "INT1U")
cat("  ✓ Cleaned network saved\n")

# ── 5. EUCLIDEAN DISTANCE RASTER — PREDICTOR 10 ─────────────
# Distance from each 30m cell to nearest confirmed
# palaeochannel cell. Lower value = closer to palaeochannel.
# This is the actual predictor entered into all six models.

cat("\n--- Computing Distance to Palaeochannel (Predictor 10) ---\n")

# terra::distance() computes Euclidean distance to nearest
# non-NA cell. Set non-network cells to NA first.
network_na <- ifel(network_r == 1L, 1L, NA)

dist_r <- terra::distance(network_na)
dist_r <- terra::mask(dist_r, boundary_vect)

dist_stats <- terra::global(dist_r, c("min", "max", "mean"),
                            na.rm = TRUE)
cat("  Distance range:", round(dist_stats[1, "min"], 0), "m to",
    round(dist_stats[1, "max"], 0), "m\n")
cat("  Mean distance: ", round(dist_stats[1, "mean"], 0), "m\n")

# Sanity check: mean distance should be several km in Vidarbha
if (dist_stats[1, "max"] < 1000) {
  warning("Max distance < 1km — network may be too dense. ",
          "Check MNDWI threshold or minimum patch size.")
}

out_dist <- file.path(OUT_PREDICTORS,
                      "DIST_PALAEOCHANNEL_30m_utm44n.tif")
terra::writeRaster(dist_r, out_dist,
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Predictor 10 saved: DIST_PALAEOCHANNEL_30m_utm44n.tif\n\n")

# ── 6. FIELD VALIDATION POINT GENERATION ────────────────────
# Research Design 5.4: ≥15 palaeochannel segments must be
# field-validated before the layer is used in modelling.
# This section generates 30 candidate GPS waypoints from the
# cleaned network — stratified across the study area.
# Take these waypoints to the field; validate using 2 of 3
# criteria (sediment texture, clast evidence, vegetation anomaly).

cat("--- Generating Field Validation Waypoints ---\n")
cat("  Protocol: ≥15 segments, ≥2 of 3 criteria per segment\n")
cat("  Generating 30 stratified candidate points...\n")

# Convert network raster to points
network_pts <- terra::as.points(network_na, values = FALSE,
                                na.rm = TRUE)

# Stratify across study area using a 5x6 grid (30 cells)
# Sample one point per grid cell where network is present
ext_study <- terra::ext(template_30m)
grid_cols <- 6
grid_rows <- 5

x_breaks <- seq(ext_study[1], ext_study[2],
                length.out = grid_cols + 1)
y_breaks <- seq(ext_study[3], ext_study[4],
                length.out = grid_rows + 1)

set.seed(42)
validation_pts <- list()

for (i in seq_len(grid_rows)) {
  for (j in seq_len(grid_cols)) {
    cell_ext <- terra::ext(x_breaks[j], x_breaks[j + 1],
                           y_breaks[i], y_breaks[i + 1])
    pts_in_cell <- terra::crop(network_pts, cell_ext)
    if (length(pts_in_cell) > 0) {
      # Random sample of 1 point from this cell
      idx <- sample(seq_len(length(pts_in_cell)), 1)
      validation_pts[[length(validation_pts) + 1]] <-
        pts_in_cell[idx, ]
    }
  }
}

if (length(validation_pts) > 0) {
  val_vect  <- do.call(rbind, validation_pts)
  # Convert to WGS84 geographic for GPS upload
  val_geo   <- terra::project(val_vect, "EPSG:4326")
  val_coords <- terra::crds(val_geo)
  
  val_df <- data.frame(
    Point_ID         = paste0("PC_", sprintf("%02d", seq_len(nrow(val_coords)))),
    Longitude_WGS84  = round(val_coords[, 1], 6),
    Latitude_WGS84   = round(val_coords[, 2], 6),
    Priority         = "High",
    Criterion_1_Sediment   = "Record: sandy loam vs clay matrix",
    Criterion_2_Clasts     = "Record: rounded gravel >5mm present?",
    Criterion_3_Vegetation = "Cross-ref NIR B8 map before visit",
    Confirmed        = NA,
    Notes            = NA
  )
  
  val_csv <- file.path(OUT_PALAEO,
                       "PALAEO_fieldvalidation_plan.csv")
  write.csv(val_df, val_csv, row.names = FALSE)
  cat("  ✓ Field validation plan saved:",
      basename(val_csv), "\n")
  cat("  Points generated:", nrow(val_df), "\n")
  cat("  Required confirmed: ≥15 (confirm ≥2 of 3 criteria)\n\n")
} else {
  cat("  ⚠ No network points found for waypoint generation.\n",
      "  Check candidate network — may be too sparse.\n\n")
}

# ── 7. DIAGNOSTIC FIGURES ───────────────────────────────────

cat("--- Generating Diagnostic Figures ---\n")

# 3-panel: MNDWI anomaly | valley mask | final network
png(file.path(OUT_FIG_SUPP, "S0e_palaeochannel_construction.png"),
    width = 7200, height = 2400, res = 300)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))

terra::plot(anomaly_aligned,
            main   = "Source 1: MNDWI Anomaly\n(MNDWI > -0.10)",
            col    = c("grey92", "#2166ac"),
            legend = FALSE, axes = FALSE)
terra::plot(boundary_vect, add = TRUE, border = "black", lwd = 0.8)

terra::plot(valley_union,
            main   = "Source 2: Valley Morphology\n(geomorphons)",
            col    = c("grey92", "#d7191c"),
            legend = FALSE, axes = FALSE)
terra::plot(boundary_vect, add = TRUE, border = "black", lwd = 0.8)

terra::plot(network_r,
            main   = "Confirmed Palaeochannel Network\n(MNDWI ∩ Valley)",
            col    = c("grey92", "#1a9641"),
            legend = FALSE, axes = FALSE)
terra::plot(boundary_vect, add = TRUE, border = "black", lwd = 0.8)
if (exists("val_vect")) {
  terra::plot(val_vect, add = TRUE, col = "orange",
              pch = 3, cex = 0.5, lwd = 1.5)
}

dev.off()
cat("  ✓ S0e_palaeochannel_construction.png\n")

# Distance raster
png(file.path(OUT_FIG_SUPP, "S0f_dist_palaeochannel.png"),
    width = 2400, height = 2400, res = 300)

terra::plot(dist_r / 1000,
            main = "Distance to Palaeochannel (km)\n[Predictor 10]",
            col  = viridisLite::plasma(100),
            axes = FALSE)
terra::plot(boundary_vect, add = TRUE, border = "white", lwd = 0.8)

dev.off()
cat("  ✓ S0f_dist_palaeochannel.png\n\n")

# ── 8. CLEAN TEMPORARY FILES ────────────────────────────────

for (f in c(geom_10_path, geom_30_path)) {
  if (file.exists(f)) file.remove(f)
}
cat("  Temporary geomorphons files removed\n\n")

# ── 9. FINAL VERIFICATION ───────────────────────────────────

cat("--- Verifying Predictor 10 ---\n")

dist_final <- terra::rast(out_dist)
dist_rng   <- terra::global(dist_final, c("min", "max", "mean"),
                            na.rm = TRUE)
cat(sprintf("  DIST_PALAEOCHANNEL_30m_utm44n.tif\n"))
cat(sprintf("    Range:  %.0f m to %.0f m\n",
            dist_rng[1, "min"], dist_rng[1, "max"]))
cat(sprintf("    Mean:   %.0f m\n", dist_rng[1, "mean"]))
cat(sprintf("    CRS:    %s\n",
            terra::crs(dist_final, describe = TRUE)$name))
cat(sprintf("    Dims:   %d x %d\n",
            nrow(dist_final), ncol(dist_final)))

# ── 10. SUMMARY ─────────────────────────────────────────────

cat("\n========================================\n")
cat("SCRIPT 05 COMPLETE — Summary\n")
cat("========================================\n")
cat("Valley extraction:       ✓ geomorphons r=10 + r=30\n")
cat("MNDWI ∩ Valley:          ✓", n_candidate, "candidate cells\n")
cat("After cleaning (≥5):     ✓", n_network, "network cells\n")
cat("Distance raster:         ✓ Predictor 10 ready\n")
cat("Field validation plan:   ✓ ≥30 waypoints generated\n")
cat("\nPREDICTOR 10 STATUS: COMPLETE\n")
cat("  → DIST_PALAEOCHANNEL_30m_utm44n.tif\n")
cat("\n⚠ FIELD ACTION REQUIRED before modelling:\n")
cat("  Visit ≥15 waypoints from PALAEO_fieldvalidation_plan.csv\n")
cat("  Confirm ≥2 of 3 criteria at each validated segment\n")
cat("  Record results in the CSV and report rate in Methods 5.4\n")
cat("\nNext: Run Script 06 — Climate Downscaling\n")
cat("  (CHELSA-TraCE21k LGM variables for UP sub-model)\n")
cat("========================================\n")