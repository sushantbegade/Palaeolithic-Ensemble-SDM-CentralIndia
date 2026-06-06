# ============================================================
# SCRIPT 10: SPATIAL BLOCK CROSS-VALIDATION DESIGN
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 10 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Designs the 5-fold spatial block cross-validation used
#   in all six algorithm scripts (11-16) and ensemble (17).
#   Spatial block CV eliminates autocorrelation artefacts
#   by ensuring training and test sites are spatially
#   separated (Roberts et al. 2017; Valavi et al. 2019).
#
#   METHOD:
#   1. Spatial autocorrelation range estimated via
#      blockCV::cv_spatial_autocor() on combined presence
#      + background dataset — range determines block size
#   2. 5-fold spatial blocks designed via blockCV::cv_spatial()
#   3. Both presence (thinned sites) and background points
#      assigned to the same block structure
#   4. Block assignments saved — used identically in all
#      algorithm scripts (consistent CV across algorithms)
#
#   FALLBACK: if autocorrelation range < 10km or > 150km,
#   a default block size of 50km is applied (appropriate
#   for a ~21,300 km² study area with 190 thinned sites).
#
# OUTPUTS (to data_processed/cv_blocks/):
#   cv_block_assignments.rds    — fold IDs for sites + bg
#   cv_block_stats.csv          — block size, range, fold N
#   Supplementary Figure S1 panel data
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp")

cat("\n========================================\n")
cat("SCRIPT 10: Spatial Block CV Design\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# Thinned pooled sites (190 sites)
sites_sf <- sf::st_read(file.path(OUT_SITES,
                                  "sites_thinned_pooled.gpkg"),
                        quiet = TRUE)
sites_sf$presence <- 1L
cat("  Thinned sites:", nrow(sites_sf), "\n")

# Primary background (N=10,000)
bg_sf <- sf::st_read(file.path(OUT_BACKGROUND,
                               "background_N10000.gpkg"),
                     quiet = TRUE)
bg_sf$presence <- 0L
cat("  Background pts:", nrow(bg_sf), "\n\n")

# ── 2. COMBINE PRESENCE + BACKGROUND ────────────────────────

cat("--- Combining Presence + Background ---\n\n")

# Align columns — keep only presence + geometry
sites_slim <- sites_sf[, "presence"]
bg_slim    <- bg_sf[, "presence"]

combined_sf <- rbind(sites_slim, bg_slim)
cat("  Combined dataset:", nrow(combined_sf),
    "records (", sum(combined_sf$presence),
    "presences +", sum(!combined_sf$presence),
    "background)\n\n")

# ── 3. SPATIAL AUTOCORRELATION RANGE ────────────────────────

cat("--- Spatial Autocorrelation Range ---\n")
cat("  Package: blockCV::cv_spatial_autocor()\n\n")

autocor_range_m <- tryCatch({
  
  autocor <- blockCV::cv_spatial_autocor(
    x      = combined_sf,
    column = "presence",
    plot   = FALSE
  )
  
  range_val <- autocor$range
  cat(sprintf("  Autocorrelation range: %.1f m (%.1f km)\n",
              range_val, range_val / 1000))
  range_val
  
}, error = function(e) {
  cat("  cv_spatial_autocor() failed:", conditionMessage(e), "\n")
  cat("  → Using default block size: 50,000 m (50 km)\n")
  50000
})

# Apply sensible bounds for study area
# Study area ~21,300 km² → reasonable block range 20-100km
if (autocor_range_m < 10000) {
  cat("  ⚠ Range < 10km — applying minimum 20km\n")
  autocor_range_m <- 20000
} else if (autocor_range_m > 150000) {
  cat("  ⚠ Range > 150km — applying maximum 100km\n")
  autocor_range_m <- 100000
}

block_size_m  <- autocor_range_m
block_size_km <- block_size_m / 1000
cat(sprintf("  Block size used: %.0f m (%.1f km)\n\n",
            block_size_m, block_size_km))

# ── 4. DESIGN 5-FOLD SPATIAL BLOCKS ─────────────────────────

cat("--- Designing 5-Fold Spatial Blocks ---\n\n")

set.seed(42)

cv_blocks <- tryCatch({
  
  blockCV::cv_spatial(
    x         = combined_sf,
    column    = "presence",
    r         = template_30m,
    k         = 5,
    size      = block_size_m,
    selection = "random",
    iteration = 100,   # tries to balance folds
    seed      = 42,
    plot      = FALSE,
    report    = FALSE
  )
  
}, error = function(e) {
  cat("  cv_spatial() error:", conditionMessage(e), "\n")
  cat("  → Trying without raster argument...\n")
  
  blockCV::cv_spatial(
    x         = combined_sf,
    column    = "presence",
    k         = 5,
    size      = block_size_m,
    selection = "random",
    iteration = 100,
    seed      = 42,
    plot      = FALSE,
    report    = FALSE
  )
})

cat("  ✓ CV block design complete\n\n")

# ── 5. EXTRACT FOLD ASSIGNMENTS ──────────────────────────────

cat("--- Fold Assignments ---\n\n")

# blockCV 3.x stores fold IDs in cv_blocks$folds_ids
fold_ids <- cv_blocks$folds_ids

cat("  Total records with fold ID:", length(fold_ids), "\n")
cat("  Fold distribution:\n")
fold_table <- table(fold_ids)
for (f in names(fold_table)) {
  n_pres <- sum(combined_sf$presence[fold_ids == as.integer(f)])
  n_bg   <- sum(!combined_sf$presence[fold_ids == as.integer(f)])
  cat(sprintf("    Fold %s: %4d total  (%3d presences + %4d background)\n",
              f, fold_table[f], n_pres, n_bg))
}
cat("\n")

# Separate fold assignments for sites and background
n_sites <- nrow(sites_sf)
site_folds <- fold_ids[seq_len(n_sites)]
bg_folds   <- fold_ids[(n_sites + 1):length(fold_ids)]

cat("  Site fold distribution:\n")
for (f in sort(unique(site_folds))) {
  cat(sprintf("    Fold %d: %d sites\n", f,
              sum(site_folds == f)))
}

# Check: every fold must have at least 3 presence sites
min_pres_per_fold <- min(table(site_folds))
if (min_pres_per_fold < 3) {
  warning("Fold with < 3 presence sites — consider larger block size")
  cat("  ⚠ Min presences per fold:", min_pres_per_fold, "\n")
} else {
  cat("  ✓ Min presences per fold:", min_pres_per_fold, "(adequate)\n")
}
cat("\n")

# ── 6. SAVE CV DESIGN ────────────────────────────────────────

cat("--- Saving CV Design ---\n\n")

# Primary output: list with fold assignments for sites and bg
cv_design <- list(
  block_size_m    = block_size_m,
  block_size_km   = block_size_km,
  n_folds         = 5L,
  autocor_range_m = autocor_range_m,
  n_sites         = n_sites,
  n_background    = nrow(bg_sf),
  site_folds      = site_folds,      # fold ID per thinned site
  bg_folds        = bg_folds,        # fold ID per background pt
  combined_folds  = fold_ids,        # full combined vector
  cv_blocks_obj   = cv_blocks        # full blockCV object
)

saveRDS(cv_design,
        file.path(OUT_CV, "cv_block_assignments.rds"))
cat("  ✓ cv_block_assignments.rds\n")

# Summary stats CSV
stats_df <- data.frame(
  parameter    = c("block_size_m","block_size_km",
                   "autocor_range_m","n_folds",
                   "n_sites","n_background",
                   "min_sites_per_fold",
                   "max_sites_per_fold"),
  value        = c(block_size_m, block_size_km,
                   autocor_range_m, 5,
                   n_sites, nrow(bg_sf),
                   min(table(site_folds)),
                   max(table(site_folds))),
  stringsAsFactors = FALSE
)

write.csv(stats_df,
          file.path(OUT_CV, "cv_block_stats.csv"),
          row.names = FALSE)
cat("  ✓ cv_block_stats.csv\n\n")

# Add fold IDs to site and background sf objects and save
sites_sf$fold_id <- site_folds
sf::st_write(sites_sf,
             file.path(OUT_CV, "sites_with_folds.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
cat("  ✓ sites_with_folds.gpkg\n")

bg_sf$fold_id <- bg_folds
sf::st_write(bg_sf,
             file.path(OUT_CV, "background_with_folds.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
cat("  ✓ background_with_folds.gpkg\n\n")

# ── 7. DIAGNOSTIC FIGURE ────────────────────────────────────

cat("--- Generating Diagnostic Figure ---\n")

fold_colours <- c("#e41a1c","#377eb8","#4daf4a",
                  "#984ea3","#ff7f00")

png(file.path(OUT_FIG_SUPP, "S10_cv_block_design.png"),
    width = 4800, height = 2400, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

# Panel 1: CV block map (sites coloured by fold)
terra::plot(template_30m, col = "grey95",
            legend = FALSE, axes = FALSE,
            main = sprintf(
              "5-Fold Spatial Block CV\nBlock size: %.0f km",
              block_size_km))
terra::plot(boundary_vect, add = TRUE,
            border = "grey40", lwd = 0.8)
for (f in 1:5) {
  pts_f <- sites_sf[site_folds == f, ]
  if (nrow(pts_f) > 0) {
    terra::plot(terra::vect(pts_f), add = TRUE,
                col = fold_colours[f], pch = 16, cex = 0.5)
  }
}
legend("bottomright",
       legend = paste("Fold", 1:5,
                      paste0("(n=", table(site_folds), ")")),
       col    = fold_colours,
       pch    = 16, cex = 0.6, bty = "n")

# Panel 2: Background coloured by fold
terra::plot(template_30m, col = "grey95",
            legend = FALSE, axes = FALSE,
            main = "Background Points by Fold\n(N=10,000)")
terra::plot(boundary_vect, add = TRUE,
            border = "grey40", lwd = 0.8)
for (f in 1:5) {
  bg_f <- bg_sf[bg_folds == f, ]
  if (nrow(bg_f) > 0) {
    terra::plot(terra::vect(bg_f), add = TRUE,
                col = paste0(
                  substr(fold_colours[f], 1, 7), "44"),
                pch = 16, cex = 0.04)
  }
}

dev.off()
cat("  ✓ S10_cv_block_design.png\n\n")

# ── 8. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 10 COMPLETE — Summary\n")
cat("========================================\n")
cat(sprintf("Autocorrelation range: %.1f km\n",
            autocor_range_m / 1000))
cat(sprintf("Block size used:       %.1f km\n",
            block_size_km))
cat("Folds:                 5\n")
cat(sprintf("Sites per fold:        %d – %d\n",
            min(table(site_folds)),
            max(table(site_folds))))
cat("\nOutputs saved to:\n ", OUT_CV, "\n")
cat("\nCV design used in ALL algorithm scripts (11-16)\n")
cat("Load with: cv_design <- readRDS(file.path(OUT_CV,\n")
cat("             'cv_block_assignments.rds'))\n")
cat("\nNext: Run Script 11 — MaxEnt Model\n")
cat("========================================\n")