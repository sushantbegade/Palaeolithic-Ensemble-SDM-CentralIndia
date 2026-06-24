# ============================================================
# SCRIPT 16b PART 2: Fix Dist_RawMat with selective keywords
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 16b2 of 25 (inserted — not in original plan)
# ============================================================
# WHAT THIS SCRIPT DOES:
# FIELD: lithologic (42 values — confirmed correct)
# STRATEGY: restrict to HIGH-QUALITY knappable rocks only
# JUSTIFICATION: Palaeolithic raw material selection was
#   quality-specific — fine-grained siliceous/hard rocks
#   preferred over coarse/weathered equivalents
#   (Petraglia & Korisettar 1998; Paddayya et al. 2002)
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)
set.seed(42)

cat("\n========================================\n")
cat("SCRIPT 16b PART 2: Fix Dist_RawMat\n")
cat("========================================\n\n")

# ── 1. LOAD LITHOLOGY ────────────────────────────────────────

lithology_sf <- sf::st_read(
  file.path(OUT_PREDICTORS, "lithology_utm44n.gpkg"),
  quiet = TRUE)

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# ── 2. SELECTIVE RAW MATERIAL KEYWORDS ───────────────────────
# TIER 1 — Primary knappable: siliceous, fine-grained, hard
# Basalt EXCLUDED: weathers to clay in Vidarbha, poor knapping
# Granite EXCLUDED: coarse-grained, not preferred
# Schist EXCLUDED: foliated, poor conchoidal fracture
# INCLUDED: chert, quartzite, quartz vein, meta-rhyolite

RAWMAT_SELECTIVE <- c(
  # Best knapping stone in Deccan context
  "chert",
  "quartzite",
  "quartz vein",
  "silicified",
  # Fine-grained volcanic — knappable
  "meta rhyolite",
  "meta basalt",      # metamorphosed = harder than fresh basalt
  # Fine-grained siliceous sedimentary
  "cherty limestone",
  "ferruginous sandstone",
  "fine grained sandstone"
)

rawmat_pattern <- paste(RAWMAT_SELECTIVE, collapse = "|")

rawmat_sf <- lithology_sf[
  grepl(rawmat_pattern,
        as.character(lithology_sf[["lithologic"]]),
        ignore.case = TRUE), ]

cat("--- Selective Keyword Match ---\n")
cat("  Field used: lithologic\n")
cat("  Keywords:", paste(RAWMAT_SELECTIVE, collapse=", "), "\n\n")
cat("  Matched polygons:", nrow(rawmat_sf), "\n")

matched_vals <- sort(unique(
  as.character(rawmat_sf[["lithologic"]])))
cat("  Matched lithologic classes:\n")
for (v in matched_vals) cat("    ✓", v, "\n")
cat("\n")

# ── 3. COMPUTE NEW DISTANCE RASTER ───────────────────────────

cat("--- Computing Selective DIST_RAWMAT ---\n\n")

if (nrow(rawmat_sf) == 0) {
  stop("No polygons matched — check keywords vs lithologic values")
}

rawmat_vect <- terra::vect(rawmat_sf)
rawmat_rast <- terra::rasterize(rawmat_vect, template_30m,
                                field      = 1,
                                background = NA)
dist_rawmat_new <- terra::distance(rawmat_rast)
dist_rawmat_new <- terra::mask(dist_rawmat_new, boundary_vect)

rng_new <- terra::global(dist_rawmat_new,
                         c("min","max","mean","sd"),
                         na.rm = TRUE)

cat("  Old DIST_RAWMAT: min=0m  max=4888m  mean=36m  SD=234m\n")
cat(sprintf("  New DIST_RAWMAT: min=%.0fm  max=%.0fm  mean=%.0fm  SD=%.0fm\n",
            rng_new[1,"min"], rng_new[1,"max"],
            rng_new[1,"mean"], rng_new[1,"sd"]))
cat("\n")

# ── 4. DIAGNOSTIC CHECK ──────────────────────────────────────

cat("--- Diagnostic ---\n")

if (rng_new[1,"max"] < 1000) {
  cat("  ⚠ Max still < 1km — selective rocks also widespread\n")
  cat("  Consider restricting to CHERT + QUARTZITE + QUARTZ VEIN only\n")
  cat("  These are rarest + highest quality in Vidarbha\n\n")
  IMPROVEMENT <- "marginal"
} else if (rng_new[1,"mean"] < 500) {
  cat("  ⚠ Mean < 500m — high-quality rocks still common\n")
  cat("  Predictor has more variation than before but still limited\n\n")
  IMPROVEMENT <- "moderate"
} else {
  cat("  ✓ Distance range substantially improved\n")
  cat("  Predictor now has genuine spatial gradient\n\n")
  IMPROVEMENT <- "substantial"
}

# ── 5. ULTRA-SELECTIVE FALLBACK (if needed) ──────────────────
# If mean still < 500m, restrict to RAREST rocks only

if (IMPROVEMENT %in% c("marginal", "moderate")) {
  
  cat("--- Ultra-Selective: Chert + Quartzite + Quartz Vein only ---\n\n")
  
  RAWMAT_ULTRA <- c("chert", "quartzite", "quartz vein",
                    "quartz vein/silicified zone")
  
  ultra_pattern <- paste(RAWMAT_ULTRA, collapse = "|")
  
  rawmat_ultra_sf <- lithology_sf[
    grepl(ultra_pattern,
          as.character(lithology_sf[["lithologic"]]),
          ignore.case = TRUE), ]
  
  cat("  Ultra-selective matched:", nrow(rawmat_ultra_sf), "polygons\n")
  ultra_vals <- sort(unique(
    as.character(rawmat_ultra_sf[["lithologic"]])))
  for (v in ultra_vals) cat("    ✓", v, "\n")
  cat("\n")
  
  if (nrow(rawmat_ultra_sf) > 0) {
    rawmat_ultra_vect <- terra::vect(rawmat_ultra_sf)
    rawmat_ultra_rast <- terra::rasterize(rawmat_ultra_vect,
                                          template_30m,
                                          field      = 1,
                                          background = NA)
    dist_ultra <- terra::distance(rawmat_ultra_rast)
    dist_ultra <- terra::mask(dist_ultra, boundary_vect)
    
    rng_ultra <- terra::global(dist_ultra,
                               c("min","max","mean","sd"),
                               na.rm = TRUE)
    cat(sprintf("  Ultra DIST_RAWMAT: min=%.0fm max=%.0fm mean=%.0fm SD=%.0fm\n",
                rng_ultra[1,"min"], rng_ultra[1,"max"],
                rng_ultra[1,"mean"], rng_ultra[1,"sd"]))
    
    # Use ultra if mean > 1km (real spatial gradient)
    if (rng_ultra[1,"mean"] > 1000) {
      dist_rawmat_new <- dist_ultra
      rng_new         <- rng_ultra
      matched_vals    <- ultra_vals
      cat("  ✓ Ultra-selective version adopted (mean > 1km)\n\n")
      IMPROVEMENT <- "ultra_adopted"
    } else {
      cat("  Mean still < 1km — chert/quartzite also widespread\n")
      cat("  Predictor dead due to genuine Vidarbha geology\n\n")
      IMPROVEMENT <- "genuinely_dead"
    }
  }
}

# ── 6. SAVE BEST VERSION ─────────────────────────────────────

cat("--- Saving Updated DIST_RAWMAT ---\n\n")

terra::writeRaster(
  dist_rawmat_new,
  file.path(OUT_PREDICTORS, "DIST_RAWMAT_30m_utm44n.tif"),
  overwrite = TRUE,
  datatype  = "FLT4S"
)

cat(sprintf("  ✓ Saved. Final stats: mean=%.0fm SD=%.0fm max=%.0fm\n\n",
            rng_new[1,"mean"], rng_new[1,"sd"], rng_new[1,"max"]))

# Save record of what was used
rawmat_record <- data.frame(
  field_used       = "lithologic",
  classes_matched  = paste(matched_vals, collapse="; "),
  n_polygons       = nrow(rawmat_sf),
  mean_dist_m      = round(rng_new[1,"mean"], 0),
  max_dist_m       = round(rng_new[1,"max"],  0),
  sd_dist_m        = round(rng_new[1,"sd"],   0),
  improvement      = IMPROVEMENT,
  stringsAsFactors = FALSE
)
write.csv(rawmat_record,
          file.path(OUT_PREDICTORS, "RAWMAT_selection_record.csv"),
          row.names = FALSE)
cat("  ✓ RAWMAT_selection_record.csv saved\n\n")

# ── 7. REGENERATE PREDICTOR STACK ────────────────────────────

cat("--- Regenerating PREDICTOR_STACK_FINAL ---\n\n")

final_names <- readRDS(file.path(OUT_PREDICTORS,
                                 "final_predictor_names.rds"))

all_pred_files <- list(
  "Elevation"          = "DEM_30m_utm44n.tif",
  "Slope"              = "SLOPE_30m_utm44n.tif",
  "Aspect"             = "ASPECT_30m_utm44n.tif",
  "TRI"                = "TRI_30m_utm44n.tif",
  "TPI"                = "TPI_30m_utm44n.tif",
  "Plan_Curvature"     = "PLANCURV_30m_utm44n.tif",
  "HAND"               = "HAND_30m_utm44n.tif",
  "Flow_Accum_log10"   = "FLOWACC_LOG10_30m_utm44n.tif",
  "Dist_River"         = "DIST_RIVER_30m_utm44n.tif",
  "Dist_Palaeochannel" = "DIST_PALAEOCHANNEL_30m_utm44n.tif",
  "Dist_RawMat"        = "DIST_RAWMAT_30m_utm44n.tif",
  "Geology"            = "GEOLOGY_30m_utm44n.tif",
  "Geomorphology"      = "GEOMORPHOLOGY_30m_utm44n.tif",
  "NDVI"               = "NDVI_30m_utm44n.tif"
)

rast_list <- lapply(final_names, function(nm) {
  r <- terra::rast(file.path(OUT_PREDICTORS,
                             all_pred_files[[nm]]))
  names(r) <- nm
  return(r)
})

new_stack <- terra::rast(rast_list)

terra::writeRaster(
  new_stack,
  file.path(OUT_PREDICTORS,
            "PREDICTOR_STACK_FINAL_30m_utm44n.tif"),
  overwrite = TRUE,
  datatype  = "FLT4S"
)

cat("  ✓ PREDICTOR_STACK_FINAL regenerated\n")
cat("  Layers:", terra::nlyr(new_stack), "\n\n")

# ── 8. FINAL DECISION ────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 16b PART 2 COMPLETE\n")
cat("========================================\n\n")

cat("RAWMAT improvement status:", IMPROVEMENT, "\n\n")

if (IMPROVEMENT == "genuinely_dead") {
  cat("CONCLUSION: Dist_RawMat dead due to real geology.\n")
  cat("  Chert/quartzite/quartz vein widespread in Vidarbha.\n")
  cat("  Archaean basement + metamorphic belt = siliceous\n")
  cat("  rocks ubiquitous. Hominins had no procurement\n")
  cat("  constraint — consistent with dense site distribution.\n\n")
  cat("ACTION: Proceed to Script 17 WITHOUT re-running 11–16.\n")
  cat("  Report in Results 6.2 and Limitations 7.5:\n")
  cat("  'Dist_RawMat received near-zero importance across\n")
  cat("   all tree-based algorithms, consistent with the\n")
  cat("   ubiquitous distribution of knappable siliceous\n")
  cat("   lithologies across the Archaean crystalline basement\n")
  cat("   of the Wainganga-Wardha basin. Raw material\n")
  cat("   availability did not constitute a binding constraint\n")
  cat("   on site location in this geological context,\n")
  cat("   consistent with the high density of Palaeolithic\n")
  cat("   sites throughout the region.'\n\n")
} else {
  cat("CONCLUSION: Dist_RawMat now has real spatial gradient.\n\n")
  cat("DECISION NEEDED — re-run Scripts 11–16?\n\n")
  cat("  IF mean_dist > 5000m: YES re-run — major change\n")
  cat("  IF mean_dist 1000–5000m: BORDERLINE\n")
  cat("    Check if Dist_RawMat now in top-5 by quick RF test\n")
  cat("  IF mean_dist < 1000m: NO — marginal change\n")
  cat("    Proceed to Script 17, note in Limitations\n\n")
  cat(sprintf("  Your mean_dist = %.0fm\n", rng_new[1,"mean"]))
  cat("  → Paste new stats here for decision\n")
}

cat("\nNext: Paste output here → me give final decision\n")
cat("========================================\n")