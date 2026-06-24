# ============================================================
# SCRIPT 07: PREDICTOR STACK ASSEMBLY + VIF SCREENING
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 07 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1. Computes Predictor 9:  Distance to Perennial River
#   2. Computes Predictor 11: Distance to Raw Material Source
#   3. Rasterizes Predictor 12: Geology (categorical)
#   4. Rasterizes Predictor 13: Geomorphology (categorical)
#   5. Assembles all 14 predictors into one raster stack
#   6. VIF screening — usdm::vifcor(), threshold VIF < 5
#      on continuous predictors only
#   7. Saves final retained stack + Table 2
#
# FIXES FROM v1:
#   - Field detection rewritten: explicit priority-ordered
#     candidate list for lithology, geology, geomorphology.
#     No longer relies on is.character() alone — catches
#     factor and numeric-coded categorical fields too.
#   - Lithology field: tries "lithologic", "formation",
#     "sld_name", "notation", "age" in order before fallback.
#   - Prints ALL unique values for every detected field so
#     user can verify correct field was chosen.
#   - geol_field / geom_field detection now length-safe.
#   - Raw material: uses correct lithology field, not the
#     centroid-code field (input_cent).
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 07: Predictor Stack + VIF\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. LOAD TEMPLATE AND BOUNDARY ───────────────────────────

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# ── HELPER: detect best class field from a data frame ────────
# Tries candidate field names in order, then falls back to
# the field with the most unique string-like values.

detect_class_field <- function(df, candidates, layer_name) {
  # Remove geometry column names
  all_fields <- names(df)[!names(df) %in% c("geom","geometry")]
  
  # Try candidate names first (case-insensitive)
  for (cand in candidates) {
    match <- all_fields[tolower(all_fields) == tolower(cand)]
    if (length(match) > 0 && length(unique(df[[match[1]]])) > 1) {
      cat(sprintf("    Field selected: '%s' (priority match)\n",
                  match[1]))
      return(match[1])
    }
  }
  
  # Fallback: field with most unique non-numeric values
  n_unique <- sapply(all_fields, function(f) {
    vals <- df[[f]]
    if (is.numeric(vals) && max(vals, na.rm=TRUE) > 1000) return(0)
    length(unique(vals[!is.na(vals)]))
  })
  best <- all_fields[which.max(n_unique)]
  
  if (length(best) == 0 || is.na(best)) {
    stop("Cannot detect class field for ", layer_name,
         ". Fields: ", paste(all_fields, collapse=", "))
  }
  cat(sprintf("    Field selected: '%s' (fallback — max unique)\n",
              best))
  return(best)
}

# ── 2. PREDICTOR 9 — DISTANCE TO PERENNIAL RIVER ────────────

cat("--- Predictor 9: Distance to Perennial River ---\n")

rivers_sf   <- sf::st_read(file.path(OUT_PREDICTORS,
                                     "rivers_utm44n.gpkg"),
                           quiet = TRUE)
rivers_vect <- terra::vect(rivers_sf)
cat("  River features:", nrow(rivers_sf), "\n")

rivers_rast <- terra::rasterize(rivers_vect, template_30m,
                                field = 1, background = NA)
dist_river  <- terra::distance(rivers_rast)
dist_river  <- terra::mask(dist_river, boundary_vect)

rng_rv <- terra::global(dist_river, c("min","max","mean"),
                        na.rm = TRUE)
cat(sprintf("  Range: %.0f m to %.0f m  (mean: %.0f m)\n",
            rng_rv[1,"min"], rng_rv[1,"max"], rng_rv[1,"mean"]))

terra::writeRaster(dist_river,
                   file.path(OUT_PREDICTORS,
                             "DIST_RIVER_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Saved: DIST_RIVER_30m_utm44n.tif\n\n")

# ── 3. PREDICTOR 11 — DISTANCE TO RAW MATERIAL ──────────────

cat("--- Predictor 11: Distance to Raw Material Source ---\n\n")

lithology_sf <- sf::st_read(file.path(OUT_PREDICTORS,
                                      "lithology_utm44n.gpkg"),
                            quiet = TRUE)
cat("  Lithology features:", nrow(lithology_sf), "\n")

# Detect lithology class field — priority order for GSI shapefiles
lith_candidates <- c("lithologic", "lithology", "lith_type",
                     "rock_type", "formation", "sld_name",
                     "notation", "age", "supergroup",
                     "group_name", "description", "litho")

lith_field <- detect_class_field(lithology_sf,
                                 lith_candidates,
                                 "Lithology")

lith_vals <- sort(unique(as.character(
  lithology_sf[[lith_field]])))
cat(sprintf("  Unique values in '%s' (%d):\n",
            lith_field, length(lith_vals)))
for (v in lith_vals) cat("    -", v, "\n")
cat("\n")

# Raw material keywords — knappable rocks for Palaeolithic
rawmat_keywords <- c(
  "quartzite", "quartz", "chert", "flint", "siliceous",
  "basalt", "trap", "deccan", "agate", "jasper",
  "silicified", "hornfels", "lydite", "archaean",
  "crystalline", "metamorphic", "schist", "gneiss",
  "granite", "sandstone", "limestone", "shale"
)

rawmat_pattern <- paste(rawmat_keywords, collapse = "|")
rawmat_sf <- lithology_sf[
  grepl(rawmat_pattern,
        as.character(lithology_sf[[lith_field]]),
        ignore.case = TRUE), ]

cat("  Keywords matched:", nrow(rawmat_sf), "polygons\n")

if (nrow(rawmat_sf) == 0) {
  cat("  ⚠ No keyword matches — checking if values are codes.\n")
  cat("  → Using ALL lithology polygons (entire study area).\n")
  cat("  → Distance will be near-zero everywhere.\n")
  cat("  → ACTION: After script, inspect lithology_utm44n.gpkg\n")
  cat("    in QGIS, identify the correct field for rock types,\n")
  cat("    update lith_field and rawmat_keywords, re-run.\n\n")
  rawmat_sf <- lithology_sf
} else {
  matched_cls <- sort(unique(as.character(rawmat_sf[[lith_field]])))
  cat("  Matched classes:\n")
  for (cls in matched_cls) cat("    ✓", cls, "\n")
  cat("\n")
}

rawmat_vect <- terra::vect(rawmat_sf)
rawmat_rast <- terra::rasterize(rawmat_vect, template_30m,
                                field = 1, background = NA)
dist_rawmat <- terra::distance(rawmat_rast)
dist_rawmat <- terra::mask(dist_rawmat, boundary_vect)

rng_rm <- terra::global(dist_rawmat, c("min","max","mean"),
                        na.rm = TRUE)
cat(sprintf("  Range: %.0f m to %.0f m  (mean: %.0f m)\n",
            rng_rm[1,"min"], rng_rm[1,"max"], rng_rm[1,"mean"]))

if (rng_rm[1,"max"] < 500) {
  cat("  ⚠ Max distance < 500m — raw material likely covers\n")
  cat("    most of study area. Check field and keywords.\n")
}

terra::writeRaster(dist_rawmat,
                   file.path(OUT_PREDICTORS,
                             "DIST_RAWMAT_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Saved: DIST_RAWMAT_30m_utm44n.tif\n\n")

# ── 4. PREDICTOR 12 — GEOLOGY (categorical raster) ──────────

cat("--- Predictor 12: Geology ---\n")

geol_sf <- sf::st_read(file.path(OUT_PREDICTORS,
                                 "geology_utm44n.gpkg"),
                       quiet = TRUE)
cat("  Geology features:", nrow(geol_sf), "\n")
cat("  Fields:", paste(names(geol_sf)[
  !names(geol_sf) %in% c("geom","geometry")],
  collapse=", "), "\n")

geol_candidates <- c("geology", "geol_type", "rock_type",
                     "lithology", "formation", "description",
                     "age", "era", "period", "supergroup",
                     "group_name", "notation", "sld_name",
                     "name", "type", "label", "class")

geol_field <- detect_class_field(geol_sf, geol_candidates,
                                 "Geology")

geol_vals <- sort(unique(as.character(geol_sf[[geol_field]])))
cat(sprintf("  Unique classes in '%s' (%d):\n",
            geol_field, length(geol_vals)))
for (v in geol_vals) cat("    -", v, "\n")
cat("\n")

# Encode as integer
geol_sf$RASTER_CODE <- as.integer(factor(
  as.character(geol_sf[[geol_field]])))

geol_lookup <- data.frame(
  code  = sort(unique(geol_sf$RASTER_CODE)),
  label = levels(factor(as.character(geol_sf[[geol_field]])))
)
write.csv(geol_lookup,
          file.path(OUT_PREDICTORS, "GEOLOGY_lookup.csv"),
          row.names = FALSE)
cat("  Lookup table saved: GEOLOGY_lookup.csv\n")

geol_vect <- terra::vect(geol_sf)
geol_rast <- terra::rasterize(geol_vect, template_30m,
                              field = "RASTER_CODE",
                              background = NA)
geol_rast <- terra::mask(geol_rast, boundary_vect)

rng_g <- terra::global(geol_rast, c("min","max"), na.rm = TRUE)
cat(sprintf("  Code range: %d to %d\n",
            as.integer(rng_g[1,1]), as.integer(rng_g[1,2])))

terra::writeRaster(geol_rast,
                   file.path(OUT_PREDICTORS,
                             "GEOLOGY_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "INT2S")
cat("  ✓ Saved: GEOLOGY_30m_utm44n.tif\n\n")

# ── 5. PREDICTOR 13 — GEOMORPHOLOGY (categorical raster) ────

cat("--- Predictor 13: Geomorphology ---\n")

geom_sf <- sf::st_read(file.path(OUT_PREDICTORS,
                                 "geomorphology_utm44n.gpkg"),
                       quiet = TRUE)
cat("  Geomorphology features:", nrow(geom_sf), "\n")
cat("  Fields:", paste(names(geom_sf)[
  !names(geom_sf) %in% c("geom","geometry")],
  collapse=", "), "\n")

geom_candidates <- c("geomorphology", "geomorph", "landform",
                     "land_form", "terrain", "form_type",
                     "type", "class", "category",
                     "description", "name", "label",
                     "unit", "morphology", "morph_type")

geom_field <- detect_class_field(geom_sf, geom_candidates,
                                 "Geomorphology")

geom_vals <- sort(unique(as.character(geom_sf[[geom_field]])))
cat(sprintf("  Unique classes in '%s' (%d):\n",
            geom_field, length(geom_vals)))
# Print first 20 only if many
for (v in head(geom_vals, 20)) cat("    -", v, "\n")
if (length(geom_vals) > 20) {
  cat("    ... and", length(geom_vals) - 20, "more\n")
}
cat("\n")

# Encode as integer
geom_sf$RASTER_CODE <- as.integer(factor(
  as.character(geom_sf[[geom_field]])))

geom_lookup <- data.frame(
  code  = sort(unique(geom_sf$RASTER_CODE)),
  label = levels(factor(as.character(geom_sf[[geom_field]])))
)
write.csv(geom_lookup,
          file.path(OUT_PREDICTORS, "GEOMORPHOLOGY_lookup.csv"),
          row.names = FALSE)
cat("  Lookup table saved: GEOMORPHOLOGY_lookup.csv\n")

geom_vect <- terra::vect(geom_sf)
geom_rast <- terra::rasterize(geom_vect, template_30m,
                              field = "RASTER_CODE",
                              background = NA)
geom_rast <- terra::mask(geom_rast, boundary_vect)

rng_gm <- terra::global(geom_rast, c("min","max"), na.rm = TRUE)
cat(sprintf("  Code range: %d to %d\n",
            as.integer(rng_gm[1,1]), as.integer(rng_gm[1,2])))

terra::writeRaster(geom_rast,
                   file.path(OUT_PREDICTORS,
                             "GEOMORPHOLOGY_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "INT2S")
cat("  ✓ Saved: GEOMORPHOLOGY_30m_utm44n.tif\n\n")

# ── 6. ASSEMBLE FULL 14-PREDICTOR STACK ─────────────────────

cat("--- Assembling 14-Predictor Stack ---\n\n")

predictor_files <- data.frame(
  num      = 1:14,
  name     = c("Elevation","Slope","Aspect","TRI","TPI",
               "Plan_Curvature","HAND","Flow_Accum_log10",
               "Dist_River","Dist_Palaeochannel",
               "Dist_RawMat","Geology","Geomorphology","NDVI"),
  filename = c("DEM_30m_utm44n.tif",
               "SLOPE_30m_utm44n.tif",
               "ASPECT_30m_utm44n.tif",
               "TRI_30m_utm44n.tif",
               "TPI_30m_utm44n.tif",
               "PLANCURV_30m_utm44n.tif",
               "HAND_30m_utm44n.tif",
               "FLOWACC_LOG10_30m_utm44n.tif",
               "DIST_RIVER_30m_utm44n.tif",
               "DIST_PALAEOCHANNEL_30m_utm44n.tif",
               "DIST_RAWMAT_30m_utm44n.tif",
               "GEOLOGY_30m_utm44n.tif",
               "GEOMORPHOLOGY_30m_utm44n.tif",
               "NDVI_30m_utm44n.tif"),
  type     = c(rep("continuous",11),
               "categorical","categorical","continuous"),
  stringsAsFactors = FALSE
)

rast_list   <- list()
all_present <- TRUE

for (i in seq_len(nrow(predictor_files))) {
  fpath <- file.path(OUT_PREDICTORS, predictor_files$filename[i])
  if (!file.exists(fpath)) {
    cat(sprintf("  ✗ P%02d %-22s MISSING\n",
                i, predictor_files$name[i]))
    all_present <- FALSE
  } else {
    r      <- terra::rast(fpath)
    names(r) <- predictor_files$name[i]
    rng    <- terra::global(r, c("min","max"), na.rm = TRUE)
    cat(sprintf("  ✓ P%02d %-22s [%9.2f, %9.2f]  %s\n",
                i, predictor_files$name[i],
                round(rng[1,1],2), round(rng[1,2],2),
                predictor_files$type[i]))
    rast_list[[predictor_files$name[i]]] <- r
  }
}

if (!all_present) stop("Missing predictors — check output above.")

pred_stack <- terra::rast(rast_list)
cat("\n  Stack:", terra::nlyr(pred_stack), "layers x",
    nrow(pred_stack), "rows x", ncol(pred_stack), "cols\n")

terra::writeRaster(pred_stack,
                   file.path(OUT_PREDICTORS,
                             "PREDICTOR_STACK_14_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Full 14-layer stack saved\n\n")

# ── 7. VIF SCREENING — CONTINUOUS PREDICTORS ONLY ───────────

cat("--- VIF Screening (usdm::vifcor, VIF < 5) ---\n\n")
cat("  Categorical predictors excluded from VIF.\n\n")

cont_names <- predictor_files$name[
  predictor_files$type == "continuous"]

cat("  Continuous predictors entered (n =",
    length(cont_names), "):\n")
for (n in cont_names) cat("   ", n, "\n")
cat("\n")

# Sample 10,000 random non-NA cells
set.seed(42)
cont_stack  <- terra::subset(pred_stack, cont_names)
sample_vals <- terra::spatSample(cont_stack,
                                 size   = 10000,
                                 method = "random",
                                 na.rm  = TRUE,
                                 as.df  = TRUE)
cat("  Sample size:", nrow(sample_vals), "cells\n\n")

cat("  Running usdm::vifcor()...\n\n")
vif_result <- usdm::vifcor(sample_vals, th = VIF_THRESHOLD)
print(vif_result)

retained_vars <- vif_result@results$Variables
excluded_vars <- cont_names[!cont_names %in% retained_vars]

cat("\n  ─────────────────────────────────────────\n")
cat("  RETAINED (VIF <", VIF_THRESHOLD, "):",
    length(retained_vars), "continuous\n")
for (v in retained_vars) cat("   ✓", v, "\n")
cat("\n  EXCLUDED (VIF >=", VIF_THRESHOLD, "):",
    length(excluded_vars), "\n")
if (length(excluded_vars) > 0) {
  for (v in excluded_vars) cat("   ✗", v, "\n")
} else {
  cat("   None\n")
}
cat("  ─────────────────────────────────────────\n\n")

# ── 8. FINAL RETAINED STACK ─────────────────────────────────

final_names <- c(retained_vars, "Geology", "Geomorphology")
final_stack <- terra::subset(pred_stack, final_names)

terra::writeRaster(final_stack,
                   file.path(OUT_PREDICTORS,
                             "PREDICTOR_STACK_FINAL_30m_utm44n.tif"),
                   overwrite = TRUE, datatype = "FLT4S")
cat("  ✓ Final retained stack saved:",
    terra::nlyr(final_stack), "layers\n\n")

# ── 9. TABLE 2 ──────────────────────────────────────────────

cat("--- Table 2: Predictor Summary ---\n\n")

# Get VIF scores from result object
vif_df <- as.data.frame(vif_result@results)
names(vif_df) <- c("name","VIF")

table2 <- predictor_files
table2$VIF_score <- NA_real_
table2$min_val   <- NA_real_
table2$max_val   <- NA_real_
table2$status    <- NA_character_

for (i in seq_len(nrow(table2))) {
  nm <- table2$name[i]
  r  <- terra::rast(file.path(OUT_PREDICTORS,
                              table2$filename[i]))
  rg <- terra::global(r, c("min","max"), na.rm = TRUE)
  table2$min_val[i] <- round(rg[1,1], 3)
  table2$max_val[i] <- round(rg[1,2], 3)
  
  if (table2$type[i] == "categorical") {
    table2$status[i]    <- "Categorical — not in VIF"
  } else if (nm %in% vif_df$name) {
    table2$VIF_score[i] <- round(vif_df$VIF[vif_df$name==nm], 3)
    table2$status[i]    <- ifelse(nm %in% retained_vars,
                                  "Retained", "Excluded (VIF >= 5)")
  } else {
    # Removed iteratively before final table
    table2$status[i]    <- "Excluded (VIF >= 5)"
  }
}

write.csv(table2,
          file.path(OUT_TABLES, "Table2_predictors_VIF.csv"),
          row.names = FALSE)
cat("  ✓ Saved: Table2_predictors_VIF.csv\n\n")

cat(sprintf("  %-2s  %-22s %-12s %6s  %s\n",
            "#","Predictor","Type","VIF","Status"))
cat("  ", paste(rep("-",68), collapse=""), "\n")
for (i in seq_len(nrow(table2))) {
  cat(sprintf("  %2d  %-22s %-12s %6s  %s\n",
              table2$num[i],
              table2$name[i],
              table2$type[i],
              ifelse(is.na(table2$VIF_score[i]), "—",
                     sprintf("%.2f", table2$VIF_score[i])),
              table2$status[i]))
}

# ── 10. SAVE RETAINED NAMES FOR DOWNSTREAM SCRIPTS ──────────

saveRDS(retained_vars,
        file.path(OUT_PREDICTORS, "retained_continuous_vars.rds"))
saveRDS(final_names,
        file.path(OUT_PREDICTORS, "final_predictor_names.rds"))
cat("\n  ✓ Retained names saved as .rds for Scripts 08-25\n\n")

# ── 11. SUMMARY ─────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 07 COMPLETE — Summary\n")
cat("========================================\n")
cat("Dist to River:       ✓ Predictor 9\n")
cat("Dist to Raw Mat:     ✓ Predictor 11\n")
cat("Geology raster:      ✓ Predictor 12\n")
cat("Geomorphology rast:  ✓ Predictor 13\n")
cat("14-layer stack:      ✓\n")
cat("VIF threshold:       <", VIF_THRESHOLD, "\n")
cat("Retained cont:      ", length(retained_vars), "\n")
cat("Excluded:           ", length(excluded_vars), "\n")
cat("Final stack layers: ", terra::nlyr(final_stack), "\n")
if (rng_rm[1,"max"] < 500) {
  cat("\n⚠ ACTION NEEDED — Raw Material Predictor:\n")
  cat("  Distance range near-zero — lithology field\n")
  cat("  '", lith_field, "' may not contain rock type names.\n",
      sep="")
  cat("  Open lithology_utm44n.gpkg in QGIS,\n")
  cat("  identify correct field, update lith_field\n")
  cat("  and rawmat_keywords in this script, re-run.\n")
}
cat("\nNext: Run Script 08 — Site Thinning (spThin)\n")
cat("========================================\n")