# ============================================================
# SCRIPT 06: CLIMATE DOWNSCALING — CHELSA-TraCE21k
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 06 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Processes CHELSA-TraCE21k palaeoclimate variables for use
#   in the UP sub-model ONLY. Reprojects, resamples to master
#   template, applies correct unit conversion, and saves.
#
#   PALAEOCLIMATE INCLUSION RULES (Research Design 5.3):
#   ┌─────────────────────────────────────────────────────┐
#   │ Pooled model:   NO palaeoclimate                    │
#   │ LP sub-model:   NO palaeoclimate                    │
#   │   Acheulian >200 ka — beyond CHELSA-TraCE21k range  │
#   │ MP sub-model:   NO palaeoclimate                    │
#   │   MIS3 ~40 ka — CHELSA-TraCE21k extends to ~22 ka   │
#   │ UP sub-model:   CHELSA-TraCE21k LGM time slice      │
#   │   BIO1  (mean annual temperature, converted to °C)  │
#   │   BIO12 (annual precipitation, mm)                  │
#   │   BIO15 (precipitation seasonality, CV)             │
#   └─────────────────────────────────────────────────────┘
#
# FILES (confirmed present and correct):
#   CHELSA_TraCE21k_bio01_-200_V.1.0.tif
#   CHELSA_TraCE21k_bio12_-200_V.1.0.tif
#   CHELSA_TraCE21k_bio15_-200_V.1.0.tif
#   Time slice: -200 (centennial index = 20,000 BP = LGM)
#   Confirmed by raster metadata: activity_id = last_glacial_period
#
# BIO01 UNIT CONVERSION (investigated and confirmed):
#   Raster metadata: variable_unit = K, Scale = 0.1
#   Terra reads stored integer values and applies Scale=0.1,
#   yielding raw terra values in KELVIN (e.g., 300.4–303.1 K).
#   Correct conversion to °C: BIO01_celsius = raw - 273.15
#   Result for study area: 27.25–29.95°C (mean 28.80°C)
#   Spatial SD = 0.42°C — meaningful gradient, predictor retained.
#   WRONG approach (used in earlier diagnostic): dividing by 10
#   gave 30.2°C, misidentified as °C×10 — incorrect.
#
# BIO12 and BIO15:
#   No unit conversion needed. Values are in mm and CV
#   respectively, loaded correctly by terra as-is.
#
# REFERENCE: Karger et al. (2023) Climate of the Past 19:439-456
# ============================================================
# HOW TO USE:
#   Run after Script 05. Script 01 sourced automatically.
#   Outputs used only in Script 20 (sub-models).
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 06: Climate Downscaling\n")
cat("========================================\n\n")

# ── 1. LOAD TEMPLATE AND BOUNDARY ───────────────────────────

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

cat("Template CRS:",
    terra::crs(template_30m, describe = TRUE)$name, "\n\n")

# ── 2. SOURCE FILES ──────────────────────────────────────────
# Hardcoded — filenames confirmed present and verified correct.

lgm_files <- list(
  "01" = file.path(PATH_CLIMATE,
                   "CHELSA_TraCE21k_bio01_-200_V.1.0.tif"),
  "12" = file.path(PATH_CLIMATE,
                   "CHELSA_TraCE21k_bio12_-200_V.1.0.tif"),
  "15" = file.path(PATH_CLIMATE,
                   "CHELSA_TraCE21k_bio15_-200_V.1.0.tif")
)

out_names <- list(
  "01" = "CHELSA_BIO1_LGM_30m_utm44n.tif",
  "12" = "CHELSA_BIO12_LGM_30m_utm44n.tif",
  "15" = "CHELSA_BIO15_LGM_30m_utm44n.tif"
)

cat("--- Verifying Source Files ---\n")
for (bio in c("01","12","15")) {
  if (file.exists(lgm_files[[bio]])) {
    cat(sprintf("  ✓ BIO%s: %s\n", bio,
                basename(lgm_files[[bio]])))
  } else {
    stop("Missing: ", lgm_files[[bio]])
  }
}
cat("\n")

# ── 3. PROCESSING PARAMETERS ────────────────────────────────

# Crop extent — WGS84, slightly larger than study area
study_bbox_geo <- terra::ext(78.0, 80.3, 19.2, 21.9)

# BIO01 Kelvin offset (confirmed from raster metadata)
KELVIN_OFFSET  <- 273.15

# ── 4. PROCESS EACH VARIABLE ────────────────────────────────

cat("--- Processing Climate Variables ---\n\n")

for (bio in c("01", "12", "15")) {
  
  cat(sprintf("  [BIO%s]\n", bio))
  cat("    Source:", basename(lgm_files[[bio]]), "\n")
  
  # Load
  r <- terra::rast(lgm_files[[bio]])
  cat("    CRS (raw):  ",
      terra::crs(r, describe = TRUE)$name, "\n")
  cat("    Res (raw):  ",
      paste(round(terra::res(r), 4), collapse = " x "), "\n")
  
  # Raw value inspection before conversion
  r_crop_raw <- terra::crop(r, study_bbox_geo)
  raw_stats  <- terra::global(r_crop_raw,
                              c("min","max","mean","sd"),
                              na.rm = TRUE)
  cat(sprintf("    Raw values: %.2f to %.2f  (mean %.2f, SD %.4f)\n",
              raw_stats[1,"min"], raw_stats[1,"max"],
              raw_stats[1,"mean"], raw_stats[1,"sd"]))
  
  # Step 1: Crop to study bbox (avoids global reproject)
  r <- terra::crop(r, study_bbox_geo)
  
  # Step 2: Unit conversion for BIO01 only
  # BIO01: stored in Kelvin (terra applies Scale=0.1 on read,
  # yielding values in K such as 300.4–303.1).
  # Convert to °C by subtracting 273.15.
  # BIO12 (mm) and BIO15 (CV): no conversion needed.
  if (bio == "01") {
    r <- r - KELVIN_OFFSET
    conv_stats <- terra::global(r, c("min","max","mean","sd"),
                                na.rm = TRUE)
    cat(sprintf("    After K→°C: %.2f to %.2f°C  (mean %.2f°C, SD %.4f°C)\n",
                conv_stats[1,"min"], conv_stats[1,"max"],
                conv_stats[1,"mean"], conv_stats[1,"sd"]))
    if (conv_stats[1,"sd"] >= 0.1) {
      cat("    ✓ Spatial variation adequate (SD >= 0.1°C)\n")
    } else {
      cat("    ⚠ Low spatial variation (SD < 0.1°C) — review\n")
    }
  }
  
  if (bio == "12") {
    cat(sprintf("    Units: mm  range %.0f to %.0f mm\n",
                raw_stats[1,"min"], raw_stats[1,"max"]))
  }
  
  if (bio == "15") {
    cat(sprintf("    Units: CV  range %.1f to %.1f\n",
                raw_stats[1,"min"], raw_stats[1,"max"]))
  }
  
  # Step 3: Reproject to UTM 44N
  r <- terra::project(r, terra::crs(template_30m),
                      method = "bilinear")
  
  # Step 4: Resample to 30m master template
  # CHELSA ~1km → 30m = spatial interpolation only.
  # Bilinear appropriate for smooth climate gradient fields.
  r <- terra::resample(r, template_30m, method = "bilinear")
  r <- terra::mask(r, boundary_vect)
  
  # Final value range after all processing
  final_stats <- terra::global(r, c("min","max","mean","sd"),
                               na.rm = TRUE)
  cat(sprintf("    Final (30m UTM): %.2f to %.2f  (SD %.4f)\n",
              final_stats[1,"min"], final_stats[1,"max"],
              final_stats[1,"sd"]))
  
  # Save
  out_path <- file.path(OUT_CLIMATE_PRO, out_names[[bio]])
  terra::writeRaster(r, out_path,
                     overwrite = TRUE, datatype = "FLT4S")
  cat(sprintf("    ✓ Saved: %s\n\n", out_names[[bio]]))
}

# ── 5. VERIFY ALL OUTPUTS ────────────────────────────────────

cat("--- Verifying Outputs ---\n\n")

all_ok <- TRUE
for (bio in c("01","12","15")) {
  fpath <- file.path(OUT_CLIMATE_PRO, out_names[[bio]])
  if (file.exists(fpath)) {
    r    <- terra::rast(fpath)
    stat <- terra::global(r, c("min","max","mean"),
                          na.rm = TRUE)
    unit_label <- switch(bio,
                         "01" = "°C",
                         "12" = "mm",
                         "15" = "CV")
    cat(sprintf("  ✓ BIO%s: %.2f to %.2f %s  (mean %.2f)\n",
                bio,
                stat[1,"min"], stat[1,"max"],
                unit_label, stat[1,"mean"]))
  } else {
    cat("  ✗ MISSING: BIO", bio, "\n")
    all_ok <- FALSE
  }
}

# ── 6. SPATIAL VARIATION SUMMARY ────────────────────────────

cat("\n--- Spatial Variation Summary (for Script 20) ---\n\n")

var_df <- data.frame(
  variable = character(),
  unit     = character(),
  min      = numeric(),
  max      = numeric(),
  mean     = numeric(),
  sd       = numeric(),
  range    = numeric(),
  use_UP   = character(),
  stringsAsFactors = FALSE
)

units_map <- list("01"="°C", "12"="mm", "15"="CV")

for (bio in c("01","12","15")) {
  r    <- terra::rast(file.path(OUT_CLIMATE_PRO, out_names[[bio]]))
  stat <- terra::global(r, c("min","max","mean","sd"),
                        na.rm = TRUE)
  rng  <- stat[1,"max"] - stat[1,"min"]
  use  <- if (stat[1,"sd"] >= 0.1) "YES" else
    "LOW VARIATION — review"
  cat(sprintf("  BIO%s (%s): range=%.3f  SD=%.4f  → %s\n",
              bio, units_map[[bio]], rng, stat[1,"sd"], use))
  var_df <- rbind(var_df, data.frame(
    variable = paste0("BIO", bio),
    unit     = units_map[[bio]],
    min      = round(stat[1,"min"],  3),
    max      = round(stat[1,"max"],  3),
    mean     = round(stat[1,"mean"], 3),
    sd       = round(stat[1,"sd"],   4),
    range    = round(rng, 3),
    use_UP   = use,
    stringsAsFactors = FALSE
  ))
}

write.csv(var_df,
          file.path(OUT_CLIMATE_PRO,
                    "S06_climate_variation_summary.csv"),
          row.names = FALSE)
cat("\n  ✓ Saved: S06_climate_variation_summary.csv\n\n")

# ── 7. DIAGNOSTIC FIGURE ────────────────────────────────────

cat("--- Generating Diagnostic Figure ---\n")

png(file.path(OUT_FIG_SUPP, "S0g_climate_LGM_check.png"),
    width = 7200, height = 2400, res = 300)
par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))

plot_titles <- list(
  "01" = "BIO1: Mean Ann. Temp (°C)\nLGM ~20 ka [K - 273.15]",
  "12" = "BIO12: Annual Precip (mm)\nLGM ~20 ka",
  "15" = "BIO15: Precip Seasonality (CV)\nLGM ~20 ka"
)

for (bio in c("01","12","15")) {
  r <- terra::rast(file.path(OUT_CLIMATE_PRO, out_names[[bio]]))
  terra::plot(r,
              main = plot_titles[[bio]],
              col  = viridisLite::viridis(100),
              axes = FALSE)
  terra::plot(boundary_vect, add = TRUE,
              border = "white", lwd = 0.8)
}

dev.off()
cat("  ✓ S0g_climate_LGM_check.png\n\n")

# ── 8. METHODS NOTE ─────────────────────────────────────────

methods_note <- paste0(
  "CHELSA-TraCE21k palaeoclimate variables (BIO1, BIO12, BIO15; ",
  "time slice -200, ~20,000 BP, confirmed as last_glacial_period ",
  "from raster metadata; Karger et al. 2023) were processed for ",
  "the Upper Palaeolithic sub-model. BIO1 raw values are stored ",
  "in Kelvin (raster metadata: variable_unit = K, Scale = 0.1; ",
  "terra applies the scale factor on read, yielding values in K). ",
  "BIO1 was converted to degrees Celsius by subtracting 273.15 ",
  "before use (BIO1_celsius = BIO1_kelvin - 273.15), yielding a ",
  "study-area range of 27.25-29.95 degrees C (SD = 0.42 degrees C). ",
  "BIO12 (annual precipitation, mm) and BIO15 (precipitation ",
  "seasonality, CV) required no unit conversion. All three variables ",
  "were reprojected to WGS84 UTM Zone 44N (EPSG:32644) and ",
  "resampled from 30 arc-second (~1 km) to 30 m using bilinear ",
  "interpolation. Palaeoclimate variables are included as spatially ",
  "variable indicators of long-term climatic gradients relevant to ",
  "the LGM occupation period, acknowledging that temporal alignment ",
  "is approximate (Karger et al. 2023). No palaeoclimate variables ",
  "are included in the pooled model, the Lower Palaeolithic sub-model ",
  "(Acheulian occupations >200 ka predate CHELSA-TraCE21k coverage), ",
  "or the Middle Palaeolithic sub-model (MIS3 ~40 ka lies outside ",
  "CHELSA-TraCE21k temporal extent)."
)

writeLines(methods_note,
           file.path(OUT_CLIMATE_PRO, "climate_methods_note.txt"))
cat("  ✓ Methods note saved: climate_methods_note.txt\n\n")

# ── 9. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 06 COMPLETE — Summary\n")
cat("========================================\n")
cat("Source: CHELSA-TraCE21k time slice -200\n")
cat("        (last_glacial_period, ~20,000 BP)\n")
cat("BIO01:  Kelvin → °C  (subtract 273.15)\n")
cat("BIO12:  mm  (no conversion)\n")
cat("BIO15:  CV  (no conversion)\n")
cat("\nFinal values (°C or native units):\n")
for (bio in c("01","12","15")) {
  r    <- terra::rast(file.path(OUT_CLIMATE_PRO, out_names[[bio]]))
  stat <- terra::global(r, c("min","max"), na.rm = TRUE)
  cat(sprintf("  BIO%s: %.2f to %.2f %s\n",
              bio, stat[1,1], stat[1,2], units_map[[bio]]))
}
cat("\nUSAGE:\n")
cat("  Pooled:  NO climate\n")
cat("  LP:      NO climate\n")
cat("  MP:      NO climate\n")
cat("  UP:      BIO1 + BIO12 + BIO15 (all retained)\n")
cat("\nOutputs: ", OUT_CLIMATE_PRO, "\n")
cat("\nNext: Run Script 07 — Predictor Stack + VIF\n")
cat("========================================\n")