# ============================================================
# SCRIPT 06: CLIMATE DOWNSCALING — CHELSA-TraCE21k
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 06 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Processes CHELSA-TraCE21k palaeoclimate variables for use
#   in cultural period sub-models ONLY (NOT the pooled model).
#
#   PALAEOCLIMATE INCLUSION RULES (Research Design 5.3):
#   ┌─────────────────────────────────────────────────────┐
#   │ LP sub-model : NO palaeoclimate                     │
#   │   Acheulian >200 ka — beyond CHELSA-TraCE21k range  │
#   │ MP sub-model : NO palaeoclimate                     │
#   │   MIS3 ~40 ka — outside CHELSA-TraCE21k extent      │
#   │   (TraCE21k only extends to ~22 ka)                 │
#   │ UP sub-model : CHELSA-TraCE21k LGM (~21 ka)         │
#   │   BIO1 (mean annual temperature)                    │
#   │   BIO12 (annual precipitation)                      │
#   │   BIO15 (precipitation seasonality)                 │
#   └─────────────────────────────────────────────────────┘
#
#   CHELSA-TraCE21k specs:
#   - Resolution: 30 arc-second (~1 km at equator)
#   - LGM time slice: -210 (= 21,000 BP) in file naming
#   - Reference: Karger et al. (2023) Climate of the Past
#   - Source: chelsa-climate.org/chelsa-trace21k/
#
#   DOWNSCALING NOTE (for Methods 5.3):
#   CHELSA-TraCE21k variables are resampled from ~1 km to 30 m
#   using bilinear interpolation. As stated in the Research
#   Design framing: "Palaeoclimate variables are included as
#   spatially variable indicators of long-term climatic
#   gradients relevant to the cultural period in question,
#   acknowledging that temporal alignment is approximate."
#   Fine-scale 30m variation in these layers reflects
#   interpolation, not genuine palaeoclimate resolution.
#
# OUTPUTS (to data_processed/climate/):
#   CHELSA_BIO1_LGM_30m_utm44n.tif   — mean annual temp (°C×10)
#   CHELSA_BIO12_LGM_30m_utm44n.tif  — annual precip (mm)
#   CHELSA_BIO15_LGM_30m_utm44n.tif  — precip seasonality (CV)
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

cat("Template loaded — CRS:",
    terra::crs(template_30m, describe = TRUE)$name, "\n\n")

# ── 2. DEFINE TARGET VARIABLES ──────────────────────────────

# LGM time slice code in CHELSA-TraCE21k file naming
# -210 = 21,000 BP (Last Glacial Maximum)
LGM_SLICE <- "-210"

# Variables needed
target_vars <- data.frame(
  bio_code  = c("bio01", "bio12", "bio15"),
  bio_num   = c("01",    "12",    "15"),
  label     = c("Mean Annual Temperature (x10 °C)",
                "Annual Precipitation (mm)",
                "Precipitation Seasonality (CV)"),
  out_name  = c("CHELSA_BIO1_LGM_30m_utm44n.tif",
                "CHELSA_BIO12_LGM_30m_utm44n.tif",
                "CHELSA_BIO15_LGM_30m_utm44n.tif"),
  stringsAsFactors = FALSE
)

cat("Target variables for UP sub-model (LGM ~21 ka):\n")
for (i in seq_len(nrow(target_vars))) {
  cat(sprintf("  BIO%s: %s\n",
              target_vars$bio_num[i], target_vars$label[i]))
}
cat("\n")

# ── 3. FIND CHELSA FILES IN CLIMATE FOLDER ──────────────────
# CHELSA-TraCE21k naming convention:
#   CHELSA_TraCE21k_bio01_-210_V1.0.tif
# Files may be nested in subdirectories — scan recursively.

cat("--- Locating CHELSA-TraCE21k Files ---\n")
cat("  Searching in:", PATH_CLIMATE, "\n\n")

all_climate_files <- list.files(PATH_CLIMATE,
                                pattern    = "\\.tif$|\\.nc$",
                                full.names = TRUE,
                                recursive  = TRUE)

cat("  Total climate files found:", length(all_climate_files), "\n")

if (length(all_climate_files) == 0) {
  stop(
    "No climate files found in: ", PATH_CLIMATE,
    "\nDownload CHELSA-TraCE21k LGM files from:",
    "\n  https://chelsa-climate.org/chelsa-trace21k/",
    "\nRequired files:",
    "\n  CHELSA_TraCE21k_bio01_-210_V1.0.tif",
    "\n  CHELSA_TraCE21k_bio12_-210_V1.0.tif",
    "\n  CHELSA_TraCE21k_bio15_-210_V1.0.tif"
  )
}

# Match each target variable
find_climate_file <- function(bio_code, time_slice, file_list) {
  # Try exact CHELSA-TraCE21k naming pattern
  patterns <- c(
    paste0("TraCE21k_", bio_code, "_", time_slice),
    paste0("trace21k_", bio_code, "_", time_slice),
    paste0(bio_code, "_", time_slice),
    paste0(bio_code, ".*", gsub("-", "", time_slice)),
    paste0("LGM.*", bio_code),
    paste0(bio_code, ".*LGM")
  )
  for (pat in patterns) {
    matches <- grep(pat, file_list,
                    value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  return(NULL)
}

found_files <- character(nrow(target_vars))
any_missing <- FALSE

for (i in seq_len(nrow(target_vars))) {
  f <- find_climate_file(target_vars$bio_code[i],
                         LGM_SLICE,
                         all_climate_files)
  if (is.null(f)) {
    cat(sprintf("  ✗ BIO%s: NOT FOUND\n", target_vars$bio_num[i]))
    any_missing <- TRUE
  } else {
    cat(sprintf("  ✓ BIO%s: %s\n",
                target_vars$bio_num[i], basename(f)))
    found_files[i] <- f
  }
}

if (any_missing) {
  cat("\n  ⚠ Some CHELSA files not found.\n")
  cat("  Files present in climate folder:\n")
  for (f in head(basename(all_climate_files), 20)) {
    cat("    ", f, "\n")
  }
  cat("\n  Attempting to use any available BIO files...\n")
  
  # Fallback: try to find BIO files without time slice filter
  for (i in seq_len(nrow(target_vars))) {
    if (found_files[i] == "") {
      fallback <- grep(target_vars$bio_code[i],
                       all_climate_files,
                       value = TRUE, ignore.case = TRUE)
      if (length(fallback) > 0) {
        cat(sprintf("  Fallback BIO%s: %s\n",
                    target_vars$bio_num[i],
                    basename(fallback[1])))
        found_files[i] <- fallback[1]
      }
    }
  }
}

cat("\n")

# ── 4. PROCESS EACH VARIABLE ────────────────────────────────

cat("--- Processing Climate Variables ---\n\n")

# Study area bounding box in WGS84 for initial crop
# (faster than reprojecting full global raster)
study_bbox_geo <- terra::ext(78.0, 80.3, 19.2, 21.9)

processed_vars <- list()

for (i in seq_len(nrow(target_vars))) {
  
  bio_num  <- target_vars$bio_num[i]
  label    <- target_vars$label[i]
  out_name <- target_vars$out_name[i]
  src_file <- found_files[i]
  
  cat(sprintf("  [BIO%s — %s]\n", bio_num, label))
  
  if (src_file == "" || !file.exists(src_file)) {
    cat(sprintf("    ✗ Skipped — file not found\n\n"))
    next
  }
  
  cat("    Source:", basename(src_file), "\n")
  
  # Load
  r_raw <- terra::rast(src_file)
  cat("    Original CRS:       ",
      terra::crs(r_raw, describe = TRUE)$name, "\n")
  cat("    Original resolution:",
      paste(round(terra::res(r_raw), 4), collapse = " x "), "\n")
  
  # Step 1: Crop to study area bounding box in native CRS
  # (avoids reprojecting entire global raster)
  r_crop <- terra::crop(r_raw, study_bbox_geo)
  cat("    ✓ Cropped to study area extent\n")
  
  # Step 2: Reproject to UTM 44N
  r_utm <- terra::project(r_crop,
                          y      = terra::crs(template_30m),
                          method = "bilinear")
  cat("    ✓ Reprojected to UTM 44N\n")
  
  # Step 3: Resample to 30m template
  # NOTE: CHELSA ~1km → 30m is interpolation, not new information.
  # Bilinear produces smooth gradients appropriate for climate fields.
  r_aligned <- terra::resample(r_utm, template_30m,
                               method = "bilinear")
  r_aligned  <- terra::mask(r_aligned, boundary_vect)
  cat("    ✓ Resampled to 30m template\n")
  
  # Step 4: Value checks
  rng <- terra::global(r_aligned, c("min", "max", "mean"),
                       na.rm = TRUE)
  
  # BIO1: stored as °C × 10 in CHELSA — check reasonable range
  if (bio_num == "01") {
    if (rng[1, "min"] < -500 || rng[1, "max"] > 500) {
      cat("    ⚠ BIO1 range unexpected:",
          round(rng[1, "min"], 1), "to",
          round(rng[1, "max"], 1), "\n")
    } else {
      cat("    BIO1 range:", round(rng[1, "min"] / 10, 1),
          "°C to", round(rng[1, "max"] / 10, 1),
          "°C (stored ×10)\n")
    }
  }
  # BIO12: annual precip in mm
  if (bio_num == "12") {
    cat("    BIO12 range:", round(rng[1, "min"], 0),
        "to", round(rng[1, "max"], 0), "mm\n")
  }
  # BIO15: coefficient of variation (0-100)
  if (bio_num == "15") {
    cat("    BIO15 range:", round(rng[1, "min"], 1),
        "to", round(rng[1, "max"], 1), "(CV)\n")
  }
  
  # Step 5: Save
  out_path <- file.path(OUT_CLIMATE_PRO, out_name)
  terra::writeRaster(r_aligned, out_path,
                     overwrite = TRUE, datatype = "FLT4S")
  cat("    ✓ Saved:", out_name, "\n\n")
  
  processed_vars[[bio_num]] <- r_aligned
}

# ── 5. STACK AND VERIFY ─────────────────────────────────────

cat("--- Verifying Climate Outputs ---\n\n")

expected_outputs <- target_vars$out_name
all_ok <- TRUE

for (f in expected_outputs) {
  fpath <- file.path(OUT_CLIMATE_PRO, f)
  if (file.exists(fpath)) {
    r   <- terra::rast(fpath)
    rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
    cat(sprintf("  ✓ %-45s  [%8.2f, %8.2f]\n",
                f,
                round(rng[1, 1], 2),
                round(rng[1, 2], 2)))
  } else {
    cat("  ✗ MISSING:", f, "\n")
    all_ok <- FALSE
  }
}

cat("\n")

# ── 6. DIAGNOSTIC FIGURE ────────────────────────────────────

n_processed <- sum(sapply(expected_outputs, function(f)
  file.exists(file.path(OUT_CLIMATE_PRO, f))))

if (n_processed > 0) {
  cat("--- Generating Climate Diagnostic Figure ---\n")
  
  png(file.path(OUT_FIG_SUPP, "S0g_climate_LGM_check.png"),
      width = 7200, height = 2400, res = 300)
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
  
  bio_labels <- c("01" = "BIO1: Mean Annual Temp (×10 °C)\nLGM ~21 ka",
                  "12" = "BIO12: Annual Precipitation (mm)\nLGM ~21 ka",
                  "15" = "BIO15: Precip Seasonality (CV)\nLGM ~21 ka")
  
  for (bio_num in c("01", "12", "15")) {
    out_name <- target_vars$out_name[target_vars$bio_num == bio_num]
    fpath    <- file.path(OUT_CLIMATE_PRO, out_name)
    if (file.exists(fpath)) {
      r <- terra::rast(fpath)
      terra::plot(r,
                  main = bio_labels[bio_num],
                  col  = viridisLite::viridis(100),
                  axes = FALSE)
      terra::plot(boundary_vect, add = TRUE,
                  border = "white", lwd = 0.8)
    } else {
      plot.new()
      title(main = paste(bio_labels[bio_num], "\n[NOT FOUND]"))
    }
  }
  
  dev.off()
  cat("  ✓ S0g_climate_LGM_check.png\n\n")
}

# ── 7. METHODS NOTE — SAVE FOR MANUSCRIPT ───────────────────

methods_note <- paste0(
  "CHELSA-TraCE21k LGM palaeoclimate variables (BIO1, BIO12, BIO15; ",
  "time slice -210, ~21,000 BP; Karger et al. 2023) were resampled ",
  "from 30 arc-second (~1 km) resolution to 30 m using bilinear ",
  "interpolation and projected to WGS84 UTM Zone 44N (EPSG:32644). ",
  "These variables are included in the Upper Palaeolithic sub-model ",
  "only, as spatially variable indicators of long-term climatic ",
  "gradients relevant to the LGM occupation period, acknowledging ",
  "that temporal alignment is approximate. No palaeoclimate variables ",
  "are included in the pooled model, the Lower Palaeolithic sub-model ",
  "(Acheulian occupations >200 ka predate CHELSA-TraCE21k coverage), ",
  "or the Middle Palaeolithic sub-model (MIS3 ~40 ka lies outside ",
  "CHELSA-TraCE21k extent, which extends only to ~22 ka)."
)

writeLines(methods_note,
           file.path(OUT_CLIMATE_PRO,
                     "climate_methods_note.txt"))
cat("  ✓ Methods note saved: climate_methods_note.txt\n\n")

# ── 8. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 06 COMPLETE — Summary\n")
cat("========================================\n")
cat("Target: CHELSA-TraCE21k LGM (~21 ka)\n")
cat("Variables: BIO1, BIO12, BIO15\n")
cat("Processed:", n_processed, "of 3 variables\n")
if (n_processed < 3) {
  cat("\n⚠ INCOMPLETE — check climate file names above.\n")
  cat("  Expected pattern: CHELSA_TraCE21k_bio01_-210_V1.0.tif\n")
  cat("  Actual files in folder listed above.\n")
}
cat("\nUSAGE REMINDER:\n")
cat("  Pooled model:    NO climate variables\n")
cat("  LP sub-model:    NO climate variables\n")
cat("  MP sub-model:    NO climate variables\n")
cat("  UP sub-model:    BIO1 + BIO12 + BIO15 (LGM)\n")
cat("\nOutputs saved to:\n  ", OUT_CLIMATE_PRO, "\n")
cat("\nNext: Run Script 07 — Predictor Stack + VIF Screening\n")
cat("  Assembles all 14 predictors into final stack\n")
cat("  Runs usdm::vifcor() — threshold VIF < 5\n")
cat("========================================\n")