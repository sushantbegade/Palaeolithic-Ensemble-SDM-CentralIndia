# ============================================================
# SCRIPT 02: LOAD ALL RAW DATA AND REPROJECT TO UTM 44N
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 02 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1. Sources Script 01 (loads all paths and packages)
#   2. Loads the district boundary shapefile
#   3. Reprojects boundary to UTM Zone 44N (EPSG:32644)
#   4. Loads and reprojects the DEM (Cartosat-1 CartoDEM)
#   5. Creates the MASTER TEMPLATE raster — 30m, UTM 44N
#      All 14 predictor rasters MUST match this template
#   6. Loads and reprojects all vector layers
#   7. Aligns Sentinel-2A rasters to master template
#   8. Loads site data and reprojects to UTM 44N
#   9. Saves all processed files to data_processed/
# ============================================================
# FIX APPLIED (v2.0):
#   terra::project() masked by other packages (likely spgwr,
#   stats, or raster). All terra spatial functions now use
#   explicit terra:: prefix to avoid namespace conflicts:
#     terra::project(), terra::resample(), terra::mask(),
#     terra::crop(), terra::writeRaster(), terra::rast(),
#     terra::vect(), terra::ext(), terra::res(), terra::crs(),
#     terra::minmax(), terra::ncell(), terra::global()
# ============================================================
# HOW TO USE:
#   Run after Script 01 completes successfully.
#   Or source directly — Script 01 is sourced automatically.
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 02: Load and Reproject\n")
cat("========================================\n\n")

# ── 1. DETECT SHAPEFILES IN EACH FOLDER ─────────────────────
# Automatically finds the .shp file in each folder.
# No need to hardcode filenames.

find_shp <- function(folder_path) {
  shp_files <- list.files(folder_path,
                          pattern     = "\\.shp$",
                          full.names  = TRUE,
                          recursive   = TRUE)
  if (length(shp_files) == 0) {
    stop(paste("No shapefile found in:", folder_path))
  }
  if (length(shp_files) > 1) {
    cat("  Multiple shapefiles found — using first:",
        basename(shp_files[1]), "\n")
  }
  return(shp_files[1])
}

# ── 2. LOAD DISTRICT BOUNDARY ────────────────────────────────

cat("--- Loading District Boundary ---\n")

shp_district <- find_shp(PATH_DISTRICT)
cat("  File found:", basename(shp_district), "\n")

district_raw <- sf::st_read(shp_district, quiet = TRUE)
cat("  CRS detected:", sf::st_crs(district_raw)$Name, "\n")
cat("  Features:", nrow(district_raw), "\n")
cat("  Fields:",
    paste(names(district_raw)[names(district_raw) != "geometry"],
          collapse = ", "), "\n")

# Reproject to UTM Zone 44N
district_utm <- sf::st_transform(district_raw, crs = CRS_PROJECT)
cat("  ✓ Reprojected to UTM Zone 44N\n")

# Merge all districts into one boundary polygon
# (Nagpur and Chandrapur are separate features — union them)
district_merged <- sf::st_union(district_utm)
district_merged <- sf::st_sf(geometry = district_merged)
cat("  ✓ Districts merged into single study area boundary\n")

# Convert to terra SpatVector for raster masking operations
boundary_vect <- terra::vect(district_merged)

# Save processed district boundary
out_path_boundary <- file.path(OUT_SITES,
                               "study_area_boundary_utm44n.gpkg")
sf::st_write(district_merged,
             out_path_boundary,
             delete_dsn = TRUE,
             quiet      = TRUE)
cat("  ✓ Saved:", basename(out_path_boundary), "\n\n")

# ── 3. LOAD AND PROCESS DEM ──────────────────────────────────

cat("--- Loading DEM (Cartosat-1 CartoDEM) ---\n")

dem_path <- file.path(PATH_DEM, "DEM.tif")
dem_raw  <- terra::rast(dem_path)

cat("  Original CRS:       ",
    terra::crs(dem_raw, describe = TRUE)$name, "\n")
cat("  Original resolution:",
    paste(round(terra::res(dem_raw), 4), collapse = " x "), "\n")
cat("  Original extent:    ",
    paste(round(as.vector(terra::ext(dem_raw)), 3),
          collapse = ", "), "\n")
cat("  Value range:         min =",
    round(terra::minmax(dem_raw)[1], 1),
    "m,  max =",
    round(terra::minmax(dem_raw)[2], 1), "m\n")

# ── KEY FIX: explicit terra::project() to avoid namespace clash ──
cat("  Reprojecting to UTM 44N at 30m...\n")

dem_utm <- terra::project(
  x      = dem_raw,
  y      = CRS_PROJECT,   # target CRS string "EPSG:32644"
  method = "bilinear",    # bilinear for continuous elevation
  res    = RESOLUTION     # 30 metres
)

cat("  ✓ Reprojected\n")
cat("  New resolution:",
    paste(round(terra::res(dem_utm), 2), collapse = " x "), "m\n")

# Clip to study area boundary
dem_clipped <- terra::crop(dem_utm, boundary_vect, mask = TRUE)
cat("  ✓ Clipped and masked to study area\n")

# Verify elevation values are sensible for Vidarbha
# Expected: ~150m to ~900m (Satpura foothills to Deccan plateau)
dem_stats <- terra::minmax(dem_clipped)
cat("  Elevation range:",
    round(dem_stats[1], 1), "m to",
    round(dem_stats[2], 1), "m\n")

if (dem_stats[1] < 0 || dem_stats[2] > 1500) {
  warning("Elevation values outside expected range for Vidarbha.",
          " Check DEM for errors.")
} else {
  cat("  ✓ Elevation values within expected range for Vidarbha\n")
}

# Save processed DEM
out_path_dem <- file.path(OUT_PREDICTORS, "DEM_30m_utm44n.tif")
terra::writeRaster(dem_clipped,
                   out_path_dem,
                   overwrite = TRUE,
                   datatype  = "FLT4S")  # Float32 for elevation
cat("  ✓ Saved:", basename(out_path_dem), "\n\n")

# ── 4. CREATE MASTER TEMPLATE RASTER ────────────────────────
# MOST IMPORTANT output of Script 02.
# Every predictor raster MUST exactly match this template:
# same extent, resolution, CRS, and origin.
# All subsequent scripts use template_30m as the alignment
# target in terra::resample().

cat("--- Creating Master Template Raster (30m) ---\n")

template_30m <- dem_clipped  # DEM is our reference grid

cat("  CRS:        ",
    terra::crs(template_30m, describe = TRUE)$name, "\n")
cat("  Resolution: ",
    paste(terra::res(template_30m), collapse = " x "), "m\n")
cat("  Extent:     ",
    paste(round(as.vector(terra::ext(template_30m)), 0),
          collapse = ", "),
    "(xmin, xmax, ymin, ymax)\n")
cat("  Dimensions: ", nrow(template_30m), "rows x",
    ncol(template_30m), "cols\n")
cat("  Total cells:", terra::ncell(template_30m), "\n")

# Count non-NA cells
n_valid <- terra::global(!is.na(template_30m), "sum",
                         na.rm = TRUE)[1, 1]
cat("  Non-NA cells:", n_valid, "\n")

out_path_template <- file.path(OUT_PREDICTORS,
                               "TEMPLATE_30m_utm44n.tif")
terra::writeRaster(template_30m,
                   out_path_template,
                   overwrite = TRUE,
                   datatype  = "FLT4S")
cat("  ✓ Master template saved:", basename(out_path_template), "\n")
cat("  → All other scripts will align to this template\n\n")

# ── 5. LOAD AND REPROJECT ALL VECTOR LAYERS ─────────────────

cat("--- Loading and Reprojecting Vector Layers ---\n")

# Helper: load shapefile, reproject to UTM 44N, save as GeoPackage
process_vector <- function(folder_path, layer_name, out_filename) {
  cat(paste0("  [", layer_name, "]\n"))
  
  shp_file <- find_shp(folder_path)
  cat("    File:", basename(shp_file), "\n")
  
  v_raw <- sf::st_read(shp_file, quiet = TRUE)
  cat("    CRS:", sf::st_crs(v_raw)$Name, "\n")
  cat("    Features:", nrow(v_raw), "\n")
  cat("    Geometry:",
      as.character(sf::st_geometry_type(v_raw)[1]), "\n")
  
  v_utm <- sf::st_transform(v_raw, crs = CRS_PROJECT)
  
  # GeoPackage — no field name truncation unlike shapefile
  out_path <- file.path(OUT_PREDICTORS, out_filename)
  sf::st_write(v_utm, out_path, delete_dsn = TRUE, quiet = TRUE)
  cat("    ✓ Saved:", out_filename, "\n\n")
  
  return(v_utm)
}

geology_utm   <- process_vector(PATH_GEOLOGY,
                                "Geology",
                                "geology_utm44n.gpkg")

geomorph_utm  <- process_vector(PATH_GEOMORPH,
                                "Geomorphology",
                                "geomorphology_utm44n.gpkg")

lithology_utm <- process_vector(PATH_LITHOLOGY,
                                "Lithology",
                                "lithology_utm44n.gpkg")

rivers_utm    <- process_vector(PATH_RIVERS,
                                "Rivers",
                                "rivers_utm44n.gpkg")

waterbody_utm <- process_vector(PATH_WATERBODY,
                                "Waterbody",
                                "waterbody_utm44n.gpkg")

# ── 6. LOAD AND ALIGN SENTINEL-2A RASTERS ───────────────────

cat("--- Loading and Aligning Sentinel-2A Outputs ---\n")

# Helper: reproject if needed, resample to master template, mask, save
# ── KEY FIX: all terra spatial calls use explicit terra:: prefix ──

align_raster <- function(file_path, layer_name,
                         out_filename, method = "bilinear") {
  cat(paste0("  [", layer_name, "]\n"))
  
  r <- terra::rast(file_path)
  cat("    Original CRS:       ",
      terra::crs(r, describe = TRUE)$name, "\n")
  cat("    Original resolution:",
      paste(round(terra::res(r), 1), collapse = " x "), "\n")
  
  # Reproject only if CRS does not already match template
  if (!identical(terra::crs(r), terra::crs(template_30m))) {
    r <- terra::project(x = r, y = template_30m, method = method)
    cat("    ✓ Reprojected to UTM 44N\n")
  } else {
    cat("    CRS already matches template — skip reproject\n")
  }
  
  # Resample to align exactly: same extent, origin, resolution
  r_aligned <- terra::resample(r, template_30m, method = method)
  r_aligned  <- terra::mask(r_aligned, boundary_vect)
  
  out_path <- file.path(OUT_PREDICTORS, out_filename)
  terra::writeRaster(r_aligned, out_path,
                     overwrite = TRUE, datatype = "FLT4S")
  cat("    ✓ Saved:", out_filename, "\n\n")
  
  return(r_aligned)
}

ndvi_aligned  <- align_raster(
  file.path(PATH_S2, "NDVI.tif"),
  "NDVI (dry season)",
  "NDVI_30m_utm44n.tif"
)

mndwi_aligned <- align_raster(
  file.path(PATH_S2, "MNDWI.tif"),
  "MNDWI (palaeochannel detection)",
  "MNDWI_30m_utm44n.tif"
)

# Also align NDBI — supplementary heritage threat analysis
ndbi_aligned  <- align_raster(
  file.path(PATH_S2, "NDBI.tif"),
  "NDBI (built-up/quarry detection)",
  "NDBI_30m_utm44n.tif"
)

# ── 7. LOAD SITES DATA ───────────────────────────────────────

cat("--- Loading Palaeolithic Sites ---\n")

sites_raw <- readxl::read_xlsx(PATH_SITES_XLSX)
cat("  Total records:", nrow(sites_raw), "\n")
cat("  Columns:", paste(names(sites_raw), collapse = ", "), "\n")

# Convert to spatial (WGS84 geographic = EPSG:4326)
# Adjust column names if your Excel uses different headers
# Expected: "Longitude" and "Latitude" — check and update if needed
lon_col <- "Longitude"
lat_col <- "Latitude"

if (!all(c(lon_col, lat_col) %in% names(sites_raw))) {
  stop(paste(
    "Longitude/Latitude columns not found in Excel.",
    "Columns present:", paste(names(sites_raw), collapse = ", ")
  ))
}

sites_geo <- sf::st_as_sf(sites_raw,
                          coords = c(lon_col, lat_col),
                          crs    = 4326,
                          remove = FALSE)  # keep lat/lon columns

# Reproject to UTM 44N
sites_utm <- sf::st_transform(sites_geo, crs = CRS_PROJECT)

# Verify all sites fall within study area boundary
sites_in_boundary <- sf::st_within(sites_utm,
                                   district_merged,
                                   sparse = FALSE)
n_outside <- sum(!sites_in_boundary)

if (n_outside > 0) {
  cat("  ⚠ WARNING:", n_outside,
      "sites fall outside the district boundary!\n")
  cat("  Check these rows:\n")
  # Use row index if no Site ID column; adjust column name if needed
  if ("Site ID" %in% names(sites_raw)) {
    outside_ids <- sites_raw[["Site ID"]][!sites_in_boundary]
  } else {
    outside_ids <- which(!sites_in_boundary)
  }
  cat("   ", paste(outside_ids, collapse = ", "), "\n")
} else {
  cat("  ✓ All", nrow(sites_utm),
      "sites confirmed within study area boundary\n")
}

# Save reprojected sites
out_sites_path <- file.path(OUT_SITES, "sites_all_utm44n.gpkg")
sf::st_write(sites_utm, out_sites_path,
             delete_dsn = TRUE, quiet = TRUE)
cat("  ✓ Sites saved:", basename(out_sites_path), "\n\n")

# ── 8. QUICK VISUAL VERIFICATION ────────────────────────────

cat("--- Quick Map Check ---\n")

png(file.path(OUT_FIG_SUPP, "S0_data_alignment_check.png"),
    width = 2400, height = 2400, res = 300)

terra::plot(dem_clipped,
            col  = gray.colors(100),
            main = "Data Alignment Check\n(DEM + Boundary + Sites)",
            axes = TRUE)
terra::plot(boundary_vect, add = TRUE,
            border = "red", lwd = 2)
terra::plot(terra::vect(rivers_utm), add = TRUE,
            col = "blue", lwd = 0.5)
terra::plot(terra::vect(sites_utm),  add = TRUE,
            col = "yellow", pch = 16, cex = 0.3)

dev.off()
cat("  ✓ Alignment check map saved to outputs/figures/supplementary/\n\n")

# ── 9. SUMMARY REPORT ───────────────────────────────────────

cat("\n========================================\n")
cat("SCRIPT 02 COMPLETE — Summary Report\n")
cat("========================================\n")
cat("Study area boundary:  ✓ Processed\n")
cat("DEM (30m UTM 44N):    ✓ Processed\n")
cat("Master template:      ✓ Created\n")
cat("Geology:              ✓ Reprojected\n")
cat("Geomorphology:        ✓ Reprojected\n")
cat("Lithology:            ✓ Reprojected\n")
cat("Rivers:               ✓ Reprojected\n")
cat("Waterbody:            ✓ Reprojected\n")
cat("NDVI (30m):           ✓ Aligned to template\n")
cat("MNDWI (30m):          ✓ Aligned to template\n")
cat("NDBI (30m):           ✓ Aligned to template\n")
cat("Sites (", nrow(sites_raw), "records):  ✓ Loaded and reprojected\n",
    sep = "")
cat("\nAll processed files saved to:\n")
cat(" ", OUT_PREDICTORS, "\n")
cat("\nNext: Run Script 03 — Terrain Derivatives from DEM\n")
cat("========================================\n")