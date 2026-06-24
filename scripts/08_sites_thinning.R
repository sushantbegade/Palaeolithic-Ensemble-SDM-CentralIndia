# ============================================================
# SCRIPT 08: SPATIAL THINNING OF SITE RECORDS
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 08 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Applies spatial thinning (spThin, 1km minimum distance)
#   to the site occurrence dataset. Cultural period pools
#   defined from the RELATIVE CHRONOLOGY column — standardized
#   period attribution based on tools and typology, NOT from
#   Cultural.Attribution (raw author-reported terminology,
#   inconsistent across literature).
#
#   Period pools (multi-period sites in multiple pools):
#     Pooled: all 197 sites
#     LP pool: sites with Lower Palaeolithic attribution
#     MP pool: sites with Middle Palaeolithic attribution
#     UP pool: sites with Upper Palaeolithic attribution
#
#   Sub-model decision rule post-thinning:
#     N >= 40:   all 6 algorithms
#     N = 25-39: MaxEnt + RF + GAM only
#     N < 25:    sub-model not viable
#
#   Attribution sensitivity (Supp S3):
#     Stage 1: GPS-confirmed sites only (field verified +
#              referenced GPS — Location.Precision column)
#     Stage 2: all sites (no probable tier in this dataset)
#
# KEY COLUMN:
#   Relative.Chronology....Lower...Middle.Palaeolithic..etc..
#   (sf reads column name spaces/brackets as dots)
#   This column contains standardized LP/MP/UP attribution
#   based on lithic assemblage analysis and typology.
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

library(spThin)

cat("\n========================================\n")
cat("SCRIPT 08: Spatial Thinning\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. LOAD SITES ────────────────────────────────────────────

cat("--- Loading Site Data ---\n")

sites_sf <- sf::st_read(file.path(OUT_SITES,
                                  "sites_all_utm44n.gpkg"),
                        quiet = TRUE)
cat("  Total sites loaded:", nrow(sites_sf), "\n")
cat("  Columns available:\n")
for (col in names(sites_sf)[names(sites_sf) != "geom"]) {
  cat("   ", col, "\n")
}
cat("\n")

sites_df <- as.data.frame(sf::st_drop_geometry(sites_sf))

# UTM coordinates (for saving back to sf)
coords_utm     <- sf::st_coordinates(sites_sf)
sites_df$UTM_X <- coords_utm[, 1]
sites_df$UTM_Y <- coords_utm[, 2]

# WGS84 coordinates (required by spThin)
sites_geo          <- sf::st_transform(sites_sf, crs = 4326)
coords_geo         <- sf::st_coordinates(sites_geo)
sites_df$Longitude_WGS84 <- coords_geo[, 1]
sites_df$Latitude_WGS84  <- coords_geo[, 2]

# ── 2. IDENTIFY RELATIVE CHRONOLOGY COLUMN ──────────────────

cat("--- Identifying Relative Chronology Column ---\n\n")

all_cols <- names(sites_df)

# Find column matching "Relative" and "Chronology"
relchron_col <- all_cols[grepl("Relative", all_cols,
                               ignore.case = TRUE) &
                           grepl("Chronology|Chrono",
                                 all_cols, ignore.case = TRUE)]

if (length(relchron_col) == 0) {
  # Fallback: any column with "chronology"
  relchron_col <- all_cols[grepl("chronol|period|chrono",
                                 all_cols, ignore.case = TRUE)]
}

if (length(relchron_col) == 0) {
  stop("Relative Chronology column not found.\n",
       "Columns available: ",
       paste(all_cols, collapse = ", "))
}

relchron_col <- relchron_col[1]
cat("  Column selected:", relchron_col, "\n\n")

chron_str <- as.character(sites_df[[relchron_col]])

cat("  All unique values in Relative Chronology:\n")
for (v in sort(unique(chron_str))) cat("    -", v, "\n")
cat("\n")

# ── 3. CULTURAL PERIOD CLASSIFICATION ───────────────────────
# Using Relative.Chronology — standardized by Begade from
# assemblage analysis. Values should be cleaner than
# Cultural.Attribution (which reflects raw author terminology).
#
# LP: contains "Lower" (Lower Palaeolithic)
# MP: contains "Middle" (Middle Palaeolithic)
# UP: contains "Upper" OR "Late Stone Age" OR "Late Palaeolithic"
#     (Upper Palaeolithic and equivalent)
#
# Multi-period attributions (e.g., "Lower and Middle
# Palaeolithic") are counted in ALL applicable period pools.

cat("--- Cultural Period Classification ---\n")
cat("  (from Relative Chronology column — typology-based)\n\n")

is_LP <- grepl("Lower", chron_str, ignore.case = TRUE)
is_MP <- grepl("Middle", chron_str, ignore.case = TRUE)
is_UP <- grepl("Upper|Late Stone Age|Late Palaeolithic",
               chron_str, ignore.case = TRUE)

cat(sprintf("  LP pool (before thinning): %d sites\n", sum(is_LP)))
cat(sprintf("  MP pool (before thinning): %d sites\n", sum(is_MP)))
cat(sprintf("  UP pool (before thinning): %d sites\n", sum(is_UP)))
cat(sprintf("  Pooled (all):              %d sites\n\n",
            nrow(sites_df)))

# Report what each pool contains
cat("  LP attribution values in pool:\n")
for (v in sort(unique(chron_str[is_LP]))) cat("    +", v, "\n")
cat("\n  MP attribution values in pool:\n")
for (v in sort(unique(chron_str[is_MP]))) cat("    +", v, "\n")
cat("\n  UP attribution values in pool:\n")
for (v in sort(unique(chron_str[is_UP]))) cat("    +", v, "\n")

# Report unclassified
unclassified <- !is_LP & !is_MP & !is_UP
if (any(unclassified)) {
  cat("\n  Unclassified (pooled model only):\n")
  for (v in sort(unique(chron_str[unclassified]))) {
    cat("    ?", v, "\n")
  }
}
cat("\n")

# ── 4. CONFIDENCE CLASSIFICATION ────────────────────────────
# Location.Precision column — confirmed present in your data:
#   "Exact ±5-10 m (Field Verified)"               → confirmed
#   "Exact ±5-10 m (Reported By Referenced Author)" → confirmed
#   "GIS Map Based"                                 → uncertain

cat("--- Confidence Classification ---\n")
cat("  (from Location.Precision column)\n\n")

prec_str <- as.character(sites_df$Location.Precision)

is_confirmed <- grepl("Field Verified|Reported By Referenced",
                      prec_str, ignore.case = TRUE)
is_uncertain <- grepl("GIS Map Based", prec_str,
                      ignore.case = TRUE)

cat(sprintf("  GPS confirmed (±5-10m): %d sites\n",
            sum(is_confirmed)))
cat(sprintf("  GIS-based (uncertain):  %d sites\n",
            sum(is_uncertain)))
cat(sprintf("  Other:                  %d sites\n\n",
            sum(!is_confirmed & !is_uncertain)))

# ── 5. THINNING HELPER ───────────────────────────────────────

thin_sites <- function(df, label, thin_dist_km = THIN_DIST_KM) {
  
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  [%-24s]   0 →   0  (skipped)\n", label))
    return(NULL)
  }
  
  thin_df <- data.frame(
    species   = "Palaeolithic",
    Longitude = df$Longitude_WGS84,
    Latitude  = df$Latitude_WGS84
  )
  
  set.seed(42)
  result <- spThin::thin(
    loc.data   = thin_df,
    lat.col    = "Latitude",
    long.col   = "Longitude",
    spec.col   = "species",
    thin.par   = thin_dist_km,
    reps       = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    verbose    = FALSE
  )
  
  n_kept <- sapply(result, nrow)
  best   <- result[[which.max(n_kept)]]
  
  thin_key  <- paste(round(best$Longitude, 5),
                     round(best$Latitude,  5))
  orig_key  <- paste(round(df$Longitude_WGS84, 5),
                     round(df$Latitude_WGS84,  5))
  kept_idx  <- which(orig_key %in% thin_key)
  thinned   <- df[kept_idx, ]
  
  cat(sprintf("  [%-24s] %3d → %3d  (%.0f%% retained)\n",
              label, nrow(df), nrow(thinned),
              100 * nrow(thinned) / nrow(df)))
  return(thinned)
}

n_thin <- function(d) if (is.null(d)) 0L else nrow(d)

df_to_sf <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sf::st_as_sf(df,
               coords = c("UTM_X", "UTM_Y"),
               crs    = 32644,
               remove = FALSE)
}

save_thinned <- function(df, filename) {
  sf_obj <- df_to_sf(df)
  if (is.null(sf_obj)) return(invisible(NULL))
  out_path <- file.path(OUT_SITES, filename)
  sf::st_write(sf_obj, out_path, delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("    → Saved: %-45s (%d sites)\n",
              filename, nrow(sf_obj)))
}

# ── 6. APPLY THINNING ────────────────────────────────────────

cat("--- Spatial Thinning (1km minimum distance) ---\n\n")

thinned_pooled <- thin_sites(sites_df,             "Pooled (all)")
thinned_LP     <- thin_sites(sites_df[is_LP, ],    "LP")
thinned_MP     <- thin_sites(sites_df[is_MP, ],    "MP")
thinned_UP     <- thin_sites(sites_df[is_UP, ],    "UP")
cat("\n")

# ── 7. SENSITIVITY ANALYSIS (Supp S3) ───────────────────────
# Stage 1: GPS-confirmed sites only (excludes GIS-estimated)
# Stage 2: all sites (no probable tier — Stage 2 = main pools)

cat("--- Attribution Sensitivity Analysis (Supp S3) ---\n\n")

thinned_LP_conf <- thin_sites(
  sites_df[is_LP & is_confirmed, ], "LP confirmed (GPS only)")
thinned_MP_conf <- thin_sites(
  sites_df[is_MP & is_confirmed, ], "MP confirmed (GPS only)")
thinned_UP_conf <- thin_sites(
  sites_df[is_UP & is_confirmed, ], "UP confirmed (GPS only)")

cat("\n  Stage 2 = all sites = same as main pools\n",
    " (no probable tier in this dataset)\n\n")

# ── 8. SAVE ALL DATASETS ─────────────────────────────────────

cat("--- Saving Thinned Datasets ---\n\n")

save_thinned(thinned_pooled,  "sites_thinned_pooled.gpkg")
save_thinned(thinned_LP,      "sites_thinned_LP.gpkg")
save_thinned(thinned_MP,      "sites_thinned_MP.gpkg")
save_thinned(thinned_UP,      "sites_thinned_UP.gpkg")
save_thinned(thinned_LP_conf, "sites_thinned_LP_confirmed.gpkg")
save_thinned(thinned_MP_conf, "sites_thinned_MP_confirmed.gpkg")
save_thinned(thinned_UP_conf, "sites_thinned_UP_confirmed.gpkg")

# ── 9. SUB-MODEL DECISION RULE ───────────────────────────────

cat("\n--- Sub-Model Decision Rule ---\n\n")
cat("  N >= 40:   all 6 algorithms\n")
cat("  N = 25-39: MaxEnt + RF + GAM only\n")
cat("  N < 25:    sub-model not viable\n\n")

apply_rule <- function(thinned_df, period) {
  n    <- n_thin(thinned_df)
  rule <- if (n >= 40) "All 6 algorithms" else
    if (n >= 25) "MaxEnt + RF + GAM only" else
      "NOT VIABLE (N < 25)"
  cat(sprintf("  %-8s N = %3d → %s\n", period, n, rule))
  return(rule)
}

r_pooled <- apply_rule(thinned_pooled, "Pooled")
r_LP     <- apply_rule(thinned_LP,     "LP")
r_MP     <- apply_rule(thinned_MP,     "MP")
r_UP     <- apply_rule(thinned_UP,     "UP")

# ── 10. TABLE 1 ──────────────────────────────────────────────

cat("\n--- Table 1: Site Count Summary ---\n\n")

table1 <- data.frame(
  Period     = c("Pooled","LP","MP","UP",
                 "LP_GPS_confirmed","MP_GPS_confirmed",
                 "UP_GPS_confirmed"),
  N_raw      = c(nrow(sites_df),
                 sum(is_LP), sum(is_MP), sum(is_UP),
                 sum(is_LP & is_confirmed),
                 sum(is_MP & is_confirmed),
                 sum(is_UP & is_confirmed)),
  N_thinned  = c(n_thin(thinned_pooled),
                 n_thin(thinned_LP),
                 n_thin(thinned_MP),
                 n_thin(thinned_UP),
                 n_thin(thinned_LP_conf),
                 n_thin(thinned_MP_conf),
                 n_thin(thinned_UP_conf)),
  Decision   = c(r_pooled, r_LP, r_MP, r_UP,
                 "Sensitivity Stage 1",
                 "Sensitivity Stage 1",
                 "Sensitivity Stage 1"),
  stringsAsFactors = FALSE
)

write.csv(table1,
          file.path(OUT_TABLES, "Table1_site_counts.csv"),
          row.names = FALSE)

cat(sprintf("  %-22s %8s %10s  %s\n",
            "Period","N raw","N thinned","Decision"))
cat("  ", paste(rep("-", 66), collapse=""), "\n")
for (i in seq_len(nrow(table1))) {
  cat(sprintf("  %-22s %8d %10d  %s\n",
              table1$Period[i],
              table1$N_raw[i],
              table1$N_thinned[i],
              table1$Decision[i]))
}
cat("\n  ✓ Table 1 saved: Table1_site_counts.csv\n\n")

# ── 11. DIAGNOSTIC MAP ───────────────────────────────────────

cat("--- Diagnostic Map ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

png(file.path(OUT_FIG_SUPP, "S0h_thinning_result.png"),
    width = 4800, height = 2400, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

terra::plot(template_30m, col = "grey95",
            legend = FALSE, axes = FALSE,
            main = sprintf("Before Thinning (N = %d)",
                           nrow(sites_df)))
terra::plot(boundary_vect, add = TRUE,
            border = "grey50", lwd = 0.8)
terra::plot(terra::vect(sites_geo), add = TRUE,
            col = "steelblue", pch = 16, cex = 0.3)

if (!is.null(thinned_pooled)) {
  tf <- df_to_sf(thinned_pooled)
  tg <- sf::st_transform(tf, 4326)
  terra::plot(template_30m, col = "grey95",
              legend = FALSE, axes = FALSE,
              main = sprintf("After Thinning 1km (N = %d)",
                             n_thin(thinned_pooled)))
  terra::plot(boundary_vect, add = TRUE,
              border = "grey50", lwd = 0.8)
  terra::plot(terra::vect(tg), add = TRUE,
              col = "firebrick", pch = 16, cex = 0.4)
}
dev.off()
cat("  ✓ S0h_thinning_result.png\n\n")

# ── 12. SUMMARY ─────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 08 COMPLETE — Summary\n")
cat("========================================\n")
cat("Column used: Relative Chronology\n")
cat("  (typology-based — NOT Cultural Attribution)\n\n")
cat(sprintf("Pooled: %3d → %3d  %s\n",
            nrow(sites_df), n_thin(thinned_pooled), r_pooled))
cat(sprintf("LP:     %3d → %3d  %s\n",
            sum(is_LP), n_thin(thinned_LP), r_LP))
cat(sprintf("MP:     %3d → %3d  %s\n",
            sum(is_MP), n_thin(thinned_MP), r_MP))
cat(sprintf("UP:     %3d → %3d  %s\n",
            sum(is_UP), n_thin(thinned_UP), r_UP))
cat("\nExpected pools (from research design):\n")
cat("  LP: ~75  MP: ~114  UP: ~79\n")
cat("If counts differ significantly, check the unique\n")
cat("values printed above for the Relative Chronology\n")
cat("column and adjust regex if needed.\n")
cat("\nNext: Run Script 09 — Bias Correction + Background\n")
cat("========================================\n")