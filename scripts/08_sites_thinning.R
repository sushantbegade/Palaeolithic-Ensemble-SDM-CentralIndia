# ============================================================
# SCRIPT 08: SPATIAL THINNING OF SITE RECORDS
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 08 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Applies spatial thinning (spThin, 1km minimum distance)
#   to reduce spatial autocorrelation bias in the site
#   occurrence dataset. Thinning applied:
#     (a) Pooled dataset — all 197 sites
#     (b) LP pool (sites with LP attribution)
#     (c) MP pool (sites with MP attribution)
#     (d) UP pool (sites with UP attribution)
#   Each thinned independently; cultural period pools include
#   all sites attributed to that period (including multi-period
#   sites, e.g., LP+MP counted in both LP and MP pools).
#
#   Sub-model decision rule applied post-thinning:
#     N >= 40 thinned: all 6 algorithms
#     N = 25-39:       MaxEnt + RF + GAM only
#     N < 25:          sub-model not viable
#
#   Attribution uncertainty sensitivity analysis (Supp S3):
#     Stage 1: CONFIRMED sites only per period
#     Stage 2: CONFIRMED + PROBABLE sites per period
#   Both thinned and saved for use in Script 20.
#
# INPUTS:
#   sites_all_utm44n.gpkg — 197 sites (from Script 02)
#   Excel: columns include Cultural Attribution,
#          Relative Chronology, Location Precision
#
# OUTPUTS (to data_processed/sites/):
#   sites_thinned_pooled.gpkg
#   sites_thinned_LP.gpkg
#   sites_thinned_MP.gpkg
#   sites_thinned_UP.gpkg
#   sites_thinned_LP_confirmed.gpkg    (sensitivity Stage 1)
#   sites_thinned_MP_confirmed.gpkg
#   sites_thinned_LP_conf_prob.gpkg    (sensitivity Stage 2)
#   sites_thinned_MP_conf_prob.gpkg
#   Table1_site_counts.csv
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
cat("  Fields:", paste(names(sites_sf)[
  !names(sites_sf) %in% c("geom","geometry")],
  collapse=", "), "\n\n")

# Extract WGS84 coordinates for spThin (requires lat/lon)
sites_geo <- sf::st_transform(sites_sf, crs = 4326)
coords    <- sf::st_coordinates(sites_geo)
sites_df  <- as.data.frame(sf::st_drop_geometry(sites_sf))
sites_df$Longitude_WGS84 <- coords[, 1]
sites_df$Latitude_WGS84  <- coords[, 2]
sites_df$UTM_X <- sf::st_coordinates(sites_sf)[, 1]
sites_df$UTM_Y <- sf::st_coordinates(sites_sf)[, 2]

# ── 2. IDENTIFY KEY COLUMNS ─────────────────────────────────
# Auto-detect cultural attribution and confidence columns

cat("--- Detecting Attribute Columns ---\n")

all_cols <- names(sites_df)

# Cultural attribution column
cult_col_candidates <- c("Cultural Attribution",
                         "cultural_attribution",
                         "Cultural.Attribution",
                         "Relative Chronology",
                         "relative_chronology",
                         "Culture", "Period", "culture",
                         "period", "attribution")
cult_col <- NULL
for (cand in cult_col_candidates) {
  if (cand %in% all_cols) { cult_col <- cand; break }
}
if (is.null(cult_col)) {
  # Try partial match
  cult_col <- all_cols[grepl("cultur|chrono|period|attrib",
                             all_cols, ignore.case = TRUE)][1]
}
cat("  Cultural attribution column:", cult_col, "\n")
cat("  Unique values:\n")
cult_vals <- sort(unique(as.character(sites_df[[cult_col]])))
for (v in cult_vals) cat("    -", v, "\n")

# Location precision / confidence column
prec_col_candidates <- c("Location Precision",
                         "location_precision",
                         "Location.Precision",
                         "Attribution Confidence",
                         "Confidence", "confidence",
                         "precision", "accuracy")
prec_col <- NULL
for (cand in prec_col_candidates) {
  if (cand %in% all_cols) { prec_col <- cand; break }
}
if (is.null(prec_col)) {
  prec_col <- all_cols[grepl("precis|confid|accur|quality",
                             all_cols, ignore.case = TRUE)][1]
}
cat("\n  Precision/confidence column:", prec_col, "\n")
if (!is.null(prec_col) && !is.na(prec_col)) {
  prec_vals <- sort(unique(as.character(sites_df[[prec_col]])))
  for (v in prec_vals) cat("    -", v, "\n")
}
cat("\n")

# ── 3. THINNING HELPER ───────────────────────────────────────

thin_sites <- function(df, label, thin_dist_km = THIN_DIST_KM) {
  
  if (nrow(df) == 0) {
    cat(sprintf("  [%s] 0 sites — skipping\n", label))
    return(NULL)
  }
  
  # spThin requires: data frame with Longitude, Latitude, species
  thin_df <- data.frame(
    species   = "Palaeolithic",
    Longitude = df$Longitude_WGS84,
    Latitude  = df$Latitude_WGS84
  )
  
  set.seed(42)
  
  # Run thinning — 100 repetitions, keep best (most sites)
  result <- spThin::thin(
    loc.data   = thin_df,
    lat.col    = "Latitude",
    long.col   = "Longitude",
    spec.col   = "species",
    thin.par   = thin_dist_km,   # minimum distance in km
    reps       = 100,
    locs.thinned.list.return = TRUE,
    write.files = FALSE,
    verbose    = FALSE
  )
  
  # Select repetition with most sites retained
  n_kept <- sapply(result, nrow)
  best   <- result[[which.max(n_kept)]]
  
  # Match thinned coordinates back to original df
  # (round to 6 decimal places to handle floating point)
  thin_key <- paste(round(best$Longitude, 5),
                    round(best$Latitude,  5))
  orig_key <- paste(round(df$Longitude_WGS84, 5),
                    round(df$Latitude_WGS84,  5))
  kept_idx <- which(orig_key %in% thin_key)
  
  thinned_df <- df[kept_idx, ]
  cat(sprintf("  [%-18s] %3d → %3d sites  (%.0f%% retained)\n",
              label, nrow(df), nrow(thinned_df),
              100 * nrow(thinned_df) / nrow(df)))
  
  return(thinned_df)
}

# Convert thinned data frame back to sf
df_to_sf <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sf::st_as_sf(df,
               coords = c("UTM_X", "UTM_Y"),
               crs    = 32644,
               remove = FALSE)
}

save_thinned <- function(sf_obj, filename) {
  if (is.null(sf_obj) || nrow(sf_obj) == 0) return(invisible(NULL))
  out_path <- file.path(OUT_SITES, filename)
  sf::st_write(sf_obj, out_path, delete_dsn = TRUE, quiet = TRUE)
  cat(sprintf("    → Saved: %s\n", filename))
}

# ── 4. DEFINE CULTURAL PERIOD POOLS ─────────────────────────
# Sites attributed to multiple periods counted in each period.
# LP pool: any site with LP in attribution string
# MP pool: any site with MP in attribution string
# UP pool: any site with UP in attribution string

cat("--- Defining Cultural Period Pools ---\n\n")

cult_str <- as.character(sites_df[[cult_col]])

# Detect attribution codes — handle various naming conventions
is_LP <- grepl("Lower|LP|Acheulian|Acheulean|LPA",
               cult_str, ignore.case = TRUE)
is_MP <- grepl("Middle|MP|MSA|MPA",
               cult_str, ignore.case = TRUE)
is_UP <- grepl("Upper|UP|Later Stone|LSA|UPA|Microlithic",
               cult_str, ignore.case = TRUE)

cat(sprintf("  LP pool (before thinning): %d sites\n", sum(is_LP)))
cat(sprintf("  MP pool (before thinning): %d sites\n", sum(is_MP)))
cat(sprintf("  UP pool (before thinning): %d sites\n", sum(is_UP)))
cat(sprintf("  Pooled (all):              %d sites\n\n",
            nrow(sites_df)))

# ── 5. APPLY THINNING ────────────────────────────────────────

cat("--- Applying Spatial Thinning (1km minimum) ---\n\n")

thinned_pooled <- thin_sites(sites_df,       "Pooled (all)")
thinned_LP     <- thin_sites(sites_df[is_LP,], "LP")
thinned_MP     <- thin_sites(sites_df[is_MP,], "MP")
thinned_UP     <- thin_sites(sites_df[is_UP,], "UP")

cat("\n")

# ── 6. ATTRIBUTION SENSITIVITY ANALYSIS (Supp S3) ───────────
# Stage 1: confirmed sites only per period
# Stage 2: confirmed + probable sites per period
# Requires a confidence/precision column

cat("--- Attribution Sensitivity Analysis (Supp S3) ---\n\n")

if (!is.null(prec_col) && !is.na(prec_col)) {
  
  prec_str <- as.character(sites_df[[prec_col]])
  
  # Detect confirmed and probable flags
  is_confirmed <- grepl("confirm|verified|certain|exact|GPS",
                        prec_str, ignore.case = TRUE)
  is_probable  <- grepl("probable|likely|moderate",
                        prec_str, ignore.case = TRUE)
  
  cat(sprintf("  Confirmed sites: %d\n", sum(is_confirmed)))
  cat(sprintf("  Probable sites:  %d\n", sum(is_probable)))
  cat(sprintf("  Other/uncertain: %d\n\n",
              sum(!is_confirmed & !is_probable)))
  
  # Stage 1: confirmed only
  thinned_LP_conf <- thin_sites(
    sites_df[is_LP & is_confirmed, ], "LP confirmed")
  thinned_MP_conf <- thin_sites(
    sites_df[is_MP & is_confirmed, ], "MP confirmed")
  
  # Stage 2: confirmed + probable
  thinned_LP_cp <- thin_sites(
    sites_df[is_LP & (is_confirmed | is_probable), ],
    "LP conf+prob")
  thinned_MP_cp <- thin_sites(
    sites_df[is_MP & (is_confirmed | is_probable), ],
    "MP conf+prob")
  
  cat("\n")
  
} else {
  cat("  ⚠ No precision/confidence column detected.\n")
  cat("  Sensitivity analysis not run.\n")
  cat("  Add attribution confidence to Excel and re-run.\n\n")
  thinned_LP_conf <- thinned_LP_cp <- NULL
  thinned_MP_conf <- thinned_MP_cp <- NULL
}

# ── 7. SAVE ALL THINNED DATASETS ────────────────────────────

cat("--- Saving Thinned Datasets ---\n\n")

save_list <- list(
  list(thinned_pooled,  "sites_thinned_pooled.gpkg"),
  list(thinned_LP,      "sites_thinned_LP.gpkg"),
  list(thinned_MP,      "sites_thinned_MP.gpkg"),
  list(thinned_UP,      "sites_thinned_UP.gpkg"),
  list(thinned_LP_conf, "sites_thinned_LP_confirmed.gpkg"),
  list(thinned_MP_conf, "sites_thinned_MP_confirmed.gpkg"),
  list(thinned_LP_cp,   "sites_thinned_LP_conf_prob.gpkg"),
  list(thinned_MP_cp,   "sites_thinned_MP_conf_prob.gpkg")
)

for (item in save_list) {
  sf_obj   <- df_to_sf(item[[1]])
  filename <- item[[2]]
  if (!is.null(sf_obj)) save_thinned(sf_obj, filename)
}

# ── 8. SUB-MODEL DECISION RULE ───────────────────────────────

cat("\n--- Sub-Model Decision Rule ---\n\n")
cat("  N >= 40:   all 6 algorithms\n")
cat("  N = 25-39: MaxEnt + RF + GAM only\n")
cat("  N < 25:    sub-model not viable\n\n")

apply_rule <- function(n, period) {
  if (is.null(n) || n == 0) {
    rule <- "NOT VIABLE (N = 0)"
  } else if (n >= 40) {
    rule <- "All 6 algorithms"
  } else if (n >= 25) {
    rule <- "MaxEnt + RF + GAM only"
  } else {
    rule <- "NOT VIABLE (N < 25)"
  }
  cat(sprintf("  %-8s N = %3d → %s\n", period,
              ifelse(is.null(n), 0, n), rule))
  return(rule)
}

n_pooled <- if (!is.null(thinned_pooled)) nrow(thinned_pooled) else 0
n_LP     <- if (!is.null(thinned_LP))     nrow(thinned_LP)     else 0
n_MP     <- if (!is.null(thinned_MP))     nrow(thinned_MP)     else 0
n_UP     <- if (!is.null(thinned_UP))     nrow(thinned_UP)     else 0

rule_pooled <- apply_rule(n_pooled, "Pooled")
rule_LP     <- apply_rule(n_LP,     "LP")
rule_MP     <- apply_rule(n_MP,     "MP")
rule_UP     <- apply_rule(n_UP,     "UP")

# ── 9. TABLE 1 — SITE COUNT SUMMARY ─────────────────────────

cat("\n--- Table 1: Site Count Summary ---\n\n")

table1 <- data.frame(
  Period         = c("Pooled","LP","MP","UP",
                     "LP_confirmed","MP_confirmed",
                     "LP_conf_prob","MP_conf_prob"),
  N_before_thin  = c(nrow(sites_df), sum(is_LP), sum(is_MP),
                     sum(is_UP),
                     if(!is.null(thinned_LP_conf))
                       nrow(sites_df[is_LP & is_confirmed,])
                     else NA,
                     if(!is.null(thinned_MP_conf))
                       nrow(sites_df[is_MP & is_confirmed,])
                     else NA,
                     if(!is.null(thinned_LP_cp))
                       nrow(sites_df[is_LP & (is_confirmed|is_probable),])
                     else NA,
                     if(!is.null(thinned_MP_cp))
                       nrow(sites_df[is_MP & (is_confirmed|is_probable),])
                     else NA),
  N_after_thin   = c(n_pooled, n_LP, n_MP, n_UP,
                     if(!is.null(thinned_LP_conf))
                       nrow(thinned_LP_conf) else NA,
                     if(!is.null(thinned_MP_conf))
                       nrow(thinned_MP_conf) else NA,
                     if(!is.null(thinned_LP_cp))
                       nrow(thinned_LP_cp)   else NA,
                     if(!is.null(thinned_MP_cp))
                       nrow(thinned_MP_cp)   else NA),
  Decision_rule  = c(rule_pooled, rule_LP, rule_MP, rule_UP,
                     rep("Sensitivity analysis", 4)),
  stringsAsFactors = FALSE
)

write.csv(table1,
          file.path(OUT_TABLES, "Table1_site_counts.csv"),
          row.names = FALSE)

cat(sprintf("  %-18s %12s %12s  %s\n",
            "Period","Before thin","After thin","Decision"))
cat("  ", paste(rep("-",65), collapse=""), "\n")
for (i in seq_len(nrow(table1))) {
  cat(sprintf("  %-18s %12s %12s  %s\n",
              table1$Period[i],
              ifelse(is.na(table1$N_before_thin[i]),
                     "—", table1$N_before_thin[i]),
              ifelse(is.na(table1$N_after_thin[i]),
                     "—", table1$N_after_thin[i]),
              table1$Decision_rule[i]))
}
cat("\n  ✓ Table 1 saved: Table1_site_counts.csv\n\n")

# ── 10. DIAGNOSTIC MAP ───────────────────────────────────────

cat("--- Generating Thinning Diagnostic Map ---\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

png(file.path(OUT_FIG_SUPP, "S0h_thinning_result.png"),
    width = 4800, height = 2400, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))

# Before thinning
terra::plot(template_30m, col="grey95", legend=FALSE, axes=FALSE,
            main=sprintf("Before Thinning\n(N = %d)", nrow(sites_df)))
terra::plot(boundary_vect, add=TRUE, border="grey50", lwd=0.8)
terra::plot(terra::vect(sites_geo), add=TRUE,
            col="steelblue", pch=16, cex=0.3)

# After thinning
if (!is.null(thinned_pooled)) {
  thinned_sf <- df_to_sf(thinned_pooled)
  thinned_geo <- sf::st_transform(thinned_sf, crs=4326)
  terra::plot(template_30m, col="grey95", legend=FALSE, axes=FALSE,
              main=sprintf("After Thinning (1km)\n(N = %d)", n_pooled))
  terra::plot(boundary_vect, add=TRUE, border="grey50", lwd=0.8)
  terra::plot(terra::vect(thinned_geo), add=TRUE,
              col="firebrick", pch=16, cex=0.4)
}

dev.off()
cat("  ✓ S0h_thinning_result.png\n\n")

# ── 11. SUMMARY ─────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 08 COMPLETE — Summary\n")
cat("========================================\n")
cat(sprintf("Pooled:  %3d → %3d  %s\n",
            nrow(sites_df), n_pooled, rule_pooled))
cat(sprintf("LP:      %3d → %3d  %s\n",
            sum(is_LP), n_LP, rule_LP))
cat(sprintf("MP:      %3d → %3d  %s\n",
            sum(is_MP), n_MP, rule_MP))
cat(sprintf("UP:      %3d → %3d  %s\n",
            sum(is_UP), n_UP, rule_UP))
cat("\nNext: Run Script 09 — Bias Correction + Background\n")
cat("========================================\n")