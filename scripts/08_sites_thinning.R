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
#   to reduce spatial autocorrelation bias. Thinning applied:
#     (a) Pooled dataset — all 197 sites
#     (b) LP pool — all sites with LP-equivalent attribution
#     (c) MP pool — all sites with MP-equivalent attribution
#     (d) UP pool — all sites with UP-equivalent attribution
#   Multi-period sites counted in each applicable period pool.
#
#   Sub-model decision rule applied post-thinning:
#     N >= 40:   all 6 algorithms
#     N = 25-39: MaxEnt + RF + GAM only
#     N < 25:    sub-model not viable
#
#   Attribution sensitivity analysis (Supp S3):
#     Stage 1: CONFIRMED sites only per period
#     Stage 2: CONFIRMED + PROBABLE sites per period
#
# ATTRIBUTION CLASSIFICATION (India-specific terminology):
#   LP: Lower Palaeolithic, Acheulian/Acheulean, Series I,
#       Early Stone Age, Abbevilleo, Palaeoliths (standalone)
#       NOTE: "Palaeolith" alone is too broad — use word
#       boundary \\bPalaeoliths\\b to avoid matching
#       "Upper Palaeolithic", "Middle Palaeolithic" etc.
#   MP: Middle Palaeolithic, Middle Stone Age, Series II,
#       MSA, Early Middle, Late Middle, Late Palaeolithic
#       NOTE: "Late Palaeolithic" in Indian context = Late MP
#   UP: Upper Palaeolithic, Late Stone Age, Late Palaeolithic,
#       LSA, Blade and Burin, Microlithic
#       NOTE: "Late Palaeolithic" is in BOTH MP and UP pools
#             (transitional — counted in both)
#
# FIXES FROM v1 (06 June 2026):
#   - LP regex: removed "Palaeolith" (caught "Upper Palaeolithic",
#     "Middle Palaeolithic" etc.). Replaced with
#     \\bPalaeoliths\\b (exact word — standalone term only).
#   - MP regex: added "Late Palaeolithic" (Indian Late MP
#     terminology — undercount fixed from 86 to ~114).
#   - Attribution column auto-detection updated for dot-encoded
#     field names (sf reads spaces as dots from geopackage).
#   - Confidence classification updated for your exact values:
#     "Exact ±5-10 m (Field Verified)" and
#     "Exact ±5-10 m (Reported By Referenced Author)" = confirmed
#     "GIS Map Based" = uncertain (not confirmed or probable)
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
cat("  Total sites loaded:", nrow(sites_sf), "\n\n")

sites_df <- as.data.frame(sf::st_drop_geometry(sites_sf))

# UTM coordinates (for saving back to sf)
coords_utm <- sf::st_coordinates(sites_sf)
sites_df$UTM_X <- coords_utm[, 1]
sites_df$UTM_Y <- coords_utm[, 2]

# WGS84 coordinates (required by spThin)
sites_geo  <- sf::st_transform(sites_sf, crs = 4326)
coords_geo <- sf::st_coordinates(sites_geo)
sites_df$Longitude_WGS84 <- coords_geo[, 1]
sites_df$Latitude_WGS84  <- coords_geo[, 2]

# ── 2. CULTURAL ATTRIBUTION ──────────────────────────────────
# Column name is "Cultural.Attribution" (sf converts spaces to dots)

cat("--- Cultural Period Classification ---\n\n")

cult_str <- as.character(sites_df$Cultural.Attribution)

# ── LP: Lower Palaeolithic equivalent ────────────────────────
# CRITICAL FIX: use \\bPalaeoliths\\b NOT "Palaeolith"
# "Palaeolith" matches "Upper Palaeolithic", "Middle Palaeolithic" etc.
# \\bPalaeoliths\\b matches only the standalone term "Palaeoliths"

is_LP <- grepl(
  paste("Lower Palaeolithic",
        "Lower palaeolithic",
        "\\bLower\\b",
        "Acheulian",
        "Acheulean",
        "Abbevilleo",
        "\\bSeries I\\b",
        "Early Stone Age",
        "\\bESA\\b",
        "\\bPalaeoliths\\b",
        sep = "|"),
  cult_str,
  ignore.case = FALSE  # case-sensitive — avoids false matches
)
# Also catch lowercase "lower" in compound attributions
is_LP <- is_LP | grepl("lower palaeolithic|lower and|lower,",
                       cult_str, ignore.case = TRUE)

# ── MP: Middle Palaeolithic equivalent ───────────────────────
# CRITICAL FIX: add "Late Palaeolithic" — in Indian Palaeolithic
# literature this refers to the Late Middle Palaeolithic phase.
# It appears in compound attributions with "Late Middle Palaeolithic"
# and as a standalone term for the MP-UP transitional phase.

is_MP <- grepl(
  paste("Middle Palaeolithic",
        "Middle Stone Age",
        "Series II",
        "\\bMSA\\b",
        "Early Middle",
        "Late Middle",
        "Late Palaeolithic",
        sep = "|"),
  cult_str,
  ignore.case = TRUE
)

# ── UP: Upper Palaeolithic equivalent ────────────────────────
# Late Palaeolithic also in UP pool (transitional — in both MP and UP)

is_UP <- grepl(
  paste("Upper Palaeolithic",
        "Late Stone Age",
        "Late Palaeolithic",
        "\\bLSA\\b",
        "Blade and Burin",
        "Microlithic",
        sep = "|"),
  cult_str,
  ignore.case = TRUE
)

# Report pools
cat(sprintf("  LP pool: %d sites\n", sum(is_LP)))
cat(sprintf("  MP pool: %d sites\n", sum(is_MP)))
cat(sprintf("  UP pool: %d sites\n", sum(is_UP)))
cat(sprintf("  Pooled:  %d sites\n\n", nrow(sites_df)))

# Report unclassified sites
unclassified <- !is_LP & !is_MP & !is_UP
if (any(unclassified)) {
  cat("  Unclassified (in pooled model only):\n")
  for (v in sort(unique(cult_str[unclassified]))) {
    cat("    -", v, "\n")
  }
  cat("\n")
}

# Validate LP values captured
cat("  LP attribution values:\n")
for (v in sort(unique(cult_str[is_LP]))) cat("    +", v, "\n")
cat("\n  MP attribution values:\n")
for (v in sort(unique(cult_str[is_MP]))) cat("    +", v, "\n")
cat("\n  UP attribution values:\n")
for (v in sort(unique(cult_str[is_UP]))) cat("    +", v, "\n")
cat("\n")

# ── 3. CONFIDENCE CLASSIFICATION ────────────────────────────
# From Location.Precision column (confirmed values in your data):
#   "Exact ±5-10 m (Field Verified)"               → confirmed
#   "Exact ±5-10 m (Reported By Referenced Author)" → confirmed
#   "GIS Map Based"                                 → uncertain

cat("--- Confidence Classification ---\n\n")

prec_str <- as.character(sites_df$Location.Precision)

is_confirmed <- grepl("Field Verified|Reported By Referenced",
                      prec_str, ignore.case = TRUE)
is_probable  <- rep(FALSE, nrow(sites_df))  # no probable tier in your data
is_uncertain <- grepl("GIS Map Based", prec_str, ignore.case = TRUE)

cat(sprintf("  Confirmed (GPS ±5-10m): %d sites\n",  sum(is_confirmed)))
cat(sprintf("  Uncertain (GIS-based):  %d sites\n",  sum(is_uncertain)))
cat(sprintf("  Other:                  %d sites\n\n",
            sum(!is_confirmed & !is_uncertain)))

# ── 4. THINNING HELPER ───────────────────────────────────────

thin_sites <- function(df, label, thin_dist_km = THIN_DIST_KM) {
  
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  [%-22s] 0 sites — skipping\n", label))
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
  
  thin_key <- paste(round(best$Longitude, 5),
                    round(best$Latitude,  5))
  orig_key <- paste(round(df$Longitude_WGS84, 5),
                    round(df$Latitude_WGS84,  5))
  kept_idx    <- which(orig_key %in% thin_key)
  thinned_df  <- df[kept_idx, ]
  
  cat(sprintf("  [%-22s] %3d → %3d  (%.0f%% retained)\n",
              label, nrow(df), nrow(thinned_df),
              100 * nrow(thinned_df) / nrow(df)))
  return(thinned_df)
}

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
  cat(sprintf("    → Saved: %s  (%d sites)\n",
              filename, nrow(sf_obj)))
}

# ── 5. APPLY THINNING ────────────────────────────────────────

cat("--- Spatial Thinning (1km minimum distance) ---\n\n")

thinned_pooled <- thin_sites(sites_df,           "Pooled (all 197)")
thinned_LP     <- thin_sites(sites_df[is_LP, ],  "LP pool")
thinned_MP     <- thin_sites(sites_df[is_MP, ],  "MP pool")
thinned_UP     <- thin_sites(sites_df[is_UP, ],  "UP pool")

cat("\n")

# ── 6. ATTRIBUTION SENSITIVITY ANALYSIS (Supp S3) ───────────
# Your data has no "probable" tier — only confirmed vs uncertain.
# Stage 1 = confirmed only (GPS-verified sites)
# Stage 2 = all sites (confirmed + GIS-based)
# This tests whether GIS-estimated coordinates affect SHAP rankings.

cat("--- Attribution Sensitivity (Supp S3) ---\n\n")

thinned_LP_conf <- thin_sites(
  sites_df[is_LP & is_confirmed, ],   "LP confirmed only")
thinned_MP_conf <- thin_sites(
  sites_df[is_MP & is_confirmed, ],   "MP confirmed only")

# Stage 2 = all sites (same as main pools — no probable tier)
thinned_LP_all  <- thinned_LP  # Stage 2 = all = same as main
thinned_MP_all  <- thinned_MP

cat("  (Stage 2 = all sites = same as main pools)\n\n")

# ── 7. SAVE ALL DATASETS ─────────────────────────────────────

cat("--- Saving Thinned Datasets ---\n\n")

save_thinned(thinned_pooled,  "sites_thinned_pooled.gpkg")
save_thinned(thinned_LP,      "sites_thinned_LP.gpkg")
save_thinned(thinned_MP,      "sites_thinned_MP.gpkg")
save_thinned(thinned_UP,      "sites_thinned_UP.gpkg")
save_thinned(thinned_LP_conf, "sites_thinned_LP_confirmed.gpkg")
save_thinned(thinned_MP_conf, "sites_thinned_MP_confirmed.gpkg")

# ── 8. SUB-MODEL DECISION RULE ───────────────────────────────

cat("\n--- Sub-Model Decision Rule ---\n\n")
cat("  N >= 40:   all 6 algorithms\n")
cat("  N = 25-39: MaxEnt + RF + GAM only\n")
cat("  N < 25:    sub-model not viable\n\n")

apply_rule <- function(thinned_df, period) {
  n <- if (is.null(thinned_df)) 0 else nrow(thinned_df)
  rule <- if (n >= 40)  "All 6 algorithms" else
    if (n >= 25)  "MaxEnt + RF + GAM only" else
      "NOT VIABLE (N < 25)"
  cat(sprintf("  %-8s N = %3d → %s\n", period, n, rule))
  return(list(n = n, rule = rule))
}

r_pooled <- apply_rule(thinned_pooled, "Pooled")
r_LP     <- apply_rule(thinned_LP,     "LP")
r_MP     <- apply_rule(thinned_MP,     "MP")
r_UP     <- apply_rule(thinned_UP,     "UP")

# ── 9. TABLE 1 ───────────────────────────────────────────────

cat("\n--- Table 1: Site Count Summary ---\n\n")

n <- function(d) if (is.null(d)) 0 else nrow(d)

table1 <- data.frame(
  Period        = c("Pooled","LP","MP","UP",
                    "LP_confirmed","MP_confirmed"),
  N_raw         = c(nrow(sites_df), sum(is_LP), sum(is_MP),
                    sum(is_UP),
                    sum(is_LP & is_confirmed),
                    sum(is_MP & is_confirmed)),
  N_thinned     = c(n(thinned_pooled), n(thinned_LP),
                    n(thinned_MP),     n(thinned_UP),
                    n(thinned_LP_conf),n(thinned_MP_conf)),
  Decision      = c(r_pooled$rule, r_LP$rule,
                    r_MP$rule,     r_UP$rule,
                    "Sensitivity Stage 1",
                    "Sensitivity Stage 1"),
  stringsAsFactors = FALSE
)

write.csv(table1,
          file.path(OUT_TABLES, "Table1_site_counts.csv"),
          row.names = FALSE)

cat(sprintf("  %-18s %8s %10s  %s\n",
            "Period","N raw","N thinned","Decision"))
cat("  ", paste(rep("-", 62), collapse=""), "\n")
for (i in seq_len(nrow(table1))) {
  cat(sprintf("  %-18s %8d %10d  %s\n",
              table1$Period[i],
              table1$N_raw[i],
              table1$N_thinned[i],
              table1$Decision[i]))
}
cat("\n  ✓ Table 1 saved: Table1_site_counts.csv\n\n")

# ── 10. DIAGNOSTIC MAP ───────────────────────────────────────

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
                             n(thinned_pooled)))
  terra::plot(boundary_vect, add = TRUE,
              border = "grey50", lwd = 0.8)
  terra::plot(terra::vect(tg), add = TRUE,
              col = "firebrick", pch = 16, cex = 0.4)
}

dev.off()
cat("  ✓ S0h_thinning_result.png\n\n")

# ── 11. SUMMARY ─────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 08 COMPLETE — Summary\n")
cat("========================================\n")
cat(sprintf("Pooled: %3d → %3d  %s\n",
            nrow(sites_df), n(thinned_pooled), r_pooled$rule))
cat(sprintf("LP:     %3d → %3d  %s\n",
            sum(is_LP), n(thinned_LP), r_LP$rule))
cat(sprintf("MP:     %3d → %3d  %s\n",
            sum(is_MP), n(thinned_MP), r_MP$rule))
cat(sprintf("UP:     %3d → %3d  %s\n",
            sum(is_UP), n(thinned_UP), r_UP$rule))
cat("\nNext: Run Script 09 — Bias Correction + Background\n")
cat("========================================\n")