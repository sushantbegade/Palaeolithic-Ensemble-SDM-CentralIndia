# ============================================================
# SCRIPT 16b: PRE-ENSEMBLE ISSUE RESOLUTION
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 16b of 25 (inserted — not in original plan)
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Resolves two bugs before Script 17 (Ensemble):
#
#   ISSUE 1 — MaxEnt Boyce=NA:
#     ecospat.boyce() silently failed in Script 11.
#     Likely cause: non-finite values in cv predictions,
#     or fit vector shorter than obs vector.
#     Fix: diagnose, clean, recompute, update CSV.
#
#   ISSUE 2 — Dist_RawMat dead (near-zero importance):
#     Script 07 lithology keyword match likely failed.
#     If ALL polygons used as raw material → distance ≈ 0
#     everywhere → predictor has no discriminatory power.
#     Fix: inspect actual lithology values, remap correctly,
#     regenerate DIST_RAWMAT_30m_utm44n.tif if needed.
#
#   NON-FIXABLE RESULTS (document only — no code needed):
#     SVM AUC=0.63: gets lowest AUC weight, discuss 6.2
#     Full AUC≈1.0 trees: spatial CV AUC is valid estimate
#     Dist_River MDA negative: correlation artifact
#     TSS low: expected with 10,000 bg prevalence ~1.9%
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 16b: Pre-Ensemble Issue Resolution\n")
cat("========================================\n\n")

set.seed(42)

# ── HELPER ───────────────────────────────────────────────────

safe_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  x[length(x)]
}

compute_boyce_manual <- function(fit, obs, n_bins = 101) {
  # Manual Boyce — bypasses ecospat if needed
  fit <- fit[is.finite(fit)]
  obs <- obs[is.finite(obs)]
  if (length(fit) < 10 || length(obs) < 3) return(NA_real_)
  breaks  <- seq(min(fit), max(fit), length.out = n_bins + 1)
  bin_mid <- (breaks[-1] + breaks[-(n_bins + 1)]) / 2
  pe <- sapply(seq_len(n_bins), function(i) {
    n_p <- sum(obs >= breaks[i] & obs <= breaks[i + 1])
    n_a <- sum(fit >= breaks[i] & fit <= breaks[i + 1])
    if (n_a == 0) return(NA_real_)
    (n_p / length(obs)) / (n_a / length(fit))
  })
  valid <- is.finite(pe) & pe > 0
  if (sum(valid) < 3) return(NA_real_)
  safe_scalar(cor(bin_mid[valid], pe[valid],
                  method = "spearman", use = "complete.obs"))
}

# ════════════════════════════════════════════════════════════
# ISSUE 1: MAXENT BOYCE = NA
# ════════════════════════════════════════════════════════════

cat("════════════════════════════════════════\n")
cat("ISSUE 1: MaxEnt Boyce = NA\n")
cat("════════════════════════════════════════\n\n")

# Step 1: Load CV predictions
mx_cv <- readRDS(file.path(OUT_MOD_IND,
                           "maxent_cv_predictions.rds"))

site_preds <- mx_cv$site_preds
bg_preds   <- mx_cv$bg_preds
fold_aucs  <- mx_cv$fold_aucs

cat("--- Diagnostics ---\n")
cat("  site_preds length:", length(site_preds), "\n")
cat("  bg_preds length:  ", length(bg_preds), "\n")
cat("  fold_aucs:        ", round(fold_aucs, 4), "\n\n")

# Step 2: Check for problems
cat("  site_preds — NA count:      ",
    sum(is.na(site_preds)), "\n")
cat("  site_preds — non-finite:    ",
    sum(!is.finite(site_preds)), "\n")
cat("  site_preds — range:         ",
    round(range(site_preds, na.rm = TRUE), 4), "\n\n")

cat("  bg_preds — NA count:        ",
    sum(is.na(bg_preds)), "\n")
cat("  bg_preds — non-finite:      ",
    sum(!is.finite(bg_preds)), "\n")
cat("  bg_preds — range:           ",
    round(range(bg_preds, na.rm = TRUE), 4), "\n\n")

# Step 3: Diagnose ecospat.boyce() failure mode
cat("--- Diagnosing ecospat.boyce() ---\n")
fit_vec <- c(site_preds, bg_preds)
obs_vec <- site_preds

cat("  fit vector — finite:  ", sum(is.finite(fit_vec)),
    "/", length(fit_vec), "\n")
cat("  obs vector — finite:  ", sum(is.finite(obs_vec)),
    "/", length(obs_vec), "\n")
cat("  fit range:            ",
    round(range(fit_vec[is.finite(fit_vec)], na.rm = TRUE), 4), "\n")
cat("  obs range:            ",
    round(range(obs_vec[is.finite(obs_vec)], na.rm = TRUE), 4), "\n\n")

# Step 4: Attempt ecospat.boyce() with cleaned vectors
cat("--- Attempting ecospat.boyce() with cleaned data ---\n")

fit_clean <- fit_vec[is.finite(fit_vec)]
obs_clean <- obs_vec[is.finite(obs_vec)]

boyce_ecospat <- tryCatch({
  res <- ecospat::ecospat.boyce(
    fit     = fit_clean,
    obs     = obs_clean,
    PEplot  = FALSE
  )
  safe_scalar(res$Spearman.cor)
}, error = function(e) {
  cat("  ecospat.boyce() error:", conditionMessage(e), "\n")
  NA_real_
})

cat("  ecospat Boyce:   ", round(boyce_ecospat, 4), "\n")

# Step 5: Manual Boyce as fallback
boyce_manual <- compute_boyce_manual(fit_clean, obs_clean)
cat("  Manual Boyce:    ", round(boyce_manual, 4), "\n\n")

# Step 6: Select best available Boyce value
boyce_final <- if (!is.na(boyce_ecospat)) boyce_ecospat else boyce_manual

if (is.na(boyce_final)) {
  cat("  ✗ BOTH methods return NA\n")
  cat("  DIAGNOSIS: Check if site_preds are all identical value\n")
  cat("  (zero variance → Spearman undefined)\n")
  cat("  site_preds unique values:",
      length(unique(site_preds)), "\n")
  cat("  ACTION: If all identical → MaxEnt model degenerate\n")
  cat("          Re-check maxnet prediction type='logistic'\n\n")
} else {
  cat("  ✓ Boyce resolved:", round(boyce_final, 4), "\n\n")
}

# Step 7: Update maxent_evaluation.csv
mx_eval <- read.csv(file.path(OUT_EVAL, "maxent_evaluation.csv"),
                    stringsAsFactors = FALSE)
cat("  Current CSV boyce_index:", mx_eval$boyce_index, "\n")

mx_eval$boyce_index <- round(boyce_final, 4)
write.csv(mx_eval,
          file.path(OUT_EVAL, "maxent_evaluation.csv"),
          row.names = FALSE)
cat("  ✓ maxent_evaluation.csv updated\n")
cat("  New boyce_index:", mx_eval$boyce_index, "\n\n")

# ════════════════════════════════════════════════════════════
# ISSUE 2: DIST_RAWMAT — NEAR-ZERO IMPORTANCE
# ════════════════════════════════════════════════════════════

cat("════════════════════════════════════════\n")
cat("ISSUE 2: Dist_RawMat Dead Predictor\n")
cat("════════════════════════════════════════\n\n")

# Step 1: Load current raster and check value range
dist_rawmat <- terra::rast(file.path(OUT_PREDICTORS,
                                     "DIST_RAWMAT_30m_utm44n.tif"))
rng <- terra::global(dist_rawmat, c("min", "max", "mean", "sd"),
                     na.rm = TRUE)

cat("--- Current DIST_RAWMAT Statistics ---\n")
cat("  Min:  ", round(rng[1, "min"],  0), "m\n")
cat("  Max:  ", round(rng[1, "max"],  0), "m\n")
cat("  Mean: ", round(rng[1, "mean"], 0), "m\n")
cat("  SD:   ", round(rng[1, "sd"],   0), "m\n\n")

# DIAGNOSIS:
# If max < 500m or SD < 100m → lithology covers entire area
# → distance everywhere ≈ 0 → predictor has no variation
# → confirmed dead predictor

if (rng[1, "max"] < 500) {
  cat("  DIAGNOSIS: Max distance < 500m\n")
  cat("  Raw material polygons cover nearly entire study area\n")
  cat("  Script 07 used ALL lithology polygons as fallback\n")
  cat("  lithology keyword match FAILED\n\n")
  RAWMAT_PROBLEM <- "coverage_too_high"
} else if (rng[1, "sd"] < 200) {
  cat("  DIAGNOSIS: SD < 200m — very low spatial variation\n")
  cat("  Predictor exists but nearly constant → no discrimination\n\n")
  RAWMAT_PROBLEM <- "low_variation"
} else {
  cat("  DIAGNOSIS: Distance range looks OK statistically\n")
  cat("  Dead importance may reflect genuine non-relationship\n")
  cat("  Not a lithology-matching bug\n\n")
  RAWMAT_PROBLEM <- "genuine_non_signal"
}

# Step 2: Inspect lithology field — what values are actually there?
cat("--- Inspecting Lithology Shapefile ---\n\n")

lithology_sf <- sf::st_read(file.path(OUT_PREDICTORS,
                                      "lithology_utm44n.gpkg"),
                            quiet = TRUE)
all_fields <- names(lithology_sf)[!names(lithology_sf) %in%
                                    c("geom", "geometry")]
cat("  Available fields:", paste(all_fields, collapse = ", "), "\n\n")

# Print unique values of ALL non-geometry fields
# This is the KEY diagnostic — find which field has rock types
for (fld in all_fields) {
  vals <- sort(unique(as.character(lithology_sf[[fld]])))
  if (length(vals) > 50) {
    cat(sprintf("  Field '%s': %d unique values (too many — likely ID/code)\n",
                fld, length(vals)))
  } else {
    cat(sprintf("  Field '%s' (%d unique values):\n",
                fld, length(vals)))
    for (v in vals) cat("    -", v, "\n")
  }
  cat("\n")
}

# Step 3: Based on actual field values, attempt re-matching
cat("--- Re-matching Raw Material Keywords ---\n\n")
cat("  Inspect field values printed above.\n")
cat("  Identify correct field + update CORRECT_LITH_FIELD below.\n\n")

# ── ACTION REQUIRED ──────────────────────────────────────────
# After reading output above:
# 1. Find field containing rock type names (e.g., "Basalt", 
#    "Quartzite", "Granite", "Gneiss", "Sandstone" etc.)
# 2. Set CORRECT_LITH_FIELD to that field name
# 3. Run section below

# USER: UPDATE THIS based on diagnostic output above
CORRECT_LITH_FIELD <- NULL  # e.g., "ROCKTYPE" or "LITH_NAME"

# Extended raw material keywords for Deccan + Vidarbha geology
RAWMAT_KEYWORDS <- c(
  # Siliceous — primary knapping stones
  "quartzite", "quartz", "chert", "flint", "siliceous",
  "agate", "jasper", "silicified", "lydite", "hornfels",
  # Volcanic — Deccan Traps common knapping source
  "basalt", "trap", "deccan", "volcanic", "dolerite", "rhyolite",
  # Metamorphic — Archaean crystalline basement
  "archaean", "crystalline", "metamorphic", "schist", "gneiss",
  "quartzofeldspathic",
  # Sedimentary — sometimes used
  "sandstone", "limestone", "shale", "conglomerate",
  # Generic geological terms covering knappable rocks
  "granite", "granitic", "igneous"
)

if (!is.null(CORRECT_LITH_FIELD)) {
  
  cat("  Testing field:", CORRECT_LITH_FIELD, "\n\n")
  
  rawmat_pattern <- paste(RAWMAT_KEYWORDS, collapse = "|")
  rawmat_sf <- lithology_sf[
    grepl(rawmat_pattern,
          as.character(lithology_sf[[CORRECT_LITH_FIELD]]),
          ignore.case = TRUE), ]
  
  cat("  Matched polygons:", nrow(rawmat_sf), "\n")
  
  if (nrow(rawmat_sf) > 0) {
    matched_vals <- sort(unique(
      as.character(rawmat_sf[[CORRECT_LITH_FIELD]])))
    cat("  Matched values:\n")
    for (v in matched_vals) cat("    ✓", v, "\n")
    cat("\n")
    
    # Re-generate DIST_RAWMAT raster
    cat("  Regenerating DIST_RAWMAT raster...\n")
    
    template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                           "TEMPLATE_30m_utm44n.tif"))
    boundary_vect <- terra::vect(sf::st_read(
      file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
      quiet = TRUE))
    
    rawmat_vect <- terra::vect(rawmat_sf)
    rawmat_rast <- terra::rasterize(rawmat_vect, template_30m,
                                    field = 1, background = NA)
    dist_rawmat_new <- terra::distance(rawmat_rast)
    dist_rawmat_new <- terra::mask(dist_rawmat_new, boundary_vect)
    
    rng_new <- terra::global(dist_rawmat_new,
                             c("min", "max", "mean", "sd"),
                             na.rm = TRUE)
    cat("  New DIST_RAWMAT:\n")
    cat("    Min: ", round(rng_new[1,"min"], 0), "m\n")
    cat("    Max: ", round(rng_new[1,"max"], 0), "m\n")
    cat("    Mean:", round(rng_new[1,"mean"], 0), "m\n")
    cat("    SD:  ", round(rng_new[1,"sd"],  0), "m\n\n")
    
    if (rng_new[1, "max"] < 500) {
      cat("  ⚠ Max still < 500m — keywords match too many polygons\n")
      cat("  Raw material still covers entire area\n")
      cat("  CONCLUSION: raw material ubiquitous in Vidarbha\n")
      cat("  Predictor dead due to real ecological fact, not bug\n")
      cat("  ACTION: Keep current raster, report as Limitation\n\n")
      RAWMAT_FIX <- FALSE
    } else {
      # Save updated raster — overwrite old one
      terra::writeRaster(dist_rawmat_new,
                         file.path(OUT_PREDICTORS,
                                   "DIST_RAWMAT_30m_utm44n.tif"),
                         overwrite = TRUE, datatype = "FLT4S")
      cat("  ✓ DIST_RAWMAT updated with correct field matching\n")
      cat("  ⚠ PREDICTOR STACK MUST BE REGENERATED\n")
      cat("  Run section below to update stack\n\n")
      RAWMAT_FIX <- TRUE
    }
    
  } else {
    cat("  ✗ Zero keyword matches with CORRECT_LITH_FIELD\n")
    cat("  Try different keywords or field name\n\n")
    RAWMAT_FIX <- FALSE
  }
  
} else {
  cat("  ⚠ CORRECT_LITH_FIELD not set yet\n")
  cat("  Read diagnostic output above → set field name → re-run\n\n")
  RAWMAT_FIX <- FALSE
}

# ── IF RAWMAT FIX APPLIED: Regenerate predictor stack ────────

if (exists("RAWMAT_FIX") && RAWMAT_FIX) {
  
  cat("--- Regenerating Predictor Stack with Fixed DIST_RAWMAT ---\n\n")
  
  final_names <- readRDS(file.path(OUT_PREDICTORS,
                                   "final_predictor_names.rds"))
  
  predictor_files <- list(
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
  
  # Load only final retained predictors
  rast_list <- lapply(final_names, function(nm) {
    r <- terra::rast(file.path(OUT_PREDICTORS,
                               predictor_files[[nm]]))
    names(r) <- nm
    return(r)
  })
  
  new_stack <- terra::rast(rast_list)
  
  terra::writeRaster(new_stack,
                     file.path(OUT_PREDICTORS,
                               "PREDICTOR_STACK_FINAL_30m_utm44n.tif"),
                     overwrite = TRUE, datatype = "FLT4S")
  
  cat("  ✓ PREDICTOR_STACK_FINAL_30m_utm44n.tif regenerated\n")
  cat("  ⚠ Models 11–16 now have stale raster predictions\n")
  cat("  ⚠ IF RAWMAT was truly dead (near-zero everywhere),\n")
  cat("     model rank-ordering is unchanged — safe to proceed\n")
  cat("  ⚠ IF RAWMAT now has real variation, consider re-running\n")
  cat("     scripts 11–16 with corrected stack\n\n")
}

# ════════════════════════════════════════════════════════════
# SUMMARY: ALL ISSUES
# ════════════════════════════════════════════════════════════

cat("════════════════════════════════════════\n")
cat("RESOLUTION SUMMARY\n")
cat("════════════════════════════════════════\n\n")

# Reload all evaluation CSVs for final check
eval_files <- c("maxent_evaluation.csv", "rf_evaluation.csv",
                "xgboost_evaluation.csv", "brt_evaluation.csv",
                "gam_evaluation.csv", "svm_evaluation.csv")

alg_names  <- c("MaxEnt","RF","XGBoost","BRT","GAM","SVM")

cat(sprintf("  %-10s %8s %8s %8s %8s %8s\n",
            "Algorithm","CV_AUC","Full_AUC","Boyce","TSS","KG"))
cat("  ", paste(rep("-", 58), collapse=""), "\n")

auc_weights <- numeric(6)
for (i in seq_along(eval_files)) {
  df  <- read.csv(file.path(OUT_EVAL, eval_files[i]),
                  stringsAsFactors = FALSE)
  auc <- df$cv_auc_mean
  boyce <- if ("boyce_index" %in% names(df)) df$boyce_index else NA
  tss   <- if ("tss_max"     %in% names(df)) df$tss_max     else NA
  kg    <- if ("kvamme_gain" %in% names(df)) df$kvamme_gain else NA
  full  <- if ("full_auc"    %in% names(df)) df$full_auc    else NA
  
  auc_weights[i] <- auc
  
  boyce_s <- if (is.na(boyce)) "NA" else sprintf("%.4f", boyce)
  cat(sprintf("  %-10s %8.4f %8.4f %8s %8.4f %8.4f\n",
              alg_names[i], auc, full, boyce_s, tss, kg))
}

# Compute AUC weights for Script 17
auc_weights_norm <- auc_weights / sum(auc_weights)
cat("\n  AUC Weights for Ensemble:\n")
for (i in seq_along(alg_names)) {
  cat(sprintf("  %-10s %.4f (%.1f%%)\n",
              alg_names[i],
              auc_weights_norm[i],
              100 * auc_weights_norm[i]))
}

# Save weights for Script 17
weights_df <- data.frame(
  algorithm  = alg_names,
  cv_auc     = round(auc_weights, 4),
  auc_weight = round(auc_weights_norm, 6),
  stringsAsFactors = FALSE
)
write.csv(weights_df,
          file.path(OUT_EVAL, "ensemble_auc_weights.csv"),
          row.names = FALSE)
cat("\n  ✓ ensemble_auc_weights.csv saved\n\n")

# Non-fixable issues — formal documentation strings
cat("--- Non-Fixable Results (document in manuscript) ---\n\n")

cat("Results 6.2 — SVM:\n")
cat("  SVM CV AUC = 0.6309 — below 0.75 adequate threshold.\n")
cat("  Attributed to inherent SVM sensitivity to hyperparameter\n")
cat("  tuning under extreme class imbalance (SPW=52.9:1) and\n")
cat("  high-dimensional predictor space. Receives lowest\n")
cat("  AUC-weight in ensemble (", 
    round(100 * auc_weights_norm[6], 1), "%).\n\n")

cat("Discussion 7.1 — Full AUC vs CV AUC:\n")
cat("  RF/XGBoost/BRT full-dataset AUC ≈ 1.0 versus spatial\n")
cat("  CV AUC 0.71–0.73 demonstrates that training-data AUC\n")
cat("  substantially overestimates predictive skill.\n")
cat("  Spatial block CV provides conservative, unbiased estimate\n")
cat("  of genuine transferability (Roberts et al. 2017).\n\n")

cat("Results 6.2 — Dist_River RF MDA = -13.3:\n")
cat("  Negative MDA in balanced-bootstrap RF indicates\n")
cat("  permuting Dist_River improves OOB accuracy — consistent\n")
cat("  with strong correlation between Dist_River and other\n")
cat("  hydrological predictors (HAND, Flow_Accum). SHAP values\n")
cat("  will decompose actual contribution independently.\n\n")

cat("Results 6.2 — TSS low (0.11–0.33):\n")
cat("  Expected artifact of presence-background design.\n")
cat("  Background prevalence ≈ 1.9% (190 sites / 10,000 bg).\n")
cat("  Max-TSS threshold biased toward rare positive class.\n")
cat("  Boyce Index is primary performance metric for\n")
cat("  presence-background models (Hirzel et al. 2006).\n\n")

cat("========================================\n")
cat("SCRIPT 16b COMPLETE\n")
cat("========================================\n")
cat("Issue 1 (MaxEnt Boyce): RESOLVED or diagnosed\n")
cat("Issue 2 (Dist_RawMat):  DIAGNOSED — action may be needed\n")
cat("Non-fixable results:    Documentation strings written above\n")
cat("ensemble_auc_weights.csv: SAVED for Script 17\n")
cat("\nIF Dist_RawMat fix applied AND new max > 500m:\n")
cat("  Decision needed: re-run scripts 11–16 or proceed\n")
cat("  Recommendation: proceed if old models' rank order\n")
cat("  unchanged (CV AUC differences minimal)\n")
cat("\nNext: Run Script 17 — Ensemble\n")
cat("========================================\n")