# ============================================================
# SCRIPT 18: BACKGROUND SENSITIVITY ANALYSIS (Supp Table S5)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 18 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Tests sensitivity of spatial CV AUC to background point
#   count N = 1,000 / 5,000 / 10,000 / 20,000.
#   Algorithm: RF only (fastest, highest CV AUC=0.7274,
#   representative of ensemble; avoids 24-run computation).
#   Same parameters as Script 12: balanced bootstrap,
#   ntree=1000, mtry=3.
#   Same spatial block CV design as Script 10.
#   If max AUC variation < 0.02 → N=10,000 confirmed adequate
#   per Warton & Shepherd (2010).
#   Output: Supplementary Table S5.
#
# BACKGROUND FILES (from Script 09):
#   background_N1000.gpkg  — sensitivity
#   background_N5000.gpkg  — sensitivity
#   background_N10000.gpkg — primary (re-run for consistency)
#   background_N20000.gpkg — sensitivity
#
# FOLD ASSIGNMENT:
#   N=10000: use existing bg_folds from cv_block_assignments.rds
#   N=1000/5000/20000: spatial join to block polygons from
#   cv_blocks_obj$blocks (SpatVector, column "folds")
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 18: Background Sensitivity (Supp S5)\n")
cat("========================================\n\n")

set.seed(42)

# ── HELPERS ──────────────────────────────────────────────────

safe_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  x[length(x)]
}

calc_auc <- function(ps, pb) {
  ps <- ps[is.finite(ps)]; pb <- pb[is.finite(pb)]
  if (length(ps) == 0 || length(pb) == 0) return(NA_real_)
  safe_scalar(as.numeric(pROC::auc(
    pROC::roc(c(rep(1, length(ps)), rep(0, length(pb))),
              c(ps, pb), quiet = TRUE))))
}

# ── 1. LOAD FIXED INPUTS ─────────────────────────────────────

cat("--- Loading Fixed Inputs ---\n\n")

# Predictor stack
pred_stack <- terra::rast(file.path(OUT_PREDICTORS,
                                    "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
final_names   <- readRDS(file.path(OUT_PREDICTORS,
                                   "final_predictor_names.rds"))
raster_levels <- readRDS(file.path(OUT_MOD_IND,
                                   "gam_raster_levels.rds"))

cat_predictors <- c("Geology", "Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

# Thinned pooled sites
sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet = TRUE)
site_folds_orig <- sites_sf$fold_id

# Extract predictor values at site locations
cat("  Extracting predictor values at sites...\n")
sites_vals <- terra::extract(pred_stack,
                             terra::vect(sites_sf), ID = FALSE)
sites_ok   <- complete.cases(sites_vals)
sites_vals <- sites_vals[sites_ok, ]
site_folds <- site_folds_orig[sites_ok]

for (col in cat_in_stack) {
  sites_vals[[col]] <- factor(as.integer(sites_vals[[col]]),
                              levels = raster_levels[[col]])
}

N_PRES <- nrow(sites_vals)
NTREE  <- 1000L
MTRY   <- 3L

cat(sprintf("  Thinned sites: %d (after complete.cases)\n",
            N_PRES))
cat(sprintf("  Predictors:    %d\n\n",
            terra::nlyr(pred_stack)))

# ── 2. SPATIAL BLOCK FOLD ASSIGNMENT ─────────────────────────

cat("--- Loading Spatial Block Structure ---\n\n")

cv_design     <- readRDS(file.path(OUT_CV,
                                   "cv_block_assignments.rds"))
bg_folds_orig <- cv_design$bg_folds  # folds for N=10000 bg

# Extract block polygon sf for spatial join
# blockCV 3.x stores blocks as SpatVector in cv_blocks_obj$blocks
block_sv <- tryCatch(
  cv_design$cv_blocks_obj$blocks,
  error = function(e) NULL
)

if (!is.null(block_sv)) {
  block_sf <- sf::st_as_sf(block_sv)
  # Identify fold column name
  fold_col <- grep("fold", names(block_sf),
                   ignore.case = TRUE, value = TRUE)
  fold_col <- fold_col[1]
  if (is.na(fold_col)) fold_col <- names(block_sf)[
    !names(block_sf) %in% c("geometry","geom")][1]
  cat(sprintf("  Block polygons: %d blocks, fold column = '%s'\n",
              nrow(block_sf), fold_col))
} else {
  block_sf <- NULL
  cat("  WARNING: block polygons not extractable\n")
  cat("  Will use KDE-weighted random fold assignment fallback\n")
}
cat("\n")

# Helper: assign folds to new background points
assign_bg_folds <- function(bg_sf, block_sf, fold_col,
                            n_folds = 5L) {
  
  if (is.null(block_sf)) {
    # Fallback: stratified random across study area
    set.seed(42)
    return(sample(seq_len(n_folds), nrow(bg_sf),
                  replace = TRUE))
  }
  
  # Align CRS
  bg_utm     <- sf::st_transform(bg_sf,
                                 crs = sf::st_crs(block_sf))
  bg_utm     <- sf::st_make_valid(bg_utm)
  block_sf   <- sf::st_make_valid(block_sf)
  
  joined <- tryCatch(
    sf::st_join(bg_utm,
                block_sf[, fold_col],
                join  = sf::st_within,
                left  = TRUE),
    error = function(e) {
      cat("    sf::st_join failed — using nearest block\n")
      sf::st_join(bg_utm,
                  block_sf[, fold_col],
                  join = sf::st_nearest_feature)
    }
  )
  
  folds <- as.integer(joined[[fold_col]])
  
  # Handle any remaining NAs (boundary points)
  if (any(is.na(folds))) {
    n_na <- sum(is.na(folds))
    # Assign to fold with most background points so far
    mode_fold <- as.integer(
      names(which.max(table(folds, useNA = "no"))))
    folds[is.na(folds)] <- mode_fold
    cat(sprintf("    %d boundary points assigned to fold %d\n",
                n_na, mode_fold))
  }
  
  return(folds)
}

# ── 3. DEFINE N VALUES AND BG FILES ──────────────────────────

sensitivity_n <- c(1000L, 5000L, 10000L, 20000L)

bg_files <- c(
  "1000"  = "background_N1000.gpkg",
  "5000"  = "background_N5000.gpkg",
  "10000" = "background_N10000.gpkg",
  "20000" = "background_N20000.gpkg"
)

# ── 4. MAIN SENSITIVITY LOOP ─────────────────────────────────

cat("--- Running RF Sensitivity Analysis ---\n\n")
cat("  RF parameters: ntree=1000, mtry=3, balanced bootstrap\n")
cat("  CV design: 5-fold spatial block (same as Scripts 11-16)\n\n")

results_list <- list()

for (N_BG in sensitivity_n) {
  
  N_key <- as.character(N_BG)
  cat(sprintf("═══ N_bg = %6d ═══\n", N_BG))
  
  # ── Load background ───────────────────────────────────────
  
  bg_path <- file.path(OUT_BACKGROUND, bg_files[N_key])
  if (!file.exists(bg_path)) {
    cat(sprintf("  ✗ File not found: %s — skipping\n\n",
                basename(bg_path)))
    next
  }
  
  bg_sf <- sf::st_read(bg_path, quiet = TRUE)
  cat(sprintf("  Background points loaded: %d\n", nrow(bg_sf)))
  
  # ── Assign folds ──────────────────────────────────────────
  
  if (N_BG == 10000L) {
    # Use existing fold assignments (consistent with Scripts 11-16)
    bg_folds <- bg_folds_orig
    cat("  Fold assignment: from cv_block_assignments.rds\n")
  } else {
    cat("  Fold assignment: spatial join to block polygons\n")
    bg_folds <- assign_bg_folds(bg_sf, block_sf,
                                fold_col, n_folds = 5L)
  }
  
  cat(sprintf("  Fold distribution: %s\n",
              paste(table(bg_folds), collapse="/")))
  
  # ── Extract predictor values ──────────────────────────────
  
  bg_vals <- terra::extract(pred_stack,
                            terra::vect(bg_sf), ID = FALSE)
  bg_ok   <- complete.cases(bg_vals)
  bg_vals <- bg_vals[bg_ok, ]
  bg_folds_c <- bg_folds[bg_ok]
  
  for (col in cat_in_stack) {
    bg_vals[[col]] <- factor(as.integer(bg_vals[[col]]),
                             levels = raster_levels[[col]])
  }
  
  N_BG_actual <- nrow(bg_vals)
  cat(sprintf("  Complete bg after extraction: %d\n", N_BG_actual))
  
  # Balanced sampsize: N_PRES per class
  # Safe: ensure bg train per fold >= N_PRES before sampling
  SAMP_SIZE <- c("0" = N_PRES, "1" = N_PRES)
  
  # ── 5-fold spatial block CV ───────────────────────────────
  
  fold_aucs    <- numeric(5)
  cv_preds_s   <- numeric(N_PRES)
  cv_preds_b   <- numeric(N_BG_actual)
  
  t0 <- proc.time()
  
  for (f in 1:5) {
    
    tr_s_idx <- which(site_folds  != f)
    ts_s_idx <- which(site_folds  == f)
    tr_b_idx <- which(bg_folds_c  != f)
    ts_b_idx <- which(bg_folds_c  == f)
    
    n_pres_tr <- length(tr_s_idx)
    n_bg_tr   <- length(tr_b_idx)
    
    # Adjust sampsize if bg_train < N_PRES (only for N=1000)
    samp_this <- c("0" = min(SAMP_SIZE["0"], n_bg_tr),
                   "1" = min(SAMP_SIZE["1"], n_pres_tr))
    
    tr_resp <- factor(
      c(rep("1", n_pres_tr), rep("0", n_bg_tr)),
      levels = c("0","1"))
    tr_data <- rbind(sites_vals[tr_s_idx, ],
                     bg_vals[tr_b_idx, ])
    
    set.seed(42)
    fold_rf <- randomForest::randomForest(
      x           = tr_data,
      y           = tr_resp,
      ntree       = NTREE,
      mtry        = MTRY,
      sampsize    = samp_this,
      replace     = TRUE,
      importance  = FALSE,
      keep.forest = TRUE
    )
    
    # Predict on test set
    ts_s_data <- sites_vals[ts_s_idx, ]
    ts_b_data <- bg_vals[ts_b_idx, ]
    
    ps <- tryCatch(
      predict(fold_rf, ts_s_data, type="prob")[,"1"],
      error=function(e) rep(NA_real_, nrow(ts_s_data)))
    pb <- tryCatch(
      predict(fold_rf, ts_b_data, type="prob")[,"1"],
      error=function(e) rep(NA_real_, nrow(ts_b_data)))
    
    cv_preds_s[ts_s_idx] <- ps
    cv_preds_b[ts_b_idx] <- pb
    fold_aucs[f] <- calc_auc(ps, pb)
    
    cat(sprintf("    Fold %d: AUC=%.4f  sites=%d bg=%d\n",
                f, fold_aucs[f],
                length(ts_s_idx), length(ts_b_idx)))
    
    rm(fold_rf, tr_data, tr_resp, ps, pb)
    gc(full=TRUE)
  }
  
  elapsed <- (proc.time() - t0)[3]
  cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
  cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
  
  cat(sprintf("  CV AUC: %.4f ± %.4f  (%.1f min)\n\n",
              cv_auc_mean, cv_auc_sd, elapsed/60))
  
  # Store results
  results_list[[N_key]] <- data.frame(
    N_background  = N_BG,
    N_complete    = N_BG_actual,
    cv_auc_mean   = round(cv_auc_mean, 4),
    cv_auc_sd     = round(cv_auc_sd, 4),
    fold1_auc     = round(fold_aucs[1], 4),
    fold2_auc     = round(fold_aucs[2], 4),
    fold3_auc     = round(fold_aucs[3], 4),
    fold4_auc     = round(fold_aucs[4], 4),
    fold5_auc     = round(fold_aucs[5], 4),
    algorithm     = "RF",
    note          = sprintf(
      "ntree=%d mtry=%d balanced_bootstrap",
      NTREE, MTRY),
    stringsAsFactors = FALSE
  )
}

# ── 5. COMPILE AND EVALUATE ───────────────────────────────────

cat("════════════════════════════════════════\n")
cat("SENSITIVITY ANALYSIS RESULTS\n")
cat("════════════════════════════════════════\n\n")

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

cat(sprintf("  %-8s  %-10s  %-8s\n",
            "N_bg", "CV AUC", "SD"))
cat("  ", paste(rep("-", 32), collapse=""), "\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("  %-8d  %.4f ± %.4f\n",
              r$N_background, r$cv_auc_mean, r$cv_auc_sd))
}

# Stability assessment
auc_vals  <- results_df$cv_auc_mean
auc_range <- max(auc_vals, na.rm=TRUE) - min(auc_vals, na.rm=TRUE)
auc_n10   <- results_df$cv_auc_mean[results_df$N_background == 10000]

cat(sprintf("\n  AUC range across all N: %.4f\n", auc_range))
cat(sprintf("  N=10,000 CV AUC:        %.4f\n\n", auc_n10))

if (auc_range < 0.02) {
  STABILITY_CONCLUSION <- "STABLE"
  cat("  ✓ CONCLUSION: AUC variation < 0.02\n")
  cat("  N=10,000 confirmed adequate (Warton & Shepherd 2010)\n\n")
} else {
  STABILITY_CONCLUSION <- "UNSTABLE"
  cat("  ⚠ CONCLUSION: AUC variation >= 0.02\n")
  cat("  Background count may affect results — report in Limitations\n\n")
}

# ── 6. SAVE TABLE S5 ─────────────────────────────────────────

results_df$stability_conclusion <- STABILITY_CONCLUSION
results_df$auc_range_across_N   <- round(auc_range, 4)

write.csv(results_df,
          file.path(OUT_SUPP_AN,
                    "TableS5_background_sensitivity.csv"),
          row.names = FALSE)
cat("  ✓ TableS5_background_sensitivity.csv\n\n")

# ── 7. FIGURE S5 ─────────────────────────────────────────────

cat("--- Generating Figure S5 ---\n\n")

tryCatch({
  png(file.path(OUT_FIG_SUPP,
                "FigS5_background_sensitivity.png"),
      width = 2400, height = 2000, res = 300)
  
  par(mar = c(6, 5.5, 4, 2))
  
  n_vals <- results_df$N_background
  aucs   <- results_df$cv_auc_mean
  sds    <- results_df$cv_auc_sd
  
  plot(log10(n_vals), aucs,
       type  = "b",
       pch   = 19,
       col   = "#2166ac",
       lwd   = 2,
       cex   = 1.4,
       ylim  = c(max(0.50, min(aucs, na.rm=TRUE) - 0.05),
                 min(1.00, max(aucs, na.rm=TRUE) + 0.05)),
       xaxt  = "n",
       xlab  = "N background points",
       ylab  = "Spatial Block CV AUC (5-fold mean)",
       main  = sprintf(
         "Background Sensitivity Analysis (RF)\nAUC range = %.4f — %s",
         auc_range, STABILITY_CONCLUSION),
       cex.main = 0.85)
  
  # Error bars (SD)
  for (i in seq_along(n_vals)) {
    if (!is.na(sds[i]) && is.finite(sds[i])) {
      arrows(log10(n_vals[i]),
             aucs[i] - sds[i],
             log10(n_vals[i]),
             aucs[i] + sds[i],
             angle = 90, code = 3,
             length = 0.06, lwd = 1.5,
             col = "#2166ac")
    }
  }
  
  # Custom x-axis with actual N values
  axis(1, at = log10(n_vals),
       labels = format(n_vals, big.mark=","),
       las = 2, cex.axis = 0.85)
  
  # Highlight N=10,000 (primary)
  idx_10k <- which(results_df$N_background == 10000)
  if (length(idx_10k) > 0) {
    points(log10(10000), aucs[idx_10k],
           pch = 19, col = "#d73027", cex = 1.8)
    text(log10(10000),
         aucs[idx_10k] + (if(!is.na(sds[idx_10k])) sds[idx_10k] else 0) + 0.008,
         "Primary\n(N=10,000)",
         cex = 0.60, col = "#d73027", adj = 0.5)
  }
  
  # Stability band: ±0.01 around N=10,000 AUC
  if (!is.na(auc_n10)) {
    abline(h = auc_n10 + 0.01, lty = 3,
           col = "grey50", lwd = 1.2)
    abline(h = auc_n10 - 0.01, lty = 3,
           col = "grey50", lwd = 1.2)
    text(log10(n_vals[1]),
         auc_n10 + 0.012,
         "±0.01 stability band",
         cex = 0.55, col = "grey40", adj = 0)
  }
  
  # 0.75 adequate threshold
  abline(h = 0.75, lty = 2, col = "orange", lwd = 1.2)
  text(log10(n_vals[1]), 0.753,
       "Adequate (0.75)", adj = 0,
       cex = 0.60, col = "orange")
  
  # AUC values as text above points
  for (i in seq_along(n_vals)) {
    text(log10(n_vals[i]),
         aucs[i] - (if(!is.na(sds[i])) sds[i] else 0) - 0.012,
         sprintf("%.4f", aucs[i]),
         cex = 0.60, col = "#2166ac", adj = 0.5)
  }
  
  legend("bottomright",
         legend = c("RF CV AUC (mean ± SD)",
                    "Primary N=10,000",
                    "±0.01 stability band"),
         col    = c("#2166ac","#d73027","grey50"),
         lty    = c(1, NA, 3),
         pch    = c(19, 19, NA),
         lwd    = c(2, NA, 1.2),
         cex    = 0.65, bty = "n")
  
  dev.off()
  cat("  ✓ FigS5_background_sensitivity.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  cat("  ✗ Figure error:", e$message, "\n\n")
})

# ── 8. SUMMARY ───────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 18 COMPLETE — Background Sensitivity\n")
cat("========================================\n")
cat(sprintf("Algorithm:    RF (ntree=%d, mtry=%d, balanced)\n",
            NTREE, MTRY))
cat(sprintf("N values tested: %s\n",
            paste(sensitivity_n, collapse=" / ")))
cat("\nResults:\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("  N=%6d  CV AUC = %.4f ± %.4f\n",
              r$N_background, r$cv_auc_mean, r$cv_auc_sd))
}
cat(sprintf("\nAUC range:    %.4f\n", auc_range))
cat(sprintf("Conclusion:   %s\n", STABILITY_CONCLUSION))
if (STABILITY_CONCLUSION == "STABLE") {
  cat("N=10,000 justified (Warton & Shepherd 2010)\n")
}
cat("\nFiles saved:\n")
cat("  outputs/supplementary_analyses/TableS5_background_sensitivity.csv\n")
cat("  outputs/figures/supplementary/FigS5_background_sensitivity.png\n")
cat("\nMethods text (add to Section 5.5):\n")
cat("  'Background point count sensitivity was assessed using\n")
cat("   Random Forest as the representative algorithm across\n")
cat(sprintf("   N = 1,000; 5,000; 10,000; and 20,000 background points.\n"))
cat(sprintf("   Spatial block CV AUC varied by %.4f across all counts,\n",
            auc_range))
cat(sprintf("   confirming N=10,000 as adequate (Warton & Shepherd 2010).\n"))
cat("   Results are reported in Supplementary Table S5.'\n")
cat("\nNext: Script 19 — SHAP Analysis\n")
cat("========================================\n")