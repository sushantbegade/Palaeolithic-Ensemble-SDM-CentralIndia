# ============================================================
# SCRIPT 17: AUC-WEIGHTED ENSEMBLE + UNCERTAINTY
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 17 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1. Loads all 6 CV prediction vectors + evaluation CSVs
#   2. Builds AUC-weighted ensemble (primary)
#   3. Builds equal-weight ensemble (comparator)
#   4. DeLong's test: which ensemble is primary?
#   5. Computes ensemble suitability surface (weighted avg
#      of 6 logistic probability rasters)
#   6. Computes uncertainty surface (SD across 6 rasters)
#   7. Classifies suitability: quartile-based 4 classes
#   8. 2x2 confidence zone map (suitability x uncertainty)
#   9. Kvamme's Gain for ensemble
#  10. All evaluation metrics for ensemble
#  11. DeLong pairwise: ensemble vs each individual algorithm
#  12. Saves all outputs for Scripts 18-25
#
# CRITICAL DESIGN NOTES:
#   All 6 algorithms output logistic probability [0,1] —
#   ensemble SD is mathematically valid (same scale).
#   MaxEnt uses type="logistic" throughout (never cloglog).
#   BRT/SVM have compressed means (0.029/0.019) but valid
#   rank-ordering — SD surface will reflect this but
#   ensemble mean is weighted so less impacted.
#
# AUC WEIGHTS (from ensemble_auc_weights.csv):
#   RF=17.6% MaxEnt=17.5% BRT=16.7%
#   XGBoost=16.6% GAM=16.4% SVM=15.1%
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 17: Ensemble + Uncertainty\n")
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

sm <- function(x, d = 4) {
  x <- safe_scalar(x)
  if (is.na(x) || !is.finite(x)) return(NA_real_)
  round(x, d)
}

compute_boyce_manual <- function(fit, obs, n_bins = 101) {
  fit <- fit[is.finite(fit)]; obs <- obs[is.finite(obs)]
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

compute_tss <- function(obs, pred,
                        thresholds = seq(0.01, 0.99, 0.01)) {
  tss_vals <- sapply(thresholds, function(t) {
    tp <- sum(obs==1 & pred>=t, na.rm=TRUE)
    tn <- sum(obs==0 & pred< t, na.rm=TRUE)
    fp <- sum(obs==0 & pred>=t, na.rm=TRUE)
    fn <- sum(obs==1 & pred< t, na.rm=TRUE)
    sens <- if (tp+fn>0) tp/(tp+fn) else NA_real_
    spec <- if (tn+fp>0) tn/(tn+fp) else NA_real_
    sens + spec - 1
  })
  safe_scalar(max(tss_vals, na.rm = TRUE))
}

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n\n")

# AUC weights
weights_df <- read.csv(file.path(OUT_EVAL,
                                 "ensemble_auc_weights.csv"),
                       stringsAsFactors = FALSE)
cat("  AUC weights loaded:\n")
for (i in seq_len(nrow(weights_df))) {
  cat(sprintf("  %-10s %.4f (%.1f%%)\n",
              weights_df$algorithm[i],
              weights_df$auc_weight[i],
              100 * weights_df$auc_weight[i]))
}
cat("\n")

# CV prediction vectors (for DeLong tests)
alg_names <- c("MaxEnt","RF","XGBoost","BRT","GAM","SVM")
cv_files  <- c("maxent_cv_predictions.rds",
               "rf_cv_predictions.rds",
               "xgboost_cv_predictions.rds",
               "brt_cv_predictions.rds",
               "gam_cv_predictions.rds",
               "svm_cv_predictions.rds")

cv_list <- lapply(cv_files, function(f) {
  readRDS(file.path(OUT_MOD_IND, f))
})
names(cv_list) <- alg_names

cat("  CV predictions loaded:\n")
for (alg in alg_names) {
  n_s <- length(cv_list[[alg]]$site_preds)
  n_b <- length(cv_list[[alg]]$bg_preds)
  cat(sprintf("  %-10s %d sites  %d bg\n", alg, n_s, n_b))
}
cat("\n")

# Sites and background spatial objects
sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet = TRUE)
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# Individual probability rasters
raster_files <- c(
  "MaxEnt"   = "maxent_pred_logistic.tif",
  "RF"       = "rf_pred_prob.tif",
  "XGBoost"  = "xgboost_pred_prob.tif",
  "BRT"      = "brt_pred_prob.tif",
  "GAM"      = "gam_pred_prob.tif",
  "SVM"      = "svm_pred_prob.tif"
)

cat("--- Loading Individual Probability Rasters ---\n\n")
rast_list <- list()
for (alg in alg_names) {
  r <- terra::rast(file.path(OUT_MOD_IND, raster_files[alg]))
  names(r) <- alg
  rng <- terra::global(r, c("min","max","mean"), na.rm=TRUE)
  cat(sprintf("  %-10s [%.4f, %.4f]  mean=%.4f\n",
              alg, rng[1,1], rng[1,2], rng[1,3]))
  rast_list[[alg]] <- r
}
pred_stack_all <- terra::rast(rast_list)
cat(sprintf("\n  Stack: %d layers\n\n", terra::nlyr(pred_stack_all)))

# ── 2. BUILD ENSEMBLE PROBABILITY SURFACES ───────────────────

cat("--- Building Ensemble Surfaces ---\n\n")

# Match weight order to raster stack order
w <- weights_df$auc_weight[
  match(alg_names, weights_df$algorithm)]
w_equal <- rep(1/6, 6)

cat("  Computing AUC-weighted ensemble...\n")
# Weighted sum across 6 layers
ens_weighted <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x * w, na.rm = FALSE)
})
names(ens_weighted) <- "ensemble_auc_weighted"

cat("  Computing equal-weight ensemble...\n")
ens_equal <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x * w_equal, na.rm = FALSE)
})
names(ens_equal) <- "ensemble_equal_weight"

cat("  Computing uncertainty surface (SD across 6 rasters)...\n")
ens_sd <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  sd(x, na.rm = FALSE)
})
names(ens_sd) <- "ensemble_sd"

# Report ranges
for (r_obj in list(ens_weighted, ens_equal, ens_sd)) {
  rng <- terra::global(r_obj, c("min","max","mean","sd"),
                       na.rm = TRUE)
  cat(sprintf("  %-28s [%.4f, %.4f]  mean=%.4f  SD=%.4f\n",
              names(r_obj),
              rng[1,1], rng[1,2], rng[1,3], rng[1,4]))
}
cat("\n")

# ── 3. CV PREDICTIONS FOR ENSEMBLE ───────────────────────────
# Weight individual CV predictions for ensemble evaluation.
# Each site/bg point gets a weighted average across algorithms.
# Note: CV predictions are spatially blocked — this is the
# honest performance estimate.

cat("--- Computing Ensemble CV Predictions ---\n\n")

# All algorithms have same site/bg counts (188/9945)
# Verify before combining
n_sites <- length(cv_list[["RF"]]$site_preds)
n_bg    <- length(cv_list[["RF"]]$bg_preds)

cat(sprintf("  Sites: %d  Background: %d\n\n", n_sites, n_bg))

# AUC-weighted ensemble CV predictions
ens_cv_sites_w <- Reduce("+", lapply(seq_along(alg_names), function(i) {
  cv_list[[alg_names[i]]]$site_preds * w[i]
}))

ens_cv_bg_w <- Reduce("+", lapply(seq_along(alg_names), function(i) {
  cv_list[[alg_names[i]]]$bg_preds * w[i]
}))

# Equal-weight ensemble CV predictions
ens_cv_sites_e <- Reduce("+", lapply(alg_names, function(a) {
  cv_list[[a]]$site_preds * (1/6)
}))

ens_cv_bg_e <- Reduce("+", lapply(alg_names, function(a) {
  cv_list[[a]]$bg_preds * (1/6)
}))

cat("  Ensemble CV prediction ranges:\n")
cat(sprintf("  AUC-weighted sites: [%.4f, %.4f]\n",
            min(ens_cv_sites_w), max(ens_cv_sites_w)))
cat(sprintf("  AUC-weighted bg:    [%.4f, %.4f]\n",
            min(ens_cv_bg_w), max(ens_cv_bg_w)))
cat(sprintf("  Equal-weight sites: [%.4f, %.4f]\n",
            min(ens_cv_sites_e), max(ens_cv_sites_e)))
cat("\n")

# ── 4. AUC COMPUTATION ───────────────────────────────────────

cat("--- AUC for Both Ensembles ---\n\n")

calc_roc <- function(site_p, bg_p) {
  pROC::roc(
    response  = c(rep(1, length(site_p)), rep(0, length(bg_p))),
    predictor = c(site_p, bg_p),
    quiet     = TRUE
  )
}

roc_weighted <- calc_roc(ens_cv_sites_w, ens_cv_bg_w)
roc_equal    <- calc_roc(ens_cv_sites_e, ens_cv_bg_e)

auc_w <- safe_scalar(as.numeric(pROC::auc(roc_weighted)))
auc_e <- safe_scalar(as.numeric(pROC::auc(roc_equal)))

cat(sprintf("  AUC-weighted ensemble CV AUC: %.4f\n", auc_w))
cat(sprintf("  Equal-weight ensemble CV AUC: %.4f\n", auc_e))
cat(sprintf("  Delta AUC (weighted - equal): %.4f\n\n",
            auc_w - auc_e))

# ── 5. DELONG'S TEST: WEIGHTED VS EQUAL ──────────────────────

cat("--- DeLong's Test: AUC-weighted vs Equal-weight ---\n\n")

delong_w_vs_e <- pROC::roc.test(roc_weighted, roc_equal,
                                method = "delong")

delong_p <- delong_w_vs_e$p.value
cat(sprintf("  DeLong p-value: %.4f\n", delong_p))

if (delong_p < 0.05) {
  PRIMARY_ENSEMBLE <- "AUC-weighted"
  PRIMARY_AUC      <- auc_w
  primary_sites    <- ens_cv_sites_w
  primary_bg       <- ens_cv_bg_w
  primary_roc      <- roc_weighted
  primary_raster   <- ens_weighted
  cat("  p < 0.05 → AUC-weighted significantly better\n")
  cat("  PRIMARY ENSEMBLE: AUC-weighted\n\n")
} else {
  PRIMARY_ENSEMBLE <- "equal-weight"
  PRIMARY_AUC      <- auc_e
  primary_sites    <- ens_cv_sites_e
  primary_bg       <- ens_cv_bg_e
  primary_roc      <- roc_equal
  primary_raster   <- ens_equal
  cat("  p >= 0.05 → no significant difference\n")
  cat("  PRIMARY ENSEMBLE: equal-weight (simpler, preferred)\n\n")
}

# ── 6. DELONG'S TEST: ENSEMBLE VS EACH INDIVIDUAL ────────────

cat("--- DeLong's Test: Ensemble vs Each Algorithm ---\n\n")

delong_results <- data.frame(
  comparison   = character(),
  auc_ensemble = numeric(),
  auc_alg      = numeric(),
  delta_auc    = numeric(),
  delong_p     = numeric(),
  significant  = logical(),
  stringsAsFactors = FALSE
)

for (alg in alg_names) {
  sp <- cv_list[[alg]]$site_preds
  bp <- cv_list[[alg]]$bg_preds
  roc_alg <- calc_roc(sp, bp)
  auc_alg <- safe_scalar(as.numeric(pROC::auc(roc_alg)))
  
  dt <- tryCatch(
    pROC::roc.test(primary_roc, roc_alg, method = "delong"),
    error = function(e) list(p.value = NA_real_)
  )
  p_val <- safe_scalar(dt$p.value)
  sig   <- !is.na(p_val) && p_val < 0.05
  
  cat(sprintf("  Ensemble vs %-10s AUC_ens=%.4f  AUC_alg=%.4f  delta=%.4f  p=%.4f  %s\n",
              alg, PRIMARY_AUC, auc_alg,
              PRIMARY_AUC - auc_alg, p_val,
              if (sig) "✓ sig" else "ns"))
  
  delong_results <- rbind(delong_results, data.frame(
    comparison   = paste0("Ensemble vs ", alg),
    auc_ensemble = round(PRIMARY_AUC, 4),
    auc_alg      = round(auc_alg, 4),
    delta_auc    = round(PRIMARY_AUC - auc_alg, 4),
    delong_p     = round(p_val, 4),
    significant  = sig,
    stringsAsFactors = FALSE
  ))
}

# Count significant outperformances
n_sig <- sum(delong_results$significant, na.rm = TRUE)
cat(sprintf("\n  Ensemble significantly outperforms %d/6 algorithms\n\n",
            n_sig))

write.csv(delong_results,
          file.path(OUT_EVAL, "delong_pairwise_results.csv"),
          row.names = FALSE)
cat("  ✓ delong_pairwise_results.csv\n\n")

# ── 7. ENSEMBLE EVALUATION METRICS ───────────────────────────

cat("--- Ensemble Evaluation Metrics ---\n\n")

# Boyce Index
boyce_ens <- compute_boyce_manual(
  c(primary_sites, primary_bg), primary_sites)

# TSS
tss_ens <- compute_tss(
  obs  = c(rep(1, length(primary_sites)),
           rep(0, length(primary_bg))),
  pred = c(primary_sites, primary_bg))

# Kvamme's Gain from raster
area_cells <- terra::global(!is.na(primary_raster),
                            "sum", na.rm=TRUE)[1,1]

# Threshold via max-TSS
obs_vec  <- c(rep(1, length(primary_sites)),
              rep(0, length(primary_bg)))
pred_vec <- c(primary_sites, primary_bg)
thresholds <- seq(0.01, 0.99, 0.01)
tss_vals   <- sapply(thresholds, function(t) {
  tp <- sum(obs_vec==1 & pred_vec>=t, na.rm=TRUE)
  tn <- sum(obs_vec==0 & pred_vec< t, na.rm=TRUE)
  fp <- sum(obs_vec==0 & pred_vec>=t, na.rm=TRUE)
  fn <- sum(obs_vec==1 & pred_vec< t, na.rm=TRUE)
  sens <- if(tp+fn>0) tp/(tp+fn) else NA_real_
  spec <- if(tn+fp>0) tn/(tn+fp) else NA_real_
  sens + spec - 1
})
optimal_threshold <- thresholds[which.max(tss_vals)]

high_cells <- terra::global(primary_raster > optimal_threshold,
                            "sum", na.rm=TRUE)[1,1]
area_pct   <- high_cells / area_cells

site_pred  <- terra::extract(primary_raster,
                             terra::vect(sites_sf))[, 2]
sites_pct  <- sum(site_pred > optimal_threshold,
                  na.rm=TRUE) / length(site_pred)
kg <- safe_scalar(1 - (area_pct / max(sites_pct, 1e-9)))

cat(sprintf("  PRIMARY ENSEMBLE:  %s\n", PRIMARY_ENSEMBLE))
cat(sprintf("  CV AUC:            %.4f\n", sm(PRIMARY_AUC)))
cat(sprintf("  DeLong p (w vs e): %.4f\n", delong_p))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_ens)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_ens)))
cat(sprintf("  Optimal threshold: %.2f\n", optimal_threshold))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Area > threshold:  %.1f%%\n", 100 * area_pct))
cat(sprintf("  Sites > threshold: %.1f%%\n", 100 * sites_pct))
cat(sprintf("  Outperforms %d/6 individuals (DeLong p<0.05)\n\n",
            n_sig))

# ── 8. UNCERTAINTY SURFACE STATISTICS ────────────────────────

cat("--- Uncertainty Surface Statistics ---\n\n")

sd_stats  <- terra::global(ens_sd, c("min","max","mean","sd"),
                           na.rm = TRUE)
cat(sprintf("  SD range:  %.4f to %.4f\n",
            sd_stats[1,"min"], sd_stats[1,"max"]))
cat(sprintf("  SD mean:   %.4f\n", sd_stats[1,"mean"]))
cat(sprintf("  SD median: %.4f\n",
            safe_scalar(terra::global(ens_sd, "median",
                                      na.rm=TRUE)[1,1])))

# High uncertainty threshold = upper quartile of SD
sd_q75 <- safe_scalar(terra::global(
  ens_sd, fun=function(x) quantile(x, 0.75, na.rm=TRUE))[1,1])
cat(sprintf("  SD Q75 (high uncertainty threshold): %.4f\n\n",
            sd_q75))

# High suitability threshold = optimal_threshold from above
# 2x2 confidence zone classification:
#   Zone 1: High suit + Low uncertainty  = HIGH CONFIDENCE PRIORITY
#   Zone 2: High suit + High uncertainty = UNCERTAIN PRIORITY
#   Zone 3: Low suit  + Low uncertainty  = CONFIRMED LOW
#   Zone 4: Low suit  + High uncertainty = UNCERTAIN LOW

cat("--- 2x2 Confidence Zone Classification ---\n\n")

conf_zone <- terra::ifel(
  primary_raster >= optimal_threshold & ens_sd < sd_q75,  1L,
  terra::ifel(
    primary_raster >= optimal_threshold & ens_sd >= sd_q75, 2L,
    terra::ifel(
      primary_raster < optimal_threshold & ens_sd < sd_q75,  3L,
      4L
    )
  )
)
conf_zone <- terra::mask(conf_zone, boundary_vect)

zone_labels <- c("High suit + Low uncert (Priority A)",
                 "High suit + High uncert (Priority B)",
                 "Low suit  + Low uncert  (Confirmed low)",
                 "Low suit  + High uncert (Uncertain low)")

zone_freq <- terra::freq(conf_zone)
zone_total <- sum(zone_freq$count)

cat(sprintf("  %-40s %8s %6s\n","Zone","Cells","Area%"))
cat("  ", paste(rep("-",58),collapse=""), "\n")
for (i in 1:4) {
  row <- zone_freq[zone_freq$value == i, ]
  if (nrow(row) > 0) {
    cat(sprintf("  %d. %-38s %8d %5.1f%%\n",
                i, zone_labels[i], row$count,
                100 * row$count / zone_total))
  }
}
cat("\n")

# ── 9. SUITABILITY CLASSIFICATION (quartile-based) ───────────

cat("--- Suitability Classification (quartiles) ---\n\n")

# Q25, Q50, Q75 of ensemble surface
q25 <- safe_scalar(terra::global(
  primary_raster, fun=function(x) quantile(x,0.25,na.rm=TRUE))[1,1])
q50 <- safe_scalar(terra::global(
  primary_raster, fun=function(x) quantile(x,0.50,na.rm=TRUE))[1,1])
q75 <- safe_scalar(terra::global(
  primary_raster, fun=function(x) quantile(x,0.75,na.rm=TRUE))[1,1])

cat(sprintf("  Q25: %.4f  Q50: %.4f  Q75: %.4f\n\n",
            q25, q50, q75))

suit_class <- terra::ifel(
  primary_raster >= q75, 4L,
  terra::ifel(
    primary_raster >= q50, 3L,
    terra::ifel(
      primary_raster >= q25, 2L, 1L
    )
  )
)
suit_class <- terra::mask(suit_class, boundary_vect)

suit_labels <- c("Low (<Q25)","Moderate (Q25-Q50)",
                 "High (Q50-Q75)","Very High (>Q75)")
suit_freq   <- terra::freq(suit_class)

cat(sprintf("  %-22s %8s %6s %10s\n",
            "Class","Cells","Area%","Sites%"))
cat("  ", paste(rep("-",52),collapse=""), "\n")

for (i in 1:4) {
  row    <- suit_freq[suit_freq$value == i, ]
  n_cell <- if (nrow(row)>0) row$count else 0
  ap     <- 100 * n_cell / zone_total
  # Sites in this class
  n_sites_cls <- sum(site_pred >= c(0,q25,q50,q75)[i] &
                       site_pred <  c(q25,q50,q75,2)[i],
                     na.rm=TRUE)
  sp <- 100 * n_sites_cls / length(site_pred)
  cat(sprintf("  %-22s %8d %5.1f%% %9.1f%%\n",
              suit_labels[i], n_cell, ap, sp))
}
cat("\n")

# ── 10. SAVE ALL RASTER OUTPUTS ──────────────────────────────

cat("--- Saving Raster Outputs ---\n\n")

# Primary ensemble
terra::writeRaster(
  primary_raster,
  file.path(OUT_MOD_ENS, "ensemble_primary.tif"),
  overwrite=TRUE, datatype="FLT4S")
cat("  ✓ ensemble_primary.tif\n")

# Both ensembles
terra::writeRaster(
  ens_weighted,
  file.path(OUT_MOD_ENS, "ensemble_auc_weighted.tif"),
  overwrite=TRUE, datatype="FLT4S")
cat("  ✓ ensemble_auc_weighted.tif\n")

terra::writeRaster(
  ens_equal,
  file.path(OUT_MOD_ENS, "ensemble_equal_weight.tif"),
  overwrite=TRUE, datatype="FLT4S")
cat("  ✓ ensemble_equal_weight.tif\n")

# Uncertainty
terra::writeRaster(
  ens_sd,
  file.path(OUT_MOD_ENS, "ensemble_uncertainty_sd.tif"),
  overwrite=TRUE, datatype="FLT4S")
cat("  ✓ ensemble_uncertainty_sd.tif\n")

# Classified surfaces
terra::writeRaster(
  suit_class,
  file.path(OUT_MOD_ENS, "ensemble_suitability_class.tif"),
  overwrite=TRUE, datatype="INT1U")
cat("  ✓ ensemble_suitability_class.tif\n")

terra::writeRaster(
  conf_zone,
  file.path(OUT_MOD_ENS, "ensemble_confidence_zones.tif"),
  overwrite=TRUE, datatype="INT1U")
cat("  ✓ ensemble_confidence_zones.tif\n\n")

# ── 11. SAVE CV PREDICTIONS ──────────────────────────────────

saveRDS(
  list(site_preds_weighted = ens_cv_sites_w,
       bg_preds_weighted   = ens_cv_bg_w,
       site_preds_equal    = ens_cv_sites_e,
       bg_preds_equal      = ens_cv_bg_e,
       primary_type        = PRIMARY_ENSEMBLE,
       primary_site_preds  = primary_sites,
       primary_bg_preds    = primary_bg,
       optimal_threshold   = optimal_threshold,
       sd_q75_threshold    = sd_q75),
  file.path(OUT_MOD_ENS, "ensemble_cv_predictions.rds")
)
cat("  ✓ ensemble_cv_predictions.rds\n\n")

# ── 12. SAVE EVALUATION TABLE ────────────────────────────────

cat("--- Saving Ensemble Evaluation ---\n\n")

ens_eval <- data.frame(
  ensemble_type      = PRIMARY_ENSEMBLE,
  cv_auc             = sm(PRIMARY_AUC),
  auc_weighted_auc   = sm(auc_w),
  equal_weight_auc   = sm(auc_e),
  delong_p_w_vs_e    = round(delong_p, 4),
  boyce_index        = sm(boyce_ens),
  tss_max            = sm(tss_ens),
  optimal_threshold  = round(optimal_threshold, 2),
  kvamme_gain        = sm(kg),
  area_pct_high      = round(100 * area_pct, 2),
  sites_pct_high     = round(100 * sites_pct, 2),
  n_sig_outperform   = n_sig,
  sd_mean            = round(sd_stats[1,"mean"], 4),
  sd_max             = round(sd_stats[1,"max"],  4),
  sd_q75             = round(sd_q75, 4),
  stringsAsFactors   = FALSE
)

write.csv(ens_eval,
          file.path(OUT_EVAL, "ensemble_evaluation.csv"),
          row.names = FALSE)
cat("  ✓ ensemble_evaluation.csv\n\n")

# ── 13. FIGURES ──────────────────────────────────────────────

cat("--- Generating Figures ---\n\n")

# Figure 8a — Primary ensemble suitability
png(file.path(OUT_FIG_MAIN, "Fig08a_ensemble_suitability.png"),
    width=2400, height=2400, res=300)
terra::plot(primary_raster,
            main = sprintf(
              "%s Ensemble — Suitability\nCV AUC=%.4f  Boyce=%.4f  KG=%.4f",
              PRIMARY_ENSEMBLE, PRIMARY_AUC,
              sm(boyce_ens), sm(kg)),
            col  = viridisLite::viridis(100),
            range = c(0, 1), axes = FALSE)
terra::plot(boundary_vect, add=TRUE, border="white", lwd=0.8)
terra::plot(terra::vect(sites_sf), add=TRUE,
            col="red", pch=16, cex=0.35)
dev.off()
cat("  ✓ Fig08a_ensemble_suitability.png\n")

# Figure 8b — Uncertainty surface
png(file.path(OUT_FIG_MAIN, "Fig08b_ensemble_uncertainty.png"),
    width=2400, height=2400, res=300)
terra::plot(ens_sd,
            main = sprintf(
              "Ensemble Uncertainty (SD)\nmean=%.4f  max=%.4f  Q75=%.4f",
              sd_stats[1,"mean"], sd_stats[1,"max"], sd_q75),
            col  = viridisLite::magma(100),
            axes = FALSE)
terra::plot(boundary_vect, add=TRUE, border="white", lwd=0.8)
# Overlay high-suitability contour
hs_contour <- terra::as.contour(
  primary_raster, levels=optimal_threshold)
terra::plot(hs_contour, add=TRUE, col="cyan", lwd=1.2)
dev.off()
cat("  ✓ Fig08b_ensemble_uncertainty.png\n")

# Figure 8c — Confidence zone map
zone_cols <- c("#2166ac","#74add1","#fee090","#f46d43")
png(file.path(OUT_FIG_MAIN, "Fig08c_confidence_zones.png"),
    width=2400, height=2400, res=300)
terra::plot(conf_zone,
            main  = "Prospection Confidence Zones\n(Suitability × Uncertainty)",
            col   = zone_cols,
            type  = "classes",
            axes  = FALSE,
            legend= FALSE)
terra::plot(boundary_vect, add=TRUE, border="white", lwd=0.8)
terra::plot(terra::vect(sites_sf), add=TRUE,
            col="white", pch=16, cex=0.3)
legend("bottomright",
       legend = c("Priority A: High suit, Low uncert",
                  "Priority B: High suit, High uncert",
                  "Confirmed low suitability",
                  "Uncertain low suitability"),
       fill   = zone_cols,
       bty    = "n", cex = 0.55)
dev.off()
cat("  ✓ Fig08c_confidence_zones.png\n")

# Figure 7 — Performance comparison bar chart
png(file.path(OUT_FIG_MAIN, "Fig07_performance_comparison.png"),
    width=3600, height=2400, res=300)

all_aucs <- c(
  sapply(alg_names, function(a) {
    safe_scalar(as.numeric(pROC::auc(
      calc_roc(cv_list[[a]]$site_preds,
               cv_list[[a]]$bg_preds))))
  }),
  PRIMARY_AUC
)
all_sds <- c(
  sapply(alg_names, function(a) {
    safe_scalar(sd(cv_list[[a]]$fold_aucs, na.rm=TRUE))
  }),
  NA_real_
)
all_labels <- c(alg_names, paste0("Ensemble\n(", PRIMARY_ENSEMBLE, ")"))
bar_cols   <- c(rep("#74add1", 6), "#d73027")

par(mar = c(6, 5, 4, 2))
bp <- barplot(all_aucs,
              names.arg = all_labels,
              col       = bar_cols,
              ylim      = c(0.5, 0.85),
              xpd       = FALSE,
              ylab      = "Spatial Block CV AUC",
              main      = "Algorithm Performance Comparison\n(5-fold Spatial Block CV AUC ± SD)",
              cex.names = 0.75,
              las       = 2)

# Error bars (SD)
for (i in seq_along(all_aucs)) {
  if (!is.na(all_sds[i])) {
    arrows(bp[i], all_aucs[i] - all_sds[i],
           bp[i], all_aucs[i] + all_sds[i],
           angle=90, code=3, length=0.05, lwd=1.5)
  }
}

# Threshold lines
abline(h=0.75, lty=2, col="orange", lwd=1.5)
abline(h=0.85, lty=2, col="darkgreen", lwd=1.5)
text(0.5, 0.755, "Adequate (0.75)", adj=0, cex=0.65, col="orange")
text(0.5, 0.855, "Strong (0.85)",   adj=0, cex=0.65, col="darkgreen")

# DeLong significance annotations
for (i in seq_along(alg_names)) {
  row <- delong_results[delong_results$comparison ==
                          paste0("Ensemble vs ", alg_names[i]), ]
  if (nrow(row) > 0 && row$significant) {
    text(bp[i], all_aucs[i] + (if(!is.na(all_sds[i])) all_sds[i] else 0) + 0.01,
         "*", cex=1.2, col="red")
  }
}

dev.off()
cat("  ✓ Fig07_performance_comparison.png\n\n")

# ── 14. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 17 COMPLETE — Ensemble\n")
cat("========================================\n")
cat(sprintf("Primary ensemble:     %s\n", PRIMARY_ENSEMBLE))
cat(sprintf("DeLong p (w vs e):    %.4f\n", delong_p))
cat(sprintf("Ensemble CV AUC:      %.4f\n", sm(PRIMARY_AUC)))
cat(sprintf("Boyce Index:          %.4f\n", sm(boyce_ens)))
cat(sprintf("TSS:                  %.4f\n", sm(tss_ens)))
cat(sprintf("Kvamme's Gain:        %.4f\n", sm(kg)))
cat(sprintf("Outperforms %d/6 individuals (DeLong p<0.05)\n",
            n_sig))
cat(sprintf("Optimal threshold:    %.2f\n", optimal_threshold))
cat(sprintf("SD mean:              %.4f\n", sd_stats[1,"mean"]))
cat(sprintf("SD Q75:               %.4f\n", sd_q75))
cat("\nOutputs saved:\n")
cat("  ensemble_primary.tif\n")
cat("  ensemble_auc_weighted.tif\n")
cat("  ensemble_equal_weight.tif\n")
cat("  ensemble_uncertainty_sd.tif\n")
cat("  ensemble_suitability_class.tif\n")
cat("  ensemble_confidence_zones.tif\n")
cat("  ensemble_cv_predictions.rds\n")
cat("  ensemble_evaluation.csv\n")
cat("  delong_pairwise_results.csv\n")
cat("  Fig07_performance_comparison.png\n")
cat("  Fig08a_ensemble_suitability.png\n")
cat("  Fig08b_ensemble_uncertainty.png\n")
cat("  Fig08c_confidence_zones.png\n")
cat("\nNext: Run Script 18 — Background Sensitivity Analysis\n")
cat("  (Supplementary Table S5 — N=1000/5000/10000/20000)\n")
cat("========================================\n")