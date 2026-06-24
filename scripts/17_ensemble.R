# ============================================================
# SCRIPT 17: AUC-WEIGHTED ENSEMBLE + UNCERTAINTY (FINAL)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 17 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   1.  Loads all 6 CV prediction vectors + evaluation CSVs
#   2.  Builds AUC-weighted ensemble (primary)
#   3.  Builds equal-weight ensemble (comparator)
#   4.  DeLong's test: weighted vs equal → primary decision
#   5.  Ensemble CV AUC (weighted CV predictions)
#   6.  Ensemble raster AUC (predictions at site/bg locations
#       extracted from full raster — additional diagnostic)
#   7.  DeLong pairwise: ensemble vs each individual algorithm
#   8.  Ensemble Boyce, TSS, Kvamme's Gain
#   9.  Uncertainty surface (SD across 6 logistic rasters)
#  10.  SD quantiles via spatSample (terra::global custom
#       fun NOT supported — fixed)
#  11.  2x2 confidence zone map (suitability × uncertainty)
#  12.  Suitability classification (quartile-based)
#  13.  Saves all rasters, RDS, CSVs
#  14.  Figures 07, 08a, 08b, 08c (all terra plotting fixed)
#
# KEY FIXES vs failed earlier attempts:
#   - terra::contour() inside PNG device → use terra::lines()
#     after extracting contour BEFORE device opens
#   - terra::global(fun=custom) → spatSample + base quantile()
#   - Fig07 SD bars → NA-safe arrows()
#   - All figure closures wrapped in tryCatch
#
# CONFIRMED RESULTS (from Script 17 partial run):
#   Primary ensemble:    AUC-weighted (DeLong p=0.0048)
#   CV AUC:             0.7239
#   Boyce Index:         0.9092
#   TSS:                 0.3238
#   Optimal threshold:   0.23
#   Kvamme's Gain:       0.6063
#   Outperforms:         4/6 individual algorithms
#   SD Q75:              0.2332
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 17: Ensemble + Uncertainty (FINAL)\n")
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
    tp <- sum(obs == 1 & pred >= t, na.rm = TRUE)
    tn <- sum(obs == 0 & pred <  t, na.rm = TRUE)
    fp <- sum(obs == 0 & pred >= t, na.rm = TRUE)
    fn <- sum(obs == 1 & pred <  t, na.rm = TRUE)
    sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    sens + spec - 1
  })
  safe_scalar(max(tss_vals, na.rm = TRUE))
}

get_tss_threshold <- function(obs, pred,
                              thresholds = seq(0.01, 0.99, 0.01)) {
  tss_vals <- sapply(thresholds, function(t) {
    tp <- sum(obs == 1 & pred >= t, na.rm = TRUE)
    tn <- sum(obs == 0 & pred <  t, na.rm = TRUE)
    fp <- sum(obs == 0 & pred >= t, na.rm = TRUE)
    fn <- sum(obs == 1 & pred <  t, na.rm = TRUE)
    sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    sens + spec - 1
  })
  thresholds[which.max(tss_vals)]
}

calc_roc <- function(site_p, bg_p) {
  pROC::roc(
    response  = c(rep(1, length(site_p)), rep(0, length(bg_p))),
    predictor = c(site_p, bg_p),
    quiet     = TRUE
  )
}

calc_kg <- function(raster, sites_vect, threshold) {
  area_cells <- terra::global(!is.na(raster),
                              "sum", na.rm = TRUE)[1, 1]
  high_cells <- terra::global(raster > threshold,
                              "sum", na.rm = TRUE)[1, 1]
  area_pct   <- high_cells / area_cells
  site_pred  <- terra::extract(raster, sites_vect)[, 2]
  sites_pct  <- sum(site_pred > threshold,
                    na.rm = TRUE) / length(site_pred)
  list(kg       = safe_scalar(1 - (area_pct /
                                     max(sites_pct, 1e-9))),
       area_pct = area_pct,
       sites_pct= sites_pct,
       site_pred= site_pred)
}

# Safe raster quantile via spatSample
rast_quantile <- function(r, probs, n = 200000) {
  set.seed(42)
  vals <- terra::spatSample(r, size = n, method = "random",
                            na.rm = TRUE, as.df = TRUE)[[1]]
  quantile(vals, probs = probs, na.rm = TRUE)
}

# ── 1. LOAD ALL INPUTS ───────────────────────────────────────

cat("--- Loading Inputs ---\n\n")

weights_df <- read.csv(file.path(OUT_EVAL,
                                 "ensemble_auc_weights.csv"),
                       stringsAsFactors = FALSE)
alg_names  <- c("MaxEnt","RF","XGBoost","BRT","GAM","SVM")
w          <- weights_df$auc_weight[
  match(alg_names, weights_df$algorithm)]

cat("  AUC weights:\n")
for (i in seq_along(alg_names)) {
  cat(sprintf("  %-10s %.4f (%.1f%%)\n",
              alg_names[i], w[i], 100 * w[i]))
}
cat("\n")

# CV predictions
cv_files <- c("maxent_cv_predictions.rds",
              "rf_cv_predictions.rds",
              "xgboost_cv_predictions.rds",
              "brt_cv_predictions.rds",
              "gam_cv_predictions.rds",
              "svm_cv_predictions.rds")
cv_list <- lapply(cv_files, function(f)
  readRDS(file.path(OUT_MOD_IND, f)))
names(cv_list) <- alg_names

# Probability rasters
raster_files <- c(
  MaxEnt  = "maxent_pred_logistic.tif",
  RF      = "rf_pred_prob.tif",
  XGBoost = "xgboost_pred_prob.tif",
  BRT     = "brt_pred_prob.tif",
  GAM     = "gam_pred_prob.tif",
  SVM     = "svm_pred_prob.tif"
)

rast_list <- lapply(alg_names, function(a) {
  r <- terra::rast(file.path(OUT_MOD_IND, raster_files[a]))
  names(r) <- a; r
})
names(rast_list) <- alg_names
pred_stack_all <- terra::rast(rast_list)

cat("  Individual rasters loaded:\n")
for (a in alg_names) {
  rng <- terra::global(rast_list[[a]],
                       c("min","max","mean"), na.rm = TRUE)
  cat(sprintf("  %-10s [%.4f, %.4f] mean=%.4f\n",
              a, rng[1,1], rng[1,2], rng[1,3]))
}
cat("\n")

sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet = TRUE)
sites_vect    <- terra::vect(sites_sf)
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# ── 2. BUILD ENSEMBLE SURFACES ───────────────────────────────

cat("--- Building Ensemble Surfaces ---\n\n")

ens_weighted <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x * w, na.rm = FALSE)
})
names(ens_weighted) <- "ensemble_auc_weighted"

ens_equal <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = FALSE)
})
names(ens_equal) <- "ensemble_equal_weight"

ens_sd <- terra::app(pred_stack_all, fun = function(x) {
  if (all(is.na(x))) return(NA_real_)
  sd(x, na.rm = FALSE)
})
names(ens_sd) <- "ensemble_sd"

for (r_obj in list(ens_weighted, ens_equal, ens_sd)) {
  rng <- terra::global(r_obj, c("min","max","mean"),
                       na.rm = TRUE)
  cat(sprintf("  %-28s [%.4f, %.4f] mean=%.4f\n",
              names(r_obj), rng[1,1], rng[1,2], rng[1,3]))
}
cat("\n")

# ── 3. ENSEMBLE CV PREDICTIONS ───────────────────────────────

cat("--- Computing Ensemble CV Predictions ---\n\n")

ens_cv_sites_w <- Reduce("+", lapply(seq_along(alg_names),
                                     function(i) cv_list[[alg_names[i]]]$site_preds * w[i]))
ens_cv_bg_w    <- Reduce("+", lapply(seq_along(alg_names),
                                     function(i) cv_list[[alg_names[i]]]$bg_preds  * w[i]))

ens_cv_sites_e <- Reduce("+", lapply(alg_names,
                                     function(a) cv_list[[a]]$site_preds * (1/6)))
ens_cv_bg_e    <- Reduce("+", lapply(alg_names,
                                     function(a) cv_list[[a]]$bg_preds   * (1/6)))

cat(sprintf("  AUC-weighted sites: [%.4f, %.4f]\n",
            min(ens_cv_sites_w), max(ens_cv_sites_w)))
cat(sprintf("  AUC-weighted bg:    [%.4f, %.4f]\n",
            min(ens_cv_bg_w),    max(ens_cv_bg_w)))
cat(sprintf("  Equal-weight sites: [%.4f, %.4f]\n\n",
            min(ens_cv_sites_e), max(ens_cv_sites_e)))

# ── 4. AUC — CV AND RASTER-BASED ─────────────────────────────

cat("--- AUC Computation ---\n\n")

roc_w  <- calc_roc(ens_cv_sites_w, ens_cv_bg_w)
roc_e  <- calc_roc(ens_cv_sites_e, ens_cv_bg_e)
auc_w  <- safe_scalar(as.numeric(pROC::auc(roc_w)))
auc_e  <- safe_scalar(as.numeric(pROC::auc(roc_e)))

cat(sprintf("  CV AUC (AUC-weighted): %.4f\n", auc_w))
cat(sprintf("  CV AUC (equal-weight): %.4f\n", auc_e))
cat(sprintf("  Delta:                 %.4f\n\n",
            auc_w - auc_e))

# Raster-based AUC: extract ensemble predictions at all
# sites and background locations from full raster
# (not CV predictions — full model, diagnostic only)
cat("  Computing raster-extracted AUC (diagnostic)...\n")
raster_pred_sites <- terra::extract(ens_weighted,
                                    sites_vect)[, 2]
# Sample 10000 bg points for raster AUC
bg_sf      <- sf::st_read(file.path(OUT_BACKGROUND,
                                    "background_N10000.gpkg"),
                          quiet = TRUE)
raster_pred_bg <- terra::extract(ens_weighted,
                                 terra::vect(bg_sf))[, 2]
roc_raster <- calc_roc(
  raster_pred_sites[is.finite(raster_pred_sites)],
  raster_pred_bg[is.finite(raster_pred_bg)])
auc_raster <- safe_scalar(as.numeric(pROC::auc(roc_raster)))

cat(sprintf("  Raster-based AUC (full model, diagnostic): %.4f\n\n",
            auc_raster))

# ── 5. DELONG'S TEST: WEIGHTED VS EQUAL ──────────────────────

cat("--- DeLong's Test: Weighted vs Equal-weight ---\n\n")

delong_w_e   <- pROC::roc.test(roc_w, roc_e, method = "delong")
delong_p     <- safe_scalar(delong_w_e$p.value)

cat(sprintf("  p-value: %.4f\n", delong_p))

if (delong_p < 0.05) {
  PRIMARY_ENSEMBLE <- "AUC-weighted"
  PRIMARY_AUC      <- auc_w
  primary_raster   <- ens_weighted
  primary_roc      <- roc_w
  primary_sites    <- ens_cv_sites_w
  primary_bg       <- ens_cv_bg_w
  cat("  AUC-weighted significantly outperforms equal-weight\n")
} else {
  PRIMARY_ENSEMBLE <- "equal-weight"
  PRIMARY_AUC      <- auc_e
  primary_raster   <- ens_equal
  primary_roc      <- roc_e
  primary_sites    <- ens_cv_sites_e
  primary_bg       <- ens_cv_bg_e
  cat("  No significant difference → equal-weight primary\n")
}
cat(sprintf("  PRIMARY: %s (CV AUC = %.4f)\n\n",
            PRIMARY_ENSEMBLE, PRIMARY_AUC))

# ── 6. DELONG'S: ENSEMBLE VS EACH INDIVIDUAL ─────────────────

cat("--- DeLong's: Ensemble vs Each Algorithm ---\n\n")

delong_results <- data.frame(
  comparison   = character(), auc_ensemble = numeric(),
  auc_alg      = numeric(),  delta_auc    = numeric(),
  delong_p     = numeric(),  significant  = logical(),
  stringsAsFactors = FALSE)

for (alg in alg_names) {
  roc_alg <- calc_roc(cv_list[[alg]]$site_preds,
                      cv_list[[alg]]$bg_preds)
  auc_alg <- safe_scalar(as.numeric(pROC::auc(roc_alg)))
  dt  <- tryCatch(
    pROC::roc.test(primary_roc, roc_alg, method="delong"),
    error = function(e) list(p.value = NA_real_))
  p_v <- safe_scalar(dt$p.value)
  sig <- !is.na(p_v) && p_v < 0.05
  cat(sprintf("  vs %-10s AUC_ens=%.4f AUC=%.4f d=%.4f p=%.4f %s\n",
              alg, PRIMARY_AUC, auc_alg,
              PRIMARY_AUC - auc_alg, p_v,
              if(sig) "✓ sig" else "ns"))
  delong_results <- rbind(delong_results, data.frame(
    comparison=paste0("Ensemble vs ",alg),
    auc_ensemble=round(PRIMARY_AUC,4),
    auc_alg=round(auc_alg,4),
    delta_auc=round(PRIMARY_AUC-auc_alg,4),
    delong_p=round(p_v,4), significant=sig,
    stringsAsFactors=FALSE))
}

n_sig <- sum(delong_results$significant, na.rm=TRUE)
cat(sprintf("\n  Ensemble significantly outperforms %d/6\n\n",
            n_sig))

write.csv(delong_results,
          file.path(OUT_EVAL,"delong_pairwise_results.csv"),
          row.names=FALSE)
cat("  ✓ delong_pairwise_results.csv\n\n")

# ── 7. ENSEMBLE EVALUATION METRICS ───────────────────────────

cat("--- Ensemble Evaluation Metrics ---\n\n")

obs_cv  <- c(rep(1, length(primary_sites)),
             rep(0, length(primary_bg)))
pred_cv <- c(primary_sites, primary_bg)

boyce_ens <- compute_boyce_manual(pred_cv, primary_sites)
tss_ens   <- compute_tss(obs_cv, pred_cv)

# TWO THRESHOLDS reported:
# (1) TSS-optimal — maximises sensitivity+specificity
# (2) Conventional 0.5 — standard probability cutoff
optimal_threshold <- get_tss_threshold(obs_cv, pred_cv)

kg_tss  <- calc_kg(primary_raster, sites_vect,
                   optimal_threshold)
kg_05   <- calc_kg(primary_raster, sites_vect, 0.5)

cat(sprintf("  PRIMARY ENSEMBLE:    %s\n", PRIMARY_ENSEMBLE))
cat(sprintf("  CV AUC:              %.4f\n", sm(PRIMARY_AUC)))
cat(sprintf("  Raster AUC (diag):   %.4f\n", sm(auc_raster)))
cat(sprintf("  DeLong p (w vs e):   %.4f\n", delong_p))
cat(sprintf("  Boyce Index:         %.4f\n", sm(boyce_ens)))
cat(sprintf("  TSS (max-TSS):       %.4f\n", sm(tss_ens)))
cat(sprintf("  TSS threshold:       %.2f\n", optimal_threshold))
cat(sprintf("  KG @ TSS threshold:  %.4f\n", sm(kg_tss$kg)))
cat(sprintf("  Area @ TSS thr:      %.1f%%\n",
            100 * kg_tss$area_pct))
cat(sprintf("  Sites @ TSS thr:     %.1f%%\n",
            100 * kg_tss$sites_pct))
cat(sprintf("  KG @ 0.5:            %.4f\n", sm(kg_05$kg)))
cat(sprintf("  Area @ 0.5:          %.1f%%\n",
            100 * kg_05$area_pct))
cat(sprintf("  Sites @ 0.5:         %.1f%%\n",
            100 * kg_05$sites_pct))
cat(sprintf("  Outperforms %d/6 individuals (DeLong p<0.05)\n\n",
            n_sig))

# ── 8. UNCERTAINTY SURFACE STATISTICS ────────────────────────

cat("--- Uncertainty Surface Statistics ---\n\n")

sd_stats <- terra::global(ens_sd, c("min","max","mean","sd"),
                          na.rm=TRUE)

# FIX: use spatSample for quantiles (terra::global custom fun unsupported)
sd_quants  <- rast_quantile(ens_sd, c(0.25, 0.50, 0.75))
sd_q25     <- sd_quants[1]
sd_median  <- sd_quants[2]
sd_q75     <- sd_quants[3]

cat(sprintf("  SD range:   %.4f to %.4f\n",
            sd_stats[1,"min"], sd_stats[1,"max"]))
cat(sprintf("  SD mean:    %.4f  SD: %.4f\n",
            sd_stats[1,"mean"], sd_stats[1,"sd"]))
cat(sprintf("  SD Q25:     %.4f\n", sd_q25))
cat(sprintf("  SD median:  %.4f\n", sd_median))
cat(sprintf("  SD Q75:     %.4f  (high uncertainty threshold)\n\n",
            sd_q75))

# ── 9. SUITABILITY CLASSIFICATION ────────────────────────────

cat("--- Suitability Classification (quartiles) ---\n\n")

suit_quants <- rast_quantile(primary_raster,
                             c(0.25, 0.50, 0.75))
q25 <- suit_quants[1]
q50 <- suit_quants[2]
q75 <- suit_quants[3]

cat(sprintf("  Q25=%.4f  Q50=%.4f  Q75=%.4f\n\n",
            q25, q50, q75))

suit_class <- terra::ifel(
  primary_raster >= q75, 4L,
  terra::ifel(
    primary_raster >= q50, 3L,
    terra::ifel(primary_raster >= q25, 2L, 1L)))
suit_class <- terra::mask(suit_class, boundary_vect)

suit_labels <- c("Low (<Q25)","Moderate (Q25-Q50)",
                 "High (Q50-Q75)","Very High (>Q75)")
suit_freq   <- terra::freq(suit_class)
zone_total  <- sum(suit_freq$count)
site_pred   <- kg_tss$site_pred

cat(sprintf("  %-22s %8s %6s %10s\n",
            "Class","Cells","Area%","Sites%"))
cat("  ", paste(rep("-",52), collapse=""), "\n")
for (i in 1:4) {
  row    <- suit_freq[suit_freq$value == i, ]
  n_cell <- if (nrow(row) > 0) row$count else 0
  ap     <- 100 * n_cell / zone_total
  # Sites in this class
  lo <- c(0, q25, q50, q75)[i]
  hi <- c(q25, q50, q75, 2)[i]
  sp <- 100 * sum(site_pred >= lo & site_pred < hi,
                  na.rm=TRUE) / sum(!is.na(site_pred))
  cat(sprintf("  %-22s %8d %5.1f%% %9.1f%%\n",
              suit_labels[i], n_cell, ap, sp))
}
cat("\n")

# ── 10. 2x2 CONFIDENCE ZONES ─────────────────────────────────

cat("--- 2x2 Confidence Zone Classification ---\n\n")
cat("  Thresholds:\n")
cat(sprintf("  Suitability: >= %.2f (TSS-optimal)\n",
            optimal_threshold))
cat(sprintf("  Uncertainty: SD >= %.4f (Q75)\n\n", sd_q75))

conf_zone <- terra::ifel(
  primary_raster >= optimal_threshold & ens_sd < sd_q75,
  1L,
  terra::ifel(
    primary_raster >= optimal_threshold & ens_sd >= sd_q75,
    2L,
    terra::ifel(
      primary_raster < optimal_threshold & ens_sd < sd_q75,
      3L, 4L)))
conf_zone <- terra::mask(conf_zone, boundary_vect)

zone_labels <- c(
  "Zone 1 — High suit, Low uncert   (Survey Priority A)",
  "Zone 2 — High suit, High uncert  (Verify needed — Priority B)",
  "Zone 3 — Low suit,  Low uncert   (Confirmed low)",
  "Zone 4 — Low suit,  High uncert  (Uncertain low)")
zone_cols_hex <- c("#2166ac","#74add1","#f7f7f7","#d1e5f0")

zone_freq  <- terra::freq(conf_zone)
zone_total2 <- sum(zone_freq$count)

cat(sprintf("  %-54s %6s\n", "Zone", "Area%"))
cat("  ", paste(rep("-",62), collapse=""), "\n")
for (i in 1:4) {
  row <- zone_freq[zone_freq$value == i, ]
  n_c <- if (nrow(row) > 0) row$count else 0
  cat(sprintf("  %s  %5.1f%%\n",
              zone_labels[i], 100 * n_c / zone_total2))
}
cat("\n")

# ── 11. SAVE RASTERS ─────────────────────────────────────────

cat("--- Saving Raster Outputs ---\n\n")

save_rast <- function(r, fname, type = "FLT4S") {
  terra::writeRaster(r,
                     file.path(OUT_MOD_ENS, fname),
                     overwrite = TRUE, datatype = type)
  cat(sprintf("  ✓ %s\n", fname))
}

save_rast(ens_weighted,  "ensemble_auc_weighted.tif")
save_rast(ens_equal,     "ensemble_equal_weight.tif")
save_rast(ens_sd,        "ensemble_uncertainty_sd.tif")
save_rast(ens_weighted,  "ensemble_primary.tif")
save_rast(suit_class,    "ensemble_suitability_class.tif", "INT1U")
save_rast(conf_zone,     "ensemble_confidence_zones.tif",  "INT1U")
cat("\n")

# ── 12. SAVE CV PREDICTIONS ──────────────────────────────────

saveRDS(
  list(site_preds_weighted = ens_cv_sites_w,
       bg_preds_weighted   = ens_cv_bg_w,
       site_preds_equal    = ens_cv_sites_e,
       bg_preds_equal      = ens_cv_bg_e,
       primary_type        = PRIMARY_ENSEMBLE,
       primary_site_preds  = primary_sites,
       primary_bg_preds    = primary_bg,
       optimal_threshold   = optimal_threshold,
       sd_q75_threshold    = sd_q75,
       suit_q25 = q25, suit_q50 = q50, suit_q75 = q75),
  file.path(OUT_MOD_ENS, "ensemble_cv_predictions.rds"))
cat("  ✓ ensemble_cv_predictions.rds\n\n")

# ── 13. SAVE EVALUATION CSV ──────────────────────────────────

ens_eval <- data.frame(
  ensemble_type      = PRIMARY_ENSEMBLE,
  cv_auc             = sm(PRIMARY_AUC),
  raster_auc_diag    = sm(auc_raster),
  auc_weighted_auc   = sm(auc_w),
  equal_weight_auc   = sm(auc_e),
  delong_p_w_vs_e    = round(delong_p, 4),
  boyce_index        = sm(boyce_ens),
  tss_max            = sm(tss_ens),
  tss_threshold      = round(optimal_threshold, 2),
  kvamme_gain_tss    = sm(kg_tss$kg),
  area_pct_tss       = round(100 * kg_tss$area_pct, 2),
  sites_pct_tss      = round(100 * kg_tss$sites_pct, 2),
  kvamme_gain_05     = sm(kg_05$kg),
  area_pct_05        = round(100 * kg_05$area_pct, 2),
  sites_pct_05       = round(100 * kg_05$sites_pct, 2),
  n_sig_outperform   = n_sig,
  sd_mean            = round(sd_stats[1,"mean"], 4),
  sd_max             = round(sd_stats[1,"max"],  4),
  sd_median          = round(sd_median, 4),
  sd_q75             = round(sd_q75, 4),
  suit_q25           = round(q25, 4),
  suit_q50           = round(q50, 4),
  suit_q75           = round(q75, 4),
  stringsAsFactors   = FALSE
)

write.csv(ens_eval,
          file.path(OUT_EVAL, "ensemble_evaluation.csv"),
          row.names = FALSE)
cat("  ✓ ensemble_evaluation.csv\n\n")

# ── 14. FIGURES ──────────────────────────────────────────────
# FIX: terra::contour() inside open device crashes.
# Solution: extract contour lines BEFORE opening device,
# then add with terra::lines() (not terra::plot()).
# All dev.off() calls wrapped in tryCatch.

cat("--- Generating Figures ---\n\n")

# Pre-extract contour at TSS threshold (BEFORE any device opens)
cat("  Pre-extracting contour lines...\n")
hs_contour <- tryCatch(
  terra::as.contour(primary_raster, levels = optimal_threshold),
  error = function(e) {
    cat("  Warning: contour extraction failed —", e$message, "\n")
    NULL
  }
)

# ── FIGURE 8a — Primary ensemble suitability ─────────────────

tryCatch({
  png(file.path(OUT_FIG_MAIN,"Fig08a_ensemble_suitability.png"),
      width = 2400, height = 2400, res = 300)
  terra::plot(primary_raster,
              main = sprintf(
                "%s Ensemble\nCV AUC=%.4f  Boyce=%.4f  KG=%.4f",
                PRIMARY_ENSEMBLE, PRIMARY_AUC,
                sm(boyce_ens), sm(kg_tss$kg)),
              col   = viridisLite::viridis(100),
              range = c(0, 1),
              axes  = FALSE)
  terra::plot(boundary_vect, add=TRUE,
              border="white", lwd=1.2)
  terra::plot(sites_vect, add=TRUE,
              col="red", pch=16, cex=0.35)
  dev.off()
  cat("  ✓ Fig08a_ensemble_suitability.png\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig08a error:", e$message, "\n")
})

# ── FIGURE 8b — Uncertainty + high-suit boundary ─────────────

tryCatch({
  png(file.path(OUT_FIG_MAIN,"Fig08b_ensemble_uncertainty.png"),
      width = 2400, height = 2400, res = 300)
  terra::plot(ens_sd,
              main = sprintf(
                "Ensemble Uncertainty (SD across 6 algorithms)\nmean=%.4f  Q75=%.4f  max=%.4f",
                sd_stats[1,"mean"], sd_q75, sd_stats[1,"max"]),
              col  = viridisLite::magma(100),
              axes = FALSE)
  terra::plot(boundary_vect, add=TRUE,
              border="white", lwd=1.2)
  
  # FIX: use terra::lines() not terra::plot() for contour
  if (!is.null(hs_contour) && length(hs_contour) > 0) {
    terra::lines(hs_contour, col="cyan", lwd=1.8)
  }
  
  terra::plot(sites_vect, add=TRUE,
              col="white", pch=16, cex=0.3)
  legend("bottomright",
         legend = c(sprintf("High suit boundary (>%.2f)",
                            optimal_threshold),
                    "Known sites"),
         col    = c("cyan","white"),
         lty    = c(1, NA), pch = c(NA, 16),
         lwd    = c(1.8, NA), cex = 0.65, bty="n")
  dev.off()
  cat("  ✓ Fig08b_ensemble_uncertainty.png\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig08b error:", e$message, "\n")
})

# ── FIGURE 8c — Confidence zones ─────────────────────────────

zone_cols <- c("#1a6faf","#74b9e0","#f5f5f5","#d6e9f7")

tryCatch({
  png(file.path(OUT_FIG_MAIN,"Fig08c_confidence_zones.png"),
      width = 2400, height = 2400, res = 300)
  terra::plot(conf_zone,
              main   = sprintf(
                "Prospection Confidence Zones\n(suit threshold=%.2f; SD threshold=%.4f)",
                optimal_threshold, sd_q75),
              col    = zone_cols,
              type   = "classes",
              axes   = FALSE,
              legend = FALSE)
  terra::plot(boundary_vect, add=TRUE,
              border="grey40", lwd=1.0)
  terra::plot(sites_vect, add=TRUE,
              col="red", pch=16, cex=0.3)
  
  # Zone area percentages in legend
  zone_pcts <- sapply(1:4, function(i) {
    row <- zone_freq[zone_freq$value == i, ]
    if (nrow(row) > 0) round(100*row$count/zone_total2,1) else 0
  })
  legend("bottomright",
         legend = c(
           sprintf("A: High suit, Low uncert  (%.1f%%)", zone_pcts[1]),
           sprintf("B: High suit, High uncert (%.1f%%)", zone_pcts[2]),
           sprintf("Confirmed low             (%.1f%%)", zone_pcts[3]),
           sprintf("Uncertain low             (%.1f%%)", zone_pcts[4])),
         fill = zone_cols,
         bty  = "n", cex = 0.60)
  dev.off()
  cat("  ✓ Fig08c_confidence_zones.png\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig08c error:", e$message, "\n")
})

# ── FIGURE 7 — Performance comparison ────────────────────────

tryCatch({
  # Compute individual ROC AUCs + SDs
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
      faucs <- cv_list[[a]]$fold_aucs
      if (!is.null(faucs) && length(faucs) > 1)
        safe_scalar(sd(faucs, na.rm=TRUE))
      else NA_real_
    }),
    NA_real_   # No SD for ensemble (composite)
  )
  all_labels <- c(alg_names, "Ensemble\n(AUC-wtd)")
  bar_cols   <- c(
    "#4393c3","#2166ac","#4dac26","#d01c8b","#f1a340","#ca0020",
    "#d73027"  # ensemble in red
  )
  
  png(file.path(OUT_FIG_MAIN,"Fig07_performance_comparison.png"),
      width = 3600, height = 2400, res = 300)
  par(mar = c(8, 5.5, 4, 2))
  
  bp <- barplot(all_aucs,
                names.arg = all_labels,
                col       = bar_cols,
                ylim      = c(0.50, max(all_aucs, na.rm=TRUE) + 0.08),
                xpd       = FALSE,
                ylab      = "Spatial Block CV AUC (5-fold mean)",
                main      = "Individual Algorithm vs Ensemble Performance\n(DeLong's test: * p<0.05)",
                cex.names = 0.80,
                las       = 2,
                border    = NA)
  
  # Error bars — NA-safe
  for (i in seq_along(all_aucs)) {
    if (!is.na(all_sds[i]) && is.finite(all_sds[i])) {
      arrows(bp[i], all_aucs[i] - all_sds[i],
             bp[i], all_aucs[i] + all_sds[i],
             angle=90, code=3, length=0.05, lwd=1.5)
    }
  }
  
  # Reference lines
  abline(h = 0.75, lty = 2, col = "orange",    lwd = 1.5)
  abline(h = 0.85, lty = 2, col = "darkgreen", lwd = 1.5)
  text(bp[1], 0.758, "Adequate (0.75)",
       adj=0, cex=0.62, col="orange")
  text(bp[1], 0.858, "Strong (0.85)",
       adj=0, cex=0.62, col="darkgreen")
  
  # AUC value labels above each bar
  for (i in seq_along(all_aucs)) {
    text(bp[i], all_aucs[i] + 0.003,
         sprintf("%.3f", all_aucs[i]),
         cex=0.60, col="black", adj=0.5)
  }
  
  # DeLong significance stars below bar labels
  for (i in seq_along(alg_names)) {
    comp <- paste0("Ensemble vs ", alg_names[i])
    row  <- delong_results[delong_results$comparison == comp, ]
    if (nrow(row) > 0 && isTRUE(row$significant)) {
      top_y <- all_aucs[i] +
        (if (!is.na(all_sds[i])) all_sds[i] else 0) + 0.015
      text(bp[i], top_y, "*", cex = 1.6, col="red", font=2)
    }
  }
  
  dev.off()
  cat("  ✓ Fig07_performance_comparison.png\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig07 error:", e$message, "\n")
})

cat("\n")

# ── 15. FINAL SUMMARY ────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 17 COMPLETE — Ensemble (FINAL)\n")
cat("========================================\n")
cat(sprintf("PRIMARY ENSEMBLE:    %s\n", PRIMARY_ENSEMBLE))
cat(sprintf("DeLong p (w vs e):   %.4f\n", delong_p))
cat(sprintf("CV AUC:              %.4f\n", sm(PRIMARY_AUC)))
cat(sprintf("Raster AUC (diag):   %.4f\n", sm(auc_raster)))
cat(sprintf("Boyce Index:         %.4f\n", sm(boyce_ens)))
cat(sprintf("TSS (max-TSS):       %.4f\n", sm(tss_ens)))
cat(sprintf("TSS threshold:       %.2f\n", optimal_threshold))
cat(sprintf("KG @ TSS threshold:  %.4f\n", sm(kg_tss$kg)))
cat(sprintf("KG @ 0.5:            %.4f\n", sm(kg_05$kg)))
cat(sprintf("Outperforms:         %d/6 individuals\n", n_sig))
cat(sprintf("SD Q75:              %.4f\n", sd_q75))
cat(sprintf("Suit Q25/Q50/Q75:    %.4f / %.4f / %.4f\n",
            q25, q50, q75))
cat("\nFiles saved:\n")
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
cat("\nNext: Script 18 — Background Sensitivity (Supp S5)\n")
cat("========================================\n")