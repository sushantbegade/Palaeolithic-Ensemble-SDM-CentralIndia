# ============================================================
# SCRIPT 21: GEOGRAPHIC TRANSFER VALIDATION
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 21 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Geographic transfer validation following Research Design
#   Section 5.13. Southern Chandrapur block (Gondpipri-
#   Korpana-Rajura talukas, ~3,500 km²) held out. Two
#   fastest top-performing algorithms (RF + MaxEnt) retrained
#   on northern ~17,800 km² only. Predictions generated for
#   held-out southern sector. Transfer AUC and Boyce Index
#   compared against spatial CV AUC from ensemble.
#
# METHOD:
#   Transfer zone: bottom 16.4% of study area y-extent
#   (3500/21300 = 0.164) = southern Chandrapur sector.
#   Override: set TALUKA_SHP path if taluka shapefile exists.
#
#   Transfer algorithms: RF + MaxEnt (top-2 CV AUC, fastest)
#   Transfer ensemble: equal-weight average of RF + MaxEnt
#
# METRICS:
#   Transfer AUC (pROC)
#   Transfer Boyce Index (manual Spearman)
#   Delta AUC = |Transfer AUC - Spatial CV AUC|
#   If delta > 0.10: flag for Discussion 7.4
#
# OUTPUTS:
#   outputs/evaluation/transfer_evaluation.csv
#   models/ensemble/transfer_rf_pred.tif
#   models/ensemble/transfer_maxent_pred.tif
#   models/ensemble/transfer_ensemble_pred.tif
#   outputs/figures/main/Fig_transfer_validation.png
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 21: Transfer Validation\n")
cat("========================================\n\n")

set.seed(42)

# ─────────────────────────────────────────────────────────────
# 1. HELPERS
# ─────────────────────────────────────────────────────────────

safe_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  x[length(x)]
}

sm <- function(x, d = 4) {
  v <- safe_scalar(x)
  if (is.na(v)) return(NA_real_)
  round(v, d)
}

calc_auc <- function(ps, pb) {
  ps <- ps[is.finite(ps)]; pb <- pb[is.finite(pb)]
  if (!length(ps) || !length(pb)) return(NA_real_)
  safe_scalar(as.numeric(pROC::auc(
    pROC::roc(c(rep(1, length(ps)), rep(0, length(pb))),
              c(ps, pb), quiet = TRUE))))
}

compute_boyce <- function(fit, obs, n_bins = 101) {
  fit <- fit[is.finite(fit)]; obs <- obs[is.finite(obs)]
  if (length(fit) < 10 || length(obs) < 3) return(NA_real_)
  brk <- seq(min(fit), max(fit), length.out = n_bins + 1)
  mid <- (brk[-1] + brk[-(n_bins + 1)]) / 2
  pe  <- sapply(seq_len(n_bins), function(i) {
    np <- sum(obs >= brk[i] & obs <= brk[i + 1])
    na <- sum(fit >= brk[i] & fit <= brk[i + 1])
    if (na == 0) return(NA_real_)
    (np / length(obs)) / (na / length(fit))
  })
  v <- is.finite(pe) & pe > 0
  if (sum(v) < 3) return(NA_real_)
  safe_scalar(cor(mid[v], pe[v], method = "spearman",
                  use = "complete.obs"))
}

compute_tss <- function(obs, pred,
                        thr = seq(0.01, 0.99, 0.01)) {
  obs  <- as.integer(obs); pred <- as.numeric(pred)
  ok   <- is.finite(obs) & is.finite(pred)
  obs  <- obs[ok]; pred <- pred[ok]
  if (!length(obs)) return(NA_real_)
  v <- sapply(thr, function(t) {
    tp <- sum(obs == 1 & pred >= t)
    tn <- sum(obs == 0 & pred <  t)
    fp <- sum(obs == 0 & pred >= t)
    fn <- sum(obs == 1 & pred <  t)
    s  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    p  <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    s + p - 1
  })
  safe_scalar(max(v, na.rm = TRUE))
}

# ─────────────────────────────────────────────────────────────
# 2. LOAD SHARED DATA
# ─────────────────────────────────────────────────────────────

cat("--- Loading Shared Data ---\n\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))
pred_stack    <- terra::rast(file.path(OUT_PREDICTORS,
                                       "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
final_names   <- readRDS(file.path(OUT_PREDICTORS,
                                   "final_predictor_names.rds"))
raster_levels <- readRDS(file.path(OUT_MOD_IND,
                                   "gam_raster_levels.rds"))
cat_predictors <- c("Geology", "Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

# All thinned pooled sites
sites_sf <- sf::st_read(file.path(OUT_SITES,
                                  "sites_thinned_pooled.gpkg"),
                        quiet = TRUE)
cat(sprintf("  Total thinned sites: %d\n", nrow(sites_sf)))

# Load CV AUC from ensemble for delta comparison
ens_eval <- tryCatch(
  read.csv(file.path(OUT_EVAL, "ensemble_evaluation.csv"),
           stringsAsFactors = FALSE),
  error = function(e) NULL)
cv_auc_ref <- if (!is.null(ens_eval) &&
                  "cv_auc_mean" %in% names(ens_eval))
  ens_eval$cv_auc_mean[1] else 0.7239
cat(sprintf("  Reference CV AUC (ensemble): %.4f\n\n",
            cv_auc_ref))

# ─────────────────────────────────────────────────────────────
# 3. DEFINE TRANSFER ZONE
# ─────────────────────────────────────────────────────────────

cat("--- Defining Transfer Zone ---\n\n")

# OPTION A: Use taluka shapefile if available
# Set TALUKA_SHP to full path of taluka-level boundary file.
# Required attribute: column with taluka names.
# If file does not exist, Option B (extent) is used.
TALUKA_SHP <- file.path(PATH_DISTRICT,
                        "Taluka_Boundaries",
                        "taluka_boundary.shp")

# Target talukas in southern Chandrapur
TRANSFER_TALUKAS <- c("Gondpipri", "Korpana", "Rajura",
                      "gondpipri", "korpana", "rajura",
                      "GONDPIPRI", "KORPANA", "RAJURA")

transfer_zone <- NULL  # will hold SpatVector

if (file.exists(TALUKA_SHP)) {
  cat("  Using taluka shapefile...\n")
  tal_sf <- sf::st_read(TALUKA_SHP, quiet = TRUE)
  tal_sf <- sf::st_transform(tal_sf, crs = 32644)
  
  # Find name column
  name_cols <- names(tal_sf)[
    !names(tal_sf) %in% c("geom", "geometry")]
  matched_col <- NULL
  for (col in name_cols) {
    vals <- as.character(tal_sf[[col]])
    if (any(tolower(vals) %in%
            c("gondpipri", "korpana", "rajura"))) {
      matched_col <- col; break
    }
  }
  
  if (!is.null(matched_col)) {
    sel <- tolower(as.character(tal_sf[[matched_col]])) %in%
      c("gondpipri", "korpana", "rajura")
    if (sum(sel) >= 1) {
      tz_sf <- sf::st_union(tal_sf[sel, ])
      transfer_zone <- terra::vect(sf::st_as_sf(tz_sf))
      area_km2 <- sum(as.numeric(
        sf::st_area(tal_sf[sel, ]))) / 1e6
      cat(sprintf("  Talukas matched: %d  Area: %.0f km²\n",
                  sum(sel), area_km2))
    }
  }
  
  if (is.null(transfer_zone))
    cat("  Taluka match failed — using extent fallback\n")
}

if (is.null(transfer_zone)) {
  # OPTION B: Bottom 16.4% of y-extent = ~3,500 km²
  # (3500 / 21300 = 0.164)
  cat("  Using spatial extent fallback (bottom 16.4% y)\n")
  ext_full <- terra::ext(boundary_vect)
  y_min    <- ext_full[3]
  y_max    <- ext_full[4]
  y_cut    <- y_min + 0.164 * (y_max - y_min)
  tz_ext   <- terra::ext(ext_full[1], ext_full[2],
                         y_min, y_cut)
  tz_rast  <- terra::crop(template_30m, tz_ext)
  tz_rast  <- terra::mask(tz_rast, boundary_vect)
  tz_pts   <- terra::as.polygons(
    !is.na(tz_rast), dissolve = TRUE)
  transfer_zone <- tz_pts[terra::values(tz_pts)[, 1] == 1]
  cat(sprintf("  y_min=%.0f  y_cut=%.0f  y_max=%.0f\n",
              y_min, y_cut, y_max))
}

cat("  Transfer zone defined ✓\n\n")

# ─────────────────────────────────────────────────────────────
# 4. SPLIT SITES INTO TRAINING / TRANSFER
# ─────────────────────────────────────────────────────────────

cat("--- Splitting Sites ---\n\n")

sites_in_tz <- tryCatch(
  as.logical(terra::relate(
    terra::vect(sites_sf), transfer_zone,
    relation = "within")),
  error = function(e)
    as.logical(terra::is.related(
      terra::vect(sites_sf), transfer_zone,
      relation = "intersects")))

# Ensure logical vector, same length as sites
if (length(sites_in_tz) != nrow(sites_sf))
  sites_in_tz <- rep(FALSE, nrow(sites_sf))
sites_in_tz[is.na(sites_in_tz)] <- FALSE

sites_transfer <- sites_sf[sites_in_tz, ]   # held-out
sites_training <- sites_sf[!sites_in_tz, ]  # training

cat(sprintf("  Training sites (north):   %d\n",
            nrow(sites_training)))
cat(sprintf("  Transfer sites (south):   %d\n",
            nrow(sites_transfer)))

if (nrow(sites_transfer) < 5) {
  cat("\n  ⚠ < 5 sites in transfer zone.\n")
  cat("  Check transfer zone definition.\n")
  cat("  Using all sites outside fold 5 as proxy.\n\n")
  
  # Fallback: use CV fold 5 as transfer zone
  cv_design <- readRDS(file.path(OUT_CV,
                                 "cv_block_assignments.rds"))
  sites_cv  <- sf::st_read(file.path(OUT_CV,
                                     "sites_with_folds.gpkg"),
                           quiet = TRUE)
  sites_transfer <- sites_cv[sites_cv$fold_id == 5, ]
  sites_training <- sites_cv[sites_cv$fold_id != 5, ]
  cat(sprintf("  Fallback — training: %d  transfer: %d\n\n",
              nrow(sites_training), nrow(sites_transfer)))
}

# ─────────────────────────────────────────────────────────────
# 5. PREPARE TRAINING DATA
# ─────────────────────────────────────────────────────────────

cat("--- Preparing Training Data ---\n\n")

# Background: exclude points in transfer zone
bg_sf <- sf::st_read(file.path(OUT_BACKGROUND,
                               "background_N10000.gpkg"),
                     quiet = TRUE)

bg_in_tz <- tryCatch(
  as.logical(terra::relate(
    terra::vect(bg_sf), transfer_zone,
    relation = "within")),
  error = function(e) rep(FALSE, nrow(bg_sf)))
bg_in_tz[is.na(bg_in_tz)] <- FALSE

bg_training <- bg_sf[!bg_in_tz, ]
cat(sprintf("  Training background: %d pts (excluded %d)\n\n",
            nrow(bg_training), sum(bg_in_tz)))

# Extract predictor values
sites_vals <- terra::extract(pred_stack,
                             terra::vect(sites_training),
                             ID = FALSE)
bg_vals    <- terra::extract(pred_stack,
                             terra::vect(bg_training),
                             ID = FALSE)
transfer_vals <- terra::extract(pred_stack,
                                terra::vect(sites_transfer),
                                ID = FALSE)

ok_s  <- complete.cases(sites_vals)
ok_b  <- complete.cases(bg_vals)
ok_t  <- complete.cases(transfer_vals)

sites_vals    <- sites_vals[ok_s, ]
bg_vals       <- bg_vals[ok_b, ]
transfer_vals <- transfer_vals[ok_t, ]

N_TR_PRES <- nrow(sites_vals)
N_TR_BG   <- nrow(bg_vals)
N_TF      <- nrow(transfer_vals)
SPW_TR    <- N_TR_BG / N_TR_PRES

cat(sprintf("  Training: %d sites + %d bg  SPW=%.1f\n",
            N_TR_PRES, N_TR_BG, SPW_TR))
cat(sprintf("  Transfer sites complete: %d\n\n", N_TF))

# Categorical factors
for (col in cat_in_stack) {
  lvls <- raster_levels[[col]]
  sites_vals[[col]]    <- factor(as.integer(sites_vals[[col]]),
                                 levels = lvls)
  bg_vals[[col]]       <- factor(as.integer(bg_vals[[col]]),
                                 levels = lvls)
  transfer_vals[[col]] <- factor(as.integer(transfer_vals[[col]]),
                                 levels = lvls)
}

# ─────────────────────────────────────────────────────────────
# 6. TRAIN RF ON NORTHERN SITES ONLY
# ─────────────────────────────────────────────────────────────

cat("--- Training RF (northern sites only) ---\n")
t0 <- proc.time()

n_preds  <- length(final_names)
NTREE    <- 1000L
MTRY     <- floor(sqrt(n_preds))
samp_sz  <- c("0" = N_TR_PRES, "1" = N_TR_PRES)

all_resp_rf <- factor(
  c(rep("1", N_TR_PRES), rep("0", N_TR_BG)),
  levels = c("0", "1"))
all_data_rf <- rbind(sites_vals, bg_vals)

set.seed(42)
rf_transfer <- randomForest::randomForest(
  x           = all_data_rf,
  y           = all_resp_rf,
  ntree       = NTREE,
  mtry        = MTRY,
  sampsize    = samp_sz,
  replace     = TRUE,
  importance  = FALSE,
  keep.forest = TRUE)

cat(sprintf("  Done in %.1f min\n",
            (proc.time() - t0)[3] / 60))
cat(sprintf("  OOB error: %.2f%%\n",
            rf_transfer$err.rate[NTREE, "OOB"] * 100))

# Predict on transfer sites
rf_pred_safe <- function(model, data) {
  d <- data
  for (col in cat_in_stack)
    d[[col]] <- factor(as.integer(d[[col]]),
                       levels = raster_levels[[col]])
  tryCatch(predict(model, d, type = "prob")[, "1"],
           error = function(e) rep(NA_real_, nrow(d)))
}

rf_tr_preds <- rf_pred_safe(rf_transfer, transfer_vals)
rf_tr_bg    <- rf_pred_safe(rf_transfer, bg_vals)

rf_tr_auc   <- calc_auc(rf_tr_preds, rf_tr_bg)
rf_tr_boyce <- compute_boyce(
  c(rf_tr_preds, rf_tr_bg), rf_tr_preds)
rf_tr_tss   <- compute_tss(
  c(rep(1L, N_TF), rep(0L, N_TR_BG)),
  c(rf_tr_preds, rf_tr_bg))

cat(sprintf("  Transfer AUC:   %.4f\n", sm(rf_tr_auc)))
cat(sprintf("  Transfer Boyce: %.4f\n", sm(rf_tr_boyce)))
cat(sprintf("  Transfer TSS:   %.4f\n\n", sm(rf_tr_tss)))

rm(all_data_rf, all_resp_rf); gc(full = TRUE)

# ─────────────────────────────────────────────────────────────
# 7. TRAIN MAXENT ON NORTHERN SITES ONLY
# ─────────────────────────────────────────────────────────────

cat("--- Training MaxEnt (northern sites only) ---\n")
t0 <- proc.time()

cats_mx    <- cat_in_stack[cat_in_stack %in% final_names]
mx_sites   <- sites_vals
mx_bg      <- bg_vals
mx_tr_site <- transfer_vals

# Convert categoricals to integer for maxnet
for (col in cat_in_stack) {
  mx_sites[[col]]   <- as.integer(mx_sites[[col]])
  mx_bg[[col]]      <- as.integer(mx_bg[[col]])
  mx_tr_site[[col]] <- as.integer(mx_tr_site[[col]])
}

all_p_mx <- c(rep(1L, N_TR_PRES), rep(0L, N_TR_BG))
all_d_mx <- rbind(mx_sites, mx_bg)

set.seed(42)
mx_transfer <- tryCatch(
  maxnet::maxnet(
    p = all_p_mx,
    data = all_d_mx,
    f = maxnet::maxnet.formula(all_p_mx, all_d_mx,
                               classes = "lqh"),
    regmult = 0.5,
    categoricals = cats_mx),
  error = function(e) {
    cat(sprintf("  MaxEnt error: %s\n", e$message))
    NULL
  })

cat(sprintf("  Done in %.1f min\n",
            (proc.time() - t0)[3] / 60))

if (!is.null(mx_transfer)) {
  mx_tr_preds <- tryCatch(
    predict(mx_transfer, mx_tr_site, type = "logistic"),
    error = function(e) rep(NA_real_, N_TF))
  mx_tr_bg_preds <- tryCatch(
    predict(mx_transfer, mx_bg, type = "logistic"),
    error = function(e) rep(NA_real_, N_TR_BG))
  
  mx_tr_auc   <- calc_auc(mx_tr_preds, mx_tr_bg_preds)
  mx_tr_boyce <- compute_boyce(
    c(mx_tr_preds, mx_tr_bg_preds), mx_tr_preds)
  mx_tr_tss   <- compute_tss(
    c(rep(1L, N_TF), rep(0L, N_TR_BG)),
    c(mx_tr_preds, mx_tr_bg_preds))
  
  cat(sprintf("  Transfer AUC:   %.4f\n", sm(mx_tr_auc)))
  cat(sprintf("  Transfer Boyce: %.4f\n", sm(mx_tr_boyce)))
  cat(sprintf("  Transfer TSS:   %.4f\n\n", sm(mx_tr_tss)))
} else {
  mx_tr_preds    <- rep(NA_real_, N_TF)
  mx_tr_bg_preds <- rep(NA_real_, N_TR_BG)
  mx_tr_auc      <- NA_real_
  mx_tr_boyce    <- NA_real_
  mx_tr_tss      <- NA_real_
  cat("  MaxEnt skipped\n\n")
}
rm(all_d_mx, all_p_mx); gc(full = TRUE)

# ─────────────────────────────────────────────────────────────
# 8. TRANSFER ENSEMBLE (equal-weight RF + MaxEnt)
# ─────────────────────────────────────────────────────────────

cat("--- Transfer Ensemble ---\n\n")

# Combine non-NA predictions
ens_tr_preds <- rowMeans(
  cbind(rf_tr_preds, mx_tr_preds), na.rm = TRUE)
ens_tr_bg    <- rowMeans(
  cbind(rf_tr_bg, mx_tr_bg_preds), na.rm = TRUE)

ens_tr_auc   <- calc_auc(ens_tr_preds, ens_tr_bg)
ens_tr_boyce <- compute_boyce(
  c(ens_tr_preds, ens_tr_bg), ens_tr_preds)
ens_tr_tss   <- compute_tss(
  c(rep(1L, N_TF), rep(0L, N_TR_BG)),
  c(ens_tr_preds, ens_tr_bg))

delta_auc <- abs(ens_tr_auc - cv_auc_ref)

cat(sprintf("  Transfer AUC:   %.4f\n", sm(ens_tr_auc)))
cat(sprintf("  Transfer Boyce: %.4f\n", sm(ens_tr_boyce)))
cat(sprintf("  Transfer TSS:   %.4f\n", sm(ens_tr_tss)))
cat(sprintf("  CV AUC (ref):   %.4f\n", cv_auc_ref))
cat(sprintf("  Delta AUC:      %.4f  %s\n",
            delta_auc,
            if (!is.na(delta_auc) && delta_auc > 0.10)
              "⚠ > 0.10 — discuss in 7.4"
            else "✓ within acceptable range"))
cat("\n")

# ─────────────────────────────────────────────────────────────
# 9. RASTER PREDICTION — TRANSFER ZONE
# ─────────────────────────────────────────────────────────────

cat("--- Generating Transfer Zone Rasters ---\n\n")

# Clip predictor stack to transfer zone + buffer
tz_buffered <- terra::buffer(transfer_zone, 5000)
pred_tz     <- terra::crop(pred_stack, tz_buffered)
pred_tz     <- terra::mask(pred_tz, boundary_vect)

# RF raster prediction
rf_pred_raster_fn <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  for (col in cat_in_stack[cat_in_stack %in% final_names]) {
    if (col %in% names(data))
      data[[col]] <- factor(as.integer(data[[col]]),
                            levels = raster_levels[[col]])
  }
  ok <- complete.cases(data)
  if (any(ok))
    result[ok] <- tryCatch(
      predict(model, data[ok, , drop = FALSE],
              type = "prob")[, "1"],
      error = function(e) rep(NA_real_, sum(ok)))
  result
}

cat("  RF raster (transfer zone)...\n")
rf_tz_path <- file.path(OUT_MOD_ENS, "transfer_rf_pred.tif")
terra::predict(pred_tz, rf_transfer,
               fun      = rf_pred_raster_fn,
               na.rm    = FALSE,
               filename = rf_tz_path,
               overwrite = TRUE,
               wopt     = list(datatype = "FLT4S"))
cat("  ✓ transfer_rf_pred.tif\n")
gc(full = TRUE)

# MaxEnt raster prediction
if (!is.null(mx_transfer)) {
  mx_pred_raster_fn <- function(model, data, ...) {
    result <- rep(NA_real_, nrow(data))
    for (col in cat_in_stack[cat_in_stack %in% final_names]) {
      if (col %in% names(data))
        data[[col]] <- as.integer(data[[col]])
    }
    ok <- complete.cases(data)
    if (any(ok))
      result[ok] <- tryCatch(
        predict(model, data[ok, , drop = FALSE],
                type = "logistic"),
        error = function(e) rep(NA_real_, sum(ok)))
    result
  }
  
  cat("  MaxEnt raster (transfer zone)...\n")
  mx_tz_path <- file.path(OUT_MOD_ENS,
                          "transfer_maxent_pred.tif")
  terra::predict(pred_tz, mx_transfer,
                 fun      = mx_pred_raster_fn,
                 na.rm    = FALSE,
                 filename = mx_tz_path,
                 overwrite = TRUE,
                 wopt     = list(datatype = "FLT4S"))
  cat("  ✓ transfer_maxent_pred.tif\n")
  gc(full = TRUE)
  
  # Ensemble raster
  cat("  Transfer ensemble raster...\n")
  rf_tz_r  <- terra::rast(rf_tz_path)
  mx_tz_r  <- terra::rast(mx_tz_path)
  ens_tz_r <- terra::app(
    terra::rast(list(rf_tz_r, mx_tz_r)),
    fun = function(x) rowMeans(x, na.rm = TRUE))
  ens_tz_path <- file.path(OUT_MOD_ENS,
                           "transfer_ensemble_pred.tif")
  terra::writeRaster(ens_tz_r, ens_tz_path,
                     overwrite = TRUE, datatype = "FLT4S")
  cat("  ✓ transfer_ensemble_pred.tif\n\n")
} else {
  ens_tz_r    <- terra::rast(rf_tz_path)
  ens_tz_path <- rf_tz_path
  cat("  (ensemble = RF only; MaxEnt failed)\n\n")
}

# ─────────────────────────────────────────────────────────────
# 10. ALSO EVALUATE FULL ENSEMBLE IN TRANSFER ZONE
#     (existing full-data ensemble from Script 17)
# ─────────────────────────────────────────────────────────────

cat("--- Full Ensemble Geographic Evaluation ---\n\n")

ens_primary_path <- file.path(OUT_MOD_ENS,
                              "ensemble_primary.tif")
full_ens_tr_auc   <- NA_real_
full_ens_tr_boyce <- NA_real_

if (file.exists(ens_primary_path)) {
  ens_full_r  <- terra::rast(ens_primary_path)
  full_tr_ext <- terra::extract(ens_full_r,
                                terra::vect(sites_transfer),
                                ID = FALSE)[, 1]
  full_bg_ext <- terra::extract(ens_full_r,
                                terra::vect(bg_sf[!bg_in_tz,]),
                                ID = FALSE)[, 1]
  full_tr_ext <- full_tr_ext[is.finite(full_tr_ext)]
  full_bg_ext <- full_bg_ext[is.finite(full_bg_ext)]
  
  full_ens_tr_auc   <- calc_auc(full_tr_ext, full_bg_ext)
  full_ens_tr_boyce <- compute_boyce(
    c(full_tr_ext, full_bg_ext), full_tr_ext)
  
  cat(sprintf("  Full-ensemble transfer AUC:   %.4f\n",
              sm(full_ens_tr_auc)))
  cat(sprintf("  Full-ensemble transfer Boyce: %.4f\n",
              sm(full_ens_tr_boyce)))
  cat("  NOTE: Full ensemble trained on all sites\n")
  cat("  (optimistic estimate; retrained above = valid)\n\n")
} else {
  cat("  ensemble_primary.tif not found — skipping\n\n")
}

# ─────────────────────────────────────────────────────────────
# 11. SAVE EVALUATION RESULTS
# ─────────────────────────────────────────────────────────────

cat("--- Saving Evaluation ---\n\n")

transfer_eval <- data.frame(
  metric                  = c(
    "n_training_sites",
    "n_transfer_sites",
    "n_training_bg",
    "n_bg_excluded_tz",
    "rf_transfer_auc",
    "rf_transfer_boyce",
    "rf_transfer_tss",
    "maxent_transfer_auc",
    "maxent_transfer_boyce",
    "maxent_transfer_tss",
    "ensemble_transfer_auc",
    "ensemble_transfer_boyce",
    "ensemble_transfer_tss",
    "cv_auc_reference",
    "delta_auc_ensemble",
    "delta_exceeds_0.10",
    "full_ensemble_tz_auc",
    "full_ensemble_tz_boyce"),
  value                   = c(
    N_TR_PRES,
    N_TF,
    N_TR_BG,
    sum(bg_in_tz),
    sm(rf_tr_auc),
    sm(rf_tr_boyce),
    sm(rf_tr_tss),
    sm(mx_tr_auc),
    sm(mx_tr_boyce),
    sm(mx_tr_tss),
    sm(ens_tr_auc),
    sm(ens_tr_boyce),
    sm(ens_tr_tss),
    round(cv_auc_ref, 4),
    sm(delta_auc),
    ifelse(!is.na(delta_auc) & delta_auc > 0.10, 1, 0),
    sm(full_ens_tr_auc),
    sm(full_ens_tr_boyce)),
  stringsAsFactors = FALSE)

write.csv(transfer_eval,
          file.path(OUT_EVAL, "transfer_evaluation.csv"),
          row.names = FALSE)
cat("  ✓ transfer_evaluation.csv\n\n")

# ─────────────────────────────────────────────────────────────
# 12. FIGURE — TRANSFER VALIDATION
# ─────────────────────────────────────────────────────────────

cat("--- Figure: Transfer Validation ---\n\n")

tryCatch({
  png(file.path(OUT_FIG_MAIN, "Fig_transfer_validation.png"),
      width = 7200, height = 2800, res = 300)
  par(mfrow = c(1, 3), mar = c(4, 4, 3.5, 2))
  
  # Panel 1: Transfer ensemble prediction map
  terra::plot(ens_tz_r,
              main  = "Transfer Ensemble\n(RF + MaxEnt, northern training)",
              col   = viridisLite::viridis(100),
              range = c(0, 1),
              axes  = FALSE, legend = TRUE,
              cex.main = 0.85)
  terra::plot(transfer_zone, add = TRUE,
              border = "red", lwd = 1.5)
  if (nrow(sites_transfer) > 0)
    terra::plot(terra::vect(sites_transfer), add = TRUE,
                col = "red", pch = 16, cex = 0.6)
  
  # Panel 2: ROC curve — transfer ensemble
  obs_vec  <- c(rep(1L, N_TF), rep(0L, length(ens_tr_bg)))
  pred_vec <- c(ens_tr_preds, ens_tr_bg)
  ok_roc   <- is.finite(obs_vec) & is.finite(pred_vec)
  
  if (sum(ok_roc) > 10 && !is.na(ens_tr_auc)) {
    roc_obj <- pROC::roc(obs_vec[ok_roc], pred_vec[ok_roc],
                         quiet = TRUE)
    plot(roc_obj, col = "#2A9D8F", lwd = 2.5,
         main = sprintf(
           "Transfer ROC\nAUC=%.4f  \u0394=%.4f",
           ens_tr_auc, delta_auc),
         cex.main = 0.85)
    abline(a = 1, b = -1, col = "grey60", lty = 2)
    legend("bottomright",
           legend = c(sprintf("Transfer (N=%d)", N_TF),
                      "Random"),
           col    = c("#2A9D8F", "grey60"),
           lwd    = c(2.5, 1), lty = c(1, 2),
           bty    = "n", cex = 0.75)
  } else {
    plot.new()
    text(0.5, 0.5, "Insufficient transfer data\nfor ROC",
         cex = 1.2, col = "grey50")
  }
  
  # Panel 3: AUC comparison bar chart
  alg_lbls <- c("RF\n(transfer)", "MaxEnt\n(transfer)",
                "Ensemble\n(transfer)", "CV AUC\n(reference)")
  auc_vals <- c(safe_scalar(rf_tr_auc),
                safe_scalar(mx_tr_auc),
                safe_scalar(ens_tr_auc),
                cv_auc_ref)
  bar_cols <- c("#E63946","#2A9D8F","#264653","#F4A261")
  valid    <- is.finite(auc_vals)
  
  if (any(valid)) {
    bp <- barplot(auc_vals[valid],
                  names.arg = alg_lbls[valid],
                  col       = bar_cols[valid],
                  ylim      = c(0, 1),
                  main      = "AUC Comparison\n(transfer vs reference)",
                  ylab      = "AUC",
                  cex.main  = 0.85, cex.names = 0.75,
                  border    = NA)
    abline(h = cv_auc_ref, col = "#F4A261",
           lty = 2, lwd = 1.5)
    abline(h = cv_auc_ref - 0.10, col = "red",
           lty = 3, lwd = 1)
    text(bp, auc_vals[valid] + 0.02,
         labels = sprintf("%.3f", auc_vals[valid]),
         cex = 0.7, col = "black")
    legend("bottomright",
           legend = c("CV AUC ref", "\u0394=0.10 threshold"),
           col    = c("#F4A261", "red"),
           lty    = c(2, 3), lwd = c(1.5, 1),
           bty    = "n", cex = 0.7)
  } else {
    plot.new()
    text(0.5, 0.5, "No valid AUC values", cex = 1.2)
  }
  
  dev.off()
  cat("  ✓ Fig_transfer_validation.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  cat(sprintf("  ✗ Figure error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 13. SUMMARY
# ─────────────────────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 21 COMPLETE — Transfer Validation\n")
cat("========================================\n\n")

cat(sprintf("Training sites:    %d (northern sector)\n",
            N_TR_PRES))
cat(sprintf("Transfer sites:    %d (southern sector)\n",
            N_TF))
cat(sprintf("Training BG pts:   %d\n\n", N_TR_BG))

cat("TRANSFER RESULTS:\n")
cat(sprintf("  RF:       AUC=%.4f  Boyce=%.4f\n",
            sm(rf_tr_auc), sm(rf_tr_boyce)))
cat(sprintf("  MaxEnt:   AUC=%.4f  Boyce=%.4f\n",
            sm(mx_tr_auc), sm(mx_tr_boyce)))
cat(sprintf("  Ensemble: AUC=%.4f  Boyce=%.4f\n",
            sm(ens_tr_auc), sm(ens_tr_boyce)))
cat(sprintf("  CV ref:   AUC=%.4f\n", cv_auc_ref))
cat(sprintf("  Delta:    %.4f  %s\n\n",
            safe_scalar(delta_auc),
            if (!is.na(delta_auc) && delta_auc > 0.10)
              "⚠ DISCUSS IN 7.4"
            else "✓ WITHIN RANGE"))

if (!is.na(delta_auc) && delta_auc <= 0.10) {
  cat("INTERPRETATION:\n")
  cat("  Delta <= 0.10 → model transfers within basin.\n")
  cat("  Supports applicability within Vidarbha region.\n")
  cat("  (Research Design Section 5.13 criterion met)\n\n")
} else {
  cat("INTERPRETATION:\n")
  cat("  Delta > 0.10 → reduced transfer performance.\n")
  cat("  Discuss predictor specificity in Discussion 7.4.\n")
  cat("  Consider SHAP for transfer zone cells separately.\n\n")
}

cat("MANUSCRIPT NOTE (Research Design 5.13 + 7.4):\n")
cat("  'Transfer validation was conducted within the\n")
cat("   Wainganga-Wardha basin and demonstrates spatial\n")
cat("   transferability at the regional scale;\n")
cat("   generalisability to other Indian Palaeolithic\n")
cat("   regions has not been tested and is not claimed.'\n\n")

cat("Files saved:\n")
cat("  outputs/evaluation/transfer_evaluation.csv\n")
cat("  models/ensemble/transfer_rf_pred.tif\n")
cat("  models/ensemble/transfer_maxent_pred.tif\n")
cat("  models/ensemble/transfer_ensemble_pred.tif\n")
cat("  outputs/figures/main/Fig_transfer_validation.png\n")
cat("\nNext: Script 22 — Evaluation Metrics Table\n")
cat("========================================\n")