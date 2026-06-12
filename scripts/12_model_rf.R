# ============================================================
# SCRIPT 12: RANDOM FOREST MODEL
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 12 of 25
# ============================================================
# PARAMETERS (Research Design 5.7.2):
#   ntree  = 1,000
#   mtry   = floor(sqrt(p))
#   Imbalance correction: BALANCED BOOTSTRAP SAMPLING
#     sampsize = c("0"=N_PRES, "1"=N_PRES) per tree
#     Draws equal numbers from each class per bootstrap
#     sample — Balanced Random Forest (Chen et al. 2004).
#     More stable than extreme classwt=52.9:1 which caused:
#       - CV AUC=0.57 (near random)
#       - Full AUC=0.996 (severe overfitting)
#       - Raster mean prediction=0.013 (all near zero)
#       - Negative variable importance scores
#   No SMOTE — balanced sampling only (per Research Design)
#
# FIX v2 vs v1:
#   - classwt=52.9 REMOVED → sampsize balanced bootstrap
#   - rf_predict_safe: factor levels from saved cat_levels
#     (not model$forest$xlevels which is unreliable)
#   - CV code: single prediction per fold (no redundant call)
#   - Wrapper: robust NA handling with complete.cases()
#
# Methods text (5.7.2): "Random Forest was fitted with
#   balanced bootstrap sampling (sampsize = N_presence per
#   class per tree) following Chen et al. (2004), correcting
#   for presence-background imbalance without SMOTE."
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 12: Random Forest Model\n")
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
  safe_scalar(max(tss_vals, na.rm=TRUE))
}

compute_boyce <- function(fit, obs, n_bins = 101) {
  fit <- fit[is.finite(fit)]; obs <- obs[is.finite(obs)]
  if (length(fit)<10 || length(obs)<3) return(NA_real_)
  breaks  <- seq(min(fit), max(fit), length.out=n_bins+1)
  bin_mid <- (breaks[-1]+breaks[-(n_bins+1)])/2
  pe <- sapply(seq_len(n_bins), function(i) {
    n_p <- sum(obs>=breaks[i] & obs<=breaks[i+1])
    n_a <- sum(fit>=breaks[i] & fit<=breaks[i+1])
    if (n_a==0) return(NA_real_)
    (n_p/length(obs))/(n_a/length(fit))
  })
  valid <- is.finite(pe) & pe>0
  if (sum(valid)<3) return(NA_real_)
  safe_scalar(cor(bin_mid[valid], pe[valid],
                  method="spearman", use="complete.obs"))
}

calc_auc <- function(ps, pb) {
  ps <- ps[is.finite(ps)]; pb <- pb[is.finite(pb)]
  if (length(ps)==0||length(pb)==0) return(NA_real_)
  as.numeric(pROC::auc(
    pROC::roc(c(rep(1,length(ps)),rep(0,length(pb))),
              c(ps,pb), quiet=TRUE)))
}

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n")

cv_design  <- readRDS(file.path(OUT_CV,
                                "cv_block_assignments.rds"))
site_folds <- cv_design$site_folds
bg_folds   <- cv_design$bg_folds

sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet=TRUE)
bg_sf    <- sf::st_read(file.path(OUT_CV,
                                  "background_with_folds.gpkg"),
                        quiet=TRUE)

final_names <- readRDS(file.path(OUT_PREDICTORS,
                                 "final_predictor_names.rds"))
pred_stack  <- terra::rast(file.path(OUT_PREDICTORS,
                                     "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))

boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES,"study_area_boundary_utm44n.gpkg"),
  quiet=TRUE))

# Load raster-derived factor levels (from GAM script)
raster_levels <- readRDS(
  file.path(OUT_MOD_IND,"gam_raster_levels.rds"))

cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]
n_predictors   <- terra::nlyr(pred_stack)

cat("  Sites:", nrow(sites_sf),
    "  Background:", nrow(bg_sf),
    "  Predictors:", n_predictors, "\n\n")

# ── 2. EXTRACT PREDICTOR VALUES ──────────────────────────────

cat("--- Extracting Predictor Values ---\n")

sites_vals <- terra::extract(pred_stack,
                             terra::vect(sites_sf), ID=FALSE)
bg_vals    <- terra::extract(pred_stack,
                             terra::vect(bg_sf),    ID=FALSE)

sites_ok     <- complete.cases(sites_vals)
bg_ok        <- complete.cases(bg_vals)
sites_vals   <- sites_vals[sites_ok, ]
bg_vals      <- bg_vals[bg_ok, ]
site_folds_c <- site_folds[sites_ok]
bg_folds_c   <- bg_folds[bg_ok]

N_PRES <- nrow(sites_vals)
N_BG   <- nrow(bg_vals)

# Apply raster-derived factor levels (all 9/23 values)
for (col in cat_in_stack) {
  lvls <- raster_levels[[col]]
  sites_vals[[col]] <- factor(as.integer(sites_vals[[col]]),
                              levels=lvls)
  bg_vals[[col]]    <- factor(as.integer(bg_vals[[col]]),
                              levels=lvls)
}

cat(sprintf("  Sites: %d  Background: %d\n\n",
            N_PRES, N_BG))

# ── 3. RF PARAMETERS ─────────────────────────────────────────

NTREE <- 1000L
MTRY  <- floor(sqrt(n_predictors))

# BALANCED BOOTSTRAP: sample N_PRES from each class per tree
# Replaces classwt=52.9:1 which caused severe overfitting
SAMP_SIZE <- c("0"=N_PRES, "1"=N_PRES)

cat("--- RF Parameters ---\n")
cat("  ntree:", NTREE, "\n")
cat("  mtry:", MTRY, "\n")
cat("  Imbalance correction: balanced bootstrap\n")
cat(sprintf("  sampsize: class 0=%d  class 1=%d\n",
            SAMP_SIZE["0"], SAMP_SIZE["1"]))
cat("  (Chen et al. 2004 — no SMOTE)\n\n")

# ── 4. NA-SAFE PREDICT WRAPPER ───────────────────────────────
# Uses raster_levels (not model$forest$xlevels which varies
# by RF version). Handles NAs from boundary cells.

rf_predict_safe <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  # Apply correct factor levels to categorical columns
  for (col in cat_in_stack) {
    if (col %in% names(data)) {
      data[[col]] <- factor(as.integer(data[[col]]),
                            levels=raster_levels[[col]])
    }
  }
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    result[complete_rows] <- predict(
      model,
      data[complete_rows,,drop=FALSE],
      type="prob"
    )[,"1"]
  }
  return(result)
}

# ── 5. FIT FINAL RF MODEL ────────────────────────────────────

cat("--- Fitting Final RF Model ---\n")

all_resp <- factor(c(rep("1",N_PRES),rep("0",N_BG)),
                   levels=c("0","1"))
all_vals <- rbind(sites_vals, bg_vals)

set.seed(42)
t0 <- proc.time()

rf_model <- randomForest::randomForest(
  x           = all_vals,
  y           = all_resp,
  ntree       = NTREE,
  mtry        = MTRY,
  sampsize    = SAMP_SIZE,
  replace     = TRUE,
  importance  = TRUE,
  keep.forest = TRUE
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))
cat("  OOB error rate:",
    round(rf_model$err.rate[NTREE,"OOB"]*100,2),"%\n")
cat("  OOB class 0 error:",
    round(rf_model$err.rate[NTREE,"0"]*100,2),"%\n")
cat("  OOB class 1 error:",
    round(rf_model$err.rate[NTREE,"1"]*100,2),"%\n")

saveRDS(rf_model,
        file.path(OUT_MOD_IND,"rf_model_final.rds"))
saveRDS(raster_levels,
        file.path(OUT_MOD_IND,"rf_cat_levels.rds"))
cat("  ✓ rf_model_final.rds\n\n")

rm(all_vals, all_resp); gc(full=TRUE)

# ── 6. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold Spatial Block CV ---\n\n")

cv_preds_sites <- numeric(N_PRES)
cv_preds_bg    <- numeric(N_BG)
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_resp <- factor(
    c(rep("1",sum(site_folds_c!=f)),
      rep("0",sum(bg_folds_c!=f))),
    levels=c("0","1"))
  tr_data <- rbind(sites_vals[site_folds_c!=f,],
                   bg_vals[bg_folds_c!=f,])
  
  # Fold-specific balanced sampsize
  n_pres_tr <- sum(site_folds_c!=f)
  samp_f    <- c("0"=n_pres_tr, "1"=n_pres_tr)
  
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_rf <- randomForest::randomForest(
    x           = tr_data,
    y           = tr_resp,
    ntree       = NTREE,
    mtry        = MTRY,
    sampsize    = samp_f,
    replace     = TRUE,
    importance  = FALSE,
    keep.forest = TRUE
  )
  
  ps <- rf_predict_safe(fold_rf, sites_vals[ts_s,])
  pb <- rf_predict_safe(fold_rf, bg_vals[ts_b,])
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f]         <- calc_auc(ps, pb)
  
  cat(sprintf("AUC=%.4f  OOB=%.1f%%  (%d sites/%d bg)\n",
              fold_aucs[f],
              rf_model$err.rate[NTREE,"OOB"]*100,
              length(ts_s), length(ts_b)))
  
  rm(fold_rf, tr_data, tr_resp, ps, pb)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. VARIABLE IMPORTANCE ───────────────────────────────────

cat("--- Variable Importance ---\n\n")

imp    <- randomForest::importance(rf_model,
                                   type=1, scale=TRUE)
imp_df <- data.frame(
  predictor            = rownames(imp),
  MeanDecreaseAccuracy = round(imp[,1],4),
  stringsAsFactors     = FALSE
)
imp_df <- imp_df[order(-imp_df$MeanDecreaseAccuracy),]

cat("  Top predictors (Mean Decrease Accuracy):\n")
for (i in seq_len(min(10,nrow(imp_df)))) {
  cat(sprintf("  %2d. %-22s %7.4f\n",
              i, imp_df$predictor[i],
              imp_df$MeanDecreaseAccuracy[i]))
}
write.csv(imp_df,
          file.path(OUT_EVAL,"rf_variable_importance.csv"),
          row.names=FALSE)
cat("  ✓ rf_variable_importance.csv\n\n")

# ── 8. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

ext_full   <- terra::ext(pred_stack)
y_step     <- (ext_full[4]-ext_full[3])/4
tile_paths <- character(4)

for (i in 1:4) {
  cat(sprintf("  Tile %d/4 ... ", i))
  te <- terra::ext(ext_full[1],ext_full[2],
                   ext_full[3]+(i-1)*y_step,
                   ext_full[3]+i*y_step)
  ts <- terra::crop(pred_stack, te)
  tile_paths[i] <- file.path("E:/R_temp",
                             sprintf("rf_tile%d.tif",i))
  t0 <- proc.time()
  
  terra::predict(ts, rf_model,
                 fun       = rf_predict_safe,
                 na.rm     = FALSE,
                 filename  = tile_paths[i],
                 overwrite = TRUE,
                 wopt      = list(datatype="FLT4S"))
  
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging tiles ... ")
out_pred <- file.path(OUT_MOD_IND,"rf_pred_prob.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

rf_raster <- terra::rast(out_pred)
rng <- terra::global(rf_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))

if (rng[1,2]>1.001||rng[1,1]< -0.001) {
  warning("RF probabilities outside [0,1]")
} else {
  cat("  ✓ Probability [0,1] confirmed\n")
}
cat("  ✓ rf_pred_prob.tif\n\n")
gc(full=TRUE)

# ── 9. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps  <- rf_predict_safe(rf_model, sites_vals)
full_pb  <- rf_predict_safe(rf_model, bg_vals)
auc_full <- safe_scalar(calc_auc(
  full_ps[is.finite(full_ps)],
  full_pb[is.finite(full_pb)]))

cv_s_v <- cv_preds_sites[is.finite(cv_preds_sites)]
cv_b_v <- cv_preds_bg[is.finite(cv_preds_bg)]

boyce_val <- compute_boyce(c(cv_s_v,cv_b_v), cv_s_v)
tss_val   <- compute_tss(
  obs  = c(rep(1,length(cv_s_v)),rep(0,length(cv_b_v))),
  pred = c(cv_s_v,cv_b_v))

area_cells <- terra::global(!is.na(rf_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(rf_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(rf_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  OOB error:         %.2f%%\n",
            rf_model$err.rate[NTREE,"OOB"]*100))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

# ── 10. SAVE OUTPUTS ─────────────────────────────────────────

metrics_df <- data.frame(
  algorithm        = "RandomForest",
  ntree            = NTREE,
  mtry             = MTRY,
  imbalance_method = "balanced_bootstrap",
  sampsize_per_cls = N_PRES,
  cv_auc_mean      = sm(cv_auc_mean),
  cv_auc_sd        = sm(cv_auc_sd),
  full_auc         = sm(auc_full),
  boyce_index      = sm(boyce_val),
  tss_max          = sm(tss_val),
  oob_error        = sm(rf_model$err.rate[NTREE,"OOB"]),
  kvamme_gain      = sm(kg),
  fold1_auc        = sm(fold_aucs[1]),
  fold2_auc        = sm(fold_aucs[2]),
  fold3_auc        = sm(fold_aucs[3]),
  fold4_auc        = sm(fold_aucs[4]),
  fold5_auc        = sm(fold_aucs[5]),
  stringsAsFactors = FALSE
)

write.csv(metrics_df,
          file.path(OUT_EVAL,"rf_evaluation.csv"),
          row.names=FALSE)
saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"rf_cv_predictions.rds"))

cat("\n  ✓ rf_evaluation.csv\n")
cat("  ✓ rf_cv_predictions.rds\n\n")

# ── 11. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_RF_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(rf_raster,
            main=sprintf(
              "Random Forest — Probability\nCV AUC=%.4f±%.4f  ntree=%d  mtry=%d  balanced",
              cv_auc_mean,cv_auc_sd,NTREE,MTRY),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_RF_prediction.png\n\n")

# ── 12. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 12 COMPLETE — Random Forest\n")
cat("========================================\n")
cat(sprintf("ntree=%d  mtry=%d\n",NTREE,MTRY))
cat("Imbalance: balanced bootstrap (Chen et al. 2004)\n")
cat(sprintf("sampsize: %d per class per tree\n",N_PRES))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("OOB error:     %.2f%%\n",
            rf_model$err.rate[NTREE,"OOB"]*100))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: probability [0,1] ✓\n")
cat("Balanced bootstrap ✓\n")
cat("NA-safe tiled prediction ✓\n")
cat("All I/O: E drive ✓\n")
cat("========================================\n")