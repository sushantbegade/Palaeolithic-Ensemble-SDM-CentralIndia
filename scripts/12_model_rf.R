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
#   ntree   = 1,000
#   mtry    = floor(sqrt(p))
#   classwt = c("0"=1, "1"=N_background/N_presence)
#   No SMOTE — class weights only (per Research Design)
#   Categorical predictors as factors
#   Output: predict(type="prob")[,"1"] → [0,1] probability
#
# FIX v2 — TILED PREDICTION NA HANDLING:
#   RF predict() fails on rows with any NA value.
#   rf_predict_safe() subsets to complete rows, predicts,
#   reconstructs full result vector with NA for incomplete rows.
#   Factor levels saved explicitly at training time and
#   reapplied in every prediction call — avoids xlevels issues.
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

# Save factor levels EXPLICITLY — used in all predict calls
# Levels = union of all values seen in sites AND background
cat_levels <- list()
for (col in cat_in_stack) {
  lvls <- sort(unique(c(
    as.integer(sites_vals[[col]]),
    as.integer(bg_vals[[col]]))))
  cat_levels[[col]] <- lvls
  sites_vals[[col]] <- factor(as.integer(sites_vals[[col]]),
                              levels=lvls)
  bg_vals[[col]]    <- factor(as.integer(bg_vals[[col]]),
                              levels=lvls)
}

cat(sprintf("  Sites: %d  Background: %d\n",
            nrow(sites_vals), nrow(bg_vals)))
cat("  Factor levels saved for:", paste(names(cat_levels),
                                        collapse=", "), "\n\n")

# ── 3. RF PARAMETERS ─────────────────────────────────────────

N_PRES <- nrow(sites_vals)
N_BG   <- nrow(bg_vals)
NTREE  <- 1000
MTRY   <- floor(sqrt(n_predictors))
CLS_WT <- c("0"=1, "1"=N_BG/N_PRES)

cat("--- RF Parameters ---\n")
cat("  ntree:", NTREE, "\n")
cat("  mtry:", MTRY, "\n")
cat(sprintf("  classwt: background=1  presence=%.1f\n",
            CLS_WT["1"]))
cat("  No SMOTE — class weights only\n\n")

# ── 4. SAFE PREDICT FUNCTION ─────────────────────────────────
# Handles NAs in newdata (cells outside study area boundary)
# and ensures categorical columns have correct factor levels.
# Called by terra::predict() during tiled raster prediction.

rf_predict_safe <- function(model, data, ...) {
  # Apply saved factor levels to categorical columns
  for (col in cat_in_stack) {
    if (col %in% names(data)) {
      data[[col]] <- factor(as.integer(data[[col]]),
                            levels = cat_levels[[col]])
    }
  }
  # Output vector — NA for incomplete rows, prediction otherwise
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    result[complete_rows] <- predict(
      model,
      data[complete_rows, , drop=FALSE],
      type = "prob"
    )[, "1"]
  }
  return(result)
}

# ── 5. FIT FINAL RF MODEL ────────────────────────────────────

cat("--- Fitting Final RF Model ---\n")

all_resp <- factor(c(rep("1",N_PRES), rep("0",N_BG)),
                   levels=c("0","1"))
all_vals <- rbind(sites_vals, bg_vals)

set.seed(42)
t0 <- proc.time()

rf_model <- randomForest::randomForest(
  x           = all_vals,
  y           = all_resp,
  ntree       = NTREE,
  mtry        = MTRY,
  classwt     = CLS_WT,
  importance  = TRUE,
  keep.forest = TRUE
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))
cat("  OOB error rate:",
    round(rf_model$err.rate[NTREE,"OOB"]*100, 2), "%\n")

saveRDS(rf_model,
        file.path(OUT_MOD_IND,"rf_model_final.rds"))
# Save factor levels alongside model for downstream scripts
saveRDS(cat_levels,
        file.path(OUT_MOD_IND,"rf_cat_levels.rds"))
cat("  ✓ rf_model_final.rds\n")
cat("  ✓ rf_cat_levels.rds\n\n")

rm(all_vals, all_resp); gc(full=TRUE)

# ── 6. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold CV ---\n\n")

cv_preds_sites <- numeric(nrow(sites_vals))
cv_preds_bg    <- numeric(nrow(bg_vals))
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_resp <- factor(
    c(rep("1",sum(site_folds_c!=f)),
      rep("0",sum(bg_folds_c!=f))),
    levels=c("0","1"))
  tr_data <- rbind(sites_vals[site_folds_c!=f,],
                   bg_vals[bg_folds_c!=f,])
  ts_s    <- which(site_folds_c==f)
  ts_b    <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_rf <- randomForest::randomForest(
    x           = tr_data,
    y           = tr_resp,
    ntree       = NTREE,
    mtry        = MTRY,
    classwt     = CLS_WT,
    importance  = FALSE,
    keep.forest = TRUE
  )
  
  # Apply rf_predict_safe for consistency (handles factors)
  ps <- rf_predict_safe(fold_rf, sites_vals[ts_s,])
  pb <- rf_predict_safe(fold_rf, bg_vals[ts_b,])
  
  # Remove NAs introduced by safe predict
  ps <- ps[is.finite(ps)]; pb <- pb[is.finite(pb)]
  
  cv_preds_sites[ts_s] <- rf_predict_safe(fold_rf,
                                          sites_vals[ts_s,])
  cv_preds_bg[ts_b]    <- rf_predict_safe(fold_rf,
                                          bg_vals[ts_b,])
  fold_aucs[f]         <- calc_auc(ps, pb)
  
  cat(sprintf("AUC=%.4f  OOB=%.1f%%  (%d sites/%d bg)\n",
              fold_aucs[f],
              fold_rf$err.rate[NTREE,"OOB"]*100,
              length(ts_s), length(ts_b)))
  
  rm(fold_rf, tr_data, tr_resp, ps, pb)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs)
cv_auc_sd   <- sd(fold_aucs)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

ext_full   <- terra::ext(pred_stack)
y_step     <- (ext_full[4]-ext_full[3])/4
tile_paths <- character(4)

for (i in 1:4) {
  cat(sprintf("  Tile %d/4 ... ", i))
  te <- terra::ext(ext_full[1], ext_full[2],
                   ext_full[3]+(i-1)*y_step,
                   ext_full[3]+i*y_step)
  ts <- terra::crop(pred_stack, te)
  tile_paths[i] <- file.path("E:/R_temp",
                             sprintf("rf_tile%d.tif",i))
  t0 <- proc.time()
  
  # rf_predict_safe handles NAs and factor levels
  terra::predict(ts, rf_model,
                 fun       = rf_predict_safe,
                 na.rm     = FALSE,  # let fun handle NAs
                 filename  = tile_paths[i],
                 overwrite = TRUE,
                 wopt      = list(datatype="FLT4S"))
  
  cat(sprintf("%.1f min\n", (proc.time()-t0)[3]/60))
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

# ── 8. VARIABLE IMPORTANCE ───────────────────────────────────

cat("--- Variable Importance ---\n\n")

imp    <- randomForest::importance(rf_model, type=1, scale=TRUE)
imp_df <- data.frame(predictor=rownames(imp),
                     MeanDecreaseAccuracy=round(imp[,1],4),
                     stringsAsFactors=FALSE)
imp_df <- imp_df[order(-imp_df$MeanDecreaseAccuracy),]

cat("  Top predictors (Mean Decrease Accuracy):\n")
for (i in seq_len(min(10,nrow(imp_df)))) {
  cat(sprintf("  %2d. %-22s %.4f\n",
              i, imp_df$predictor[i],
              imp_df$MeanDecreaseAccuracy[i]))
}
write.csv(imp_df,
          file.path(OUT_EVAL,"rf_variable_importance.csv"),
          row.names=FALSE)
cat("  ✓ rf_variable_importance.csv\n\n")

# ── 9. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

auc_full  <- safe_scalar(calc_auc(
  rf_predict_safe(rf_model, sites_vals),
  rf_predict_safe(rf_model, bg_vals)))

boyce_val <- compute_boyce(
  fit = c(cv_preds_sites, cv_preds_bg),
  obs = cv_preds_sites)

tss_val <- compute_tss(
  obs  = c(rep(1,length(cv_preds_sites)),
           rep(0,length(cv_preds_bg))),
  pred = c(cv_preds_sites, cv_preds_bg))

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
            cv_auc_mean, cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  OOB error:         %.2f%%\n",
            rf_model$err.rate[NTREE,"OOB"]*100))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Area > 0.5:        %.1f%%\n", 100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n", 100*sites_pct))

# ── 10. SAVE OUTPUTS ─────────────────────────────────────────

metrics_df <- data.frame(
  algorithm   = "RandomForest",
  ntree       = NTREE, mtry = MTRY,
  cv_auc_mean = sm(cv_auc_mean),
  cv_auc_sd   = sm(cv_auc_sd),
  full_auc    = sm(auc_full),
  boyce_index = sm(boyce_val),
  tss_max     = sm(tss_val),
  oob_error   = sm(rf_model$err.rate[NTREE,"OOB"]),
  kvamme_gain = sm(kg),
  fold1_auc   = sm(fold_aucs[1]),
  fold2_auc   = sm(fold_aucs[2]),
  fold3_auc   = sm(fold_aucs[3]),
  fold4_auc   = sm(fold_aucs[4]),
  fold5_auc   = sm(fold_aucs[5]),
  stringsAsFactors=FALSE
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
    width=2400, height=2400, res=300)
terra::plot(rf_raster,
            main=sprintf(
              "Random Forest — Probability\nCV AUC=%.4f±%.4f  ntree=%d  mtry=%d",
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
cat(sprintf("ntree=%d  mtry=%d  classwt 1:%.0f\n",
            NTREE,MTRY,CLS_WT["1"]))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("OOB error:     %.2f%%\n",
            rf_model$err.rate[NTREE,"OOB"]*100))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: probability [0,1] ✓\n")
cat("NA-safe tiled prediction ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 13 — XGBoost\n")
cat("========================================\n")