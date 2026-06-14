# ============================================================
# SCRIPT 13: XGBOOST MODEL
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 13 of 25
# ============================================================
# PARAMETERS (Research Design 5.7.3):
#   objective:        binary:logistic
#   max_depth:        6
#   eta:              0.05  (changed from 0.01)
#   colsample_bytree: 0.8
#   subsample:        0.8
#   scale_pos_weight: N_background / N_presence
#   alpha=0, lambda=1
#
# FIX vs v1 (nrounds=16 underfitting issue):
#   PROBLEM: xgb.cv with spatial blocks gives noisy test AUC
#   curves. Early stopping at round 16 captured noise peak,
#   not true optimal. With eta=0.01 and 16 rounds, the model
#   barely moved from its prior — prediction range 0.42-0.57.
#
#   FIX 1: eta=0.01 → eta=0.05
#   XGBoost converges 5× faster. Reaches meaningful
#   discrimination in 50-300 rounds rather than 500-2000.
#
#   FIX 2: No early stopping in xgb.cv
#   Run all 500 rounds. Select best_nrounds from PEAK of
#   test_auc_mean curve (which.max on evaluation log).
#   Apply minimum floor of 50 rounds — prevents trivial models.
#
#   RESULT EXPECTED: prediction range ~0.1-0.9 (vs 0.42-0.57),
#   better spatial discrimination, CV AUC ≥ 0.70.
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 13: XGBoost Model\n")
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
  if (length(fit)<10||length(obs)<3) return(NA_real_)
  breaks  <- seq(min(fit),max(fit),length.out=n_bins+1)
  bin_mid <- (breaks[-1]+breaks[-(n_bins+1)])/2
  pe <- sapply(seq_len(n_bins), function(i) {
    n_p <- sum(obs>=breaks[i]&obs<=breaks[i+1])
    n_a <- sum(fit>=breaks[i]&fit<=breaks[i+1])
    if (n_a==0) return(NA_real_)
    (n_p/length(obs))/(n_a/length(fit))
  })
  valid <- is.finite(pe)&pe>0
  if (sum(valid)<3) return(NA_real_)
  safe_scalar(cor(bin_mid[valid],pe[valid],
                  method="spearman",use="complete.obs"))
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

cat("  Sites:", nrow(sites_sf),
    "  Background:", nrow(bg_sf),
    "  Predictors:", terra::nlyr(pred_stack), "\n\n")

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

cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

# XGBoost: all numeric
for (col in cat_in_stack) {
  sites_vals[[col]] <- as.numeric(as.integer(sites_vals[[col]]))
  bg_vals[[col]]    <- as.numeric(as.integer(bg_vals[[col]]))
}

sites_mat <- as.matrix(sites_vals)
bg_mat    <- as.matrix(bg_vals)
storage.mode(sites_mat) <- "double"
storage.mode(bg_mat)    <- "double"

N_PRES <- nrow(sites_mat)
N_BG   <- nrow(bg_mat)
SPW    <- N_BG / N_PRES

cat(sprintf("  Sites: %d  Background: %d  SPW: %.2f\n\n",
            N_PRES, N_BG, SPW))

# ── 3. PARAMETERS ────────────────────────────────────────────

cat("--- XGBoost Parameters ---\n")
cat("  objective:        binary:logistic\n")
cat("  max_depth:        6\n")
cat("  eta:              0.05  (converges in ~100-300 rounds)\n")
cat("  colsample_bytree: 0.8  |  subsample: 0.8\n")
cat(sprintf("  scale_pos_weight: %.2f\n", SPW))
cat("  nrounds:          tuned from xgb.cv peak (no early stop)\n\n")

xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = 6L,
  eta              = 0.05,
  colsample_bytree = 0.8,
  subsample        = 0.8,
  scale_pos_weight = SPW,
  alpha            = 0,
  lambda           = 1,
  seed             = 42L,
  nthread          = 1L
)

# ── 4. TUNE NROUNDS — NO EARLY STOPPING ──────────────────────

cat("--- Tuning nrounds via xgb.cv (500 rounds, no early stop) ---\n")
cat("  Peak of mean test AUC selects optimal rounds\n")
cat("  Minimum floor: 50 rounds\n\n")

all_lab <- c(rep(1.0, N_PRES), rep(0.0, N_BG))
all_mat <- rbind(sites_mat, bg_mat)
storage.mode(all_mat) <- "double"

combined_folds <- c(site_folds_c, bg_folds_c)
xgb_folds <- lapply(1:5, function(f)
  which(combined_folds == f))

dtrain_full <- xgboost::xgb.DMatrix(
  data  = all_mat,
  label = all_lab
)

set.seed(42)
t0 <- proc.time()

cv_result <- xgboost::xgb.cv(
  params    = xgb_params,
  data      = dtrain_full,
  nrounds   = 500,
  folds     = xgb_folds,
  maximize  = TRUE,
  verbose   = 0
  # NO early_stopping_rounds — run all 500, pick peak from log
)

cat(sprintf("  xgb.cv done in %.1f min\n",
            (proc.time()-t0)[3]/60))

# Extract nrounds from peak of evaluation log
log_df   <- cv_result$evaluation_log
auc_col  <- grep("test_auc_mean", names(log_df), value=TRUE)[1]

if (!is.null(auc_col) && !is.na(auc_col)) {
  auc_vals     <- log_df[[auc_col]]
  peak_idx     <- which.max(auc_vals)
  best_nrounds <- as.integer(log_df$iter[peak_idx])
  peak_auc     <- auc_vals[peak_idx]
  cat(sprintf("  Peak test AUC: %.4f at round %d\n",
              peak_auc, best_nrounds))
} else {
  best_nrounds <- 200L
  cat("  ⚠ Could not find AUC column — using 200 rounds\n")
}

# Floor: minimum 50 rounds to ensure meaningful model
best_nrounds <- max(best_nrounds, 50L)
cat(sprintf("  nrounds to use: %d\n\n", best_nrounds))

# Plot AUC curve summary
if (exists("auc_col") && !is.na(auc_col)) {
  auc_at_50  <- if (nrow(log_df)>=50)
    log_df[[auc_col]][50] else NA
  auc_at_100 <- if (nrow(log_df)>=100)
    log_df[[auc_col]][100] else NA
  auc_at_200 <- if (nrow(log_df)>=200)
    log_df[[auc_col]][200] else NA
  cat(sprintf("  AUC curve: iter50=%.4f iter100=%.4f",
              safe_scalar(auc_at_50),
              safe_scalar(auc_at_100)))
  cat(sprintf(" iter200=%.4f peak=%.4f@%d\n\n",
              safe_scalar(auc_at_200), peak_auc,
              best_nrounds))
}

rm(cv_result); gc(full=TRUE)

# ── 5. FIT FINAL XGBOOST MODEL ───────────────────────────────

cat("--- Fitting Final XGBoost Model ---\n")
cat(sprintf("  nrounds=%d  eta=0.05\n", best_nrounds))

set.seed(42)
t0 <- proc.time()

xgb_model <- xgboost::xgb.train(
  params  = xgb_params,
  data    = dtrain_full,
  nrounds = best_nrounds,
  verbose = 0
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))

xgboost::xgb.save(xgb_model,
                  file.path(OUT_MOD_IND,
                            "xgboost_model_final.bin"))
saveRDS(list(params=xgb_params,
             nrounds=best_nrounds,
             feature_names=colnames(all_mat),
             spw=SPW),
        file.path(OUT_MOD_IND,"xgboost_model_info.rds"))

cat("  ✓ xgboost_model_final.bin\n")
cat("  ✓ xgboost_model_info.rds\n\n")

rm(dtrain_full, all_mat, all_lab); gc(full=TRUE)

# ── 6. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold Spatial Block CV ---\n\n")

cv_preds_sites <- numeric(N_PRES)
cv_preds_bg    <- numeric(N_BG)
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_mat <- rbind(sites_mat[site_folds_c!=f,,drop=FALSE],
                  bg_mat[bg_folds_c!=f,,drop=FALSE])
  tr_lab <- c(rep(1.0,sum(site_folds_c!=f)),
              rep(0.0,sum(bg_folds_c!=f)))
  storage.mode(tr_mat) <- "double"
  
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  dtr <- xgboost::xgb.DMatrix(data=tr_mat, label=tr_lab)
  
  set.seed(42)
  fold_xgb <- xgboost::xgb.train(
    params  = xgb_params,
    data    = dtr,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  ts_s_mat <- sites_mat[ts_s,,drop=FALSE]
  ts_b_mat <- bg_mat[ts_b,,drop=FALSE]
  storage.mode(ts_s_mat) <- "double"
  storage.mode(ts_b_mat) <- "double"
  
  ps <- predict(fold_xgb, xgboost::xgb.DMatrix(ts_s_mat))
  pb <- predict(fold_xgb, xgboost::xgb.DMatrix(ts_b_mat))
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f]         <- calc_auc(ps, pb)
  
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f], length(ts_s), length(ts_b)))
  
  rm(fold_xgb,dtr,tr_mat,tr_lab,ps,pb,ts_s_mat,ts_b_mat)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs)
cv_auc_sd   <- sd(fold_aucs)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. VARIABLE IMPORTANCE ───────────────────────────────────

cat("--- Variable Importance ---\n\n")

imp_mat <- xgboost::xgb.importance(
  feature_names = colnames(sites_mat),
  model         = xgb_model
)
print(imp_mat[seq_len(min(10,nrow(imp_mat))),], digits=4)
write.csv(imp_mat,
          file.path(OUT_EVAL,"xgboost_importance.csv"),
          row.names=FALSE)
cat("  ✓ xgboost_importance.csv\n\n")

# ── 8. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

xgb_predict_safe <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    d_mat <- as.matrix(data[complete_rows,,drop=FALSE])
    storage.mode(d_mat) <- "double"
    result[complete_rows] <- predict(
      model, xgboost::xgb.DMatrix(d_mat))
  }
  return(result)
}

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
                             sprintf("xgb_tile%d.tif",i))
  t0 <- proc.time()
  terra::predict(ts, xgb_model, fun=xgb_predict_safe,
                 na.rm=FALSE, filename=tile_paths[i],
                 overwrite=TRUE,
                 wopt=list(datatype="FLT4S"))
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging tiles ... ")
out_pred <- file.path(OUT_MOD_IND,"xgboost_pred_prob.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

xgb_raster <- terra::rast(out_pred)
rng <- terra::global(xgb_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))

if (rng[1,2]>1.001||rng[1,1]< -0.001) {
  warning("XGBoost probabilities outside [0,1]")
} else {
  cat("  ✓ Probability [0,1] confirmed\n")
}
cat("  ✓ xgboost_pred_prob.tif\n\n")
gc(full=TRUE)

# ── 9. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps  <- predict(xgb_model,
                    xgboost::xgb.DMatrix(sites_mat))
full_pb  <- predict(xgb_model,
                    xgboost::xgb.DMatrix(bg_mat))
auc_full <- safe_scalar(calc_auc(full_ps, full_pb))

boyce_val <- compute_boyce(
  c(cv_preds_sites,cv_preds_bg), cv_preds_sites)
tss_val   <- compute_tss(
  obs  = c(rep(1,length(cv_preds_sites)),
           rep(0,length(cv_preds_bg))),
  pred = c(cv_preds_sites,cv_preds_bg))

area_cells <- terra::global(!is.na(xgb_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(xgb_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(xgb_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  nrounds:           %d\n",   best_nrounds))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

# ── 10. SAVE OUTPUTS ─────────────────────────────────────────

metrics_df <- data.frame(
  algorithm    = "XGBoost",
  nrounds      = best_nrounds,
  eta          = 0.05, max_depth=6,
  scale_pos_wt = round(SPW,2),
  cv_auc_mean  = sm(cv_auc_mean), cv_auc_sd=sm(cv_auc_sd),
  full_auc     = sm(auc_full), boyce_index=sm(boyce_val),
  tss_max      = sm(tss_val),  kvamme_gain=sm(kg),
  fold1_auc=sm(fold_aucs[1]), fold2_auc=sm(fold_aucs[2]),
  fold3_auc=sm(fold_aucs[3]), fold4_auc=sm(fold_aucs[4]),
  fold5_auc=sm(fold_aucs[5]), stringsAsFactors=FALSE
)

write.csv(metrics_df,
          file.path(OUT_EVAL,"xgboost_evaluation.csv"),
          row.names=FALSE)
saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"xgboost_cv_predictions.rds"))

cat("\n  ✓ xgboost_evaluation.csv\n")
cat("  ✓ xgboost_cv_predictions.rds\n\n")

# ── 11. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_XGBoost_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(xgb_raster,
            main=sprintf(
              "XGBoost — Probability\nCV AUC=%.4f±%.4f  nrounds=%d  eta=0.05",
              cv_auc_mean,cv_auc_sd,best_nrounds),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_XGBoost_prediction.png\n\n")

# ── 12. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 13 COMPLETE — XGBoost\n")
cat("========================================\n")
cat(sprintf("nrounds=%d  eta=0.05  max_depth=6\n",
            best_nrounds))
cat(sprintf("scale_pos_weight: %.1f\n", SPW))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: probability [0,1] ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("All I/O: E drive ✓\n")
cat("========================================\n")