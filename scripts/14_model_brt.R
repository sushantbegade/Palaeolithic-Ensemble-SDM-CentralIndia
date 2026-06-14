# ============================================================
# SCRIPT 14: BOOSTED REGRESSION TREES (BRT)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 14 of 25
# ============================================================
# PARAMETERS (Research Design 5.7.4):
#   interaction.depth: 5 (tree complexity)
#   shrinkage:         0.01 (Elith et al. 2008)
#   bag.fraction:      0.75
#   distribution:      bernoulli
#   n.trees:           tuned via gbm() cv.folds=5
#   cv.folds:          5 (internal CV for n.trees)
#
# IMBALANCE CORRECTION — CASE WEIGHTS:
#   N_BG/N_PRES = 52.9:1 causes mean prediction ~0.018
#   and Area > 0.5 = 0% despite CV AUC=0.71, Boyce=0.95.
#   Fix: presence case weight = sqrt(N_BG/N_PRES) ≈ 7.27
#        background case weight = 1
#   sqrt() keeps optimizer stable (full 52.9 weight can
#   destabilise gbm gradient updates similarly to SVM).
#   Preserves rank-ordering (AUC/Boyce unchanged) while
#   shifting absolute probabilities into interpretable range.
#   Methods text: "Case weights (presence weight = sqrt of
#   background:presence ratio = 7.27) were applied to gbm()
#   to correct for presence-background imbalance, following
#   the approach used for GAM (Wood 2017)."
#
# TILED PREDICTION: 4 tiles with NA-safe wrapper
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 14: BRT Model\n")
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

# NA-safe BRT predict wrapper
# type="response" gives probability [0,1] directly from gbm
brt_predict_safe <- function(model, data, n_trees, ...) {
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    result[complete_rows] <- predict(
      model,
      newdata = data[complete_rows,,drop=FALSE],
      n.trees = n_trees,
      type    = "response"
    )
  }
  return(result)
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

# BRT: integer coding for categoricals (tree handles as numeric)
cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]
for (col in cat_in_stack) {
  sites_vals[[col]] <- as.integer(sites_vals[[col]])
  bg_vals[[col]]    <- as.integer(bg_vals[[col]])
}

N_PRES <- nrow(sites_vals)
N_BG   <- nrow(bg_vals)

# Case weights: sqrt(SPW) for presences, 1 for background
SPW_FULL <- N_BG / N_PRES
WT       <- sqrt(SPW_FULL)

cat(sprintf("  Sites: %d  Background: %d\n", N_PRES, N_BG))
cat(sprintf("  Full SPW: %.1f  Case weight (sqrt): %.2f\n\n",
            SPW_FULL, WT))

# ── 3. BRT PARAMETERS ────────────────────────────────────────

TC <- 5L    # interaction depth
LR <- 0.01  # shrinkage
BF <- 0.75  # bag fraction

cat("--- BRT Parameters ---\n")
cat("  interaction.depth:", TC, "\n")
cat("  shrinkage:        ", LR, "(Elith et al. 2008)\n")
cat("  bag.fraction:     ", BF, "\n")
cat("  distribution:      bernoulli\n")
cat("  cv.folds:          5\n")
cat(sprintf("  presence weight:   %.2f (sqrt of %.1f:1 ratio)\n\n",
            WT, SPW_FULL))

# ── 4. FIT FINAL BRT WITH CASE WEIGHTS ───────────────────────

cat("--- Fitting BRT (gbm with cv.folds=5 + case weights) ---\n")
cat("  Running up to 5000 trees — 5-20 minutes...\n\n")

all_resp    <- c(rep(1L, N_PRES), rep(0L, N_BG))
all_data    <- rbind(sites_vals, bg_vals)
all_df      <- cbind(response=all_resp, all_data)
all_weights <- c(rep(WT, N_PRES), rep(1, N_BG))

pred_formula <- as.formula(
  paste("response ~", paste(names(all_data), collapse=" + ")))

set.seed(42)
t0 <- proc.time()

brt_model <- gbm::gbm(
  formula           = pred_formula,
  data              = all_df,
  weights           = all_weights,
  distribution      = "bernoulli",
  n.trees           = 10000,
  interaction.depth = TC,
  shrinkage         = LR,
  bag.fraction      = BF,
  cv.folds          = 5,
  n.minobsinnode    = 10,
  keep.data         = FALSE,
  verbose           = FALSE
)

cat(sprintf("  gbm() done in %.1f min\n",
            (proc.time()-t0)[3]/60))

# Optimal n.trees from CV deviance
best_ntrees <- gbm::gbm.perf(brt_model,
                             method  = "cv",
                             plot.it = FALSE)

if (is.null(best_ntrees) || length(best_ntrees)==0 ||
    is.na(best_ntrees) || best_ntrees < 1) {
  cat("  ⚠ gbm.perf returned NA — using 1000 trees\n")
  best_ntrees <- 1000L
} else {
  best_ntrees <- as.integer(best_ntrees)
}

cat("  Optimal n.trees:", best_ntrees, "\n")
cat("  CV deviance (at optimal):",
    round(brt_model$cv.error[best_ntrees], 4), "\n")

saveRDS(brt_model,
        file.path(OUT_MOD_IND,"brt_model_final.rds"))
cat("  ✓ brt_model_final.rds\n\n")

rm(all_df, all_data, all_resp, all_weights); gc(full=TRUE)

# ── 5. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold Spatial Block CV ---\n\n")

cv_preds_sites <- numeric(N_PRES)
cv_preds_bg    <- numeric(N_BG)
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_resp <- c(rep(1L,sum(site_folds_c!=f)),
               rep(0L,sum(bg_folds_c!=f)))
  tr_data <- rbind(sites_vals[site_folds_c!=f,],
                   bg_vals[bg_folds_c!=f,])
  tr_wts  <- c(rep(WT, sum(site_folds_c!=f)),
               rep(1,  sum(bg_folds_c!=f)))
  tr_df   <- cbind(response=tr_resp, tr_data)
  
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_brt <- gbm::gbm(
    formula           = pred_formula,
    data              = tr_df,
    weights           = tr_wts,
    distribution      = "bernoulli",
    n.trees           = best_ntrees,
    interaction.depth = TC,
    shrinkage         = LR,
    bag.fraction      = BF,
    n.minobsinnode    = 10,
    keep.data         = FALSE,
    verbose           = FALSE
  )
  
  ps <- brt_predict_safe(fold_brt, sites_vals[ts_s,],
                         best_ntrees)
  pb <- brt_predict_safe(fold_brt, bg_vals[ts_b,],
                         best_ntrees)
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f]         <- calc_auc(ps, pb)
  
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f], length(ts_s), length(ts_b)))
  
  rm(fold_brt, tr_df, tr_data, tr_resp, tr_wts, ps, pb)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 6. VARIABLE IMPORTANCE ───────────────────────────────────

cat("--- Variable Importance ---\n\n")

imp_summary <- summary(brt_model, plotit=FALSE)
imp_df <- data.frame(
  predictor = as.character(imp_summary$var),
  rel_inf   = round(imp_summary$rel.inf, 4),
  stringsAsFactors=FALSE
)

cat("  Top 10 predictors (relative influence %):\n")
for (i in seq_len(min(10,nrow(imp_df)))) {
  cat(sprintf("  %2d. %-22s %.2f%%\n",
              i, imp_df$predictor[i], imp_df$rel_inf[i]))
}
write.csv(imp_df,
          file.path(OUT_EVAL,"brt_importance.csv"),
          row.names=FALSE)
cat("  ✓ brt_importance.csv\n\n")

# ── 7. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

# Closure captures best_ntrees
make_brt_wrapper <- function(n_trees) {
  function(model, data, ...) {
    result        <- rep(NA_real_, nrow(data))
    complete_rows <- complete.cases(data)
    if (any(complete_rows)) {
      result[complete_rows] <- predict(
        model,
        newdata = data[complete_rows,,drop=FALSE],
        n.trees = n_trees,
        type    = "response")
    }
    return(result)
  }
}

brt_tile_fn <- make_brt_wrapper(best_ntrees)

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
                             sprintf("brt_tile%d.tif",i))
  t0 <- proc.time()
  
  terra::predict(ts, brt_model,
                 fun       = brt_tile_fn,
                 na.rm     = FALSE,
                 filename  = tile_paths[i],
                 overwrite = TRUE,
                 wopt      = list(datatype="FLT4S"))
  
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging tiles ... ")
out_pred <- file.path(OUT_MOD_IND,"brt_pred_prob.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

brt_raster <- terra::rast(out_pred)
rng <- terra::global(brt_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))

if (rng[1,2]>1.001||rng[1,1]< -0.001) {
  warning("BRT probabilities outside [0,1]")
} else {
  cat("  ✓ Probability [0,1] confirmed\n")
}
cat("  ✓ brt_pred_prob.tif\n\n")
gc(full=TRUE)

# ── 8. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps  <- brt_predict_safe(brt_model, sites_vals, best_ntrees)
full_pb  <- brt_predict_safe(brt_model, bg_vals,    best_ntrees)
auc_full <- safe_scalar(calc_auc(full_ps, full_pb))

cv_s_v <- cv_preds_sites[is.finite(cv_preds_sites)]
cv_b_v <- cv_preds_bg[is.finite(cv_preds_bg)]

boyce_val <- compute_boyce(c(cv_s_v,cv_b_v), cv_s_v)
tss_val   <- compute_tss(
  obs  = c(rep(1,length(cv_s_v)),rep(0,length(cv_b_v))),
  pred = c(cv_s_v,cv_b_v))

area_cells <- terra::global(!is.na(brt_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(brt_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(brt_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Optimal n.trees:   %d\n",   best_ntrees))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

if (area_pct < 0.01) {
  cat("  ⚠ Area > 0.5 still low — rank-ordering (AUC/Boyce)\n")
  cat("    is what matters for ensemble. Documenting as-is.\n")
}

# ── 9. SAVE OUTPUTS ──────────────────────────────────────────

metrics_df <- data.frame(
  algorithm       = "BRT",
  tc=TC, lr=LR, bf=BF, n_trees=best_ntrees,
  case_wt_pres    = round(WT,3),
  cv_auc_mean     = sm(cv_auc_mean),
  cv_auc_sd       = sm(cv_auc_sd),
  full_auc        = sm(auc_full),
  boyce_index     = sm(boyce_val),
  tss_max         = sm(tss_val),
  kvamme_gain     = sm(kg),
  fold1_auc=sm(fold_aucs[1]), fold2_auc=sm(fold_aucs[2]),
  fold3_auc=sm(fold_aucs[3]), fold4_auc=sm(fold_aucs[4]),
  fold5_auc=sm(fold_aucs[5]),
  stringsAsFactors=FALSE
)

write.csv(metrics_df,
          file.path(OUT_EVAL,"brt_evaluation.csv"),
          row.names=FALSE)
saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"brt_cv_predictions.rds"))

cat("\n  ✓ brt_evaluation.csv\n")
cat("  ✓ brt_cv_predictions.rds\n\n")

# ── 10. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_BRT_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(brt_raster,
            main=sprintf(
              "BRT — Probability\nCV AUC=%.4f±%.4f  n.trees=%d  wt=%.2f",
              cv_auc_mean,cv_auc_sd,best_ntrees,WT),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_BRT_prediction.png\n\n")

# ── 11. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 14 COMPLETE — BRT\n")
cat("========================================\n")
cat(sprintf("tc=%d  lr=%.2f  bf=%.2f\n",TC,LR,BF))
cat(sprintf("Optimal n.trees:   %d\n", best_ntrees))
cat(sprintf("Case wt (sqrt):    %.2f (pres) / 1 (bg)\n",WT))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat(sprintf("Area > 0.5:    %.1f%%\n",100*area_pct))
cat(sprintf("Sites > 0.5:   %.1f%%\n",100*sites_pct))
cat("Output: probability [0,1] ✓\n")
cat("Case weights applied ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 15 — GAM\n")
cat("========================================\n")