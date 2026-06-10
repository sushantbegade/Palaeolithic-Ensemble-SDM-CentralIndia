# ============================================================
# SCRIPT 15: GENERALISED ADDITIVE MODEL (GAM)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 15 of 25
# ============================================================
# SPECIFICATION (Research Design 5.7.5):
#   family:  binomial(link="logit")
#   method:  REML (preferred for smooth term selection)
#   select:  TRUE (shrinkage penalty — automatic variable
#            selection; uninformative smooths shrunk to zero)
#   k:       5 basis dimensions for all continuous smooths
#            (sufficient for N=188 presence records)
#   Continuous predictors: s(x, k=5)
#   Categorical predictors: s(x, bs="re") random effects
#   Formula built from final_predictor_names.rds
#
# METHODS NOTE (verbatim for manuscript 5.7.5):
#   "GAM model selection used the shrinkage penalty
#    (select=TRUE in mgcv) which drives uninformative smooth
#    terms toward zero, providing automatic variable selection
#    within a single model fit (Wood 2017). method='REML' is
#    preferred for smooth term selection (Wood 2011)."
#
# OUTPUT: predict(type="response") → probability [0,1]
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
cat("SCRIPT 15: GAM Model\n")
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

# NA-safe GAM predict wrapper
gam_predict_safe <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    result[complete_rows] <- predict(
      model,
      newdata = data[complete_rows,,drop=FALSE],
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

cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]
cont_names     <- final_names[!final_names %in% cat_predictors]

cat("  Sites:", nrow(sites_sf),
    "  Background:", nrow(bg_sf), "\n")
cat("  Continuous predictors:", length(cont_names), "\n")
cat("  Categorical predictors:",
    paste(cat_in_stack, collapse=", "), "\n\n")

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

# Categoricals as factors for GAM bs="re"
for (col in cat_in_stack) {
  all_levels <- sort(unique(c(
    as.integer(sites_vals[[col]]),
    as.integer(bg_vals[[col]]))))
  sites_vals[[col]] <- factor(as.integer(sites_vals[[col]]),
                              levels=all_levels)
  bg_vals[[col]]    <- factor(as.integer(bg_vals[[col]]),
                              levels=all_levels)
}

N_PRES <- nrow(sites_vals)
N_BG   <- nrow(bg_vals)
cat(sprintf("  Sites: %d  Background: %d\n\n", N_PRES, N_BG))

# ── 3. BUILD GAM FORMULA ─────────────────────────────────────

cat("--- Building GAM Formula ---\n\n")

# Continuous predictors: s(x, k=5)
# Categorical predictors: s(x, bs="re") random effects
smooth_terms <- c(
  paste0("s(", cont_names, ", k=5)"),
  paste0("s(", cat_in_stack, ", bs='re')")
)

gam_formula <- as.formula(
  paste("presence ~", paste(smooth_terms, collapse=" + ")))

cat("  Formula:\n  presence ~\n")
for (term in smooth_terms) cat("   +", term, "\n")
cat("\n")

# ── 4. FIT FINAL GAM MODEL ───────────────────────────────────

cat("--- Fitting Final GAM ---\n")
cat("  family: binomial(link=logit)\n")
cat("  method: REML  |  select: TRUE\n")
cat("  Training: all", N_PRES+N_BG, "records\n\n")

all_resp <- c(rep(1L, N_PRES), rep(0L, N_BG))
all_data <- rbind(sites_vals, bg_vals)
all_df   <- cbind(presence=all_resp, all_data)

set.seed(42)
t0 <- proc.time()

gam_model <- mgcv::gam(
  formula = gam_formula,
  family  = binomial(link="logit"),
  method  = "REML",
  select  = TRUE,
  data    = all_df
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))

# Summary of smooth terms
gam_sum <- summary(gam_model)
cat("\n  Smooth term EDFs (effective degrees of freedom):\n")
cat("  (EDF ≈ 0 = shrunk to zero by select=TRUE)\n\n")
edf_df <- data.frame(
  term    = rownames(gam_sum$s.table),
  edf     = round(gam_sum$s.table[,"edf"], 3),
  p_value = round(gam_sum$s.table[,"p-value"], 4),
  stringsAsFactors=FALSE
)
for (i in seq_len(nrow(edf_df))) {
  shrunk <- if (edf_df$edf[i] < 0.05) " ← shrunk to zero" else ""
  cat(sprintf("  %-30s EDF=%.3f  p=%.4f%s\n",
              edf_df$term[i], edf_df$edf[i],
              edf_df$p_value[i], shrunk))
}

cat(sprintf("\n  Deviance explained: %.1f%%\n",
            100 * (1 - gam_model$deviance/gam_model$null.deviance)))
cat(sprintf("  GCV score: %.4f\n", gam_model$gcv.ubre))

saveRDS(gam_model,
        file.path(OUT_MOD_IND,"gam_model_final.rds"))
cat("  ✓ gam_model_final.rds\n\n")

rm(all_df, all_data, all_resp); gc(full=TRUE)

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
  tr_df   <- cbind(presence=tr_resp, tr_data)
  
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_gam <- tryCatch(
    mgcv::gam(
      formula = gam_formula,
      family  = binomial(link="logit"),
      method  = "REML",
      select  = TRUE,
      data    = tr_df
    ),
    error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(fold_gam)) {
    fold_aucs[f] <- NA_real_
    next
  }
  
  # Ensure test data factor levels match training
  ts_data_s <- sites_vals[ts_s,]
  ts_data_b <- bg_vals[ts_b,]
  for (col in cat_in_stack) {
    ts_data_s[[col]] <- factor(
      as.integer(ts_data_s[[col]]),
      levels=levels(tr_data[[col]]))
    ts_data_b[[col]] <- factor(
      as.integer(ts_data_b[[col]]),
      levels=levels(tr_data[[col]]))
  }
  
  ps <- gam_predict_safe(fold_gam, ts_data_s)
  pb <- gam_predict_safe(fold_gam, ts_data_b)
  
  # Remove NAs before AUC
  ps_valid <- ps[is.finite(ps)]
  pb_valid <- pb[is.finite(pb)]
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f]         <- calc_auc(ps_valid, pb_valid)
  
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f], length(ts_s), length(ts_b)))
  
  rm(fold_gam, tr_df, tr_data, tr_resp,
     ts_data_s, ts_data_b, ps, pb)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 6. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

# Wrapper: ensure factor levels match training data
# and handle NAs
gam_tile_wrapper <- function(model, data, ...) {
  result <- rep(NA_real_, nrow(data))
  # Apply factor levels for categorical columns
  for (col in cat_in_stack) {
    if (col %in% names(data)) {
      train_levels <- levels(gam_model$model[[col]])
      if (!is.null(train_levels)) {
        data[[col]] <- factor(as.integer(data[[col]]),
                              levels=as.integer(train_levels))
      } else {
        data[[col]] <- factor(as.integer(data[[col]]))
      }
    }
  }
  complete_rows <- complete.cases(data)
  if (any(complete_rows)) {
    result[complete_rows] <- tryCatch(
      predict(model,
              newdata = data[complete_rows,,drop=FALSE],
              type    = "response"),
      error = function(e) rep(NA_real_, sum(complete_rows))
    )
  }
  return(result)
}

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
                             sprintf("gam_tile%d.tif",i))
  t0 <- proc.time()
  
  terra::predict(ts, gam_model,
                 fun       = gam_tile_wrapper,
                 na.rm     = FALSE,
                 filename  = tile_paths[i],
                 overwrite = TRUE,
                 wopt      = list(datatype="FLT4S"))
  
  cat(sprintf("%.1f min\n", (proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging tiles ... ")
out_pred <- file.path(OUT_MOD_IND,"gam_pred_prob.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

gam_raster <- terra::rast(out_pred)
rng <- terra::global(gam_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))
cat("  ✓ gam_pred_prob.tif\n\n")
gc(full=TRUE)

# ── 7. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

all_resp2 <- c(rep(1L,N_PRES), rep(0L,N_BG))
all_data2 <- rbind(sites_vals, bg_vals)
for (col in cat_in_stack) {
  lvls <- levels(gam_model$model[[col]])
  if (!is.null(lvls)) {
    all_data2[[col]] <- factor(as.integer(all_data2[[col]]),
                               levels=as.integer(lvls))
  }
}
full_ps  <- gam_predict_safe(gam_model, sites_vals)
full_pb  <- gam_predict_safe(gam_model, bg_vals)
auc_full <- safe_scalar(calc_auc(
  full_ps[is.finite(full_ps)],
  full_pb[is.finite(full_pb)]))

cv_s_valid <- cv_preds_sites[is.finite(cv_preds_sites)]
cv_b_valid <- cv_preds_bg[is.finite(cv_preds_bg)]

boyce_val <- compute_boyce(
  fit = c(cv_s_valid, cv_b_valid),
  obs = cv_s_valid)

tss_val <- compute_tss(
  obs  = c(rep(1,length(cv_s_valid)),
           rep(0,length(cv_b_valid))),
  pred = c(cv_s_valid, cv_b_valid))

area_cells <- terra::global(!is.na(gam_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(gam_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(gam_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Deviance expl.:    %.1f%%\n",
            100*(1-gam_model$deviance/gam_model$null.deviance)))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

# ── 8. SAVE OUTPUTS ──────────────────────────────────────────

# Predictors shrunk to zero by select=TRUE
shrunk_vars <- edf_df$term[edf_df$edf < 0.05]
cat("\n  Shrunk to zero by select=TRUE:",
    if (length(shrunk_vars)==0) "none" else
      paste(shrunk_vars, collapse=", "), "\n")

metrics_df <- data.frame(
  algorithm       = "GAM",
  family          = "binomial(logit)",
  method          = "REML",
  select          = TRUE,
  deviance_expl   = round(
    100*(1-gam_model$deviance/gam_model$null.deviance),2),
  shrunk_terms    = paste(shrunk_vars, collapse=";"),
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
          file.path(OUT_EVAL,"gam_evaluation.csv"),
          row.names=FALSE)
write.csv(edf_df,
          file.path(OUT_EVAL,"gam_smooth_edfs.csv"),
          row.names=FALSE)
saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"gam_cv_predictions.rds"))

cat("  ✓ gam_evaluation.csv\n")
cat("  ✓ gam_smooth_edfs.csv\n")
cat("  ✓ gam_cv_predictions.rds\n\n")

# ── 9. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_GAM_prediction.png"),
    width=2400, height=2400, res=300)
terra::plot(gam_raster,
            main=sprintf(
              "GAM — Probability\nCV AUC=%.4f±%.4f  Dev.expl=%.1f%%",
              cv_auc_mean,cv_auc_sd,
              100*(1-gam_model$deviance/gam_model$null.deviance)),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_GAM_prediction.png\n\n")

# ── 10. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 15 COMPLETE — GAM\n")
cat("========================================\n")
cat("family=binomial  method=REML  select=TRUE\n")
cat(sprintf("Deviance expl.: %.1f%%\n",
            100*(1-gam_model$deviance/gam_model$null.deviance)))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: probability [0,1] ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 16 — SVM\n")
cat("========================================\n")