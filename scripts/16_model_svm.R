# ============================================================
# SCRIPT 16: SUPPORT VECTOR MACHINE (SVM)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 16 of 25
# ============================================================
# PARAMETERS (Research Design 5.7.6):
#   Kernel:    RBF (Gaussian)
#   type:      C-svc (C-classification)
#   Grid:      C ∈ {0.01, 0.1, 1, 10, 100}
#              sigma ∈ {0.001, 0.01, 0.1, 1} → 20 combinations
#   Tuning:    fold-1 AUC screening → top 5 → full 5-fold CV
#   Probs:     prob.model=TRUE (Platt scaling, internal)
#   Scaling:   all continuous predictors z-standardised
#              (mean=0, SD=1) before fitting — REQUIRED for SVM
#              Scaling parameters saved for prediction
#   Categoricals: one-hot encoded (model.matrix) — SVM needs
#                 all-numeric input
#
# TILED PREDICTION: 4 tiles with NA-safe wrapper.
#   Wrapper applies same z-standardisation and one-hot
#   encoding as training before calling predict().
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 16: SVM Model\n")
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
cont_names     <- final_names[!final_names %in% cat_predictors]

cat("  Sites:", nrow(sites_sf),
    "  Background:", nrow(bg_sf), "\n\n")

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

N_PRES <- nrow(sites_vals); N_BG <- nrow(bg_vals)
cat(sprintf("  Sites: %d  Background: %d\n\n", N_PRES, N_BG))

# ── 3. PREPROCESSING ─────────────────────────────────────────
# (a) Compute scaling params from TRAINING data
# (b) z-standardise continuous predictors (mean=0, SD=1)
# (c) One-hot encode categorical predictors

cat("--- Preprocessing ---\n")

# Load raster-derived factor levels (from GAM script)
raster_levels <- readRDS(
  file.path(OUT_MOD_IND,"gam_raster_levels.rds"))

# Scaling parameters from combined sites+background
all_raw <- rbind(sites_vals, bg_vals)

scale_params <- setNames(
  lapply(cont_names, function(col) {
    list(mean = mean(all_raw[[col]], na.rm=TRUE),
         sd   = max(sd(all_raw[[col]], na.rm=TRUE), 1e-10))
  }),
  cont_names
)

# Function: preprocess a data frame → numeric matrix
preprocess_svm <- function(df, scale_p, cat_cols,
                           rast_lvls) {
  df2 <- df
  # Scale continuous
  for (col in names(scale_p)) {
    df2[[col]] <- (df2[[col]] - scale_p[[col]]$mean) /
      scale_p[[col]]$sd
  }
  # One-hot encode categorical
  for (col in cat_cols) {
    df2[[col]] <- factor(as.integer(df2[[col]]),
                         levels=rast_lvls[[col]])
  }
  # model.matrix drops intercept
  mm <- model.matrix(~ . - 1, data=df2)
  return(mm)
}

# Preprocess all data
sites_mat <- preprocess_svm(sites_vals, scale_params,
                            cat_in_stack, raster_levels)
bg_mat    <- preprocess_svm(bg_vals,    scale_params,
                            cat_in_stack, raster_levels)

cat("  Continuous: z-standardised\n")
cat("  Categorical: one-hot encoded\n")
cat("  Feature matrix columns:", ncol(sites_mat), "\n\n")

# Save preprocessing parameters for prediction
saveRDS(list(scale_params=scale_params,
             cat_levels=raster_levels,
             cat_cols=cat_in_stack,
             cont_cols=cont_names,
             feature_names=colnames(sites_mat)),
        file.path(OUT_MOD_IND,"svm_preprocess_params.rds"))

# ── 4. GRID SEARCH — FOLD-1 SCREENING ────────────────────────

cat("--- Grid Search: Fold-1 Screening (20 combos) ---\n\n")

C_grid     <- c(0.01, 0.1, 1, 10, 100)
sigma_grid <- c(0.001, 0.01, 0.1, 1)

# Fold-1 training/test split
tr_idx_s <- which(site_folds_c != 1)
tr_idx_b <- which(bg_folds_c   != 1)
ts_idx_s <- which(site_folds_c == 1)
ts_idx_b <- which(bg_folds_c   == 1)

tr_mat  <- rbind(sites_mat[tr_idx_s,], bg_mat[tr_idx_b,])
tr_resp <- factor(c(rep("pos",length(tr_idx_s)),
                    rep("neg",length(tr_idx_b))),
                  levels=c("neg","pos"))
ts_mat_s <- sites_mat[ts_idx_s,]
ts_mat_b <- bg_mat[ts_idx_b,]

screen_res <- data.frame(
  C=numeric(), sigma=numeric(), fold1_auc=numeric(),
  stringsAsFactors=FALSE)

n_combo <- 0
t0 <- proc.time()

for (C in C_grid) {
  for (sigma in sigma_grid) {
    n_combo <- n_combo + 1
    cat(sprintf("  [%02d/20] C=%-6.3f sigma=%-6.4f ... ",
                n_combo, C, sigma))
    
    m <- tryCatch(
      kernlab::ksvm(tr_mat, tr_resp,
                    type       = "C-svc",
                    kernel     = "rbfdot",
                    kpar       = list(sigma=sigma),
                    C          = C,
                    prob.model = TRUE,
                    scaled     = FALSE),
      error=function(e) NULL
    )
    
    if (is.null(m)) {
      cat("FAILED\n")
      screen_res <- rbind(screen_res,
                          data.frame(C=C,sigma=sigma,fold1_auc=NA_real_,
                                     stringsAsFactors=FALSE))
    } else {
      ps <- tryCatch(
        kernlab::predict(m, ts_mat_s, type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_mat_s)))
      pb <- tryCatch(
        kernlab::predict(m, ts_mat_b, type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_mat_b)))
      auc_f1 <- calc_auc(ps[is.finite(ps)],
                         pb[is.finite(pb)])
      cat(sprintf("AUC=%.4f\n", auc_f1))
      screen_res <- rbind(screen_res,
                          data.frame(C=C,sigma=sigma,fold1_auc=auc_f1,
                                     stringsAsFactors=FALSE))
      rm(m, ps, pb)
    }
    gc(full=TRUE)
  }
}

cat(sprintf("\n  Phase 1 done in %.1f min\n",
            (proc.time()-t0)[3]/60))
rm(tr_mat, tr_resp, ts_mat_s, ts_mat_b); gc(full=TRUE)

screen_res <- screen_res[
  order(-screen_res$fold1_auc, na.last=TRUE),]
cat("\n  Top 10 by fold-1 AUC:\n")
cat(sprintf("  %-8s %-8s  %s\n","C","sigma","AUC"))
cat("  ", paste(rep("-",25),collapse=""), "\n")
for (i in seq_len(min(10,nrow(screen_res)))) {
  cat(sprintf("  %-8.3f %-8.4f  %.4f\n",
              screen_res$C[i],screen_res$sigma[i],
              screen_res$fold1_auc[i]))
}

# ── 5. FULL CV ON TOP 5 ──────────────────────────────────────

cat("\n--- Full 5-Fold CV (top 5 candidates) ---\n\n")

top5 <- head(screen_res[!is.na(screen_res$fold1_auc),],5)

cv_grid <- data.frame(C=numeric(),sigma=numeric(),
                      cv_auc_mean=numeric(),cv_auc_sd=numeric(),
                      stringsAsFactors=FALSE)

for (row_i in seq_len(nrow(top5))) {
  C_i     <- top5$C[row_i]
  sig_i   <- top5$sigma[row_i]
  cat(sprintf("  C=%-6.3f sigma=%-6.4f :", C_i, sig_i))
  fa <- numeric(5)
  
  for (f in 1:5) {
    tr_mat_f <- rbind(sites_mat[site_folds_c!=f,],
                      bg_mat[bg_folds_c!=f,])
    tr_r_f   <- factor(
      c(rep("pos",sum(site_folds_c!=f)),
        rep("neg",sum(bg_folds_c!=f))),
      levels=c("neg","pos"))
    ts_s_f <- sites_mat[site_folds_c==f,]
    ts_b_f <- bg_mat[bg_folds_c==f,]
    
    m_f <- tryCatch(
      kernlab::ksvm(tr_mat_f, tr_r_f,
                    type="C-svc", kernel="rbfdot",
                    kpar=list(sigma=sig_i), C=C_i,
                    prob.model=TRUE, scaled=FALSE),
      error=function(e) NULL)
    
    if (is.null(m_f)) { fa[f] <- NA_real_ } else {
      ps_f <- tryCatch(
        kernlab::predict(m_f,ts_s_f,type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_s_f)))
      pb_f <- tryCatch(
        kernlab::predict(m_f,ts_b_f,type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_b_f)))
      fa[f] <- calc_auc(ps_f[is.finite(ps_f)],
                        pb_f[is.finite(pb_f)])
      rm(m_f,ps_f,pb_f)
    }
    rm(tr_mat_f,tr_r_f,ts_s_f,ts_b_f); gc(full=TRUE)
  }
  
  m_auc <- mean(fa,na.rm=TRUE); s_auc <- sd(fa,na.rm=TRUE)
  cat(sprintf(" AUC=%.4f±%.4f\n",m_auc,s_auc))
  cv_grid <- rbind(cv_grid,
                   data.frame(C=C_i,sigma=sig_i,
                              cv_auc_mean=m_auc,cv_auc_sd=s_auc,
                              stringsAsFactors=FALSE))
}

cv_grid  <- cv_grid[order(-cv_grid$cv_auc_mean),]
best_C   <- cv_grid$C[1]
best_sig <- cv_grid$sigma[1]

cat(sprintf("\n  BEST: C=%.3f  sigma=%.4f  CV AUC=%.4f±%.4f\n\n",
            best_C, best_sig,
            cv_grid$cv_auc_mean[1], cv_grid$cv_auc_sd[1]))

write.csv(merge(screen_res, cv_grid,
                by=c("C","sigma"), all.x=TRUE),
          file.path(OUT_TABLES,"TableS_SVM_tuning.csv"),
          row.names=FALSE)

# ── 6. FINAL SVM MODEL ───────────────────────────────────────

cat("--- Fitting Final SVM Model ---\n")
cat(sprintf("  C=%.3f  sigma=%.4f  prob.model=TRUE\n\n",
            best_C, best_sig))

all_mat  <- rbind(sites_mat, bg_mat)
all_resp <- factor(c(rep("pos",N_PRES),rep("neg",N_BG)),
                   levels=c("neg","pos"))

set.seed(42)
t0 <- proc.time()

svm_model <- kernlab::ksvm(
  all_mat, all_resp,
  type       = "C-svc",
  kernel     = "rbfdot",
  kpar       = list(sigma=best_sig),
  C          = best_C,
  prob.model = TRUE,
  scaled     = FALSE   # already scaled manually
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))
cat("  Support vectors:", kernlab::nSV(svm_model), "\n")

saveRDS(svm_model,
        file.path(OUT_MOD_IND,"svm_model_final.rds"))
cat("  ✓ svm_model_final.rds\n\n")

rm(all_mat, all_resp); gc(full=TRUE)

# ── 7. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold CV ---\n\n")

cv_preds_sites <- numeric(N_PRES)
cv_preds_bg    <- numeric(N_BG)
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_mat_f <- rbind(sites_mat[site_folds_c!=f,],
                    bg_mat[bg_folds_c!=f,])
  tr_r_f   <- factor(
    c(rep("pos",sum(site_folds_c!=f)),
      rep("neg",sum(bg_folds_c!=f))),
    levels=c("neg","pos"))
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_svm <- kernlab::ksvm(
    tr_mat_f, tr_r_f,
    type="C-svc", kernel="rbfdot",
    kpar=list(sigma=best_sig), C=best_C,
    prob.model=TRUE, scaled=FALSE)
  
  ps <- kernlab::predict(fold_svm,
                         sites_mat[ts_s,],
                         type="probabilities")[,"pos"]
  pb <- kernlab::predict(fold_svm,
                         bg_mat[ts_b,],
                         type="probabilities")[,"pos"]
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f]         <- calc_auc(ps,pb)
  
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f], length(ts_s), length(ts_b)))
  
  rm(fold_svm,tr_mat_f,tr_r_f,ps,pb); gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs)
cv_auc_sd   <- sd(fold_aucs)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 8. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

# Preprocess parameters saved earlier — captured in closure
svm_tile_fn <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (!any(complete_rows)) return(result)
  
  d_sub <- data[complete_rows,,drop=FALSE]
  
  # Apply z-scaling to continuous columns
  for (col in names(scale_params)) {
    if (col %in% names(d_sub)) {
      d_sub[[col]] <- (d_sub[[col]] -
                         scale_params[[col]]$mean) /
        scale_params[[col]]$sd
    }
  }
  
  # One-hot encode categoricals
  for (col in cat_in_stack) {
    if (col %in% names(d_sub)) {
      d_sub[[col]] <- factor(as.integer(d_sub[[col]]),
                             levels=raster_levels[[col]])
    }
  }
  
  # model.matrix — same encoding as training
  mm <- tryCatch(
    model.matrix(~.-1, data=d_sub),
    error=function(e) NULL)
  
  if (is.null(mm) || ncol(mm)==0) return(result)
  
  # Align columns to training matrix
  train_cols <- readRDS(
    file.path(OUT_MOD_IND,"svm_preprocess_params.rds"))$feature_names
  missing_cols <- setdiff(train_cols, colnames(mm))
  for (mc in missing_cols) mm <- cbind(mm, setNames(
    data.frame(rep(0,nrow(mm))), mc))
  mm <- mm[, train_cols, drop=FALSE]
  
  preds <- tryCatch(
    kernlab::predict(model, mm,
                     type="probabilities")[,"pos"],
    error=function(e) rep(NA_real_,nrow(mm)))
  
  result[complete_rows] <- preds
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
                             sprintf("svm_tile%d.tif",i))
  t0 <- proc.time()
  terra::predict(ts, svm_model, fun=svm_tile_fn,
                 na.rm=FALSE, filename=tile_paths[i],
                 overwrite=TRUE,
                 wopt=list(datatype="FLT4S"))
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging tiles ... ")
out_pred <- file.path(OUT_MOD_IND,"svm_pred_prob.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

svm_raster <- terra::rast(out_pred)
rng <- terra::global(svm_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))
cat("  ✓ svm_pred_prob.tif\n\n")
gc(full=TRUE)

# ── 9. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps  <- kernlab::predict(svm_model, sites_mat,
                             type="probabilities")[,"pos"]
full_pb  <- kernlab::predict(svm_model, bg_mat,
                             type="probabilities")[,"pos"]
auc_full <- safe_scalar(calc_auc(full_ps,full_pb))

boyce_val <- compute_boyce(
  c(cv_preds_sites,cv_preds_bg), cv_preds_sites)
tss_val   <- compute_tss(
  obs  = c(rep(1,length(cv_preds_sites)),
           rep(0,length(cv_preds_bg))),
  pred = c(cv_preds_sites,cv_preds_bg))

area_cells <- terra::global(!is.na(svm_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(svm_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(svm_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", sm(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", sm(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", sm(tss_val)))
cat(sprintf("  Kvamme's Gain:     %.4f\n", sm(kg)))
cat(sprintf("  Support vectors:   %d\n",
            kernlab::nSV(svm_model)))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

# ── 10. SAVE OUTPUTS ─────────────────────────────────────────

metrics_df <- data.frame(
  algorithm="SVM", kernel="RBF",
  C=best_C, sigma=best_sig,
  n_sv=kernlab::nSV(svm_model),
  cv_auc_mean=sm(cv_auc_mean), cv_auc_sd=sm(cv_auc_sd),
  full_auc=sm(auc_full), boyce_index=sm(boyce_val),
  tss_max=sm(tss_val), kvamme_gain=sm(kg),
  fold1_auc=sm(fold_aucs[1]), fold2_auc=sm(fold_aucs[2]),
  fold3_auc=sm(fold_aucs[3]), fold4_auc=sm(fold_aucs[4]),
  fold5_auc=sm(fold_aucs[5]), stringsAsFactors=FALSE
)

write.csv(metrics_df,
          file.path(OUT_EVAL,"svm_evaluation.csv"),
          row.names=FALSE)
saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"svm_cv_predictions.rds"))

cat("\n  ✓ svm_evaluation.csv\n")
cat("  ✓ svm_cv_predictions.rds\n\n")

# ── 11. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_SVM_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(svm_raster,
            main=sprintf(
              "SVM — Probability (Platt)\nCV AUC=%.4f±%.4f  C=%.3f  σ=%.4f",
              cv_auc_mean,cv_auc_sd,best_C,best_sig),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_SVM_prediction.png\n\n")

# ── 12. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 16 COMPLETE — SVM\n")
cat("========================================\n")
cat(sprintf("C=%.3f  sigma=%.4f\n", best_C, best_sig))
cat(sprintf("Support vectors: %d\n",
            kernlab::nSV(svm_model)))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: Platt-calibrated probability [0,1] ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 17 — Ensemble\n")
cat("========================================\n")