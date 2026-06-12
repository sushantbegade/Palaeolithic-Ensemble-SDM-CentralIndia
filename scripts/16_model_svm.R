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
#   type:      C-svc
#   Grid:      C ∈ {0.1, 1, 10, 50, 100, 500}
#              sigma ∈ {0.01, 0.05, 0.1, 0.5, 1.0} → 30 combos
#   Tuning:    fold-1 AUC screening → top 5 → full 5-fold CV
#              PREFERENCE: converged models over non-converged
#   class.weights: c(neg=1, pos=sqrt(N_BG/N_PRES))
#              Moderate weighting — corrects imbalance without
#              overwhelming the SVM optimizer (full SPW=52.9
#              causes convergence failures; sqrt ≈ 7.3 stable)
#   Probs:     prob.model=TRUE (Platt scaling, internal)
#   Scaling:   continuous predictors z-standardised
#   Categoricals: one-hot encoded (model.matrix)
#
# CONVERGENCE TRACKING:
#   "maximum number of iterations reached" = non-converged.
#   Models are flagged; converged models preferred in selection.
#   If all top-5 non-converged, select by AUC among them.
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
  ps <- ps[is.finite(ps)]; pb <- pb[is.finite(pb)]
  if (length(ps)==0 || length(pb)==0) return(NA_real_)
  as.numeric(pROC::auc(
    pROC::roc(c(rep(1,length(ps)),rep(0,length(pb))),
              c(ps,pb), quiet=TRUE)))
}

# Fit SVM and track convergence warnings
fit_svm_tracked <- function(x, y, C_val, sig_val, cw) {
  converged  <- TRUE
  warn_msgs  <- character(0)
  
  m <- withCallingHandlers(
    tryCatch(
      kernlab::ksvm(x, y,
                    type          = "C-svc",
                    kernel        = "rbfdot",
                    kpar          = list(sigma=sig_val),
                    C             = C_val,
                    class.weights = cw,
                    prob.model    = TRUE,
                    scaled        = FALSE),
      error=function(e) NULL
    ),
    warning=function(w) {
      if (grepl("maximum number of iterations",
                conditionMessage(w))) {
        converged  <<- FALSE
        warn_msgs  <<- c(warn_msgs, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  list(model=m, converged=converged)
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

raster_levels  <- readRDS(
  file.path(OUT_MOD_IND,"gam_raster_levels.rds"))

cat("  Sites:", nrow(sites_sf),
    "  Background:", nrow(bg_sf), "\n\n")

# ── 2. EXTRACT AND PREPROCESS ────────────────────────────────

cat("--- Extracting and Preprocessing ---\n")

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

# Class weight — sqrt(ratio) for stability
# Full ratio (52.9) causes convergence failures
# sqrt(52.9) ≈ 7.3 corrects imbalance without overwhelming optimizer
SPW_FULL <- N_BG / N_PRES
SPW_SVM  <- sqrt(SPW_FULL)
CW <- c("neg"=1, "pos"=SPW_SVM)

cat(sprintf("  Sites: %d  Background: %d\n", N_PRES, N_BG))
cat(sprintf("  Full SPW=%.1f  SVM class weight=%.2f (sqrt)\n\n",
            SPW_FULL, SPW_SVM))

# Scaling parameters from full dataset
all_raw <- rbind(sites_vals, bg_vals)
scale_params <- setNames(
  lapply(cont_names, function(col) {
    list(mean=mean(all_raw[[col]],na.rm=TRUE),
         sd  =max(sd(all_raw[[col]],na.rm=TRUE),1e-10))
  }), cont_names)

# Preprocess: z-scale + one-hot encode
preprocess_svm <- function(df, sp, cat_cols, rl) {
  df2 <- df
  for (col in names(sp)) {
    df2[[col]] <- (df2[[col]]-sp[[col]]$mean)/sp[[col]]$sd
  }
  for (col in cat_cols) {
    df2[[col]] <- factor(as.integer(df2[[col]]),
                         levels=rl[[col]])
  }
  model.matrix(~.-1, data=df2)
}

sites_mat <- preprocess_svm(sites_vals, scale_params,
                            cat_in_stack, raster_levels)
bg_mat    <- preprocess_svm(bg_vals,    scale_params,
                            cat_in_stack, raster_levels)

feat_names <- colnames(sites_mat)
cat(sprintf("  Feature matrix: %d columns\n",
            length(feat_names)))
cat("  Continuous: z-standardised ✓\n")
cat("  Categorical: one-hot encoded ✓\n\n")

# Save preprocessing params
saveRDS(list(scale_params=scale_params,
             cat_levels=raster_levels,
             cat_cols=cat_in_stack,
             cont_cols=cont_names,
             feature_names=feat_names),
        file.path(OUT_MOD_IND,"svm_preprocess_params.rds"))

# ── 3. GRID SEARCH — FOLD-1 SCREENING ────────────────────────

cat("--- Grid Search: Fold-1 Screening (30 combos) ---\n")
cat("  Convergence tracked — converged models preferred\n\n")

# Refined grid: focus on sigma=0.05-0.5 (productive region)
C_grid     <- c(0.1, 1, 10, 50, 100, 500)
sigma_grid <- c(0.01, 0.05, 0.1, 0.5, 1.0)

# Fold-1 split
tr_s <- which(site_folds_c!=1); tr_b <- which(bg_folds_c!=1)
ts_s <- which(site_folds_c==1); ts_b <- which(bg_folds_c==1)

tr_mat <- rbind(sites_mat[tr_s,], bg_mat[tr_b,])
tr_lab <- factor(c(rep("pos",length(tr_s)),
                   rep("neg",length(tr_b))),
                 levels=c("neg","pos"))

screen_res <- data.frame(C=numeric(),sigma=numeric(),
                         fold1_auc=numeric(),
                         converged=logical(),
                         stringsAsFactors=FALSE)
t0 <- proc.time(); n <- 0

for (C in C_grid) {
  for (sig in sigma_grid) {
    n <- n+1
    cat(sprintf("  [%02d/30] C=%6.1f sigma=%5.3f ... ",
                n, C, sig))
    
    res <- fit_svm_tracked(tr_mat, tr_lab, C, sig, CW)
    conv_sym <- if (res$converged) "✓" else "!"
    
    if (is.null(res$model)) {
      cat("FAILED\n")
      screen_res <- rbind(screen_res,
                          data.frame(C=C,sigma=sig,fold1_auc=NA_real_,
                                     converged=FALSE,stringsAsFactors=FALSE))
    } else {
      ps <- tryCatch(
        kernlab::predict(res$model, sites_mat[ts_s,],
                         type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,length(ts_s)))
      pb <- tryCatch(
        kernlab::predict(res$model, bg_mat[ts_b,],
                         type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,length(ts_b)))
      auc_f1 <- calc_auc(ps, pb)
      cat(sprintf("%s AUC=%.4f\n", conv_sym, auc_f1))
      screen_res <- rbind(screen_res,
                          data.frame(C=C,sigma=sig,fold1_auc=auc_f1,
                                     converged=res$converged,
                                     stringsAsFactors=FALSE))
      rm(res,ps,pb)
    }
    gc(full=TRUE)
  }
}

cat(sprintf("\n  Phase 1: %.1f min\n",
            (proc.time()-t0)[3]/60))

rm(tr_mat,tr_lab); gc(full=TRUE)

# Sort: converged first, then by AUC
screen_res <- screen_res[order(
  !screen_res$converged,
  -screen_res$fold1_auc,
  na.last=TRUE),]

cat("\n  Top 10 (✓=converged, !=not converged):\n")
cat(sprintf("  %-7s %-7s  %-8s  %s\n",
            "C","sigma","AUC","Conv"))
cat("  ", paste(rep("-",32),collapse=""), "\n")
for (i in seq_len(min(10,nrow(screen_res)))) {
  r     <- screen_res[i,]
  auc_s <- if (is.na(r$fold1_auc)) "—" else
    sprintf("%.4f", r$fold1_auc)
  cat(sprintf("  %-7.2f %-7.4f  %-8s  %s\n",
              r$C, r$sigma, auc_s,
              if(r$converged) "✓" else "!"))
}

# ── 4. FULL CV ON TOP 5 (CONVERGED PREFERRED) ───────────────

cat("\n--- Full 5-Fold CV (top 5 candidates) ---\n\n")

top5 <- head(screen_res[!is.na(screen_res$fold1_auc),], 5)

cv_grid <- data.frame(C=numeric(),sigma=numeric(),
                      cv_auc_mean=numeric(),
                      cv_auc_sd=numeric(),
                      n_converged=integer(),
                      stringsAsFactors=FALSE)

for (r in seq_len(nrow(top5))) {
  C_i   <- top5$C[r]; sig_i <- top5$sigma[r]
  cat(sprintf("  C=%6.1f sigma=%.3f :", C_i, sig_i))
  fa <- numeric(5); n_conv <- 0L
  
  for (f in 1:5) {
    tr_f <- rbind(sites_mat[site_folds_c!=f,],
                  bg_mat[bg_folds_c!=f,])
    lb_f <- factor(
      c(rep("pos",sum(site_folds_c!=f)),
        rep("neg",sum(bg_folds_c!=f))),
      levels=c("neg","pos"))
    ts_sf <- sites_mat[site_folds_c==f,]
    ts_bf <- bg_mat[bg_folds_c==f,]
    
    res_f <- fit_svm_tracked(tr_f, lb_f, C_i, sig_i, CW)
    if (res_f$converged) n_conv <- n_conv+1L
    
    if (!is.null(res_f$model)) {
      ps <- tryCatch(
        kernlab::predict(res_f$model,ts_sf,
                         type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_sf)))
      pb <- tryCatch(
        kernlab::predict(res_f$model,ts_bf,
                         type="probabilities")[,"pos"],
        error=function(e) rep(NA_real_,nrow(ts_bf)))
      fa[f] <- calc_auc(ps,pb)
      rm(res_f,ps,pb)
    } else {
      fa[f] <- NA_real_
    }
    rm(tr_f,lb_f,ts_sf,ts_bf); gc(full=TRUE)
  }
  
  m_auc <- mean(fa,na.rm=TRUE); s_auc <- sd(fa,na.rm=TRUE)
  cat(sprintf(" AUC=%.4f±%.4f  (%d/5 conv)\n",
              m_auc, s_auc, n_conv))
  cv_grid <- rbind(cv_grid,
                   data.frame(C=C_i,sigma=sig_i,
                              cv_auc_mean=m_auc,cv_auc_sd=s_auc,
                              n_converged=n_conv,
                              stringsAsFactors=FALSE))
}

# Select best — prefer models where ≥3/5 folds converged
cv_grid <- cv_grid[order(
  -(cv_grid$n_converged >= 3),  # converged first
  -cv_grid$cv_auc_mean),]

best_C   <- cv_grid$C[1]
best_sig <- cv_grid$sigma[1]

cat(sprintf("\n  BEST: C=%.1f  sigma=%.3f",
            best_C, best_sig))
cat(sprintf("  CV AUC=%.4f±%.4f  (%d/5 converged)\n\n",
            cv_grid$cv_auc_mean[1], cv_grid$cv_auc_sd[1],
            cv_grid$n_converged[1]))

write.csv(merge(screen_res, cv_grid,
                by=c("C","sigma"), all.x=TRUE),
          file.path(OUT_TABLES,"TableS_SVM_tuning.csv"),
          row.names=FALSE)
cat("  ✓ TableS_SVM_tuning.csv\n\n")

# ── 5. FINAL SVM MODEL ───────────────────────────────────────

cat("--- Fitting Final SVM Model ---\n")
cat(sprintf("  C=%.1f  sigma=%.3f  class.weights: neg=1 pos=%.2f\n\n",
            best_C, best_sig, SPW_SVM))

all_mat  <- rbind(sites_mat, bg_mat)
all_resp <- factor(c(rep("pos",N_PRES),rep("neg",N_BG)),
                   levels=c("neg","pos"))

set.seed(42)
t0 <- proc.time()

final_res <- fit_svm_tracked(all_mat, all_resp,
                             best_C, best_sig, CW)
svm_model <- final_res$model

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))
cat("  Converged:", final_res$converged, "\n")
cat("  Support vectors:", kernlab::nSV(svm_model), "\n")

saveRDS(svm_model,
        file.path(OUT_MOD_IND,"svm_model_final.rds"))
cat("  ✓ svm_model_final.rds\n\n")

rm(all_mat, all_resp); gc(full=TRUE)

# ── 6. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

cat("--- 5-Fold CV ---\n\n")

cv_preds_sites <- numeric(N_PRES)
cv_preds_bg    <- numeric(N_BG)
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  
  tr_f <- rbind(sites_mat[site_folds_c!=f,],
                bg_mat[bg_folds_c!=f,])
  lb_f <- factor(
    c(rep("pos",sum(site_folds_c!=f)),
      rep("neg",sum(bg_folds_c!=f))),
    levels=c("neg","pos"))
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  res_f <- fit_svm_tracked(tr_f, lb_f, best_C, best_sig, CW)
  
  if (!is.null(res_f$model)) {
    ps <- tryCatch(
      kernlab::predict(res_f$model, sites_mat[ts_s,],
                       type="probabilities")[,"pos"],
      error=function(e) rep(NA_real_,length(ts_s)))
    pb <- tryCatch(
      kernlab::predict(res_f$model, bg_mat[ts_b,],
                       type="probabilities")[,"pos"],
      error=function(e) rep(NA_real_,length(ts_b)))
    cv_preds_sites[ts_s] <- ps
    cv_preds_bg[ts_b]    <- pb
    fold_aucs[f] <- calc_auc(ps,pb)
    conv_sym <- if(res_f$converged)"✓" else "!"
    cat(sprintf("AUC=%.4f %s  (%d sites/%d bg)\n",
                fold_aucs[f], conv_sym,
                length(ts_s), length(ts_b)))
    rm(res_f,ps,pb)
  } else {
    cat("FAILED\n"); fold_aucs[f] <- NA_real_
  }
  rm(tr_f,lb_f); gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

svm_tile_fn <- function(model, data, ...) {
  result        <- rep(NA_real_, nrow(data))
  complete_rows <- complete.cases(data)
  if (!any(complete_rows)) return(result)
  
  d_sub <- data[complete_rows,,drop=FALSE]
  for (col in names(scale_params)) {
    if (col %in% names(d_sub))
      d_sub[[col]] <- (d_sub[[col]] -
                         scale_params[[col]]$mean) /
        scale_params[[col]]$sd
  }
  for (col in cat_in_stack) {
    if (col %in% names(d_sub))
      d_sub[[col]] <- factor(as.integer(d_sub[[col]]),
                             levels=raster_levels[[col]])
  }
  mm <- tryCatch(model.matrix(~.-1, data=d_sub),
                 error=function(e) NULL)
  if (is.null(mm)) return(result)
  
  # Align columns
  missing_c <- setdiff(feat_names, colnames(mm))
  for (mc in missing_c)
    mm <- cbind(mm, matrix(0, nrow=nrow(mm), ncol=1,
                           dimnames=list(NULL,mc)))
  mm <- mm[, feat_names, drop=FALSE]
  
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

# ── 8. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps  <- tryCatch(
  kernlab::predict(svm_model, sites_mat,
                   type="probabilities")[,"pos"],
  error=function(e) rep(NA_real_,N_PRES))
full_pb  <- tryCatch(
  kernlab::predict(svm_model, bg_mat,
                   type="probabilities")[,"pos"],
  error=function(e) rep(NA_real_,N_BG))
auc_full <- safe_scalar(calc_auc(full_ps,full_pb))

boyce_val <- compute_boyce(
  c(cv_preds_sites,cv_preds_bg), cv_preds_sites)
tss_val <- compute_tss(
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

# ── 9. SAVE OUTPUTS ──────────────────────────────────────────

metrics_df <- data.frame(
  algorithm="SVM", kernel="RBF",
  C=best_C, sigma=best_sig,
  class_wt_pos=round(SPW_SVM,3),
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

# ── 10. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_SVM_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(svm_raster,
            main=sprintf(
              "SVM (Platt) — Probability\nCV AUC=%.4f±%.4f  C=%.1f  σ=%.3f  wt=%.1f",
              cv_auc_mean,cv_auc_sd,best_C,best_sig,SPW_SVM),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_SVM_prediction.png\n\n")

# ── 11. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 16 COMPLETE — SVM\n")
cat("========================================\n")
cat(sprintf("C=%.1f  sigma=%.3f\n", best_C, best_sig))
cat(sprintf("class.weights: neg=1  pos=%.2f (sqrt ratio)\n",
            SPW_SVM))
cat(sprintf("Support vectors: %d\n",
            kernlab::nSV(svm_model)))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: Platt-calibrated probability [0,1] ✓\n")
cat("Convergence tracked ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 17 — Ensemble\n")
cat("========================================\n")