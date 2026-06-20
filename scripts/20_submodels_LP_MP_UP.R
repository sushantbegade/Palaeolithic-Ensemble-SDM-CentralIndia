# ============================================================
# SCRIPT 20: CULTURAL PERIOD SUB-MODELS — COMPLETE FINAL
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 20 of 25
# ============================================================
# ALL FIXES APPLIED:
#   FIX 1 — GAM INFINITE HANG (PRIMARY FIX):
#     setTimeLimit(elapsed=60) around every GAM fold fit.
#     Throws error if REML optimizer exceeds 60 seconds.
#     tryCatch catches it → NULL → fold skipped gracefully.
#     Root cause: REML + select=TRUE + extreme SPW + tiny
#     fold sizes (fold 5 = 4 sites) → infinite loop.
#   FIX 2 — GAM k reduced: k=3 (was 4/5). With N=70 and
#     11 continuous predictors, k=5 = 44 spline params,
#     too many for smallest training folds (~60 sites).
#     k=3 = 22 max params, stable.
#   FIX 3 — GAM gam.control: maxit=100, mgcv.tol=1e-3,
#     nthreads=1. Caps optimizer iterations hard.
#   FIX 4 — GAM weight: sqrt(SPW) prevents logit overflow.
#     SPW=142 → sqrt=11.9 (stable). Full 142 → NaN preds.
#   FIX 5 — GAM/SVM: suppressWarnings(tryCatch()) only.
#     withCallingHandlers(invokeRestart("muffleWarning"))
#     crashes when muffleWarning restart not available.
#   FIX 6 — Duplicate fgam removed. Single definition.
#   FIX 7 — UP background: bg_folds_orig[bg_ok_UP] correct.
#     bg_folds_base is length 9945 (post-NA removal).
#     bg_ok_UP indexes into full 10000-length vector.
#   FIX 8 — safe_mean/safe_sd: NaN-safe display helpers.
#   FIX 9 — SVM: same timeout pattern as GAM for safety.
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 20: Cultural Period Sub-Models\n")
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

safe_mean <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  mean(v)
}

safe_sd <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  v <- v[is.finite(v)]
  if (length(v) < 2) return(NA_real_)
  sd(v)
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

get_tss_thr <- function(obs, pred,
                        thr = seq(0.01, 0.99, 0.01)) {
  obs  <- as.integer(obs); pred <- as.numeric(pred)
  ok   <- is.finite(obs) & is.finite(pred)
  obs  <- obs[ok]; pred <- pred[ok]
  v <- sapply(thr, function(t) {
    tp <- sum(obs == 1 & pred >= t)
    tn <- sum(obs == 0 & pred <  t)
    fp <- sum(obs == 0 & pred >= t)
    fn <- sum(obs == 1 & pred <  t)
    s  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    p  <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    s + p - 1
  })
  thr[which.max(v)]
}

compute_kg <- function(sp, bp) {
  sp <- sp[is.finite(sp)]; bp <- bp[is.finite(bp)]
  if (!length(sp) || !length(bp)) return(NA_real_)
  thr <- get_tss_thr(c(rep(1, length(sp)),
                       rep(0, length(bp))),
                     c(sp, bp))
  bga <- mean(bp > thr, na.rm = TRUE)
  sa  <- mean(sp > thr, na.rm = TRUE)
  if (sa < 1e-9) return(NA_real_)
  1 - bga / sa
}

assign_folds <- function(pts_sf, block_sf, fold_col) {
  pts_utm <- sf::st_transform(pts_sf,
                              crs = sf::st_crs(block_sf))
  pts_utm <- sf::st_make_valid(pts_utm)
  blk     <- sf::st_make_valid(block_sf)
  joined  <- tryCatch(
    sf::st_join(pts_utm, blk[, fold_col],
                join = sf::st_within, left = TRUE),
    error = function(e)
      sf::st_join(pts_utm, blk[, fold_col],
                  join = sf::st_nearest_feature))
  folds <- as.integer(joined[[fold_col]])
  if (any(is.na(folds))) {
    mode_f <- as.integer(names(which.max(
      table(folds, useNA = "no"))))
    folds[is.na(folds)] <- mode_f
  }
  folds
}

make_xgb_mat <- function(df, fn, cats) {
  d <- df[, fn, drop = FALSE]
  for (col in cats[cats %in% fn])
    d[[col]] <- as.numeric(as.integer(d[[col]]))
  m <- as.matrix(d); storage.mode(m) <- "double"
  colnames(m) <- fn; m
}

make_rf_df <- function(df, fn, cats, rl) {
  d <- df[, fn, drop = FALSE]
  for (col in cats[cats %in% fn])
    d[[col]] <- factor(as.integer(d[[col]]),
                       levels = rl[[col]])
  d
}

PRED_COLORS <- c(
  Elevation          = "E63946", Aspect         = "F4A261",
  TRI                = "2A9D8F", TPI            = "264653",
  Plan_Curvature     = "E9C46A", HAND           = "06D6A0",
  Flow_Accum_log10   = "118AB2", Dist_River     = "073B4C",
  Dist_Palaeochannel = "8338EC", Dist_RawMat    = "FF6B6B",
  NDVI               = "56CFE1", Geology        = "FF9F1C",
  Geomorphology      = "4CC9F0", BIO1           = "B5838D",
  BIO12              = "6D6875", BIO15          = "E2A0FF")
to_col <- function(nm) {
  hex <- PRED_COLORS[nm]
  ifelse(is.na(hex), "#AAAAAA", paste0("#", hex))
}

# ─────────────────────────────────────────────────────────────
# 2. SHARED RESOURCES
# ─────────────────────────────────────────────────────────────

cat("--- Loading Shared Resources ---\n\n")

cv_design     <- readRDS(file.path(OUT_CV,
                                   "cv_block_assignments.rds"))
bg_folds_orig <- cv_design$bg_folds  # length 10000 FULL

block_sf <- sf::st_as_sf(cv_design$cv_blocks_obj$blocks)
fold_col <- grep("fold", names(block_sf),
                 ignore.case = TRUE, value = TRUE)[1]
cat(sprintf("  Blocks: %d  fold col='%s'\n",
            nrow(block_sf), fold_col))

bg_sf <- sf::st_read(file.path(OUT_BACKGROUND,
                               "background_N10000.gpkg"),
                     quiet = TRUE)

pred_stack_base  <- terra::rast(file.path(OUT_PREDICTORS,
                                          "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
final_names_base <- readRDS(file.path(OUT_PREDICTORS,
                                      "final_predictor_names.rds"))
raster_levels    <- readRDS(file.path(OUT_MOD_IND,
                                      "gam_raster_levels.rds"))
cat_predictors   <- c("Geology", "Geomorphology")
cat_base         <- cat_predictors[
  cat_predictors %in% final_names_base]

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

cat("  Extracting base bg predictor values...\n")
bg_raw_base   <- terra::extract(pred_stack_base,
                                terra::vect(bg_sf), ID = FALSE)
bg_ok_base    <- complete.cases(bg_raw_base)
bg_vals_base  <- bg_raw_base[bg_ok_base, ]
bg_folds_base <- bg_folds_orig[bg_ok_base]  # 9945
N_BG_BASE     <- nrow(bg_vals_base)
cat(sprintf("  Background (base stack): %d points\n\n",
            N_BG_BASE))

# ─────────────────────────────────────────────────────────────
# 3. CORE SUB-MODEL FUNCTION
# ─────────────────────────────────────────────────────────────

run_period <- function(period_name,
                       sites_path,
                       pred_stack,
                       final_names,
                       cat_in_stack,
                       rl          = raster_levels,
                       bg_vals_p   = bg_vals_base,
                       bg_folds_pp = bg_folds_base) {
  
  cat(sprintf(
    "\n════════════════════════════════════\nPERIOD: %s\n════════════════════════════════════\n\n",
    period_name))
  
  # ── sites + fold assignment ──────────────────────────────
  sites_sf_p   <- sf::st_read(sites_path, quiet = TRUE)
  sites_sf_p$presence <- 1L
  site_folds_p <- assign_folds(sites_sf_p, block_sf, fold_col)
  cat(sprintf("  Sites: %d  fold dist: %s\n",
              nrow(sites_sf_p),
              paste(table(site_folds_p), collapse = "/")))
  
  sites_raw    <- terra::extract(pred_stack,
                                 terra::vect(sites_sf_p),
                                 ID = FALSE)
  ok_s         <- complete.cases(sites_raw)
  sites_raw    <- sites_raw[ok_s, ]
  site_folds_c <- site_folds_p[ok_s]
  N_PRES       <- nrow(sites_raw)
  N_BG         <- nrow(bg_vals_p)
  SPW          <- N_BG / N_PRES
  
  cat(sprintf("  Complete sites: %d  BG: %d  SPW: %.1f\n\n",
              N_PRES, N_BG, SPW))
  
  # ── data variants ───────────────────────────────────────
  sites_xgb <- make_xgb_mat(sites_raw, final_names,
                            cat_in_stack)
  bg_xgb    <- make_xgb_mat(bg_vals_p, final_names,
                            cat_in_stack)
  sites_rf  <- make_rf_df(sites_raw, final_names,
                          cat_in_stack, rl)
  bg_rf     <- make_rf_df(bg_vals_p, final_names,
                          cat_in_stack, rl)
  
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
    nthread          = 1L,
    seed             = 42L)
  
  rf_pred_safe <- function(model, data) {
    d <- data
    for (col in cat_in_stack[cat_in_stack %in% final_names])
      d[[col]] <- factor(as.integer(d[[col]]),
                         levels = rl[[col]])
    tryCatch(predict(model, d, type = "prob")[, "1"],
             error = function(e) rep(NA_real_, nrow(d)))
  }
  
  NTREE <- 1000L
  MTRY  <- floor(sqrt(length(final_names)))
  
  # ── CV storage ──────────────────────────────────────────
  alg_codes <- c("MX", "RF", "XG", "BR", "GA", "SV")
  alg_names <- c("MaxEnt", "RF", "XGBoost",
                 "BRT", "GAM", "SVM")
  cv_s      <- matrix(NA_real_, N_PRES, 6,
                      dimnames = list(NULL, alg_codes))
  cv_b      <- matrix(NA_real_, N_BG,   6,
                      dimnames = list(NULL, alg_codes))
  fold_aucs <- matrix(NA_real_, 5, 6,
                      dimnames = list(NULL, alg_codes))
  
  cat("  Running 5-fold CV:\n")
  
  # ── [1] MaxEnt ──────────────────────────────────────────
  cat("  [MaxEnt] ")
  t0 <- proc.time()
  cats_mx  <- cat_in_stack[cat_in_stack %in% final_names]
  all_p_mx <- c(rep(1L, N_PRES), rep(0L, N_BG))
  all_d_mx <- rbind(sites_raw, bg_vals_p)
  
  for (f in 1:5) {
    idx_all <- c(site_folds_c, bg_folds_pp)
    tr_p    <- all_p_mx[idx_all != f]
    tr_d    <- all_d_mx[idx_all != f, , drop = FALSE]
    ts_s    <- which(site_folds_c == f)
    ts_b    <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    
    m <- tryCatch(
      maxnet::maxnet(
        p = tr_p, data = tr_d,
        f = maxnet::maxnet.formula(tr_p, tr_d,
                                   classes = "lqh"),
        regmult = 0.5, categoricals = cats_mx),
      error = function(e) NULL)
    
    if (!is.null(m)) {
      ps <- tryCatch(
        predict(m, sites_raw[ts_s, , drop = FALSE],
                type = "logistic"),
        error = function(e) rep(NA_real_, length(ts_s)))
      pb <- tryCatch(
        predict(m, bg_vals_p[ts_b, , drop = FALSE],
                type = "logistic"),
        error = function(e) rep(NA_real_, length(ts_b)))
      cv_s[ts_s, "MX"] <- ps
      cv_b[ts_b, "MX"] <- pb
      fold_aucs[f, "MX"] <- calc_auc(ps, pb)
    }
    rm(m); gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f (%.1fmin)\n",
              safe_mean(fold_aucs[, "MX"]),
              safe_sd(fold_aucs[,   "MX"]),
              (proc.time() - t0)[3] / 60))
  rm(all_d_mx, all_p_mx); gc(full = TRUE)
  
  # ── [2] Random Forest ───────────────────────────────────
  cat("  [RF]     ")
  t0 <- proc.time()
  
  for (f in 1:5) {
    ts_s   <- which(site_folds_c == f)
    ts_b   <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    n_p_tr <- sum(site_folds_c != f)
    n_b_tr <- sum(bg_folds_pp  != f)
    if (n_p_tr < 5 || n_b_tr < 5) next
    
    samp_f  <- c("0" = min(n_p_tr, n_b_tr),
                 "1" = min(n_p_tr, n_b_tr))
    tr_resp <- factor(
      c(rep("1", n_p_tr), rep("0", n_b_tr)),
      levels = c("0", "1"))
    tr_data <- rbind(
      sites_rf[site_folds_c != f, , drop = FALSE],
      bg_rf[bg_folds_pp  != f,    , drop = FALSE])
    
    set.seed(42)
    frf <- randomForest::randomForest(
      x = tr_data, y = tr_resp,
      ntree    = NTREE, mtry = MTRY,
      sampsize = samp_f, replace = TRUE,
      importance = FALSE, keep.forest = TRUE)
    
    ps <- rf_pred_safe(frf, sites_rf[ts_s, , drop = FALSE])
    pb <- rf_pred_safe(frf, bg_rf[ts_b,    , drop = FALSE])
    cv_s[ts_s, "RF"] <- ps
    cv_b[ts_b, "RF"] <- pb
    fold_aucs[f, "RF"] <- calc_auc(ps, pb)
    rm(frf, tr_data, tr_resp); gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f (%.1fmin)\n",
              safe_mean(fold_aucs[, "RF"]),
              safe_sd(fold_aucs[,   "RF"]),
              (proc.time() - t0)[3] / 60))
  
  # ── [3] XGBoost ─────────────────────────────────────────
  cat("  [XGB]    ")
  t0 <- proc.time()
  
  for (f in 1:5) {
    ts_s   <- which(site_folds_c == f)
    ts_b   <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    tr_mat <- rbind(
      sites_xgb[site_folds_c != f, , drop = FALSE],
      bg_xgb[bg_folds_pp  != f,    , drop = FALSE])
    tr_lab <- c(rep(1, sum(site_folds_c != f)),
                rep(0, sum(bg_folds_pp  != f)))
    dtr    <- xgboost::xgb.DMatrix(tr_mat, label = tr_lab)
    set.seed(42)
    fxgb <- xgboost::xgb.train(params = xgb_params,
                               data   = dtr,
                               nrounds = 100L,
                               verbose = 0)
    ps <- predict(fxgb, xgboost::xgb.DMatrix(
      sites_xgb[ts_s, , drop = FALSE]))
    pb <- predict(fxgb, xgboost::xgb.DMatrix(
      bg_xgb[ts_b,    , drop = FALSE]))
    cv_s[ts_s, "XG"] <- ps
    cv_b[ts_b, "XG"] <- pb
    fold_aucs[f, "XG"] <- calc_auc(ps, pb)
    rm(fxgb, dtr, tr_mat, tr_lab); gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f (%.1fmin)\n",
              safe_mean(fold_aucs[, "XG"]),
              safe_sd(fold_aucs[,   "XG"]),
              (proc.time() - t0)[3] / 60))
  
  # ── [4] BRT ─────────────────────────────────────────────
  cat("  [BRT]    ")
  t0 <- proc.time()
  WT_brt <- sqrt(SPW)
  
  pred_formula_brt <- as.formula(
    paste("presence ~",
          paste(final_names, collapse = " + ")))
  
  # One initial fit to get optimal n.trees
  all_r  <- c(rep(1L, N_PRES), rep(0L, N_BG))
  all_d  <- rbind(sites_raw, bg_vals_p)
  all_df <- data.frame(presence = all_r, all_d,
                       check.names = FALSE)
  all_w  <- c(rep(WT_brt, N_PRES), rep(1, N_BG))
  
  set.seed(42)
  brt_init <- tryCatch(
    gbm::gbm(pred_formula_brt,
             data = all_df, weights = all_w,
             distribution = "bernoulli",
             n.trees = 5000L, interaction.depth = 5L,
             shrinkage = 0.01, bag.fraction = 0.75,
             cv.folds = 5L, n.minobsinnode = 5L,
             keep.data = FALSE, verbose = FALSE),
    error = function(e) NULL)
  
  best_nt <- if (!is.null(brt_init)) {
    nt <- tryCatch(
      as.integer(gbm::gbm.perf(brt_init, method = "cv",
                               plot.it = FALSE)),
      error = function(e) NA_integer_)
    if (is.null(nt) || is.na(nt) || nt < 100L) 500L
    else nt
  } else 500L
  
  rm(brt_init, all_df, all_d, all_r, all_w); gc(full = TRUE)
  
  for (f in 1:5) {
    ts_s   <- which(site_folds_c == f)
    ts_b   <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    n_p_tr <- sum(site_folds_c != f)
    n_b_tr <- sum(bg_folds_pp  != f)
    tr_r   <- c(rep(1L, n_p_tr), rep(0L, n_b_tr))
    tr_d   <- rbind(
      sites_raw[site_folds_c != f, , drop = FALSE],
      bg_vals_p[bg_folds_pp  != f, , drop = FALSE])
    tr_w   <- c(rep(WT_brt, n_p_tr), rep(1, n_b_tr))
    tr_df  <- data.frame(presence = tr_r, tr_d,
                         check.names = FALSE)
    
    set.seed(42)
    fbrt <- tryCatch(
      gbm::gbm(pred_formula_brt,
               data = tr_df, weights = tr_w,
               distribution = "bernoulli",
               n.trees = best_nt, interaction.depth = 5L,
               shrinkage = 0.01, bag.fraction = 0.75,
               n.minobsinnode = 5L,
               keep.data = FALSE, verbose = FALSE),
      error = function(e) NULL)
    
    if (!is.null(fbrt)) {
      ps <- tryCatch(
        predict(fbrt, sites_raw[ts_s, , drop = FALSE],
                n.trees = best_nt, type = "response"),
        error = function(e) rep(NA_real_, length(ts_s)))
      pb <- tryCatch(
        predict(fbrt, bg_vals_p[ts_b, , drop = FALSE],
                n.trees = best_nt, type = "response"),
        error = function(e) rep(NA_real_, length(ts_b)))
      cv_s[ts_s, "BR"] <- ps
      cv_b[ts_b, "BR"] <- pb
      fold_aucs[f, "BR"] <- calc_auc(ps, pb)
      rm(fbrt)
    }
    gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f  n.trees=%d (%.1fmin)\n",
              safe_mean(fold_aucs[, "BR"]),
              safe_sd(fold_aucs[,   "BR"]),
              best_nt,
              (proc.time() - t0)[3] / 60))
  
  # ── [5] GAM — FIX 1: setTimeLimit prevents infinite hang
  #             FIX 2: k=3 (safe for small N)
  #             FIX 3: gam.control caps optimizer
  #             FIX 4: sqrt(SPW) weight
  #             FIX 5: suppressWarnings(tryCatch()) only
  #             FIX 6: single fgam definition
  cat("  [GAM]    ")
  t0 <- proc.time()
  
  WT_gam <- sqrt(SPW)           # FIX 4: prevent overflow
  
  # FIX 2: k=3 safe for N~60-70 sites per training fold
  # FIX 3: hard cap on REML iterations
  k_val    <- 3L
  gam_ctrl <- mgcv::gam.control(
    maxit     = 100,    # max REML outer iterations
    mgcv.tol  = 1e-3,   # looser tolerance → faster convergence
    nthreads  = 1L)     # single thread → stable
  
  cont_p <- final_names[!final_names %in% cat_in_stack]
  cats_p <- cat_in_stack[cat_in_stack %in% final_names]
  
  # Parametric categorical terms (not bs="re") — more
  # robust with small N and unequal fold sizes
  gam_terms     <- c(paste0("s(", cont_p,
                            ", k=", k_val, ")"),
                     cats_p)
  gam_formula_p <- as.formula(paste("presence ~",
                                    paste(gam_terms, collapse = " + ")))
  
  for (f in 1:5) {
    ts_s   <- which(site_folds_c == f)
    ts_b   <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    n_p_tr <- sum(site_folds_c != f)
    n_b_tr <- sum(bg_folds_pp  != f)
    if (n_p_tr < 5 || n_b_tr < 5) next
    
    tr_r  <- c(rep(1L, n_p_tr), rep(0L, n_b_tr))
    tr_d  <- rbind(
      sites_rf[site_folds_c != f, , drop = FALSE],
      bg_rf[bg_folds_pp  != f,    , drop = FALSE])
    tr_w  <- c(rep(WT_gam, n_p_tr), rep(1, n_b_tr))
    tr_df <- data.frame(presence = as.numeric(tr_r),
                        tr_d, check.names = FALSE)
    
    # FIX 1: setTimeLimit(60) — throws error if GAM
    # optimizer exceeds 60 seconds per fold.
    # tryCatch catches timeout → NULL → fold skipped.
    # Reset with Inf after each fit (CRITICAL).
    set.seed(42)
    setTimeLimit(elapsed = 60, transient = TRUE)
    fgam <- tryCatch(
      suppressWarnings(
        mgcv::gam(gam_formula_p,
                  family  = binomial(link = "logit"),
                  method  = "REML",
                  select  = TRUE,
                  weights = tr_w,
                  data    = tr_df,
                  control = gam_ctrl)),
      error = function(e) NULL)
    setTimeLimit(elapsed = Inf, transient = TRUE)
    
    if (!is.null(fgam)) {
      ts_sd <- sites_rf[ts_s, , drop = FALSE]
      ts_bd <- bg_rf[ts_b,    , drop = FALSE]
      # Ensure factor levels align with training data
      for (col in cats_p) {
        ts_sd[[col]] <- factor(as.integer(ts_sd[[col]]),
                               levels = rl[[col]])
        ts_bd[[col]] <- factor(as.integer(ts_bd[[col]]),
                               levels = rl[[col]])
      }
      ps <- tryCatch(
        suppressWarnings(
          predict(fgam, ts_sd, type = "response")),
        error = function(e) rep(NA_real_, nrow(ts_sd)))
      pb <- tryCatch(
        suppressWarnings(
          predict(fgam, ts_bd, type = "response")),
        error = function(e) rep(NA_real_, nrow(ts_bd)))
      ps[!is.finite(ps)] <- NA_real_
      pb[!is.finite(pb)] <- NA_real_
      cv_s[ts_s, "GA"] <- ps
      cv_b[ts_b, "GA"] <- pb
      fold_aucs[f, "GA"] <- calc_auc(ps, pb)
      rm(fgam)
    }
    gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f (%.1fmin)\n",
              safe_mean(fold_aucs[, "GA"]),
              safe_sd(fold_aucs[,   "GA"]),
              (proc.time() - t0)[3] / 60))
  
  # ── [6] SVM — FIX 5/9: suppressWarnings(tryCatch())
  #              single definition + 90s timeout
  cat("  [SVM]    ")
  t0 <- proc.time()
  WT_svm <- sqrt(SPW)
  CW_svm <- c("neg" = 1, "pos" = WT_svm)
  
  # Scaling params from full combined data
  cont_cols <- final_names[!final_names %in% cat_in_stack]
  all_cont  <- rbind(sites_xgb, bg_xgb)[, cont_cols,
                                        drop = FALSE]
  scale_p <- setNames(
    lapply(cont_cols, function(col)
      list(mean = mean(all_cont[, col], na.rm = TRUE),
           sd   = max(sd(all_cont[, col], na.rm = TRUE),
                      1e-10))),
    cont_cols)
  rm(all_cont)
  
  svm_prep <- function(mat_num) {
    d <- as.data.frame(mat_num, stringsAsFactors = FALSE)
    for (col in cont_cols)
      d[[col]] <- (d[[col]] - scale_p[[col]]$mean) /
        scale_p[[col]]$sd
    cats_svm <- cat_in_stack[cat_in_stack %in% final_names]
    for (col in cats_svm)
      d[[col]] <- factor(as.integer(d[[col]]),
                         levels = rl[[col]])
    model.matrix(~. - 1, data = d)
  }
  
  sites_svm <- svm_prep(sites_xgb)
  bg_svm    <- svm_prep(bg_xgb)
  
  for (f in 1:5) {
    ts_s   <- which(site_folds_c == f)
    ts_b   <- which(bg_folds_pp  == f)
    if (!length(ts_s) || !length(ts_b)) next
    n_p_tr <- sum(site_folds_c != f)
    n_b_tr <- sum(bg_folds_pp  != f)
    if (n_p_tr < 5 || n_b_tr < 5) next
    
    tr_mat <- rbind(
      sites_svm[site_folds_c != f, , drop = FALSE],
      bg_svm[bg_folds_pp  != f,    , drop = FALSE])
    tr_lab <- factor(
      c(rep("pos", n_p_tr), rep("neg", n_b_tr)),
      levels = c("neg", "pos"))
    
    # FIX 9: 90s timeout on SVM (can also hang under
    # convergence failure) + suppressWarnings(tryCatch())
    set.seed(42)
    setTimeLimit(elapsed = 90, transient = TRUE)
    fsvm <- suppressWarnings(tryCatch(
      kernlab::ksvm(tr_mat, tr_lab,
                    type          = "C-svc",
                    kernel        = "rbfdot",
                    kpar          = list(sigma = 0.01),
                    C             = 1.0,
                    class.weights = CW_svm,
                    prob.model    = TRUE,
                    scaled        = FALSE),
      error = function(e) NULL))
    setTimeLimit(elapsed = Inf, transient = TRUE)
    
    if (!is.null(fsvm)) {
      ps <- tryCatch(
        kernlab::predict(
          fsvm, sites_svm[ts_s, , drop = FALSE],
          type = "probabilities")[, "pos"],
        error = function(e) rep(NA_real_, length(ts_s)))
      pb <- tryCatch(
        kernlab::predict(
          fsvm, bg_svm[ts_b, , drop = FALSE],
          type = "probabilities")[, "pos"],
        error = function(e) rep(NA_real_, length(ts_b)))
      cv_s[ts_s, "SV"] <- ps
      cv_b[ts_b, "SV"] <- pb
      fold_aucs[f, "SV"] <- calc_auc(ps, pb)
      rm(fsvm)
    }
    gc(full = TRUE)
  }
  cat(sprintf("AUC=%.4f\u00b1%.4f (%.1fmin)\n",
              safe_mean(fold_aucs[, "SV"]),
              safe_sd(fold_aucs[,   "SV"]),
              (proc.time() - t0)[3] / 60))
  
  rm(sites_svm, bg_svm, scale_p); gc(full = TRUE)
  
  # ── Evaluation summary ──────────────────────────────────
  cat("\n  Evaluation Summary:\n")
  
  auc_vec  <- numeric(6)
  eval_list <- list()
  
  for (i in seq_along(alg_codes)) {
    sp  <- cv_s[, alg_codes[i]]
    bp  <- cv_b[, alg_codes[i]]
    fa  <- fold_aucs[, alg_codes[i]]
    mau <- safe_mean(fa)
    sau <- safe_sd(fa)
    bo  <- compute_boyce(c(sp, bp), sp)
    ts  <- compute_tss(c(rep(1L, N_PRES), rep(0L, N_BG)),
                       c(sp, bp))
    kg  <- compute_kg(sp, bp)
    auc_vec[i] <- if (is.na(mau)) 0 else mau
    
    cat(sprintf(
      "  %-10s AUC=%.4f\u00b1%.4f B=%.3f T=%.3f K=%.3f\n",
      alg_names[i],
      sm(mau, 4), sm(sau, 4),
      sm(bo, 3),  sm(ts, 3),  sm(kg, 3)))
    
    eval_list[[alg_names[i]]] <- data.frame(
      period      = period_name,
      algorithm   = alg_names[i],
      cv_auc_mean = sm(mau), cv_auc_sd = sm(sau),
      boyce_index = sm(bo),  tss_max   = sm(ts),
      kvamme_gain = sm(kg),  n_sites   = N_PRES,
      n_bg        = N_BG,    spw       = round(SPW, 1),
      stringsAsFactors = FALSE)
  }
  
  # AUC-weighted sub-model ensemble
  w_sub <- auc_vec / max(sum(auc_vec, na.rm = TRUE), 1e-9)
  ens_s <- vapply(seq_len(N_PRES), function(i) {
    v <- cv_s[i, ]; ok <- is.finite(v)
    if (!any(ok)) NA_real_ else sum(v[ok] * w_sub[ok])
  }, numeric(1))
  ens_b <- vapply(seq_len(N_BG), function(i) {
    v <- cv_b[i, ]; ok <- is.finite(v)
    if (!any(ok)) NA_real_ else sum(v[ok] * w_sub[ok])
  }, numeric(1))
  
  ens_auc   <- calc_auc(ens_s, ens_b)
  ens_boyce <- compute_boyce(c(ens_s, ens_b), ens_s)
  ens_tss   <- compute_tss(c(rep(1L, N_PRES), rep(0L, N_BG)),
                           c(ens_s, ens_b))
  ens_kg    <- compute_kg(ens_s, ens_b)
  
  cat(sprintf(
    "  %-10s AUC=%.4f       B=%.3f T=%.3f K=%.3f\n",
    "Ensemble",
    sm(ens_auc, 4), sm(ens_boyce, 3),
    sm(ens_tss, 3), sm(ens_kg, 3)))
  
  eval_list[["Ensemble"]] <- data.frame(
    period = period_name, algorithm = "Ensemble_AUC_wtd",
    cv_auc_mean = sm(ens_auc), cv_auc_sd = NA_real_,
    boyce_index = sm(ens_boyce), tss_max  = sm(ens_tss),
    kvamme_gain = sm(ens_kg),   n_sites  = N_PRES,
    n_bg = N_BG, spw = round(SPW, 1),
    stringsAsFactors = FALSE)
  eval_df <- do.call(rbind, eval_list)
  
  # ── Final XGBoost + 250m SHAP ───────────────────────────
  cat(sprintf("\n  Fitting XGBoost + 250m SHAP...\n"))
  
  dtrain <- xgboost::xgb.DMatrix(
    rbind(sites_xgb, bg_xgb),
    label = c(rep(1, N_PRES), rep(0, N_BG)))
  set.seed(42)
  xgb_fin <- xgboost::xgb.train(params = xgb_params,
                                data   = dtrain,
                                nrounds = 100L,
                                verbose = 0)
  rm(dtrain); gc(full = TRUE)
  
  # 250m prediction grid
  tmpl_250 <- terra::rast(
    ext = terra::ext(template_30m), res = 250,
    crs = terra::crs(template_30m))
  p250m <- terra::resample(pred_stack, tmpl_250,
                           method = "bilinear")
  p250m <- terra::mask(p250m, boundary_vect)
  gdf   <- terra::as.data.frame(p250m, xy = TRUE,
                                cells = TRUE, na.rm = TRUE)
  cids  <- gdf$cell
  nc    <- nrow(gdf)
  
  gmat <- as.matrix(gdf[, final_names, drop = FALSE])
  for (col in cat_in_stack[cat_in_stack %in% final_names])
    gmat[, col] <- as.numeric(as.integer(gmat[, col]))
  storage.mode(gmat) <- "double"
  colnames(gmat)     <- final_names
  
  CHUNK  <- 50000L
  n_ch   <- ceiling(nc / CHUNK)
  S_list <- vector("list", n_ch)
  
  for (ch in seq_len(n_ch)) {
    i1 <- (ch - 1L) * CHUNK + 1L
    i2 <- min(ch * CHUNK, nc)
    cm  <- gmat[i1:i2, , drop = FALSE]
    raw <- predict(xgb_fin,
                   xgboost::xgb.DMatrix(cm),
                   predcontrib = TRUE)
    S_list[[ch]] <- raw[, seq_len(length(final_names)),
                        drop = FALSE]
    rm(raw, cm); gc(full = TRUE)
  }
  S_grid <- do.call(rbind, S_list)
  colnames(S_grid) <- final_names
  rm(S_list, xgb_fin, p250m, gdf, gmat); gc(full = TRUE)
  
  dom_idx  <- apply(abs(S_grid), 1, which.max)
  dom_freq <- sort(table(final_names[dom_idx]),
                   decreasing = TRUE)
  
  cat("  Top 3 dominant drivers:\n")
  for (dr in names(dom_freq)[1:min(3, length(dom_freq))])
    cat(sprintf("    %-22s %.1f%%\n",
                dr, 100 * as.numeric(dom_freq[dr]) / nc))
  
  dom_rast <- tmpl_250
  terra::values(dom_rast) <- NA_integer_
  dom_rast[cids] <- as.integer(dom_idx)
  dom_rast <- terra::mask(dom_rast, boundary_vect)
  rm(S_grid); gc(full = TRUE)
  
  list(period = period_name, eval_df = eval_df,
       dom_rast = dom_rast, dom_freq = dom_freq,
       n_sites = N_PRES, ens_auc = ens_auc,
       final_names = final_names)
}

# ─────────────────────────────────────────────────────────────
# 4. RUN LP
# ─────────────────────────────────────────────────────────────

result_LP <- run_period(
  period_name  = "LP",
  sites_path   = file.path(OUT_SITES,
                           "sites_thinned_LP.gpkg"),
  pred_stack   = pred_stack_base,
  final_names  = final_names_base,
  cat_in_stack = cat_base)

# ─────────────────────────────────────────────────────────────
# 5. RUN MP
# ─────────────────────────────────────────────────────────────

result_MP <- run_period(
  period_name  = "MP",
  sites_path   = file.path(OUT_SITES,
                           "sites_thinned_MP.gpkg"),
  pred_stack   = pred_stack_base,
  final_names  = final_names_base,
  cat_in_stack = cat_base)

# ─────────────────────────────────────────────────────────────
# 6. RUN UP + CHELSA
# ─────────────────────────────────────────────────────────────

cat("\n════════════════════════════════════\n")
cat("UP: Adding CHELSA-TraCE21k LGM\n")
cat("════════════════════════════════════\n\n")

chelsa_files <- c(
  BIO1  = file.path(OUT_CLIMATE_PRO,
                    "CHELSA_BIO1_LGM_30m_utm44n.tif"),
  BIO12 = file.path(OUT_CLIMATE_PRO,
                    "CHELSA_BIO12_LGM_30m_utm44n.tif"),
  BIO15 = file.path(OUT_CLIMATE_PRO,
                    "CHELSA_BIO15_LGM_30m_utm44n.tif"))
chelsa_ok <- all(file.exists(chelsa_files))

if (chelsa_ok) {
  chelsa_rasts <- lapply(names(chelsa_files), function(v) {
    r <- terra::rast(chelsa_files[v]); names(r) <- v; r})
  pred_stack_UP <- c(pred_stack_base,
                     terra::rast(chelsa_rasts))
  cat(sprintf("  Stack: %d layers\n",
              terra::nlyr(pred_stack_UP)))
  
  # VIF on augmented continuous predictors
  cont_base <- final_names_base[
    !final_names_base %in% cat_base]
  cont_aug  <- c(cont_base, "BIO1", "BIO12", "BIO15")
  aug_ok    <- cont_aug[cont_aug %in% names(pred_stack_UP)]
  aug_sub   <- terra::subset(pred_stack_UP, aug_ok)
  set.seed(42)
  vif_samp <- terra::spatSample(aug_sub, size = 10000,
                                method = "random",
                                na.rm = TRUE, as.df = TRUE)
  vif_res <- tryCatch(usdm::vifcor(vif_samp, th = 5),
                      error = function(e) NULL)
  if (!is.null(vif_res)) {
    ret_cont <- vif_res@results$Variables
    excl     <- cont_aug[!cont_aug %in% ret_cont]
    cat(sprintf("  Retained: %s\n",
                paste(ret_cont, collapse = ", ")))
    if (length(excl))
      cat(sprintf("  Excluded VIF>5: %s\n",
                  paste(excl, collapse = ", ")))
  } else {
    ret_cont <- cont_aug
  }
  final_names_UP  <- c(ret_cont, cat_base)
  cat_UP          <- cat_base
  
  # FIX 7: bg_folds_orig[bg_ok_UP] — correct length match
  cat("  Extracting UP bg...\n")
  bg_raw_UP   <- terra::extract(pred_stack_UP,
                                terra::vect(bg_sf), ID = FALSE)
  bg_ok_UP    <- complete.cases(bg_raw_UP)
  bg_vals_UP  <- bg_raw_UP[bg_ok_UP, ]
  bg_folds_UP <- bg_folds_orig[bg_ok_UP]  # indexes 10000-vec
  cat(sprintf("  UP background: %d points\n\n",
              nrow(bg_vals_UP)))
  
  result_UP <- run_period(
    period_name  = "UP",
    sites_path   = file.path(OUT_SITES,
                             "sites_thinned_UP.gpkg"),
    pred_stack   = pred_stack_UP,
    final_names  = final_names_UP,
    cat_in_stack = cat_UP,
    rl           = raster_levels,
    bg_vals_p    = bg_vals_UP,
    bg_folds_pp  = bg_folds_UP)
} else {
  cat("  CHELSA missing — UP uses base stack\n\n")
  result_UP <- run_period(
    period_name  = "UP",
    sites_path   = file.path(OUT_SITES,
                             "sites_thinned_UP.gpkg"),
    pred_stack   = pred_stack_base,
    final_names  = final_names_base,
    cat_in_stack = cat_base)
}

# ─────────────────────────────────────────────────────────────
# 7. ATTRIBUTION SENSITIVITY (Supp S3)
# ─────────────────────────────────────────────────────────────

cat("\n════════════════════════════════════\n")
cat("Attribution Sensitivity (Supp S3)\n")
cat("════════════════════════════════════\n\n")

sens_list <- list()

for (per in c("LP", "MP")) {
  cat(sprintf("  [%s] ", per))
  s1 <- file.path(OUT_SITES,
                  sprintf("sites_thinned_%s_confirmed.gpkg", per))
  s2 <- file.path(OUT_SITES,
                  sprintf("sites_thinned_%s.gpkg",           per))
  stages <- list()
  
  for (stg in c("S1", "S2")) {
    path <- if (stg == "S1") s1 else s2
    if (!file.exists(path)) next
    sf_s <- sf::st_read(path, quiet = TRUE)
    N_s  <- nrow(sf_s)
    if (N_s < 15) {
      cat(sprintf("%s N=%d<15 ", stg, N_s)); next
    }
    v  <- terra::extract(pred_stack_base,
                         terra::vect(sf_s), ID = FALSE)
    ok <- complete.cases(v); v <- v[ok, ]
    N_v <- nrow(v)
    
    xp <- list(objective = "binary:logistic",
               max_depth = 6L, eta = 0.05,
               colsample_bytree = 0.8, subsample = 0.8,
               scale_pos_weight = N_BG_BASE / N_v,
               nthread = 1L, seed = 42L)
    
    xm    <- make_xgb_mat(v, final_names_base, cat_base)
    bg_xm <- make_xgb_mat(bg_vals_base, final_names_base,
                          cat_base)
    dm    <- xgboost::xgb.DMatrix(
      rbind(xm, bg_xm),
      label = c(rep(1, N_v), rep(0, N_BG_BASE)))
    
    set.seed(42)
    xm_fit <- xgboost::xgb.train(params = xp, data = dm,
                                 nrounds = 100L,
                                 verbose = 0)
    rm(dm)
    
    dm2  <- xgboost::xgb.DMatrix(rbind(xm, bg_xm))
    shap <- predict(xm_fit, dm2, predcontrib = TRUE)
    shap <- shap[, seq_len(length(final_names_base)),
                 drop = FALSE]
    colnames(shap) <- final_names_base
    imp  <- sort(colMeans(abs(shap)), decreasing = TRUE)
    stages[[stg]] <- list(n = N_v, top3 = names(imp)[1:3])
    rm(xm_fit, shap, dm2, xm, bg_xm); gc(full = TRUE)
  }
  
  if (length(stages) == 2) {
    t1  <- stages$S1$top3; t2 <- stages$S2$top3
    stb <- setequal(t1, t2)
    cat(sprintf(
      "S1(N=%d)[%s] vs S2(N=%d)[%s] → %s\n",
      stages$S1$n, paste(t1, collapse = ","),
      stages$S2$n, paste(t2, collapse = ","),
      if (stb) "STABLE" else "CHANGED"))
    sens_list[[per]] <- data.frame(
      period = per,
      n_s1   = stages$S1$n, n_s2 = stages$S2$n,
      t1_s1 = t1[1], t2_s1 = t1[2], t3_s1 = t1[3],
      t1_s2 = t2[1], t2_s2 = t2[2], t3_s2 = t2[3],
      stable = stb, stringsAsFactors = FALSE)
  } else {
    cat("insufficient stages\n")
  }
}

if (length(sens_list) > 0) {
  write.csv(do.call(rbind, sens_list),
            file.path(OUT_SUPP_AN,
                      "SuppS3_attribution_sensitivity.csv"),
            row.names = FALSE)
  cat("\n  ✓ SuppS3_attribution_sensitivity.csv\n")
}

# ─────────────────────────────────────────────────────────────
# 8. SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────

cat("\n--- Saving Outputs ---\n\n")

all_eval <- do.call(rbind, list(
  result_LP$eval_df,
  result_MP$eval_df,
  result_UP$eval_df))
write.csv(all_eval,
          file.path(OUT_EVAL, "submodel_evaluation.csv"),
          row.names = FALSE)
cat("  ✓ submodel_evaluation.csv\n")

for (res in list(result_LP, result_MP, result_UP)) {
  pd <- file.path(OUT_MOD_SUB, res$period)
  dir.create(pd, recursive = TRUE, showWarnings = FALSE)
  outf <- file.path(pd, sprintf(
    "%s_dominant_driver_250m.tif", res$period))
  terra::writeRaster(res$dom_rast, outf,
                     overwrite = TRUE, datatype = "INT1U")
  cat(sprintf("  ✓ %s dominant driver raster\n",
              res$period))
}

# ─────────────────────────────────────────────────────────────
# 9. FIGURE 12
# ─────────────────────────────────────────────────────────────

cat("\n--- Figure 12: Diachronic Driver Maps ---\n\n")

tryCatch({
  png(file.path(OUT_FIG_MAIN,
                "Fig12_diachronic_driver_maps.png"),
      width = 7200, height = 2800, res = 300)
  par(mfrow = c(1, 3), mar = c(2, 2, 3.5, 1))
  
  for (res in list(result_LP, result_MP, result_UP)) {
    fn_p   <- res$final_names
    freq   <- res$dom_freq
    active <- names(freq)[as.numeric(freq) > 0]
    a_codes <- match(active, fn_p)
    a_cols  <- vapply(active, to_col, character(1))
    
    rcl    <- cbind(from  = a_codes,
                    to    = seq_along(a_codes))
    dr_plt <- suppressWarnings(
      terra::classify(res$dom_rast, rcl, others = NA))
    
    terra::plot(dr_plt,
                main   = sprintf("%s Sub-Model\nTop: %s",
                                 res$period, names(freq)[1]),
                col    = a_cols, type = "classes",
                axes   = FALSE, legend = FALSE, cex.main = 0.80)
    terra::plot(boundary_vect, add = TRUE,
                border = "black", lwd = 0.8)
    
    tot  <- sum(as.numeric(freq[active]))
    leg  <- sprintf("%s (%.1f%%)", active,
                    100 * as.numeric(freq[active]) / tot)
    legend("bottomright", legend = leg, fill = a_cols,
           bty = "n", cex = 0.45,
           ncol = if (length(active) > 7) 2L else 1L)
  }
  dev.off()
  cat("  ✓ Fig12_diachronic_driver_maps.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  cat("  ✗ Fig12 error:", e$message, "\n\n")
})

# ─────────────────────────────────────────────────────────────
# 10. SUMMARY
# ─────────────────────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 20 COMPLETE — Sub-Models\n")
cat("========================================\n")

for (res in list(result_LP, result_MP, result_UP)) {
  er  <- res$eval_df
  ens <- er[er$algorithm == "Ensemble_AUC_wtd", ]
  rf  <- er[er$algorithm == "RF", ]
  cat(sprintf("\n%s (N=%d):\n", res$period, res$n_sites))
  cat(sprintf("  RF:       AUC=%.4f  B=%.4f\n",
              ens$cv_auc_mean, ens$boyce_index))
  cat(sprintf("  Ensemble: AUC=%.4f  B=%.4f\n",
              rf$cv_auc_mean,  rf$boyce_index))
  top3 <- names(res$dom_freq)[
    1:min(3, length(res$dom_freq))]
  cat(sprintf("  Drivers:  %s\n",
              paste(top3, collapse = " > ")))
}

cat("\nFixes applied:\n")
cat("  FIX 1 — GAM setTimeLimit(60s/fold) ✓\n")
cat("  FIX 2 — GAM k=3 (safe for N~60-70) ✓\n")
cat("  FIX 3 — GAM gam.control(maxit=100) ✓\n")
cat("  FIX 4 — GAM/BRT/SVM sqrt(SPW) wt   ✓\n")
cat("  FIX 5 — suppressWarnings(tryCatch)  ✓\n")
cat("  FIX 6 — Single fgam definition      ✓\n")
cat("  FIX 7 — bg_folds_orig[bg_ok_UP]     ✓\n")
cat("  FIX 8 — safe_mean/safe_sd           ✓\n")
cat("  FIX 9 — SVM setTimeLimit(90s/fold)  ✓\n")
cat("\nFiles:\n")
cat("  outputs/evaluation/submodel_evaluation.csv\n")
cat("  models/submodels/LP|MP|UP/*_dominant_driver.tif\n")
cat("  outputs/supplementary_analyses/SuppS3_*.csv\n")
cat("  outputs/figures/main/Fig12_*.png\n")
cat("\nNext: Script 21 — Transfer Validation\n")
cat("========================================\n")