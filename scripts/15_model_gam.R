# ============================================================
# SCRIPT 15: GENERALISED ADDITIVE MODEL (GAM)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 15 of 25
# ============================================================
# SPECIFICATION (Research Design 5.7.5):
#   family: binomial(link="logit")
#   method: REML  |  select: TRUE
#   k: 5 basis dimensions for all continuous smooths
#   Continuous: s(x, k=5)
#   Categorical: s(x, bs="re") random effects smooths
#   Case weights: presence = N_BG/N_PRES, background = 1
#
# FIX v3 — RASTER UNIQUE VALUES:
#   terra::global(x,"unique") not available in all versions.
#   Replaced with terra::freq() which is memory-efficient
#   and universally supported across terra versions.
#   terra::freq(raster) returns a data frame of unique values
#   and their counts without loading all values into RAM.
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
cat("  Continuous:", length(cont_names),
    "  Categorical:", paste(cat_in_stack, collapse=", "),
    "\n\n")

# ── 2. RASTER FACTOR LEVELS VIA terra::freq() ────────────────
# terra::freq() computes unique value frequency table
# without loading full raster into memory — safe for large
# rasters and universally supported across terra versions.

cat("--- Extracting Raster Factor Levels (terra::freq) ---\n")

raster_levels <- list()
for (col in cat_in_stack) {
  r_layer  <- terra::subset(pred_stack,
                            which(names(pred_stack)==col))
  freq_tbl <- terra::freq(r_layer)
  # freq() returns data frame with columns: layer, value, count
  all_vals <- sort(unique(
    as.integer(freq_tbl$value[!is.na(freq_tbl$value)])))
  raster_levels[[col]] <- all_vals
  cat(sprintf("  %s: %d unique values (%d to %d)\n",
              col, length(all_vals),
              min(all_vals), max(all_vals)))
}
cat("\n")

# ── 3. EXTRACT PREDICTOR VALUES ──────────────────────────────

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

# Apply raster-derived factor levels
for (col in cat_in_stack) {
  lvls <- raster_levels[[col]]
  sites_vals[[col]] <- factor(as.integer(sites_vals[[col]]),
                              levels=lvls)
  bg_vals[[col]]    <- factor(as.integer(bg_vals[[col]]),
                              levels=lvls)
}

N_PRES <- nrow(sites_vals)
N_BG   <- nrow(bg_vals)
WT     <- N_BG / N_PRES

cat(sprintf("  Sites: %d  Background: %d\n",
            N_PRES, N_BG))
cat(sprintf("  Presence case weight: %.2f\n\n", WT))

# ── 4. BUILD GAM FORMULA ─────────────────────────────────────

cat("--- Building GAM Formula ---\n\n")

smooth_terms <- c(
  paste0("s(", cont_names, ", k=5)"),
  paste0("s(", cat_in_stack, ", bs='re')")
)

gam_formula <- as.formula(
  paste("presence ~", paste(smooth_terms, collapse=" + ")))

for (term in smooth_terms) cat("  +", term, "\n")
cat("\n")

# ── 5. FIT FINAL GAM ─────────────────────────────────────────

cat("--- Fitting Final GAM ---\n")
cat("  method: REML  select: TRUE\n")
cat(sprintf("  Presence weight: %.1f  Background: 1\n\n", WT))

all_resp    <- c(rep(1L,N_PRES), rep(0L,N_BG))
all_data    <- rbind(sites_vals, bg_vals)
all_df      <- cbind(presence=all_resp, all_data)
all_weights <- c(rep(WT,N_PRES), rep(1,N_BG))

set.seed(42)
t0 <- proc.time()

gam_model <- mgcv::gam(
  formula = gam_formula,
  family  = binomial(link="logit"),
  method  = "REML",
  select  = TRUE,
  weights = all_weights,
  data    = all_df
)

cat(sprintf("  Done in %.1f min\n",
            (proc.time()-t0)[3]/60))

gam_sum <- summary(gam_model)
cat("\n  Smooth term EDFs:\n\n")

edf_df <- data.frame(
  term    = rownames(gam_sum$s.table),
  edf     = round(gam_sum$s.table[,"edf"], 3),
  p_value = round(gam_sum$s.table[,"p-value"], 4),
  stringsAsFactors=FALSE
)
for (i in seq_len(nrow(edf_df))) {
  note <- if (edf_df$edf[i] < 0.05) " ← shrunk" else ""
  cat(sprintf("  %-32s EDF=%.3f  p=%.4f%s\n",
              edf_df$term[i], edf_df$edf[i],
              edf_df$p_value[i], note))
}

dev_expl <- 100*(1-gam_model$deviance/gam_model$null.deviance)
cat(sprintf("\n  Deviance explained: %.1f%%\n", dev_expl))

saveRDS(gam_model,
        file.path(OUT_MOD_IND,"gam_model_final.rds"))
saveRDS(raster_levels,
        file.path(OUT_MOD_IND,"gam_raster_levels.rds"))
cat("  ✓ gam_model_final.rds\n")
cat("  ✓ gam_raster_levels.rds\n\n")

rm(all_df, all_data, all_resp, all_weights); gc(full=TRUE)

# ── 6. 5-FOLD SPATIAL BLOCK CV ───────────────────────────────

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
  tr_wts  <- c(rep(WT,sum(site_folds_c!=f)),
               rep(1, sum(bg_folds_c!=f)))
  tr_df   <- cbind(presence=tr_resp, tr_data)
  
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  
  set.seed(42)
  fold_gam <- tryCatch(
    mgcv::gam(gam_formula,
              family=binomial(link="logit"),
              method="REML", select=TRUE,
              weights=tr_wts, data=tr_df),
    error=function(e) {
      cat("ERROR:", conditionMessage(e),"\n"); NULL
    }
  )
  
  if (is.null(fold_gam)) {
    fold_aucs[f] <- NA_real_; next
  }
  
  ts_s_d <- sites_vals[ts_s,]
  ts_b_d <- bg_vals[ts_b,]
  for (col in cat_in_stack) {
    ts_s_d[[col]] <- factor(as.integer(ts_s_d[[col]]),
                            levels=raster_levels[[col]])
    ts_b_d[[col]] <- factor(as.integer(ts_b_d[[col]]),
                            levels=raster_levels[[col]])
  }
  
  ps <- suppressWarnings(predict(fold_gam, ts_s_d,
                                 type="response"))
  pb <- suppressWarnings(predict(fold_gam, ts_b_d,
                                 type="response"))
  
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f] <- calc_auc(ps[is.finite(ps)],
                           pb[is.finite(pb)])
  
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f], length(ts_s), length(ts_b)))
  
  rm(fold_gam, tr_df, tr_data, tr_resp, tr_wts,
     ts_s_d, ts_b_d, ps, pb)
  gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs, na.rm=TRUE)
cv_auc_sd   <- sd(fold_aucs, na.rm=TRUE)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

gam_tile_fn <- function(model, data, ...) {
  result <- rep(NA_real_, nrow(data))
  for (col in cat_in_stack) {
    if (col %in% names(data)) {
      data[[col]] <- factor(as.integer(data[[col]]),
                            levels=raster_levels[[col]])
    }
  }
  ok <- complete.cases(data)
  if (any(ok)) {
    result[ok] <- suppressWarnings(
      predict(model, newdata=data[ok,,drop=FALSE],
              type="response"))
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
  terra::predict(ts, gam_model, fun=gam_tile_fn,
                 na.rm=FALSE, filename=tile_paths[i],
                 overwrite=TRUE,
                 wopt=list(datatype="FLT4S"))
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
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

# ── 8. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_ps <- suppressWarnings(
  predict(gam_model, sites_vals, type="response"))
full_pb <- suppressWarnings(
  predict(gam_model, bg_vals, type="response"))
auc_full <- safe_scalar(calc_auc(
  full_ps[is.finite(full_ps)],
  full_pb[is.finite(full_pb)]))

cv_s_v <- cv_preds_sites[is.finite(cv_preds_sites)]
cv_b_v <- cv_preds_bg[is.finite(cv_preds_bg)]

boyce_val <- compute_boyce(c(cv_s_v,cv_b_v), cv_s_v)
tss_val   <- compute_tss(
  obs  = c(rep(1,length(cv_s_v)),rep(0,length(cv_b_v))),
  pred = c(cv_s_v,cv_b_v))

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
cat(sprintf("  Deviance expl.:    %.1f%%\n", dev_expl))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

shrunk <- edf_df$term[edf_df$edf < 0.05]
cat("  Shrunk to zero:", if(length(shrunk)==0) "none" else
  paste(shrunk, collapse=", "), "\n")

# ── 9. SAVE OUTPUTS ──────────────────────────────────────────

metrics_df <- data.frame(
  algorithm="GAM", method="REML", select=TRUE,
  case_wt=round(WT,2), deviance_expl=round(dev_expl,2),
  shrunk_terms=paste(shrunk,collapse=";"),
  cv_auc_mean=sm(cv_auc_mean), cv_auc_sd=sm(cv_auc_sd),
  full_auc=sm(auc_full), boyce_index=sm(boyce_val),
  tss_max=sm(tss_val), kvamme_gain=sm(kg),
  fold1_auc=sm(fold_aucs[1]), fold2_auc=sm(fold_aucs[2]),
  fold3_auc=sm(fold_aucs[3]), fold4_auc=sm(fold_aucs[4]),
  fold5_auc=sm(fold_aucs[5]), stringsAsFactors=FALSE
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

cat("\n  ✓ gam_evaluation.csv\n")
cat("  ✓ gam_smooth_edfs.csv\n")
cat("  ✓ gam_cv_predictions.rds\n\n")

# ── 10. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")
png(file.path(OUT_FIG_MAIN,"Fig_GAM_prediction.png"),
    width=2400,height=2400,res=300)
terra::plot(gam_raster,
            main=sprintf("GAM — Probability\nCV AUC=%.4f±%.4f  Dev=%.1f%%",
                         cv_auc_mean,cv_auc_sd,dev_expl),
            col=viridisLite::viridis(100),range=c(0,1),axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_GAM_prediction.png\n\n")

# ── 11. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 15 COMPLETE — GAM\n")
cat("========================================\n")
cat("method=REML  select=TRUE\n")
cat(sprintf("Case weight: presence=%.1f  background=1\n",WT))
cat(sprintf("Factor levels from terra::freq() ✓\n"))
cat(sprintf("Deviance expl.: %.1f%%\n", dev_expl))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n", sm(boyce_val)))
cat(sprintf("TSS:           %.4f\n", sm(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n", sm(kg)))
cat("Output: probability [0,1] ✓\n")
cat("Tiled prediction (4 tiles) ✓\n")
cat("\nNext: Run Script 16 — SVM\n")
cat("========================================\n")