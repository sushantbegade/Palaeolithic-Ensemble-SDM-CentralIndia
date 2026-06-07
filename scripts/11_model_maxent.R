# ============================================================
# SCRIPT 11: MAXENT MODEL (maxnet)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 11 of 25
# ============================================================
# PARAMETER TUNING — SEQUENTIAL MANUAL GRID SEARCH:
#   Full grid: 5 FC × 7 RM = 35 combinations
#   Phase 1: Fold-1 AUC screening (one model at a time)
#   Phase 2: Full 5-fold CV on top 5 candidates
#   No ENMeval — avoids simultaneous RAM accumulation crash
#
# TILED PREDICTION: 4 horizontal tiles (~6M cells each)
# ALL predictions type="logistic" — never cloglog
# BOYCE FIX: safe_scalar() ensures all metrics are length 1
# ============================================================

# ── 0. SOURCE SCRIPT 01 ─────────────────────────────────────

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 11: MaxEnt Model\n")
cat("========================================\n\n")

set.seed(42)
cat("R tempdir:   ", tempdir(), "\n")
cat("terra tmpdir:", terra::terraOptions()$tempdir, "\n\n")

# ── HELPERS ──────────────────────────────────────────────────

# Ensure any metric is exactly length-1 numeric
safe_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  x[length(x)]
}

# Manual TSS — no PresenceAbsence dependency
compute_tss <- function(obs, pred,
                        thresholds = seq(0.01, 0.99, 0.01)) {
  tss_vals <- sapply(thresholds, function(t) {
    tp   <- sum(obs==1 & pred>=t, na.rm=TRUE)
    tn   <- sum(obs==0 & pred< t, na.rm=TRUE)
    fp   <- sum(obs==0 & pred>=t, na.rm=TRUE)
    fn   <- sum(obs==1 & pred< t, na.rm=TRUE)
    sens <- if (tp+fn>0) tp/(tp+fn) else NA_real_
    spec <- if (tn+fp>0) tn/(tn+fp) else NA_real_
    sens + spec - 1
  })
  safe_scalar(max(tss_vals, na.rm=TRUE))
}

# Fit one maxnet model
fit_maxnet <- function(pres, data, fc, rm, cats) {
  set.seed(42)
  maxnet::maxnet(
    p            = pres,
    data         = data,
    f            = maxnet::maxnet.formula(pres, data,
                                          classes = tolower(fc)),
    regmult      = rm,
    categoricals = cats
  )
}

# Compute AUC from presence/background predictions
calc_auc <- function(ps, pb) {
  as.numeric(pROC::auc(
    pROC::roc(c(rep(1,length(ps)), rep(0,length(pb))),
              c(ps, pb), quiet=TRUE)))
}

# ── 1. LOAD INPUTS ───────────────────────────────────────────

cat("--- Loading Inputs ---\n")

cv_design  <- readRDS(file.path(OUT_CV,
                                "cv_block_assignments.rds"))
site_folds <- cv_design$site_folds
bg_folds   <- cv_design$bg_folds

sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet = TRUE)
bg_sf    <- sf::st_read(file.path(OUT_CV,
                                  "background_with_folds.gpkg"),
                        quiet = TRUE)

final_names <- readRDS(file.path(OUT_PREDICTORS,
                                 "final_predictor_names.rds"))
pred_stack  <- terra::rast(file.path(OUT_PREDICTORS,
                                     "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))

boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

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

for (col in cat_in_stack) {
  sites_vals[[col]] <- as.integer(sites_vals[[col]])
  bg_vals[[col]]    <- as.integer(bg_vals[[col]])
}

cat(sprintf("  Sites: %d  Background: %d\n\n",
            nrow(sites_vals), nrow(bg_vals)))
gc(full=TRUE)

# ── 3. PHASE 1 — FOLD-1 SCREENING (35 combos) ───────────────

cat("--- Phase 1: Fold-1 Screening (35 combinations) ---\n\n")

fc_list <- c("L","LQ","LQH","LQHP","LQHPT")
rm_list <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)

# Build fold-1 train/test split
scr_tr_pres <- c(rep(1L,sum(site_folds_c!=1)),
                 rep(0L,sum(bg_folds_c!=1)))
scr_tr_data <- rbind(sites_vals[site_folds_c!=1,],
                     bg_vals[bg_folds_c!=1,])
scr_ts_s    <- which(site_folds_c==1)
scr_ts_b    <- which(bg_folds_c==1)

screen_results <- data.frame(fc=character(), rm=numeric(),
                             fold1_auc=numeric(),
                             stringsAsFactors=FALSE)
t0 <- proc.time(); n <- 0

for (fc in fc_list) {
  for (rm in rm_list) {
    n <- n + 1
    cat(sprintf("  [%02d/35] FC=%-5s RM=%.1f ... ", n, fc, rm))
    m <- tryCatch(
      fit_maxnet(scr_tr_pres, scr_tr_data,
                 fc, rm, cat_in_stack),
      error = function(e) NULL)
    if (is.null(m)) {
      cat("FAILED\n")
      auc_f1 <- NA_real_
    } else {
      ps <- predict(m, sites_vals[scr_ts_s,], type="logistic")
      pb <- predict(m, bg_vals[scr_ts_b,],    type="logistic")
      auc_f1 <- calc_auc(ps, pb)
      cat(sprintf("AUC=%.4f\n", auc_f1))
      rm(m, ps, pb)
    }
    screen_results <- rbind(screen_results,
                            data.frame(fc=fc, rm=rm, fold1_auc=auc_f1,
                                       stringsAsFactors=FALSE))
    gc(full=TRUE)
  }
}

cat(sprintf("\n  Phase 1 done in %.1f min\n",
            (proc.time()-t0)[3]/60))

rm(scr_tr_pres, scr_tr_data); gc(full=TRUE)

screen_results <- screen_results[
  order(-screen_results$fold1_auc, na.last=TRUE),]

cat("\n  Top 10 by fold-1 AUC:\n")
cat(sprintf("  %-6s %-5s  %s\n","FC","RM","AUC"))
cat("  ", paste(rep("-",22),collapse=""), "\n")
for (i in seq_len(min(10,nrow(screen_results)))) {
  cat(sprintf("  %-6s %-5.1f  %.4f\n",
              screen_results$fc[i], screen_results$rm[i],
              screen_results$fold1_auc[i]))
}

# ── 4. PHASE 2 — FULL CV ON TOP 5 ───────────────────────────

cat("\n--- Phase 2: Full 5-Fold CV (top 5 candidates) ---\n\n")

top5 <- head(screen_results[!is.na(screen_results$fold1_auc),], 5)
full_cv <- data.frame(fc=character(), rm=numeric(),
                      cv_auc_mean=numeric(), cv_auc_sd=numeric(),
                      stringsAsFactors=FALSE)

for (r in seq_len(nrow(top5))) {
  fc_i <- top5$fc[r]; rm_i <- top5$rm[r]
  cat(sprintf("  FC=%-5s RM=%.1f :", fc_i, rm_i))
  fa <- numeric(5)
  for (f in 1:5) {
    tr_p <- c(rep(1L,sum(site_folds_c!=f)),
              rep(0L,sum(bg_folds_c!=f)))
    tr_d <- rbind(sites_vals[site_folds_c!=f,],
                  bg_vals[bg_folds_c!=f,])
    ts_s <- which(site_folds_c==f)
    ts_b <- which(bg_folds_c==f)
    m <- tryCatch(fit_maxnet(tr_p,tr_d,fc_i,rm_i,cat_in_stack),
                  error=function(e) NULL)
    if (is.null(m)) { fa[f] <- NA_real_ } else {
      ps <- predict(m,sites_vals[ts_s,],type="logistic")
      pb <- predict(m,bg_vals[ts_b,],   type="logistic")
      fa[f] <- calc_auc(ps,pb)
      rm(m,ps,pb)
    }
    rm(tr_d,tr_p); gc(full=TRUE)
  }
  m_auc <- mean(fa,na.rm=TRUE); s_auc <- sd(fa,na.rm=TRUE)
  cat(sprintf(" AUC=%.4f±%.4f\n", m_auc, s_auc))
  full_cv <- rbind(full_cv,
                   data.frame(fc=fc_i,rm=rm_i,cv_auc_mean=m_auc,
                              cv_auc_sd=s_auc,stringsAsFactors=FALSE))
}

full_cv  <- full_cv[order(-full_cv$cv_auc_mean),]
best_fc  <- full_cv$fc[1]
best_rm  <- full_cv$rm[1]

cat(sprintf("\n  BEST: FC=%s  RM=%.1f  CV AUC=%.4f±%.4f\n\n",
            best_fc, best_rm,
            full_cv$cv_auc_mean[1], full_cv$cv_auc_sd[1]))

tuning_tbl <- merge(screen_results, full_cv,
                    by=c("fc","rm"), all.x=TRUE)
write.csv(tuning_tbl,
          file.path(OUT_TABLES,"TableS2_MaxEnt_tuning.csv"),
          row.names=FALSE)
cat("  ✓ TableS2_MaxEnt_tuning.csv\n\n")

# ── 5. FINAL MODEL ───────────────────────────────────────────

cat("--- Final Model (all data, best params) ---\n")

all_pres <- c(rep(1L,nrow(sites_vals)), rep(0L,nrow(bg_vals)))
all_vals <- rbind(sites_vals, bg_vals)

maxnet_model <- fit_maxnet(all_pres, all_vals,
                           best_fc, best_rm, cat_in_stack)

cat("  Features:", length(maxnet_model$betas), "\n")
saveRDS(maxnet_model,
        file.path(OUT_MOD_IND,"maxent_model_final.rds"))
cat("  ✓ maxent_model_final.rds\n\n")

rm(all_vals, all_pres); gc(full=TRUE)

# ── 6. FULL 5-FOLD CV (BEST PARAMS) ─────────────────────────

cat("--- Full 5-Fold CV (best params) ---\n\n")

cv_preds_sites <- numeric(nrow(sites_vals))
cv_preds_bg    <- numeric(nrow(bg_vals))
fold_aucs      <- numeric(5)

for (f in 1:5) {
  cat(sprintf("  Fold %d: ", f))
  tr_p <- c(rep(1L,sum(site_folds_c!=f)),
            rep(0L,sum(bg_folds_c!=f)))
  tr_d <- rbind(sites_vals[site_folds_c!=f,],
                bg_vals[bg_folds_c!=f,])
  ts_s <- which(site_folds_c==f)
  ts_b <- which(bg_folds_c==f)
  fm   <- fit_maxnet(tr_p,tr_d,best_fc,best_rm,cat_in_stack)
  ps   <- predict(fm,sites_vals[ts_s,],type="logistic")
  pb   <- predict(fm,bg_vals[ts_b,],   type="logistic")
  cv_preds_sites[ts_s] <- ps
  cv_preds_bg[ts_b]    <- pb
  fold_aucs[f] <- calc_auc(ps,pb)
  cat(sprintf("AUC=%.4f  (%d sites/%d bg)\n",
              fold_aucs[f],length(ts_s),length(ts_b)))
  rm(fm,tr_d,tr_p,ps,pb); gc(full=TRUE)
}

cv_auc_mean <- mean(fold_aucs)
cv_auc_sd   <- sd(fold_aucs)
cat(sprintf("\n  CV AUC: %.4f ± %.4f\n\n",
            cv_auc_mean, cv_auc_sd))

# ── 7. TILED RASTER PREDICTION ───────────────────────────────

cat("--- Tiled Raster Prediction (4 tiles) ---\n\n")

ext_full   <- terra::ext(pred_stack)
y_step     <- (ext_full[4]-ext_full[3]) / 4
tile_paths <- character(4)

for (i in 1:4) {
  cat(sprintf("  Tile %d/4 ... ", i))
  te <- terra::ext(ext_full[1], ext_full[2],
                   ext_full[3]+(i-1)*y_step,
                   ext_full[3]+i*y_step)
  ts <- terra::crop(pred_stack, te)
  tile_paths[i] <- file.path("E:/R_temp",
                             sprintf("mx_tile%d.tif",i))
  t0 <- proc.time()
  terra::predict(ts, maxnet_model, type="logistic",
                 na.rm=TRUE, filename=tile_paths[i],
                 overwrite=TRUE, wopt=list(datatype="FLT4S"))
  cat(sprintf("%.1f min\n",(proc.time()-t0)[3]/60))
  rm(ts); gc(full=TRUE)
}

cat("  Merging ... ")
out_pred <- file.path(OUT_MOD_IND,"maxent_pred_logistic.tif")
terra::merge(terra::sprc(lapply(tile_paths,terra::rast)),
             filename=out_pred, overwrite=TRUE,
             wopt=list(datatype="FLT4S"))
file.remove(tile_paths)
cat("done\n")

maxent_raster <- terra::rast(out_pred)
rng <- terra::global(maxent_raster,c("min","max","mean"),
                     na.rm=TRUE)
cat(sprintf("  Range: %.4f to %.4f (mean %.4f)\n",
            rng[1,1],rng[1,2],rng[1,3]))
if (rng[1,2]>1.001||rng[1,1]< -0.001) {
  warning("Values outside [0,1]")
} else { cat("  ✓ logistic [0,1] confirmed\n") }
cat("  ✓ maxent_pred_logistic.tif\n\n")
gc(full=TRUE)

# ── 8. EVALUATION METRICS ────────────────────────────────────

cat("--- Evaluation Metrics ---\n\n")

full_p   <- c(predict(maxnet_model,sites_vals,type="logistic"),
              predict(maxnet_model,bg_vals,    type="logistic"))
full_t   <- c(rep(1L,nrow(sites_vals)),rep(0L,nrow(bg_vals)))
auc_full <- safe_scalar(as.numeric(pROC::auc(
  pROC::roc(full_t,full_p,quiet=TRUE))))

# Boyce Index — safe_scalar ensures length-1 output
boyce_val <- tryCatch({
  res <- ecospat::ecospat.boyce(
    fit=c(cv_preds_sites,cv_preds_bg),
    obs=cv_preds_sites, PEplot=FALSE)
  safe_scalar(res$Spearman.cor)
}, error=function(e) NA_real_)

# TSS — manual
tss_val <- tryCatch(
  compute_tss(
    obs  = c(rep(1,length(cv_preds_sites)),
             rep(0,length(cv_preds_bg))),
    pred = c(cv_preds_sites,cv_preds_bg)),
  error=function(e) NA_real_)

# Kvamme's Gain
area_cells <- terra::global(!is.na(maxent_raster),
                            "sum",na.rm=TRUE)[1,1]
high_cells <- terra::global(maxent_raster>0.5,
                            "sum",na.rm=TRUE)[1,1]
area_pct   <- high_cells/area_cells
site_pred  <- terra::extract(maxent_raster,
                             terra::vect(sites_sf))[,2]
sites_pct  <- sum(site_pred>0.5,na.rm=TRUE)/length(site_pred)
kg         <- safe_scalar(1-(area_pct/max(sites_pct,1e-9)))

cat(sprintf("  CV AUC (primary):  %.4f ± %.4f\n",
            cv_auc_mean,cv_auc_sd))
cat(sprintf("  Full AUC (diag):   %.4f\n", safe_scalar(auc_full)))
cat(sprintf("  Boyce Index:       %.4f\n", safe_scalar(boyce_val)))
cat(sprintf("  TSS (max-TSS):     %.4f\n", safe_scalar(tss_val)))
cat(sprintf("  Kvamme's Gain:     %.4f\n", safe_scalar(kg)))
cat(sprintf("  Area > 0.5:        %.1f%%\n",100*area_pct))
cat(sprintf("  Sites > 0.5:       %.1f%%\n",100*sites_pct))

# ── 9. SAVE ALL OUTPUTS ──────────────────────────────────────

# safe_metric: round a scalar, return NA if not finite
sm <- function(x, d=4) {
  x <- safe_scalar(x)
  if (is.na(x) || !is.finite(x)) return(NA_real_)
  round(x, d)
}

metrics_df <- data.frame(
  algorithm   = "MaxEnt",
  best_fc     = best_fc,
  best_rm     = best_rm,
  cv_auc_mean = sm(cv_auc_mean),
  cv_auc_sd   = sm(cv_auc_sd),
  full_auc    = sm(auc_full),
  boyce_index = sm(boyce_val),
  tss_max     = sm(tss_val),
  kvamme_gain = sm(kg),
  fold1_auc   = sm(fold_aucs[1]),
  fold2_auc   = sm(fold_aucs[2]),
  fold3_auc   = sm(fold_aucs[3]),
  fold4_auc   = sm(fold_aucs[4]),
  fold5_auc   = sm(fold_aucs[5]),
  stringsAsFactors = FALSE
)

write.csv(metrics_df,
          file.path(OUT_EVAL,"maxent_evaluation.csv"),
          row.names=FALSE)

saveRDS(list(site_preds=cv_preds_sites,
             bg_preds=cv_preds_bg,
             fold_aucs=fold_aucs),
        file.path(OUT_MOD_IND,"maxent_cv_predictions.rds"))

cat("\n  ✓ maxent_evaluation.csv\n")
cat("  ✓ maxent_cv_predictions.rds\n\n")

# ── 10. FIGURE ───────────────────────────────────────────────

cat("--- Diagnostic Figure ---\n")

png(file.path(OUT_FIG_MAIN,"Fig_MaxEnt_prediction.png"),
    width=2400, height=2400, res=300)
terra::plot(maxent_raster,
            main=sprintf(
              "MaxEnt — Logistic Suitability\nCV AUC=%.4f±%.4f  FC=%s  RM=%.1f",
              cv_auc_mean,cv_auc_sd,best_fc,best_rm),
            col=viridisLite::viridis(100), range=c(0,1), axes=FALSE)
terra::plot(boundary_vect,add=TRUE,border="white",lwd=0.8)
terra::plot(terra::vect(sites_sf),add=TRUE,
            col="red",pch=16,cex=0.3)
dev.off()
cat("  ✓ Fig_MaxEnt_prediction.png\n\n")

# ── 11. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 11 COMPLETE — MaxEnt\n")
cat("========================================\n")
cat(sprintf("Best FC: %s  |  RM: %.1f\n",best_fc,best_rm))
cat(sprintf("CV AUC:        %.4f ± %.4f\n",cv_auc_mean,cv_auc_sd))
cat(sprintf("Boyce Index:   %.4f\n",safe_scalar(boyce_val)))
cat(sprintf("TSS:           %.4f\n",safe_scalar(tss_val)))
cat(sprintf("Kvamme's Gain: %.4f\n",safe_scalar(kg)))
cat("Output: logistic [0,1] ✓\n")
cat("Tiled prediction (4×6M cells) ✓\n")
cat("All I/O: E drive ✓\n")
cat("\nNext: Run Script 12 — Random Forest\n")
cat("========================================\n")