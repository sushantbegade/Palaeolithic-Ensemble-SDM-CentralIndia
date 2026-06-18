# ============================================================
# SCRIPT 19: SHAP ANALYSIS
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 19 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   TreeSHAP via shapviz for XGBoost (exact, native, fast)
#   TreeSHAP via treeshap for RF (exact, fast subsample)
#   kernelshap fallback for RF if treeshap unavailable
#
#   Outputs:
#   Table 4: global SHAP importance (mean |SHAP|) both models
#   Figure 9: SHAP beeswarm (XGBoost + RF side by side)
#   Figure 10: SHAP dependence plots (top 5 XGBoost predictors)
#   Figure 11: SPATIAL dominant driver map — CENTRAL FIGURE
#   Cohen's kappa: RF vs XGBoost dominant drivers at subsample
#   Chi-squared: spatial heterogeneity across landscape zones
#   shap_global_importance.csv + shap_dominant_driver_map.tif
#
# SHAP SPATIAL GRID: 250m resolution
#   Created by resampling 30m predictor stack to 250m.
#   Within-site catchment scale (Whallon 2006).
#   ~280k-340k valid cells; computed in 50k-cell chunks.
#   XGBoost SHAP on full grid; RF SHAP on 500-row subsample.
#
# LANDSCAPE ZONES (elevation-based, Section 7.2):
#   Zone 1 — Nagpur pediplain:   DEM 270–380 m
#   Zone 2 — Chandrapur basin:   DEM < 270 m
#   Zone 3 — Satpura foothills:  DEM > 380 m
#   Verify thresholds against study area DEM statistics.
#   Adjust before manuscript writing if needed.
#
# COHEN'S KAPPA:
#   Computed at 500-row training subsample locations.
#   If kappa > 0.65 → single averaged importance discussed.
#   If kappa < 0.65 → RF and XGBoost reported separately.
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

cat("\n========================================\n")
cat("SCRIPT 19: SHAP Analysis\n")
cat("========================================\n\n")

set.seed(42)

# ── 1. INSTALL TREESHAP IF NEEDED ────────────────────────────

cat("--- Package Check ---\n")
if (!requireNamespace("treeshap", quietly = TRUE)) {
  cat("  Installing treeshap...\n")
  install.packages("treeshap")
}
library(treeshap)
library(shapviz)
cat("  treeshap:", as.character(packageVersion("treeshap")), "✓\n")
cat("  shapviz: ", as.character(packageVersion("shapviz")),  "✓\n\n")

# ── 2. HELPERS ───────────────────────────────────────────────

# Extract SHAP matrix from shapviz object (version-safe)
get_shap_matrix <- function(sv) {
  if (!is.null(sv$S)) return(as.matrix(sv$S))
  stop("Cannot extract SHAP matrix from shapviz object")
}

# Cohen's kappa (no extra package)
compute_kappa <- function(x, y) {
  lvls <- union(unique(x), unique(y))
  x    <- factor(x, levels = lvls)
  y    <- factor(y, levels = lvls)
  n    <- length(x)
  k    <- length(lvls)
  conf <- table(x, y)
  Po   <- sum(diag(conf)) / n
  Pe   <- sum((rowSums(conf)/n) * (colSums(conf)/n))
  if (abs(1 - Pe) < 1e-10) return(NA_real_)
  (Po - Pe) / (1 - Pe)
}

# 13-predictor colour palette (visually distinct)
PRED_COLORS <- c(
  "Elevation"          = "#E63946",
  "Aspect"             = "#F4A261",
  "TRI"                = "#2A9D8F",
  "TPI"                = "#264653",
  "Plan_Curvature"     = "#E9C46A",
  "HAND"               = "#06D6A0",
  "Flow_Accum_log10"   = "#118AB2",
  "Dist_River"         = "#073B4C",
  "Dist_Palaeochannel" = "#8B2FC9",
  "Dist_RawMat"        = "#FF6B6B",
  "NDVI"               = "#56CFE1",
  "Geology"            = "#FF9F1C",
  "Geomorphology"      = "#4CC9F0"
)

# ── 3. LOAD MODELS + METADATA ────────────────────────────────

cat("--- Loading Models ---\n\n")

# XGBoost
xgb_model <- xgboost::xgb.load(
  file.path(OUT_MOD_IND, "xgboost_model_final.bin"))
xgb_info  <- readRDS(
  file.path(OUT_MOD_IND, "xgboost_model_info.rds"))
# Set feature names (may not be stored in .bin file)
xgb_model$feature_names <- xgb_info$feature_names
cat("  XGBoost loaded. nrounds =", xgb_info$nrounds, "\n")
cat("  Feature names:", paste(xgb_info$feature_names, collapse=", "), "\n\n")

# Random Forest
rf_model <- readRDS(file.path(OUT_MOD_IND, "rf_model_final.rds"))
cat("  RF loaded. ntree =", rf_model$ntree, "\n\n")

# Predictor metadata
final_names   <- readRDS(file.path(OUT_PREDICTORS,
                                   "final_predictor_names.rds"))
raster_levels <- readRDS(file.path(OUT_MOD_IND,
                                   "gam_raster_levels.rds"))
cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

pred_stack <- terra::rast(file.path(OUT_PREDICTORS,
                                    "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))

# ── 4. EXTRACT TRAINING DATA ─────────────────────────────────

cat("--- Extracting Training Data ---\n\n")

sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet = TRUE)
bg_sf    <- sf::st_read(file.path(OUT_BACKGROUND,
                                  "background_N10000.gpkg"),
                        quiet = TRUE)

sites_raw <- terra::extract(pred_stack,
                            terra::vect(sites_sf), ID = FALSE)
bg_raw    <- terra::extract(pred_stack,
                            terra::vect(bg_sf), ID = FALSE)

sites_ok  <- complete.cases(sites_raw)
bg_ok     <- complete.cases(bg_raw)
sites_raw <- sites_raw[sites_ok, ]
bg_raw    <- bg_raw[bg_ok, ]

N_PRES <- nrow(sites_raw)
N_BG   <- nrow(bg_raw)
cat(sprintf("  Sites: %d   Background: %d\n\n", N_PRES, N_BG))

# ── TWO DATA VERSIONS ────────────────────────────────────────
# Version A (XGBoost): all numeric, cats as integer
# Version B (RF):      cats as factor

make_xgb_matrix <- function(df) {
  d <- df
  for (col in cat_in_stack) {
    d[[col]] <- as.numeric(as.integer(d[[col]]))
  }
  m <- as.matrix(d[, final_names, drop = FALSE])
  storage.mode(m) <- "double"
  colnames(m) <- final_names
  m
}

make_rf_df <- function(df) {
  d <- df
  for (col in cat_in_stack) {
    d[[col]] <- factor(as.integer(d[[col]]),
                       levels = raster_levels[[col]])
  }
  d[, final_names, drop = FALSE]
}

sites_xgb <- make_xgb_matrix(sites_raw)
bg_xgb    <- make_xgb_matrix(bg_raw)
all_xgb   <- rbind(sites_xgb, bg_xgb)
all_label <- c(rep(1L, N_PRES), rep(0L, N_BG))

sites_rf  <- make_rf_df(sites_raw)
bg_rf     <- make_rf_df(bg_raw)
all_rf    <- rbind(sites_rf, bg_rf)

# Subsample for RF SHAP (500 rows — treeshap fast enough)
set.seed(42)
# Ensure all sites included + random bg
pres_idx <- seq_len(N_PRES)
bg_sub_idx <- N_PRES + sample(N_BG, min(312L, N_BG))
sub_idx  <- c(pres_idx, bg_sub_idx)  # 188 + 312 = 500

X_rf_sub  <- all_rf[sub_idx, , drop = FALSE]
X_xgb_sub <- all_xgb[sub_idx, , drop = FALSE]
# Numeric version for treeshap
X_rf_num_sub <- as.data.frame(X_xgb_sub)  # all numeric, same as xgb

cat(sprintf("  XGBoost matrix: %d x %d\n",
            nrow(all_xgb), ncol(all_xgb)))
cat(sprintf("  RF subsample:   %d rows (all sites + 312 bg)\n\n",
            nrow(X_rf_sub)))

# ── 5. XGBOOST SHAP ON FULL TRAINING DATA ────────────────────

cat("--- XGBoost TreeSHAP (full training data) ---\n\n")

t0 <- proc.time()
dmat_all <- xgboost::xgb.DMatrix(all_xgb, label = all_label)
sv_xgb   <- shapviz::shapviz(xgb_model,
                             X_pred = dmat_all,
                             X      = as.data.frame(all_xgb))
S_xgb    <- get_shap_matrix(sv_xgb)
colnames(S_xgb) <- final_names

cat(sprintf("  SHAP matrix: %d x %d\n",
            nrow(S_xgb), ncol(S_xgb)))
cat(sprintf("  Done in %.1f min\n\n",
            (proc.time()-t0)[3]/60))

# Global importance: mean |SHAP| across all rows
imp_xgb <- colMeans(abs(S_xgb))
imp_xgb_ranked <- sort(imp_xgb, decreasing = TRUE)

cat("  XGBoost mean |SHAP| ranking:\n")
for (i in seq_along(imp_xgb_ranked)) {
  cat(sprintf("  %2d. %-22s %.4f\n", i,
              names(imp_xgb_ranked)[i],
              imp_xgb_ranked[i]))
}
cat("\n")

# ── 6. RF SHAP VIA TREESHAP (subsample) ──────────────────────

cat("--- RF TreeSHAP (500-row subsample via treeshap) ---\n\n")

RF_SHAP_SUCCESS <- FALSE

sv_rf <- tryCatch({
  
  cat("  Unifying RF model for treeshap...\n")
  # treeshap needs numeric data (no factors)
  unified_rf <- treeshap::randomForest.unify(
    rf_model, X_rf_num_sub)
  
  cat(sprintf("  Computing TreeSHAP on %d rows...\n",
              nrow(X_rf_num_sub)))
  t0 <- proc.time()
  shap_rf_result <- treeshap::treeshap(
    unified_rf, X_rf_num_sub, verbose = FALSE)
  cat(sprintf("  Done in %.1f min\n\n",
              (proc.time()-t0)[3]/60))
  
  sv_obj <- shapviz::shapviz(shap_rf_result)
  RF_SHAP_SUCCESS <- TRUE
  sv_obj
  
}, error = function(e) {
  
  cat("  treeshap failed:", e$message, "\n")
  cat("  Trying kernelshap fallback (200 rows)...\n\n")
  
  if (!requireNamespace("kernelshap", quietly=TRUE)) {
    install.packages("kernelshap")
  }
  library(kernelshap)
  
  # Tiny subsample for kernelshap
  set.seed(42)
  tiny_idx <- sample(nrow(X_rf_sub), min(200L, nrow(X_rf_sub)))
  X_tiny   <- X_rf_sub[tiny_idx, , drop=FALSE]
  bg_tiny  <- X_rf_sub[sample(nrow(X_rf_sub), 50L), , drop=FALSE]
  
  rf_pred_fn <- function(model, newdata) {
    nd <- newdata
    for (col in cat_in_stack) {
      nd[[col]] <- factor(as.integer(nd[[col]]),
                          levels = raster_levels[[col]])
    }
    as.numeric(predict(model, nd, type="prob")[,"1"])
  }
  
  t0 <- proc.time()
  ks_rf <- kernelshap::kernelshap(
    object   = rf_model,
    X        = X_tiny,
    bg_X     = bg_tiny,
    pred_fun = rf_pred_fn,
    verbose  = FALSE
  )
  cat(sprintf("  kernelshap done in %.1f min\n\n",
              (proc.time()-t0)[3]/60))
  
  RF_SHAP_SUCCESS <<- TRUE
  shapviz::shapviz(ks_rf)
})

if (RF_SHAP_SUCCESS) {
  S_rf <- get_shap_matrix(sv_rf)
  colnames(S_rf) <- final_names[seq_len(ncol(S_rf))]
  imp_rf <- colMeans(abs(S_rf))
  imp_rf_ranked <- sort(imp_rf, decreasing = TRUE)
  cat("  RF mean |SHAP| ranking:\n")
  for (i in seq_along(imp_rf_ranked)) {
    cat(sprintf("  %2d. %-22s %.4f\n", i,
                names(imp_rf_ranked)[i], imp_rf_ranked[i]))
  }
  cat("\n")
} else {
  cat("  ⚠ RF SHAP unavailable — using MDA from Script 12\n\n")
  imp_df_rf <- read.csv(
    file.path(OUT_EVAL, "rf_variable_importance.csv"),
    stringsAsFactors = FALSE)
  imp_rf <- setNames(imp_df_rf$MeanDecreaseAccuracy,
                     imp_df_rf$predictor)
  S_rf <- NULL
}

# ── 7. TABLE 4 — GLOBAL SHAP IMPORTANCE ──────────────────────

cat("--- Table 4: Global SHAP Importance ---\n\n")

# Align by predictor name
pred_order <- names(imp_xgb_ranked)  # XGBoost rank order

table4 <- data.frame(
  predictor          = pred_order,
  xgb_mean_abs_shap  = round(imp_xgb[pred_order], 4),
  xgb_rank           = seq_along(pred_order),
  rf_importance      = round(imp_rf[pred_order], 4),
  rf_importance_type = if (RF_SHAP_SUCCESS) "mean_abs_SHAP"
  else "MeanDecreaseAccuracy",
  stringsAsFactors   = FALSE
)
# Add rf rank
table4$rf_rank <- rank(-table4$rf_importance)

write.csv(table4,
          file.path(OUT_TABLES, "Table4_SHAP_importance.csv"),
          row.names = FALSE)
cat("  ✓ Table4_SHAP_importance.csv\n\n")
print(table4[, c("predictor","xgb_rank","xgb_mean_abs_shap",
                 "rf_rank","rf_importance")])
cat("\n")

# ── 8. FIGURE 9 — SHAP BEESWARM ──────────────────────────────

cat("--- Figure 9: SHAP Beeswarm ---\n\n")

tryCatch({
  if (RF_SHAP_SUCCESS && !is.null(S_rf)) {
    png(file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
        width=4800, height=2800, res=300)
    par(mfrow=c(1,2), mar=c(4,12,3,1))
    
    # XGBoost beeswarm (manual — no ggplot dependency)
    # Order by mean |SHAP|
    ord <- order(imp_xgb)
    boxplot(as.data.frame(S_xgb[, ord]),
            horizontal = TRUE,
            las        = 2,
            col        = viridisLite::viridis(ncol(S_xgb)),
            main       = "XGBoost SHAP Values\n(full training data)",
            xlab       = "SHAP value",
            cex.axis   = 0.7,
            outline    = FALSE)
    abline(v=0, lty=2, col="grey40")
    
    # RF beeswarm
    ord_rf <- order(imp_rf[colnames(S_rf)])
    boxplot(as.data.frame(S_rf[, ord_rf]),
            horizontal = TRUE,
            las        = 2,
            col        = viridisLite::viridis(ncol(S_rf)),
            main       = "RF SHAP Values\n(training subsample)",
            xlab       = "SHAP value",
            cex.axis   = 0.7,
            outline    = FALSE)
    abline(v=0, lty=2, col="grey40")
    
    dev.off()
    cat("  ✓ Fig09_SHAP_beeswarm.png (2-panel)\n\n")
    
  } else {
    # XGBoost only
    png(file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
        width=2400, height=2800, res=300)
    par(mar=c(4,12,3,1))
    ord <- order(imp_xgb)
    boxplot(as.data.frame(S_xgb[, ord]),
            horizontal=TRUE, las=2,
            col  = viridisLite::viridis(ncol(S_xgb)),
            main = "XGBoost SHAP Values (full training data)",
            xlab = "SHAP value",
            cex.axis=0.7, outline=FALSE)
    abline(v=0, lty=2, col="grey40")
    dev.off()
    cat("  ✓ Fig09_SHAP_beeswarm.png (XGBoost only)\n\n")
  }
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig09 error:", e$message, "\n\n")
})

# ── 9. FIGURE 10 — SHAP DEPENDENCE PLOTS ─────────────────────

cat("--- Figure 10: SHAP Dependence Plots (top 5) ---\n\n")

top5_preds <- names(imp_xgb_ranked)[1:5]
cat("  Top 5 predictors:", paste(top5_preds, collapse=", "), "\n\n")

tryCatch({
  png(file.path(OUT_FIG_MAIN,"Fig10_SHAP_dependence.png"),
      width=4800, height=4000, res=300)
  par(mfrow=c(2,3), mar=c(4,4,3,1))
  
  for (pred in top5_preds) {
    pred_vals <- all_xgb[, pred]
    shap_vals <- S_xgb[, pred]
    
    # Color by second highest interacting predictor
    # (predictor with highest |SHAP| correlation with residuals)
    shap_resid <- S_xgb[, pred] - mean(S_xgb[, pred])
    # Find most correlated OTHER predictor's values
    cors <- sapply(setdiff(final_names, pred), function(p) {
      tryCatch(abs(cor(all_xgb[,p], shap_resid, use="complete.obs")),
               error=function(e) 0)
    })
    color_pred <- setdiff(final_names, pred)[which.max(cors)]
    color_vals <- all_xgb[, color_pred]
    
    # Scale color_vals to 0-1 for colour gradient
    cv_range <- range(color_vals, na.rm=TRUE)
    cv_scaled <- (color_vals - cv_range[1]) /
      max(cv_range[2] - cv_range[1], 1e-10)
    pt_cols <- viridisLite::viridis(100)[
      pmax(1, pmin(100, round(cv_scaled * 99) + 1))]
    
    plot(pred_vals, shap_vals,
         pch  = 16, cex = 0.3,
         col  = pt_cols,
         xlab = pred,
         ylab = paste("SHAP value for", pred),
         main = sprintf("%s\n(coloured by %s)", pred, color_pred),
         cex.main = 0.75)
    abline(h=0, lty=2, col="grey60")
    # Smooth trend
    tryCatch({
      lo <- loess(shap_vals ~ pred_vals, span=0.5)
      xseq <- seq(min(pred_vals,na.rm=TRUE),
                  max(pred_vals,na.rm=TRUE), length=100)
      lines(xseq, predict(lo, xseq), col="red", lwd=2)
    }, error=function(e) NULL)
  }
  
  # 6th panel: global importance bar
  par(mar=c(4,10,3,1))
  barplot(rev(imp_xgb_ranked),
          horiz   = TRUE,
          las     = 2,
          col     = PRED_COLORS[names(imp_xgb_ranked)],
          main    = "Mean |SHAP| (XGBoost)",
          xlab    = "Mean |SHAP value|",
          cex.names=0.65)
  
  dev.off()
  cat("  ✓ Fig10_SHAP_dependence.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig10 error:", e$message, "\n\n")
})

# ── 10. COHEN'S KAPPA ─────────────────────────────────────────

cat("--- Cohen's Kappa: RF vs XGBoost Dominant Drivers ---\n\n")

if (RF_SHAP_SUCCESS && !is.null(S_rf)) {
  
  # XGBoost dominant driver at subsample locations
  S_xgb_sub <- S_xgb[sub_idx, , drop=FALSE]
  dom_xgb_sub <- colnames(S_xgb_sub)[
    apply(abs(S_xgb_sub), 1, which.max)]
  
  # RF dominant driver at same subsample
  dom_rf_sub <- colnames(S_rf)[
    apply(abs(S_rf), 1, which.max)]
  
  kappa_val <- compute_kappa(dom_xgb_sub, dom_rf_sub)
  
  cat(sprintf("  Agreement table (top predictors):\n"))
  agree_tab <- table(XGBoost=dom_xgb_sub, RF=dom_rf_sub)
  print(agree_tab)
  
  cat(sprintf("\n  Cohen's kappa: %.4f\n", kappa_val))
  
  if (!is.na(kappa_val) && kappa_val > 0.65) {
    KAPPA_DECISION <- "AVERAGE"
    cat("  kappa > 0.65 → maps are concordant\n")
    cat("  → Discuss averaged dominant driver interpretation\n\n")
  } else {
    KAPPA_DECISION <- "SEPARATE"
    cat("  kappa <= 0.65 → maps discordant\n")
    cat("  → Report RF and XGBoost maps separately\n\n")
  }
  
  saveRDS(list(kappa_val=kappa_val,
               decision=KAPPA_DECISION,
               dom_xgb=dom_xgb_sub,
               dom_rf =dom_rf_sub),
          file.path(OUT_SHAP, "cohen_kappa_result.rds"))
  
} else {
  kappa_val   <- NA_real_
  KAPPA_DECISION <- "NOT_COMPUTED"
  cat("  RF SHAP unavailable — kappa not computed\n\n")
}

# ── 11. 250m PREDICTION GRID ──────────────────────────────────

cat("--- Setting Up 250m Prediction Grid ---\n\n")

# Create exact 250m template
template_250m <- terra::rast(
  ext = terra::ext(template_30m),
  res = 250,
  crs = terra::crs(template_30m)
)

# Resample predictor stack to 250m (bilinear for continuous)
cat("  Resampling predictor stack to 250m...\n")
pred_250m <- terra::resample(pred_stack, template_250m,
                             method="bilinear")
pred_250m <- terra::mask(pred_250m, boundary_vect)

# Extract values with cell IDs for mapping back
grid_df <- terra::as.data.frame(pred_250m, xy=TRUE,
                                cells=TRUE, na.rm=TRUE)
cell_ids <- grid_df$cell
n_cells  <- nrow(grid_df)

cat(sprintf("  250m grid: %d valid cells\n", n_cells))

# Prepare XGBoost matrix from 250m grid
grid_mat <- as.matrix(grid_df[, final_names, drop=FALSE])
# Ensure all numeric (cats already integer in raster)
for (col in cat_in_stack) {
  grid_mat[, col] <- as.numeric(as.integer(grid_mat[, col]))
}
storage.mode(grid_mat) <- "double"
colnames(grid_mat) <- final_names

cat(sprintf("  Grid matrix: %d x %d\n\n",
            nrow(grid_mat), ncol(grid_mat)))

# ── 12. XGBOOST SHAP ON 250m GRID (chunked) ──────────────────

cat("--- XGBoost SHAP on 250m Grid (chunked) ---\n\n")

CHUNK_SIZE <- 50000L
n_chunks   <- ceiling(n_cells / CHUNK_SIZE)
cat(sprintf("  Processing %d cells in %d chunks of %d...\n\n",
            n_cells, n_chunks, CHUNK_SIZE))

S_grid_list <- vector("list", n_chunks)
t0_grid     <- proc.time()

for (ch in seq_len(n_chunks)) {
  idx_s <- (ch - 1L) * CHUNK_SIZE + 1L
  idx_e <- min(ch * CHUNK_SIZE, n_cells)
  chunk_mat  <- grid_mat[idx_s:idx_e, , drop=FALSE]
  chunk_dmat <- xgboost::xgb.DMatrix(chunk_mat)
  
  sv_chunk <- shapviz::shapviz(
    xgb_model,
    X_pred = chunk_dmat,
    X      = as.data.frame(chunk_mat))
  
  S_grid_list[[ch]] <- get_shap_matrix(sv_chunk)
  rm(sv_chunk, chunk_dmat); gc(full=TRUE)
  
  cat(sprintf("  Chunk %d/%d: cells %d-%d  (%.1f min elapsed)\n",
              ch, n_chunks, idx_s, idx_e,
              (proc.time()-t0_grid)[3]/60))
}

S_grid <- do.call(rbind, S_grid_list)
colnames(S_grid) <- final_names
rm(S_grid_list); gc(full=TRUE)

cat(sprintf("\n  SHAP grid complete: %d x %d  (%.1f min total)\n\n",
            nrow(S_grid), ncol(S_grid),
            (proc.time()-t0_grid)[3]/60))

# ── 13. DOMINANT DRIVER MAP ───────────────────────────────────

cat("--- Computing Dominant Driver per Cell ---\n\n")

# Dominant driver = predictor with highest |SHAP| per cell
dom_idx  <- apply(abs(S_grid), 1, which.max)
dom_name <- final_names[dom_idx]
dom_code <- as.integer(dom_idx)

# Dominant driver frequency
dom_freq <- sort(table(dom_name), decreasing=TRUE)
cat("  Dominant driver distribution (250m cells):\n")
total_cells <- length(dom_name)
for (pred in names(dom_freq)) {
  pct <- 100 * dom_freq[pred] / total_cells
  cat(sprintf("  %-22s %7d cells  %5.1f%%\n",
              pred, dom_freq[pred], pct))
}
cat("\n")

# Create dominant driver raster
dom_rast <- template_250m
terra::values(dom_rast) <- NA_integer_
dom_rast[cell_ids] <- dom_code
dom_rast <- terra::mask(dom_rast, boundary_vect)

# Save as integer raster
terra::writeRaster(
  dom_rast,
  file.path(OUT_SHAP, "shap_dominant_driver_map.tif"),
  overwrite=TRUE, datatype="INT1U")
cat("  ✓ shap_dominant_driver_map.tif\n\n")

# Save lookup table (code → predictor name)
driver_lookup <- data.frame(
  code  = seq_along(final_names),
  predictor = final_names,
  color     = PRED_COLORS[final_names],
  stringsAsFactors=FALSE)
write.csv(driver_lookup,
          file.path(OUT_SHAP, "dominant_driver_lookup.csv"),
          row.names=FALSE)
cat("  ✓ dominant_driver_lookup.csv\n\n")

# ── 14. LANDSCAPE ZONES ───────────────────────────────────────

cat("--- Landscape Zone Classification ---\n\n")
cat("  Method: DEM elevation thresholds\n")
cat("  Zone 1 — Nagpur pediplain:   270–380 m\n")
cat("  Zone 2 — Chandrapur basin:   < 270 m\n")
cat("  Zone 3 — Satpura foothills:  > 380 m\n\n")

dem_r    <- terra::rast(file.path(OUT_PREDICTORS,
                                  "DEM_30m_utm44n.tif"))
dem_250m <- terra::resample(dem_r, template_250m,
                            method="bilinear")
dem_250m <- terra::mask(dem_250m, boundary_vect)

dem_stats <- terra::global(dem_250m, c("min","max","mean","sd"),
                           na.rm=TRUE)
cat(sprintf("  DEM stats (250m): min=%.0f max=%.0f mean=%.0f SD=%.0f\n\n",
            dem_stats[1,1], dem_stats[1,2],
            dem_stats[1,3], dem_stats[1,4]))

# NOTE: Adjust thresholds 270 and 380 based on DEM stats above
zone_rast <- terra::ifel(
  dem_250m > 380, 3L,
  terra::ifel(dem_250m > 270, 1L, 2L))
zone_rast <- terra::mask(zone_rast, boundary_vect)

zone_names <- c("1"="Nagpur pediplain",
                "2"="Chandrapur basin",
                "3"="Satpura foothills")

# Extract zone values at grid cells
zone_vals <- terra::extract(zone_rast, 
                            terra::xyFromCell(template_250m, cell_ids),
                            ID=FALSE)[,1]

zone_freq_tbl <- table(zone=zone_vals)
cat("  Zone distribution (250m cells):\n")
for (z in names(zone_freq_tbl)) {
  cat(sprintf("  Zone %s %-20s %7d cells %5.1f%%\n",
              z, zone_names[z], zone_freq_tbl[z],
              100*zone_freq_tbl[z]/sum(zone_freq_tbl)))
}
cat("\n")

# ── 15. CHI-SQUARED TEST ──────────────────────────────────────

cat("--- Chi-Squared Test: Spatial Heterogeneity ---\n\n")

# Use cells where both zone and dominant driver are available
valid_chi <- !is.na(zone_vals) & !is.na(dom_name)
chi_table  <- table(
  Zone   = zone_vals[valid_chi],
  Driver = dom_name[valid_chi]
)

cat("  Contingency table (zones x dominant drivers):\n")
print(chi_table)
cat("\n")

# Check expected cell frequencies
n_chi    <- sum(chi_table)
n_zone   <- rowSums(chi_table)
n_driver <- colSums(chi_table)
min_expected <- min(outer(n_zone, n_driver, "*")) / n_chi
cat(sprintf("  Min expected frequency: %.2f\n", min_expected))

if (min_expected < 5) {
  cat("  ⚠ Some expected frequencies < 5 — using simulation\n\n")
  chi_result <- chisq.test(chi_table,
                           simulate.p.value = TRUE, B = 10000)
} else {
  chi_result <- chisq.test(chi_table)
}

cat(sprintf("  Chi-squared: %.2f\n", chi_result$statistic))
cat(sprintf("  df:          %d\n",   chi_result$parameter))
cat(sprintf("  p-value:     %.6f\n", chi_result$p.value))

if (chi_result$p.value < 0.001) {
  cat("  ✓ CONCLUSION: Dominant driver distribution differs\n")
  cat("    significantly across landscape zones (p < 0.001)\n\n")
} else if (chi_result$p.value < 0.05) {
  cat("  ✓ CONCLUSION: Significant spatial heterogeneity\n\n")
} else {
  cat("  NS: No significant difference across zones\n\n")
}

# Zone-specific dominant driver (which predictor tops each zone)
cat("  Top dominant driver per zone:\n")
zone_driver_df <- data.frame(
  zone        = character(),
  zone_label  = character(),
  top_driver  = character(),
  pct_cells   = numeric(),
  stringsAsFactors=FALSE)

for (z in rownames(chi_table)) {
  zone_row <- chi_table[z, ]
  top_drv  <- names(which.max(zone_row))
  top_pct  <- 100 * max(zone_row) / sum(zone_row)
  cat(sprintf("  Zone %s %-20s → %-22s %.1f%% of zone cells\n",
              z, zone_names[z], top_drv, top_pct))
  zone_driver_df <- rbind(zone_driver_df, data.frame(
    zone=z, zone_label=zone_names[z],
    top_driver=top_drv, pct_cells=round(top_pct,1),
    stringsAsFactors=FALSE))
}
cat("\n")

# Save chi-squared results
chi_summary <- data.frame(
  chi_statistic  = round(chi_result$statistic, 3),
  df             = chi_result$parameter,
  p_value        = round(chi_result$p.value, 6),
  simulated_p    = (min_expected < 5),
  n_cells        = n_chi,
  conclusion     = if(chi_result$p.value<0.05)
    "Significant spatial heterogeneity" else "NS",
  stringsAsFactors=FALSE)

write.csv(chi_summary,
          file.path(OUT_SHAP,"shap_chi_squared_result.csv"),
          row.names=FALSE)
write.csv(zone_driver_df,
          file.path(OUT_SHAP,"shap_zone_dominant_drivers.csv"),
          row.names=FALSE)
cat("  ✓ shap_chi_squared_result.csv\n")
cat("  ✓ shap_zone_dominant_drivers.csv\n\n")

# ── 16. FIGURE 11 — DOMINANT DRIVER MAP ──────────────────────

cat("--- Figure 11: Dominant Driver Map (Central Figure) ---\n\n")

tryCatch({
  # Pre-compute zone boundaries for overlay
  zone_contours <- tryCatch(
    terra::as.contour(zone_rast, nlevels=3),
    error=function(e) NULL)
  
  # Subset colours to predictors actually dominant somewhere
  active_drivers <- names(dom_freq)[dom_freq > 0]
  active_cols    <- PRED_COLORS[active_drivers]
  active_codes   <- match(active_drivers, final_names)
  
  png(file.path(OUT_FIG_MAIN,"Fig11_dominant_driver_map.png"),
      width=3600, height=3600, res=300)
  
  # Reclassify to sequential codes for plotting (only active drivers)
  dom_plot <- terra::classify(
    dom_rast,
    cbind(active_codes,
          seq_along(active_codes)))
  
  terra::plot(dom_plot,
              main   = sprintf(
                "SHAP Dominant Predictor per Cell (250m grid)\nChi²=%.1f, p=%.4f",
                chi_result$statistic, chi_result$p.value),
              col    = active_cols,
              type   = "classes",
              axes   = FALSE,
              legend = FALSE,
              cex.main = 0.80)
  
  terra::plot(boundary_vect, add=TRUE,
              border="black", lwd=1.0)
  
  # Zone boundary overlay (white dashed lines)
  if (!is.null(zone_contours) && length(zone_contours) > 0) {
    terra::lines(zone_contours, col="white",
                 lwd=1.2, lty=2)
  }
  
  # Site locations
  terra::plot(terra::vect(sites_sf), add=TRUE,
              col="white", pch=16, cex=0.25)
  
  # Zone labels (approximate centres)
  text(x=c(740000, 810000, 720000),
       y=c(2360000, 2240000, 2200000),
       labels=c("Nagpur\npediplain",
                "Chandrapur\nbasin",
                "Satpura\nfoothills"),
       col="white", cex=0.65, font=2)
  
  # Legend — active drivers with % coverage
  legend_labels <- sprintf(
    "%s (%.1f%%)",
    active_drivers,
    100 * dom_freq[active_drivers] / total_cells)
  
  legend("bottomright",
         legend = legend_labels,
         fill   = active_cols,
         title  = "Dominant Predictor",
         bty    = "n",
         cex    = 0.52,
         ncol   = if(length(active_drivers)>8) 2L else 1L)
  
  dev.off()
  cat("  ✓ Fig11_dominant_driver_map.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig11 error:", e$message, "\n\n")
})

# ── 17. SUPPLEMENTARY FIGURE — ZONE × DRIVER HEATMAP ─────────

cat("--- Supplementary: Zone x Driver Heatmap ---\n\n")

tryCatch({
  chi_prop <- prop.table(chi_table, margin=1) * 100
  
  png(file.path(OUT_FIG_SUPP,
                "FigS_shap_zone_driver_heatmap.png"),
      width=3600, height=1600, res=300)
  
  n_drivers <- ncol(chi_prop)
  n_zones   <- nrow(chi_prop)
  
  image(t(chi_prop),
        col    = viridisLite::plasma(100),
        axes   = FALSE,
        main   = "% Cells with Dominant Driver by Zone",
        xlab   = "Driver", ylab="Zone")
  
  axis(1, at=seq(0,1,length.out=n_drivers),
       labels=colnames(chi_prop), las=2, cex.axis=0.55)
  axis(2, at=seq(0,1,length.out=n_zones),
       labels=paste0("Z",rownames(chi_prop),
                     " ",zone_names[rownames(chi_prop)]),
       las=2, cex.axis=0.65)
  
  # Cell values
  for (i in seq_len(n_zones)) {
    for (j in seq_len(n_drivers)) {
      val <- chi_prop[i,j]
      if (val > 1) {
        text((j-1)/(n_drivers-1),
             (i-1)/(max(n_zones-1,1)),
             sprintf("%.1f", val),
             cex=0.45, col="white")
      }
    }
  }
  dev.off()
  cat("  ✓ FigS_shap_zone_driver_heatmap.png\n\n")
}, error=function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Heatmap error:", e$message, "\n\n")
})

# ── 18. SAVE ALL SHAP DATA ────────────────────────────────────

cat("--- Saving SHAP Data ---\n\n")

# Global importance CSV
imp_csv <- data.frame(
  predictor         = names(imp_xgb_ranked),
  xgb_mean_abs_shap = round(imp_xgb_ranked, 4),
  xgb_rank          = seq_along(imp_xgb_ranked),
  stringsAsFactors  = FALSE
)
if (RF_SHAP_SUCCESS) {
  imp_csv$rf_mean_abs_shap <- round(imp_rf[imp_csv$predictor], 4)
} else {
  imp_csv$rf_mda <- round(imp_rf[imp_csv$predictor], 4)
}
write.csv(imp_csv,
          file.path(OUT_SHAP,"shap_global_importance.csv"),
          row.names=FALSE)
cat("  ✓ shap_global_importance.csv\n")

# Dominant driver frequency
dom_freq_csv <- data.frame(
  predictor  = names(dom_freq),
  n_cells    = as.integer(dom_freq),
  pct_cells  = round(100*as.numeric(dom_freq)/total_cells, 2),
  stringsAsFactors=FALSE
)
write.csv(dom_freq_csv,
          file.path(OUT_SHAP,"dominant_driver_frequency.csv"),
          row.names=FALSE)
cat("  ✓ dominant_driver_frequency.csv\n\n")

# ── 19. SUMMARY ──────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 19 COMPLETE — SHAP Analysis\n")
cat("========================================\n")
cat(sprintf("RF SHAP:     %s\n",
            if(RF_SHAP_SUCCESS) "✓ computed (subsample)" else
              "⚠ MDA fallback"))
cat(sprintf("Cohen kappa: %.4f (%s)\n",
            kappa_val, KAPPA_DECISION))
cat(sprintf("Chi-squared: %.2f  p=%.6f\n",
            chi_result$statistic, chi_result$p.value))
cat(sprintf("Grid cells:  %d (250m resolution)\n", n_cells))
cat(sprintf("Active dominant drivers: %d/%d predictors\n",
            length(active_drivers), length(final_names)))
cat("\nTop dominant driver per zone:\n")
for (i in seq_len(nrow(zone_driver_df))) {
  cat(sprintf("  Zone %s %-20s → %s (%.1f%%)\n",
              zone_driver_df$zone[i],
              zone_driver_df$zone_label[i],
              zone_driver_df$top_driver[i],
              zone_driver_df$pct_cells[i]))
}
cat("\nXGBoost SHAP importance (top 5):\n")
for (i in 1:min(5,length(imp_xgb_ranked))) {
  cat(sprintf("  %d. %-22s %.4f\n", i,
              names(imp_xgb_ranked)[i],
              imp_xgb_ranked[i]))
}
cat("\nFiles saved to outputs/shap/:\n")
cat("  Table4_SHAP_importance.csv\n")
cat("  shap_global_importance.csv\n")
cat("  shap_dominant_driver_map.tif\n")
cat("  dominant_driver_lookup.csv\n")
cat("  dominant_driver_frequency.csv\n")
cat("  cohen_kappa_result.rds\n")
cat("  shap_chi_squared_result.csv\n")
cat("  shap_zone_dominant_drivers.csv\n")
cat("  Fig09_SHAP_beeswarm.png\n")
cat("  Fig10_SHAP_dependence.png\n")
cat("  Fig11_dominant_driver_map.png\n")
cat("  FigS_shap_zone_driver_heatmap.png\n")
cat("\nNOTE: Verify DEM zone thresholds (270m, 380m)\n")
cat("against study area DEM statistics printed above.\n")
cat("Adjust if needed before manuscript writing.\n")
cat("\nNext: Script 20 — Cultural Period Sub-Models (LP/MP/UP)\n")
cat("========================================\n")