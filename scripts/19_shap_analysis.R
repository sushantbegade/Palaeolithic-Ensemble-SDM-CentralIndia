# ============================================================
# SCRIPT 19: SHAP ANALYSIS (FINAL — xgboost 3.2.x compatible)
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 19 of 25
# ============================================================
# XGBOOST 3.2.x FIX:
#   xgb_model$feature_names <- ... BLOCKED (ALTREP object).
#   Solution: predict(model, dmat, predcontrib=TRUE) returns
#   exact TreeSHAP matrix directly (n × p+1; drop last col
#   = bias term). Create shapviz from matrix, not model.
#   Feature names assigned from xgb_info$feature_names.
#
# SHAP SPATIAL GRID: 250m (~300k cells), chunked 50k/chunk.
# RF SHAP: treeshap 0.4.0 on 500-row subsample.
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

# ── 1. PACKAGES ───────────────────────────────────────────────

cat("--- Package Check ---\n")
for (pkg in c("treeshap","shapviz")) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    cat("  Installing", pkg, "...\n")
    install.packages(pkg)
  }
  library(pkg, character.only=TRUE)
  cat(sprintf("  %-10s %s ✓\n", pkg,
              as.character(packageVersion(pkg))))
}
cat("\n")

# ── 2. HELPERS ────────────────────────────────────────────────

# Cohen's kappa — no extra package
compute_kappa <- function(x, y) {
  lvls <- union(unique(x), unique(y))
  x <- factor(x, levels=lvls); y <- factor(y, levels=lvls)
  n <- length(x); k <- length(lvls)
  conf <- table(x, y)
  Po <- sum(diag(conf)) / n
  Pe <- sum((rowSums(conf)/n) * (colSums(conf)/n))
  if (abs(1-Pe) < 1e-10) return(NA_real_)
  (Po - Pe) / (1 - Pe)
}

# 13-predictor colour palette
PRED_COLORS <- c(
  "Elevation"          = "#E63946",
  "Aspect"             = "#F4A261",
  "TRI"                = "#2A9D8F",
  "TPI"                = "#264653",
  "Plan_Curvature"     = "#E9C46A",
  "HAND"               = "#06D6A0",
  "Flow_Accum_log10"   = "#118AB2",
  "Dist_River"         = "#073B4C",
  "Dist_Palaeochannel" = "#8338EC",
  "Dist_RawMat"        = "#FF6B6B",
  "NDVI"               = "#56CFE1",
  "Geology"            = "#FF9F1C",
  "Geomorphology"      = "#4CC9F0"
)

# ── 3. LOAD MODELS + METADATA ─────────────────────────────────

cat("--- Loading Models ---\n\n")

# XGBoost — load model (warning about format = non-fatal)
xgb_model  <- suppressWarnings(
  xgboost::xgb.load(file.path(OUT_MOD_IND,
                              "xgboost_model_final.bin")))
xgb_info   <- readRDS(file.path(OUT_MOD_IND,
                                "xgboost_model_info.rds"))
feat_names <- xgb_info$feature_names   # use independently
# DO NOT assign to xgb_model$feature_names (ALTREP error in v3.2)

cat("  XGBoost loaded (v3.2.x compatible)\n")
cat("  nrounds:", xgb_info$nrounds, "\n")
cat("  Features:", paste(feat_names, collapse=", "), "\n\n")

# RF
rf_model <- readRDS(file.path(OUT_MOD_IND, "rf_model_final.rds"))
cat("  RF loaded. ntree =", rf_model$ntree, "\n\n")

# Predictor metadata
final_names   <- readRDS(file.path(OUT_PREDICTORS,
                                   "final_predictor_names.rds"))
raster_levels <- readRDS(file.path(OUT_MOD_IND,
                                   "gam_raster_levels.rds"))
cat_predictors <- c("Geology","Geomorphology")
cat_in_stack   <- cat_predictors[cat_predictors %in% final_names]

pred_stack    <- terra::rast(file.path(OUT_PREDICTORS,
                                       "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES,"study_area_boundary_utm44n.gpkg"),
  quiet=TRUE))

# ── 4. EXTRACT TRAINING DATA ──────────────────────────────────

cat("--- Extracting Training Data ---\n\n")

sites_sf <- sf::st_read(file.path(OUT_CV,
                                  "sites_with_folds.gpkg"),
                        quiet=TRUE)
bg_sf    <- sf::st_read(file.path(OUT_BACKGROUND,
                                  "background_N10000.gpkg"),
                        quiet=TRUE)

sites_raw <- terra::extract(pred_stack,
                            terra::vect(sites_sf), ID=FALSE)
bg_raw    <- terra::extract(pred_stack,
                            terra::vect(bg_sf),    ID=FALSE)

sites_ok <- complete.cases(sites_raw)
bg_ok    <- complete.cases(bg_raw)
sites_raw <- sites_raw[sites_ok, ]
bg_raw    <- bg_raw[bg_ok, ]

N_PRES <- nrow(sites_raw); N_BG <- nrow(bg_raw)
cat(sprintf("  Sites: %d   Background: %d\n\n", N_PRES, N_BG))

# XGBoost version: all numeric, cats as integer
make_xgb_mat <- function(df) {
  d <- df[, final_names, drop=FALSE]
  for (col in cat_in_stack)
    d[[col]] <- as.numeric(as.integer(d[[col]]))
  m <- as.matrix(d)
  storage.mode(m) <- "double"
  colnames(m) <- final_names
  m
}

# RF version: cats as factor
make_rf_df <- function(df) {
  d <- df[, final_names, drop=FALSE]
  for (col in cat_in_stack)
    d[[col]] <- factor(as.integer(d[[col]]),
                       levels=raster_levels[[col]])
  d
}

sites_xgb <- make_xgb_mat(sites_raw)
bg_xgb    <- make_xgb_mat(bg_raw)
all_xgb   <- rbind(sites_xgb, bg_xgb)

sites_rf  <- make_rf_df(sites_raw)
bg_rf     <- make_rf_df(bg_raw)

# 500-row subsample for RF SHAP (all sites + 312 bg)
set.seed(42)
sub_bg_idx <- N_PRES + sample(N_BG, min(312L, N_BG))
sub_idx    <- c(seq_len(N_PRES), sub_bg_idx)
X_rf_sub   <- make_rf_df(rbind(sites_raw, bg_raw)[sub_idx, ])
X_xgb_sub  <- all_xgb[sub_idx, , drop=FALSE]

cat(sprintf("  XGBoost matrix:  %d x %d\n",
            nrow(all_xgb), ncol(all_xgb)))
cat(sprintf("  RF subsample:    %d rows\n\n", nrow(X_rf_sub)))

# ── 5. XGBOOST SHAP — predict(predcontrib=TRUE) ───────────────
# FIX: bypasses xgb_model$feature_names assignment (ALTREP error)
# predcontrib=TRUE → exact TreeSHAP, output = n × (p+1) matrix
# Last column = SHAP bias/intercept → DROP IT

cat("--- XGBoost TreeSHAP (predcontrib=TRUE) ---\n\n")

t0 <- proc.time()

dmat_all  <- xgboost::xgb.DMatrix(all_xgb)
shap_raw  <- predict(xgb_model, dmat_all, predcontrib=TRUE)

cat(sprintf("  Raw SHAP output: %d x %d (last col = bias)\n",
            nrow(shap_raw), ncol(shap_raw)))

# Drop bias column (last), assign feature names
S_xgb <- shap_raw[, seq_len(length(feat_names)), drop=FALSE]
colnames(S_xgb) <- feat_names

cat(sprintf("  SHAP matrix:     %d x %d\n",
            nrow(S_xgb), ncol(S_xgb)))
cat(sprintf("  Done in %.1f min\n\n",
            (proc.time()-t0)[3]/60))

# Create shapviz object FROM MATRIX (not from model)
# shapviz(matrix, X=data.frame) is fully supported in v0.10.x
sv_xgb <- shapviz::shapviz(S_xgb,
                           X = as.data.frame(all_xgb))
cat("  shapviz object created from SHAP matrix ✓\n\n")

# Global importance
imp_xgb        <- colMeans(abs(S_xgb))
imp_xgb_ranked <- sort(imp_xgb, decreasing=TRUE)

cat("  XGBoost mean |SHAP| ranking:\n")
for (i in seq_along(imp_xgb_ranked)) {
  cat(sprintf("  %2d. %-22s %.4f\n", i,
              names(imp_xgb_ranked)[i],
              imp_xgb_ranked[i]))
}
cat("\n")

# ── 6. RF SHAP VIA TREESHAP ───────────────────────────────────

cat("--- RF TreeSHAP (treeshap 0.4.0, 500-row subsample) ---\n\n")

RF_SHAP_SUCCESS <- FALSE
S_rf <- NULL
sv_rf <- NULL
imp_rf <- NULL

# treeshap needs numeric data (factors → integer)
X_rf_num_sub <- make_xgb_mat(
  rbind(sites_raw, bg_raw)[sub_idx, ])

sv_rf <- tryCatch({
  
  cat("  Unifying RF for treeshap...\n")
  unified_rf <- treeshap::randomForest.unify(
    rf_model, as.data.frame(X_rf_num_sub))
  
  cat(sprintf("  Computing RF TreeSHAP on %d rows...\n",
              nrow(X_rf_num_sub)))
  t0 <- proc.time()
  shap_rf_obj <- treeshap::treeshap(
    unified_rf,
    as.data.frame(X_rf_num_sub),
    verbose = FALSE)
  cat(sprintf("  Done in %.1f min\n\n",
              (proc.time()-t0)[3]/60))
  
  sv_obj <- shapviz::shapviz(shap_rf_obj)
  RF_SHAP_SUCCESS <<- TRUE
  sv_obj
  
}, error = function(e) {
  
  cat("  treeshap error:", e$message, "\n")
  cat("  Fallback: kernelshap on 150 rows...\n\n")
  
  if (!requireNamespace("kernelshap", quietly=TRUE))
    install.packages("kernelshap")
  library(kernelshap)
  
  set.seed(42)
  tiny_idx <- sample(nrow(X_rf_sub), min(150L, nrow(X_rf_sub)))
  X_tiny   <- X_rf_sub[tiny_idx, , drop=FALSE]
  bg_tiny  <- X_rf_sub[sample(nrow(X_rf_sub), 50L), , drop=FALSE]
  
  rf_pred_fn <- function(model, newdata) {
    nd <- newdata
    for (col in cat_in_stack)
      nd[[col]] <- factor(as.integer(nd[[col]]),
                          levels=raster_levels[[col]])
    as.numeric(predict(model, nd, type="prob")[,"1"])
  }
  
  t0 <- proc.time()
  ks <- kernelshap::kernelshap(
    object=rf_model, X=X_tiny, bg_X=bg_tiny,
    pred_fun=rf_pred_fn, verbose=FALSE)
  cat(sprintf("  kernelshap done in %.1f min\n\n",
              (proc.time()-t0)[3]/60))
  
  RF_SHAP_SUCCESS <<- TRUE
  shapviz::shapviz(ks)
})

if (RF_SHAP_SUCCESS && !is.null(sv_rf)) {
  # Extract SHAP matrix — works for both treeshap and kernelshap
  S_rf <- tryCatch(sv_rf$S, error=function(e) NULL)
  if (is.null(S_rf)) S_rf <- tryCatch(
    as.matrix(sv_rf[["S"]]), error=function(e) NULL)
  
  if (!is.null(S_rf)) {
    colnames(S_rf) <- feat_names[seq_len(ncol(S_rf))]
    imp_rf        <- colMeans(abs(S_rf))
    imp_rf_ranked <- sort(imp_rf, decreasing=TRUE)
    
    cat("  RF mean |SHAP| ranking:\n")
    for (i in seq_along(imp_rf_ranked)) {
      cat(sprintf("  %2d. %-22s %.4f\n", i,
                  names(imp_rf_ranked)[i],
                  imp_rf_ranked[i]))
    }
    cat("\n")
  } else {
    RF_SHAP_SUCCESS <- FALSE
    cat("  ⚠ Could not extract RF SHAP matrix\n\n")
  }
}

# Fallback importance: RF MDA from Script 12
if (!RF_SHAP_SUCCESS || is.null(imp_rf)) {
  cat("  Using RF MDA from Script 12 as fallback\n\n")
  imp_df_rf <- read.csv(
    file.path(OUT_EVAL,"rf_variable_importance.csv"),
    stringsAsFactors=FALSE)
  imp_rf <- setNames(imp_df_rf$MeanDecreaseAccuracy,
                     imp_df_rf$predictor)
}

# ── 7. TABLE 4 ────────────────────────────────────────────────

cat("--- Table 4: Global SHAP Importance ---\n\n")

pred_ord <- names(imp_xgb_ranked)
table4   <- data.frame(
  predictor         = pred_ord,
  xgb_mean_abs_shap = round(imp_xgb[pred_ord], 4),
  xgb_rank          = seq_along(pred_ord),
  rf_importance     = round(imp_rf[pred_ord], 4),
  rf_type           = if (RF_SHAP_SUCCESS && !is.null(S_rf))
    "mean_abs_SHAP" else "MeanDecreaseAccuracy",
  stringsAsFactors  = FALSE
)
table4$rf_rank <- rank(-table4$rf_importance, ties.method="first")

write.csv(table4,
          file.path(OUT_TABLES,"Table4_SHAP_importance.csv"),
          row.names=FALSE)
cat("  ✓ Table4_SHAP_importance.csv\n\n")
print(table4[, c("predictor","xgb_rank","xgb_mean_abs_shap",
                 "rf_rank","rf_importance")])
cat("\n")

# ── 8. FIGURE 9 — BEESWARM ────────────────────────────────────

cat("--- Figure 9: Beeswarm ---\n\n")

tryCatch({
  
  has_rf_sv <- RF_SHAP_SUCCESS && !is.null(sv_rf) &&
    !is.null(S_rf)
  
  p9 <- if (has_rf_sv) {
    # Combined plot: XGBoost + RF via shapviz
    p_xgb <- shapviz::sv_importance(sv_xgb, kind="beeswarm",
                                    max_display=13L) +
      ggplot2::ggtitle("XGBoost SHAP (full training data)") +
      ggplot2::theme(plot.title=ggplot2::element_text(size=9))
    p_rf  <- shapviz::sv_importance(sv_rf,  kind="beeswarm",
                                    max_display=13L) +
      ggplot2::ggtitle("RF SHAP (500-row subsample)") +
      ggplot2::theme(plot.title=ggplot2::element_text(size=9))
    patchwork::wrap_plots(p_xgb, p_rf, ncol=2)
  } else {
    shapviz::sv_importance(sv_xgb, kind="beeswarm",
                           max_display=13L) +
      ggplot2::ggtitle("XGBoost SHAP (full training data)")
  }
  
  ggplot2::ggsave(
    file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
    plot   = p9,
    width  = if (has_rf_sv) 16 else 9,
    height = 8,
    dpi    = 300,
    units  = "in")
  cat("  ✓ Fig09_SHAP_beeswarm.png\n\n")
  
}, error = function(e) {
  cat("  ✗ Fig09 error:", e$message, "\n")
  cat("  Trying base-R fallback...\n")
  
  tryCatch({
    ord <- order(imp_xgb)
    n_p <- ncol(S_xgb)
    # Sample for plotting speed
    set.seed(42)
    row_sub <- sample(nrow(S_xgb), min(3000, nrow(S_xgb)))
    
    png(file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
        width=2400, height=3000, res=300)
    par(mar=c(4,12,3,1))
    boxplot(as.data.frame(S_xgb[row_sub, ord]),
            horizontal=TRUE, las=2, outline=FALSE,
            col=viridisLite::viridis(n_p),
            main="XGBoost SHAP Values\n(full training data)",
            xlab="SHAP value", cex.axis=0.65)
    abline(v=0, lty=2, col="grey40")
    dev.off()
    cat("  ✓ Fig09_SHAP_beeswarm.png (base-R fallback)\n\n")
  }, error=function(e2) {
    tryCatch(dev.off(), error=function(x) NULL)
    cat("  ✗ Fig09 fallback also failed:", e2$message, "\n\n")
  })
})

# ── 9. FIGURE 10 — DEPENDENCE PLOTS ──────────────────────────

cat("--- Figure 10: Dependence Plots (top 5) ---\n\n")

top5 <- names(imp_xgb_ranked)[1:5]
cat("  Top 5:", paste(top5, collapse=", "), "\n\n")

tryCatch({
  
  # Find best interaction predictor for each top-5
  get_interact <- function(focal, S, X) {
    others <- setdiff(colnames(S), focal)
    res_shap <- S[,focal] - mean(S[,focal])
    cors <- sapply(others, function(p)
      tryCatch(abs(cor(X[,p], res_shap, use="complete.obs")),
               error=function(e) 0))
    others[which.max(cors)]
  }
  
  plots10 <- lapply(top5, function(pred) {
    int_pred <- get_interact(pred, S_xgb, all_xgb)
    shapviz::sv_dependence(sv_xgb,
                           v         = pred,
                           color_var = int_pred) +
      ggplot2::theme(plot.title=ggplot2::element_text(size=8))
  })
  
  p10 <- patchwork::wrap_plots(plots10, ncol=3,
                               nrow=2) +
    patchwork::plot_annotation(
      title="SHAP Dependence Plots — Top 5 Predictors (XGBoost)")
  
  ggplot2::ggsave(
    file.path(OUT_FIG_MAIN,"Fig10_SHAP_dependence.png"),
    plot=p10, width=16, height=11,
    dpi=300, units="in")
  cat("  ✓ Fig10_SHAP_dependence.png\n\n")
  
}, error=function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig10 error:", e$message, "\n\n")
})

# ── 10. COHEN'S KAPPA ────────────────────────────────────────

cat("--- Cohen's Kappa (RF vs XGBoost at subsample) ---\n\n")

kappa_val      <- NA_real_
KAPPA_DECISION <- "NOT_COMPUTED"

if (!is.null(S_rf) && RF_SHAP_SUCCESS) {
  
  # XGBoost dominant at same 500 subsample rows
  S_xgb_sub    <- S_xgb[sub_idx, , drop=FALSE]
  dom_xgb_sub  <- feat_names[apply(abs(S_xgb_sub), 1, which.max)]
  dom_rf_sub   <- feat_names[apply(abs(S_rf),       1, which.max)]
  
  kappa_val <- compute_kappa(dom_xgb_sub, dom_rf_sub)
  
  cat("  Dominant driver agreement (top 5 shown):\n")
  agree_tab <- sort(table(XGB=dom_xgb_sub, RF=dom_rf_sub),
                    decreasing=TRUE)
  print(head(as.data.frame(agree_tab), 10))
  
  cat(sprintf("\n  Cohen's kappa: %.4f\n", kappa_val))
  KAPPA_DECISION <- if (!is.na(kappa_val) && kappa_val > 0.65)
    "AVERAGE" else "SEPARATE"
  cat(sprintf("  Decision: %s\n\n", KAPPA_DECISION))
  
  saveRDS(list(kappa=kappa_val, decision=KAPPA_DECISION,
               dom_xgb=dom_xgb_sub, dom_rf=dom_rf_sub),
          file.path(OUT_SHAP,"cohen_kappa_result.rds"))
  cat("  ✓ cohen_kappa_result.rds\n\n")
} else {
  cat("  RF SHAP unavailable — kappa skipped\n\n")
}

# ── 11. 250m GRID SETUP ───────────────────────────────────────

cat("--- Setting Up 250m Prediction Grid ---\n\n")

template_250m <- terra::rast(
  ext = terra::ext(template_30m),
  res = 250,
  crs = terra::crs(template_30m))

cat("  Resampling predictor stack to 250m...\n")
pred_250m <- terra::resample(pred_stack, template_250m,
                             method="bilinear")
pred_250m <- terra::mask(pred_250m, boundary_vect)

grid_df  <- terra::as.data.frame(pred_250m, xy=TRUE,
                                 cells=TRUE, na.rm=TRUE)
cell_ids <- grid_df$cell
n_cells  <- nrow(grid_df)
cat(sprintf("  250m grid: %d valid cells\n\n", n_cells))

# Prepare grid matrix
grid_mat <- as.matrix(grid_df[, final_names, drop=FALSE])
for (col in cat_in_stack)
  grid_mat[,col] <- as.numeric(as.integer(grid_mat[,col]))
storage.mode(grid_mat) <- "double"
colnames(grid_mat) <- final_names

# ── 12. XGBOOST SHAP ON 250m GRID (chunked) ──────────────────

cat("--- XGBoost SHAP on 250m Grid ---\n\n")

CHUNK_SIZE  <- 50000L
n_chunks    <- ceiling(n_cells / CHUNK_SIZE)
S_grid_list <- vector("list", n_chunks)
t0_grid     <- proc.time()

cat(sprintf("  %d cells → %d chunks of %d\n\n",
            n_cells, n_chunks, CHUNK_SIZE))

for (ch in seq_len(n_chunks)) {
  idx_s <- (ch-1L)*CHUNK_SIZE + 1L
  idx_e <- min(ch*CHUNK_SIZE, n_cells)
  
  chunk_mat  <- grid_mat[idx_s:idx_e, , drop=FALSE]
  chunk_dmat <- xgboost::xgb.DMatrix(chunk_mat)
  
  # predict(..., predcontrib=TRUE) = exact TreeSHAP
  # Returns n × (p+1); drop last col (bias)
  shap_chunk <- predict(xgb_model, chunk_dmat,
                        predcontrib=TRUE)
  S_grid_list[[ch]] <- shap_chunk[,
                                  seq_len(length(feat_names)), drop=FALSE]
  
  rm(chunk_dmat, shap_chunk); gc(full=TRUE)
  
  cat(sprintf("  Chunk %d/%d: cells %d-%d  (%.1f min)\n",
              ch, n_chunks, idx_s, idx_e,
              (proc.time()-t0_grid)[3]/60))
}

S_grid <- do.call(rbind, S_grid_list)
colnames(S_grid) <- feat_names
rm(S_grid_list); gc(full=TRUE)

cat(sprintf("\n  SHAP grid: %d x %d  total %.1f min\n\n",
            nrow(S_grid), ncol(S_grid),
            (proc.time()-t0_grid)[3]/60))

# ── 13. DOMINANT DRIVER ───────────────────────────────────────

cat("--- Dominant Driver per Cell ---\n\n")

dom_idx  <- apply(abs(S_grid), 1, which.max)
dom_name <- feat_names[dom_idx]
dom_code <- as.integer(dom_idx)

dom_freq   <- sort(table(dom_name), decreasing=TRUE)
total_cells <- length(dom_name)

cat("  Driver frequency (250m grid):\n")
for (pred in names(dom_freq)) {
  cat(sprintf("  %-22s %7d cells  %5.1f%%\n",
              pred, dom_freq[pred],
              100*dom_freq[pred]/total_cells))
}
cat("\n")

# Create dominant driver raster
dom_rast <- template_250m
terra::values(dom_rast) <- NA_integer_
dom_rast[cell_ids] <- dom_code
dom_rast <- terra::mask(dom_rast, boundary_vect)

terra::writeRaster(dom_rast,
                   file.path(OUT_SHAP,"shap_dominant_driver_map.tif"),
                   overwrite=TRUE, datatype="INT1U")

driver_lookup <- data.frame(
  code=seq_along(feat_names), predictor=feat_names,
  stringsAsFactors=FALSE)
write.csv(driver_lookup,
          file.path(OUT_SHAP,"dominant_driver_lookup.csv"),
          row.names=FALSE)

dom_freq_df <- data.frame(
  predictor=names(dom_freq),
  n_cells=as.integer(dom_freq),
  pct_cells=round(100*as.numeric(dom_freq)/total_cells,2),
  stringsAsFactors=FALSE)
write.csv(dom_freq_df,
          file.path(OUT_SHAP,"dominant_driver_frequency.csv"),
          row.names=FALSE)

cat("  ✓ shap_dominant_driver_map.tif\n")
cat("  ✓ dominant_driver_lookup.csv\n")
cat("  ✓ dominant_driver_frequency.csv\n\n")

# ── 14. LANDSCAPE ZONES ───────────────────────────────────────

cat("--- Landscape Zone Classification ---\n\n")

dem_r    <- terra::rast(file.path(OUT_PREDICTORS,
                                  "DEM_30m_utm44n.tif"))
dem_250m <- terra::resample(dem_r, template_250m,
                            method="bilinear")
dem_250m <- terra::mask(dem_250m, boundary_vect)

dem_stats <- terra::global(dem_250m,
                           c("min","max","mean","sd"),
                           na.rm=TRUE)
cat(sprintf("  DEM statistics (250m):\n"))
cat(sprintf("  min=%.0f  max=%.0f  mean=%.0f  SD=%.0f\n\n",
            dem_stats[1,1], dem_stats[1,2],
            dem_stats[1,3], dem_stats[1,4]))
cat("  Zone thresholds: <270m = Chandrapur basin,\n")
cat("  270-380m = Nagpur pediplain, >380m = Satpura foothills\n")
cat("  Verify against DEM statistics above and adjust if needed\n\n")

zone_rast <- terra::ifel(
  dem_250m > 380, 3L,
  terra::ifel(dem_250m > 270, 1L, 2L))
zone_rast <- terra::mask(zone_rast, boundary_vect)

zone_names <- c("1"="Nagpur pediplain",
                "2"="Chandrapur basin",
                "3"="Satpura foothills")

# Extract zone at grid cell centres
zone_vals <- terra::extract(
  zone_rast,
  terra::xyFromCell(template_250m, cell_ids))[,1]

zone_freq <- table(zone=zone_vals)
cat("  Zone distribution:\n")
for (z in names(zone_freq)) {
  cat(sprintf("  Zone %s %-20s %7d cells %5.1f%%\n",
              z, zone_names[z], zone_freq[z],
              100*zone_freq[z]/sum(zone_freq)))
}
cat("\n")

# ── 15. CHI-SQUARED TEST ──────────────────────────────────────

cat("--- Chi-Squared: Spatial Heterogeneity ---\n\n")

valid_chi <- !is.na(zone_vals) & !is.na(dom_name)
chi_table <- table(
  Zone   = zone_vals[valid_chi],
  Driver = dom_name[valid_chi])

# Min expected frequency check
n_chi  <- sum(chi_table)
min_ex <- min(outer(rowSums(chi_table),
                    colSums(chi_table), "*")) / n_chi
cat(sprintf("  Min expected frequency: %.2f\n", min_ex))

if (min_ex < 5) {
  chi_result <- chisq.test(chi_table,
                           simulate.p.value=TRUE, B=10000)
  cat("  Using simulation (min expected < 5)\n\n")
} else {
  chi_result <- chisq.test(chi_table)
  cat("\n")
}

cat(sprintf("  Chi-squared: %.2f\n",  chi_result$statistic))
cat(sprintf("  p-value:     %.6f\n",  chi_result$p.value))
cat(sprintf("  Conclusion:  %s\n\n",
            if(chi_result$p.value<0.05)
              "Significant spatial heterogeneity (p<0.05)"
            else "No significant difference across zones"))

# Zone-specific dominant driver
cat("  Top dominant driver per zone:\n")
zone_dr_df <- data.frame(
  zone=character(), zone_label=character(),
  top_driver=character(), pct_cells=numeric(),
  stringsAsFactors=FALSE)

for (z in rownames(chi_table)) {
  zrow    <- chi_table[z,]
  top_drv <- names(which.max(zrow))
  top_pct <- 100 * max(zrow) / sum(zrow)
  cat(sprintf("  Zone %s %-20s → %-22s %.1f%%\n",
              z, zone_names[z], top_drv, top_pct))
  zone_dr_df <- rbind(zone_dr_df, data.frame(
    zone=z, zone_label=zone_names[z],
    top_driver=top_drv,
    pct_cells=round(top_pct,1),
    stringsAsFactors=FALSE))
}
cat("\n")

write.csv(data.frame(
  chi_stat   = round(chi_result$statistic,3),
  p_value    = round(chi_result$p.value,6),
  simulated  = (min_ex<5),
  n_cells    = n_chi,
  conclusion = if(chi_result$p.value<0.05)
    "Significant" else "NS",
  stringsAsFactors=FALSE),
  file.path(OUT_SHAP,"shap_chi_squared_result.csv"),
  row.names=FALSE)
write.csv(zone_dr_df,
          file.path(OUT_SHAP,"shap_zone_dominant_drivers.csv"),
          row.names=FALSE)
cat("  ✓ shap_chi_squared_result.csv\n")
cat("  ✓ shap_zone_dominant_drivers.csv\n\n")

# ── 16. FIGURE 11 — DOMINANT DRIVER MAP ──────────────────────

cat("--- Figure 11: Dominant Driver Map ---\n\n")

tryCatch({
  
  # Extract contour BEFORE opening device
  zone_lines <- tryCatch(
    terra::as.contour(zone_rast, levels=c(1.5, 2.5)),
    error=function(e) NULL)
  
  active_drivers <- names(dom_freq)[dom_freq > 0]
  active_cols    <- PRED_COLORS[active_drivers]
  active_codes   <- match(active_drivers, feat_names)
  
  # Reclassify to sequential codes for plotting
  rcl_mat  <- cbind(active_codes, seq_along(active_codes))
  dom_plot <- terra::classify(dom_rast, rcl_mat,
                              others=NA)
  
  png(file.path(OUT_FIG_MAIN,
                "Fig11_dominant_driver_map.png"),
      width=3600, height=3600, res=300)
  
  terra::plot(dom_plot,
              main   = sprintf(
                "SHAP Dominant Predictor per Cell (250m grid)\nChi²=%.1f  p=%.4f",
                chi_result$statistic, chi_result$p.value),
              col    = active_cols,
              type   = "classes",
              axes   = FALSE,
              legend = FALSE,
              cex.main = 0.80)
  
  terra::plot(boundary_vect,
              add=TRUE, border="black", lwd=0.8)
  
  # Zone boundaries
  if (!is.null(zone_lines) && length(zone_lines) > 0) {
    terra::lines(zone_lines, col="white",
                 lwd=1.8, lty=2)
  }
  
  # Site locations
  terra::plot(terra::vect(sites_sf),
              add=TRUE, col="white", pch=16, cex=0.25)
  
  # Zone labels
  zone_label_df <- data.frame(
    x = c(745000, 815000, 720000),
    y = c(2365000, 2240000, 2195000),
    label = c("Nagpur\npediplain",
              "Chandrapur\nbasin",
              "Satpura\nfoothills"))
  for (i in seq_len(nrow(zone_label_df))) {
    text(zone_label_df$x[i], zone_label_df$y[i],
         zone_label_df$label[i],
         col="white", cex=0.65, font=2)
  }
  
  # Legend with % coverage
  leg_labels <- sprintf("%s\n(%.1f%%)",
                        active_drivers,
                        100*dom_freq[active_drivers]/total_cells)
  
  legend("bottomright",
         legend = leg_labels,
         fill   = active_cols,
         title  = "Dominant Predictor",
         bty    = "n",
         cex    = 0.50,
         ncol   = if(length(active_drivers)>7) 2L else 1L)
  
  dev.off()
  cat("  ✓ Fig11_dominant_driver_map.png\n\n")
}, error=function(e) {
  tryCatch(dev.off(), error=function(x) NULL)
  cat("  ✗ Fig11 error:", e$message, "\n\n")
})

# ── 17. GLOBAL IMPORTANCE CSV ─────────────────────────────────

imp_csv <- data.frame(
  predictor         = names(imp_xgb_ranked),
  xgb_mean_abs_shap = round(imp_xgb_ranked, 4),
  xgb_rank          = seq_along(imp_xgb_ranked),
  rf_importance     = round(imp_rf[names(imp_xgb_ranked)], 4),
  rf_type           = if(RF_SHAP_SUCCESS && !is.null(S_rf))
    "SHAP" else "MDA",
  stringsAsFactors  = FALSE
)
write.csv(imp_csv,
          file.path(OUT_SHAP,"shap_global_importance.csv"),
          row.names=FALSE)
cat("  ✓ shap_global_importance.csv\n\n")

# ── 18. SUMMARY ───────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 19 COMPLETE — SHAP Analysis\n")
cat("========================================\n")
cat(sprintf("XGBoost SHAP: predcontrib=TRUE (v3.2.x fix) ✓\n"))
cat(sprintf("RF SHAP:      %s\n",
            if(RF_SHAP_SUCCESS && !is.null(S_rf))
              "treeshap ✓" else "MDA fallback"))
cat(sprintf("Grid cells:   %d (250m resolution)\n", n_cells))
cat(sprintf("Active drivers: %d/%d predictors\n",
            length(active_drivers), length(feat_names)))
cat(sprintf("Cohen kappa:  %.4f (%s)\n",
            kappa_val, KAPPA_DECISION))
cat(sprintf("Chi-squared:  %.2f  p=%.6f\n",
            chi_result$statistic, chi_result$p.value))
cat("\nTop driver per zone:\n")
for (i in seq_len(nrow(zone_dr_df))) {
  cat(sprintf("  Zone %s %-20s → %s (%.1f%%)\n",
              zone_dr_df$zone[i], zone_dr_df$zone_label[i],
              zone_dr_df$top_driver[i],
              zone_dr_df$pct_cells[i]))
}
cat("\nXGBoost top 5 predictors (mean |SHAP|):\n")
for (i in 1:5) {
  cat(sprintf("  %d. %-22s %.4f\n", i,
              names(imp_xgb_ranked)[i],
              imp_xgb_ranked[i]))
}
cat("\nFiles saved:\n")
cat("  outputs/tables/Table4_SHAP_importance.csv\n")
cat("  outputs/shap/shap_global_importance.csv\n")
cat("  outputs/shap/shap_dominant_driver_map.tif\n")
cat("  outputs/shap/dominant_driver_lookup.csv\n")
cat("  outputs/shap/dominant_driver_frequency.csv\n")
cat("  outputs/shap/cohen_kappa_result.rds\n")
cat("  outputs/shap/shap_chi_squared_result.csv\n")
cat("  outputs/shap/shap_zone_dominant_drivers.csv\n")
cat("  outputs/figures/main/Fig09_SHAP_beeswarm.png\n")
cat("  outputs/figures/main/Fig10_SHAP_dependence.png\n")
cat("  outputs/figures/main/Fig11_dominant_driver_map.png\n")
cat("\nNOTE: Verify zone thresholds (270m, 380m) against\n")
cat("DEM statistics printed above. Adjust before manuscript.\n")
cat("\nNext: Script 20 — Cultural Period Sub-Models\n")
cat("========================================\n")