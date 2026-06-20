# ============================================================
# SCRIPT 22: EVALUATION METRICS — TABLES 3 AND 5
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 22 of 25
# ============================================================
# WHAT THIS SCRIPT DOES:
#   Compiles all model evaluation results into publication-
#   ready Tables 3 and 5 (Research Design Sections 10.2).
#
#   TABLE 3 — Model performance:
#     Rows: 6 algorithms + Ensemble (primary) +
#           LP/MP/UP ensemble sub-models + Transfer
#     Cols: Spatial CV AUC (mean±SD), Full AUC,
#           Boyce Index, TSS (max-TSS), Kvamme's Gain
#
#   TABLE 5 — DeLong pairwise AUC comparison:
#     All algorithm pairs + ensemble vs each individual
#     + AUC-weighted vs equal-weight
#     Cols: Algorithm A, Algorithm B, AUC_A, AUC_B,
#           Delta AUC, DeLong p-value, Significant
#
#   Also produces:
#     Figure 7: CV AUC bar chart with 95% CI across
#               all algorithms + ensemble (±1.96×SD)
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

cat("\n========================================\n")
cat("SCRIPT 22: Evaluation Metrics Tables\n")
cat("========================================\n\n")

# ─────────────────────────────────────────────────────────────
# 1. LOAD ALL EVALUATION FILES
# ─────────────────────────────────────────────────────────────

cat("--- Loading Evaluation Files ---\n\n")

safe_read <- function(path, label) {
  if (file.exists(path)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    cat(sprintf("  ✓ %-45s %d rows\n",
                basename(path), nrow(df)))
    return(df)
  }
  cat(sprintf("  ✗ MISSING: %s\n", basename(path)))
  return(NULL)
}

eval_mx  <- safe_read(file.path(OUT_EVAL,
                                "maxent_evaluation.csv"),
                      "MaxEnt")
eval_rf  <- safe_read(file.path(OUT_EVAL,
                                "rf_evaluation.csv"),
                      "RF")
eval_xgb <- safe_read(file.path(OUT_EVAL,
                                "xgboost_evaluation.csv"),
                      "XGBoost")
eval_brt <- safe_read(file.path(OUT_EVAL,
                                "brt_evaluation.csv"),
                      "BRT")
eval_gam <- safe_read(file.path(OUT_EVAL,
                                "gam_evaluation.csv"),
                      "GAM")
eval_svm <- safe_read(file.path(OUT_EVAL,
                                "svm_evaluation.csv"),
                      "SVM")
eval_ens <- safe_read(file.path(OUT_EVAL,
                                "ensemble_evaluation.csv"),
                      "Ensemble")
eval_sub <- safe_read(file.path(OUT_EVAL,
                                "submodel_evaluation.csv"),
                      "Sub-models")
eval_trf <- safe_read(file.path(OUT_EVAL,
                                "transfer_evaluation.csv"),
                      "Transfer")
delong   <- safe_read(file.path(OUT_EVAL,
                                "delong_pairwise_results.csv"),
                      "DeLong")
cat("\n")

# ─────────────────────────────────────────────────────────────
# 2. HELPER: safe extract from eval row
# ─────────────────────────────────────────────────────────────

get_val <- function(df, col, default = NA_real_,
                    row = 1L) {
  if (is.null(df) || !col %in% names(df)) return(default)
  v <- suppressWarnings(as.numeric(df[[col]][row]))
  if (!length(v) || is.na(v)) return(default)
  v
}

fmt <- function(x, d = 4) {
  if (is.na(x)) return(NA_character_)
  formatC(round(x, d), format = "f", digits = d)
}

fmt_auc_sd <- function(m, s) {
  if (is.na(m)) return(NA_character_)
  if (is.na(s)) return(fmt(m, 4))
  sprintf("%.4f \u00b1 %.4f", m, s)
}

# ─────────────────────────────────────────────────────────────
# 3. TABLE 3 — MAIN MODEL PERFORMANCE
# ─────────────────────────────────────────────────────────────

cat("--- Building Table 3: Model Performance ---\n\n")

# Helper to build one row from an eval data frame
make_row <- function(label, df, row = 1L) {
  data.frame(
    Model          = label,
    CV_AUC_mean    = get_val(df, "cv_auc_mean", row = row),
    CV_AUC_SD      = get_val(df, "cv_auc_sd",   row = row),
    CV_AUC_display = fmt_auc_sd(
      get_val(df, "cv_auc_mean", row = row),
      get_val(df, "cv_auc_sd",   row = row)),
    Full_AUC       = get_val(df, "full_auc",     row = row),
    Boyce_Index    = get_val(df, "boyce_index",  row = row),
    TSS_max        = get_val(df, "tss_max",      row = row),
    Kvamme_Gain    = get_val(df, "kvamme_gain",  row = row),
    Fold1_AUC      = get_val(df, "fold1_auc",    row = row),
    Fold2_AUC      = get_val(df, "fold2_auc",    row = row),
    Fold3_AUC      = get_val(df, "fold3_auc",    row = row),
    Fold4_AUC      = get_val(df, "fold4_auc",    row = row),
    Fold5_AUC      = get_val(df, "fold5_auc",    row = row),
    stringsAsFactors = FALSE)
}

# Sub-model ensemble rows (period = LP/MP/UP, alg = Ensemble)
sub_row <- function(period_label, sub_df) {
  if (is.null(sub_df)) return(NULL)
  r <- sub_df[sub_df$algorithm == "Ensemble_AUC_wtd" &
                sub_df$period  == period_label, ]
  if (nrow(r) == 0) return(NULL)
  data.frame(
    Model          = paste0(period_label, " Ensemble"),
    CV_AUC_mean    = get_val(r, "cv_auc_mean"),
    CV_AUC_SD      = NA_real_,
    CV_AUC_display = fmt(get_val(r, "cv_auc_mean"), 4),
    Full_AUC       = NA_real_,
    Boyce_Index    = get_val(r, "boyce_index"),
    TSS_max        = get_val(r, "tss_max"),
    Kvamme_Gain    = get_val(r, "kvamme_gain"),
    Fold1_AUC = NA_real_, Fold2_AUC = NA_real_,
    Fold3_AUC = NA_real_, Fold4_AUC = NA_real_,
    Fold5_AUC = NA_real_,
    stringsAsFactors = FALSE)
}

# Transfer row
trf_val <- function(metric) {
  if (is.null(eval_trf)) return(NA_real_)
  r <- eval_trf[eval_trf$metric == metric, ]
  if (nrow(r) == 0) return(NA_real_)
  suppressWarnings(as.numeric(r$value[1]))
}

trf_row <- data.frame(
  Model          = "Transfer (RF+MaxEnt, southern sector)",
  CV_AUC_mean    = trf_val("ensemble_transfer_auc"),
  CV_AUC_SD      = NA_real_,
  CV_AUC_display = fmt(trf_val("ensemble_transfer_auc"), 4),
  Full_AUC       = NA_real_,
  Boyce_Index    = trf_val("ensemble_transfer_boyce"),
  TSS_max        = trf_val("ensemble_transfer_tss"),
  Kvamme_Gain    = NA_real_,
  Fold1_AUC = NA_real_, Fold2_AUC = NA_real_,
  Fold3_AUC = NA_real_, Fold4_AUC = NA_real_,
  Fold5_AUC = NA_real_,
  stringsAsFactors = FALSE)

table3 <- do.call(rbind, Filter(Negate(is.null), list(
  make_row("MaxEnt",   eval_mx),
  make_row("RF",       eval_rf),
  make_row("XGBoost",  eval_xgb),
  make_row("BRT",      eval_brt),
  make_row("GAM",      eval_gam),
  make_row("SVM",      eval_svm),
  make_row("Ensemble (AUC-weighted, primary)", eval_ens),
  sub_row("LP", eval_sub),
  sub_row("MP", eval_sub),
  sub_row("UP", eval_sub),
  trf_row
)))

rownames(table3) <- NULL

# Print Table 3
cat(sprintf("  %-40s  %18s  %7s  %7s  %7s  %7s\n",
            "Model", "CV AUC (mean\u00b1SD)",
            "Full", "Boyce", "TSS", "KG"))
cat("  ", paste(rep("-", 100), collapse = ""), "\n")
for (i in seq_len(nrow(table3))) {
  r <- table3[i, ]
  cat(sprintf("  %-40s  %18s  %7s  %7s  %7s  %7s\n",
              r$Model,
              ifelse(is.na(r$CV_AUC_display), "—",
                     r$CV_AUC_display),
              ifelse(is.na(r$Full_AUC),    "—",
                     fmt(r$Full_AUC, 4)),
              ifelse(is.na(r$Boyce_Index), "—",
                     fmt(r$Boyce_Index, 3)),
              ifelse(is.na(r$TSS_max),     "—",
                     fmt(r$TSS_max, 3)),
              ifelse(is.na(r$Kvamme_Gain), "—",
                     fmt(r$Kvamme_Gain, 3))))
}

write.csv(table3,
          file.path(OUT_TABLES, "Table3_model_performance.csv"),
          row.names = FALSE)
cat("\n  ✓ Table3_model_performance.csv\n\n")

# ─────────────────────────────────────────────────────────────
# 4. TABLE 5 — DELONG PAIRWISE AUC COMPARISONS
# ─────────────────────────────────────────────────────────────

cat("--- Building Table 5: DeLong Pairwise AUC ---\n\n")

if (!is.null(delong) && nrow(delong) > 0) {
  # Clean and standardise DeLong table
  table5 <- delong
  
  # Ensure required columns exist
  req_cols <- c("algorithm_A", "algorithm_B",
                "auc_A", "auc_B", "delta_auc",
                "delong_p", "significant_p05")
  missing_cols <- req_cols[!req_cols %in% names(table5)]
  
  if (length(missing_cols) > 0) {
    cat(sprintf("  ⚠ Missing columns: %s\n",
                paste(missing_cols, collapse = ", ")))
    cat("  Attempting column name inference...\n")
    # Try common alternatives
    alt_map <- list(
      algorithm_A    = c("alg_A","algo_A","model_A","pair_A"),
      algorithm_B    = c("alg_B","algo_B","model_B","pair_B"),
      auc_A          = c("AUC_A","auc1","AUC1"),
      auc_B          = c("AUC_B","auc2","AUC2"),
      delta_auc      = c("delta","diff","AUC_diff"),
      delong_p       = c("p_value","pvalue","p.value","p"),
      significant_p05 = c("significant","sig","p_sig"))
    for (std in missing_cols) {
      for (alt in alt_map[[std]]) {
        if (alt %in% names(table5)) {
          names(table5)[names(table5) == alt] <- std
          break
        }
      }
    }
  }
  
  # Add significance stars
  if ("delong_p" %in% names(table5)) {
    table5$sig_stars <- ifelse(
      is.na(table5$delong_p), "",
      ifelse(table5$delong_p < 0.001, "***",
             ifelse(table5$delong_p < 0.01, "**",
                    ifelse(table5$delong_p < 0.05,
                           "*", "ns"))))
  }
  
  # Print summary
  cat(sprintf("  %-22s  %-22s  %6s  %6s  %8s  %4s\n",
              "Algorithm A", "Algorithm B",
              "AUC_A", "AUC_B", "p-value", "Sig"))
  cat("  ", paste(rep("-", 80), collapse = ""), "\n")
  for (i in seq_len(nrow(table5))) {
    r <- table5[i, ]
    cat(sprintf(
      "  %-22s  %-22s  %6.4f  %6.4f  %8.4f  %s\n",
      ifelse("algorithm_A" %in% names(r),
             as.character(r$algorithm_A), "?"),
      ifelse("algorithm_B" %in% names(r),
             as.character(r$algorithm_B), "?"),
      ifelse("auc_A" %in% names(r) &&
               is.finite(as.numeric(r$auc_A)),
             as.numeric(r$auc_A), NA),
      ifelse("auc_B" %in% names(r) &&
               is.finite(as.numeric(r$auc_B)),
             as.numeric(r$auc_B), NA),
      ifelse("delong_p" %in% names(r) &&
               is.finite(as.numeric(r$delong_p)),
             as.numeric(r$delong_p), NA),
      ifelse("sig_stars" %in% names(r),
             as.character(r$sig_stars), "")))
  }
  
  write.csv(table5,
            file.path(OUT_TABLES,
                      "Table5_DeLong_pairwise.csv"),
            row.names = FALSE)
  cat("\n  ✓ Table5_DeLong_pairwise.csv\n\n")
  
  # Count significant comparisons
  if ("sig_stars" %in% names(table5)) {
    n_sig <- sum(table5$sig_stars %in% c("*","**","***"),
                 na.rm = TRUE)
    n_tot <- nrow(table5)
    cat(sprintf("  Significant (p<0.05): %d / %d pairs\n\n",
                n_sig, n_tot))
  }
  
} else {
  cat("  DeLong results not found — recomputing...\n\n")
  
  # Recompute from CV predictions
  alg_names <- c("MaxEnt","RF","XGBoost","BRT","GAM","SVM")
  alg_files <- c("maxent","rf","xgboost","brt","gam","svm")
  
  cv_preds <- list()
  for (i in seq_along(alg_names)) {
    fp <- file.path(OUT_MOD_IND,
                    paste0(alg_files[i],
                           "_cv_predictions.rds"))
    if (file.exists(fp)) {
      p <- readRDS(fp)
      cv_preds[[alg_names[i]]] <- p
    }
  }
  
  rows <- list()
  algs_loaded <- names(cv_preds)
  for (i in seq_along(algs_loaded)) {
    for (j in seq_along(algs_loaded)) {
      if (j <= i) next
      a <- algs_loaded[i]; b <- algs_loaded[j]
      pa <- cv_preds[[a]]; pb_a <- cv_preds[[a]]
      pa_b <- cv_preds[[b]]
      
      sp_a <- pa$site_preds; bg_a <- pa$bg_preds
      sp_b <- pa_b$site_preds; bg_b <- pa_b$bg_preds
      n_s  <- min(length(sp_a), length(sp_b))
      n_b  <- min(length(bg_a), length(bg_b))
      
      roc_a <- tryCatch(pROC::roc(
        c(rep(1,n_s),rep(0,n_b)),
        c(sp_a[1:n_s], bg_a[1:n_b]), quiet=TRUE),
        error=function(e) NULL)
      roc_b <- tryCatch(pROC::roc(
        c(rep(1,n_s),rep(0,n_b)),
        c(sp_b[1:n_s], bg_b[1:n_b]), quiet=TRUE),
        error=function(e) NULL)
      
      if (!is.null(roc_a) && !is.null(roc_b)) {
        tst <- tryCatch(
          pROC::roc.test(roc_a, roc_b,
                         method = "delong"),
          error = function(e) NULL)
        p_val <- if (!is.null(tst)) tst$p.value else NA_real_
        rows[[length(rows)+1]] <- data.frame(
          algorithm_A     = a,
          algorithm_B     = b,
          auc_A           = round(as.numeric(
            pROC::auc(roc_a)), 4),
          auc_B           = round(as.numeric(
            pROC::auc(roc_b)), 4),
          delta_auc       = round(abs(
            as.numeric(pROC::auc(roc_a)) -
              as.numeric(pROC::auc(roc_b))),4),
          delong_p        = round(p_val, 4),
          significant_p05 = ifelse(!is.na(p_val) &
                                     p_val < 0.05, 1L, 0L),
          sig_stars       = ifelse(is.na(p_val), "",
                                   ifelse(p_val<0.001,"***",
                                          ifelse(p_val<0.01, "**",
                                                 ifelse(p_val<0.05, "*","ns")))),
          stringsAsFactors = FALSE)
      }
    }
  }
  
  if (length(rows) > 0) {
    table5 <- do.call(rbind, rows)
    write.csv(table5,
              file.path(OUT_TABLES,
                        "Table5_DeLong_pairwise.csv"),
              row.names = FALSE)
    cat(sprintf("  ✓ Recomputed %d pairs\n\n",
                nrow(table5)))
  } else {
    cat("  ✗ Could not compute DeLong table\n\n")
    table5 <- NULL
  }
}

# ─────────────────────────────────────────────────────────────
# 5. FIGURE 7 — CV AUC BAR CHART WITH 95% CI
# ─────────────────────────────────────────────────────────────

cat("--- Figure 7: CV AUC Comparison Bar Chart ---\n\n")

tryCatch({
  # Collect AUC values for all 7 models
  alg_order <- c("MaxEnt","RF","XGBoost","BRT",
                 "GAM","SVM","Ensemble")
  eval_list_all <- list(
    MaxEnt   = eval_mx,
    RF       = eval_rf,
    XGBoost  = eval_xgb,
    BRT      = eval_brt,
    GAM      = eval_gam,
    SVM      = eval_svm,
    Ensemble = eval_ens)
  
  auc_m <- sapply(alg_order, function(a)
    get_val(eval_list_all[[a]], "cv_auc_mean"))
  auc_s <- sapply(alg_order, function(a)
    get_val(eval_list_all[[a]], "cv_auc_sd"))
  
  # 95% CI = mean ± 1.96 × SD
  ci_hi <- auc_m + 1.96 * auc_s
  ci_lo <- auc_m - 1.96 * auc_s
  ci_hi[ci_hi > 1] <- 1
  ci_lo[ci_lo < 0] <- 0
  
  bar_cols <- c(
    MaxEnt   = "#2A9D8F", RF      = "#E63946",
    XGBoost  = "#F4A261", BRT     = "#264653",
    GAM      = "#E9C46A", SVM     = "#8338EC",
    Ensemble = "#06D6A0")
  
  valid_idx <- which(is.finite(auc_m))
  
  png(file.path(OUT_FIG_MAIN,
                "Fig07_AUC_comparison.png"),
      width = 3600, height = 2800, res = 300)
  
  par(mar = c(5, 5, 4, 2))
  bp <- barplot(
    auc_m[valid_idx],
    names.arg = alg_order[valid_idx],
    col       = bar_cols[valid_idx],
    ylim      = c(0.50, 1.00),
    xpd       = FALSE,
    main      = "Spatial Block CV AUC — All Models\n(±1.96 SD = 95% CI)",
    ylab      = "Spatial Block CV AUC (5-fold mean)",
    xlab      = "",
    border    = NA,
    cex.main  = 1.1,
    cex.lab   = 0.95,
    cex.names = 0.85,
    las       = 1)
  
  # Performance threshold lines
  abline(h = 0.75, col = "grey60",
         lty = 2, lwd = 1.2)
  abline(h = 0.85, col = "grey40",
         lty = 3, lwd = 1.2)
  
  # Error bars (95% CI)
  for (i in seq_along(valid_idx)) {
    idx <- valid_idx[i]
    if (!is.na(ci_hi[idx]) && !is.na(ci_lo[idx]))
      arrows(bp[i], ci_lo[idx], bp[i], ci_hi[idx],
             code = 3, angle = 90, length = 0.05,
             col = "black", lwd = 1.5)
  }
  
  # Value labels above bars
  text(bp, auc_m[valid_idx] + 0.008,
       labels = sprintf("%.4f", auc_m[valid_idx]),
       cex = 0.75, col = "black")
  
  # Threshold legend
  legend("topright",
         legend = c("AUC=0.75 (adequate)",
                    "AUC=0.85 (strong)"),
         col    = c("grey60", "grey40"),
         lty    = c(2, 3), lwd = c(1.2, 1.2),
         bty    = "n", cex = 0.80)
  
  # Ensemble star
  ens_idx <- which(valid_idx == which(
    alg_order == "Ensemble"))
  if (length(ens_idx) > 0 &&
      !is.na(auc_m[which(alg_order == "Ensemble")])) {
    text(bp[ens_idx],
         auc_m[which(alg_order == "Ensemble")] + 0.022,
         labels = "\u2605", col = "#06D6A0", cex = 1.2)
  }
  
  dev.off()
  cat("  ✓ Fig07_AUC_comparison.png\n\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  cat(sprintf("  ✗ Fig07 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 6. SUMMARY STATISTICS FOR MANUSCRIPT
# ─────────────────────────────────────────────────────────────

cat("--- Summary Statistics for Manuscript ---\n\n")

# Best individual algorithm
best_idx <- which.max(table3$CV_AUC_mean[1:6])
best_alg <- table3$Model[best_idx]
best_auc <- table3$CV_AUC_mean[best_idx]

# Ensemble vs best individual
ens_row_idx <- which(grepl("Ensemble.*primary",
                           table3$Model,
                           ignore.case = TRUE))
ens_auc <- if (length(ens_row_idx) > 0)
  table3$CV_AUC_mean[ens_row_idx[1]] else NA_real_

# DeLong significant pairs
if (!is.null(table5) && "sig_stars" %in% names(table5)) {
  n_sig_ens <- if ("algorithm_A" %in% names(table5))
    sum(grepl("Ensemble", table5$algorithm_A) &
          table5$sig_stars %in% c("*","**","***"),
        na.rm = TRUE) +
    sum(grepl("Ensemble", table5$algorithm_B) &
          table5$sig_stars %in% c("*","**","***"),
        na.rm = TRUE)
  else 0L
} else {
  n_sig_ens <- NA_integer_
}

cat(sprintf("  Best individual algorithm: %s  AUC=%.4f\n",
            best_alg, best_auc))
cat(sprintf("  Ensemble CV AUC:           %.4f\n",
            safe_scalar(ens_auc)))
cat(sprintf("  Ensemble Boyce Index:      %.4f\n",
            get_val(eval_ens, "boyce_index")))
cat(sprintf("  Transfer AUC:              %.4f\n",
            trf_val("ensemble_transfer_auc")))
cat(sprintf("  Transfer Boyce:            %.4f\n",
            trf_val("ensemble_transfer_boyce")))
cat(sprintf("  Transfer Delta AUC:        %.4f\n",
            trf_val("delta_auc_ensemble")))

if (!is.na(n_sig_ens))
  cat(sprintf("  DeLong sig vs ensemble:    %d pairs\n",
              n_sig_ens))

# AUC range across 6 algorithms
auc_6 <- table3$CV_AUC_mean[1:6]
auc_6 <- auc_6[is.finite(auc_6)]
if (length(auc_6) > 1)
  cat(sprintf("  AUC range (6 alg):  %.4f – %.4f\n",
              min(auc_6), max(auc_6)))

cat("\n")

# ─────────────────────────────────────────────────────────────
# 7. EXPORT FOLD-LEVEL AUC TABLE (Supplementary)
# ─────────────────────────────────────────────────────────────

cat("--- Fold-Level AUC Table ---\n\n")

fold_cols <- paste0("fold", 1:5, "_auc")
fold_rows <- list()

for (i in seq_len(min(7, nrow(table3)))) {
  r <- table3[i, ]
  fold_vals <- sapply(fold_cols, function(c)
    ifelse(c %in% names(r) && !is.na(r[[c]]),
           r[[c]], NA_real_))
  fold_rows[[i]] <- data.frame(
    Model  = r$Model,
    Fold1  = fold_vals[1], Fold2 = fold_vals[2],
    Fold3  = fold_vals[3], Fold4 = fold_vals[4],
    Fold5  = fold_vals[5],
    Mean   = r$CV_AUC_mean, SD = r$CV_AUC_SD,
    stringsAsFactors = FALSE)
}

fold_table <- do.call(rbind, fold_rows)

# Print fold table
cat(sprintf("  %-30s  %6s  %6s  %6s  %6s  %6s  %6s  %6s\n",
            "Model", "F1","F2","F3","F4","F5",
            "Mean","SD"))
cat("  ", paste(rep("-", 80), collapse=""), "\n")
for (i in seq_len(nrow(fold_table))) {
  r <- fold_table[i, ]
  cat(sprintf("  %-30s  %6s  %6s  %6s  %6s  %6s  %6.4f  %6.4f\n",
              r$Model,
              ifelse(is.na(r$Fold1),"—",sprintf("%.4f",r$Fold1)),
              ifelse(is.na(r$Fold2),"—",sprintf("%.4f",r$Fold2)),
              ifelse(is.na(r$Fold3),"—",sprintf("%.4f",r$Fold3)),
              ifelse(is.na(r$Fold4),"—",sprintf("%.4f",r$Fold4)),
              ifelse(is.na(r$Fold5),"—",sprintf("%.4f",r$Fold5)),
              ifelse(is.na(r$Mean), 0, r$Mean),
              ifelse(is.na(r$SD),   0, r$SD)))
}

write.csv(fold_table,
          file.path(OUT_TABLES,
                    "TableS_fold_level_AUC.csv"),
          row.names = FALSE)
cat("\n  ✓ TableS_fold_level_AUC.csv\n\n")

# ─────────────────────────────────────────────────────────────
# 8. SUMMARY
# ─────────────────────────────────────────────────────────────

cat("========================================\n")
cat("SCRIPT 22 COMPLETE — Evaluation Tables\n")
cat("========================================\n\n")

cat(sprintf("Table 3 rows:     %d\n", nrow(table3)))
if (!is.null(table5))
  cat(sprintf("Table 5 pairs:    %d\n", nrow(table5)))

cat("\nFiles saved:\n")
cat("  outputs/tables/Table3_model_performance.csv\n")
cat("  outputs/tables/Table5_DeLong_pairwise.csv\n")
cat("  outputs/tables/TableS_fold_level_AUC.csv\n")
cat("  outputs/figures/main/Fig07_AUC_comparison.png\n")
cat("\nNext: Script 23 — Main Figures\n")
cat("========================================\n")