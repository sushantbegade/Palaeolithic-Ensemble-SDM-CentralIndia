# Palaeolithic Ensemble SDM — Central India

**Paper:** A Multi-Model Ensemble Framework with Spatial Explainability for Predicting Open-Air Palaeolithic Site Distribution in Central India: Bias Correction, Algorithmic Diversity, and Spatially Heterogeneous Predictor Importance

**Authors:** Sushant Begade*¹ and Prabash Sahu¹

**Affiliation:** ¹ Department of Ancient Indian History, Culture and Archaeology, Rashtrasant Tukadoji Maharaj Nagpur University, Nagpur – 440033, India

**ORCIDs:** Sushant Begade: [0009-0003-0804-1763](https://orcid.org/0009-0003-0804-1763) | Prabash Sahu: [0000-0003-0691-0403](https://orcid.org/0000-0003-0691-0403)

**Target Journal:** Journal of Archaeological Science 

**Status:** Under Review

---

## Study Overview

This repository contains all R analysis scripts for a six-algorithm ensemble species distribution model (SDM) applied to open-air Palaeolithic site prediction across the Wainganga-Wardha basin, Nagpur and Chandrapur Districts, central India (~21,300 km²).

**Key methodological contributions:**
- Six-algorithm AUC-weighted ensemble (MaxEnt, RF, XGBoost, BRT, GAM, SVM) with all outputs unified on logistic probability scale
- Survey-effort bias correction via kernel density-weighted background sampling (Phillips et al. 2009)
- Spatial block cross-validation (blockCV 3.1; 50 km blocks; variogram range = 50 km)
- SHAP (SHapley Additive exPlanations) spatial decomposition — first application to Indian Palaeolithic predictive modelling, to our knowledge
- Novel palaeochannel proximity predictor (Sentinel-2A MNDWI ∩ DEM valley morphology)
- Geographic transfer validation on held-out southern Chandrapur sector
- Cultural period sub-models (LP, MP, UP) with attribution uncertainty sensitivity analysis

---

## Key Results

| Algorithm | CV AUC ± SD | Full AUC | Boyce Index | TSS | Kvamme's Gain |
|---|---|---|---|---|---|
| MaxEnt | 0.7254 ± 0.038 | 0.794 | 0.867 | 0.336 | 0.716 |
| RF | 0.7274 ± 0.028 | 0.998 | 0.876 | 0.329 | 0.842 |
| XGBoost | 0.6867 ± 0.019 | 0.999 | 0.836 | 0.279 | 0.960 |
| BRT | 0.6921 ± 0.009 | 1.000 | 0.908 | 0.284 | 0.998 |
| GAM | 0.6805 ± 0.036 | 0.794 | 0.735 | 0.280 | 0.593 |
| SVM | 0.6265 ± 0.019 | 0.860 | 0.852 | 0.138 | 0.970 |
| **Ensemble (AUC-wtd)** | **0.7239** | — | **0.909** | **0.324** | **0.606** |

**Ensemble outperforms:** XGBoost (p=0.003), BRT (p=0.010), GAM (p<0.001), SVM (p<0.001) — DeLong's test

**SHAP dominant drivers (XGBoost, 250m grid, N=339,975 cells):**
1. TRI — 34.4% of cells
2. Elevation — 17.7%
3. Geomorphology — 15.8%
4. NDVI — 13.8%

Chi-squared spatial heterogeneity: χ²=54,352, p<0.001 — predictor importance is significantly non-stationary across the landscape.

**Transfer validation (southern Chandrapur, N=7 independent sites):** AUC=0.823, Boyce=0.929, delta AUC=0.099 < 0.10 threshold → transferable within basin.

**Background sensitivity (N=1,000–20,000):** AUC range=0.0099 < 0.02 → N=10,000 confirmed adequate (Warton & Shepherd 2010).

---

## Repository Structure

```
Palaeolithic-Ensemble-SDM-CentralIndia/
├── scripts/
│   ├── 01_setup_packages.R
│   ├── 02_load_and_reproject.R
│   ├── 03_terrain_derivatives.R
│   ├── 04_sentinel_processing.R
│   ├── 05_palaeochannel.R
│   ├── 06_climate_downscaling.R
│   ├── 07_predictor_stack_vif.R
│   ├── 08_sites_thinning.R
│   ├── 09_bias_correction_background.R
│   ├── 10_spatial_cv_design.R
│   ├── 11_model_maxent.R
│   ├── 12_model_rf.R
│   ├── 13_model_xgboost.R
│   ├── 14_model_brt.R
│   ├── 15_model_gam.R
│   ├── 16_model_svm.R
│   ├── 16b_pre_ensemble_issue_resolution.R
│   ├── 16b2_pre_ensemble_issue_resolution.R
│   ├── 17_ensemble.R
│   ├── 18_background_sensitivity.R
│   ├── 19_shap_analysis.R
│   ├── 20_submodels_LP_MP_UP.R
│   ├── 21_transfer_validation.R
│   ├── 22_evaluation_metrics.R
│   ├── 23_figures_main.R
│   ├── 24_figures_supplementary.R
│   └── 25_tables_export.R
├── renv.lock
├── Nagpur_Chandrapur_Ensemble.Rproj (Create this file using RStudio)
└── README.md
```

---

## Script Descriptions

| Script | Purpose |
|---|---|
| `01_setup_packages.R` | Install and load all required R packages; set global options |
| `02_load_and_reproject.R` | Load all raw data, reproject to UTM Zone 44N (EPSG:32644) |
| `03_terrain_derivatives.R` | DEM-derived predictors: Aspect, TRI, TPI, Plan Curvature, HAND, Flow Accumulation (WhiteboxTools v2.3) |
| `04_sentinel_processing.R` | NDVI (dry season) and MNDWI from Sentinel-2A Level-2A, May 2025 composite |
| `05_palaeochannel.R` | Palaeochannel layer: MNDWI anomaly (>−0.10) ∩ DEM valley morphology; distance raster |
| `06_climate_downscaling.R` | CHELSA-TraCE21k delta downscaling to 30m UTM 44N; Kelvin→Celsius correction for BIO1 |
| `07_predictor_stack_vif.R` | Assemble 14-variable predictor stack; VIF screening (usdm::vifcor, threshold < 5); Slope excluded (VIF=74.82); 13 variables retained |
| `08_sites_thinning.R` | Spatial thinning at 1 km minimum distance (spThin); cultural period flags (LP/MP/UP) |
| `09_bias_correction_background.R` | KDE survey-effort bias surface (ks package, Hscv bandwidth); 10,000 bias-weighted background points within district boundary |
| `10_spatial_cv_design.R` | Spatial block CV design (blockCV 3.1); variogram-based block size = 50 km; 5-fold assignment for sites and background |
| `11_model_maxent.R` | MaxEnt via maxnet 0.1.4; ENMeval 2.0 tuning (35 FC×RM combinations, AICc); **always type="logistic"**; best params: FC=LQH, RM=0.5 |
| `12_model_rf.R` | Random Forest; balanced bootstrap (NOT classwt); ntree=1000, mtry=3 |
| `13_model_xgboost.R` | XGBoost; objective="binary:logistic"; nrounds=100 floor; eta=0.05; scale_pos_weight=52.9 |
| `14_model_brt.R` | Boosted Regression Trees (gbm); interaction.depth=5, shrinkage=0.01, bag.fraction=0.75; n.trees=6058 |
| `15_model_gam.R` | GAM (mgcv); **select=TRUE + method="REML"** (never dredge()); k=5 all continuous; all 13 terms retained (none shrunk to zero) |
| `16_model_svm.R` | SVM (kernlab); RBF kernel; C=1.0, sigma=0.010; Platt scaling for probabilities; z-standardised continuous predictors |
| `16b_pre_ensemble_issue_resolution.R` | Issue resolution: MaxEnt Boyce NA fix (manual Spearman); XGBoost nrounds floor 50→100 |
| `16b2_pre_ensemble_issue_resolution.R` | Issue resolution: Dist_RawMat predictor fix — selective knappable lithologies only; mean distance 36m→6742m |
| `17_ensemble.R` | AUC-weighted ensemble (primary; DeLong p=0.005 vs equal-weight); equal-weight ensemble; all outputs on logistic [0,1] scale |
| `18_background_sensitivity.R` | Ensemble SD uncertainty surface; confidence zones (2×2 suitability×uncertainty); background sensitivity (N=1k/5k/10k/20k) |
| `19_shap_analysis.R` | TreeSHAP via shapviz; XGBoost (full data) + RF (500-row subsample); 250m grid (N=339,975 cells); dominant driver map; Cohen's kappa |
| `20_submodels_LP_MP_UP.R` | Cultural period sub-models (LP N=70, MP N=110, UP N=74); 5-algorithm ensemble (GAM excluded — REML timeout); attribution sensitivity analysis |
| `21_transfer_validation.R` | Geographic hold-out: southern Chandrapur sector (Gondpipri–Korpana–Rajura; N=7 sites); transfer AUC=0.823, delta=0.099 |
| `22_evaluation_metrics.R` | AUC, Boyce Index, TSS (max-TSS threshold), Kvamme's Gain, DeLong's pairwise tests; all algorithms + ensemble + sub-models |
| `23_figures_main.R` | Main manuscript Figures 1–12 |
| `24_figures_supplementary.R` | Supplementary Figures S1–S6 |
| `25_tables_export.R` | All tables exported as CSV and XLSX |

---

## Critical Methodological Notes

These fixes were applied during analysis and are documented in scripts 16b and 16b2:

| Item | Incorrect (initial) | Correct (final) |
|---|---|---|
| MaxEnt output | `type="cloglog"` (default) | `type="logistic"` — mandatory for valid ensemble |
| GAM selection | `MuMIn::dredge()` | `select=TRUE + method="REML"` |
| RF imbalance | `classwt` | Balanced bootstrap (`sampsize`) |
| XGBoost nrounds | Floor at 50 | Floor at 100 (noisy CV peak at round 28) |
| Dist_RawMat | All lithological units (mean=36m) | Selective knappable lithologies only (mean=6742m) |
| Ensemble comparison | Arbitrary 0.02 threshold | DeLong's test (pROC::roc.test) |
| Background extent | Unbounded | District boundary only |
| ODMAP citation | Feng et al. (2019) | Zurell et al. (2020) |

---

## Environmental Predictor Stack (Final — 13 variables)

| # | Variable | Type | Source | VIF |
|---|---|---|---|---|
| 1 | Elevation | Continuous | Cartosat-1 CartoDEM 30m | 1.22 |
| 2 | Aspect | Continuous | Derived: DEM | 1.00 |
| 3 | TRI | Continuous | Wilson et al. (2007) | 1.19* |
| 4 | TPI | Continuous | Weiss (2001) | 1.19 |
| 5 | Plan Curvature | Continuous | Derived: DEM | 1.18 |
| 6 | HAND | Continuous | WhiteboxTools v2.3 | 1.58 |
| 7 | Flow Accumulation (log10) | Continuous | D8 algorithm | 1.25 |
| 8 | Distance to Perennial River | Continuous | SOI 1:50,000 | 1.05 |
| 9 | Distance to Palaeochannel | Continuous | MNDWI ∩ DEM valley | 1.02 |
| 10 | Distance to Raw Material Source | Continuous | NGDR + field GPS | 1.04 |
| 11 | NDVI (dry season) | Continuous | Sentinel-2A May 2025 | 1.05 |
| 12 | Geology | Categorical | NGDR 1:50,000 | — |
| 13 | Geomorphology | Categorical | NGDR 1:50,000 | — |

*Slope excluded: VIF=74.82 (collinear with TRI). TRI retained as mechanistically more informative.

---

## How to Reproduce the Analysis

```r
# 1. Clone repository
# git clone https://github.com/sushantbegade/Palaeolithic-Ensemble-SDM-CentralIndia

# 2. Open project in RStudio
# Open: Nagpur_Chandrapur_Ensemble.Rproj

# 3. Restore exact package versions
renv::restore()

# 4. Download data from Zenodo
# Zenodo DOI: [INSERT AT SUBMISSION]
# Place in data_raw/ directory as specified in Script 02

# 5. Run scripts in order
source("scripts/01_setup_packages.R")
source("scripts/02_load_and_reproject.R")
# ... through to ...
source("scripts/25_tables_export.R")
```

> **Note:** All scripts use `set.seed(42)` before every stochastic operation. All rasters are aligned to `TEMPLATE_30m_utm44n.tif`. Runtime for full analysis: approximately 6–10 hours depending on hardware.

---

## Data Availability

All environmental predictor rasters, processed site coordinates, model prediction surfaces, SHAP rasters, and evaluation outputs are openly available at:

**Zenodo DOI:** [INSERT AT SUBMISSION]

**Raw data sources:**

| Dataset | Portal |
|---|---|
| Geological, Geomorphological, Lithological data | [National Geoscience Data Repository (NGDR)](https://geodataindia.gov.in) |
| River and Waterbody data | [India-WRIS](https://indiawris.gov.in) |
| Cartosat-1 CartoDEM 30m | [Bhoonidhi](https://bhoonidhi.nrsc.gov.in) |
| Sentinel-2A Level-2A imagery | [Bhoonidhi](https://bhoonidhi.nrsc.gov.in) |
| CHELSA-TraCE21k palaeoclimate | [CHELSA](https://chelsa-climate.org/chelsa-trace21k/) |

> Raw GPS field data are available from the corresponding author upon reasonable request, subject to site confidentiality constraints.

---

## Sub-Model Results

| Sub-model | N sites | Algorithms | Ens CV AUC | Boyce | Top SHAP driver |
|---|---|---|---|---|---|
| Lower Palaeolithic | 70 | 5 (GAM excluded) | 0.815 | 0.929 | Elevation (42.9%) |
| Middle Palaeolithic | 110 | 5 (GAM excluded) | 0.728 | 0.855 | TRI (43.4%) |
| Upper Palaeolithic | 74 | 5 (GAM excluded) | 0.785 | 0.946 | TRI (31.7%) |

GAM excluded from all sub-models due to REML timeout under extreme class imbalance at sub-model sample sizes.

Attribution sensitivity: LP rankings STABLE across confirmed-only vs confirmed+probable sites. MP Rank 3 changed (HAND→NDVI) — interpretive claims for MP limited to confirmed-site results.

---

## Preceding Paper

This ensemble paper is the direct methodological successor to:

>	Begade, S. 2026. A Predictive Model for the Identification of Open-Air Palaeolithic Sites in Nagpur and Chandrapur Districts, Maharashtra, India. Journal of Archaeological Science: Reports, 69, 105541. https://doi.org/10.1016/j.jasrep.2025.105541

---

## Citation

```
Begade, S. and Sahu, P., 2026. A Multi-Model Ensemble Framework with Spatial 
Explainability for Predicting Open-Air Palaeolithic Site Distribution in Central India: 
Bias Correction, Algorithmic Diversity, and Spatially Heterogeneous Predictor Importance. 
Journal of Archaeological Science. (Under Review)
```

---

## License

**Code:** [MIT License](LICENSE)

**Data:** [Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

---

## Reproducibility

- R version: 4.5.x (see Supplementary Table S3 for full `sessionInfo()`)
- Package versions locked in `renv.lock`
- ODMAP reporting checklist: Zurell et al. (2020) — Supplementary S1
- Pre-analysis design: Begade & Sahu (2026) Research Design v4.0 (available on request)