# Palaeolithic Ensemble SDM — Central India

**Paper:** A Multi-Model Ensemble Framework with Spatial Explainability 
for Predicting Open-Air Palaeolithic Site Distribution in Central India: 
Bias Correction, Algorithmic Diversity, and Spatially Heterogeneous 
Predictor Importance

**Author:** Sushant Begade and Prabash Sahu
**Affiliation:** Department of Ancient Indian History, Culture and Archaeology, Rashtrasant Tukadoji Maharaj Nagpur University, Nagpur – 440033, India.   
**ORCID:** Sushant Begade: 0009-0003-0804-1763; Prabash Sahu: 0000-0003-0691-0403
**Target journal:** Journal of Archaeological Science  

---

## Repository Contents

| Folder/File | Description |
|---|---|
| `scripts/01_setup_packages.R` | Install and load all required R packages |
| `scripts/02_load_and_reproject.R` | Load all raw data, reproject to UTM 44N |
| `scripts/03_terrain_derivatives.R` | DEM-derived predictors (slope, TRI, HAND etc.) |
| `scripts/04_sentinel_processing.R` | NDVI and MNDWI from Sentinel-2A |
| `scripts/05_palaeochannel.R` | Palaeochannel layer construction |
| `scripts/06_climate_downscaling.R` | CHELSA-TraCE21k delta downscaling to 30m |
| `scripts/07_predictor_stack_vif.R` | Assemble 14-variable predictor stack, VIF screen |
| `scripts/08_sites_thinning.R` | Spatial thinning at 1km, cultural period flags |
| `scripts/09_bias_correction_background.R` | KDE bias surface, 10,000 background points |
| `scripts/10_spatial_cv_design.R` | Spatial block cross-validation design |
| `scripts/11_model_maxent.R` | MaxEnt (maxnet, logistic output) |
| `scripts/12_model_rf.R` | Random Forest (class-weighted) |
| `scripts/13_model_xgboost.R` | XGBoost (regularised gradient boosting) |
| `scripts/14_model_brt.R` | Boosted Regression Trees |
| `scripts/15_model_gam.R` | Generalised Additive Models (select=TRUE) |
| `scripts/16_model_svm.R` | Support Vector Machine (Platt scaling) |
| `scripts/17_ensemble.R` | AUC-weighted ensemble + uncertainty surface |
| `scripts/18_uncertainty.R` | Ensemble SD surface and confidence zones |
| `scripts/19_shap_analysis.R` | SHAP values, dominant driver map |
| `scripts/20_submodels_LP_MP_UP.R` | Cultural period sub-models |
| `scripts/21_transfer_validation.R` | Geographic hold-out validation |
| `scripts/22_evaluation_metrics.R` | AUC, TSS, Boyce, Kvamme's Gain, DeLong |
| `scripts/23_figures_main.R` | Main manuscript figures 1–12 |
| `scripts/24_figures_supplementary.R` | Supplementary figures S1–S5 |
| `scripts/25_tables_export.R` | All tables exported as CSV |
| `renv.lock` | Exact R package versions for reproducibility |

---

## How to Reproduce the Analysis

1. Clone this repository
2. Open `Nagpur_Chandrapur_Ensemble.Rproj` in RStudio
3. Run `renv::restore()` to install exact package versions
4. Download data from Zenodo DOI: [INSERT DOI AT SUBMISSION]
5. Run scripts 01 through 25 in numbered order

---

## Data Availability

All environmental predictor rasters, processed site coordinates, 
model outputs, and SHAP surfaces are available at:  
**Zenodo DOI:** [INSERT AT SUBMISSION]

Raw data sources:
- Geological, Geomorphological and Lithological Data: geodataindia.gov.in
- River and Waterbody Data: indiawris.gov.in
- DEM and Sentinel-2A: bhoonidhi.nrsc.gov.in
- CHELSA-TraCE21k: chelsa-climate.org/chelsa-trace21k/
- WorldClim v2.1: worldclim.org

---

## Citation

Begade, S., [Year]. [Full paper title]. 
Journal of Archaeological Science. 
https://doi.org/[INSERT AT PUBLICATION]

---

## License

Code: MIT License  
Data: Creative Commons Attribution 4.0 (CC-BY 4.0)
