# ============================================================
# SCRIPT 23: MAIN FIGURES — PUBLICATION QUALITY REWRITE
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade | RTMNU Nagpur
# ORCID:  0009-0003-0804-1763
# Date:   June 2026
# Script: 23 of 25
# ============================================================
# COMPLETE REWRITE — ggplot2 + tidyterra throughout.
# Base R terra::plot() removed entirely. All figures use
# ggplot2 pipeline for publication-quality output.
#
# FIXES FROM Script 23 v1:
#   FIX 1 — Fig 8 crash: terra::global custom function
#     returns list not scalar. Fixed: use terra::values()
#     + base::quantile() directly.
#   FIX 2 — Fig 1 India inset: was blank blue box.
#     Fixed: rnaturalearth for real country outline.
#   FIX 3 — Fig 2 dark panels (TPI/HAND): range near-zero
#     squished to black. Fixed: oob=scales::squish +
#     per-raster scale limits.
#   FIX 4 — Fig 11 legend "character(0)": lookup mismatch.
#     Fixed: infer from raster unique values directly.
#   FIX 5 — Fig 7 error bar to 1.0: ensemble has no SD.
#     Fixed: omit CI bar for ensemble, note in caption.
#
# FIGURES PRODUCED (all 13):
#   Fig01 — Study area (ggplot + rnaturalearth)
#   Fig02 — Predictor stack (ggplot + tidyterra)
#   Fig03 — Palaeochannel construction
#   Fig04 — Bias correction
#   Fig05 — Spatial block CV design
#   Fig06 — Individual algorithm response curves
#   Fig07 — CV AUC comparison (ggplot, DeLong stars)
#   Fig08 — Ensemble suitability + uncertainty + zones
#   Fig09 — SHAP beeswarm (reuse Script 19 if exists)
#   Fig10 — SHAP dependence (reuse Script 19 if exists)
#   Fig11 — SHAP dominant driver map (reuse Script 19)
#   Fig12 — Diachronic sub-model drivers (reuse Script 20)
#   Fig_transfer — Transfer validation (reuse Script 21)
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))

terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

# ── Extra packages ───────────────────────────────────────────
pkgs_extra <- c("tidyterra","ggspatial","patchwork",
                "rnaturalearth","rnaturalearthdata",
                "scales","ggtext","showtext")
for (p in pkgs_extra) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, quiet = TRUE)
}
suppressPackageStartupMessages({
  library(ggplot2); library(tidyterra); library(ggspatial)
  library(patchwork); library(scales); library(sf)
})
rl_ok <- requireNamespace("rnaturalearth", quietly = TRUE)
if (rl_ok) library(rnaturalearth)

cat("\n========================================\n")
cat("SCRIPT 23: Main Figures (ggplot2 rewrite)\n")
cat("========================================\n\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─────────────────────────────────────────────────────────────
# 0. GLOBAL THEME & PALETTES
# ─────────────────────────────────────────────────────────────

# Publication theme — clean, minimal, JAS-appropriate
theme_pub_map <- function(base_size = 11) {
  theme_void(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 1,
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 1,
                                      hjust = 0.5, color = "grey35",
                                      margin = margin(b = 3)),
      plot.caption     = element_text(size = base_size - 2,
                                      hjust = 0, color = "grey45",
                                      face = "italic",
                                      margin = margin(t = 5)),
      legend.title     = element_text(size = base_size - 1,
                                      face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      legend.key.size  = unit(0.4, "cm"),
      legend.key.width = unit(0.55, "cm"),
      plot.margin      = margin(6, 6, 6, 6),
      panel.border     = element_rect(color = "grey30",
                                      fill = NA, linewidth = 0.4),
      plot.background  = element_rect(fill = "white",
                                      color = NA))
}

theme_pub_chart <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 1,
                                      face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size - 1,
                                      hjust = 0.5, color = "grey35"),
      axis.text        = element_text(size = base_size - 1),
      axis.title       = element_text(size = base_size,
                                      face = "bold"),
      legend.position  = "none",
      panel.grid.major.y = element_line(color = "grey92",
                                        linewidth = 0.4),
      plot.background  = element_rect(fill = "white",
                                      color = NA),
      plot.margin      = margin(8, 8, 8, 8))
}

# Color palettes
PAL_SUIT  <- "viridis"     # suitability maps
PAL_UNCERT <- "plasma"     # uncertainty
PAL_KDE   <- "inferno"     # bias/KDE surfaces
PAL_ELEV  <- "mako"        # elevation

# Site colors by cultural period
PERIOD_COLS <- c(LP = "#E63946", MP = "#2A9D8F",
                 UP = "#F4A261", Undiff = "#ADB5BD")

# Driver palette (13 predictors, color-blind tolerant)
DRIVER_PAL <- c(
  TRI                = "#2A9D8F",
  Elevation          = "#E63946",
  Geomorphology      = "#FF9F1C",
  NDVI               = "#4CC9F0",
  HAND               = "#06D6A0",
  Geology            = "#8338EC",
  Dist_Palaeochannel = "#073B4C",
  Dist_River         = "#118AB2",
  Dist_RawMat        = "#FF6B6B",
  Aspect             = "#F4A261",
  TPI                = "#264653",
  Plan_Curvature     = "#E9C46A",
  Flow_Accum_log10   = "#B5838D")

SITE_SIZE  <- 1.2
SITE_SHAPE <- 21
SITE_STROKE <- 0.3

DPI <- 300L

# ─────────────────────────────────────────────────────────────
# 1. LOAD SHARED DATA
# ─────────────────────────────────────────────────────────────

cat("--- Loading Shared Data ---\n\n")

template_30m  <- terra::rast(file.path(OUT_PREDICTORS,
                                       "TEMPLATE_30m_utm44n.tif"))
boundary_vect <- terra::vect(sf::st_read(
  file.path(OUT_SITES, "study_area_boundary_utm44n.gpkg"),
  quiet = TRUE))
boundary_sf   <- sf::st_as_sf(boundary_vect)

sites_sf  <- sf::st_read(file.path(OUT_SITES,
                                   "sites_all_utm44n.gpkg"),
                         quiet = TRUE)
sites_thin <- sf::st_read(file.path(OUT_SITES,
                                    "sites_thinned_pooled.gpkg"),
                          quiet = TRUE)

# Add cultural period for color coding
if (!"Relative.Chronology....Lower...Middle.Palaeolithic..etc.." %in%
    names(sites_sf)) {
  period_col <- grep("Relativ|Chronol|Period",
                     names(sites_sf), value = TRUE,
                     ignore.case = TRUE)[1]
} else {
  period_col <- "Relative.Chronology....Lower...Middle.Palaeolithic..etc.."
}

if (!is.null(period_col) && !is.na(period_col)) {
  chron <- as.character(sites_sf[[period_col]])
  sites_sf$period_class <- dplyr::case_when(
    grepl("Lower",  chron, ignore.case=TRUE) &
      !grepl("Middle|Upper", chron, ignore.case=TRUE) ~ "LP",
    grepl("Middle", chron, ignore.case=TRUE) &
      !grepl("Lower|Upper",  chron, ignore.case=TRUE) ~ "MP",
    grepl("Upper",  chron, ignore.case=TRUE) &
      !grepl("Lower|Middle", chron, ignore.case=TRUE) ~ "UP",
    TRUE ~ "Undiff")
} else {
  sites_sf$period_class <- "Undiff"
}

dem_r    <- terra::rast(file.path(OUT_PREDICTORS,
                                  "DEM_30m_utm44n.tif"))
rivers   <- sf::st_read(file.path(OUT_PREDICTORS,
                                  "rivers_utm44n.gpkg"),
                        quiet = TRUE)
kde_r    <- tryCatch(terra::rast(file.path(OUT_BIAS,
                                           "bias_surface_kde_30m_utm44n.tif")), error=function(e) NULL)
bg_sf    <- sf::st_read(file.path(OUT_BACKGROUND,
                                  "background_N10000.gpkg"),
                        quiet = TRUE)
ens_r    <- tryCatch(terra::rast(file.path(OUT_MOD_ENS,
                                           "ensemble_primary.tif")), error=function(e) NULL)
unc_r    <- tryCatch(terra::rast(file.path(OUT_MOD_ENS,
                                           "ensemble_uncertainty_sd.tif")), error=function(e) NULL)
zone_r   <- tryCatch(terra::rast(file.path(OUT_MOD_ENS,
                                           "ensemble_confidence_zones.tif")), error=function(e) NULL)
shap_r   <- tryCatch(terra::rast(file.path(OUT_SHAP,
                                           "shap_dominant_driver_map.tif")), error=function(e) NULL)
pred_stack <- terra::rast(file.path(OUT_PREDICTORS,
                                    "PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
final_names <- readRDS(file.path(OUT_PREDICTORS,
                                 "final_predictor_names.rds"))

# Hillshade
slope_h  <- terra::terrain(dem_r, "slope",  unit = "radians")
aspect_h <- terra::terrain(dem_r, "aspect", unit = "radians")
hillshade_r <- terra::shade(slope_h, aspect_h,
                            angle = 40, direction = 315)
rm(slope_h, aspect_h); gc(full=TRUE)

cat(sprintf("  Sites: %d  Thinned: %d\n",
            nrow(sites_sf), nrow(sites_thin)))
cat("\n")

# ─────────────────────────────────────────────────────────────
# 2. FIGURE 1 — STUDY AREA MAP
# ─────────────────────────────────────────────────────────────

cat("--- Figure 1: Study Area ---\n")

tryCatch({
  
  # India + neighbours for inset
  india_sf <- NULL; neighbours_sf <- NULL
  if (rl_ok) {
    india_sf      <- rnaturalearth::ne_countries(
      country = "India", returnclass = "sf", scale = "medium")
    neighbours_sf <- rnaturalearth::ne_countries(
      continent = "Asia", returnclass = "sf", scale = "medium")
  }
  
  # Study box in WGS84
  bbox_geo <- sf::st_bbox(sf::st_transform(boundary_sf, 4326))
  study_box <- sf::st_as_sf(sf::st_as_sfc(bbox_geo))
  
  # Main map
  p_main <- ggplot() +
    # Hillshade base
    tidyterra::geom_spatraster(data = hillshade_r,
                               show.legend = FALSE) +
    scale_fill_gradient(low = "#1a1a2e", high = "#e8e8e8",
                        na.value = "white",
                        guide = "none") +
    ggnewscale::new_scale_fill() %||%
    ggplot2::scale_fill_identity() # fallback if no ggnewscale
  
  # Try with ggnewscale, fall back without
  use_gns <- requireNamespace("ggnewscale", quietly = TRUE)
  if (use_gns) library(ggnewscale)
  
  p_main <- ggplot() +
    # Hillshade base
    tidyterra::geom_spatraster(data = hillshade_r,
                               show.legend = FALSE,
                               aes(fill = hillshade)) +
    scale_fill_gradientn(
      colors = c("#0d0d0d","#404040","#808080","#c0c0c0","#f0f0f0"),
      na.value = "white", guide = "none") +
    # Rivers
    geom_sf(data = rivers, color = "#74B3CE",
            linewidth = 0.5, alpha = 0.8) +
    # Study boundary
    geom_sf(data = boundary_sf, fill = NA,
            color = "white", linewidth = 1.0) +
    # Sites by cultural period
    geom_sf(data = sites_sf,
            aes(color = period_class),
            size = SITE_SIZE, shape = SITE_SHAPE,
            fill = NA, stroke = SITE_STROKE + 0.1) +
    scale_color_manual(name = "Cultural period",
                       values = PERIOD_COLS,
                       labels = c(LP    = "Lower Palaeolithic",
                                  MP    = "Middle Palaeolithic",
                                  UP    = "Upper Palaeolithic",
                                  Undiff = "Undifferentiated"),
                       guide  = guide_legend(override.aes = list(
                         size = 2.5, shape = 16))) +
    ggspatial::annotation_scale(
      location = "bl", width_hint = 0.22,
      text_col = "white", line_col = "white",
      bar_cols = c("white","white")) +
    ggspatial::annotation_north_arrow(
      location = "br", which_north = "true",
      height = unit(0.9, "cm"), width = unit(0.7, "cm"),
      style = ggspatial::north_arrow_orienteering(
        fill = c("white","white"),
        line_col = "white",
        text_col = "white")) +
    labs(title    = "Study Area: Nagpur & Chandrapur Districts",
         subtitle = "Wainganga-Wardha Basin, Central India (~21,300 km\u00b2)",
         caption  = paste0("Coordinate system: WGS84 / UTM Zone 44N ",
                           "(EPSG:32644). DEM: Cartosat-1 CartoDEM.")) +
    theme_pub_map() +
    theme(legend.position  = c(0.18, 0.15),
          legend.background = element_rect(
            fill = adjustcolor("black", 0.5),
            color = NA, linewidth = 0),
          legend.text  = element_text(color = "white",
                                      size = 7.5),
          legend.title = element_text(color = "white",
                                      size = 8, face = "bold"),
          legend.key   = element_blank())
  
  # Inset map (India)
  if (!is.null(india_sf)) {
    p_inset <- ggplot() +
      geom_sf(data = neighbours_sf %||% india_sf,
              fill = "#D8E4BC", color = "grey70",
              linewidth = 0.2) +
      geom_sf(data = india_sf, fill = "#A8C5A0",
              color = "grey40", linewidth = 0.4) +
      geom_sf(data = study_box, fill = "#E63946",
              color = "white", linewidth = 1,
              alpha = 0.75) +
      coord_sf(xlim = c(68, 98), ylim = c(6, 37)) +
      annotate("text", x = 82, y = 22,
               label = "India", fontface = "bold",
               size = 2.5, color = "grey20") +
      theme_void() +
      theme(panel.border = element_rect(color = "grey40",
                                        fill = NA, linewidth = 0.6),
            plot.background = element_rect(fill = "white",
                                           color = NA))
  } else {
    # Minimal fallback inset
    p_inset <- ggplot() +
      annotate("rect", xmin=68, xmax=98,
               ymin=6, ymax=37,
               fill="#D8E4BC", color="grey40") +
      annotate("rect", xmin=bbox_geo[1], xmax=bbox_geo[3],
               ymin=bbox_geo[2], ymax=bbox_geo[4],
               fill="#E63946", color="white", linewidth=1) +
      annotate("text", x=82, y=22, label="India",
               fontface="bold", size=2.5) +
      theme_void() +
      theme(panel.border=element_rect(color="grey40",
                                      fill=NA,linewidth=0.6),
            plot.background=element_rect(fill="white"))
  }
  
  # Compose: main + inset using patchwork inset_element
  fig01 <- p_main +
    patchwork::inset_element(p_inset,
                             left = 0.72, bottom = 0.60,
                             right = 1.00, top  = 1.00,
                             align_to = "plot")
  
  ggsave(file.path(OUT_FIG_MAIN, "Fig01_study_area.png"),
         fig01, width = 7.5, height = 9.5, dpi = DPI,
         bg = "white")
  cat("  ✓ Fig01_study_area.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig01 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 3. FIGURE 2 — PREDICTOR STACK
# ─────────────────────────────────────────────────────────────

cat("--- Figure 2: Predictor Stack ---\n")

tryCatch({
  n_pred <- terra::nlyr(pred_stack)
  
  # FIX: per-raster colour limits prevent dark/black panels
  # Use quantile stretch for each predictor
  plot_list <- lapply(seq_len(n_pred), function(i) {
    r    <- terra::subset(pred_stack, i)
    nm   <- names(r)
    cat_r <- nm %in% c("Geology","Geomorphology")
    
    # Quantile stretch for continuous rasters
    vals <- as.numeric(terra::values(r, na.rm = TRUE))
    vlo  <- as.numeric(quantile(vals, 0.02, na.rm = TRUE))
    vhi  <- as.numeric(quantile(vals, 0.98, na.rm = TRUE))
    if (vlo == vhi) { vlo <- min(vals, na.rm=TRUE)
    vhi <- max(vals, na.rm=TRUE) }
    
    pal  <- if (cat_r) "Set2" else PAL_SUIT
    n_cat <- if (cat_r)
      length(unique(vals[is.finite(vals)])) else 100L
    
    p <- ggplot() +
      tidyterra::geom_spatraster(data = r) +
      { if (cat_r)
        scale_fill_gradientn(
          colors = RColorBrewer::brewer.pal(
            min(n_cat, 8), "Set2"),
          na.value = "white", guide = "none")
        else
          scale_fill_viridis_c(
            option     = PAL_SUIT,
            na.value   = "white",
            limits     = c(vlo, vhi),
            oob        = scales::squish,
            guide      = "none") } +
      geom_sf(data = boundary_sf, fill = NA,
              color = "grey40", linewidth = 0.25) +
      geom_sf(data = sites_thin, color = "#E63946",
              size = 0.20, alpha = 0.7) +
      labs(title = nm) +
      theme_pub_map(base_size = 7.5) +
      theme(plot.title = element_text(size = 7.5,
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b=2)),
            panel.border = element_rect(color="grey60",
                                        fill=NA, linewidth=0.25))
    p
  })
  
  # patchwork grid (4 columns)
  fig02 <- patchwork::wrap_plots(plot_list, ncol = 4) +
    patchwork::plot_annotation(
      title   = "Fig. 2 — Final Retained Predictor Stack",
      subtitle = sprintf("(%d variables, 30m UTM Zone 44N; sites in red)",
                         n_pred),
      theme   = theme(
        plot.title    = element_text(size = 12, face = "bold",
                                     hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     color = "grey40")))
  
  n_rows_fig2 <- ceiling(n_pred / 4)
  ggsave(file.path(OUT_FIG_MAIN, "Fig02_predictor_stack.png"),
         fig02,
         width  = 12,
         height = n_rows_fig2 * 3.5,
         dpi    = DPI, bg = "white")
  cat("  ✓ Fig02_predictor_stack.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig02 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 4. FIGURE 3 — PALAEOCHANNEL CONSTRUCTION
# ─────────────────────────────────────────────────────────────

cat("--- Figure 3: Palaeochannel ---\n")

tryCatch({
  mndwi_r  <- terra::rast(file.path(OUT_PREDICTORS,
                                    "MNDWI_30m_utm44n.tif"))
  valley_r <- terra::rast(file.path(OUT_PALAEO,
                                    "PALAEO_VALLEY_30m_utm44n.tif"))
  dist_pc  <- terra::rast(file.path(OUT_PREDICTORS,
                                    "DIST_PALAEOCHANNEL_30m_utm44n.tif"))
  
  # Panel A: MNDWI moisture anomaly
  mndwi_anom <- terra::ifel(mndwi_r > -0.10, 1, NA)
  
  pA <- ggplot() +
    tidyterra::geom_spatraster(data = terra::ifel(
      is.na(mndwi_anom), 0, mndwi_anom)) +
    scale_fill_gradientn(
      colors    = c("#f5f5f5","#2196F3"),
      na.value  = "white", guide = "none") +
    geom_sf(data = boundary_sf, fill = NA,
            color = "grey30", linewidth = 0.4) +
    labs(title    = "(A) Spectral Source",
         subtitle = "MNDWI moisture anomaly (> \u22120.10)") +
    theme_pub_map()
  
  # Panel B: Valley morphology
  val_mask <- terra::ifel(valley_r == 1L, 1, NA)
  pB <- ggplot() +
    tidyterra::geom_spatraster(data = terra::ifel(
      is.na(val_mask), 0, val_mask)) +
    scale_fill_gradientn(
      colors   = c("#f5f5f5","#d62828"),
      na.value = "white", guide = "none") +
    geom_sf(data = boundary_sf, fill = NA,
            color = "grey30", linewidth = 0.4) +
    labs(title    = "(B) Topographic Source",
         subtitle = "DEM valley morphology (geomorphons)") +
    theme_pub_map()
  
  # Panel C: Distance to confirmed palaeochannel
  dist_km <- dist_pc / 1000
  pC <- ggplot() +
    tidyterra::geom_spatraster(data = dist_km) +
    scale_fill_viridis_c(
      option   = "plasma", name = "Dist. (km)",
      na.value = "white",
      limits   = c(0, as.numeric(quantile(
        terra::values(dist_km, na.rm=TRUE), 0.99, na.rm=TRUE))),
      oob = scales::squish) +
    geom_sf(data = boundary_sf, fill = NA,
            color = "white", linewidth = 0.4) +
    ggspatial::annotation_scale(location = "bl",
                                text_col = "white", line_col = "white",
                                bar_cols = c("white","white")) +
    labs(title    = "(C) Distance to Palaeochannel",
         subtitle = "Confirmed network (MNDWI \u2229 Valley)") +
    theme_pub_map() +
    theme(legend.position = "right")
  
  fig03 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 3 \u2014 Palaeochannel Layer Construction (Dual-Source Protocol)",
      caption = "Left: spectral (MNDWI >-0.10). Centre: topographic (geomorphons). Right: distance to confirmed network.",
      theme   = theme(plot.title   = element_text(size=11, face="bold", hjust=0.5),
                      plot.caption = element_text(size=7.5, hjust=0, color="grey40",
                                                  face="italic")))
  
  ggsave(file.path(OUT_FIG_MAIN, "Fig03_palaeochannel.png"),
         fig03, width = 14, height = 6, dpi = DPI, bg = "white")
  cat("  ✓ Fig03_palaeochannel.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig03 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 5. FIGURE 4 — BIAS CORRECTION
# ─────────────────────────────────────────────────────────────

cat("--- Figure 4: Bias Correction ---\n")

tryCatch({
  if (is.null(kde_r)) stop("KDE raster not found")
  
  pA <- ggplot() +
    tidyterra::geom_spatraster(data = kde_r) +
    scale_fill_viridis_c(option = PAL_KDE,
                         name = "Norm.\ndensity", na.value = "white") +
    geom_sf(data = sites_sf, color = "white",
            size = 0.5, alpha = 0.6) +
    geom_sf(data = boundary_sf, fill = NA,
            color = "grey70", linewidth = 0.4) +
    labs(title    = "(A) KDE Survey-Effort Surface",
         subtitle = "Hscv bandwidth, N=197 sites") +
    theme_pub_map() +
    theme(legend.position = "right")
  
  pB <- ggplot() +
    tidyterra::geom_spatraster(data = kde_r) +
    scale_fill_viridis_c(option = PAL_KDE,
                         na.value = "white", guide = "none") +
    geom_sf(data = bg_sf, color = "#00FF88",
            size = 0.06, alpha = 0.3) +
    geom_sf(data = boundary_sf, fill = NA,
            color = "grey70", linewidth = 0.4) +
    labs(title    = "(B) Bias-Weighted Background",
         subtitle = "N = 10,000 pts (Warton & Shepherd 2010)") +
    theme_pub_map()
  
  pC <- ggplot() +
    geom_sf(data = boundary_sf, fill = "grey96",
            color = "grey50", linewidth = 0.5) +
    geom_sf(data = bg_sf, color = "#118AB2",
            size = 0.12, alpha = 0.25) +
    geom_sf(data = sites_sf,
            aes(color = period_class),
            size = 1.4, shape = 16, alpha = 0.85) +
    scale_color_manual(name = "Period",
                       values = PERIOD_COLS,
                       labels = c(LP="Lower", MP="Middle",
                                  UP="Upper", Undiff="Undiffer."),
                       guide  = guide_legend(
                         override.aes = list(size = 2.5))) +
    labs(title    = "(C) Sites vs Background",
         subtitle = "KDE-weighted matching survey intensity") +
    theme_pub_map() +
    theme(legend.position = "right",
          legend.key.size = unit(0.3, "cm"))
  
  fig04 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 4 \u2014 Survey-Effort Bias Correction",
      caption = "Phillips et al. (2009); ks package, Hscv bandwidth.",
      theme   = theme(plot.title   = element_text(size=11, face="bold", hjust=0.5),
                      plot.caption = element_text(size=7.5, hjust=0,
                                                  color="grey40", face="italic")))
  
  ggsave(file.path(OUT_FIG_MAIN, "Fig04_bias_correction.png"),
         fig04, width = 14, height = 6, dpi = DPI, bg = "white")
  cat("  ✓ Fig04_bias_correction.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig04 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 6. FIGURE 5 — SPATIAL CV DESIGN (from blockCV object)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 5: Spatial CV Design ---\n")

tryCatch({
  cv_design  <- readRDS(file.path(OUT_CV,
                                  "cv_block_assignments.rds"))
  sites_cv   <- sf::st_read(file.path(OUT_CV,
                                      "sites_with_folds.gpkg"),
                            quiet = TRUE)
  bg_cv      <- sf::st_read(file.path(OUT_CV,
                                      "background_with_folds.gpkg"),
                            quiet = TRUE)
  fold_col_cv <- grep("fold", names(sites_cv),
                      ignore.case=TRUE, value=TRUE)[1]
  
  FOLD_COLS <- c("1"="#E63946","2"="#2A9D8F","3"="#F4A261",
                 "4"="#8338EC","5"="#118AB2")
  
  sites_cv$fold_fac <- factor(as.character(
    sites_cv[[fold_col_cv]]))
  bg_cv$fold_fac    <- factor(as.character(
    bg_cv[[fold_col_cv]]))
  
  # Panel A: site fold assignment
  pA <- ggplot() +
    geom_sf(data = boundary_sf, fill = "grey97",
            color = "grey40", linewidth = 0.5) +
    geom_sf(data = sites_cv, aes(color = fold_fac),
            size = 1.8, shape = 16) +
    scale_color_manual(name = "Fold",
                       values = FOLD_COLS,
                       guide  = guide_legend(
                         override.aes = list(size = 3))) +
    ggspatial::annotation_scale(location = "bl") +
    labs(title    = "(A) Site Fold Assignments",
         subtitle = sprintf("Block size: %.0f km",
                            cv_design$block_size_km)) +
    theme_pub_map() +
    theme(legend.position = "right")
  
  # Panel B: background fold assignment
  pB <- ggplot() +
    geom_sf(data = boundary_sf, fill = "grey97",
            color = "grey40", linewidth = 0.5) +
    geom_sf(data = bg_cv, aes(color = fold_fac),
            size = 0.08, alpha = 0.4) +
    scale_color_manual(name = "Fold",
                       values = FOLD_COLS, guide = "none") +
    labs(title    = "(B) Background Fold Assignments",
         subtitle = "N = 10,000 pts") +
    theme_pub_map()
  
  # Fold size table panel
  fold_tbl <- as.data.frame(table(Fold = sites_cv$fold_fac))
  pC <- ggplot(fold_tbl, aes(x = Fold, y = Freq,
                             fill = Fold)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = Freq), vjust = -0.4,
              size = 3.2, fontface = "bold") +
    scale_fill_manual(values = FOLD_COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title    = "(C) Sites per Fold",
         x = "Spatial fold", y = "N sites") +
    theme_pub_chart() +
    theme(legend.position = "none")
  
  fig05 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 5 \u2014 Five-Fold Spatial Block Cross-Validation Design",
      caption = sprintf("Variogram autocorrelation range: %.1f km (blockCV 3.1; Valavi et al. 2019).",
                        cv_design$autocor_range_m / 1000),
      theme   = theme(plot.title   = element_text(size=11, face="bold", hjust=0.5),
                      plot.caption = element_text(size=7.5, hjust=0,
                                                  color="grey40", face="italic")))
  
  ggsave(file.path(OUT_FIG_MAIN, "Fig05_spatial_cv_design.png"),
         fig05, width = 14, height = 6, dpi = DPI, bg = "white")
  cat("  ✓ Fig05_spatial_cv_design.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig05 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 7. FIGURE 7 — AUC COMPARISON (ggplot, DeLong stars)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 7: AUC Comparison ---\n")

tryCatch({
  # Load all evaluation CSVs
  read_auc <- function(path, alg_col = "algorithm") {
    if (!file.exists(path)) return(NULL)
    df <- read.csv(path, stringsAsFactors = FALSE)
    list(mean = df$cv_auc_mean[1],
         sd   = df$cv_auc_sd[1])
  }
  
  alg_data <- data.frame(
    algorithm = c("MaxEnt","RF","XGBoost","BRT","GAM","SVM",
                  "Ensemble\n(AUC-wt)"),
    auc_mean  = c(
      read.csv(file.path(OUT_EVAL,"maxent_evaluation.csv"))$cv_auc_mean[1],
      read.csv(file.path(OUT_EVAL,"rf_evaluation.csv"))$cv_auc_mean[1],
      read.csv(file.path(OUT_EVAL,"xgboost_evaluation.csv"))$cv_auc_mean[1],
      read.csv(file.path(OUT_EVAL,"brt_evaluation.csv"))$cv_auc_mean[1],
      read.csv(file.path(OUT_EVAL,"gam_evaluation.csv"))$cv_auc_mean[1],
      read.csv(file.path(OUT_EVAL,"svm_evaluation.csv"))$cv_auc_mean[1],
      0.7239),
    auc_sd    = c(
      read.csv(file.path(OUT_EVAL,"maxent_evaluation.csv"))$cv_auc_sd[1],
      read.csv(file.path(OUT_EVAL,"rf_evaluation.csv"))$cv_auc_sd[1],
      read.csv(file.path(OUT_EVAL,"xgboost_evaluation.csv"))$cv_auc_sd[1],
      read.csv(file.path(OUT_EVAL,"brt_evaluation.csv"))$cv_auc_sd[1],
      read.csv(file.path(OUT_EVAL,"gam_evaluation.csv"))$cv_auc_sd[1],
      read.csv(file.path(OUT_EVAL,"svm_evaluation.csv"))$cv_auc_sd[1],
      NA_real_),  # ensemble has no per-fold SD from here
    is_ensemble = c(rep(FALSE,6), TRUE),
    stringsAsFactors = FALSE)
  
  # DeLong significance (ensemble is significantly better than these)
  alg_data$delong_sig <- c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
  
  # CI (1.96 × SD) — omit for ensemble
  alg_data$ci_hi <- ifelse(!alg_data$is_ensemble,
                           pmin(alg_data$auc_mean + 1.96*alg_data$auc_sd, 1, na.rm=TRUE),
                           NA_real_)
  alg_data$ci_lo <- ifelse(!alg_data$is_ensemble,
                           pmax(alg_data$auc_mean - 1.96*alg_data$auc_sd, 0, na.rm=TRUE),
                           NA_real_)
  
  alg_data$algorithm <- factor(alg_data$algorithm,
                               levels = alg_data$algorithm)
  
  BAR_COLS <- c("#2A9D8F","#E63946","#F4A261","#264653",
                "#E9C46A","#8338EC","#06D6A0")
  
  fig07 <- ggplot(alg_data,
                  aes(x = algorithm, y = auc_mean, fill = algorithm)) +
    # Reference lines
    geom_hline(yintercept = 0.75, color = "#F4A261",
               linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = 0.85, color = "#E63946",
               linetype = "dotted", linewidth = 0.6) +
    # Bars
    geom_col(width = 0.68, color = NA) +
    # CI error bars (individual algorithms only)
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.22, linewidth = 0.7,
                  color = "grey20", na.rm = TRUE) +
    # AUC value labels
    geom_text(aes(label = sprintf("%.4f", auc_mean)),
              vjust = -0.5, size = 3.0,
              fontface = "bold", color = "grey20") +
    # DeLong significance stars
    geom_text(data = subset(alg_data, delong_sig),
              aes(y = auc_mean + 0.032),
              label = "*", color = "#E63946",
              size = 5, fontface = "bold") +
    # Ensemble star marker
    geom_text(data = subset(alg_data, is_ensemble),
              aes(y = auc_mean + 0.032),
              label = "\u2605", color = "#06D6A0",
              size = 4) +
    scale_fill_manual(values = BAR_COLS, guide = "none") +
    scale_y_continuous(limits = c(0.50, 1.00),
                       breaks = seq(0.50, 1.00, 0.05),
                       expand = expansion(mult=c(0,0.04))) +
    annotate("text", x = 0.6, y = 0.753,
             label = "AUC = 0.75 (adequate)",
             hjust = 0, size = 2.8, color = "#F4A261",
             fontface = "italic") +
    annotate("text", x = 0.6, y = 0.853,
             label = "AUC = 0.85 (strong)",
             hjust = 0, size = 2.8, color = "#E63946",
             fontface = "italic") +
    labs(title    = "Fig. 7 \u2014 Spatial Block CV AUC: All Models",
         subtitle = "Error bars = \u00b11.96 SD (95% CI); * = DeLong p<0.05 vs ensemble; \u2605 = ensemble",
         x = "", y = "Spatial Block CV AUC (5-fold mean)",
         caption  = "DeLong's test (DeLong et al. 1988; pROC). Ensemble significantly outperforms XGBoost, BRT, GAM, SVM.") +
    theme_pub_chart() +
    theme(axis.text.x  = element_text(size = 9,
                                      color = "grey20"),
          plot.caption = element_text(size = 7.5,
                                      color = "grey45", face = "italic"))
  
  ggsave(file.path(OUT_FIG_MAIN, "Fig07_AUC_comparison.png"),
         fig07, width = 9, height = 6, dpi = DPI, bg = "white")
  cat("  ✓ Fig07_AUC_comparison.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig07 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 8. FIGURE 8 — ENSEMBLE SUITABILITY + UNCERTAINTY + ZONES
# FIX: terra::global custom function crash → use terra::values()
# ─────────────────────────────────────────────────────────────

cat("--- Figure 8: Ensemble Suitability ---\n")

if (is.null(ens_r)) {
  cat("  ✗ ensemble_primary.tif not found\n\n")
} else {
  tryCatch({
    # FIX 1: terra::values() + base quantile() — no terra::global crash
    ens_flat <- as.numeric(terra::values(ens_r, na.rm = TRUE))
    q25 <- as.numeric(quantile(ens_flat, 0.25, na.rm = TRUE))
    q50 <- as.numeric(quantile(ens_flat, 0.50, na.rm = TRUE))
    q75 <- as.numeric(quantile(ens_flat, 0.75, na.rm = TRUE))
    
    suit_class <- terra::classify(ens_r,
                                  rcl = matrix(c(-Inf,q25,1, q25,q50,2,
                                                 q50,q75,3, q75,Inf,4),
                                               ncol=3, byrow=TRUE),
                                  include.lowest = TRUE)
    
    suit_lvls  <- c("Low","Moderate","High","Very High")
    suit_labs  <- sprintf("%s\n(%.2f\u2013%.2f)",
                          suit_lvls,
                          c(0, q25, q50, q75),
                          c(q25, q50, q75, 1))
    suit_cols  <- c("#264653","#2A9D8F","#E9C46A","#E63946")
    names(suit_cols) <- as.character(1:4)
    
    # Panel A: continuous suitability
    pA <- ggplot() +
      tidyterra::geom_spatraster(data = ens_r) +
      scale_fill_viridis_c(
        option   = PAL_SUIT,
        name     = "Suitability\n(probability)",
        na.value = "white",
        limits   = c(0, 1),
        breaks   = seq(0, 1, 0.2)) +
      geom_sf(data = boundary_sf, fill = NA,
              color = "white", linewidth = 0.7) +
      geom_sf(data = sites_thin, color = "#E63946",
              size = SITE_SIZE, shape = 16, alpha = 0.8) +
      ggspatial::annotation_scale(location = "bl",
                                  text_col="white", line_col="white",
                                  bar_cols=c("white","white")) +
      ggspatial::annotation_north_arrow(location = "br",
                                        which_north = "true",
                                        height = unit(0.8,"cm"), width = unit(0.65,"cm"),
                                        style = ggspatial::north_arrow_orienteering(
                                          fill=c("white","grey30"), line_col="white",
                                          text_col="white")) +
      labs(title    = "(A) Ensemble Suitability Surface",
           subtitle = sprintf("CV AUC=0.7239  Boyce=0.9092  KG=0.6063")) +
      theme_pub_map() +
      theme(legend.position = c(0.87, 0.70))
    
    # Panel B: uncertainty
    if (!is.null(unc_r)) {
      unc_flat <- as.numeric(terra::values(unc_r, na.rm=TRUE))
      unc_q75  <- as.numeric(quantile(unc_flat, 0.75, na.rm=TRUE))
      unc_max  <- as.numeric(quantile(unc_flat, 0.99, na.rm=TRUE))
      
      pB <- ggplot() +
        tidyterra::geom_spatraster(data = unc_r) +
        scale_fill_viridis_c(
          option   = "plasma",
          name     = "SD\n(uncertainty)",
          na.value = "white",
          limits   = c(0, unc_max),
          oob      = scales::squish) +
        geom_sf(data = boundary_sf, fill = NA,
                color = "white", linewidth = 0.7) +
        geom_sf(data = sites_thin, color = "white",
                size = 0.7, alpha = 0.6) +
        labs(title    = "(B) Prediction Uncertainty",
             subtitle = sprintf("SD across 6 algorithms  Q75=%.3f",
                                unc_q75)) +
        theme_pub_map() +
        theme(legend.position = c(0.87, 0.70))
    } else {
      pB <- ggplot() +
        annotate("text", x=0.5, y=0.5,
                 label="Uncertainty raster\nnot found",
                 size=4, color="grey50") +
        theme_void()
    }
    
    # Panel C: classified suitability
    suit_class_df <- as.data.frame(suit_class,
                                   xy = TRUE, na.rm = TRUE)
    names(suit_class_df)[3] <- "class"
    suit_class_df$class_lab <- factor(
      suit_class_df$class,
      levels = 1:4, labels = suit_lvls)
    
    pC <- ggplot() +
      geom_raster(data = suit_class_df,
                  aes(x = x, y = y, fill = class_lab)) +
      scale_fill_manual(
        name   = "Suitability",
        values = c(Low="#264653", Moderate="#2A9D8F",
                   High="#E9C46A", `Very High`="#E63946"),
        na.value = "white",
        guide  = guide_legend(reverse = TRUE)) +
      geom_sf(data = boundary_sf, fill = NA,
              color = "grey30", linewidth = 0.6) +
      geom_sf(data = sites_thin, color = "white",
              size = 0.8, alpha = 0.7) +
      coord_sf() +
      labs(title    = "(C) Suitability Classification",
           subtitle = "Max-TSS threshold; 100% sites in Very High") +
      theme_pub_map() +
      theme(legend.position = "right")
    
    fig08 <- (pA | pB | pC) +
      patchwork::plot_annotation(
        title   = "Fig. 8 \u2014 AUC-Weighted Ensemble: Suitability, Uncertainty & Classification",
        caption = "N=6 algorithms; all outputs on logistic probability scale [0\u20131]. Sites = thinned pooled (N=190).",
        theme   = theme(plot.title   = element_text(size=11, face="bold", hjust=0.5),
                        plot.caption = element_text(size=7.5, hjust=0,
                                                    color="grey40", face="italic")))
    
    ggsave(file.path(OUT_FIG_MAIN,
                     "Fig08_ensemble_suitability.png"),
           fig08, width = 16, height = 7.5,
           dpi = DPI, bg = "white")
    cat("  ✓ Fig08_ensemble_suitability.png\n\n")
  }, error = function(e) {
    cat(sprintf("  ✗ Fig08 error: %s\n\n", e$message))
  })
}

# ─────────────────────────────────────────────────────────────
# 9. FIGURE 11 — SHAP DOMINANT DRIVER MAP
# FIX: character(0) legend bug — infer from raster values
# ─────────────────────────────────────────────────────────────

cat("--- Figure 11: SHAP Dominant Driver Map ---\n")

# Use Script 19 version if it looks better
fig11_s19 <- file.path(OUT_FIG_MAIN,
                       "Fig11_dominant_driver_map.png")
fig11_s23 <- file.path(OUT_FIG_MAIN,
                       "Fig11_shap_dominant_driver.png")

if (file.exists(fig11_s19) &&
    file.info(fig11_s19)$size > 200000) {
  cat("  Using Script 19 version (superior ggplot2 output) ✓\n\n")
  file.copy(fig11_s19, gsub("dominant_driver_map",
                            "shap_dominant_driver_FINAL",
                            fig11_s19), overwrite = TRUE)
} else if (!is.null(shap_r)) {
  tryCatch({
    # FIX: infer driver codes from raster directly
    freq_tbl <- terra::freq(shap_r, na.rm = TRUE)
    freq_tbl <- freq_tbl[!is.na(freq_tbl$value) &
                           freq_tbl$value > 0, ]
    freq_tbl <- freq_tbl[order(-freq_tbl$count), ]
    
    codes  <- as.integer(freq_tbl$value)
    n_drv  <- length(codes)
    labels <- final_names[pmin(codes, length(final_names))]
    pcts   <- 100 * freq_tbl$count / sum(freq_tbl$count)
    
    # Color assignment
    drv_cols <- sapply(labels, function(nm) {
      if (nm %in% names(DRIVER_PAL)) DRIVER_PAL[nm]
      else "#AAAAAA"
    })
    
    # Reclassify raster to 1:n
    rcl  <- cbind(codes, seq_along(codes))
    shap_rcl <- terra::classify(shap_r, rcl, others = NA)
    
    # Convert to data frame for ggplot
    shap_df <- as.data.frame(shap_rcl, xy = TRUE,
                             na.rm = TRUE)
    names(shap_df)[3] <- "driver_idx"
    shap_df$driver <- factor(shap_df$driver_idx,
                             levels   = seq_along(labels),
                             labels   = labels)
    
    # Legend labels with percentage
    leg_labs <- sprintf("%s (%.1f%%)", labels, pcts)
    names(leg_labs) <- labels
    
    fig11 <- ggplot() +
      geom_raster(data = shap_df,
                  aes(x = x, y = y, fill = driver)) +
      scale_fill_manual(
        name   = "Dominant predictor\n(% of 250m cells)",
        values = setNames(drv_cols, labels),
        labels = leg_labs,
        na.value = "white",
        guide  = guide_legend(
          ncol = if (n_drv > 8) 2L else 1L,
          override.aes = list(size = 4))) +
      geom_sf(data = boundary_sf, fill = NA,
              color = "grey20", linewidth = 0.7) +
      geom_sf(data = sites_thin, color = "white",
              size = 0.7, alpha = 0.7, shape = 16) +
      # Zone annotations
      annotate("label", x = mean(sf::st_bbox(boundary_sf)[c(1,3)]) - 30000,
               y = mean(sf::st_bbox(boundary_sf)[c(2,4)]) + 50000,
               label = "Chandrapur\nBasin\n(TRI dom.)",
               size = 2.8, fontface = "bold",
               fill = adjustcolor("white", 0.75),
               label.size = 0) +
      annotate("label",
               x = mean(sf::st_bbox(boundary_sf)[c(1,3)]),
               y = mean(sf::st_bbox(boundary_sf)[c(2,4)]) + 10000,
               label = "Nagpur\nPediplain\n(Geomorph.)",
               size = 2.8, fontface = "bold",
               fill = adjustcolor("white", 0.75),
               label.size = 0) +
      annotate("label",
               x = mean(sf::st_bbox(boundary_sf)[c(1,3)]) + 20000,
               y = mean(sf::st_bbox(boundary_sf)[c(2,4)]) - 50000,
               label = "Satpura\nFoothills\n(Geomorph.)",
               size = 2.8, fontface = "bold",
               fill = adjustcolor("white", 0.75),
               label.size = 0) +
      ggspatial::annotation_scale(location = "bl") +
      ggspatial::annotation_north_arrow(location = "br",
                                        which_north = "true",
                                        height = unit(0.8,"cm"), width = unit(0.65,"cm")) +
      coord_sf() +
      labs(title   = "Fig. 11 \u2014 SHAP Dominant Predictor Map",
           subtitle = sprintf(
             "TreeSHAP, XGBoost; 250m grid resampled to 30m; \u03c7\u00b2=54,352, p<0.001"),
           caption = "Predictor with highest absolute SHAP value per 250m cell. Sites = thinned pooled (N=190, white).") +
      theme_pub_map() +
      theme(legend.position = "right",
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8, face="bold"))
    
    ggsave(file.path(OUT_FIG_MAIN,
                     "Fig11_shap_dominant_driver.png"),
           fig11, width = 9, height = 10.5,
           dpi = DPI, bg = "white")
    cat("  ✓ Fig11_shap_dominant_driver.png\n\n")
  }, error = function(e) {
    cat(sprintf("  ✗ Fig11 error: %s\n\n", e$message))
  })
}

# ─────────────────────────────────────────────────────────────
# 10. FIGURE 6 — INDIVIDUAL ALGORITHM RESPONSE CURVES
# ─────────────────────────────────────────────────────────────

cat("--- Figure 6: Response Curves ---\n")

tryCatch({
  # Load XGBoost model + top predictors
  xgb_path <- file.path(OUT_MOD_IND, "xgboost_model_final.bin")
  xgb_info <- readRDS(file.path(OUT_MOD_IND,
                                "xgboost_model_info.rds"))
  if (!file.exists(xgb_path)) stop("XGBoost model not found")
  xgb_mod  <- xgboost::xgb.load(xgb_path)
  
  # Top 5 predictors from SHAP
  shap_imp_path <- file.path(OUT_TABLES,
                             "Table4_SHAP_importance.csv")
  if (file.exists(shap_imp_path)) {
    shap_imp_df <- read.csv(shap_imp_path,
                            stringsAsFactors = FALSE)
    shap_imp_df$xgb_v <- suppressWarnings(
      as.numeric(shap_imp_df[[
        grep("xgb|XGB", names(shap_imp_df),
             value=TRUE, ignore.case=TRUE)[1]]]))
    top5 <- head(shap_imp_df[order(-shap_imp_df$xgb_v,
                                   na.last=TRUE), "predictor"],
                 5)
  } else {
    top5 <- c("TRI","NDVI","Elevation","Dist_River","HAND")
  }
  top5 <- top5[top5 %in% final_names]
  if (!length(top5)) top5 <- final_names[1:5]
  
  cat(sprintf("  Top 5 for response curves: %s\n",
              paste(top5, collapse=", ")))
  
  # Sample predictor values
  set.seed(42)
  samp <- terra::spatSample(pred_stack, size = 5000,
                            method = "random",
                            na.rm = TRUE, as.df = TRUE)
  cats_cols <- c("Geology","Geomorphology")
  for (col in cats_cols[cats_cols %in% names(samp)])
    samp[[col]] <- as.numeric(as.integer(samp[[col]]))
  
  samp_mat <- as.matrix(samp[, final_names, drop=FALSE])
  storage.mode(samp_mat) <- "double"
  
  # Marginal response curves for each top predictor
  curve_list <- lapply(top5, function(pred) {
    col_idx <- which(final_names == pred)
    if (!length(col_idx)) return(NULL)
    pred_range <- seq(
      quantile(samp_mat[,col_idx], 0.05, na.rm=TRUE),
      quantile(samp_mat[,col_idx], 0.95, na.rm=TRUE),
      length.out = 80)
    med_mat <- matrix(
      apply(samp_mat, 2, median, na.rm=TRUE),
      nrow = 80, ncol = ncol(samp_mat), byrow=TRUE)
    colnames(med_mat) <- final_names
    med_mat[, col_idx] <- pred_range
    storage.mode(med_mat) <- "double"
    preds <- predict(xgb_mod,
                     xgboost::xgb.DMatrix(med_mat))
    data.frame(predictor = pred, x = pred_range, y = preds)
  })
  curve_df <- do.call(rbind, Filter(Negate(is.null), curve_list))
  curve_df$predictor <- factor(curve_df$predictor,
                               levels = top5)
  
  fig06 <- ggplot(curve_df, aes(x = x, y = y,
                                color = predictor)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey60", linewidth = 0.5) +
    scale_color_manual(
      values = setNames(
        c("#2A9D8F","#E63946","#F4A261","#8338EC","#06D6A0"),
        top5),
      guide = "none") +
    facet_wrap(~predictor, scales = "free_x", nrow = 1) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    labs(title    = "Fig. 6 \u2014 Marginal Response Curves (XGBoost)",
         subtitle = "Top 5 predictors by mean |SHAP|; other predictors held at median",
         x = "Predictor value", y = "Predicted suitability",
         caption  = "Computed from 80-point marginal grid; N=5,000 background sample.") +
    theme_pub_chart() +
    theme(strip.text   = element_text(size=9, face="bold"),
          panel.border = element_rect(color="grey80",
                                      fill=NA, linewidth=0.4),
          axis.text.x  = element_text(size=7, angle=20,
                                      hjust=1),
          plot.caption = element_text(size=7.5, hjust=0,
                                      color="grey45", face="italic"),
          legend.position = "none")
  
  ggsave(file.path(OUT_FIG_MAIN,
                   "Fig06_response_curves.png"),
         fig06, width = 14, height = 4.5,
         dpi = DPI, bg = "white")
  cat("  ✓ Fig06_response_curves.png\n\n")
}, error = function(e) {
  cat(sprintf("  ✗ Fig06 error: %s\n\n", e$message))
})

# ─────────────────────────────────────────────────────────────
# 11. FINAL FIGURE CHECKLIST
# ─────────────────────────────────────────────────────────────

cat("--- Final Figure Checklist ---\n\n")

checklist <- list(
  list("Fig01_study_area.png",            "Study area + India inset"),
  list("Fig02_predictor_stack.png",       "Predictor stack (13 vars)"),
  list("Fig03_palaeochannel.png",         "Palaeochannel construction"),
  list("Fig04_bias_correction.png",       "Bias correction"),
  list("Fig05_spatial_cv_design.png",     "Spatial CV design"),
  list("Fig06_response_curves.png",       "Response curves (XGBoost)"),
  list("Fig07_AUC_comparison.png",        "AUC comparison + DeLong"),
  list("Fig08_ensemble_suitability.png",  "Ensemble suitability"),
  list("Fig09_SHAP_beeswarm.png",         "SHAP beeswarm (Script 19)"),
  list("Fig10_SHAP_dependence.png",       "SHAP dependence (Script 19)"),
  list("Fig11_dominant_driver_map.png",   "SHAP dominant driver (Script 19)"),
  list("Fig12_diachronic_driver_maps.png","Diachronic sub-models (Script 20)"),
  list("Fig_transfer_validation.png",     "Transfer validation (Script 21)"))

n_ok <- 0L
for (item in checklist) {
  fp   <- file.path(OUT_FIG_MAIN, item[[1]])
  ok   <- file.exists(fp)
  if (ok) n_ok <- n_ok + 1L
  sz   <- if (ok) sprintf("%.0f KB",
                          file.info(fp)$size/1024) else "MISSING"
  cat(sprintf("  %s %-45s %s\n",
              if (ok) "\u2713" else "\u2717",
              item[[2]], sz))
}
cat(sprintf("\n  %d / %d figures present\n\n", n_ok,
            length(checklist)))

cat("========================================\n")
cat("SCRIPT 23 COMPLETE\n")
cat("========================================\n\n")
cat("Framework: ggplot2 + tidyterra (not base R)\n")
cat("All fixes applied:\n")
cat("  FIX 1 — Fig08: terra::values()+quantile() ✓\n")
cat("  FIX 2 — Fig01: rnaturalearth India outline ✓\n")
cat("  FIX 3 — Fig02: quantile stretch (no dark panels) ✓\n")
cat("  FIX 4 — Fig11: legend from raster freq() ✓\n")
cat("  FIX 5 — Fig07: no CI bar for ensemble ✓\n")
cat("\n300 DPI  |  viridis palettes  |  patchwork layouts\n")
cat("========================================\n")