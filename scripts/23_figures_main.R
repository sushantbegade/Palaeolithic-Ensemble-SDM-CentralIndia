# ============================================================
# SCRIPT 23: ALL 12 MAIN FIGURES — BUG-FIXED
# ============================================================
# Paper: A Multi-Model Ensemble Framework with Spatial
#        Explainability for Predicting Open-Air Palaeolithic
#        Site Distribution in Central India
# Author: Sushant Begade* and Prabash Sahu
# ORCID:  0009-0003-0804-1763; 0000-0003-0691-0403
# Script: 23 of 25
# ============================================================
# Fixes applied vs previous version:
#  [1] Fig01 fill conflict (DEM continuous + sites discrete)
#      → sites now use aes(COLOR=period) NOT fill; DEM keeps fill
#  [2] Fig01 India inset REMOVED (user request)
#  [3] ALL maps: UTM 44N retained (NO WGS84 transform; user request)
#  [4] Fig07 "Discrete to continuous": geom_hline() replaces
#      annotate("segment",..,x=-Inf..); x=integer positions in annotate
#  [5] Fig07 type safety: as.numeric() on all adf columns
#  [6] Fig11/12: terra::freq() — na.rm=TRUE removed (not valid arg)
#  [7] Fig09/10: SHAP RDS absent → rebuild from XGBoost model inline
# ============================================================

source(file.path(
  "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset",
  "R_Analysis/Palaeolithic-Ensemble-SDM-CentralIndia/scripts",
  "01_setup_packages.R"
))
terra::terraOptions(tempdir = "E:/R_temp", memmax = 4)

for (p in c("tidyterra","ggspatial","patchwork","scales","shapviz")) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, quiet = TRUE)
}
suppressPackageStartupMessages({
  library(ggplot2); library(tidyterra); library(ggspatial)
  library(patchwork); library(scales); library(sf); library(dplyr)
})
HAS_SHAP <- requireNamespace("shapviz", quietly = TRUE)

cat("\n========================================\n")
cat("SCRIPT 23: All 12 Main Figures\n")
cat("========================================\n\n")
`%||%` <- function(a, b) if (!is.null(a)) a else b
DPI <- 300L

# ─────────────────────────────────────────────────────────────
# 0. DESIGN SYSTEM
# ─────────────────────────────────────────────────────────────

TERRAIN_PAL <- c(
  "#004529","#006d2c","#238b45","#41ae76","#66c2a4",
  "#99d8c9","#ccece6","#ffffcc","#ffeda0","#fed976",
  "#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#800026")

# Cultural period: use COLOUR (not fill) on map to avoid DEM fill conflict
PERIOD_COLS <- c(LP="#1f77b4",MP="#d62728",UP="#2ca02c",Undiff="#9467bd")
PERIOD_LABS <- c(LP="Lower Palaeolithic",MP="Middle Palaeolithic",
                 UP="Upper Palaeolithic",Undiff="Undifferentiated")

FOLD_COLS  <- c("1"="#E63946","2"="#2ca02c","3"="#F4A261","4"="#8338EC","5"="#1f77b4")
ALG_COLS   <- c(MaxEnt="#2A9D8F",RF="#E63946",XGBoost="#F4A261",
                BRT="#264653",GAM="#E9C46A",SVM="#8338EC",Ensemble="#06D6A0")
DRIVER_PAL <- c(
  TRI="#2A9D8F",Elevation="#E63946",Geomorphology="#FF9F1C",
  NDVI="#4CC9F0",HAND="#06D6A0",Geology="#8338EC",
  Dist_Palaeochannel="#073B4C",Dist_River="#118AB2",
  Dist_RawMat="#FF6B6B",Aspect="#F4A261",TPI="#264653",
  Plan_Curvature="#E9C46A",Flow_Accum_log10="#B5838D",
  BIO1="#C77DFF",BIO12="#9B5DE5",BIO15="#F15BB5")

# Map theme: bw base, coord axis labels in UTM
theme_map_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace% theme(
    plot.title       = element_text(size=base_size+1, face="bold",
                                    hjust=0.5, margin=margin(b=4)),
    plot.subtitle    = element_text(size=base_size-1, hjust=0.5,
                                    color="grey35", margin=margin(b=3)),
    plot.caption     = element_text(size=7.5, hjust=0, color="grey45",
                                    face="italic", margin=margin(t=5)),
    axis.title       = element_text(size=base_size-1, face="bold"),
    axis.text        = element_text(size=base_size-2, color="grey25"),
    panel.grid.major = element_line(color="white", linewidth=0.5),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size=base_size-1, face="bold"),
    legend.text      = element_text(size=base_size-2),
    legend.key.size  = unit(0.38,"cm"),
    legend.background= element_rect(fill="white",color="grey70",linewidth=0.3),
    plot.background  = element_rect(fill="white",color=NA),
    plot.margin      = margin(5,5,5,5))
}

theme_chart_pub <- function(base_size = 11) {
  theme_bw(base_size=base_size) %+replace% theme(
    plot.title       = element_text(size=base_size+1,face="bold",hjust=0.5),
    plot.subtitle    = element_text(size=base_size-1,hjust=0.5,
                                    color="grey35",margin=margin(b=4)),
    plot.caption     = element_text(size=7.5,hjust=0,color="grey45",
                                    face="italic",margin=margin(t=5)),
    axis.title       = element_text(size=base_size-1,face="bold"),
    axis.text        = element_text(size=base_size-2,color="grey20"),
    panel.grid.major = element_line(color="grey93",linewidth=0.4),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size=base_size-1,face="bold"),
    strip.background = element_rect(fill="grey96",color="grey70"),
    legend.title     = element_text(size=base_size-1,face="bold"),
    legend.text      = element_text(size=base_size-2),
    plot.background  = element_rect(fill="white",color=NA),
    plot.margin      = margin(6,8,6,8))
}

# UTM axis label formatter (m → readable)
utm_fmt <- function(x) {
  ifelse(x >= 1e6, paste0(round(x/1e3), "k"), format(round(x), big.mark=","))
}

# ─────────────────────────────────────────────────────────────
# 1. LOAD SHARED DATA (all in UTM 44N — NO WGS84 transform)
# ─────────────────────────────────────────────────────────────

cat("--- Loading Shared Data ---\n\n")

# Boundaries — keep in UTM
boundary_utm <- sf::st_as_sf(terra::vect(sf::st_read(
  file.path(OUT_SITES,"study_area_boundary_utm44n.gpkg"),quiet=TRUE)))

# Sites
sites_sf   <- sf::st_read(file.path(OUT_SITES,"sites_all_utm44n.gpkg"),quiet=TRUE)
sites_thin <- sf::st_read(file.path(OUT_SITES,"sites_thinned_pooled.gpkg"),quiet=TRUE)

# Cultural period
per_col <- grep("Relativ|Chronol",names(sites_sf),
                value=TRUE,ignore.case=TRUE)[1]
if (!is.na(per_col %||% NA)) {
  ch <- as.character(sites_sf[[per_col]])
  sites_sf$period <- dplyr::case_when(
    grepl("Lower",ch,ignore.case=T)&!grepl("Middle|Upper",ch,ignore.case=T)~"LP",
    grepl("Middle",ch,ignore.case=T)&!grepl("Lower|Upper",ch,ignore.case=T)~"MP",
    grepl("Upper",ch,ignore.case=T)&!grepl("Lower|Middle",ch,ignore.case=T)~"UP",
    TRUE~"Undiff")
} else sites_sf$period <- "Undiff"
sites_sf$period <- factor(sites_sf$period, levels=names(PERIOD_COLS))

# Rivers — UTM, major only (>5 km)
rivers_utm_all <- sf::st_read(file.path(OUT_PREDICTORS,"rivers_utm44n.gpkg"),quiet=TRUE)
river_len      <- as.numeric(sf::st_length(rivers_utm_all))
rivers_utm     <- rivers_utm_all[river_len > 5000, ]
cat(sprintf("  Rivers (major >5 km): %d\n", nrow(rivers_utm)))

# DEM — aggregate to ~120m for display speed, keep in UTM
dem_r   <- terra::rast(file.path(OUT_PREDICTORS,"DEM_30m_utm44n.tif"))
dem_agg <- terra::aggregate(dem_r, fact=4, fun="mean")
names(dem_agg) <- "elevation"
dem_range <- as.numeric(terra::global(dem_agg, c("min","max"), na.rm=TRUE))

# Predictor stack
pred_stack  <- terra::rast(file.path(OUT_PREDICTORS,"PREDICTOR_STACK_FINAL_30m_utm44n.tif"))
final_names <- readRDS(file.path(OUT_PREDICTORS,"final_predictor_names.rds"))

# Ensemble rasters — aggregate for display
ens_r   <- tryCatch(terra::rast(file.path(OUT_MOD_ENS,"ensemble_primary.tif")),error=function(e)NULL)
unc_r   <- tryCatch(terra::rast(file.path(OUT_MOD_ENS,"ensemble_uncertainty_sd.tif")),error=function(e)NULL)
shap_r  <- tryCatch(terra::rast(file.path(OUT_SHAP,"shap_dominant_driver_map.tif")),error=function(e)NULL)
kde_r   <- tryCatch(terra::rast(file.path(OUT_BIAS,"bias_surface_kde_30m_utm44n.tif")),error=function(e)NULL)
bg_sf   <- sf::st_read(file.path(OUT_BACKGROUND,"background_N10000.gpkg"),quiet=TRUE)

ens_agg <- if (!is.null(ens_r)) {
  r <- terra::aggregate(ens_r, fact=3, fun="mean"); names(r) <- "suitability"; r
} else NULL
unc_agg <- if (!is.null(unc_r)) {
  r <- terra::aggregate(unc_r, fact=3, fun="mean"); names(r) <- "uncertainty"; r
} else NULL

cat(sprintf("  Sites: %d | Thinned: %d\n",nrow(sites_sf),nrow(sites_thin)))
cat("  Data loaded ✓\n\n")

# ─────────────────────────────────────────────────────────────
# 2. FIGURE 1 — STUDY AREA (UTM 44N, NO India inset)
# FIX: sites use aes(COLOR) not fill → no conflict with DEM fill
# ─────────────────────────────────────────────────────────────

cat("--- Figure 1: Study Area ---\n")

tryCatch({
  fig01 <- ggplot() +
    # DEM terrain — uses FILL scale (continuous)
    tidyterra::geom_spatraster(data = dem_agg,
                               show.legend = TRUE) +
    scale_fill_gradientn(
      name   = "Elevation\n(m a.s.l.)",
      colors = TERRAIN_PAL,
      limits = c(dem_range[1], dem_range[2]),
      na.value = "white",
      breaks = seq(150, 550, 100),
      guide  = guide_colorbar(
        barwidth  = unit(0.5,"cm"),
        barheight = unit(4.0,"cm"),
        ticks.colour = "grey30")) +
    # Major rivers
    geom_sf(data  = rivers_utm,
            color = "#2171b5",
            linewidth = 0.30,
            alpha = 0.70) +
    # Boundary
    geom_sf(data  = boundary_utm,
            fill  = NA,
            color = "grey15",
            linewidth = 0.80) +
    # Sites — use COLOR aesthetic (NOT fill) → NO CONFLICT with DEM fill
    geom_sf(data = sites_sf,
            aes(color = period),
            shape = 19,
            size  = 2.2,
            alpha = 0.88) +
    scale_color_manual(
      name   = "Cultural period",
      values = PERIOD_COLS,
      labels = PERIOD_LABS,
      guide  = guide_legend(
        override.aes = list(size=3.0,shape=19,alpha=1))) +
    ggspatial::annotation_scale(
      location   = "bl",
      width_hint = 0.22,
      text_cex   = 0.72,
      pad_x      = unit(0.4,"cm"),
      pad_y      = unit(0.4,"cm")) +
    ggspatial::annotation_north_arrow(
      location    = "tr",
      which_north = "true",
      pad_x       = unit(0.4,"cm"),
      pad_y       = unit(0.5,"cm"),
      height      = unit(1.0,"cm"),
      width       = unit(0.80,"cm"),
      style       = ggspatial::north_arrow_fancy_orienteering(
        text_size = 8)) +
    coord_sf(expand = FALSE) +   # UTM 44N — no transform
    scale_x_continuous(labels = utm_fmt) +
    scale_y_continuous(labels = utm_fmt) +
    labs(x = "Easting (m)", y = "Northing (m)",
         title    = "Study Area: Nagpur & Chandrapur Districts",
         subtitle = "Wainganga-Wardha Basin, Central India (~21,300 km\u00b2)  |  UTM Zone 44N",
         caption  = "DEM: Cartosat-1 CartoDEM (30m). Rivers: SOI 1:50,000 (major only, >5 km). N = 197 sites.") +
    theme_map_pub() +
    theme(legend.position = "right",
          legend.key      = element_blank(),
          legend.spacing.y = unit(0.4,"cm"))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig01_study_area.png"),
         fig01, width=8.5, height=10.5, dpi=DPI, bg="white")
  cat("  ✓ Fig01_study_area.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig01: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 3. FIGURE 2 — PREDICTOR STACK (13 retained, UTM)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 2: Predictor Stack ---\n")

tryCatch({
  n_pred   <- terra::nlyr(pred_stack)
  cat_pred <- c("Geology","Geomorphology")
  
  plist <- lapply(seq_len(n_pred), function(i) {
    r_agg <- terra::aggregate(terra::subset(pred_stack,i),
                              fact=4, fun="mean")
    nm  <- names(terra::subset(pred_stack,i))
    names(r_agg) <- "value"
    is_cat <- nm %in% cat_pred
    vals <- as.numeric(terra::values(r_agg, na.rm=TRUE))
    vlo  <- as.numeric(quantile(vals,0.01,na.rm=TRUE))
    vhi  <- as.numeric(quantile(vals,0.99,na.rm=TRUE))
    if (!is.finite(vlo)||vlo==vhi) { vlo<-min(vals,na.rm=T); vhi<-max(vals,na.rm=T) }
    
    ggplot() +
      tidyterra::geom_spatraster(data=r_agg) +
      { if (is_cat)
        scale_fill_viridis_d(na.value="white", guide="none")
        else
          scale_fill_viridis_c(
            option = if(nm %in% c("Elevation","HAND","TPI")) "mako"
            else if(nm=="NDVI") "viridis"
            else if(grepl("Dist",nm)) "plasma" else "viridis",
            limits=c(vlo,vhi), oob=scales::squish,
            na.value="white", guide="none") } +
      geom_sf(data=boundary_utm, fill=NA,
              color="grey30", linewidth=0.30) +
      # Sites as color (no fill conflict)
      geom_sf(data=sites_thin, color="#d62728",
              shape=19, size=0.40, alpha=0.70) +
      coord_sf(expand=FALSE) +
      labs(title=nm) +
      theme_void() +
      theme(
        plot.title   = element_text(size=8.5,face="bold",hjust=0.5,margin=margin(b=1)),
        panel.border = element_rect(fill=NA,color="grey55",linewidth=0.35),
        plot.background = element_rect(fill="white",color=NA),
        plot.margin  = margin(2,2,2,2))
  })
  
  fig02 <- patchwork::wrap_plots(plist, ncol=4L) +
    patchwork::plot_annotation(
      title    = "Fig. 2 \u2014 Final Retained Predictor Stack",
      subtitle = sprintf("%d predictors (post-VIF screening, VIF<5); 30m UTM 44N; red = sites (N=190)", n_pred),
      theme=theme(
        plot.title    = element_text(size=13,face="bold",hjust=0.5),
        plot.subtitle = element_text(size=9,hjust=0.5,color="grey40"),
        plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig02_predictor_stack.png"),
         fig02, width=13, height=ceiling(n_pred/4)*3.6,
         dpi=DPI, bg="white")
  cat("  ✓ Fig02_predictor_stack.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig02: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 4. FIGURE 3 — PALAEOCHANNEL (UTM, 3 panels)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 3: Palaeochannel ---\n")

tryCatch({
  mndwi_r  <- terra::rast(file.path(OUT_PREDICTORS,"MNDWI_30m_utm44n.tif"))
  valley_r <- terra::rast(file.path(OUT_PALAEO,"PALAEO_VALLEY_30m_utm44n.tif"))
  dist_pc  <- terra::rast(file.path(OUT_PREDICTORS,"DIST_PALAEOCHANNEL_30m_utm44n.tif"))
  
  agg3 <- function(r, nm, fun="mean") {
    rr <- terra::aggregate(r, fact=3, fun=fun)
    names(rr) <- nm; rr
  }
  mndwi_agg <- agg3(mndwi_r,  "MNDWI")
  valley_agg <- agg3(valley_r, "valley")
  dist_agg   <- agg3(dist_pc/1000, "dist_km")
  mv <- as.numeric(terra::values(mndwi_agg, na.rm=TRUE))
  
  base_layers <- list(
    coord_sf(expand=FALSE),
    scale_x_continuous(labels=utm_fmt),
    scale_y_continuous(labels=utm_fmt),
    labs(x=NULL,y=NULL))
  
  pA <- ggplot() +
    geom_sf(data=boundary_utm, fill="grey98",
            color="grey30", linewidth=0.5) +
    tidyterra::geom_spatraster(data=mndwi_agg) +
    scale_fill_gradient2(
      name="MNDWI",low="#8B0000",mid="grey96",high="#08519c",
      midpoint=-0.10,
      limits=c(as.numeric(quantile(mv,.01,na.rm=T)),
               as.numeric(quantile(mv,.99,na.rm=T))),
      oob=scales::squish, na.value="white") +
    geom_sf(data=boundary_utm, fill=NA, color="grey25", linewidth=0.5) +
    base_layers +
    labs(title="(A) Spectral Source",
         subtitle="MNDWI (Sentinel-2A, May 2025)") +
    theme_map_pub() + theme(legend.position="right")
  
  pB <- ggplot() +
    tidyterra::geom_spatraster(data=dem_agg) +
    scale_fill_gradientn(colors=TERRAIN_PAL, name=NULL,
                         guide="none", na.value="white") +
    # Valley overlay — use NEW fill only if valley is binary 0/1
    tidyterra::geom_spatraster(
      data=terra::ifel(valley_agg > 0.5, 1L, NA_integer_),
      aes(), show.legend=FALSE) +
    scale_fill_gradientn(colors="#d62728", na.value="transparent",
                         guide="none") +
    geom_sf(data=boundary_utm, fill=NA, color="grey25", linewidth=0.5) +
    base_layers +
    labs(title="(B) Topographic Source",
         subtitle="Valley morphology (geomorphons; DEM background)") +
    theme_map_pub()
  
  dv <- as.numeric(terra::values(dist_agg, na.rm=TRUE))
  pC <- ggplot() +
    tidyterra::geom_spatraster(data=dist_agg) +
    scale_fill_viridis_c(
      option="inferno", direction=-1, name="Dist. (km)",
      limits=c(0, as.numeric(quantile(dv,.99,na.rm=T))),
      oob=scales::squish, na.value="white") +
    geom_sf(data=boundary_utm, fill=NA, color="grey70", linewidth=0.5) +
    ggspatial::annotation_scale(location="bl") +
    ggspatial::annotation_north_arrow(location="tr",
                                      which_north="true",height=unit(0.75,"cm"),width=unit(0.60,"cm")) +
    base_layers +
    labs(title="(C) Distance to Palaeochannel",
         subtitle="Confirmed network: MNDWI \u2229 Valley morphology") +
    theme_map_pub() + theme(legend.position="right")
  
  fig03 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 3 \u2014 Palaeochannel Layer Construction (Dual-Source Protocol)",
      caption = "Sentinel-2A L2A, May 2025 (dry season). WhiteboxTools v2.3 geomorphons (search radius 10\u201330 cells). Threshold: MNDWI > \u22120.10.",
      theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5),
                  plot.caption=element_text(size=7.5,hjust=0,color="grey45",face="italic"),
                  plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig03_palaeochannel.png"),
         fig03, width=15, height=7, dpi=DPI, bg="white")
  cat("  ✓ Fig03_palaeochannel.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig03: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 5. FIGURE 4 — BIAS CORRECTION (UTM, 3 panels)
# FIX: sites use color, bg points = steel blue not neon
# ─────────────────────────────────────────────────────────────

cat("--- Figure 4: Bias Correction ---\n")

tryCatch({
  if (is.null(kde_r)) stop("KDE raster not found")
  kde_agg <- terra::aggregate(kde_r, fact=4, fun="mean")
  names(kde_agg) <- "kde"
  base_map <- list(
    geom_sf(data=boundary_utm,fill=NA,color="grey50",linewidth=0.5),
    coord_sf(expand=FALSE),
    scale_x_continuous(labels=utm_fmt),
    scale_y_continuous(labels=utm_fmt),
    labs(x=NULL,y=NULL))
  
  pA <- ggplot() +
    tidyterra::geom_spatraster(data=kde_agg) +
    scale_fill_viridis_c(option="inferno",name="Norm.\ndensity",
                         na.value="white") +
    geom_sf(data=sites_sf, color="white",
            shape=19, size=0.9, alpha=0.75) +
    base_map +
    ggspatial::annotation_scale(location="bl") +
    labs(title="(A) KDE Survey-Effort Surface",
         subtitle="Hscv bandwidth; N = 197 sites") +
    theme_map_pub() + theme(legend.position="right")
  
  pB <- ggplot() +
    tidyterra::geom_spatraster(data=kde_agg) +
    scale_fill_viridis_c(option="inferno",guide="none",na.value="white") +
    geom_sf(data=bg_sf, color="#bdc9d6",
            shape=19, size=0.10, alpha=0.40) +
    base_map +
    labs(title="(B) Bias-Weighted Background",
         subtitle="N = 10,000 pts (Warton & Shepherd 2010)") +
    theme_map_pub()
  
  pC <- ggplot() +
    geom_sf(data=boundary_utm,fill="grey97",color="grey30",linewidth=0.5) +
    geom_sf(data=bg_sf, color="#abd9e9",
            shape=19, size=0.12, alpha=0.30) +
    geom_sf(data=sites_sf, aes(color=period),
            shape=19, size=2.0, alpha=0.88) +
    scale_color_manual(name="Period",values=PERIOD_COLS,
                       labels=PERIOD_LABS,
                       guide=guide_legend(override.aes=list(size=2.5,shape=19))) +
    base_map +
    labs(title="(C) Sites vs Background",
         subtitle="Background intensity tracks survey effort") +
    theme_map_pub() +
    theme(legend.position="right",
          legend.text=element_text(size=7.5),
          legend.key.size=unit(0.32,"cm"))
  
  fig04 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 4 \u2014 Survey-Effort Bias Correction",
      caption = "Target-group approach (Phillips et al. 2009). ks package, Hscv bandwidth. Background sampled within district boundary.",
      theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5),
                  plot.caption=element_text(size=7.5,hjust=0,color="grey45",face="italic"),
                  plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig04_bias_correction.png"),
         fig04, width=15, height=7, dpi=DPI, bg="white")
  cat("  ✓ Fig04_bias_correction.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig04: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 6. FIGURE 5 — SPATIAL CV DESIGN (UTM)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 5: Spatial CV Design ---\n")

tryCatch({
  cv_design <- readRDS(file.path(OUT_CV,"cv_block_assignments.rds"))
  sites_cv  <- sf::st_read(file.path(OUT_CV,"sites_with_folds.gpkg"),quiet=TRUE)
  bg_cv     <- sf::st_read(file.path(OUT_CV,"background_with_folds.gpkg"),quiet=TRUE)
  fold_col  <- grep("fold",names(sites_cv),ignore.case=TRUE,value=TRUE)[1]
  sites_cv$fold_f <- factor(as.character(sites_cv[[fold_col]]))
  bg_cv$fold_f    <- factor(as.character(bg_cv[[fold_col]]))
  fold_tbl  <- as.data.frame(table(Fold=sites_cv$fold_f))
  
  pA <- ggplot() +
    geom_sf(data=boundary_utm, fill="grey98",
            color="grey30", linewidth=0.5) +
    geom_sf(data=bg_cv, aes(color=fold_f),
            size=0.08, alpha=0.35) +
    geom_sf(data=sites_cv, aes(color=fold_f),
            shape=19, size=2.2, alpha=0.90) +
    scale_color_manual(name="Fold",values=FOLD_COLS) +
    ggspatial::annotation_scale(location="bl") +
    ggspatial::annotation_north_arrow(location="tr",
                                      which_north="true",height=unit(0.8,"cm"),width=unit(0.65,"cm")) +
    guides(color=guide_legend(
      override.aes=list(size=2.8,shape=19,alpha=1))) +
    coord_sf(expand=FALSE) +
    scale_x_continuous(labels=utm_fmt) +
    scale_y_continuous(labels=utm_fmt) +
    labs(title="(A) Five-Fold Spatial Block Assignment",
         subtitle=sprintf("Block size: %.0f km | blockCV 3.1",
                          cv_design$block_size_km),
         x=NULL, y=NULL) +
    theme_map_pub() + theme(legend.position="right")
  
  pB <- ggplot(fold_tbl, aes(x=Fold, y=Freq, fill=Fold)) +
    geom_col(width=0.65, color=NA) +
    geom_text(aes(label=Freq), vjust=-0.4,
              size=3.5, fontface="bold", color="grey20") +
    scale_fill_manual(values=FOLD_COLS, guide="none") +
    scale_y_continuous(limits=c(0,max(fold_tbl$Freq)*1.18),
                       expand=expansion(mult=c(0,0))) +
    labs(title="(B) Sites per Fold",
         x="Spatial fold", y="N sites") +
    theme_chart_pub()
  
  fig05 <- (pA | pB) + patchwork::plot_layout(widths=c(2.5,1)) +
    patchwork::plot_annotation(
      title   = "Fig. 5 \u2014 Five-Fold Spatial Block Cross-Validation Design",
      caption = sprintf("Variogram range: %.0f km. blockCV 3.1 (Valavi et al. 2019).",
                        cv_design$autocor_range_m/1000),
      theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5),
                  plot.caption=element_text(size=7.5,hjust=0,color="grey45",face="italic"),
                  plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig05_spatial_cv_design.png"),
         fig05, width=12, height=7, dpi=DPI, bg="white")
  cat("  ✓ Fig05_spatial_cv_design.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig05: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 7. FIGURE 6 — RESPONSE CURVES (LOESS-smoothed)
# ─────────────────────────────────────────────────────────────

cat("--- Figure 6: Response Curves ---\n")

tryCatch({
  xgb_path <- file.path(OUT_MOD_IND,"xgboost_model_final.bin")
  if (!file.exists(xgb_path)) stop("XGBoost model not found")
  xgb_mod <- xgboost::xgb.load(xgb_path)
  
  xgb_imp  <- read.csv(file.path(OUT_EVAL,"xgboost_importance.csv"),
                       stringsAsFactors=FALSE)
  names(xgb_imp)[1] <- "predictor"
  gain_col <- grep("Gain|gain",names(xgb_imp),value=TRUE)[1]
  if (is.na(gain_col)) gain_col <- names(xgb_imp)[2]
  top5 <- head(xgb_imp$predictor[order(-xgb_imp[[gain_col]])], 5)
  top5 <- top5[top5 %in% final_names]
  if (!length(top5)) top5 <- final_names[1:5]
  
  set.seed(42)
  samp <- terra::spatSample(pred_stack, size=5000,
                            method="random", na.rm=TRUE, as.df=TRUE)
  for (col in c("Geology","Geomorphology")[c("Geology","Geomorphology") %in% names(samp)])
    samp[[col]] <- as.numeric(as.integer(samp[[col]]))
  samp_mat <- as.matrix(samp[,final_names,drop=FALSE])
  storage.mode(samp_mat) <- "double"
  
  units_map <- c(TRI="TRI",NDVI="NDVI",Elevation="Elevation (m)",
                 Dist_River="Dist. river (m)",HAND="HAND (m)",
                 Geomorphology="Geomorphology",Dist_Palaeochannel="Dist. palaeochannel (m)",
                 Dist_RawMat="Dist. raw material (m)",Aspect="Aspect (\u00b0)",
                 Plan_Curvature="Plan curvature",TPI="TPI",
                 Flow_Accum_log10="Flow accumulation (log\u2081\u2080)",Geology="Geology")
  
  curve_df <- do.call(rbind, lapply(top5, function(pred) {
    ci <- which(final_names==pred); if (!length(ci)) return(NULL)
    x_seq <- seq(as.numeric(quantile(samp_mat[,ci],0.02,na.rm=T)),
                 as.numeric(quantile(samp_mat[,ci],0.98,na.rm=T)),
                 length.out=200)
    med_m <- matrix(apply(samp_mat,2,median,na.rm=T),
                    nrow=200,ncol=ncol(samp_mat),byrow=TRUE)
    colnames(med_m) <- final_names
    med_m[,ci] <- x_seq; storage.mode(med_m) <- "double"
    preds <- predict(xgb_mod, xgboost::xgb.DMatrix(med_m))
    data.frame(predictor=pred, x=x_seq, y=preds)
  }))
  curve_df$pred_label <- sapply(as.character(curve_df$predictor),
                                function(p) units_map[p] %||% p)
  ordered_labs <- unique(curve_df$pred_label[
    match(top5, curve_df$predictor)])
  curve_df$pred_label <- factor(curve_df$pred_label, levels=ordered_labs)
  curve_cols <- setNames(
    c("#2A9D8F","#E63946","#F4A261","#8338EC","#06D6A0")[seq_along(top5)],
    levels(curve_df$pred_label))
  
  fig06 <- ggplot(curve_df, aes(x=x,y=y,color=pred_label)) +
    geom_line(color="grey80", linewidth=0.5) +
    geom_smooth(method="loess", formula=y~x, span=0.35,
                se=TRUE, linewidth=1.4, alpha=0.15) +
    geom_hline(yintercept=0.5, linetype="dashed",
               color="grey40", linewidth=0.6) +
    scale_color_manual(values=curve_cols, guide="none") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
    facet_wrap(~pred_label, scales="free_x", nrow=1) +
    labs(title    = "Fig. 6 \u2014 Marginal Response Curves (XGBoost)",
         subtitle = "Top 5 by Gain importance; others at median; shading = 95% CI (LOESS, span=0.35)",
         x="Predictor value", y="Predicted suitability",
         caption  = "Grey = raw XGBoost step-function. Smooth = LOESS. Dashed = 0.5 decision boundary.") +
    theme_chart_pub() +
    theme(strip.text=element_text(size=9,face="bold"),
          axis.text.x=element_text(size=8,angle=25,hjust=1))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig06_response_curves.png"),
         fig06, width=14, height=5.5, dpi=DPI, bg="white")
  cat("  ✓ Fig06_response_curves.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig06: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 8. FIGURE 7 — AUC COMPARISON
# FIX: geom_hline() for ref lines; as.numeric() on all cols
# ─────────────────────────────────────────────────────────────

cat("--- Figure 7: AUC Comparison ---\n")

tryCatch({
  eval_files <- c(MaxEnt="maxent",RF="rf",XGBoost="xgboost",
                  BRT="brt",GAM="gam",SVM="svm")
  adf <- do.call(rbind, lapply(names(eval_files), function(nm) {
    fp <- file.path(OUT_EVAL, paste0(eval_files[nm],"_evaluation.csv"))
    if (!file.exists(fp)) return(NULL)
    d  <- read.csv(fp, stringsAsFactors=FALSE)
    mc <- grep("cv_auc_mean|cv_auc$",names(d),value=TRUE,ignore.case=TRUE)[1]
    sc <- grep("cv_auc_sd|auc_sd",  names(d),value=TRUE,ignore.case=TRUE)[1]
    data.frame(
      algorithm  = nm,
      auc_mean   = as.numeric(d[[mc]][1]),
      auc_sd     = if(!is.na(sc %||% NA)) as.numeric(d[[sc]][1]) else NA_real_,
      is_ens     = FALSE,
      delong_sig = nm %in% c("XGBoost","BRT","GAM","SVM"),
      stringsAsFactors=FALSE)
  }))
  
  # Known results — ensure numeric types explicitly
  adf <- rbind(adf,
               data.frame(algorithm="Ensemble (AUC-wt)",
                          auc_mean=0.7239, auc_sd=NA_real_,
                          is_ens=TRUE, delong_sig=FALSE,
                          stringsAsFactors=FALSE))
  
  # Force numeric (prevents "discrete to continuous" error)
  adf$auc_mean <- as.numeric(adf$auc_mean)
  adf$auc_sd   <- as.numeric(adf$auc_sd)
  
  adf$ci_lo <- ifelse(!adf$is_ens,
                      pmax(adf$auc_mean - 1.96*adf$auc_sd, 0.5), NA_real_)
  adf$ci_hi <- ifelse(!adf$is_ens,
                      pmin(adf$auc_mean + 1.96*adf$auc_sd, 1.0), NA_real_)
  adf$algorithm <- factor(adf$algorithm, levels=adf$algorithm)
  adf$label_y   <- ifelse(is.na(adf$ci_hi),
                          adf$auc_mean + 0.025,
                          adf$ci_hi    + 0.025)
  
  # Colour map — must match factor levels exactly
  bar_cols <- setNames(
    c("#2A9D8F","#E63946","#F4A261","#264653",
      "#E9C46A","#8338EC","#06D6A0"),
    levels(adf$algorithm))
  
  # Significance data subsets
  sig_df  <- subset(adf,  delong_sig)
  ens_df  <- subset(adf,  is_ens)
  err_df  <- subset(adf, !is_ens)
  
  fig07 <- ggplot(adf, aes(x=algorithm, y=auc_mean)) +
    # FIX: geom_hline (NOT annotate segment) → no discrete/continuous conflict
    geom_hline(yintercept=0.75, linetype="dashed",
               color="#F4A261", linewidth=0.7) +
    geom_hline(yintercept=0.85, linetype="dotted",
               color="#E63946", linewidth=0.7) +
    # Reference labels at integer x positions
    annotate("text", x=1L, y=0.754,
             label="AUC = 0.75 (adequate)",
             hjust=0, size=2.9, color="#F4A261", fontface="italic") +
    annotate("text", x=1L, y=0.854,
             label="AUC = 0.85 (strong)",
             hjust=0, size=2.9, color="#E63946", fontface="italic") +
    # Error bars — only non-ensemble rows
    geom_errorbar(data=err_df,
                  aes(ymin=ci_lo, ymax=ci_hi),
                  width=0.22, linewidth=0.9, color="grey20") +
    # Points
    geom_point(aes(fill=algorithm), shape=21,
               size=5.5, stroke=0.5, color="white") +
    scale_fill_manual(values=bar_cols, guide="none") +
    # AUC value labels
    geom_text(aes(y=label_y, label=sprintf("%.4f",auc_mean)),
              size=3.1, fontface="bold", color="grey20") +
    # DeLong significance stars (* above CI bar)
    { if (nrow(sig_df)>0)
      geom_text(data=sig_df, aes(y=ci_hi+0.055), label="*",
                color="#E63946", size=7, fontface="bold")
      else geom_blank() } +
    # Ensemble star
    { if (nrow(ens_df)>0)
      geom_text(data=ens_df, aes(y=auc_mean+0.043),
                label="\u2605", color="#06D6A0", size=5)
      else geom_blank() } +
    scale_y_continuous(limits=c(0.50,1.00),
                       breaks=seq(0.50,1.00,0.05),
                       expand=expansion(mult=c(0.02,0.05))) +
    labs(title    = "Fig. 7 \u2014 Spatial Block CV AUC: All Models",
         subtitle = "\u00b11.96 SD (95% CI shown for individual algorithms); \u2605 = ensemble; * = DeLong p<0.05 vs ensemble",
         x=NULL, y="Spatial Block CV AUC (5-fold mean)",
         caption  = "DeLong et al. (1988); pROC package. Ensemble significantly outperforms XGBoost, BRT, GAM, SVM (p<0.05).") +
    theme_chart_pub() +
    theme(axis.text.x=element_text(size=10,color="grey15"))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig07_AUC_comparison.png"),
         fig07, width=9, height=7, dpi=DPI, bg="white")
  cat("  ✓ Fig07_AUC_comparison.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig07: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 9. FIGURE 8 — ENSEMBLE SUITABILITY + UNCERTAINTY + ZONES
# ─────────────────────────────────────────────────────────────

cat("--- Figure 8: Ensemble Suitability ---\n")

if (is.null(ens_agg)) {
  cat("  ✗ Ensemble raster not found\n\n")
} else tryCatch({
  ens_vals <- as.numeric(terra::values(ens_agg, na.rm=TRUE))
  e_hi <- as.numeric(quantile(ens_vals, 0.995, na.rm=TRUE))
  q25  <- as.numeric(quantile(ens_vals, 0.25,  na.rm=TRUE))
  q50  <- as.numeric(quantile(ens_vals, 0.50,  na.rm=TRUE))
  q75  <- as.numeric(quantile(ens_vals, 0.75,  na.rm=TRUE))
  cat(sprintf("  Suit range: %.4f\u2013%.4f | Q75=%.4f\n",
              min(ens_vals,na.rm=T), max(ens_vals,na.rm=T), q75))
  
  base_ens <- list(
    geom_sf(data=boundary_utm,fill=NA,color="white",linewidth=0.7),
    coord_sf(expand=FALSE),
    scale_x_continuous(labels=utm_fmt),
    scale_y_continuous(labels=utm_fmt),
    labs(x=NULL,y=NULL))
  
  pA <- ggplot() +
    tidyterra::geom_spatraster(data=ens_agg) +
    scale_fill_viridis_c(option="viridis", name="Suit.",
                         limits=c(0,e_hi), oob=scales::squish, na.value="white",
                         breaks=round(seq(0,e_hi,length.out=5),3)) +
    base_ens +
    geom_sf(data=sites_thin, color="#d62728",
            shape=19, size=1.5, alpha=0.85) +
    ggspatial::annotation_scale(location="bl") +
    ggspatial::annotation_north_arrow(location="tr",
                                      which_north="true",height=unit(0.75,"cm"),width=unit(0.60,"cm")) +
    labs(title="(A) Ensemble Suitability",
         subtitle="CV AUC=0.7239 | Boyce=0.9092 | KG=0.6063") +
    theme_map_pub() + theme(legend.position="right")
  
  if (!is.null(unc_agg)) {
    unc_vals <- as.numeric(terra::values(unc_agg, na.rm=TRUE))
    unc_q75  <- as.numeric(quantile(unc_vals, 0.75, na.rm=TRUE))
    unc_hi   <- as.numeric(quantile(unc_vals, 0.98, na.rm=TRUE))
    pB <- ggplot() +
      tidyterra::geom_spatraster(data=unc_agg) +
      scale_fill_viridis_c(option="plasma", name="SD",
                           limits=c(0,unc_hi), oob=scales::squish, na.value="white") +
      base_ens +
      geom_sf(data=sites_thin, color="white",
              shape=19, size=0.9, alpha=0.75) +
      labs(title="(B) Prediction Uncertainty",
           subtitle=sprintf("SD across 6 algorithms | Q75=%.3f",unc_q75)) +
      theme_map_pub() + theme(legend.position="right")
  } else pB <- ggplot()+annotate("text",x=0.5,y=0.5,label="Uncertainty\nnot found")+theme_void()
  
  # Classified suitability
  suit_rcl <- terra::classify(ens_agg,
                              matrix(c(-Inf,q25,1,q25,q50,2,q50,q75,3,q75,Inf,4),ncol=3,byrow=TRUE),
                              include.lowest=TRUE)
  names(suit_rcl) <- "class"
  sc_df <- as.data.frame(suit_rcl, xy=TRUE, na.rm=TRUE)
  sc_df$class_f <- factor(sc_df$class, levels=1:4,
                          labels=c("Low","Moderate","High","Very High"))
  
  pC <- ggplot() +
    geom_raster(data=sc_df, aes(x=x,y=y,fill=class_f)) +
    scale_fill_manual(name="Suitability",
                      values=c(Low="#264653",Moderate="#2A9D8F",
                               High="#E9C46A",`Very High`="#E63946"),
                      guide=guide_legend(reverse=TRUE), na.value="white") +
    geom_sf(data=boundary_utm,fill=NA,color="grey30",linewidth=0.6) +
    geom_sf(data=sites_thin, color="white",
            shape=19, size=1.0, alpha=0.80) +
    coord_sf(expand=FALSE) +
    scale_x_continuous(labels=utm_fmt) +
    scale_y_continuous(labels=utm_fmt) +
    labs(title="(C) Suitability Classification",
         subtitle="Quartile thresholds; 100% of sites in Very High class",
         x=NULL, y=NULL) +
    theme_map_pub() + theme(legend.position="right")
  
  fig08 <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title   = "Fig. 8 \u2014 AUC-Weighted Ensemble: Suitability, Uncertainty & Classification",
      caption = "All 6 algorithms on logistic probability scale [0\u20131]. Red/white dots = thinned pooled sites (N=190).",
      theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5),
                  plot.caption=element_text(size=7.5,hjust=0,color="grey45",face="italic"),
                  plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig08_ensemble_suitability.png"),
         fig08, width=16, height=7.5, dpi=DPI, bg="white")
  cat("  ✓ Fig08_ensemble_suitability.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig08: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 10. FIGURE 9 — SHAP BEESWARM
# FIX: search multiple paths; rebuild from model if RDS absent
# ─────────────────────────────────────────────────────────────

cat("--- Figure 9: SHAP Beeswarm ---\n")

tryCatch({
  # Step 1: Check if Script 19 already saved high-quality version
  s19_paths <- c(
    file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
    file.path(OUT_FIG_MAIN,"Fig09_shap_beeswarm.png"))
  s19_ok <- Find(function(p) file.exists(p) && file.info(p)$size>300000, s19_paths)
  if (!is.null(s19_ok)) {
    cat("  Using Script 19 high-quality version ✓\n\n")
  } else {
    # Step 2: Load or rebuild SHAP values
    shap_rds_paths <- c(
      file.path(OUT_SHAP,"xgboost_shap_values.rds"),
      file.path(OUT_SHAP,"shap_values_xgboost.rds"),
      file.path(OUT_SHAP,"shap_data.rds"),
      file.path(OUT_MOD_IND,"xgboost_shap_values.rds"))
    shap_rds <- Find(file.exists, shap_rds_paths)
    
    if (!is.null(shap_rds)) {
      shap_data <- readRDS(shap_rds)
      shap_mat  <- shap_data$shap_matrix
      X_mat     <- shap_data$predictor_matrix
    } else {
      # Rebuild inline from XGBoost model
      cat("  RDS not found — rebuilding SHAP from model...\n")
      xgb_path <- file.path(OUT_MOD_IND,"xgboost_model_final.bin")
      if (!file.exists(xgb_path)) stop("XGBoost model not found")
      xgb_mod_s <- xgboost::xgb.load(xgb_path)
      set.seed(42)
      samp_s <- terra::spatSample(pred_stack, size=2000,
                                  method="random", na.rm=TRUE, as.df=TRUE)
      for (col in c("Geology","Geomorphology")[c("Geology","Geomorphology") %in% names(samp_s)])
        samp_s[[col]] <- as.numeric(as.integer(samp_s[[col]]))
      X_mat    <- as.matrix(samp_s[,final_names,drop=FALSE])
      storage.mode(X_mat) <- "double"
      shap_raw <- predict(xgb_mod_s,
                          xgboost::xgb.DMatrix(X_mat), predcontrib=TRUE)
      # Remove BIAS column (last col)
      shap_mat <- shap_raw[, final_names, drop=FALSE]
      # Save for reuse
      saveRDS(list(shap_matrix=shap_mat, predictor_matrix=X_mat),
              file.path(OUT_SHAP,"xgboost_shap_values.rds"))
      cat("  SHAP rebuilt and saved ✓\n")
    }
    
    if (HAS_SHAP) {
      sv <- shapviz::shapviz(shap_mat, X=X_mat)
      p9 <- shapviz::sv_importance(sv, kind="beeswarm",
                                   max_display=length(final_names), alpha=0.55, size=1.2) +
        labs(title    = "Fig. 9 \u2014 SHAP Global Importance (XGBoost)",
             subtitle = "TreeSHAP; 2,000-cell random sample; colour = predictor value") +
        theme_chart_pub() +
        theme(legend.position="right")
    } else {
      # Fallback bar chart
      imp_df <- data.frame(
        predictor = final_names,
        mean_abs  = colMeans(abs(shap_mat), na.rm=TRUE))
      imp_df <- imp_df[order(-imp_df$mean_abs),]
      imp_df$predictor <- factor(imp_df$predictor,
                                 levels=rev(imp_df$predictor))
      p9 <- ggplot(imp_df, aes(x=predictor, y=mean_abs, fill=predictor)) +
        geom_col(width=0.7) +
        coord_flip() +
        scale_fill_manual(
          values=colorRampPalette(c("#264653","#2A9D8F","#E9C46A","#E63946"))(nrow(imp_df)),
          guide="none") +
        labs(title="Fig. 9 \u2014 SHAP Global Importance (XGBoost)",
             x=NULL, y="Mean |SHAP value|") +
        theme_chart_pub()
    }
    
    ggsave(file.path(OUT_FIG_MAIN,"Fig09_SHAP_beeswarm.png"),
           p9, width=9, height=7, dpi=DPI, bg="white")
    cat("  ✓ Fig09_SHAP_beeswarm.png\n\n")
  }
}, error=function(e) cat(sprintf("  ✗ Fig09: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 11. FIGURE 10 — SHAP DEPENDENCE PLOTS
# FIX: same SHAP rebuild path as Fig 9
# ─────────────────────────────────────────────────────────────

cat("--- Figure 10: SHAP Dependence ---\n")

tryCatch({
  s19_paths10 <- c(
    file.path(OUT_FIG_MAIN,"Fig10_SHAP_dependence.png"),
    file.path(OUT_FIG_MAIN,"Fig10_shap_dependence.png"))
  s19_ok10 <- Find(function(p) file.exists(p) && file.info(p)$size>400000, s19_paths10)
  
  if (!is.null(s19_ok10)) {
    cat("  Using Script 19 high-quality version ✓\n\n")
  } else {
    # Load SHAP (just rebuilt / already saved above)
    shap_rds <- file.path(OUT_SHAP,"xgboost_shap_values.rds")
    if (!file.exists(shap_rds)) stop("Run Fig09 block first")
    shap_data <- readRDS(shap_rds)
    shap_mat  <- shap_data$shap_matrix
    X_mat     <- shap_data$predictor_matrix
    
    if (HAS_SHAP) {
      sv <- shapviz::shapviz(shap_mat, X=X_mat)
      top5_s <- names(sort(colMeans(abs(shap_mat),na.rm=T),
                           decreasing=T))[1:5]
      
      dep_list <- lapply(top5_s, function(pred) {
        shapviz::sv_dependence(sv, v=pred, color_var="auto",
                               alpha=0.50, size=0.9) +
          labs(title=pred) +
          theme_chart_pub(base_size=9) +
          theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
                legend.key.size=unit(0.30,"cm"))
      })
    } else {
      # Fallback scatter plots
      top5_s <- names(sort(colMeans(abs(shap_mat),na.rm=T),
                           decreasing=T))[1:5]
      dep_list <- lapply(top5_s, function(pred) {
        if (!pred %in% colnames(shap_mat)) return(NULL)
        df <- data.frame(x=X_mat[,pred], y=shap_mat[,pred])
        df <- df[is.finite(df$x)&is.finite(df$y),]
        ggplot(df,aes(x,y)) +
          geom_point(size=0.5,alpha=0.4,color="#2A9D8F") +
          geom_smooth(method="loess",formula=y~x,
                      color="#E63946",linewidth=1.1,se=FALSE) +
          geom_hline(yintercept=0,linetype="dashed",color="grey50") +
          labs(title=pred,x=pred,y="SHAP value") +
          theme_chart_pub(base_size=9)
      })
    }
    
    fig10 <- patchwork::wrap_plots(Filter(Negate(is.null),dep_list),
                                   nrow=2L) +
      patchwork::plot_annotation(
        title   = "Fig. 10 \u2014 SHAP Dependence Plots \u2014 Top 5 Predictors (XGBoost)",
        caption = "TreeSHAP. Colour = interacting predictor with highest |SHAP| correlation (auto-selected).",
        theme=theme(
          plot.title=element_text(size=12,face="bold",hjust=0.5),
          plot.caption=element_text(size=7.5,hjust=0,color="grey45",face="italic"),
          plot.background=element_rect(fill="white",color=NA)))
    
    ggsave(file.path(OUT_FIG_MAIN,"Fig10_SHAP_dependence.png"),
           fig10, width=14, height=9, dpi=DPI, bg="white")
    cat("  ✓ Fig10_SHAP_dependence.png\n\n")
  }
}, error=function(e) cat(sprintf("  ✗ Fig10: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 12. FIGURE 11 — SHAP DOMINANT DRIVER MAP (UTM)
# FIX: terra::freq() WITHOUT na.rm=TRUE argument
# ─────────────────────────────────────────────────────────────

cat("--- Figure 11: SHAP Dominant Driver Map ---\n")

tryCatch({
  s19_11 <- file.path(OUT_FIG_MAIN,"Fig11_dominant_driver_map.png")
  if (file.exists(s19_11) && file.info(s19_11)$size > 300000) {
    cat("  Using Script 19 version ✓\n\n")
  } else if (!is.null(shap_r)) {
    # FIX: terra::freq() — NO na.rm argument
    freq_t  <- terra::freq(shap_r)
    freq_t  <- freq_t[!is.na(freq_t$value) & freq_t$value > 0, ]
    freq_t  <- freq_t[order(-freq_t$count), ]
    codes   <- as.integer(freq_t$value)
    labels  <- final_names[pmin(codes, length(final_names))]
    pcts    <- 100 * freq_t$count / sum(freq_t$count)
    drv_c   <- sapply(labels, function(nm)
      if (nm %in% names(DRIVER_PAL)) DRIVER_PAL[nm] else "#AAAAAA")
    
    shap_rcl <- terra::classify(shap_r,
                                cbind(codes, seq_along(codes)), others=NA)
    names(shap_rcl) <- "driver_idx"
    shap_df <- as.data.frame(shap_rcl, xy=TRUE, na.rm=FALSE)
    shap_df <- shap_df[!is.na(shap_df$driver_idx), ]
    shap_df$driver <- factor(shap_df$driver_idx,
                             levels=seq_along(labels), labels=labels)
    
    leg_labs <- sprintf("%s (%.1f%%)", labels, pcts)
    names(leg_labs) <- labels
    
    fig11 <- ggplot() +
      geom_raster(data=shap_df, aes(x=x,y=y,fill=driver)) +
      scale_fill_manual(
        name   = "Dominant predictor\n(% of 250m cells)",
        values = setNames(drv_c, labels),
        labels = leg_labs,
        na.value = "white",
        guide  = guide_legend(
          ncol = if(length(labels)>8) 2L else 1L,
          override.aes = list(size=3.5))) +
      geom_sf(data=boundary_utm, fill=NA,
              color="grey20", linewidth=0.75) +
      geom_sf(data=sites_thin, color="white",
              shape=19, size=0.9, alpha=0.75) +
      ggspatial::annotation_scale(location="bl") +
      ggspatial::annotation_north_arrow(location="tr",
                                        which_north="true", height=unit(0.8,"cm"),
                                        width=unit(0.65,"cm")) +
      coord_sf(expand=FALSE) +
      scale_x_continuous(labels=utm_fmt) +
      scale_y_continuous(labels=utm_fmt) +
      labs(title   = "Fig. 11 \u2014 SHAP Dominant Predictor per 250m Cell",
           subtitle = "TreeSHAP (XGBoost); 250m grid; \u03c7\u00b2=54,352, p<0.001",
           caption  = "White dots = thinned sites (N=190). Predictor with highest |SHAP value| per 250m cell.",
           x="Easting (m)", y="Northing (m)") +
      theme_map_pub() +
      theme(legend.position="right",
            legend.text=element_text(size=7.5),
            legend.title=element_text(size=8.5,face="bold"))
    
    ggsave(file.path(OUT_FIG_MAIN,"Fig11_shap_dominant_driver.png"),
           fig11, width=9, height=11.5, dpi=DPI, bg="white")
    cat("  ✓ Fig11_shap_dominant_driver.png\n\n")
  } else cat("  ✗ Fig11: SHAP raster not found\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig11: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 13. FIGURE 12 — DIACHRONIC SUB-MODEL DRIVER MAPS (UTM)
# FIX: terra::freq() WITHOUT na.rm=TRUE; handle UP CHELSA preds
# ─────────────────────────────────────────────────────────────

cat("--- Figure 12: Diachronic Driver Maps ---\n")

tryCatch({
  periods    <- c("LP","MP","UP")
  per_titles <- c(
    LP = "Lower Palaeolithic\n(N=70, Ens AUC=0.8153)",
    MP = "Middle Palaeolithic\n(N=110, Ens AUC=0.7278)",
    UP = "Upper Palaeolithic\n(N=74, Ens AUC=0.7845)")
  # Extended pool: main + potential CHELSA predictors
  pred_pool <- c(final_names, "BIO1","BIO12","BIO15")
  
  plist12 <- lapply(periods, function(per) {
    rp <- tryCatch(
      terra::rast(file.path(OUT_MOD_SUB, per,
                            sprintf("%s_dominant_driver_250m.tif",per))),
      error=function(e) NULL)
    if (is.null(rp)) return(
      ggplot()+annotate("text",x=0.5,y=0.5,
                        label=sprintf("%s raster\nnot found",per),
                        size=4,color="grey50")+
        theme_void()+theme(plot.background=element_rect(fill="white",color=NA)))
    
    # FIX: terra::freq() without na.rm
    freq_p  <- terra::freq(rp)
    freq_p  <- freq_p[!is.na(freq_p$value) & freq_p$value > 0, ]
    freq_p  <- freq_p[order(-freq_p$count), ]
    codes_p <- as.integer(freq_p$value)
    labs_p  <- pred_pool[pmin(codes_p, length(pred_pool))]
    pcts_p  <- 100 * freq_p$count / sum(freq_p$count)
    cols_p  <- sapply(labs_p, function(nm)
      if (nm %in% names(DRIVER_PAL)) DRIVER_PAL[nm] else "#AAAAAA")
    
    rp_rcl <- terra::classify(rp,
                              cbind(codes_p, seq_along(codes_p)), others=NA)
    names(rp_rcl) <- "di"
    df_p <- as.data.frame(rp_rcl, xy=TRUE, na.rm=FALSE)
    df_p <- df_p[!is.na(df_p$di), ]
    df_p$driver <- factor(df_p$di,
                          levels=seq_along(labs_p), labels=labs_p)
    
    leg_p <- sprintf("%s\n%.1f%%", labs_p, pcts_p)
    names(leg_p) <- labs_p
    
    ggplot() +
      geom_raster(data=df_p, aes(x=x,y=y,fill=driver)) +
      scale_fill_manual(
        name="Dominant\npredictor",
        values=setNames(cols_p, labs_p),
        labels=leg_p, na.value="white",
        guide=guide_legend(ncol=1,
                           override.aes=list(size=3.2))) +
      geom_sf(data=boundary_utm, fill=NA,
              color="grey20", linewidth=0.6) +
      coord_sf(expand=FALSE) +
      scale_x_continuous(labels=utm_fmt) +
      scale_y_continuous(labels=utm_fmt) +
      labs(title    = per_titles[per],
           subtitle = sprintf("Top driver: %s (%.1f%%)", labs_p[1], pcts_p[1]),
           x=NULL, y=NULL) +
      theme_map_pub(base_size=9) +
      theme(legend.position="right",
            legend.text=element_text(size=6.5),
            legend.title=element_text(size=7.5,face="bold"),
            legend.key.size=unit(0.30,"cm"),
            plot.title=element_text(size=9.5,face="bold",lineheight=1.2))
  })
  
  fig12 <- patchwork::wrap_plots(plist12, nrow=1L) +
    patchwork::plot_annotation(
      title   = "Fig. 12 \u2014 Diachronic SHAP Dominant Predictor Maps",
      subtitle = "XGBoost TreeSHAP sub-models; 250m grid; GAM excluded (REML timeout); 5-algorithm ensemble",
      caption  = "LP=Lower; MP=Middle; UP=Upper Palaeolithic. Colours consistent with Fig. 11.",
      theme=theme(
        plot.title    = element_text(size=12,face="bold",hjust=0.5),
        plot.subtitle = element_text(size=8.5,hjust=0.5,color="grey40"),
        plot.caption  = element_text(size=7.5,hjust=0,color="grey45",face="italic"),
        plot.background=element_rect(fill="white",color=NA)))
  
  ggsave(file.path(OUT_FIG_MAIN,"Fig12_diachronic_driver_maps.png"),
         fig12, width=16, height=9, dpi=DPI, bg="white")
  cat("  ✓ Fig12_diachronic_driver_maps.png\n\n")
}, error=function(e) cat(sprintf("  ✗ Fig12: %s\n\n",e$message)))

# ─────────────────────────────────────────────────────────────
# 14. FINAL CHECKLIST
# ─────────────────────────────────────────────────────────────

cat("--- Final Checklist ---\n\n")
items <- list(
  c("Fig01_study_area.png",            "Fig 1:  Study area (UTM 44N)"),
  c("Fig02_predictor_stack.png",       "Fig 2:  Predictor stack (13)"),
  c("Fig03_palaeochannel.png",         "Fig 3:  Palaeochannel"),
  c("Fig04_bias_correction.png",       "Fig 4:  Bias correction"),
  c("Fig05_spatial_cv_design.png",     "Fig 5:  Spatial CV design"),
  c("Fig06_response_curves.png",       "Fig 6:  Response curves (LOESS)"),
  c("Fig07_AUC_comparison.png",        "Fig 7:  AUC comparison"),
  c("Fig08_ensemble_suitability.png",  "Fig 8:  Ensemble suitability"),
  c("Fig09_SHAP_beeswarm.png",         "Fig 9:  SHAP beeswarm"),
  c("Fig10_SHAP_dependence.png",       "Fig 10: SHAP dependence"),
  c("Fig11_shap_dominant_driver.png",  "Fig 11: SHAP dominant driver"),
  c("Fig12_diachronic_driver_maps.png","Fig 12: Diachronic LP/MP/UP"))

n_ok <- 0L
for (it in items) {
  fp <- file.path(OUT_FIG_MAIN, it[1])
  ok <- file.exists(fp)
  if (ok) n_ok <- n_ok + 1L
  sz <- if(ok) sprintf("%.0f KB",file.info(fp)$size/1024) else "MISSING"
  cat(sprintf("  %s %-40s %s\n",
              if(ok)"\u2713" else "\u2717", it[2], sz))
}
cat(sprintf("\n  %d / %d figures present\n\n", n_ok, length(items)))
cat("========================================\n")
cat("SCRIPT 23 COMPLETE\n")
cat("========================================\n\n")
cat("All bugs fixed:\n")
cat("  [1] Fill conflict Fig01 → sites use COLOR ✓\n")
cat("  [2] India inset removed (UTM 44N only) ✓\n")
cat("  [3] All maps stay in UTM 44N ✓\n")
cat("  [4] Fig07 geom_hline() + as.numeric() ✓\n")
cat("  [5] terra::freq() without na.rm ✓\n")
cat("  [6] SHAP rebuilt from model if RDS absent ✓\n")
cat("========================================\n")