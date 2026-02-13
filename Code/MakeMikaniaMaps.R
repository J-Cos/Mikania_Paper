# =============================================================================
# 03_MikaniaMaps.R
# Produces publication-ready maps of Mikania classification results
# for each band of Data/Mikania_Results_10m_3Bands.tif
#
# Bands:
#   1 = mikania_prob_mean       (continuous, 0–1)
#   2 = mikania_uncertainty_std (continuous, 0–0.2)
#   3 = mikania_binary_optimized (binary, 0/1)
# =============================================================================

library(tidyverse)
library(terra)
library(tidyterra)
library(cowplot)
library(rnaturalearth)

# ---- Data paths (relative to project root) ----------------------------------
raster_path  <- "Data/Mikania_Results_10m_3Bands.tif"
cpc_path     <- "Data/CPC.shp"
chitwan_path <- "Data/WDPA_WDOECM_Apr2024_Public_10905_shp-polygons.shp"
parsa_path   <- "Data/WDPA_WDOECM_Apr2024_Public_10089_shp-polygons.shp"

# ---- Load data once ---------------------------------------------------------
mikania_rast <- rast(raster_path)

# PA boundaries (Chitwan NP + Parsa WR), reproject to match raster (UTM 45N)
cpc_crs  <- vect(cpc_path)                                    # just for CRS
chitwan  <- project(vect(chitwan_path), cpc_crs)
parsa    <- project(vect(parsa_path),  cpc_crs)
mgmt_areas <- rbind(chitwan, parsa)

# Crop/mask raster to PA boundaries only (exclude buffer zone)
mikania_rast <- mask(crop(mikania_rast, mgmt_areas), mgmt_areas)

# Nepal country boundary for the locator inset
nepal <- ne_countries(scale = "medium", country = "Nepal", returnclass = "sf") |>
  vect() |>
  project(cpc_crs)

# ---- Define 3 zoom inset extents -------------------------------------------
# A: Western Chitwan — shifted south for maximum coverage
# B: Central Chitwan — core area
# C: Parsa — zoomed in tighter on the main cluster
zoom_list <- list(
  A = ext(200000, 225000, 3040000, 3055000),
  B = ext(245000, 270000, 3038000, 3053000),
  C = ext(280000, 292000, 3017400, 3024600)
)

# ---- Mapping function -------------------------------------------------------
make_mikania_map <- function(band = 1,
                             out_dir = "Figures",
                             width = 12, height = 8, dpi = 300) {

  # Select the band
  r <- mikania_rast[[band]]
  is_binary <- (band == 3)

  # For binary: make categorical with both levels
  if (is_binary) {
    levels(r) <- data.frame(id = c(0, 1), label = c("Absent", "Present"))
  }

  # For continuous bands: get global range for consistent colour scale
  if (!is_binary) {
    val_range <- range(values(r, na.rm = TRUE))
  }

  # --- Helper: build fill scale (ensures consistency) -----------------------
  add_fill <- function(p, show_legend = TRUE) {
    if (is_binary) {
      p + scale_fill_manual(
        values = c("Absent" = "#440154", "Present" = "#FDE725"),
        na.value = NA, na.translate = FALSE, name = NULL,
        guide = if (show_legend) "legend" else "none"
      )
    } else {
      p + scale_fill_viridis_c(
        option = "viridis", na.value = NA,
        limits = val_range,
        name = if (band == 1) "Probability" else "Uncertainty (SD)",
        guide = if (show_legend) "colourbar" else "none"
      )
    }
  }

  # --- Build the main map ----------------------------------------------------
  p_main <- ggplot() +
    geom_spatraster(data = r)

  p_main <- add_fill(p_main, show_legend = TRUE)

  # PA boundaries only (bold), no buffer zone
  pa_bb <- ext(mgmt_areas)
  p_main <- p_main +
    geom_spatvector(data = mgmt_areas, fill = NA, colour = "black",
                    linewidth = 0.8)

  # Add labelled zoom boxes
  for (i in seq_along(zoom_list)) {
    ze  <- zoom_list[[i]]
    lbl <- names(zoom_list)[i]
    p_main <- p_main +
      annotate("rect",
               xmin = ze[1], xmax = ze[2],
               ymin = ze[3], ymax = ze[4],
               fill = NA, colour = "black", linewidth = 0.4) +
      annotate("text",
               x = ze[1] + (ze[2] - ze[1]) * 0.5,
               y = ze[4] + (ze[4] - ze[3]) * 0.08,
               label = lbl, size = 3.5, fontface = "bold")
  }

  p_main <- p_main +
    theme_void() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.90, 0.80),
      legend.background = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    coord_sf(
      xlim = c(pa_bb[1] - 3000, pa_bb[2] + 3000),
      ylim = c(pa_bb[3] - 3000, pa_bb[4] + 3000),
      expand = FALSE
    )

  # --- Build the 3 zoom insets -----------------------------------------------
  zoom_plots <- list()
  for (i in seq_along(zoom_list)) {
    ze  <- zoom_list[[i]]
    lbl <- names(zoom_list)[i]

    r_zoom    <- crop(r, ze)
    mgmt_zoom <- crop(mgmt_areas, ze)

    p_z <- ggplot() +
      geom_spatraster(data = r_zoom)

    p_z <- add_fill(p_z, show_legend = FALSE)

    p_z <- p_z +
      geom_spatvector(data = mgmt_zoom, fill = NA, colour = "black",
                      linewidth = 0.8) +
      # Panel label
      annotate("text", x = ze[1] + (ze[2] - ze[1]) * 0.06,
               y = ze[4] - (ze[4] - ze[3]) * 0.08,
               label = lbl, size = 4, fontface = "bold", hjust = 0) +
      theme_void() +
      theme(
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black",
                                    linewidth = 0.5),
        plot.margin = margin(2, 2, 2, 2)
      ) +
      coord_sf(
        xlim = c(ze[1], ze[2]),
        ylim = c(ze[3], ze[4]),
        expand = FALSE
      )

    zoom_plots[[i]] <- p_z
  }

  # --- Build the Nepal locator inset -----------------------------------------
  # PA outlines filled solid black, no buffer zone
  p_nepal <- ggplot() +
    geom_spatvector(data = nepal, fill = "grey80", colour = "grey30",
                    linewidth = 0.3) +
    geom_spatvector(data = mgmt_areas, fill = "black", colour = "black",
                    linewidth = 0.4) +
    theme_void() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black",
                                  linewidth = 0.5),
      plot.margin = margin(2, 2, 2, 2)
    )

  # --- Assemble with cowplot --------------------------------------------------
  # Layout: main map upper portion, 3 zoom insets along bottom,
  #         Nepal inset positioned above zoom inset A
  inset_w <- 0.30
  inset_h <- 0.26
  gap     <- 0.02

  p_assembled <- ggdraw() +
    draw_plot(p_main,  x = 0, y = 0.22, width = 1, height = 0.78) +
    # Nepal inset — just above inset A
    draw_plot(p_nepal, x = gap, y = inset_h + gap,
              width = 0.22, height = 0.18) +
    # Zoom insets A, B, C along the bottom
    draw_plot(zoom_plots[[1]], x = gap,                      y = 0.0,
              width = inset_w, height = inset_h) +
    draw_plot(zoom_plots[[2]], x = gap + inset_w + gap,      y = 0.0,
              width = inset_w, height = inset_h) +
    draw_plot(zoom_plots[[3]], x = gap + 2*(inset_w + gap),  y = 0.0,
              width = inset_w, height = inset_h)

  # --- Save -------------------------------------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, paste0("MikaniaMap_Band", band, ".png"))
  ggsave(out_file, p_assembled, width = width, height = height,
         dpi = dpi, bg = "white")
  message("Saved: ", out_file)

  invisible(p_assembled)
}

# ---- Generate all three maps ------------------------------------------------
for (b in 1:3) {
  make_mikania_map(band = b)
}

message("Done — all maps saved to Figures/")
