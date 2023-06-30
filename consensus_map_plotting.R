africa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                          returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")

mean_afr_2010 = terra::rast("consensus_maps/africa/2010/mean/africa_2010.grd")
crs(mean_afr_2010) = "EPSG:4326"
stdev_afr_2010 = terra::rast("consensus_maps/africa/2010/stdev/africa_2010_sd.grd")
crs(stdev_afr_2010) = "EPSG:4326"

mean_afr_2050_RCP4.5 = terra::rast("consensus_maps/africa/2050/RCP4.5/mean/africa_2050_RCP4.5.grd")
crs(mean_afr_2050_RCP4.5) = "EPSG:4326"
stdev_afr_2050_RCP4.5 = terra::rast("consensus_maps/africa/2050/RCP4.5/stdev/africa_2050_RCP4.5_stdev.grd")
crs(stdev_afr_2050_RCP4.5) = "EPSG:4326"

mean_afr_2070_RCP4.5 = terra::rast("consensus_maps/africa/2070/RCP4.5/mean/africa_2070_RCP4.5.grd")
crs(mean_afr_2070_RCP4.5) = "EPSG:4326"
stdev_afr_2070_RCP4.5 = terra::rast("consensus_maps/africa/2070/RCP4.5/stdev/africa_2070_RCP4.5_stdev.grd")
crs(stdev_afr_2070_RCP4.5) = "EPSG:4326"

mean_afr_2050_RCP8.5 = terra::rast("consensus_maps/africa/2050/RCP8.5/mean/africa_2050_RCP8.5.grd")
crs(mean_afr_2050_RCP8.5) = "EPSG:4326"
stdev_afr_2050_RCP8.5 = terra::rast("consensus_maps/africa/2050/RCP8.5/stdev/africa_2050_RCP8.5_stdev.grd")
crs(stdev_afr_2050_RCP8.5) = "EPSG:4326"

mean_afr_2070_RCP8.5 = terra::rast("consensus_maps/africa/2070/RCP8.5/mean/africa_2070_RCP8.5.grd")
crs(mean_afr_2070_RCP8.5) = "EPSG:4326"
stdev_afr_2070_RCP8.5 = terra::rast("consensus_maps/africa/2070/RCP8.5/stdev/africa_2070_RCP8.5_stdev.grd")
crs(stdev_afr_2070_RCP8.5) = "EPSG:4326"

# Calculate the difference rasters
diff_2010_2050 <- mean_afr_2050_RCP8.5 - mean_afr_2010
diff_2050_2070 <- mean_afr_2070_RCP8.5 - mean_afr_2050_RCP8.5
diff_2010_2070 = mean_afr_2070_RCP8.5 - mean_afr_2010

combined_diff = diff_2010_2050 + diff_2050_2070

# using the tidyterra package to plot the map
mean_afr_2010_2070_plot = ggplot() +
  
  tidyterra::geom_spatraster(data = diff_2010_2070) +
  scale_fill_whitebox_c(
    palette = "muted",
    #breaks = seq(0, 1, 0.2),
    limits = c(-0.4, 0.4)
  ) +
  
  # Add country borders
  geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
  
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "P(suitability)",
    subtitle = "Mean difference 2010 - 2070"
  ) +
  # Create title for the legend
  theme(legend.position = "right") +
  # Add scale bar to bottom-right of map
  ggspatial::annotation_scale(
    location = "bl",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Change appearance of the legend
  guides(
    fill = guide_colorbar(ticks = FALSE)
  ) +
  coord_sf(ylim = c(-35, 40))

ggsave(plot = mean_afr_2010_2070_plot, filename = "consensus_maps/africa/mean_afr_2010_2070.png", dpi=400,height=6,width = 6)

##############################################################
# SOUTH AFRICA ONLY
##############################################################

##############
# 2010 - 2070
##############

# Get map to project our model over
sa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                      returnclass = "sf") %>%
  # dplyr::filter(name == c("South Africa", "Lesotho", "Swaziland"))
  dplyr::filter(name %in% c("South Africa", "Lesotho", "Swaziland"))

sa_ext = st_set_crs(sa_ext, 4326)

# Load provincial borders data
provincial_borders <- ne_states(country = "South Africa", returnclass = "sf")

# Clip the raster to the extent of South Africa
diff_2010_2070_sa <- raster::mask(diff_2010_2070, sa_ext)


# Plot mean map  
mean_sa_diff_2010_2070 <- ggplot() +
  # Plot world boundary
  geom_sf(data = sa_ext, fill = NA) +
  # Plot MaxEnt prediction raster
  geom_spatraster(
    data = diff_2010_2070_sa,
    maxcell = 5e+6         # maxcell = Inf
  ) +
  geom_sf(data = sa_ext, fill = NA) +
  # Plot provincial borders
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  # Control raster colour and legend
  scale_fill_whitebox_c(
    palette = "muted",
    #breaks = seq(0, 1, 0.2),
    limits = c(-0.4, 0.4)
  ) +
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Suitability",
    subtitle = "Mean difference 2010 - 2070"
  ) +
  # Crops map to just the geographic extent of SA
  coord_sf(
    xlim = c(15, 34),
    ylim = c(-21, -36),
    crs = 4326,
    expand = FALSE
  ) +
  theme(
    legend.position = "right"
  ) +
  # Add scale bar to bottom-right of map
  ggspatial::annotation_scale(
    location = "br",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) 

mean_sa_diff_2010_2070
ggsave(plot = mean_sa_diff_2010_2070, filename = "consensus_maps/africa/mean_southafr_2010_2070.png", dpi=400,height=6,width = 6)



##############
# 2010 - 2050
##############

# using the tidyterra package to plot the map
mean_afr_2010_2050_plot = ggplot() +
  
  tidyterra::geom_spatraster(data = diff_2010_2050) +
  scale_fill_whitebox_c(
    palette = "muted",
    #breaks = seq(0, 1, 0.2),
    limits = c(-0.4, 0.4)
  ) +
  
  # Add country borders
  geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
  
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "P(suitability)",
    subtitle = "Mean difference 2010 - 2050"
  ) +
  # Create title for the legend
  theme(legend.position = "right") +
  # Add scale bar to bottom-right of map
  ggspatial::annotation_scale(
    location = "bl",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Change appearance of the legend
  guides(
    fill = guide_colorbar(ticks = FALSE)
  ) +
  coord_sf(ylim = c(-35, 40))

ggsave(plot = mean_afr_2010_2050_plot, filename = "consensus_maps/africa/mean_afr_2010_2050.png", dpi=400,height=6,width = 6)


# SOUTH AFRICA ONLY

# Get map to project our model over
sa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                      returnclass = "sf") %>%
  # dplyr::filter(name == c("South Africa", "Lesotho", "Swaziland"))
  dplyr::filter(name %in% c("South Africa", "Lesotho", "Swaziland"))

sa_ext = st_set_crs(sa_ext, 4326)

# Load provincial borders data
provincial_borders <- ne_states(country = "South Africa", returnclass = "sf")

# Clip the raster to the extent of South Africa
diff_2010_2050_sa <- raster::mask(diff_2010_2050, sa_ext)


# Plot mean map  
mean_sa_diff_2010_2050 <- ggplot() +
  # Plot world boundary
  geom_sf(data = sa_ext, fill = NA) +
  # Plot MaxEnt prediction raster
  geom_spatraster(
    data = diff_2010_2050_sa,
    maxcell = 5e+6         # maxcell = Inf
  ) +
  geom_sf(data = sa_ext, fill = NA) +
  # Plot provincial borders
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  # Control raster colour and legend
  scale_fill_whitebox_c(
    palette = "muted",
    #breaks = seq(0, 1, 0.2),
    limits = c(-0.4, 0.4)
  ) +
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Suitability",
    subtitle = "Mean difference 2010 - 2050"
  ) +
  # Crops map to just the geographic extent of SA
  coord_sf(
    xlim = c(15, 34),
    ylim = c(-21, -36),
    crs = 4326,
    expand = FALSE
  ) +
  theme(
    legend.position = "right"
  ) +
  # Add scale bar to bottom-right of map
  ggspatial::annotation_scale(
    location = "br",          # 'bl' = bottom left
    style = "ticks",
    width_hint = 0.2
  ) +
  # Add north arrow
  ggspatial::annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.175, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) 

mean_sa_diff_2010_2050

ggsave(plot = mean_sa_diff_2010_2050, filename = "consensus_maps/africa/mean_southafr_2010_2050.png", dpi=400,height=6,width = 6)

