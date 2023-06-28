# Script: Demon on how to calculate summary statistics across multiple raster layers 

# -----------------------------------------------------------------------------
# Session setup
# -----------------------------------------------------------------------------

# Load required packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyverse,
  raster,
  terra,
  tidyterra
)


# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# -----------------------------------------------------------------------------
# Load and clean data 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# USA
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUTURE CLIMATE RCP4.5 (2050)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# A-PRIORI
# -----------------------------------------------------------------------------

# DEFAULT SETTINGS

# Import the desired raster layers as 'spatRast' objects 
usa.apriori.all.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.covselcombined.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.covselglm.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")

# TUNED SETTINGS

usa.apriori.all.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.covselcombined.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.covselglm.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.apriori.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")

# -----------------------------------------------------------------------------
# ALL 19
# -----------------------------------------------------------------------------

# DEFAULT

usa.all19.covsel.combined.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.covsel.glm.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.all19.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.R2VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")

# TUNED

usa.all19.covsel.combined.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.covsel.glm.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.all19.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.R2VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.all19.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")

# -----------------------------------------------------------------------------
# LITERATURE
# -----------------------------------------------------------------------------

# DEFAULT

usa.literature.aidoo.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.fordjour.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.narouikhandan.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.wang.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/USA_results/default_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")

# TUNED

usa.literature.aidoo.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.fordjour.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.narouikhandan.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")
usa.literature.wang.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/USA_results/tuned_maxent/Diaphorina_citri/RCP4.5/2050_RCP4.5_ensembled.grd")





# Get map to project our model over
wld_ext <- rnaturalearth::ne_countries(scale = "medium",
                                       returnclass = "sf") %>%
  dplyr::filter(continent == "North America")

# Crop each of the four rasters to extent of usa

usa.apriori.all.default <- terra::mask(usa.apriori.all.default, wld_ext)
usa.apriori.covselcombined.default <- terra::mask(usa.apriori.covselcombined.default, wld_ext)
usa.apriori.covselglm.default <- terra::mask(usa.apriori.covselglm.default, wld_ext)
usa.apriori.pcmvR2.default <- terra::mask(usa.apriori.pcmvR2.default, wld_ext)
usa.apriori.R2.default <- terra::mask(usa.apriori.R2.default, wld_ext)
usa.apriori.VIF.default <- terra::mask(usa.apriori.VIF.default, wld_ext)

usa.apriori.all.tuned <- terra::mask(usa.apriori.all.tuned, wld_ext)
usa.apriori.covselcombined.tuned <- terra::mask(usa.apriori.covselcombined.tuned, wld_ext)
usa.apriori.covselglm.tuned <- terra::mask(usa.apriori.covselglm.tuned, wld_ext)
usa.apriori.pcmvR2.tuned <- terra::mask(usa.apriori.pcmvR2.tuned, wld_ext)
usa.apriori.R2.tuned <- terra::mask(usa.apriori.R2.tuned, wld_ext)
usa.apriori.VIF.tuned <- terra::mask(usa.apriori.VIF.tuned, wld_ext)

usa.all19.covsel.combined.default <- terra::mask(usa.all19.covsel.combined.default, wld_ext)
usa.all19.covsel.glm.default <- terra::mask(usa.all19.covsel.glm.default, wld_ext)
usa.all19.all19.default <- terra::mask(usa.all19.all19.default, wld_ext)     
usa.all19.pcmvR2.default <- terra::mask(usa.all19.pcmvR2.default, wld_ext)        
usa.all19.R2.default = terra::mask(usa.all19.R2.default, wld_ext)             
usa.all19.R2VIF.default = terra::mask(usa.all19.R2VIF.default, wld_ext)            
usa.all19.VIF.default = terra::mask(usa.all19.VIF.default, wld_ext)            

usa.all19.covsel.combined.tuned <- terra::mask(usa.all19.covsel.combined.tuned, wld_ext)           
usa.all19.covsel.glm.tuned <- terra::mask(usa.all19.covsel.glm.tuned, wld_ext)
usa.all19.all19.tuned <-   terra::mask(usa.all19.all19.tuned, wld_ext)   
usa.all19.pcmvR2.tuned <- terra::mask(usa.all19.pcmvR2.tuned, wld_ext)        
usa.all19.R2.tuned =      terra::mask(usa.all19.R2.tuned, wld_ext)        
usa.all19.R2VIF.tuned =   terra::mask(usa.all19.R2VIF.tuned, wld_ext)          
usa.all19.VIF.tuned =     terra::mask(usa.all19.VIF.tuned, wld_ext)        

usa.literature.aidoo.default <- terra::mask(usa.literature.aidoo.default, wld_ext)
usa.literature.fordjour.default <- terra::mask(usa.literature.fordjour.default, wld_ext)
usa.literature.narouikhandan.default <- terra::mask(usa.literature.narouikhandan.default, wld_ext)
usa.literature.wang.default <- terra::mask(usa.literature.wang.default, wld_ext)

usa.literature.aidoo.tuned <- terra::mask(usa.literature.aidoo.tuned, wld_ext)
usa.literature.fordjour.tuned <- terra::mask(usa.literature.fordjour.tuned, wld_ext)
usa.literature.narouikhandan.tuned <- terra::mask(usa.literature.narouikhandan.tuned, wld_ext)
usa.literature.wang.tuned <- terra::mask(usa.literature.wang.tuned, wld_ext)




# Set the CRS projection
terra::crs(usa.apriori.all.default) <- "epsg:4326"
terra::crs(usa.apriori.covselcombined.default) <- "epsg:4326"
terra::crs(usa.apriori.covselglm.default) <- "epsg:4326"
terra::crs(usa.apriori.pcmvR2.default) <- "epsg:4326"
terra::crs(usa.apriori.R2.default) <- "epsg:4326"
terra::crs(usa.apriori.VIF.default) <- "epsg:4326"
terra::crs(usa.apriori.all.tuned) <- "epsg:4326"
terra::crs(usa.apriori.covselcombined.tuned) <- "epsg:4326"
terra::crs(usa.apriori.covselglm.tuned) <- "epsg:4326"
terra::crs(usa.apriori.pcmvR2.tuned) <- "epsg:4326"
terra::crs(usa.apriori.R2.tuned) <- "epsg:4326"
terra::crs(usa.apriori.VIF.tuned) <- "epsg:4326"
terra::crs(usa.all19.covsel.combined.default) <- "epsg:4326"
terra::crs(usa.all19.covsel.glm.default) <- "epsg:4326"
terra::crs(usa.all19.all19.default) <- "epsg:4326"
terra::crs(usa.all19.pcmvR2.default) <- "epsg:4326"
terra::crs(usa.all19.R2.default) <- "epsg:4326"
terra::crs(usa.all19.R2VIF.default) <- "epsg:4326"
terra::crs(usa.all19.VIF.default) <- "epsg:4326"
terra::crs(usa.all19.covsel.combined.tuned) <- "epsg:4326"
terra::crs(usa.all19.covsel.glm.tuned) <- "epsg:4326"
terra::crs(usa.all19.all19.tuned) <- "epsg:4326"
terra::crs(usa.all19.pcmvR2.tuned) <- "epsg:4326"
terra::crs(usa.all19.R2.tuned) <- "epsg:4326"
terra::crs(usa.all19.R2VIF.tuned) <- "epsg:4326"
terra::crs(usa.all19.VIF.tuned) <- "epsg:4326"
terra::crs(usa.literature.aidoo.default) <- "epsg:4326"
terra::crs(usa.literature.fordjour.default) <- "epsg:4326"
terra::crs(usa.literature.narouikhandan.default) <- "epsg:4326"
terra::crs(usa.literature.wang.default) <- "epsg:4326"
terra::crs(usa.literature.aidoo.tuned) <- "epsg:4326"
terra::crs(usa.literature.fordjour.tuned) <- "epsg:4326"
terra::crs(usa.literature.narouikhandan.tuned) <- "epsg:4326"
terra::crs(usa.literature.wang.tuned) <- "epsg:4326"


# Plot some maps 
terra::plot(usa.apriori.all.default)
terra::plot(usa.apriori.covselcombined.default)
terra::plot(usa.apriori.covselglm.default)
terra::plot(usa.apriori.pcmvR2.default)

# -----------------------------------------------------------------------------
# Perform raster operations 
# -----------------------------------------------------------------------------

# To perform operations across multiple rasters, they need to be stored
# in the same variable 
r_stack <- c(
  usa.apriori.all.default, 
  usa.apriori.covselcombined.default, 
  usa.apriori.covselglm.default, 
  usa.apriori.pcmvR2.default, 
  usa.apriori.R2.default, 
  usa.apriori.VIF.default, 
  usa.apriori.all.tuned,
  usa.apriori.covselcombined.tuned, 
  usa.apriori.covselglm.tuned, 
  usa.apriori.pcmvR2.tuned, 
  usa.apriori.R2.tuned, 
  usa.apriori.VIF.tuned, 
  usa.all19.covsel.combined.default, 
  usa.all19.covsel.glm.default, 
  usa.all19.all19.default,      
  usa.all19.pcmvR2.default,         
  usa.all19.R2.default,              
  usa.all19.R2VIF.default,           
  usa.all19.VIF.default,      
  usa.all19.covsel.combined.tuned,    
  usa.all19.covsel.glm.tuned, 
  usa.all19.all19.tuned,    
  usa.all19.pcmvR2.tuned,       
  usa.all19.R2.tuned,      
  usa.all19.R2VIF.tuned,           
  usa.all19.VIF.tuned,        
  usa.literature.aidoo.default, 
  usa.literature.fordjour.default, 
  usa.literature.narouikhandan.default, 
  usa.literature.wang.default, 
  usa.literature.aidoo.tuned, 
  usa.literature.fordjour.tuned, 
  usa.literature.narouikhandan.tuned, 
  usa.literature.wang.tuned
)

class(r_stack) # Note, the new object is still a spatial raster object 

# Calculate mean prediction for each raster cell
r_mean <- terra::app(r_stack, fun = mean)

# Plot the new layer
terra::plot(r_mean)

# Calculate standard deviation prediction for each raster cell
r_sd <- terra::app(r_stack, fun = sd)

# Plot the new layer
terra::plot(r_sd)

# Plot mean map  
plot_mean <- ggplot() +
  # Plot world boundary
  geom_sf(data = wld_ext, fill = NA) +
  # Plot MaxEnt prediction raster
  geom_spatraster(
    data = r_mean,
    maxcell = 5e+6         # maxcell = Inf
  ) +
  geom_sf(data = wld_ext, fill = NA) +
  # geom_point(data = sp_data_inv, 
  #           aes(x = lon, y = lat)) + 
  # Control raster colour and legend
  scale_fill_whitebox_c(
    palette = "muted",
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1)
  ) +
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Suitability",
    subtitle = "(a) Mean"
  ) +
  coord_sf(
    xlim = c(-127, -65),
    ylim = c(10, 50),
    crs = 4326,
    expand = FALSE
  ) +
  theme(
    legend.position = "right"
  ) +
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
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(text_size = 10)
  ) 
plot_mean

# Plot sd map  
plot_sd <- ggplot() +
  # Plot world boundary
  geom_sf(data = wld_ext, fill = NA) +
  # Plot MaxEnt prediction raster
  geom_spatraster(
    data = r_sd,
    maxcell = 5e+6         # maxcell = Inf
  ) +
  geom_sf(data = wld_ext, fill = NA) +
  # geom_point(data = sp_data_inv, 
  #           aes(x = lon, y = lat)) + 
  # Control raster colour and legend
  scale_fill_whitebox_c(
    palette = "muted",
    breaks = seq(0, 0.3, 0.05),
    limits = c(0, 0.3)
  ) +
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Suitability",
    subtitle = "(b) Standard deviation"
  ) +
  coord_sf(
    xlim = c(-127, -65),
    ylim = c(10, 50),
    crs = 4326,
    expand = FALSE
  ) +
  theme(
    legend.position = "right"
  ) +
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
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering(text_size = 10)
  ) 
plot_sd

# Combine plots 
plots <- cowplot::plot_grid(
  plot_mean,
  plot_sd,
  nrow = 1
)
plots

# Combine plots 
ggsave(
  here::here("./figures/usa_RCP4.5_2050.png"),
  dpi = 600,
  height = 10,
  width = 10
)

ggsave(
  here::here("./figures/usa_RCP4.5_2050.svg"),
  height = 10,
  width = 10
)

