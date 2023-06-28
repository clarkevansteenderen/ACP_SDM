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
# BRAZIL
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CURRENT CLIAMTE (2010)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# A-PRIORI
# -----------------------------------------------------------------------------

# DEFAULT SETTINGS

# Import the desired raster layers as 'spatRast' objects 
brazil.apriori.all.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.covselcombined.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.covselglm.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")

# TUNED SETTINGS

brazil.apriori.all.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.covselcombined.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.covselglm.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.apriori.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")

# -----------------------------------------------------------------------------
# ALL 19
# -----------------------------------------------------------------------------

# DEFAULT

brazil.all19.covsel.combined.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.covsel.glm.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.all19.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.R2VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")

# TUNED

brazil.all19.covsel.combined.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.covsel.glm.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.all19.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.R2VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.all19.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")

# -----------------------------------------------------------------------------
# LITERATURE
# -----------------------------------------------------------------------------

# DEFAULT

brazil.literature.aidoo.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.fordjour.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.narouikhandan.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.wang.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/Brazil_results/default_maxent/Diaphorina_citri/2010_ensembled.grd")

# TUNED

brazil.literature.aidoo.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.fordjour.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.narouikhandan.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")
brazil.literature.wang.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/Brazil_results/tuned_maxent/Diaphorina_citri/2010_ensembled.grd")





# Get map to project our model over
wld_ext <- rnaturalearth::ne_countries(scale = "medium",
                                       returnclass = "sf") %>%
  dplyr::filter(name == "Brazil")

# Crop each of the four rasters to extent of brazil

brazil.apriori.all.default <- terra::mask(brazil.apriori.all.default, wld_ext)
brazil.apriori.covselcombined.default <- terra::mask(brazil.apriori.covselcombined.default, wld_ext)
brazil.apriori.covselglm.default <- terra::mask(brazil.apriori.covselglm.default, wld_ext)
brazil.apriori.pcmvR2.default <- terra::mask(brazil.apriori.pcmvR2.default, wld_ext)
brazil.apriori.R2.default <- terra::mask(brazil.apriori.R2.default, wld_ext)
brazil.apriori.VIF.default <- terra::mask(brazil.apriori.VIF.default, wld_ext)

brazil.apriori.all.tuned <- terra::mask(brazil.apriori.all.tuned, wld_ext)
brazil.apriori.covselcombined.tuned <- terra::mask(brazil.apriori.covselcombined.tuned, wld_ext)
brazil.apriori.covselglm.tuned <- terra::mask(brazil.apriori.covselglm.tuned, wld_ext)
brazil.apriori.pcmvR2.tuned <- terra::mask(brazil.apriori.pcmvR2.tuned, wld_ext)
brazil.apriori.R2.tuned <- terra::mask(brazil.apriori.R2.tuned, wld_ext)
brazil.apriori.VIF.tuned <- terra::mask(brazil.apriori.VIF.tuned, wld_ext)

brazil.all19.covsel.combined.default <- terra::mask(brazil.all19.covsel.combined.default, wld_ext)
brazil.all19.covsel.glm.default <- terra::mask(brazil.all19.covsel.glm.default, wld_ext)
brazil.all19.all19.default <- terra::mask(brazil.all19.all19.default, wld_ext)     
brazil.all19.pcmvR2.default <- terra::mask(brazil.all19.pcmvR2.default, wld_ext)        
brazil.all19.R2.default = terra::mask(brazil.all19.R2.default, wld_ext)             
brazil.all19.R2VIF.default = terra::mask(brazil.all19.R2VIF.default, wld_ext)            
brazil.all19.VIF.default = terra::mask(brazil.all19.VIF.default, wld_ext)            

brazil.all19.covsel.combined.tuned <- terra::mask(brazil.all19.covsel.combined.tuned, wld_ext)           
brazil.all19.covsel.glm.tuned <- terra::mask(brazil.all19.covsel.glm.tuned, wld_ext)
brazil.all19.all19.tuned <-   terra::mask(brazil.all19.all19.tuned, wld_ext)   
brazil.all19.pcmvR2.tuned <- terra::mask(brazil.all19.pcmvR2.tuned, wld_ext)        
brazil.all19.R2.tuned =      terra::mask(brazil.all19.R2.tuned, wld_ext)        
brazil.all19.R2VIF.tuned =   terra::mask(brazil.all19.R2VIF.tuned, wld_ext)          
brazil.all19.VIF.tuned =     terra::mask(brazil.all19.VIF.tuned, wld_ext)        

brazil.literature.aidoo.default <- terra::mask(brazil.literature.aidoo.default, wld_ext)
brazil.literature.fordjour.default <- terra::mask(brazil.literature.fordjour.default, wld_ext)
brazil.literature.narouikhandan.default <- terra::mask(brazil.literature.narouikhandan.default, wld_ext)
brazil.literature.wang.default <- terra::mask(brazil.literature.wang.default, wld_ext)

brazil.literature.aidoo.tuned <- terra::mask(brazil.literature.aidoo.tuned, wld_ext)
brazil.literature.fordjour.tuned <- terra::mask(brazil.literature.fordjour.tuned, wld_ext)
brazil.literature.narouikhandan.tuned <- terra::mask(brazil.literature.narouikhandan.tuned, wld_ext)
brazil.literature.wang.tuned <- terra::mask(brazil.literature.wang.tuned, wld_ext)




# Set the CRS projection
terra::crs(brazil.apriori.all.default) <- "epsg:4326"
terra::crs(brazil.apriori.covselcombined.default) <- "epsg:4326"
terra::crs(brazil.apriori.covselglm.default) <- "epsg:4326"
terra::crs(brazil.apriori.pcmvR2.default) <- "epsg:4326"
terra::crs(brazil.apriori.R2.default) <- "epsg:4326"
terra::crs(brazil.apriori.VIF.default) <- "epsg:4326"
terra::crs(brazil.apriori.all.tuned) <- "epsg:4326"
terra::crs(brazil.apriori.covselcombined.tuned) <- "epsg:4326"
terra::crs(brazil.apriori.covselglm.tuned) <- "epsg:4326"
terra::crs(brazil.apriori.pcmvR2.tuned) <- "epsg:4326"
terra::crs(brazil.apriori.R2.tuned) <- "epsg:4326"
terra::crs(brazil.apriori.VIF.tuned) <- "epsg:4326"
terra::crs(brazil.all19.covsel.combined.default) <- "epsg:4326"
terra::crs(brazil.all19.covsel.glm.default) <- "epsg:4326"
terra::crs(brazil.all19.all19.default) <- "epsg:4326"
terra::crs(brazil.all19.pcmvR2.default) <- "epsg:4326"
terra::crs(brazil.all19.R2.default) <- "epsg:4326"
terra::crs(brazil.all19.R2VIF.default) <- "epsg:4326"
terra::crs(brazil.all19.VIF.default) <- "epsg:4326"
terra::crs(brazil.all19.covsel.combined.tuned) <- "epsg:4326"
terra::crs(brazil.all19.covsel.glm.tuned) <- "epsg:4326"
terra::crs(brazil.all19.all19.tuned) <- "epsg:4326"
terra::crs(brazil.all19.pcmvR2.tuned) <- "epsg:4326"
terra::crs(brazil.all19.R2.tuned) <- "epsg:4326"
terra::crs(brazil.all19.R2VIF.tuned) <- "epsg:4326"
terra::crs(brazil.all19.VIF.tuned) <- "epsg:4326"
terra::crs(brazil.literature.aidoo.default) <- "epsg:4326"
terra::crs(brazil.literature.fordjour.default) <- "epsg:4326"
terra::crs(brazil.literature.narouikhandan.default) <- "epsg:4326"
terra::crs(brazil.literature.wang.default) <- "epsg:4326"
terra::crs(brazil.literature.aidoo.tuned) <- "epsg:4326"
terra::crs(brazil.literature.fordjour.tuned) <- "epsg:4326"
terra::crs(brazil.literature.narouikhandan.tuned) <- "epsg:4326"
terra::crs(brazil.literature.wang.tuned) <- "epsg:4326"


# Plot some maps 
terra::plot(brazil.apriori.all.default)
terra::plot(brazil.apriori.covselcombined.default)
terra::plot(brazil.apriori.covselglm.default)
terra::plot(brazil.apriori.pcmvR2.default)

# -----------------------------------------------------------------------------
# Perform raster operations 
# -----------------------------------------------------------------------------

# To perform operations across multiple rasters, they need to be stored
# in the same variable 
r_stack <- c(
  brazil.apriori.all.default, 
  brazil.apriori.covselcombined.default, 
  brazil.apriori.covselglm.default, 
  brazil.apriori.pcmvR2.default, 
  brazil.apriori.R2.default, 
  brazil.apriori.VIF.default, 
  brazil.apriori.all.tuned,
  brazil.apriori.covselcombined.tuned, 
  brazil.apriori.covselglm.tuned, 
  brazil.apriori.pcmvR2.tuned, 
  brazil.apriori.R2.tuned, 
  brazil.apriori.VIF.tuned, 
  brazil.all19.covsel.combined.default, 
  brazil.all19.covsel.glm.default, 
  brazil.all19.all19.default,      
  brazil.all19.pcmvR2.default,         
  brazil.all19.R2.default,              
  brazil.all19.R2VIF.default,           
  brazil.all19.VIF.default,      
  brazil.all19.covsel.combined.tuned,    
  brazil.all19.covsel.glm.tuned, 
  brazil.all19.all19.tuned,    
  brazil.all19.pcmvR2.tuned,       
  brazil.all19.R2.tuned,      
  brazil.all19.R2VIF.tuned,           
  brazil.all19.VIF.tuned,        
  brazil.literature.aidoo.default, 
  brazil.literature.fordjour.default, 
  brazil.literature.narouikhandan.default, 
  brazil.literature.wang.default, 
  brazil.literature.aidoo.tuned, 
  brazil.literature.fordjour.tuned, 
  brazil.literature.narouikhandan.tuned, 
  brazil.literature.wang.tuned
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
  # Crops map to just the geographic extent of SA
  #coord_sf(
  # xlim = c(15, 33.5),
  #  ylim = c(-22, -35),
  #  crs = 4326,
  #  expand = FALSE
  # ) + 
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
  # Crops map to just the geographic extent of brazil
  # coord_sf(
  # xlim = c(15, 33.5),
  # ylim = c(-22, -35),
  # crs = 4326,
  # expand = FALSE
  #) + 
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
  here::here("./figures/brazil_2010.png"),
  dpi = 600,
  height = 10,
  width = 10
)

