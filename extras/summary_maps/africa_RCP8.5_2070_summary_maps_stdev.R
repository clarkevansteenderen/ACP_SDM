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
# AFRICA
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# FUTURE CLIAMTE RCP 8.5 (2070)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# A-PRIORI
# -----------------------------------------------------------------------------

# DEFAULT SETTINGS

# Import the desired raster layers as 'spatRast' objects 
africa.apriori.all.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.covselcombined.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.covselglm.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")

# TUNED SETTINGS

africa.apriori.all.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_all_predictors/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.covselcombined.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_combined/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.covselglm.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_covsel_glm/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_PCMV_R2/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_R2/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.apriori.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/A_PRIORI/a_priori_reduced_VIF/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")

# -----------------------------------------------------------------------------
# ALL 19
# -----------------------------------------------------------------------------

# DEFAULT

africa.all19.covsel.combined.default <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.covsel.glm.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.all19.default <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.pcmvR2.default <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.R2.default =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.R2VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.VIF.default =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")

# TUNED

africa.all19.covsel.combined.tuned <-            terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_combined/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.covsel.glm.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_covsel_glm/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.all19.tuned <-      terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_predictors/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.pcmvR2.tuned <-         terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_PCMV_R2/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.R2.tuned =              terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.R2VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_R2_VIF/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.all19.VIF.tuned =             terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/ALL_19/all_19_reduced_VIF/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")

# -----------------------------------------------------------------------------
# LITERATURE
# -----------------------------------------------------------------------------

# DEFAULT

africa.literature.aidoo.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.fordjour.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.narouikhandan.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.wang.default <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/Africa_results/default_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")

# TUNED

africa.literature.aidoo.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/aidoo_set/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.fordjour.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/fordjour_set/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.narouikhandan.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/naroui_khandan_set/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")
africa.literature.wang.tuned <- terra::rast("D:/POST_DOC/Updated_psyllid_SDM/FINAL_RESULTS/LITERATURE/wang_set/Africa_results/tuned_maxent/Diaphorina_citri/RCP8.5/2070_RCP8.5_ensembled.grd")





# Get map to project our model over
wld_ext <- rnaturalearth::ne_countries(scale = "medium",
                                       returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")

# Crop each of the four rasters to extent of Africa

africa.apriori.all.default <- terra::mask(africa.apriori.all.default, wld_ext)
africa.apriori.covselcombined.default <- terra::mask(africa.apriori.covselcombined.default, wld_ext)
africa.apriori.covselglm.default <- terra::mask(africa.apriori.covselglm.default, wld_ext)
africa.apriori.pcmvR2.default <- terra::mask(africa.apriori.pcmvR2.default, wld_ext)
africa.apriori.R2.default <- terra::mask(africa.apriori.R2.default, wld_ext)
africa.apriori.VIF.default <- terra::mask(africa.apriori.VIF.default, wld_ext)

africa.apriori.all.tuned <- terra::mask(africa.apriori.all.tuned, wld_ext)
africa.apriori.covselcombined.tuned <- terra::mask(africa.apriori.covselcombined.tuned, wld_ext)
africa.apriori.covselglm.tuned <- terra::mask(africa.apriori.covselglm.tuned, wld_ext)
africa.apriori.pcmvR2.tuned <- terra::mask(africa.apriori.pcmvR2.tuned, wld_ext)
africa.apriori.R2.tuned <- terra::mask(africa.apriori.R2.tuned, wld_ext)
africa.apriori.VIF.tuned <- terra::mask(africa.apriori.VIF.tuned, wld_ext)

africa.all19.covsel.combined.default <- terra::mask(africa.all19.covsel.combined.default, wld_ext)
africa.all19.covsel.glm.default <- terra::mask(africa.all19.covsel.glm.default, wld_ext)
africa.all19.all19.default <- terra::mask(africa.all19.all19.default, wld_ext)     
africa.all19.pcmvR2.default <- terra::mask(africa.all19.pcmvR2.default, wld_ext)        
africa.all19.R2.default = terra::mask(africa.all19.R2.default, wld_ext)             
africa.all19.R2VIF.default = terra::mask(africa.all19.R2VIF.default, wld_ext)            
africa.all19.VIF.default = terra::mask(africa.all19.VIF.default, wld_ext)            

africa.all19.covsel.combined.tuned <- terra::mask(africa.all19.covsel.combined.tuned, wld_ext)           
africa.all19.covsel.glm.tuned <- terra::mask(africa.all19.covsel.glm.tuned, wld_ext)
africa.all19.all19.tuned <-   terra::mask(africa.all19.all19.tuned, wld_ext)   
africa.all19.pcmvR2.tuned <- terra::mask(africa.all19.pcmvR2.tuned, wld_ext)        
africa.all19.R2.tuned =      terra::mask(africa.all19.R2.tuned, wld_ext)        
africa.all19.R2VIF.tuned =   terra::mask(africa.all19.R2VIF.tuned, wld_ext)          
africa.all19.VIF.tuned =     terra::mask(africa.all19.VIF.tuned, wld_ext)        

africa.literature.aidoo.default <- terra::mask(africa.literature.aidoo.default, wld_ext)
africa.literature.fordjour.default <- terra::mask(africa.literature.fordjour.default, wld_ext)
africa.literature.narouikhandan.default <- terra::mask(africa.literature.narouikhandan.default, wld_ext)
africa.literature.wang.default <- terra::mask(africa.literature.wang.default, wld_ext)

africa.literature.aidoo.tuned <- terra::mask(africa.literature.aidoo.tuned, wld_ext)
africa.literature.fordjour.tuned <- terra::mask(africa.literature.fordjour.tuned, wld_ext)
africa.literature.narouikhandan.tuned <- terra::mask(africa.literature.narouikhandan.tuned, wld_ext)
africa.literature.wang.tuned <- terra::mask(africa.literature.wang.tuned, wld_ext)




# Set the CRS projection
terra::crs(africa.apriori.all.default) <- "epsg:4326"
terra::crs(africa.apriori.covselcombined.default) <- "epsg:4326"
terra::crs(africa.apriori.covselglm.default) <- "epsg:4326"
terra::crs(africa.apriori.pcmvR2.default) <- "epsg:4326"
terra::crs(africa.apriori.R2.default) <- "epsg:4326"
terra::crs(africa.apriori.VIF.default) <- "epsg:4326"
terra::crs(africa.apriori.all.tuned) <- "epsg:4326"
terra::crs(africa.apriori.covselcombined.tuned) <- "epsg:4326"
terra::crs(africa.apriori.covselglm.tuned) <- "epsg:4326"
terra::crs(africa.apriori.pcmvR2.tuned) <- "epsg:4326"
terra::crs(africa.apriori.R2.tuned) <- "epsg:4326"
terra::crs(africa.apriori.VIF.tuned) <- "epsg:4326"
terra::crs(africa.all19.covsel.combined.default) <- "epsg:4326"
terra::crs(africa.all19.covsel.glm.default) <- "epsg:4326"
terra::crs(africa.all19.all19.default) <- "epsg:4326"
terra::crs(africa.all19.pcmvR2.default) <- "epsg:4326"
terra::crs(africa.all19.R2.default) <- "epsg:4326"
terra::crs(africa.all19.R2VIF.default) <- "epsg:4326"
terra::crs(africa.all19.VIF.default) <- "epsg:4326"
terra::crs(africa.all19.covsel.combined.tuned) <- "epsg:4326"
terra::crs(africa.all19.covsel.glm.tuned) <- "epsg:4326"
terra::crs(africa.all19.all19.tuned) <- "epsg:4326"
terra::crs(africa.all19.pcmvR2.tuned) <- "epsg:4326"
terra::crs(africa.all19.R2.tuned) <- "epsg:4326"
terra::crs(africa.all19.R2VIF.tuned) <- "epsg:4326"
terra::crs(africa.all19.VIF.tuned) <- "epsg:4326"
terra::crs(africa.literature.aidoo.default) <- "epsg:4326"
terra::crs(africa.literature.fordjour.default) <- "epsg:4326"
terra::crs(africa.literature.narouikhandan.default) <- "epsg:4326"
terra::crs(africa.literature.wang.default) <- "epsg:4326"
terra::crs(africa.literature.aidoo.tuned) <- "epsg:4326"
terra::crs(africa.literature.fordjour.tuned) <- "epsg:4326"
terra::crs(africa.literature.narouikhandan.tuned) <- "epsg:4326"
terra::crs(africa.literature.wang.tuned) <- "epsg:4326"


# Plot some maps 
terra::plot(africa.apriori.all.default)
terra::plot(africa.apriori.covselcombined.default)
terra::plot(africa.apriori.covselglm.default)
terra::plot(africa.apriori.pcmvR2.default)

# -----------------------------------------------------------------------------
# Perform raster operations 
# -----------------------------------------------------------------------------

# To perform operations across multiple rasters, they need to be stored
# in the same variable 
r_stack <- c(
  africa.apriori.all.default, 
  africa.apriori.covselcombined.default, 
  africa.apriori.covselglm.default, 
  africa.apriori.pcmvR2.default, 
  africa.apriori.R2.default, 
  africa.apriori.VIF.default, 
  africa.apriori.all.tuned,
  africa.apriori.covselcombined.tuned, 
  africa.apriori.covselglm.tuned, 
  africa.apriori.pcmvR2.tuned, 
  africa.apriori.R2.tuned, 
  africa.apriori.VIF.tuned, 
  africa.all19.covsel.combined.default, 
  africa.all19.covsel.glm.default, 
  africa.all19.all19.default,      
  africa.all19.pcmvR2.default,         
  africa.all19.R2.default,              
  africa.all19.R2VIF.default,           
  africa.all19.VIF.default,      
  africa.all19.covsel.combined.tuned,    
  africa.all19.covsel.glm.tuned, 
  africa.all19.all19.tuned,    
  africa.all19.pcmvR2.tuned,       
  africa.all19.R2.tuned,      
  africa.all19.R2VIF.tuned,           
  africa.all19.VIF.tuned,        
  africa.literature.aidoo.default, 
  africa.literature.fordjour.default, 
  africa.literature.narouikhandan.default, 
  africa.literature.wang.default, 
  africa.literature.aidoo.tuned, 
  africa.literature.fordjour.tuned, 
  africa.literature.narouikhandan.tuned, 
  africa.literature.wang.tuned
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
  coord_sf(
    xlim = c(-20, 52),
    ylim = c(-36, 40),
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
  #Crops map to just the geographic extent of Africa
  coord_sf(
    xlim = c(-20, 52),
    ylim = c(-36, 40),
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
  here::here("./figures/africa_RCP8.5_2070.png"),
  dpi = 600,
  height = 10,
  width = 10
)

