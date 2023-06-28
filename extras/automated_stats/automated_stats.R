####################################################################################
# AUTOMATED MODEL VALIDATIONS
####################################################################################

source("predictor_sets.R")

predictors <- terra::rast( list.files(
  here::here("./data/environmental_layers/current/wc2.1_2.5m/") ,
  full.names = TRUE,
  pattern = '.tif'
))

# Import and crop KG shapefile 
kg_layer <- rgdal::readOGR(here::here("./data/shapefiles/koppen_geiger"), 
                           "WC05_1975H_Koppen", 
                           verbose = FALSE)


# Reproject KG-layer
geo_proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kg_layer <- sp::spTransform(kg_layer, geo_proj)


##########################################
# Set up the files
##########################################

categories = c("A_PRIORI", "ALL_19", "LITERATURE")

a_priori_result_files = c("a_priori_all_predictors", "a_priori_reduced_PCMV_R2", "a_priori_reduced_R2",
                          "a_priori_reduced_VIF", "a_priori_covsel_combined", "a_priori_covsel_glm")

all_19_result_files = c("all_19_predictors",
                        "all_19_reduced_PCMV_R2", "all_19_reduced_R2",
                        "all_19_reduced_R2_VIF", "all_19_reduced_VIF", "all_19_covsel_combined",
                        "all_19_covsel_glm")

literature_result_files = c("aidoo_set", "fordjour_set", "naroui_khandan_set", "wang_set")

locations = c("Africa_results", "Brazil_results", "USA_results")

model = c("default_maxent", "tuned_maxent")

##########################################
# A PRIORI TARGET FILES
##########################################

a_priori_filepaths = list()

  for(m in 1:length(a_priori_result_files)){
        for(q in 1:length(locations)){
          for(r in 1:length(model)){
  
        target_files = paste("FINAL_RESULTS/", 
                    categories[1],
                    "/",
                    a_priori_result_files[m],
                    "/",
                    locations[q],
                    "/",
                    model[r],
                    "/",
                    "Diaphorina_citri/2010_ensembled.grd",
                    sep ="")
        #print(target_files)
        a_priori_filepaths[[length(a_priori_filepaths) + 1]] = target_files # add the filepath to the list

    } # for(r in 1:length(model))
  } # for(q in 1:length(locations))
}# for(m in 1:length(a_priori_result_files))

length(a_priori_filepaths)

a_priori_filepaths_africa = a_priori_filepaths[grep("Africa", a_priori_filepaths)]  
a_priori_filepaths_brazil = a_priori_filepaths[grep("Brazil", a_priori_filepaths)] 
a_priori_filepaths_usa = a_priori_filepaths[grep("USA", a_priori_filepaths)] 

##########################################
# ALL 19 TARGET FILES
##########################################

all_19_filepaths = list()

for(m in 1:length(all_19_result_files)){
  for(q in 1:length(locations)){
    for(r in 1:length(model)){
      
      target_files = paste("FINAL_RESULTS/", 
                           categories[2],
                           "/",
                           all_19_result_files[m],
                           "/",
                           locations[q],
                           "/",
                           model[r],
                           "/",
                           "Diaphorina_citri/2010_ensembled.grd",
                           sep ="")
      #print(target_files)
      all_19_filepaths[[length(all_19_filepaths) + 1]] = target_files # add the filepath to the list
      
    } # for(r in 1:length(model))
  } # for(q in 1:length(locations))
}# for(m in 1:length(all_19_result_files))

length(all_19_filepaths)

all_19_filepaths_africa = all_19_filepaths[grep("Africa", all_19_filepaths)]  
all_19_filepaths_brazil = all_19_filepaths[grep("Brazil", all_19_filepaths)] 
all_19_filepaths_usa = all_19_filepaths[grep("USA", all_19_filepaths)] 

##########################################
# LITERATURE TARGET FILES
##########################################

literature_filepaths = list()

for(m in 1:length(literature_result_files)){
  for(q in 1:length(locations)){
    for(r in 1:length(model)){
      
      target_files = paste("FINAL_RESULTS/", 
                           categories[3],
                           "/",
                           literature_result_files[m],
                           "/",
                           locations[q],
                           "/",
                           model[r],
                           "/",
                           "Diaphorina_citri/2010_ensembled.grd",
                           sep ="")
      #print(target_files)
      literature_filepaths[[length(literature_filepaths) + 1]] = target_files # add the filepath to the list
      
    } # for(r in 1:length(model))
  } # for(q in 1:length(locations))
}# for(m in 1:length(literature_result_files))

length(literature_filepaths)

literature_filepaths_africa = literature_filepaths[grep("Africa", literature_filepaths)]  
literature_filepaths_brazil = literature_filepaths[grep("Brazil", literature_filepaths)] 
literature_filepaths_usa = literature_filepaths[grep("USA", literature_filepaths)] 


##############################################################################################################################
##########################################
# AFRICA VALIDATIONS
##########################################
##############################################################################################################################

a_priori_filepaths_africa
all_19_filepaths_africa
literature_filepaths_africa

##############################################################################################################################
##########################################
# A PRIORI SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_a_priori_africa <- data.frame(matrix(nrow=length(a_priori_filepaths_africa), ncol=10))
colnames(results_df_a_priori_africa) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                          "lr", "df", "pr", "file", "predictor_set", "model")

##########################################
# For loop that runs through each result
# file, and extracts the relevant
# statistical values and saves the map
# for current suitability, and the probability curve
##########################################

for(k in 1:length(a_priori_filepaths_africa)){

  #print(a_priori_filepaths_africa[[k]])
  file_name_clean <- stringr::str_remove(stringr::str_remove(a_priori_filepaths_africa[[k]], 
                                           "FINAL_RESULTS/.*?/"), 
                                "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_a_priori_africa$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_a_priori_africa$model[k] = model_type
  results_df_a_priori_africa$locality[k] = "africa"
  
  # check which set of predictors you need, based on the current file at hand
  if(length( a_priori_filepaths_africa[k][grep("a_priori_all", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_all
    results_df_a_priori_africa$predictor_set[k] = "a_priori_all"
  } else if(length( a_priori_filepaths_africa[k][grep("a_priori_reduced_PCMV_R2", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_pcmv_r2
    results_df_a_priori_africa$predictor_set[k] = "a_priori_reduced_PCMV_R2"
  } else if(length( a_priori_filepaths_africa[k][grep("a_priori_reduced_R2", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_r2
    results_df_a_priori_africa$predictor_set[k] = "a_priori_reduced_R2"
  } else if(length( a_priori_filepaths_africa[k][grep("a_priori_reduced_VIF", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_VIF
    results_df_a_priori_africa$predictor_set[k] = "a_priori_reduced_VIF"
  } else if(length( a_priori_filepaths_africa[k][grep("a_priori_covsel_combined", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined_apriori
    results_df_a_priori_africa$predictor_set[k] = "a_priori_covsel_combined"
  } else if(length( a_priori_filepaths_africa[k][grep("a_priori_covsel_glm", a_priori_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm_apriori
    results_df_a_priori_africa$predictor_set[k] = "a_priori_covsel_glm"
  } else {
    # handle the case where none of the conditions are met
    predictor_subset = NULL
  }
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  africa_pred <- terra::rast(
    here::here(
      a_priori_filepaths_africa[k]
    )
  )
  
  # Set the CRS
  crs(africa_pred) <- "EPSG:4326"
  
  # Get map of Africa to project our model over
  africa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "Africa")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = africa_pred) +
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    coord_sf(
      xlim = c(-20, 54),
      ylim = c(-40, 38),
      crs = 4326,
      expand = FALSE
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
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", a_priori_filepaths_africa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in Africa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -48 & lat < 40 & lon > -26 & lon < 55
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(africa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = africa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(africa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = africa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_a_priori_africa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_a_priori_africa$lr[k] = c$`LR Chisq`
  results_df_a_priori_africa$df[k] = c$Df
  results_df_a_priori_africa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
 print(paste("OPTIMAL VALUE = ", optimal)) 
 
 results_df_a_priori_africa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_a_priori_africa$sensitivity[k] = sensitivity
  
#######################################################
# PLOT MAP, THRESHOLDED BY OPTIMAL VALUE
#######################################################
  
africa_binary_optimal = africa_pred >= optimal
  
values_at_points <- terra::extract(x = africa_binary_optimal, # thresholded map
                    y = pts_species, 
                    xy = TRUE,
                    na.rm = TRUE) %>%
                    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
values_at_points = na.omit(values_at_points) 

# Plot using this optimal threshold, and add the species GPS points as an overlay
thresholded_map = ggplot() +
  tidyterra::geom_spatraster(data = africa_pred, aes(fill = median)) +
  
  geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
  
  scale_fill_whitebox_c(
    palette = "muted",
    breaks = seq(0, 1, 0.2),
    limits = c(optimal, 1) # here we set the optimal threshold
  ) +
  # Control axis and legend labels 
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "P(suitability)"
  ) +
  
  coord_sf(
    xlim = c(-20, 54),
    ylim = c(-40, 38),
    crs = 4326,
    expand = FALSE
  ) +
  
  # Create title for the legend
  theme(legend.position = "right" ) +
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
  geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
  scale_color_manual(values = c("red", "black")) +
  labs(color = "Present")

ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
       plot = thresholded_map, width = 8, height = 10)
  
}#for(k in a_priori_filepaths_africa){



##############################################################################################################################
##########################################
# ALL 19 PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_all_19_africa <- data.frame(matrix(nrow=length(all_19_filepaths_africa), ncol=10))
colnames(results_df_all_19_africa) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                        "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(all_19_filepaths_africa)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(all_19_filepaths_africa[[k]], 
                                           "FINAL_RESULTS/.*?/"), 
                                "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_all_19_africa$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_all_19_africa$model[k] = model_type
  results_df_all_19_africa$locality[k] = "africa"
  
  # check which set of predictors you need, based on the current file at hand
  # all 19 files
  if(length( all_19_filepaths_africa[k][grep("all_19_predictors", all_19_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$all_19_all
    results_df_all_19_africa$predictor_set[k] = "all_19"
  }else if(length( all_19_filepaths_africa[k][grep("all_19_reduced_PCMV_R2", all_19_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_pcmv_r2 
    results_df_all_19_africa$predictor_set[k] = "all_19_reduced_PCMV_R2"
  }else if(length( all_19_filepaths_africa[k][grepl("all_19_reduced_R2_VIF", all_19_filepaths_africa[k], fixed = TRUE)]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2_VIF
    results_df_all_19_africa$predictor_set[k] = "all_19_reduced_R2_VIF"
  }else if(length( all_19_filepaths_africa[k][grepl("all_19_reduced_R2", all_19_filepaths_africa[k], fixed =TRUE)]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2
    results_df_all_19_africa$predictor_set[k] = "all_19_reduced_R2"
  }else if(length( all_19_filepaths_africa[k][grep("all_19_reduced_VIF", all_19_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_VIF
    results_df_all_19_africa$predictor_set[k] = "all_19_reduced_VIF"
  }else if(length( all_19_filepaths_africa[k][grep("all_19_covsel_combined", all_19_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined
    results_df_all_19_africa$predictor_set[k] = "all_19_covsel_combined"
  }else if(length( all_19_filepaths_africa[k][grep("all_19_covsel_glm", all_19_filepaths_africa[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm
    results_df_all_19_africa$predictor_set[k] = "all_19_covsel_glm"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  africa_pred <- terra::rast(
    here::here(
      all_19_filepaths_africa[k]
    )
  )
  
  # Set the CRS
  crs(africa_pred) <- "EPSG:4326"
  
  # Get map of Africa to project our model over
  africa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "Africa")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = africa_pred) +
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    coord_sf(
      xlim = c(-20, 54),
      ylim = c(-40, 38),
      crs = 4326,
      expand = FALSE
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
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", all_19_filepaths_africa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in Africa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -48 & lat < 40 & lon > -26 & lon < 55
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(africa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = africa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(africa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = africa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_all_19_africa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_all_19_africa$lr[k] = c$`LR Chisq`
  results_df_all_19_africa$df[k] = c$Df
  results_df_all_19_africa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_all_19_africa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_all_19_africa$sensitivity[k] = sensitivity
  
  #######################################################
  # PLOT MAP, THRESHOLDED BY OPTIMAL VALUE
  #######################################################
  
  africa_binary_optimal = africa_pred >= optimal
  
  values_at_points <- terra::extract(x = africa_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points) 
  
  # Plot using this optimal threshold, and add the species GPS points as an overlay
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = africa_pred, aes(fill = median)) +
    
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
    
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-20, 54),
      ylim = c(-40, 38),
      crs = 4326,
      expand = FALSE
    ) +
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in all_19_filepaths_africa){



##############################################################################################################################
##########################################
# LITERATURE PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_literature_africa <- data.frame(matrix(nrow=length(literature_filepaths_africa), ncol=10))
colnames(results_df_literature_africa) <- c("locality","optimal", "sensitivity", "suit_score_coef", 
                                            "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(literature_filepaths_africa)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(literature_filepaths_africa[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_literature_africa$file[k] = file_name_clean
  model_type = sub(".*/", "", file_name_clean)
  results_df_literature_africa$model[k] = model_type
  results_df_literature_africa$locality[k] = "africa"
  
  # check which set of predictors you need, based on the current file at hand
  
  if(length( literature_filepaths_africa[k][grep("wang", literature_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$wang
    results_df_literature_africa$predictor_set[k] = "wang"
  }else if(length( literature_filepaths_africa[k][grep("aidoo", literature_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$aidoo
    results_df_literature_africa$predictor_set[k] = "aidoo"
  }else if(length( literature_filepaths_africa[k][grep("fordjour", literature_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$fordjour
    results_df_literature_africa$predictor_set[k] = "fordjour"
  }else if(length( literature_filepaths_africa[k][grep("naroui", literature_filepaths_africa[k])]) ){
    predictor_subset = predictor_sets$naroui_khandan
    results_df_literature_africa$predictor_set[k] = "naroui_khandan"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  africa_pred <- terra::rast(
    here::here(
      literature_filepaths_africa[k]
    )
  )
  
  # Set the CRS
  crs(africa_pred) <- "EPSG:4326"
  
  # Get map of Africa to project our model over
  africa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "Africa")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = africa_pred) +
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    coord_sf(
      xlim = c(-20, 54),
      ylim = c(-40, 38),
      crs = 4326,
      expand = FALSE
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
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", literature_filepaths_africa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in Africa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -48 & lat < 40 & lon > -26 & lon < 55
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(africa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = africa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(africa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = africa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_literature_africa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_literature_africa$lr[k] = c$`LR Chisq`
  results_df_literature_africa$df[k] = c$Df
  results_df_literature_africa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_literature_africa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_literature_africa$sensitivity[k] = sensitivity
  
  #######################################################
  # PLOT MAP, THRESHOLDED BY OPTIMAL VALUE
  #######################################################
  
  africa_binary_optimal = africa_pred >= optimal
  
  values_at_points <- terra::extract(x = africa_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points) 
  
  # Plot using this optimal threshold, and add the species GPS points as an overlay
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = africa_pred, aes(fill = median)) +
    
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.2) +
    
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-20, 54),
      ylim = c(-40, 38),
      crs = 4326,
      expand = FALSE
    ) +
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in literature_filepaths_africa){


results_df_a_priori_africa
results_df_all_19_africa
results_df_literature_africa

africa_results = rbind(results_df_a_priori_africa,
      results_df_all_19_africa,
      results_df_literature_africa)


##############################################################################################################################
##########################################
# BRAZIL VALIDATIONS
##########################################
##############################################################################################################################

a_priori_filepaths_brazil
all_19_filepaths_brazil
literature_filepaths_brazil

##############################################################################################################################
##########################################
# A PRIORI SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_a_priori_brazil <- data.frame(matrix(nrow=length(a_priori_filepaths_brazil), ncol=10))
colnames(results_df_a_priori_brazil) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                          "lr", "df", "pr", "file", "predictor_set", "model")

##########################################
# For loop that runs through each result
# file, and extracts the relevant
# statistical values and saves the map
# for current suitability, and the probability curve
##########################################

for(k in 1:length(a_priori_filepaths_brazil)){
  
  #print(a_priori_filepaths_brazil[[k]])
  file_name_clean <- stringr::str_remove(stringr::str_remove(a_priori_filepaths_brazil[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_a_priori_brazil$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_a_priori_brazil$model[k] = model_type
  results_df_a_priori_brazil$locality[k] = "brazil"
  
  # check which set of predictors you need, based on the current file at hand
  if(length( a_priori_filepaths_brazil[k][grep("a_priori_all", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$a_priori_all
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_all"
  } else if(length( a_priori_filepaths_brazil[k][grep("a_priori_reduced_PCMV_R2", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_pcmv_r2
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_reduced_PCMV_R2"
  } else if(length( a_priori_filepaths_brazil[k][grep("a_priori_reduced_R2", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_r2
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_reduced_R2"
  } else if(length( a_priori_filepaths_brazil[k][grep("a_priori_reduced_VIF", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_VIF
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_reduced_VIF"
  } else if(length( a_priori_filepaths_brazil[k][grep("a_priori_covsel_combined", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined_apriori
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_covsel_combined"
  } else if(length( a_priori_filepaths_brazil[k][grep("a_priori_covsel_glm", a_priori_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm_apriori
    results_df_a_priori_brazil$predictor_set[k] = "a_priori_covsel_glm"
  } else {
    # handle the case where none of the conditions are met
    predictor_subset = NULL
  }
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  brazil_pred <- terra::rast(
    here::here(
      a_priori_filepaths_brazil[k]
    )
  )
  
  # Set the CRS
  crs(brazil_pred) <- "EPSG:4326"
  
  # Get map of brazil to project our model over
  brazil_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "South America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # Create title for the legend
    theme(legend.position = "right") +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", a_priori_filepaths_brazil[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in brazil -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -33.74219 & lat < 5.257959 & lon > -74.00205 & lon < -34.80547
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(brazil_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = brazil_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(brazil_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = brazil_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_a_priori_brazil$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_a_priori_brazil$lr[k] = c$`LR Chisq`
  results_df_a_priori_brazil$df[k] = c$Df
  results_df_a_priori_brazil$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_a_priori_brazil$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_a_priori_brazil$sensitivity[k] = sensitivity
  
  ##################################################
  # THRESHOLDED MAP
  ##################################################
  
  brazil_binary_optimal = brazil_pred >= optimal
  
  values_at_points <- terra::extract(x = brazil_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred, aes(fill = median)) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    ) +
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in a_priori_filepaths_brazil){



##############################################################################################################################
##########################################
# ALL 19 PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_all_19_brazil <- data.frame(matrix(nrow=length(all_19_filepaths_brazil), ncol=10))
colnames(results_df_all_19_brazil) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                        "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(all_19_filepaths_brazil)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(all_19_filepaths_brazil[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_all_19_brazil$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_all_19_brazil$model[k] = model_type
  results_df_all_19_brazil$locality[k] = "brazil"
  
  # check which set of predictors you need, based on the current file at hand
  # all 19 files
  if(length( all_19_filepaths_brazil[k][grep("all_19_predictors", all_19_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$all_19_all
    results_df_all_19_brazil$predictor_set[k] = "all_19"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_reduced_PCMV_R2", all_19_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_pcmv_r2 
    results_df_all_19_brazil$predictor_set[k] = "all_19_reduced_PCMV_R2"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_reduced_R2_VIF", all_19_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2_VIF
    results_df_all_19_brazil$predictor_set[k] = "all_19_reduced_R2_VIF"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_reduced_R2", all_19_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2
    results_df_all_19_brazil$predictor_set[k] = "all_19_reduced_R2"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_reduced_VIF", all_19_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_VIF
    results_df_all_19_brazil$predictor_set[k] = "all_19_reduced_VIF"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_covsel_combined", all_19_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined
    results_df_all_19_brazil$predictor_set[k] = "all_19_covsel_combined"
  }else if(length( all_19_filepaths_brazil[k][grep("all_19_covsel_glm", all_19_filepaths_brazil[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm
    results_df_all_19_brazil$predictor_set[k] = "all_19_covsel_glm"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  brazil_pred <- terra::rast(
    here::here(
      all_19_filepaths_brazil[k]
    )
  )
  
  # Set the CRS
  crs(brazil_pred) <- "EPSG:4326"
  
  # Get map of brazil to project our model over
  brazil_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "South America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # Create title for the legend
    theme(legend.position = "right") +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", all_19_filepaths_brazil[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in brazil -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -33.74219 & lat < 5.257959 & lon > -74.00205 & lon < -34.80547
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(brazil_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = brazil_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(brazil_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = brazil_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_all_19_brazil$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_all_19_brazil$lr[k] = c$`LR Chisq`
  results_df_all_19_brazil$df[k] = c$Df
  results_df_all_19_brazil$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_all_19_brazil$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_all_19_brazil$sensitivity[k] = sensitivity
  
  ##################################################
  # THRESHOLDED MAP
  ##################################################
  
  brazil_binary_optimal = brazil_pred >= optimal
  
  values_at_points <- terra::extract(x = brazil_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred, aes(fill = median)) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 

    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in all_19_filepaths_brazil){



##############################################################################################################################
##########################################
# LITERATURE PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_literature_brazil <- data.frame(matrix(nrow=length(literature_filepaths_brazil), ncol=10))
colnames(results_df_literature_brazil) <- c("locality","optimal", "sensitivity", "suit_score_coef", 
                                            "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(literature_filepaths_brazil)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(literature_filepaths_brazil[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_literature_brazil$file[k] = file_name_clean
  model_type = sub(".*/", "", file_name_clean)
  results_df_literature_brazil$model[k] = model_type
  results_df_literature_brazil$locality[k] = "brazil"
  
  # check which set of predictors you need, based on the current file at hand
  
  if(length( literature_filepaths_brazil[k][grep("wang", literature_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$wang
    results_df_literature_brazil$predictor_set[k] = "wang"
  }else if(length( literature_filepaths_brazil[k][grep("aidoo", literature_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$aidoo
    results_df_literature_brazil$predictor_set[k] = "aidoo"
  }else if(length( literature_filepaths_brazil[k][grep("fordjour", literature_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$fordjour
    results_df_literature_brazil$predictor_set[k] = "fordjour"
  }else if(length( literature_filepaths_brazil[k][grep("naroui", literature_filepaths_brazil[k])]) ){
    predictor_subset = predictor_sets$naroui_khandan
    results_df_literature_brazil$predictor_set[k] = "naroui_khandan"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  brazil_pred <- terra::rast(
    here::here(
      literature_filepaths_brazil[k]
    )
  )
  
  # Set the CRS
  crs(brazil_pred) <- "EPSG:4326"
  
  # Get map of brazil to project our model over
  brazil_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "South America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # Create title for the legend
    theme(legend.position = "right") +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", literature_filepaths_brazil[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in brazil -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > -33.74219 & lat < 5.257959 & lon > -74.00205 & lon < -34.80547
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(brazil_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = brazil_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(brazil_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = brazil_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_literature_brazil$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_literature_brazil$lr[k] = c$`LR Chisq`
  results_df_literature_brazil$df[k] = c$Df
  results_df_literature_brazil$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_literature_brazil$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_literature_brazil$sensitivity[k] = sensitivity
  
  ##################################################
  # THRESHOLDED MAP
  ##################################################
  
  brazil_binary_optimal = brazil_pred >= optimal
  
  values_at_points <- terra::extract(x = brazil_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = brazil_pred, aes(fill = median)) +
    geom_sf(data = brazil_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-74, -34.81),
      ylim = c(-33.74, 5.26),
      crs = 4326,
      expand = FALSE
    ) + 
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in literature_filepaths_brazil){


results_df_a_priori_brazil
results_df_all_19_brazil
results_df_literature_brazil

brazil_results = rbind(results_df_a_priori_brazil,
                       results_df_all_19_brazil,
                       results_df_literature_brazil)



##############################################################################################################################
##########################################
# USA VALIDATIONS
##########################################
##############################################################################################################################

a_priori_filepaths_usa
all_19_filepaths_usa
literature_filepaths_usa

##############################################################################################################################
##########################################
# A PRIORI SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_a_priori_usa <- data.frame(matrix(nrow=length(a_priori_filepaths_usa), ncol=10))
colnames(results_df_a_priori_usa) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                          "lr", "df", "pr", "file", "predictor_set", "model")

##########################################
# For loop that runs through each result
# file, and extracts the relevant
# statistical values and saves the map
# for current suitability, and the probability curve
##########################################

for(k in 1:length(a_priori_filepaths_usa)){

  #print(a_priori_filepaths_usa[[k]])
  file_name_clean <- stringr::str_remove(stringr::str_remove(a_priori_filepaths_usa[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_a_priori_usa$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_a_priori_usa$model[k] = model_type
  results_df_a_priori_usa$locality[k] = "usa"
  
  # check which set of predictors you need, based on the current file at hand
  if(length( a_priori_filepaths_usa[k][grep("a_priori_all", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_all
    results_df_a_priori_usa$predictor_set[k] = "a_priori_all"
  } else if(length( a_priori_filepaths_usa[k][grep("a_priori_reduced_PCMV_R2", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_pcmv_r2
    results_df_a_priori_usa$predictor_set[k] = "a_priori_reduced_PCMV_R2"
  } else if(length( a_priori_filepaths_usa[k][grep("a_priori_reduced_R2", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_r2
    results_df_a_priori_usa$predictor_set[k] = "a_priori_reduced_R2"
  } else if(length( a_priori_filepaths_usa[k][grep("a_priori_reduced_VIF", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$a_priori_reduced_VIF
    results_df_a_priori_usa$predictor_set[k] = "a_priori_reduced_VIF"
  } else if(length( a_priori_filepaths_usa[k][grep("a_priori_covsel_combined", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined_apriori
    results_df_a_priori_usa$predictor_set[k] = "a_priori_covsel_combined"
  } else if(length( a_priori_filepaths_usa[k][grep("a_priori_covsel_glm", a_priori_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm_apriori
    results_df_a_priori_usa$predictor_set[k] = "a_priori_covsel_glm"
  } else {
    # handle the case where none of the conditions are met
    predictor_subset = NULL
  }
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  usa_pred <- terra::rast(
    here::here(
      a_priori_filepaths_usa[k]
    )
  )
  
  # Set the CRS
  crs(usa_pred) <- "EPSG:4326"
  
  # Get map of usa to project our model over
  usa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "North America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred) +
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
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
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", a_priori_filepaths_usa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in usa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > 10 & lon < -60
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(usa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = usa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(usa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = usa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_a_priori_usa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_a_priori_usa$lr[k] = c$`LR Chisq`
  results_df_a_priori_usa$df[k] = c$Df
  results_df_a_priori_usa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_a_priori_usa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_a_priori_usa$sensitivity[k] = sensitivity
  
  ############################################
  # THRESHOLDED MAP
  ############################################
  
  usa_binary_optimal = usa_pred >= optimal
  
  values_at_points <- terra::extract(x = usa_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred, aes(fill = median)) +
    
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
  
}#for(k in a_priori_filepaths_usa){



##############################################################################################################################
##########################################
# ALL 19 PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_all_19_usa <- data.frame(matrix(nrow=length(all_19_filepaths_usa), ncol=10))
colnames(results_df_all_19_usa) <- c("locality", "optimal", "sensitivity", "suit_score_coef", 
                                        "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(all_19_filepaths_usa)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(all_19_filepaths_usa[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_all_19_usa$file[k] = file_name_clean
  
  model_type = sub(".*/", "", file_name_clean)
  results_df_all_19_usa$model[k] = model_type
  results_df_all_19_usa$locality[k] = "usa"
  
  # check which set of predictors you need, based on the current file at hand
  # all 19 files
  if(length( all_19_filepaths_usa[k][grep("all_19_predictors", all_19_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$all_19_all
    results_df_all_19_usa$predictor_set[k] = "all_19"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_reduced_PCMV_R2", all_19_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_pcmv_r2 
    results_df_all_19_usa$predictor_set[k] = "all_19_reduced_PCMV_R2"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_reduced_R2_VIF", all_19_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2_VIF
    results_df_all_19_usa$predictor_set[k] = "all_19_reduced_R2_VIF"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_reduced_R2", all_19_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_r2
    results_df_all_19_usa$predictor_set[k] = "all_19_reduced_R2"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_reduced_VIF", all_19_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$all_19_reduced_VIF
    results_df_all_19_usa$predictor_set[k] = "all_19_reduced_VIF"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_covsel_combined", all_19_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$covsel_combined
    results_df_all_19_usa$predictor_set[k] = "all_19_covsel_combined"
  }else if(length( all_19_filepaths_usa[k][grep("all_19_covsel_glm", all_19_filepaths_usa[k])]) ) {
    predictor_subset = predictor_sets$covsel_glm
    results_df_all_19_usa$predictor_set[k] = "all_19_covsel_glm"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  usa_pred <- terra::rast(
    here::here(
      all_19_filepaths_usa[k]
    )
  )
  
  # Set the CRS
  crs(usa_pred) <- "EPSG:4326"
  
  # Get map of usa to project our model over
  usa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "North America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred) +
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # Create title for the legend
    theme(legend.position = "right") +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", all_19_filepaths_usa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in usa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > 10 & lon < -60
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(usa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = usa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(usa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = usa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_all_19_usa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_all_19_usa$lr[k] = c$`LR Chisq`
  results_df_all_19_usa$df[k] = c$Df
  results_df_all_19_usa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_all_19_usa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_all_19_usa$sensitivity[k] = sensitivity
  
  ############################################
  # THRESHOLDED MAP
  ############################################
  
  usa_binary_optimal = usa_pred >= optimal
  
  values_at_points <- terra::extract(x = usa_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred, aes(fill = median)) +
    
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in all_19_filepaths_usa){



##############################################################################################################################
##########################################
# LITERATURE PREDICTOR SET OF RESULTS
##########################################
##############################################################################################################################

predictor_subset = NULL

results_df_literature_usa <- data.frame(matrix(nrow=length(literature_filepaths_usa), ncol=10))
colnames(results_df_literature_usa) <- c("locality","optimal", "sensitivity", "suit_score_coef", 
                                            "lr", "df", "pr", "file", "predictor_set", "model")

for(k in 1:length(literature_filepaths_usa)){
  
  file_name_clean <- stringr::str_remove(stringr::str_remove(literature_filepaths_usa[[k]], 
                                                             "FINAL_RESULTS/.*?/"), 
                                         "/Diaphorina_citri/2010_ensembled\\.grd$")
  
  results_df_literature_usa$file[k] = file_name_clean
  model_type = sub(".*/", "", file_name_clean)
  results_df_literature_usa$model[k] = model_type
  results_df_literature_usa$locality[k] = "usa"
  
  # check which set of predictors you need, based on the current file at hand
  
  if(length( literature_filepaths_usa[k][grep("wang", literature_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$wang
    results_df_literature_usa$predictor_set[k] = "wang"
  }else if(length( literature_filepaths_usa[k][grep("aidoo", literature_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$aidoo
    results_df_literature_usa$predictor_set[k] = "aidoo"
  }else if(length( literature_filepaths_usa[k][grep("fordjour", literature_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$fordjour
    results_df_literature_usa$predictor_set[k] = "fordjour"
  }else if(length( literature_filepaths_usa[k][grep("naroui", literature_filepaths_usa[k])]) ){
    predictor_subset = predictor_sets$naroui_khandan
    results_df_literature_usa$predictor_set[k] = "naroui_khandan"
  }else predictor_subset = NULL 
  
  # specify the reduced predictor set
  reduced_pred = terra::subset(
    x = predictors,               # SpatRast containing WORLDCLIM layers 
    subset = predictor_subset
  )
  
  climPred <- predictor_subset
  
  # Import the ensemble MaxEnt raster prediction
  usa_pred <- terra::rast(
    here::here(
      literature_filepaths_usa[k]
    )
  )
  
  # Set the CRS
  crs(usa_pred) <- "EPSG:4326"
  
  # Get map of usa to project our model over
  usa_ext <- rnaturalearth::ne_countries(scale = "medium",
                                            returnclass = "sf") %>%
    dplyr::filter(continent == "North America")
  
  
  #######################################################
  # PLOT THE MAP STRAIGHT UP
  #######################################################
  
  # using the tidyterra package to plot the map
  current_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred) +
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # Create title for the legend
    theme(legend.position = "right") +
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
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  image_address <- gsub("Diaphorina_citri/2010_ensembled.grd", "", literature_filepaths_usa[[k]])
  
  ggsave(paste(image_address, "current_suitability.svg", sep = ""), 
         plot = current_map, width = 8, height = 10)
  
  #print(ggplot())  # Render the plot
  print("Map plotted")
  
  #######################################################
  # GET GPS RECORDS
  #######################################################
  
  # Import the newly downloaded GPS records -> all occurrences across the world,
  # then filter to the USA
  species <- 
    readr::read_csv(
      here::here("./data/gps/Diaphorina_citri.csv")
    ) %>%
    dplyr::select(
      species,
      lat = decimalLatitude,
      lon = decimalLongitude
    ) %>%
    # Keep only the records in usa -> SPECIFY COUNTRY HERE
    dplyr::filter(
      lat > 10 & lon < -60
    )
  
  #######################################################
  # Convert GPS coords to a 'terra' SpatialVector 
  #######################################################
  
  pts_species <- terra::vect(species[, 2:3])
  
  #######################################################
  # Set CRS for GPS coords 
  #######################################################
  
  crs(pts_species) <- "EPSG:4326"
  
  #######################################################
  # Check the CRS for the points = CRS for the raster 
  #######################################################
  
  crs(pts_species) == crs(usa_pred)
  
  #######################################################
  # Extract MaxEnt suitability scores at each 
  # recorded GPS point
  #######################################################
  
  clim_pred_gps <- 
    terra::extract(
      x = usa_pred,
      y = pts_species,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2 # assign the name of the second column ("median") to "suit_score"
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 1)
  
  #######################################################
  # Coerce GPS records into SPDF
  #######################################################
  
  recordsSpatialInv <- sp::SpatialPointsDataFrame(
    coords = cbind(species$lon, species$lat),
    data = species,
    proj4string = CRS(
      '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    )
  )
  
  # Select KG ecoregions in which there is at least one GPS record
  ecoContainInv <- kg_layer[recordsSpatialInv, ]
  
  ecoContainInv <- sf::st_as_sf(ecoContainInv)
  
  # get shape files for the specified predictors
  bg_area_inv <- terra::mask(reduced_pred, ecoContainInv)
  
  set.seed(2023)
  
  bg_points_inv <- terra::spatSample(
    x = bg_area_inv,        # Raster of background area to sample points from 
    size = 1000,        # How many background points do we want?
    method = "random",  # Random points
    replace = FALSE,    # Sample without replacement
    na.rm = TRUE,       # Remove background points that have NA climate data
    as.df = TRUE,       # Return background points as data.frame object
    xy = TRUE           # Return lat/lon values for each background point
    #cells = TRUE       # Return the cell numbers in which the background points fall
  ) %>%
    # Rename lon and lat columns to be consistent with GPS data for focal species 
    dplyr::rename(
      lon = x,
      lat = y
    )
  
  # Convert GPS coords to a 'terra' SpatialVector 
  bgptsInv <- terra::vect(bg_points_inv)
  
  # Set CRS for GPS coords 
  crs(bgptsInv) <- "EPSG:4326"
  
  # Check the CRS for the points = CRS for the raster 
  crs(bgptsInv) == crs(usa_pred)
  
  # Extract MaxEnt suitability scores at each background GPS point
  bg_clim_predInv <- 
    terra::extract(
      x = usa_pred,
      y = bgptsInv,
      xy = TRUE,
      na.rm = TRUE
    ) %>%
    # Clean df
    dplyr::select(
      lat = y,
      lon = x,
      suit_score = 2
    ) %>%
    tidyr::drop_na(suit_score) %>%
    dplyr::mutate(pres = 0)
  
  ##############################################################
  # Combine background and GPS points for occurrences  
  ##############################################################
  
  # Combine dataframes 
  clim_data <- 
    dplyr::bind_rows(
      clim_pred_gps,
      bg_clim_predInv
    )
  
  ##############################################################
  # Fit statistical model  
  ##############################################################
  
  # Run a simple binomial GLM
  # Does the MaxEnt suitability score significantly affect the presence of the psyllid?
  # This compares the suitability scores at the locations where the psyillid was recorded, to
  # background points where it shouldn't occur
  
  mod1 <- glm(
    # Response variable
    pres ~ 
      # Fixed effects 
      suit_score, 
    data = clim_data,
    family = binomial(link = "logit")
  )
  
  print(paste("GLM coefficient = ", mod1$coefficients[2]) ) # suit score
  #print( summary(mod1) )
  results_df_literature_usa$suit_score_coef[k] = mod1$coefficients[2]
  
  # Check model diagnostics 
  # DHARMa::simulateResiduals(fittedModel = mod1, plot = TRUE)
  
  # Test parameter significance 
  c = car::Anova(
    mod1, 
    test = "LR",
    type = "II"
  )
  
  print(paste("LR chisq = ", c$`LR Chisq`) )
  print(paste("DF = ", c$Df) )
  print(paste("Pr = ", c$`Pr(>Chisq)`) )
  
  results_df_literature_usa$lr[k] = c$`LR Chisq`
  results_df_literature_usa$df[k] = c$Df
  results_df_literature_usa$pr[k] = c$`Pr(>Chisq)`
  
  ##############################################################
  # Plot statistical model  
  ##############################################################
  
  # Extract model predictions 
  preds <- ggeffects::ggpredict(
    mod1, 
    terms = c("suit_score [0:1 by = 0.01]")) %>%
    as.data.frame() %>%
    dplyr::rename(
      suit_score = x
    )
  head(preds)
  
  clim_data_presences = dplyr::filter(clim_data, pres == 1)
  clim_data_absences = dplyr::filter(clim_data, pres == 0)
  
  # Plot model predictions 
  prob_curve = preds %>%
    ggplot2::ggplot(data = ., 
                    aes(x = suit_score,
                        y = predicted)) +
    geom_rug(data = clim_data_presences, aes(x= suit_score, y = pres), sides = "t") + 
    geom_rug(data = clim_data_absences, aes(x= suit_score, y = pres), sides = "b") +
    geom_line() +
    geom_ribbon(aes(ymin = conf.low,
                    ymax = conf.high),
                alpha = 0.1) +
    labs(
      x = "MaxEnt suitability score",
      y = "Probability of establishment"
    )
  
  ggsave(paste(image_address, "prob_curve.svg", sep = ""), 
         plot = prob_curve, width = 10, height = 7)
  
  ##############################################################
  # Calculate model accuracy metrics 
  ##############################################################
  
  #devtools::install_github("selva86/InformationValue")
  
  # Use model to predict probability of default
  predicted <- terra::predict(
    mod1, 
    clim_data, 
    type = "response"
  ) 
  #head(predicted)
  
  # Find optimal cutoff probability to use to maximize accuracy
  optimal <- InformationValue::optimalCutoff(
    clim_data$pres, 
    #optimiseFor = "Ones",
    predicted)[1]
  
  print(paste("OPTIMAL VALUE = ", optimal)) 
  
  results_df_literature_usa$optimal[k] = optimal
  
  # Calculate sensitivity -> MOST IMPORTANT METRIC
  # - Percentage of true positives 
  sensitivity = InformationValue::sensitivity(
    actuals = as.factor(clim_data$pres),
    predicted = predicted, 
    threshold = optimal
    #threshold = 0.5
  )
  
  print(paste("SENSITIVITY = ", sensitivity)) 
  
  results_df_literature_usa$sensitivity[k] = sensitivity
  
  ############################################
  # THRESHOLDED MAP
  ############################################
  
  usa_binary_optimal = usa_pred >= optimal
  
  values_at_points <- terra::extract(x = usa_binary_optimal, # thresholded map
                                     y = pts_species, 
                                     xy = TRUE,
                                     na.rm = TRUE) %>%
    dplyr::mutate(median = dplyr::if_else(median == TRUE, 1, 0))
  
  values_at_points = na.omit(values_at_points)
  
  thresholded_map = ggplot() +
    tidyterra::geom_spatraster(data = usa_pred, aes(fill = median)) +
    
    geom_sf(data = usa_ext, fill = NA, color = "black", size = 0.2) +
    
    scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(optimal, 1) # here we set the optimal threshold
    ) +
    # Control axis and legend labels 
    labs(
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    
    coord_sf(
      xlim = c(-160, -55),
      ylim = c(0, 60),
      crs = 4326,
      expand = FALSE
    ) +
    
    # Create title for the legend
    theme(legend.position = "right" ) +
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
    geom_point(data = values_at_points, aes(x = x, y = y, color = as.factor(median) ), cex = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(color = "Present")
  
  ggsave(paste(image_address, "thresholded_map.svg", sep = ""), 
         plot = thresholded_map, width = 8, height = 10)
  
}#for(k in literature_filepaths_usa){


results_df_a_priori_usa
results_df_all_19_usa
results_df_literature_usa

usa_results = rbind(results_df_a_priori_usa,
                       results_df_all_19_usa,
                       results_df_literature_usa)


##########################################################
# POST PROCESSING
##########################################################

all_results = rbind(africa_results, 
                    brazil_results, 
                    usa_results)

all_results

# fill out the model setting parameters in a new column

all_results$model_settings = ""

for(k in 1:nrow(all_results)){
  
  if( (all_results$predictor_set[k] == "all_19_reduced_R2" || 
       all_results$predictor_set[k] == "all_19_reduced_R2_VIF") && 
       all_results$model[k] == "tuned_maxent" ) #if
       all_results$model_settings[k] = "LQH1"
  
  if(all_results$predictor_set[k] == "all_19_reduced_VIF" && 
     all_results$model[k] == "tuned_maxent") all_results$model_settings[k] = "H1"
  
  if( (all_results$predictor_set[k] == "a_priori_reduced_PCMV_R2" || 
       all_results$predictor_set[k] == "a_priori_covsel_glm" ||
       all_results$predictor_set[k] == "all_19_covsel_glm") && 
     all_results$model[k] == "tuned_maxent") all_results$model_settings[k] = "LQH1"
  
  if(all_results$model_settings[k] == "") all_results$model_settings[k] = "LQH2"
  
}#for

all_results$locality = as.factor(all_results$locality)
all_results$predictor_set = factor(all_results$predictor_set, 
                          levels = c("a_priori_all", "a_priori_reduced_PCMV_R2", 
                                     "a_priori_reduced_R2", 
                                     "a_priori_reduced_VIF",
                                     "a_priori_covsel_combined",
                                     "a_priori_covsel_glm",
                                     "all_19", "all_19_reduced_PCMV_R2",
                                     "all_19_reduced_R2",
                                     "all_19_reduced_R2_VIF", "all_19_reduced_VIF",
                                     "all_19_covsel_combined",
                                     "all_19_covsel_glm",
                                     "aidoo", "fordjour", "naroui_khandan", "wang"))
all_results$model = as.factor(all_results$model)

str(all_results)

write.csv(all_results, file = "all_results_output_table.csv", quote = FALSE, row.names = FALSE)

levels(all_results$predictor_set)

################
# optimal values
################

ggplot2::ggplot(data = subset(all_results, locality == "africa"), 
                aes(x = predictor_set, y = optimal,  shape = model)) +
  geom_point(size = 2) + 
  #scale_color_manual(values = c("black", "steelblue")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Optimal value", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")

##################
# all localities
##################

# AS BARS
ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = optimal,  fill = model)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_fill_manual(values = c("lightgrey", "forestgreen")) +
  #scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Optimal value", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality, labeller = labeller(locality = c(africa = "Africa", brazil = "Brazil", usa = "USA")))

# AS POINTS
ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = optimal,  shape = model)) +
  geom_point(size = 2) + 
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Optimal value", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality, labeller = labeller(locality = c(africa = "Africa", brazil = "Brazil", usa = "USA")))


################
# sensitivity
################

ggplot2::ggplot(data = subset(all_results, locality == "africa"), 
                aes(x = predictor_set, y = sensitivity,  shape = model)) +
  geom_point(size = 2) + 
  #scale_color_manual(values = c("black", "steelblue")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Sensitivity", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")

##################
# all localities
##################

# PLOT AS POINTS
ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = sensitivity, shape = model)) +
  geom_point(size = 2) + 
 # scale_colour_manual(values = c("black", "forestgreen")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Sensitivity", x = "Predictor set", shape = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality, labeller = labeller(locality = c(africa = "Africa", brazil = "Brazil", usa = "USA")))

# PLOT AS BARS
ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = sensitivity, fill = model)) +
  #geom_point(size = 2) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_fill_manual(values = c("lightgrey", "forestgreen")) +
  #scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Sensitivity", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality, labeller = labeller(locality = c(africa = "Africa", brazil = "Brazil", usa = "USA")))


################
# suit_score coefficient
################

ggplot2::ggplot(data = subset(all_results, locality == "africa"), 
                aes(x = predictor_set, y = suit_score_coef,  shape = model)) +
  geom_point(size = 2) + 
  #scale_color_manual(values = c("black", "steelblue")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Suitability score coefficient", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")

##################
# all localities
##################

ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = suit_score_coef,  shape = model)) +
  geom_point(size = 2) + 
  #scale_color_manual(values = c("black", "steelblue")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Suitability score coefficient", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality)

################
# likelihood ratio
################

ggplot2::ggplot(data = subset(all_results, locality == "africa"), 
                aes(x = predictor_set, y = lr,  shape = model)) +
  geom_point(size = 2) + 
  #scale_color_manual(values = c("black", "steelblue")) +
  scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Likelihood ratio", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right")

##################
# all localities
##################

ggplot2::ggplot(data = all_results, 
                aes(x = predictor_set, y = lr,  fill = model)) +
  #geom_point(size = 2) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("black", "steelblue")) +
  #scale_shape_manual(values = c(16, 3)) +
  labs(title = "",
       y = "Likelihood ratio", x = "Predictor set", colour = "MaxEnt model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
  facet_wrap(~locality, labeller = labeller(locality = c(africa = "Africa", brazil = "Brazil", usa = "USA")) )
