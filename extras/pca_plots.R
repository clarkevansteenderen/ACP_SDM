source("predictor_sets.R")

# How many points in the USA?

species_usa <- sp_data %>%
  dplyr::filter(lon > -127, lon< -65, lat > 10, lat < 50)

nrow(species_usa)

# How many points in Brazil?

species_braz <- sp_data %>%
  dplyr::filter(lon > -73.99, lon < -34.80, lat > -33.75, lat < 5.27)

nrow(species_braz)

# How many points in Africa?

species_afr <- sp_data %>%
  dplyr::filter(lon > -26, lon < 55, lat > -48, lat < 40)

nrow(species_afr)


#############################
# USA climate data - current
#############################

USA_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_usa$lon, species_usa$lat) )
)

head(USA_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
USA_speciesSwd <- cbind(species_usa, USA_speciesEnv)


# Process data 
USA_gps_mega <- as.data.frame(USA_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(USA_gps_mega)

# Reproject CRS 
USA_gpsSpatial <- USA_gps_mega 
coordinates(USA_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(USA_gpsSpatial, "proj4string") <- crs_wgs84

#############################
# BRAZIL climate data - current
#############################

BRAZ_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_braz$lon, species_braz$lat) )
)

head(BRAZ_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
BRAZ_speciesSwd <- cbind(species_braz, BRAZ_speciesEnv)


# Process data 
BRAZ_gps_mega <- as.data.frame(BRAZ_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(BRAZ_gps_mega)

# Reproject CRS 
BRAZ_gpsSpatial <- BRAZ_gps_mega 
coordinates(BRAZ_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(BRAZ_gpsSpatial, "proj4string") <- crs_wgs84


#############################
# AFRICA climate data - current
#############################

AFR_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_afr$lon, species_afr$lat) )
)

head(AFR_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
AFR_speciesSwd <- cbind(species_afr, AFR_speciesEnv)


# Process data 
AFR_gps_mega <- as.data.frame(AFR_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(AFR_gps_mega)

# Reproject CRS 
AFR_gpsSpatial <- AFR_gps_mega 
coordinates(AFR_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(AFR_gpsSpatial, "proj4string") <- crs_wgs84


#############################
# USA climate data - current
#############################

USA_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_usa$lon, species_usa$lat) )
)

head(USA_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
USA_speciesSwd <- cbind(species_usa, USA_speciesEnv)


# Process data 
USA_gps_mega <- as.data.frame(USA_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(USA_gps_mega)

# Reproject CRS 
USA_gpsSpatial <- USA_gps_mega 
coordinates(USA_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(USA_gpsSpatial, "proj4string") <- crs_wgs84

#############################
# BRAZIL climate data - current
#############################

BRAZ_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_braz$lon, species_braz$lat) )
)

head(BRAZ_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
BRAZ_speciesSwd <- cbind(species_braz, BRAZ_speciesEnv)


# Process data 
BRAZ_gps_mega <- as.data.frame(BRAZ_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(BRAZ_gps_mega)

# Reproject CRS 
BRAZ_gpsSpatial <- BRAZ_gps_mega 
coordinates(BRAZ_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(BRAZ_gpsSpatial, "proj4string") <- crs_wgs84


#############################
# NATIVE RANGE climate data - current
#############################

NATIVE_speciesEnv <- base::data.frame(
  raster::extract(predictors, cbind(species_native$lon, species_native$lat) )
)

head(NATIVE_speciesEnv)

# Combine thinned records with environmental data Swd = samples with data
NATIVE_speciesSwd <- cbind(species_native, NATIVE_speciesEnv)


# Process data 
NATIVE_gps_mega <- as.data.frame(NATIVE_speciesSwd) %>%
  # Keep only the lat/lon columns and key environmental variables
  dplyr::select(
    decimalLatitude = lat,
    decimalLongitude = lon,
    predictor_sets$all_19_all
  ) %>%
  dplyr::mutate(species = "Diaphorina_citri") %>%
  dplyr::select(
    species,
    everything()
  ) %>%
  # Drop rows with NA 
  tidyr::drop_na()
head(NATIVE_gps_mega)

# Reproject CRS 
NATIVE_gpsSpatial <- NATIVE_gps_mega 
coordinates(NATIVE_gpsSpatial) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(NATIVE_gpsSpatial, "proj4string") <- crs_wgs84


########################################
# PCA PLOTS
########################################

NATIVE_gpsSpatial = as.data.frame(NATIVE_gpsSpatial)
NATIVE_gpsSpatial$region = "native"
USA_gpsSpatial = as.data.frame(USA_gpsSpatial)
USA_gpsSpatial$region = "USA"
AFR_gpsSpatial = as.data.frame(AFR_gpsSpatial)
AFR_gpsSpatial$region = "africa"
BRAZ_gpsSpatial = as.data.frame(BRAZ_gpsSpatial)
BRAZ_gpsSpatial$region = "brazil"

climate_data_current = rbind(NATIVE_gpsSpatial, USA_gpsSpatial, AFR_gpsSpatial, BRAZ_gpsSpatial)

climate_df = climate_data_current[4:22]

# Get the current column names
current_names <- colnames(climate_df)

# Create new column names
new_names <- paste0("bio_", seq_along(current_names))

# Set the new column names
colnames(climate_df) <- new_names

pca_res = prcomp(climate_df, scale. = TRUE)
loadings = pca_res$rotation
loadings

summary(pca_res)
eig.val = factoextra::get_eigenvalue(pca_res)
factoextra::fviz_eig(pca_res, col.var="blue")
var = factoextra::get_pca_var(pca_res)
corrplot(var$cos2, is.corr = FALSE)
factoextra::fviz_cos2(pca_res, choice = "var", axes = 1:2)

factoextra::fviz_pca_var(pca_res,
                         col.var = "cos2", # Color by the quality of representation
                         gradient.cols = c("darkorchid4", "gold", "darkorange"),
                         repel = TRUE)

factoextra::fviz_pca_var(pca_res,
                         col.var = "contrib", # Color by contributions to the PC
                         gradient.cols = c("darkorchid4", "gold", "darkorange"),
                         repel = TRUE)


factoextra::fviz_pca(pca_res, label = "var", col.var = "black", habillage = climate_data_current$region,
                     palette = c("#00AFBB", "#E7B800", "#FC4E07", "forestgreen"),
                     addEllipses = TRUE  ) + theme_classic()

# Contributions of variables to PC1
a<-factoextra::fviz_contrib(pca_res, choice = "var", axes = 1)
# Contributions of variables to PC2
b<-factoextra::fviz_contrib(pca_res, choice = "var", axes = 2)
grid.arrange(a,b, ncol=2, top='Contribution of the variables to the first two PCs')

pca_res_dataframe = data.frame("region" = climate_data_current$region, pca_res$x[,1:2])

clim_pca_plot = ggplot2::ggplot(data = pca_res_dataframe) +
  geom_point(aes(x = PC1, y = PC2, col = region)) + 
  theme_classic() 

clim_pca_plot_2 = clim_pca_plot +
  ggConvexHull::geom_convexhull(aes(x = PC1, y = PC2, fill = region, colour = region), alpha = 0.08) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept=0, linetype="dashed", colour= "black", linewidth = 0.8) +
  scale_fill_manual(values = c(USA = "black", africa = "darkorange1", brazil = "red", native = "forestgreen")) + 
  scale_color_manual(values = c(USA = "black", africa = "darkorange1", brazil = "red", native = "forestgreen")) 

clim_pca_plot_2




kmeans = factoextra::eclust(climate_df, k = 4) 

ggplot2::autoplot(cluster::fanny(climate_df, 4), frame = TRUE) +
  theme_classic()

ggplot2::autoplot(cluster::pam(climate_df, 4), frame = TRUE, frame.type = 'norm') +
  theme_classic()
