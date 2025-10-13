
rm(list = ls())          # remove all variables of the work space

# Import packages
library(terra)
library(tidyverse)
library(sf)
library(stars)

# Area
study_area <- ext(2, 6, 42.3, 43.7) 

# Functions

combine_ncfile <- function(path_non_int, path_interim){
  # Combine non-interim and interim periods
  ## Open the .nc files
  raster <- rast(path_non_int)
  raster_int <- rast(path_interim)
  
  ## Give a number to each month to rename the columns
  i_max <- nlyr(raster)
  i_max_int <- nlyr(raster_int)
  new_names <- as.character(1:i_max)
  new_names_int <- as.character(1:i_max_int + i_max)
  
  names(raster) <- new_names
  names(raster_int) <- new_names_int
  ## Combine the two rasters 
  raster_combine <- c(raster, raster_int)
  return(raster_combine)
}


#  Dynamic covariates ----

## SST ----

# Import Sea Surface Temperature (SST) and average by season 
raster_sst <- combine_ncfile("1.data/covariates/copernicus/SST.nc", "1.data/covariates/copernicus/SST_int.nc") %>% 
  terra::crop(study_area) %>% 
  terra::extend(study_area) 

# Average seasonal SST across all years
mean_season <- raster_sst %>% as_tibble() %>% 
  bind_cols(crds(raster_sst)) %>% 
  pivot_longer(-c(x, y), names_to = "date_nb", values_to = "value") %>% 
  mutate(month_nb = (as.numeric(date_nb)-1)%%12 + 1) %>% 
  mutate(season_nb =  ceiling(month_nb/3)) %>% 
  group_by(x, y, season_nb) %>%  
  summarize(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE)
  ) %>% 
  ungroup()


# Create one col for each season and compute all season mean SST
raster_mean_SST <- mean_season %>% 
  pivot_wider(names_from = season_nb, values_from = c(mean_value, sd_value)) %>% 
  rename(mean_winter_SST = mean_value_1,
         mean_spring_SST = mean_value_2,
         mean_summer_SST = mean_value_3,
         mean_autumn_SST = mean_value_4,
         sd_winter_SST = sd_value_1,
         sd_spring_SST = sd_value_2,
         sd_summer_SST = sd_value_3,
         sd_autumn_SST = sd_value_4) %>% 
  terra::rast(type="xyz", crs="", digits=6, extent=NULL) 

all_SST <- c(app(raster_sst, mean),app(raster_sst, sd)) %>% 
  resample(raster_mean_SST)

names(all_SST) <- c("mean_SST", "sd_SST")
raster_mean_SST <- c(raster_mean_SST, all_SST)

plot(raster_mean_SST)

rm(mean_season, raster_sst, all_SST)


## Chloropyll a ----

CHL <- combine_ncfile("1.data/covariates/copernicus/CHL.nc", "1.data/covariates/copernicus/CHL_int.nc") 

## Salinity ----
SAL <- combine_ncfile("1.data/covariates/copernicus/SAL.nc", "1.data/covariates/copernicus/SAL_int.nc") 

## Sea surface height ----
SSH <- combine_ncfile("1.data/covariates/copernicus/SSH.nc", "1.data/covariates/copernicus/SSH_int.nc") 

## Velocity ---- 
VEL_North <- combine_ncfile("1.data/covariates/copernicus/VEL_North.nc", "1.data/covariates/copernicus/VEL_North_int.nc")
VEL_East <- combine_ncfile("1.data/covariates/copernicus/VEL_East.nc", "1.data/covariates/copernicus/VEL_East_int.nc")


raster_dynamic_cov <- c(app(CHL, mean), app(CHL, sd),
                        app(SAL, mean), app(SAL, sd),
                        app(SSH, mean), app(SSH, sd),
                        app(VEL_North, mean), app(VEL_North, sd),
                        app(VEL_East, mean), app(VEL_East, sd)
) %>% 
  terra::crop(raster_mean_SST) %>% 
  terra::extend(raster_mean_SST) 


names(raster_dynamic_cov) <- c("mean_CHL", "sd_CHL",
                               "mean_SAL", "sd_SAL",
                               "mean_SSH", "sd_SSH",
                               "mean_VEL_North", "sd_VEL_North",
                               "mean_VEL_East", "sd_VEL_East"
)
plot(raster_dynamic_cov)

rm(VEL_East, VEL_North, SSH, SAL, CHL)

# Static covariates ---- 

## Bathymetry in m ----
depth_raster <- rast("1.data/covariates/marspec/bathy_30s") %>% 
  raster::crop(study_area) %>% 
  as.numeric()

## Distance to shore, slope and concavity ----

# - distance to shore in km (biogeo_05), 
# - slope in degrees (biogeo_06), 
# - concavity in degrees (biogeo_07)


file_nb <- c(5, 6, 7)
for (i in file_nb){
  depth_raster <- rast(paste0("1.data/covariates/marspec/biogeo01_07_30s/biogeo0",i,"_30s")) %>% 
    crop(study_area) %>% 
    as.numeric() %>% 
    c(depth_raster)
}

rm(i, file_nb)

names(depth_raster) <- c("concavity", "slope", "dist_to_shore", "bathymetry")
varnames(depth_raster) <- c("concavity", "slope", "dist_to_shore", "bathymetry")


crs(raster_mean_SST) <- crs(depth_raster)

# Use resample to have same extent and resolution
raster_static_cov <- depth_raster %>% 
  terra::resample(raster_mean_SST) 

plot(raster_static_cov)

rm(depth_raster)

# Combine all the rasters ----
combined_rasters <- c(raster_mean_SST, raster_static_cov, raster_dynamic_cov) 
plot(combined_rasters)

rm(raster_mean_SST, raster_dynamic_cov, raster_static_cov)

## Add log value ----

log_raster <- log(combined_rasters)
names(log_raster) <- paste0("log_", names(log_raster))

combined_rasters_log <- log_raster[[c("log_mean_CHL")]]

layers_to_keep <- c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
                    "mean_SST", "slope", "dist_to_shore", "bathymetry", 
                    "mean_SAL", "mean_SSH", "mean_VEL_North", "mean_VEL_East")

combined_rasters_red <- combined_rasters[[c(layers_to_keep)]]


final_raster <- c(combined_rasters_log, combined_rasters_red) 

rm(combined_rasters, layers_to_keep, combined_rasters_log, combined_rasters_red, log_raster)

terra::writeRaster(final_raster, filename = "1.data/covariates/all_covariates.tif", overwrite=TRUE)

# Create covariates grid ----
rm(list = ls())

load("1.data/spatial_objects/grid.rdata")
raster_cov <- read_stars("1.data/covariates/all_covariates.tif")

grid <- grid %>% st_transform(st_crs(raster_cov))

polygon_values <- split(raster_cov, "band") %>% 
  st_extract(at = st_centroid(grid), mean) %>% 
  split("band") 

covariates <- st_sf(grid) %>%
  mutate(polygon_values$band) %>%
  st_drop_geometry() %>% 
  left_join(grid, by = "grid_id") %>%
  st_as_sf() %>%
  st_transform(2154)


save(covariates, file = "1.data/covariates.rdata")


