rm(list= ls())  

library(dplyr)
library(sf)
library(lubridate)

`%!in%` = Negate(`%in%`)

# DATA
load("1.data/spatial_objects/grid.rdata")

telemetry <- readRDS("1.data/telemetry/telemetry_brut.rds")

telemetry <- telemetry %>%
  rename(ind_ID = "individual_local_identifier",
         species = "individual_taxon_canonical_name",
         X = "location_long",
         Y = "location_lat",
         speed = "ground_speed",
         time = timestamp,
         satellite = "gps_satellite_count",
         hdop = "gps_hdop",
         altitude = "height_above_msl") %>%
  dplyr::select(ind_ID, species, X, Y, time, speed, satellite, hdop, altitude) %>%
  mutate_at("time", as_datetime) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = st_crs(grid)) %>% 
  mutate(lon = st_coordinates(.)[, 1],  
         lat = st_coordinates(.)[, 2]) %>%
  filter(time > as.Date("2016-12-31")) %>%
  st_intersection(grid %>% select(geometry)) %>%
  mutate(species = case_when(species == "Thalasseus sandvicensis" ~ "sterne caugek",
                             species == "Puffinus yelkouan" ~ "puffin yelkouan",
                             species == "Calonectris diomedea" ~ "puffin de scopoli",
                             species == "Larus michahellis" ~ "goeland leucophee"))
  


save(telemetry, file = "1.data/all_seabirds_telemetry.rdata")




