rm(list = ls())  

library(dplyr)
library(sf)
library(purrr)
library(ggplot2)
`%!in%` = Negate(`%in%`)


# Covariates
load("1.data/covariates.rdata")

xlim <- st_bbox(covariates)[c(1,3)]
ylim <- st_bbox(covariates)[c(2,4)]

# Coutour france
europe <- st_read("1.data/spatial_objects/Europe/Europe_merged.shp", quiet = TRUE) %>% 
  dplyr::filter(COUNTRY %in% c("France",  "Spain", "Italy"))%>%
  st_transform(2154) %>%
  st_union()

grid_buff_select <- covariates %>%
  st_union() %>%
  st_buffer(20000)

## Load Count data 
colonies_fr <- readxl::read_excel("1.data/covariates/colonies/ROMN 2020-2022_POUR CEFE.xlsx") %>%
  mutate_at(c("Colonie_longitude (degrés décimaux)", "Colonie_latitude (degrés décimaux)"), as.numeric) %>%
  rename( Lat = "Colonie_latitude (degrés décimaux)", 
          Lon = "Colonie_longitude (degrés décimaux)",
          EFF_MIN_synthese = "EFF_MIN synthèse",
          EFF_MAX_synthese = "EFF_MAX synthèse",
          Espece = ESPECE,
          Annee = ANNEE,
          Commune = COMMUNE,
          Secteur = SECTEUR,
          Dept = DEPARTEMENT) %>%
  filter(is.na(Lon) == FALSE) %>%
  st_as_sf(coords = c("Lon", "Lat"), remove = FALSE, crs = 4326) %>% 
  st_transform(crs = st_crs(grid_buff_select)) %>%
  st_intersection(grid_buff_select) %>%
  group_by(Lat, Lon) %>%
  mutate(ID_colonie = cur_group_id()) %>%
  ungroup() %>%
  mutate(EFF = case_when(is.na(EFF_MAX_synthese) == TRUE ~ EFF_MIN_synthese,
                         TRUE ~ EFF_MAX_synthese)) %>%
  mutate(EFF = case_when(is.na(EFF) == TRUE ~ as.numeric(EFF_MAX),
                         TRUE ~ as.numeric(EFF))) %>%
  filter(is.na(EFF) == FALSE) %>%
  filter(EFF != 0) %>%
  mutate(Espece = case_when(
    Espece == "Sterne pierregarin" ~ "sterne pierregarin",
    Espece == "Goéland leucophée"  ~ "goeland leucophee",
    Espece == "Mouette mélanocéphale" ~ "mouette melanocephale",
    Espece == "Mouette rieuse" ~ "mouette rieuse",
    Espece == "Sterne caugek" ~ "sterne caugek",
    Espece == "Océanite tempête" ~ "oceanite tempete",
    Espece == "Puffin de Scopoli" ~ "puffin de scopoli",
    Espece == "Puffin yelkouan" ~ "puffin yelkouan",
    TRUE ~ as.character(Espece))) %>%
  filter(Espece %in% c("sterne pierregarin", "oceanite tempete", "puffin de scopoli", "puffin yelkouan", "sterne caugek", "goeland leucophee", "mouette melanocephale","mouette rieuse")) %>%
  select(ID_colonie, Dept, Commune, Secteur, Lon, Lat, Espece, Annee, EFF) %>%
  group_by(ID_colonie, Espece) %>%
  filter(EFF == max(EFF)) %>%
  distinct(ID_colonie, Espece, .keep_all = TRUE) %>%
  ungroup(.)

rm(grid_buff_select)

save(colonies_fr, file = "1.data/colonies.rdata")
