rm(list = ls())

library(dplyr)
library(stringi)
library(stringr)
library(lubridate)
library(sf)
library(ggplot2)
library(readxl)
library(readr)
`%!in%` = Negate(`%in%`)

load("1.data/spatial_objects/grid.rdata")


# Pelmed  --------------------------

# 2017 -> 2020
load("1.data/count/PELMED/pelmed2017_2021.rdata")

## Observations -----------------------------

pelmed_obs <- pelobs %>%  
  select(lat, lon, podSize, transect, routeType, sighting, date, 
         hhmmss, famille_fr, nom_fr, geometry) %>% 
  mutate(nom_fr = stri_trans_general(str = nom_fr, id = "Latin-ASCII")) %>%  # remove accents
  mutate(nom_fr = tolower(nom_fr)) %>% # in lower case
  rename(effectif = podSize) %>% 
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) %>%
  st_transform(st_crs(grid)) %>%
  filter(nom_fr %!in% c("anatide ind.")) %>%
  mutate(nom_fr = case_when(nom_fr %in% c("alcide ind.") ~ "alcide ind",
                            nom_fr %in% c("cormoran ind.") ~ "cormoran ind",
                            nom_fr %in% c("oceanite ind.") ~ "oceanite ind",
                            nom_fr %in% c("grand goeland ind.") ~ "goeland ind",
                            nom_fr %in% c("labbe ind.") ~ "labbe ind",
                            nom_fr %in% c("petit puffin ind.") ~ "petit puffin ind",
                            nom_fr %in% c("puffin ind.") ~ "puffin ind",
                            nom_fr %in% c("puffin cendre / de scopoli", "grand puffin ind.") ~ "grand puffin ind",
                            nom_fr %in% c("laride ind.") ~ "laride ind",
                            nom_fr %in% c("mouette ind.") ~ "mouette ind",
                            nom_fr %in% c("plongeon ind.") ~ "plongeon ind",
                            nom_fr %in% c("sterne ind.") ~ "sterne ind",
                            TRUE ~ as.character(nom_fr)))

## Efforts -------------------

pelmed_eff <- peleff %>% 
  select(effort, date, hhmmss, seaState, lat, lon, legLengKm, geometry) %>% 
  mutate(year = year(date),
         transect_name = 1:nrow(.)) %>%
  st_transform(st_crs(grid))

rm(peleff, pelobs)

# Parc Naturel Marin  -----------------------

## Observations -------------

obs <- st_read(dsn = '1.data/count/PNM/obs_oiseaux_pt_4326.shp', quiet = TRUE) %>%
  mutate(date = dmy(date))  

# New data to add
obs1 <- read.csv("1.data/count/PNM/MegaObs_pri23_pt_wgs84.csv", sep = ";") %>%
  rename(age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs2 <- read.csv("1.data/count/PNM/MegaObs_ete23_pt_wgs84.csv", sep = ";") %>%
  rename(age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs3 <- read.csv("1.data/count/PNM/MegaObs_aut23_pt_wgs84.csv", sep = ";") %>%
  rename(transect = num.T,
         age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs4 <- read.csv("1.data/count/PNM/MegaObs_hiv24_pt_wgs84_vf.csv", sep = ";") %>%
  rename(transect = num.T) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs5 <- read.csv("1.data/count/PNM/MegaObs_pri24_pt_wgs84.csv", sep = ";") %>%
  rename(transect = num_T,
         age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs6 <- read.csv("1.data/count/PNM/MegaObs_ete24_pt_wgs84.csv", sep = ";") %>%
  rename(transect = num_T,
         age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  

obs7 <- read.csv("1.data/count/PNM/Megaobs_aut24_pt_wgs84.csv", sep = ";") %>%
  rename(transect = num_T,
         age = Oiseaux.Age) %>%
  select(any_of(names(obs))) %>%
  mutate(date = dmy(date))  


# Combine all datasets 
pnm_obs <- obs %>%
  select(any_of(names(obs1))) %>%
  st_drop_geometry() %>%
  rbind(obs1, obs2, obs3, obs4, obs5, obs6, obs7) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(st_crs(grid))

rm(obs1, obs3, obs4, obs5, obs6, obs7, obs)

## Efforts ---------------------------

effort <- st_read(dsn = '1.data/count/PNM/Megaobs_transect_ln_2154.shp', quiet = TRUE) %>%
  rename(date = date_tr,
         transect = id_trace,
         LengKm = long_km) 

# New data to add
eff1 <- st_read(dsn = '1.data/count/PNM/Megaobs_pri23_transect_ln_l93.shp', quiet = TRUE) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  select(any_of(names(effort))) 


eff2 <- st_read(dsn = '1.data/count/PNM/Megaobs_ete23_transect_ln_l93v3.shp', quiet = TRUE) %>%
  left_join(obs2 %>% select(transect, date) %>% distinct(transect, .keep_all = TRUE), by = "transect") %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  mutate(camp = "2023_ete") %>%
  select(any_of(names(effort)))

eff3 <- st_read(dsn = '1.data/count/PNM/Megaobs_aut23_transect_ln_l93.shp', quiet = TRUE) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  rename(transect = num.T,
         date = date_tr) %>%
  select(any_of(names(effort)))

eff4 <- st_read(dsn = '1.data/count/PNM/Megaobs_hiv24_transect_ln_l93.shp', quiet = TRUE) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  rename(transect = num.T) %>%
  select(any_of(names(effort)))

eff5 <- st_read(dsn = '1.data/count/PNM/Megaobs_pri24_transect_ln_l93.shp', quiet = TRUE) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  rename(transect = num_T) %>%
  select(any_of(names(effort)))

eff6 <- st_read(dsn = '1.data/count/PNM/Megaobs_ete24_transect_ln_l93.shp', quiet = TRUE) %>%
  rename(date = transect) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  rename(transect = num_T) %>%
  select(any_of(names(effort)))

eff7 <- st_read(dsn = '1.data/count/PNM/Megaobs_aut24_transect_ln_l93.shp', quiet = TRUE) %>%
  st_transform(st_crs(effort)) %>%
  mutate(LengKm = as.numeric(st_length(.)/1000)) %>%
  rename(transect = num_T) %>%
  select(any_of(names(effort)))

# on combine tout 
pnm_eff <- effort %>%
  st_transform(st_crs(effort)) %>%
  select(any_of(names(eff1))) %>%
  rbind(eff1, eff2, eff3, eff4, eff5, eff6, eff7) %>%
  mutate(year = year(date),
         session = as.factor(camp %>% str_remove("_"))) %>%
  st_transform(st_crs(grid)) 

rm(effort, eff1, eff2, eff3, eff4, eff5, eff6, eff7, obs2)

# Certaines obs n'ont pas de transect !

pnm_obs <- pnm_obs %>%
  group_split(date) %>%
  map_dfr(~ {
    vdate <- unique(.x$date)  # Récupère la date
    veff <- pnm_eff %>% filter(date == vdate)  # Filtre les transects du jour
    
    if (nrow(veff) > 0) {
      .x %>%
        mutate(transect_id = veff$transect[st_nearest_feature(.x, veff)])
    } else {
      .x %>%
        mutate(transect_id = NA)  # Si pas de transect, on met NA
    }
  }) %>%
  select(!transect) %>%
  rename(transect = transect_id) 


pnm_obs <- pnm_obs %>%
  filter(type == "oiseau") %>% 
  rename(nom_fr = espece,
         effectif = nb) %>% 
  mutate(nom_fr = stri_trans_general(str = nom_fr, id = "Latin-ASCII")) %>%  # remove accents
  mutate(nom_fr = tolower(nom_fr)) %>% # in lower case 
  filter(nom_fr %!in% c("grand dauphin", "limicole sp", "dauphin sp", "chevalier sylvain", "chevalier gambette", "poisson lune", "dauphin bleu et blanc","rorqual commun", "meduse oeuf au plat", "dauphin sp.", "exocet","chauve-souris sp.", "chiroptere sp.","vulcain", "morosphynx", "cordulie sp.", "espadon", "autre espece","echasses blanche", "anatide ind", "aigrette garzette","avocette elegante","faucon crecerelle","hirondelle rustique","gravelot a collier interrompu", "tadorne de belon","foulque macroule","alouette sp.", "fringille sp.", "chardonneret elegant", "spatule blanche", "ibis falcinelle", "huitres pie","passereau indetermine","bergeronnette grise", "rouge-queue a front blanc","rouge-queue noir","milan noir","pigeon biset domestique" ,"hirondelle des rochers","martinet noir", "hirondelle indeterminee","chechoier gambette","engoulevent d'europe","pipit rousseline","hirondelle de fenetre","gobemouche noir","pouillot sp.","tourterelle turque","bergeronnette printaniere","serin cini","traquet motteux","pouillot de bonelli","limicoles sp.","busard des roseaux","busard cendre","huppe fasciee","heron cendre","butor etoile","aigrette sp","balbuzard pecheur","epervier d'europe","bergeronette printaniere","fauvette a tete noire","rougegorge familier","becasseau variable","faucon sp","heron pourpe","hirondelle de rocher","linotte melodieuse","bergeronnette des ruisseaux","bergeronnette sp","etourneau sansonnet","tarin des aulnes","heron garde-boeuf", "rouge-gorge familier","grive musicienne","pigeon ramier","alouette des champs","canard colvert","pouillot veloce","canard souchet","heron sp","traquet motteux/oreillard","heron pourpre","traquet sp","hirondelle de rivage","chechoier sylvain","limicoles sp","bondree apivore","faucon kobez","canard sp","sarcelle sp","pinson des arbres","faucon crecerelle/crecerellette","oedicneme criard","caille des bles","oiseau indetermine","pipit indetermine","epervier","hirondelle sp","rougequeue noire","pipit farlouse","passereau sp" ,"martinet sp.","busard sp.","rapace indetermine","cygne noir","grue cendree","faucon","bergeronnette","pipit","flamand rose", "rougequeue","ibis falcinelle ","grive","rougequeue noir","petit limicole indetermine ","bergeronnette printaniere ","pouillot","rouge gorge","rouge queue noir","martinet","passereau ind�termin�" ,"grande aigrette") )  %>%
  mutate(nom_fr = case_when(nom_fr %in% c("alcide", "alcide indetermine") ~ "alcide ind",
                            nom_fr %in% c("cormoran ind.") ~ "cormoran ind",
                            nom_fr %in% c("oceanite ind.") ~ "oceanite ind",
                            nom_fr %in% c("goeland", "goeland indetermine") ~ "goeland ind",
                            nom_fr %in% c("labbe ind.") ~ "labbe ind",
                            nom_fr %in% c("petit puffin") ~ "petit puffin ind",
                            nom_fr %in% c("puffin", "puffin indetermine") ~ "puffin ind",
                            nom_fr %in% c("puffin cendre / de scopoli", "grand puffin ind.") ~ "grand puffin ind",
                            nom_fr %in% c("laride ind.") ~ "laride ind",
                            nom_fr %in% c("mouette", "mouette indeterminee") ~ "mouette ind",
                            nom_fr %in% c("plongeon ind.") ~ "plongeon ind",
                            nom_fr %in% c("sterne", "sterne indeterminee") ~ "sterne ind",
                            nom_fr %in% c("puffin yelkouan de mediterranee") ~ "puffin yelkouan",
                            nom_fr %in% c("cormoran huppe de mediterranee") ~ "cormoran huppe" ,
                            nom_fr %in% c("goeland d audouin") ~ "goeland d'audouin" ,
                            
                            TRUE ~ as.character(nom_fr)))


# Migralion (biotope) --------------------

## Observations ----
load("1.data/count/Biotope/prenup2022_v3.rdata")

prenup22_obs <-  prenup22_obs %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect)

load("1.data/count/Biotope/postnup2022.rdata")

postnup22_obs <-  postnup2022_obs %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect) %>%
  st_transform(st_crs(prenup22_obs)) %>%
  filter(Transect != "Hors-transect")

load("1.data/count/Biotope/prenup2023.rdata")

prenup23_obs <-  prenup2023_obs %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect) %>%
  st_transform(st_crs(prenup22_obs)) %>%
  filter(Transect != "Hors-transect") %>%
  filter(Transect != "OFF effort")

load("1.data/count/Biotope/postnup2023.rdata")

postnup23_obs <-  postnup2023_obs %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect) %>%
  st_transform(st_crs(prenup22_obs)) %>%
  filter(Transect != "Hors-transect") %>%
  filter(Transect != "OFF effort")


postnup24_obs <- read_excel("1.data/count/Biotope/Export_DATA_bateau_v1.7.2_Migralion_2024-08.xlsx", 
                            sheet = "DATA distance sampling") %>%
  filter(is.na(X_tablette) == FALSE) %>%
  st_as_sf(coords= c("X_tablette", "Y_tablette"), crs = 2154, remove = FALSE) %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect) %>%
  st_transform(st_crs(prenup22_obs)) %>%
  filter(Transect != "Hors-transect") %>%
  filter(Transect != "OFF effort")


prenup24_obs <- read_excel("1.data/count/Biotope/Export_DATA_bateau_v1.7.2_Migralion_2024-03.xlsx", 
                           sheet = "DATA distance sampling") %>%
  st_as_sf(coords= c("X_tablette", "Y_tablette"), crs = 2154, remove = FALSE) %>%
  select(Espece, Effectif, Date_UTC, Time_UTC, Gisement, Distance, Transect) %>%
  st_transform(st_crs(prenup22_obs)) %>%
  filter(Transect != "Hors-transect") %>%
  filter(Transect != "OFF effort")



migralion_obs <- rbind(prenup22_obs, postnup22_obs, prenup23_obs, postnup23_obs, prenup24_obs, postnup24_obs) %>%
  rename(nom_fr = Espece,
         date = Date_UTC,
         effectif = Effectif) %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(year = lubridate::year(date)) %>%
  mutate(nom_fr = stri_trans_general(str = nom_fr, id = "Latin-ASCII")) %>%  # remove accents
  mutate(nom_fr = tolower(nom_fr)) %>%
  st_transform(st_crs(grid)) %>%
  filter(nom_fr %!in% c("grand dauphin", "limicole sp", "dauphin sp", "chevalier sylvain", "chevalier gambette", "poisson lune", "dauphin bleu et blanc","rorqual commun", "meduse oeuf au plat", "dauphin sp.", "exocet","chauve-souris sp.", "chiroptere sp.","vulcain", "morosphynx", "cordulie sp.", "espadon", "autre espece","echasses blanche", "anatide ind", "aigrette garzette","avocette elegante","faucon crecerelle","hirondelle rustique","gravelot a collier interrompu", "tadorne de belon","foulque macroule","alouette sp.", "fringille sp.", "chardonneret elegant", "spatule blanche", "ibis falcinelle", "huitres pie","passereau sp.","bergeronnette grise", "rouge-queue a front blanc","rouge-queue noir","milan noir","pigeon biset domestique" ,"hirondelle des rochers","martinet noir", "hirondelle sp.","chechoier gambette","engoulevent d'europe","pipit rousseline","hirondelle de fenetre","gobemouche noir","pouillot sp.","tourterelle turque","bergeronnette printaniere","serin cini","traquet motteux","pouillot de bonelli","limicoles sp.","busard des roseaux","busard cendre","huppe fasciee","heron cendre","butor etoile","aigrette sp","balbuzard pecheur","epervier d'europe","bergeronette printaniere","fauvette a tete noire","rougegorge familier","becasseau variable","faucon sp","heron pourpe","hirondelle de rocher","linotte melodieuse","bergeronnette des ruisseaux","bergeronnette sp","etourneau sansonnet","tarin des aulnes","heron garde-boeuf", "rouge-gorge familier","grive musicienne","pigeon ramier","alouette des champs","canard colvert","pouillot veloce","canard souchet","heron sp","traquet motteux/oreillard","heron pourpre","traquet sp","hirondelle de rivage","chechoier sylvain","limicoles sp","bondree apivore","faucon kobez","canard sp","sarcelle sp","pinson des arbres","faucon crecerelle/crecerellette","oedicneme criard","caille des bles","oiseau indetermine","pipit indetermine","epervier","hirondelle sp","rougequeue noire","pipit farlouse","passereau sp" ,"martinet ind","busard ind","rapace indtermine","cygne noir","grue cendree","faucon","bergeronnette","pipit","flamand rose", "rougequeue","ibis falcinelle ","grive","rougequeue noir","petit limicole indtermine ","bergeronnette printaniere ","pouillot","rouge gorge","rouge queue noir","martinet","passereau indtermin�" ,"grande aigrette") )  %>%
  mutate(nom_fr = case_when(nom_fr %in% c("alcide sp.", "alcide sp") ~ "alcide ind",
                            nom_fr %in% c("type cormoran") ~ "cormoran ind",
                            nom_fr %in% c("labbe parasite/pomarin", "labbe sp") ~ "labbe ind",
                            nom_fr %in% c("puffin yelkouan/baleare", "puffin yelkouan/baleares") ~ "petit puffin ind",
                            nom_fr %in% c("petit larides", "larides sp", "sterne/guifette sp", "petit laridee type pygmee") ~ "laride ind",
                            nom_fr %in% c("mouette sp.", "mouette sp") ~ "mouette ind",
                            nom_fr %in% c("sterne sp") ~ "sterne ind",
                            nom_fr %in% c("puffin sp") ~ "puffin ind",
                            TRUE ~ as.character(nom_fr)))


rm(postnup2022_obs, postnup2023_obs, postnup22_obs, postnup23_obs, postnup24_obs, prenup2023_obs, prenup22_obs, prenup23_obs, prenup24_obs)


## Efforts ----

postnup24_eff <- read_excel("1.data/count/Biotope/Export_DATA_bateau_v1.7.2_Migralion_2024-08.xlsx", 
                            sheet = "DATA Effort") %>%
  rename(Transect = Transect_effort) %>%
  mutate(
    debut_geom = st_as_sfc(WKT_debut_effort, crs = 2154),  # Convertir WKT en sf
    fin_geom   = st_as_sfc(WKT_fin_effort, crs = 2154)
  ) %>%
  select(Transect, Time_UTC_debut_effort, debut_geom, fin_geom) %>%
  st_sf() %>%
  mutate(geometry = purrr::map2(debut_geom, fin_geom, ~ st_linestring(rbind(.x, .y))) %>% st_sfc()) %>%
  st_as_sf() %>%
  mutate(date = as.Date(Time_UTC_debut_effort)) %>%
  st_drop_geometry() %>%
  select(Transect, date, geometry) %>%
  st_as_sf(crs = 2154) %>%
  st_transform((st_crs(prenup22_eff)))


prenup24_eff <- read_excel("1.data/count/Biotope/Export_DATA_bateau_v1.7.2_Migralion_2024-03.xlsx", 
                           sheet = "DATA Effort") %>%
  rename(Transect = Transect_effort) %>%
  mutate(
    debut_geom = st_as_sfc(WKT_debut_effort, crs = 2154),  # Convertir WKT en sf
    fin_geom   = st_as_sfc(WKT_fin_effort, crs = 2154)
  ) %>%
  select(Transect, Time_UTC_debut_effort, debut_geom, fin_geom) %>%
  st_sf() %>%
  mutate(geometry = purrr::map2(debut_geom, fin_geom, ~ st_linestring(rbind(.x, .y))) %>% st_sfc()) %>%
  st_as_sf() %>%
  mutate(date = as.Date(Time_UTC_debut_effort)) %>%
  st_drop_geometry() %>%
  select(Transect, date, geometry) %>%
  st_as_sf(crs = 2154) %>%
  st_transform((st_crs(prenup22_eff)))

prenup22_eff <- prenup22_eff %>%
  rename(date = day) %>%
  mutate(Transect = NA)

postnup22_eff <- postnup2022_eff %>%
  rename(date = day,
         Transect = Transect_effort) %>%
  st_transform((st_crs(prenup22_eff)))

prenup23_eff <- prenup2023_eff %>%
  rename(date = day,
         Transect = Transect_effort) %>%
  st_transform((st_crs(prenup22_eff)))

postnup23_eff <- postnup2023_eff %>%
  rename(date = day,
         Transect = Transect_effort) %>%
  st_transform((st_crs(prenup22_eff)))


migralion_eff <- rbind(prenup22_eff, postnup22_eff, prenup23_eff, postnup23_eff, prenup24_eff, postnup24_eff) %>%
  mutate(year = year(date)) %>%
  st_transform(st_crs(grid))

rm(postnup2022_eff, postnup2023_eff, postnup22_eff, postnup23_eff, postnup24_eff, prenup2023_eff, prenup22_eff, prenup23_eff, prenup24_eff)

# SAMM ----

load("1.data/count/SAMM/birdSAMM.rdata")

# Observtions
samm_obs <- birdsamm %>% 
  select(transect, passage, flight, date, hhmmss, podSize, observer, 
         lat, lon, speed, altitude, famille_fr, nom_fr, geometry) %>% 
  rename(effectif = podSize) %>% 
  mutate(nom_fr = stri_trans_general(str = nom_fr, id = "Latin-ASCII")) %>%  # remove accents
  mutate(nom_fr = tolower(nom_fr)) %>% 
  mutate(year = year(date)) %>%
  filter(year > 2016) %>%
  st_transform(st_crs(grid)) %>% 
  mutate_at("date", as.Date) %>%
  mutate(nom_fr = case_when(nom_fr %in% c("alcide ind.") ~ "alcide ind",
                            nom_fr %in% c("cormoran ind.") ~ "cormoran ind",
                            nom_fr %in% c("grand goeland ind.") ~ "goeland ind",
                            nom_fr %in% c("labbe ind.") ~ "labbe ind",
                            nom_fr %in% c("petit puffin ind.") ~ "petit puffin ind",
                            nom_fr %in% c("puffin cendre / de scopoli") ~ "grand puffin ind",
                            nom_fr %in% c("laride ind.") ~ "laride ind",
                            nom_fr %in% c("mouette ind.") ~ "mouette ind",
                            nom_fr %in% c("plongeon ind.") ~ "plongeon ind",
                            nom_fr %in% c("sterne ind.") ~ "sterne ind",
                            TRUE ~ as.character(nom_fr)))

# Efforts
samm_eff <- effsamm %>% 
  select(nom_suivi, flight, date, hhmmss, DATE_TIME1, seaState, swell, turbidity, skyGlint, 
         glareFrom, glareTo, glareSever, cloudCover, lat, lon, speed, altitude,
         aircraft, seaStIndex, geometry) %>% 
  mutate(year = lubridate::year(date),
         Time = scale(as.Date(DATE_TIME1))) %>%
  filter(year > 2016) %>%
  st_transform(st_crs(grid)) %>% 
  mutate_at("date", as.Date)

save(pelmed_obs, pelmed_eff, 
     samm_obs, samm_eff, 
     pnm_obs, pnm_eff,
     migralion_obs, migralion_eff,
     file = "1.data/all_seabirds_counts.rdata")


# data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
#                   migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
#                   pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
#                   samm = list(obs_data = samm_obs, effort_data = samm_eff))
# 
# obs <- map_dfr(
#   .x = names(data_list),
#   .f = ~ data_list[[.x]]$obs %>%
#     mutate(dataset = .x)  # Ajoute le nom de l'élément comme colonne
# ) %>%
#   #filter(nom_fr %!in% c("grand dauphin", "chevalier sylvain", "chevalier gambette", "poisson lune", "dauphin bleu et blanc","rorqual commun", "meduse oeuf au plat", "dauphin ind", "exocet","chauve-souris ind", "chiroptere ind","vulcain", "moroindynx", "cordulie ind", "einddon", "autre eindce","echasses blanche", "anatide ind", "aigrette garzette","avocette elegante","faucon crecerelle","hirondelle rustique","gravelot a collier interrompu", "tadorne de belon","foulque macroule","alouette ind", "fringille ind", "chardonneret elegant", "indtule blanche", "ibis falcinelle", "huitres pie","passereau ind","bergeronnette grise", "rouge-queue a front blanc","rouge-queue noir","milan noir","pigeon biset domestique" ,"hirondelle des rochers","martinet noir", "hirondelle ind","chechoier gambette","engoulevent d'europe","pipit rousseline","hirondelle de fenetre","gobemouche noir","pouillot ind","tourterelle turque","bergeronnette printaniere","serin cini","traquet motteux","pouillot de bonelli","limicoles ind","busard des roseaux","busard cendre","huppe fasciee","heron cendre","butor etoile","aigrette ind","balbuzard pecheur","epervier d'europe","bergeronette printaniere","fauvette a tete noire","rougegorge familier","becasseau variable","faucon ind","heron pourpe","hirondelle de rocher","linotte melodieuse","bergeronnette des ruisseaux","bergeronnette ind","etourneau sansonnet","tarin des aulnes","heron garde-boeuf", "rouge-gorge familier","grive musicienne","pigeon ramier","alouette des champs","canard colvert","pouillot veloce","canard souchet","heron ind","traquet motteux/oreillard","heron pourpre","traquet ind","hirondelle de rivage","chechoier sylvain","limicole ind","bondree apivore","faucon kobez","canard ind","sarcelle ind","pinson des arbres","faucon crecerelle/crecerellette","oedicneme criard","caille des bles","oiseau indtermine","pipit indtermine","epervier","hirondelle indterminee","rougequeue noire","pipit farlouse","passereau indtermine" ,"martinet ind","busard ind","rapace indtermine","cygne noir","grue cendree","faucon","bergeronnette","pipit","flamand rose", "rougequeue","ibis falcinelle ","grive","rougequeue noir","petit limicole indtermine ","bergeronnette printaniere ","pouillot","rouge gorge","rouge queue noir","martinet","passereau indtermin�" ,"grande aigrette") )     %>%
#   
#   
#   
#   
#   mutate(nom_fr = case_when(nom_fr %in% c("alcide ind") ~ "alcidés indéterminés",
#                             nom_fr %in% c("cormoran ind") ~ "cormorans indéterminés",
#                             nom_fr %in% c("goeland ind") ~ "goélands indéterminés",
#                             nom_fr %in% c("labbe ind") ~ "labbes indéterminés",
#                             nom_fr %in% c("puffin ind") ~ "puffins indéterminés",
#                             nom_fr %in% c("petit puffin ind", "puffin yelkouan/baleare", "puffin yelkouan/baleares", "petit puffin") ~ "petits puffins indéterminés",
#                             nom_fr %in% c("grand puffin ind") ~ "grands puffins indéterminés",
#                             nom_fr %in% c("laride ind") ~ "laridés indéterminés",
#                             nom_fr %in% c("mouette ind") ~ "mouettes indéterminées",
#                             nom_fr %in% c("oceanite ind") ~ "océanites indéterminés",
#                             nom_fr %in% c("plongeon ind") ~ "plongeons indéterminés",
#                             nom_fr %in% c("sterne ind",  "sterne",  "sterne indterminee") ~ "sternes indéterminées",
# 
#                             TRUE ~ as.character(nom_fr))) %>%
# 
#   
#   
#   mutate(nom_fr = case_when(nom_fr %in% c("cormoran huppe") ~ "cormoran huppé",
#                             nom_fr %in% c("goeland d'audouin") ~ "goéland d'Audouin",
#                             nom_fr %in% c("grebe huppe") ~ "grèbe huppé",
#                             nom_fr %in% c("goeland leucophee") ~ "goéland leucophée",
#                             nom_fr %in% c("guillemot de troil") ~ "guillemot de Troïl",
#                             nom_fr %in% c("goeland brun") ~ "goéland brun",
#                             nom_fr %in% c("goeland railleur") ~ "goéland railleur",
#                             nom_fr %in% c("labbe a longue queue") ~ "labbe à longue queue",
#                             nom_fr %in% c("fou de bassan") ~ "fou de Bassan",
#                             TRUE ~ as.character(nom_fr))) %>%
#   
#   group_by(nom_fr) %>%
#   summarize(effectif_total = sum(effectif, na.rm = TRUE), .groups = "drop",
#             nb_obs= n()) %>%
#   
#   mutate(sp = case_when(nom_fr %in% c("macareux moine", "pingouin torda", "guillemot de Troïl", "alcidés indéterminés") ~ "alcidés",
#                         nom_fr %in% c("cormoran huppé", "grand cormoran","cormorans indéterminés") ~ "cormorans",
#                         nom_fr %in% c("goéland leucophée", "goéland brun", "goéland d'Audouin", "goéland railleur", "goélands indéterminés") ~ "goélands",
#                         nom_fr %in% c("labbe parasite", "grand labbe", "labbe pomarin", "labbe à longue queue", "labbes indéterminés") ~ "labbes",
#                         nom_fr %in% c("fou de Bassan") ~ "fous",
#                         nom_fr %in% c("puffin des baleares", "puffin yelkouan", "puffin fuligineux", "puffin de scopoli", "grand puffins indéterminés", "petits puffins indéterminés", "puffins indéterminés") ~ "puffins",
#                         nom_fr %in% c("grèbe huppé") ~ "grèbes",
#                         nom_fr %in% c("guifette moustac", "guifette noire") ~ "guifettes",
#                         nom_fr %in% c("laridés indéterminés") ~ "laridés",
#                         nom_fr %in% c("mouette melanocephale", "mouette rieuse","mouette tridactyle", "mouette pygmee", "mouettes indéterminées") ~ "mouettes",
#                         nom_fr %in% c("oceanite tempete", "océanites indéterminés") ~ "océanites",
#                         nom_fr %in% c("plongeon arctique",  "plongeon imbrin",  "plongeons indéterminés") ~ "plongeons",
#                         nom_fr %in% c("sterne caugek", "sterne pierregarin", "sterne naine", "sterne caspienne", "sternes indéterminées") ~ "sternes",
# 
# 
#                         TRUE ~ as.character(NA)))
# 
# 
# pastel_colors_13 <- c("#8f5d5d","#81b29a","#f2cc8f","#FBB4AE", "#babf95","#3d405b", "#f8a07e", "#fac484" ,"#005b6e", "#b1c7b3", "#FAD02E", "#eb7f86", "#F28D35")
# 
# 
# df <- obs %>%
#   mutate(nb_obs_trunc = pmin(nb_obs, 2000))
# 
# ggplot(df, aes(x = reorder(nom_fr, nb_obs), y = nb_obs_trunc, fill = sp)) +
#   geom_col(width = 0.7) +
#   
#   # Valeurs dans les barres
#   geom_text(
#     aes(label = case_when(
#       nb_obs <= 30 ~ "",
#       nb_obs > 2000 ~ paste0(nb_obs, "+"),
#       TRUE ~ as.character(nb_obs)
#     )),
#     color = "white",
#     hjust = 1.1,
#     size = 3.5
#   ) +
#   scale_fill_manual(values = pastel_colors_13) +
#   # # Nombre d'obsevervations
#   geom_text(
#     aes(label = paste0("(",effectif_total,")")),
#     hjust = -0.2,
#     color = "black",
#     size = 3.5
#   ) +
#   ylim(0,2100) +
#   #scale_y_log10(labels = scales::label_number()) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  coord_flip() +
#   facet_grid(sp ~ ., scales = "free_y", space = "free") +
#   theme_minimal() +
#   labs(x = NULL, y = "Nombre d'observations") +
#   theme(
#     strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),  # titre du groupe
#     strip.placement = "outside",
#     legend.position = "none"
#   )
