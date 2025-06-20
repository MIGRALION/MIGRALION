# RUN RSF ---------

rm(list = ls())  

library(dplyr)
library(sf)
library(purrr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(nimble)
library(MCMCvis)
library(ggspatial)
library(stringr)
library(colorspace)
library(mgcv)
`%!in%` = Negate(`%in%`)

# Path ----
local_path <- "C:/Users/queroue/Documents/Migratlane/Git/Migralion_/"

# Functions ----
source(paste0(local_path, "2.code/functions.R"))

# Data ----
load(paste0(local_path, "1.data/all_seabirds_telemetry.rdata"))
load("1.data/colonies.rdata")

# Covariates ----
load(paste0(local_path, "1.data/covariates.rdata"))
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

# Code ----

## Without individuals effects ----

code_rsf <- nimbleCode({
  
  
  for(i in 1:n.occ.cov){
    
    beta[i] ~ dnorm(0, sd = 3)
    sd_pop[i] ~ dunif(0, 5)

     for(j in 1:nindividual){
       # add individual heterogeneity
       beta_ind[j, i] ~ dnorm(beta[i], sd = sd_pop[i])
     }
  }
  
  # Likelihood
  for(t in 1:npoints){
    
    logit(omega[t]) <- inprod(beta_ind[idind[t], 1:n.occ.cov], XN[site_id[t], 1:n.occ.cov])
    #logit(omega[t]) <- inprod(beta[1:n.occ.cov], XN[site_id[t], 1:n.occ.cov])
    kase[t] ~ dbinom(omega[t], w[t])
    
    
  }
  
})

### Goéland leucophée  ------------------------------------------------------------------------------------------

#### Hors reproduction -----------------------------------------------------


species <- "goeland leucophee"
month_to_keep <- c("09","10","11","12","01","02")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)

results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_HIV_gps_rdm.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)

#### Reproduction (sans colonies) -----------------------------------------------------

species <- "goeland leucophee"
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_ETE_gps_rdm_without_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)


#### Reproduction (avec colonies) -----------------------------------------------------

species <- "goeland leucophee"
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_ETE_gps_rdm_with_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep, covariates)


### Sterne caugek  ------------------------------------------------------------------------------------------

#### Hors reproduction -----------------------------------------------------


species <- "sterne caugek"
month_to_keep <- c("09","10","11","12","01","02")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)

results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/caugek/caugek_HIV_gps_rdm.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)

#### Reproduction (sans colonies) -----------------------------------------------------

species <- "sterne caugek"
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/caugek/caugek_ETE_gps_rdm_without_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)


#### Reproduction (avec colonies) -----------------------------------------------------

species <- "sterne caugek"
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/caugek/caugek_ETE_gps_rdm_with_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep, covariates)

### Puffin yelkouan  ------------------------------------------------------------------------------------------

#### Hors reproduction -----------------------------------------------------


species <- "puffin yelkouan"
month_to_keep <- c("09","10","11","12","01")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)

results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_HIV_gps_rdm.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)

#### Reproduction (sans colonies) -----------------------------------------------------

species <- "puffin yelkouan"
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_gps_rdm_without_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)


#### Reproduction (avec colonies) -----------------------------------------------------

species <- "puffin yelkouan"
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_gps_rdm_with_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep, covariates)


### Puffin yelkouan  ------------------------------------------------------------------------------------------

#### Hors reproduction -----------------------------------------------------


species <- "puffin yelkouan"
month_to_keep <- c("09","10","11","12","01")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = FALSE)

results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_HIV_gps_fix.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)

#### Reproduction (sans colonies) -----------------------------------------------------

species <- "puffin yelkouan"
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = FALSE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_gps_fix_without_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)


#### Reproduction (avec colonies) -----------------------------------------------------

species <- "puffin yelkouan"
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = FALSE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_gps_fix_with_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep, covariates)

### Puffin de Scopoli  ------------------------------------------------------------------------------------------

#### Reproduction (sans colonies) -----------------------------------------------------

species <- "puffin de scopoli"
month_to_keep <- c("04","05","06","07","08","09")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))

data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = TRUE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/scopoli/scopoli_ETE_gps_rdm_without_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep)


#### Reproduction (avec colonies) -----------------------------------------------------

species <- "puffin de scopoli"
month_to_keep <- c("04","05","06","07","08","09")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)

data_gps_aug <- simulate_random_point(data_gps, covariates)
cov_rsf <- get_cov_rsf(data_gps_aug, selected_cov)
data_rsf <- get_data_rsf(data_gps_aug, cov_rsf)
out_rsf <- run_rsf(data_rsf, code_rsf, random_effect = FALSE)


results <- list(out_rsf = out_rsf,
                data_rsf = data_rsf,
                cov_rsf = cov_rsf,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)

save(results, file = paste0(local_path,"3.results/scopoli/scopoli_ETE_gps_fix_with_col.rdata"))
rm(out_rsf, data_rsf, cov_rsf, results, data_gps_aug, data_gps, species, month_to_keep, covariates)
