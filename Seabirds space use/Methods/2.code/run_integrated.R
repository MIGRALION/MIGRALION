# RUN INTEGRATED ---------

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

# PATH
local_path <- "C:/Users/queroue/Documents/Migratlane/Git/Migralion_/"

# FUNCTIONS
source(paste0(local_path, "2.code/functions.R"))

# DATA
load(paste0(local_path, "1.data/all_seabirds_telemetry.rdata"))
load(paste0(local_path, "1.data/all_seabirds_counts.rdata"))

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
                  samm = list(obs_data = samm_obs, effort_data = samm_eff))

rm(pelmed_eff, pelmed_obs, migralion_eff, migralion_obs, pnm_eff, pnm_obs, samm_eff, samm_obs)

# Covariates
load(paste0(local_path, "1.data/covariates.rdata"))
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

# Modèle ---------
code_int <- nimbleCode({
  
  
  for(j in 1: ndatasets) {
    
    p[j] ~ dunif(0, 1)
    
  }
  
  beta0_nmix ~ dnorm(0, sd = 5)
  beta0_gps ~ dnorm(0, sd = 5)
  
  for(i in 1:ncov){
    
    beta[i] ~ dnorm(0, sd = 3)
    #   sd_pop[i] ~ dunif(0, 5)
    #   
    #   for(j in 1:nindividual){
    #     # add individual heterogeneity
    #     beta_ind[j, i] ~ dnorm(beta[i], sd = sd_pop[i])
    #   }
  }
  
  # Likelihood GPS 
  for(t in 1:npoints_gps){
    
    #logit(omega[t]) <- inprod(beta_ind[idind[t], 1:n.occ.cov], XN[site_id[t], 1:n.occ.cov])
    logit(omega[t]) <- beta0_gps + inprod(beta[1:ncov], XN[site_id_gps[t], 1:ncov])
    kase[t] ~ dbinom(omega[t], w[t])
    
    
  }
  
  
  
  # Likelihood Nmix 
  
  # State process
  for(i in 1:nsite_nmix){
    
    log(lambda[site_nmix[i]]) <- beta0_nmix + inprod(beta[1:ncov], XN[site_nmix[i],1:ncov])
    
  }
  
  
  # Observation process
  for (i in 1:npoints_nmix){
    
    N[i] ~ dpois(lambda[site_id_nmix[i]] * surface[i])
    nobs[i] ~ dbin(p[dataset_id[i]], N[i])
    
    
    # Posterior Predictive checks 
    
    # Compute discrepancy for real and simulated data
    # Expected count at this site and this survey
    exp_count[i] <- N[i] * p[dataset_id[i]]
    
    # Discrepancy for the real data
    # (small value added to denominator to avoid potential divide by zero)
    E[i] <- pow((nobs[i] - exp_count[i]), 2) / (exp_count[i] + 0.5)
    
    # Simulate new count from model
    nobs.rep[i] ~ dbin(p[dataset_id[i]], N[i])
    
    # Discrepancy for the simulated data
    E.rep[i] <- pow((nobs.rep[i] - exp_count[i]), 2) / (exp_count[i] + 0.5)
    
  }
  
  
  # chi-squared test statistics
  
  for (nd in 1:ndatasets){
    
    fit[nd] <- sum(E[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1])])
    fit.rep[nd] <- sum(E.rep[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1])])
    
  }
  
})


# Leucophée ----- 

## Hors reproduction ----
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

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/leucophee/leucophee_HIV_int_fix.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)

## Reproduction ----
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

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/leucophee/leucophee_ETE_int_fix_without_col.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)

# Caugek ----- 

## Hors reproduction ----
species <- "sterne caugek"
month_to_keep <- c("09","10","11","12","01","02")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)


data_gps_aug <- simulate_random_point(data_gps, covariates)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/caugek/caugek_HIV_int_fix.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)

## Reproduction ----
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

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/caugek/caugek_ETE_int_fix_without_col.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)

# Yekouan ----- 

## Hors reproduction ----
species <- "puffin yelkouan"
month_to_keep <- c("09","10","11","12","01")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)


data_gps_aug <- simulate_random_point(data_gps, covariates)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_HIV_int_fix.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)

## Reproduction ----
species <- "puffin yelkouan"
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)


data_gps_aug <- simulate_random_point(data_gps, covariates)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_int_fix_without_col.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)


# Puffin de scopoli ----

## Reproduction ----
species <- "puffin de scopoli"
month_to_keep <- c("04","05","06","07","08","09")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


data_gps <- telemetry %>%
  filter_sp_rsf(species) %>%
  filter_month_rsf(month_to_keep) %>%
  st_intersection(covariates) %>%
  filter_per_hour() %>%
  get_new_id_rsf(ind_ID)


data_gps_aug <- simulate_random_point(data_gps, covariates)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)


cov_int <- get_cov_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = selected_cov) 
data_int <- get_data_int(data_nmix_obj = nmix_list, data_gps_obj = data_gps_aug, cov_obj = cov_int) 
out_int <- run_int(data_int, code_int, random_effect = FALSE)


results <- list(out_int = out_int,
                data_int = data_int,
                cov_int = cov_int,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list,
                data_gps = data_gps,
                data_gps_aug = data_gps_aug)


save(results, file = paste0(local_path,"3.results/scopoli/scopoli_ETE_int_fix_without_col.rdata"))

rm(out_int, data_int, cov_int, results, data_gps_aug, data_gps, species_gps, species_nmix, month_to_keep,  nmix_list, nmix_tibble_list)
