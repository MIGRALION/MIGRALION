# RUN NMIX ---------

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
load(paste0(local_path, "1.data/all_seabirds_counts.rdata"))

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
                  samm = list(obs_data = samm_obs, effort_data = samm_eff))

rm(pelmed_eff, pelmed_obs, migralion_eff, migralion_obs, pnm_eff, pnm_obs, samm_eff, samm_obs)

# Covariates
load(paste0(local_path, "1.data/covariates.rdata"))
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

# Code 


code_nmix <- nimbleCode({
  
  # Priors 
  
  
  kappa ~ dunif(min = 0.01, max = 100)
  
  
  for(j in 1: ndatasets) {
    
    p[j] ~ dunif(0, 1)
    
  }
  
  
  for(i in 1:n.occ.cov){
    
    beta[i] ~ dnorm(0, sd = 3)
    
  }
  
  # Likelihood 
  
  # State process
  for(i in 1:nsites_total){
    
    log(lambda[i]) <- inprod(beta[1:n.occ.cov], XN[i,1:n.occ.cov])
    
    
  }
  
  
  # Observation process
  for (i in 1:nsampled_points){
    
    
    prob[i] <- kappa / (kappa + lambda[site_id[i]] * surface[i])
    N[i] ~ dnegbin(prob = prob[i], size = kappa)    # Modèle Negative Binomial
    
    
    nobs[i] ~ dbin(p[dataset_nb[i]], N[i])
    
    
    
    # Posterior Predictive checks 
    
    # Compute discrepancy for real and simulated data
    # Expected count at this site and this survey
    exp_count[i] <- N[i] * p[dataset_nb[i]]
    
    # Discrepancy for the real data
    # (small value added to denominator to avoid potential divide by zero)
    E[i] <- pow((nobs[i] - exp_count[i]), 2) / (exp_count[i] + 0.5)
    
    # Simulate new count from model
    nobs.rep[i] ~ dbin(p[dataset_nb[i]], N[i])
    
    # Discrepancy for the simulated data
    E.rep[i] <- pow((nobs.rep[i] - exp_count[i]), 2) / (exp_count[i] + 0.5)
    
  }
  
  
  # chi-squared test statistics
  
  for (nd in 1:ndatasets){
    
    fit[nd] <- sum(E[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1])])
    fit.rep[nd] <- sum(E.rep[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1])])
    
  }
  
})

## Hors reproduction -----------------------------------------------------

### Macareux moine  ------------------------------------------------------------------------------------------

species <- c("macareux moine")
month_to_keep <- c("12","01","02","03")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/macareux/macareux_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Pingouin torda  ------------------------------------------------------------------------------------------

species <- c("pingouin torda")
month_to_keep <- c("12","01","02","03")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/torda/torda_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Fou de Bassan  ------------------------------------------------------------------------------------------

species <- c("fou de bassan")
month_to_keep <- c("10","11","12","01","02","03","04")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/fou/fou_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Goéland leucophée ------------------------------------------------------------------------------------------

species <- c("goeland leucophee")
month_to_keep <- c("09","10","11","12","01","02")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette tridactyle ------------------------------------------------------------------------------------------

species <- c("mouette tridactyle")
month_to_keep <- c("11","12","01","02","03","04")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/tridactyle/tridactyle_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Labbes ------------------------------------------------------------------------------------------

species <- c("labbe parasite", "grand labbe","labbe ind","labbe pomarin")
month_to_keep <- c("09","10","11","12","01","02","03","04","05")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/labbes/labbes_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette rieuse ------------------------------------------------------------------------------------------

species <- c("mouette rieuse")
month_to_keep <- c("09","10","11","12","01","02","03")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/rieuse/rieuse_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette mélanocéphale ------------------------------------------------------------------------------------------

species <- c("mouette melanocephale")
month_to_keep <- c("09","10","11","12","01","02","03")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/melanocephale/melanocephale_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette pygmée  ------------------------------------------------------------------------------------------

species <- c("mouette pygmee")
month_to_keep <- c("10","11","12","01","02","03","04","05")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/pygmee/pygmee_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Sterne caugek  ------------------------------------------------------------------------------------------

species <- c("sterne caugek")
month_to_keep <- c("09","10","11","12","01","02")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/caugek/caugek_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Puffin des Baléares ------------------------------------------------------------------------------------------

species <- c("puffin des baleares")
month_to_keep <- c("08","09","10","11","12")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/baleares/baleares_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Puffin yelkouan  ------------------------------------------------------------------------------------------

species <- c("puffin yelkouan")
month_to_keep <- c("09","10","11","12","01")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_HIV_nmix.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

## Reproduction (sans colonies) -----------------------------------------------------

### Goéland leucophée  ------------------------------------------------------------------------------------------

species <- c("goeland leucophee")
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette rieuse -------------------------

species <- c("mouette rieuse")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/rieuse/rieuse_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette mélanocéphale -------------------------

species <- c("mouette melanocephale")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/melanocephale/melanocephale_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Puffin des baléares -------------------------

species <- c("puffin des baleares")
month_to_keep <- c("01","02","03","04","05","06")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/baleares/baleares_ETE_nmix_red_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Puffin yelkouan -------------------------

species <- c("puffin yelkouan")
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)


### Puffin de Scopoli -------------------------

species <- c("grand puffin ind", "puffin de scopoli")
month_to_keep <- c("04","05","06","07","08","09")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/scopoli/scopoli_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)


### Océanites  -------------------------

species <- c("oceanite tempete",  "oceanite ind")
month_to_keep <- c("04","05","06","07","08")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/oceanite/oceanite_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)


### Sterne caugek -------------------------

species <- c("sterne caugek")
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/caugek/caugek_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)


### Sterne pierregarin -------------------------

species <- c("sterne pierregarin")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/pierregarin/pierregarin_ETE_nmix_without_col.rdata"))
rm(out_nmix, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

## Reproduction (avec colonies) -----------------------------------------------------

load("1.data/colonies.rdata")

### Goéland leucophée  ------------------------------------------------------------------------------------------

species <- c("goeland leucophee")
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")

covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/leucophee/leucophee_ETE_nmix_with_col.rdata"))
rm(out_nmix, covariates, data_nmix, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mouette rieuse -------------------------

species <- c("mouette rieuse")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/rieuse/rieuse_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)

### Mélanocéphale -------------------------

species <- c("mouette melanocephale")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/melanocephale/melanocephale_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)


### Puffin yelkouan -------------------------

species <- c("puffin yelkouan")
month_to_keep <- c("02","03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/yelkouan/yelkouan_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)


### Puffin de Scopoli -------------------------

species <- c("grand puffin ind", "puffin de scopoli")
month_to_keep <- c("04","05","06","07","08","09")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry", "prox_colonies")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/scopoli/scopoli_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)


### Océanites  -------------------------

species <- c("oceanite tempete",  "oceanite ind")
month_to_keep <- c("04","05","06","07","08")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/oceanite/oceanite_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)


### Sterne caugek -------------------------

species <- c("sterne caugek")
month_to_keep <- c("03","04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)

nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/caugek/caugek_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)


### Sterne pierregarin -------------------------

species <- c("sterne pierregarin")
month_to_keep <- c("04","05","06","07")
selected_cov <- c("log_mean_CHL", "mean_SST", "dist_to_shore", "bathymetry")
load(paste0(local_path, "1.data/covariates.rdata"))
covariates <- get_proximity_colony_score(colonies_obj = colonies_fr, grid_obj = covariates, species_obj = species)


nmix_tibble_list <- map(data_list,
                        ~ {filter_month(.x, month_to_keep) %>%
                            filter_one_species(species)}) %>%
  keep(~ nrow(.x$eff) > 0)


nmix_list <- map(nmix_tibble_list,
                 ~ {get_count_and_effort(.x, covariates)
                 }) %>%
  Filter(f = function(x) nrow(x) != 0)

cov_nmix <- get_occ_cov(nmix_list, selected_cov)
data_nmix <- get_data_nmix(nmix_list, cov_nmix, gam = TRUE)
out_nmix <- run_nmix(data_nmix, code_nmix)

results <- list(out_nmix = out_nmix,
                data_nmix = data_nmix,
                cov_nmix = cov_nmix,
                nmix_list = nmix_list,
                nmix_tibble_list = nmix_tibble_list)

save(results, file = paste0(local_path,"3.results/pierregarin/pierregarin_ETE_nmix_with_col.rdata"))
rm(out_nmix, data_nmix, covariates, cov_nmix, nmix_list, nmix_tibble_list, results)
