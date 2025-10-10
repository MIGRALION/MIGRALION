#################################
#
# Integrated model of migratory intensity across the gulf of Lion 
#
# Integrated GAM with radar echos treated as counts in an IPP 
# and telemetry locations used as presence in an RSF/IPP
#
# Author : Coline Canonne & Sebastien Roques
# Date : 26/06/2025
#
#################################

# ---- Setup ----
rm(list=ls())
library(tidyverse); library(sf);library(nimble);library(MCMCvis);library(mgcv)
library(ggplot2);library(viridis); library(ggspatial); library(ggpattern)
library(rstudioapi);setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# ---  Load shapes of the study area ----
load("../1.data/grid_GOL.Rdata")
load("../1.data/cost_line_shapes.Rdata")
load("../1.data/cities_names.Rdata")

# ---  Load data ----
# For reasons of confidentiality, some of the data is not available here. 
# load("../1.data/private/Echos_grid.Rdata")
# telem_locs <- readRDS('../1.data/private/gps_migration_4.rds')


# ------ Data formatting  --------

# Create study area grid 
poly <- as.data.frame(matrix(c(2, 44, 7, 44, 7, 41, 2, 41), 4,2, byrow = TRUE)) %>% 
  st_as_sf(coords = c("V1", "V2"), crs = 4326) %>%
  summarise((geometry = st_combine(geometry))) %>%
  st_cast("POLYGON")

grid_sf <- st_make_grid(poly, cellsize = 0.05, what="polygons", square=T) %>% 
  st_sf()

# Radar boat data
effort_R <- as.vector(scale(Echos_grid$nb_pix, center = FALSE, scale = TRUE))  # Number of cells sampled by the radar at each scan
mtr_R <- as.vector(scale(Echos_grid$mean_mtr, center = FALSE, scale = TRUE)) # Mean traffic rate of the night measured at the coast (La palissade)

ggplot() +
  geom_sf(data = coast_line) +
  geom_sf(data = Echos_grid,aes(fill= nb_counts/(mean_mtr*nb_pix))) + 
  facet_wrap(~PassageID) +
  coord_sf(xlim = c(2, 7), ylim = c(41, 44)) 


# Telemetry data and generate pseudo-absences
telem_locs_sf <- st_as_sf(telem_locs, coords = c("x", "y"), remove=F, crs = 4326)

background_pts <- grid_sf %>%
  amt::random_points(n = nrow(telem_locs) * 10) %>%
  st_as_sf(coords = c("x_", "y_"), crs = 4326) %>% 
  mutate(x=st_coordinates(.)[,1], y=st_coordinates(.)[,2])

dataRSF <- bind_rows(telem_locs, background_pts) %>%
  mutate(case_ = case_when(case_ == FALSE ~ 0, TRUE ~ 1)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)

dfRSF <- as.data.frame(dataRSF)
w = dfRSF$case_
w[w == 0] <- 1000

ggplot() +
  geom_sf(data = telem_locs_sf, aes(color = as.factor(species)), show.legend = F) +
  geom_sf(data = coast_line, lwd=1, alpha=0.6)+
  coord_sf(xlim = c(2,7),ylim=c(41,44))


## --- Build spline matrices ----

coords_echos_grid <- Echos_grid %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  as_tibble()

coords_telemetry <- dataRSF %>% 
  st_coordinates() %>% 
  as_tibble() 

X <- c(coords_echos_grid$X, coords_telemetry$X)
Y <- c(coords_echos_grid$Y, coords_telemetry$Y)
coords_tot <- data.frame(X,Y)

# Create X-Y splines for the RSF-GAM
k = 29 # Number of nodes for the spline
spline_matrix <- smoothCon(s(X, Y, k = k + 1), data = coords_tot, absorb.cons = T)


## --- Bundle data ----

cst <- list(ntot = nrow(dfRSF), w = w, k = k,
            nsites = nrow(Echos_grid), effort = effort_R, mtr = mtr_R)

IM_data <-  list(kase = dfRSF$case_,
                 X = cbind(rep(1, length(X)), spline_matrix[[1]]$X),
                 nechosBoat =  Echos_grid$nb_counts)

initial_values <-  list(beta = rep(0, k), intIPP = 0, intRSF = 0)

# ---------------------------------------------------------------------------- |
# ---- NIMBLE MODEL ----

integrated <- nimbleCode({
  # PRIORS --
  intRSF ~ dnorm(0, sd = 50)
  intIPP ~ dnorm(0, sd = 50)
  for (i in 1:k){beta[i] ~ dnorm(0, sd = 2)}
  
  # LIKELIHOOD --
  # Telemetry RSF
  for (i in 1:ntot){
    logit(omega[i]) <- inprod(c(intRSF, beta[]), X[nsites + i, ])
    kase[i] ~ dbinom(omega[i], w[i])
  }
  # Radar IPP 
  for(j in 1:nsites){
    nechosBoat[j] ~ dpois(lambda[j] * effort[j] * mtr[j])
    log(lambda[j]) <- inprod(c(intIPP,beta[]), X[j,])
  }
})

## --- Run MCMC model ----

Rmod <- nimbleModel(code = integrated, constants = cst, data = IM_data, inits = initial_values)
Rmod$initializeInfo()
Rmod$calculate()
conf <- configureMCMC(Rmod)
conf$removeSamplers('beta', print = FALSE)
conf$addSampler(target = paste0("beta[1:",k , "]"), type = 'AF_slice')
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmod)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 2000, nchains = 2, samplesAsCodaMCMC = TRUE, summary = TRUE)
save.image(file = "results.Rdata")

# ---------------------------------------------------------------------------- |

# --- Results ----

MCMCsummary(samples)
MCMCtrace(samples$samples, pdf=FALSE)
MCMCplot(object = samples$samples)

### Make predictions from model results with coefficient of variation (CV)

GOL_grid <- st_intersection(grid_sf, st_union(shape_GOL))

coord_GOL_pred <- GOL_grid %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  as_tibble() 
Xp_xy <- PredictMat(sm_xy[[1]], coord_GOL_pred)

GOL_predict <- matrix(NA, dim(coord_GOL_pred)[1], nb_it)
for (it in 1:nb_it){
  GOL_predict[,it] <- exp(Xp_xy %*% 
                            samples$samples$chain1[it, c(paste0("beta[", 1:k, "]"))])
}
GOL_grid$predicted_values  <- apply(GOL_predict, 1, function(x) median(x))
GOL_grid$predicted_CV  <- apply(GOL_predict, 1, function(x) sd(x) / mean(x) * 100)

# Scale relative intensity between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
GOL_grid$predicted_values <- range01(GOL_grid$predicted_values)



