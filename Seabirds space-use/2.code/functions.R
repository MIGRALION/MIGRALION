#                              -  Functions -

#------------------------------ Covariates -----------------------
get_occ_cov <- function(nmix_tibble, selected_cov_obj){
  
  
  cov_no_scale <- map(nmix_tibble, function(x) select(x, all_of(c(selected_cov_obj, "grid_id")))) %>% 
    bind_rows() %>% 
    unique() %>% 
    arrange(by = grid_id) %>% 
    select(all_of(selected_cov_obj)) 
  
  cov_scaled <- cov_no_scale %>% 
    rename_with(~ paste0(., "_scaled")) %>%
    mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
  
  all_cov <- cbind(cov_no_scale, cov_scaled)
  
  return(list(cov_scaled = cov_scaled, cov_no_scale = cov_no_scale))
}

get_cov_const <- function(all_cov, scaled = T) {
  
  if(scaled == T) {
    cov_const <- all_cov$cov_scaled %>%
      select(ends_with("_scaled")) %>%
      mutate(intersect = 1) %>% 
      relocate(intersect) %>% 
      as.matrix()
  }
  
  if(scaled == F) {
    cov_const <- all_cov$cov_no_scale %>%
      select(!ends_with("_scaled")) %>%
      mutate(intersect = 1) %>% 
      relocate(intersect) %>% 
      as.matrix()
  }
  
  return(cov_const)
}

get_proximity_colony_score <- function(colonies_obj, grid_obj, species_obj){
  
  colonies_sp <- colonies_obj %>%
    filter(Espece %in% species_obj)
  
  colonies_size <- colonies_sp %>% pull(EFF) %>% as.vector()
  
  dist <- st_distance(grid_obj, colonies_sp) / 1000 
  
  # Convertir en data.frame et retirer les unités
  dist_df <- dist %>%
    as.data.frame() %>%
    units::drop_units()
  
  if(length(colonies_size) > 1){
    # Appliquer la division colonne par colonne
    resultat <- map2_dfc(dist_df, colonies_size, ~ .y / (.x + 1)) %>%
      mutate(colonies_prox_score = rowSums(., na.rm = TRUE)) %>%
      pull(colonies_prox_score)}
  
  if(length(colonies_size) == 1){
    
    a <- dist_df/ colonies_size
    
    resultat <- dist_df %>%
      mutate(colonies_prox_score = a$.) %>%
      pull(colonies_prox_score)
  }
  
  grid_obj <- grid_obj %>%
    mutate(prox_colonies = resultat) 
  
  return(grid_obj)
  
  
}

get_grid_scaled <- function(grid_obj, cov_obj){
  
  selected_cov <- colnames(cov_obj$cov_no_scale)
  
  stats <- cov_obj$cov_no_scale %>%
    select(all_of(selected_cov)) %>%
    summarise(across(everything(), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))
  
  means <- stats %>% select(ends_with("_mean")) %>% as.list()
  sds <- stats %>% select(ends_with("_sd")) %>% as.list()
  
  grid_scaled <- grid_obj %>%
    select(all_of(selected_cov)) %>%
    mutate(across(selected_cov, 
                  ~ (. - means[[paste0(cur_column(), "_mean")]]) / 
                    sds[[paste0(cur_column(), "_sd")]]))  %>%
    rename_with(~ paste0(., "_scaled")) %>%
    rename(geometry = geometry_scaled)
  
  
  return(grid_scaled)
  
}

get_cov_rsf <- function(data_obj, selected_cov_obj){
  
  cov_no_scale <- data_obj %>%
    distinct(grid_id, .keep_all = TRUE) %>% 
    arrange(by = grid_id) %>%
    select(all_of(selected_cov_obj)) %>% 
    st_drop_geometry() 
  
  
  cov_scaled <- cov_no_scale %>% 
    rename_with(~ paste0(., "_scaled")) %>%
    mutate(across(everything(), ~ (. - mean(.)) / sd(.)))
  
  all_cov <- cbind(cov_no_scale, cov_scaled)
  
  return(list(cov_scaled = cov_scaled, cov_no_scale = cov_no_scale))
  
}

get_cov_int <-  function(data_nmix_obj, data_gps_obj, cov_obj){
  # Gestion des covariables poru les GAM  
  cov_no_scale_nmix <- map(data_nmix_obj, function(x) select(x, all_of(c(cov_obj, "grid_id")))) %>% 
    bind_rows() %>% 
    unique() %>% 
    arrange(by = grid_id) %>% 
    select(grid_id, all_of(cov_obj)) 
  
  cov_no_scale_rsf <- data_gps_obj %>%
    distinct(grid_id, .keep_all = TRUE) %>% 
    arrange(by = grid_id) %>%
    select(grid_id, all_of(cov_obj)) %>% 
    st_drop_geometry() 
  
  cov_no_scale <- rbind(cov_no_scale_nmix, cov_no_scale_rsf) %>%
    distinct(grid_id, .keep_all = TRUE)
  
  cov_scaled <- cov_no_scale %>% 
    mutate(across(any_of(cov_obj), ~ (. - mean(.)) / sd(.))) %>%
    rename_with(~ paste0(., "_scaled"), any_of(cov_obj))
  
  all_cov <- list(cov_scaled = cov_scaled, cov_no_scale = cov_no_scale)
}



#------------------------------ Nmixture -------------------------------------

## Prepare data -----

filter_month <- function(data, months){
  
  obs <- data$obs %>% filter(format(date, "%m") %in% months)
  eff <- data$eff %>% filter(format(date, "%m") %in% months)
  
  return(list(obs = obs, eff = eff))
}


filter_one_species <- function(data, species){
  #obs <- data$obs %>% filter(species_name == species)
  obs <- data$obs %>% filter(nom_fr %in% species)
  
  return(list(obs = obs, eff = data$eff))
  
}

get_count_and_effort <- function(data, grid){
  
  #' Get the count data and compute transect length for all sampled cells
  #' @data: a list containing two sf data frame, obs with the observed counts and
  #' eff which contains the transects. Both data frame should have a column
  #' "session" indicating the corresponding session for each transect or observation
  #' @grid: a sf data frame containing the grid, with a column "id_cells"
  #' @Output: same format as grid with additionnal columns (transect_length and 
  #' effectif for each session)
  
  intersect_grid_obs <- st_intersection(data$obs, grid) 
  intersect_grid_eff <- st_intersection(data$eff, grid) 
  
  id_sampled_cells <- intersect_grid_eff %>% pull(grid_id) %>% unique() %>% sort()
  
  sampled_gridcells <- grid %>% 
    filter(grid_id %in% id_sampled_cells) %>% 
    as_tibble() %>% 
    st_drop_geometry()
  
  # # Get session names (only when the species is observed)

  session_names <- data$eff %>% pull(date) %>% unique()
  
  
  for (k in seq_along(session_names)){
    sessionK <- session_names[k]
    # transect length
    intersect_eff_sessionK <- intersect_grid_eff %>% 
      #filter(session == sessionK) 
      filter(date == sessionK) 
    
    df_transect_length <- tibble(grid_id = intersect_eff_sessionK$grid_id, 
                                 len = st_length(intersect_eff_sessionK)) %>% 
      group_by(grid_id) %>% 
      summarise(transectLength = sum(len))
    
    # Retrieve effectif 
    effectif_df <- intersect_grid_obs %>%
      as_tibble() %>% 
      filter(date == sessionK) %>% 
      select(grid_id, effectif) %>% 
      group_by(grid_id) %>% 
      summarise(effectif = sum(effectif)) %>% 
      arrange(grid_id)
    
    # Add in the grid
    sampled_gridcells <- sampled_gridcells %>% 
      left_join(df_transect_length, by = join_by(grid_id)) %>% 
      left_join(effectif_df, by = join_by(grid_id)) %>% 
      mutate(effectif = ifelse(is.na(effectif) & !is.na(transectLength), 0, effectif)) %>%
      rename_with(~ paste0(.x, "_", as.character(sessionK), recycle0 = TRUE), 
                  .cols = c(effectif, transectLength))
  }
  return(sampled_gridcells)
}



## Run model ----

get_data_nmix <- function(data_list, cov_obj, gam = FALSE) {
  
  result_list <- map(data_list, ~ {
    .x %>%
      select(starts_with("transectLength"), starts_with("effectif"), grid_id) %>%
      relocate(grid_id) %>%
      pivot_longer(-c(grid_id), names_to = c(".value", "session"), names_pattern = "(.+)_(.+)") %>%
      # Keep only cells that have been sampled
      filter(!is.na(transectLength)) %>% 
      # Scale transect length (mean = 0 and sd = 1)
      mutate(transectLength_std = (transectLength - mean(transectLength))/sd(transectLength)) %>%
      mutate(session_nb = as.numeric(as.factor(session)))
  }) %>% 
    bind_rows(.id = "dataset_nb") %>% 
    mutate(dataset_nb = as.numeric(as.factor(dataset_nb))) %>%
    arrange(dataset_nb)
  
  
  # Get the observation data
  data <- list(nobs = result_list$effectif)
  
  ndatasets <- max(result_list$dataset_nb)
  
  
  # Get the constants
  
  XN <- get_cov_const(cov_obj)
  
  if(gam == TRUE){ 
    
    data_spline <- as.data.frame(XN)
    
    XN <- rep(1, n_distinct(result_list$grid_id)) %>%
      as.data.frame()
    
    for(i in 2:ncol(data_spline)){
      
      a <- data_spline[,i] %>%
        as.data.frame()
      
      colnames(a)[1] <- "cov"
      
      cov_spline <- smoothCon(s(cov, k = 3), data = a, absorb.cons = T)
      
      XN <- cbind(XN, cov_spline[[1]]$X)
      
      rm(cov_spline, a)
      
    }
  }
  
  
  constants <- list(XN = XN,
                    surface = result_list$transectLength,
                    site_id = get_new_id(result_list, id_column = "grid_id"),
                    n.occ.cov = ncol(XN),
                    dataset_nb = result_list$dataset_nb,
                    nsites_total = n_distinct(result_list$grid_id),
                    nsampled_points = nrow(result_list),
                    ndatasets = ndatasets)
  
  # Set initial values
  # N0 <- result_list %>% 
  #   group_by(grid_id) %>% 
  #   summarise(N0 = max(effectif)) %>% 
  #   arrange(by = "grid_id") %>% 
  #   pull(N0)
  # 
  initial.values <- list(beta = rnorm(constants$n.occ.cov, 0, 1),
                         p = rep(0.5, constants$ndatasets))#, 
  #N = N0)
  
  return(list(data = data, constants = constants, inits = initial.values))
  
  
}


get_new_id <- function(df_w_id, id_column){
  # Pull the id and create a data frame with the current and the new id
  sampled_sites_id <- df_w_id %>% pull(id_column) %>% unique() %>% sort()
  id_correspondance <- tibble(old_id = sampled_sites_id, new_id = 1:length(sampled_sites_id))
  
  # Put the new id in the good row by correspondence with old id and pull them
  new_sites_id <- df_w_id %>% select(all_of(id_column)) %>% 
    left_join(id_correspondance, by = join_by({{ id_column }} == "old_id")) %>%
    pull(new_id) 
  return(new_sites_id)
}


run_nmix <- function(data_obj, code){
  
  data_obj$inits$N <- data_nmix$data$nobs +1
  data_obj$inits$p <- runif(data_nmix$constants$ndatasets, 0, 1)
  data_obj$inits$kappa <- 1
  
  # Parameters 
  parameters.to.save <- c("beta", "fit", "fit.rep", "p","kappa")
  
  # Add info consts
  data_obj$constants$npoints_dataset <- cumsum(c(0, table(data_obj$constants$dataset_nb) %>% unname()))
  # Add info inits
  #data_obj$inits$kappa <- 0.5
  
  # Define model
  Nmix.model <- nimbleModel(code = code, 
                            constants = data_obj$constants,
                            data = data_obj$data, 
                            inits = data_obj$inits)
  
  Nmix.model$initializeInfo()
  Nmix.model$calculate()
  
  
  # Configure model
  confo <- configureMCMC(Nmix.model, monitors = parameters.to.save, thin = 10)
  
  
  confo$removeSampler(paste0("beta[", 1:data_obj$constants$n.occ.cov, "]"))
  confo$addSampler(paste0("beta[", 1:data_obj$constants$n.occ.cov, "]"), 
                   type = "AF_slice")
  
  ## Build and compile MCMC
  Rmcmco <- buildMCMC(confo)
  Cmodelo <- compileNimble(Nmix.model)
  Cmcmco <- compileNimble(Rmcmco, project = Cmodelo)
  
  # Run
  mcmc.output <- runMCMC(Cmcmco, 
                         niter = 100000, 
                         nburnin = 50000, 
                         nchains = 3, 
                         thin = 10, 
                         samplesAsCodaMCMC = TRUE)
  
  return(mcmc.output)
}

#---------------------------------- RSF ---------------------------------------


## Prepare data ---------------------------------------------------------------

filter_per_hour <- function(data_obj){
  
  data_obj %>%
    mutate(hourly_time = floor_date(time, "hour")) %>% 
    group_by(ind_ID, hourly_time) %>%  
    arrange(time) %>%
    slice(1) %>%
    ungroup()
}


filter_month_rsf <- function(data_obj, months){
  data_obj %>%
    filter(format(time, "%m") %in% months)
}


filter_sp_rsf <- function(data_obj, species_name){
  
  data_obj %>%
    filter(species %in% species_name)
  
}



get_new_id_rsf <- function(data_obj, id_column){
  
  data_obj %>%
    mutate(id = as.numeric(factor({{id_column}}, levels = unique({{id_column}}), labels = seq_along(unique({{id_column}})))))
  
}

simulate_random_point <- function(data_obj, grid_obj, n_aug = 10, w = 1000, buffer_size = 5000) {
  
  
  all_data <- data_obj %>%
    select(all_of(colnames(grid_obj %>% st_drop_geometry)), id) %>%
    mutate(weight = 1,
           presence = 1)
  
  for(i in unique(data_obj$id)) {
    
    # area <- data_gps %>%
    #   filter(id == i) %>%  # Filtrer les données
    #   st_union() %>% 
    #   st_convex_hull() %>% # Fusionner les points en une seule géométrie
    #   st_buffer(buffer_size) %>%
    #   st_intersection(grid_obj)
    
    area <- grid_obj %>%
      st_union()
    
    random_points <- st_sample(area, size = nrow(data_obj %>% filter(id == i)) * n_aug, type = "random") %>%
      st_as_sf() %>%
      st_join(grid_obj) %>%
      mutate(id = rep(i, nrow(data_obj %>% filter(id == i)) * n_aug)) %>%
      rename(geometry = x) %>%
      mutate(weight = w,
             presence = 0) %>%
      select(all_of(colnames(all_data)))
    
    
    all_data <- bind_rows(all_data, random_points)
    
    rm(random_points, area)
  }
  
  return(all_data)
  
}

## Run model ----

get_data_rsf <- function(data_obj, cov_obj){
  
  XN <- get_cov_const(cov_obj)
  
  data_spline <- as.data.frame(XN)
  
  XN <- rep(1, nrow(XN)) %>%
    as.data.frame()
  
  for(i in 2:ncol(data_spline)){
    
    a <- data_spline[,i] %>%
      as.data.frame()
    
    colnames(a)[1] <- "cov"
    
    cov_spline <- smoothCon(s(cov, k = 3), data = a, absorb.cons = T)
    
    XN <- cbind(XN, cov_spline[[1]]$X)
    
    rm(cov_spline, a)
    
  }
  
  # Get the observation data
  data <- list(kase = data_obj$presence)
  
  
  constants <-  list(npoints = nrow(data_obj),
                     idind = as.numeric(data_obj$id),
                     XN = XN,
                     site_id = get_new_id(data_obj, id_column = "grid_id"),
                     nindividual = n_distinct(data_obj$id),
                     n.occ.cov = ncol(XN),
                     w = data_obj$weight)
  
  #inits
  inits <-  list(beta = rep(0, constants$n.occ.cov),
                 sd_pop = runif(constants$n.occ.cov, 0, 5),
                 beta_ind = matrix(0, ncol = constants$n.occ.cov , nrow = constants$nindividual))
  
  return(list(data = data, constants = constants, inits = inits))
  
}


run_rsf <- function(data_obj, code, random_effect = FALSE) {
  
  
  # Nimble pre run
  rsf.model <- nimbleModel(code = code, 
                           constants = data_obj$constants, 
                           data = data_obj$data, 
                           inits = data_obj$inits)
  
  
  # rsf.model$initializeInfo()
  # rsf.model$calculate()
  
  if(random_effect == FALSE){
    
    parameters <- c("beta")
    
  }
  
  
  if(random_effect == TRUE){
    
    parameters <- c("beta", "beta_ind", "sd_pop")
    
  }
  
  
  # configure model
  model_configuration <- configureMCMC(rsf.model, monitors = parameters, thin = 1)
  
  
  if(random_effect == FALSE){
    
      
      model_configuration$removeSampler(c(paste0("beta[",1:data_obj$constants$n.occ.cov, "]")))
      model_configuration$addSampler(c(paste0("beta[",1:data_obj$constants$n.occ.cov, "]")),
                                     type = "AF_slice")
    
  }
  
  
  
  if(random_effect == TRUE){
    
    model_configuration$removeSampler(c(paste0("beta[",1:data_obj$constants$n.occ.cov, "]")))
    model_configuration$addSampler(c(paste0("beta[",1:data_obj$constants$n.occ.cov, "]")),
                                   type = "AF_slice")
    
    for(i in 1:data_obj$constants$nindividual){
      
      model_configuration$removeSampler(c(paste0("beta_ind[",i,", ", 1:data_obj$constants$n.occ.cov, "]")))
      model_configuration$addSampler(c(paste0("beta_ind[",i,", ", 1:data_obj$constants$n.occ.cov, "]")),
                                     type = "AF_slice")
      
    }
    
  }
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(model_configuration)
  Cmodel <- compileNimble(rsf.model)
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
  
  # Run
  samplesRSF <- runMCMC(Cmcmc, 
                        niter = 5000, 
                        nburnin = 1000, 
                        thin = 1, 
                        nchains = 3)
  
  return(samplesRSF)
}

#------------------------------ Integrated model ----------

get_id_int <- function(id_column, data_obj){
  # Pull the id and create a data frame with the current and the new id
  new_id <- id_column %>% mutate(new_id = 1:nrow(.))
  data_obj <- data_obj %>%
    left_join(new_id, by = "grid_id")
  return(data_obj$new_id)
}

get_data_int <- function(data_nmix_obj, data_gps_obj, cov_obj) {
  
  XN <- get_cov_const(cov_obj)
  
  data_spline <- as.data.frame(XN) 
  
  XN <- cov_obj$cov_scaled %>%
    select(grid_id)
  
  
  for(i in 2:ncol(data_spline)){
    
    a <- data_spline[,i] %>%
      as.data.frame()
    
    colnames(a)[1] <- "cov"
    
    cov_spline <- smoothCon(s(cov, k = 3), data = a, absorb.cons = T)
    
    XN <- cbind(XN, cov_spline[[1]]$X)
    
    rm(cov_spline, a)
    
  }
  
  
  # Data Nmixture
  result_list <- map(data_nmix_obj, ~ {
    .x %>%
      select(starts_with("transectLength"), starts_with("effectif"), grid_id) %>%
      relocate(grid_id) %>%
      pivot_longer(-c(grid_id), names_to = c(".value", "session"), names_pattern = "(.+)_(.+)") %>%
      # Keep only cells that have been sampled
      filter(!is.na(transectLength)) %>% 
      # Scale transect length (mean = 0 and sd = 1)
      mutate(transectLength_std = (transectLength - mean(transectLength))/sd(transectLength)) %>%
      mutate(session_nb = as.numeric(as.factor(session)))
  }) %>% 
    bind_rows(.id = "dataset_nb") %>% 
    mutate(dataset_nb = as.numeric(as.factor(dataset_nb))) %>%
    arrange(dataset_nb)
  
  
  # Get the observation data
  data <- list(nobs = result_list$effectif,
               kase = data_gps_obj$presence)
  
  
  constants <-  list(XN = XN[,2:ncol(XN)],
                     ncov = ncol(XN) - 1,
                     site_id_gps = get_id_int(XN %>% select(grid_id), data_gps_obj),
                     nind_gps = n_distinct(data_gps_aug$id),
                     npoints_gps = nrow(data_gps_obj),
                     id_ind_gps = as.numeric(data_gps_aug$id),
                     w = data_gps_obj$weight,
                     surface = result_list$transectLength,
                     dataset_id = result_list$dataset_nb,
                     npoints_nmix = nrow(result_list),
                     site_id_nmix = get_id_int(XN %>% select(grid_id), result_list),
                     site_nmix = unique(get_id_int(XN %>% select(grid_id), result_list)),
                     nsite_nmix = length(unique(get_id_int(XN %>% select(grid_id), result_list))),
                     ndatasets = max(result_list$dataset_nb),
                     npoints_dataset = cumsum(c(0, table(result_list$dataset_nb) %>% unname())))
  
  #inits
  
  inits <-  list(beta = rep(0, constants$ncov),
                 beta0_nmix = 0,
                 beta0_gps = 0,
                 #sd_pop = runif(constants$ncov, 0, 5),
                 #beta_ind = matrix(0, ncol = constants$ncov , nrow = constants$nind_gps),
                 p = rep(0.5, constants$ndatasets),
                 N = data$nobs + 1)
  
  return(list(data = data, constants = constants, inits = inits))
  
}

run_int <- function (data_obj, code_obj, random_effect = FALSE){ 
  
  # Parameters 
  if(random_effect == FALSE){
  
    parameters.to.save <- c("beta", "fit", "fit.rep", "p", "beta0_nmix", "beta0_gps")
  
  }
  
  if(random_effect == TRUE){
    
    parameters.to.save  <- c("beta", "beta_ind", "sd_pop", "fit", "fit.rep", "p", "beta0_nmix", "beta0_gps")
    
  }
  
  # Nimble pre run
  int.model <- nimbleModel(code = code_int, 
                           constants = data_int$constants, 
                           data = data_int$data, 
                           inits = data_int$inits)
  
  # Configure model
  confo <- configureMCMC(int.model, monitors = parameters.to.save, thin = 10)
  
  if(random_effect == FALSE){
    
    
    confo$removeSampler(c(paste0("beta[", 1:data_int$constants$ncov, "]"),"beta0_nmix", "beta0_gps"))
    confo$addSampler(c(paste0("beta[", 1:data_int$constants$ncov, "]"),"beta0_nmix", "beta0_gps"),
                     type = "AF_slice")
    
  }
  
  
  
  if(random_effect == TRUE){
    
    confo$removeSampler(c(paste0("beta[", 1:data_int$constants$ncov, "]"),"beta0_nmix", "beta0_gps"))
    confo$addSampler(c(paste0("beta[", 1:data_int$constants$ncov, "]"),"beta0_nmix", "beta0_gps"),
                     type = "AF_slice")
    
    for(i in 1:data_obj$constants$nindividual){
      
      confo$removeSampler(c(paste0("beta_ind[",i,", ", 1:data_int$constants$ncov, "]")))
      confo$addSampler(c(paste0("beta_ind[",i,", ", 1:data_int$constants$ncov, "]")),
                                     type = "AF_slice")
      
    }
    
  }
  
  

  
  
  ## Build and compile MCMC
  Rmcmco <- buildMCMC(confo)
  Cmodelo <- compileNimble(int.model)
  Cmcmco <- compileNimble(Rmcmco, project = Cmodelo)
  
  # Run
  mcmc.output <- runMCMC(Cmcmco, 
                         niter = 5000, 
                         nburnin = 1000, 
                         nchains = 3, 
                         thin = 5, 
                         samplesAsCodaMCMC = TRUE) 
  
  return(mcmc.output)
  
}

#------------------------------ Calculate results ----------

func_cv <- function(x) { sd(x) / mean(x) * 100 }

make_spatial_prediction_gam <- function(cov_obj, grid_obj, out_obj){ 
  
  beta <- out_obj$chain1 %>% as.data.frame() %>%
    select(starts_with("beta")) %>%
    rbind(out_obj$chain2 %>%
            as.data.frame() %>% select(starts_with("beta")))
  
  # Get the constants
  XN <- get_cov_const(cov_obj)
  
  data_spline <- as.data.frame(XN)
  
  grid_scaled <- get_grid_scaled(grid_obj, cov_obj) 
  
  selected_cov <- colnames(cov_obj$cov_no_scale)
  
  XN_pred <- rep(1, nrow(grid_scaled)) %>%
    as.data.frame()
  
  for(i in 2:ncol(data_spline)){
    
    a <- data_spline[,i] %>%
      as.data.frame()
    
    colnames(a)[1] <- "cov"
    
    cov_spline <- smoothCon(s(cov, k = 3), data = a, absorb.cons = T)
    
    grid_pred <- grid_scaled %>% select(colnames(data_spline[i])) %>%
      st_drop_geometry() %>%
      rename(cov = colnames(data_spline[i]))
    
    pred <- PredictMat(object = cov_spline[[1]], 
                       data = grid_pred)
    
    XN_pred <- cbind(XN_pred, pred)
    
    rm(cov_spline, a, pred, grid_pred)
    
    
  }
  
  XN_pred <- as.matrix(XN_pred)
  
  a <- matrix(NA, nrow = nrow(XN_pred), ncol = nrow(beta))
  
  for(i in 1:nrow(beta)){
    a[,i] <- exp(as.vector(XN_pred[,2:ncol(XN_pred)] %*% as.vector(as.matrix(beta[i, 2:ncol(XN_pred)])) ))
  }
  
  
  grid_pred <- grid_scaled %>% 
    mutate(mean = apply(a, 1, mean)) %>%
    mutate(cv = apply(a, 1, func_cv))
  
  return(grid_pred)
  
}

make_spatial_prediction_int <- function(cov_obj, grid_obj, out_obj){
  
  
  
  beta <- out_obj$chain1 %>% as.data.frame() %>%
    select(starts_with("beta[")) %>%
    rbind(out_obj$chain2 %>%
            as.data.frame() %>% select(starts_with("beta[")))
  
  # Get the constants
  XN <- get_cov_const(cov_obj)
  
  data_spline <- as.data.frame(XN)
  
  grid_scaled <- get_grid_scaled(grid_obj, cov_obj) 
  
  selected_cov <- colnames(cov_obj$cov_no_scale)
  
  XN_pred <- rep(1, nrow(grid_scaled)) %>%
    as.data.frame()
  
  for(i in 2:ncol(data_spline)){
    
    a <- data_spline[,i] %>%
      as.data.frame()
    
    colnames(a)[1] <- "cov"
    
    cov_spline <- smoothCon(s(cov, k = 3), data = a, absorb.cons = T)
    
    grid_pred <- grid_scaled %>% select(colnames(data_spline[i])) %>%
      st_drop_geometry() %>%
      rename(cov = colnames(data_spline[i]))
    
    pred <- PredictMat(object = cov_spline[[1]], 
                       data = grid_pred)
    
    XN_pred <- cbind(XN_pred, pred)
    
    rm(cov_spline, a, pred, grid_pred)
    
    
  }
  
  XN_pred <- as.matrix(XN_pred)[,2: ncol(XN_pred)] 
  
  a <- matrix(NA, nrow = nrow(XN_pred), ncol = nrow(beta))
  
  for(i in 1:nrow(beta)){
    a[,i] <- exp(as.vector(XN_pred[,1:ncol(XN_pred)] %*% as.vector(as.matrix(beta[i, 1:ncol(XN_pred)])) ))
  }
  
  grid_pred <- grid_scaled %>% 
    mutate(mean = apply(a, 1, mean)) %>%
    mutate(cv = apply(a, 1, func_cv))
  
  return(grid_pred)
  
}

make_cov_prediction_gam <- function(cov_obj, const_obj, out_obj, gps = FALSE) {
  
  
  beta <- out_obj$chain1 %>% as.data.frame() %>%
    select(starts_with("beta")) %>%
    rbind(out_obj$chain2 %>%
            as.data.frame() %>% select(starts_with("beta")))
  
  if(gps == TRUE){
    
    df <- data.frame(w = const_obj$w, site_id = const_obj$site_id) %>%
      filter(w == 1) %>%
      distinct(site_id) %>%
      pull(site_id)
    
    selected_cov <- colnames(cov_obj$cov_no_scale) 
    
    names(const_obj$XN) <- make.unique(names(const_obj$XN))
    
    XN <- const_obj$XN %>%
      slice(df) %>% as.matrix()
    
    cov_effect <-  cov_obj$cov_no_scale %>%
      slice(df)
    
    a <- array(NA, c(nrow(XN),nrow(beta),length(selected_cov)))
    
  }
  
  if(gps == FALSE){
    
    selected_cov <- colnames(cov_obj$cov_no_scale) 
    
    XN <- const_obj$XN %>% as.matrix()
    
    cov_effect <-  cov_obj$cov_no_scale
    
    a <- array(NA, c(nrow(XN),nrow(beta),length(selected_cov)))
  }
  
  for(j in 1:length(selected_cov)){
    
    
    
    for(i in 1:nrow(beta)){
      
      a[,i,j] <- exp(as.vector(XN[,(2*j):(2*j+1)] %*% as.vector(as.matrix(beta[i,(2*j):(2*j+1)]))))
      
    }
    
    cov_effect <- cov_effect %>%
      mutate(!!paste0("pred_mean_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.5))) %>%
      mutate(!!paste0("pred_025_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.025))) %>%
      mutate(!!paste0("pred_975_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.975))) %>%
      mutate(!!paste0("pred_25_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.25))) %>%
      mutate(!!paste0("pred_75_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.75)))
    
  }
  
  df_mean <- cov_effect %>%
    select(contains("pred_mean")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(mean = value) %>%
    select(-variable)
  
  df_025 <- cov_effect %>%
    select(contains("pred_025")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_025 = value) %>%
    select(-variable)
  
  df_975 <- cov_effect %>%
    select(contains("pred_975")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_975 = value) %>%
    select(-variable)
  
  df_25 <- cov_effect %>%
    select(contains("pred_25")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_25 = value) %>%
    select(-variable)
  
  df_75 <- cov_effect %>%
    select(contains("pred_75")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_75 = value) %>%
    select(-variable)
  
  color <- c("#8f5d5d","#81b29a","#f2cc8f","#babf95","#3d405b", "#f8a07e", "#fac484" ,"#005b6e", "#b1c7b3", "#f3e79b", "#eb7f86")
  
  out_cov <- cov_effect %>%
    select(-contains("pred")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    cbind(df_mean, df_025, df_975, df_25, df_75)
  
  a <- ggplot(data = out_cov, aes(x = value, y = mean)) +
    geom_point(aes(color = variable), size = 0.5) +
    geom_ribbon(aes(ymin = q_025, ymax = q_975, fill = variable), alpha = 0.3) +  # Intervalle de confiance
    geom_ribbon(aes(ymin = q_25, ymax = q_75, fill = variable), alpha = 0.4) +  # Intervalle de confiance
    facet_wrap(~ variable, scales = "free" ) +
    scale_color_manual(values = color[1:length(selected_cov)])+
    scale_fill_manual(values = color[1:length(selected_cov)])+
    #ylim(0,20) +
    theme_minimal() +
    ylab("Estimated effect") +
    xlab("Covariates value") +
    theme(legend.position = 'none') 
  
  return(a)
  
}

make_cov_prediction_int <- function(cov_obj, const_obj, out_obj, gps = FALSE) {
  
  
  beta <- out_obj$chain1 %>% as.data.frame() %>%
    select(starts_with("beta[")) %>%
    rbind(out_obj$chain2 %>%
            as.data.frame() %>% select(starts_with("beta[")))
  
  selected_cov <- cov_obj$cov_scaled %>%
    select(!grid_id) %>%
    colnames()
  
  XN <- const_obj$XN %>% as.matrix()
  
  cov_effect <-  cov_obj$cov_no_scale %>% select(!grid_id)
  
  a <- array(NA, c(nrow(XN),nrow(beta),length(selected_cov)))
  
  for(j in 1:length(selected_cov)){
    
    
    for(i in 1:nrow(beta)){
      
      a[,i,j] <- exp(as.vector(XN[,(2*j - 1):(2*j)] %*% as.vector(as.matrix(beta[i,(2*j -1):(2*j)]))))
      
    }
    
    cov_effect <- cov_effect %>%
      mutate(!!paste0("pred_mean_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.5))) %>%
      mutate(!!paste0("pred_025_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.025))) %>%
      mutate(!!paste0("pred_975_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.975))) %>%
      mutate(!!paste0("pred_25_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.25))) %>%
      mutate(!!paste0("pred_75_", selected_cov[j]) := apply(a[,,j], 1, function(x) quantile(x, probs = 0.75)))
    
  }
  
  df_mean <- cov_effect %>%
    select(contains("pred_mean")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(mean = value) %>%
    select(-variable)
  
  df_025 <- cov_effect %>%
    select(contains("pred_025")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_025 = value) %>%
    select(-variable)
  
  df_975 <- cov_effect %>%
    select(contains("pred_975")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_975 = value) %>%
    select(-variable)
  
  df_25 <- cov_effect %>%
    select(contains("pred_25")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_25 = value) %>%
    select(-variable)
  
  df_75 <- cov_effect %>%
    select(contains("pred_75")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(variable = str_remove(variable, "pred_")) %>%
    rename(q_75 = value) %>%
    select(-variable)
  
  color <- c("#8f5d5d","#81b29a","#f2cc8f","#babf95","#3d405b", "#f8a07e", "#fac484" ,"#005b6e", "#b1c7b3", "#f3e79b", "#eb7f86")
  
  out_cov <- cov_effect %>%
    select(-contains("pred")) %>%  # Sélectionne les colonnes contenant "pred"
    pivot_longer(cols = everything(), 
                 names_to = "variable", 
                 values_to = "value") %>%
    cbind(df_mean, df_025, df_975, df_25, df_75)
  
  a <- ggplot(data = out_cov, aes(x = value, y = mean)) +
    geom_point(aes(color = variable), size = 0.5) +
    geom_ribbon(aes(ymin = q_025, ymax = q_975, fill = variable), alpha = 0.3) +  # Intervalle de confiance
    geom_ribbon(aes(ymin = q_25, ymax = q_75, fill = variable), alpha = 0.4) +  # Intervalle de confiance
    facet_wrap(~ variable, scales = "free" ) +
    scale_color_manual(values = color[1:length(selected_cov)])+
    scale_fill_manual(values = color[1:length(selected_cov)])+
    #ylim(0,20) +
    theme_minimal() +
    ylab("Estimated effect") +
    xlab("Covariates value") +
    theme(legend.position = 'none') 
  
  return(a)
}

# ----------------------- Visualisation ---------------------

## Models results ----


plot_spatial_prediction <- function(cov_obj, grid_obj, out_obj, xlim_obj = xlim, ylim_obj = ylim, world_obj = europe, integrated = FALSE){
  
  if(integrated == FALSE) {
     grid_pred_obj <- make_spatial_prediction_gam(cov_obj, grid_obj, out_obj)
   }

   if(integrated == TRUE) {
    grid_pred_obj <- make_spatial_prediction_int(cov_obj, grid_obj, out_obj)
  }
  
  grid_pred_obj <- grid_pred_obj %>%
    mutate(cv = case_when(cv > 100 ~ 100,
                          TRUE ~ cv)) %>%
    mutate(mean = (mean/max(mean)))
  
  legend_position = "bottom"
  plot_title = ""
  plot_title_size = 10
  legend_title_size = 9
  axis_text_size = 6
  remove_x_axis = F
  
  
  map_theme <- theme_bw() +
    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, 
                                face = "bold", 
                                size = plot_title_size,
                                margin = margin(0,0,0,0)),
      plot.title.position = 'plot',
      axis.text = element_text(size = axis_text_size),
      axis.line = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      panel.border = element_blank()
    ) 
  
  
  
  p1 <- ggplot() +
    geom_sf(data = grid_pred_obj, aes(fill = mean), color = NA, lwd = 0) +
    geom_sf(data = world_obj) +
    xlim(xlim_obj) +
    ylim(ylim_obj) + 
    scale_fill_gradientn(colours= c("#1b263b","#009392", "#fbc99d", "#ffa86b","#ff9286","#d17b88") ,
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.50, 0.75, 1),
                         guide = guide_colourbar(
                           title = "Relative space-use",
                           title.theme = element_text(
                             face = "bold",
                             size = legend_title_size,
                             hjust = 0.5),
                           title.position = 'top',
                           ticks = FALSE,
                           title.hjust = .5,
                           barwidth = unit(6, 'lines'),
                           barheight = unit(.4, 'lines'),
                           label.theme = element_text(size = legend_title_size * 0.65))) +
    labs(title = plot_title) +
    ggspatial::annotation_scale(location = "tl", 
                                pad_y = unit(0.12, "cm"),
                                pad_x = unit(0.1, "cm"),
                                width_hint = 0.22,
                                height = unit(0.07, "cm"),
                                text_cex = 0.4) +
    map_theme
  
  
  
  p2 <- ggplot() +
    geom_sf(data = grid_pred_obj, aes(fill = cv), color = NA, lwd = 0) +
    geom_sf(data = world_obj) +
    xlim(xlim_obj) +
    ylim(ylim_obj) + 
    scale_fill_gradientn(colours= c("#f8edeb","#e8ddb5","#f7af9d","#e8a598","#b56576","#6d597a","#355070") ,
                         limits = c(0, 100),
                         breaks = c(0, 25, 50, 75, 100),
                         labels = c("0", "25", "50", "75", ">100"),
                         guide = guide_colourbar(
                           title = "",
                           title.theme = element_text(
                             face = "bold",
                             size = legend_title_size,
                             hjust = 0.5),
                           title.position = 'top',
                           ticks = FALSE,
                           title.hjust = .5,
                           barwidth = unit(0.2, 'lines'),
                           barheight = unit(2, 'lines'),
                           label.theme = element_text(size = legend_title_size * 0.5))) +
    labs(title = "Coefficient of variation") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, 
                                face = "bold", 
                                size = plot_title_size,
                                margin = margin(0,0,0,0)),
      plot.title.position = 'plot',
      axis.line = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      panel.border = element_blank()
    ) 
  
  cowplot::ggdraw() +
    cowplot::draw_plot(p1) +
    cowplot::draw_plot(p2, x = 0.6, y = 0.2, width = 0.3, height = 0.3)
  
  
}

plot_p <- function(sum_out_obj, nmix_list_obj, title = ""){ 
  
  vec <- names(nmix_list_obj)
  fact_vec <- as.factor(vec)
  cor_campagne <- data.frame(dataset_id = as.numeric(fact_vec), campagne = vec)
  
  detect <- sum_out_obj %>%
    mutate(parameter = rownames(.)) %>%
    filter(str_sub(parameter, 1, 2) == "p[") %>%
    mutate(dataset_id = as.numeric(str_extract(str_sub(parameter, 3), "\\d+"))) %>%
    left_join(cor_campagne, by = "dataset_id") %>%
    rename(IC025 = "2.5%",
           IC975 = "97.5%") 
  
  ggplot(detect, aes(x=campagne, y=mean)) +
    geom_point(size = 2) +
    geom_pointrange(aes(ymin=IC025, ymax=IC975)) +
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title,
         x = "",
         y = "Detection probability")
}


##  Posterior Predictive Checks ------------------------------------------------

retrieve_fit_values <- function(mcmc.output_obj){
  mcmc.output_obj %>% map( ~ {
    .x %>% as_tibble() %>% janitor::clean_names()
  }) %>%
    bind_rows() %>%
    select(starts_with("fit")) %>%
    mutate(id = 1:nrow(.)) %>%
    pivot_longer(-id) %>%
    mutate(type = ifelse(str_detect(string = name, pattern = "fit_rep"), "fit_rep", "fit"),
           data_source = str_sub(name,-1)) %>%
    select(-name) %>%
    pivot_wider(names_from = type, values_from = value) 
}


compute_pvalue_nmix <- function(tibble_Nmix){
  
  pval_df <- tibble_Nmix %>% 
    group_by(data_source) %>% 
    summarise(pvalue = round(sum(fit_rep > fit) / length(fit), 2))
  
  return(pval_df)
}


plot_PPC <- function(mcmc.output_obj, nmix_list_obj, title = ""){ 
  
  vec <- names(nmix_list_obj)
  fact_vec <- as.factor(vec)
  cor_campagne <- data.frame(data_source = as.character(as.numeric(fact_vec)), campagne = vec)
  
  tibble_Nmix <- retrieve_fit_values(mcmc.output_obj) %>%
    left_join(cor_campagne, by = "data_source")
  
  pvalues_df <- tibble_Nmix %>%
    group_by(data_source) %>%
    summarise(x_max = max(fit),
              x_min = min(fit),
              y_max = max(fit_rep),
              y_min = min(fit_rep)) %>%
    mutate(x = 0.1 * (x_max - x_min) + x_min,
           y = 0.9 * (y_max - y_min) + y_min) %>%
    full_join(compute_pvalue_nmix(tibble_Nmix = tibble_Nmix),
              by = join_by(data_source)) %>%
    left_join(cor_campagne, by = "data_source")
  
  
  ggplot(data = tibble_Nmix, aes(x = fit, y = fit_rep, color = fit<fit_rep)) +
    geom_point(shape = 20) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.75) +
    facet_wrap(~campagne, scales = "free") +
    geom_text(data = pvalues_df, mapping = aes(x = x, y = y, label = pvalue), color = "black") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = title,
         x = "Observed values",
         y = "Simulated values")
}

## Covariates ----
plot_cov <- function(grid_obj, selected_cov_obj, xlim_obj = xlim, ylim_obj = ylim, world_obj = europe){
  
  
  df_pred <- grid_obj %>%
    mutate(across(all_of(selected_cov_obj), ~ (. - mean(.)) / sd(.))) %>%
    pivot_longer(cols = selected_cov_obj, 
                 names_to = "variable", 
                 values_to = "value") %>%
    mutate(fill_scaled = value) %>%
    group_by(variable) 
  
  legend_position = "bottom"
  plot_title = ""
  plot_title_size = 10
  legend_title_size = 9
  axis_text_size = 6
  remove_x_axis = F
  
  
  map_theme <- theme_bw() +
    theme(
      legend.position = legend_position,
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, 
                                face = "bold", 
                                size = plot_title_size,
                                margin = margin(0,0,0,0)),
      plot.title.position = 'plot',
      axis.text = element_text(size = axis_text_size),
      axis.line = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      panel.border = element_blank()
    ) 
  
  ggplot() +
    geom_sf(data = df_pred, aes(fill = fill_scaled), color = NA, lwd = 0) +
    geom_sf(data = world_obj) +
    xlim(xlim_obj) +
    ylim(ylim_obj) + 
    facet_wrap(~ variable) +
    scale_fill_gradientn(colours= c("#1b263b","#005b6e","#009392", "#72aaa1", "#b1c7b3", "#f3e79b", "#fac484", "#f8a07e", "#eb7f86") ,
                         limits = c(-5,5),
                         guide = guide_colourbar(
                           title = "Covariate intensity",
                           title.theme = element_text(
                             face = "bold",
                             size = legend_title_size,
                             hjust = 0.5),
                           title.position = 'top',
                           ticks = FALSE,
                           title.hjust = .5,
                           barwidth = unit(6, 'lines'),
                           barheight = unit(.4, 'lines'),
                           label.theme = element_text(size = legend_title_size * 0.65))) +
    labs(title = plot_title) +
    ggspatial::annotation_scale(location = "br", 
                                pad_y = unit(0.12, "cm"),
                                pad_x = unit(0.1, "cm"),
                                width_hint = 0.22,
                                height = unit(0.07, "cm"),
                                text_cex = 0.4) +
    map_theme
  
}

plot_prox_colonies <- function(colonies, species, grid_obj, world_obj = europe, xlim_obj = xlim, ylim_obj = ylim) {
  
  
  colonies_sp <- colonies %>%
    filter(Espece %in% species)
  
  
  ggplot() +
    geom_sf(data = world_obj, fill = "#DADDD8") +
    geom_sf(data = grid_obj, aes(fill = prox_colonies)) +
    geom_sf(data = colonies_sp, aes(size = EFF), color = "black", alpha = 0.5) +
    labs(size = "Nombre de couples par colonie") +
    xlim(xlim_obj) +
    ylim(6133800, 6300000) +
    scale_fill_gradientn(colours= c("#1b263b","#005b6e","#009392", "#72aaa1", "#b1c7b3", "#f3e79b", "#fac484", "#f8a07e", "#eb7f86") ,
                         guide = guide_colourbar(
                           title = "Covariate intensity",
                           title.theme = element_text(
                             face = "bold",
                             hjust = 0.5),
                           title.position = 'top',
                           ticks = FALSE,
                           title.hjust = .5,
                           barwidth = unit(2, 'lines'),
                           barheight = unit(6, 'lines'))) +
    theme(#legend.position = "none",
      panel.grid = element_line(colour = "transparent"),
      plot.title = element_text(lineheight = 0.8, face = "bold"),
      axis.text = element_text(size = 10),
      strip.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "azure"),
      panel.border = element_rect(fill = NA))
  
}


## Efforts -----

plot_effort <- function(data, title = "", xlim_obj = xlim, ylim_obj = ylim, world_obj = europe) {
  ggplot() +
    geom_sf(data = world_obj) +
    geom_sf(data = data$eff, color = "black", size = 0.1) +
    geom_sf(data = data$obs, color = "#e07a5f", size = 1) +
    xlim(xlim_obj) +
    ylim(ylim_obj) + 
    theme_minimal() +
    labs(caption = paste("Nombre de jours de suivi = ", length(unique(data$eff$date)), sep = " ")) +
    labs(title = title)
}

plot_all_effort <- function(data_list, effort = FALSE, world_obj = europe){
  
  plots_list <- imap(data_list, ~ plot_effort(.x, title = .y, world_obj = world_obj))
  
  # Combiner les graphiques en une seule figure
  cowplot::plot_grid(plotlist = plots_list)
  
}

plot_colonies <- function(colonies, species, world_obj = europe, xlim_obj = xlim, ylim_obj = ylim) {
  
  
  colonies_sp <- colonies %>%
    filter(Espece %in% species)
  
  ggplot() +
    geom_sf(data = world_obj, fill = "#DADDD8") +
    geom_sf(data = colonies_sp, aes(size = EFF), color = "#006d77", alpha = 0.5) +
    labs(size = "Nombre de couples par colonie") +
    xlim(xlim_obj) +
    ylim(6133800, 6300000) +
    theme(#legend.position = "none",
      panel.grid = element_line(colour = "transparent"),
      plot.title = element_text(lineheight = 0.8, face = "bold"),
      axis.text = element_text(size = 10),
      strip.background = element_rect(fill = "white"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "azure"),
      panel.border = element_rect(fill = NA))
  
}

## Count data ------

plot_count <- function(data, title = "", xlim_obj = xlim, ylim_obj = ylim, world_obj = europe) {
  
  ggplot() +
    geom_sf(data = world_obj) +
    geom_sf(data = data$eff, color = "black", size = 0.1) +
    geom_sf(data = data$obs, color = "#e07a5f", size = 1) +
    xlim(xlim_obj) +
    ylim(ylim_obj) + 
    theme_minimal() +
    labs(caption = paste("Sightings = ", sum(data$obs$effectif), sep = " ")) +
    labs(title = title)
  
  
}

plot_all_count <- function(data_list, effort = FALSE){
  
  plots_list <- imap(data_list, ~ plot_count(.x, title = .y))
  
  # Combiner les graphiques en une seule figure
  cowplot::plot_grid(plotlist = plots_list)
  
}

plot_obs <- function(data_list, grid_obj, data_nmix, title = "", xlim_obj = xlim, ylim_obj = ylim, world_obj = europe, in.log = TRUE) {
  
  # select geomtry et effectif et concerver le nom de la liste, tou tça un seul tableau.
  
  intersect_grid_obs <- map_dfr(names(data_list), ~ data_list[[.x]]$obs %>% 
                                  select(effectif, date, geometry) %>% 
                                  mutate(data = .x)) %>%
    st_intersection(grid_obj %>% select(grid_id)) %>%
    group_by(grid_id, date, data) %>%
    mutate(n = sum(effectif)) %>%
    ungroup(.) %>%
    distinct(grid_id, date, data, .keep_all = TRUE) 
  
  
  intersect_grid_eff <- map_dfr(names(data_list), ~ data_list[[.x]]$eff %>% 
                                  select(date,geometry) %>% 
                                  mutate(data = .x)) %>%
    st_intersection(grid_obj %>% select(grid_id)) 
  
  # Filtrer les observations en fonction des cellules échantillonnées
  df_final <- intersect_grid_obs %>%
    semi_join(intersect_grid_eff %>% st_drop_geometry(), by = c("date", "data", "grid_id"))
  
  
  
  
  
  if(in.log == TRUE){
    
    # Définir les 5 valeurs réparties de manière logarithmique
    breaks_values <- exp(seq(log(min(df_final$n)), log(max(df_final$n)), length.out = 5))
    
    plot <- ggplot() +
      geom_sf(data = world_obj) +
      geom_sf(data= covariates, fill = "white") +
      geom_sf(data = intersect_grid_eff) +
      geom_sf(data = df_final, aes(color = log(n), size = log(n)), alpha = 0.9) +
      xlim(xlim_obj) +
      ylim(ylim_obj) + 
      scale_color_gradientn(colours = c("#1b263b","#005b6e","#009392", "#72aaa1", "#b1c7b3", "#f3e79b", "#fac484", "#f8a07e", "#eb7f86"), ,
                            breaks = log(breaks_values), # Choisir les valeurs affichées
                            labels = round(breaks_values),
                            name = "Sightings")  +
      scale_size_continuous(range = c(1, 10), guide = 'none') +
      labs(caption = paste("Nombre d'individus observés = ", sum(df_final$n), "\n Nombre de données pour la modélisation =", length(data_nmix$data$nobs), sep = " ")) +
      theme_minimal() 
  }
  
  if(in.log == FALSE){
    
    # Définir les 5 valeurs réparties de manière logarithmique
    breaks_values <- seq(min(df_final$n), max(df_final$n), length.out = 5)
    
    plot <- ggplot() +
      geom_sf(data = world_obj) +
      geom_sf(data= covariates, fill = "white") +
      geom_sf(data = intersect_grid_eff) +
      geom_sf(data = df_final, aes(color = n, size = n), alpha = 0.9) +
      xlim(xlim_obj) +
      ylim(ylim_obj) + 
      scale_color_gradientn(colours = c("#1b263b","#005b6e","#009392", "#72aaa1", "#b1c7b3", "#f3e79b", "#fac484", "#f8a07e", "#eb7f86"), ,
                            breaks = breaks_values, # Choisir les valeurs affichées
                            labels = round(breaks_values),
                            name = "Sightings")  +
      scale_size_continuous(range = c(1, 10), guide = 'none') +
      labs(caption = paste("Nombre d'individus observés = ", sum(df_final$n), "\n Nombre de données pour la modélisation =", length(data_nmix$data$nobs), sep = " ")) +
      theme_minimal() 
  }
  
  return(plot)
}

get_obs <- function(data, grid_obj) {
  
  
  intersect_grid_obs <- map_dfr(names(data), ~ data[[.x]]$obs %>% 
                                  select(effectif, date, nom_fr) %>% 
                                  mutate(data = .x)) %>%
    st_intersection(grid_obj %>% select(grid_id))
  
  
  intersect_grid_eff <- map_dfr(names(data), ~ data[[.x]]$eff %>% 
                                  select(date,geometry) %>% 
                                  mutate(data = .x)) %>%
    st_intersection(grid_obj %>% select(grid_id)) 
  
  # Filtrer les observations en fonction des cellules échantillonnées
  df_final <- intersect_grid_obs %>%
    semi_join(intersect_grid_eff %>% st_drop_geometry(), by = c("date", "data", "grid_id")) %>%
    rename(Campagne = data) %>%
    rename(Nom = nom_fr) %>%
    rename(Date = date) %>%
    rename(Effectif = effectif) %>%
    select(Campagne, Date, Nom, Effectif) %>%
    st_drop_geometry()
  
}

## GPS data ------
plot_gps_data <- function(data_obj, data_rsf, xlim_obj = xlim, ylim_obj = ylim, world_obj = europe){
  
  ggplot() +
    geom_sf(data = world_obj, fill = "#EDEDE9") +
    geom_path(data = data_obj, aes(x = lon, y = lat, group = ind_ID, color = ind_ID), 
              linewidth = 0.2) +
    geom_point(data = data_obj, aes(x = lon, y = lat, color = ind_ID), 
               linewidth = 0.2) +
    xlim(xlim_obj) +
    ylim(ylim_obj) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_line(colour = "transparent"),
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "azure"),
          panel.border = element_rect(fill = NA))  +
    labs(caption = paste("Nombre d'individus = ", length(unique(data_obj$ind_ID)), " \n  Nombre de points gps utilisés pour la modélisation =", sum(data_rsf$data$kase), sep = " ")) 
  
}

plot_gps_data_aug <- function(data_obj, xlim_obj = xlim, ylim_obj = ylim, world_obj = europe){
  
  ggplot() +
    geom_sf(data = world_obj) +
    geom_sf(data = data_obj %>% filter(presence == 0), color = "#81b29a") +
    geom_sf(data = data_obj %>% filter(presence == 1), color = "#8f5d5d") +
    xlim(xlim_obj) +
    ylim(ylim_obj) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_line(colour = "transparent"),
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "azure"),
          panel.border = element_rect(fill = NA))  +
    labs(caption = paste("Nombre de points GPS = ", nrow(data_obj %>% filter(presence == 1)), sep = " ")) 
}
