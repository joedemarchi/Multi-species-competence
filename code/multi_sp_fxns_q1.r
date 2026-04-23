## Multi-species Functions for Q1


### set discretized values
set_discretized_values = function(min_size, max_size, bins){
  
  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model
  
  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix
  
  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)
  
  
  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins
  
  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])
  
  # Width of cells
  h = y[2] - y[1]
  
  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
              h=h))
  
}
### This works

# # Example usage
# 
# lower = -5
# upper = 10
# bins = 100
# 
# disc <- set_discretized_values(lower, upper, bins)


##____________________________________________________________________________##
### Parameter functions

# Growth Function --------------------------------------------------------------
growth_fxn = function(xlower, xupper, xnow, params){
  # The Bd growth function on an infected individual as a CDF. 
  #
  # Uses the CDF to calculate the probability of the growth function 
  #
  # Parameters
  # ----------
  # xlower, xupper : float
  # 	Specify the interval into which growth occurs
  # xnow : The midpoint of the current load class
  # params : list
  # 	A list of parameters
  #
  # Returns 
  # -------
  # : probability of transitioning from x_now to [x_next - h/2, x_next + h/2]
  
  mu = params$growth_int + params$growth_slope * xnow
  
  sigma = params$growth_sigma
  prob = pnorm(xupper, mean=mu, sd=sigma) - pnorm(xlower, mean=mu, sd=sigma)
  return(prob)
  
}

# Infection Loss Function --------------------------------------------------------------
loss_fxn = function(xnow, params) {
  
  # Parameters
  # -----------
  # xnow: float
  #The (log) infection load at time t
  # params : list
  #The parameters for the IPM
  #
  # Returns
  # -------
  # : probability of losing infection over a time step
  
  logit = params$loss_int + params$loss_slope*xnow
  prob = 1 / (1 + exp(-logit))
  
  # # scale to 6 days if R.clamitans
  # if(species == "R.clamitans"){
  #   prob = 1 - (1 - prob)^(6/4)
  #   
  # }
  return(prob)
  
}

loss_sp_fxn = function(xnow, params, species) {
  
  # Parameters
  # -----------
  # xnow: float
  #The (log) infection load at time t
  # params : list
  #The parameters for the IPM
  #
  # Returns
  # -------
  # : probability of losing infection over a time step
  
  logit = params$loss_int + params$loss_slope*xnow
  prob = 1 / (1 + exp(-logit))
  
  # scale to 6 days if R.clamitans
  if(species == "R.clamitans"){
    prob = 1 - (1 - prob)^(6/4)

  }
  return(prob)
  
}

## Initial Infection Function --------------------------------------------------------------
init_inf_fxn = function(xlower, xupper, params){
  # The initial infection function for an individual (independent of dose)
  # 
  # Parameters
  # ----------
  # xlower, xupper : float or array
  # 	The lower and upper loads of the bins
  # params : list
  #	A list of parameters for the IPM
  # 
  # Returns
  # -------
  # : probability of gaining an initial infection between xlower and xupper
  
  mu = params$init_inf_int
  sigma = params$init_inf_sigma
  prob = pnorm(xupper, mean=mu, sd=sigma) - pnorm(xlower, mean=mu, sd=sigma)
  return(prob)
  
}

# Survival Function (Don't need) --------------------------------------------------------------
surv_fxn = function(xnow, params){
  # The survival probability given infection load xnow
  #
  # Parameters
  # -----------
  # xnow: float
  #	The (log) infection load at time t
  # params : list
  #	The parameters for the IPM
  #
  # Returns
  # -------
  # : probability of survival over a time step
  
  logit = params$surv_int + params$surv_slope*xnow
  prob = 1 / (1 + exp(-logit))
  return(prob)
}


## Transmission Functions --------------------------------------------------------------
trans_fxn = function(x, params){
  # The infection given contact function
  #
  # Parameters
  # ----------
  # x : float
  # 	Concentration/dose of zoospores in the environment.
  #   On the log_e scale
  # params : list
  # 	List of IPM model parameters
  #
  # Return
  # ------
  # : probability of infection given contact
  
  logit = params$trans_int + params$trans_slope*x
  prob = 1 / (1 + exp(-logit))
  
  # # scale to 6 days if R.clamitans
  # if(species == "R.clamitans"){
  #   prob = 1 - (1 - prob)^(6/4) 
  #   
  # }
  # 
  return(prob)
}

trans_sp_fxn = function(x, params, species){
  # The infection given contact function
  #
  # Parameters
  # ----------
  # x : float
  # 	Concentration/dose of zoospores in the environment.
  #   On the log_e scale
  # params : list
  # 	List of IPM model parameters
  #
  # Return
  # ------
  # : probability of infection given contact
  
  logit = params$trans_int + params$trans_slope*x
  prob = 1 / (1 + exp(-logit))
  
  # scale to 6 days if R.clamitans
  if(species == "R.clamitans"){
    prob = 1 - (1 - prob)^(6/4)

  }

  return(prob)
}


## Exponential link
trans_exp_fxn = function(x, params, species){
  # The infection given contact function
  #
  # Parameters
  # ----------
  # x : float
  # 	Concentration/dose of zoospores in the environment

  # params : list
  # 	List of IPM model parameters
  #
  # Return
  # ------
  # : probability of infection per zoospore contact
  # with theta per zoospore(params$trans_exp) = exp(theta) / 1e6
    prob = 1 -exp(- params$trans_exp * x)
  
  # scale to 6 days if R.clamitans
  if(species == "R.clamitans"){
    prob = 1 - (1 - prob)^(6/4)

  }

  return(prob)
}

### ---------------------------------------------------------------------------
# standardize pond volume
# Calculate pond volume
calculate_pond_volume <- function(diameter, avg_depth) {
  # Calculate pond volume (m^3) assuming circular pond
  radius <- diameter / 2
  volume <- pi * (radius^2) * avg_depth
  volume_liters <- volume * 1000 # convert to liters
  return(volume_liters)
}

##_______________________________________________________________________
# Movement by volume
movement_volume <- function(ipm_params, contact_depth){
  
  #convert area moved per day (m^2) in volume searched (m^3)
  movement_m3 <- ipm_params$movement_distance * contact_depth

  #convert to Liters
  movement_m3 <- movement_m3 * 1000 # convert to liters
  
  return(movement_m3)
}

### density-dependent metamorphosis function

dd_meta_fxn <- function(N_adult, K_meta){
  1/(1 + N_adult/K_meta)
}

### ___________________________________________________________________________##
## Build IPM kernels
#### Multi stage NGM #####################################################

build_FU_bullfrog_2stage <- function(params_A, params_T,
                                     pond_diameter,
                                     lower, upper, bins,
                                     contact_depth = 0.5,
                                     avg_depth = 0.5,
                                     species_A ,
                                     species_T ,
                                     tad_to_adult_map = NULL
) {
  
  disc <- set_discretized_values(lower, upper, bins)
  xlower <- disc$bnd[1:bins]
  xupper <- disc$bnd[2:(bins + 1)]
  y      <- disc$y
  
  # Pond volume (L) and contact volume per step (L/step)
  pond_volume_L <- calculate_pond_volume(pond_diameter, avg_depth = avg_depth)
  movement_L_A  <- movement_volume(params_A, contact_depth = contact_depth)
  movement_L_T  <- movement_volume(params_T, contact_depth = contact_depth)
  
  # Linearized "beta0" term per step (your exponential link version)
  beta0_A <- params_A$trans_exp * (movement_L_A / pond_volume_L)
  beta0_T <- params_T$trans_exp * (movement_L_T / pond_volume_L)
  
  # Susceptible densities (DFE)
  S_A <- if (is.null(params_A$df_density)) 1 else params_A$df_density
  S_T <- if (is.null(params_T$df_density)) 1 else params_T$df_density
  
  #---------------------------------
  # No density dependence + emigration at metamorphosis
  #---------------------------------
  
  # emigration fraction of *metamorphosing* tadpoles that leave the focal population
  emig_meta <- params_A$emig_meta
  if (is.null(emig_meta) || length(emig_meta) == 0 || !is.finite(emig_meta)) emig_meta <- 0
  emig_meta <- max(0, min(0.999, emig_meta))
  stay_meta <- 1 - emig_meta
  
  # per-step metamorphosis probability (already on 6-day step scale in your params)
  meta_prob <- if (!is.null(params_T$meta_successB)) params_T$meta_successB else params_T$meta_successG
  
  # NO DD multiplier
  meta_eff <- meta_prob
  
  # NO DD tadpole survival (just per-step survival)
  tad_surv <- params_T$s0
  
  # “metamorphosing infected tadpoles survive & metamorphose” term (scaled by emigration)
  meta_tad <- tad_surv * meta_eff * stay_meta
  
  ##--------------------------------------------------------

  ##--------------------------------------------------------
  
  # Initial infection load distribution into bins
  G0_A <- init_inf_fxn(xlower, xupper, params_A)
  G0_T <- init_inf_fxn(xlower, xupper, params_T)

  # I -> Z (shedding) per step (keep your scaling choices)
  ItoZ_A <- exp(params_A$lambda_int) * exp(y)^params_A$lambda_slope
  ItoZ_T <- exp(params_T$lambda_int) * exp(y)^params_T$lambda_slope
  
  # Convert to "total shed captured" scaling you used:
  ItoZ_A <- 3 * ItoZ_A          # 0.3L -> 1L scaling 
  ItoZ_T <- 3 * ItoZ_T
  
  Zvect_A <- ItoZ_A * 6         # 1 day -> 6 day step scaling 
  Zvect_T <- ItoZ_T * 6
  
  # Dimensions: [IA (bins), IT (bins), Z (1)]  => total = 2*bins + 1
  n <- 2 * bins + 1
  idx_A <- 1:bins
  idx_T <- (bins + 1):(2 * bins)
  idx_Z <- 2 * bins + 1
  
  # --------------------------
  # F matrix (new infection edges)
  # --------------------------
  F <- matrix(0, n, n)
  
  
  # Z -> I (new infections)
  # adult infections
  StoI_A <- params_A$s0 * S_A * beta0_A * G0_A
  F[idx_A, idx_Z] <- StoI_A 
  # tadpole infections Z -> infected tadpoles
  StoI_T <- tad_surv * (1-meta_eff) * S_T * beta0_T * G0_T
  F[idx_T, idx_Z] <- StoI_T
  # tadpole infects and metamorphose in same step
  if(is.null(tad_to_adult_map)){
    tad_to_adult_map <- diag(bins)
  }
  # --- tadpoles that metamorphose AND get infected in the same step become infected adults ---
  if (is.null(tad_to_adult_map)) tad_to_adult_map <- diag(bins)
  
  StoI_T_meta <- tad_surv * meta_eff * stay_meta * S_T * beta0_T * G0_T   # length = bins
  F[idx_A, idx_Z] <- F[idx_A, idx_Z] + (tad_to_adult_map %*% StoI_T_meta)
  
  
  # infected adults / tadpoles -> Z
  F[idx_Z, idx_A] <- Zvect_A
  F[idx_Z, idx_T] <- Zvect_T
  
  ### Expand out f matrix for right size
  
  # --------------------------
  # U matrix (within-infected transitions)
  # --------------------------
  U <- matrix(0, n, n)
  
  # Adult ItoI block
  G_A <- sapply(y, function(x) growth_fxn(xlower, xupper, x, params_A))
  L_A <- loss_sp_fxn(y, params_A, species_A)
  sI_A <- params_A$s0
  ItoI_A <- G_A %*% diag(sI_A * (1 - L_A))
  
  # Tadpole ItoI block
  G_T <- sapply(y, function(x) growth_fxn(xlower, xupper, x, params_T))
  L_T <- loss_sp_fxn(y, params_T, species_T)
  sI_T <- params_T$s0
  ItoI_T <- G_T %*% diag(tad_surv * (1-meta_eff) * (1 - L_T))
  
  # TI -> AI metamorphose infected
  if (is.null(tad_to_adult_map)) tad_to_adult_map <- diag(bins)
  U[idx_A, idx_T] <- tad_to_adult_map %*% (G_T %*% diag(meta_tad * (1 - L_T)))
  
  ### ADD in a TI to AI block and a block for AI to TI of all zeros
  
  U[idx_A, idx_A] <- ItoI_A
  U[idx_T, idx_T] <- ItoI_T
  
  # Z -> Z decay (within-environment)
  zsurv <- exp(-0.3 * 6)   # zoospore decay
  
  # Z to Z decay
  U[idx_Z, idx_Z] <- zsurv
  
  list(F = F, U = U,
       indexing = list(idx_A = idx_A, idx_T = idx_T, idx_Z = idx_Z),
       disc = disc,
       dd = list(meta_eff = meta_eff, tad_surv = tad_surv, meta_tad = meta_tad))
}


compute_R0_from_FU <- function(F, U) {
  I <- diag(nrow(F))
  K <- F %*% solve(I - U)
  # vals <- max(Re(eigen(K)$values))
  # R0 <- vals^2  # square because shedding is in F matrix
  # return(R0)
  
  vals <- eigen(K, only.values = TRUE)$values
  
  r0 <- max(Mod(vals))
  
  R0 <- r0^2
  
  
  return(R0)
}

compute_bullfrog_2stage_R0_current <- function(ipm_d,
                                               adult_name,
                                               tadpole_name,
                                               pond_diameter,
                                               lower, upper, bins) {
  
  params_A <- ipm_d[[adult_name]]
  params_T <- ipm_d[[tadpole_name]]
  
  mats <- build_FU_bullfrog_2stage(
    params_A = params_A,
    params_T = params_T,
    pond_diameter = pond_diameter,
    lower = lower, upper = upper, bins = bins,
    species_A = adult_name,
    species_T = tadpole_name
  )
  
  compute_R0_from_FU(mats$F, mats$U)
}

## Map tadpole to corresponding adult

is_tad_label <- function(sp) grepl("tadpole", sp, ignore.case = TRUE)

tad_parent <- function(sp) {
  # tadpole labels to correct adult label
  if (sp == "Racl.tadpole") return("R.clamitans")
  if (sp == "Raca.tadpole") return("R.catesbeianus")
  stop("Unknown tadpole label: ", sp)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Build a community FU builder:
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

build_FU_community <- function(species_list, ipm_dists,
                                       pond_diameter, lower, upper, bins,
                                       contact_depth = 0.5, avg_depth = 0.5) {
  
  disc <- set_discretized_values(lower, upper, bins)
  xlower <- disc$bnd[1:bins]
  xupper <- disc$bnd[2:(bins + 1)]
  y      <- disc$y
  
  pond_volume_L <- calculate_pond_volume(pond_diameter, avg_depth = avg_depth)
  zsurv <- exp(-0.3 * 6)
  
  n_blocks <- length(species_list)
  n <- n_blocks * bins + 1
  idx_Z <- n
  
  # block indices
  idx_I <- setNames(vector("list", n_blocks), species_list)
  for (i in seq_along(species_list)) {
    idx_I[[i]] <- ((i - 1) * bins + 1):(i * bins)
  }
  
  F <- matrix(0, n, n)
  U <- matrix(0, n, n)
  U[idx_Z, idx_Z] <- zsurv
  
  for (sp in species_list) {
    
    ii <- idx_I[[sp]]
    
    # Choose the correct parameter set
    if (!is_tad_label(sp)) {
      pars <- ipm_dists[[sp]]
      if (is.null(pars)) stop("No params for species: ", sp)
      
      move_L <- movement_volume(pars, contact_depth = contact_depth)
      beta0 <- pars$trans_exp * (move_L / pond_volume_L)
      S_star <- pars$df_density  # naive DFE
      if (is.null(S_star)) stop("df_density missing for ", sp)
      
      
      # Z->I
      G0 <- init_inf_fxn(xlower, xupper, pars)
      StoI <- pars$s0 * S_star * beta0 * G0
      
      # I->Z
      ItoZ <- 3 * (exp(pars$lambda_int) * exp(y)^pars$lambda_slope) * 6
      
      
      # ItoI
      G <- sapply(y, function(x) growth_fxn(xlower, xupper, x, pars))
      L <- loss_sp_fxn(y, pars, sp)
      ItoI <- G %*% diag(pars$s0 * (1 - L))
      
      
    } else {
      # Tadpole block
      parent <- tad_parent(sp)
      parsA <- ipm_dists[[parent]]     # for DFE density
      parsT <- ipm_dists[[sp]]         # tad parameters
      if (is.null(parsA)) stop("No adult params for parent: ", parent)
      if (is.null(parsT)) stop("No tadpole params for: ", sp)
      
      move_L <- movement_volume(parsT, contact_depth = contact_depth)
      beta0 <- parsT$trans_exp * (move_L / pond_volume_L)
      
      # --- DFE densities---
      N_ad_DFE <- parsA$df_density
      tad_DFE <- parsT$df_density
      if (is.null(N_ad_DFE)) stop("df_density missing for parent adult: ", parent)
      if (is.null(tad_DFE))  stop("df_density missing for tadpole: ", sp)
      
      
      #------ Density dependnet metamorphosis --------
      meta_mult <- dd_meta_fxn(N_ad_DFE, parsA$K_meta)
      
      # species-specific tad K + metamorphosis prob (per step)
      if (parent == "R.clamitans") {
        tad_K    <- parsT$K_g
        meta_prob <- parsT$meta_successG
      } else if (parent == "R.catesbeianus") {
        tad_K    <- parsT$K_b
        meta_prob <- parsT$meta_successB
      } else {
        stop("Unknown tad parent: ", parent)
      }
      
      meta_eff <- meta_prob * meta_mult
      
      # ------- density dependent tadpole survival ------
      
      tad_surv <- dd_surv(parsT$s0, tad_DFE, tad_K)
      
      # Naive invasion population size
      S_star <- tad_DFE
      
      #ZtoI scaled by metamorphosis and survival
      G0 <- init_inf_fxn(xlower, xupper, parsT)
      StoI <- (tad_surv * (1 - meta_eff)) * S_star * beta0 * G0
      
    
      # ItoZ shedding
      ItoZ <- (tad_surv * (1 - meta_eff)) * 3 * (exp(parsT$lambda_int) * exp(y)^parsT$lambda_slope) * 6
      
      # ItoI
      G <- sapply(y, function(x) growth_fxn(xlower, xupper, x, parsT))
      L <- loss_sp_fxn(y, parsT, sp)
      ItoI <- G %*% diag((tad_surv * (1 - meta_eff)) * (1 - L))
    }
    
    # fill community matrices
    F[ii, idx_Z] <- StoI
    F[idx_Z, ii] <- ItoZ
    U[ii, ii] <- ItoI
  }
  
  list(F = F, U = U)
}


compute_community_R0_multi <- function(species_list, ipm_dists, lower, upper, bins, pond_diameter) {
  
  mats <- build_FU_community(
    species_list = species_list,
    ipm_dists = ipm_dists,
    pond_diameter = pond_diameter,
    lower = lower, upper = upper, bins = bins
  )
  
  F <- mats$F; U <- mats$U
  K <- F %*% solve(diag(nrow(F)) - U)
  vals <- eigen(K, only.values = TRUE)$values
  r0 <- max(Mod(vals))
  R0 <- r0^2
}

########----------Adult only FU----------------###################


build_FU_adult_only <- function(params_A,
                                pond_diameter,
                                lower, upper, bins,
                                contact_depth = 0.5,
                                avg_depth = 0.5,
                                species_A) {
  
  disc <- set_discretized_values(lower, upper, bins)
  xlower <- disc$bnd[1:bins]
  xupper <- disc$bnd[2:(bins + 1)]
  y      <- disc$y
  
  pond_volume_L <- calculate_pond_volume(pond_diameter, avg_depth = avg_depth)
  movement_L_A  <- movement_volume(params_A, contact_depth = contact_depth)
  beta0_A <- params_A$trans_exp * (movement_L_A / pond_volume_L)
  
  # DFE susceptible density
  S_A <- if (is.null(params_A$df_density)) 1 else params_A$df_density
  
  # initial infection load distribution
  G0_A <- init_inf_fxn(xlower, xupper, params_A)   # length bins
  
  # shedding
  ItoZ_A <- exp(params_A$lambda_int) * exp(y)^params_A$lambda_slope
  ItoZ_A <- 3 * ItoZ_A
  Zvect_A <- ItoZ_A * 6
  
  # Dimensions: [IA (bins), Z (1)] => bins + 1
  n <- bins + 1
  idx_A <- 1:bins
  idx_Z <- bins + 1
  
  # --------------------------
  # F matrix
  # --------------------------
  F <- matrix(0, n, n)
  
  StoI_A <- params_A$s0 * S_A * beta0_A * G0_A
  F[idx_A, idx_Z] <- StoI_A
  F[idx_Z, idx_A] <- Zvect_A
  
  # --------------------------
  # U matrix
  # --------------------------
  U <- matrix(0, n, n)
  
  G_A <- sapply(y, function(x) growth_fxn(xlower, xupper, x, params_A))
  L_A <- loss_sp_fxn(y, params_A, species_A)
  ItoI_A <- G_A %*% diag(params_A$s0 * (1 - L_A))
  
  zsurv <- exp(-0.3 * 6)
  
  U[idx_A, idx_A] <- ItoI_A
  U[idx_Z, idx_Z] <- zsurv
  
  list(F = F, U = U,
       indexing = list(idx_I = idx_A, idx_Z = idx_Z),
       disc = disc)
}


#-------------------------------Build Tadpole only FU------------#####

build_FU_tadpole_only <- function(params_T,
                                  pond_diameter,
                                  lower, upper, bins,
                                  contact_depth = 0.5,
                                  avg_depth = 0.5,
                                  species_T,
                                  include_metamorphosis_loss = TRUE) {
  
  disc <- set_discretized_values(lower, upper, bins)
  xlower <- disc$bnd[1:bins]
  xupper <- disc$bnd[2:(bins + 1)]
  y      <- disc$y
  
  pond_volume_L <- calculate_pond_volume(pond_diameter, avg_depth = avg_depth)
  movement_L_T  <- movement_volume(params_T, contact_depth = contact_depth)
  beta0_T <- params_T$trans_exp * (movement_L_T / pond_volume_L)
  
  # DFE susceptible density (tadpoles)
  S_T <- if (is.null(params_T$df_density)) 1 else params_T$df_density
  
  # per-step metamorphosis probability (NO density dependence)
  if (species_T == "Raca.tadpole") meta_prob <- params_T$meta_successB
  if (species_T == "Racl.tadpole") meta_prob <- params_T$meta_successG
  # meta_prob <- if (is.null(meta_prob) || !is.finite(meta_prob)) 0 else meta_prob
  # meta_prob <- max(0, min(0.999, meta_prob))
  meta_eff <- if (isTRUE(include_metamorphosis_loss)) meta_prob else 0
  
  # per-step tadpole survival
  tad_surv <- params_T$s0
  
  # initial infection load distribution
  G0_T <- init_inf_fxn(xlower, xupper, params_T)   # length bins
  
  # shedding
  ItoZ_T <- exp(params_T$lambda_int) * exp(y)^params_T$lambda_slope
  ItoZ_T <- 3 * ItoZ_T
  Zvect_T <- ItoZ_T * 6
  
  # Dimensions: [IT (bins), Z (1)] => bins + 1
  n <- bins + 1
  idx_T <- 1:bins
  idx_Z <- bins + 1
  
  # --------------------------
  # F matrix
  # --------------------------
  F <- matrix(0, n, n)
  
  StoI_T <- tad_surv * (1 - meta_eff) * S_T * beta0_T * G0_T
  F[idx_T, idx_Z] <- StoI_T
  F[idx_Z, idx_T] <- Zvect_T
  
  # --------------------------
  # U matrix
  # --------------------------
  U <- matrix(0, n, n)
  
  G_T <- sapply(y, function(x) growth_fxn(xlower, xupper, x, params_T))
  L_T <- loss_sp_fxn(y, params_T, species_T)
  ItoI_T <- G_T %*% diag(tad_surv * (1 - meta_eff) * (1 - L_T))
  
  zsurv <- exp(-0.3 * 6)
  
  U[idx_T, idx_T] <- ItoI_T
  U[idx_Z, idx_Z] <- zsurv
  
  list(F = F, U = U,
       indexing = list(idx_I = idx_T, idx_Z = idx_Z),
       disc = disc,
       dd = list(meta_eff = meta_eff, tad_surv = tad_surv))
}


###---------------------- Compute R0 toggle---------------------#

compute_R0_stage <- function(ipm_d,
                             stage = c("adult", "tadpole", "full"),
                             adult_name,
                             tadpole_name,
                             pond_diameter,
                             lower, upper, bins,
                             include_metamorphosis_loss = TRUE) {
  
  stage <- match.arg(stage)
  
  if (stage == "adult") {
    mats <- build_FU_adult_only(
      params_A = ipm_d[[adult_name]],
      pond_diameter = pond_diameter,
      lower = lower, upper = upper, bins = bins,
      species_A = adult_name
    )
    return(compute_R0_from_FU(mats$F, mats$U))
  }
  
  if (stage == "tadpole") {
    mats <- build_FU_tadpole_only(
      params_T = ipm_d[[tadpole_name]],
      pond_diameter = pond_diameter,
      lower = lower, upper = upper, bins = bins,
      species_T = tadpole_name,
      include_metamorphosis_loss = include_metamorphosis_loss
    )
    return(compute_R0_from_FU(mats$F, mats$U))
  }
  
  if(stage == "full") {
    mats <- build_FU_bullfrog_2stage(
      params_A = ipm_d[[adult_name]],
      params_T = ipm_d[[tadpole_name]],
      pond_diameter = pond_diameter,
      lower = lower, upper = upper, bins = bins,
      species_A = adult_name,
      species_T = tadpole_name
    )
    return(compute_R0_from_FU(mats$F, mats$U))
  }
}

##------------------------------------------------------------------------------------##
#-------------------------Build full community R0 -----------------------------------##
##-----------------------------------------------------------------------------------##

build_FU_stage_community <- function(species_list, ipm_dists,
                                     pond_diameter, lower, upper, bins,
                                     contact_depth = 0.5, avg_depth = 0.5,
                                     tad_to_adult_map = NULL,
                                     include_infect_then_meta = TRUE,
                                     include_emigration = TRUE,
                                     z_decay = 0.3) {
  
  disc <- set_discretized_values(lower, upper, bins)
  xlower <- disc$bnd[1:bins]
  xupper <- disc$bnd[2:(bins + 1)]
  y      <- disc$y
  
  pond_volume_L <- calculate_pond_volume(pond_diameter, avg_depth = avg_depth)
  zsurv <- exp(-z_decay * 6)
  
  # Build block indices (each listed "species" gets bins slots)
  n_blocks <- length(species_list)
  n <- n_blocks * bins + 1
  idx_Z <- n
  
  idx_I <- setNames(vector("list", n_blocks), species_list)
  for (i in seq_along(species_list)) {
    idx_I[[ species_list[[i]] ]] <- ((i - 1) * bins + 1):(i * bins)
  }
  
  F <- matrix(0, n, n)
  U <- matrix(0, n, n)
  U[idx_Z, idx_Z] <- zsurv
  
  # default mapping: same bin -> same bin
  if (is.null(tad_to_adult_map)) tad_to_adult_map <- diag(bins)
  
  for (sp in species_list) {
    
    ii <- idx_I[[sp]]
    
    if (!is_tad_label(sp)) {
      # -------------------------
      # ADULT BLOCK
      # -------------------------
      parsA <- ipm_dists[[sp]]
      if (is.null(parsA)) stop("No params for adult: ", sp)
      
      move_L <- movement_volume(parsA, contact_depth = contact_depth)
      beta0  <- parsA$trans_exp * (move_L / pond_volume_L)
      S_star <- if (is.null(parsA$df_density)) 1 else parsA$df_density
      
      # Z -> I
      G0 <- init_inf_fxn(xlower, xupper, parsA)
      StoI <- parsA$s0 * S_star * beta0 * G0
      
      # I -> Z
      ItoZ <- 3 * (exp(parsA$lambda_int) * exp(y)^parsA$lambda_slope) * 6
      
      # I -> I (within infected)
      G <- sapply(y, function(x) growth_fxn(xlower, xupper, x, parsA))
      L <- loss_sp_fxn(y, parsA, sp)
      ItoI <- G %*% diag(parsA$s0 * (1 - L))
      
      # fill matrices
      F[ii, idx_Z] <- StoI
      F[idx_Z, ii] <- ItoZ
      U[ii, ii]    <- ItoI
      
    } else {
      # -------------------------
      # TADPOLE BLOCK (plus coupling to its adult parent if present)
      # -------------------------
      parent <- tad_parent(sp)
      parsT  <- ipm_dists[[sp]]
      parsA  <- ipm_dists[[parent]]
      if (is.null(parsT)) stop("No params for tadpole: ", sp)
      if (is.null(parsA)) stop("No params for tadpole parent adult: ", parent)
      
      move_L <- movement_volume(parsT, contact_depth = contact_depth)
      beta0  <- parsT$trans_exp * (move_L / pond_volume_L)
      S_star <- if (is.null(parsT$df_density)) 1 else parsT$df_density
      
      # metamorphosis probability (per step)
      meta_prob <- if (!is.null(parsT$meta_successB)) parsT$meta_successB else parsT$meta_successG
      if (is.null(meta_prob) || !is.finite(meta_prob)) meta_prob <- 0
      meta_prob <- max(0, min(0.999, meta_prob))
      meta_eff <- meta_prob
      
      #emigration during metamorphosis 
      emig_meta <- 0
      if (isTRUE(include_emigration)) {
        emig_meta <- parsA$emig_meta
        if (is.null(emig_meta) || !is.finite(emig_meta)) emig_meta <- 0
        emig_meta <- max(0, min(0.999, emig_meta))
      }
      stay_meta <- 1 - emig_meta
      
      tad_surv <- parsT$s0
      
      # Z -> IT (new infections that remain tadpoles this step)
      G0 <- init_inf_fxn(xlower, xupper, parsT)
      StoI_T <- tad_surv * (1 - meta_eff) * S_star * beta0 * G0
      
      # I -> Z (shedding from infected tadpoles that remain tadpoles this step)
      ItoZ_T <- (tad_surv * (1 - meta_eff)) * 3 * (exp(parsT$lambda_int) * exp(y)^parsT$lambda_slope) * 6
      
      # IT -> IT (within infected, conditioned on staying tadpole)
      G <- sapply(y, function(x) growth_fxn(xlower, xupper, x, parsT))
      L <- loss_sp_fxn(y, parsT, sp)
      ItoI_T <- G %*% diag(tad_surv * (1 - meta_eff) * (1 - L))
      
      # fill tadpole block pieces
      F[ii, idx_Z] <- StoI_T
      F[idx_Z, ii] <- ItoZ_T
      U[ii, ii]    <- ItoI_T
      
      # -------------------------
      # Coupling: infected tadpole -> infected adult via metamorphosis (within infected transition)
      # Only if the adult parent is included in species_list
      # -------------------------
      if (parent %in% species_list) {
        jjA <- idx_I[[parent]]
        
        meta_tad <- tad_surv * meta_eff * stay_meta  # fraction that survive + metamorphose + stay
        
        # IT -> IA (metamorphose while infected)
        U[jjA, ii] <- U[jjA, ii] + (tad_to_adult_map %*% (G %*% diag(meta_tad * (1 - L))))
        
        # get infected + metamorphose in same step becomes infected adult (F)
        if (isTRUE(include_infect_then_meta)) {
          StoI_T_meta <- tad_surv * meta_eff * stay_meta * S_star * beta0 * G0
          F[jjA, idx_Z] <- F[jjA, idx_Z] + (tad_to_adult_map %*% StoI_T_meta)
        }
      }
    }
  }
  
  list(F = F, U = U, indexing = list(idx_I = idx_I, idx_Z = idx_Z), disc = disc)
}


compute_R0_stage_community <- function(species_list, ipm_dists,
                                       pond_diameter, lower, upper, bins,
                                       contact_depth = 0.5, avg_depth = 0.5,
                                       tad_to_adult_map = NULL,
                                       include_infect_then_meta = TRUE,
                                       include_emigration = TRUE,
                                       z_decay = 0.3) {
  
  mats <- build_FU_stage_community(
    species_list = species_list,
    ipm_dists = ipm_dists,
    pond_diameter = pond_diameter,
    lower = lower, upper = upper, bins = bins,
    contact_depth = contact_depth, avg_depth = avg_depth,
    tad_to_adult_map = tad_to_adult_map,
    include_infect_then_meta = include_infect_then_meta,
    include_emigration = include_emigration,
    z_decay = z_decay
  )
  
  compute_R0_from_FU(mats$F, mats$U)
}


#---------------------------------------------------------------------------
# Elasticity fxns
#----------------------------------------------------------------------------

perturb_params <- function(ipm_params, sp, par_name, delta, x0) {
  # Make a copy of the full param list
  ipm_pert <- ipm_params
  x0_raw <- exp(x0)
  
  ## ----------------- Intercepts: simple multiplicative ----------------- ##
  if (par_name %in% c("growth_int", "lambda_int", "loss_int", "trans_exp")) {
    
    old_val <- ipm_pert[[sp]][[par_name]]
    new_val <- old_val * (1 + delta)
    ipm_pert[[sp]][[par_name]] <- new_val
    
    ## ----------------- Slopes: anchored so mu(x0) is fixed ---------------- ##
  } else if (par_name == "growth_slope") {
    
    slope0 <- ipm_pert[[sp]]$growth_slope
    int0   <- ipm_pert[[sp]]$growth_int
    slope1 <- slope0 * (1 + delta)
    int1   <- int0 + (slope0 - slope1) * x0
    
    ipm_pert[[sp]]$growth_slope <- slope1
    ipm_pert[[sp]]$growth_int   <- int1
    
  } else if (par_name == "lambda_slope") {
    
    slope0 <- ipm_pert[[sp]]$lambda_slope
    int0   <- ipm_pert[[sp]]$lambda_int
    slope1 <- slope0 * (1 + delta)
    int1   <- int0 + (slope0 - slope1) * x0
    
    ipm_pert[[sp]]$lambda_slope <- slope1
    ipm_pert[[sp]]$lambda_int   <- int1
    
  } else if (par_name == "loss_slope") {
    
    slope0 <- ipm_pert[[sp]]$loss_slope
    int0   <- ipm_pert[[sp]]$loss_int
    slope1 <- slope0 * (1 + delta)
    int1   <- int0 + (slope0 - slope1) * x0
    
    ipm_pert[[sp]]$loss_slope <- slope1
    ipm_pert[[sp]]$loss_int   <- int1
    
    # } else if (par_name == "trans_slope") {
    #   
    #   slope0 <- ipm_pert[[sp]]$trans_slope
    #   int0   <- ipm_pert[[sp]]$trans_int
    #   slope1 <- slope0 * (1 + delta)
    #   int1   <- int0 + (slope0 - slope1) * x0
    #   
    #   ipm_pert[[sp]]$trans_slope <- slope1
    #   ipm_pert[[sp]]$trans_int   <- int1
    
  } else {
    stop("Unknown parameter name: ", par_name)
  }
  
  ipm_pert
}



### Elasticity and sensitvity calculation

calculate_R0_elasticity <- function(ipm_params, species_list,
                                    delta) {
  
  # discretization for x0 (used in anchored slope perturbation)
  disc <- set_discretized_values(min_size = -15, max_size = 18, bins = 150)
  x0   <- mean(disc$y)
  
  # the parameters you want elasticities for
  species_list <- names(ipm_params)
  par_names <- c("growth_int", "lambda_int", "loss_int", "trans_exp",
                 "growth_slope", "lambda_slope", "loss_slope")
  
  res_list <- list()
  
  for (sp in species_list) {
    # base R0 for this species
    R0_base <- compute_community_R0(sp, ipm_params,
                                    lower = -15, upper = 18,
                                    bins = 150, pond_diameter = 22)
    
    for (par in par_names) {
      # base parameter value p0
      is_slope <-grepl("_slope", par)
      p0 <- ipm_params[[sp]][[par]]
      
      den_p <- if(is_slope) abs(p0) else p0
      
      # upper and lower perturbed parameter sets
      ipm_upper <- perturb_params(ipm_params, sp, par, +delta, x0)
      ipm_lower <- perturb_params(ipm_params, sp, par, -delta, x0)
      
      R0_up  <- compute_community_R0(sp, ipm_upper,
                                     lower = -15, upper = 18,
                                     bins = 150, pond_diameter = 22)
      R0_low <- compute_community_R0(sp, ipm_lower,
                                     lower = -15, upper = 18,
                                     bins = 150, pond_diameter = 22)
      
      # central-difference sensitivity
      # derivative wrt the *absolute* parameter value:
      # dR0/dp ≈ (R0_up - R0_low) / (2 * delta * p0)
      sens <- (R0_up - R0_low) / (2 * delta * den_p)
      
      # elasticity = (dR0/dp) * (p0 / R0_base)
      elas <- sens * (den_p / R0_base)
      
      res_list[[length(res_list) + 1]] <- data.frame(
        species     = sp,
        parameter   = par,
        p0          = p0,
        R0_base     = R0_base,
        sensitivity = sens,
        elasticity  = elas
      )
    }
  }
  
  do.call(rbind, res_list)
}






################################################################################
