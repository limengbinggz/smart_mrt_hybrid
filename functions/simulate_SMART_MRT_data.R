#######################################################
#' Simulate data from a SMART-MRT hybrid design study.
#######################################################

library(data.table)
# library(foreach)
# library(doMC)


expit = function(x){
  return(1 / (1 + exp(-x)))
}


### simulate data --------------------
#' Simulation functions
#' @description Generate a dataset representing a possible hybrid SMART-MRT study 
#' with hypothesized parameters
#' @param N number of subjects to be generated
#' @param T number of time points (e.g., days or weeks) in the observation process
#' @param T_SMART the time point at which re-randomization of ADI options occur
#' @param p_responder proportion of responders, which we assume is a constant across ADIs
#' @param eta_A a length-4 numeric coefficient vector representing the effect of
#'  A_{t-1}, state_t, Z_1, Z_2 on the treatment probability at time t
#' @param eta_state a length-3 numeric coefficient vector representing the effect of
#'  A_{t-1}, Z_1, Z_2 on the state at time t
#' @param beta_interaction a length-6 numeric coefficient vector representing the interaction
#'  between treatment A_t and ADI options: intercept, Z_1, Z_2, Z_1 * Z_2, state_t, state_t * Z_1
#' @param gamma_Z a length-5 numeric coefficient vector representing the interaction
#'  between treatment A_t and ADI options: Z_1, Z_2, Z_1 * Z_2, state_t * Z_1
#'  responder status
#' @param error_cor autocorrelation between consecutive error terms in proximal outcomes
#' @param error_sd standard deviation of the residual errors in proximal outcomes
#' @param mc.cores number of cores to use in data simulation in parallel
#' @import parallel  
#' @import doParallel  
#' @example 
#' N = 100
#' T = 112
#' T_SMART = 28
#' p_SMART = rep(0.5, length(T_SMART) + 1)
#' p_responder = 0.5
#' set.seed(50)
#' sim_data <- generate_random_data_hybrid(N)
generate_random_data_hybrid <- function(N, T = 112, T_SMART = 28, p_SMART = rep(0.5, length(T_SMART) + 1), 
                                        pt_setting = "Aconstant",
                                        responder_setting = c("none", "constant", "Hdependent"), 
                                        eta_A = c(-0.8, 0.8, 0.1, 0.1),
                                        eta_state = c(0.2, 0.05, 0.08),
                                        beta_interaction = c(0.3, -0.3, 0.3, -0.3, 0.2, 0.2),
                                        gamma_Z = c(0.1, -0.1, -0.1, 0.1, 0.08),
                                        error_cor = sqrt(0.2), error_sd = 1,
                                        seed = 1, mc.cores = 4L){ 

  # SZ1_centered <- function(eta_state, A_t_previous, CZ2){
  #   # P(Xt = 1 | Z1 = 1)
  #   p_Xt1_Z11 <- expit(eta_state %*% c(A_t_previous, 1, CZ2))
  #   # P(Xt = 1 | Z1 = -1)
  #   p_Xt1_Z1minus1 <- expit(eta_state %*% c(A_t_previous, -1, CZ2))
  #   0.5 * (p_Xt1_Z11 - 1 + p_Xt1_Z11 - p_Xt1_Z1minus1 + 1 - p_Xt1_Z1minus1)
  # }
  SZ1_centered <- function(eta_state, A_t_previous, CZ2){
    # P(Xt = 1 | Z1 = 1)
    p_Xt1_Z11 <- expit(eta_state %*% c(A_t_previous, 1, CZ2))
    # P(Xt = 1 | Z1 = -1)
    p_Xt1_Z1minus1 <- expit(eta_state %*% c(A_t_previous, -1, CZ2))
    0.5 * (p_Xt1_Z11 - p_Xt1_Z1minus1)
  }
  
  Z_all_numeric <- matrix(c(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)), ncol = 2, byrow = T)
  Z_all <- sapply(list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1), c(1, 0), c(-1, 0)), function(x) paste0(x, collapse = ""))
  Y_counterfactual_name_all <- c("Y_counterfactual_Z11", "Y_counterfactual_Z1m1", "Y_counterfactual_Zm11", "Y_counterfactual_Zm1m1",
                                 "Y_counterfactual_Z11", "Y_counterfactual_Zm11")
  # for responders only
  Z_all_numeric_R1 <- matrix(c(c(1, 0), c(-1, 0)), ncol = 2, byrow = T)
  Z_all_R1 <- sapply(list(c(1, 0), c(-1, 0)), function(x) paste0(x, collapse = ""))
  Y_counterfactual_name_all_R1 <- c("Y_counterfactual_Z11", "Y_counterfactual_Zm11")
  
  simulate_data_one <- function(i){
    # stage-1 ADI option
    Z_1 <- sample(c(1, -1), 1, prob = c(p_SMART[1], 1 - p_SMART[1]))
    Z_2 <- 0
    R_c <- 0
    data_i <- data.frame(id = i, time = 1:T, 
                         Z = 0, Z1 = Z_1, Z2 = 0,
                         responder = 0, A = NA, state = NA, Y = NA, trt_prob_t = NA, state_prob_t = NA,
                         residual_error = NA, Y_counterfactual_A = NA)
    
    # randomly generate the initial state
    # state_t <- sample(c(1, -1), 1)
    A_t_previous <- 0
    A_t_previous_c <- A_t_previous
    error_t_previous <- 0
    R <- 0
    X1_previous <- 0
    for (t in 1:T) {
      if (t == T_SMART){
        # responder status
        if (responder_setting == "none"){
          p_responder <- 0
        } else if (responder_setting == "constant"){
          p_responder <- 0.5
        } else if (responder_setting == "Z1dependent"){
          p_responder <- 0.6 * (Z_1 == 1) + 0.45 * (Z_1 == -1)  #expit(0.5 * Z_1)
        } else if (responder_setting == "Hdependent"){
          p_responder <- expit(-0.62 + data_i[1, 'state_t_c'] + data_i[t-1, 'A_t_c'] + 0.5 * Z_1)
        }
        R <- rbinom(1, 1, p_responder)
        # centered responder status
        R_c <- R - p_responder #(p_responder-1)
        # stage-2 ADI option
        # Z_2 <- (1 - R) * sample(c(1, -1), 1, prob = c(p_SMART[2], 1 - p_SMART[2])) #R * Z_1 + 
        Z_2 <- R * Z_1 + (1 - R) * sample(c(1, -1), 1, prob = c(p_SMART[2], 1 - p_SMART[2])) #
      }
      
      C_t <- (t >= T_SMART) * 1
      CZ2 <- C_t * (1 - R) * Z_2
      # calculate the probability of state being +1
      state_prob_t <- expit(eta_state %*% c(A_t_previous, Z_1, CZ2))#t/T + 
      # generate state at time t
      state_t <- sample(c(1, -1)*2, 1, prob = c(state_prob_t, 1 - state_prob_t))
      # state_t <- sample(c(2, -2), 1, prob = c(state_prob_t, 1 - state_prob_t))
      # state_t <- sample(c(1, 0), 1, prob = c(state_prob_t, 1 - state_prob_t))
      # state_t <- sample(c(1, 2), 1, prob = c(state_prob_t, 1 - state_prob_t))
      # centered state
      # state_t_c <- state_t - 2 * state_prob_t + 1
      state_t_c <- state_t + (- 4 * state_prob_t + 2)
      # state_t_c <- state_t - state_prob_t
      # state_t_c <- state_t - 2 + state_prob_t
      
      ## generate another time-varying covariate X1
      X1_prob <- expit(c(0.1, 0.2) %*% c(Z_1, CZ2))
      X1 <- sample(c(1, -1)*2, 1, prob = c(X1_prob, 1 - X1_prob))
      # X1 <- -0.2 * X1_previous + Z_1 + CZ2 + rnorm(1, 0, 0.5)
      # X1_previous <- X1
      
      # calculate the probability of action being 1
      if (pt_setting == "Aconstant"){
        trt_prob_t <- 0.5
      } else if (pt_setting == "ADependOnZ1"){
        trt_prob_t <- c(0.6, 0.4) %*% c(Z_1==1, Z_1==-1)
      } else if (pt_setting == "ADependOnZ1Z2"){
        trt_prob_t <- c(0.6, 0.4, -0.2, 0.2) %*% c(Z_1==1, Z_1==-1, CZ2==1, CZ2==-1)
      } else if (pt_setting == "ADependOnX1"){
        trt_prob_t <- expit(eta_A %*% c(A_t_previous, X1, Z_1, CZ2))
      } else if (pt_setting == "ADependOnX1X2"){
        trt_prob_t <- expit(eta_A %*% c(A_t_previous, state_t, Z_1, CZ2))
      }
      # generate action at time t
      A_t <- rbinom(1, 1, trt_prob_t)
      # counter factual At
      A_t_counterfactual <- 1 - A_t
      # centered action at time t
      A_t_c <- A_t - trt_prob_t
      # generate residual error
      error_t <- rnorm(1, mean = 0, sd = error_sd)
      error_total <- error_t + error_cor * error_t_previous
      SZ1_c <- state_t * Z_1 - SZ1_centered(eta_state, A_t_previous, CZ2)
      
      # if (t == 1){
      #   cat("\nstate_t_c =", state_t_c)
      #   cat("\nA_t_previous_c =", A_t_previous_c)
      #   cat("\nA_t_c =", A_t_c)
      #   cat("\nc(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1) =", c(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1))
      #   cat("\nc(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c) =", c(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c))
      #   cat("\nerror_total =", error_total)
      #   cat("\n0.5 * state_t_c - 0.1 * A_t_previous_c =", 0.5 * state_t_c - 0.1 * A_t_previous_c)
      #   cat("\nA_t_c * (beta_interaction %*% c(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1)) =", A_t_c * (beta_interaction %*% c(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1)))
      #   cat("\ngamma_Z %*% c(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c) =", gamma_Z %*% c(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c))
      # }
      
      # generate response Y_{t+1}
      Y_t_next <- 0.5 * state_t_c + 
        - 0.1 * A_t_previous_c +
        A_t_c * (beta_interaction %*% c(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1)) + #A_t_c
        gamma_Z %*% c(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c) + #SZ1_c
        error_total
      
      # generate counterfactual response Y_{t+1} if the other A_t is taken
      Y_counterfactual_A <- 0.5 * state_t_c + 
        - 0.1 * A_t_previous_c +
        (A_t_counterfactual - trt_prob_t) * (beta_interaction %*% c(1, Z_1, CZ2, CZ2 * Z_1, state_t_c, state_t_c * Z_1)) + #SZ1_c
        gamma_Z %*% c(Z_1, CZ2, CZ2 * Z_1, state_t_c * Z_1, C_t * R_c) + #SZ1_c
        error_total

      # record variables to data
      data_i[t, c('A', 'A_t_c', 'A_t_previous_c', 'state', 'state_t_c', 'trt_prob_t', 'state_prob_t', 'SZ1_c', 'X1')] <- 
        c(A_t, A_t_c, A_t_previous_c, state_t, state_t_c, trt_prob_t, state_prob_t, state_t_c * Z_1, X1)
      data_i[t, 'Y'] <- Y_t_next
      data_i[t, 'Y_counterfactual_A'] <- Y_counterfactual_A
      data_i[t, 'residual_error'] <- error_total
      
      # update the variables at t-1
      error_t_previous <- error_total
      A_t_previous <- A_t
      A_t_previous_c <- A_t_c
    }
    data_i <- data.table(data_i)
    data_i$Z <- c(rep(Z_1, T_SMART-1), rep(Z_2, T - T_SMART+1))
    data_i$responder <- R
    data_i$p_responder <- p_responder
    data_i$Z2 <- Z_2
    data_i[, Ct := (time >= T_SMART)]
    data_i[, CZ2 := (time >= T_SMART) * Z_2 * (1 - responder)]


    # compute all counterfactuals
    for (cf_index in 1:nrow(Z_all_numeric)) {
      Z1_cf <- Z_all_numeric[cf_index, 1]; Z2_cf <- Z_all_numeric[cf_index, 2]
      Z_cf_index <- which(Z_all == paste0(c(Z1_cf, Z2_cf), collapse = ""))
      data_i[, paste(Y_counterfactual_name_all[Z_cf_index]) := 
               0.5 * state_t_c + 
               - 0.1 * A_t_previous_c +
               A_t_c * (as.matrix(data_i[, .(1, Z1_cf, Ct * Z2_cf, Ct * Z2_cf * Z1_cf, state_t_c, state_t_c * Z1_cf)]) %*% beta_interaction) + 
               as.matrix(data_i[, .(Z1_cf, Ct * Z2_cf, Ct * Z2_cf * Z1_cf, state_t_c * Z1_cf, Ct * R_c)]) %*% gamma_Z + 
               residual_error]
    }
    # For a nonresponder, we have all counterfactuals. So no need to do anything here.
    # For a responder, we only have the counterfactual Y(Z1, Z_R=0, Z_NR=0)
    for (cf_index in 1:nrow(Z_all_numeric_R1)) {
      Z1_cf <- Z_all_numeric_R1[cf_index, 1]; Z2_cf <- Z_all_numeric_R1[cf_index, 2]
      Z_cf_index <- which(Z_all_R1 == paste0(c(Z1_cf, Z2_cf), collapse = ""))
      data_i[(time >= T_SMART) & (responder == 1), paste(Y_counterfactual_name_all_R1[Z_cf_index]) := 
               0.5 * state_t_c + 
               - 0.1 * A_t_previous_c +
               A_t_c * (as.matrix(data_i[(time >= T_SMART) & (responder == 1), .(1, Z1_cf, Ct * Z2_cf, Ct * Z2_cf * Z1_cf, state_t_c, state_t_c * Z1_cf)]) %*% beta_interaction) + 
               as.matrix(data_i[(time >= T_SMART) & (responder == 1), .(Z1_cf, Ct * Z2_cf, Ct * Z2_cf * Z1_cf, state_t_c * Z1_cf, Ct * R_c)]) %*% gamma_Z + 
               residual_error]
    }
    data_i[(time >= T_SMART) & (responder == 1), Y_counterfactual_Z1m1 := Y_counterfactual_Z11]
    data_i[(time >= T_SMART) & (responder == 1), Y_counterfactual_Zm1m1 := Y_counterfactual_Zm11]
    
    # # if this person is a responder, we only need to generate the counterfactuals in the first stage
    # # the counterfactuals in the second stage are assumed to be equal to the observed
    # if (R == 1){
    #   # for the second stage, replace the counterfactuals with the observed values
    #   data_i[time >= T_SMART, paste(Y_counterfactual_name_all[2]) := get(Y_counterfactual_name_all[1])]
    #   data_i[time >= T_SMART, paste(Y_counterfactual_name_all[3]) := get(Y_counterfactual_name_all[4])]
    # }
    
    return(data_i)
    # sim_data[[i]] <- data_i
  }
  # cl <- makeCluster(mc.cores)
  # registerDoParallel(cl)
  # sim_data <- foreach(i=1:N) %dopar% {
  #   set.seed(seed+i)
  #   simulate_data_one(i)
  # }
  # 
  
  # cl <- makeCluster(mc.cores)
  # clusterSetRNGStream(cl, seed)  ## set seed
  # sim_data <- parSapply(cl, 1:N, simulate_data_one)
  # stopCluster(cl)
  # 
  # cl <- makeCluster(mc.cores)
  # doParallel::registerDoParallel(cl)
  # # registerDoMC(mc.cores)  # Register a parallel backend with two workers
  # sim_data <- foreach(i = 1:N) %dopar% {
  #   set.seed(seed*i)
  #   simulate_data_one(i)
  # }
  # stopCluster(cl)
  # # # sim_data <- do.call(rbind, sim_data)
  # # unregister <- function() {
  # #   env <- foreach:::.foreachGlobals
  # #   rm(list=ls(name=env), pos=env)
  # # }
  # # unregister()
  
  set.seed(seed)
  # RNGkind("L'Ecuyer-CMRG")
  # sim_data <- parallel::mclapply(1:N, simulate_data_one, mc.cores = mc.cores)
  sim_data <-lapply(1:N, simulate_data_one)
  sim_data <- do.call(rbind, sim_data)
  sim_data <- data.table(sim_data)
  sim_data$stage <- (sim_data$time >= T_SMART) + 1
  return(sim_data)
}





