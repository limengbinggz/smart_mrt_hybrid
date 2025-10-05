#######################################################
#' Simulate data and fit with WCLS.
#######################################################
if(Sys.info()["sysname"] %in% c("Darwin")){
  curr_dir <- "/Users/mengbing/Dropbox (University of Michigan)/from_box/research/GSRA_walter"
  setwd(curr_dir)
} else if (substr(Sys.info()["nodename"],1, 2) == "cn" | Sys.info()["user"] == "mengbing"){ # biostat cluster
  curr_dir <- "/home/mengbing/research/GSRA_walter"
  setwd(curr_dir)
} else{ # greatlakes
  curr_dir <- "/home/mengbing/research/GSRA_walter"
  setwd(curr_dir)
}
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
source("functions/simulate_SMART_MRT_data.R")
source("functions/utils.R")
library(geepack)
# library(wgeesel)
library(rootSolve) #multiroot
library(parallel)


### Simulate data -------------------------------------------
seed <- as.integer(args[1])
pt_setting <- as.character(args[2])
responder_setting = as.character(args[3]) 
N <- as.integer(args[4])
# seed = 83
# pt_setting = "ADependOnZ1Z2" "ADependOnX1" "ADependOnZ1" "Aconstant"
# responder_setting = "Hdependent" "Z1dependent" "constant" "none"
# N = 200
cat('\n====================== pt_setting:', pt_setting, '======================\n')
cat('\n====================== responder_setting:', responder_setting, '======================\n')
cat('\n Seed', seed, '\n')
T = 50
T_SMART = 14
p_SMART = rep(0.5, length(T_SMART) + 1)

set.seed(seed)
eta_A = c(-0.2, 0.5, 0.1, 0.2)
eta_state = c(-1, 0.1, 0.2) #c(-1, 0.1, 0.2)
beta_interaction = c(0.4, -0.3, 0.2, -0.1, 0.4, 0.2) #c(0.4, -0.3, 0.3, -0.3, 0.2, 0.2)
gamma_Z = c(0.2, -0.1, -0.1, 0.2, 0.2)


contrasts_geeglm <- function(gee_coef, gee_cov, contrast){
  contrast_est = contrast %*% gee_coef 
  contrast_se = t(contrast) %*% gee_cov %*% contrast
  output = data.table(Estimate = c(contrast_est),
                      SE = sqrt(diag(contrast_se)))
  output[, lci := Estimate - 1.96*SE]
  output[, uci := Estimate + 1.96*SE]
  return(output)
}



run_ee <- function(data, p_tilde){
  
  ## next we need to augment data to accommodate for SMART design
  rows_to_replicate <- data[responder == 1, ]
  rows_to_replicate[, weight_SMART := 1 / p_SMART[1]]
  rows_to_replicate_observed <- rows_to_replicate
  rows_to_replicate_observed$replicant <- 1L
  rows_to_replicate_observed[, Z2 := Z1]
  rows_to_replicate_pseudo <- rows_to_replicate
  rows_to_replicate_pseudo$replicant <- 2L
  rows_to_replicate_pseudo[, Z2 := Z1 * (-1)]
  rows_not_to_replicate <- data[responder == 0, ]
  rows_not_to_replicate$replicant <- 1
  rows_not_to_replicate[, weight_SMART := 1 / (p_SMART[1] * p_SMART[2])]
  data_augmented <- rbind(rows_not_to_replicate,
                          rows_to_replicate_observed,
                          rows_to_replicate_pseudo)
  rm(rows_to_replicate, rows_to_replicate_observed, rows_to_replicate_pseudo, rows_not_to_replicate)
  data_augmented <- data_augmented[order(id, time, replicant)]
  
  # calculate weight of the MRT actions
  # data_augmented <- data
  # data_augmented[, weight_SMART := 1]
  # data_augmented[, p_tilde := trt_prob_t]
  # data_augmented[time < T_SMART, p_tilde := mean(A), by = .(time)]
  # data_augmented[time >= T_SMART, p_tilde := mean(A), by = .(time)]
  # data_augmented[, p_tilde := mean(A), by = .(Z1)]#rho
  data_augmented[, weight_MRT := (p_tilde / trt_prob_t)^A * ((1 - p_tilde) / (1 - trt_prob_t))^(1-A)]
  data_augmented[, weight := weight_MRT * weight_SMART]
  data_augmented[, A_tilde_c := A - p_tilde]
  
  # compute the centered control variables by Z
  data_augmented[, state_mean := sum(weight_SMART*state) / sum(weight_SMART), by = .(time, Z1, Z2)]
  data_augmented[, state_c := state - state_mean]
  data_augmented[, stateZ1_mean := sum(weight_SMART*state*Z1) / sum(weight_SMART), by  = .(time, Z1, Z2)]
  data_augmented[, stateZ1_c := state*Z1 - stateZ1_mean]
  data_augmented[, Z2_after := (time >= T_SMART)*Z2]
  data_augmented[, C_t_star := (time >= T_SMART)*1]
  data_augmented[, before := (time < T_SMART)*1]
  data_augmented[, after := (time >= T_SMART)*1]
  
  
  ## run GEE
  model1 <- geeglm(Y ~ -1 + state_c + stateZ1_c +
                  A_tilde_c + A_tilde_c:Z1 +
                  A_tilde_c:Z2:after + A_tilde_c:Z1:Z2:after +
                  before + before:Z1 +
                  after + after:(Z1*Z2),
                data = data_augmented, weights = weight, id = id)
  select_var1 <- c("A_tilde_c", "A_tilde_c:Z1",
                  "A_tilde_c:Z2:after", "A_tilde_c:Z1:Z2:after",
                  "before", "Z1:before",
                  "after", "Z1:after", "Z2:after", "Z1:Z2:after")
  select_idx <- match(select_var1, names(model1$coefficients))
  solution <- model1$coefficients[select_idx]
  cov_mat <- model1$geese$vbeta[select_idx,]
  cov_mat <- cov_mat[, select_idx]
  names(solution) <- gsub("A_tilde_c", "(At - pt_tilde)", names(solution))
  # solution
  
  
  
  
  ### solve for the effect of Z under MRT randomization probability
  # # get the predicted Y under the alternation MRT randomization probability
  # data_augmented[, Y_predicted := model1$linear.predictors]
  # ## run GEE
  # model2 <- geeglm(Y_predicted ~ 1 + state_c + stateZ1_c +
  #                    Z1*Z2_after,
  #                  data = data_augmented, weights = weight_SMART, id = id)

  # get the variances of the predicted Y
  design_mat_model1 <- model.matrix( ~ -1 +
                                      A_tilde_c + A_tilde_c:Z1 +
                                      A_tilde_c:Z2:after + A_tilde_c:Z1:Z2:after +
                                      before + before:Z1 +
                                      after + after:(Z1*Z2), data = data_augmented)
  design_mat_model1 <- design_mat_model1[, select_var1]

  data_augmented[, Y_predicted := design_mat_model1 %*% solution]
  ## run GEE
  model2 <- geeglm(Y ~ state_c + stateZ1_c + Z1*Z2_after,
                   data = data_augmented, weights = weight_SMART, id = id, corstr = "independence")
  select_var2 <- c("(Intercept)", "Z1", "Z2_after", "Z1:Z2_after")
  select_idx <- match(select_var2, names(model2$coefficients))
  solution_mrt_prob <- model2$coefficients[select_idx]
  cov_mat_mrt_prob <- model2$geese$vbeta[select_idx,]
  cov_mat_mrt_prob <- cov_mat_mrt_prob[, select_idx]
  names(solution_mrt_prob) <- gsub("A_tilde_c:", "(At - pt_tilde):", names(solution_mrt_prob))
  
  return(list(solution = solution, cov_mat = cov_mat,
              solution_mrt_prob = solution_mrt_prob, cov_mat_mrt_prob = cov_mat_mrt_prob,
              data_augmented = data_augmented))
}




calculate_marginal_effects <- function(solution, cov_mat, solution_mrt_prob, cov_mat_mrt_prob, data_augmented){
  ### calculate the marginal effects ----------------------------
  p_beta <- length(grep("At - pt_tilde", names(solution)))
  p_gamma <- length(solution) - p_beta
  p_beta_true <- 4
  p_gamma_true <- 4
  # true_coef <- c(beta_interaction[1:p_beta], gamma_Z[1:p_gamma])
  true_coef <- c(beta_interaction[1:p_beta_true], 0, gamma_Z[1:(p_gamma_true-1)])
  
  coef_names <- c(paste0("\\beta_", 1:p_beta-1), paste0("\\gamma_", 1:p_gamma-1))
  
  Z_all_numeric <- matrix(c(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)), ncol = 2, byrow = T)
  Z_all <- sapply(list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)), function(x) paste0(x, collapse = ""))
  Y_counterfactual_name_all <- c("Y_counterfactual_Z11", "Y_counterfactual_Z1m1", "Y_counterfactual_Zm11", "Y_counterfactual_Zm1m1")
  
  # randomization probabilities in 2 stages  
  if (pt_setting == "Aconstant"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.5, 0.5) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.5, 0.5) %*% c(Z2==1, Z2==-1)
    }
  } else if (pt_setting == "ADependOnZ1"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1)
    }
  } else if (pt_setting == "ADependOnZ1Z2"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.6, 0.4, -0.2, 0.2) %*% c(Z1==1, Z1==-1, (1-R)*Z2==1, (1-R)*Z2==-1)
    }
  }
  
  p_responder_empirical <- mean(unique(data[, .(id, responder)])$responder)
  if (responder_setting == "none"){
    p_responder <- c(0, 0)
  } else if (responder_setting == "constant"){
    p_responder <- c(0.5, 0.5)
  } else if (responder_setting == "Z1dependent"){
    p_responder <- c(0.45, 0.6)  #expit(0.5 * Z_1)
  } else if (responder_setting == "Hdependent"){
    p_responder <- c(mean(unique(data[Z1 == -1, .(id, responder)])$responder), 
                     mean(unique(data[Z1 == 1, .(id, responder)])$responder))  
  }
  dat_marginal_effects <- c()
  
  
  
  ## 1. comparing At = 1 vs At = 0 given Z ----------------
  
  ## first look at stage 1
  before <- 1; after <- 0
  ## Effect of At for a given Z = (1, 0) at t < t^*
  Z_ <- c(1, 0)
  contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  # se <- output$SE
  # lci <- output$lci
  # uci <- output$uci
  contrast_true <- c(1, c(Z_, prod(Z_)), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  # truth_across_Z_stage1 <- c(truth_across_Z_stage1, truth)
  mc_estimate <- data[(Z1 == Z_[1] & time < T_SMART), mean((Y - Y_counterfactual_A) * (2*A-1))]
  label <- paste0('t < t^*: E[Y(A = 1, Z_1 = ', Z_[1], ")] - E[Y(A = 0, Z_1 = ", Z_[1], ")]")
  coef_labels <- paste_coef(contrast, coef_names)
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  ## Effect of At for a given Z = (-1, 0) at t < t^*
  Z_ <- c(-1, 0)
  contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
  # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  contrast_true <- c(1, c(Z_, prod(Z_)), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  # truth_across_Z_stage1 <- c(truth_across_Z_stage1, truth)
  mc_estimate <- data[(Z1 == Z_[1] & time < T_SMART), mean((Y - Y_counterfactual_A) * (2*A-1))]
  label <- paste0('t < t^*: E[Y(A = 1, Z_1 = ', Z_[1], ")] - E[Y(A = 0, Z_1 = ", Z_[1], ")]")
  coef_labels <- paste_coef(contrast, coef_names)
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  ## first look at stage 2
  before <- 0; after <- 1
  ## Effect of At for a given Z2 = 1 at t >= t^*
  # store the effect given each Z
  for (Z_index in 1:nrow(Z_all_numeric)) {
    Z_ <- Z_all_numeric[Z_index, ]
    contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
    # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
    contrast <- c(contrast_beta, rep(0, p_gamma))
    output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
    estimate <- output$Estimate
    
    # truth <- contrast %*% true_coef
    # responders only have Z1 effects
    truth_r <- c(1, Z_[1], 0, 0) %*% true_coef[1:p_beta_true]
    # nonresponders only have Z1 and Z2 effects
    truth_nr <- c(1, c(Z_, prod(Z_))) %*% true_coef[1:p_beta_true] 
    truth <- p_responder[(Z_[1]+1)/2+1] * truth_r + (1 - p_responder[(Z_[1]+1)/2+1]) * truth_nr
    mc_estimate <- data_augmented[(Z1 == Z_[1] & Z2 == Z_[2] & time >= T_SMART), 
                                  sum(weight_SMART * (Y - Y_counterfactual_A) * (2*A-1)) / sum(weight_SMART)]
    label <- paste0('t >= t^*: E[Y(A = 1, Z = (', Z_[1], ", ",  Z_[2], ")] - E[Y(A = 0, Z = (", Z_[1], ", ",  Z_[2], ")]")
    coef_labels <- paste_coef(contrast, coef_names)
    dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  }
  
  
  ## Effect of At averaging over Z
  # stage 1
  before <- 1; after <- 0
  Z_ <- c(0, 0)
  contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
  # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  mc_estimate <- mean(data[(time < T_SMART), (Y - Y_counterfactual_A)*(2*A - 1)])
  contrast_true <- c(1, rep(0, p_beta_true-1), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  label <- paste0("t < t^*: E[Y(A = 1, Z1)] - E[Y(A = 0, Z1)]")
  coef_labels <- paste_coef(contrast, coef_names)
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  
  # stage 2
  before <- 0; after <- 1
  Z_ <- c(0, 0)
  contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
  # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  # estimate <- contrast %*% solution
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  contrast_true <- c(1, rep(0, p_beta_true-1), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  mc_estimate <- data_augmented[(time >= T_SMART), sum(weight_SMART * (Y - Y_counterfactual_A) * (2*A-1)) / sum(weight_SMART)]
  label <- paste0("t >= t^*: E[Y(A = 1, Z1, Z2)] - E[Y(A = 0, Z1, Z2)]")
  coef_labels <- paste_coef(contrast, coef_names)
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  
  
  
  ## 2. comparing Z's at specific At ----------------
  At_all <- 0:1
  At <- At_all[1]
  for (At in At_all) {
    
    # first look at stage 1 
    before <- 1; after <- 0
    pt_all <- unique(data[(time < T_SMART), .(trt_prob_t, Z1)])
    mc_estimate_all <- data[(time < T_SMART) & A == At, mean(Y), by = Z1]
    
    Z_ <- c(1, 0)
    pt_At <- A_stage1_prob_calculate(Z_[1])
    contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
    contrast_gamma <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
    contrast1 <- c(contrast_beta*(At-p_tilde), contrast_gamma)
    contrast_beta_true <- c(1, c(Z_, prod(Z_)))
    contrast_gamma_true <- c(1, Z_, prod(Z_))
    truth1 <- (At-pt_At) * contrast_beta_true %*% true_coef[1:p_beta_true] + contrast_gamma_true %*% true_coef[1:p_gamma_true+p_beta_true]
    
    Z_ <- c(-1, 0)
    pt_At <- A_stage1_prob_calculate(Z_[1])
    contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
    contrast_gamma <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
    contrast2 <- c(contrast_beta*(At-p_tilde), contrast_gamma)
    contrast_beta_true <- c(1, c(Z_, prod(Z_)))
    contrast_gamma_true <- c(1, Z_, prod(Z_))
    truth2 <- (At-pt_At) * contrast_beta_true %*% true_coef[1:p_beta_true] + contrast_gamma_true %*% true_coef[1:p_gamma_true+p_beta_true]
    
    output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast1 - contrast2)
    estimate <- output$Estimate
    mc_estimate <- mc_estimate_all[Z1 == 1, V1] - mc_estimate_all[Z1 == -1, V1]
    truth <- truth1 - truth2
    label <- paste0("t < t^*: E[Y(A = ", At, ", Z_1 = 1) - Y(A = ", At, ", Z_1 = -1)]")
    dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - truth))
    
    
    
    
    ### then look at stage 2 
    before <- 0; after <- 1
    pt_all <- unique(data[(time > T_SMART), .(trt_prob_t, Z1, Z2, responder)])
    mc_estimate_all <- data_augmented[(time >= T_SMART) & A == At, sum(weight_SMART*Y) / sum(weight_SMART), by = .(Z1, Z2)]
    
    all_combn_index <- combn(1:4, 2)
    for (combn_index in 1:ncol(all_combn_index)) {
      Z_ <- Z_all_numeric[all_combn_index[,combn_index][1], ]
      pt_At <- c(A_stage2_prob_calculate(Z_[1], 0, R = 1),
                 A_stage2_prob_calculate(Z_[1], Z_[2], R = 0))
      contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
      contrast_gamma <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
      contrast1 <- c(contrast_beta*(At-p_tilde), contrast_gamma)
      # responders only have Z1 effects
      truth_r <- (At-pt_At[1]) * c(1, Z_[1], 0, 0) %*% true_coef[1:p_beta_true] + c(1, Z_[1], 0, 0) %*% true_coef[1:p_gamma_true+p_beta_true]
      # nonresponders have Z1 and Z2 effects
      contrast_beta_true <- c(1, c(Z_, prod(Z_)))
      contrast_gamma_true <- c(1, Z_, prod(Z_))
      truth_nr <- (At-pt_At[2]) * c(1, c(Z_, prod(Z_))) %*% true_coef[1:p_beta_true] + contrast_gamma_true %*% true_coef[1:p_gamma_true+p_beta_true]
      truth1 <- p_responder[(Z_[1]+1)/2+1] * truth_r + (1 - p_responder[(Z_[1]+1)/2+1]) * truth_nr
      
      Z_prime <- Z_all_numeric[all_combn_index[,combn_index][2], ]
      pt_At <- c(A_stage2_prob_calculate(Z_prime[1], 0, R = 1),
                 A_stage2_prob_calculate(Z_prime[1], Z_prime[2], R = 0))
      contrast_beta <- c(1, Z_prime[1], after*Z_prime[2], after*prod(Z_prime))
      contrast_gamma <- c(before, before*Z_prime[1], after, Z_prime[1]*after, after*Z_prime[2], after*prod(Z_prime))
      contrast2 <- c(contrast_beta*(At-p_tilde), contrast_gamma)
      # responders only have Z1 effects
      truth_r <- (At-pt_At[1]) * c(1, Z_prime[1], 0, 0) %*% true_coef[1:p_beta_true] + c(1, Z_prime[1], 0, 0) %*% true_coef[1:p_gamma_true+p_beta_true]
      # nonresponders have Z1 and Z2 effects
      contrast_beta_true <- c(1, c(Z_prime, prod(Z_prime)))
      contrast_gamma_true <- c(1, Z_prime, prod(Z_prime))
      truth_nr <- (At-pt_At[2]) * c(1, c(Z_prime, prod(Z_prime))) %*% true_coef[1:p_beta_true] + contrast_gamma_true %*% true_coef[1:p_gamma_true+p_beta_true]
      truth2 <- p_responder[(Z_prime[1]+1)/2+1] * truth_r + (1 - p_responder[(Z_prime[1]+1)/2+1]) * truth_nr
      
      
      output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast1 - contrast2)
      estimate <- output$Estimate
      mc_estimate <- mc_estimate_all[Z1 == Z_[1] & Z2 == Z_[2], V1] - mc_estimate_all[Z1 == Z_prime[1] & Z2 == Z_prime[2], V1]
      truth <- truth1 - truth2
      label <- paste0("t >= t^*: E[Y(A = ", At, ", Z = (", Z_[1], ", ", Z_[2], ")) - Y(A = ", At, ", Z = (", Z_prime[1], ", ", Z_prime[2], "))]")
      dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - truth))
    }
  }
  ## END fixed At
  
  
    
    
    
  ## 3. comparing Z's averaging over At ----------------
  solution_combined <- c(solution[1:p_beta], solution_mrt_prob)
  coef_names <- paste0("\\gamma_", 1:length(solution_mrt_prob)-1)
  
  # first look at stage 1 
  mc_estimate_all <- data_augmented[(time < T_SMART), sum(weight_SMART * Y) / sum(weight_SMART), by = .(Z1)]
  before <- 1; after <- 0
  Z_ <- c(1, 0)
  contrast1 <- c(1, Z_, prod(Z_))
  contrast_true <- c(rep(0, p_beta_true), c(1, Z_, prod(Z_)))
  truth1 <- contrast_true %*% true_coef
  
  Z_ <- c(-1, 0)
  contrast2 <- c(1, Z_, prod(Z_))
  contrast_true <- c(rep(0, p_beta_true), c(1, Z_, prod(Z_)))
  truth2 <- contrast_true %*% true_coef
  
  output <- contrasts_geeglm(gee_coef = solution_mrt_prob, gee_cov = cov_mat_mrt_prob, contrast1 - contrast2)
  estimate <- output$Estimate
  mc_estimate <- mc_estimate_all[Z1 == 1, V1] - mc_estimate_all[Z1 == -1, V1]
  truth <- truth1 - truth2
  label <- paste0("t < t^*: E[Y(A, Z_1 = 1) - Y(A, Z_1 = -1)]")
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  
  
  # next look at stage 2
  mc_estimate_all <- data_augmented[(time >= T_SMART), sum(weight_SMART * Y) / sum(weight_SMART), by = .(Z1, Z2)]
  all_combn_index <- combn(1:4, 2)
  for (combn_index in 1:ncol(all_combn_index)) {
    Z_ <- Z_all_numeric[all_combn_index[,combn_index][1], ]
    contrast1 <- c(1, Z_, prod(Z_))
    # responders only have Z1 effects
    contrast_r1 <- c(rep(0, p_beta_true), (c(1, Z_[1], 0, 0)))
    truth_r1 <- contrast_r1 %*% true_coef
    contrast_r0 <- c(rep(0, p_beta_true), (c(1, Z_, prod(Z_))))
    truth_r0 <- contrast_r0 %*% true_coef
    truth1 <- p_responder[(Z_[1]+1)/2+1] * truth_r1 + (1 - p_responder[(Z_[1]+1)/2+1]) * truth_r0
    
    Z_prime <- Z_all_numeric[all_combn_index[,combn_index][2], ]
    contrast2 <- c(1, Z_prime, prod(Z_prime))
    # responders only have Z1 effects
    contrast_r1 <- c(rep(0, p_beta_true), (c(1, Z_prime[1], 0, 0)))
    truth_r1 <- contrast_r1 %*% true_coef
    contrast_r0 <- c(rep(0, p_beta_true), (c(1, Z_prime, prod(Z_prime))))
    truth_r0 <- contrast_r0 %*% true_coef
    truth2 <- p_responder[(Z_prime[1]+1)/2+1] * truth_r1 + (1 - p_responder[(Z_prime[1]+1)/2+1]) * truth_r0
    
    output <- contrasts_geeglm(gee_coef = solution_mrt_prob, gee_cov = cov_mat_mrt_prob, contrast1 - contrast2)
    estimate <- output$Estimate
    mc_estimate <- mc_estimate_all[Z1 == Z_[1] & Z2 == Z_[2], V1] - mc_estimate_all[Z1 == Z_prime[1] & Z2 == Z_prime[2], V1]
    truth <- truth1 - truth2
    label <- paste0("t >= t^*: E[Y(A, Z = (", Z_[1], ", ", Z_[2], ")) - Y(A, Z = (", Z_prime[1], ", ", Z_prime[2], "))]")
    dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - truth))
    
  }
  
  
  
  dat_marginal_effects <- data.table(dat_marginal_effects)
  colnames(dat_marginal_effects)[1] <- c("label")
  colnames(dat_marginal_effects)[(ncol(dat_marginal_effects)-6):ncol(dat_marginal_effects)] <- c("estimate", "se", "lci", "uci",
                                                                                                 "mc_estimate", "truth", "bias")
  rownames(dat_marginal_effects) <- NULL
  dat_marginal_effects <- data.table(dat_marginal_effects)
  for (col in colnames(dat_marginal_effects)[-c(1)]){
    set(dat_marginal_effects, j=col, value =as.numeric(dat_marginal_effects[[col]]))
  }
  dat_marginal_effects[, cover := (truth >= lci & truth <= uci)*1]
  dat_marginal_effects[, contrast_index := 1:nrow(dat_marginal_effects)]
  return(dat_marginal_effects)
}



p_tilde = 0.5
data <- generate_random_data_hybrid(N, T = T, T_SMART = T_SMART, pt_setting = pt_setting, responder_setting = responder_setting, 
                                    eta_A = eta_A, eta_state = eta_state,
                                    beta_interaction = beta_interaction, gamma_Z = gamma_Z, seed = seed, mc.cores = 1L)
out <- run_ee(data, p_tilde)
solution <- out$solution
cov_mat <- out$cov_mat
solution_mrt_prob <- out$solution_mrt_prob
cov_mat_mrt_prob <- out$cov_mat_mrt_prob
data_augmented <- out$data_augmented

dat_marginal_effects <- calculate_marginal_effects(solution, cov_mat, solution_mrt_prob, cov_mat_mrt_prob, data_augmented)
dat_marginal_effects
dat_marginal_effects


path_name <- paste0("simulations/data/hdwcls")
if(!dir.exists(path_name)) dir.create(path_name)
file_name <- paste0(path_name, "/N", N, "_", pt_setting, "_R", responder_setting, "_seed", seed, ".RData")
save(solution, dat_marginal_effects,  file = file_name)















### alternatively, use predicted values in the first step as outcomes in the second step --------------
run_ee <- function(data, p_tilde){
  
  ## next we need to augment data to accommodate for SMART design
  rows_to_replicate <- data[responder == 1, ]
  rows_to_replicate[, weight_SMART := 1 / p_SMART[1]]
  rows_to_replicate_observed <- rows_to_replicate
  rows_to_replicate_observed$replicant <- 1L
  rows_to_replicate_observed[, Z2 := Z1]
  rows_to_replicate_pseudo <- rows_to_replicate
  rows_to_replicate_pseudo$replicant <- 2L
  rows_to_replicate_pseudo[, Z2 := Z1 * (-1)]
  rows_not_to_replicate <- data[responder == 0, ]
  rows_not_to_replicate$replicant <- 1
  rows_not_to_replicate[, weight_SMART := 1 / (p_SMART[1] * p_SMART[2])]
  data_augmented <- rbind(rows_not_to_replicate,
                          rows_to_replicate_observed,
                          rows_to_replicate_pseudo)
  rm(rows_to_replicate, rows_to_replicate_observed, rows_to_replicate_pseudo, rows_not_to_replicate)
  data_augmented <- data_augmented[order(id, time, replicant)]
  
  # calculate weight of the MRT actions
  data_augmented[, weight_MRT := (p_tilde / trt_prob_t)^A * ((1 - p_tilde) / (1 - trt_prob_t))^(1-A)]
  data_augmented[, weight := weight_MRT * weight_SMART]
  data_augmented[, A_tilde_c := A - p_tilde]
  
  # compute the centered control variables by Z
  data_augmented[, state_mean := sum(weight_SMART*state) / sum(weight_SMART), by = .(time, Z1, Z2)]
  data_augmented[, state_c := state - state_mean]
  data_augmented[, stateZ1_mean := sum(weight_SMART*state*Z1) / sum(weight_SMART), by  = .(time, Z1, Z2)]
  data_augmented[, stateZ1_c := state*Z1 - stateZ1_mean]
  data_augmented[, Z2_after := (time >= T_SMART)*Z2]
  data_augmented[, before := (time < T_SMART)*1]
  data_augmented[, after := (time >= T_SMART)*1]
  
  
  ## run GEE
  model1 <- geeglm(Y ~ -1 + state_c + stateZ1_c +
                     A_tilde_c + A_tilde_c:Z1 +
                     A_tilde_c:Z2:after + A_tilde_c:Z1:Z2:after +
                     before + before:Z1 +
                     after + after:(Z1*Z2),
                   data = data_augmented, weights = weight, id = id)
  select_var1 <- c("A_tilde_c", "A_tilde_c:Z1",
                   "A_tilde_c:Z2:after", "A_tilde_c:Z1:Z2:after",
                   "before", "Z1:before",
                   "after", "Z1:after", "Z2:after", "Z1:Z2:after")
  select_idx <- match(select_var1, names(model1$coefficients))
  solution <- model1$coefficients[select_idx]
  cov_mat <- model1$geese$vbeta[select_idx,]
  cov_mat <- cov_mat[, select_idx]
  names(solution) <- gsub("A_tilde_c", "(At - pt_tilde)", names(solution))
  
  
  
  
  ### solve for the effect of Z under MRT randomization probability
  # get the predicted Y under the alternation MRT randomization probability
  design_mat_model1 <- model.matrix( ~ -1 + state_c + stateZ1_c +
                                       A_tilde_c + A_tilde_c:Z1 +
                                       A_tilde_c:Z2:after + A_tilde_c:Z1:Z2:after +
                                       before + before:Z1 +
                                       after + after:(Z1*Z2), data = data_augmented)
  design_mat_model1 <- design_mat_model1[, select_var1]
  data_augmented[, Y_predicted := design_mat_model1 %*% solution]
  ## run GEE
  # model2 <- geeglm(Y ~ Z1*Z2_after, 
  #               data = data_augmented, weights = weight_SMART, id = id)
  model2 <- geeglm(Y_predcited ~ Z1*Z2_after,
                   data = data_augmented, weights = weight_SMART, id = id)
  select_var2 <- c("(Intercept)", "Z1", "Z2_after", "Z1:Z2_after")
  select_idx <- match(select_var2, names(model2$coefficients))
  solution_mrt_prob <- model2$coefficients[select_idx]
  cov_mat_mrt_prob <- model2$geese$vbeta[select_idx,]
  cov_mat_mrt_prob <- cov_mat_mrt_prob[, select_idx]
  names(solution_mrt_prob) <- gsub("A_tilde_c:", "(At - pt_tilde):", names(solution_mrt_prob))
  
  
  
  # get the variances of the predicted Y
  # get the design matrix without A
  design_mat_m <- model.matrix(Y ~ Z1*Z2_after, data = data_augmented)
  design_mat_m <- design_mat_m[, select_var2]
  design_mat_m_weighted <- design_mat_m * data_augmented$weight_SMART
  design_mat_m_list <- split(data.table(design_mat_m_weighted), data_augmented$id)
  design_mat_m_list <- lapply(design_mat_m_list, as.matrix)
  N <- length(design_mat_m_list)
  
  design_mat_m_weighted_sqrt <- design_mat_m * sqrt(data_augmented$weight_SMART)
  mat_multiplier <- solve(crossprod(design_mat_m_weighted_sqrt))
  
  # 1. \sum_i \sum_t W_i m_it t(m_it) var(Y_it)
  # initialize \sum_i \sum_t m_it t(m_it) var(Y_it)
  mm_var_Y <- 0
  for (i in 1:nrow(data_augmented)) {
    # mm_var_Y <- mm_var_Y + tcrossprod(design_mat_m_weighted[i,]) * data_augmented$Y_predicted[i]
    mm_var_Y <- mm_var_Y + tcrossprod(design_mat_m_weighted[i,]) * (t(design_mat_model1[i,]) %*% cov_mat %*% design_mat_model1[i,])[1]
  }
  
  
  # mm_var_Y <- 0
  # for (i in 1:nrow(data_augmented)) {
  #   mm_var_Y <- mm_var_Y + design_mat_m_weighted[i,] * data_augmented$Y_predicted[i]
  # }
  # mat_multiplier %*% mm_var_Y
  
  
  
  # 2. covariance part
  calculate_cov <- function(i){
    mm_cov_t_tprime <- 0
    design_mat_m_i <- design_mat_m_list[[i]]
    for (t in 2:nrow(design_mat_m_i)) {
      for (tprime in 1:(t-1)) {
        mm_cov_t_tprime <- mm_cov_t_tprime + tcrossprod(design_mat_m_i[t,], design_mat_m_i[tprime,]) * (t(design_mat_model1[t,]) %*% cov_mat %*% design_mat_model1[tprime,])[1]
        mm_cov_t_tprime <- mm_cov_t_tprime + tcrossprod(design_mat_m_i[tprime,], design_mat_m_i[t,]) * (t(design_mat_model1[tprime,]) %*% cov_mat %*% design_mat_model1[t,])[1]
      }
    }
    return(mm_cov_t_tprime)
  }
  mm_cov_t_tprime <- mclapply(1:N, calculate_cov, mc.cores = 5L)
  mm_cov_t_tprime <- Reduce('+', mm_cov_t_tprime)
  
  cov_mat_mrt_prob <- mat_multiplier %*% (mm_var_Y + mm_cov_t_tprime) %*% mat_multiplier
  
  
  
  return(list(solution = solution, cov_mat = cov_mat,
              solution_mrt_prob = solution_mrt_prob, cov_mat_mrt_prob = cov_mat_mrt_prob,
              data_augmented = data_augmented))
}




out <- run_ee(data, p_tilde)
solution <- out$solution
cov_mat <- out$cov_mat
solution_mrt_prob <- out$solution_mrt_prob
cov_mat_mrt_prob <- out$cov_mat_mrt_prob
data_augmented <- out$data_augmented

dat_marginal_effects <- calculate_marginal_effects(solution, cov_mat, solution_mrt_prob, cov_mat_mrt_prob, data_augmented)
dat_marginal_effects
dat_marginal_effects


# dat_marginal_effects <- calculate_marginal_effects(solution)
path_name <- paste0("simulations/data/hdwcls_predicted")
if(!dir.exists(path_name)) dir.create(path_name)
file_name <- paste0(path_name, "/N", N, "_", pt_setting, "_R", responder_setting, "_seed", seed, ".RData")
save(solution, dat_marginal_effects,  file = file_name)


