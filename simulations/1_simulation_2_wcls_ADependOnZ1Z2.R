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
library(rootSolve) #multiroot
source("functions/simulate_SMART_MRT_data.R")
source("functions/utils.R")
library(geepack)



### Simulate data -------------------------------------------
seed <- as.integer(args[1])
pt_setting <- as.character(args[2])
responder_setting = as.character(args[3]) 
N <- as.integer(args[4])
# seed = 83
# pt_setting = "ADependOnZ1Z2" "ADependOnX1X2" "ADependOnX1"  "ADependOnX1" "ADependOnZ1" "Aconstant"  
# responder_setting = "Hdependent" "Z1dependent" "constant" "none" 
# N <- 100
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


calculate_marginal_effects <- function(solution, cov_mat){
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
  
  # randomization probabilities in 2 stages    # randomization probabilities in 2 stages  
  if (pt_setting == "Aconstant"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.5, 0.5) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.5, 0.5) %*% c(Z2==1, Z2==-1) * (1-R) + A_stage1_prob_calculate(Z1) * R
    }
  } else if (pt_setting == "ADependOnZ1"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1) * (1-R) + A_stage1_prob_calculate(Z1) * R
    }
  } else if (pt_setting == "ADependOnZ1Z2"){
    A_stage1_prob_calculate <- function(Z1){
      c(0.6, 0.4) %*% c(Z1==1, Z1==-1)
    }
    A_stage2_prob_calculate <- function(Z1, Z2, R){
      c(0.6, 0.4, -0.2, 0.2) %*% c(Z1==1, Z1==-1, (1-R)*Z2==1, (1-R)*Z2==-1)
    }
  } else if (pt_setting == "ADependOnX1X2"){
    A_stage1_prob_calculate <- function(Z1_){
      data[time < T_SMART & Z1 == Z1_, mean(trt_prob_t), ]
    }
    A_stage2_prob_calculate <- function(Z1_, Z2_, R){
      data[time >= T_SMART & Z1 == Z1_ & CZ2 == Z2_ & responder == R, mean(trt_prob_t), ]
    }
  }
  
  
  ## compute marginal responder probabilities given Z1
  if (responder_setting == "none"){
    p_responder <- c(0, 0)
  } else if (responder_setting == "constant"){
    p_responder <- c(0.5, 0.5)
  } else if (responder_setting == "Z1dependent"){
    p_responder <- c(0.45, 0.6) 
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
  # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  contrast_true <- c(1, c(Z_, prod(Z_)), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  mc_estimate <- data[(Z1 == Z_[1] & time < T_SMART), mean((Y - Y_counterfactual_A) * (2*A-1))]
  label <- paste0('t < t^*: E[Y(A = 1, Z_1 = ', Z_[1], ")] - E[Y(A = 0, Z_1 = ", Z_[1], ")]")
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
  mc_estimate <- data[(Z1 == Z_[1] & time < T_SMART), mean((Y - Y_counterfactual_A) * (2*A-1))]
  label <- paste0('t < t^*: E[Y(A = 1, Z_1 = ', Z_[1], ")] - E[Y(A = 0, Z_1 = ", Z_[1], ")]")
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  ## first look at stage 2
  before <- 0; after <- 1
  ## Effect of At for a given Z2 = 1 at t >= t^*
  # Z_index <- 2
  for (Z_index in 1:nrow(Z_all_numeric)) {
    Z_ <- Z_all_numeric[Z_index, ]
    contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
    # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
    contrast <- c(contrast_beta, rep(0, p_gamma))
    output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
    estimate <- output$Estimate
    
    # responders only have Z1 effects
    truth_r <- c(1, Z_[1], 0, 0) %*% true_coef[1:p_beta_true]
    # nonresponders only have Z1 and Z2 effects
    truth_nr <- c(1, c(Z_, prod(Z_))) %*% true_coef[1:p_beta_true] 
    truth <- p_responder[(Z_[1]+1)/2+1] * truth_r + (1 - p_responder[(Z_[1]+1)/2+1]) * truth_nr
    mc_estimate <- data_augmented[(Z1 == Z_[1] & Z2 == Z_[2] & time >= T_SMART), 
                                  sum(weight_SMART * (Y - Y_counterfactual_A) * (2*A-1)) / sum(weight_SMART)]
    label <- paste0('t >= t^*: E[Y(A = 1, Z = (', Z_[1], ", ",  Z_[2], ")] - E[Y(A = 0, Z = (", Z_[1], ", ",  Z_[2], ")]")
    dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  }
  
  
  ## Effect of At averaging over Z
  # stage 1
  before <- 1; after <- 0
  Z_ <- c(0, 0)
  contrast_beta <- c(1, Z_[1], after*Z_[2], after*prod(Z_))
  # contrast_beta <- c(before, before*Z_[1], after, Z_[1]*after, after*Z_[2], after*prod(Z_))
  contrast <- c(contrast_beta, rep(0, p_gamma))
  # estimate <- contrast %*% solution
  output <- contrasts_geeglm(gee_coef = solution, gee_cov = cov_mat, contrast)
  estimate <- output$Estimate
  mc_estimate <- mean(data[(time < T_SMART), (Y - Y_counterfactual_A)*(2*A - 1)])
  contrast_true <- c(1, rep(0, p_beta_true-1), rep(0, p_gamma_true))
  truth <- contrast_true %*% true_coef
  label <- paste0("t < t^*: E[Y(A = 1, Z1)] - E[Y(A = 0, Z1)]")
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
  dat_marginal_effects <- rbind(dat_marginal_effects, c(label, unlist(output), mc_estimate, truth, estimate - mc_estimate))
  
  
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


data <- generate_random_data_hybrid(N, T = T, T_SMART = T_SMART, pt_setting = pt_setting, responder_setting = responder_setting, 
                                    eta_A = eta_A, eta_state = eta_state,
                                    beta_interaction = beta_interaction, gamma_Z = gamma_Z, seed = seed, mc.cores = 1L) 




## create augmented data to calculate the truth
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

# compute the centered control variables by Z
data_augmented[, Z2_after := (time >= T_SMART)*Z2]
data_augmented[, before := (time < T_SMART)*1]
data_augmented[, after := (time >= T_SMART)*1]



### fit WCLS -------------------------------------------
p_tilde <- 0.5
# calculate weight of the MRT actions
data[, weight_MRT := (p_tilde / trt_prob_t)^A * ((1 - p_tilde) / (1 - trt_prob_t))^(1-A)]
data[, A_tilde_c := A - p_tilde]
data[, Z2_after := (time >= T_SMART)*Z2]
data[, before := (time < T_SMART)*1]
data[, after := (time >= T_SMART)*1]

## run GEE
out <- geeglm(Y ~ state + state:Z1 +
                A_tilde_c + A_tilde_c:(Z1*Z2_after),
              data = data, weights = weight_MRT, id = id)
select_var <- c("A_tilde_c", "Z1:A_tilde_c",
                "A_tilde_c:Z2_after", "Z1:A_tilde_c:Z2_after")

# out <- geeglm(Y ~ -1 + state + state:Z1 +
#                 A_tilde_c + A_tilde_c:Z1 +
#                 A_tilde_c:Z2:after + A_tilde_c:Z1:Z2:after +
#                 before + before:Z1 +
#                 after + after:(Z1*Z2),
#               data = data, weights = weight_MRT, id = id)
# select_var <- c("A_tilde_c", "Z1:A_tilde_c",
#                 "A_tilde_c:Z2:after", "Z1:A_tilde_c:Z2:after",
#                 "before", "Z1:before",
#                 "after", "Z1:after", "Z2:after", "Z1:Z2:after")
select_idx <- match(select_var, names(out$coefficients))
solution <- out$coefficients[select_idx]
cov_mat <- out$geese$vbeta[select_idx,]
cov_mat <- cov_mat[, select_idx]
names(solution) <- gsub("A_tilde_c", "(At - pt_tilde)", names(solution))


dat_marginal_effects <- calculate_marginal_effects(solution, cov_mat)
dat_marginal_effects
dat_marginal_effects

path_name <- "simulations/data/wcls"
if(!dir.exists(path_name)) dir.create(path_name)
file_name <- paste0(path_name, "/N", N, "_", pt_setting, "_R", responder_setting, "_seed", seed, ".RData")
save(solution, dat_marginal_effects, file = file_name)






