
paste_coef <- function(coef_values, coef_names){
  coef_labels <- c()
  for (j in 1:length(coef_names)) {
    label_j <- NULL
    if (coef_values[j] > 1){
      if (length(coef_labels) != 0){
        label_j <- paste0("+", as.integer(coef_values[j]), coef_names[j])
      } else{ # if the first element
        label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
      }
    } else if (coef_values[j] == 1){
      if (length(coef_labels) != 0){
        label_j <- paste0("+", coef_names[j])
      } else{ # if the first element
        label_j <- paste0(coef_names[j])
      }
    } else if (coef_values[j] == -1){
      label_j <- paste0("-", coef_names[j])
    } else if (coef_values[j] < -1){
      label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
    }
    coef_labels <- c(coef_labels, label_j)
  }
  coef_labels <- paste0(coef_labels, collapse = "")
  coef_labels <- paste0("$", coef_labels, "$")
  return(paste0(coef_labels, collapse = ""))
}


# At_prob = p_tilde
# coef_beta_values = contrast_beta; coef_beta_names = coef_names[1:p_beta]
# coef_gamma_values = contrast_gamma; coef_gamma_names = coef_names[1:p_gamma+p_beta]
paste_coef_fixedAt <- function(At, At_prob, coef_beta_values, coef_beta_names, coef_gamma_values, coef_gamma_names){
  # coef_labels <- paste_coef(coef_beta_values, coef_beta_names[1:p_beta])
  # coef_labels <- paste0("(A_t - ", p_tilde, ") (", coef_labels, ")")
  # coef_labels <- paste0(coef_labels, "+", paste_coef(contrast_gamma, coef_names[1:p_gamma+p_beta]))
  
  coef_labels <- paste0("(", At, "- ", p_tilde, ") (")
  coef_names <- coef_beta_names
  coef_values <- coef_beta_values
  for (j in 1:length(coef_names)) {
    label_j <- NULL
    if (coef_values[j] > 1){
      if (length(coef_labels) == 1){
        label_j <- paste0("+", as.integer(coef_values[j]), coef_names[j])
      } else{ # if the first element
        label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
      }
    } else if (coef_values[j] == 1){
      if (length(coef_labels) > 1){
        label_j <- paste0("+", coef_names[j])
      } else{ # if the first element
        label_j <- paste0(coef_names[j])
      }
    } else if (coef_values[j] == -1){
      label_j <- paste0("-", coef_names[j])
    } else if (coef_values[j] < -1){
      label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
    }
    coef_labels <- c(coef_labels, label_j)
  }
  coef_labels <- paste0(coef_labels, collapse = "")
  coef_labels <- paste0(coef_labels, ") ")
  
  coef_names <- coef_gamma_names
  coef_values <- coef_gamma_values
  for (j in 1:length(coef_names)) {
    label_j <- NULL
    if (coef_values[j] > 1){
      if (length(coef_labels) == 1){
        label_j <- paste0("+", as.integer(coef_values[j]), coef_names[j])
      } else{ # if the first element
        label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
      }
    } else if (coef_values[j] == 1){
      # if (length(coef_labels) > 1){
      label_j <- paste0("+", coef_names[j])
      # } else{ # if the first element
      #   label_j <- paste0(coef_names[j])
      # }
    } else if (coef_values[j] == -1){
      label_j <- paste0("-", coef_names[j])
    } else if (coef_values[j] < -1){
      label_j <- paste0(as.integer(coef_values[j]), coef_names[j])
    }
    coef_labels <- c(coef_labels, label_j)
  }
  coef_labels <- paste0(coef_labels, collapse = "")
  
  coef_labels <- paste0("$", coef_labels, "$")
  return(coef_labels)
}

