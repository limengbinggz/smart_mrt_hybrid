
if(Sys.info()["sysname"] %in% c("Darwin")){
  curr_dir <- "/Users/mengbing/Dropbox (University of Michigan)/from_box/research/GSRA_walter/simulations/"
  setwd(curr_dir)
} else{ # biostat cluster
  curr_dir <- "/home/mengbing/research/GSRA_walter/simulations/"
  setwd(curr_dir)
}
library(data.table)
library(knitr)
library(kableExtra)
library(dplyr)
options(knitr.kable.NA = '')

# helper paste to remove NA
paste3 <- function(..., sep=", ") {
   L <- list(...)
   L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
   ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
               gsub(paste0(sep,sep),sep,
                    do.call(paste,c(L,list(sep=sep)))))
   is.na(ret) <- ret==""
   ret
   }

pt_settings <- c("ADependOnZ1Z2") 
pt_setting_labels <- c("$p_t$ depends on $Z_1$ and $Z_2$. ") 
responder_setting <- c("Hdependent")
responder_setting_label <- c("Responder probability depends on history. ") 
methods <- c("hdwcls", "wcls", "wr") 
method_labels <- c("Hybrid", "WCLS", "WR")
# settings <- c("ADependOnZ1") #,"Aconstant", "ADependOnZ1Z2"
# setting_labels <- c("$p_t = 0.5$ is constant. ") #, "$p_t$ depends on $Z_1$. ", "$p_t$ depends on $Z_1$ and $Z_2$. "
# setting_labels <- paste(setting_labels, "$P_R$ depends on $Z_1$. ")
# methods <- c("hdwcls", "wcls", "wr") 
# method_labels <- c("Hybrid", "WCLS", "WR")

Ns <- c(100, 400)
seeds <- 1000
# seeds <- 100
# N = 200
T = 50
T_SMART = 14
p_SMART = rep(0.5, length(T_SMART) + 1)
p_responder = 0.5
num_digits_rounding <- 3

for (N in Ns) {
  
  # seed = 1
  # setting = "noSMART"
  # m = 2
  p_beta <- 4
  p_gamma <- 4
  dat <- list()
  dat_marginal <- list()
  for (s in 1:length(pt_settings)) {
    pt_setting <- pt_settings[s]
    print(pt_setting)
    pt_setting_label <- pt_setting_labels[s]
    # cat('\n======================Setting:', pt_setting, '======================\n')
    
    eta_A = c(0, 0, 0.1, 0.2)
    eta_state = c(-1, 0.1, 0.2)
    beta_interaction = c(0.4, -0.3, 0.3, -0.3, 0.2, 0.2)
    gamma_Z = c(0.2, -0.1, -0.1, 0.2, 0.2)
    
    pt_setting_labels[s] <- ""
    dat_marginal_effects_setting <- c()
    for (m in 1:length(methods)) {
      method <- methods[m]
      method_label <- method_labels[m]
      print(method)
      path_name <- paste0("data/", method) #, "/", pt_setting
      
      p <- p_beta + p_gamma
      estmates <- c()# matrix(0, nrow = p, ncol = seeds)
      dat_marginal_effects_setting_method <- c()
      for (seed in 1:seeds) {
        # cat('\n Seed', seed, '\n')
        file_name <- paste0(path_name, "/N", N, "_", pt_setting, "_R", responder_setting, "_seed", seed, ".RData")
        # tryCatch({
        load(file_name)
        estmates <- rbind(estmates, solution)
        # estmates[, seed] <- solution
        dat_marginal_effects[, contrast_index := 1:nrow(dat_marginal_effects)]
        dat_marginal_effects[, cover := truth >= lci & estimate <= uci]
        dat_marginal_effects <- dat_marginal_effects[, .(seed, contrast_index, label, truth, mc_estimate, estimate,
                                                         se, cover)]
        dat_marginal_effects_setting_method <- rbind(dat_marginal_effects_setting_method, dat_marginal_effects)
        # }, error = function(cond) {message("No file")})
      }
      if (!is.null(dat_marginal_effects_setting_method)){
        dat_marginal_effects_setting_method[, pt_setting := pt_setting]
        dat_marginal_effects_setting_method[, method_label := method_label]
        dat_marginal_effects_setting_method[, method_index := m]
        # dat_marginal_effects_setting_method <- rbind(dat_marginal_effects_setting_method, dat_marginal_effects)
        
      }
      estimates_point <- round(rowMeans(estmates), 2)
      estimates_sd <- round(apply(estmates, 1, sd), 2)
      result <- mapply(function(x, y) paste0(x, " (", y, ")"), x = estimates_point, y = estimates_sd)
      result <- mapply(function(x, y) if_else(is.na(x) | is.na(y), "", paste0("\\makecell{", round(x, 2), "\\\\(", round(y, 2), ")}")), x = estimates_point, y = estimates_sd)
      dat_marginal_effects_setting <- rbind(dat_marginal_effects_setting, dat_marginal_effects_setting_method)
    }
    dat_marginal_effects_setting <- data.table(dat_marginal_effects_setting)
    # to_replace <- which(dat_marginal_effects_setting$label == "t <= t^*: E[Y(A = 1, Z = (1, -1)] - E[Y(A = 0, Z = (1, -1)]")
    # if (length(to_replace) > 0) {
    #   dat_marginal_effects_setting$label[to_replace] <- "t < t^*: E[Y(A = 1, Z = (1, -1)] - E[Y(A = 0, Z = (1, -1)]"
    # }
    dat_marginal[[pt_setting]] <- dat_marginal_effects_setting
    
    ## add ratio of asymptotic variance between hybrid and WR
    dat_wr <- dat_marginal[[pt_setting]][method_label == "WR", ]
    dat_se_hybrid <- dat_marginal[[pt_setting]][method_label == "Hybrid", .(seed, label, se)]
    colnames(dat_se_hybrid)[3] <- "se_hybrid"
    dat_wr <- merge(x = dat_wr, y = dat_se_hybrid, by = c("seed", "label"), all.x = TRUE)
    dat_wr[, se_hybrid_ratio := (se / se_hybrid)^2]
    dat_marginal[[pt_setting]] <- merge(x = dat_marginal[[pt_setting]], y = dat_wr[, .(seed, label, method_label, se_hybrid_ratio)], by = c("seed", "label", "method_label"), all.x = TRUE)
    dat_marginal[[pt_setting]] <- dat_marginal[[pt_setting]][order(seed, method_label, contrast_index), ]
    
    
    ## add MC truth
    dat_marginal_truth <- dat_marginal[[pt_setting]][method_label == "Hybrid", mean(truth), by = c('label')]
    colnames(dat_marginal_truth)[2] <- c("mean_truth")
    dat_marginal[[pt_setting]] <- merge(x = dat_marginal[[pt_setting]], y = dat_marginal_truth, by = c('label'), all.x = TRUE)
    dat_marginal[[pt_setting]][, uci := estimate + 1.96 * se]
    dat_marginal[[pt_setting]][, lci := estimate - 1.96 * se]
    dat_marginal[[pt_setting]][, cover_mc := (truth >= lci & truth <= uci) * 1]
  }
  
  
  
  
  ## summarize results ------------------------------------------
  dat_marginal_summary <- list()
  for (s in 1:length(pt_settings)) {
    pt_setting <- pt_settings[s]
    dat_marginal_summary[[pt_setting]] <- dat_marginal[[pt_setting]][, mean(estimate - mean_truth), by = c('contrast_index', 'method_index', 'method_label', 'label')] 
    colnames(dat_marginal_summary[[pt_setting]])[ncol(dat_marginal_summary[[pt_setting]])] <- c("mean_est")
    dat_marginal_truth <- dat_marginal[[pt_setting]][method_label == "Hybrid", mean(truth), by = c('label')]
    colnames(dat_marginal_truth)[ncol(dat_marginal_truth)] <- c("mean_truth")
    dat_marginal_se <- dat_marginal[[pt_setting]][, mean(se), by = c('contrast_index', 'method_index', 'method_label', 'label')]
    colnames(dat_marginal_se)[ncol(dat_marginal_se)] <- c("mean_se")
    dat_marginal_cover <- dat_marginal[[pt_setting]][, mean(cover_mc), by = c('contrast_index', 'method_index', 'method_label', 'label')]
    colnames(dat_marginal_cover)[ncol(dat_marginal_cover)] <- c("mean_cover")
    dat_marginal_summary[[pt_setting]] <- merge(x = dat_marginal_truth, y = dat_marginal_summary[[pt_setting]], by = c('label'))
    dat_marginal_summary[[pt_setting]] <- merge(x = dat_marginal_se, y = dat_marginal_summary[[pt_setting]], by = c('contrast_index', 'method_index', 'method_label', 'label'))
    dat_marginal_summary[[pt_setting]] <- merge(x = dat_marginal_cover, y = dat_marginal_summary[[pt_setting]], by = c('contrast_index', 'method_index', 'method_label', 'label'))
    dat_marginal_summary[[pt_setting]][, mean_truth := if_else(is.na(mean_truth), "", paste0(round(mean_truth, 2)))]
    dat_marginal_summary[[pt_setting]][, est := if_else(is.na(mean_est), "", paste0(round(mean_est, 2)))]
    dat_marginal_summary[[pt_setting]][, se_est := if_else(is.na(mean_se), "", paste0(round(mean_se, 2)))]
    dat_marginal_summary[[pt_setting]][, cover_est := if_else(is.na(mean_cover), "", paste0(round(mean_cover, 2)))]
    # summarize ratio of se's
    dat_se_ratio <- dat_marginal[[pt_setting]][method_label =="WR", .(seed, method_label, label, se_hybrid_ratio)]
    dat_se_ratio[, mean_se_ratio := round(mean(se_hybrid_ratio), 2), by = label]
    dat_se_ratio[, sd_se_ratio := round(sd(se_hybrid_ratio), 2), by = label]
    dat_se_ratio <- unique(dat_se_ratio[, .(method_label, label, mean_se_ratio, sd_se_ratio)])
    dat_marginal_summary[[pt_setting]] <- merge(x = dat_marginal_summary[[pt_setting]], y = dat_se_ratio, by = c('method_label', 'label'), all.x = TRUE)
    
    # select columns to the final table
    dat_marginal_summary[[pt_setting]] <- dat_marginal_summary[[pt_setting]][, .(method_label, method_index, contrast_index, label, mean_truth, est, se_est, cover_est, mean_se_ratio, sd_se_ratio)]
    
    dat_marginal_summary[[pt_setting]]$label <- paste0("$", dat_marginal_summary[[pt_setting]]$label, "$")
    dat_marginal_summary[[pt_setting]]$label <- gsub("E", "\\\\EE", dat_marginal_summary[[pt_setting]]$label)
    dat_marginal_summary[[pt_setting]]$label <- gsub("A", "A_t", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]]$label <- gsub("Z1", "Z_1", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]]$label <- gsub("Z2", "Z_2", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]]$label <- gsub("Z =", "\\\\bar Z =", dat_marginal_summary[[pt_setting]]$label)
    # dat_marginal_summary[[pt_setting]]$label <- gsub("Z", "\\\\bar{Z}", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]]$label <- gsub("<=", "\\\\leq", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]]$label <- gsub(">=", "\\\\geq", dat_marginal_summary[[pt_setting]]$label) 
    dat_marginal_summary[[pt_setting]] <- dat_marginal_summary[[pt_setting]][order(method_index, contrast_index), ]
    dat_marginal_summary[[pt_setting]] <- dat_marginal_summary[[pt_setting]][method_label == "WR", contrast_index := contrast_index + 22]
    dat_marginal_summary[[pt_setting]] <- reshape(dat_marginal_summary[[pt_setting]][,-2], idvar = c("contrast_index", "label", "mean_truth"), timevar = "method_label", direction = "wide")
    dat_marginal_summary[[pt_setting]] <- dat_marginal_summary[[pt_setting]][, -1]
    colnames(dat_marginal_summary[[pt_setting]]) <- gsub("est.", "", colnames(dat_marginal_summary[[pt_setting]]))
    dat_marginal_summary[[pt_setting]][, mean_se_ratio.WCLS := NULL]
    dat_marginal_summary[[pt_setting]][, sd_se_ratio.WCLS := NULL]
    dat_marginal_summary[[pt_setting]][, mean_se_ratio.Hybrid := NULL]
    dat_marginal_summary[[pt_setting]][, sd_se_ratio.Hybrid := NULL]
    
    # colnames(dat_marginal_summary[[pt_setting]]) <- c("Contrast", "Coefficients", "True", "Estimated", "SE", "Coverage", "MC Estimated")
    colnames(dat_marginal_summary[[pt_setting]])[1:2] <- c("Contrast", "True")
    # colnames(dat_marginal_summary[[pt_setting]])[3:5] <- c("Mean", "True")
    dat_marginal_summary[[pt_setting]] <- cbind(paste0("(", 1:nrow(dat_marginal_summary[[pt_setting]]), ")"), dat_marginal_summary[[pt_setting]])
  }
  
  
  
  
  ## Split into 3 tables for each scenario ------------------------------------------
  for (s in 1:length(pt_settings)) {
    pt_setting <- pt_settings[s]
    # remove the linear combination of coefficients in printing table
    tab_subset <- dat_marginal_summary[[pt_setting]]#[, -c("mc_Hybrid", "mc_WCLS", "mc_WR")] #, "MC Estimated"
    
    ###  the first table: contrasts between At = 1 vs At = 0 -----------
    table_1 <- tab_subset[1:8,]
    # remove WR method from here
    table_1 <- table_1[, -grep("WR", colnames(table_1), value = TRUE), with = FALSE]
    fileConn <- file(paste0("output/table_N", N, "_", pt_setting, "_R", responder_setting, "_contrastsA.txt"))
    table_1 %>%
      kbl(caption = "Comparison of effects on the proximal outcome between $A_t = 1$ and 0, for a fixed DTR or averaging over DTRs. ", 
          col.names = c("", "Contrast", "True", rep(c("Bias", "SE", "CP"), 2)),
          label = "sim:results:contrastsA",
          format = "latex", booktabs = TRUE, position = 'H', 
          align = c('c', 'X', rep('l', ncol(tab_subset)-2)),
          table.envir = "subtable",
          linesep = "", escape = FALSE) %>%
      add_header_above(c(" ", " ", " ",
                         "Hybrid" = 3, "WCLS" = 3)) %>%
      row_spec(c(6), hline_after = T) %>%
      kable_styling(latex_options = "repeat_header", latex_table_env = "tabularx") %>%
      gsub("\\{cX", "\\{\\\\linewidth\\}\\{cX", .) %>%
      gsub("\\caption", "\\scriptsize\\\\caption\n", .) %>%
      writeLines(con = fileConn)
    close(fileConn)
    cat("\n\n")
    
    
    
    ###  the second table: contrasts between DTRs averaging over At -----------
    table_2 <- tab_subset[23:29,]
    # remove WR method from here
    table_2 <- table_2[, -grep("WCLS", colnames(table_2), value = TRUE), with = FALSE]
    table_2[, V1 := paste0("(", 1:nrow(table_2), ")")]
    fileConn <- file(paste0("output/table_N", N, "_", pt_setting, "_R", responder_setting, "_contrastsZ.txt"))
    table_2 %>%
      kbl(caption = "Comparison of effects on the proximal outcome between DTRs averaging over MRT treatments. ", 
          col.names = c("", "Contrast", "True", rep(c("Bias", "SE", "CP"), 2), "mRE", "sdRE"),
          label = "sim:results:contrastsZ",
          format = "latex", booktabs = TRUE, position = 'H', 
          align = c('c', 'X', rep('l', ncol(tab_subset)-2)),
          table.envir = "subtable",
          linesep = "", escape = FALSE) %>%
      add_header_above(c(" ", " ", " ",
                         "Hybrid" = 3, "WR" = 3)) %>%
      row_spec(c(1), hline_after = T) %>%
      kable_styling(latex_options = "repeat_header", latex_table_env = "tabularx") %>%
      gsub("\\{cX", "\\{\\\\linewidth\\}\\{cX", .) %>%
      gsub("\\caption", "\\scriptsize\\\\caption\n", .) %>%
      writeLines(con = fileConn)
    close(fileConn)
    cat("\n\n")
    
    
    ###  the third table: contrasts between DTRs averaging over At -----------
    table_3 <- tab_subset[9:22,]
    # remove WR method from here
    table_3 <- table_3[, -grep("WCLS", colnames(table_3), value = TRUE), with = FALSE]
    table_3 <- table_3[, -grep("WR", colnames(table_3), value = TRUE), with = FALSE]
    table_3[, V1 := paste0("(", 1:nrow(table_3), ")")]
    fileConn <- file(paste0("output/table_N", N, "_", pt_setting, "_R", responder_setting, "_contrastsZFixedA.txt"))
    table_3 %>%
      kbl(caption = "Comparison of effects on the proximal outcome between DTRs for a fixed MRT treatment $A_t$. ", 
          col.names = c("", "Contrast", "True", rep(c("Bias", "SE", "CP"), 1)),
          label = "sim:results:contrastsZ:fixedA",
          format = "latex", booktabs = TRUE, position = 'H', 
          align = c('c', 'X', rep('l', ncol(tab_subset)-2)),
          table.envir = "subtable",
          linesep = "", escape = FALSE) %>%
      add_header_above(c(" ", " ", " ",
                         "Hybrid" = 3)) %>%
      row_spec(c(7), hline_after = T) %>%
      kable_styling(latex_options = "repeat_header", latex_table_env = "tabularx") %>%
      gsub("\\{cX", "\\{\\\\linewidth\\}\\{cX", .) %>%
      gsub("\\caption", "\\scriptsize\\\\caption\n", .) %>%
      writeLines(con = fileConn)
    close(fileConn)
    cat("\n\n")
  }
  
}



