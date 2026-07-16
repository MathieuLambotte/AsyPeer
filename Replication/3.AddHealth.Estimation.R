###########################################################################
################# Houndetoungan and Lambotte (202x) #######################
######## Asymmetries in Peer Effects for Continuous Outcomes ##############
########################### Add Health Estimation #########################
###########################################################################
# This script replicates our empirical results using Add Health data.
# Before using this script, use first the script '2.AddHealth.data.R' to 
# prepare the data set.
# The data set should be saved in the folder with the path `OutDataPath`

# Last update: 2026-07-15

rm(list = ls())

library(PartialNetwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(AsyPeer)
library(openxlsx)


OutDataPath <- "PATH/TO/WHERE/PREPARED/DATA/WILL/BE/SAVED" # Where prepared data will be saved
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED" # Where results will be saved

# List of outcome variables
depvar  <- c("smoke","fight", "optimism", "drink")

# Remove fake isolated
rmfiso  <- FALSE

# This function estimates the symmetric and asymmetric model for each outcome
festim  <- function(outcome) {
  set.seed(123)
  cat("Outcome: ", outcome, "\n", sep = "")
  
  ########################################################
  ##################### Data Loading #####################
  ########################################################
  load(paste0(OutDataPath, "/", outcome, ".Rda"))
  
  exovar  <- c("age", "female", "grade", "hispanic", "racewhite", 
               "raceblack", "raceasian", "melhigh", "memhigh", "memiss", 
               "mjprof", "mjother", "mjmiss", "match", "nmatch")
  
  G       <- norm.network(G)
  match   <- data$match
  nmatch  <- data$nmatch
  y       <- data$y
  X       <- as.matrix(data[,exovar])
  GX      <- peer.avg(G, X)
  drop    <- NULL
  if (rmfiso)
    drop  <- (match == 0) & (nmatch > 0) # Observations to be dropped: false isolated
  
  ########################################################
  ###################### Estimation ######################
  ########################################################
  ## Arguments for gen.inst
  gen.inst.arg <- list(estimator = "rf",
                       full = TRUE,
                       nfold = 5,
                       num.trees = 1500,
                       mtry = 2 * ncol(X),
                       replace = FALSE,
                       sample.fraction = 0.8,
                       num.threads = 8,
                       seed = 123)
  
  ## Estimation assuming symmetry
  ssym <- summary(asypeer.estim(formula = y ~ X + GX, 
                                Glist = G, 
                                fixed.effects = TRUE,
                                asymmetry = FALSE,
                                HAC = "cluster",
                                drop = drop,
                                gen.inst.arg = gen.inst.arg), 
                  diagnostic = TRUE, KP = TRUE) 
  
  print(ssym)
  
  ## Estimation assuming asymmetry
  # Results seems good but weak instruments for future expectation
  sasym <- summary(asypeer.estim(formula = y ~ X + GX, 
                                 Glist = G, 
                                 fixed.effects = TRUE,
                                 spillover = FALSE,
                                 asymmetry = TRUE,
                                 HAC = "cluster",
                                 drop = drop,
                                 gen.inst.arg = gen.inst.arg), 
                   diagnostic = TRUE, KP = TRUE)
  
  print(sasym)
  
  ## Saving results
  save(ssym, sasym, file = paste0(OutResPath, "/", outcome, 
                                  ifelse(rmfiso, ".noFakeIso", ""), ".Rda"))
  invisible(NULL)
}

# We apply the function to each outcome
sapply(depvar, function(k) festim(k))

# Full name of the outcome
OUTCOME <- c("Smoking","Fighting","Optimism","Drinking")

# This function gathers the estimates to create 
# Argument k is the index of the outcome in depvar
fres <- function(k) {
  # Load results
  load(paste0(OutResPath, "/", depvar[k], 
              ifelse(rmfiso, ".noFakeIso", ""), ".Rda"))
  
  # First column
  coef_names <- rownames(sasym$coefficients)
  coef_names <- c(coef_names[1:2], "beta", coef_names[-(1:2)])
  coef_names <- unlist(lapply(coef_names, \(s) paste0(c("est:", "se:"), s)))
  
  outest     <- data.frame(coef = coef_names)
  
  # Extract estimates and SEs
  coef_asym  <-  unlist(lapply(rownames(sasym$coefficients), \(s) paste0(c("est:", "se:"), s)))
  coef_sym   <-  unlist(lapply(rownames(ssym$coefficients), \(s) paste0(c("est:", "se:"), s)))
  outest     <- outest %>%
    left_join(data.frame(coef = coef_asym,
                         value1 = c(apply(sasym$coefficients, 1, \(s) {
                           c(sprintf("%.3f", s[1]), 
                             paste0("(", sprintf("%.3f", s[2]), ")"))
                         }))), by = "coef") %>%
    left_join(data.frame(coef = coef_sym,
                         value2 = c(apply(ssym$coefficients, 1, \(s) {
                           c(sprintf("%.3f", s[1]), 
                             paste0("(", sprintf("%.3f", s[2]), ")"))
                         }))), by = "coef")%>%
    mutate(across(c(value1, value2), ~ ifelse(is.na(.), "-", .)))
  
  # Test
  test  <- data.frame(coef  = c("Test Symmetry",
                                "",
                                "1st Stage F (ybar)", 
                                "1st Stage F (ycheck)",
                                "KP LM test p-value"),
                      value1 = c(sasym$gmm$diffbeta[1],
                                 sasym$gmm$diffbeta[2],
                                 sasym$diagnostics[1:2, "statistic"],
                                 sasym$diagnostics[3, "p-value"]),
                      value2 = c(NA,NA,ssym$diagnostics[1, "statistic"], NA, 
                                 ssym$diagnostics[2, "p-value"])) %>%
    mutate(across(c(value1, value2), ~ ifelse(is.na(.), "-", sprintf("%.3f", .))))
  test[2,2] <-  paste0("(", test[2,2], ")")
  rownames(test) <- NULL
  out    <- outest %>% bind_rows(test)
  
  # Total Peer Effect range
  minbeta <-  min(sasym$gmm$Estimate[c("betal", "betah")])
  maxbeta <- max(sasym$gmm$Estimate[c("betal", "betah")])
  boundl <- unname(minbeta/ (1 + minbeta))
  boundh <- unname(maxbeta / (1 + maxbeta))
  bound <- unname(ssym$gmm$Estimate["beta"] / (1 + ssym$gmm$Estimate["beta"]))
  PE_range <- data.frame(cbind("Total Peer Effects Range",paste0("[", deparse(round(boundl,3)),", ", deparse(round(boundh,3)),"]", sep = ""),
                               round(bound,3)))
  colnames(PE_range) <- c("coef","value1","value2")
  out    <- out %>% bind_rows(PE_range)
  
  # Colnames
  colnames(out) <- c("coef", paste0(OUTCOME[k], c("_asym", "_sym")))
  
  return(out)
}
out    <- lapply(1:length(depvar), fres)
out    <- out[[1]] %>% left_join(out[[2]], by = "coef") %>% 
  left_join(out[[3]], by = "coef") %>% left_join(out[[4]], by = "coef")

# Pretty coeff names
name_col <- rep("", nrow(out)-6) # 6 last rows are the tests and total peer effect
name_col[seq(1, length(name_col), by = 2)] <- c("$\\beta^l$", "$\\beta^h$","$\beta$","Age","Female", "Grade","Hispanic","White","Black","Asian","melhigh","memhigh","memiss",
                                                "mjprof","mjother","mjmiss","degree_matched","degree_nmatch","Age","Female", "Grade","Hispanic","White","Black","Asian","melhigh","memhigh","memiss",
                                                "mjprof","mjother","mjmiss","degree_matched","degree_nmatch")  # 1st, 3rd, ... rows

out$coef[1:length(name_col)] <- name_col

wb <- createWorkbook()
addWorksheet(wb, "Estimates")

writeData(wb, "Estimates", out, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, paste0(OutResPath, "/estimation", 
                        ifelse(rmfiso, ".noFakeIso", ""), ".xlsx"), overwrite = TRUE)

