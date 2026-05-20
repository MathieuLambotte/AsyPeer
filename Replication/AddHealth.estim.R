##############################################################################################################
##############################################################################################################
#################### Asymmetric Peer Effect Models by A. Houndetoungan and M. Lambotte #######################
##############################################################################################################
# This script replicates our empirical study using Add Health data
# Please, prior use the script AddHealth.data.R to prepare the data set for
# each outcome in am expect form for this scrip.
# This scripts assumes that the prepared data are saved in the folder with path `OutDataPath`

# Last update: 2026-02-14

rm(list = ls())

library(PartialNetwork)
library(dplyr)
library(AsyPeer)
library(openxlsx)

# This defines the path to the folder where prepared data using the script 
# AddHealth.data.R are saved.
# Note that the symbol / at the end of the path is necessary.
OutDataPath <- "PATH/TO/WHERE/PREPARED/DATA/IS/SAVED/"

# This defines the path to the folder where results will be saved.
# Note that the symbol / at the end of the path is necessary.
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED/" 

OutDataPath <- "~/Dropbox/Academy/1.Papers/AsymmetricPeer/AsyPeerCode/SOULD_BE_DELETED_LATER/"
OutResPath  <- "~/Dropbox/Academy/1.Papers/AsymmetricPeer/AsyPeerCode/Application/"

# OutDataPath <- "D:/Home/mlambotte/Documents/Dropbox/AsymmetricPeer/AsyPeerCode/SOULD_BE_DELETED_LATER/"
# OutResPath <- "D:/Home/mlambotte/Documents/Dropbox/AsymmetricPeer/AsyPeerCode/Application/"

# List of outcome variable names
depvar  <- c("smoke", "drink", "fight", "futureperception")

# Results
# smoke 1.92 and 2.65
# drink 1.03 et 0.89
# fight 0.187 and 0.672
# futureperception 2.679 and -0.461

## This function prepares the data and estimate the models for each outcome
festim <-function(outcome) {
  cat("Outcome: ", outcome, "\n", sep = "")
  
  ########################################################
  ################### Data Preparation ###################
  ########################################################
  
  load(file = paste0(OutDataPath, outcome, ".Rda"))
  exovar  <- c("age", "female", "grade", "hispanic", "racewhite", 
               "raceblack", "raceasian", "melhigh", "memhigh", "memiss", 
               "mjprof", "mjother", "mjmiss")
  G       <- norm.network(G)
  match   <- data$match
  nmatch  <- data$nmatch
  y       <- data$y
  X       <- as.matrix(data[,exovar])
  GX      <- peer.avg(G, X)
  drop    <- (match == 0) & (nmatch > 0)
  
  ########################################################
  ###################### Estimation ######################
  ########################################################
  ## Estimation assuming symmetry
  ssym <- summary(asypeer.estim(formula = y ~ X + GX, 
                                Glist = G, 
                                power = 2,
                                fixed.effects = TRUE,
                                asymmetry = FALSE,
                                HAC = "hetero",
                                drop = drop), 
                  diagnostic = TRUE, KPtest = TRUE)
  
  ## Estimation assuming asymmetry
  sasym <- summary(asypeer.estim(formula = y ~ X + GX, 
                                 Glist = G, 
                                 estimator = "glm", 
                                 power = 2,
                                 nfold = 5, 
                                 fixed.effects = TRUE,
                                 spillover = FALSE,
                                 asymmetry = TRUE,
                                 HAC = "hetero",
                                 drop = drop), 
                   diagnostic = TRUE, KPtest = TRUE)
  
  ## Saving results
  save(ssym, sasym, file = paste0(OutResPath, outcome, ".Rda"))
  invisible(NULL)
}

# We apply the function to each outcome
sapply(depvar, festim) 

# Name of the outcome
OUTCOME <- c("Smoking", "Drinking","Fighting","Future Perception")# This function puts results together

# This function gathers the estimates to create 
# Argument k is the index of the outcome in depvar
fres <- function(k) {
  # Load results
  load(paste0(OutResPath, depvar[k], ".Rda"))
  
  # Coefficient names
  coef_names <- rownames(sasym$coefficients)
  
  # Extract estimates and SEs
  est <- round(sasym$coefficients[coef_names, "Estimate"], 3)
  se  <- round(sasym$coefficients[coef_names, "Std. Error"], 3)
  
  # Format SEs with parentheses
  se_fmt <- paste0("(", se, ")")
  
  # Interleave estimates and SEs
  values <- as.vector(rbind(est, se_fmt))
  
  # Create a 'Name' column: only for estimate rows, blank for SE rows
  name_col <- rep("", length(values))
  name_col[seq(1, length(values), by = 2)] <- c("$\\beta^l", "$\\beta^h$","Age","Female", "Grade","Hispanic","White","Black","Asian","melhigh","memhigh","memiss",
                                                "mjprof","mjother","mjmiss","Age","Female", "Grade","Hispanic","White","Black","Asian","melhigh","memhigh","memiss",
                                                "mjprof","mjother","mjmiss")  # 1st, 3rd, ... rows
  name_test   <-c("SymTest","")
  values_test <- as.vector(rbind(round(sasym$gmm$diffbeta[1],3), 
                                 paste0("(", round(sasym$gmm$diffbeta[2],3), ")")))
  test_table1 <-data.frame(
    Coefficient = name_test,
    Asymmetric = values_test,
    stringsAsFactors = FALSE
  )
  
  # Put in a data.frame
  pub_table1   <- data.frame(
    Coefficient = name_col,
    Asymmetric = values,
    stringsAsFactors = FALSE
  )
  diag <- sasym$diagnostics[,3]
  
  diag_df <- data.frame(
    Coefficient   = names(diag),
    Asymmetric = round(as.numeric(diag), 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  pub_table1 <-rbind(pub_table1[1:4,],c("$\\beta$",NA),c("",NA),
                     pub_table1[5:nrow(pub_table1),],diag_df,test_table1)
  
  # Coefficient names
  coef_names <- rownames(ssym$coefficients)
  
  # Extract estimates and SEs
  est <- round(ssym$coefficients[coef_names, "Estimate"], 3)
  se  <- round(ssym$coefficients[coef_names, "Std. Error"], 3)
  
  # Format SEs with parentheses
  se_fmt <- paste0("(", se, ")")
  
  # Interleave estimates and SEs
  values <- as.vector(rbind(est, se_fmt))
  
  # Put in a data.frame
  pub_table2 <- data.frame(
    Symmetric = values,
    stringsAsFactors = FALSE
  )
  diag <- ssym$diagnostics[,3]
  
  diag_df2 <- data.frame(
    Symmetric = round(as.numeric(diag), 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  pub_table2<-rbind(NA,NA,NA,NA,pub_table2,diag_df2[1,],NA,NA,diag_df2[2,],NA,NA)
  cbind(pub_table1,pub_table2)
}
out    <- lapply(1:length(depvar), fres)

names(out) <- OUTCOME

# Extract the Coefficient column (same for all)
coef_col <- out[[1]]$Coefficient

# Combine all value columns into a single wide table
wide_table <- data.frame(
  Coefficient = coef_col,
  Smoking_Asymmetric  = out$Smoking$Asymmetric,
  Smoking_Symmetric   = out$Smoking$Symmetric,
  Drinking_Asymmetric = out$Drinking$Asymmetric,
  Drinking_Symmetric  = out$Drinking$Symmetric,
  Fighting_Asymmetric    = out$Fighting$Asymmetric,
  Fighting_Symmetric     = out$Fighting$Symmetric,
  Perception_Asymmetric    = out$`Future Perception`$Asymmetric,
  Perception_Symmetric     = out$`Future Perception`$Symmetric,
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "Estimates")

writeData(wb, "Estimates", wide_table, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, paste0(OutResPath, "/estimation.xlsx"), overwrite = TRUE)

