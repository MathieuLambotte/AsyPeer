###########################################################################
################# Houndetoungan and Lambotte (202x) #######################
######## Asymmetries in Peer Effects for Continuous Outcomes ##############
###################### Add Health Data Extraction #########################
###########################################################################

# Last updated: 2026-07-15

## This script imports raw Add Health data and prepares them for estimation.
## The outputs of this script include an .Rda file for each outcome.
## Each .Rda file includes a data set (y and X variables) and the network matrix G.
## Add Health data are imported from the path InDataPath.
## The .Rda files will be saved at the path OutDataPath.

rm(list = ls())
library(dplyr)
library(haven)
library(openxlsx)

InDataPath  <- "PATH/TO/RAW/DATA" # Where raw data are located
OutDataPath <- "PATH/TO/WHERE/PREPARED/DATA/WILL/BE/SAVED" # Where prepared data will be saved
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED" # Where results will be saved

# Importing data sets
# Friendship data set (WAVE I)
sfriend  <- read_xpt(paste0(InDataPath, "/sfriend.xpt")) %>% arrange(SQID) %>% 
  filter(!(SQID %in% c("999999", ""))) %>%  # Remove student with missing questionnaire ID (cannot be matched)
  mutate(across(ends_with("AID"), as.character)) 

# Inschool data set (WAVE I)
Inschool <- read_xpt(paste0(InDataPath, "/Inschool.xpt")) %>% arrange(SQID) %>% 
  filter(!(SQID %in% c("999999", "")), AID != "") %>% # Remove student with missing questionnaire ID and missing student ID (cannot be matched)
  filter(!(SSCHLCDE %in% c("080", "175"))) %>% # Schools with many missing values. Only 2 students remain after the data cleaning process
  mutate(across(c("SQID", "AID", "SSCHLCDE"), as.character))

# Variable to keep in the Inschool data set
to_keep  <- c("SQID", "AID", "SSCHLCDE", "S1", "S2", "S3", "S4", "S6A", "S6B", "S6C", "S6E", "S10A", "S10B", "S10C", "S10D", 
              "S12", "S14", paste0("S44A", 1:33), grep("S45", names(Inschool), value = TRUE), "S46A", "S46B", "S46C", 
              "S46D", "S48", grep("S59", names(Inschool), value = TRUE), "S62H", "S62K", "S62M", "S62N", "S62O", "S62P", 
              "S63", "S64")
Inschool <- Inschool %>% select(all_of(to_keep)) 

# Add the friendship data set to the Inschool data set
Inschool <- Inschool %>% left_join(sfriend, by = "SQID")

# Create variables that will be used
Inschool <- Inschool %>% 
  # We recoded some variables before constructing the ones that will be used in the models.
  # Replace 99 with NA
  mutate(across(c(all_of("S1"), starts_with("S45"), starts_with("S59"), starts_with("S62")), ~ ifelse(. == 99, NA, .))) %>%
  # Replace 9 with NA
  mutate(across(c(all_of(c("S2", "S48", "S63", "S64")), starts_with("S46")), ~ ifelse(. == 9, NA, .))) %>%
  # Replace values out of 1:4 with NA (this is for the GPA)
  mutate(across(all_of(c("S10A", "S10B", "S10C", "S10D")), ~ ifelse(. %in% (1:4), ., NA))) %>%
  # Variable about perception of future vary from zero to eight. 
  # But S45C (be killed) and S45D (get HIV) capture negative perception. 
  # Thus, we do 8 - variable so that they vary from zero to eight, where higher values indicate a more positive perception of the future.
  mutate(across(all_of(c("S45C", "S45D")), ~ 8 - .)) %>% 
  # Get trouble is converted to the number of times per week.
  mutate(across(all_of(c("S46A", "S46B", "S46C", "S46D")), ~ ifelse(. == 1, 0, ifelse(. == 2, 1, ifelse(. == 3, 2.5, ifelse(. == 4, 5, .)))))) %>%
  # Dangerous activities (we convert in week)
  mutate(across(starts_with("S59"), ~ ifelse(. == 1, 0.029, #1.5 per year / 52 week
                                             ifelse(. == 2, 0.25, # 0.5 per month / 4 weeks
                                                    ifelse(. == 3, 0.625, # 2.5 per month / 4 weeks
                                                           ifelse(. == 4, 1.5, # 1.5 per week
                                                                  ifelse(. == 5, 4, # 4 per weeks
                                                                         ifelse(. == 6, 7, .)))))))) %>% # everyday
  # Self esteem
  mutate(across(starts_with("S62"), ~ ifelse(. == 5, 0, ifelse(. == 4, 0.25, ifelse(. == 3, 0.5, ifelse(. == 2, 0.75, .)))))) %>%
  # Model variable creation
  mutate(age = S1, age2 = age^2, S1 = NULL, 
         male = as.integer(S2 == 1), female = as.integer(S2 == 2), S2 = NULL,
         grade = ifelse(S3 %in% 6:12, S3, NA), grade2 = grade^2, S3 = NULL,
         hispanic = as.integer(S4 == 1), S4 = NULL,
         racewhite = as.integer(S6A == 1), S6A = NULL,
         raceblack = as.integer(S6B == 1), S6B = NULL,
         raceasian = as.integer(S6C == 1), S6C = NULL, S6E = NULL,
         gpa = 5 - rowMeans(select(., c("S10A", "S10B", "S10C", "S10D")), na.rm = TRUE),
         across(all_of(c("S10A", "S10B", "S10C", "S10D")), ~ NULL),
         melhigh = as.integer(S12 %in% c(1, 2, 9, 10)), # less than HS
         mehigh = as.integer(S12 %in% c(3, 4)), # HS
         memhigh = as.integer(S12 %in% (5:8)), # more than HS
         memiss = as.integer(S12 >= 11 | is.na(S12)),
         S12 = NULL,
         mjhome = as.integer(S14 %in% c(1, 16:18)), # at hme
         mjprof = as.integer(S14 %in% c(2, 3)), # professional
         mjother = as.integer(S14 %in% c(4:15, 19)), # other
         mjmiss = as.integer(S14 >= 20 | is.na(S14)), # missing
         S14 = NULL,
         nclubs = ifelse(is.na(rowMeans(select(., all_of(paste0("S44A", 1:33))), na.rm = TRUE)), NA, 
                         rowSums(select(., all_of(paste0("S44A", 1:33))), na.rm = TRUE)), # This is because rowMeans returns NA when all variables are NA but not rowSums
         across(all_of(paste0("S44A", 1:33)), ~ NULL),
         optimism = rowMeans(select(., all_of(c("S45A", "S45E"))), na.rm = TRUE), 
         # 45.a WILL LIVE TO AGE 35
         # 45.e WILL GRADUATE FROM COLLEGE
         # High values mean that it will happen.
         across(starts_with("S45"), ~ NULL),
         trouble = rowMeans(select(., all_of(c("S46A", "S46B", "S46C", "S46D"))), na.rm = TRUE),
         across(all_of(c("S46A", "S46B", "S46C", "S46D")), ~ NULL),
         academiceffort = ifelse(S48 == 2, 0.66, ifelse(S48 == 3, 0.33, ifelse(S48 == 4, 0, S48))),
         S48 = NULL,
         smoke = S59A,
         drink = S59B,
         risky = rowMeans(select(., starts_with("S59")), na.rm = TRUE),
         across(starts_with("S59"), ~ NULL),
         selfesteem = rowMeans(select(., starts_with("S62")), na.rm = TRUE),
         across(starts_with("S62"), ~ NULL),
         physicalexercise = ifelse(S63 == 1, 1.5, ifelse(S63 == 2, 4.5, ifelse(S63 == 3, 6.5, ifelse(S63 == 4, 7.5, S63)))),
         S63 = NULL,
         fight = ifelse(S64 == 1, 1.5, # 1.5 per year
                        ifelse(S64 == 2, 4, # 4 per year 
                               ifelse(S64 == 3, 6.5, # 6.5 per year
                                      ifelse(S64 == 4, 8, S64)))), # 8 per year
         S64 = NULL) %>%
  arrange(SSCHLCDE, AID) %>% mutate(SCID = SSCHLCDE)

# Creating network data
mislist <- c(55555555, 77777777, 88888888, 99999999, 99959995)
# Check if mislist is ID
if (sum(Inschool$AID %in% mislist) > 0) {
  stop("mislist is an ID")
} else {
  cat("mislist is not an ID: OK", "\n")
}

# Exogenous characteristics (excluding reference variables for identification)
exovar       <- c("age", "age2", "female", "grade", "grade2", "hispanic", "racewhite", "raceblack", 
                  "raceasian", "melhigh", "memhigh", "memiss", "mjprof", "mjother",  "mjmiss")

# List of outcomes 
depvar       <- c("smoke","fight", "optimism", "drink")
# This loop creates data for each outcome
# The data consist of a data set of y and X as well as the network.
for (outcome in depvar) {
  # Exogenous characteristics and outcome
  va.names   <- c(exovar, outcome)
  
  # Remove observations with missing values
  data       <- Inschool %>% filter(if_all(all_of(va.names), ~ !is.na(.)))
  Sch        <- unique(data$SCID) # School ID
  M          <- length(Sch)       # Number of schools
  
  # remove friends from different groups
  # remove self friendship
  # remove friends non found
  N       <- nrow(data)
  dscf    <- rep(0, N)
  sfre    <- rep(0, N)
  nffr    <- rep(0, N)
  f_coln  <- paste0(rep(c("MF", "FF"), each = 5), rep(1:5, 2), "AID")
  for (i in 1:N) {
    for (j in f_coln) {
      k   <- which(data$AID == data[[j]][i])
      # remove if different school
      if (length(k) != 0) { # if friend found
        if(data$SCID[i] != data$SCID[k]) {
          data[[j]][i]   <- -1
          dscf[i]        <- dscf[i] + 1
        }
        # remove if self friendship
        if(data$AID[i] == data$AID[k]) {
          data[[j]][i]   <- -2
          sfre[i]        <- sfre[i] + 1
        }
      } else { # if friend not found
        if (!((data[[j]][i] %in% mislist) | is.na(data[[j]][i]))) { # if not missing or error code
          data[[j]][i]   <- -3
          nffr[i]        <- nffr[i] + 1
        }
      }
    }
  }
  cat("remove", sum(dscf), "link(s) because students from different schools: their code are recoded as -1", "\n")
  cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
  cat("remove", sum(nffr), "non-found friends (and not error code) : their code are recoded as -3", "\n")
  
  # dependent variable
  data[["y"]] <- data[[outcome]]
  mislistmis  <- c(55555555, 99999999, 99959995) # within school error codes
  
  # Network 
  G           <- vector("list", M)
  # Number of unmatched links
  nmatch      <- vector("list", M)
  # In this loop, G and nmatch will be filled
  for (m in 1:M) {
    cat("Outcome: ", outcome, " ** School: ", m, "/", M, "\n", sep = "")
    datam          <- data %>% filter(SCID == Sch[m])
    Nm             <- nrow(datam)
    Gm             <- matrix(0, Nm, Nm)
    nmatchm        <- numeric() 
    for (s in 1:Nm) {
      idx          <- which(datam[["AID"]] %in% unlist(datam[f_coln][s,]))
      Gm[s, idx]   <- 1
      nmatchm[s]   <- length(which(datam[f_coln][s,] %in% c(mislistmis, -3))) # Number of missing links
    }
    diag(Gm)       <- 0
    G[[m]]         <- Gm
    nmatch[[m]]    <- nmatchm
  }
  data             <- data %>% select("AID", "SSCHLCDE", "SCID", "y", starts_with("S45"), 
                                      all_of(exovar))
  data[["match"]]  <- unlist(lapply(G, rowSums))
  data[["nmatch"]] <- unlist(nmatch)
  save(data, G, exovar, file = paste0(OutDataPath, "/", outcome, ".Rda")) # Saving data
}

# Data summary
results <- list()
for(outcome in depvar){
  load(paste0(OutDataPath, "/", outcome, ".Rda"))
  # data <- data %>%  filter(match != 0 | nmatch <= 0) Remove fake isolated
  exovar  <- c("age", "female", "grade", "hispanic", "racewhite", 
               "raceblack", "raceasian", "melhigh", "memhigh", "memiss", 
               "mjprof", "mjother", "mjmiss","match","nmatch")
  YX       <- cbind(data$y, as.matrix(data[,exovar]))
  namesXY  <- c("Outcome", "Age", "Female", "Grade", "Hispanic", 
                "\\quad White", "\\quad Black", "\\quad Asian",
                "\\quad $<$ High", "\\quad $>$ High", "\\quad Missing",
                "\\quad Professional", "\\quad Other", "\\quad Missing",
                "\\quad Matched", "\\quad Unmatched", 
                "$n$ (students)", "$M$ (schools)")
  XY <- data.frame(Id   = 1:length(namesXY),
                   Var  = namesXY,
                   Mean = c(round(colMeans(YX, na.rm = TRUE), 3), nrow(data), 
                            length(unique(data$SSCHLCDE))),
                   SD   = c(round(apply(YX, 2, sd, na.rm = TRUE), 3), NA, NA),
                   Min  = c(round(apply(YX, 2, min, na.rm = TRUE), 3), NA, NA),
                   Max  = c(round(apply(YX, 2, max, na.rm = TRUE), 3), NA, NA))
  
  results[[outcome]] <- XY
}

stat_des <- results[[1]] %>% 
  left_join(results[[2]] %>% select(-"Var"), by = "Id") %>%
  left_join(results[[3]] %>% select(-"Var"), by = "Id") %>%
  left_join(results[[4]] %>% select(-"Var"), by = "Id") %>% select(-"Id")

wb       <- createWorkbook()
addWorksheet(wb, "Summary")

## First header row (outcome names)
header1   <- c("")
for(outcome in depvar)
  header1 <- c(header1, outcome, "", "", "")

header2 <- c("")
for(i in seq_along(depvar))
  header2 <- c(header2, "Mean", "SD", "Min", "Max")

## Write headers
writeData(wb, "Summary", t(header1), startRow = 1, colNames = FALSE)
writeData(wb, "Summary", t(header2), startRow = 2, colNames = FALSE)

## Merge outcome headers
for(i in seq_along(depvar)){
  first_col <- 2 + (i - 1) * 4
  mergeCells(wb, "Summary",
             cols = first_col:(first_col + 3),
             rows = 1)
}

## Write data
writeData(wb, "Summary",
          stat_des,
          startRow = 3,
          colNames = FALSE)

saveWorkbook(wb, paste0(OutResPath,"/","stats.des.xlsx"), overwrite = TRUE)
