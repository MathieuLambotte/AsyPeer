##############################################################################################################
##############################################################################################################
#################### Asymmetric Peer Effect Models by A. Houndetoungan and M. Lambotte #######################
##############################################################################################################


rm(list = ls())

library(PartialNetwork)
library(dplyr)
library(AsyPeer)

OutDataPath <- "PATH/TO/WHERE/PREPARED/DATA/IS/SAVED/" # Where prepared data for each outcome are saved (/ at the end is important)
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED/" # Where results should be saved

OutDataPath <- "~/Dropbox/Academy/1.Papers/AsymmetricPeer/AsyPeerCode/SOULD_BE_DELETED_LATER/"
OutResPath  <- "~/Dropbox/Academy/1.Papers/AsymmetricPeer/AsyPeerCode/Application"
  
# List of outcome variables
depvar  <- c( "smoke", "drink", "risky")

outcome <- "smoke"
cat("Outcome: ", outcome, "\n", sep = "")

########################################################
################### Data Preparation ###################
########################################################

load(file = paste0(OutDataPath, outcome, ".Rda"))

match   <- data$match
nmatch  <- data$nmatch
y       <- data$y
X       <- as.matrix(data[,exovar])
GX      <- peer.avg(norm.network(G), data[,exovar])
drop    <- (match == 0) & (nmatch > 0)
G       <- norm.network(G)
rm(list = "data")
invisible(gc())

########################################################
###################### Estimation ######################
########################################################
## Estimation without asymetry
est_sym    <- asypeer.estim(formula = y ~ X + GX, 
                            Glist = G, 
                            estimator = "glm", 
                            power = 3,
                            nfold = 5, 
                            nthread = 1, 
                            fixed.effects = TRUE,
                            asymmetry = FALSE)

summary(est_sym)


est_asy    <- asypeer.estim(formula = y ~ X + GX, 
                            Glist = G, 
                            estimator = "glm", 
                            power = 3,
                            nfold = 5, 
                            nthread = 1, 
                            fixed.effects = TRUE,
                            asymmetry = TRUE)
summary(est_asy)

est_sym    <- asypeer.estim(form, Glist = norm.network(G), estimator = "glm", power = 3,
                        nfold = 5, nthread = 1, data = data, fixed.effects = T,
                        spillover = FALSE, asymmetry = FALSE)
summary(est_sym)

data_iso<-data%>%filter(match==0)
ols_iso <- lm(paste0("y ~ ", paste0(c(exovar), collapse = " + ")),data=data_iso)
summary(ols_iso)



