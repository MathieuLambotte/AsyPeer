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
  
OutDataPath <- "D:/Home/mlambotte/Documents/Dropbox/AsymmetricPeer/AsyPeerCode/SOULD_BE_DELETED_LATER/"
OutResPath <- "D:/Home/mlambotte/Documents/Dropbox/AsymmetricPeer/AsyPeerCode/Application/"
# List of outcome variables
depvar  <- c( "smoke", "drink", "risky")

outcome <- "smoke"
cat("Outcome: ", outcome, "\n", sep = "")
#smoke 1.92 et 2.65
#drink 1.03 et 0.89
#risky 0.97 1.0; 1.55
#risky
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
#rm(list = "data") 
#invisible(gc())

########################################################
###################### Estimation ######################
########################################################
## Estimation assuming symmetry
est_sym    <- asypeer.estim(formula = y ~ X + GX, 
                            Glist = G, 
                            estimator = "glm", 
                            power = 2,
                            nfold = 5, 
                            fixed.effects = TRUE,
                            asymmetry = FALSE,
                            HAC = "hetero",
                            drop = drop)
summary(est_sym, diagnostic = TRUE, KPtest = FALSE)


est_asy    <- asypeer.estim(formula = y ~ X + GX, 
                            Glist = G, 
                            estimator = "glm", 
                            power = 2,
                            nfold = 5, 
                            fixed.effects = TRUE,
                            spillover=FALSE,
                            asymmetry = TRUE,
                            HAC = "group-iid",
                            drop = drop)
summary(est_asy, diagnostic = TRUE, KPtest = FALSE)



## Other tests
if(any(data$y<=0)){ # Because y cannot be <= 0 for the CES-based model
  y <- data$y +1
} else {
    y <- data$y
    }
z    <- fitted(lm(y ~ X + GX + as.factor(data$SCID)))
z[z <= 0] <- min(z[z > 0]) # Because z cannot be <= 0 for the CES-based model

# Test with different weight.rho
# The plot is completely different
plotcespeer(formula = y ~ X+GX, instrument = ~ z, Glist = G, fixed.effect=TRUE,
            grid.rho = seq(-100, 100, 1),drop=drop, nthread = 7,
            weight.rho = 5)

est <- cesconfpeer(formula = y ~ X+GX, instrument = ~ z, Glist = G,
                   interval = c(-50, 50),fixed.effect = TRUE,drop=drop,
                   nthread = 7, n.rho0 = 6, n.beta0 = 1, weight.rho = 5)
summary(est)

