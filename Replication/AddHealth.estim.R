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

outcome <- "drink"
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
y       <- data$y + 1 
Z    <- fitted(lm(y ~ X + GX + as.factor(data$SCID)))
plotcespeer(formula = y ~ X, instrument = ~ z, Glist = G, 
            grid.rho = seq(-20, 20, 0.5))
est <- cesconfpeer(formula = y ~ X, instrument = ~ zces, Glist = G)
summary(est)

Zces[Zces <= 0] <- min(Zces[Zces > 0]) # Because z cannot be <= 0 for the CES-based model

summary(cesconfpeer(formula = y ~ X + GX, 
        instrument = ~ Zces, 
        Glist = G, 
        fixed.effects = TRUE,
        drop = drop, # only isolated students
        radius = 5,
        grid.rho = seq(-100, 100, 1)))

GGX <- peer.avg(G, GX)
genpeer(formula = y ~ X + GX, 
        excluded.instruments = ~ GGX,
        Glist = G,
        fixed.effects = TRUE,
        structural = TRUE)

asypeer.estim(formula = y ~ X + GX, 
              Glist = G, 
              estimator = "glm", 
              power = 1,
              nfold = 5, 
              nthread = 1, 
              fixed.effects = TRUE,
              asymmetry = FALSE,
              drop = match == 0)



data_iso<-data%>%filter(match==0)
ols_iso <- lm(paste0("y ~ ", paste0(c(exovar), collapse = " + ")),data=data_iso)
summary(ols_iso)



