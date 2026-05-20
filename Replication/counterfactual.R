##############################################################################################################
##############################################################################################################
#################### Asymmetric Peer Effect Models by A. Houndetoungan and M. Lambotte #######################
##############################################################################################################
# This script replicates our counterfactual exercices using Add Health data
# Please, prior use the script AddHealth.data.R to prepare the data set for
# each outcome in am expect form for this scrip.
# This scripts assumes that the prepared data are saved in the folder with path `OutDataPath`

# Last update: 2026-03-17

rm(list = ls())

library(PartialNetwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(AsyPeer)
library(openxlsx)

#Simulate data
ngr  <- 50  # Number of subnets
nvec <- rep(100, ngr)  # Size of subnets
n    <- sum(nvec)

### Simulating Data
#fixed parameters values
gamma   <- c(2, -0.5, 1, -0.2, 0.6)
sigma <- 1 #variance of eps
delta <- 0
beta<- c(2,1) #two parameters even if symmetric
X       <- cbind(runif(n, 0, 5), rpois(n, 2))
eps     <- rnorm(n, 0, sigma)
## Network matrix
# Network degree distribution: fitting the empirical distribution of AddHealth
# Degree possible values (support):
dvalues    <- 0:10  # Agents can have up to 10 friends
# Degree distribution:
ddvalues   <- c(0.22175143, 0.09047220, 0.10325461, 0.11262459, 0.11128805, 0.10039670,
                0.08578010, 0.07272753, 0.05633362, 0.03411014, 0.01126104)
degree   <- lapply(1:ngr, function(s) sample(dvalues, nvec[s], replace = TRUE, prob = ddvalues))

# Peer group for each agent in each subnetwork
peerg    <- lapply(1:ngr, function(s){
  lapply(1:nvec[s], function(i) sample((1:nvec[s])[-i], degree[[s]][i]))
})    

# Network matrix
G        <- list()
for (s in 1:ngr) {
  Gs     <- matrix(0, nvec[s], nvec[s])
  for (i in 1:nvec[s]) {
    Gs[i, peerg[[s]][[i]]] <- 1
  }
  G[[s]] <- Gs
}

Gnorm   <- norm.network(G)
GX      <- peer.avg(Gnorm, X) 
## Generating `y`
y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta=delta,beta = beta,
                 gamma = gamma, epsilon = eps, nthread = 5)
y <- y$y
## Estimation assuming symmetry
sym <- asypeer.estim(formula = y ~ X + GX, 
                     Glist = Gnorm, 
                     power = 2,
                     fixed.effects = TRUE,
                     asymmetry = FALSE,
                     HAC = "hetero")

## Estimation assuming asymmetry
asym <- asypeer.estim(formula = y ~ X + GX, 
                      Glist =Gnorm, 
                      estimator = "glm", 
                      power = 2,
                      nfold = 5, 
                      fixed.effects = TRUE,
                      spillover = FALSE,
                      asymmetry = TRUE,
                      HAC = "hetero")

########################################################
################### Targeted policy ##################
########################################################
#list_school <- split(seq_len(nrow(data)), data[,"SCID"])
cs          <- c(0, cumsum(nvec))
list_school <- lapply(1:ngr, function(i) (cs[i] + 1):cs[i + 1])
sch_id      <- 2

## Symmetric model
beta        <- sym$gmm$Estimate["beta"]
target_y0   <- y[list_school[[sch_id]]]
target_X    <- X[list_school[[sch_id]],]
target_GX   <- GX[list_school[[sch_id]],]
target_G    <- Gnorm[[sch_id]]
target_n    <- length(target_y0)
target_niso <- rowSums(target_G) > 0
I           <- diag(target_n)

# Order in which agents should be targeted 
D         <- diag(ifelse(target_niso, 1 / (1 + beta), 1))   
target_or <- order(D %*% solve(I - beta * t(target_G) / (1 + beta), rep(1, target_n)), 
                   decreasing = TRUE)

# Spillover
## For the symmetric model
E         <- matrix(0, target_n, target_n)
for (i in 1:target_n) {
  E[i, head(target_or, i)] <- ifelse(target_niso[target_or[1:i]], 1 / (1 + beta), 1)
}
spill_sym <- c(E %*% solve(I - beta * t(target_G) / (1 + beta), rep(1, target_n)))


#######################
## Asymmetric model
# we need the residuals at the baseline
target_avg <- peer.asyavg(formula = ~target_y0, Glist=target_G)
target_ybar <- target_avg$peer.avg 
target_ycheck <- target_avg$hpeer.avg-target_avg$hdegree*target_y0 
# predicted y at baseline
target_Xtemp<-target_X
target_Xtemp[target_niso,]<-target_Xtemp[target_niso,]/(1+asym$gmm$Estimate["betal"])
target_GXtemp<-target_GX
target_GXtemp[target_niso,]<-target_GXtemp[target_niso,]/(1+asym$gmm$Estimate["betal"])
target_yhat<-cbind(target_ybar,target_ycheck,target_Xtemp,target_GXtemp)%*%
  (c(asym$gmm$Estimate["betal"]/(1+asym$gmm$Estimate["betal"]),
     (asym$gmm$Estimate["betah"]-asym$gmm$Estimate["betal"])/(1+asym$gmm$Estimate["betal"]),
     asym$gmm$Estimate[-c(1,2)])) #eq 3.3 paper
target_res <-  (target_y0-target_yhat) 
target_res[target_niso]=target_res[target_niso]*(1+asym$gmm$Estimate["betal"])

spill_asym<-cbind(numeric(target_n),numeric(target_n))
for(i in 1:target_n){
  shock<-numeric(target_n)
  shock[head(target_or, i)]=1
  spill_asym[i,1]<-i
  spill_asym[i,2]<-sum(asypeer.sim(formula = ~  -1 + shock + target_X + target_GX, Glist = target_G, delta=0,
                                       beta = asym$gmm$Estimate[c(1,2)], gamma = c(1,asym$gmm$Estimate[-c(1,2)]), 
                                      epsilon = target_res, nthread = 5)$y-target_y0)
}

df_plot <- data.frame(x = spill_asym[,1],
  sym = diff(c(0, spill_sym)),
  asym = diff(c(0, spill_asym[,2])))

df_long <- reshape2::melt(df_plot, id.vars = "x",
                          variable.name = "model",
                          value.name = "spillover")
ggplot(df_long, aes(x = x, y = spillover, color = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("sym" = "#1f77b4", "asym" = "#d62728"),
    labels = c("sym" = "Symmetric", "asym" = "Asymmetric")) +
  labs(x = "Number of targeted agents",
    y = "Marginal Spillover",
    color = "Model") +
  theme_minimal() 

