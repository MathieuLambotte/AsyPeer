###########################################################################
################# Houndetoungan and Lambotte (202x) #######################
######## Asymmetries in Peer Effects for Continuous Outcomes ##############
######################## Monte Carlo Simulations ##########################
###########################################################################

# Last update: 2026-05-20

set.seed(123)
rm(list = ls())

library(PartialNetwork)
library(openxlsx)
library(dplyr)
library(parallel)
library(AsyPeer)

OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED/" # Where results should be saved
OutResPath  <- "~/Dropbox/Academy/1.Papers/AsymmetricPeer/AsyPeerCode/Replication"
OutResPath  <- "D:/Home/mlambotte/Documents/Dropbox/AsymmetricPeer/AsyPeerCode/package/AsyPeer/Replication"
ngr  <- 50  # Number of subnets
nvec <- rep(50, ngr)  # Size of subnets
n    <- sum(nvec)

### Simulating Data
#fixed parameters values
gamma   <- c(2, -0.5, 1, -0.2, 0.6)
sigma <- 1 #variance of eps
delta <- 0

## Network matrix
# Network degree distribution: fitting the empirical distribution of AddHealth
# Degree possible values (support):
dvalues    <- 0:10  # Agents can have up to 10 friends
# Degree distribution:
ddvalues   <- c(0.22175143, 0.09047220, 0.10325461, 0.11262459, 0.11128805, 0.10039670,
                0.08578010, 0.07272753, 0.05633362, 0.03411014, 0.01126104)  # Degree distribution

festim<-function(beta, fixed.effects = TRUE, HAC="group-iid", power = 3, nfold = 2, weight = "IV"){
  X       <- cbind(runif(n, 0, 5), rpois(n, 2))
  eps     <- rnorm(n, 0, sigma)
  
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

avg <- peer.asyavg(formula = ~y, Glist=Gnorm)
ybar <- avg$peer.avg
ycheck <- avg$hpeer.avg-avg$hdegree*y

est1 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                     estimator = "ols", power = power, nfold = nfold, nthread = 5,
                     weight = weight, spillover = FALSE, asymmetry = TRUE)

est2 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "glm", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = TRUE)

est3 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "rf", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = TRUE)

est4 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "ols", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = FALSE)

est5  <- asypeer.estim(y ~ X + GX, excluded.instruments = ~ ybar+ycheck,
                      Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "ols", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = TRUE)

coef1<-summary(est1)$coefficients[,1]
names(coef1)<-paste0(names(coef1),"_LM")
coef2<-summary(est2)$coefficients[,1]
names(coef2)<-paste0(names(coef2),"_GLM")
coef3<-summary(est3)$coefficients[,1]
names(coef3)<-paste0(names(coef3),"_RF")
coef4<-summary(est4)$coefficients[,1]
names(coef4)<-paste0(names(coef4),"_SYM")
coef5<-summary(est5)$coefficients[,1]
names(coef5)<-paste0(names(coef5),"_OLS")

c(coef1,coef2,coef3,coef4,coef5)
}


festim_homo<-function(beta, fixed.effects = TRUE, HAC="group-iid", power = 3, nfold = 2, weight = "IV"){
  X       <- cbind(runif(n, 0, 5), rpois(n, 2))
  eps     <- rnorm(n, 0, sigma)
  
degree   <- lapply(1:ngr, function(s) sample(dvalues, nvec[s], replace = TRUE, prob = ddvalues))


#homophilous network formation
#parameters of the network formation model
rho <- c(-2, 0.5, -0.2)
# compute distance between individuals
tmp <- c(0, cumsum(nvec))
X1l <- lapply(1:ngr, function(x) X[c(tmp[x] + 1):tmp[x+1],1])
X2l <- lapply(1:ngr, function(x) X[c(tmp[x] + 1):tmp[x+1],2])
dist.net <- function(x, y) abs(x - y)
X1.mat <- lapply(1:ngr, function(m) {
  matrix(kronecker(X1l[[m]], X1l[[m]], FUN = dist.net), nvec[m])})
X2.mat <- lapply(1:ngr, function(m) {
  matrix(kronecker(X2l[[m]], X2l[[m]], FUN = dist.net), nvec[m])})
# true network
Xnet <- as.matrix(cbind("Const" = 1,
                        "dX1" = mat.to.vec(X1.mat),
                        "dX2" = mat.to.vec(X2.mat)))
U <- Xnet %*% rho + rlogis(n)
Umat <- vec.to.mat(U, nvec, normalise = FALSE)

G <- lapply(1:ngr, function(g) {
  Ug <- Umat[[g]]
  ng <- nvec[g]
  dg <- degree[[g]]   # vector of target degrees for group g
  
  Gg <- matrix(0, ng, ng)
  
  for (i in 1:ng) {
    # exclude self-links
    Ui <- Ug[i, ]
    Ui[i] <- -Inf
    
    # pick top d_i neighbors
    idx <- order(Ui, decreasing = TRUE)[0:dg[i]]
    Gg[i, idx] <- 1
  }
  Gg
})

Gnorm   <- norm.network(G)
GX      <- peer.avg(Gnorm, X) 
## Generating `y`
y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta=delta,beta = beta,
                 gamma = gamma, epsilon = eps, nthread = 5)
y <- y$y

avg <- peer.asyavg(formula = ~y, Glist=Gnorm)
ybar <- avg$peer.avg
ycheck <- avg$hpeer.avg-avg$hdegree*y

est1 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                     estimator = "ols", power = power, nfold = nfold, nthread = 5,
                     weight = weight, spillover = FALSE, asymmetry = TRUE)

est2 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "glm", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = TRUE)

est3 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "rf", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = TRUE)

est4 <- asypeer.estim(y ~ X + GX, Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                      estimator = "ols", power = power, nfold = nfold, nthread = 5,
                      weight = weight, spillover = FALSE, asymmetry = FALSE)

est5  <- asypeer.estim(y ~ X + GX, excluded.instruments = ~ ybar+ycheck,
                       Glist = Gnorm, fixed.effects = fixed.effects, HAC=HAC,
                       estimator = "ols", power = power, nfold = nfold, nthread = 5,
                       weight = weight, spillover = FALSE, asymmetry = TRUE)

coef1<-summary(est1)$coefficients[,1]
names(coef1)<-paste0(names(coef1),"_OLS")
coef2<-summary(est2)$coefficients[,1]
names(coef2)<-paste0(names(coef2),"_GLM")
coef3<-summary(est3)$coefficients[,1]
names(coef3)<-paste0(names(coef3),"_RF")
coef4<-summary(est4)$coefficients[,1]
names(coef4)<-paste0(names(coef4),"_SYM")
coef5<-summary(est5)$coefficients[,1]
names(coef5)<-paste0(names(coef5),"_OLS")

c(coef1, coef2, coef3, coef4, coef5)
}

# Number of simulations
nsim      <- 1e3
RNGkind("L'Ecuyer-CMRG")

cl <- makeCluster(5) # Use one less than the number of cores

# Export the function to the cluster
clusterExport(cl, varlist = c("festim","festim_homo","ngr","dvalues","ddvalues","nvec","n","gamma","sigma","delta"))

# Load necessary libraries in the cluster
clusterEvalQ(cl, {
  library(PartialNetwork)
  library(AsyPeer)
})
#DGP A
betaA=c(1,1)
clusterExport(cl, varlist = c("betaA"))
EstA <- parLapply(cl, 1:nsim, function(i) festim(beta=betaA))

#DGP B
betaB=c(2,1)
clusterExport(cl, varlist = c("betaB"))
EstB <- parLapply(cl, 1:nsim, function(i) festim(beta=betaB))

#DGP C
betaC=c(1,2)
clusterExport(cl, varlist = c("betaC"))
EstC <- parLapply(cl, 1:nsim, function(i) festim(beta=betaC))

#DGP D
betaD=c(1,-0.4)
clusterExport(cl, varlist = c("betaD"))
EstD <- parLapply(cl, 1:nsim, function(i) festim(beta=betaD))



save(EstA,EstB,EstC,EstD,
     file = paste0(OutResPath, "/Simulations.Rda"))


# Summary
# Summary function
Sumfunc    <- function(x) {
  c(mean = mean(x), sd = sd(x))
}

Est  <- as.data.frame(rbind(apply(do.call(rbind, EstA), 2, Sumfunc),
                                  apply(do.call(rbind, EstB), 2, Sumfunc),
                                  apply(do.call(rbind, EstC), 2, Sumfunc),
                                  apply(do.call(rbind, EstD), 2, Sumfunc)))

rownames(Est)<-c("mean_DGPA","sd_DGPA","mean_DGPB","sd_DGPB","mean_DGPC","sd_DGPC",
                 "mean_DGPD","sd_DGPD")
Est <- Est %>% mutate(across(everything(), ~ sprintf("%.3f", .))) %>%  # Ensures exactly 3 decimals
  mutate(across(everything(), ~ ifelse(row_number() %% 2 == 0, paste0("(", trimws(.), ")"), .)))  # Parentheses without spaces 

Est <- Est %>%
  tibble::rownames_to_column("row") %>%
  mutate(
    stat = ifelse(str_detect(row, "mean"), "mean", "sd"),
    DGP  = case_when(
      str_detect(row, "DGPA") ~ "DGP A",
      str_detect(row, "DGPB") ~ "DGP B",
      str_detect(row, "DGPC") ~ "DGP C",
      str_detect(row, "DGPD") ~ "DGP D"
    )
  ) %>%
  select(-row) %>%
  pivot_longer(-c(stat, DGP), names_to = "param", values_to = "value") %>%
  arrange(param, stat) %>%
  pivot_wider(names_from = DGP, values_from = value)

colnames(Est)<-c("stat","param",paste0("DGP A: $\\boldsymbol\\beta^l = (", paste0(betaA[1], collapse = ", "), ")$, ", 
         "$\\beta^h = ", betaA[2], "$"),
  paste0("DGP B: $\\boldsymbol\\beta^l = (", paste0(betaB[1], collapse = ", "), ")$, ", 
         "$\\beta^h = ", betaB[2], "$"),
  paste0("DGP C: $\\boldsymbol\\beta^l = (", paste0(betaC[1], collapse = ", "), ")$, ", 
         "$\\beta^h = ", betaC[2], "$"),
  paste0("DGP D: $\\boldsymbol\\beta^l = (", paste0(betaD[1], collapse = ", "), ")$, ", 
         "$\\beta^h = ", betaD[2], "$"))
# Export to Excel
# Create the workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "Simu results")

writeData(wb, "Simu results", Est, keepNA = TRUE, na.string = "", startRow = 1, startCol = 1)

# Save the workbook
saveWorkbook(wb, paste0(OutResPath, "/simulations.xlsx"), overwrite = TRUE)
