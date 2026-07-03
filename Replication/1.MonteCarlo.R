###########################################################################
################# Houndetoungan and Lambotte (202x) #######################
######## Asymmetries in Peer Effects for Continuous Outcomes ##############
######################## Monte Carlo Simulations ##########################
###########################################################################

# Last update: 2026-06-30

rm(list = ls())

library(AsyPeer)
library(PartialNetwork)
library(openxlsx)
library(dplyr)

OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED" # Where results should be saved

### Sample size
ngr   <- 50  # Number of subnets
nvec  <- rep(50, ngr)  # Size of subnets
n     <- sum(nvec)
cumsn <- c(0, cumsum(nvec))

### variance of eps
sigma <- 1

### This function simulates the data and estimates the model
# beta takes the conformity parameter value
# segregated is a Boolean indicating whether the network is segregated
festim  <- function(beta, segregated){
  
  ### gamma
  gamma <- c(1.4, -0.8, 0.7, -0.5)
  if (segregated) {
    gamma[3:4] <- 0 # No contextual effects for the segregated network
  }
  
  ### Possible degree values (support)
  dvalues  <- 0:10  # Agents can have up to 10 friends
  
  ### Degree distribution
  ddvalues <- c(0.22175143, 0.09047220, 0.10325461, 0.11262459, 0.11128805, 0.10039670,
                0.08578010, 0.07272753, 0.05633362, 0.03411014, 0.01126104)  
  
  ### Fixed effects, correlated with X1
  eff   <- rep(rnorm(ngr, 10, 1), nvec)
  
  ### X and Network simulation
  G     <- list()
  X     <- cbind(NA, rpois(n, 2))
  if (segregated) {
    
    # X1 is simulated to reflect the segregated grouping structure
    for (m in 1:ngr) {
      effm <- eff[(cumsn[m] + 1):cumsn[m + 1]]
      X[(cumsn[m] + 1):cumsn[m + 1], 1] <- 
        c(runif(20, effm[1:20] - 2, effm[1:20] + 2), 
          runif(30, effm[21:50] + 10, effm[21:50] + 20))
    }
    
    for (m in 1:ngr) {
      
      Gm    <- matrix(0, nvec[m], nvec[m])
      
      # The three groups
      n1    <- 10
      n2    <- 10
      n3    <- nvec[m] - n1 - n2
      gr1   <- 1:n1
      gr2   <- (n1 + 1):(n1 + n2)
      gr3   <- (n1 + n2 + 1):nvec[m]
      
      # Nodes from gr2 choose friends in gr3; these links are not reciprocal
      deg1  <- sample(dvalues, n1, replace = TRUE, prob = ddvalues)
      ifr1  <- lapply(1:n1, \(i) sample(gr2, deg1[i]))
      
      for (i in 1:n1) {
        Gm[gr1[i], ifr1[[i]]] <- Gm[ifr1[[i]], gr1[i]] <- 1
      }
      
      # Nodes from gr1 choose friends in gr2; these links are reciprocal
      deg2  <- sample(dvalues, n2, replace = TRUE, prob = ddvalues)
      ifr2  <- lapply(1:n2, \(i) {
        sample(gr3, deg2[i])
      })
      
      for (i in 1:n2) {
        Gm[gr2[i], ifr2[[i]]] <- 1
      }
      
      G[[m]] <- Gm
      
    }
    
  } else {
    
    # X
    X[,1]    <- runif(n, eff - 2, eff + 2) 
    
    # Degree simulation
    degree   <- lapply(1:ngr, \(s) sample(dvalues, nvec[s], replace = TRUE, 
                                          prob = ddvalues))
    
    # Peer group for each agent in each subnetwork
    peerg    <- lapply(1:ngr, \(s){
      lapply(1:nvec[s], function(i) sample((1:nvec[s])[-i], degree[[s]][i]))
    })    
    
    # Network
    G        <- list()
    for (m in 1:ngr) {
      Gm     <- matrix(0, nvec[m], nvec[m])
      for (i in 1:nvec[m]) {
        Gm[i, peerg[[m]][[i]]] <- 1
      }
      G[[m]] <- Gm
    }
    
  }
  
  ### GX 
  Gnorm <- norm.network(G) # Row-normalized network
  GX    <- peer.avg(Gnorm, X)
  
  ### Error term
  eps   <- rnorm(n, 0, sigma)
  
  ### Generating y
  y     <- asypeer.sim(formula = ~ -1 + eff + X + GX, Glist = Gnorm, delta = 0, 
                       beta = beta, gamma = c(1, gamma), epsilon = eps)$y
  
  ### Arguments for instrument generation (gen.inst)
  gen.inst.arg <- list(estimator = "rf",
                       full = TRUE,
                       nfold = 5,
                       num.trees = 1000,
                       mtry = 2 * ncol(X),
                       replace = FALSE,
                       sample.fraction = 0.8,
                       num.threads = 8,
                       power = 2,
                       seed = as.integer(runif(1, 0, 1e9)))
  
  # # using xgboost to generate instruments
  # gen.inst.arg <- list(estimator = "xgboost",
  #                      full = TRUE,
  #                      nfold = 5,
  #                      nrounds = 100,
  #                      nthread = 8,
  #                      subsample = 0.8,
  #                      power = 2,
  #                      seed = as.integer(runif(1, 0, 1e9)))
  
  fm    <- NULL
  if (segregated) {
    fm  <- y ~ X # No contextual effects for segregated networks
  } else {
    fm  <- y ~ X + GX
  }
  
  ### Estimation of the asymmetric model
  est1  <- asypeer.estim(fm, Glist = Gnorm, fixed.effects = TRUE, 
                         asymmetry = TRUE, gen.inst.arg = gen.inst.arg)
  
  ### Estimation of the symmetric model
  est2  <- asypeer.estim(fm, Glist = Gnorm, fixed.effects = TRUE, 
                         asymmetry = FALSE, gen.inst.arg = gen.inst.arg)
  
  ### Output
  coef1 <- est1$gmm$Estimate
  names(coef1) <- paste0("ASYM_", names(coef1))
  coef2 <- est2$gmm$Estimate
  names(coef2) <- paste0("SYM_", names(coef2))
  print(c(coef1, coef2))
  c(coef1, coef2)
}

### Number of simulations
nsim   <- 1000

### This function simulates and estimates nsim times for each specification
run_sim <- function(beta, segregated) {
  do.call(cbind, lapply(1:nsim, \(i) {
    message(sprintf("Iteration: %d", i))
    festim(beta = beta, segregated = segregated)
  }))
}

### DGP 1
set.seed(123)
beta1 <- c(1, 1)
print("DGP 1 - random")
Est1  <- run_sim(beta = beta1, segregated = FALSE)
print("DGP 1 - segregated")
Est1H <- run_sim(beta1, segregated = TRUE)

### DGP 2
set.seed(123)
beta2 <- c(0.4, 2.6)
print("DGP 2 - random")
Est2  <- run_sim(beta2, segregated = FALSE)
print("DGP 2 - segregated")
Est2H <- run_sim(beta2, segregated = TRUE)

### DGP 3
set.seed(123)
beta3 <- c(2.6, 0.4)
print("DGP 3 - random")
Est3  <- run_sim(beta3, segregated = FALSE)
print("DGP 3 - segregated")
Est3H <- run_sim(beta3, segregated = TRUE)

### DGP 4
set.seed(123)
beta4 <- c(-0.4, 2.6)
print("DGP 4 - random")
Est4  <- run_sim(beta4, segregated = FALSE)
print("DGP 4 - segregated")
Est4H <- run_sim(beta4, segregated = TRUE)

save(Est1, Est1H, Est2, Est2H, Est3, Est3H, Est4, Est4H, 
     file = paste0(OutResPath, "/Simulations.Rda"))

### Export in an Excel file
load(paste0(OutResPath, "/Simulations.Rda"))

# Summary function
Sumfunc    <- function(x) {
  c(mean = mean(x), sd = sd(x))
}

# Gather results to construct tables to be exported
Est <- lapply(1:4, \(k) {
  # True values of beta
  b    <- get(paste0("beta", k))
  bl   <- b[1]
  bh   <- b[2]
  
  # Compute mean and sd
  tp   <- as.data.frame(t(apply(get(paste0("Est", k)), 1, Sumfunc)))
  tpH  <- as.data.frame(t(apply(get(paste0("Est", k, "H")), 1, Sumfunc)))
  
  # Construct tables
  outk <- data.frame(Model = c("betal", "betah", "beta",  
                               paste0("gamma1_", 1:2),
                               paste0("gamma2_", 1:2))) %>% 
    ## Adding columns for asymmetry and random network
    left_join(data.frame(Model = c("betal", "betah", 
                                   paste0("gamma1_", 1:2),
                                   paste0("gamma2_", 1:2)),
                         tp["ASY" == substr(rownames(tp), 1, 3),]) %>%
                rename(mean_asym = mean, sd_asym = sd), by = "Model") %>%
    ## Adding columns for symmetry and random network
    left_join(data.frame(Model = c("beta", 
                                   paste0("gamma1_", 1:2),
                                   paste0("gamma2_", 1:2)),
                         tp["SYM" == substr(rownames(tp), 1, 3),]) %>%
                rename(mean_sym = mean, sd_sym = sd), by = "Model") %>% 
    # Adding columns for asymmetry and segregated network
    left_join(data.frame(Model = c("betal", "betah", 
                                   paste0("gamma1_", 1:2)),
                         tpH["ASY" == substr(rownames(tpH), 1, 3),]) %>%
                rename(mean_asymH = mean, sd_asymH = sd), by = "Model") %>%
    # Adding columns for symmetry and segregated network
    left_join(data.frame(Model = c("beta", 
                                   paste0("gamma1_", 1:2)),
                         tpH["SYM" == substr(rownames(tpH), 1, 3),]) %>%
                rename(mean_symH = mean, sd_symH = sd), by = "Model") %>%
    # Formatting
    mutate(across(starts_with("sd_"), 
                  ~ ifelse(is.na(.), NA, paste0("(", sprintf("%.3f", .), ")"))),
           across(starts_with("mean_"), 
                  ~ ifelse(is.na(.), NA, sprintf("%.3f", .))))
  
  addrow <- data.frame(matrix(NA, 1, ncol(outk), dimnames =list(NULL, colnames(outk))))
  outk   <- addrow %>% bind_rows(outk) 
  outk$Model <- c(paste0("DGP ", k,
                         ": $\\beta^l = ", sprintf("%.2f", bl),
                         "$, $\\beta^h = ", sprintf("%.2f", bh)),
                  "$\\beta^l$", "$\\beta^h$", "$\\beta$",
                  "$\\gamma_{1,1}$", "$\\gamma_{1,2}$",
                  "$\\gamma_{2,1}$", "$\\gamma_{2,2}$")
  outk
})

Est <- Est[[1]] %>% bind_rows(Est[[2]]) %>% 
  bind_rows(Est[[3]]) %>% bind_rows(Est[[4]])


# Create an Excel file
wb       <- createWorkbook()

# Table 1 for the paper
shname   <- "Table1"
rstart   <- 3
addWorksheet(wb, shname)
writeData(wb, shname, Est, 
          keepNA = TRUE, na.string = "", startRow = rstart)

# First row
val      <- c("Random Network", "Segregated Network")
pos      <- c(2, 6)
for (k in 1:length(pos)) {
  writeData(wb, shname, val[k], startCol = pos[k], startRow = rstart - 1)
  mergeCells(wb, shname, cols = pos[k]:(pos[k] + 3), rows = rstart - 1)
}

# Second row
val      <- rep(c("Asymmetry", "Symmetry"), 2)
pos      <- seq(2, 8, 2)
for (k in 1:length(pos)) {
  writeData(wb, shname, val[k], startCol = pos[k], startRow = rstart)
  mergeCells(wb, shname, cols = pos[k]:(pos[k] + 1), rows = rstart)
}

# DGP details
pos      <- seq(1, 25, 8)
for (k in 1:length(pos)) {
  mergeCells(wb, shname, cols = 1:9, rows = rstart + pos[k])
}

# Export
saveWorkbook(wb, paste0(OutResPath, "/Simulations.xlsx"), overwrite = TRUE)

