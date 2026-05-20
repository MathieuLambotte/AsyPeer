##############################################################################################################
##############################################################################################################
#################### Asymmetric Peer Effect Models by A. Houndetoungan and M. Lambotte #######################
##############################################################################################################
# This script replicates the targeted policy simulated in Section 2.3.

# Last update: 2026-02-12
rm(list = ls())
library(AsyPeer)

# Network
G   <- matrix(c(0, rep(1/2, 2), rep(0, 6)), 3, byrow = TRUE)
rownames(G) <- colnames(G) <- c("i", "j", "k")

# beta
beta    <- c(1.5, 0.5)

# alpha
alpha   <- c(1.2, 1.5, 1)
eps     <- rep(0, 3)

# Outcome before the intervention
y1 <- asypeer.sim(formula = ~ -1 + alpha, Glist = G, delta = 0, beta = beta,
                  gamma = 1, epsilon = eps)$y
round(y1, 1) # j is high and k is low

# Increasing alphaj
alpha2    <- alpha
alpha2[2] <- alpha2[2] + 1
y2 <- asypeer.sim(formula = ~ -1 + alpha2, Glist = G, delta = 0, beta = beta,
                  gamma = 1, epsilon = eps)$y
round(y2, 1) # j is high and k is low
round(sum(y2 - y1), 1) # 1.1. There are social multiplier effects

# Increasing alphak
alpha3    <- alpha
alpha3[3] <- alpha3[3] + 1
y3 <- asypeer.sim(formula = ~ -1 + alpha3, Glist = G, delta = 0, beta = beta,
                  gamma = 1, epsilon = eps)$y
round(y3, 1) # j and k are high
round(sum(y3 - y1), 1) # 1.2: There are social multiplier effects but less efficient than increasing j

# Increasing alpha1
alpha4    <- alpha
alpha4[1] <- alpha4[1] + 1
y4 <- asypeer.sim(formula = ~ -1 + alpha4, Glist = G, delta = 0, beta = beta,
                  gamma = 1, epsilon = eps)$y
round(y4, 1) # j and k are low
round(sum(y4 - y1), 1) # 0.5: There are no social multiplier effects
