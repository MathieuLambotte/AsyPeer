###########################################################################
################# Houndetoungan and Lambotte (202x) #######################
######## Asymmetries in Peer Effects for Continuous Outcomes ##############
########################### Targeting policy ##############################
###########################################################################
# This script replicates our empirical results using Add Health data.
# Before using this script, use first the script '2.AddHealth.data.R' to 
# prepare the data set, and '3.AddHealth.Estimation.R' to estimate the 
# This script assumes that the estimated models are saved in the folder with 
# path `OutResPath`

# Last update: 2026-07-25

library(ggplot2)
library(ggh4x)
library(tidyr)
library(dplyr)
library(AsyPeer)
library(PartialNetwork)
set.seed(123)

OutDataPath <- "PATH/TO/WHERE/PREPARED/DATA/WILL/BE/SAVED"
OutResPath  <- "PATH/TO/WHERE/RESULTS/WILL/BE/SAVED"

# List of outcome variables
depvar  <- c("smoke", "fight", "optimism", "drink")

# Remove fake isolated
rmfiso  <- FALSE

# This function computes spillover for each outcome
spillovers <- function(outcome){
  
  cat("Outcome: ", outcome, "\n", sep = "")
  
  ### Data and result loading
  load(paste0(OutDataPath, "/", outcome, ".Rda"))
  load(paste0(OutResPath, "/", outcome, 
              ifelse(rmfiso, ".noFakeIso", ""), ".Rda"))
  
  ## Variables and networks
  exovar  <- c("age", "female", "grade", "hispanic", "racewhite", 
               "raceblack", "raceasian", "melhigh", "memhigh", "memiss", 
               "mjprof", "mjother", "mjmiss", "match", "nmatch")
  
  G       <- norm.network(G)
  y       <- data$y
  
  # Find 2 networks around the average size
  schools <- unique(data$SSCHLCDE)
  sizes   <- sapply(schools, \(m) sum(data$SSCHLCDE == m))
  avgsize <- mean(sizes)
  pick    <- which(schools %in% sample(schools[abs(sizes - avgsize) <= 100], 2))    
  schools <- schools[pick]
  
  # Estimation
  spill1  <- spillover(asymodel = sasym,
                       symodel = ssym,
                       Glist = G,
                       targ.net = pick[1],
                       treatment = ifelse(outcome == "optimism", 1, -1),
                       nthread = 25,
                       print = TRUE,
                       tol = 1e-9) 
  
  spill2   <- spillover(asymodel = sasym,
                        symodel = ssym,
                        Glist = G,
                        targ.net = pick[2],
                        treatment = ifelse(outcome == "optimism", 1, -1),
                        nthread = 25,
                        print = TRUE,
                        tol = 1e-9) 
  
  save(spill1, spill2, schools, 
       file = paste0(OutResPath, "/", outcome, 
                     ifelse(rmfiso, ".noFakeIso", ""), ".spillover.Rda"))
  invisible(NULL)
}

lapply(depvar, spillovers)

# This function creates the data set for ggplot
data_plot <- function(outcome){
  
  ## Load estimation
  load(paste0(OutResPath, "/", outcome, 
              ifelse(rmfiso, ".noFakeIso", ""), ".spillover.Rda"))
  
  ## Round spillover to zero when it is lower that 1e-8
  ## Because it lead to loss = -1000%. This is only due to numerical precisions
  spill1$spillover <- as.data.frame(spill1$spillover) %>% 
    mutate_all(~ifelse(abs(.) < 1e-8, 0, .))
  spill2$spillover <- as.data.frame(spill2$spillover) %>% 
    mutate_all(~ifelse(abs(.) < 1e-8, 0, .))
  
  ## Compute loss with respect to symmetric
  loss1 <- data.frame(fwd = 1 -  spill1$spillover[,"sym"] / spill1$spillover[,"asym.forw"],
                      bwd = 1 -  spill1$spillover[,"sym"] / spill1$spillover[,"asym.backw"])
  
  loss2 <- data.frame(fwd = 1 -  spill2$spillover[,"sym"] / spill2$spillover[,"asym.forw"],
                      bwd = 1 -  spill2$spillover[,"sym"] / spill2$spillover[,"asym.backw"])
  
  ## data school 1
  n             <- length(spill1$budget)
  spillover_df1 <- data.frame(Budget = rep(spill1$budget, 5), 
                              Series = c(rep("Symmetric", n),
                                         rep("Asymmetric (Backward)", n),
                                         rep("Asymmetric (Forward)", n),
                                         rep("Asymmetric (Backward)", n),
                                         rep("Asymmetric (Forward)", n)),
                              Metric = c(rep("Spillover", 3 * n),
                                         rep("Loss", 2 * n)),
                              Value = c(spill1$spillover$sym, 
                                        spill1$spillover$asym.backw,
                                        spill1$spillover$asym.forw,
                                        loss1$bwd,
                                        loss1$fwd) * 100,
                              SCHID = schools[1])
  
  ## data school 2
  n             <- length(spill2$budget)
  spillover_df2 <- data.frame(Budget = rep(spill2$budget, 5), 
                              Series = c(rep("Symmetric", n),
                                         rep("Asymmetric (Backward)", n),
                                         rep("Asymmetric (Forward)", n),
                                         rep("Asymmetric (Backward)", n),
                                         rep("Asymmetric (Forward)", n)),
                              Metric = c(rep("Spillover", 3 * n),
                                         rep("Loss", 2 * n)),
                              Value = c(spill2$spillover$sym, 
                                        spill2$spillover$asym.backw,
                                        spill2$spillover$asym.forw,
                                        loss2$bwd,
                                        loss2$fwd) * 100,
                              SCHID = schools[2])
  
  spillover_df1 %>% bind_rows(spillover_df2) %>%
    mutate(Outcome = outcome,
           Metric = factor(Metric, levels = c("Spillover", "Loss")),
           Value = ifelse(is.na(Value), 0, Value)) # because 0 spillover / 0 spillover = NA
}

df <- do.call(rbind, lapply(depvar, data_plot)) %>% 
  ## change outcome in factor
  mutate(Outcome = factor(Outcome, levels = c("smoke", "fight", "optimism", "drink"),
                          label = c("Smoke", "Fight","Optimism","Drink"))) %>%
  ## School index
  group_by(Outcome) %>%
  mutate(SchoolSlot = factor(as.numeric(as.factor(SCHID)),
                             labels = c("School 1", "School 2"))) %>%
  ungroup() %>%
  mutate(FacetCol = factor(paste(SchoolSlot, Metric, sep = "_"),
                           levels = c("School 1_Spillover", "School 1_Loss",
                                      "School 2_Spillover", "School 2_Loss")))

### Data set for school label
schid_labels <- df %>%
  group_by(Outcome, FacetCol, SCHID) %>%
  summarise(
    x = 0.54 * max(Budget, na.rm = TRUE),
    .groups = "drop"
  ) 

### Plot
(graph <- ggplot(df, aes(Budget, Value, colour = Series)) +
    geom_line(linewidth = 0.8) +
    facet_grid2(rows = vars(Outcome),
                cols = vars(FacetCol),
                scales = "free",
                independent = "all",
                labeller = labeller(FacetCol = as_labeller(function(x) {
                  ifelse(grepl("Spillover", x), "Panel A: Spillover (%)", "Panel B: Symmetric Loss (%)")}))) +
    geom_text(data = schid_labels,
              aes(x = x, y = Inf, label = paste("School:", SCHID)),
              hjust = -0.1, vjust = 1.5,
              size = 5, colour = "#334455", fontface = "italic",
              family = "Palatino",
              inherit.aes = FALSE) +
    scale_colour_manual(values = c("Symmetric" = "red",
                                   "Asymmetric (Forward)" = "#22AA99",
                                   "Asymmetric (Backward)" = "blue"),
                        breaks = c("Symmetric",
                                   "Asymmetric (Forward)",
                                   "Asymmetric (Backward)")) +
    labs(x = expression("Budget (" * kappa * ")"), y = NULL) +
    theme_bw() +
    
    theme(
      legend.position = "bottom",
      text = element_text(family = "Palatino"),
      legend.title = element_blank(),
      legend.text = element_text(size = 14, margin = margin(r = 20)),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
)


ggsave(paste0("spillovers", ifelse(rmfiso, ".noFakeIso", ""), ".pdf"), 
       path = OutResPath, plot = graph, device = "pdf", 
       width = 15, height = 9)

## Separate figures for outcome
## This is used in our slides

for (k in 1:4) {
  dfk <- df %>% filter(Outcome == levels(Outcome)[k]) %>%
    mutate(school = factor(paste0("School ", SCHID)))
  
  (graph <- ggplot(dfk, aes(Budget, Value, colour = Series)) +
      geom_line(linewidth = 0.8) +
      facet_grid2(rows = vars(school),
                  cols = vars(Metric),
                  scales = "free",
                  independent = "all",
                  labeller = labeller(Metric = as_labeller(function(x) {
                    ifelse(grepl("Spillover", x), "Panel A: Spillover (%)", "Panel B: Symmetric Loss (%)")}))) +
      scale_colour_manual(values = c("Symmetric" = "red",
                                     "Asymmetric (Forward)" = "#22AA99",
                                     "Asymmetric (Backward)" = "blue"),
                          breaks = c("Symmetric",
                                     "Asymmetric (Forward)",
                                     "Asymmetric (Backward)"))+
      labs(y = NULL) +
      theme_bw() +
      
      theme(
        legend.position = "bottom",
        text = element_text(family = "Palatino"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, margin = margin(r = 20)),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  )
  
  
  ggsave(paste0("spillovers.", depvar[k], ifelse(rmfiso, ".noFakeIso", ""), ".pdf"), 
         path = OutResPath, plot = graph, device = "pdf", width = 7.6, height = 4.6)
}
