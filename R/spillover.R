#' @title Find the optimal targeting order in an intervention
#' @description
#' `spillover` identifies the optimal targeting order in an intervention aiming at 
#' maximizing the utilitarian welfare and computes the associated individual spillovers. 
#' The targeting order and the spillovers can be generated via different targeting methods:
#' using the symmetric specification or the asymmetric specification. In the latter case,
#' the optimal order is computed by forward or backward optimization.
#' 
#' @param asymodel An object of class \code{\link[Asypeer]{asypeer.estim}} or 
#' \code{\link[Asypeer]{summary.asypeer.estim}}, 
#' estimated using the asymmetric specification (\code{asymmetry=TRUE}).
#' @param symodel An object of class \code{\link[Asypeer]{asypeer.estim}} or 
#' \code{\link[Asypeer]{summary.asypeer.estim}},
#'  estimated using the symmetric specification (\code{asymmetry=FALSE}).
#' 
#' @param Glist The adjacency matrix or list of adjacency matrices. For networks 
#'   consisting of multiple subnets (e.g., schools), \code{Glist} must be a list, 
#'   where the \code{s}-th element is an \eqn{n_s x n_s} adjacency matrix and 
#'   \eqn{n_s} is the number of nodes in subnet \code{s}.
#'   
#' @param targ.net A scalar indicating the index in \code{Glist} of the network
#' for which the targeting order and spillovers must be computed.
#' 
#' @param data An optional data frame, list, or environment (or an object that can 
#'   be coerced to a data frame via \link[base]{as.data.frame}) containing the 
#'   variables in the model. If a variable is not found in \code{data}, it is 
#'   retrieved from the environment of \code{asymodel}.
#'   
#' @param treatment A scalar, or a vector whose length equals the 
#' number of agents in the targeted network, indicating the value of the treatment.
#' 
#' @param nthread The number of CPU cores (threads) used for parallel computation.
#' 
#' @param print A logical value indicating whether a progress bar should be printed. 
#' 
#' @param tol A numeric tolerance used to assess the convergence, post-intervention,
#'  to a Nash equilibrium
#'   
#' @return A list containing:
#'     \item{targeted}{A matrix or list containing the indices 
#'     of the intervention's targets, ordered by their influence.}
#'     \item{budget}{A matrix indicating the intervention's budget for any number 
#'     of targeted agents.}
#'     \item{diff.y}{A matrix indicating the difference between the baseline welfare 
#'     and the welfare post-intervention, for any number of targeted agents.}
#'     \item{spillover}{A matrix giving the spillover associated
#'     with each targeted agent.}
#' @examples
#' \donttest{
#' if (requireNamespace("PartialNetwork", quietly = TRUE)) {
#'   library(PartialNetwork)
#'   ngr  <- 50  # Number of subnets
#'   nvec <- rep(30, ngr)  # Size of subnets
#'   n    <- sum(nvec)
#'   
#'   ### Simulating Data
#'   ## Network matrix
#'   G   <- lapply(1:ngr, function(z) {
#'     Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
#'     diag(Gz) <- 0
#'     # Adding isolated nodes (important for the structural model)
#'     niso <- sample(0:nvec[z], 1, prob = ((nvec[z] + 1):1)^5 / sum(((nvec[z] + 1):1)^5))
#'     if (niso > 0) {
#'       Gz[sample(1:nvec[z], niso), ] <- 0
#'     }
#'     Gz
#'   })
#'   
#'   Gnorm   <- norm.network(G)
#'   X       <- cbind(rnorm(n, 0, 2), rpois(n, 2))
#'   GX      <- peer.avg(Gnorm, X)
#'   delta   <- 0.25
#'   beta    <- c(0.3, 1.6)
#'   gamma   <- c(4, 1, -0.7, 0, -0.5)
#'   eps     <- rnorm(n, 0, 0.5)
#'   
#'   ## Generating `y`
#'   y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta = delta, beta = beta,
#'                    gamma = gamma, epsilon = eps)
#'   y <- y$y
#'   
#'   ### Estimating a symmetric peer effects model
#'   mod1 <- asypeer.estim(formula = y ~ X + GX, Glist = Gnorm, spillover = TRUE,
#'                         asymmetry = FALSE)
#'   summary(mod1)
#'   
#'   ### Estimating an asymmetric peer effects model
#'   mod2 <- asypeer.estim(formula = y ~ X + GX, Glist = Gnorm, spillover = TRUE)
#'   summary(mod2)
#'   
#'   ## treatment
#'   treat <- spillover(asymodel = mod2, symodel = mod1, Glist = Gnorm, targ.net = 1,
#'                      treatment = 1)
#'   
#'   ## Plot
#'   plot(treat)
#' }}
#' @export
#' @importFrom stats as.formula
spillover <- function(asymodel,
                      symodel,
                      Glist,
                      targ.net,
                      data,
                      treatment = 1,
                      nthread = 1,
                      print = TRUE,
                      tol = 1e-12) {
  
  tp        <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread <- tp
  }
  
  # Check models
  stopifnot(inherits(asymodel, c("asypeer.estim", "summary.asypeer.estim")))
  if (asymodel$model.info$asymmetry == FALSE) {
    stop("`asymodel` is not an asymmetric model.")
  }
  
  sym <- FALSE
  if (!missing(symodel)) {
    if (!is.null(symodel)) {
      sym <- TRUE
      stopifnot(inherits(symodel, c("asypeer.estim", "summary.asypeer.estim")))
      if (symodel$model.info$asymmetry == TRUE) {
        stop("`symodel` is not an asymmetric model.")
      }
      if (asymodel$model.info$yname != symodel$model.info$yname) {
        stop("The dependent variable is not the same in `asymodel` and `symodel`.")
      }
    }
  }
  
  # Targeted network
  targ.net <- unlist(targ.net)
  if(length(targ.net) > 1) {
    stop("Only a single network can be targed.")
  }
  
  # Network
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist  <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list.")
    }
  }
  
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  cumsn    <- c(0, cumsum(nvec))
  
  if (targ.net > M) {
    stop("`Glist` is a list of ", M, " network(s). Thus `targ.net` cannot be larger.")
  }
  Glist    <- Glist[[targ.net]]
  n        <- nrow(Glist)
  idpeer   <- lapply(1:n, \(i) which(Glist[i,] > tol) - 1)
  degree   <- rowSums(Glist)
  isolates <- sapply(idpeer, \(s) ifelse(length(s) == 0, 1, 0))
  
  # treatment
  if (length(treatment == 1)) {
    treatment <- rep(treatment, n)
  }
  if (length(treatment) != n) {
    stop("`treatment` must be either a scalar or a vector whose length equals the number of individuals in the treated network.")
  }
  
  # Outcome
  form     <- as.formula(paste0(asymodel$model.info$yname, " ~ 1"))
  if (missing(data)) {
    data   <- environment(asymodel)
  }
  y        <- formula2data(formula = form, data = data, fixed.effects = FALSE,
                           simulations = FALSE)$y
  if (length(y) != sum(nvec)) {
    stop("The size of `y` does not match `Glist`.")
  }
  y        <- y[(cumsn[targ.net] + 1):cumsn[targ.net + 1]]
  
  # Parameters
  betal    <- asymodel$gmm$Estimate["betal"]
  betah    <- asymodel$gmm$Estimate["betah"]
  delta    <- 0
  if (asymodel$model.info$spillover) {
    delta  <- asymodel$gmm$Estimate["delta"]
  }
  betdel   <- c(betal, betah, delta)
  if ((abs((delta + betal) / (1 + betal)) >= 1) | (abs((delta + betah) / (1 + betah)) >= 1)) {
    warning("The total peer effect is outside (-1, 1) for the asymmetric model. ",
            "The best response dynamics may not converge.")
  }
      
  # Estimate alpha
  alpha    <- falpha(betadelta = betdel, G = Glist, y = y, isolates = isolates)
  
  # Spillover
  # Asymmetric model
  out      <- list()
  seed     <- round(runif(1, 0, 1e9))
  RASym    <- fRankASym(y = y, alpha = alpha, betadelta = betdel,
                        G = Glist, idpeer = idpeer, treat = treatment, 
                        d = degree, tol = tol, nthread = nthread, seed = seed,
                        print = print)
  
  # Targeted index
  targeted <- list(asym.forw  = RASym$Rank[,1] + 1,
                   asym.backw = RASym$Rank[,2] + 1)
  
  streat   <- data.frame(asym.forw  = RASym$sum.treat[,1],
                         asym.backw = RASym$sum.treat[,2])
  
  diff.y   <- data.frame(asym.forw = RASym$Diff.sumy[,1],
                         asym.backw = RASym$Diff.sumy[,2])
  
  spill    <- data.frame(asym.forw = RASym$Spillover[,1] - 1,
                         asym.backw = RASym$Spillover[,2] - 1)
  
  # Symmetric model if provided
  RSym     <- NULL
  if (sym) {
    sbet   <- symodel$gmm$Estimate["beta"]
    sdel   <- 0
    if (symodel$model.info$spillover) {
      sdel <- symodel$gmm$Estimate["delta"]
    }
    if (abs((sdel + sbet) / (1 + sbet)) >= 1) {
      warning("The total peer effect is outside (-1, 1) for the symmetric model. ",
              "The best response dynamics may not converge.")
    }
    
    RSym   <- fRankSym(beta = sbet, delta = sdel, betadelta = betdel,
                       y = y, alpha = alpha, G = Glist, idpeer = idpeer,
                       treat = treatment, d = degree, isolates = isolates, 
                       tol = tol, nthread = nthread)
    
    targeted$sym <- lapply(RSym$setRank, \(k) k + 1)
    streat$sym   <- RSym$sum.treat
    diff.y$sym   <- RSym$Diff.sumy
    spill$sym    <- RSym$Spillover - 1
  }
  
  out        <- list(targeted  = targeted, 
                     budget    = streat,
                     diff.y    = diff.y, 
                     spillover = spill)
  class(out) <- "spillover"
  return(out)
  
}

#' @title Plotting the effect of a targeted intervention
#' @description
#'  `plot.spillover` creates plots that illustrate the spillovers across different
#'  intervention's budgets and targeting methods.
#' @param x An object of class \code{\link{spillover}} as returned by the function \link{spillover}.
#' @param metric A character string specifying which plots should be generated.
#' Available options include \code{"spillover"}, \code{"gain"}, \code{"outcome"}
#' and any combination of these options.
#' @param ... Further arguments passed to or from other methods.
#' @export
#' @importFrom graphics par layout legend lines plot.new
plot.spillover <- function(x, metric = "all", ...) {
  
  # Graph
  metric <- unique(tolower(metric))
  if (length(metric) == 1) {
    if (metric == "all") {
      metric <- c("spillover", "gain", "outcome")
    }
  }
  if (!all(metric %in% c("spillover", "gain", "outcome"))) {
    stop("Expected metrics include 'spillover', 'gain', and 'outcome'.")
  }
  
  # Plot
  plotspillover <- "spillover" %in% metric
  plotgain <- "gain" %in% metric
  plotoutcome   <- "outcome" %in% metric

  ## layout 
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  if (length(metric) == 1) {
    layout(matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE), 
           heights = c(4, 0.8))
    par(mar = c(5, 5, 3, 2))
  } else if (length(metric) == 2) {
    layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE), 
           heights = c(4, 0.8))
    par(mar = c(4, 5, 3, 1))
  } else {
    layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, ncol = 3, byrow = TRUE), 
           heights = c(4, 0.6))
    par(mar = c(4, 4, 3, 1))
  }
  
  # plots
  iplot <- c("spillover" = 1, "gain" = 2, "outcome" = 3)
  for (k in 1:length(metric)) {
    fplot(x = x, which = iplot[metric[k]])
  }
  
  # Legend
  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new() 
  if (any(c("spillover", "outcome") %in% metric)) {
    legend("bottom",
           xpd = NA,
           horiz = FALSE,
           bty = "n",
           legend = c("Symmetric", "Asymmetric (Forward)", "Asymmetric (Backward)"),
           col = c("red", "#2a9", "blue"),
           lty = 1)
  } else {
    legend("bottom",
           xpd = NA,
           horiz = FALSE,
           bty = "n",
           legend = c("Asymmetric (Forward)", "Asymmetric (Backward)"),
           col = c("#2a9", "blue"),
           lty = 1)
  }

  invisible(NULL)
}


fplot <- function(x, which) {
  
  if (which == 1) {
    
    plot(x$budget[,"sym"], x$spillover$sym * 100, type = "l",
         xlab = "Budget", ylab = "Spillover (%)", col = "red", 
         ylim = c(min(unlist(x$spillover)), max(unlist(x$spillover))) * 100)
    lines(x$budget[,"asym.backw"], x$spillover$asym.backw * 100, col = "blue")
    lines(x$budget[,"asym.forw"], x$spillover$asym.forw * 100, col = "#2a9")
    
  } else if (which == 2) {
    
    tp   <- as.data.frame(apply(x$spillover, 2, \(s) ifelse(abs(s) < 1e-8, 0, s)))
    gain <- (tp[,1:2] / tp[,3] - 1) * 100
    
    plot(x$budget[,"sym"], gain$asym.backw, type = "l",
         xlab = "Budget", ylab = "Gain (%)", col = "blue", 
         ylim = c(min(unlist(gain), na.rm = TRUE), 
                  max(unlist(gain), na.rm = TRUE)))
    lines(x$budget[,"asym.backw"], gain$asym.forw, col = "#2a9")
    
  } else {
    
    plot(x$budget[,"sym"], x$diff.y$sym, type = "l",
         xlab = "Budget", ylab = "Outcome increase", col = "red", 
         ylim = c(min(unlist(x$diff.y)), max(unlist(x$diff.y))))
    lines(x$budget[,"asym.backw"], x$diff.y$asym.backw, col = "blue")
    lines(x$budget[,"asym.forw"], x$diff.y$asym.forw, col = "#2a9")
    
  }
  
}

