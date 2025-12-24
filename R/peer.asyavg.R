#' @title Computing peer (asymmetric) average values
#' 
#' @param formula An object of class \link[stats]{formula}, which should be specified as \code{y ~ X}, 
#' where `y` is the outcome used to characterized which friends have a higher or lower outcome and `X` is the vector 
#'   or matrix of variables whose average values among peers should be computed.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), 
#'   `Glist` must be a list of subnets, with the \code{m}-th element being an \eqn{n_m \times n_m} 
#'   adjacency matrix, where \eqn{n_m} is the number of nodes in the \code{m}-th subnet.  
#' @param data An optional data frame, list, or environment (or an object that can be coerced to a 
#'   data frame via \link[base]{as.data.frame}) containing the variables in the model. If a variable 
#'   is not found in `data`, it is taken from \code{environment(formula)}, typically the environment 
#'   from which `peer.asyavg` is called.
#' @param nthread A strictly positive integer specifying the number of threads used in 
#'   computationally intensive steps.
#'   
#' @description
#' `peer.asyavg` computes the average values of a matrix `X` among friends, 
#' differentiating between friends with a higher or lower outcome than the agent.
#' 
#' @examples
#' if (requireNamespace("PartialNetwork", quietly = TRUE)) {
#' library(PartialNetwork)
#' ngr  <- 50  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#' ### Simulating Data
#' ## Network matrix
#' G <- lapply(1:ngr, function(z) {
#'   Gz       <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
#'   diag(Gz) <- 0
#'   # Adding isolated nodes (important for the structural model)
#'   niso <- sample(0:nvec[z], 1, prob = ((nvec[z] + 1):1)^5 / sum(((nvec[z] + 1):1)^5))
#'   if (niso > 0) {
#'     Gz[sample(1:nvec[z], niso), ] <- 0
#'   }
#'   Gz
#' })
#' Gnorm   <- norm.network(G)
#' X       <- cbind(rnorm(n, 0, 2), rpois(n, 2))
#' GX      <- peer.avg(Gnorm, X)
#' colnames(X)<-c("X1","X2") 
#' delta   <- 0.25
#' beta    <- c(0.3, 0.6)
#' gamma   <- c(4, 1, -0.7, 0, -0.5) 
#' eps     <- rnorm(n, 0, 0.5) 
#' 
#' ## Generating `y`
#' y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta = delta, beta = beta, 
#'                 gamma = gamma, epsilon = eps, nthread = 5)
#' y <- y$y
#' 
#' ### Computing averages X among peers
#' peeravg <- peer.asyavg(formula = y ~ X, Glist = Gnorm, nthread = 1)
#' }
#' @export
peer.asyavg <- function(formula, Glist, data, nthread = 1) {
  
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  Glist      <- fGnormalise(Glist, nthread)
  
  formula    <- as.formula(formula)
  hasy       <- (length(as.formula(formula)) == 3)
  if (missing(data)) {
    data     <- env(formula)
  }
  f.t.data   <- formula2data(formula = formula, data = data, fixed.effects = TRUE,
                                simulations = !hasy) #Intercept is not necessary
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  
  S          <- length(Glist)
  nvec       <- unlist(lapply(Glist, ncol))
  
  ## Thread
  tp         <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread  <- tp
  }
  
  if (sum(nvec) != nrow(X)) {
    stop("Glist and V do not match")
  }
  cumsn      <- c(0, cumsum(nvec))
  tp         <- NULL
  if (hasy) {
    tp <- highlowstat2(y = y, X = X, G = Glist, cumsn = cumsn, nvec = nvec, 
                       ngroup = S, nthread = nthread)
  } else {
    tp <- highlowstat1(X = X, G = Glist, cumsn = cumsn, nvec = nvec, 
                       ngroup = S, nthread = nthread)
  }
  g     <- NULL
  pavg  <- NULL
  hg    <- NULL
  hpavg <- NULL
  lg    <- NULL
  lpavg <- NULL
  if (hasy) {
    g     <- c(tp$g)
    hg    <- c(tp$gh)
    lg    <- c(tp$gl)
    pavg  <- cbind(tp$ybar, tp$Xbar)
    hpavg <- cbind(tp$ybh, tp$Xbh)
    lpavg <- cbind(tp$ybl, tp$Xbl)
  } else {
    if (Kx == 1) {
      g     <- c(tp$g)
      hg    <- c(tp$gh)
      lg    <- c(tp$gl)
      pavg  <- c(tp$Xbar)
      hpavg <- c(tp$Xbh)
      lpavg <- c(tp$Xbl)
    } else {
      g     <- tp$g
      hg    <- tp$gh
      lg    <- tp$gl
      pavg  <- tp$Xbar
      hpavg <- tp$Xbh
      lpavg <- tp$Xbl
    }
  }
  
  list(degree     = g,
       peer.avg   = pavg,
       hdegree    = hg,
       hpeer.avg  = hpavg,
       ldegree    = lg,
       lpeer.avg  = lpavg)
}