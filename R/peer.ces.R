##################################
### ACKNOWLEDGMENT
### Taken from QuantilePeer and modified according to the MIT License; QuantilePeer is cited.
### Houndetoungan A (2025). QuantilePeer: Quantile Peer Effect Models. 
### doi:10.32614/CRAN.package.QuantilePeer, R package version 0.0.1,
##################################


#' @title Estimation of CES-Based Peer Effects Models With Conformity
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified as \code{y ~ x1 + x2}, 
#' where `y` is the outcome and `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param instrument An object of class \link[stats]{formula} indicating the excluded instrument. It should be specified as \code{~ z},  
#' where `z` is the excluded instrument for the outcome. Following Boucher et al. (2024), it can be an OLS exogenous prediction of `y`.  
#' This prediction is used to compute instruments for the CES function of peer outcomes.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet.
#' 
#' @param data An optional data frame, list, or environment (or an object that can be coerced by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `cesconfpeer` is called.
#' 
#' @param tol A tolerance value used in the QR factorization to identify columns of explanatory variable and instrument matrices that ensure a full-rank matrix (see the \link[base]{qr} function).
#' The same tolerance is also used in the to minimize the concentrated GMM objective function (see \link[stats]{optimise}).
#' 
#' @param drop A dummy vector of the same length as the sample, indicating whether an observation should be dropped. 
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user because they may still be friends with other agents.
#' 
#' @param fixed.effects A logical value or string specifying whether the model includes subnet fixed effects. 
#' 
#' @param HAC A character string specifying the correlation structure of the idiosyncratic errors 
#'   for covariance computation. Options are `"iid"` for independent errors; `"group iid"` for 
#'   independence within the groups of isolated and non-isolated players; `"hetero"` for 
#'   heteroskedastic but non-autocorrelated errors; and `"cluster"` for heteroskedastic errors with 
#'   potential within-subnetwork correlation.
#'   
#' @param set.rho A fixed value for the CES substitution parameter to estimate a constrained model.  
#' 
#' @param interval A numerical vector of two elements indicating support for \eqn{\rho}.  
#' 
#' @param nthread Number of CPU cores (threads) used to run parts of the estimation in parallel.
#' 
#' @param arg.optim A list of additional arguments passed to \link[stats]{optim},
#'   such as `method` and `control`, for models with a flexible \eqn{\rho}, or 
#'   to \link[stats]{optimize}, such as `tol`, for models with a fixed (or grid) 
#'   value of \eqn{\rho}. 
#'   For models with a flexible \eqn{\rho}, the objective function is concentrated 
#'   in a function of \eqn{\rho} and \eqn{\beta}, which is optimized using 
#'   \link[stats]{optim}. For models with a fixed \eqn{\rho}, the objective 
#'   function is concentrated in a function of \eqn{\beta} alone, which is 
#'   optimized using \link[stats]{optimize}.
#'
#' @param n.rho0,n.beta0 Integer values indicating the number of starting points 
#'   to test for \eqn{\rho} and \eqn{\beta}, respectively. Because the objective 
#'   function may exhibit local optima, testing multiple starting values is useful, 
#'   especially for \eqn{\rho}. On a reasonable parameter set, testing ten starting 
#'   values often suffices.
#' 
#' @param grid.rho A finite grid of values for the CES substitution parameter \eqn{\rho} (see Details).
#' This grid is used to plot the objective function after maximizing over the other parameters.
#' It is helpful for determining a reasonable interval for \eqn{\rho}.
#'
#' @param weight.rho The \eqn{\rho} value used to compute the weighting matrix for
#'   the first GMM estimation or to plot the objective function. The weight is 
#'   computed as proportional to \eqn{(Z'Z)^{-1}}, where \eqn{Z} is the matrix of 
#'   instruments. If `weight.rho` is omitted, the identity matrix is used.
#'   
#' @param weight.optim A Boolean indicating whether the GMM with the optimal 
#'   weighting matrix should be computed.
#'   
#' @param ... Further arguments passed to or from other methods.
#' 
#' @description
#' `cesconfpeer` estimates the CES-based peer effects model introduced by Boucher et al. (2024) with conformity only. See Details.
#' @details
#' Let \eqn{\mathcal{N}} denote a set of \eqn{n} agents indexed by \eqn{i \in [1, n]}.  
#' Agents are connected through a network represented by an adjacency matrix \eqn{\mathbf{G} = [g_{ij}]} of dimension \eqn{n \times n}, where \eqn{g_{ij} = 1} if agent \eqn{j} is a friend of agent \eqn{i}, and \eqn{g_{ij} = 0} otherwise.  
#' In weighted networks, \eqn{g_{ij}} can be a nonnegative value (not necessarily binary) that measures the intensity of the outgoing link from \eqn{i} to \eqn{j}. The model accommodates such networks.  
#' Note that the network may consist of multiple independent subnets (e.g., schools).  
#' The `Glist` argument is the list of subnets. For a single subnet, `Glist` should be a list containing one matrix.\cr
#'
#' The specification of the CES-based peer effects model differs for isolated and non-isolated individuals.  
#' For an **isolated** agent \eqn{i}, the specification is similar to a standard linear-in-means model without social interactions:
#' \deqn{y_i = \mathbf{x}_i^{\prime}\gamma + \varepsilon_i,}
#' where \eqn{\varepsilon_i} is an idiosyncratic error term and \eqn{\gamma} captures the effect of \eqn{\mathbf{x}_i} on \eqn{y_i}.\cr  
#'
#' If agent \eqn{i} is **non-isolated**, the specification is:
#' \deqn{y_i = \dfrac{\beta}{1 + \beta}\left(\sum_{j = 1}^n g_{ij}y_j^{\rho}\right)^{1/\rho} + (1 - \lambda_2)\mathbf{x}_i^{\prime}\dfrac{\gamma}{1 + \beta} + \varepsilon_i,}
#' where \eqn{\gamma} captures conformity effects—how agent \eqn{i} tends to conform to friends through the social norm \eqn{\left(\sum_{j = 1}^n g_{ij}y_j^{\rho}\right)^{1/\rho}}.\cr  
#'
#' The parameter \eqn{\rho} determines the form of the social norm:
#' - When \eqn{\rho > 1}, individuals are more sensitive to peers with high outcomes.  
#' - When \eqn{\rho < 1}, individuals are more sensitive to peers with low outcomes.  
#' - When \eqn{\rho = 1}, peer effects are uniform across peer outcome values.
#' 
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, including the number of subnets, the number of observations, and other key details.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#'     \item{first.gmm}{Initial GMM estimation results using the identity matrix as the weighting matrix.}
#' @examples 
#' \donttest{
#' if (requireNamespace("PartialNetwork", quietly = TRUE)) {
#' library(PartialNetwork)
#' set.seed(123)
#' ngr  <- 30  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#' ### Simulating Data
#' ## Network matrix
#' G   <- lapply(1:ngr, function(z) {
#'   Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
#'   diag(Gz) <- 0
#'   # Adding isolated nodes 
#'   niso <- sample(0:nvec[z], 1, prob = ((nvec[z] + 1):1)^5 / sum(((nvec[z] + 1):1)^5))
#'   if (niso > 0) {
#'     Gz[sample(1:nvec[z], niso), ] <- 0
#'   }
#'   Gz
#' })
#' 
#' Gnorm   <- norm.network(G)
#' X       <- cbind(rnorm(n, 0, 2), rpois(n, 2))
#' delta   <- 0
#' beta    <- c(1, 1)
#' gamma   <- c(6, 1, 0.7)
#' eps     <- rnorm(n, 0, 0.5)
#' 
#' ## Generating `y`
#' y <- asypeer.sim(formula = ~ X, Glist = Gnorm, delta = delta, beta = beta,
#'                  gamma = gamma, epsilon = eps)
#' y <- y$y
#' 
#' ### Estimating the asymmetric peer effects model
#' z   <- predict(lm(y ~ X))
#' plotcespeer(formula = y ~ X, instrument = ~ z, Glist = Gnorm, 
#'             grid.rho = seq(-20, 20, 1))
#' est <- cesconfpeer(formula = y ~ X, instrument = ~ z, Glist = Gnorm,
#'                    n.rho0 = 2, n.beta0 = 1)
#' summary(est)
#' }
#' }
#' @export
#' @importFrom stats optimise optim 
cesconfpeer <- function(formula, 
                        instrument, 
                        Glist, 
                        fixed.effects = FALSE,
                        set.rho       = NULL, 
                        interval      = c(-100, 100), 
                        tol           = 1e-8, 
                        drop          = NULL, 
                        HAC           = "group-iid", 
                        nthread       = 1,
                        data,
                        n.rho0        = 10,
                        n.beta0       = 5,
                        weight.rho    = NULL,
                        weight.optim  = TRUE,
                        arg.optim     = list(),
                        ...) {
  ## Thread
  tp        <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread <- tp
  }
  
  # set.rho
  if (!is.null(set.rho)) {
    if(set.rho == 0) {
      stop("`rho` cannot be zero.")
    }
  }
  
  ## Network
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist  <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  Glist      <- fGnormalise(Glist, nthread)
  
  ## HAC
  HACn    <- NULL
  if (tolower(HAC) %in% c("iid","i.i.d")) {
    HAC   <- "iid"
    HACn  <- 0
  } else if (tolower(HAC) %in% c("group-iid", "g-iid", "group iid", "giid", "groupiid")) {
    HAC   <- "group-iid"
    HACn  <- 1
  } else if (tolower(HAC) %in% c("hetero","het","heter","heteroskedasticity","heteroscedasticity")) {
    HAC   <- "hetero"
    HACn  <- 2
  } else if (tolower(HAC) %in% c("cluster", "clustered","clus")) {
    HAC   <- "cluster"
    HACn  <- 3
  } else {
    stop("This HAC option is not available.")
  }
  if (HACn == 2 & fixed.effects > 0) {
    HAC   <- "cluster"
    HACn  <- 3
  }
  
  ## sizes
  dg       <- fnetwork(Glist = Glist)
  S        <- dg$S
  ldg      <- dg$ldg
  SIs      <- dg$SIs
  SnIs     <- dg$SnIs
  nvec     <- dg$nvec
  n        <- dg$n
  cumsn    <- dg$cumsn
  idpeer   <- dg$idpeer
  Is       <- dg$Is
  lIs      <- dg$lIs
  nIs      <- dg$nIs
  lnIs     <- dg$lnIs
  n_iso    <- length(Is)
  n_niso   <- n - n_iso
  dg       <- dg$dg
  if (!all(dg %in% c(0, 1))) {
    stop("G is not row-normalized.")
  }
  
  ## Formula to data
  formula  <- as.formula(formula)
  f.t.data <- formula2data(formula = formula, data = data, fixed.effects = fixed.effects,
                           simulations = FALSE)
  
  # X, exogenous variables
  X      <- f.t.data$X
  xname  <- f.t.data$xname
  Kx     <- length(xname)
  
  # outcome
  y      <- f.t.data$y
  yname  <- f.t.data$yname
  
  # Instrument
  inst     <- as.formula(instrument); instrument <- inst
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format `~ z` (a single instrument).")
  f.t.data <- formula2data(formula = inst, data = data, fixed.effects = FALSE, 
                           simulations = TRUE)
  z       <- f.t.data$X
  zename  <- f.t.data$xname
  z       <- z[, zename != "(Intercept)"]
  zename  <- zename[zename != "(Intercept)"]
  if (length(zename) != 1) {
    stop("Only one excluded instrument can be used.")
  }
  
  # Create additional variables
  hasIs    <- fCESdatainit(y = y, z = z, G = Glist, nvec = nvec, S = S, 
                           ldg = ldg, lIs = lIs, lnIs = lnIs, drop = drop)
  frindex  <- hasIs$friendindex
  frzeroy  <- hasIs$frzeroy
  frzeroz  <- hasIs$frzeroz
  ldg_st   <- hasIs$ldg
  dg_st    <- hasIs$dg
  S_st     <- hasIs$S
  SIs_st   <- hasIs$SIs
  SnIs_st  <- hasIs$SnIs
  yFMiMa   <- cbind(hasIs$yFmin, hasIs$yFmax)
  zFMiMa   <- cbind(hasIs$zFmin, hasIs$zFmax)
  lIs      <- hasIs$lIs # After selection
  Is       <- hasIs$Is # After selection
  lnIs     <- hasIs$lnIs # After selection
  nIs      <- hasIs$nIs # After selection
  hasIs    <- hasIs$hasIs
  niso     <- length(Is)
  nniso    <- length(nIs)
  n_st     <- niso + nniso
  sel      <- sort(c(Is, nIs))
  if ((1 %in% frzeroy) | (1 %in% frzeroz)) {
    stop("The outcome `y` and the instrument `z` must be strictly positive.")
  }
  
  # idXiso and idXniso
  if (!is.null(weight.rho)) {
    if(weight.rho == 0) {
      stop("`rho` cannot be zero.")
    }
  }
  CESd     <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                       cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                       lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                       n = n, Kx = Kx, S = S, rho = ifnullset(weight.rho, 1), 
                       FE = fixed.effects, deriv = FALSE, nthread = nthread)
  
  idXiso   <- fcheckrank(X = CESd[Is + 1, 1:Kx, drop = FALSE], tol = tol)
  idXniso  <- fcheckrank(X = CESd[nIs + 1, 1:Kx, drop = FALSE], tol = tol)
  if (length(unique(c(idXiso, idXniso))) != Kx) {
    stop("The regressor matrix X is not full rank.")
  }
  if (niso > 0) {
    if (length(idXiso) == 0) {
      stop("The regressor matrix X is not full rank for isolated nodes")
    }
  }
  
  ## degree of freedom
  dfiso    <- NULL
  dfniso   <- NULL
  if (fixed.effects) {
    Kxiso  <- length(idXiso)
    dfiso  <- niso - SIs_st - Kxiso
    dfniso <- nniso - SnIs_st - Kx + Kxiso - 1 - is.null(set.rho)
  } else {
    Kxiso  <- length(idXiso)
    dfiso  <- niso - Kxiso
    dfniso <- nniso - Kx + Kxiso - 1 - is.null(set.rho)
  }
  if ((dfiso < 30) && HACn == 1) {
    HAC    <- "iid"
    HACn   <- 0
  }
  if (dfniso <= 0) {
    stop("Insufficient number of observations for non-isolated nodes.")
  }
  Kiso  <- length(idXiso)
  Kniso <- length(idXniso)
  
  # Interval for rho
  rhomin <- -100
  rhomax <- 100
  
  if (is.null(set.rho)) {
    if (!all(is.finite(interval))) {
      stop("`interval` is not finite.")
    }
    if (length(interval) != 2) {
      stop("`interval` must be a valid two-dimensional vector.")
    }
    if (interval[1] >= interval[2]) {
      stop("`interval` must have its first element smaller than its second.")
    }
    
    rhomin <- interval[1]
    rhomax <- interval[2]
  }
  
  # Additional arguments for optim or optimize
  optmethod  <- ifnullset(arg.optim$method, "Nelder-Mead")
  optcontrol <- ifnullset(arg.optim$control, list(maxit = 1e3))
  opttol     <- ifnullset(arg.optim$tol, .Machine$double.eps^0.25)
  
  # Optimization
  gmm1       <- NULL
  gmm2       <- NULL
  k          <- NULL  # for R CMD check   
  if (is.null(set.rho)) {
    
    cl       <- NULL
    on.exit({
      registerDoSEQ()
      try(stopCluster(cl), silent = TRUE)
    }, add = TRUE)
    
    cl       <- makeCluster(nthread)
    registerDoParallel(cl)
    
    # First estimation when W = I
    Kz       <- Kiso + Kniso + 2
    est1     <- NULL
    W        <- S * diag(Kz)
    if (!is.null(weight.rho)) {
      Z   <- matrix(0, n, Kz)
      Z[nIs + 1, 1] <- CESd[nIs + 1, Kx + 3]
      Z[nIs + 1, 2] <- CESd[nIs + 1, Kx + 4]
      if (Kiso > 0) { # if there are isolated
        Z[Is + 1, 3:(2 + Kiso)] <- CESd[Is + 1, idXiso + 1]
      }
      Z[nIs + 1, (3 + Kiso):(Kz)] <- CESd[nIs + 1, idXniso + 1]
      
      # first Optimization
      W     <- solve(crossprod(Z[sel + 1,]) / S) 
    }
    
    if (rhomax * rhomin < 0) { # estimation for positive rho and for negative rho
      
      n.rho0 <- ceiling(n.rho0 / 2) * 2
      rho0   <- c(seq(rhomin, -1e-4, length.out = n.rho0 / 2), 
                  seq(1e-4, rhomax, length.out = n.rho0/ 2))
      beta0  <- seq(-0.9, 0.9, length.out = n.beta0) / (1 - seq(-0.9, 0.9, length.out = n.beta0))
      init   <- cbind(rep(rho0, each = n.beta0), rep(beta0, n.rho0))
      
      est1   <- foreach(k         = 1:nrow(init),
                         .packages = "AsyPeer"
      ) %dorng% {
        # Optimization
        ARG   <- list(par = init[k,], fn = fCESobj, 
                      method = optmethod, control = optcontrol, 
                      X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                      cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                      lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, 
                      zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                      sel = sel,  FE = fixed.effects, 
                      rhomin = ifelse(init[k, 1] < 0, rhomin, 0),
                      rhomax = ifelse(init[k, 1] < 0, 0, rhomax))
        do.call(optim, ARG)
      }
      
    } else {
      
      rho0   <- seq(rhomin, rhomax, length.out = n.rho0)
      beta0  <- seq(-0.9, 0.9, length.out = n.beta0) / (1 - seq(-0.9, 0.9, length.out = n.beta0))
      init   <- cbind(rep(rho0, each = n.beta0), rep(beta0, n.rho0))
      
      est1   <- foreach(k         = 1:nrow(init),
                        .packages = "AsyPeer"
      ) %dorng% {
        # Optimization
        ARG   <- list(par = init[k,], fn = fCESobj, 
                      method = optmethod, control = optcontrol, 
                      X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                      cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                      lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, 
                      zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                      sel = sel,  FE = fixed.effects, 
                      rhomin = rhomin, rhomax = rhomax)
        do.call(optim, ARG)
      }
      
    }
    est1       <- est1[[which.min(sapply(est1, \(k) k$value))]]

    ## Full parameters
    est1$theta <- fCESparm(theta01 = est1$par, 
                           X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                           cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, 
                           lIso = lIs, lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, 
                           yFMiMa = yFMiMa, zFMiMa = zFMiMa, W = W, idXiso = idXiso, 
                           idXniso = idXniso, sel = sel, FE = fixed.effects)
  
    # Second estimation
    est2 <- est1
    if (weight.optim) {
      W    <- fCESWeight_2ins(theta = est1$theta,
                              X = X, y = y, z = z, G = Glist, friendindex = frindex,
                              cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs,
                              lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa,
                              zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                              sel = sel, HACn = HACn, FE = fixed.effects, dfiso = dfiso, 
                              dfniso = dfniso, nthread = nthread)
      
      if (rhomax * rhomin < 0) { # estimation for positive rho and for negative rho
        
        rho0   <- c(est1$par[1], # former value and we add many others
                    seq(rhomin, -1e-4, length.out = n.rho0 / 2), 
                    seq(1e-4, rhomax, length.out = n.rho0/ 2))
        beta0  <- est1$par[2] # only former value
        init   <- cbind(rho0, rep(beta0, n.rho0 + 1))
        
        est2   <- foreach(k         = 1:nrow(init),
                          .packages = "AsyPeer"
        ) %dorng% {
          # Optimization
          ARG   <- list(par = init[k,], fn = fCESobj, 
                        method = optmethod, control = optcontrol, 
                        X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                        cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                        lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, 
                        zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                        sel = sel,  FE = fixed.effects, 
                        rhomin = ifelse(init[k, 1] < 0, rhomin, 0),
                        rhomax = ifelse(init[k, 1] < 0, 0, rhomax))
          do.call(optim, ARG)
        }
        
      } else {
        
        rho0   <- c(est1$par[1], seq(rhomin, rhomax, length.out = n.rho0))
        beta0  <- est1$par[2]
        init   <- cbind(rho0, rep(beta0, n.rho0 + 1))
        
        est2   <- foreach(k         = 1:nrow(init),
                          .packages = "AsyPeer"
        ) %dorng% {
          # Optimization
          ARG   <- list(par = init[k,], fn = fCESobj, 
                        method = optmethod, control = optcontrol, 
                        X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                        cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                        lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, 
                        zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                        sel = sel,  FE = fixed.effects, 
                        rhomin = rhomin, rhomax = rhomax)
          do.call(optim, ARG)
        }
        
      }
      est2       <- est2[[which.min(sapply(est2, \(k) k$value))]]
      
      ## Compute all parameters
      est2$theta <- fCESparm(theta01 = est2$par, 
                             X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                             cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, 
                             lIso = lIs, lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, 
                             yFMiMa = yFMiMa, zFMiMa = zFMiMa, W = W, idXiso = idXiso, 
                             idXniso = idXniso, sel = sel, FE = fixed.effects)
    }
    
    ## Covariance matrix
    Cov <- fCESparmCov(theta = est2$theta,
                        X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                        cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                        lnIso = lnIs, Iso = Is, nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, 
                        zFMiMa = zFMiMa, W = W, idXiso = idXiso, idXniso = idXniso,
                        sel = sel, FE = fixed.effects, HACn = HACn, dfiso = dfiso, 
                        dfniso = dfniso, nthread = nthread)
  
    gmm1 <- list(Estimate = est1$theta,
                 optim    = c(list("objective" = est1$value),
                              est1[c("counts", "convergence", "message")]))
    gmm2 <- list(Estimate      = est2$theta,
                 cov           = Cov$cov,
                 sigma         = c(overall = Cov$serr, isolates = Cov$serr_iso, nonisolates = Cov$serr_niso),
                 testlinear    = c("stat" = Cov$Testrho, "p-value" = 1 - pchisq(Cov$Testrho, df = 1)),
                 unscale.resid = Cov$unscale.resid,
                 optim         = c(list("objective" = est2$value),
                                   est2[c("counts", "convergence", "message")]))
  } else {
    rho  <- set.rho
    CESd <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                        cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                        lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                        n = n, Kx = Kx, S = S, rho = rho, FE = fixed.effects, 
                        deriv = FALSE, nthread = nthread)
    
    # this is how data are organized in data
    # X: 1 to Kx
    # y: Kx + 1
    # Gy: Kx + 2
    # Gz: Kx + 3
    # dGz: Kx + 4
    
    # Instrument
    Kz  <- 1 + Kiso + Kniso
    Z   <- matrix(0, n, Kz)
    Z[nIs + 1, 1] <- CESd[nIs + 1, Kx + 3]
    if (Kiso > 0) { # if there are isolated
      Z[Is + 1, 2:(1 + Kiso)] <- CESd[Is + 1, idXiso + 1]
    }
    Z[nIs + 1, (2 + Kiso):(Kz)] <- CESd[nIs + 1, idXniso + 1]
    
    # first Optimization
    W     <- solve(crossprod(Z[sel + 1,]) / S) 
    ARG   <- list(f = fCESobjrho, interval = c(-0.5, 100), tol = opttol,
                  y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx],
                  Z = Z, nIso = nIs, W = W, sel = sel, n = n, S = S)
    est1  <- do.call(optimize, ARG)
    
    ## Full parameter
    est1$theta <- fCESparmrho(beta = est1$minimum, rho = rho,
                              y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx],
                              Z = Z, nIso = nIs, W = W, sel = sel, n = n, S = S)
    
    # Second estimation
    est2  <- est1
    if (weight.optim) {
      ARG$W <- fCESWeight_1ins(theta = est1$theta, y = CESd[,Kx + 1], Gy = CESd[,Kx + 2],
                               X = CESd[,1:Kx], Z = Z, Iso = Is, nIso = nIs, 
                               lIso = lIs, lnIso = lnIs, sel = sel, n = n, S = S,
                               Kx = Kx, Kz = Kz, HACn = HACn, dfiso = dfiso, 
                               dfniso = dfniso)
      est2  <- do.call(optimize, ARG)
      
      ## Full parameter
      est2$theta <- fCESparmrho(beta = est2$minimum, rho = rho,
                                y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx],
                                Z = Z, nIso = nIs, W = W, sel = sel, n = n, S = S)
    }
    
    ## Covariance matrix
    Cov <- fCESparmCovrho(theta = est2$theta,
                          y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx], Z = Z,
                          Iso = Is, nIso = nIs, lIso = lIs, lnIso = lnIs, W = W, 
                          sel = sel, n = n, S = S, Kx = Kx, Kz = Kz, HACn = HACn,
                          dfiso = dfiso, dfniso = dfniso)
    
    gmm1 <- list(Estimate = est1$theta,
                 optim    = list("objective" = est1$objective))
    gmm2 <- list(Estimate      = est2$theta,
                 cov           = Cov$cov,
                 sigma         = c(overall = Cov$serr, isolates = Cov$serr_iso, nonisolates = Cov$serr_niso),
                 unscale.resid = Cov$unscale.resid,
                 optim         = list("objective" = est2$objective))
  }

  names(gmm1$Estimate) <- names(gmm2$Estimate) <-
    colnames(gmm2$cov) <- rownames(gmm2$cov) <- c("rho", "beta", xname)
  
  if (!weight.optim) {
    gmm1      <- NULL
  }
  
  out         <- list(model.info = list(n = n_st, ngroup = S, nvec = nvec, formula = formula, 
                                          instrument = instrument, fixed.effects = fixed.effects, 
                                          HAC = HAC, set.rho = set.rho, yname = yname, xnames = xname, 
                                          zname = zename),
                      gmm        = gmm2,
                      first.gmm  = gmm1)
  class(out)  <- "cesconfpeer"
  out
}


#' @title Summary for the Estimation of CES-based Peer Effects Models
#' @param object An object of class \code{\link{cesconfpeer}}.
#' @param ... Further arguments passed to or from other methods.
#' @param x An object of class \code{\link{summary.cesconfpeer}} or \code{\link{cesconfpeer}}.
#' @param fullparameters A logical value indicating whether all parameters should be summarized (may be useful for the structural model).
#' @description Summary and print methods for the class \code{\link{cesconfpeer}}.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the number of subnets, number of observations, and other key details.}
#'     \item{coefficients}{A summary of the estimates, standard errors, and p-values.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#' @export
#' @method summary cesconfpeer
summary.cesconfpeer <- function(object, fullparameters = TRUE, ...) {
  stopifnot(inherits(object, "cesconfpeer"))
  coef           <- fcoef(Estimate = object$gmm$Estimate, cov = object$gmm$cov)
  if (!is.null(object$model.info$set.rho)) {
    coef[1, -c(1, 2)]  <- NA
  }
  out            <- c(object["model.info"], 
                      list(coefficients = coef),
                      object[c("gmm", "first.search")], list(...))
  class(out)     <- "summary.cesconfpeer"
  out
}

#' @rdname summary.cesconfpeer
#' @export
print.summary.cesconfpeer <- function(x, ...) {
  hete <- x$model.info$HAC
  hete <- ifelse(hete %in% c("iid", "group-iid"), hete,
                 ifelse(hete == "hetero", "Individual", "Cluster"))
  sig_overall  <- x$gmm$sigma["overall"]
  sig_iso      <- x$gmm$sigma["isolates"]
  sig_niso     <- x$gmm$sigma["nonisolates"]
  FE           <- x$model.info$fixed.effects
  
  cat("Formula: ", deparse(x$model.info$formula),
      "\nExcluded instrument: ", deparse(x$model.info$instrument), 
      "\nFixed effects: ", paste0(toupper(substr(FE, 1, 1)), tolower(substr(FE, 2, nchar(FE)))), "\n", sep = "")
  
  coef       <- x$coefficients
  coef[,1:2] <- round(coef[,1:2], 7)
  coef[,3]   <- round(coef[,3], 5)
  
  cat("\nCoefficients:\n")
  fprintcoeft(coef)
  
  cat("---\nSignif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n\n")
  cat("HAC: ", hete, sep = "")
  if (x$model.info$HAC == "iid") {
    cat(", sigma : ", format(sig_overall, digits = 5), sep = "")
  } else if (x$model.info$HAC == "group-iid") {
    cat(", sigma (isolated): ", format(sig_iso, digits = 5), ", (non-isolated): ", format(sig_niso, digits = 5), sep = "")
  }
  cat("\nCES parameter -- testing whether rho = 1: prob = ", 
      ifelse(is.null(x$model.info$set.rho), round(x$gmm$testlinear["p-value"], 5), NA), "\n")
  class(x) <- "print.summary.cesconfpeer"
  invisible(x)
}

#' @rdname summary.cesconfpeer
#' @export
print.cesconfpeer <- function(x, ...) {
  print(summary(x))
}



#' @rdname cesconfpeer
#' @export
plotcespeer <- function(formula, 
                        instrument, 
                        Glist, 
                        fixed.effects = FALSE,
                        grid.rho      = seq(-50, 50, 0.1), 
                        tol           = 1e-8, 
                        drop          = NULL,
                        nthread       = 1,
                        data,
                        weight.rho    = NULL,
                        arg.optim     = list(),
                        ...) {
  # Grid for rho
  grid.rho <- grid.rho[grid.rho != 0]
  if (length(grid.rho) < 2) {
    stop("At least two nonzero values for rho is needed.")
  }
  
  ## Thread
  tp        <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread <- tp
  }
  
  ## Network
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist  <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  Glist      <- fGnormalise(Glist, nthread)
  
  ## sizes
  dg       <- fnetwork(Glist = Glist)
  S        <- dg$S
  ldg      <- dg$ldg
  SIs      <- dg$SIs
  SnIs     <- dg$SnIs
  nvec     <- dg$nvec
  n        <- dg$n
  cumsn    <- dg$cumsn
  idpeer   <- dg$idpeer
  Is       <- dg$Is
  lIs      <- dg$lIs
  nIs      <- dg$nIs
  lnIs     <- dg$lnIs
  n_iso    <- length(Is)
  n_niso   <- n - n_iso
  dg       <- dg$dg
  if (!all(dg %in% c(0, 1))) {
    stop("G is not row-normalized.")
  }
  
  ## Formula to data
  formula  <- as.formula(formula)
  f.t.data <- formula2data(formula = formula, data = data, fixed.effects = fixed.effects,
                           simulations = FALSE)
  
  # X, exogenous variables
  X      <- f.t.data$X
  xname  <- f.t.data$xname
  Kx     <- length(xname)
  
  # outcome
  y      <- f.t.data$y
  yname  <- f.t.data$yname
  
  # Instrument
  inst     <- as.formula(instrument); instrument <- inst
  if(length(inst) != 2) stop("The `excluded.instruments` argument must be in the format `~ z` (a single instrument).")
  f.t.data <- formula2data(formula = inst, data = data, fixed.effects = FALSE, 
                           simulations = TRUE)
  z       <- f.t.data$X
  zename  <- f.t.data$xname
  z       <- z[, zename != "(Intercept)"]
  zename  <- zename[zename != "(Intercept)"]
  if (length(zename) != 1) {
    stop("Only one excluded instrument can be used.")
  }
  
  # Create additional variables
  hasIs    <- fCESdatainit(y = y, z = z, G = Glist, nvec = nvec, S = S, 
                           ldg = ldg, lIs = lIs, lnIs = lnIs, drop = drop)
  frindex  <- hasIs$friendindex
  frzeroy  <- hasIs$frzeroy
  frzeroz  <- hasIs$frzeroz
  ldg_st   <- hasIs$ldg
  dg_st    <- hasIs$dg
  S_st     <- hasIs$S
  SIs_st   <- hasIs$SIs
  SnIs_st  <- hasIs$SnIs
  yFMiMa   <- cbind(hasIs$yFmin, hasIs$yFmax)
  zFMiMa   <- cbind(hasIs$zFmin, hasIs$zFmax)
  lIs      <- hasIs$lIs # After selection
  Is       <- hasIs$Is # After selection
  lnIs     <- hasIs$lnIs # After selection
  nIs      <- hasIs$nIs # After selection
  hasIs    <- hasIs$hasIs
  niso     <- length(Is)
  nniso    <- length(nIs)
  n_st     <- niso + nniso
  sel      <- sort(c(Is, nIs))
  if ((1 %in% frzeroy) | (1 %in% frzeroz)) {
    stop("The outcome `y` and the instrument `z` must be strictly positive.")
  }
  
  # idXiso and idXniso
  if (!is.null(weight.rho)) {
    if(weight.rho == 0) {
      stop("`rho` cannot be zero.")
    }
  }
  CESd     <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                       cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                       lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                       n = n, Kx = Kx, S = S, rho = ifnullset(weight.rho, 1), 
                       FE = fixed.effects, deriv = FALSE, nthread = nthread)
  
  idXiso   <- fcheckrank(X = CESd[Is + 1, 1:Kx, drop = FALSE], tol = tol)
  idXniso  <- fcheckrank(X = CESd[nIs + 1, 1:Kx, drop = FALSE], tol = tol)
  if (length(unique(c(idXiso, idXniso))) != Kx) {
    stop("The regressor matrix X is not full rank.")
  }
  if (niso > 0) {
    if (length(idXiso) == 0) {
      stop("The regressor matrix X is not full rank for isolated nodes")
    }
  }
  
  ## degree of freedom
  dfiso    <- NULL
  dfniso   <- NULL
  if (fixed.effects) {
    Kxiso  <- length(idXiso)
    dfiso  <- niso - SIs_st - Kxiso
    dfniso <- nniso - SnIs_st - Kx + Kxiso - 1
  } else {
    Kxiso  <- length(idXiso)
    dfiso  <- niso - Kxiso
    dfniso <- nniso - Kx + Kxiso - 1
  }
  if (dfniso <= 0) {
    stop("Insufficient number of observations for non-isolated nodes.")
  }
  Kiso  <- length(idXiso)
  Kniso <- length(idXniso)
  Kz    <- Kiso + Kniso + 2
  W     <- S * diag(Kz)
  if (!is.null(weight.rho)) {
    Z   <- matrix(0, n, Kz)
    Z[nIs + 1, 1] <- CESd[nIs + 1, Kx + 3]
    Z[nIs + 1, 2] <- CESd[nIs + 1, Kx + 4]
    if (Kiso > 0) { # if there are isolated
      Z[Is + 1, 3:(2 + Kiso)] <- CESd[Is + 1, idXiso + 1]
    }
    Z[nIs + 1, (3 + Kiso):(Kz)] <- CESd[nIs + 1, idXniso + 1]
    
    # first Optimization
    W     <- solve(crossprod(Z[sel + 1,]) / S) 
  }

  # Tolerence argument in optimize
  opttol <- ifnullset(arg.optim$tol, .Machine$double.eps^0.25)
  
  # Estimation
  cl     <- NULL
  on.exit({
    registerDoSEQ()
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  cl     <- makeCluster(nthread)
  registerDoParallel(cl)
  
  nrho   <- length(grid.rho)
  rho    <- NULL  # for R CMD check    
  estim  <- foreach(rho       = grid.rho,
                     .packages = "AsyPeer",
                     .combine  = "cbind"
  ) %dorng% {
    # Data
    CESd <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex,
                     cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs,
                     lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa,
                     n = n, Kx = Kx, S = S, rho = rho, FE = fixed.effects,
                     deriv = FALSE, nthread = 1)

    # this is how data are organized in data
    # X: 1 to Kx
    # y: Kx + 1
    # Gy: Kx + 2
    # Gz: Kx + 3
    # dGz: Kx + 4

    # Instrument
    Z   <- matrix(0, n, 2 + Kiso + Kniso)
    Z[nIs + 1, 1] <- CESd[nIs + 1, Kx + 3]
    Z[nIs + 1, 2] <- CESd[nIs + 1, Kx + 4]
    if (Kiso > 0) { # if there are isolated
      Z[Is + 1, 3:(2 + Kiso)] <- CESd[Is + 1, idXiso + 1]
    }
    Z[nIs + 1, (3 + Kiso):(2 + Kiso + Kniso)] <- CESd[nIs + 1, idXniso + 1]

    # Optimization
    ARG   <- list(f = fCESobjrho, interval = c(-0.5, 100), tol = opttol,
                  y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx],
                  Z = Z, nIso = nIs, W = W, sel = sel, n = n, S = S)
    opt   <- do.call(optimize, ARG)
    
    # Full parameter
    theta <- fCESparmrho(beta = opt$minimum, rho = rho,
                         y = CESd[,Kx + 1], Gy = CESd[,Kx + 2], X = CESd[,1:Kx],
                         Z = Z, nIso = nIs, W = W, sel = sel, n = n, S = S)
    c(opt$objective, theta)
  }
  
  # Plotdata
  plotdata  <- list(estimate = t(estim[-1,]), objective = estim[1,])
  colnames(plotdata$estimate) <- c("rho", "beta", xname)
  rownames(plotdata$estimate) <- names(plotdata$objective) <- NULL
  
  ARG <- list(x    = plotdata$estimate[,"rho"], 
              y    = plotdata$objective, 
              xlab = "rho", ylab = "objective", type = "l",
              ...)
  do.call(plot, ARG)
  
  invisible(plotdata)
}

