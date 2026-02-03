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
#' @param ctl.lbfgs A list of control parameters, such as `maxit`, `eps_f`, and `eps_g`, 
#' to be passed to `optim_lbfgs()` from the \pkg{RcppNumerical} package.
#' 
#' @param grid.rho A finite grid of values for the CES substitution parameter \eqn{\rho} (see Details).
#' This grid is used to plot the objective function after maximizing over the other parameters.
#' It is helpful for determining a reasonable interval for \eqn{\rho}.
#'
#' @param details A Boolean indicating whether optimization details should be printed when plotting the objective function.
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
#' @examples 
#' \donttest{
#' if (requireNamespace("PartialNetwork", quietly = TRUE)) {
#' library(PartialNetwork)
#' ngr  <- 50  # Number of subnets
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
#'             grid.rho = seq(-20, 20, 0.5))
#' est <- cesconfpeer(formula = y ~ X, instrument = ~ z, Glist = Gnorm)
#' summary(est)
#' }
#' }
#' @export
#' @importFrom stats optimise 
cesconfpeer <- function(formula, 
                        instrument, 
                        Glist, 
                        fixed.effects = FALSE,
                        set.rho  = NULL, 
                        interval = c(-100, 100), 
                        tol = 1e-8, 
                        drop = NULL, 
                        HAC = "group-iid", 
                        nthread = 1,
                        data,
                        ctl.lbfgs = list(),
                        ...) {
  ## Thread
  tp        <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread <- tp
  }
  
  # set.rho
  if (!is.null(set.rho)) {
    stopifnot(set.rho != 0)
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
  CESd     <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                       cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                       lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                       n = n, Kx = Kx, S = S, rho = 1, FE = fixed.effects, 
                       deriv = FALSE, nthread = nthread)
  
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
  
  # control
  maxit  <- ifnullset(ctl.lbfgs$maxit, 1e6)
  eps_f  <- ifnullset(ctl.lbfgs$eps_f, 1e-9)
  eps_g  <- ifnullset(ctl.lbfgs$eps_g, 1e-9)
  
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
  
  # Optimization
  est <- suppressWarnings(
    fCESMain(setrho = set.rho, X = X, y = y, z = z, G = Glist, 
                  friendindex = frindex, cumsn = cumsn, frzeroy = frzeroy, 
                  frzeroz = frzeroz, lIso = lIs, lnIso = lnIs, Iso = Is, nIso = nIs,
                  nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, idXiso = idXiso, 
                  idXniso = idXniso, sel = sel, n = n, Kx = Kx, S = S, HACn = HACn, 
                  dfiso = dfiso, dfniso = dfniso, FE = fixed.effects, rhomin = rhomin,
                  rhomax = rhomax, maxit = maxit, eps_f = eps_f, eps_g = eps_g, 
                  nthread = nthread)
  )
  
  names(est$estimate) <- c("rho", "beta", xname)
  colnames(est$cov)   <- rownames(est$cov) <- c("rho", "beta", xname)
  sigma                 <- c(overall = est$serr, isolates = est$serr_iso, 
                         nonisolates = est$serr_niso)
  gmm        <- list(Estimate = est$estimate, 
                     cov      = est$cov, 
                     sigma    = sigma)
  
  if (is.null(set.rho)) {
    gmm   <- c(gmm, list("testlinear" = c("stat"   = est$Testrho,
                                         "p-value" = 1 - pchisq(est$Testrho, df = 1))))
  }
  
  gmm     <- c(gmm, list(unscale.resid = est$unscale.resid,
                         lbfgs = est[c("objective", "gradient", "status")]))
  
  out         <- list(model.info   = list(n = n_st, ngroup = S, nvec = nvec, formula = formula, 
                                          instrument = instrument, fixed.effects = fixed.effects, 
                                          HAC = HAC, set.rho = set.rho, yname = yname, xnames = xname, 
                                          zname = zename),
                      gmm          = gmm)
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
                        grid.rho = seq(-50, 50, 0.1), 
                        tol = 1e-8, 
                        details = TRUE,
                        drop = NULL, 
                        nthread = 1,
                        data,
                        ctl.lbfgs = list(),
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
  CESd     <- fCESdata(X = X, y = y, z = z, G = Glist, friendindex = frindex, 
                       cumsn = cumsn, frzeroy = frzeroy, frzeroz = frzeroz, lIso = lIs, 
                       lnIso = lnIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                       n = n, Kx = Kx, S = S, rho = 1, FE = fixed.effects, 
                       deriv = FALSE, nthread = nthread)
  
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
  
  # control
  maxit  <- ifnullset(ctl.lbfgs$maxit, 1e6)
  eps_f  <- ifnullset(ctl.lbfgs$eps_f, 1e-9)
  eps_g  <- ifnullset(ctl.lbfgs$eps_g, 1e-9)
  
  # plot data
  plotdat <- suppressWarnings(
    fCESplotdata(gridrho = grid.rho, X = X, y = y, z = z, G = Glist, 
                          friendindex = frindex, cumsn = cumsn, frzeroy = frzeroy, 
                          frzeroz = frzeroz, lIso = lIs, lnIso = lnIs, Iso = Is, 
                          nIso = nIs, nvec = nvec, yFMiMa = yFMiMa, zFMiMa = zFMiMa, 
                          idXiso = idXiso, idXniso = idXniso, sel = sel, 
                          FE = fixed.effects, maxit = maxit, eps_f = eps_f, 
                          eps_g = eps_g, nthread = nthread)
  )

  
  ARG <- list(x    = plotdat$estimate[,1], 
              y    = plotdat$objective, 
              xlab = "rho", ylab = "objective", type = "l",
              ...)
  do.call(plot, ARG)
  
  if (details) {
    cat("Summary of the gradient of rho:\n")
    print(summary(plotdat$gradient))
    
    cat("\nFrequency table of the estimation status:")
    print(table(plotdat$status))
  }
  
  invisible(plotdat)
}

