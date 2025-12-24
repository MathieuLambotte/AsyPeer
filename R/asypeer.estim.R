#' @title Estimating the Asymmetric Peer Effects Model
#'
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. 
#'   `formula` should be specified as \code{y ~ x1 + x2 + ....}, where `y` is the outcome and `x1`, `x2`, .... 
#'   are the vectors or matrices of control variables, which may include contextual variables such as averages among peers.
#'
#' @param excluded.instruments A \link[stats]{formula} specifying the excluded
#'   instruments. It should be written as \code{~ z1 + z2 + ...}, where \code{z1},
#'   \code{z2}, ... are the vectors or matrices of excluded instruments for the two
#'   endogenous variables in the asymmetric model: the average peers' outcomes and
#'   the average difference in outcomes between an agent and her peers who have a
#'   higher outcome. If omitted, default instruments are generated using
#'   \link{gen.instrument}.
#'
#' @param common.gamma A character vector indicating the control variables assumed to have the same marginal effect on marginal utility for isolated and non-isolated agents. 
#' If omitted, all control variables are included. A value of \code{NULL} or \code{character(0)} indicates that none of the variables.
#'
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), 
#'   `Glist` must be a list of subnets, with the \code{s}-th element being an \eqn{n_s \times n_s} 
#'   adjacency matrix, where \eqn{n_s} is the number of nodes in the \code{s}-th subnet.
#'
#' @param data An optional data frame, list, or environment (or an object that can be coerced to a 
#'   data frame via \link[base]{as.data.frame}) containing the variables in the model. If a variable 
#'   is not found in `data`, it is taken from \code{environment(formula)}, typically the environment 
#'   from which `asypeer.estim` is called.
#'
#' @param weight A character string specifying the weighting matrix used in the GMM estimation. 
#'   Available options are: `"identity"` for GMM with the identity matrix as the weighting matrix; 
#'   `"IV"` for the standard instrumental-variable GMM estimator; and `"optimal"` for GMM with the 
#'   optimal weighting matrix.
#'
#' @param tol A tolerance value used in the QR factorization to detect collinear columns in the 
#'   matrices of explanatory variables and instruments, ensuring a full-rank matrix (see 
#'   \link[base]{qr}).
#'
#' @param HAC A character string specifying the correlation structure of the idiosyncratic errors 
#'   for covariance computation. Options are `"iid"` for independent errors; `"group iid"` for 
#'   independence within the groups of isolated and non-isolated players; `"hetero"` for 
#'   heteroskedastic but non-autocorrelated errors; and `"cluster"` for heteroskedastic errors with 
#'   potential within-subnetwork correlation.
#'
#' @param fixed.effects A logical value indicating whether the model includes subnetwork fixed effects.
#'
#' @param nthread Number of CPU cores (threads) used to run parts of the estimation in parallel.
#' 
#' @param drop A dummy vector of the same length as the sample, indicating whether an observation should be dropped.
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user because they may still be friends with other agents.
#'   
#' @param spillover A logical value indicating if the model included a spillover component in additional to the asymmetric conformity effect.
#'    
#' @param asymmetry A logical value indicating if the preference for conformity is asymmetric or not.
#'    
#' @param ... Further arguments passed to or from other methods.   
#'
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the number of subnets, number of observations, and other key details.}
#'     \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix, and related statistics.}
#'     \item{data}{A list including the original data used to estimate the model, as well as the endogenous variables and their instruments.}
#' 
#' @description
#' `asypeer.estim` estimates the asymmetric peer effects model. The instruments are generated automatically following the procedure described 
#'   in the paper and implemented in \link{gen.instruments}. The user can also supply her own instruments.
#' @examples
#' if (requireNamespace("PartialNetwork", quietly = TRUE)) {
#' library(PartialNetwork)
#' ngr  <- 50  # Number of subnets
#' nvec <- rep(30, ngr)  # Size of subnets
#' n    <- sum(nvec)
#' 
#'### Simulating Data
#' ## Network matrix
#' G   <- lapply(1:ngr, function(z) {
#'  Gz <- matrix(rbinom(nvec[z]^2, 1, 0.3), nvec[z], nvec[z])
#'  diag(Gz) <- 0
#'  # Adding isolated nodes (important for the structural model)
#'  niso <- sample(0:nvec[z], 1, prob = ((nvec[z] + 1):1)^5 / sum(((nvec[z] + 1):1)^5))
#'  if (niso > 0) {
#'    Gz[sample(1:nvec[z], niso), ] <- 0
#'  }
#'  Gz
#' })
#' 
#' Gnorm   <- norm.network(G)
#' X       <- cbind(rnorm(n, 0, 2), rpois(n, 2))
#' GX      <- peer.avg(Gnorm, X)
#' delta   <- 0.25
#' beta    <- c(0.3, 0.6)
#' gamma   <- c(4, 1, -0.7, 0, -0.5) 
#' eps     <- rnorm(n, 0, 0.5) 
#' 
#' ## Generating `y`
#' y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta = delta, beta = beta, 
#'                 gamma = gamma, epsilon = eps)
#' y <- y$y
#' 
#' ### Estimating the asymmetric peer effects model
#' est <- asypeer.estim(formula=y ~ X + GX, Glist = Gnorm)
#' summary(est, diagnostic = TRUE)
#' }                  
#' @export
#' @importFrom stats pchisq
#' @importFrom stats optimize
asypeer.estim <- function(formula,
                          excluded.instruments,
                          Glist,
                          data,
                          common.gamma,
                          spillover = TRUE,
                          asymmetry = TRUE,
                          weight = "IV",
                          HAC = "group-iid",
                          fixed.effects = FALSE,
                          tol = 1e-10,
                          nthread = 1,
                          drop = NULL,
                          ...){
  ## Thread
  tp        <- fnthreads(nthread = nthread)
  if ((tp == 1) & (nthread != 1)) {
    warning("OpenMP is not available. Sequential processing is used.")
    nthread <- tp
  }
  
  ## Instrument
  Z        <- NULL
  detInst  <- list()
  zname    <- NULL
  if (missing(excluded.instruments)) {
    if (missing(data)) {
      data <- env(formula)
    }
    excluded.instruments <- NULL
    ARG    <-  list(formula = formula, Glist = Glist, data = data, 
                    nthread = nthread, drop = drop, tol = tol, ...) ## Remember to add ... when the function is completed
    Z      <- do.call(gen.inst, ARG)
    detInst$estimator <- Z$model.info$estimator
    detInst$power     <- Z$model.info$power
    Z                 <- Z$instruments
    zname             <- colnames(Z)
  }  else {
    if(length(excluded.instruments) != 2) stop("The `excluded.instruments` argument must be in the format ~ z1 + z2 + ...")
    if (missing(data)) {
      data      <- env(formula)
    }
    f.t.data    <- formula2data(formula = excluded.instruments, 
                                data = data, fixed.effects = FALSE, 
                                simulations = TRUE)
    Z           <- f.t.data$X
    zname       <- f.t.data$xname
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
  ## Weight
  if (tolower(weight) %in% c("i","identity")) {
    weight   <- "I"
  } else if (tolower(weight) %in% c("iv")) {
    weight   <- "IV"
  } else if (tolower(weight) %in% c("opti", "optimal","optim")) {
    weight   <- "optimal"
  } else {
    stop("This weighting option is not available.")
  }
  
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
  SIso     <- dg$SIs
  SnIso    <- dg$SnIs
  nvec     <- dg$nvec
  n        <- dg$n
  cumsn    <- dg$cumsn
  idpeer   <- dg$idpeer
  Iso      <- dg$Is
  lIso     <- dg$lIs
  nIso     <- dg$nIs
  lnIso    <- dg$lnIs
  n_iso    <- length(Iso)
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
  
  # endogenous variables
  endo   <- highlowstat1(X = as.matrix(y), G = Glist, cumsn = cumsn, nvec = nvec, 
                         ngroup = S, nthread = nthread)
  if(asymmetry){
  endo   <- cbind(yb = endo$Xbar, ydot = endo$Xbh - endo$gh * y)
  colnames(endo) <- paste0(yname, c("_bar", "_dot"))
  } else {
    endo   <- endo$Xbar
    colnames(endo) <- paste0(yname, "_bar")
  }
  #drop
  if (!is.null(drop)) {
    dg       <- fdrop(drop = drop, ldg = ldg, S = S, nvec = nvec, y = y, 
                      X = X, Z = Z, endo = endo)
    S        <- dg$S
    SIso     <- dg$SIs
    SnIso    <- dg$SnIs
    nvec     <- dg$nvec
    n        <- dg$n
    cumsn    <- c(0, cumsum(nvec))
    Iso      <- dg$Is
    lIso     <- dg$lIs
    nIso     <- dg$nIs
    lnIso    <- dg$lnIs
    n_iso    <- length(Iso)
    n_niso   <- n - n_iso
    y        <- dg$y
    X        <- dg$X
    Z        <- dg$Z
    endo     <- dg$endo
    dg       <- dg$dg
  }
  invisible(gc())
  if (spillover && SnIso == 0) {
    stop("At least one subnet containing isolated nodes is required to estimate spillover effects.")
  }
  if (SnIso < 2) {
    stop("At least two subnets containing non-isolated nodes are required.")
  }
  
  # Demean
  if (fixed.effects) {
    y      <- c(Demean_separate(X = as.matrix(y), cumsn = cumsn, lIso = lIso, 
                                lnIso = lnIso, nthread = nthread))
    endo   <- Demean_separate(X = endo, cumsn = cumsn, lIso = lIso, 
                              lnIso = lnIso, nthread = nthread)
    X      <- Demean_separate(X = X, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                              nthread = nthread)
    Z      <- Demean_separate(X = Z, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                              nthread = nthread)
  }
  
  # X variables 
  colnames(X) <- xname
  X_iso       <- X * (1 - dg)
  X_niso      <- X * dg
  
  # common.gamma 
  if(missing(common.gamma)) {
    common.gamma <- xname
  }
  if(!all(common.gamma %in% xname)){
    stop("At least one variable name in `common.gamma` is not written correctly.")
  }
  xnameiso  <- xname[fcheckrank(X = X_iso[Iso + 1,], tol = tol) + 1]
  xnameniso <- xname[fcheckrank(X = X_niso[nIso + 1,], tol = tol) + 1]
  cgamma    <- xnameiso[(xnameiso %in% common.gamma) & (xnameiso %in% xnameniso)]
  ncgamma   <- setdiff(xnameiso, cgamma)
  cgamma    <- which(xname %in% cgamma) - 1
  ncgamma   <- which(xname %in% ncgamma) - 1
  if ((length(cgamma) == 0) && spillover) {
    stop("`common.gamma` is needed to estimate spillover effects.")
  }
  
  ## [endo, X_iso, X_niso] should be full rank
  if (length(fcheckrank(X = cbind(endo, X_iso, X_niso), tol = tol)) < (Kx + 1L + asymmetry + spillover)) {
    stop("The design matrix is not full rank.")
  }
  
  # Instruments
  Z           <- cbind(X_iso, X_niso, Z)
  zname       <- c(paste0("iso_", xname), paste0("niso_", xname), zname)
  colnames(Z) <- zname
  ## Remove unecessary instruments
  keepZ <- fcheckrank(X = Z, tol = tol) + 1
  Z     <- Z[, keepZ, drop = FALSE]
  zname <- zname[keepZ]  
  
  ## Check the number of instruments
  Kz     <- ncol(Z)
  if (Kz < (Kx + 1L + asymmetry + spillover)) {
    stop("Insufficient number of instruments: the model is not identified.")
  }
  if ((HAC == "cluster") & (Kz >= S)) {
    stop("Errors cannot be clustered because the number of instruments is larger than the number of subnets.")
  }
  
  ## degree of freedom
  dfiso    <- NULL
  dfniso   <- NULL
  if (fixed.effects) {
    Kxiso  <- length(xnameiso)
    dfiso  <- n_iso - SIso - Kxiso
    dfniso <- n_niso - SnIso - asymmetry -1  - spillover - Kx + Kxiso
  } else {
    Kxiso  <- length(xnameiso)
    dfiso  <- n_iso - Kxiso
    dfniso <- n_niso - 1 - asymmetry - spillover - Kx + Kxiso
  }
  if ((dfiso <= 0) && HACn == 1) {
    HAC    <- "iid"
    HACn   <- 0
  }
  if (spillover && (dfiso <= 0)) {
    stop("Insufficient number of observations for isolated nodes.")
  }
  if (dfniso <= 0) {
    stop("Insufficient number of observations for non-isolated nodes.")
  }
  
  ## Weight
  if (weight == "I"){
    W    <- diag(nrow = ncol(Z))
  } else if(weight %in% c("IV", "optimal")){ 
    W    <- solve(crossprod(Z) / n)
  }
  
  ## Optimization
  fGMM     <- gmm_obj_nospil
  fWopt    <- W_optimal_nospil
  fest     <- compute_estimate_nospil
  if (spillover){
    fGMM   <- gmm_obj
    fWopt  <- W_optimal
    fest   <- compute_estimate
  }
  
  # solve with the first step weighting matrix
  gmm    <- optimize(f = fGMM, Z = Z, y = y, endo = endo, X_iso = X_iso,
                     c_gamma = cgamma, nc_gamma = ncgamma,
                     X_niso = X_niso, W = W, S = S, lower = -0.999, upper =  20)
  
  #get the estimate of beta_l
  betal  <- gmm$minimum
  
  # get the optimal weighting matrice given the estimated beta_l
  if(weight == "optimal"){
    W    <- fWopt(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                  c_gamma = cgamma, nc_gamma = ncgamma,
                  X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, cumsn = cumsn,
                  dfiso = dfiso, dfniso = dfniso, HAC = HACn, S = S)
    
    #gmm with optimal W
    gmm  <- optimize(f = fGMM, Z = Z, y = y, endo = endo, X_iso = X_iso, 
                     c_gamma = cgamma, nc_gamma = ncgamma,
                     X_niso = X_niso, W = W, S = S, lower = -0.999, upper =  20)
    
    #get the optimal estimate of beta_l
    betal <- gmm$minimum
  }
  
  # and the associated estimates of phi
  estim   <- fest(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                  c_gamma = cgamma, nc_gamma = ncgamma,
                  X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, 
                  cumsn = cumsn, dfiso = dfiso, dfniso = dfniso, HAC = HACn, S = S)
  
  ## Shape the output
  gmm     <- list(redparms = estim$redparm,
                  strparms = estim$strparm,
                  redcov   = estim$redcov,
                  strcov   = estim$strcov,
                  s2       = c(overall = estim$s2, isolates = estim$s2iso, nonisolates = estim$s2niso),
                  Sargan   = list(stat   = ifelse(estim$Jdf > 0, estim$JStat, NA),
                                  df     = estim$Jdf,
                                  pvalue = ifelse(estim$Jdf > 0, 1 - pchisq(estim$JStat, estim$Jdf), NA)))
  
  gname   <- paste0("gamma:", xname)
  if (length(ncgamma) > 0) {
    gname <- c(gname, paste0("isogamma:", xname[ncgamma + 1]))
  }
  if(spillover){
    if(asymmetry){
    names(gmm$redparms) <- colnames(gmm$redcov) <- rownames(gmm$redcov) <- 
      c("betal", "theta1", "theta2", gname) 
    names(gmm$strparms) <- colnames(gmm$strcov) <- rownames(gmm$strcov) <- 
      c("betal", "betah", "delta", gname)
    } else {
      names(gmm$redparms) <- colnames(gmm$redcov) <- rownames(gmm$redcov) <- 
        c("beta", "theta1", gname) 
      names(gmm$strparms) <- colnames(gmm$strcov) <- rownames(gmm$strcov) <- 
        c("beta", "delta", gname)
    }
  } else {
    if(asymmetry){
    names(gmm$redparms) <- colnames(gmm$redcov) <- rownames(gmm$redcov) <- 
      c("theta1","theta2", gname) 
    names(gmm$strparms) <- colnames(gmm$strcov) <- rownames(gmm$strcov) <- 
      c("betal", "betah", gname)
    } else {
      names(gmm$redparms) <- colnames(gmm$redcov) <- rownames(gmm$redcov) <- 
        c("theta1", gname) 
      names(gmm$strparms) <- colnames(gmm$strcov) <- rownames(gmm$strcov) <- 
        c("beta", gname)
    }
  }
  
  out     <- list(model.info  = c(list(n = n, n_iso=n_iso, ngroup = S, nvec = nvec, formula = formula, 
                                       excluded.instruments = excluded.instruments, spillover = spillover,
                                       weight = weight, HAC = HAC, fixed.effects = fixed.effects, asymmetry = asymmetry,
                                       tol = tol, xname = xname, yname = yname, zname = zname, 
                                       dfiso = dfiso, dfniso = dfniso, common.gamma = xname[cgamma + 1],
                                       ncommon.gamma = xname[ncgamma + 1]), detInst),
                  gmm         = gmm,
                  data        = list(dependent = y, exogenous = X,
                                     endogenous.variables = endo, instruments = Z, 
                                     degree = dg, isolates = lapply(lIso, \(x) x + 1),
                                     non.isolates = lapply(lnIso, \(x) x + 1)))
  class(out) <- "asypeer.estim"
  out
}

#' @title Summary and Print Methods for the Asymmetric Peer Effects Model
#' @param object An object of class \code{\link{asypeer.estim}} as returned by the function \link{asypeer.estim}.
#' @param structural A logical value indicating whether the summary should display the structural parameters (\code{structural = TRUE}) or the reduced-form parameters (\code{structural = FALSE}).
#' @param diagnostics,diagnostic A logical value indicating whether diagnostic tests for the IV GMM should be performed.
#'   These include an F-test of the first-stage regression for weak instruments, a Wu-Hausman test
#'   for endogeneity, and a Hansen's J-test for overidentifying restrictions (the latter only when
#'   the number of instruments exceeds the number of regressors).
#' @param KPtest A logical value indicating whether a Kleibergen-Paap Wald test (5% level) should be performed in addition to the standard F test 
#'   of the first-stage regression for weak instruments.
#' @param nthread A strictly positive integer specifying the number of threads used in
#'   computationally intensive steps of the estimation procedure.
#' @param x An object of class \code{\link{summary.asypeer.estim}} or \code{\link{asypeer.estim}} as returned by the function \link{summary.asypeer.estim} or \link{asypeer.estim}, respectively.
#' @param ... Further arguments passed to or from other methods.
#'
#' @description
#' Summary and print methods for objects of class \code{\link{asypeer.estim}}.
#'
#' @return A list containing:
#'   \item{model.info}{A list containing information about the model, such as the number of subnets,
#'     number of observations, and other key details.}
#'   \item{coefficients}{A summary of coefficient estimates, standard errors, and p-values.}
#'   \item{diagnostics}{A summary of the diagnostic tests for the instrumental-variable regression,
#'     if requested.}
#'   \item{gmm}{A list of GMM estimation results, including parameter estimates, the covariance matrix,
#'     and related statistics.}
#'
#' @export

summary.asypeer.estim <- function(object, 
                                  structural  = TRUE, 
                                  diagnostic  = FALSE, 
                                  diagnostics = FALSE,
                                  KPtest      = diagnostics || diagnostic,  
                                  nthread     = 1L,
                                  ...) {
  stopifnot(inherits(object, "asypeer.estim"))
  diagn   <- NULL
  if (diagnostic || diagnostics) {
    tp        <- fnthreads(nthread = nthread)
    if ((tp == 1) & (nthread != 1)) {
      warning("OpenMP is not available. Sequential processing is used.")
      nthread <- tp
    }
    diagn   <- fdiagnostic(object, KPtest, nthread)
  }
  yname   <- object$model.info$yname
  xnames  <- object$model.info$xnames
  if (structural) {
    est   <- object$gmm$strparms
    covt  <- object$gmm$strcov
  } else{
    est   <- object$gmm$redparms
    covt  <- object$gmm$redcov
  }
  
  coef    <- fcoef(Estimate = est, cov = covt)
  out     <- c(object["model.info"], 
               list(coefficients = coef, 
                    diagnostics  = diagn),
               object["gmm"], list(...))
  class(out) <- "summary.asypeer.estim"
  out
}

#' @rdname summary.asypeer.estim
#' @export
print.summary.asypeer.estim <- function(x, ...) {
  esti <- ifelse(x$model.info$weight == "identity", "GMM (Weight: Identity Matrix)", 
                 ifelse(x$model.info$weight == "optimal", "GMM (Weight: Optimal)",
                        ifelse(x$model.info$weight == "IV", "GMM (Weight: IV)" )))
  hete <- x$model.info$HAC
  hete <- ifelse(hete %in% c("iid", "group-iid"), hete,
                 ifelse(hete == "hetero", "Individual", "Cluster"))
  sig_overall  <- x$gmm$s2["overall"]
  sig_iso      <- x$gmm$s2["isolates"]
  sig_niso     <- x$gmm$s2["nonisolates"]
  FE   <- x$model.info$fixed.effects
  cat("Formula: ", deparse(x$model.info$formula),
      "\nExcluded instruments: ", ifelse(!is.null(x$model.info$excluded.instruments), 
                                         deparse(x$model.info$excluded.instruments),
                                         paste("(G^p)X with max(p) =", max(x$model.info$power),"and", 
                                               ifelse(x$model.info$estimator == "rf", "Random Forest", 
                                                      toupper(x$model.info$estimator)), "predictions")),
      "\n\n# Subnetworks:", x$model.info$ngroup,
      "\n# Isolates:", x$model.info$n_iso,
      "\n# Non-isolates:", x$model.info$n - x$model.info$n_iso,
      "\n\nEstimator: ", esti,
      "\nFixed effects: ", ifelse(FE, "Yes", "No"), "\n", sep = "")
  
  coef       <- x$coefficients
  coef[,1:2] <- round(coef[,1:2], 7)
  coef[,3]   <- round(coef[,3], 5)
  cat("\nCoefficients:\n")
  fprintcoeft(coef)
  
  if (!is.null(x$diagnostics)) {
    coef       <- x$diagnostics
    coef[,3]   <- round(coef[,3], 5)
    cat("\nDiagnostic tests:\n")
    fprintcoeft(coef) 
  }
  cat("---\nSignif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01 \u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n")
  
  cat("\nHAC: ", hete, sep = "")
  if(hete %in% c("group-iid")){
    cat(", sigma (isolates): ", format(sig_iso, digits = 5), ", (non-isolates): ", format(sig_niso, digits = 5), sep = "")
  } else if (hete == "iid") {
    cat(", sigma: ", format(sig_overall, digits = 5), sep = "")
  }
  
  ## range
  delta   <- 0
  if (x$model.info$spillover) {
    delta <- x$gmm$strparms["delta"]
  }
  minbeta <- ifelse(x$model.info$asymmetry, min(x$gmm$strparms[c("betal", "betah")]), x$gmm$strparms["beta"])
  maxbeta <- ifelse(x$model.info$asymmetry, max(x$gmm$strparms[c("betal", "betah")]), x$gmm$strparms["beta"])
  if (minbeta > -1) {
    if(x$model.info$asymmetry){ 
      boundl <- unname((delta + minbeta) / (1 + maxbeta))
      boundh <- unname((delta + maxbeta) / (1 + minbeta))
      if ((abs(boundl) < 1) && (abs(boundh) < 1)) {
        cat("\nTotal Peer Effects Range:", " [", deparse(round(boundl,4)),", ", deparse(round(boundh,4)),"]\n", sep = "")
      } else{
        warning("Total Peer effects are outside the [-1, 1] interval, there might be multiple equilibria.")
      }
    } else {
      bound  <- unname((delta + minbeta) / (1 + minbeta))
      if (abs(bound) < 1) {
        cat("\nTotal Peer Effects: ",deparse(round(bound,4)), "\n", sep = "")
      } else{
        warning("Total Peer effects are outside the [-1, 1] interval, there might be multiple equilibria.")
      }
    }
  } else{
    warning("Conformity parameter lower than -1")
  }
  
  class(x) <- "print.summary.asypeer.estim"
  invisible(x)
}

#' @rdname summary.asypeer.estim
#' @export
print.asypeer.estim <- function(x, ...) {
  print(summary(x))
}
