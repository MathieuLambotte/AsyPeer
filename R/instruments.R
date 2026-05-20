#' @title Generating Instruments for the Asymmetric Peer Effects Model
#'
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. 
#'   The formula should be specified as \code{y ~ x1 + x2}, where \code{y} is the outcome and `x1`, `x2`, .... 
#'   are the vectors or matrices of control variables, which may include contextual variables such as 
#'   peer averages.
#' @param ... Further arguments passed to or from other methods.  
#' @param Glist The adjacency matrix or list of adjacency matrices. For networks composed of 
#'   multiple subnets (e.g., schools), \code{Glist} must be a list, where the \code{s}-th element 
#'   is an \eqn{n_s \times n_s} adjacency matrix, and \eqn{n_s} is the number of nodes in the 
#'   \code{s}-th subnet.
#'
#' @param data An optional data frame, list, or environment (or an object that can be coerced to a 
#'   data frame via \link[base]{as.data.frame}) containing the variables in the model. If a variable 
#'   is not found in \code{data}, it is searched for in \code{environment(formula)}, typically the 
#'   environment from which \code{gen.instrument} is called.
#'

#' @param estimator A vector of two character strings specifying the estimators 
#' used for the extensive margin \eqn{P(y_j > y_i)} and the intensive margin 
#' \eqn{E(y_j - y_i \mid y_j - y_i > 0)}.
#'
#' Supported parametric options include `"ols"` and `"lasso"` for linear models, 
#' and `"glm"` for generalized linear models (only for the extensive margin). 
#' Nonparametric options include `"rf"` (random forest) and `"xgboost"` 
#' (gradient boosting).
#'
#' For example, `estimator = c("glm", "xgboost")` uses a generalized linear model 
#' for the extensive margin and gradient boosting for the intensive margin.
#'
#' @param power A numeric vector of length 2 indicating the maximum walk lengths to be used (typically k in \eqn{G^kX}). 
#' The two entries allow specifying a different value of k for each endogenous variable.
#'
#' @param nfold A strictly positive integer specifying the number of folds used when estimating the probability model via cross-fitting.
#' @param drop A dummy vector of the same length as the sample, indicating whether an observation should be dropped.
#' This can be used, for example, to remove false isolates or to estimate the model only on non-isolated agents.
#' These observations cannot be directly removed from the network by the user because they may still be friends with other agents.
#'
#' @param checkrank A logical value indicating whether the linearly dependent columns in the matrix of generated instruments should be drooped.
#' @param nthread Number of CPU cores (threads) used to run parts of the estimation in parallel.
#' @param tol A numeric tolerance used in QR factorization to detect linearly dependent columns in the 
#'   matrices of explanatory variables and instruments, ensuring a full-rank matrix 
#'   (see \link[base]{qr}).
#' @param asymmetry A logical value indicating if the preference for conformity is asymmetric or not.
#' 
#' @description
#' `gen.instrument` generates instruments for the endogenous variables in asymmetric peer effects models.
#' @return A list containing:
#'     \item{model.info}{A list with information about the model, such as the estimator, the number of folds, and other key details.}
#'     \item{instruments}{A matrix of generated instruments.}
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
#' delta   <- 0.25
#' beta    <- c(0.3, 0.6)
#' gamma   <- c(4, 1, -0.7, 0, -0.5) 
#' eps     <- rnorm(n, 0, 0.5) 
#' 
#' ## Generating `y`
#' y <- asypeer.sim(formula = ~ X + GX, Glist = Gnorm, delta = delta, 
#'                  beta = beta, gamma = gamma, epsilon = eps)
#' y <- y$y
#' 
#' ### Generating instruments
#' ins <- gen.instrument(formula = y ~ X, Glist = Gnorm, 
#'                       estimator = c("logit", "ols"))}
#'  
#' @export
gen.instrument <- function(formula,
                           Glist, 
                           data,
                           asymmetry = TRUE,
                           estimator = c("logit", "ols"),
                           power     = c(1, 1),
                           nfold     = 2,
                           checkrank = TRUE,
                           tol       = 1e-10,
                           nthread   = 1,
                           drop = NULL,
                           ...) {
  ## power for G
  power   <- as.integer(power)
  if (length(power) == 1) {
    power <- rep(power, 2)
  } else if (length(power) != 2) {
    stop("`power` is expected to be a vector containing 2 integers: The maximal power of G to compute instruments for ybar and ycheck, respectively.")
  }
  if (any(power < 1)) {
    stop("`power` cannot be negative or zero.")
  }
  
  ## estimator
  if (length(estimator) == 1) {
    estimator    <- rep(estimator, 2)
  }
  
  if (tolower(estimator[1]) %in% c("lin", "linear", "ols", "lm")) {
    estimator[1] <- "OLS"
  } else if (tolower(estimator[1]) %in% c("logit", "logistic")) {
    estimator[1] <- "Logit"
  } else if (tolower(estimator[1]) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
    estimator[1] <- "Random Forest"
  } else if (tolower(estimator[1]) %in% c("xgboost")) {
    estimator[1] <- "XGBoost"
  } else if (tolower(estimator[1]) %in% c("lasso")) {
    estimator[1] <- "LASSO"
  } else {
    stop("This estimator is not available.")
  }
  
  
  if (tolower(estimator[2]) %in% c("lin", "linear", "ols", "lm")) {
    estimator[2] <- "OLS"
  }  else if (tolower(estimator[2]) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
    estimator[2] <- "Random Forest"
  } else if (tolower(estimator[2]) %in% c("xgboost")) {
    estimator[2] <- "XGBoost"
  } else if (tolower(estimator[2]) %in% c("lasso")) {
    estimator[2] <- "LASSO"
  } else {
    stop("This estimator is not available.")
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
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  Glist    <- fGnormalise(Glist, nthread)
  
  ## sizes
  dg       <- fnetwork(Glist = Glist)
  S        <- dg$S
  if (S < 2) {
    stop("At least two subnets are required.")
  }
  nvec     <- dg$nvec
  n        <- dg$n
  cumsn    <- dg$cumsn
  idpeer   <- dg$idpeer
  dg       <- dg$dg
  if (!all(dg %in% c(0, 1))) {
    stop("G is not row-normalized.")
  }
  
  ## Formula to data
  formula  <- as.formula(formula)
  if (missing(data)) {
    data  <- env(formula)
  }
  f.t.data <- formula2data(formula = formula, data = data, fixed.effects = TRUE,
                           simulations = FALSE)
  y        <- f.t.data$y
  yname    <- f.t.data$yname
  
  X        <- f.t.data$X
  xname    <- f.t.data$xname
  X        <- cbind(1 - dg, X * (1 - dg), dg, X * dg)
  xname    <- c("iso", paste0("iso_", xname), "niso", paste0("niso_", xname))
  
  ### Drop
  if (is.null(drop)){
    drop  <- rep(0, n)
  }
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  keep    <- !as.logical(drop)
  
  ### Instrument for ybar
  insBary <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                          power = power[1], nthread = nthread)
  
  ### Instrument for ycheck
  out     <- NULL
  if(asymmetry){
    Xtp      <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                             power = power[2], nthread = nthread)
    ## Dyadic dada
    group    <- rep(0:(S - 1), nvec)
    IDi      <- unlist(lapply(1:S, \(s) 0:(nvec[s] - 1)))
    gij      <- lapply(1:n, \(i) Glist[[group[i] + 1]][IDi[i] + 1, idpeer[[i]] + 1])
    ddni     <- sapply(gij, length)
    ddncs    <- c(0, cumsum(ddni))
    ddkeep   <- rep(keep, ddni)
    
    tp       <- fdataML(y = y, X = Xtp, group = group, IDi = IDi, gij = gij,
                        idpeer = idpeer, ddni = ddni, ddncs = ddncs, ncs = cumsn,
                        nthread = nthread)
    
    ddyext   <- as.integer(tp$ddy[ddkeep, 7]) # extensive margin
    ddyint   <- tp$ddy[ddkeep, 6] - tp$ddy[ddkeep, 5] #intensive margin
    ddX      <- tp$ddXj - tp$ddXi
    ddX      <- ddX[ddkeep, ,drop = FALSE]
    
    ## Check rank of ddX
    ddX      <- as.data.frame(ddX[, fcheckrank(X = ddX, tol = tol) + 1, drop = FALSE])
    
    ## Fold construction
    nfold    <- as.integer(nfold)
    if (nfold > S) {
      nfold  <- S
      warning("The number of folds exceeds the number of subnets; it has been reset to the number of subnets.")
    } 
    if (nfold == 1) {
      stop("At lead two folds is required.")
    }
    id_fold  <- fassignfold(tp$ddy[ddkeep, 1], nfold = nfold)
    
    ## Prediction
    ARG      <- list(ddyext = ddyext, ddyint = ddyint, ddX = ddX, id_fold = id_fold, 
                     estimatorext = estimator[1], estimatorint = estimator[2], 
                     nthread = nthread, ...)
    
    insChey  <- fInstChecky(rhoddX = as.matrix(do.call(mpredict, ARG)), ddni = ddni, 
                            nthread = nthread)
    
    GinsChey <- peeravgpower(G = Glist, V = insChey, cumsn = cumsn, nvec = nvec, 
                             power = power[1], nthread = nthread)
    
    out      <- cbind(insBary, insChey, GinsChey)
    colnames(out) <- c(sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, xname)),
                       "y_check_hat", sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, "y_check_hat")))
  } else {
    out      <- insBary
    colnames(out) <- c(sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, xname)))
    
  }
  if (checkrank) {
    keepcol    <- fcheckrank(X = out[keep, , drop = FALSE], tol = tol) + 1
    out        <- out[, keepcol, drop = FALSE]
  }
  out[!keep, ] <- NA
  list(model.info  = list(power = power, estimator = estimator,
                          nfold = nfold, tol = tol),
       instruments = out)
}

#' @rdname gen.instrument
#' @export
gen.instruments <- function(formula,
                            Glist, 
                            data,
                            asymmetry = TRUE,
                            estimator = c("logit", "ols"),
                            power     = c(1, 1),
                            nfold     = 2,
                            checkrank = TRUE,
                            tol       = 1e-10,
                            nthread   = 1,
                            drop      = NULL,
                            ...) { 
  if (missing(data)) {
    data  <- env(formula)
  }
  ARG <- list(formula = formula, Glist = Glist, data = data, asymmetry = asymmetry, 
              estimator = estimator, power = power, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread,
              drop = drop, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.insts <- function(formula,
                      Glist, 
                      data,
                      asymmetry = TRUE,
                      estimator = c("logit", "ols"),
                      power     = c(1, 1),
                      nfold     = 2,
                      checkrank = TRUE,
                      tol       = 1e-10,
                      nthread   = 1,
                      drop      = NULL,
                      ...) { 
  if (missing(data)) {
    data  <- env(formula)
  }
  ARG <- list(formula = formula, Glist = Glist, data = data, asymmetry = asymmetry, 
              estimator = estimator, power = power, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread,
              drop = drop, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.inst <- function(formula,
                     Glist, 
                     data,
                     asymmetry = TRUE,
                     estimator = c("logit", "ols"),
                     power     = c(1, 1),
                     nfold     = 2,
                     checkrank = TRUE,
                     tol       = 1e-10,
                     nthread   = 1,
                     drop      = NULL,
                     ...) { 
  if (missing(data)) {
    data  <- env(formula)
  }
  ARG <- list(formula = formula, Glist = Glist, data = data, asymmetry = asymmetry, 
              estimator = estimator, power = power, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread,
              drop = drop, ...)
  do.call(gen.instrument, ARG)
}