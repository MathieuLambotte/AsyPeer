#' @title Generate Instruments for the Asymmetric Peer Effects Model
#'
#' @param formula An object of class \link[stats]{formula}: a symbolic description 
#'   of the model. The formula should be specified as \code{y ~ x1 + x2}, where 
#'   \code{y} is the outcome variable and \code{x1}, \code{x2}, ... are control 
#'   variables, which may include contextual variables such as peer averages.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @param Glist The adjacency matrix or a list of adjacency matrices. For networks 
#'   composed of multiple subnets (e.g., schools), \code{Glist} must be a list, 
#'   where the \code{s}-th element is an \eqn{n_s \times n_s} adjacency matrix, 
#'   and \eqn{n_s} is the number of nodes in subnet \code{s}.
#'
#' @param data An optional data frame, list, or environment (or an object that can 
#'   be coerced to a data frame via \link[base]{as.data.frame}) containing the 
#'   variables in the model. If a variable is not found in \code{data}, it is 
#'   searched for in \code{environment(formula)}, typically the environment from 
#'   which \code{gen.instrument} is called.
#'
#' @param drop A logical vector of the same length as the sample, indicating which 
#'   observations should be dropped. This can be used, for example, to remove false 
#'   isolates or to estimate the model only on non-isolated agents. These 
#'   observations cannot be physically removed from the network structure because 
#'   they may still be connected to other agents.
#'
#' @param estimator A character vector of length 2 specifying the estimators used 
#'   for the intensive margin \eqn{E(y_j - y_i \mid y_j - y_i > 0)} (or \eqn{E(\bar{y}_i})) 
#'   and the extensive margin \eqn{P(y_j > y_i)}.
#'
#'   Supported parametric options include \code{"ols"} and \code{"lasso"} for linear 
#'   models, and \code{"glm"} for generalized linear models (only for the extensive 
#'   margin). Nonparametric options include \code{"rf"} (random forest) and 
#'   \code{"xgboost"} (gradient boosting).
#'
#'   For example, \code{estimator = c("xgboost", "glm")} uses gradient boosting for 
#'   the intensive margin and a generalized linear model for the extensive margin.
#'
#' @param power A numeric vector of length 2 indicating the maximum walk lengths 
#'   (typically \eqn{k} in \eqn{G^k X}) used to construct or estimate instruments for 
#'   \eqn{\bar{y}_i} and \eqn{\check{y}_i}. The two entries allow different values 
#'   of \eqn{k} for each endogenous variable.
#'
#' @param nfold A strictly positive integer specifying the number of folds used for 
#'   cross-fitting in the estimation of the probability model.
#'
#' @param checkrank A logical value indicating whether linearly dependent columns 
#'   in the matrix of generated instruments should be dropped.
#'
#' @param nthread Number of CPU cores (threads) used to run parts of the estimation 
#'   in parallel.
#'
#' @param tol A numeric tolerance used in QR factorization to detect linearly 
#'   dependent columns in the matrices of explanatory variables and instruments, 
#'   ensuring a full-rank matrix (see \link[base]{qr}).
#'
#' @param asymmetry A logical value indicating whether preferences for conformity 
#'   are asymmetric.
#'
#' @param full A logical value. If \code{TRUE}, predictions for both 
#'   \eqn{\bar{y}_i} and \eqn{\check{y}_i} are used as instruments. If \code{FALSE}, 
#'   only the prediction for \eqn{\check{y}_i} is used as an instrument, while 
#'   \eqn{\bar{y}_i} is instrumented using the usual instruments, \eqn{G^k X}.
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
#'                       estimator = c("ols", "logit"))}
#'  
#' @export
gen.instrument <- function(formula,
                           Glist, 
                           data,
                           asymmetry = TRUE,
                           estimator = "ols",
                           power     = c(1, 1),
                           full      = FALSE,
                           nfold     = 2,
                           checkrank = TRUE,
                           tol       = 1e-10,
                           nthread   = 1,
                           drop      = NULL,
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
  if (!asymmetry & full){
    if (length(estimator) != 1) {
      stop("A single prediction method is required for the symmetric model.")
    }
  } else if (asymmetry) {
    if (length(estimator) == 1) {
      estimator    <- rep(estimator, 2)
    } 
    if (length(estimator) != 2) {
      stop("Two prediction methods (for intensive and extensive margins) are required for the asymmetric model.")
    }
  }
  
  estimatorint <- NULL
  estimatorext <- NULL
  if (length(estimator) >= 1) {
    if (tolower(estimator[1]) %in% c("lin", "linear", "ols", "lm")) {
      estimatorint <- "OLS"
    }  else if (tolower(estimator[1]) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
      estimatorint <- "Random Forest"
    } else if (tolower(estimator[1]) %in% c("xgboost")) {
      estimatorint <- "XGBoost"
    } else if (tolower(estimator[1]) %in% c("lasso")) {
      estimatorint <- "LASSO"
    } else {
      stop("This estimator is not available.")
    }
  }
  
  if (length(estimator) == 2) {
    if (tolower(estimator[2]) %in% c("lin", "linear", "ols", "lm")) {
      estimatorext <- "OLS"
    } else if (tolower(estimator[2]) %in% c("logit", "logistic")) {
      estimatorext <- "Logit"
    } else if (tolower(estimator[2]) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
      estimatorext <- "Random Forest"
    } else if (tolower(estimator[2]) %in% c("xgboost")) {
      estimatorext <- "XGBoost"
    } else if (tolower(estimator[2]) %in% c("lasso")) {
      estimatorext <- "LASSO"
    } else {
      stop("This estimator is not available.")
    }
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
  
  ### ybar
  endo   <- highlowstat1(X = as.matrix(y), G = Glist, cumsn = cumsn, nvec = nvec, 
                         ngroup = S, nthread = nthread)
  ybar   <- as.numeric(endo$Xbar)
  
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
  insyBar    <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                             power = power[1], nthread = nthread)
  insyBar_cn <- sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, xname))
  
  ### Folds construction for individuals
  if (nfold > S) {
    nfold  <- S
    warning("The number of folds exceeds the number of subnets; it has been reset to the number of subnets.")
  } 
  if (nfold == 1) {
    stop("At lead two folds is required.")
  }
  group      <- rep(0:(S - 1), nvec)
  # Only for non isolated but group should take value from 0, 1, 2, ... without jump
  id_fold_i  <- fassignfold(group, nfold = nfold)
  
  ### Instrument for ycheck
  out        <- NULL
  if(asymmetry){
    Xtp      <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                             power = power[2], nthread = nthread)
    
    ## Dyadic dada
    IDi      <- unlist(lapply(1:S, \(s) 0:(nvec[s] - 1)))
    gij      <- lapply(1:n, \(i) Glist[[group[i] + 1]][IDi[i] + 1, idpeer[[i]] + 1])
    ddni     <- sapply(gij, length)
    ddncs    <- c(0, cumsum(ddni))
    tp       <- fdataML(y = y, X = Xtp, group = group, IDi = IDi, gij = gij,
                        idpeer = idpeer, ddni = ddni, ddncs = ddncs, ncs = cumsn)
    ddyext   <- as.integer(tp$ddy[, 7]) # extensive margin
    ddyint   <- tp$ddy[, 6] - tp$ddy[, 5] #intensive margin
    ddX      <- tp$ddXj - tp$ddXi
    
    ## Check rank of ddX
    ddX      <- as.data.frame(ddX[, fcheckrank(X = ddX, tol = tol) + 1, drop = FALSE])
    
    ## Fold construction for dyadic
    nfold      <- as.integer(nfold)
    id_fold_d  <- fassignfold(tp$ddy[, 1], nfold = nfold)
    
    ## Prediction
    ARG      <- list(ddyext = ddyext, ddyint = ddyint, ddX = ddX, id_fold = id_fold_d, 
                     estimatorext = estimatorext, estimatorint = estimatorint, 
                     nthread = nthread, ...)
    
    insChey  <- fInstChecky(rhoddX = as.matrix(do.call(mpredict_ch, ARG)), 
                            ddni = ddni, nthread = nthread)
    
    GinsChey    <- peeravgpower(G = Glist, V = insChey, cumsn = cumsn, nvec = nvec, 
                                power = power[1], nthread = nthread)
    
    GinsChey_cn <- sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, "y_check_hat"))
    
    if (full){ # In this case we also predict Gy exogenously
      X_for_yb  <- cbind(insyBar, GinsChey)
      X_for_yb  <- as.data.frame(X_for_yb[, fcheckrank(X = X_for_yb, tol = tol) + 1, drop = FALSE])
      colnames(X_for_yb) <- paste0("X", 1:ncol(X_for_yb))
      
      ARG       <- list(Gy = ybar, X = X_for_yb, id_fold = id_fold_i, 
                        estimatorint = estimatorint, nthread = nthread, ...)
      insyBar   <- do.call(mpredict_bar, ARG)
      out       <- cbind(insyBar, insChey)
      colnames(out) <- c("y_bar_hat", "y_check_hat")
    } else {
      out       <- cbind(insyBar, insChey, GinsChey)
      colnames(out) <- c(insyBar_cn, "y_check_hat",  GinsChey_cn)
    }
    
  } else {
    
    if (full){ # In this case we also predict Gy exogenously
      X_for_yb  <- insyBar
      X_for_yb  <- as.data.frame(X_for_yb[, fcheckrank(X = X_for_yb, tol = tol) + 1, drop = FALSE])
      colnames(X_for_yb) <- paste0("X", 1:ncol(X_for_yb))
      
      ARG       <- list(Gy = ybar, X = X_for_yb, id_fold = id_fold_i, 
                        estimatorint = estimatorint, nthread = nthread, ...)
      insyBar   <- do.call(mpredict_bar, ARG)
      out       <- as.matrix(insyBar)
      colnames(out) <- "y_bar_hat"
    } else {
      out       <- insyBar
      colnames(out) <- insyBar_cn
    }
    
  }
  
  if (checkrank) {
    keepcol    <- fcheckrank(X = out[keep, , drop = FALSE], tol = tol) + 1
    out        <- out[, keepcol, drop = FALSE]
  }
  out[!keep, ] <- 0
  list(model.info  = list(power = power, estimator = c(estimatorint, estimatorext), 
                          nfold = nfold, tol = tol, full = full),
       instruments = out)
}

#' @rdname gen.instrument
#' @export
gen.instruments <- function(formula,
                            Glist, 
                            data,
                            asymmetry = TRUE,
                            estimator = c("ols", "logit"),
                            power     = c(1, 1),
                            full      = FALSE,
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
              estimator = estimator, power = power, full = full, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread, drop = drop, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.insts <- function(formula,
                      Glist, 
                      data,
                      asymmetry = TRUE,
                      estimator = c("ols", "logit"),
                      power     = c(1, 1),
                      full      = FALSE,
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
              estimator = estimator, power = power, full = full, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread, drop = drop, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.inst <- function(formula,
                     Glist, 
                     data,
                     asymmetry = TRUE,
                     estimator = c("ols", "logit"),
                     power     = c(1, 1),
                     full      = FALSE,
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
              estimator = estimator, power = power, full = full, nfold = nfold, 
              checkrank = checkrank, tol = tol, nthread = nthread, drop = drop,...)
  do.call(gen.instrument, ARG)
}