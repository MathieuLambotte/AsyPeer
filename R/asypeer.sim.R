#' @title Simulating Asymmetric Peer Effect Models with Continuous Outcome
#' @param formula A formula object (\link[stats]{formula}): a symbolic description of the model. `formula` should be specified as, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are control variables, which can include contextual variables such as averages or quantiles among peers.
#' @param Glist The adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` must be a list of subnets, with the `s`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `s`-th subnet.
#' @param parms A vector defining the true values of \eqn{(\delta, \boldsymbol{\beta}', \boldsymbol{\gamma}')'}, where \eqn{\delta} measure spillover effects, 
#' \eqn{\boldsymbol{\beta} = (\beta^-, \beta^+)'} captures asymmetric conformity effects, and \eqn{\boldsymbol{\gamma}} is the parameter vector associated with 
#' the exogenous characteristics in agent types. 
#' The parameters \eqn{\delta}, \eqn{\boldsymbol{\beta}}, and \eqn{\boldsymbol{\gamma}} can also be specified separately using the arguments `delta`, `beta`, and `gamma`.
#' @param delta The true value of the spillover effect parameter.
#' @param beta The true value of the asymmetric conformity parameters.
#' @param gamma The true value of the vector \eqn{\gamma}.
#' @param epsilon A vector of idiosyncratic error terms. If not specified, it will be simulated from a standard normal distribution (see the model specification in the Details section of \code{\link{cespeer}}). 
#' @param maxit The maximum number of iterations for the Fixed Point Iteration Method.
#' @param data An optional data frame, list, or environment containing the model variables. If a variable is not found in `data`, it is retrieved from \code{environment(formula)}, typically the environment from which `cespeer.sim` is called.
#' @param tol The tolerance value used in the Fixed Point Iteration Method to compute the outcome `y`. The process stops if the \eqn{\ell_1}-distance 
#' between two consecutive values of `y` is less than `tol`.
#' @param init An optional initial guess for the equilibrium.
#' @param nthread Number of CPU cores (threads) used to run parts of the simulation in parallel.
#' @description
#' `asypeer.sim` simulates asymmetric peer effect models with continuous responses.
#' @return A list containing:
#'     \item{y}{The simulated variable.}
#'     \item{epsilon}{The idiosyncratic error.}
#'     \item{init}{The initial guess.}
#'     \item{iteration}{The number of iterations before convergence.}
#' @examples 
#' @export
asypeer.sim <- function(formula, Glist, parms, beta, delta, gamma, epsilon, init, 
                        tol = 1e-10, maxit = 500, nthread = 1, data){
  # Network
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  dg       <- fnetwork(Glist = Glist)
  S        <- dg$S
  nvec     <- dg$nvec
  n        <- dg$n
  igr      <- dg$igr
  cumsn    <- dg$cumsn
  Is       <- dg$Is
  nIs      <- dg$nIs
  ldg      <- dg$ldg
  idpeer   <- dg$idpeer
  dg       <- dg$dg
  
  # Data
  f.t.data <- formula.to.data(formula = formula, data = data, simulations = TRUE, fixed.effects = FALSE)
  formula  <- f.t.data$formula
  X        <- f.t.data$X
  if (nrow(X) != n) stop("The number of observations does not match the number of nodes in the network.")
  Kx       <- ncol(X)
  eps      <- NULL
  if(missing(epsilon)){
    eps    <- rnorm(n)
  } else{
    eps    <- c(epsilon)
    if (!(length(eps) %in% c(1, n))) stop("`epsilon` must be either a scalar or an n-dimensional vector.")
    if (length(eps) == 1) eps <- rep(eps, n)
  }
  
  # parameters
  del    <- NULL
  bet    <- NULL
  gam    <- NULL
  if (missing(parms)) {
    if (missing(delta) | missing(beta) | missing(gamma)) {
      stop("Define either `parms` or `delta`, `beta`, and `gamma`.")
    }
    if (length(beta) != 2){
      stop("length(beta) is different from 2. See details on the model.")
    }
    del <- delta
    bet <- beta
    if (length(gamma) != Kx) stop("length(gamma) is different from ncol(X).")
    gam <- gamma
  } else{
    if (!missing(delta) | !missing(beta) | !missing(gamma)) {
      stop("Define either `parms` or `delta`, `beta`, and `gamma`.")
    }
    if (length(parms) != (3 + Kx)) stop("length(parms) is different from 3 + ncol(X). See details on the model.")
    bet   <- parms[1:2]
    del   <- parms[3]
    gam   <- tail(parms, Kx)
  }
  lambda  <- (del + bet)/(1 + beta)
  if (any(abs(lambda) >= max(dg))) {
    warning("The absolute value of the total peer effects is greater than or equal to one, which may lead to multiple or no equilibria.")
  }
  
  # Solving the game
  ## alpha
  alpha   <- X %*% gam + eps
  
  ## init
  if (missing(init)) {
    init   <- c(alpha)
  }
  if (length(init) == 1){
    init   <- rep(init, n)
  } else if (length(init) != n) {
    stop("`init` is not an n-vector.")
  }
  y        <- init + 0 # copy so that y not linked to init
  t        <- fNashE(y = y, alpha = alpha, G = Glist, peffects = c(bet, del), 
                     cumsn = cumsn, nvec = nvec, idpeer = idpeer, d = dg, 
                     ngroup = S, tol = tol, maxit = maxit, nthread = nthread)
  # Output
  list("y"         = c(y),
       "epsilon"   = eps,
       "init"      = init,
       "iteration" = t)
} 