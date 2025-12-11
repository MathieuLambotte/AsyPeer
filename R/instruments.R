#' @export
gen.instrument <- function(formula,
                           Glist, 
                           data,
                           model     = "linear",
                           power     = c(1, 1),
                           sepiso    = TRUE,
                           diffX     = TRUE,
                           nfold     = 2,
                           checkrank = TRUE,
                           tol       = 1e-10,
                           nthread   = 1,
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
  
  ## Model
  if (tolower(model) %in% c("lin", "linear", "ols", "lm")) {
    model   <- "ols"
  } else if (tolower(model) %in% c("logit", "logistic", "glm")) {
    model   <- "glm"
  } else if (tolower(model) %in% c("rf", "r-f", "random forest", "random-forest", "randomforest")) {
    model   <- "RF"
  } else {
    stop("This model is not available.")
  }
  
  ## Network
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  Glist = fGnormalise(Glist, nthread)
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
  if (length(as.formula(formula)) != 3) {
    stop("formula is expected to be y ~ x1 + x2 + ...")
  } 
  f.t.data <- formula.to.data(formula = formula, data = data, fixed.effects = TRUE,
                              simulations = FALSE)
  ### We sta rt with instruments fro check y
  ## Row data
  y        <- f.t.data$y
  X        <- f.t.data$X
  Kx       <- ncol(X)
  xname    <- f.t.data$xname
  yname    <- f.t.data$yname
  if (sepiso) {
    X      <- cbind(1 - dg, X * (1 - dg), dg, X * dg)
  } else {
    X      <- cbind(1, X)
  }
  xname    <- c("iso", paste0("iso_", xname), "niso", paste0("niso_", xname))
  Xtp      <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                           power = power[2], nthread = nthread)
  ## Dyadic dada
  group    <- rep(0:(S - 1), nvec)
  IDi      <- unlist(lapply(1:S, \(s) 0:(nvec[s] - 1)))
  gij      <- lapply(1:n, \(i) Glist[[group[i] + 1]][IDi[i] + 1, idpeer[[i]] + 1])
  ddni     <- sapply(gij, length)
  ddncs    <- c(0, cumsum(ddni))
  
  tp       <- fdataML(y = y, X = Xtp, group = group, IDi = IDi, gij = gij,
                      idpeer = idpeer, ddni = ddni, ddncs = ddncs, ncs = cumsn,
                      nthread = nthread)
  
  ddy      <- tp$ddy[,7]
  ddX      <- NULL
  if (diffX) {
    ddX    <- tp$ddXj - tp$ddXi
  } else {
    ddX    <- cbind(tp$ddXi, tp$ddXj)
  }
  
  ## Check rank of ddX
  ddX      <- ddX[, fcheckrank(X = ddX, tol = tol) + 1, drop = FALSE]
  
  ## Fold construction
  if (nfold > S) {
    nfold  <- S
    warning("The number of folds exceeds the number of subnets; it has been reset to the number of subnets.")
  }
  id_fold  <- fassignfold(tp$ddy[,1], nfold = nfold)
  
  ## Prediction
  ARG      <- list(ddy = ddy, ddX = ddX, id_fold = id_fold, model = model, 
                   nthread = nthread, ...)
  rhoddX   <- tp$ddy[,4] * do.call(mpredict, ARG) * (tp$ddXj - tp$ddXi)
  insChey  <- fInstChecky(rhoddX = rhoddX, ddni = ddni, nthread = nthread)
  
  ### Instrument for 
  insBary  <- peeravgpower(G = Glist, V = X, cumsn = cumsn, nvec = nvec, 
                           power = power[1], nthread = nthread)
  
  out      <- cbind(insBary, insChey)
  colnames(out) <- c(c(sapply(paste0("G", ifelse(1:power[1] == 1, "", 1:power[1]), "_"), \(x) paste0(x, xname))),
                     c(sapply(paste0("rhoG", ifelse(1:power[2] == 1, "", 1:power[2]), "_"), \(x) paste0(x, xname))))
  
  out[, fcheckrank(X = out, tol = tol) + 1, drop = FALSE]
}

#' @rdname gen.instrument
#' @export
gen.instruments <- function(formula,
                            Glist, 
                            data,
                            model     = "linear",
                            power     = c(1, 1),
                            sepiso    = TRUE,
                            diffX     = TRUE,
                            nfold     = 2,
                            checkrank = TRUE,
                            tol       = 1e-10,
                            nthread   = 1,
                            ...) { 
  ARG <- list(formula = formula, Glist = Glist, data = data,  model = model,
              power = power, sepiso = sepiso, diffX = diffX, nfold = nfold,
              checkrank = checkrank, tol = tol, nthread = nthread, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.insts <- function(formula,
                      Glist, 
                      data,
                      model     = "linear",
                      power     = c(1, 1),
                      sepiso    = TRUE,
                      diffX     = TRUE,
                      nfold     = 2,
                      checkrank = TRUE,
                      tol       = 1e-10,
                      nthread   = 1,
                      ...) { 
  ARG <- list(formula = formula, Glist = Glist, data = data,  model = model,
              power = power, sepiso = sepiso, diffX = diffX, nfold = nfold,
              checkrank = checkrank, tol = tol, nthread = nthread, ...)
  do.call(gen.instrument, ARG)
}

#' @rdname gen.instrument
#' @export
gen.inst <- function(formula,
                     Glist, 
                     data,
                     model     = "linear",
                     power     = c(1, 1),
                     sepiso    = TRUE,
                     diffX     = TRUE,
                     nfold     = 2,
                     checkrank = TRUE,
                     tol       = 1e-10,
                     nthread   = 1,
                     ...) { 
  ARG <- list(formula = formula, Glist = Glist, data = data,  model = model,
              power = power, sepiso = sepiso, diffX = diffX, nfold = nfold,
              checkrank = checkrank, tol = tol, nthread = nthread, ...)
  do.call(gen.instrument, ARG)
}