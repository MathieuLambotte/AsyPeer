#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats as.formula
formula.to.data <- function(formula,
                            data, 
                            simulations   = FALSE,
                            fixed.effects = FALSE) {
  
  ## Extract data from the formula
  if (missing(data)) {
    data      <- env(formula)
  }
  formula     <- as.formula(formula)
  yname       <- NULL
  
  if (simulations) {
    if(length(formula) != 2) stop("The `formula` argument is invalid. For data simulation, the expected format is `~ X1 + X2 + ...`.")
  } else {
    if(length(formula) != 3) stop("The `formula` argument is invalid. For estimation, the expected format is `y ~ X1 + X2 + ...`.")
    yname     <- all.vars(formula)[1]
  }
  
  ## call model.frame()
  mf          <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y           <- model.response(mf, "numeric")
  X           <- model.matrix(terms(formula, data = data, rhs = 1), mf)
  xname       <- colnames(X)
  intercept   <- "(Intercept)" %in% xname
  if(fixed.effects & intercept){
    X         <- X[, xname != "(Intercept)", drop = FALSE]
    xname    <- xname[xname != "(Intercept)"]
    intercept <- FALSE
  }
  
  list("formula"   = formula, 
       "X"         = X, 
       "y"         = y,
       "intercept" = intercept,
       "yname"     = yname,
       "xname"     = xname)
}

#' @importFrom utils head
#' @importFrom utils tail
fnetwork   <- function(Glist) {
  # Isol is true isolated than can be removed. But this argument is no longer used. 
  # See the more general argument which is now drop
  S        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  cumsn    <- c(0, cumsum(nvec))
  
  ldg      <- lapply(Glist, function(g) round(apply(g, 1, sum), 7))
  idpeer   <- do.call(c, lapply(1:S, function(s) lapply(1:nvec[s], function(i) which(G[[s]][i,] > 0) - 1)))
  dg       <- unlist(ldg)
  SIs      <- round(sum(sapply(ldg, function(s) any(s == 0))))
  SnIs     <- round(sum(sapply(ldg, function(s) any(s != 0))))
  lIs      <- lapply(1:S, function(s) which(ldg[[s]] == 0) - 1 + cumsn[s])
  Is       <- unlist(lIs)
  lnIs     <- lapply(1:S, function(s) which(ldg[[s]] != 0) - 1 + cumsn[s])
  nIs      <- unlist(lnIs)
  
  list(dg = dg, ldg = ldg, S = S, nvec = nvec, n = n, cumsn = cumsn, Is = Is, nIs = nIs, 
       lIs = lIs, lnIs = lnIs, SIs = SIs, SnIs = SnIs, idpeer = idpeer)
}

fcheckrank <- function(X, tol = 1e-10) {
  which(fcheckrankEigen(X, tol)) - 1
}

fcoef           <- function(Estimate, cov) {
  coef           <- cbind(Estimate, sqrt(diag(cov)), 0, 0)
  coef[,3]       <- coef[,1]/coef[,2]
  coef[,4]       <- 2*(1 - pnorm(abs(coef[,3])))
  colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coef
}

fprintcoeft <- function(coef) {
  pval      <- coef[,ncol(coef)]
  pval_pt   <- sapply(pval, function(s){ifelse(is.na(s), "NA", ifelse(s < 2e-16, "<2e-16", format(s, digit = 4)))})
  refprob   <- c(0.001, 0.01, 0.05, 0.1)
  refstr    <- c("***",  "**", "*", ".", "")
  str       <- sapply(pval, function(s) ifelse(is.na(s), "", refstr[1 + sum(s > refprob)]))
  out       <- data.frame(coef[,-ncol(coef), drop = FALSE], "P" = pval_pt, "S" = str); 
  colnames(out) <- c(colnames(coef), "")
  print(out)
}

testStargan <-function(object, y, X_iso, X_niso, endo, Z,
                      S, cumsn, lIso, lnIso, weight, HACn){
  ## degree of freedom
  Iso    <- unlist(lIso)
  nIso   <- unlist(lnIso)
  dfiso  <- object$model.info$dfiso
  dfniso <- object$model.info$dfniso

  # GMM OLS
  Z      <- cbind(endo, Z)
  Z      <- Z[, fcheckrank(X = Z, tol = object$model.info$tol) + 1, 
              drop = FALSE]
  
  if (weight == "I"){
    W    <- diag(nrow = ncol(Z))
  } else if(weight %in% c("IV", "optimal")){ 
    W    <- solve(crossprod(Z) / S)
  }
  gmm    <- optimize(f = gmm_obj, Z = Z, y = y, endo = endo, X_iso = X_iso, 
                     X_niso = X_niso, W = W, S = S, lower = -0.999, upper =  20)
  betal  <- gmm$minimum
  
  if(weight == "optimal"){
    W    <- W_optimal(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                      X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, cumsn = cumsn,
                      dfiso = dfiso, dfniso = dfniso, HAC = HACn, S = S)
    
    #gmm with optimal W
    gmm  <- optimize(f = gmm_obj, Z = Z, y = y, endo = endo, X_iso = X_iso, 
                     X_niso = X_niso, W = W, S = S, lower = -0.999, upper =  20)
    
    #get the optimal estimate of beta_l
    betal <- gmm$minimum
  }
  ols     <- compute_estimate(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                          X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, 
                          cumsn = cumsn, dfiso = dfiso, dfniso = dfniso, HAC = HACn, S = S)
  
  #J stat
  Jstat      <- ols$JStat - object$gmm$Sargan$stat
  df         <- ols$Jdf - object$gmm$Sargan$df
  pvalue     <- 1 - pchisq(Jstat, df)
  out        <- c(Jstat, df, pvalue)
  names(out) <- c("Jstat", "df", "pvalue")
  return(out)
}

#' @importFrom stats pf
fdiagnostic <- function(object, KPtest, nthread) {
  fixed.effects <- object$model.info$fixed.effects
  nvec      <- object$model.info$nvec
  lIso      <- lapply(object$data$isolates, \(x) x - 1)
  lnIso     <- lapply(object$data$non.isolates, \(x) x - 1)
  HAC       <- object$model.info$HAC
  HACn      <- (0:3)[HAC == c("iid","group-iid", "hetero", "cluster")]
  cumsn     <- c(0, cumsum(nvec))
  y         <- as.matrix(object$data$dependent)
  endo      <- object$data$endogenous.variables
  dg        <- object$data$degree
  X_iso     <- object$data$exogenous * (1 - dg)
  X_niso    <- object$data$exogenous * dg
  Kx        <- length(object$model.info$xname)
  Z         <- object$data$instruments
  weight    <- object$model.info$weight
  S         <- object$model.info$ngroup
  if (fixed.effects) {
    y       <- c(Demean(X = y, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                        nthread = nthread))
    endo    <- Demean(X = endo, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                      nthread = nthread)
    X_iso   <- Demean(X = X_iso, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                      nthread = nthread)
    X_niso  <- Demean(X = X_niso, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                      nthread = nthread)
    Z       <- Demean(X = Z, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                      nthread = nthread)
  }
  index  <- which(!(object$model.info$zname %in% 
                      c(paste0("iso_", object$model.info$xname), 
                        paste0("niso_", object$model.info$xname)))) - 1
  
  ## Weak instrument test
  tpF    <- fFstat(y = endo, X = Z, index = index, cumsn = cumsn, HAC = HACn, 
                   nthread = nthread)
  tpKP   <- NULL
  if (KPtest) {
    tpKP <- fKPstat(endo_ = endo, X = X, Z_ = Z, index = index, cumsn = cumsn, 
                      HAC = HACn)
  }

  ## Endogeneity test
  tpend  <- testStargan(object = object, y = y, X_iso = X_iso, X_niso = X_niso, 
                        endo = endo, Z = Z, S = S, cumsn = cumsn, lIso = lIso,
                        lnIso = lnIso, weight = weight, HACn = HACn)
  out    <- cbind(df1        = unlist(c(rep(tpF$df1, 2), tpKP$df, tpend["df"], 
                                        object$gmm$Sargan["df"])),
                  df2        = c(rep(tpF$df2, 2), rep(NA, 2 + KPtest)),
                  statistic  = unlist(c(tpF$F, tpKP$stat, tpend["Jstat"], object$gmm$Sargan["stat"])),
                  "p-value"  = unlist(c(rep(NA, 2 + KPtest), tpend["pvalue"], object$gmm$Sargan["pvalue"])))
  out[1:2, 4]   <- pf(out[1:2, 3], out[1:2, 1], out[1:2, 2], lower.tail = FALSE)
  out[3, 4]     <- pchisq(out[3, 3], out[3, 1], lower.tail = FALSE)
  rn            <- paste0("Weak instruments (", c("ybar", "ydot"), ")")
  if (KPtest) {
    rn          <- c(rn, "Kleibergen-Paap rk Wald", "Hausman", "Sargan J")
  } else {
    rn          <- c(rn, "Hausman", "Sargan J")
  }
  rownames(out) <- rn
  out
}

