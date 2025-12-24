#' @importFrom formula.tools env
#' @importFrom stats model.frame terms model.response model.matrix as.formula lm glm runif predict binomial
#' @importFrom ranger ranger
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel  
#' @importFrom doRNG registerDoRNG "%dorng%"
#' @importFrom foreach foreach
fdrop <- function(drop, ldg, nvec, S, y, X, Z, endo) {
  n        <- sum(nvec)
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  ncs      <- c(0, cumsum(nvec))
  lkeep    <- lapply(1:S, function(s) drop[(ncs[s] + 1):ncs[s + 1]] != 1)
  keep     <- unlist(lkeep)
  gkeep    <- sapply(1:S, function(s) sum(lkeep[[s]]) >= 1) # Groups We keep
  ldg      <- lapply(1:S, function(s) ldg[[s]][lkeep[[s]]])[gkeep]
  dg       <- unlist(ldg)
  S        <- length(ldg)
  nvec     <- sapply(ldg, length)
  n        <- sum(nvec)
  ncs      <- c(0, cumsum(nvec))
  SIs      <- sum(sapply(ldg, function(s) any(s == 0)))
  SnIs     <- sum(sapply(ldg, function(s) any(s != 0)))
  lIs      <- lapply(1:S, function(s) which(ldg[[s]] == 0) - 1 + ncs[s])
  Is       <- unlist(lIs)
  lnIs     <- lapply(1:S, function(s) which(ldg[[s]] != 0) - 1 + ncs[s])
  nIs      <- unlist(lnIs)
  y        <- y[keep]
  X        <- X[keep, , drop = FALSE]
  Z        <- Z[keep, , drop = FALSE]
  endo     <- endo[keep, , drop = FALSE]
  list(dg = dg, S = S, nvec = nvec, n = n,
       Is = Is, nIs = nIs, lIs = lIs, lnIs = lnIs, SIs = SIs, 
       SnIs = SnIs,  y = y, X = X, Z = Z, endo = endo)
}

fcheckrank <- function(X, tol = 1e-10) {
  which(fcheckrankEigen(X, tol)) - 1
}

mpredict  <- function(ddy, ddX, id_fold, estimator, nthread, ...){
  #Given a vector of fold id, create a list of the corresponding row of each ddyad
  # belonging in each fold
  id_list <- split(seq_along(id_fold), id_fold)
  seed    <- as.integer(runif(1, 0, 1e9))
  
  cl      <- makeCluster(nthread)
  registerDoParallel(cl)
  registerDoRNG(seed)
  lrho    <- foreach(k         = id_list, 
                     .export   = "mpredict_fold", #comment out
                     .packages = c("ranger", "AsyPeer") #Remember to add "NameOfThePackage"
  ) %dorng% {
    #each observation in fold k is predicted using a model trained
    #on the observations of the other folds
    ARG <- list(ddX = ddX, ddy = ddy, id_listk = k, estimator = estimator, ...)
    do.call(mpredict_fold, ARG) 
  }
  stopCluster(cl)
  
  rho   <- numeric(nrow(ddX))
  for (k in 1:length(id_list)) {
    rho[id_list[[k]]] <- lrho[[k]]
  }
  
  return(rho)
}


mpredict_fold <-function(ddX, ddy, id_listk, estimator, ...){
  #gather the observations from the other folds, expect id_listk
  ddX_train <- data.frame(ddX[-id_listk, ,drop = FALSE])
  ddy_train <- ddy[-id_listk]
  
  #gather the observations from the fold k
  ddX_k <- data.frame(ddX[id_listk, , drop = FALSE])
  if (estimator == "ols") {
    ARG         <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    model_train <- do.call(lm, ARG) 
    rho_k       <- predict(model_train, newdata = ddX_k)
  } else if (estimator == "glm") {
    ARG           <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    if (is.null(ARG$family)) { # If the user does not set link, use logit
      ARG$family  <- binomial(link = "logit")
    }
    model_train <- do.call(glm, ARG)  
    rho_k       <- predict(model_train, newdata = ddX_k, type = "response")
  } else if (estimator == "RF") {
    ddy_train   <- as.factor(ddy_train)
    ARG         <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    model_train <- do.call(ranger, ARG)  
    rho_k       <- predict(model_train, data=ddX_k)$predictions
  }
  return(rho_k)
}

formula2data <- function(formula,
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
    if(length(formula) != 2) stop("The `formula` argument is invalid. The expected format is `~ x1 + x2 + ...`.")
  } else {
    if(length(formula) != 3) stop("The `formula` argument is invalid. The expected format is `y ~ x1 + x2 + ...`.")
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
  idpeer   <- do.call(c, lapply(1:S, function(s) lapply(1:nvec[s], function(i) which(Glist[[s]][i,] > 0) - 1)))
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

#' @importFrom stats pnorm
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

testSargan <-function(object, y, X_iso, X_niso, endo, Z,
                      S, cumsn, lIso, lnIso, weight, HACn){
  ## degree of freedom
  Iso    <- unlist(lIso)
  nIso   <- unlist(lnIso)
  dfiso  <- object$model.info$dfiso
  dfniso <- object$model.info$dfniso
  spillover <- object$model.info$spillover
  cgamma    <- which(object$model.info$xname %in% object$model.info$common.gamma) - 1
  ncgamma   <- which(object$model.info$xname %in% object$model.info$ncommon.gamma) - 1
  # GMM OLS
  Z      <- cbind(endo, Z)
  Z      <- Z[, fcheckrank(X = Z, tol = object$model.info$tol) + 1, 
              drop = FALSE]
  
  if (weight == "I"){
    W    <- diag(nrow = ncol(Z))
  } else if(weight %in% c("IV", "optimal")){ 
    W    <- solve(crossprod(Z) / S)
  }
  
  fGMM     <- gmm_obj_nospil
  fWopt    <- W_optimal_nospil
  fest     <- compute_estimate_nospil
  if (spillover){
    fGMM   <- gmm_obj
    fWopt  <- W_optimal
    fest   <- compute_estimate
  }
  
  gmm    <- optimize(f = fGMM, Z = Z, y = y, endo = endo, X_iso = X_iso, 
                     X_niso = X_niso, W = W, S = S, c_gamma = cgamma, 
                     nc_gamma = ncgamma, lower = -0.999, upper =  20)
  betal  <- gmm$minimum
  
  if(weight == "optimal"){
    W    <- fWopt(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                  X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, cumsn = cumsn,
                  dfiso = dfiso, dfniso = dfniso, HAC = HACn, S = S,
                  c_gamma = cgamma, nc_gamma = ncgamma)
    
    #gmm with optimal W
    gmm  <- optimize(f = fGMM, Z = Z, y = y, endo = endo, X_iso = X_iso, 
                     X_niso = X_niso, W = W, S = S, c_gamma = cgamma, 
                     nc_gamma = ncgamma, lower = -0.999, upper =  20)
    
    #get the optimal estimate of beta_l
    betal <- gmm$minimum
  }
  ols     <- fest(betal = betal, Z = Z, y = y, endo = endo, X_iso = X_iso,
                  X_niso = X_niso, W = W, Iso = Iso, nIso = nIso, 
                  cumsn = cumsn, dfiso = dfiso, dfniso = dfniso, HAC = HACn, 
                  c_gamma = cgamma, nc_gamma = ncgamma, S = S)
  
  #J stat
  Jstat      <- ols$JStat - object$gmm$Sargan$stat
  df         <- ols$Jdf - object$gmm$Sargan$df
  pvalue     <- 1 - pchisq(Jstat, df)
  out        <- c(Jstat, df, pvalue)
  print(ols)
  names(out) <- c("Jstat", "df", "pvalue")
  return(out)
}

#' @importFrom stats pf
fdiagnostic <- function(object, KPtest, nthread) {
  fixed.effects <- object$model.info$fixed.effects
  spillover <- object$model.info$spillover
  asymmetry <- object$model.info$asymmetry
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
    y       <- c(Demean_separate(X = y, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                                 nthread = nthread))
    endo    <- Demean_separate(X = endo, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                               nthread = nthread)
    X_iso   <- Demean_separate(X = X_iso, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                               nthread = nthread)
    X_niso  <- Demean_separate(X = X_niso, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                               nthread = nthread)
    Z       <- Demean_separate(X = Z, cumsn = cumsn, lIso = lIso, lnIso = lnIso, 
                               nthread = nthread)
  }
  index  <- which(!(object$model.info$zname %in% 
                      c(paste0("iso_", object$model.info$xname), 
                        paste0("niso_", object$model.info$xname)))) - 1
  
  ## Weak instrument test
  tpF    <- fFstat(y = endo, X = Z, index = index, cumsn = cumsn, HAC = HACn, 
                   nthread = nthread)
  tpKP    <- NULL
  if (KPtest) {
    xname <- object$model.info$xname
    zname <- object$model.info$zname
    X     <- cbind(X_iso, X_niso)[, c(paste0("iso_", xname), paste0("niso_", xname)) %in% zname, drop = FALSE]
    tpKP  <- fKPstat(endo_ = endo, X = X, Z_ = Z, index = index, cumsn = cumsn, 
                     HAC = HACn)
    tpKP$pvalue <- pchisq(tpKP$stat, tpKP$df, lower.tail = FALSE)
  }
  
  ## Endogeneity test
  # tpend  <- testSargan(object = object, y = y, X_iso = X_iso, X_niso = X_niso, 
  #                      endo = endo, Z = Z, S = S, cumsn = cumsn, lIso = lIso,
  #                      lnIso = lnIso, weight = weight, HACn = HACn)
  out    <- cbind(df1        = unlist(c(rep(tpF$df1, 1 + asymmetry), tpKP$df, 
                                        object$gmm$Sargan["df"])),
                  df2        = c(rep(tpF$df2, 1 + asymmetry), rep(NA, 1 + KPtest)),
                  statistic  = unlist(c(tpF$F, tpKP$stat, object$gmm$Sargan["stat"])),
                  "p-value"  = unlist(c(rep(NA, 1 + asymmetry), tpKP["pvalue"], object$gmm$Sargan["pvalue"])))
  out[1:(1 + asymmetry), 4]   <- pf(out[1:(1 + asymmetry), 3], out[1:(1 + asymmetry), 1], out[1:(1 + asymmetry), 2], lower.tail = FALSE)
  rn            <- if (asymmetry) {
    paste0("Weak instruments (", c("ybar", "ycheck"), ")")
  } else {
    "Weak instruments"
  }
  if (KPtest) {
    rn          <- c(rn, "Kleibergen-Paap rk Wald", "Sargan J")
  } else {
    rn          <- c(rn, "Sargan J")
  }
  rownames(out) <- rn
  out
}
