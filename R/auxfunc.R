#' @importFrom formula.tools env
#' @importFrom stats model.frame terms model.response model.matrix as.formula lm glm runif predict binomial
#' @importFrom ranger ranger
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel  
#' @importFrom doRNG registerDoRNG "%dorng%"
#' @importFrom foreach foreach registerDoSEQ
#' @importFrom xgboost xgb.train xgb.DMatrix xgb.params
#' @importFrom glmnet cv.glmnet
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

############ Function to predict ych
mpredict_ch  <- function(ddyext, ddyint, ddX, id_fold, estimatorext, 
                         estimatorint, nthread, DoRNGseed, ...){
  #Given a vector of fold id, create a list of the corresponding row of each ddyad
  # belonging in each fold
  id_list <- split(seq_along(id_fold), id_fold)
  
  # Arguments
  ARG    <- list(ddX = ddX, ddyext = ddyext, ddyint = ddyint, 
                 estimatorext = estimatorext, estimatorint = estimatorint, ...)
  
  # Prediction
  cl      <- NULL
  on.exit({
    registerDoSEQ()
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  cl      <- makeCluster(nthread)
  registerDoParallel(cl)
  registerDoRNG(DoRNGseed)
  
  lycheckh <- foreach(k         = id_list, 
                      .export   = "mpredict_fold_ch", #comment out
                      .packages = c("ranger", "xgboost", "glmnet", "AsyPeer") 
  ) %dorng% {
    #each observation in fold k is predicted using a model trained
    do.call(mpredict_fold_ch, c(ARG, list(id_listk = k)))
  }
  
  ycheckh  <- numeric(length(ddyext))
  for (k in 1:length(id_list)) {
    ycheckh[id_list[[k]]] <- lycheckh[[k]]
  }
  return(ycheckh)
}


mpredict_fold_ch <-function(ddX, ddyext, ddyint, id_listk, 
                            estimatorext, estimatorint, ...){
  exthat <- NULL
  inthat <- NULL
  dots   <- list(...)
  
  # gather the observations from the fold k
  ddX_k  <- ddX[id_listk, , drop = FALSE]
  
  ### Extensive margin
  # gather the observations from the other folds, expect id_listk
  ddX_train <- ddX[-id_listk, ,drop = FALSE]
  ddy_train <- ddyext[-id_listk]
  if (!all(c(0, 1) %in% ddy_train)) {
    stop("The response variable for the instrument construction must contain both y_j > y_i and y_j <= y_i in all folds.")
  }
  
  if (estimatorext == "OLS") {
    
    dotname     <- setdiff(names(formals(lm)), c("formula", "data", "..."))
    ARG         <- c(list(formula = ddy_train ~ ., data = ddX_train),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(lm, ARG)  
    exthat      <- predict(model_train, newdata = ddX_k)
    
  } else if (estimatorext == "Logit") {
    
    dotname     <- setdiff(names(formals(glm)), 
                           c("formula", "data", "family", "..."))
    ARG         <- c(list(formula = ddy_train ~ ., data = ddX_train, 
                          family = binomial(link = "logit")),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(glm, ARG)  
    exthat      <- predict(model_train, newdata = ddX_k, type = "response")
    
  } else if (estimatorext == "Random Forest") {
    
    ddy_train   <- factor(ddy_train, levels = c(0, 1))
    dotname     <- setdiff(names(formals(ranger)), 
                           c("formula", "data", "probability", "..."))
    ARG         <- c(list(formula = ddy_train ~ ., data = ddX_train, 
                          probability = TRUE), dots[names(dots) %in% dotname])
    model_train <- do.call(ranger, ARG)  
    exthat      <- predict(model_train, data = ddX_k)$predictions[,"1"]
    
  } else if (estimatorext == "XGBoost") {
    
    # Training data
    ddtrain     <- xgb.DMatrix(data = as.matrix(ddX_train),
                               label = ddy_train)
    # prediction
    ddpred      <- xgb.DMatrix(data = as.matrix(ddX_k))
    # parameters
    dotname     <- setdiff(names(formals(xgb.params)),  
                           c("objective", "nthread", "..."))
    ddpar       <- c(list(objective = "binary:logistic",
                          nthread = 1), 
                     dots[names(dots) %in% dotname])
    # Training argument
    dotname <- setdiff(names(formals(xgb.train)),  
                       c("params", "data", "objective", 
                         "verbose", "xgb_model", "evals", "..."))
    ARG     <- c(list(params = ddpar, data = ddtrain, verbose = 0),
                 dots[names(dots) %in% dotname])
    ARG$nrounds <- fassignnull(ARG$nrounds, 200)
    # Training
    model_train   <- do.call(xgb.train, ARG)
    # Prediction 
    exthat  <- predict(model_train, newdata = ddpred)
    
  } else if (estimatorext == "LASSO") {
    
    dotname     <- setdiff(names(formals(cv.glmnet)), 
                           c("y", "x", "family", "alpha", "..."))
    ARG         <- c(list(y = ddy_train, x = as.matrix(ddX_train), 
                          family = "binomial", alpha = 1),
                     dots[names(dots) %in% dotname])
    fitlasso    <- do.call(cv.glmnet, ARG) 
    exthat      <- as.numeric(predict(fitlasso, newx = as.matrix(ddX_k), s = "lambda.min", 
                                      type = "response"))
  }
  
  ### Intensive margin
  idx_pos   <- (ddyext[-id_listk] == 1)
  ddX_train <- ddX[-id_listk, , drop = FALSE][idx_pos, , drop = FALSE]
  ddy_train <- ddyint[-id_listk][idx_pos]
  
  if (estimatorint == "OLS") {
    
    dotname     <- setdiff(names(formals(lm)), c("formula", "data", "..."))
    ARG         <- c(list(formula = ddy_train ~ ., data = ddX_train),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(lm, ARG) 
    inthat      <- predict(model_train, newdata = ddX_k)
    
  } else if (estimatorint == "Random Forest") {
    
    dotname     <- setdiff(names(formals(ranger)), c("formula", "data", "..."))
    ARG         <- c(list(formula = ddy_train ~ ., data = ddX_train),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(ranger, ARG)  
    inthat      <- predict(model_train, data = ddX_k)$predictions
    
  } else if (estimatorint == "XGBoost") {
    
    # Training data
    ddtrain     <- xgb.DMatrix(data = as.matrix(ddX_train),
                               label = ddy_train)
    # prediction
    ddpred      <- xgb.DMatrix(data = as.matrix(ddX_k))
    # parameters
    dotname     <- setdiff(names(formals(xgb.params)),  
                           c("objective", "nthread", "..."))
    ddpar       <- c(list(objective = "reg:squarederror",
                          nthread = 1), 
                     dots[names(dots) %in% dotname])
    # Training arguments
    dotname <- setdiff(names(formals(xgb.train)),  
                       c("params", "data", "objective", 
                         "verbose", "xgb_model", "evals", "..."))
    ARG     <- c(list(params = ddpar, data = ddtrain, verbose = 0),
                 dots[names(dots) %in% dotname])
    ARG$nrounds <- fassignnull(ARG$nrounds, 200)
    # Training
    model_train   <- do.call(xgb.train, ARG)
    # Prediction 
    inthat <- predict(model_train, newdata = ddpred)

  } else if (estimatorint == "LASSO") {
    
    dotname     <- setdiff(names(formals(cv.glmnet)), c("y", "x", "alpha", "..."))
    ARG         <- c(list(y = ddy_train, x = as.matrix(ddX_train), alpha = 1),
                     dots[names(dots) %in% dotname])
    fitlasso    <- do.call(cv.glmnet, ARG) 
    inthat      <- as.numeric(predict(fitlasso, newx = as.matrix(ddX_k), 
                                      s = "lambda.min"))
  }
  
  return(exthat * inthat)
}


############ Function to predict ych
mpredict_bar  <- function(Gy, X, id_fold, estimatorint, nthread, DoRNGseed, ...){
  # Given a vector of fold id, create a list of the corresponding row of each ddyad
  # belonging in each fold
  id_list <- split(seq_along(id_fold), id_fold)
  
  # Prediction
  cl      <- NULL
  on.exit({
    registerDoSEQ()
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  cl      <- makeCluster(nthread)
  registerDoParallel(cl)
  registerDoRNG(DoRNGseed)
  
  ARG     <- list(Gy = Gy, X = X, estimatorint = estimatorint, ...)
  lybarh  <-  foreach(k         = id_list,
                      .export   = "mpredict_fold_bar", #comment out
                      .packages = c("ranger", "xgboost", "glmnet", "AsyPeer")
  ) %dorng% {
    #each observation in fold k is predicted using a model trained
    do.call(mpredict_fold_bar, c(ARG, list(id_listk = k)))
  }
  
  # Prediction
  ybarh   <- numeric(length(Gy))
  for (k in 1:length(id_list)) {
    ybarh[id_list[[k]]] <- lybarh[[k]]
  }

  return(ybarh)
}


mpredict_fold_bar <-function(Gy, X, id_listk, estimatorint, ...){
  ybarh  <- NULL
  dots   <- list(...)
  
  # Y for the train sample
  Gy_train      <- Gy[-id_listk]
  # X for fold k
  X_k           <- data.frame(X[id_listk, , drop = FALSE])
  # X for other folds
  X_train       <- data.frame(X[-id_listk, ,drop = FALSE])
  
  # Estimation 
  if (estimatorint == "OLS") {
    
    dotname     <- setdiff(names(formals(lm)), c("formula", "data", "..."))
    ARG         <- c(list(formula = Gy_train ~ ., data = X_train),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(lm, ARG) 
    ybarh       <- predict(model_train, newdata = X_k)
    
  } else if (estimatorint == "Random Forest") {
    
    dotname     <- setdiff(names(formals(ranger)), c("formula", "data", "..."))
    ARG         <- c(list(formula = Gy_train ~ ., data = X_train),
                     dots[names(dots) %in% dotname])
    model_train <- do.call(ranger, ARG)  
    ybarh       <- predict(model_train, data = X_k)$predictions
    
  } else if (estimatorint == "XGBoost") {
    
    # Training data
    dtrain      <- xgb.DMatrix(data = as.matrix(X_train),
                               label = Gy_train)
    # prediction
    dpred       <- xgb.DMatrix(data = as.matrix(X_k))
    # parameters
    dotname     <- setdiff(names(formals(xgb.params)),  
                           c("objective", "nthread", "..."))
    dpar       <- c(list(objective = "reg:squarederror",
                          nthread = 1), 
                     dots[names(dots) %in% dotname])
    # Training arguments
    dotname <- setdiff(names(formals(xgb.train)),  
                       c("params", "data", "objective", 
                         "verbose", "xgb_model", "evals", "..."))
    ARG     <- c(list(params = dpar, data = dtrain, verbose = 0),
                 dots[names(dots) %in% dotname])
    ARG$nrounds <- fassignnull(ARG$nrounds, 200)
    # Training
    model_train  <- do.call(xgb.train, ARG)
    # Prediction 
    ybarh <- predict(model_train, newdata = dpred)
    
  } else if (estimatorint == "LASSO") {
    
    dotname     <- setdiff(names(formals(cv.glmnet)), c("y", "x", "alpha", "..."))
    ARG         <- c(list(y = Gy_train, x = as.matrix(X_train), alpha = 1),
                     dots[names(dots) %in% dotname])
    fitlasso    <- do.call(cv.glmnet, ARG) 
    ybarh       <- as.numeric(predict(fitlasso, newx = as.matrix(X_k), 
                                      s = "lambda.min"))
  }
  
  return(ybarh)
}


########## Transform formulat to data
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

#' @importFrom stats pf
fdiagnostic <- function(object, KP, SW, nthread) {
  fixed.effects <- object$model.info$fixed.effects
  spillover <- object$model.info$spillover
  asymmetry <- object$model.info$asymmetry
  nvec      <- object$model.info$nvec
  lIso      <- lapply(object$data$isolates, \(x) x - 1)
  lnIso     <- lapply(object$data$non.isolates, \(x) x - 1)
  HAC       <- object$model.info$HAC
  HACn      <- c(0:3, 3)[HAC == c("iid","group-iid", "hetero", "cluster", "cluster (bootstrap)")]
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
  
  # Number of included instruments
  Kinc   <- ncol(Z) - length(index)
  
  ## Weak instrument test
  if (!asymmetry) {
    SW  <- FALSE
  }
  
  tpF   <- NULL
  
  if (SW) {
    tpF <- fswstat(endo = endo, Z = Z, K = Kinc, nthread = nthread)
    rownames(tpF) <- paste0("Weak instruments - Cond. F (", c("ybar", "ycheck"), ")")
  } else {
    tpF <- fFstat(endo = endo, Z = Z, K = Kinc, cumsn = cumsn, 
                  cluster = (HACn == 3), nthread = nthread)
    rn  <- if (asymmetry) {
      paste0("Weak instruments F (", c("ybar", "ycheck"), ")")
    } else {
      "Weak instruments F"
    }
    rownames(tpF) <- rn
  }
  out   <- tpF
  
  # LM test KP
  tpKP  <- NULL
  if (KP) {
    tpKP <- fKPstat(endo = endo, Z = Z, K = Kinc, cumsn = cumsn, 
                    cluster = (HACn == 3))
    rownames(tpKP) <- "Kleibergen-Paap rk Wald (LM)"
    out  <- rbind(out, tpKP)
  }
  
  ### Sargan
  tpS    <- matrix(c(object$gmm$Sargan$df, NA, 
                     object$gmm$Sargan$stat, object$gmm$Sargan$pvalue), 1)
  rownames(tpS) <- "Sargan J"
  out    <- rbind(out, tpS)
  
  colnames(out) <- c("df1", "df2", "statistic", "p-value")
  out
}


fCESdatainit  <- function (y, z, G, nvec, S, ldg, lIs, lnIs, drop) {
  n           <- sum(nvec)
  if (length(drop) == 0) {
    drop      <- rep(0, n)
  }
  if (any(!(drop %in% 0:1) | !is.finite(drop))) {
    stop("`drop` must be a binary (0/1) variable.")
  }
  if (length(drop) != n) {
    stop("`drop` must be a vector of length n.")
  }
  ncs         <- c(0, cumsum(nvec))
  friendindex <- lapply(1:S, function(m) {
    lapply(1:nvec[m], function(s) {
      which(G[[m]][s,] > 0) - 1
    })})
  frzeroy     <- as.integer(unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      any(y[friendindex[[m]][[s]] + ncs[m] + 1] <= 0)
    })})))
  frzeroz     <- as.integer(unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      any(z[friendindex[[m]][[s]] + ncs[m] + 1] <= 0)
    })})))
  lsel        <- lapply(1:S, function(m) drop[(ncs[m] + 1):ncs[m + 1]] != 1)
  
  # Max and Min of friend y and z
  yFmax       <- unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, max(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  yFmin       <- unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, min(y[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  zFmax       <- unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, max(z[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  zFmin       <- unlist(lapply(1:S, function(m){
    lapply(1:nvec[m], function(s){
      ifelse(ldg[[m]][s] > 0, min(z[friendindex[[m]][[s]] + ncs[m] + 1]), NA)
    })
  }))
  
  # In selection variables
  ldg         <- lapply(1:S, function(m) ldg[[m]][lsel[[m]]])
  lIs         <- lapply(1:S, function(m) lIs[[m]][lsel[[m]][lIs[[m]] - ncs[m] + 1]])
  lnIs        <- lapply(1:S, function(m) lnIs[[m]][lsel[[m]][lnIs[[m]] - ncs[m] + 1]])
  Is          <- unlist(lIs)
  nIs         <- unlist(lnIs)
  
  # In selection variables if empty groups are removed
  keepg       <- sapply(1:S, function(m) length(ldg[[m]]) > 0)
  ldg         <- ldg[keepg]
  S           <- length(ldg)
  SIs         <- sum(sapply(lIs, function(s) length(s) > 0))
  SnIs        <- sum(sapply(lnIs, function(s) length(s) > 0))
  
  list(friendindex = do.call(c, friendindex), frzeroy = frzeroy, frzeroz = frzeroz, 
       S = S, SIs = SIs, SnIs = SnIs, ldg = ldg, dg = unlist(ldg), lIs = lIs, 
       Is = Is, lnIs = lnIs, nIs = nIs, hasIso = (length(Is) > 0), yFmax = yFmax, 
       yFmin = yFmin, zFmax = zFmax, zFmin = zFmin)
}

fassignnull <- function(x, val){
  if (is.null(x)) {
    return(val)
  } else {
    return(x)
  }
} 

