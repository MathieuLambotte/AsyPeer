#' @importFrom stats lm
#' @importFrom stats glm
#' @importFrom stats runif
#' @importFrom randomForest randomForest
#' 
mpredict  <- function(ddy, ddX, id_fold, model, nthread, ...){
  #Given a vector of fold id, create a list of the corresponding row of each ddyad
  # belonging in each fold
  id_list <- split(seq_along(id_fold), id_fold)
  seed    <- as.integer(runif(1, 0, 1e9))
  
  cl      <- makeCluster(nthread)
  registerDoParallel(cl)
  registerDoRNG(seed)
  lrho    <- foreach(k         = id_list, 
                     .export   = "mpredict_fold",
                     .packages = c("randomForest") #Remember to add "NameOfThePackage"
  ) %dorng% {
    #each observation in fold k is predicted using a model trained
    #on the observations of the other folds
    ARG <- list(ddX = ddX, ddy = ddy, id_listk = k, model = model, ...)
    do.call(mpredict_fold, ARG) 
  }
  stopCluster(cl)
  
  rho   <- numeric(nrow(ddX))
  for (k in 1:length(id_list)) {
    rho[id_list[[k]]] <- lrho[[k]]
  }

  return(rho)
}


mpredict_fold <-function(ddX, ddy, id_listk, model, ...){
  #gather the observations from the other folds, expect id_listk
  ddX_train <- data.frame(ddX[-id_listk, ,drop = FALSE])
  ddy_train <- ddy[-id_listk]
  
  #gather the observations from the fold k
  ddX_k <- data.frame(ddX[id_listk, , drop = FALSE])
  if (model == "ols") {
    ARG         <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    model_train <- do.call(lm, ARG) 
    rho_k       <- predict(model_train, newdata = ddX_k)
  } else if (model == "glm") {
    ARG           <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    if (is.null(ARG$family)) { # If the user does not set link, use logit
      ARG$family  <- binomial(link = logit)
    }
    model_train <- do.call(glm, ARG)  
    rho_k       <- predict(model_train, newdata = ddX_k, type = "response")
  } else if (model=="RF") {
    ddy_train   <- as.factor(ddy_train)
    ARG         <- list(formula = ddy_train ~ ., data = ddX_train, ...)
    model_train <- do.call(randomForest, ARG)  
    rho_k       <- predict(model_train, newdata=ddX_k,type="prob")[,"1"]
  }
  return(rho_k)
}
