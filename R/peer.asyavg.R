#' @title Computing peer (asymetric) average values
#' @export
peer.asyavg <- function(formula, Glist, data, nthread = 1) {
  
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  
  formula    <- as.formula(formula)
  hasy       <- (length(as.formula(formula)) == 3)
  f.t.data   <- formula.to.data(formula = formula, data = data, fixed.effects = TRUE,
                                simulations = !hasy) #Intercept is not necessary
  y          <- f.t.data$y
  X          <- f.t.data$X
  Kx         <- ncol(X)
  xname      <- f.t.data$xname
  yname      <- f.t.data$yname
  
  S          <- length(Glist)
  nvec       <- unlist(lapply(Glist, ncol))
  
  if (sum(nvec) != nrow(X)) {
    stop("Glist and V do not match")
  }
  cumsn      <- c(0, cumsum(nvec))
  tp         <- NULL
  if (hasy) {
    tp <- highlowstat2(y = y, X = X, G = Glist, cumsn = cumsn, nvec = nvec, 
                       ngroup = S, nthread = nthread)
  } else {
    tp <- highlowstat1(X = X, G = Glist, cumsn = cumsn, nvec = nvec, 
                       ngroup = S, nthread = nthread)
  }
  g     <- NULL
  pavg  <- NULL
  hg    <- NULL
  hpavg <- NULL
  lg    <- NULL
  lpavg <- NULL
  if (hasy) {
    g     <- c(tp$g)
    hg    <- c(tp$gh)
    lg    <- c(tp$gl)
    pavg  <- cbind(tp$ybar, tp$Xbar)
    hpavg <- cbind(tp$ybh, tp$Xbh)
    lpavg <- cbind(tp$ybl, tp$Xbl)
  } else {
    if (Kx == 1) {
      g     <- c(tp$g)
      hg    <- c(tp$gh)
      lg    <- c(tp$gl)
      pavg  <- c(tp$Xbar)
      hpavg <- c(tp$Xbh)
      lpavg <- c(tp$Xbl)
    } else {
      g     <- tp$g
      hg    <- tp$gh
      lg    <- tp$gl
      pavg  <- tp$Xbar
      hpavg <- tp$Xbh
      lpavg <- tp$Xbl
    }
  }
  
  list(degree     = g,
       peer.avg   = pavg,
       hdegree    = hg,
       hpeer.avg  = hpavg,
       ldegree    = lg,
       lpeer.avg  = lpavg)
}