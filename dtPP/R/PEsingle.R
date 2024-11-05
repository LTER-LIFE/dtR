## =============================================================================
## Single Fits of FRRF data and the PE models
## =============================================================================

  #Define PE Models
  # Model.I refers to PE model normalized to irradiance, 
  # useful for fluorescence measures
  # note: double quotes so that it also works with data.frames of parameters
  # (as long as they are as long as I)
  # Webb (1974); 

  fWebb   <- function(I, p = c(alpha=0.2, ek=200)) 
    return(data.frame(x = I, 
                      y = p[[1]]*p[[2]]*(1-exp(-1*I/p[[2]]))))  
  
  fWebb.I <- function(I, p = c(alpha=0.2, ek=200)) 
    return(data.frame(x = I, 
                      y = p[[1]]*p[[2]]*(1-exp(-1*I/p[[2]]))/I))
  
  # Jassby and Platt (1976)
  fJP     <- function(I, p = c(alpha=0.15, ek=250))
    return(data.frame(x = I, 
                      y = p[[1]]*p[[2]]*tanh(I/p[[2]])))  
  
  fJP.I   <- function(I, p = c(alpha=0.15, ek=250)) 
    return(data.frame(x = I, 
                      y = p[[1]]*p[[2]]*tanh(I/p[[2]])/I))
  
  #  Platt, Galegos and Harrison (1980)
  fPG     <- function(I, p = c(alpha=0.2, beta=0, ps=40)){
    pg <- p[[3]]*(1-exp(-1*p[[1]]*I/p[[3]]))*exp(-1*p[[2]]*I/p[[3]])
    pg[is.nan(pg)] <- 0
    return(data.frame(x = I, 
                      y = pg))
    }
  
  fPG.I   <- function(I, p = c(alpha=0.2, beta=0, ps=40)) 
    return(data.frame(x = I, 
                      y = p[[3]]*(1-exp(-1*p[[1]]*I/p[[3]]))*exp(-1*p[[2]]*I/p[[3]])/I))
  
  # Eilers and Peeters (1988);
  fEP     <- function(I, p = c(alpha=0.2, eopt=830, ps=40)){
    ep <- I/((1/(p[[1]]*p[[2]]^2))*I^2+(1/p[[3]]-2/(p[[1]]*p[[2]]))*I+(1/p[[1]]))
    ep[is.nan(ep)] <- 0
    return(data.frame(x = I, 
                      y = ep))
    }
  
  fEP.I   <- function(I, p = c(alpha=0.2, eopt=830, ps=40)) 
    return(data.frame(x = I, 
                      y =   1/((1/(p[[1]]*p[[2]]^2))*I^2+(1/p[[3]]-2/(p[[1]]*p[[2]]))*I+(1/p[[1]]))))
  
## -----------------------------------------------------------
## Generate predictions based on an FRRF fit.
## -----------------------------------------------------------

predict.fitPI <- function(object, I, modargs = NULL, normalized = FALSE, ...) {
  x     <- I
  p     <- object$par
  model <- object$model
  
  if (is.function(model) & is.null(modargs))
    Pred <- model(x, p)$y

  else if (is.function(model))
    Pred <- do.call("model", c(alist(x, p), modargs))$y
 
  else if (model == "Webb") # Webb et al. (1974) 
    Pred <- fWebb(x, p)$y
  
  else if (model == "JP")   # Jassby - Platt
    Pred <- fJP(x, p)$y
 
  else if (model == "PG") 
    Pred <- fPG(x, p)$y
   
  else if (model == "EP")
    Pred <- fEP(x, p)$y
  
  else
    stop ("model '", model, "' not found")

  if (normalized) 
    Pred <- Pred/x
  
  return(Pred)
}

## -----------------------------------------------------------
## Plot fitted model(s) versus the data 
## -----------------------------------------------------------

plot.fitPI <- function(x, modargs = NULL, normalized = FALSE, ...) {

   if (! inherits (x, c("fitmultiPI", "fitPI")))
     stop ("'x' not of correct class, should be created with fitPI or fitmultiPI")

   dots <- list(...)
   if (is.null(dots$xlab)) 
     dots$xlab <- "I"
   
   ylab <- dots$ylab
   if (is.null(dots$ylab))
     dots$ylab <- "response"
       
   Plotit <- function(x) {
     if (normalized){
       ii <- which(x$I != 0)
       do.call ("plot", c(alist(x = x$I[ii], y = x$response[ii]/x$I[ii]), dots))
     }     
     else      
       do.call ("plot", c(alist(x = x$I, y = x$response), dots))
    I <- seq(min(x$I), max(x$I), length.out = 100) 
    pp  <- predict (x, I = I, modargs = modargs)
    lines(I, pp)
   }
   
   if ( inherits (x, "FRRF"))
    Plotit(x)
   
   else {
    nv <- length(x$ssr)
    nc <- min(ceiling(sqrt(nv)), 3)
    nr <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
    mf <- par(mfrow = mfrow)
    Ask <- prod(par("mfrow")) < nv && dev.interactive()
    ask <- par(ask = Ask)
    A <- attributes(x)$modFit
    obnames <- names(x$ids)
    for (i in 1:nv) {
      dots$main <- obnames[i]
      Plotit(A[[i]])
    }
    par(ask = ask) 
  }
}

## -----------------------------------------------------------
## Estimate r squared (uncorrected) of a fit
## -----------------------------------------------------------

r.squared <- function (object) {      # mss/(mss + rss)
  
  rsq.single <- function(x) {
    pred <- predict(x, I = x$I)
    mss <- sum((pred - mean(pred))^2)
    rss <- sum(x$residuals^2)
    mss/(mss + rss)
  }
  
  if (inherits(object , "FRRF"))
    return(rsq.single(object))
  else {
    vec <- vector()
    for (i in 1:nrow(object$par))
      vec[i] <- rsq.single(attributes(object)$modFit[[i]])
    return(vec)
  }  
}

## ======================================================
## ======================================================
## Find best-fit parameters
## ======================================================
## ======================================================
  
    
fitPI <- function(model, I, response, normalized = FALSE, pini = NULL, 
  modargs = NULL, lower = 0, upper = Inf, 
  checkoutlier = FALSE, min.rsq = 0.995, ...){
  
  fit <- fitOne(model = model, I = I, response = response, 
                normalized = normalized, pini = pini, 
                modargs = modargs, lower = lower, upper = upper, ...)
  
  fit$removed <- NA
  
  if (checkoutlier) {  # see if removing an outlier improves the fit
    rsq  <- r.squared(fit)
    
    if (rsq < min.rsq){  # only if the fit was not very good to start with
     N    <- length(response)
     AIC0 <- AIC(fit)  # akaike information criterion
     pini <- fit$par
    
     for (i in 1:N) {
      FF <- fitOne(model = model, I = I[-i], response = response[-i], 
                   normalized = normalized, pini = pini, 
                   modargs = modargs, lower = lower, upper = upper, ...)
      AICi <- AIC(FF)
      if (AICi < AIC0){
        AIC0 <- AICi
        fit  <- FF
        fit$removed <- c(I = I[i], response = response[i])
      }
     }
    }
  }
  return(fit)
}

fitOne <- function(model, I, response, normalized = FALSE, pini = NULL, 
                    modargs = NULL, lower = 0, upper = Inf, ...){
    
    x <- I
    y <- response
    
    if (normalized)  {
      ii <- x > 0
      y  <- y/x
    }      
    
    fit <- list()
    
    #Initial Parameter Estimates  - in case pini = NULL
    if (is.null(pini)) {
      beta  <- 0.05
      eopt  <- max(I, na.rm=TRUE)
      ek    <- eopt/2 # mean(range(x))
      ps    <- max(response, na.rm=TRUE)
      alpha <- ps/ek   # max(y)
    }
    
    if (is.function(model)) {
      if (!normalized) 
        model.1 <- function(p) (y - do.call("model", c(alist(x, p), modargs))$y)
      else
        model.1 <- function(p) (y[ii] - do.call("model", c(alist(x[ii], p), modargs))$y/x[ii])
      
      if (is.null(pini))
        stop("'pini' should be specified if the 'model' is a function")
      FIT <- try(modFit(f = model.1, p = pini,
                        ...), silent = TRUE)
      if (inherits(FIT, "try-error"))
        fit <- list(ms = NA, ssr = NA, 
                    residuals = rep(NA, length(x)), par = NA)
      else  {
        fit[1:length(pini)] <- FIT$par
        names(fit)[1:length(pini)] <- names(pini)
        fit <- c(fit, FIT)
      }  
    } else if (model == "Webb"){#Call Webb et al. (1974) Model Normalized to E
      
      if (!normalized) 
        model.1 <- function(p) (y - fWebb(x, p)$y)
      else 
        model.1 <- function(p) (y[ii] - fWebb.I(x[ii], p)$y)
      
      if (is.null(pini)) {
        pini <- c(alpha = alpha, ek = ek)
      }  
      
      upper <- rep(upper, length.out = 2)
      lower <- rep(lower, length.out = 2)
      FIT <- try(modFit(f = model.1, p = pini,
                        lower = lower, upper =upper, # c(Inf, 2500), 
                        ...), silent = TRUE)
      if (inherits(FIT, "try-error"))
        fit <- list(alpha = NA, ek = NA, ms = NA, ssr = NA, 
                    residuals = rep(NA, length(x)), par = NA)
      else
        fit <- c(fit, alpha = FIT$par[1], ek = FIT$par[2],FIT)
    }
    
    else if (model == "JP"){
      
      if (!normalized) 
        model.1 <- function(p) (y - fJP(x, p)$y)
      else 
        model.1 <- function(p) (y[ii] - fJP.I(x[ii], p)$y)
      
      if (is.null(pini))
        pini <- c(alpha = alpha, ek = ek)
      
      upper <- rep(upper, length.out = 2)
      lower <- rep(lower, length.out = 2)
      
      FIT <- try(modFit(f = model.1, p = pini, 
                        lower = lower, upper = upper, #c(Inf, 2500), 
                        ...), silent = TRUE)
      
      if (inherits(FIT, "try-error"))
        fit <- list(alpha = NA, ek = NA, ms = NA, 
                    ssr = NA, residuals = rep(NA, length(x)), par = NA)
      else
        fit <- c(fit, alpha = FIT$par[1], ek = FIT$par[2], FIT)
    }
    
    else if (model == "PG"){
      
      if (!normalized) 
        model.1 <- function(p) (y - fPG(x, p)$y)
      else 
        model.1 <- function(p) (y[ii] - fPG.I(x[ii], p)$y)
      if (is.null(pini))
        pini <- c(alpha = alpha, beta = beta, ps = ps)
      
      upper <- rep(upper, length.out = 3)
      lower <- rep(lower, length.out = 3)
      
      FIT <- try(modFit(f = model.1, p = pini, 
                        lower = lower, upper = upper, #c(Inf, 10, 1000),
                        ...), silent = TRUE)
      if (inherits(FIT, "try-error"))
        fit <- list(alpha = NA, beta = NA, ps = NA, 
                    ms = NA, ssr = NA, residuals = rep(NA, length(x)), par = NA)
      else
        fit <- c(fit, alpha = FIT$par[1], beta = FIT$par[2], ps = FIT$par[3], FIT)
    }
    
    else if (model == "EP"){
      
      if (!normalized) 
        model.1 <- function(p) (y - fEP(x, p)$y)
      else 
        model.1 <- function(p) (y[ii] - fEP.I(x[ii], p)$y)
      
      if (is.null(pini))
        pini <- c(alpha = alpha, eopt = eopt, ps = ps)
      
      upper <- rep(upper, length.out = 3)
      lower <- rep(lower, length.out = 3)
      
      FIT <- try(modFit(f = model.1, p = pini, 
                        lower = lower, upper = upper, #upper = c(Inf, 2500, 1000), 
                        ...), silent = TRUE)
      if (inherits(FIT, "try-error"))
        fit <- list(alpha = NA, eopt = NA, ps = NA, 
                    ms = NA, ssr = NA, residuals = rep(NA, length(x)), par = NA)
      else
        fit <- c(fit, alpha = FIT$par[[1]], eopt = FIT$par[[2]], ps = FIT$par[[3]], FIT)    
    } else
      stop ("model '", model, "' not found")
    
    fit$I          <- I
    fit$response   <- response
    fit$model      <- model
    fit$normalized <- normalized
    
    class(fit)     <- c("fitPI", "FRRF", "modFit")
    return(fit)
    
  }
  
