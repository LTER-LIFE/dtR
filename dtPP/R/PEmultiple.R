
## -----------------------------------------------------------
## Get data, 3 columns on which to fit a model - NOT USED
## -----------------------------------------------------------

getPIdata <- function(data, which = NULL, id = NULL, I = NULL, minPAR = 0.1, 
  environment = NULL) {
  P <- ifelse(is.null(I), "I", I)
  if (! P %in% colnames(data))
    stop("column with I data, ", P, " not found - supply name of column with I data" )
  if (is.null(I))
    I <- pmax(minPAR, data$I)
  else if (is.character(I))
    I <- pmax(minPAR, data[,I])
  
    
  T <- LAT <- LON <- NA
  if (is.null(id)) {
    ii <- c(0, which(I[-1]/I[-nrow(data)] < 0.1), nrow(data))
    id <- as.vector(unlist(sapply(2:length(ii), FUN = function(i) rep(i-1, ii[i] - ii[i-1]))))
    if (! is.null(data$Date) & ! is.null(data$Time)) {
      # this is all too complex
      TT <- as.POSIXlt(paste(data$Date, data$Time), format = "%d/%m/%Y %T")
      T1 <- TT[ii[-length(ii)]+1]
      T2 <- c(T1[-1], TT[nrow(data)])
      
      T <- lapply(2:(length(ii)), FUN = function(i) mean(TT[(ii[i-1]+1) : ii[i]]))
      Tdiff <- (as.POSIXlt(T2) - as.POSIXlt(T1))
      T <- as.vector(t(as.data.frame(T)))
    } else   if (! is.null(data$Date)) {
      # this is all too complex
      TT <- as.POSIXlt(data$Date, format = "%d/%m/%Y")
      T1 <- TT[ii[-length(ii)]+1]
      T2 <- c(T1[-1], TT[nrow(data)])

      T <- lapply(2:(length(ii)), FUN = function(i) mean(TT[(ii[i-1]+1) : ii[i]]))
      Tdiff <- (as.POSIXlt(T2) - as.POSIXlt(T1))
      T <- as.vector(t(as.data.frame(T)))
    } else {
      T1 <- T2 <- T <- Tdiff <- NA
    }

    envir <- data.frame(ini_time = T1, end_time = T2, mean_time = T, duration_min = Tdiff)
    cn <- names(data)
    for (e in environment) {
      if (e %in% cn) 
        envir[,e] <- sapply(2:(length(ii)-1), FUN = function(i) mean(data[(ii[i-1]+1) : ii[i], e]))
      else
        warning("environment variable ", e, " not found")  
    }
  }
  if (is.null(which))
    Dat <- data.frame(id = id, I = I, data)
  
  else  {
    Dat <- data.frame(id = id, I = I, data[, which])
    Dat <- Dat[! is.na(Dat[,3]), ]
    names(Dat)[-(1:2)] <- which
  }
   
  attr (Dat, "environment") <- envir
  return(Dat)
}  

coef.fitPI <- function(object, ...) {
  if (inherits(object, "FRRF")) return (summary(object))
  object$par
}

coef.fitmultiPI <- function(object, ...) {
  if (inherits(object, "FRRF"))  return (summary(object))
  object$par
}

predict.fitmultiPI <- function (object, I, ...) {
  y <- NULL
  for (i in 1:nrow(object$par)) {
    FIT <- object
    class(FIT) <- c("fitPI", "modFit")
    FIT$par <- object$par[i,]
    y <- cbind(y, predict(FIT, I, ...))
  } 
  y   
}

summary.fitmultiPI <- function (object, ...) {
  n <- ncol(object$par)
  m <- nrow(object$par)
  ids <- names(object$ids)
  PP <- NULL
  sigma <- vector(length = m)
  df <- vector(length = m)
  residualVariance <- vector(length = m)
  modVariance <- vector(length = m)
  niter <- vector(length = m)
  rsq <- r.squared(object)
  for (i in 1:nrow(object$par)) {
    FIT <- attributes(object)$modFit[[i]]
    if (! is.na(FIT$ssr)) {
      SM <- summary(FIT)
      PP <- rbind(PP, data.frame(ids = rep(ids[i], n), SM$par))
      sigma[i] <- SM$sigma
      df[i] <- SM$df[2]
      residualVariance[i] <- SM$residualVariance
      modVariance[i] <- SM$modVariance
      niter[i] <- SM$niter
    } else {
      SMp <- matrix(nrow = n, ncol = 4, data = NA)
      rownames(SMp) <- colnames(object$par)
      colnames(SMp) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)") 
      PP <- rbind(PP, data.frame(ids = rep(ids[i], n), SMp))
    }    
  } 
  ans <- list(residuals = object$residuals, residualVariance = residualVariance, 
        sigma = sigma, modVariance = modVariance, df = df, rsq = rsq, 
        niter = niter, par = PP, ids = object$ids)
  class(ans) <- "summary.fitmultiPI"      
  ans 
}

print.summary.fitmultiPI <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{

    cat("\nParameters:\n")
    printCoefmat(x$par, digits = digits, ...)

    cat("\nResidual standard error and degrees of freedom:\n")
    stderrdf <- as.data.frame(rbind(sigma = x$sigma, df = x$df, rsq = x$rsq ))
    names(stderrdf) <- names(x$ids)
    printCoefmat(stderrdf, digits = digits, ...)

#    cat("\nResidual standard error:", format(signif(x$sigma, 
#        digits)))
#    cat("\nDegrees of freedom:", df, "\n")
    invisible(x)
}


fitmultiPI <- function(model, id, I, response, data = NULL, 
  normalized = FALSE, pini = NULL, 
 ...) {

  if (is.character(id))
    if (length(id) == 1) id <- data[,id]
  if (is.character(I))
    I <- data[,I]
  if (is.character(response))
    response <- data[,response]
    
  iun       <- unique(id)
  n         <- length(iun)
  ssr       <- rep(NA, n)
  ms        <- rep(NA, n)
  residuals <- rep(NA, length(id))

  pars <- NULL       
  ids  <- table(id)
  mFIT <- list()
  
  for (i in 1:n){
    
    #Identify & isolate data
    pos   <- which(id == iun[i])
    myfit <- fitPI(model, I[pos], response[pos], normalized, pini = pini, ...) 
    pars  <- rbind(pars, myfit$par)    

    ssr[i]         <- myfit$ssr
    ms[i]          <- myfit$ms
    residuals[pos] <- myfit$residuals
    
    mFIT[[i]] <- myfit
    
  }
    
   FIT <- list(par = pars, ssr = ssr, ms = ms, residuals = residuals, ids = ids)
   FIT$model <- model
   FIT$normalized <- normalized
   class(FIT) <- c("fitmultiPI", "fitPI", "modFit")
   attr(FIT, "modFit") <- mFIT 
   if (! is.null(data))
     attr(FIT, "environment") <- attributes(data)$environment 
   return(FIT)  
}
