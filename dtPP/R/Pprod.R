## =============================================================================
## Depth-integrated photosynthesis function
## =============================================================================

integratedPP <- function(
    z0=0, zn=100,             # max water depth to take into account
    nbox  = 100,
    zi = NULL, 
    
    times = 0:10,  # time for which PP needs to be calculated
    convFac = 1.,  # conversion to other units
    
    # photosynthesis-irradiance function and parameters
    PI.par = c(alpha=0.05, eopt=1000, ps=50),
    PI.fun = NULL,   # default PS function: Eilers-Peeters
    
    # Light as a function of time (t)
    It.data = data.frame(time=0:10, I=(0:10)*10),
    It.fun = NULL,
    
    # light as a function of water depth (z)
    kz   = 0.1,  # extinction coefficient (1 or 3 values)
    Iz.data = NULL,
    Iz.fun = NULL, # light as function of depth: function(z, ...) exp(-z*kz)
    
    # morphology (used for depth integration) - NOT VARIABLE IN TIME
    Mz.data = 1, #data.frame(Depth=c(0, zn), Surf=c(1,1))
    Mz.fun = NULL,
    
    Ht.data=NULL,
    Ht.fun=NULL,
    
    avgOver = NULL, avgTime = 1
    )                    
  
{  
  variableH <- FALSE  # water depth is constant by default
  
  # water depth
  
  if (! is.null(Ht.data) | ! is.null(Ht.fun)){ # water depth is variable
    
    if (!is.null(Ht.data)){
      Ht.fun <- approxfun(x=Ht.data, rule=2, ties=mean)
    } else if (!is.function(Ht.fun)) 
      stop ("'Ht.data' or 'Ht.fun' should be present")
    
    Ht.data   <- Ht.fun(times)
    variableH <- TRUE
    Hmax      <- max(Ht.data)  # maximum elevation
  } else Ht.data <- 0
  Hmax <- max(Ht.data)
  
  # --------------------------
  # Calculation grid
  # --------------------------
  
  if (is.null(zi)) {
    L    <- (zn+Hmax-z0)
    dx.1 <- L/nbox/10
  
    setup.grid <- function(dx.1){
      f.root <- function(p, dx.1, nbox, L) dx.1 * (p^(nbox)-1)/(p-1)- L
      p.est <- uniroot(f=f.root, dx=dx.1, nbox=nbox, L=L, 
                       lower=1.001, upper=10, tol=1e-20)$root
      DX <- dx.1*p.est^(0:(nbox-1))
      z0+c(0, cumsum(DX))
    }
    zi <- setup.grid(dx.1)
  } else zn <- max(zi)
  
  # thickness of the layers (required for integrating)
  dz <- diff(zi)
  
  # depth in center of the layers (required for output)
  z <- 0.5*(zi[-1] + zi[-length(zi)])  # depth in middle of box
  nbox  <- length(z)
  
  Dt.data <- zn+Ht.data  # height of water over all times
  
  # --------------------------
  # check input - TO DO
  # --------------------------
  
  # zn is positive?
  # M.fun OK? (i.e. estimate M.fun for range in z - should not have NAs)  
  
  # times?
  # It.fun OK? (i.e. estimate It.fun for range in times - should not have NAs)  
  
  # photosynthesis parameters PI.par
  # relative light profile for all z layers - stays constant
  
  # ----------------------------------
  # Light as function of time
  # ----------------------------------
  
  if (!is.null(It.data)){
    
    if (length(It.data) == 1)
      It.fun <- function(t) It.data
    else if (ncol(It.data) == 2)
      It.fun <- approxfun(x=It.data, rule=2, ties=mean)
    else
      stop ("'It.data' should be a data.frame or matrix with 2 columns, or one value")
  } else if (is.null(It.fun)) 
    stop ("'It.data' or 'It.fun' should be present")
  
  It.data <- It.fun(times)

  # ----------------------------------
  # vertical light extinction profile
  # ----------------------------------
  Ibot    <- 0
  # No light profile, but extinction coefficient(s)
  KZ <- NULL
  if (is.null(Iz.data) & ! is.null(kz)){
    if (length(kz) == 1){
        Iz.data <- exp(-kz*z)  # it is a vector
        Ibot <- It.data*exp(-kz*zn)
        KZ <- kz
    } else if (length(kz) == 3){
        Iz.data <- kz[1]*exp(-kz[2]*z)+(1-kz[1])*exp(-kz[3]*z)  
        Ibot <- It.data*(kz[1]*exp(-kz[2]*zn)+(1-kz[1])*exp(-kz[3]*zn))  
        KZ <- kz[1]*kz[2] + (1-kz[1])*kz[3]
      
    } else if (ncol(kz) == 2){  # time, kz
        kz <- approx(kz, xout=times, rule=2, ties=mean)$y
        # Iz.data is a matrix
        Iz.data <- outer(kz, z, FUN=function(x,y) exp(-x*y))
        Ibot    <- It.data*exp(-kz*zn) 
        KZ <- kz
        
    } else if (ncol(kz) == 4){  # time, p, k1, k2 (2 fractions)
        KK <- data.frame(
          p  = approx(kz[,  1:2 ], xout=times, rule=2, ties=mean)$y,
          k1 = approx(kz[,c(1,3)], xout=times, rule=2, ties=mean)$y,
          k2 = approx(kz[,c(1,4)], xout=times, rule=2, ties=mean)$y
        )
        KZ <- KK$p*KK$k1 + (1-KK$p)*KK$k2
        # Iz.data is a matrix
        Iz.data <- outer(1:nrow(KK), z, FUN=function(x, y) 
              KK$p[x]*exp(-KK$k1[x]*y) + (1-KK$p[x])*exp(-KK$k2[x]*y))
        Ibot <- It.data*(KK[,1]*exp(-KK[,2]*zn)+
                        (1-KK[,1])*exp(-KK[,3]*Dt.data))  
    } else
      stop ("'kz' not correct")
  } else {
    if (!is.null(Iz.data)){
      if (length(Iz.data) == 1)
         Iz.fun <- function(z) Iz.data
      else if (ncol(Iz.data) == 2)  # depth, z
         Iz.fun <- approxfun(x=Iz.data, rule=2, ties=mean)
      else
         stop ("'Iz.data' should be a data.frame or matrix with 2 columns, or one value")
    } else if (is.null(Iz.fun)) 
      stop ("'Iz.data' or 'Iz.fun' should be present")
  
    Iz.data <- Iz.fun(z)
    Ibot <- It.data*Iz.fun(zn)  
    
  }
  names(Ibot) <- NULL
  
  if (max(Iz.data) > 1)
    stop ("'Iz.data' is RELATIVE light per depth and should be at most 1")

  # ----------------------------------
  # morphology - volume for each layer
  # ----------------------------------
  
  if (!is.null(Mz.data)){
    if (length(Mz.data) == 1) 
      Mz.fun <- function(z) Mz.data
    else if (ncol(Mz.data) == 2)  # Mz.data is two columns
      Mz.fun <- approxfun(x=Mz.data, rule=2, ties=mean)
  } else if (is.null(Mz.fun)) 
    stop ("'Mz.data' or 'Mz.fun' should be present")
  
  dVol <- Mz.fun(z)*dz 
  meanSurf <- sum(dVol) / diff(range(zi))

  # ----------------------------------
  # PI parameters
  # ----------------------------------
  # if the parameters are a data.frame, the first column should be depth or time
  PI.type <- "depth"  # no need to transpose the data before calling PI.fun
  
  if (is.data.frame(PI.par) | is.matrix(PI.par)){
    if (nrow(PI.par) > 1) {
      nc     <- ncol(PI.par)
      zp     <- PI.par[,1]  # depth or time
      pnames <- colnames(PI.par)
      ppar   <- PI.par
      
      knownnames <- sum(c("depth", "time", "z", "t") %in% tolower(pnames))
      if (knownnames == 0 ) 
        stop("'PI.par' should have at least 't', and/or 'z' as column names")
      if (knownnames > 2 ) 
        stop("'PI.par' columnnames ambiguous: only one of 'z' or 'depth' or of 't' and 'time' should be present")
      
      if ( knownnames == 2){
        PI.type <- "zt"
        stop ( "method not yet written for variable time and space inputs of parameters")
      } else {
      
        if (tolower(pnames[1]) %in% c("depth", "z"))
          xout <- z
        else if (tolower(pnames[1] %in% c("time", "t")))  {
          xout     <- times
          PI.type <- "time"
        }
      
        # for all parameters:
        PI.par <- NULL
        for (i in 2:nc){
          PI.par <- cbind(PI.par, 
                          approx(zp, ppar[,i], xout=xout, rule=2, ties=mean)$y)
        }
        colnames(PI.par) <- pnames[-1]
        PI.par <- as.data.frame(PI.par)
      }
    } else   # only one value or a vector
      PI.par <- as.vector(PI.par)
  }
  
  if (is.null(PI.fun)) 
    PI.fun <- function(I) {
      ep <- I/((1/(PI.par[[1]]*PI.par[[2]]^2))*I^2+
                  (1/PI.par[[3]]-2/(PI.par[[1]]*PI.par[[2]]))*I+
                  (1/PI.par[[1]]))
      ep[is.na(ep)] <- 0
      
      return(ep)
  }

  # THIS IS RATHER SLOW
  # total photosynthesis, integrated over all layers (thickness dz) 
  # Function to estimate the production per layer; 
  #  Ptime <- function(irr)
  #    sum(PI.fun(Irel*irr)*dVol)
  # integrated photosynthesis for all light intensities 
  #  PP <- sapply(I, FUN=Ptime)
  
  # Memory intensive but faster - estimate TOTAL pp per layer and time
  
  # light at all depth intervals x times
  if (! is.matrix(Iz.data))                        # kz=1 value
    P_D <- outer(It.data, Iz.data, FUN="*")        
  else                                             # kz=timeseries
    P_D <- sweep(Iz.data, MARGIN=1, STATS = It.data,  FUN="*")  
  
  if (PI.type == "time")
    RR  <- PI.fun(P_D)*convFac                     # PS at these t x D  
  else if (PI.type == "depth")
    RR  <- t(PI.fun(t(P_D)))*convFac               # PS at these t x D  
  
  if (any (dim(RR)-dim(P_D) !=0))
    stop ("PI.fun(I) should return a matrix of same dimensions as light matrix")
  
  # correct for changing water depth  
  if (variableH) {
    Z_water <- outer(Dt.data, z, FUN=function(x,y)y < x)
    P_D     <- P_D*Z_water  
    RR      <- RR*Z_water  
  } 
  PP  <- sweep(RR, MARGIN=2, STATS=dVol, FUN="*")  # multiply depth with Volume

  # timeseries
  ts <- data.frame(times=times, PP=rowSums(PP)/meanSurf, TotalDepth=Dt.data, 
                    Isurf=It.data, Ibot=Ibot)  # PP per m2
  if (! is.null(KZ)) ts$kz <- KZ
  if (is.vector(PI.par))
    PI.par <- matrix(nrow=length(z), ncol=length(PI.par), byrow=TRUE, data=PI.par,
                     dimnames=list(NULL,names(PI.par)))
  
  # profile
  if (is.matrix(Iz.data)) 
    Iz.data <- apply(Iz.data, MARGIN=2, FUN=mean, na.rm=TRUE)
  
  PP_v    <- colMeans(RR)
  profile <- data.frame(z=z, PP_v=PP_v, Iz_I0=Iz.data, surface=Mz.data)
  
  if (PI.type == "depth") 
    profile <- data.frame(profile, PI.par)
  else
    ts <- data.frame(ts, PI.par)
  
  profile$dz <- dz
  
  if (! is.null(avgOver)){
    if (inherits(what = "POSIXt", x=ts$times) | inherits(what = "POSIXct", x=ts$times)){
      ts <- average_timeseries(ts, avgOver=avgOver, avgTime=avgTime, datetime="times",
                               value=colnames(ts)[colnames(ts)!="times"])
      LL <- data.frame(times=times, P_D)
      P_D <- as.matrix(average_timeseries(LL, avgOver=avgOver, avgTime=avgTime, datetime="times",
                                value=colnames(LL)[colnames(LL)!="times"])[,-1])
      LL <- data.frame(times=times, RR)
      RR <- as.matrix(average_timeseries(LL, avgOver=avgOver, avgTime=avgTime, datetime="times",
                                          value=colnames(LL)[colnames(LL)!="times"])[,-1])
    }
  }
  
  RES <- list(ts=ts, profile=profile, light=P_D, prod=RR)
  
  attr(RES, "Description") <- data.frame(
    variable   =c("ts.PP", "ts.TotalDepth",  "ts.Isurf", "ts.Ibot", 
               "profile.z", "profile.PP_v", "profile.Iz_I0", "profile.surface",
               "light", "prod"),
    description=c("integrated primary production", 
                  "Total water depth over which PP is estimated",
                  "Light at surface (I0)", 
                  "Light at the lowest grid point", 
                  "mean depth of the box at which photosynthesis is estimated", 
                  "Volumetric photosynthesis rate at depth z",
                  "fraction of light intensity at depth z (Iz/I0)", 
                  "surface area at middle of each box",
                  "light intensity for all times and depths",
                  "volumetric photosynthesis for all times and depths"),
    units     = c("mass/surface/time", "length, e.g. m", 
                  "Light intensity, as used in PI data, e.g. uEinst/m2/s",
                  "Light intensity",
                  "length unit, e.g. m",
                  "mass/volume/time", "(fraction)", "length^2",
                  "Light intensity", "mass/volume/time"))
  class(RES) <- c("integratedPP", class(RES))
  return(RES)
}

## =============================================================================
## General plotting functions for depth-integrated photosynthesis
## =============================================================================

plot.integratedPP <- function(x, mass="ugC", length="m", time="s", 
                              light="uEinst/m2/s", type="l", 
                              las=1, which="ts", ...){
  # for labelling the plots
  surface <- paste( length, "^2", sep="")
  vol     <- paste( length, "^3", sep="")
  
  pplab <- paste(mass, surface, time, sep="/")
  ppvol <- paste(mass, vol, time, sep="/")
  
  # number of figures to plot
  
  ldots <- list(...)
  
  nv <- 0  # columns to plot
  
  if (length(which) == 0)
    stop("'which' cannot be NULL; should be 'profile', 'ts' or a selection of names in these")
  
  which.ts <- which.prof <- NULL
  
  colname.ts <- colnames(x$ts)
  if ("ts" %in% which)   
    which.ts <- colname.ts[-1]  # first one is time
  else
    which.ts <- colname.ts[which(colname.ts %in% which)]
  
  colname.prof <- colnames(x$profile)
  if ("profile" %in% which) 
    which.prof <- colname.prof[-c(1, ncol(x$profile))]
  else
    which.prof <- colname.prof[which(colname.prof %in% which)]

  nv <- length(which.ts) + length(which.prof)
  
  if (is.null(ldots$mfrow)) {
    nc <- min(ceiling(sqrt(nv)), 3)
    nr <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
    pm <- par(mfrow=mfrow)
  }  else pm <- par(mfrow=ldots$mfrow)

  xlab <- "times"
  units <- c(PP  = pplab, TotalDepth=length, Isurf=light, Ibot=light, 
             z = length, PP.v=ppvol, Iz_I0="-", surface=surface)
  
  unitunknown <- c(which.ts[!which.ts %in% names(units) ],
                   which.prof[!which.prof %in% names(units) ])
  newunit <- rep("-", times=length(unitunknown))
  names(newunit) <- unitunknown
  units <- c(units, newunit)
  
  main <- c(PP  = "Integrated production", TotalDepth="Total depth", 
            Isurf= "Light at surface", Ibot="Light at bottom", 
             z = "depth", PP.v="mean volumetric production",
            Iz_I0="mean light penetration (relative)", 
            surface="horizontal surface area")
  mainunknown <- c(which.ts[!which.ts %in% names(main) ],
                   which.prof[!which.prof %in% names(main) ])
  newmain <- mainunknown
  names(newmain) <- newmain
  main <- c(main, newmain)
  
  if (length(which.ts)){

    for (var in which.ts)
     plot(x$ts$times, x$ts[,var], xlab=xlab, ylab=units[var], 
        main=main[var],type=type, las=las,...)
  }
  
  if (length(which.prof)){
      
    for (var in which.prof)
      plot (x$profile[,var], x$profile$z, ylab= paste("depth", length), xlab=units[var], 
                         main=main[var], 
                         ylim=rev(range(x$profile$z)), type=type,  las=las, ...)
  }
  
  par(mfrow=pm)
}

## =============================================================================

image2D.integratedPP <- function(z, mass="ugC", length="m", time="s", 
                              light="uEinst/m2/s",  las=1, 
                              ...){
  if (is.null(list(...)$mfrow)) {
    mfrow <- c(1,2)
    pm <- par(mfrow=mfrow)
  }
  depth <- z$profile$z
  times <- z$ts$times
  xlab  <- "times"
  ylab  <- paste("depth,", length)
  ylim  <- rev(range(depth))  
  vol   <- paste( length, "^3", sep="")
  
  image2D(x=times, y=depth, z=z$light, ylim=ylim, ylab=ylab, 
          clab=light, main="light", las=las, ...)
  image2D(x=times, y=depth, z=z$prod, ylim=ylim, ylab=ylab, 
          clab=paste(mass, vol, time, sep="/"),
          main="volumetric production", las=las, ...)
  if (! is.null(mfrow)) 
    par(mfrow=pm)
}
