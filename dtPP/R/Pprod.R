## =============================================================================
## =============================================================================
## Depth-integrated photosynthesis function
## =============================================================================
## =============================================================================

integratedPP <- function(
    zmin  = 0,          # shallowest water depth
    zmax  = 100,        # max water depth to take into account
    nbox  = 100,        # number of grid cells
    zi = NULL,          # depths over which to estimate the PP
    
    times = 0:10,       # time for which PP needs to be calculated
    convFac = 1.,       # conversion to other units
    
    # photosynthesis-irradiance function and parameters
    PI.par = c(alpha=0.05, eopt=1000, ps=50),
    PI.fun = NULL,      # default PS function: Eilers-Peeters
    
    # Light as a function of time (t)
    It.data = data.frame(time=0:10, I=(0:10)*10),
    It.fun = NULL,
    
    # light as a function of water depth (z)
    kz   = 0.1,         # extinction coefficient (1 or 3 values)
    Iz.data = NULL,
    Iz.fun = NULL,      # light ~ depth function(z, ...) exp(-z*kz)
    
    # morphology (used for depth integration) - NOT VARIABLE IN TIME
    Mz.data = 1,        # data.frame(Depth=c(0, zmax), Surf=c(1,1))
    Mz.fun = NULL,
    
    Ht.data=NULL,
    Ht.fun=NULL,
    
    avgOver = NULL, avgTime = 1,
    
    verbose = TRUE,
    unit = list(mass = "mgC", length = "m", time = "h", light = "uEinst/m2/s")
    
    )                    
  
{  
  variableH <- FALSE  # water depth is constant by default
  
  classTime <- class(times)
  rTime     <- range(times)
  
  if (any(classTime %in% c("numeric", "integer"))) 
    classTime <- c("numeric", "integer")
  
  # ----------------------------------------------------------------------------
  # Functions to check the data
  # ----------------------------------------------------------------------------
  
  # check compatibility of "times" and data
  
  CheckTime <- function(Data, var){  
    
    # first column should be of same class as "times"
    if (! inherits(Data[,1], classTime))
      stop ("cannot use timeseries `", var, "`first column is not compatible with `times`")
    
    if (verbose)
      if (max(Data[,1], na.rm=TRUE) < rTime[1] |
          min(Data[,1], na.rm=TRUE) > rTime[2]) 
        warning("no overlap between output 'times' and dataseries ", var)
  }
  
  # --------------------------
  # check general data series structure -  return function
  # --------------------------
  
  DataSeries <- function(Data, var, indep = "time")  {  
    
    nc   <- ncol(Data)
    TFUN <- NULL
    err  <- paste ("`", var, 
                  ".data' should be a data.frame or matrix with 2 columns ('", 
                  indep, "', value), a vector of same length as",
                  indep, ", or one value", sep="")
    
    
    if (length(Data) == 1 & indep == "time")   # one value
      return(rep(Data, times=length(times)))
    
    else if (length(Data) == 1 & indep == "depth")
      return(rep(Data, times=nbox))
    
    else if (is.null(nc)) {   # a vector 
      if (indep == "time" & length(Data) == length(times) |
         (indep == "depth" & length(Data) == nbox)) 
        return(Data)
      
      else 
        stop (err)
    
    } else if (nc == 2 ){
      if (indep == "time") CheckTime (Data, var)
      TFUN  <- approxfun(x = Data, rule=2, ties=mean)
    }
    else stop (err)    
    
    if (indep == "time") 
      return(TFUN (times))
    
    else    # function of depth (z)
      return(TFUN (z))

  } 

  # --------------------------
  # data or function - return data
  # --------------------------
  
  getdata <- function(Data, fun, var, indep="time"){
    
    if (! is.null(Data))
      return(DataSeries(Data, var, indep))
    
    else if (is.null(fun)) 
      stop ("`" , var, ".data` or `", var, ".fun` should be present")
  
    if (indep == "time") 
      return(fun (times))
    
    else    # function of depth (z)
      return(fun (z))
  }

  
  # ----------------------------------------------------------------------------
  # water depth as a function of time
  # ----------------------------------------------------------------------------
  
  if (! is.null(Ht.data) | 
      ! is.null(Ht.fun)) { # water depth is variable
    
    Ht.data <- getdata(Ht.data, Ht.fun, "Ht", indep="time")
      
    variableH <- TRUE

  } else Ht.data <- 0
  
  Hmax <- max(Ht.data)
  
  # ----------------------------------------------------------------------------
  # Calculation grid
  # ----------------------------------------------------------------------------
  
  if (is.null(zi)) {
    
    L    <- zmax + Hmax - zmin
    if ( L < 0) L <- 0.01
    dz.1 <- L/nbox/10       # size of first box
  
    setup.grid <- function(dx.1){
      f.root <- function(p, dx.1, nbox, L) 
        dx.1 * (p^(nbox)-1)/(p-1)- L
      
      p.est <- uniroot(f = f.root, dx = dx.1, nbox=nbox, L=L, 
                       lower = 1.001, upper = 10, tol = 1e-20)$root
      DX <- dx.1*p.est^(0:(nbox-1))
      return(zmin + c(0, cumsum(DX)))
    }
    zi <- setup.grid(dz.1)
  
  } else {
    zmax <- max(zi)
    if (length(zi) < 2) stop ("'zi' should contain at least two numbers: uppen and lower boundary")
  }
  
  # thickness of the layers (required for integrating)
  dz    <- diff(zi)
  
  # depth in center of the layers (required for output)
  z     <- 0.5*(zi[-1] + zi[-length(zi)])  
  nbox  <- length(z)
  
  Dt.data <- pmax(0, zmax + Ht.data)  # height of water over all times
  # if (max(Dt.data) == 0 ){
  #
  #}
  
  # ----------------------------------------------------------------------------
  # check input - TO DO
  # ----------------------------------------------------------------------------
  
  # zmax is positive?
  # M.fun OK? (i.e. estimate M.fun for range in z - should not have NAs)  
  
  # times?
  # It.fun OK? (i.e. estimate It.fun for range in times - should not have NAs)  
  
  # photosynthesis parameters PI.par
  # relative light profile for all z layers - stays constant
  
  # ----------------------------------------------------------------------------
  # Light as function of time
  # ----------------------------------------------------------------------------
  
  It.data <- getdata(It.data, It.fun, "It", indep="time")

  # ----------------------------------------------------------------------------
  # vertical light profile (given or defined by extinction coeff)
  # ----------------------------------------------------------------------------
  
  Ibot    <- 0
  KZ      <- NULL
  
  # light profile or light function  

  if (! is.null(Iz.data) | ! is.null(Iz.fun)) {
    Iz.data <- getdata(Iz.data, Iz.fun, "Iz", indep="depth")
  
    Ibot    <- It.data * Iz.data[length(Iz.data)]  # light in last box
  
  } else if (! is.null(kz)) {
    
    # No light profile, but extinction coefficient(s)
    
    ncol_kz <- ncol(kz)
    
    if (length(kz) == 1){
        
        Iz.data <- exp(-kz*z) 
        Ibot    <- It.data*exp(-kz*Dt.data)   # light at the bottom
        KZ      <- kz                         # for output
    
    } else if (length(kz) == 3){      # two extinction coefficients
        
        Iz.data <-    kz[1]  * exp(-kz[2]*z) + 
                   (1-kz[1]) * exp(-kz[3]*z)  
        
        Ibot    <- It.data* (kz[1]     * exp(-kz[2]*Dt.data) +
                             (1-kz[1]) * exp(-kz[3]*Dt.data))  
        
        # average exctinction coeff
        KZ      <- kz[1]*kz[2] + (1-kz[1])*kz[3]   
      
    } else if (is.null(ncol_kz)) {  # kz already matches the timeseries
       if (length(kz) != length(times))
         stop ("'kz' not correct")
       Kz <- kz
       # Iz.data is a matrix
       Iz.data <- outer(kz, z, 
                        FUN=function(x,y) exp(-x*y))
       
       Ibot    <- It.data * exp(-kz*Dt.data) 
       
    } else if (ncol_kz == 2) {  # time, kz values
        
        CheckTime (kz, "kz")
        kz      <- approx(kz, xout=times, rule=2, ties=mean)$y
        
        # Iz.data is a matrix
        Iz.data <- outer(kz, z, 
                         FUN=function(x,y) exp(-x*y))
        Ibot    <- It.data * exp(-kz*Dt.data) 
        KZ      <- kz
        
    } else if (ncol(kz) %in% 3:4) {  # time, p, k1, k2 (2 fractions)
        
       if (ncol(kz) == 4){
        CheckTime (kz, "kz")
        KK <- data.frame(
          p  = approx(kz[,  1:2 ], xout=times, rule=2, ties=mean)$y,
          k1 = approx(kz[,c(1,3)], xout=times, rule=2, ties=mean)$y,
          k2 = approx(kz[,c(1,4)], xout=times, rule=2, ties=mean)$y
        )
       } else {
         p   = kz[,  1] 
         k1  = kz[,  2]
         k2  = kz[,  3]
         
       }
        
        KZ <- KK$p * KK$k1 + (1-KK$p) * KK$k2
        
        Iz.data <- outer(1:nrow(KK), z, 
                         FUN = function(x, y) 
                            KK$p[x]  * exp(-KK$k1[x]*y) + 
                         (1-KK$p[x]) * exp(-KK$k2[x]*y))
        
        Ibot    <- It.data*(    KK[,1]  * exp(-KK[,2]*Dt.data)+
                             (1-KK[,1]) * exp(-KK[,3]*Dt.data))  
        
    } else
      
      stop ("'kz' not correct")
    
  } else stop ("either `kz`, or `Iz.data` or `Iz.fun` should be defined") 
  
  names(Ibot) <- NULL
  
  if (max(Iz.data) > 1)
        stop ("'Iz.data' is RELATIVE light per depth and should be at most 1")

  # ----------------------------------------------------------------------------
  # morphology - volume for each layer, a function of depth
  # ----------------------------------------------------------------------------
  
  Mz.data  <- getdata(Mz.data, Mz.fun, "Mz", indep="depth")

  dVol     <- Mz.data*dz 
  meanSurf <- sum(dVol) / diff(range(zi))

  # ----------------------------------------------------------------------------
  # PI parameters
  # ----------------------------------------------------------------------------
  
  PI.type <- "depth"  # no need to transpose the data before calling PI.fun
  
  if (is.data.frame(PI.par) | is.matrix(PI.par)){
    err <- "'PI.par' columnnames ambiguous: only one of 'z' ('depth') and/or one of 't' ('time') should be present"
    
    # if the parameters are a data.frame or matrix, 
    # the first column should be depth or time
    
    if (nrow(PI.par) > 1) {
      
      nc     <- ncol(PI.par)
      zp     <- PI.par[,1]        # depth or time
      pnames <- colnames(PI.par)
      ppar   <- PI.par
      
      knownnames <- sum(c("depth", "time", "z", "t") %in% tolower(pnames))
      
      if (knownnames == 0 ) 
        stop("'PI.par' should have at least 't' (time), and/or 'z' (depth) as column names")
      
      if (knownnames > 2 ) stop(err)
      
      if ( knownnames == 2){ # both a function of z and t
        
        PI.type <- "zt"
        
        # interpolate over t and z for all parameters, and put in a list
        
        zz <- which (tolower(pnames) %in% c("depth", "z"))
        tt <- which (tolower(pnames) %in% c("time" , "t"))
        
        if (length(zz) > 1 | length (tt) > 1) stop(err)
        
        # all parameters are a matrix
        PI.par <- list()
        
        np     <- (1:nc)[-c(zz, tt)]  # position of the parameters
        
        CheckTime (ppar[ , c(tt, zz, 1)], "PI.par")
        
        for (i in 1: length(np) ) {  # loop over all parameters
          
          ii          <- np[i]
          # a matrix with dimension (times, depth)
          PI.par[[i]] <- map_xy(input_xyv = ppar[,c(tt, zz, i)],
                                output_x  = times, 
                                output_y  = z)$v
        } 
        names(PI.par) <- pnames[np]
        
      } else {  # nrow(PI.par) == 1
      
        if (tolower(pnames[1]) %in% c("depth", "z"))
          xout <- z
        
        else if (tolower(pnames[1] %in% c("time", "t")))  {
          CheckTime (ppar, "PI.par")
          xout     <- times
          PI.type <- "time"
        }
      
        # for all parameters:
        PI.par <- NULL
        
        for (i in 2:nc){
          PI.par <- cbind(PI.par, 
                          approx(zp, ppar[,i], 
                                 xout = xout, rule = 2, ties = mean)$y)
        }
        colnames(PI.par) <- pnames[-1]
        PI.par <- as.data.frame(PI.par)
      }
      
    } else # only one value or a vector
      
      PI.par <- as.vector(PI.par)

  } else if (is.list(PI.par)){
    
    nc     <- length(PI.par)
    
    pnames <- names(PI.par)
    ppar   <- PI.par
    
    knownnames <- sum(c("depth", "time", "z", "t") %in% tolower(pnames))
    
    if (knownnames != 2 ) stop (err)

    PI.type <- "zt"
    
    # interpolate over t and z for all parameters, and put in a list
    
    zz   <- which(tolower(pnames) %in% c("depth", "z"))
    tt   <- which(tolower(pnames) %in% c("time" , "t"))
    
    if (length(zz) > 1 | length (tt) > 1) stop(err)
    
    # for all parameters:
    PI.par <- list()
    np     <- (1:nc) [-c(zz, tt)]  # position of the parameters in 
    
    CheckTime (data.frame(ppar[[tt]], 1), "PI.par")
    
    for (i in 1: length(np)){  # loop over all parameters
      
      ii <- np[i]
      PI.par[[i]] <- map_tx(input_t  = ppar[[tt]], 
                            input_x  = ppar[[zz]], 
                            input_2D = ppar[[ii]],
                            output_t = times, 
                            output_x = z)$v
    } 
    names(PI.par) <- pnames[np]
  } 
    
  # ----------------------------------------------------------------------------
  # PI function
  # ----------------------------------------------------------------------------
  
  if (is.null(PI.fun)) {
    
    # Eilers-Peeters model
    PI.fun <- function(I) {
      
      ep <- I/((1/(PI.par[[1]]*PI.par[[2]]^2))*I^2+
                  (1/PI.par[[3]]-2/(PI.par[[1]]*PI.par[[2]]))*I+
                  (1/PI.par[[1]]))
      ep[is.na(ep)] <- 0
      
      return(ep)
    } 
  } else if (is.character(PI.fun)) {
    
    if (PI.fun == "fWebb") 
      
      PI.fun <- function(I) {
        pif <- PI.par[[1]]*PI.par[[2]]*(1-exp(-1*I/PI.par[[2]]))
        pif[is.na(pif)] <- 0
    
        return(pif)
      }
    
    else if (PI.fun == "fJP") 
      
      PI.fun <- function(I) {
        pif <- PI.par[[1]]*PI.par[[2]]*tanh(I/PI.par[[2]])
        pif[is.na(pif)] <- 0
        
        return(pif)
      }
    
    else if (PI.fun == "fPG") 
      
      PI.fun <- function(I) {
        pif <- PI.par[[3]]*(1-exp(-1*PI.par[[1]]*I/PI.par[[3]]))*
                              exp(-1*PI.par[[2]]*I/PI.par[[3]])
        pif[is.na(pif)] <- 0
        
        return(pif)
      }
    
    else stop("PI.fun not known, choose one of fWebb, fJP, fPG, fEP ")    
  
  } else if (! is.function(PI.fun))
      stop("PI.fun should either be NULL, the name of a function or a function")   

  # ----------------------------------------------------------------------------
  # Light profile at all depth intervals x times
  # ----------------------------------------------------------------------------
  
  if (! is.matrix(Iz.data))                        # Iz is a vector
    I_tz <- outer(It.data, Iz.data, 
                 FUN = "*")        
  
  else                                             # Iz is a matrix ([0-1])
    I_tz <- sweep(Iz.data, MARGIN = 1, STATS = It.data,  
                 FUN = "*")  
  
  # ----------------------------------------------------------------------------
  # Photosynthesis rates per volume PP_vol (mass/volume/time)
  # ----------------------------------------------------------------------------
  
  if (PI.type %in% c("time", "zt"))
    PP_vol  <- PI.fun(I_tz)*convFac                     # PS at these t x D  
  
  else if (PI.type == "depth")
    PP_vol  <- t(PI.fun(t(I_tz)))*convFac               # PS at these t x D  
  
  if (any (dim(PP_vol) - dim(I_tz) !=0))
    stop ("PI.fun(I) should return a matrix of same dimensions as light matrix")
  
  # ----------------------------------------------------------------------------
  # correct for changing water depth 
  # ----------------------------------------------------------------------------
  
  
  if (variableH) {
    
    # find where layers > max depth - NO LIGHT THERE 
    Z_water <- outer(Dt.data, z, 
                     FUN = function(x, y) y < x)
    
    I_tz    <- I_tz   * Z_water  # will 0 the layers deeper than max depth
    PP_vol  <- PP_vol * Z_water  # photosynthesis rates = 0
  } 
  
  # ----------------------------------------------------------------------------
  # take into account the volume -> total PP in the box, (mass/time)
  # ----------------------------------------------------------------------------
  
  PP_tot <- sweep(PP_vol,            
                  MARGIN = 2, 
                  STATS  = dVol, # dVol will be = dz*1 if Mz.data=1 
                  FUN    = "*")  # multiply PS with Volume

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # Output: timeseries
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ts <- data.frame(times      = times, 
                   PP         = rowSums(PP_tot)/meanSurf,  # PP per m2
                   TotalDepth = Dt.data, 
                   Isurf      = It.data, 
                   Ibot       = Ibot) 
  
  if (! is.null(KZ)) ts$kz <- KZ
  
  # ----------------------------------------------------------------------------
  # Output: profile
  # ----------------------------------------------------------------------------
  
  if (is.matrix(Iz.data)) 
    Iz.data <- apply(Iz.data, 
                     MARGIN = 2, 
                     FUN    = mean, 
                     na.rm  = TRUE)
  
  PP_v    <- colMeans(PP_vol)
  
  profile <- data.frame(z       = z, 
                        PP_v    = PP_v, 
                        Iz_I0   = Iz.data, 
                        surface = Mz.data)
  
  if (is.vector(PI.par) & ! PI.type =="zt")
    
    PI.par <- matrix(nrow     = length(z), 
                     ncol     = length(PI.par), 
                     byrow    = TRUE, 
                     data     = PI.par,
                     dimnames = list(NULL, names(PI.par)))
  
  if (PI.type == "depth") 
    profile <- data.frame(profile, PI.par)
  
  else if (PI.type == "time") 
    ts      <- data.frame(ts, PI.par)
  
  else if (PI.type == "zt"){ 
    profile <- data.frame(profile, 
                          sapply(PI.par, 
                                 FUN = colMeans))
    
    ts      <- data.frame(ts, 
                          sapply(PI.par, 
                                 FUN = rowMeans))
  }
  
  profile$dz <- dz
  
  # ----------------------------------------------------------------------------
  # Output: average over time
  # ----------------------------------------------------------------------------
  
  if (! is.null(avgOver)){
    
    if (inherits(what = "POSIXt",  x=ts$times) | 
        inherits(what = "POSIXct", x=ts$times)){
      
      # average time series data.frame
      ts <- average_timeseries(ts, 
                               avgOver  = avgOver, 
                               avgTime  = avgTime, 
                               datetime = "times",
                               value    = colnames(ts)[colnames(ts) != "times"])
      
      # average photosynthesis data I_tz
      LL  <- data.frame(times = times, 
                        I_tz)
      
      I_tz <- as.matrix(
           average_timeseries(LL, 
                              avgOver  = avgOver, 
                              avgTime  = avgTime, 
                              datetime = "times",
                              value    = colnames(LL)[colnames(LL) != "times"])[,-1])
      LL  <- data.frame(times = times, 
                        PP_vol)
      
      PP_vol <- as.matrix(
           average_timeseries(LL, 
                              avgOver  = avgOver, 
                              avgTime  = avgTime, 
                              datetime = "times",
                              value    = colnames(LL)[colnames(LL)!="times"])[,-1])
    }
  }
  
  # ----------------------------------------------------------------------------
  # All output combined
  # ----------------------------------------------------------------------------

  RES <- list(ts      = ts, 
              profile = profile, 
              light   = I_tz, 
              prod    = PP_vol)
  
  unit.default <- list(mass = "mgC", length = "m", time = "h", light = "uEinst/m2/s")
  unit.default[names(unit)] <- unit

  unit    <- as.list(unit.default)
  
  surface <- paste( unit$length, "^2", sep="")
  vol     <- paste( unit$length, "^3", sep="")
    
  pplab   <- paste(unit$mass, surface, unit$time, sep="/")
  ppvol   <- paste(unit$mass, vol,     unit$time, sep="/")
    
  varunits <-  c(PP           = pplab, 
              TotalDepth   = unit$length,  
              Isurf        = unit$light, 
              Ibot         = unit$light, 
              z            = unit$length, 
              PP_v         = ppvol, 
              Iz_I0        = "fraction", 
              surface      = surface,
              light        = unit$light, 
              prod         = ppvol)
    
  
  attr(RES, "unit")        <- unlist(unit)
  attr(RES, "varunits")    <- varunits
  attr(RES, "description") <- data.frame(
    
    variable    = c("ts.PP", "ts.TotalDepth",  "ts.Isurf", "ts.Ibot", 
                    "profile.z", "profile.PP_v", 
                    "profile.Iz_I0", "profile.surface",
                    "light", "prod"),
    description = c("integrated primary production", 
                  "Total water depth over which PP is estimated",
                  "Light at surface (I0)", 
                  "Light at the lowest grid point", 
                  "mean depth of the box at which photosynthesis is estimated", 
                  "Volumetric photosynthesis rate at depth z",
                  "fraction of light intensity at depth z (Iz/I0)", 
                  "surface area at middle of each box",
                  "light intensity for all times and depths",
                  "volumetric photosynthesis for all times and depths"),
    unit        = varunits)
  
  class(RES) <- c("integratedPP", class(RES))
  return(RES)
}

## =============================================================================
## =============================================================================
## General plotting functions for depth-integrated photosynthesis
## =============================================================================
## =============================================================================

## =============================================================================

plot.integratedPP <- function(x, ..., 
                              type = "l", 
                              las = 1, which = "ts"){

  # number of PP data series to plot
  
  ldots <- list(...)
  ndots <- names(ldots)
  nx <- 0
  x2 <- list()
  
  lld <- ldots
  
  if (length(lld)) {
    
    for (i in 1:length(lld)) 
      
       if (inherits(lld[[i]], "integratedPP")) {
          
          x2[[nx <- nx + 1]] <- ldots[[i]]
          ldots[[i]]         <- NULL
          names(x2)[nx]      <- ndots[i]
       }
  }
  
  # which variables to plot

  if (length(which) == 0)
    stop("'which' cannot be NULL; should be 'profile', 'ts' or a selection of names in these")
  
  which.ts <- which.prof <- NULL
  
  colname.ts <- colnames(x$ts)
  
  if ("ts" %in% which)   
    which.ts <- colname.ts[-1]  # first one is time
  else
    which.ts <- colname.ts[which(colname.ts %in% sub("ts.", "", which))]
  
  colname.prof <- colnames(x$profile)
  
  if ("profile" %in% which) 
    which.prof <- colname.prof[-c(1, ncol(x$profile))]
  else
    which.prof <- colname.prof[which(colname.prof %in% sub("profile.", "", which))]

  # number of figures to plot
  mfrow <- ldots$mfrow
  if (is.null(names(ldots)["mfrow"])) {
    nv    <- length(which.ts) + length(which.prof)
    nc    <- min(ceiling(sqrt(nv)), 3)
    nr    <- min(ceiling(nv/nc), 3)
    mfrow <- c(nr, nc)
    pm <- par(mfrow=mfrow)
    
  }  else pm <- par(mfrow=ldots$mfrow)

  # for labelling the plots

  lab    <- "times"
  un     <- attributes(x)$unit
  units  <- attributes(x)$varunits

  unitunknown <- c(which.ts  [!which.ts   %in% names(units) ],
                   which.prof[!which.prof %in% names(units) ])
  newunit <- rep("-", times = length(unitunknown))
  names(newunit) <- unitunknown
  units <- c(units, newunit)
  
  main <- c(PP         = "Integrated production", 
            TotalDepth = "Total depth", 
            Isurf      = "Light at surface", 
            Ibot       = "Light at bottom", 
            z          = "depth", 
            PP.v       = "mean volumetric production",
            Iz_I0      = "mean light penetration (relative)", 
            surface    = "horizontal surface area")
  
  newmain <- c(which.ts   [!which.ts   %in% names(main) ],
               which.prof [!which.prof %in% names(main) ])
  
  names(newmain) <- newmain
  main <- c(main, newmain)
  
  # plot the time series data
  
  if (length(which.ts)){
    xlab  <- "times"
    
    for (var in which.ts){
      
      y <- x$ts[,var]
      if (nx > 0){ 
        for (j in 1:nx) 
          y <- cbind(y, x2[[j]]$ts[,var])
      }
       do.call("matplot", 
               c(alist(x = x$ts$times, y = y, 
                 xlab = xlab, ylab = units[var], 
                 main = main[var], type = type, las = las), 
               ldots))
    }
  }
  
  if (length(which.prof)){
      
    for (var in which.prof){
      xx <- x$profile[,var]
      y  <- x$profile$z
      if (nx > 0) { 
        for (j in 1:nx) {
          xx <- cbind(xx, x2[[j]]$profile[,var])
          y  <- cbind(y, x2[[j]]$profile$z)
        }
      }
      do.call("matplot", 
              c(alist(x = xx, y = y, 
                ylab = paste("depth", un["length"]), xlab = units[var], 
                main = main[var], ylim = rev(range(x$profile$z)), 
                type = type, las = las), 
              ldots))
     }
  }
  par(mfrow=pm)
}

## =============================================================================
## Two-dimensional plot of photosynthesis result
## =============================================================================

image2D.integratedPP <- function(z, las = 1, 
                              ...){
  
  mfrow <- ldots$mfrow
  if (is.null(names(list(...))['mfrow'])) {
    mfrow <- c(1,2)
    pm <- par(mfrow=mfrow)
  }
  depth <- z$profile$z
  times <- z$ts$times
  xlab  <- "times"
  un       <- attributes(z)$unit
  varunits <- attributes(z)$varunits
  ylab  <- varunits["TotalDepth"]
  ylim  <- rev(range(depth))  

  image2D(x    = times, 
          y    = depth, 
          z    = z$light, 
          ylim = ylim, 
          ylab = ylab, 
          clab = varunits["light"], 
          main = "light", las = las, ...)
  image2D(x    = times, 
          y    = depth, 
          z    = z$prod, 
          ylim = ylim, 
          ylab = ylab, 
          clab = varunits["prod"],
          main = "volumetric production", las = las, ...)
  if (! is.null(mfrow)) 
    par(mfrow=pm)
}
