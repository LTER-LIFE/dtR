# ==============================================================================
# ==============================================================================
# Interpolate spatial time series defined in (xy) points to other (xy) points 
# ==============================================================================
# ==============================================================================

interpolate_xyt <-  function(input.xytv,  # longitude (x), latitude (y), time, value
                             output.xy, 
                             output.t, 
                             wide = TRUE){  # wide or long format
  
  attrs <- NULL
  
  if (inherits(input.xytv, "dtLife"))
    attrs <- attributes(input.xytv)
  
  # input.xytv: 4 columns: longitude (x), latitude (y), time, value
  if (ncol(input.xytv) != 4) 
    stop ("'input.xytv' should have 4 columns: longitude (x), latitude (y), time, value")

  # different xy values

  input.xy  <- unique(input.xytv[, 1:2])
  ni <- nrow(input.xy)
  
  if (is.vector(output.xy))
    output.xy <- matrix(nrow=1, data=output.xy)
  else 
    output.xy <- unique(output.xy[,1:2])
  coord.xy  <- output.xy
  no <- nrow(output.xy)
  
  ##### step 1: weighting points #####
  # distance between two points
  # euclidean_distance
  
  DD1 <- outer(output.xy [, 1], input.xy[, 1], FUN="-")
  DD2 <- outer(output.xy [, 2], input.xy[, 2], FUN="-")
  Distance <- sqrt(DD1^2+DD2^2)
  rm(list=c("DD1", "DD2"))
  
  # find three closest points 
  # which data sets to use for each output station, 
  imin <- min(3, ni)  # at most 3 columns selected
  i.close <- matrix(nrow=no, byrow=TRUE, 
                    data=unlist(apply(Distance, 
                                      MARGIN = 1, 
                                      FUN    = order, 
                                      simplify=FALSE)))[,1:imin]

  # weighing values ~ inverse distance
  w.close <- NULL
  for (i in 1:ncol(i.close))
     w.close <- cbind(w.close, 
                      Distance[cbind(1:nrow(i.close), i.close[,i])])
  ii <- which(w.close == 0)
  w.close <- 1/w.close
  
  # rescale so that sum=1
  w.close   <- w.close/rowSums(w.close, na.rm=TRUE)
  
  if (length(ii)) w.close[ii] <- 1
  w.close[is.nan(w.close)] <- 0
  
  # columns that are used in weighing
  i.unique <- unique(as.vector(i.close))
  
  ##### step2: interpolate to times for all inputs #####
  xytin <- matrix(nrow=length(output.t), ncol=max(i.unique))
  
  for (i in i.unique){
    # interpolate to required times
    ii <- which(input.xytv[,1] == input.xy[i,1] & 
                input.xytv[,2] == input.xy[i,2])
    V.out  <- approx(x=input.xytv[ii,3], y=input.xytv[ii,4], xout=output.t)$y
    xytin[,i] <- V.out
  }
  
  stnames <- paste("st", 1:no, sep="")
  Stations <- data.frame(stations=stnames,
                         x=coord.xy[,1], 
                         y=coord.xy[,2])
  row.names(Stations) <- NULL
  
  if (wide){
    result     <- matrix(nrow = length(output.t), 
                         ncol = 1+nrow(output.xy))
    result[,1] <- output.t

    for (i in 1:no){
      iu <- xytin[ ,i.close[i,]]  # data that need to be averaged
      wu <- w.close[i,]           # averaging weights
      result[,i+1] <- rowSums(sweep(iu, MARGIN=2, STATS=wu, FUN="*"))
    }
  
    colnames(result) <- c("time", stnames)
  } else {
    result     <- NULL
    for (i in 1:no){
      iu  <- xytin[ ,i.close[i,]]  # data that need to be averaged
      wu  <- w.close[i,]           # averaging weights
      res <- rowSums(sweep(iu, MARGIN=2, STATS=wu, FUN="*"))
      result <- rbind(result,  
                      cbind(Stations[i,], output.t, res, row.names=NULL))
    }
    
  }
  atout <- attributes(result)
  if (! is.null(attrs))
    attributes(result) <- c(atout, attrs[!names(attrs) %in% names(atout)])
  
  attr(result, "stations") <- Stations
  attr(result, "processing") <- c(
    atout$processing, paste("interpolated from 2D-time input, at:", Sys.time()))
  result
}

# ==============================================================================
# ==============================================================================
# Maps xy points to other xy points - uses least distance weighing 
# (3 closest pts) 
# ==============================================================================
# ==============================================================================

map_xy <-  function(input.xyv, # latitude (x), longitude (y), value
                    input.x, input.y, input.2D,
                    output.xy, 
                    output.x, output.y){
  # check input
  if (missing (input.xyv)){
    if (missing(input.x) | missing(input.y) | missing(input.2D) )
      stop("either input.xyv should be input or (input.x, input.y & input.2D)")
  }
    
  # check output
  if (missing (output.xy)){
    if (missing(output.x) | missing(output.y) )
      stop("either output.xy should be input or (output.x & output.y)")
  }
  
  # call correct function
  if (! missing(input.xyv) & 
      ! missing (output.xy))
    result <- interpolate_xy_xy(input.xyv = input.xyv, # latitude (x), longitude (y), value
                                output.xy = output.xy) 
  
  else if (! missing(input.x) & ! missing(input.y) & ! missing(input.2D) & 
           ! missing(output.xy))
    result <- interpolate_2D_xy (input.x   = input.x, 
                                 input.y   = input.y, # latitude (x), longitude (y)
                                 input.2D  = input.2D,         #  depth
                                 output.xy = output.xy)
  
  else if (!missing(input.xyv) & ! missing(output.x) & 
           !missing(output.y))
    result <- interpolate_xy_2D(input.xyz = input.xyv, 
                               output.x  = output.x, 
                               output.y  = output.y)
  
  else if (! missing(input.x) & ! missing(input.y) & ! missing(input.2D) & 
           ! missing(output.x) & !missing(output.y))
    result <- interpolate_2D_2D (input.x  = input.x, 
                                 input.y  = input.y, 
                                 input.2D = input.2D, 
                                 output.x = output.x, 
                                 output.y = output.y)
  
  else 
    stop ("one of the required arguments is missing")
  
  attr(result, "processing") <- paste("remapped , at:", 
                                      Sys.time())
  
  return(result)
      
}

# ==============================================================================

interpolate_xy_xy <-  function(input.xyv, # latitude (x), longitude (y), value
                            output.xy){
  
  if (ncol(input.xyv) != 3) 
    stop ("'input.xyv' should have 3 columns: latitude (x), longitude (y), value")

  ##### step 1: weighting points #####
  # distance between two points 

  input.xy  <- input.xyv[, 1:2]
  ni <- nrow(input.xy)
  if (ni == 0) stop("input.xyv is empty")
  
  output.xy <- unique(output.xy)
  coord.xy <- output.xy[,1:2]
  no <- nrow(output.xy)
  if (no == 0) stop("output.xy is empty")

  # euclidean_distance
  
  DD1 <- outer(output.xy [, 1], input.xy[, 1], FUN="-")
  DD2 <- outer(output.xy [, 2], input.xy[, 2], FUN="-")
  Distance <- sqrt(DD1^2+DD2^2)
  rm(list=c("DD1", "DD2"))
  
  
  # which data sets to use for each output station, 
  imin <- min(3, ni)  # at most 3 columns selected
  
  # find three closest points - this is slightly faster, but less clear
  #  ir   <- 1: no
  #  C1 <-apply(Distance, MARGIN=1, FUN=which.min)
  #  Distance[cbind(ir, C1)] <- NA
  
  #  C2 <- apply(Distance, MARGIN=1, FUN=which.min)
  #  if (! length(C2)) C2 <- NA else Distance[cbind(ir, C2)] <- NA
  
  #  C3 <- apply(Distance, MARGIN=1, FUN=which.min)
  #  if (! length(C3)) C3 <- NA 
  
  #    i.close   <- cbind(C1, C2, C3)
  #    w.close   <- 1/cbind(Distance[cbind(ir,C1)], 
  #                         Distance[cbind(ir,C2)], 
  #                         Distance[cbind(ir,C3)]) 

  i.close <- matrix(nrow=no, byrow=TRUE, 
                    data=unlist(apply(Distance, 
                                      MARGIN   = 1, 
                                      FUN      = order, 
                                      simplify = FALSE)))[,1:imin]
# weighing values ~ inverse distance
  w.close <- NULL
  for (i in 1:ncol(i.close))
    w.close <- cbind(w.close, 
                     Distance[cbind(1:nrow(i.close), i.close[,i])])
  
  ii <- which(w.close == 0)
  
  # weighing values ~ inverse distance
  w.close <- 1/w.close
  
  # rescale so that sum=1
  w.close   <- w.close/rowSums(w.close, na.rm=TRUE)
  if (length(ii)) w.close[ii] <- 1
  w.close[is.nan(w.close)] <- 0
  
  
  # columns that are used in weighing
  i.unique <- unique(as.vector(i.close))

  ##### step2: interpolate #####

  result <- sapply(1:no, FUN=function(i) sum(w.close[i,]*input.xyv[i.close[i,],3]))
  result <- data.frame(coord.xy, value=result)
  colnames(result)[3]<- colnames(input.xyv)[3]
  
  result
}

# ==============================================================================
# ==============================================================================
# interpolate 2D gridded points to xy points
# ==============================================================================
# ==============================================================================

interpolate_2D_xy <-  function(
                        input.x, input.y, # latitude (x), longitude (y)
                        input.2D,         #  depth
                        output.xy){

  dd <- dim(input.2D)
  if (dd[1] != length(input.x)) stop ("length of 'input.x' should be = nrow(input.2D)")
  if (dd[2] != length(input.y)) stop ("length of 'input.y' should be = ncol(input.2D)")
  dlon <- diff (input.x)  # diff(input.2D$longitude)
  dlat <- diff (input.y)  # diff(input.2D$latitude)
  
  equidistant <- diff(range(dlon)) == 0 &
                 diff(range(dlat)) == 0
  
  if (! equidistant){  # SLOW!
   input.2D.grid        <- expand.grid(input.x, input.y)
   names(input.2D.grid) <- c("longitude", "latitude")
   input.2D.grid$depth  <- as.vector(input.2D)
   input.2D.grid        <- na.omit(input.2D.grid)

   gg    <- na.omit(output.xy)
  
   Result <- map_xy(input.xyv = input.2D.grid, 
                    output.xy = gg[,1:2])
  } else {
    lonmin <- min(input.x)
    latmin <- min(input.y)
    ilon   <- function(x) 1+(x-lonmin)/dlon[1]
    ilat   <- function(x) 1+(x-latmin)/dlat[1]
    ix  <- sapply(output.xy[,1], FUN=ilon)
    iy  <- sapply(output.xy[,2], FUN=ilat)
    
    ### CHECK IF THIS IS WHAT WE WANT !!!! - EXTRAPOLATION
    iix <- pmin(length(input.x) , as.integer(ix))  # index to pt on right
    iix <- pmax(1 , iix)  # index to pt on left
    iiy <- pmin(length(input.y) , as.integer(iy))  # index to pt above
    iiy <- pmax(1, iiy)  # index to pt below

    dx  <- (output.xy[,1]-input.x[iix])  /dlon[1]
    dy  <- (output.xy[,2]-input.y[iiy] ) /dlat[1]
    dxvalue <- input.2D[,-1] - input.2D[,-ncol(input.2D)]
    dxvalue <- cbind(0, dxvalue, 0)
    dyvalue <- input.2D[-1,] - input.2D[-nrow(input.2D),]
    dyvalue <- rbind(0, dyvalue, 0)
    
    newz <- input.2D[cbind(iix, iiy)]+
                0.5*(dxvalue[cbind(iix, iiy)]*dx +
                     dyvalue[cbind(iix, iiy)]*dy) 
    row.names(output.xy) <- NULL
    Result <- data.frame(output.xy, z=newz)
    
  }
  Result 
}

# ==============================================================================
# ==============================================================================
# maps xy values on a 2D regular grid
# ==============================================================================
# ==============================================================================

interpolate_xy_2D <- function(input.xyz, 
                   output.x=NULL, output.y=NULL){
  
  if (is.null(output.x)) 
    output.x <- sort(unique(input.xyz[,1]))
  if (is.null(output.y)) 
    output.y <- sort(unique(input.xyz[,2]))
  
  xo  <- unique(output.x)
  xor <- range(xo)
  dxo <- diff(xor)/(length(xo)-1)
  
  yo <- unique(output.y)
  yor <- range(yo)
  dyo <- diff(yor)/(length(yo)-1)
  
  dxi <- min(diff(sort(unique(input.xyz[,1]))))
  dyi <- min(diff(sort(unique(input.xyz[,2]))))
  
  if (dxi > dxo | dyi > dyo) {
    Ux <- diff(range(diff(sort(unique(input.xyz[,1])))))
    Uy <- diff(range(diff(sort(unique(input.xyz[,2])))))
    if (Ux == 0 & Uy == 0) {  # uniform grid
       OO <- interpolate_xy_2D(input.xyz, output.x=NULL, output.y=NULL)
       return(interpolate_2D_2D(OO$x, OO$y, OO$z, output.x, output.y))  # RECURSIVE!!!??
    }
    else stop("cannot map: input resolution too low and nonuniform input grid")
  } else {
  # map input to output grid cells
  ix <- as.integer(1+(input.xyz[,1]-xor[1])/dxo)
  iy <- as.integer(1+(input.xyz[,2]-yor[1])/dyo)
  
  ix.lack <- which(!1:length(output.x) %in% unique(ix))
  iy.lack <- which(!1:length(output.y) %in% unique(iy))
  
  lx <- length(ix.lack)
  ly <- length(iy.lack)
  if (lx+ly) {
    if (!lx) ix.lack <- 1
    if (!ly) iy.lack <- 1
    ix.lack <- rep(ix.lack, length.out=lx+ly)
    iy.lack <- rep(iy.lack, length.out=lx+ly)
    ix <- c(ix, ix.lack)
    iy <- c(iy, iy.lack)
    toadd <- input.xyz[1:(lx+ly), ]
    toadd[,] <- 0
    input.xyz <- rbind(input.xyz, toadd)
    }
  
  # sum all variables within a grid cell, and count the number of elements
  zsum <- tapply(input.xyz[,3], 
                 INDEX = list(ix, iy), 
                 FUN   = sum)
  zlen <- tapply(input.xyz[,3], 
                 INDEX = list(ix, iy), 
                 FUN   = length)
  
  # create x and y values - select only values 
  x  <- xor[1] + as.integer(rownames(zsum))* dxo
  ii <- which(x >= xor[1] & x<= xor[2])
  
  y  <- yor[1] + as.integer(colnames(zsum))* dyo
  jj <- which(y >= yor[1] & y<= yor[2])
  
  out <- list(x[ii], y[jj], zsum[ii,jj] / zlen[ii,jj])
  names(out) <- colnames(input.xyz)[1:3]
  return(out)
  }
}

# ==============================================================================
# ==============================================================================
# maps gridded values on a regular grid
# ==============================================================================
# ==============================================================================

interpolate_2D_2D <- function(input.x, input.y, input.2D, 
                              output.x=NULL, output.y=NULL){
  
  if (is.null(output.x)) 
    output.x <- sort(unique(input.x))
  if (is.null(output.y)) 
    output.y <- sort(unique(input.y))
  
 # output grid size
  xo  <- unique(output.x)
  xor <- range(xo)

  yo <- unique(output.y)
  yor <- range(yo)

  xi  <- unique(input.x)
  xir <- range(xi)
  
  yi <- unique(input.y)
  yir <- range(yi)
  
  if (diff(range(diff(xi))) != 0 | diff(range(diff(yi))) != 0)
    stop("cannot interpolate - grid sizes are not constant")
  
  # input grid size
  dxi <- min(diff(sort(unique(input.x))))
  dyi <- min(diff(sort(unique(input.y))))
  
  # map output to input grid cells
  ix <- pmin(length(input.x), as.integer(1+(output.x-xir[1])/dxi))
  iy <- pmin(length(input.y), as.integer(1+(output.y-yir[1])/dyi))

  ix <- pmax(1, ix)
  iy <- pmax(1, iy)
  
  in2D <- cbind(#input.2D[,1], 
                input.2D, 
                input.2D[,ncol(input.2D)])
  in2D <- rbind(#in2D[1,],
                in2D,     
                in2D    [nrow(in2D),])
  
  IG   <- as.matrix(expand.grid(ix, iy))
  IGx1 <- t(t(IG) + c(1,0))
  OX   <- rep(output.x, times=length(output.y))
  zx   <- in2D[IG] + (in2D[IGx1]-in2D[IG])*(OX-input.x[IG[,1]])/dxi  
  zx   <- matrix(nrow=length(output.x), ncol=length(output.y), data=zx)
#  ZX<<- zx
  
  IG   <- as.matrix(expand.grid(iy, ix))[, c(2,1)]
  IGy1 <- t(t(IG) + c(0,1))
  OY   <- rep(output.y, times=length(output.x))
  zy   <- in2D[IG] + (in2D[IGy1]-in2D[IG])*(OY-input.y[IG[,2]])/dyi  
  zy   <- matrix(nrow=length(output.x), ncol=length(output.y), data=zy, byrow=TRUE)
  zz   <- 0.5*(zx+zy)
#  ZY<<- zy
  
  out <- list(x=output.x, y=output.y, z=(zx+zy)/2)
  out
}


# ==============================================================================

interpolate_xt <-  function(input.xtv, # position (x), time (t), value
                            output.x, output.t){
  
  if (ncol(input.xtv) != 3) 
    stop ("'input.xtv' should have 3 columns: position (x), time (t), value")
 
   ############## TO DO ############### 
 }



######################## WORK IN PROGRESS  ###################################

# ==============================================================================
# ==============================================================================
# pick time series for a specific location
# ==============================================================================
# ==============================================================================

# ==============================================================================
# Check if a point falls in a triangle
# ==============================================================================

ptintri <- function(P = c(1,  3), 
                    Triangle.x,    # x-points
                    Triangle.y){ 
  
  # A =[ux  vx  wx  # Trianglex
  #     uy  vy  wy
  #      1   1   1]
  A <- rbind(Triangle.x, 
             Triangle.y, 
             rep(1, 3))
  B <- c(P, 1)
  
  all(solve(A, B) > 0)
}

# ==============================================================================
# calculates coefficients for triangular interpolation
# ==============================================================================

ptcoeff <- function(P = c(1,  3), 
                    Triangle.x, 
                    Triangle.y){   
  # A =[ux  vx  wx
  #     uy  vy  wy
  #      1   1   1]
  A <- rbind( Triangle.x, Triangle.y, rep(1, 3))
  B <- c(P, 1)
  
  solve(A, B)
}

# ==============================================================================
# Creates all unique triangles from a set of points
# ==============================================================================

CreateTriangles <- function(n = 5) {   
  Triangles <- expand.grid(1:n, 2:n, 3:n)
  
  # check for same points- in each set of 3
  TT        <- apply(Triangles, 
                     MARGIN = 1, 
                     FUN    = anyDuplicated)  
  Triangles <- Triangles[TT==0, ]       # remove duplicates                      
  
  # check for same rows
  Triangles <- t(apply(Triangles, 
                       MARGIN = 1, 
                       FUN    = function(x) x[order(x)]))
  
  TT        <- unique(apply(Triangles, 
                            MARGIN = 1, 
                            FUN    = function(x) paste(x, collapse="_")))
  
  Triangles <- matrix(ncol  = 3, 
                      data  = as.integer(unlist(strsplit(TT, "_"))), 
                      byrow = TRUE)
  
  return(Triangles)
}

# ==============================================================================
# https://codeplea.com/triangular-interpolation
# Barycentric Coordinates
# as not enough evenly spread weather stations to cover the entire region.
# ==============================================================================


#pick_xyt <- function(longitude=5.2, latitude=53.2, 
#                     input=NULL,
#                     output.t=unique(input$datetime),
#                     what="waterheight",
#                     shape=Wadshape,
#                     plotit = FALSE, 
#                     n      = NULL,
#                     verbose = FALSE,
#                     ...){
  # 1. check if longitude, latitude fall in the Waddensea. 
  #    Check if 'station' attribute is present in object
#  if (what == "waterheight"){
#    if (is.null(input)) input <- HeightHR
#    pickdata(longitude=longitude, latitude=latitude, 
#             output.t=output.t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
    
#  } else if (what == "temperature") {
#    if (is.null(input)) input <- TempHR
#    pickdata(longitude=longitude, latitude=latitude, 
#             output.t=output.t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
    
#  } else if (what %in% colnames(RWSbiogeo2021)){
#    if (is.null(input)) {
#      select <- c("station", "longitude", "latitude",  "datetime", what) 
#      input <- RWSbiogeo2021[, select]
#      attributes(input)$format <- "long"
#      input <- dtreshape(input, swap="station")
#    }
#    pickdata(longitude=longitude, latitude=latitude, 
#             output.t=output.t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
#  } else {
#    if (is.null(input)) stop (" 'input' should be provided for variable ", what)
#    pickdata(longitude=longitude, latitude=latitude, 
#             output.t=output.t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
#  }
#}

# ==============================================================================
# Interpolation using barycentric coordinates
# ==============================================================================

pickdata <- function(longitude=5.2, latitude=53.2, 
                     depth= NA,  # depth versus NAP
                     from="2021-01-01", to="2021-02-01", by=3600,
                     output.t=seq(as.POSIXct(from),     # default = hourly
                                  as.POSIXct(to), 
                                  by=by),
                     what = "watertemperature",
                     input  = WadTempLR,
                     stats = attributes(input)$stations,
                     n      = NULL,
                     shape  = Shape,
                     verbose = FALSE, ...){
  
  if (!attributes(input)$format == "wide"){
    
    stop ("can only do this in wide format")
    
  }
  
  output.t <-  output.t [order(output.t)]
  input    <- input[order(input$datetime), ]
  
  # Make a selection of the input based on times
  if (is.numeric(output.t)) # index to the datetime
    input <- input[output.t,]
  else {
#    ii <- which(input$datetime >= min(output.t) &     
#                input$datetime <= max(output.t))
#    
#    if (! length(ii))
#      stop ("cannot interpolate: no overlap in time between input and output.t")
#    
#    ir <- range(ii) + c(-1,1)   # one more before and after
#    ir <- pmin(nrow(input), ir)
#    ir <- pmax(ir, 1)
#    input  <- input[ir[1]:ir[2], ]
    ist     <- which(colnames(input) %in% stats$station)
    notNAs  <- apply(input[,ist], MARGIN=2, FUN=function(x) sum(!is.na(x))) 
    iremove <- names(notNAs)[which(notNAs <= 2)]
    if (length(iremove)){
      if (verbose) warning("removing station ", iremove, " as this has too many NAs")
      input <- input[, -which(colnames(input) %in% iremove)]
      stats    <- stats[! stats$station %in% iremove, ]   # station information
    }
  }

  # distance between the point and all stations - order stations
  asp      <- 1/cos((mean(latitude)*pi)/180)        # y/x aspect ratio
  Dist2    <- (stats$longitude - longitude)^2+      # distance^2
              (stats$latitude  - latitude)^2*asp^2
  
  # create surrounding triangles based on closest points
  if (is.null(n)) n <- 20
  ntri     <- min(n, nrow(stats))
  
  Dsort    <- order(Dist2)[1:ntri]  # selection of closest points
  stats    <- stats[Dsort,]
  
  # Creates all unique triangles from a set of points
  Triangles   <- CreateTriangles(ntri)
  Triangles.x <- matrix(ncol=3, 
                        data=stats$longitude[Triangles])
  Triangles.y <- matrix(ncol=3, 
                        data=stats$latitude [Triangles])
  
  # from triangles that embrace point: find one with smallest mean distance
  P    <- c(longitude, latitude)
  Keep <- meanD <-  NULL
  ic   <- 0
  
  for (i in 1:nrow(Triangles)){
    if (ptintri(P, 
                Triangles.x[i,], Triangles.y[i,])) {  # point falls in triangle
      
      Keep  <- c(Keep, i)
      meanD <- c(meanD,       # mean distance of point to triangle edges
                 mean(sqrt((Triangles.x[i,]-longitude)^2 +
                          ((Triangles.y[i,]-latitude)*asp)^2)))
    }
  }
  
  # point is enclosed by triangle(s) 
  if (! is.null(Keep)){
    
    i          <- Keep[which.min(meanD)]
    statselect <- stats[Triangles[i, ],]
    
    # weighing coefficients
    wu       <- ptcoeff(P, Triangles.x[i,], Triangles.y[i,])
    legtitle <- "surrounding triangle"
    
  } else {                    # point is NOT enclosed by triangle(s) 
    if (verbose) 
      warning("no surrounding points found - used the 3 closest points")
    
    statselect <- stats[1:3, ]  # 3 closest points
    
    # weighing coefficients : inverse distance 
    distance   <- sqrt((statselect$longitude-longitude)^2 +
                         ((statselect$latitude-latitude)*asp)^2)
    wu <- 1/distance
    wu <- wu/sum(wu)
    legtitle <- "3 closest points"
  }
  
  statnames  <- statselect$station
  
  # select only input stations we need 
  input <- input[,c("datetime", statnames)]         # stations
  
  # interpolate time series to the requested output times
  res <- apply(input[,-1], 
               MARGIN = 2,
               FUN    = function(v) 
                 approx(x    = input[ ,1], 
                        y    = v, 
                        xout = output.t, 
                              rule = 2)$y)     
  
  # take weighted average of input data: multiply all columns (margin=2)
  # with correct weighing factor (wu) and sum over the rows
  res <- data.frame(datetime = output.t, 
                    rowSums(sweep(res, 
                                  MARGIN = 2, 
                                  STATS  = wu, 
                                  FUN    = "*"), 
                            na.rm = TRUE))
  colnames(res)[2] <- what
  
  if (! is.na(depth))
    res[,2] <- pmax(res[,2], -depth)
  
  res <- res[order(res[,1]), ]
  
  attributes(res)$inputcol  <- colnames(input)[-1]
  attributes(res)$selected  <- statselect
  attributes(res)$weights   <- wu
  
  
  res
}

# pick_xyt(plotit=TRUE)
# A <- pick_xyt(plotit=TRUE, input=TempHR, output.t=1:40)
# A <- pick_xyt(plotit=TRUE, input=TempHR, output.t=seq(from=as.POSIXct("2021-03-01"),to= as.POSIXct("2021-05-01"), by=60))
# A <- pick_xyt(plotit=TRUE, input=TempHR)
# A <- pick_xyt(plotit=TRUE, input=TempLR, output.t = seq("15-01-2021", "15-12-2021", by=3600*24))
#  A <- pick_xyt(plotit=TRUE, input=TempHR, output.t=seq(from=as.POSIXct("2021-03-01"),to= as.POSIXct("2021-03-31"), by=60), verbose=TRUE, ylab="m", latitude=53.15)

#  A <- pick_xyt(plotit=TRUE, what="NH4", input=RWSbiogeo2021, verbose=TRUE, ylab="m", latitude=53.15)

