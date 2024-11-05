# ==============================================================================
# ==============================================================================
# Maps xy points to other xy points - large input datasets
# ==============================================================================
# ==============================================================================

map_xy <-  function(input_xyv, # latitude (x), longitude (y), value
                    input_x, input_y, input_2D,
                    output_xy, 
                    output_x, output_y){
  # check input
  if (missing (input_xyv)){
    if (missing(input_x) | missing(input_y) | missing(input_2D) )
      stop("either input_xyv should be input or (input_x, input_y & input_2D)")
  }
    
  # check output
  if (missing (output_xy)){
    if (missing(output_x) | missing(output_y) )
      stop("either output_xy should be input or (output_x & output_y)")
  }
  
  # call correct function
  if (! missing(input_xyv) & 
      ! missing (output_xy))
    result <- interpolate_xy_xy(input_xyv = input_xyv, # latitude (x), longitude (y), value
                                output_xy = output_xy) 
  
  else if (! missing(input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing(output_xy))
    result <- interpolate_2D_xy (input_x   = input_x, 
                                 input_y   = input_y, # latitude (x), longitude (y)
                                 input_2D  = input_2D,         #  depth
                                 output_xy = output_xy)
  
  else if (!missing(input_xyv) & ! missing(output_x) & 
           !missing(output_y))
    result <- interpolate_xy_2D(input_xyz = input_xyv, 
                               output_x  = output_x, 
                               output_y  = output_y)
  
  else if (! missing(input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing(output_x) & !missing(output_y))
    result <- interpolate_2D_2D (input_x  = input_x, 
                                 input_y  = input_y, 
                                 input_2D = input_2D, 
                                 output_x = output_x, 
                                 output_y = output_y)
  
  else 
    stop ("one of the required arguments is missing")
  
  attr(result, "processing") <- paste("remapped , at:", 
                                      Sys.time())
  
  return(result)
      
}

# ==============================================================================

interpolate_xy_xy <-  function(input_xyv, # latitude (x), longitude (y), value
                            output_xy){
  
  if (ncol(input_xyv) != 3) 
    stop ("'input_xyv' should have 3 columns: latitude (x), longitude (y), value")

  ##### step 1: weighting points #####
  # distance between two points 

  input_xy  <- input_xyv[, 1:2]
  ni <- nrow(input_xy)
  if (ni == 0) stop("input_xyv is empty")
  
  output_xy <- unique(output_xy)
  coord.xy <- output_xy[,1:2]
  no <- nrow(output_xy)
  if (no == 0) stop("output_xy is empty")

  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xy[,1]) < min(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  if (min(input_xy[,1]) > max(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xy[,2]) < min(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  if (min(input_xy[,2]) > max(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  # euclidean_distance
  
  DD1 <- outer(output_xy [, 1], input_xy[, 1], FUN="-")
  DD2 <- outer(output_xy [, 2], input_xy[, 2], FUN="-")
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

  result <- sapply(1:no, FUN=function(i) sum(w.close[i,]*input_xyv[i.close[i,],3]))
  result <- data.frame(coord.xy, value=result)
  colnames(result)[3]<- colnames(input_xyv)[3]
  
  result
}

# ==============================================================================
# ==============================================================================
# interpolate 2D gridded points to xy points
# ==============================================================================
# ==============================================================================

interpolate_2D_xy <-  function(
                        input_x, input_y, # latitude (x), longitude (y)
                        input_2D,         #  depth
                        output_xy){

  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_x) < min(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  if (min(input_x) > max(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_y) < min(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  if (min(input_y) > max(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  
  dd <- dim(input_2D)
  if (dd[1] != length(input_x)) stop ("length of 'input_x' should be = nrow(input_2D)")
  if (dd[2] != length(input_y)) stop ("length of 'input_y' should be = ncol(input_2D)")
  dlon <- diff (input_x)  # diff(input_2D$longitude)
  dlat <- diff (input_y)  # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)) == 0 &
                 diff(range(dlat)) == 0
  
  if (! equidistant){  # SLOW!
   input_2D.grid        <- expand.grid(input_x, input_y)
   names(input_2D.grid) <- c("longitude", "latitude")
   input_2D.grid$depth  <- as.vector(input_2D)
   input_2D.grid        <- na.omit(input_2D.grid)

   gg    <- na.omit(output_xy)
  
   Result <- map_xy(input_xyv = input_2D.grid, 
                    output_xy = gg[,1:2])
  } else {
    lonmin <- min(input_x)
    latmin <- min(input_y)
    ilon   <- function(x) 1+(x-lonmin)/dlon[1]
    ilat   <- function(x) 1+(x-latmin)/dlat[1]
    ix  <- sapply(output_xy[,1], FUN=ilon)
    iy  <- sapply(output_xy[,2], FUN=ilat)
    
    ### CHECK IF THIS IS WHAT WE WANT !!!! - EXTRAPOLATION
    iix <- pmin(length(input_x) , as.integer(ix))  # index to pt on right
    iix <- pmax(1 , iix)  # index to pt on left
    iiy <- pmin(length(input_y) , as.integer(iy))  # index to pt above
    iiy <- pmax(1, iiy)  # index to pt below

    dx  <- (output_xy[,1]-input_x[iix])  /dlon[1]
    dy  <- (output_xy[,2]-input_y[iiy] ) /dlat[1]
    dxvalue <- input_2D[,-1] - input_2D[,-ncol(input_2D)]
    dxvalue <- cbind(0, dxvalue, 0)
    dyvalue <- input_2D[-1,] - input_2D[-nrow(input_2D),]
    dyvalue <- rbind(0, dyvalue, 0)
    
    newz <- input_2D[cbind(iix, iiy)]+
                0.5*(dxvalue[cbind(iix, iiy)]*dx +
                     dyvalue[cbind(iix, iiy)]*dy) 
    row.names(output_xy) <- NULL
    Result <- data.frame(output_xy, v=newz)
    
  }
  Result 
}

# ==============================================================================
# ==============================================================================
# maps xy values on a 2D regular grid
# ==============================================================================
# ==============================================================================

interpolate_xy_2D <- function(input_xyz, 
                   output_x=NULL, output_y=NULL){
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_xyz[,1]))
  if (is.null(output_y)) 
    output_y <- sort(unique(input_xyz[,2]))
  
  xo  <- unique(output_x)
  xor <- range(xo)
  dxo <- diff(xor)/(length(xo)-1)
  
  yo <- unique(output_y)
  yor <- range(yo)
  dyo <- diff(yor)/(length(yo)-1)
  
  dxi <- min(diff(sort(unique(input_xyz[,1]))))
  dyi <- min(diff(sort(unique(input_xyz[,2]))))
  
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyz[,1]) < min(output_x))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  if (min(input_xyz[,1]) > max(output_x))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyz[,2]) < min(output_y))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  if (min(input_xyz[,2]) > max(output_y))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")

  # increase resolution
  if (dxi > dxo | dyi > dyo) {
#    Ux <- diff(range(diff(sort(unique(input_xyz[,1])))))
#    Uy <- diff(range(diff(sort(unique(input_xyz[,2])))))
#    if (Ux == 0 & Uy == 0) {  # uniform grid
       OO <- interpolate_xy_2D(input_xyz, output_x=NULL, output_y=NULL)
       return(interpolate_2D_2D(OO$x, OO$y, OO$v, output_x, output_y))  # RECURSIVE!!!??
#    }
#    else stop("cannot map: input resolution too low and nonuniform input grid")
  } else {
    
  # map input to output grid cells
  # use approx to estimate the ix and iy values
   app  <- approx(x = xo, y=1:length(xo), xout=input_xyz[,1], rule=2)
   ix   <- as.integer(app$y)
   app  <- approx(x = yo, y=1:length(yo), xout=input_xyz[,2], rule=2)
   iy   <- as.integer(app$y)
    
#  ix <- as.integer(1+(input_xyz[,1]-xor[1])/dxo)
  iy <- as.integer(1+(input_xyz[,2]-yor[1])/dyo)

  # No points to map in these output intervals.
  ix.lack <- which(!1:length(xo) %in% unique(ix))
  iy.lack <- which(!1:length(yo) %in% unique(iy))
  
  lx <- length(ix.lack)
  ly <- length(iy.lack)
  
  if (lx+ly) {
    
    if (!lx) ix.lack <- 1
    if (!ly) iy.lack <- 1
    
    ix.lack <- rep(ix.lack, length.out=lx+ly)
    iy.lack <- rep(iy.lack, length.out=lx+ly)
    ix      <- c(ix, ix.lack)
    iy      <- c(iy, iy.lack)
    
    # create empty vector with all lacking values
    toadd   <- input_xyz[1:(lx+ly), ]
    toadd[,] <- 0
    input_xyz <- rbind(input_xyz, toadd)
  }
  
  # sum all variables within a grid cell, and count the number of elements
  zsum <- tapply(input_xyz[,3], 
                 INDEX = list(ix, iy), 
                 FUN   = sum)
  zlen <- tapply(input_xyz[,3], 
                 INDEX = list(ix, iy), 
                 FUN   = length)
  
  # create x and y values - select only values within ranges
#  x  <- xor[1] + as.integer(rownames(zsum))* dxo
#  ii <- which(x >= xor[1] & x<= xor[2])
# x <- x[ii]  
  
#  y  <- yor[1] + as.integer(colnames(zsum))* dyo
#  jj <- which(y >= yor[1] & y<= yor[2])
#  y <- y[jj]

#  val <- zsum[ii,jj] / zlen[ii,jj]
  
  x <- xo[as.integer(rownames(zsum))]
  y <- yo[as.integer(colnames(zsum))]
  val <- zsum/zlen
  out <- list(x, y, val)
  names(out) <- colnames(input_xyz)[1:3]
  return(out)
  }
}

# ==============================================================================
# ==============================================================================
# maps gridded values on a regular grid
# ==============================================================================
# ==============================================================================

interpolate_2D_2Dold <- function(input_x, input_y, input_2D, 
                              output_x=NULL, output_y=NULL){
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_x))
  if (is.null(output_y)) 
    output_y <- sort(unique(input_y))
  
 # output grid size
  xo  <- unique(output_x)
  xor <- range(xo)

  yo <- unique(output_y)
  yor <- range(yo)

  xi  <- unique(input_x)
  xir <- range(xi)
  
  yi <- unique(input_y)
  yir <- range(yi)
  
  if (diff(range(diff(xi))) != 0 | diff(range(diff(yi))) != 0)
    stop("cannot interpolate - grid sizes are not constant")
  
  # input grid size
  dxi <- min(diff(sort(unique(input_x))))
  dyi <- min(diff(sort(unique(input_y))))
  
  # map output to input grid cells
  ix <- pmin(length(input_x), as.integer(1+(output_x-xir[1])/dxi))
  iy <- pmin(length(input_y), as.integer(1+(output_y-yir[1])/dyi))

  ix <- pmax(1, ix)
  iy <- pmax(1, iy)
  
  in2D <- cbind(input_2D, 
                input_2D[,ncol(input_2D)])
  in2D <- rbind(in2D,     
                in2D    [nrow(in2D),])
  
  IG   <- as.matrix(expand.grid(ix, iy))
  IGx1 <- t(t(IG) + c(1,0))
  OX   <- rep(output_x, times=length(output_y))
  vx   <- in2D[IG] + (in2D[IGx1]-in2D[IG])*(OX-input_x[IG[,1]])/dxi  
  vx   <- matrix(nrow=length(output_x), ncol=length(output_y), data=vx)
#  ZX<<- zx
  
  IG   <- as.matrix(expand.grid(iy, ix))[, c(2,1)]
  IGy1 <- t(t(IG) + c(0,1))
  OY   <- rep(output_y, times=length(output_x))
  vy   <- in2D[IG] + (in2D[IGy1]-in2D[IG])*(OY-input_y[IG[,2]])/dyi  
  vy   <- matrix(nrow=length(output_x), ncol=length(output_y), data=vy, byrow=TRUE)
  vv   <- 0.5*(vx+vy)
#  ZY<<- zy
  
  out <- list(x=output_x, y=output_y, v=vv)
  out
}

# ==============================================================================
# ==============================================================================

interpolate_2D_2D <- function(input_x, input_y, input_2D, 
                              output_x=NULL, output_y=NULL){
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_x))
  if (is.null(output_y)) 
    output_y <- sort(unique(input_y))
  
  # output grid size
  xo  <- unique(output_x)
  yo  <- unique(output_y)

  # interpolate in x direction  
  # yn = y[ix] + fac*(y[ix+1]-y[ix])

  # use approx to estimate the ix and fac values
  app  <- approx(x = input_x, y=1:length(input_x), xout=xo, rule=2)

  ix   <- as.integer(app$y)
  fac  <- app$y- ix
  ixp1 <- pmin(ix+1, max(ix, na.rm=TRUE))
  no <- length(fac) # number of output values

  OO <- sapply(1:no, 
               FUN=function(x) input_2D[ix[x], ]+    
             (input_2D[ixp1[x], ]- input_2D[ix[x], ])*fac[x])
  
  input_2D <- t(OO)
  
  # interpolate in y direction
  app  <- approx(x = input_y, y=1:length(input_y), xout=yo, rule=2)
  ix   <- as.integer(app$y)
  fac  <- app$y- ix
  ixp1 <- pmin(ix+1, max(ix, na.rm=TRUE))
  no   <- length(fac) # number of output values
  
  output_2D <- sapply(1:no, 
                FUN=function(x) input_2D[, ix[x]] +    
                 (input_2D[, ixp1[x]]-input_2D[, ix[x]])*fac[x])

  out <- list(x=output_x, y=output_y, v=output_2D)
  out
}

# ==============================================================================

  map_tx <-  function(input_txv, # time (t), position (x), value
                      input_t, input_x, input_2D,
                      output_tx,
                      output_t, output_x){

    # convert time to numeric...
    if (!missing (input_txv)){
      
    if (! is.numeric(input_txv[,1]) )
      input_txv[,1] <- as.numeric(input_txv[,1])
    } 
    
    if (!missing (input_t)){
      if (! is.numeric(input_t))
        input_t <- as.numeric(input_t)
    }        
    if (!missing (output_t)){
      ot <- output_t
      if (!is.numeric(output_t))
        output_t <- as.numeric(output_t)
    }        
    if (!missing (output_tx)){
      ot <- output_tx[,1]
      if (!is.numeric(output_tx[,1]))
        output_tx[,1] <- as.numeric(output_tx[,1])
    }        
    
    RES <- map_xy (input_txv, # latitude (x), longitude (y), value
                   input_t, input_x, input_2D,
                   output_tx, 
                   output_t, output_x)
    
    if (!missing (output_t)){
      RES <- list(t=ot, x=RES$y, v=RES$v)
    }        
    
    if (!missing (output_tx)){
      RES[,1] <- ot
    }        
    
    RES
    
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
#                     output_t=unique(input$datetime),
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
#             output_t=output_t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
    
#  } else if (what == "temperature") {
#    if (is.null(input)) input <- TempHR
#    pickdata(longitude=longitude, latitude=latitude, 
#             output_t=output_t, input=input, 
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
#             output_t=output_t, input=input, 
#             shape=shape, plotit=plotit, n=n, verbose=verbose, 
#             what=what, ...)
#  } else {
#    if (is.null(input)) stop (" 'input' should be provided for variable ", what)
#    pickdata(longitude=longitude, latitude=latitude, 
#             output_t=output_t, input=input, 
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
                     output_t=seq(as.POSIXct(from),     # default = hourly
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
  
  output_t <-  output_t [order(output_t)]
  input    <- input[order(input$datetime), ]
  
  # Make a selection of the input based on times
  if (is.numeric(output_t)) # index to the datetime
    input <- input[output_t,]
  else {
#    ii <- which(input$datetime >= min(output_t) &     
#                input$datetime <= max(output_t))
#    
#    if (! length(ii))
#      stop ("cannot interpolate: no overlap in time between input and output_t")
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
                        xout = output_t, 
                              rule = 2)$y)     
  
  # take weighted average of input data: multiply all columns (margin=2)
  # with correct weighing factor (wu) and sum over the rows
  res <- data.frame(datetime = output_t, 
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
  attributes(res)$ps   <- wu
  
  
  res
}

# pick_xyt(plotit=TRUE)
# A <- pick_xyt(plotit=TRUE, input=TempHR, output_t=1:40)
# A <- pick_xyt(plotit=TRUE, input=TempHR, output_t=seq(from=as.POSIXct("2021-03-01"),to= as.POSIXct("2021-05-01"), by=60))
# A <- pick_xyt(plotit=TRUE, input=TempHR)
# A <- pick_xyt(plotit=TRUE, input=TempLR, output_t = seq("15-01-2021", "15-12-2021", by=3600*24))
#  A <- pick_xyt(plotit=TRUE, input=TempHR, output_t=seq(from=as.POSIXct("2021-03-01"),to= as.POSIXct("2021-03-31"), by=60), verbose=TRUE, ylab="m", latitude=53.15)

#  A <- pick_xyt(plotit=TRUE, what="NH4", input=RWSbiogeo2021, verbose=TRUE, ylab="m", latitude=53.15)

