# ==============================================================================
# ==============================================================================
# Maps xy points to other xy points - to be used for sparse input datasets
# ==============================================================================
# ==============================================================================

map_sparse <- function(input_xyv,  # longitude (x), latitude (y), value
                       output_xy, 
                       output_x, output_y, 
                       barycentric = FALSE,
                       nmean      = 3,
                       geographic = TRUE,   # longitude and latitude
                       rule = 2,
                       all.output = FALSE)  {   
  
  # check input
  if (missing (output_xy)){
    if (missing(output_x) | missing(output_y) )
      stop("either output_xy should be input or (output_x & output_y)")
    output_xy <- expand.grid(output_x, output_y)
  }
  
  attrs <- NULL
  
  if (inherits(input_xyv, "dtLife"))
    attrs <- attributes(input_xyv)
  
  # input_xyv: 3 columns: longitude (x), latitude (y), value
  if (ncol(input_xyv) != 3) 
    stop ("'input_xyv' should have 3 columns: longitude (x), latitude (y), value")
  
  # input xy values
  input_xy  <- unique(input_xyv[, 1:2])
  
  # output xy values
  if (is.vector(output_xy))
    output_xy <- matrix(nrow = 1, ncol = 2, data = output_xy)
  else 
    output_xy <- unique(output_xy[,1:2])
  
  no        <- nrow(output_xy)
  
  ##### step 1: weighting points #####
  # distance between two points
  # euclidean_distance
  
  if (geographic) {
    if (min(output_xy[,2], na.rm=TRUE) < -90 |
        max(output_xy[,2], na.rm=TRUE) >  90)
      stop ("second column of output_xytv should contain valid latitude, -90:90")
    
    if (min(input_xy[,2], na.rm=TRUE) < -90 |
        max(input_xy[,2], na.rm=TRUE) >  90)
      stop ("second column of input_xy should contain valid latitude, -90:90")
    asp      <- 1/cos((mean(output_xy[,2], na.rm=TRUE)*pi)/180)  # y/x aspect ratio
  } else asp = 1
  
  if (barycentric == 1)  # barycentric; triangulation fixed
    W <- findWeights.bc(input_xy, output_xy, asp = asp)
  
  else if (barycentric == 2)  # barycentric; triangulation smallest
    W <- findWeights.bc.2(input_xy, output_xy, asp = asp)
  
  else if (barycentric == 0)   # inverse distance
    W <- findWeights.id(input_xy, output_xy, asp = asp, nmean = nmean)   
  
  # columns that are used in weighing
  i.unique <- unique(as.vector(W$idx))
  
  if (barycentric != 0)
    nmean <- min(3, nrow(input_xy))  # at most 3 columns selected
  else      
    nmean <- min(nmean, nrow(input_xy))  # at most 3 columns selected
  
  result <- output_xy
  
  ##### step 3: create matrix equation and solve
  
  lNA    <- no
  
  f      <- input_xyv[,3]
  M      <- sparseMatrix(i  = rep(1:lNA, each = nmean),
                         j    = as.numeric(t(W$idx)),
                         x    = as.numeric(t(W$p)),
                         dims = c(lNA, length(f)))
  
  result$value  <- as.vector(M %*% as.matrix(f))
  colnames(result) <- colnames(input_xyv)
  atout <- attributes(result)
  
  if (! is.null(attrs))
    attributes(result) <- c(atout, attrs[!names(attrs) %in% names(atout)])
  
  attributes(result)$format <- "long"
  attr(result, "processing") <- c(
    atout$processing, paste("interpolated from 2D input, at:", Sys.time()))
  
  if (all.output){
    for (i in 1:no){
      attributes(result)$interpolation[[i]] <- data.frame(
        x = input_xy[W$idx[i,],1],
        y = input_xy[W$idx[i,],2],
        stat.index  = W$idx[i,],
        stat.weight = W$p[i,]
      ) 
    }
    
    attributes(result)$barycentric <- W$bc
  }
  
  result
}


# ==============================================================================
# ==============================================================================
# Maps xy points to other xy points - large input datasets
# ==============================================================================
# ==============================================================================

map_dense <- function(input_xyv, # 3 columns: longitude (x), latitude (y), value (v)
                      input_x, input_y, input_2D,
                      output_xy, 
                      output_x, output_y,
                      geographic = TRUE){  # coordinates are longitude/latitude  ?
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
  if (! missing (input_xyv) & 
      ! missing (output_xy))
    
    result <- interpolate_xy_xy(input_xyv = input_xyv, # latitude (x), longitude (y), value
                                output_xy = output_xy,
                                geographic = geographic) 
  
  else if (! missing (input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing (output_xy))
    
    result <- interpolate_2D_xy (input_x   = input_x, 
                                 input_y   = input_y, 
                                 input_2D  = input_2D,        
                                 output_xy = output_xy,
                                 geographic = geographic)
  
  else if (! missing(input_xyv) & 
           ! missing(output_x)  & ! missing(output_y))
    
    result <- interpolate_xy_2D(input_xyz = input_xyv, 
                               output_x  = output_x, 
                               output_y  = output_y,
                               geographic = geographic)
  
  else if (! missing(input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing(output_x) & !missing(output_y))
    
    result <- interpolate_2D_2D (input_x  = input_x, 
                                 input_y  = input_y, 
                                 input_2D = input_2D, 
                                 output_x = output_x, 
                                 output_y = output_y,
                                 geographic = geographic)
  
  else 
    stop ("one of the required arguments is missing")
  
  attr(result, "processing") <- paste("remapped , at:", 
                                      Sys.time())
  
  return(result)
}

# ==============================================================================

interpolate_xy_xy <- function(input_xyv, # latitude (x), longitude (y), value
                              output_xy,
                              geographic){
  
  if (ncol(input_xyv) != 3) 
    stop ("'input_xyv' should have 3 columns: latitude (x), longitude (y), value")

  ##### step 1: weighting points #####
  # distance between two points 

  input_xy  <- input_xyv[, 1:2]
  ni        <- nrow(input_xy)
  if (ni == 0) stop("input_xyv is empty")
  
  output_xy <- unique(output_xy)
  coord.xy  <- output_xy[,1:2]
  
  if (geographic) asp <- aspectratio(coord.xy[,2]) else asp <- 1
  
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
  DD2 <- outer(output_xy [, 2], input_xy[, 2], FUN="-")*asp 
  Distance <- sqrt(DD1^2+DD2^2)
  rm(list=c("DD1", "DD2"))
  
  
  # which data sets to use for each output station, 
  imin <- min(3, ni)  # at most 3 columns selected
  
  # find imin closest points 
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
  
  # weighing values of the closest points ~ inverse distance
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
                        output_xy,
                        geographic){

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
  
  if (geographic) asp <- aspectratio(input_y) else asp <- 1
  
  dlon <- diff (input_x)      # diff(input_2D$longitude)
  dlat <- diff (input_y)*asp  # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)) == 0 &
                 diff(range(dlat)) == 0
  
  if (! equidistant){  # SLOW!
   input_2D.grid        <- expand.grid(input_x, input_y)
   names(input_2D.grid) <- c("longitude", "latitude")
   input_2D.grid$depth  <- as.vector(input_2D)
   input_2D.grid        <- na.omit(input_2D.grid)

   gg    <- na.omit(output_xy)
  
   Result <- map_dense(input_xyv = input_2D.grid, 
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
                              output_x=NULL, output_y=NULL,
                              geographic){
  
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
       return(interpolate_2D_2D(OO$x, OO$y, OO$v, output_x, output_y, 
                                geographic=geographic))  # RECURSIVE!!!??
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

interpolate_2D_2D <- function(input_x, input_y, input_2D, 
                              output_x=NULL, output_y=NULL, 
                              geographic){
  
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
  no   <- length(fac) # number of output values

  OO   <- sapply(1:no, 
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
                      output_t, output_x,
                      geographic = TRUE){

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
    
    RES <- map_dense (input_txv, # latitude (x), longitude (y), value
                      input_t, input_x, input_2D,
                      output_tx, 
                      output_t, output_x, geographic=geographic)
    
    if (!missing (output_t)){
      RES <- list(t=ot, x=RES$y, v=RES$v)
    }        
    
    if (!missing (output_tx)){
      RES[,1] <- ot
    }        
    
    RES
    
   ############## TO DO ############### 
 }
