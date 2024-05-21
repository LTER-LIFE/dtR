# ==============================================================================
# distance between two points
# ==============================================================================

getDistance <- function(input.xy, 
                        output.xy, 
                        asp = 1) {   # y/x aspect ratio
  
  DD1 <- outer(output.xy [, 1], input.xy[, 1], FUN="-")
  DD2 <- outer(output.xy [, 2], input.xy[, 2], FUN="-")
  
  Distance <- sqrt(DD1^2 + DD2^2*asp^2)
  
  rm(list=c("DD1", "DD2"))
  Distance
}

# ==============================================================================
# finds weights based on inverse distances
# ==============================================================================

findWeights.id <- function(input.xy, output.xy, asp=1, nmean=3){

  Distance <- getDistance(input.xy, output.xy, asp=asp)
  
  # find three closest points 
  # which data sets to use for each output station, 
  imin <- min(nmean, nrow(input.xy))  # at most 3 columns selected
  
  i.close <- matrix(nrow  = nrow(output.xy), 
                    byrow = TRUE, 
                    data  = unlist(apply(Distance, 
                                         MARGIN = 1, 
                                         FUN    = order, 
                                         simplify=FALSE)))[,1:imin]
  i.close <- t(apply(X      = Distance, 
                   MARGIN = 1, 
                   FUN    = order))[,1:imin]
  
  # weighing values ~ inverse distance
  w.close <- NULL
  
  for (i in 1:ncol(i.close))
    w.close <- cbind(w.close, 
                     Distance[cbind(1:nrow(i.close), i.close[,i])])
  
  ii      <- which(w.close == 0)
  w.close <- 1/w.close
  
  # rescale so that sum=1
  w.close   <- w.close/rowSums(w.close, na.rm=TRUE)
  
  if (length(ii)) 
    w.close[ii] <- 1
  w.close[is.nan(w.close)] <- 0
  
  bc  <- rep(FALSE, times=nrow(output.xy))
  
  return(list(idx = i.close, p = w.close, bc=bc))
}

# ==============================================================================
# finds barycentric weights
# ==============================================================================

findWeights.bc <- function(input.xy, output.xy, asp=1, dn=NULL){
  
  if (is.null(dn))
    dn  <- delaunayn(input.xy)
  
  tri <- tsearch(input.xy[,1],
                 input.xy[,2],
                 dn, 
                 output.xy[,1],
                 output.xy[,2],
                 bary = TRUE)
  idx <- dn[tri$idx,]
  p   <- tri$p
  bc  <- rep(TRUE, times=nrow(output.xy))
  
  isNA  <- which(is.na(idx[, 1])) # points NOT embraced by triangles

  
  if (length(isNA)) {
    WW <- findWeights.id(input.xy, output.xy[isNA, ], asp=asp, nmean=3)
    
    idx[isNA, ] <- WW$idx
    p[isNA, ]   <- WW$p
    
    bc[isNA] <- FALSE
   
  }
  return(list(idx = idx, p = p, bc=bc))
}

### https://dahtah.wordpress.com/2013/03/06/barycentric-interpolation-fast-interpolation-on-arbitrary-grids/
#2D barycentric interpolation at points Xi for a function with values f measured at locations X
#For N-D interpolation simply replace tsearch with tsearchn and modify the sparse matrix definition to have non-zero values in the right spots.
#interp.barycentric <- function(X,f,Xi)
#{
#  require(geometry)
#  require(Matrix)
#  dn <- delaunayn(X)
#  tri <- tsearch(X[,1],X[,2],dn,Xi[,1],Xi[,2],bary=T)
  #For each line in Xi, defines which points in X contribute to the interpolation
 # active <- dn[tri$idx,]
  #Define the interpolation as a sparse matrix operation. Faster than using apply, probably slower than a C implementation
#  M <- sparseMatrix(i=rep(1:nrow(Xi),each=3),j=as.numeric(t(active)),x=as.numeric(t(tri$p)),dims=c(nrow(Xi),length(f)))
#  as.numeric(M%*%f)
  
#}

# ============================================================================
# A function to create the barycentric matrix
# ============================================================================
# Q: given 3 points with coordinated (x1, y1), (x2, y2), (x3, y3)
# and a point (x, y)
#    what are the barycentric coefficients (p1, p2, p3) so that x
#    is a linear combination of the 3 points
# if p1, p2, p3 are > 0 then the point x is embraced by the 3 pts

# System of equations: 
# x1*p1 + x2*p2 + x3*p3 = x
# y1*p1 + y2*p2 + y3*p3 = y
# 1* p1 + 1 *p2 + 1 *p3 = 1

CreateA <- function(X,   # 3 x-values of 
                    Y)   # 3 y-values of points
  matrix(nrow  = 3, 
         byrow = TRUE, 
         data  = c(X, Y, 1, 1, 1))


# ============================================================================
# A function to find all different triangles of n points
# ============================================================================

CreateTriangles <- function(n = 5) {   
  
  # all possible combinations, includes double combinations
  Triangles <- expand.grid(1:n, 2:n, 3:n)
  
  # check for same points in each set of 3
  TT        <- apply(Triangles, 
                     MARGIN = 1, 
                     FUN    = anyDuplicated)  
  
  Triangles <- Triangles[TT==0, ]       # remove duplicates where TT=1                    
  
  # check for same rows - first order the points 
  Triangles <- t(apply(Triangles, 
                       MARGIN = 1, 
                       FUN    = function(x) x[order(x)]))
  
  # make a string for each row, 3 points collated, and find unique strings
  TT        <- unique(apply(Triangles, 
                            MARGIN = 1, 
                            FUN    = function(x) paste(x, collapse="_")))
  # split it again
  Triangles <- matrix(ncol  = 3, 
                      data  = as.integer(unlist(strsplit(TT, "_"))), 
                      byrow = TRUE)
  
  # order it by 1, 2, 3 (this assumes that the point 1 is closer than 2 than 3)
  return(Triangles[order(Triangles[,3], Triangles[,2], Triangles[,1]), ])
}


# ============================================================================
# Function to apply barycentric interpolation
# Typically, input.xy contains few points 
# ============================================================================

findWeights.bc.2 <- function(input.xy, output.xy, asp=1, nmax=10){
  
  nmax <- min(nmax, nrow(input.xy))
  NS <- CreateTriangles(n = nmax)
  
  output.xy.B <- cbind(output.xy, 1)
  
  # distance between data points and all points to map -
  # note: for latitude and longitude: apply the aspect ratio
  Dist <- getDistance(input.xy, output.xy, asp=asp)
  
  # for each point to map, sort datapoints according to closeness with data points
  i.close <- t(apply(X      = Dist, 
                     MARGIN = 1, 
                     FUN    = order))
  
  # Apply barycentric interpolation
  
  # Results will be stored here: index to data points and the factors
  IRES <- FAC <- matrix(nrow = nrow(output.xy), ncol=3, NA)
  
  # All the points that still need to be mapped
  iA     <- 1:nrow(output.xy)
  
  for (TT in 1: nrow(NS)){   # loop over possible triangles; closest points first
    
    isel <- NS[TT, ]   # points to select
    
    A    <- NULL       # for all points to be mapped: the interpolation factors
    
    for (i in iA)      # for all indices still in iA : solve the equations
      A <- rbind(A, 
                 solve(CreateA(input.xy[i.close[i, isel],1], 
                               input.xy[i.close[i, isel],2]), 
                       b = output.xy.B[i,]))
    
    # select results that show that data points embrace the b-points
    ii <- apply(X      = A, 
                MARGIN = 1, 
                FUN    = function(x) all(x >=0))
    
    if (sum(ii) > 0) { 
      IRES[iA[ii], ] <- i.close[iA[ii], isel]
      FAC [iA[ii], ] <- A[ii, ]
      
      iA <- iA[!ii]    # remove those from the points to be mapped
    }
    
    if (length(iA) == 0) break()
    # plot(c(A[ii,] %*% c(x1, x2, x3)), c(A[ii,]%*% c(y1, y2, y3)))
  }
  Embrace <- rep(TRUE, times = nrow(output.xy))
  
  # not all have been mapped - for those that have no embracing points, 
  #  we interpolate with the inverse distances
  
  if (length(iA)){
    
    Embrace[iA] <- FALSE
    IRES[iA, ]  <- i.close[iA, 1:3]    # 3 closest points
    ff   <- 1/cbind(Dist[cbind(iA, IRES[iA, 1])],
                    Dist[cbind(iA, IRES[iA, 2])],
                    Dist[cbind(iA, IRES[iA, 3])])  
    ii   <- which(ff == 0)
    ff   <- ff/rowSums(ff, na.rm=TRUE) 
    if (length(ii)) ff[ii] <- 1
    ff[is.nan(ff)] <- 0
    
    FAC [iA, ] <- ff
  }
  
  return(list(idx = IRES, p=FAC, bc=Embrace))
  
}


# ==============================================================================
# ==============================================================================
# Interpolate spatial time series defined in (xy) points to other (xy) points 
# ==============================================================================
# ==============================================================================

interpolate_xyt <-  function(input.xytv,  # longitude (x), latitude (y), time, value
                             output.xy, 
                             output.t, 
                             barycentric = FALSE,
                             nmean      = 3,
                             ID         = NULL,   #unique identifier for the output
                             geographic = TRUE)   # longitude and latitude
                              {   
  
  attrs <- NULL
  
  if (inherits(input.xytv, "dtLife"))
    attrs <- attributes(input.xytv)
  
  # input.xytv: 4 columns: longitude (x), latitude (y), time, value
  if (ncol(input.xytv) != 4) 
    stop ("'input.xytv' should have 4 columns: longitude (x), latitude (y), time, value")

  tname <- names(input.xytv)[3] # to label the output
  
  # input xy values
  input.xy  <- unique(input.xytv[, 1:2])
  ni        <- nrow(input.xy)
  
  # output xy values
  if (is.vector(output.xy))
    output.xy <- matrix(nrow = 1, ncol = 2, data = output.xy)
  else 
    output.xy <- unique(output.xy[,1:2])
  
  no        <- nrow(output.xy)
  
  ##### step 1: weighting points #####
  # distance between two points
  # euclidean_distance
  
  if (geographic) {
    if (min(output.xy[,2], na.rm=TRUE) < -90 |
        max(output.xy[,2], na.rm=TRUE) >  90)
      stop ("second column of output.xytv should contain valid latitude, -90:90")
  
    if (min(input.xy[,2], na.rm=TRUE) < -90 |
        max(input.xy[,2], na.rm=TRUE) >  90)
      stop ("second column of input.xy should contain valid latitude, -90:90")
    asp      <- 1/cos((mean(output.xy[,2], na.rm=TRUE)*pi)/180)  # y/x aspect ratio
  } else asp = 1
  
  if (barycentric == 1)  # barycentric; triangulation fixed
    W <- findWeights.bc(input.xy, output.xy, asp = asp)
  
  else if (barycentric == 2)  # barycentric; triangulation smallest
    W <- findWeights.bc.2(input.xy, output.xy, asp = asp)
  
  else if (barycentric == 0)   # inverse distance
    W <- findWeights.id(input.xy, output.xy, asp = asp, nmean = nmean)   
  
  # columns that are used in weighing
  i.unique <- unique(as.vector(W$idx))
  
  ##### step2: interpolate to output times for all inputs #####
  xytin <- matrix(nrow = length(output.t), 
                  ncol = max(i.unique), 
                  data = 0)
  
  for (i in i.unique){
    # interpolate to required times
    ii <- which(input.xytv[,1] == input.xy[i, 1] & 
                input.xytv[,2] == input.xy[i, 2])
    V.out  <- approx(x    = input.xytv[ii, 3],   # time
                     y    = input.xytv[ii, 4],   # value
                     xout = output.t, rule=2)$y
    xytin[,i] <- V.out
  }
  
  if (is.null(ID))
    stnames <- paste("st", 1:no, sep="")
  else {
    stnames <- ID
  
    if (length(stnames) != no) 
      stop("'ID' should be of length = output")
  }
  
  Stations <- data.frame(ID = stnames,
                         x  = output.xy[,1], 
                         y  = output.xy[,2])
  row.names(Stations) <- NULL
  
  ##### step 3: create matrix equation and solve
 
  lNA    <- no

  f      <- t(xytin)
  M      <- sparseMatrix(i  = rep(1:lNA, each = nmean),
                         j    = as.numeric(t(W$idx)),
                         x    = as.numeric(t(W$p)),
                         dims = c(lNA, nrow(f)))
  result  <- M %*% as.matrix(f)
  result <- data.frame(datetime = output.t, 
                       t(as.matrix(result)))
  colnames(result) <- c(tname,  stnames)

  atout <- attributes(result)
  
  if (! is.null(attrs))
    attributes(result) <- c(atout, attrs[!names(attrs) %in% names(atout)])
  
  attributes(result)$format <- "wide"
  if (! is.null(ID)) attr(result, "ID") <- Stations
  
  names(Stations)[1] <- "stations"
  attr(result, "stations") <- Stations
  
  attr(result, "processing") <- c(
    atout$processing, paste("interpolated from 2D-time input, at:", Sys.time()))

  attributes(result)$barycentric <- W$bc
  if (! is.null(ID)) attr(result, "ID") <- Stations
  
  result
}

# ==============================================================================

match_timeseries <- function(input.xytv, # timeseries in latitude, longitude, time, value
                             data       # should have ID with (ID, latitude, longitude) in its attributes
                             ) {  
  # data on which to map the timeseries
  dname    <- deparse(substitute(data))
  tsname   <- deparse(substitute(timeseries))

  timeseries <- as.data.frame(input.xytv)
  
  #  if (is.character(timeseries[,dtname]))
  #    timeseries[,dtname] <- as.Date(timeseries[,dtname])
  
  output.xy <- meta(data)$ID  [, c("ID", "longitude", "latitude")]                                    
  output.t  <- unique(timeseries[, 3])
  
  matched   <- interpolate_xyt(timeseries, 
                               output.xy = output.xy[,-1], 
                               output.t  = output.t, 
                               ID        = output.xy$ID)
  
  attributes(matched)$ID         <- output.xy
  attributes(matched)$variables  <- meta(input.xytv)$variables
  attributes(matched)$processing <- paste("matched timeseries named '", tsname, 
                                          "'with", dname, "at", Sys.time())
  matched
  
}

