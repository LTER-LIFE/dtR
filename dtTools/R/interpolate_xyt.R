# ==============================================================================
# distance between two sets of points (input_xy, and output_xy)
# ==============================================================================

getDistance <- function(input_xy, 
                        output_xy, 
                        asp = 1) {   # y/x aspect ratio (e.g. longitude)
  
  DD1 <- outer(output_xy [, 1], input_xy[, 1], FUN = "-")
  DD2 <- outer(output_xy [, 2], input_xy[, 2], FUN = "-")
  
  Distance <- sqrt(DD1^2 + DD2^2*asp^2)  # take into account aspect ratio
  
  rm(list=c("DD1", "DD2"))
  Distance
}

# ==============================================================================
# finds weights based on inverse distances
# ==============================================================================

findWeights.id <- function(input_xy, 
                           output_xy, 
                           asp   = 1, 
                           nmean = 3){  # number of points

  # distance between input and output points
  Distance <- getDistance(input_xy, output_xy, asp = asp)
  
  # find nmean closest points 
  imin <- min(nmean, nrow(input_xy))  # at most nmean columns selected
  
  i.close <- t(apply(X      = Distance, 
                     MARGIN = 1, 
                     FUN    = order))[, 1:imin]
  
  if (!is.matrix(i.close)) 
    i.close <- matrix(nrow = 1, data = i.close)
  
  # weighing values 
  w.close <- NULL
  
  ii <- 1:nrow(i.close)
  
  # estimate distances
  for (i in 1:ncol(i.close))
    w.close <- cbind(w.close, 
                     Distance[cbind(ii, i.close[,i])])
  
  # keep track of points that are exactly matched
  ii      <- which(w.close == 0)
  
  # inverse distance
  w.close <- 1/w.close   
  
  # rescale so that sum=1
  w.close   <- w.close / rowSums(w.close, na.rm = TRUE)
  
  if (length(ii)) 
    w.close[ii] <- 1
  
  w.close[is.nan(w.close)] <- 0
  
  # bc is a logical denoting barycentric interpolation (FALSE here)
  bc  <- rep(FALSE, times = nrow(output_xy))
  
  return(list(idx = i.close, 
              p   = w.close, 
              bc  = bc))
}

# ==============================================================================
# finds barycentric weights - based on delaunay triangulation
# ==============================================================================

findWeights.bc <- function(input_xy, 
                           output_xy, 
                           asp = 1, 
                           dn  = NULL){  # delaunay triangulation
  
  # perform delaunay triangulation
  if (is.null(dn))
    dn  <- delaunayn(input_xy)
  
  # find interpolation weights and indices
  tri <- tsearch(input_xy[ ,1],
                 input_xy[ ,2],
                 dn, 
                 as.vector(output_xy[,1]),
                 as.vector(output_xy[,2]),
                 bary = TRUE)
  
  # 3 points : indices
  idx <- dn[tri$idx, ]
  
  if (!is.matrix(idx)) 
    idx <- matrix(nrow = 1, data = idx)
  
  # corresponding weights
  p   <- tri$p
  
  # bc is a logical denoting barycentric interpolation
  bc  <- rep(TRUE, times = nrow(output_xy))
  
  # points NOT embraced by triangles
  isNA  <- which(is.na(idx[, 1]))

  if (length(isNA)) {
    # use inverse distance for these points
    WW <- findWeights.id(input_xy, output_xy[isNA, ], 
                         asp = asp, nmean = 3)
    
    idx[isNA, ] <- WW$idx
    p  [isNA, ] <- WW$p
    bc [isNA  ] <- FALSE
   
  }
  return(list(idx = idx, p = p, bc = bc))
}

### https://dahtah.wordpress.com/2013/03/06/barycentric-interpolation-fast-interpolation-on-arbitrary-grids/

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

CreateA <- function(X,   # 3 x-values 
                    Y)   # 3 y-values 
  matrix(nrow  = 3, 
         byrow = TRUE, 
         data  = c(X, 
                   Y, 
                   1, 1, 1))


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
  
  # remove duplicates where TT=1   
  Triangles <- Triangles[TT == 0, ]                    
  
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
  return(Triangles[order(Triangles[,3], 
                         Triangles[,2], 
                         Triangles[,1]), ])
}

# ==============================================================================
# finds barycentric weights - triangulation based on per point closest triangle
# ==============================================================================

findWeights.bc.2 <- function(input_xy, 
                             output_xy, 
                             asp  = 1, 
                             nmax = 10){  # max number of triangles
  
  nmax <- min(nmax, nrow(input_xy))
  NS   <- CreateTriangles(n = nmax)
  
  output_xy.B <- cbind(output_xy, 1)
  
  # distance between data points and all points to map -
  # for latitude and longitude: apply the aspect ratio
  Dist <- getDistance(input_xy, 
                      output_xy, 
                      asp = asp)
  
  # for each point to map, sort datapoints according to closeness with data points
  i.close <- t(apply(X      = Dist, 
                     MARGIN = 1, 
                     FUN    = order))
  
  # Apply barycentric interpolation
  
  # Results will be stored here: index to data points and the factors
  IRES <- FAC <- matrix(nrow = nrow(output_xy), ncol = 3, NA)
  
  # All the points that still need to be mapped
  isNA     <- 1:nrow(output_xy)
  
  for (TT in 1: nrow(NS)){ # loop over possible triangles; closest points first
    
    isel <- NS[TT, ]   # points to select
    
    A    <- NULL       # for all points to be mapped: the interpolation factors
    
    for (i in isNA)      # for all indices still in isNA : solve the equations
      A <- rbind(A, 
                 solve(CreateA(input_xy[i.close[i, isel], 1], 
                               input_xy[i.close[i, isel], 2]), 
                       b = output_xy.B[i,]))
    
    # select results that show that data points embrace the b-points
    ii <- apply(X      = A, 
                MARGIN = 1, 
                FUN    = function(x) all(x >=0))
    
    if (sum(ii) > 0) { 
      IRES[isNA[ii], ] <- i.close[isNA[ii], isel]
      FAC [isNA[ii], ] <- A[ii, ]
      
      isNA <- isNA[!ii]    # remove those from the points to be mapped
    }
    
    if (length(isNA) == 0) break()
  }
  
  Embrace <- rep(TRUE, times = nrow(output_xy))

  if (length(isNA)){
  
  # not all have been mapped - for those that have nout embracing points, 
  #  we interpolate with the inverse distances
    
    Embrace[isNA] <- FALSE
    IRES[isNA, ]  <- i.close[isNA, 1:3]    # 3 closest points
    ff   <- 1/cbind(Dist[cbind(isNA, IRES[isNA, 1])],
                    Dist[cbind(isNA, IRES[isNA, 2])],
                    Dist[cbind(isNA, IRES[isNA, 3])])  
    ii   <- which(ff == 0)
    ff   <- ff/rowSums(ff, na.rm=TRUE) 
    
    if (length(ii)) 
      ff[ii] <- 1
    
    ff[is.nan(ff)] <- 0
    
    FAC [isNA, ] <- ff
  }
  
  return(list(idx = IRES, p = FAC, bc = Embrace))
  
}

# ==============================================================================
# ==============================================================================
# Interpolate spatial time series defined in (xy) points to other (xy) points 
# Typically, input_xy contains few points 
# ==============================================================================
# ==============================================================================

interpolate_sparse <-  function(input_xytv,  # longitude (x), latitude (y), time, value
                             output_xy, 
                             output_t, 
                             barycentric = FALSE,  # can also be 1 or 2
                             nmean      = 3,
                             geographic = TRUE,   # longitude and latitude
                             rule       = 2,      # how to use approx
                             ID         = NULL,   # unique identifier for the output
                             all.output = FALSE)  {   
  
  ##### step 0: prepare/check the inputs #####
  # attributes of input saved
  attrs <- NULL
  
  if (inherits(input_xytv, "dtLife"))
    attrs <- attributes(input_xytv)
  
  # check input_xytv: 4 columns: longitude (x), latitude (y), time, value
  if (ncol(input_xytv) != 4) 
    stop ("'input_xytv' should have 4 columns: longitude (x), latitude (y), time, value")

  # name of the time column
  tname <- names(input_xytv)[3] # to label the output
  if (is.null(tname)) tname <- "datetime"
  
  # input xy values: unique set 0f first two columns
  input_xy  <- unique(input_xytv[, 1:2])

  # output xy values: 
  if (is.vector(output_xy))
    output_xy <- matrix(nrow = 1, ncol = 2, data = output_xy)
  else 
    output_xy <- unique(output_xy[,1:2])
  
  nout <- nrow(output_xy)
  
  # euclidean_distance aspect ratio
  if (geographic) {  # latitude/longitude coordinates
    
    if (min(output_xy[,2], na.rm=TRUE) < -90 |
        max(output_xy[,2], na.rm=TRUE) >  90)
      stop ("second column of output_xytv should contain valid latitude, -90:90")
    
    if (min(input_xy[,2], na.rm=TRUE) < -90 |
        max(input_xy[,2], na.rm=TRUE) >  90)
      stop ("second column of input_xy should contain valid latitude, -90:90")
    asp      <- aspectratio(output_xy[,2])  # y/x aspect ratio
    
  } else asp = 1
  
  ##### step 1: weighting points #####
  # distance between two points
  
  if (barycentric == 1)  # barycentric; triangulation fixed
    W <- findWeights.bc(input_xy, output_xy, asp = asp)
  
  else if (barycentric == 2)  # barycentric; triangulation smallest
    W <- findWeights.bc.2(input_xy, output_xy, asp = asp)
  
  else if (barycentric == 0)   # inverse distance
    W <- findWeights.id(input_xy, output_xy, asp = asp, nmean = nmean)   
  
  # columns that are used in weighing
  i.unique <- unique(as.vector(W$idx))
  
  if (barycentric != 0)
    nmean <- min(3, nrow(input_xy))      # at most 3 columns selected
  
  else      
    nmean <- min(nmean, nrow(input_xy))  # at most nmean columns selected

  ##### step2: interpolate to output times (output_t) for all inputs #####
  xytin <- matrix(nrow = length(output_t), 
                  ncol = max(i.unique), 
                  data = 0)
  
  for (i in i.unique){
    # interpolate to required times
    ii <- which(input_xytv[,1] == input_xy[i, 1] & 
                input_xytv[,2] == input_xy[i, 2])
    V.out  <- approx(x    = input_xytv[ii, 3],   # time
                     y    = input_xytv[ii, 4],   # value
                     xout = output_t, 
                     rule = rule)$y
    xytin[,i] <- V.out
  }
  
  
  # names of the output
  if (is.null(ID))
    stnames <- paste("st", 1:nout, sep="")
  else {
    stnames <- ID
  
    if (length(stnames) != nout) 
      stop("'ID' should be of length = output")
  }
  
  Stations <- data.frame(ID = stnames,
                         x  = output_xy[, 1], 
                         y  = output_xy[, 2])
  row.names(Stations) <- NULL
  
  ##### step 3: create matrix equation and solve
 
  if (length(i.unique) == 1){  
  # only one input dataset - needs to be repeated for all outputs
    result = matrix(ncol = length(output_t), 
                    data = rep(xytin, each = nout))    
  } else { 
  # sparse matrix algebra to solve in one go  
   f      <- t(xytin)
   M      <- sparseMatrix(i    = rep(1:nout, each = nmean),
                          j    = as.numeric(t(W$idx)),
                          x    = as.numeric(t(W$p)),
                          dims = c(nout, nrow(f)))
   result <- M %*% as.matrix(f)
  }  
  
  # colnames 
  c.names <- c("longitude", "latitude", as.character(output_t))
  result <- cbind(output_xy, 
                  as.matrix(result))      # STAYS A MATRIX
                       
  colnames(result) <- c.names # c(tname,  stnames)
  class(result) <- c("dtLife", class(result))

  atout <- attributes(result)
  
  if (! is.null(attrs))
    attributes(result) <- c(atout, attrs[!names(attrs) %in% c(names(atout), "class")])
  
  attr(result, "datetime_type") <- class(output_t)
  attr(result, "datetime") <- output_t
  attributes(result)$format <- "wide"
  
  if (! is.null(ID)) attr(result, "ID") <- Stations
  
  names(Stations)[1] <- "stations"
  attr(result, "stations") <- Stations
  
  attr(result, "processing") <- c(
    atout$processing, paste("interpolated from 2D-time input, at:", Sys.time()))

  if (all.output){
    
    for (i in 1:nout){
      attributes(result)$interpolation[[i]] <- data.frame(
        x           = input_xy[W$idx[i,],1],
        y           = input_xy[W$idx[i,],2],
        stat.index  = W$idx[i,],
        stat.weight = W$p[i,]
      ) 
   }
  names(attributes(result)$interpolation) <- stnames
  if (! is.null(ID)) attr(result, "ID") <- Stations
  
  attributes(result)$barycentric <- W$bc
  if (! is.null(ID)) attr(result, "ID") <- Stations
  }
  result
}


# ==============================================================================

match_timeseries <- function(input_xytv, # timeseries in latitude, longitude, time, value
                             data       # should have ID with (ID, latitude, longitude) in its attributes
                             ) {  
  # data on which to map the timeseries
  dname    <- deparse(substitute(data))
  tsname   <- deparse(substitute(timeseries))

  timeseries <- as.data.frame(input_xytv)
  
  #  if (is.character(timeseries[,dtname]))
  #    timeseries[,dtname] <- as.Date(timeseries[,dtname])
  
  output_xy <- meta(data)$ID  [, c("ID", "longitude", "latitude")]                                    
  output_t  <- unique(timeseries[, 3])
  
  matched   <- interpolate_sparse(timeseries, 
                                  output_xy = output_xy[,-1], 
                                  output_t  = output_t, 
                                  ID        = output_xy$ID)
  
  attributes(matched)$ID         <- output_xy
  attributes(matched)$variables  <- meta(input_xytv)$variables
  attributes(matched)$processing <- paste("matched timeseries named '", tsname, 
                                          "'with", dname, "at", Sys.time())
  matched
  
}

