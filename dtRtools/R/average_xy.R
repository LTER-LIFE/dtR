
# ==============================================================================
# ==============================================================================
# averages xy points from a dense grid to xy points at lower resolution
# ==============================================================================
# ==============================================================================

average_xy <- function(
                      input_xyv, # 3 columns: longitude (x), latitude (y), value (v)
                      input_x, input_y, input_2D,
                      output_x, output_y){  # coordinates are longitude/latitude  ?
  # check input
  if (missing (input_xyv)){
    if (missing(input_x) | missing(input_y) | missing(input_2D) )
      stop("either input_xyv should be input or (input_x, input_y & input_2D)")
  }
    
  # call correct function
  if (! missing(input_xyv) )
    
    result <- average_xy_2D(input_xyv = input_xyv, 
                            output_x   = output_x, 
                            output_y   = output_y)
  
  else if (! missing(input_x) & ! missing(input_y) & ! missing(input_2D))
    
    result <- average_2D_2D (input_x  = input_x, 
                             input_y  = input_y, 
                             input_2D = input_2D, 
                             output_x = output_x, 
                             output_y = output_y)
  
  else 
    stop ("one of the required arguments is missing")
  
  attr(result, "processing") <- paste("averaged in grid , at:", 
                                      Sys.time())
  
  return(result)
}

# ==============================================================================
# ==============================================================================
# maps xy values on a 2D regular grid
# ==============================================================================
# ==============================================================================

average_xy_2D <- function(input_xyv, 
                          output_x = NULL, 
                          output_y = NULL){
  
  ##### check input 

  input_xyv <- as.matrix(input_xyv[, 1:3])
  ni        <- nrow(input_xyv)
  if (ni == 0) stop("input_xyv is empty")
  if (any (is.na(input_xyv[, 1:2])))
    stop ("cannot proceed: some of the input coordinates are NA")
  
  output_x <- unlist(output_x)
  output_y <- unlist(output_y)
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_xyv[,1]))
  
  if (is.null(output_y)) 
    output_y <- sort(unique(input_xyv[,2]))

  if (any (is.na(output_x)))
    stop ("cannot proceed: some elements of 'output_x' are NA")
  
  if (any (is.na(output_y)))
    stop ("cannot proceed: some elements of 'output_y' are NA")
  
  ll <- length(output_x)
  output_x <- sort(unique(output_x))
  if (length(output_x) != ll) 
    warning ("some of values in output_x were duplicated - have created unique values")
  
  ll <- length(output_y)
  output_y <- sort(unique(output_y))
  if (length(output_y) != ll) 
    warning ("some of values in output_y were duplicated - have created unique values")
  
  dlon <- diff (output_x)      # diff(input_2D$longitude)
  dlat <- diff (output_y)      # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)/mean(dlon)) <= 1e-10 &
    diff(range(dlat)/mean(dlat)) <= 1e-10
  
  if (! equidistant) 
    stop ("cannot proceed -  output grid should be evenly spaced")
  
  range_x <- c(min(output_x), max(output_x), dlon)
  range_y <- c(min(output_y), max(output_y), dlat)
  
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyv[,1]) < min(output_x))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  if (min(input_xyv[,1]) > max(output_x))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyv[,2]) < min(output_y))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  if (min(input_xyv[,2]) > max(output_y))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  storage.mode(input_xyv) <-  "double" 
  storage.mode(range_x) <- storage.mode(range_y) <- "double"
  
  val <- average_xy_2D_r_cpp(input_xyv = input_xyv, 
                               range_x = range_x, 
                               range_y = range_y)

  out <- list(output_x, output_y, val)
  names(out) <- colnames(input_xyv)[1:3]
  return(out)

}

# ==============================================================================
# ==============================================================================

average_2D_2D <- function(input_x, input_y, input_2D, 
                          output_x = NULL, output_y = NULL){
  
  dd <- dim(input_2D)
  
  if (dd[1] != length(input_x)) 
    stop ("length of 'input_x' should be = nrow(input_2D)")
  if (dd[2] != length(input_y)) 
    stop ("length of 'input_y' should be = ncol(input_2D)")
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_x))
  
  if (is.null(output_y)) 
    output_y <- sort(unique(input_y))

  if (any (is.na(input_x)))
    stop ("cannot proceed: some elements of 'input_x' are NA")
  
  if (any (is.na(input_y)))
    stop ("cannot proceed: some elements of 'input_y' are NA")

  # input_x and input_y should contain unique values - if not: stop
  ll <- length(input_x)
  input_x <- unique(input_x)
  if (length(input_x) != ll) 
    stop ("some of values in input_x were duplicated - input_x should contain unique values")
  
  ll <- length(input_y)
  input_y <- unique(input_y)
  if (length(input_y) != ll) 
    stop ("some of values in input_y were duplicated - input_y should contain unique values")
  
  
  if (any (is.na(output_x)))
    stop ("cannot proceed: some elements of 'output_x' are NA")
  
  if (any (is.na(output_y)))
    stop ("cannot proceed: some elements of 'output_y' are NA")
  
  ll <- length(output_x)
  output_x <- sort(unique(output_x))
  if (length(output_x) != ll) 
    warning ("some of values in output_x were duplicated - have created unique values")
  
  ll <- length(output_y)
  output_y <- sort(unique(output_y))
  if (length(output_y) != ll) 
    warning ("some of values in output_y were duplicated - have created unique values")
  
  dlon <- diff (output_x)      # diff(input_2D$longitude)
  dlat <- diff (output_y)      # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)/mean(dlon)) <= 1e-10 &
    diff(range(dlat)/mean(dlat)) <= 1e-10
  
  if (! equidistant) 
    stop ("cannot proceed -  output grid should be evenly spaced")
  
  range_x <- c(min(output_x), max(output_x), dlon)
  range_y <- c(min(output_y), max(output_y), dlat)
  
  storage.mode(input_x)  <- storage.mode(input_y) <- "double"
  storage.mode(input_2D) <- "double"
  
  storage.mode(range_x) <- storage.mode(range_y) <- "double"
  
    Result <- average_2D_2D_r_cpp(
      input_x = input_x,
      input_y = input_y,
      input_2D = input_2D,
      range_x = range_x, 
      range_y = range_y)
    
  out <- list(x=output_x, y=output_y, v=Result)
  out
}

# ==============================================================================
# ==============================================================================
# AVERAGING TIME SERIES
# ==============================================================================
# ==============================================================================

average_xyt <- function(input_xyv,   # matrix: x, y, v1, v2, ....
                        input_t,     # vector:       t1, t2, ...  
                        output_xy,
                        output_t){
  
  ll <- nrow(input_xyv)
  input_xyv <- na.omit(input_xyv)  # cannot work with NAs
  if (nrow(input_xyv)!= ll) 
    warning ("some of values in input_xyv were NAs - wherever that happened, the entire row was removed")
  
  input_t   <- as.numeric(input_t)
  output_t  <- as.numeric(output_t)
  
  input_xy  <- input_xyv[, 1:2]
  input_2Dt <- input_xyv[, -(1:2)]
  
  if (ncol(input_2Dt) !=  length(input_t))
    stop ("input_t and input_xyv are not compatible: ncol(input_xyv)-2 != length(input_t)")
  
  if (min(input_t) > min(output_t) |
      max(input_t) < max(output_t)) 
    stop ("cannot interpolate: times in 'input_t' do not embrace 'output_t'")
  
  # check the output x and y grid
  xx <- sort(unique(output_xy[,1]))  # different x-s
  yy <- sort(unique(output_xy[,2]))  # different y-s
  
  dx <- range(diff(xx))
  
  if (diff(dx) != 0) 
    stop ("dx is not the same everywere") else 
      dx <- dx[1]
  
  dy <- range(diff(yy))
  if (diff(dy) != 0) 
    stop ("dy is not the same everywere") else 
      dy <- dy[1]
  
  # the C++ code requires input of the minimum x, maximum x and dx
  outx <- c(range(xx), dx)
  outy <- c(range(yy), dy)
  
  output_xy <- as.matrix(output_xy)
  
  mode(input_xy) <- mode(input_2Dt) <- "double"
  mode(output_xy) <- mode(outx) <- mode(outy) <- "double"
  mode(input_t)   <- mode(output_t)  <- "double"
  
  out <-  average_grid_t(input_xy,
                         input_2Dt,
                         input_t, 
                         output_xy,
                         outx, 
                         outy, 
                         output_t)
  return(out)
}
