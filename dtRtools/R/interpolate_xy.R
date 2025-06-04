
# ==============================================================================
# ==============================================================================
# Maps xy points to other xy points 
# ==============================================================================
# ==============================================================================

interpolate_xy <- function(
                      input_xyv, # 3 columns: longitude (x), latitude (y), value (v)
                      input_x, input_y, input_2D,
                      output_xy, 
                      output_x, output_y,
                      nmean = 3,
                      asp = c("geographic", "mean", "none")){  # y/x aspect ratio
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
                                nmean     = nmean,
                                asp       = asp) 
  
  else if (! missing (input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing (output_xy))
    
    result <- interpolate_2D_xy (input_x    = input_x, 
                                 input_y    = input_y, 
                                 input_2D   = input_2D,        
                                 output_xy  = output_xy,
                                 asp        = asp)
  
  else if (! missing(input_xyv) & 
           ! missing(output_x)  & ! missing(output_y))
    
    result <- interpolate_xy_2D(input_xyv = input_xyv, 
                                output_x   = output_x, 
                                output_y   = output_y,
                                nmean      = nmean,
                                asp        = asp)
  
  else if (! missing(input_x) & ! missing(input_y) & ! missing(input_2D) & 
           ! missing(output_x) & !missing(output_y))
    
    result <- interpolate_2D_2D (input_x  = input_x, 
                                 input_y  = input_y, 
                                 input_2D = input_2D, 
                                 output_x = output_x, 
                                 output_y = output_y,
                                 asp      = asp)
  
  else 
    stop ("one of the required arguments is missing")
  
  attr(result, "processing") <- paste("remapped , at:", 
                                      Sys.time())
  
  return(result)
}

# ==============================================================================

interpolate_xy_xy <- function(input_xyv, # latitude (x), longitude (y), value
                              output_xy, 
                              nmean = 3,
                              asp = c("geographic", "mean", "none")){
  
  if (ncol(input_xyv) != 3) 
    stop ("'input_xyv' should have 3 columns: latitude (x), longitude (y), value")

  ##### check input 
  input_xyv <- as.matrix(input_xyv[, 1:3])
  ni        <- nrow(input_xyv)
  if (ni == 0) stop("input_xyv is empty")
  
  if (any (is.na(input_xyv[, 1:2])))
    stop ("cannot proceed: some of the input coordinates are NA")
  
  #
  output_xy <- as.matrix(unique(output_xy[, 1:2]))
  output_xy <- na.omit(output_xy)

  if (is.numeric(asp)) {
    asp <- asp[1]
  } else {
    asp <- match.arg (asp, c("geographic", "mean", "none"))
    if (asp == "geographic")
      asp <- aspect_coord(output_xy[,2]) 
    else if (asp == "mean")
      asp = mean(output_xy[,1], na.rm = TRUE)/mean(output_xy[,2], na.rm = TRUE)
    else  
      asp <- 1
  }
  no <- nrow(output_xy)
  if (no == 0) stop("output_xy is empty")

  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyv[,1]) < min(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  
  if (min(input_xyv[,1]) > max(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_xyv[,2]) < min(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  if (min(input_xyv[,2]) > max(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  storage.mode(input_xyv) <- storage.mode(output_xy) <- "double"
  result <- interpolate_xy_xy_cpp(input_xyv = input_xyv, 
                                  output_xy = output_xy, 
                                  nmean     = as.integer(nmean), 
                                  asp       = as.double(asp))
  
  result <- data.frame(output_xy, value=result)
  colnames(result) <- colnames(input_xyv)
  
  result
}

# ==============================================================================
# ==============================================================================
# interpolate 2D gridded points to xy points
# ==============================================================================
# ==============================================================================

interpolate_2D_xy <-  function(
                        input_x, input_y, # latitude (x), longitude (y)
                        input_2D,         #  z-values
                        output_xy,
                        asp = c("geographic", "mean", "none")){

  if (any (is.na(input_x)))
    stop ("cannot proceed: some elements of 'input_x' are NA")
  
  if (any (is.na(input_y)))
    stop ("cannot proceed: some elements of 'input_y' are NA")
  
  # Check overlap of x and of y - IF NO overlap: stop
  if (max(input_x) < min(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  
  if (min(input_x) > max(output_xy[,1]))
    stop( "cannot perform mapping: x-variables of in-output do not overlap")
  
  # Check overlap of x and of y
  if (max(input_y) < min(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  if (min(input_y) > max(output_xy[,2]))
    stop( "cannot perform mapping: y-variables of in-output do not overlap")
  
  ll <- length(input_x)
  input_x <- unique(input_x)
  if (length(input_x) != ll) 
    stop ("some of values in input_x were duplicated - input_x should contain unique values")
  
  ll <- length(input_y)
  input_y <- unique(input_y)
  if (length(input_y) != ll) 
    stop ("some of values in input_y were duplicated - input_y should contain unique values")
  
  dd <- dim(input_2D)
  if (dd[1] != length(input_x)) stop ("length of 'input_x' should be = nrow(input_2D)")
  if (dd[2] != length(input_y)) stop ("length of 'input_y' should be = ncol(input_2D)")
  
  if (is.numeric(asp)) {
    asp <- asp[1]
  } else {
    asp <- match.arg (asp, c("geographic", "mean", "none"))
    if (asp == "geographic")
      asp <- aspect_coord(output_xy[,2]) 
    else if (asp == "mean")
      asp = mean(output_xy[,1], na.rm = TRUE)/mean(output_xy[,2], na.rm = TRUE)
    else  
      asp <- 1
  }
  
  dlon <- diff (input_x)      # diff(input_2D$longitude)
  dlat <- diff (input_y)      # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)/mean(dlon)) <= 1e-10 &
                 diff(range(dlat)/mean(dlat)) <= 1e-10
  
  output_xy <- as.matrix(na.omit(output_xy[, 1:2]))
  input_2D  <- as.matrix(input_2D)
  
  input_x   <- unlist(input_x)  
  input_y   <- unlist(input_y)  
  
  storage.mode(input_x) <- storage.mode(input_y) <- storage.mode(input_2D) <- "double" 
  storage.mode(output_xy) <- "double"

  if (! equidistant){  # SLOWer!
    Result <- interpolate_2D_xy_cpp(input_x   = input_x, 
                                    input_y   = input_y,
                                    input_2D  = input_2D, 
                                    output_xy = output_xy, 
                                    asp       = as.double(asp))
    Result <- data.frame(output_xy, v=Result)
    
  } else {  # equidistant input data points
    
    Result <- interpolate_2D_xy_r_cpp(
                range_x   = as.double(c(min(input_x), max(input_x), dlon[1])), 
                range_y   = as.double(c(min(input_y), max(input_y), dlat[1])),
                input_2D  = input_2D, 
                output_xy = output_xy, 
                asp       = as.double(asp))
    Result <- data.frame(output_xy, v=Result)

  }
  Result 
}

# ==============================================================================
# ==============================================================================
# maps xy values on a 2D regular grid
# ==============================================================================
# ==============================================================================

interpolate_xy_2D <- function(input_xyv, 
                              output_x = NULL, 
                              output_y = NULL,
                              nmean = 3,
                              asp = c("geographic", "mean", "none")){
  
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
  
  if (is.numeric(asp)) {
    asp <- asp[1]
  } else {
    asp <- match.arg (asp, c("geographic", "mean", "none"))
    if (asp == "geographic")
      asp <- aspect_coord(output_y) 
    else if (asp == "mean")
      asp = mean(output_x, na.rm = TRUE)/mean(output_y, na.rm = TRUE)
    else  
      asp <- 1
  }
  
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
  
  storage.mode(input_xyv) <- storage.mode(output_x) <- storage.mode(output_y) <- "double" 

  val <- interpolate_xy_2D_cpp(input_xyv = input_xyv, 
                               output_x = output_x, 
                               output_y = output_y, 
                               nmean = as.integer(nmean), 
                               asp = as.double(asp))

  out <- list(output_x, output_y, val)
  names(out) <- colnames(input_xyv)[1:3]
  return(out)

}

# ==============================================================================
# ==============================================================================

interpolate_2D_2D <- function(input_x, input_y, input_2D, 
                              output_x = NULL, output_y = NULL, 
                              asp = c("geographic", "mean", "none")){
  
  dd <- dim(input_2D)
  
  if (dd[1] != length(input_x)) 
    stop ("length of 'input_x' should be = nrow(input_2D)")
  if (dd[2] != length(input_y)) 
    stop ("length of 'input_y' should be = ncol(input_2D)")
  
  # input_x and input_y should contain unique values - if not: stop
  ll <- length(input_x)
  input_x <- unique(input_x)
  if (length(input_x) != ll) 
    stop ("some of values in input_x were duplicated - input_x should contain unique values")
  
  ll <- length(input_y)
  input_y <- unique(input_y)
  if (length(input_y) != ll) 
    stop ("some of values in input_y were duplicated - input_y should contain unique values")
  
  # output_x and output_y should contain unique values - warning if not
  ll <- length(output_x)
  output_x <- sort(unique(output_x))
  if (length(output_x) != ll) 
    warning ("some of values in output_x were duplicated - have created unique values")
  
  ll <- length(output_y)
  output_y <- sort(unique(output_y))
  if (length(output_y) != ll) 
    warning ("some of values in output_y were duplicated - have created unique values")
  
  if (is.numeric(asp)) {
    asp <- asp[1]
  } else {
    asp <- match.arg (asp, c("geographic", "mean", "none"))
    if (asp == "geographic")
      asp <- aspect_coord(output_y) 
    else if (asp == "mean")
      asp = mean(output_x, na.rm = TRUE)/mean(output_y, na.rm = TRUE)
    else  
      asp <- 1
  }
  
  if (is.null(output_x)) 
    output_x <- sort(unique(input_x))
  
  if (is.null(output_y)) 
    output_y <- sort(unique(input_y))

  if (any (is.na(input_x)))
    stop ("cannot proceed: some elements of 'input_x' are NA")
  
  if (any (is.na(input_y)))
    stop ("cannot proceed: some elements of 'input_y' are NA")

  if (any (is.na(output_x)))
    stop ("cannot proceed: some elements of 'output_x' are NA")
  
  if (any (is.na(output_y)))
    stop ("cannot proceed: some elements of 'output_y' are NA")
  
  dlon <- diff (input_x)      # diff(input_2D$longitude)
  dlat <- diff (input_y)      # diff(input_2D$latitude)
  
  equidistant <- diff(range(dlon)/mean(dlon)) <= 1e-10 &
    diff(range(dlat)/mean(dlat)) <= 1e-10
  
  storage.mode(input_x)  <- storage.mode(input_y) <- "double"
  storage.mode(input_2D) <- "double"
  storage.mode(output_x) <- storage.mode(output_y) <- "double"
  
  if (! equidistant){  # SLOWer!
    Result <- interpolate_2D_2D_cpp (input_x = input_x, 
                                     input_y = input_y, 
                                     input_2D = input_2D, 
                                     output_x = output_x, 
                                     output_y = output_y, 
                                     asp = as.double(asp))
  } else {
    Result <- interpolate_2D_2D_r_cpp(
      range_x   = as.double(c(min(input_x), max(input_x), dlon[1])), 
      range_y   = as.double(c(min(input_y), max(input_y), dlat[1])),
      input_2D  = input_2D, 
      output_x  = output_x, 
      output_y  = output_y, 
      asp       = as.double(asp))
    
  }
  out <- list(x=output_x, y=output_y, v=Result)
  out
}

