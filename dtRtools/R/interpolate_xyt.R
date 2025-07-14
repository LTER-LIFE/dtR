# ==============================================================================
# ==============================================================================
# Maps xy timeseries to other xy timeseries
# ==============================================================================
# ==============================================================================

interpolate_xyt <-  function(
              input_xytv,  # longitude (x), latitude (y), time, value
              output_xy, 
              output_x, output_y,
              output_t, 
              nmean = 3,
              asp   = c("geographic", "mean", "none"),  # y/x aspect ratio
              rule  = 2, 
              ID    = NULL)    # unique identifier for the output
                   {   
  
  ##### step 0: prepare/check the inputs #####
  # attributes of input saved
  attrs <- NULL
  
  if (inherits(input_xytv, "dtLife"))
    attrs <- attributes(input_xytv)
  
  # check output
  if (missing (output_xy)){
    if (missing(output_x) | missing(output_y) )
      stop("either output_xy should be input or (output_x & output_y)")
  }
  
  # check input_xytv: 4 columns: longitude (x), latitude (y), time, value
  if (ncol(input_xytv) != 4) 
    stop ("'input_xytv' should have 4 columns: longitude (x), latitude (y), time, value")
  
  # name of the time column
  tname <- names(input_xytv)[3] # to label the output
  if (is.null(tname)) tname <- "datetime"
  
  # input xy values: unique set 0f first two columns
  input_xy  <- as.matrix(unique(input_xytv[, 1:2]))
  
  nmean <- min(nmean, nrow(input_xy))  # at most nmean columns selected
  
  ##### interpolate to output times (output_t) for all inputs #####
  n.unique <- nrow(input_xy)
  input_t <- matrix(nrow = length(output_t), 
                    ncol = n.unique, 
                    data = 0)

  for (i in 1:n.unique){
    # interpolate input to required times
    ii <- which(input_xytv[,1] == input_xy[i, 1] & 
                input_xytv[,2] == input_xy[i, 2])
    V.out  <- approx(x    = input_xytv[ii, 3],   # time
                     y    = input_xytv[ii, 4],   # value
                     xout = output_t, 
                     rule = rule)$y
    input_t[,i] <- V.out
  }
  
  # output xy values: 
  if (! missing (output_xy)){
    
    if (is.vector(output_xy))
      output_xy <- matrix(nrow = 1, ncol = 2, data = output_xy)
    
    else 
      output_xy <- as.matrix(unique(output_xy[,1:2]))
  
    nout <- nrow(output_xy)
  
    # euclidean_distance aspect ratio
    
    if (is.numeric(asp)) {
      asp <- asp[1]
    } else {
      asp <- match.arg (asp, c("geographic", "mean", "none"))
      if (asp == "geographic"){
        if (min(output_xy[,2], na.rm = TRUE) < -90 |
            max(output_xy[,2], na.rm = TRUE) >  90)
          stop ("second column of output_xytv should contain valid latitude, -90:90")
        
        if (min(input_xy[,2], na.rm = TRUE) < -90 |
            max(input_xy[,2], na.rm = TRUE) >  90)
          stop ("second column of input_xy should contain valid latitude, -90:90")
        
        asp <- aspect_coord(output_xy[,2])
      } else if (asp == "mean")
        asp <- mean(output_xy[,1], na.rm = TRUE)/mean(output_xy[,2], na.rm = TRUE)
      else  
        asp <- 1
    }

    storage.mode(input_xy)  <- "double"
    storage.mode(input_t)   <- "double"
    storage.mode(output_xy) <- "double"
  
    Result <- interpolate_ts_cpp (input_xy  = input_xy, 
                                  input_t   = input_t, 
                                  output_xy = output_xy, 
                                  nmean     = as.integer(nmean),
                                  asp       = as.double(asp))

  } else {  #output_x and output_y values
    output_x <- unlist(output_x)
    output_y <- unlist(output_y)
    
    if (is.null(output_x)) 
      output_x <- sort(unique(input_xytv[,1]))
    
    if (is.null(output_y)) 
      output_y <- sort(unique(input_xytv[,2]))
    
    if (any (is.na(output_x)))
      stop ("cannot proceed: some elements of 'output_x' are NA")
    
    if (any (is.na(output_y)))
      stop ("cannot proceed: some elements of 'output_y' are NA")
    
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
    if (max(input_xytv[,1]) < min(output_x))
      stop( "cannot perform mapping: x-variables of in-output do not overlap")
    
    if (min(input_xytv[,1]) > max(output_x))
      stop( "cannot perform mapping: x-variables of in-output do not overlap")
    
    # Check overlap of x and of y - IF NO overlap: stop
    if (max(input_xytv[,2]) < min(output_y))
      stop( "cannot perform mapping: y-variables of in-output do not overlap")
    
    if (min(input_xytv[,2]) > max(output_y))
      stop( "cannot perform mapping: y-variables of in-output do not overlap")
    
    storage.mode(input_xy) <- storage.mode(output_x) <- storage.mode(output_y) <- "double" 
    
    Result <- interpolate_ts_xy_2D_cpp(input_xy = input_xy, 
                                    input_t     = input_t, 
                                    output_x    = output_x, 
                                    output_y    = output_y, 
                                    nmean       = as.integer(nmean), 
                                    asp         = as.double(asp))
    output_xy <- expand.grid(latitude = output_y, longitude = output_x)[, 2:1]  
    nout      <- nrow(output_xy)
    
  }  
  
  c.names <- c("longitude", "latitude", as.character(output_t))
  result <- cbind(output_xy, 
                  as.matrix(Result))      # STAYS A MATRIX
  
  colnames(result) <- c.names # c(tname,  stnames)
  class(result) <- c("dtLife", class(result))
  
  atout <- attributes(result)
  
  if (! is.null(attrs))
    attributes(result) <- c(atout, attrs[!names(attrs) %in% c(names(atout), "class")])
  
  attr(result, "datetime_type") <- class(output_t)
  attr(result, "datetime") <- output_t
  attributes(result)$format <- "wide"
  
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
  if (! is.null(ID)) attr(result, "ID") <- Stations
  
  names(Stations)[1] <- "stations"
  attr(result, "stations") <- Stations
  
  attr(result, "processing") <- c(
    atout$processing, paste("interpolated from 2D-time input, at:", Sys.time()))
  
  # names of the output
  if (is.null(ID))
    stnames <- paste("st", 1:nout, sep="")
  else {
    stnames <- ID
    
    if (length(stnames) != nout) 
      stop("'ID' should be of length = output")
  }
  
  return(result)
}  

