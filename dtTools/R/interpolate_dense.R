
interpolate_dense <- function(input_xyv,   # matrix: x, y, v1, v2, ....
                              input_t,     # vector:       t1, t2, ...  
                              output_xy,
                              output_t){
  
  input_xyv <- na.omit(input_xyv)  # cannot work with NAs
  input_t   <- as.numeric(input_t)
  output_t  <- as.numeric(output_t)
  
  if (ncol(input_xyv)-2 !=  length(input_t))
    stop ("input_t and input_xyv are not compatible")

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
  
  mode(input_xyv) <- mode(output_xy) <- mode(outx) <- mode(outy) <- "double"
  mode(input_t)   <- mode(output_t)  <- "double"
  
  out <-  average_grid_t(input_xyv, 
                         input_t, 
                         output_xy,
                         outx, 
                         outy, 
                         output_t)
  return(out)
}
