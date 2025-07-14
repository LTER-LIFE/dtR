# ==============================================================================
# ==============================================================================
## Plotting dtLife dataframes
# ==============================================================================
# ==============================================================================

plot.dtLife <- function(x, mfrow = NULL, ylab = NULL, main = NULL, 
                        use_names = FALSE,  # names to label the plot
                        cvar = NULL,        # color variable
                        type = "p", ...){
  
    attn <- attributes(x)
    format <- attributes(x)$format
    
    if (is.null(format))
       if (!"variable" %in% colnames(x)) format <- "wide" else format <- "long"
    
    cn    <- colnames(x)
    ntime <- which(cn == "datetime")  # everything preceding this ignored
    
    if (format != "wide") { # rearrange x
      ii <- 1:(ncol(x)-1)
      ii <- ii[-ntime]
      ii <- c(ii, ntime, ncol(x))
      x <- x[ ,ii]
      cn <- colnames(x)
      ntime <- which(cn == "datetime")  # everything preceding this ignored
    }
    #      stop ("dataset needs to be in 'wide' format before it can be plotted")
    
    # best arrangement of the plots
    nv    <- ncol(x) - ntime          # number of variables
    var   <- attn$variables
    
    if (is.null(mfrow)){
      nc <- min(ceiling(sqrt(nv)), 3)
      nr <- min(ceiling(nv/nc), 3)
      mfrow <- c(nr, nc)
    }
    if (! is.na(mfrow[1])) mf <- par(mfrow = mfrow)
    
    Main <- main
    Ylab <- ylab
    
    if (! is.null(cvar)) {    # variable used for coloring the plots
      cv.unique <- unique(x[,cvar])
      if (is.null(cv.unique))
        stop ("cannot plot - 'cvar' unknown")
    }
    
    for (i in (ntime+1) : ncol(x)){
      
      # units and names of the variables
      varname <- cn[i]
      iv <- which(var$variable ==  varname)
     
      if (length(iv)){
        if (is.null(Ylab)) ylab <- var[iv[1], "unit"]
        if (is.null(Main)) {
          if (use_names )
            main <- varname
          else  
            main <- var[iv[1], "description"]
        }
      } else {
        if (is.null(Ylab)) ylab <- "-"
        if (is.null(Main)) main <- varname
      } 
      
      # color variable is present
      if (! is.null (cvar)){
        remove <- unique(c(which (is.na(x$datetime)),
                           which(is.na(x[,i])))      )
        if (length(remove)) xx     <- x[-remove, ] else xx <- x # remove NAs
        plot(xx$datetime, xx[,i], 
             xlab = "time", ylab = ylab, 
             main = main, type = "n", ...)
        
        for (j in 1:length(cv.unique)){
          ii <- which(xx[,cvar] == cv.unique[j])
          
          if (length(ii)){
            ss <- xx[ii, ]
        
          if (! is.null(ss)){
            ss <- ss[order(ss$datetime),]
            lines(ss$datetime, ss[,i], col = j, ...)
          }
          }
        }
      } else {
        xx <- na.omit(data.frame(x$datetime, 
                                 x[,i]))  # remove NAs
        plot(xx[,1], xx[,2], 
             xlab = "time", ylab = ylab, 
             main = main, type = type, ...)
        
      }
    }  
    if (! is.na(mfrow[1])) par(mfrow = mf)  # restore
}
