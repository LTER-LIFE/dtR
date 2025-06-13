## ==========================================
## ==========================================
## Plotting dtLife dataframes
## ==========================================
## ==========================================

plot.dtLife <- function(x, mfrow=NULL, ylab=NULL, main=NULL, 
                       cvar=NULL, type="p", ...){
  
    format <- attributes(x)$format
    
    if (is.null(format))
    if (!"variable" %in% colnames(x)) format <- "wide" else format <- "long"
    
    if (format != "wide") 
      stop ("dataset needs to be in 'wide' format before it can be plotted")
    
    cn    <- colnames(x)
    ntime <- which(cn == "datetime")  # everything preceding this ignored
    
    # best arrangement of the plots
    nv    <- ncol(x) - ntime          # number of variables
    var   <- attributes(x)$variables
    
    if (is.null(mfrow)){
      nc <- min(ceiling(sqrt(nv)), 3)
      nr <- min(ceiling(nv/nc), 3)
      mfrow <- c(nr, nc)
    }
    mf <- par(mfrow=mfrow)
    
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
        if (is.null(Main)) main <- var[iv[1], "description"]
      } else {
        if (is.null(Ylab)) ylab <- "-"
        if (is.null(Main)) main <- varname
      } 
      
      # color variable is present
      if (! is.null (cvar)){
        remove <- unique(c(which (is.na(x$datetime)),which(is.na(x[,i]))))
        if (length(remove)) xx     <- x[-remove, ] else xx <- x # remove NAs
        plot(xx$datetime, xx[,i], xlab="time", ylab=ylab, 
             main=main, type="n", ...)
        
        for (j in 1:length(cv.unique)){
          ii <- which(xx[,cvar]==cv.unique[j])
          
          if (length(ii)){
            ss <- xx[ii, ]
        
          if (! is.null(ss)){
            ss <- ss[order(ss$datetime),]
            lines(ss$datetime, ss[,i], col=j, ...)
          }
          }
        }
      } else {
        xx <- na.omit(data.frame(x$datetime, x[,i]))  # remove NAs
        plot(xx[,1], xx[,2], xlab="time", ylab=ylab, 
             main=main, type=type, ...)
        
      }
    }  
    par(mfrow=mf)
}

## ==========================================
## ==========================================
## Plotting bathymetry data, and points
## ==========================================
## ==========================================

plot_bathymetry <- function(bat, pts=NULL, type=NULL, ptlist=NULL, negativeDepth=FALSE, ...){

  ldots <- list(...)
  
  # type of plot
  if (is.null(pts) & is.null(type))
    type <- "image"
  else if (! is.null(pts) & is.null(type))
    type <- "contour"
  
  # axes
  if (is.null(ldots$las)) 
    ldots$las <- 1
    
  if (is.null(ldots$xlab))
    ldots$xlab=expression(""^o~E)
  
  if (is.null(ldots$ylab))
    ldots$ylab=expression(""^o~N)

  if (is.null(ldots$xlim))
    ldots$xlim <- range(bat$longitude)  
  if (is.null(ldots$ylim))
    ldots$ylim <- range(bat$latitude) 
  
  # aspect ratio
  asp <- ldots$asp
  
  if (is.null(asp))
    asp <- bat$asp
  if (is.null(asp))
    asp <- 1/cos((mean(ldots$ylim) * pi)/180)
  ldots$asp <- NULL
  
  # colorkey
  colkey <- ldots$colkey
  ldots$colkey <- NULL
  
  if (is.null(colkey)) 
    colkey <- list()
#  if (! is.null(pts))
#     colkey$plot <- FALSE
#  else
     colkey$plot <- TRUE
  
  if (is.null(colkey$width))
    colkey$width <- 0.5
  if (is.null(colkey$length))
    colkey$length <- 0.5
  
  if (colkey$plot & is.null(ldots$clab) & is.null(pts))
    ldots$clab <- "m"

  if (! is.null(pts)){

#   if (is.list(pts) & ! is.data.frame(pts)){
#        ptlist <- pts
#        ptlist$longitude <- ptlist$latitude <- NULL
#        pts <- data.frame(longitude=pts$longitude, latitude=pts$latitude)
#        if (! is.null(ptlist$colvar)) 
#          pts$colvar <- ptlist$colvar
#        ptlist$colvar <- NULL
#    }
    if (is.data.frame(pts)) 
      pts <- as.matrix(pts)
    if (is.matrix(pts)) 
    pts <- pts[which(pts[,1] >= ldots$xlim[1]     & 
                     pts[,1] <= ldots$xlim[2]     &
                     pts[,2]  >= ldots$ylim[1]     & 
                     pts[,2]  <= ldots$ylim[2]),  ]

    if (is.null(ldots$pch))
      ldots$pch <- ptlist$pch
    if (is.null(ldots$pch))
      ldots$pch <- 18
    ptlist$pch <- NULL
   }
  if (negativeDepth) 
    depth <- -bat$depth
  else
    depth <- bat$depth
     
  add <- FALSE
    if (type %in% c("image", "both")) {
      do.call("image2D", 
                      c(alist(x=bat$longitude, y=bat$latitude, z=depth, 
                              asp=asp, colkey=colkey), ldots))
      colkey <- list(plot=FALSE)
      add <- TRUE
    }  
    
  if (! is.null(pts)){
    if (is.matrix(pts)){
      x <- pts[,1]
      y <- pts[,2]
      if (ncol(pts) == 2) {
        colvar <- NULL
        ck <- list(plot=FALSE)
      } else{
        colvar <- pts[,3]
        ck <- colkey
      }
    } else {
      colvar <- NULL
      ck <- list(plot=FALSE)
      x <- pts[[1]]
      y <- pts[[2]]
    }
        do.call("points2D", 
                c(alist(x=x, y=y, colvar=colvar,
                        colkey=ck, add=add, asp=asp), 
                c(ldots, ptlist)))
        colkey <- list(plot=FALSE)
        add <- TRUE
    }     
    if (type %in% c("both", "contour")) {
      cnt <- ldots$contour
      ldots$contour <- NULL
      
      ldots      <- c(ldots, cnt)
      levels     <- ldots$levels
      nlevels    <- ldots$nlevels
      drawlabels <- ldots$drawlabels
      if (is.null(drawlabels)) 
        drawlabels <- FALSE
      
      if (is.null(nlevels))
         nlevels <- length(levels)
       
      if (is.null(nlevels)) 
        nlevels <- 10
      
      if (nlevels == 0)
        nlevels <- 10
      ldots$levels <- ldots$nlevels <- ldots$drawlabels <- NULL

      if (is.null(cnt$col)){
         if (type == "both")  ldots$col <- "black"
       } else if (! add){

         if (! is.null(nlevels))
           ldots$col <- ldots$col[as.integer(seq(from=1, to=length(ldots$col), length.out=nlevels))]
       }
        do.call("contour2D",
                      c(alist(x=bat$longitude, y=bat$latitude, z=depth, 
                              levels=levels, nlevels=nlevels,
                              drawlabels=drawlabels, asp=asp,  
                              colkey=colkey, add=add), ldots))
    }
}
