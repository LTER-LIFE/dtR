# ==============================================================================
# ==============================================================================
## Plotting bathymetry data, and points
# ==============================================================================
# ==============================================================================

plot_bathymetry <- function(bat, 
                            type = NULL, 
                            pts = NULL, ptlist = NULL, 
                            negativeDepth = FALSE, 
                            ...)
      plot_map(bat, type = type, 
               pts = pts, ptlist = ptlist, 
               negativeDepth = negativeDepth, ...) 

# ==============================================================================
# ==============================================================================
## Plotting 2D data, and points
# ==============================================================================
# ==============================================================================

plot_dtmap      <- function(bat, # typically a list with latitude, longitude, matrix
                            varname = NULL, # name of matrix to plot
                            type = NULL, # image or contour
                            pts = NULL, ptlist = NULL,  ...){
  
  ldots <- list(...)
  
  use_points <- FALSE
  
  if (! missing(bat)) {
    
    if (is.data.frame(bat) | is.matrix(bat)) {
      use_points <- TRUE
      if (! is.null(pts)) 
        stop ("cannot proceed with both 'bat' and 'pts' being a data.frame or matrix")
      pts <- bat
    }

  } else   # missing (bat)
    use_points <- TRUE
  
  if (use_points){
    
    if (is.null(pts)) 
      stop ("at least one of 'bat' or 'pts' should be specified")
    
    # ... 
    ldots  <- update_ldots(ldots, y = pts[, 2])
    
    colkey <- get_colkey(ldots)
    ldots$colkey <- NULL
    
    # color variable
    if (ncol(pts) >= 3) 
      colvar <- pts[, 3] 
    else 
      colvar <- NULL
    
    do.call("points2D", 
            c(alist(x = pts[,1], y = pts[, 2], colvar = colvar,
                    colkey = colkey), 
                    c(ldots, ptlist)))
  } else {

    # find matrix to plot
    if (is.null(varname)){  
      ism <- lapply(bat, 
                    FUN = function(x) is.matrix(x) | is.data.frame(x))
      
      varname <- names(bat)[which(unlist(ism))][1]  # take the first
    }
  
    plot_map(bat, type = type, 
             pts = pts, ptlist = ptlist, 
             name = varname, ...) 
  }
}

# ==============================================================================
# ==============================================================================
## Main plotting function of 2D maps
# ==============================================================================
# ==============================================================================

# bat is a list with latitude, longitude and a 2D matrix
plot_map <- function(bat, type = NULL, 
                     pts = NULL, ptlist = NULL, negativeDepth = FALSE, 
                     name = "depth", ...){

  ldots <- list(...)
  
  # type of plot
  if (is.null(pts) & is.null(type))
    type <- "image"
  
  else if (! is.null(pts) & is.null(type))
    type <- "contour"
  
  if (negativeDepth) 
    depth <- -bat[[name]]
  else
    depth <- bat[[name]]
  
  bx <- bat$longitude
  if (is.null(bx)) 
    bx <- bat$x
  if (is.null(bx)) 
    stop("cannot find 'longitude' or 'x' in 'bat'")
  
  by <- bat$latitude
  if (is.null(by)) 
    by <- bat$y
  if (is.null(by)) 
    stop("cannot find 'latitude' or 'y' in 'bat'")
  
  # axes
  if (is.null(ldots$xlim))
    ldots$xlim <- range(bx)  
  
  if (is.null(ldots$ylim))
    ldots$ylim <- range(by) 
  
  ldots  <- update_ldots(ldots, y = ldots$ylim)
  
  colkey <- get_colkey(ldots)
  ldots$colkey <- NULL
  
  if (colkey$plot          & 
      is.null(ldots$clab)  & 
      is.null(pts))   ldots$clab <- "m"
  
  if (! is.null(pts)){

    if (is.data.frame(pts)) 
      pts <- as.matrix(pts)
    
    if (is.matrix(pts)) 
    
      pts <- pts[which(pts[,1]  >= ldots$xlim[1]     & 
                       pts[,1]  <= ldots$xlim[2]     &
                       pts[,2]  >= ldots$ylim[1]     & 
                       pts[,2]  <= ldots$ylim[2]),  ]

    if (is.null(ldots$pch))
      ldots$pch <- ptlist$pch
    if (is.null(ldots$pch))
      ldots$pch <- 18
    ptlist$pch <- NULL
  }
     
  add <- FALSE
  
  if (type %in% c("image", "both")) {
      do.call("image2D", 
              c(alist(x = bx, y = by, z = depth, 
                      colkey = colkey), ldots))
      colkey <- list(plot = FALSE)
      add <- TRUE
  }  
    
  if (! is.null(pts)){
    
    if (is.matrix(pts) | is.data.frame(pts)){
      x <- pts[ ,1]
      y <- pts[ ,2]
      
      if (ncol(pts) == 2) {
        colvar <- NULL
        ck <- list(plot = FALSE)
        
        if (is.null(ptlist))
          ptlist <- list(col = "black")
        else
          if (is.null(ptlist$col)) ptlist$col <- "black"
      } else{
        colvar <- pts[,3]
        ck <- colkey
      }
      
    } else {
      colvar <- NULL
      ck <- list(plot = FALSE)
      x <- pts[[1]]
      y <- pts[[2]]
    }
    do.call("points2D", 
            c(alist(x = x, y = y, colvar = colvar,
                    colkey = ck, add = add), 
            c(ldots, ptlist)))
    
    colkey <- list(plot = FALSE)
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
           ldots$col <- ldots$col[as.integer(seq(from = 1, 
                                                 to   = length(ldots$col), 
                                                 length.out = nlevels))]
       }
        do.call("contour2D",
                c(alist(x = bx, y = by, z = depth, 
                        levels = levels, nlevels = nlevels,
                        drawlabels = drawlabels,   
                        colkey = colkey, add = add), 
                  ldots))
  }
}

# ==============================================================================
# ==============================================================================
## Accessory functions
# ==============================================================================
# ==============================================================================

## ==========================================
## update ellipsis
## ==========================================

update_ldots <- function(ldots, y){
  if (is.null(ldots$las)) 
    ldots$las <- 1
  
  if (is.null(ldots$xlab))
    ldots$xlab <- expression(""^o~E)
  
  if (is.null(ldots$ylab))
    ldots$ylab <- expression(""^o~N)
  
  # aspect ratio
  if (is.null(ldots$asp))
    ldots$asp <- 1/cos((mean(y, na.rm = TRUE) * pi)/180)
  
  ldots
}

## ==========================================
## update the color key
## ==========================================

get_colkey <- function(ldots){
  colkey <- ldots$colkey
  
  
  if (is.null(colkey)) 
    colkey <- list()
  
  else if (is.logical(colkey)){
    colkey <- list(plot = colkey)
  }
  
  if (is.null(colkey$plot)) 
    colkey$plot <- TRUE
  
  if (is.null(colkey$width))
    colkey$width <- 0.5
  
  if (is.null(colkey$length))
    colkey$length <- 0.5
  colkey
}

