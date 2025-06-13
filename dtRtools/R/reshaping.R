# ==============================================================================
# ==============================================================================
# From wide to long and vice versa
# ==============================================================================
# ==============================================================================

# rough estimate of the format, if not recorded in its attributes
getformat <- function(atts){
  
  if (! is.null(atts$format)) 
    format <- atts$format
  
  else if ("variable" %in% atts$names | 
           "value"    %in% atts$names) 
    format <- "long" 
  
  else 
    format <- "wide"  # can also be long!
  
  return(format)
}

# ==============================================================================

dt_reshape <- function(x,                # data to reshape
                      swap = "station", ...){  # the name of the column to expand or create
  
  atts   <- attributes(x)
  format <- getformat(atts)[1]  # estimate the format
  
  if (format == "wide") {
    xx <- dt_tolong(x, swap = swap, ...)

  } else {
    xx <- dt_towide(x, swap = swap)
  }
  xx
}

# ==============================================================================

dt_towide <- function(x, 
                      swap   = "station", 
                      newcol = NULL){  # name of the column that will be expanded to columns
  
  atts <- attributes(x)
  atts <- atts[! names(atts) %in% c( "names", "row.names")]
  
  cn <- colnames(x)
  
  if (is.numeric(swap)){
    iswap <- swap
    swap  <- cn[iswap]
  } else {
    if (!swap %in% cn)
      stop("cannot go to wide format based on ", swap, 
           " as column does not exist")
    iswap <- which(cn %in% swap)
  }
  # it is assumed that datetime, date or time marks the end of the id section
  if (! is.null(newcol)){
    if (is.numeric(newcol)) 
      nt <- newcol[1]
    else
      nt <- which(cn == newcol)
  } else{
    nt                   <- which(cn == "datetime")
    if (!length(nt)) nt  <- which(cn == "date")
    if (!length(nt)) nt  <- which(cn == "time")
  } 
  
  idvar <- cn   [1:nt]
  idvar <- idvar[-iswap]
  
  xx <- reshape(x, 
                direction="wide", 
                idvar   = idvar,   
                timevar = swap, 
                sep=".")
  
  cn  <- colnames(xx)
  cnx <- cn[ncol(xx)]
  zz <- gregexpr(".", cnx, fixed=TRUE)[[1]]
  hd <- substr(cnx, 1, zz)
  
  cn <- gsub(hd, "", cn, fixed=TRUE)
  colnames(xx) <- cn  
  
  atts$processing <- c(atts$processing , 
                       paste("reshaped to 'wide' at", Sys.time()))
  atts$format <- "wide"
  attributes(xx) <- c(attributes(xx)[c("names", "row.names")],
                      atts)
  xx
}

# ==============================================================================

dt_tolong <- function (x, 
                      swap  = "station",
                      vname = "value",
                      na.rm = TRUE){
  
  atts <- attributes(x)
  atts <- atts[! names(atts) %in% c( "names", "row.names")]
  
  cn <- colnames(x)
  if (swap %in% cn)
    stop("cannot create a column named ", swap, " as this already exists")
  
  # it is assumed that datetime marks the beginning of the varying section
  nt    <- which(cn == "datetime")
  if (!length(nt)) nt <- which(cn == "date")
  if (!length(nt)) nt  <- which(cn == "time")
  
  nc    <- (nt+1):ncol(x)
  xx    <- reshape(x, 
                   direction= "long", 
                   timevar  = swap, 
                   varying  = list(nc) , 
                   idvar    = "value")
  cnxx   <- cn[nc]
  lookup <- data.frame(v=1:length(cnxx), name=cnxx)
  nn     <- lookup$name[match(xx[,swap], lookup$v)]
  xx[,swap] <- nn
  xx            <- xx[,-ncol(xx)]
  colnames(xx)[ncol(xx)] <- vname
  
  atts$processing <- c(atts$processing , 
                       paste("reshaped to 'long' at", Sys.time()))
  if (na.rm) {
    ina <- which(is.na(xx[,ncol(xx)]))
    if (length(ina)){
    xx <- xx[-ina,]
      atts$processing <- c(atts$processing , 
                           paste("removed ", length(ina), "NA values from long format at", Sys.time()))
    }
  }
  atts$format <- "long"
  attributes(xx) <- c(attributes(xx)[c("names", "row.names")],
                      atts)
  xx
}

