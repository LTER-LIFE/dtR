# ==============================================================================
# ==============================================================================
# make datasets globally available
# ==============================================================================
# ==============================================================================

if (getRversion() >= "2.15.1")  
  utils::globalVariables(c("Marsdiep", "Wad_waterTempLR", "Shape",
                           "Wad_weather", "Sediment",
                           "KNMIstations", "RWSstations",
                           "Veluwe_mean_tree_temperature",
                           "Veluwe_tree_budstage"))

# ==============================================================================
# ==============================================================================
# Extracting metadata from an object
# ==============================================================================
# ==============================================================================

meta <- function(x){
  
  if (inherits(x, "dtBathymetry"))
    Meta <- x[!names(x) %in% c("longitude", "latitude", "depth", "contours")]
  
  else{
    att <- attributes(x)
    natt <- names(att)
    Meta <- att[natt[!natt %in% c("names", "row.names", "class", "reshapeWide")]]
  }
  
  Meta
}

# ==============================================================================
# ==============================================================================
# Subsetting without loosing the attributes
# ==============================================================================
# ==============================================================================

`[.dtLife` <- function(x, i, j, ...) {
  attrs <- attributes(x)
  cls   <- class(x)
  out   <- x
  
  class(out)   <- class(out)[2] # to avoid recursion

  out   <- out[i, j, ...]       # subset it
  atout <- attributes(out)
  
  if (is.data.frame(out)){
    
    attributes(out) <- c(attributes(out), 
                         attrs[!names(attrs) %in% c(names(atout), "row.names", "names")])
    class(out) <- c("dtLife", class(out))
    
    # adapt stations and variables  
    nout <- names(out)
    
    if ("station" %in% nout & ! is.null(attributes(out)$stations)){
      stations <- unique(out$station)
      attstat  <- attributes(out)$stations
      attstat  <- subset(attstat, 
                         subset = attstat$station %in% stations)
      attributes(out)$stations <- attstat
      
    } else if (! is.null(attributes(out)$stations)) {
      attstat <- attributes(out)$stations
      ii <- which(!attstat$station  %in% colnames(out))
      if (length(ii)) attstat <- attstat[-ii, ]
      attributes(out)$stations <- attstat
      
    }
    
    if ("variable" %in% nout & ! is.null(attributes(out)$variables)){
      variables <- unique(out$variable)
      attvar <- attributes(out)$variables
      attvar <- subset(attvar, 
                       subset = attstat$variable %in% variables)
      attributes(out)$variables <- attvar
    }
    
    new <- "Selected" 
    if (! missing(i) & ! missing(j)) 
       new <- paste(new, length(i), "rows and", length(j), "columns") 
    else if (! missing(i)) 
      new <- paste(new, length(i), "rows") 
    else if (! missing(j)) 
      new <- paste(new, length(j), "columns")
  
    attributes(out)$processing <- c(attributes(out)$processing,
                                  paste(new, "at", Sys.time()))
  } 
  out
}

# ==============================================================================
# ==============================================================================

subset.dtLife <- function(x, subset, select, drop = FALSE, ..., attr = NULL){
  attrs <- attributes(x)
  cls   <- class(x)
  out   <- x
  class(out) <- class(out)[2]

    
  e <- substitute(subset)
  if (! missing(e)) {
  r <- eval(e, out, parent.frame())
  } else {r <- rep(TRUE, length = nrow(out))}
  if (!is.logical(r)) stop("'subset' must be logical")
  r <- r & !is.na(r)
  if (length(r) != nrow(out)) stop ("'subset' evaluation did not provide a selection?")
  if (sum(r) == 0) stop ("'subset' evaluation did not provide a selection?")
  
  vars <- if (missing(select))
    rep_len(TRUE, ncol(out))
  else {        
    nl <- as.list(seq_along(out))
    names(nl) <- names(out)
    eval(substitute(select), nl, parent.frame())
  }  
  
  out <- out[r, vars, drop = drop]    # take subset

  if (is.data.frame(out)){
    atout <- attributes(out)
    attributes(out) <- c(atout, attrs[!names(attrs) %in% names(atout)])

    attributes(out)$processing <- c(attributes(out)$processing,
            paste("Subsetted", sum(r), "from", nrow(x), "rows at", Sys.time())
                                  )
    if (length(attr))
      attributes(out) <- c(attributes(out), attr)

    # adapt stations and variables  
    nout <- names(out)
    if ("station" %in% nout & ! is.null(attributes(out)$stations)){
      stations <- unique(out$station)
      attstat <- attributes(out)$stations
      attstat <- attstat[attstat$station %in% stations,]
      attributes(out)$stations <- attstat
    }
    if ("variable" %in% nout & ! is.null(attributes(out)$variables)){
      variables <- unique(out$variable)
      attvar <- attributes(out)$variables
      attvar <- attvar[attvar$variable %in% variables,]
      attributes(out)$variables <- attvar
    }
    
    class(out) <- cls
  }
  out
}
