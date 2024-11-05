# ==============================================================================
# ==============================================================================
# average time series over time
# ==============================================================================
# ==============================================================================

average_timeseries <- function(input, 
                               avgOver  = "hour", 
                               avgTime  = 1, 
                               datetime = "datetime", 
                               value    = "value",
                               by       = NULL,
                               ... ){
  
   if ( avgOver == "year" & ! avgTime%%1 == 0) 
     avgOver <- "yy" 
   
   if ( avgOver %in% c("mday", "wday", "yday", "mon", "year")) {
    
#    if (avgTime %% 1 != 0 ) 
#      stop("'avgTime' should be an integer") 
     Date <- input[,datetime]
     
     if (is.numeric(Date)){  # see if it is defined in the object's attributes
       metadate <- meta(input)$datetime
       if (! is.null(metadate))
         Date <- metadate
     }
       
     DD      <- as.POSIXlt(Date)
  
     DDnames <- names(unlist(DD[1]))
     ii      <- which(DDnames == avgOver)
     
     if (! length(ii)) 
        stop("'avgOver' should be one of: ", paste(DDnames, collapse=", ")) 
    
     DD <- as.integer((DD[, avgOver] %/% avgTime) * avgTime + avgTime/2)
     
     if (length(unique(na.omit(DD))) < 1)
       stop("averaging time unit 'avgOver' and time 'avgTime' does not create suitable averages")
     
  } else  if (! avgOver %in% c("sec", "min", "hour", "day")) 
    stop ("avgOver should be one of sec, min, hour, day, year, mday, wday, yday, mon")
  
  else {  # in seconds
    avgSecond <- switch(avgOver,
                    "sec"  = 1,
                    "min"  = 60,
                    "hour" = 3600,
                    "day"  = 86400,
                    "yy" = 86400*365.25 )
    fac <- avgTime* avgSecond                  
    DD <- as.numeric(input[,datetime]) %/% fac
    DD <- as.POSIXct(DD*fac, origin="1970-01-01")
  }
#  ignore <- which(colnames(input) %in% c(datetime, value))
#    DI <- data.frame(input[, -ignore], DD)
  if (is.null(by)) 
    DI <- data.frame(DD)
  else
    DI <- data.frame(input[, by], DD)
  
  HH <- aggregate(input[,value], by= as.list(DI), FUN=mean, ...)
  
  colnames(HH) <- c(by, datetime, value)
  if (is.null(by))
     ii <- order(HH[, datetime])
  else
     ii <- order(HH[, by[1]], HH[,datetime])
   
  HH <- HH[ii,]
  
  if (inherits(input, "dtLife")){
    attrs <- attributes(input)
    avoid <- unique(c(names(attributes(HH)), "dim", "dimnames"))
    attributes(HH) <- c(attributes(HH), 
                        attrs[!names(attrs) %in% avoid])
    attributes(HH)$processing <- c(attributes(HH)$processing,
                                   paste("averaged over", avgTime, avgOver, 
                                         "at", Sys.time()))
                        
    class(HH) <- class(input)    
  }
  HH  
}
