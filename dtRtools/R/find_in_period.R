# ==============================================================================
# ==============================================================================
# find first, last or number of occurrences in a period
# ==============================================================================
# ==============================================================================

first_in_period <- function(x,          # dataframe
                            datetime,   # name of data column
                            value,      # name of value column 
                            period,     # name of period column
                            target = 0, # target value that value should have
                            by     = NULL, # group by
                            type   = c("both", "ascending", "descending")) 
  
  find_in_period(x =,x, datetime = {{datetime}}, 
                 value = {{value}}, period = {{period}}, 
                 target = {{target}}, by = {{by}}, 
                 type = type, which = "first")

# ==============================================================================

last_in_period <- function(x, 
                           datetime, 
                           value, 
                           period, 
                           target = 0, 
                           by     = NULL, 
                           type   = c("both", "ascending", "descending")) 
  
  find_in_period(x=,x, datetime = {{datetime}}, 
                 value = {{value}}, period = {{period}}, 
                 target = {{target}}, by = {{by}}, 
                 type = type, which = "last")

# ==============================================================================

count_in_period <- function(x, 
                            datetime, 
                            value, 
                            period, 
                            target = 0, 
                            by     = NULL, 
                            type   = c("both", "ascending", "descending")) 
  
  find_in_period(x=,x, datetime = {{datetime}}, 
                 value = {{value}}, period = {{period}}, 
                 target = {{target}}, by = {{by}}, 
                 type = type, which = "count")

# ==============================================================================

all_in_period <- function(x, 
                          datetime, 
                          value, 
                          period, 
                          target = 0, 
                          by     = NULL, 
                          type   = c("both", "ascending", "descending")) 
  
  find_in_period(x=,x, datetime = {{datetime}}, 
                 value = {{value}}, period = {{period}}, 
                 target = {{target}}, by = {{by}}, 
                 type = type, which = "all")

# ==============================================================================
# ==============================================================================
# main function
# ==============================================================================
# ==============================================================================

find_in_period <- function(x, 
         datetime, 
         value, 
         period, 
         target = 0, 
         by     = NULL, 
         type   = c("both", "ascending", "descending"),
         which = "first"           # first, last, count, all, 
){
  
  att  <- attributes(x)
  natt <- names(att)
  
  # argument checking
  type <- match.arg(type)
  
  if (missing(x))      
    stop("'x' should be a data.frame, tibble or matrix")
  
  if (missing(datetime))   
    stop("'datetime' should be a numeric vector or a vector with POSIXct values")
  
  if (missing(period)) 
    stop("'period' should be given")
  
  # get the names of the time, period, value, and by argument
  ntime   <- deparse(substitute(datetime, env = parent.frame()))
  ntime   <- gsub(x = ntime, "\"", "")  # in case it was a string
  
  nperiod <- deparse(substitute(period,   env = parent.frame())) 
  nperiod <- gsub(x = nperiod, "\"", "")
  
  nvalue  <- deparse(substitute(value,    env = parent.frame()))
  nvalue  <- gsub(x = nvalue, "\"", "")
  
  nby     <- deparse(substitute(by,       env = parent.frame()))
  nby     <- gsub(x = nby, "\"", "")
  
  if (nby == "NULL")  
    nby   <- NULL
  
  namesby    <- c(nby, nperiod)
  namesorder <- c(ntime, nby)
  
  if (! all(c(namesby, ntime, nvalue) %in% colnames(x)))
    stop ("not all arguments are present in 'x'")
  
  # remove NAs and 
  # calculate whether value is smaller than target, grouped by (by and datetime)
  x  <- x[! is.na(x[ , nvalue]), ]    
  
  vv <- alist(x[, namesorder[1]])
  if (length(namesorder) > 1)
    for (i in 2: length(namesorder))
      vv[[i]] <- x[, namesorder[i]]
  
  x         <- x[do.call("order", vv),]
  x$smaller <- x[, nvalue] <= target 

  # if cumsum: make cumulative sums per factor
  # Note; more complex
#  if (which == "cumsum"){
#     as.vector(aggregate(x[, value], 
#                         by = list(x[, nby]), 
#                         FUN  = function(x) cumsum(x))[, -1])
#  }
  
  # find where consecutive values change from smaller <-> larger 
  # and that belong to same period
  embrace <- which(x$smaller[-nrow(x)]  != x$smaller[-1]                       &
                   x[-nrow(x), nperiod] == x[-1, nperiod])                
  
  # extract time and values from these
  zz <- data.frame(      x[embrace,  namesby],
                    t1 = x[embrace,    ntime],
                    t2 = x[embrace+1,  ntime],
                    v1 = x[embrace,   nvalue],
                    v2 = x[embrace+1, nvalue])
  
  if (is.character(zz$t1)){
    zz$t1 <- as.Date(zz$t1)
    zz$t2 <- as.Date(zz$t2)
  }
  
  zz$asc <- (zz$v2 > zz$v1)
  
  if (type == "ascending") 
    zz <- zz[  zz$asc, ]
  
  else if (type == "descending") 
    zz <- zz[ !zz$asc, ]
  
  if (nrow(zz) == 0) return (NULL)
  
  names(zz)[1:length(namesby)] <- namesby
  
  
  if (which == "count")  {# number of occurrences
    Z <- aggregate(1:nrow(zz), 
                   by  = data.frame(zz[,namesby]), 
                   FUN = function(x) length(x))
    result <- Z
    names(result)[1:length(namesby)] <- namesby
    names(result)[ncol(result)] <- "count"
    
  } else {
    
    if (which == "all")   # all occurrences
      result <- zz
    
    else if (which == "last")   # last occurrence
      Z <- aggregate(1:nrow(zz), 
                     by  = data.frame(zz[, namesby]), 
                     FUN = function(x) x[length(x)])
    
    else            # first occurrence
      Z <- aggregate(1:nrow(zz), 
                     by  = data.frame(zz[, namesby]), 
                     FUN = function(x) x[1])
    
    # extract it
    if (which != "all") 
      result <- zz[Z[,ncol(Z)] , ]
    
    # estimate time of occurrence
    result$first     <- with(result, t1 + (target-v1) * (t2-t1)/(v2-v1) ) 
    
    ifirst <- which(names(result) == "first")
    
    if (which == "last") 
      names(result)[ifirst] <- "last"
    
    else if (which == "all") 
      names(result)[ifirst] <- "match"
    
    result$ascending <- result$asc
    
    # remove columns
    result$t1 <- result$t2 <- result$v1 <- result$v2 <- result$asc <- NULL
    
  }  
  result$target    <- target                                                  
  
  result <- as.data.frame(result[order(result[,nperiod]),])
  
  # add attributes
  Meta <- att[natt[!natt %in% names(attributes(result))]]
  attributes(result) <- c(attributes(result), Meta)
  attributes(result)$processing <- c(attributes(result)$processing,
            paste("estimated first", ntime, "in period ", nperiod, "by", nby, 
                  "setting the target for ", nvalue, "equal to", target, 
                   "at", Sys.time()))
  
  result
}
