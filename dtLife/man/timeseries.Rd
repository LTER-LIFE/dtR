\name{average_timeseries}
\alias{average_timeseries}

\title{
  averages time series over time
}

\description{
  \code{average_timeseries} takes time averages of time series values. 
}

\usage{
average_timeseries(input, avgOver = "hour",  avgTime  = 1, 
                   datetime = "datetime", value = "value", by = NULL, ...)
}

\arguments{
  \item{input }{2-columned data.frame with time series input data, one of its columns should have a datetime value (either of POSIXct, POSIXlt, or Date. The other column should have the values.
  }
  \item{datetime, value }{Name of the column with the times, and values respectively.   
  }
 \item{avgOver }{Time instance over which has to be averaged, i.e. one of "sec", "min", "hour", "mon", "year", 
   "mday" (day of the month), "wday" (day of the week), "yday" (day of the year).  
   }
 \item{avgTime }{The length of time units over which need to be averaged.
 }
 \item{by }{The name(s) of grouping elements by which the timeseries is aggregated.
 }
 \item{... }{Any other elements passed to the aggregation function.
 }
}

\value{
  Function \code{average_timeseries} returns a \code{data.frame} with the same columns as the original timeseries.
}

\examples{

# ========================================================
# time averages of weather time series - long format
# ========================================================

data(Weather2021)

head(Weather2021)  # timeseries every hour for several stations
dim(Weather2021)

# ================================================
# time averages of several time series at once - wide format
# ================================================

# take averages over all the stations

toselect <- c("windspeed", "temperature", "radiation")
average_timeseries(Weather2021, 
                   avgOver = "mon",  
                   value   = toselect,  
                   na.rm   = TRUE)

# take monthly averages for all stations separately

toselect <- c("windspeed", "temperature", "radiation")
AT <- average_timeseries(
                   Weather2021, 
                   avgOver = "mon",  
                   value   = toselect, 
                   by      = "station",
                   na.rm   = TRUE)

head(AT)

# ================================================
# time averages of one time series in wide format
# ================================================
 
head(WadTempLR)

toaverage <- colnames(WadTempLR)[-1]
average_timeseries(WadTempLR,
                   avgOver = "mon",  
                   value   = toaverage,
                   na.rm   = TRUE
                   )

}

\seealso{
  \code{\link{Marsdiep}} bathymetric data.
  \code{\link{Sediment}} data set with sediment composition.
  \code{\link{WadTempLR}} data set with water temperature .
  \code{\link{Weather2021}} data set with weather information.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
