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
  \item{input }{data.frame with time series input data, one of its columns should have a datetime value (either of POSIXct, POSIXlt, or Date. The other column(s) should contain the values to average.
  }
  \item{datetime }{Name of the column with the times (only one name)  
  }
  \item{value }{Name of the column(s) with the values.   
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

data(Wad_weather)

head(Wad_weather)  # timeseries every hour for several stations
dim(Wad_weather)

# ========================================================
# daily averages, subset of the data
# ========================================================

WW208 <- subset(Wad_weather, 
                subset = (station == 208 & 
                          datetime < "2021-01-15 00:00:00"))

ZZ <- average_timeseries(WW208, avgOver = "day", value = "temperature")

plot(WW208$datetime, WW208$temperature, 
     xlab = " time", ylab = "dgC", main = "temperature timeseries",
     type = "l")
     
T0 <- WW208$datetime[1]     
abline(v = seq(from = T0, to = T0 + 86400*15, by = 86400),
       lty = 2)     

points(ZZ, col = "red", pch = 16)

legend("bottomleft", col = 1:2, lty = c(1, NA), pch = c(NA, 16),
       legend = c("hourly (original)", "daily (averaged)"), cex = 0.7)

# ========================================================
# time averages of several time series at once - wide format
# ========================================================

# take averages over all the stations

toselect <- c("windspeed", "temperature", "radiation")
average_timeseries(Wad_weather, 
                   avgOver = "mon",  
                   value   = toselect,  
                   na.rm   = TRUE)

# take monthly averages for all stations separately

toselect <- c("windspeed", "temperature", "radiation")
AT <- average_timeseries(
                   Wad_weather, 
                   avgOver = "mon",  
                   value   = toselect, 
                   by      = "station",
                   na.rm   = TRUE)

head(AT)

# ========================================================
# time averages of one time series in wide format
# ========================================================
 
head(Wad_watertemp_LR)

toaverage <- colnames(Wad_watertemp_LR)[-1]
average_timeseries(Wad_watertemp_LR,
                   avgOver = "mon",  
                   value   = toaverage,
                   na.rm   = TRUE
                   )

}

\seealso{
  \code{\link{Marsdiep}} bathymetric data.
  \code{\link{Sediment}} data set with sediment composition.
  \code{\link{Wad_watertemp_LR}} data set with (low resolution) water temperature from the Waddensea.
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
