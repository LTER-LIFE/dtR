\name{reshape timeseries}
\alias{dtreshape}
\alias{dttowide}
\alias{dttolong}

\title{
 Reshapes timeseries from long to wide format and vice versa - only usable for simple conversions.
}

\description{
  \code{dtreshape} is a general function to cast a dataset in another format. 
  
  \code{dttowide} transforms a dataset from long to wide format. 

  \code{dttolong} transforms a dataset from wide to long format. 

  All these functions require that a column called "datetime" is present.
}
\usage{
dtreshape (x, swap = "station", ...)
dttolong  (x, swap = "station", vname = "value", na.rm=TRUE)
dttowide  (x, swap = "station")

}

\arguments{
  \item{x }{Data.frame with data to reshape. It should contain a column called "datetime".  
  }
  \item{swap }{The name of the column to expand (when going to wide format) or to collapse (when going to long format). 
  }
  \item{vname }{The name of the column that contains the values in the long format.
  }
  \item{na.rm }{To remove the NA values in the long format.
  }
  \item{... }{Any other arguments passed to dttolong.
  }
}

\details{
These functions are based on the function \link[stats]{reshape}. 

In both dttowide and dttolong, the column called \code{datetime} plays an important role:

in \code{dttowide}, \code{datetime} marks the end of the \code{id} section. 
The \code{idvar} are all columns prior to and including \code{datetime}, except for the column name passed via argument \code{swap}. 

for \code{dttolong}, it is assumed that \code{datetime} marks the beginning of the \code{varying} section.

\code{swap} is passed as the timevar

}

\value{
  Function \code{dttowide} returns a \code{data.frame} in \code{wide} format.
  
  Function \code{dttolong} returns a \code{data.frame} in \code{long} format.
  
  Function \code{dtreshape} tries to derive the original format of the dataset and returns the data in the other format.

}


\examples{

## =============================================================
## a simple dataset in wide format
## =============================================================

# low resolution temperature data

head(WadTempLR)
meta(WadTempLR)$format

# reshape to long format, 
# keep the na values
# ----------------------------

WadTempL <- dtreshape(x     = WadTempLR, 
                      na.rm = FALSE)

head(WadTempL)
meta(WadTempL)$format

# reshape to long format, 
# remove the na values
# ----------------------------

WadTempL2 <- dtreshape(x    = WadTempLR)

head(WadTempL2)

# back to original format
# ----------------------------

WW       <- dtreshape(x     = WadTempL2)

head(WW)  # ordering is different

head(WW[order(WW$datetime), ])   # original ordering

## =============================================================
## convert between two wide formats
## =============================================================

# three weather variables, multiple stations
# columns = weather variables

head(Weather2021)
meta(Weather2021)$format

# Convert to long format
# ----------------------------

# does not work with default arguments (swaps according to station):
# head(dttolong(Weather2021))

WeatherL <- dttolong(x    = Weather2021, 
                     swap = "variable")   # name of the long variable

head(WeatherL)

# rearrange it so that the variable name comes before the datetime
#                      the station name after the datetime

WeatherL2 <- WeatherL[ ,c(3,2,1,4)]
head(WeatherL2)

# reshape to wide format; columns are now the stations 

WeatherS <- dtreshape(x    = WeatherL2, 
                      swap = "station")  # long name to expand

head(WeatherS)
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
