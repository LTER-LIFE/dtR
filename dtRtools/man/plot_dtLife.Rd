\name{plotting data}
\alias{plot.dtLife}
\title{
  plots data.frames in long or wide format, that are of dtLife class
}

\description{
  \code{plot.dtLife} makes a set of plots for all variables in an object, against column datetime.
}

\usage{
\method{plot}{dtLife}(x, mfrow = NULL, ylab = NULL, main = NULL, 
      use_names = FALSE, cvar = NULL, type = "p", ...)
}

\arguments{
  \item{x }{The object of class \code{dtLife} that is to be plotted.
  }
  \item{ylab, main }{arguments passed to the plotting functions. 
  If NULL, the names will be extracted from the metadata. For ylab, this will be the column \code{unit}; for main, this will be \code{description}, unless \code{use_names} is TRUE.
  }  
  \item{use_names }{when \code{FALSE}, and \code{main} is \code{NULL}, then the \code{description} of the variables metadata will be used to label the plot, if  \code{use_names = FALSE}, the variable names will be used.
  } 
  \item{mfrow, type }{arguments passed to the plotting functions. If \code{mfrow} = NULL, the default, then the best arrangement will be estimated, based on the number of plots; upon returning from the function, the original mfrow will be reset. If  \code{mfrow} = NA, then the settings imposed outside of the function will not be altered.
  }  
  \item{cvar }{ Name of the variable that is used to color the plots.
  }
  \item{... }{Any other arguments passed to the plotting functions.
  }
}

\value{
  Function \code{meta} returns a \code{data.frame} with the objects attributes, except the names, row.names and class.
}

\examples{

# ---------------------------------------------
# plotting weather data in long format
# ---------------------------------------------

# The plotting function makes use of the attributes of the dtLife object
# attributes are kept when subsetting the object

head(Wad_weather)

# columns are the variables; stations are in column "station"
# during plotting, the station column will be ignored
# default plots the description of the variables (too long here)

plot(subset(Wad_weather, 
            subset = station == 235,
            select = c("datetime", "radiation", "temperature")),
     type = "l")

# how to color the datapoints can be specified
plot(Wad_weather, type = "l", cvar = "station", 
     use_names = TRUE)

# to plot the temperature for different stations, 
# the data need to be converted to wide format first

select   <- c("station", "datetime", "temperature")
Wad_temp <- dt_towide(Wad_weather[, select])

plot(Wad_temp, type = "l", ylim = c(-15, 30)) 

# -----------------------------------------------
# plotting water temperature data in wide format
# -----------------------------------------------

head(Wad_watertemp_LR)

# this data is in wide format
# columns are stations; data are temperatures
# xlim is set to ensure that all plots have the same x-axis

plot(Wad_watertemp_LR[, 1:5], ylab = "dgC", las = 1, 
     xlim = range(Wad_watertemp_LR$datetime))
     
}

\seealso{
  \code{\link{Marsdiep}} bathymetric data.
  \code{\link{Wad_watertemp_LR}} data set with (low resolution) water temperature data from the Waddensea.
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
