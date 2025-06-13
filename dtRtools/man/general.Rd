\name{meta}
\alias{meta}
\alias{plot.dtLife}
\title{
  extract metadata from an object - keeps the attributes for [] and \code{subset}
}

\description{
  \code{meta} returns the metadata from an object; generally this is stored in the object's attributes. 
}

\usage{
meta(x)

\method{plot}{dtLife}(x, mfrow = NULL, ylab = NULL, main = NULL, 
                      cvar = NULL, type = "p", ...)
}

\arguments{
  \item{x }{The object whose metadata have to be extracted or that is to be plotted.
  }
  \item{mfrow, ylab, main, type }{arguments passed to the plotting functions.
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

meta(Wad_watertemp_LR)
meta(Marsdiep)

meta(Wad_weather)

# subsetting keeps the attributes
meta(Wad_weather[, 1:3])
meta(Wad_weather[1:3, ])
meta(subset(Wad_weather, subset = station == 208))

# but not if it is a matrix

names(meta(as.matrix(Wad_weather[, 1:3])))

# The plotting function makes use of the attributes
plot(Wad_weather, type = "l", cvar = "station")
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
