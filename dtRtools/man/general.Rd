\name{general tools}
\alias{meta}
\title{
  extract metadata from an object - keeps the attributes for [] and \code{subset}
}

\description{
  \code{meta} returns the metadata from an object; generally this is stored in the object's attributes. 
}

\usage{
meta(x)
}

\arguments{
  \item{x }{The object whose metadata have to be extracted .
  }
}

\value{
  Function \code{meta} returns a \code{data.frame} with the objects attributes, except the names, row.names and class.
  
  Subsetting an object keeps the metadata.
}

\examples{

# ------------------------
# the metadata
# ------------------------

meta(Wad_watertemp_LR)
meta(Marsdiep)

meta(Wad_weather)

# subsetting keeps the attributes
meta(Wad_weather[, 1:3])
meta(Wad_weather[1:3, ])

meta(subset(Wad_weather, subset = station == 208))

# but not if it is a matrix

names(meta(as.matrix(Wad_weather[, 1:3])))

# ------------------------
# plotting
# ------------------------

# The plotting function makes use of the attributes

# The default uses the description of the variables as title (too long here)
plot(subset(Wad_weather, 
            subset = station == 235,
            select = c("datetime", "radiation", "temperature")),
     type = "l")

# how to color the datapoints can be specified
plot(Wad_weather, type = "l", cvar = "station", 
     use_names = TRUE)    # uses names for labelling rather than description
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
