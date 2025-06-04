\name{meta}
\alias{meta}

\title{
  extract metadata from an object
}

\description{
  \code{meta} returns the metadata from an object; generally this is stored in the object's attributes. 
}

\usage{
meta(x)
}

\arguments{
  \item{x }{The object whose metdata have to be extracted.
  }
}

\value{
  Function \code{meta} returns a \code{data.frame} with the objects attributes, except the names, row.names and class.
}

\examples{

meta(Wad_watertemp_LR)
meta(Marsdiep)
meta(Wad_weather)

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
