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

meta(WadTempLR)
meta(Marsdiep)
meta(Weather2021)

}

\seealso{
  \code{\link{Marsdiep}} bathymetric data.
  \code{\link{WadTempLR}} data set with water temperature data from the Waddensea.
  \code{\link{Weather2021}} data set with weather information.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
