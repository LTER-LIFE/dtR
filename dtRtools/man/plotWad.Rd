\name{plot_bathymetry}
\alias{plot_bathymetry}
\title{
  Simple plot function for bathymetric data.
}

\usage{
plot_bathymetry(bat, pts = NULL, type = NULL, ptlist = NULL, 
                negativeDepth = FALSE, ...)
}

\value{
returns nothing
}
\arguments{
  \item{bat }{A list with the bathymetric information. Should at least contain the vectors named \code{longitude} and \code{latitude} and the matrix named \code{depth} (of dimension \code{(length(longtiude), length(latitude))}. 
  }
  \item{pts }{A data.frame with 3 columns, latitude, longitude, value, describing the points whose values should be plotted.  
  }
  \item{ptlist }{A list with settings for the points.  
  }
  \item{type }{How the bathymetry will be depicted, one of "image", "contour" or "both". 
  The default is to use a colored image in the case of absence of pts, or contours in the case of presence of pts.
  }
  \item{negativeDepth }{If \code{FALSE}, then depth will be positive numbers; if \code{FALSE} then depth will be depicted as negative numbers. The latter is very useful if the image is shaded. 
  }
  \item{... }{Any other arguments passed to the plotting functions.
}
}

\description{
 Simple plotting of dataset obtained from EMODnet. 
 On behalf of the EMODnet project, http://www.emodnet-bathymetry.eu/.
}

\examples{

# ---------------------------------------------------
# Plot the bathymetry only
# ---------------------------------------------------

plot_bathymetry(Marsdiep)
plot_bathymetry(Marsdiep, type = "both")
plot_bathymetry(Marsdiep, type = "contour")

# add shade to emphasise topography - use negative depths.
plot_bathymetry(Marsdiep, 
               shade = 0.15, 
               NAcol = grey(0.7),
               col   = ramp.col(col=c("darkblue", "lightblue", "white")),
               negativeDepth = TRUE)

# ---------------------------------------------------
# Plot the bathymetry with points on top
# ---------------------------------------------------

plot_bathymetry(Marsdiep, 
               pts    = RWSstations[,3:2], 
               ptlist = list(cex=3), 
               main   = "Position of RWS stations")

}
\seealso{
  \code{\link{read_bathymetry}} to read bathymetric data.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}

