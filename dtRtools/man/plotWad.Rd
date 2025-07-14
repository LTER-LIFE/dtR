\name{plot_bathymetry}
\alias{plot_bathymetry}
\alias{plot_dtmap}
\title{
  Simple plot functions for bathymetric or 2D data.
}

\usage{
plot_bathymetry(bat, type = NULL, pts = NULL, ptlist = NULL, 
                negativeDepth = FALSE, ...)

plot_dtmap     (bat, varname = NULL, type = NULL, 
                pts = NULL, ptlist = NULL,  ...)

}

\value{
returns nothing
}

\arguments{
  \item{bat }{A list with the mapping or bathymetric information. Should at least contain the vectors named \code{longitude} and \code{latitude} and, for \code{plot_bathymetry}, the matrix named \code{depth} (of dimension \code{(length(longitude), length(latitude))}. 
  }
  \item{varname }{The name of the 2D variable to be plotted.
  }
  \item{pts }{A data.frame with 3 columns, latitude, longitude, value, describing the points whose values should be plotted.  If \code{bat} is missing, but \code{pts} is passed, then the plot will be performed with points2D.
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
 Simple plotting of maps. 
 An aspect ratio is estimated based on the average of the latitudes, so this plot will not produce faithful maps if the latitude range is very large.
 
 \code{plot_dtmap} can be used to plot 2D matrices, or a set of points on coordinate systems.
 
 If the first argument to \code{plot_dtmap} is not named, it will be tested if it is a list (in which case, it will be assumed to represent \code{bat}), or a data.frame or matrix (assumed to represent \code{pts}).
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

# ---------------------------------------------------
# Plot other maps
# ---------------------------------------------------

# mask_shape returns a list with longitude, latitude, and mask
mask <- mask_shape(bat   = Marsdiep, 
                   shape = Shape)
                   
plot_dtmap(mask, colkey = FALSE, col = "grey")

# use plot_dtmap to plot points

plot_dtmap(Sediment, pch = 18,
           main   = "Sediment data")

}
\seealso{
  \code{\link{read_bathymetry}} to read bathymetric data.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}

