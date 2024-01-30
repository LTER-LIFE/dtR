\name{mask_2D}
\alias{mask_2D}
\alias{mask_xy}
\alias{mask_bat}
\title{
  Creates a mask for spatial data, indicating whether they located in an area.
}

\description{
  \code{mask_2D} creates a mask for all combinations of longitude, latitude data (both vectors). 

  \code{mask_xy} creates a mask for given combinations of longitude, latitude data, input as a two-columned matrix or data.frame. 

  \code{mask_bat} applies a mask to a bathymetric data set }

\usage{
mask_2D(longitude, latitude, shape, EPSG=4326)
mask_xy(coordinates, shape, EPSG=4326)
mask_bat(bat, shape, EPSG=4326)
}

\arguments{
  \item{longitude, latitude }{Vectors with coordinates.
  }
 \item{coordinates }{Two-columned matrix or data.frame with the latitude (1st columne), and 
  longitude (second column) of the positions that need to be masked.
  }
  \item{bat }{A list with the bathymetric information. Should at least contain the vectors named \code{longitude} and \code{latitude} and the matrix named \code{depth} (of dimension \code{(length(longtiude), length(latitude))}. 
  }
 \item{shape }{The boundary shape to be used for the selection. 
  }
  \item{EPSG }{EPSG number of the coordinates. The default is the WGS84. 
  }
}

\value{
  Functions \code{mask_2D} and \code{mask_bat} return a \code{list} with 3 items: 
  
  \itemize{
     \item latitude, longitude: the vectors as input to the function,
     \item mask, a matrix with dimensions (length(latitude), length(longitude)), with the value \code{1} for coordinates falling within the shape, \code{NA} for coordinates outside of the shape. 
  }

  Function \code{mask_xy} returns a \code{data.frame} with 3 columns: 
  
  \itemize{
     \item first and second columns contain the latitude- and longitude-coordinates, input to the function
     \item third column, called mask, contains the value \code{1} for coordinates falling within the shape, \code{NA} for coordinates outside of the shape.
  }
}


\examples{

## =============================================================
## 2D mask
## =============================================================

 WW <- with(Marsdiep, mask_2D(longitude = longitude, 
                              latitude  = latitude,
                              shape     = Shape))

 pm <- par(mfrow=c(2,2))
 asp <- Marsdiep$asp
 
 with(Marsdiep, 
   image2D(x = longitude, 
           y = latitude, 
           z = depth, asp=asp,
           main = "Bathymetric file"))
 with(WW, 
   image2D(x = longitude, 
           y = latitude, 
           z = mask, asp=asp,
           main = "Waddensea mask"))
 with(Marsdiep, 
   image2D(x = longitude, 
           y = latitude, 
           z = depth*WW$mask, asp=asp, 
           main = "Masked bathymetry"))

 par(mfrow=pm)

## =============================================================
## xy mask
## =============================================================

 WW <- mask_xy(coordinates = Sediment[, c("longitude", "latitude")],
               shape       = Shape)

 pm <- par(mfrow=c(2,2))
 asp <- Marsdiep$asp
 with(Sediment, 
    points2D(x      = longitude, 
             y      = latitude, 
             colvar = mdGrain,
             main   = "Grain data, unmasked", 
             asp=asp, log="c", pch=".", cex=2))
 
 # make the NA in the mask = 1
 MM            <- WW$mask
 MM[is.na(MM)] <- 0
 
 with(WW, plot(x    = longitude, 
               y    = latitude, 
               col  = MM+1,
               main = "Waddensea mask", 
               pch = ".", cex = 2, asp = asp))

 legend("bottomright", pch = 18, col = 1:2, 
        legend=c("external", "Waddensea"))
 
 with(Sediment, 
    points2D(x      = longitude, 
             y      = latitude, 
             colvar = mdGrain*WW$mask,
             main   = "Grain data, masked", 
             log = "c", pch = ".", asp = asp, 
             cex = 2, NAcol = "black"))     

 par(mfrow=pm)
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
