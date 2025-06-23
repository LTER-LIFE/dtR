\name{mask_shape}
\alias{mask_shape}
\title{
  Creates a mask for spatial data, indicating whether they located in an area.
}

\description{
  \code{mask_shape} creates a mask for (combinations of) longitude, latitude data.
}

\usage{
mask_shape(coordinates, longitude, latitude, shape, bat, EPSG = 4326)
}

\arguments{
  \item{longitude, latitude }{Vectors with coordinates. Will be ignored if \code{bat} is present.
  }
 \item{coordinates }{Two-columned matrix or data.frame with the longitude (1st column), and 
  latitude (second column) of the positions that need to be masked. Will be ignored if \code{bat} is present.
  }
 \item{bat }{A list with vectors \code{longitude} and \code{latitude}; if \code{asp}, and \code{processing} is present in the list, it will be updated.
 
 }  
 \item{shape }{The boundary shapefile to be used for the selection. 
  }
  \item{EPSG }{EPSG number of the coordinates. The default is the WGS84. 
  }
}

\value{
  When \code{bat}, or \code{longitude} and \code{latitude} is used for input, the function will return a \code{list} with 3 items: 
  
  \itemize{
     \item latitude, longitude: the vectors as input to the function,
     \item mask, a matrix with dimensions (length(latitude), length(longitude)), with the value \code{1} for coordinates falling within the shape, \code{NA} for coordinates outside of the shape. 
     \item asp, the aspect ratio that can be used for plotting
     \item processing, strings with the processing that occurred on the data set
  }

    When \code{coordinates} is input, \code{mask_shape}  returns a \code{data.frame} with 3 columns: 
  
  \itemize{
     \item first and second columns contain the latitude- and longitude-coordinates, as in \code{coordinates}.
     \item third column, called mask, contains the value \code{1} for coordinates falling within the shape, and \code{NA} for coordinates outside of the shape.
  }
}


\examples{

## =============================================================
## 2D mask
## =============================================================

 WW <- with(Marsdiep, mask_shape(longitude = longitude, 
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

 WW <- mask_shape(coordinates = Sediment[, c("longitude", "latitude")],
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
  \code{\link{Wad_watertemp_LR}} data set with (low resolution) water temperature from the Waddensea.
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
