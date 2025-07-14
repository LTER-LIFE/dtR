\name{Bathymetry}
\alias{Bathymetry}
\alias{Marsdiep}
\alias{Shape}

\title{
  Detailed bathymetry of the Marsdiep area in the Dutch Waddensea.
  
  A shape file
}

\usage{
data(Marsdiep)
data(Shape)
}


\value{
Dataset Marsdiep is a \code{list} with the following items:

\itemize{
  \item \code{longitude}, longitude at the middle of the grid in x [dgE], length = 693.
  \item \code{latitude}, latitude at the middle of the grid in y [dgN], length = 259.
  \item \code{depth}, sea-floor depth below NAP (normaal amsterdams peil) [m], (693, 258).
  \item \code{source}, reference to the website from which the data were extracted.
  \item \code{asp}, the aspect ratio that can be used for plotting.
  \item \code{contours}, a list with contours that can be used for plotting.
  }
}

\description{
 Dataset obtained from EMODnet. 
 The data were downloaded from  http://www.emodnet-bathymetry.eu/.
 
 They were created in the EMODnet project from data obtained from https://waterinfo.rws.nl/  and https://www.seadatanet.org/
 
 The water depth are estimated as - elevation, where elevation is defined as:
 \itemize{
   \item sdn_P02_urn: SDN:P02::MBAN
   \item sdn_parameter_name: Sea-floor height (above Lowest Astronomical Tide datum)   \item sdn_parameter_urn: SDN:P01::HGHTALAT
   \item sdn_uom_name: Metres
   \item sdn_uom_urn: SDN:P06::ULAA
   \item standard_name: geoid_height_above_reference_ellipsoid
  }
  
  grid_mapping_epsg_code: EPSG:4326 (WGS 84)
  
  The original data from EMODnet bathymetry were expressed in Lowest Astronomical Tide datum (LAT)

  For consistency with the water height data, they were converted to NAP using a data set from 
https://english.defensie.nl/topics/hydrography/documents/applications/2020/06/12/nllat2018)

  \code{Marsdiep} has the original resolution of the data.
}

\examples{

# ---------------------------------------------------
# print the metadata
# ---------------------------------------------------

meta(Marsdiep)

# ---------------------------------------------------
# Plot the high-resolution data
# ---------------------------------------------------

with(Marsdiep, 
    image2D(x   = longitude, 
            y   = latitude, 
            z   = depth, 
            asp = asp, las = 1))
with(Marsdiep, 
     invisible(lapply(contours, 
                      FUN = lines, col = "white")))


# =============================
# Convert this to xyz format
# =============================

ZZ <- with(Marsdiep, 
     interpolate_xy(input_x   = longitude, 
                    input_y   = latitude, 
                    input_2D  = depth, 
                    output_xy = expand.grid(x = longitude, 
                                            y = latitude)))

head(ZZ)

points2D(ZZ[,1], ZZ[,2], colvar = ZZ[,3])

}
\seealso{
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

