\name{Bathymetry}
\alias{Bathymetry}
\alias{WadDepth}
\alias{WadShape}

\title{
  Coarse bathymetry of the Dutch Waddensea area wadensea shape file.
}

\usage{
data(WadDepth)
}


\value{
WadDepth is a \code{list} with the following items:

\itemize{
  \item \code{longitude}, longitude at the middle of the grid in x [dgE], length=693.
  \item \code{latitude}, latitude at the middle of the grid in y [dgN], length=259.
  \item \code{depth}, sea-floor depth below NAP (normaal amsterdams peil) [m], (693, 258).
  \item \code{source}, reference to the website from which the data were extracted.
  \item \code{asp}, the aspect ratio that can be used for plotting.
  \item \code{contours}, a list with contours that can be used for plotting.
  \item \code{coast}, a data.frame with the coastline.
  }
}

\description{
 Dataset obtained from EMODnet. 
 The data were downloaded from  http://www.emodnet-bathymetry.eu/.
 
 They were created in the EMODnet project from data obtained from https://waterinfo.rws.nl/  and https://www.seadatanet.org/
 
 The water depth are estimated as - elevation, where elevation is defined as:
 \itemize{
   \item sdn_P02_urn: SDN:P02::MBAN
   \item sdn_parameter_name: Sea-floor height (above Lowest Astronomical Tide datum) {bathymetric height}
   \item sdn_parameter_urn: SDN:P01::HGHTALAT
   \item sdn_uom_name: Metres
   \item sdn_uom_urn: SDN:P06::ULAA
   \item standard_name: geoid_height_above_reference_ellipsoid
  }
  
  grid_mapping_epsg_code: EPSG:4326 (WGS 84)
  
  The original data from EMODnet bathymetry were expressed in Lowest Astronomical Tide datum (LAT)

  For consistency with the water height data, they were converted to NAP using a data set from 
https://english.defensie.nl/topics/hydrography/documents/applications/2020/06/12/nllat2018)

  For \code{WadDepth}, the resolution of the original dataset was reduced 
  with a factor 4.
  
}

\examples{

# ---------------------------------------------------
# Plot the high-resolution data
# ---------------------------------------------------

mf <- par(mfrow=c(2,1))

with(WadDepth, {
   image2D(longitude, latitude, z=depth, asp=asp, las=1)
   invisible(lapply(contours, FUN=lines))
})

with(Marsdiep, image2D(longitude, latitude, z=depth, asp=asp, las=1))
with(Marsdiep, invisible(lapply(contours, FUN=lines, col="white")))

par(mfrow=mf)   

# ---------------------------------------------------
# Plot the bathymetry only
# ---------------------------------------------------

plotBathymetry(WadDepth)
plotBathymetry(WadDepth, type="both")
plotBathymetry(WadDepth, type="contour")

# add shade to emphasise topography - use negative depths.
plotBathymetry(WadDepth, 
               shade = 0.15, 
               NAcol = grey(0.7),
               col   = ramp.col(col=c("darkblue", "lightblue", "white")),
               negativeDepth=TRUE)

# ---------------------------------------------------
# Plot the bathymetry with points on top
# ---------------------------------------------------

plotBathymetry(WadDepth, 
               pts    = WadGrain[,-1], 
               ptlist = list(cex=0.5), 
               log    = "c", clab="um",
               main   = "Grain size")

}
\seealso{
  \code{\link{WadGrain}} data set with sediment composition.
  \code{\link{Water2021}} data set with water temperature and height.
  \code{\link{Weather2021}} data set with weather information.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

