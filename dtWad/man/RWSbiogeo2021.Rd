\name{WadBiogeo}
\alias{WadBiogeo}
\title{
  Biogeochemical data for all monitoring stations in the Waddensea, for 2021, at 0.5 to 1 monthly intervals.
}

\usage{
data(WadBiogeo)
}


\value{
WadBiogeo is a  data.frame of class \code{dtLife}, with the following items: 
\itemize{
  \item \code{station, longitude, latitude}, name and position (in EPSG 25831) of station - see \code{attributes(object)$stations}
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{Chl, HCO3, NH4, NO2, NO3, Ntot, O2, SPM, PO4, Ptot, SiO2, Ext, O2_percent, Secchi, pH, S}, the measured variables.
  }

 and the following attributes:
\itemize{
  \item \code{stations}, description of the stations, with their position
  \item \code{variables}, description of the measured variables, and their units
  \item \code{datasource}, data originator
  \item \code{file}, name of the files from which the data came
  \item \code{processing}, track of manipulations performed on the data
  \item \code{fun}, name of the function used to generate the data
  \item \code{format}, wide format
  }  
}

\description{
 Dataset obtained from the Rijkswaterstaat. 
 
 Downloaded from https://waterinfo.rws.nl/

}

\examples{

# ---------------------------------------------------
# Plot the position of the data
# ---------------------------------------------------

with(WadDepth, 
   image2D(x   = longitude, 
           y   = latitude, 
           z   = depth, clab = "m", 
           main = "Position of RWS biogeochemical monitoring stations"),
           asp  = asp, las=1, 
           colkey = list(length=0.5, width=0.5)
          )


stats  <- attributes(WadBiogeo)$stations
points(x = stats$longitude, 
       y = stats$latitude, 
       pch = 18, cex=3)

# ---------------------------------------------------
# A quick look at the data:
# ---------------------------------------------------

plot(WadBiogeo, pch=18)

# different color per station
plot(WadBiogeo, 
     pch = 18, cvar = "station")


# ---------------------------------------------------
# Select one variable
# ---------------------------------------------------
BG <- WadBiogeo

# use as.Date to get rid of the times - keep only sampling days
BG$datetime <- as.Date(BG$datetime) 

# reshape to wide format (columns = stations)
CHL <- reshape(BG[,c(1,4,5)], 
               direction = "wide", 
               timevar   = "station", 
               idvar     = "datetime")
CHL <- CHL[order(CHL$datetime),]

}
\seealso{
  \code{\link{Bathymetry}} bathymetric data.
  \code{\link{WadGrain}} data set with sediment composition.
  \code{\link{Weather2021}} data set with weather information.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

