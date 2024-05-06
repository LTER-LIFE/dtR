\name{Grainsize}
\alias{Wad_sediment}
\title{
  Sediment grain size data
}

\usage{
data(Wad_sediment)
}


\value{
  a data.frame with the following columns: 
\itemize{
  \item \code{station}, a unique identifier for a station
  \item \code{longitude}, [dgE]
  \item \code{latitude}, [dgN]
  \item \code{mdGrain}, the median grain size, [micrometer]
  \item \code{gravel}, fraction of sediment (in weight) > 2 mm
  \item \code{silt}, fraction of sediment (in weight) < 64 micrometer
  \item \code{sand}, fraction of sediment (in weight) inbetween 64 micrometer and 2 mm
  }
}

\description{
 Dataset obtained from rijkswaterstaat, as transformed to netcdf in
 tranformation to NetCDF: https://opendap.deltares.nl/thredds/catalog/opendap/rijkswaterstaat/sedimentatlas_waddenzee/catalog.html.  
 
From the cumulative data (cumphi) versus grain size, the mdGrain was estimated as the median grain size, i.e. where half is larger and smaller, in terms of weight.

Gravel was estimated as the weight fration > 2000 micrometer;
Silt was estimated as the weight fration < 64 micrometer;
Sand was estimated as the weight fration < 2000 micrometer - silt fraction
}

\note{
a suitable aspect ratio can be estimated as:
1/cos((mean(latitude) * pi)/180)
}
\examples{

## =============================================================
## Plot the median grain size
## =============================================================

mf <- par(mfrow=c(2,2))

asp <- 1/cos((mean(Wad_sediment$latitude)*pi)/180)
with(Wad_sediment, points2D(longitude, latitude, colvar=mdGrain, 
       log="c", clab="um", main="Median grain size",
       pch=".", las=1, cex=2, asp=asp))         
with(Wad_sediment, points2D(longitude, latitude, colvar=silt, 
       log="c", clab="", main="percent silt",
       pch=".", las=1, cex=2, asp=asp))         
with(Wad_sediment, points2D(longitude, latitude, colvar=sand, 
       log="", clab="", main="percent sand",
       pch=".", las=1, cex=2, asp=asp)) 
with(Wad_sediment, points2D(longitude, latitude, colvar=sand, 
       log="", clab="", main="percent sand",
       xlim=c(4.5, 5), ylim=c(52.8, 53.4),
       pch=".", las=1, cex=2, asp=asp))        
par(mfrow=mf) 

# ---------------------------------------------------
# Plot the bathymetry with Grainsize points on top
# ---------------------------------------------------

plotBathymetry(Wad_depth, 
               pts    = Wad_sediment[,c("longitude", "latitude", "mdGrain")], 
               ptlist = list(cex=0.5), 
               log    = "c", clab = "um",
               main   = "Grain size")


## =============================================================
## find corresponding depths for grain size data
## =============================================================

# ----------------------------------
# use only data from the Waddensea
# ----------------------------------

# bathymetry
batmasked        <- mask_bat(Wad_depth,
                             shape = Wad_shape)  
batmasked$depth  <- Wad_depth$depth * batmasked$mask 

# grain size
grainmasked <- na.omit(mask_xy(Wad_sediment[ ,c("longitude", "latitude")],
                               shape = Wad_shape))

system.time(
  Grain.ext <- map_xy(input.x   = batmasked$longitude, 
                      input.y   = batmasked$latitude, 
                      input.2D  = batmasked$depth,
                      output.xy = grainmasked)
)

Grain.ext <- merge(Grain.ext, Wad_sediment)
with(Grain.ext, plot(z, mdGrain))

# ----------------------------------
# full dataset
# ----------------------------------

system.time(
  Grain.all <- map_xy(
                  input.x   = Wad_depth$longitude, 
                  input.y   = Wad_depth$latitude, 
                  input.2D  = Wad_depth$depth,
                  output.xy = Wad_sediment[, c("longitude", "latitude")])
)

GG <- mask_xy(Grain.all[ ,c("longitude", "latitude")],
              shape = Wad_shape)
              
GG$mask[is.na(GG$mask)] <- 0

Grain.all <- merge(Grain.all, Wad_sediment)
Grain.all <- merge(Grain.all, GG)

with(Grain.all, plot(z, mdGrain, col=mask+1, pch=".", cex=2,
     main="sediment-depth relationships"))
legend("topleft", col=1:2, pch=18, legend=c("northsea", "waddensea"))

}
\seealso{
TO BE FILLED
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

