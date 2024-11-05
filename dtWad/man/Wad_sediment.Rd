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
a suitable aspect ratio for plotting can be estimated as:
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

# ---------------------------------------------------
# Create image rather than points
# ---------------------------------------------------

# this is requested at too high resolution - a lot of NA, where there 
# no values in the grid cell
Mdgrain2D <- map_dense(
       input_xyv = Wad_sediment[,c("longitude", "latitude", "mdGrain")],
       output_x = seq(from       = min(Wad_sediment$longitude), 
                      to         = max(Wad_sediment$longitude),
                      length.out = 500),
       output_y = seq(from       = min(Wad_sediment$latitude), 
                      to         = max(Wad_sediment$latitude),
                      length.out = 500))

# find NAs and select only the ones in the Wadden:

#NAs     <- which(is.na(Mdgrain2D$mdGrain), arr.ind = TRUE)
#NAcoord <- data.frame(longitude = Mdgrain2D$longitude[NAs[,1]], 
#                      latitude  = Mdgrain2D$latitude [NAs[,2]])
#NAs     <- mask_shape(coordinates = NAcoord, shape=Wad_shape)
#NAcoord <- subset(NAcoord, subset = ! is.na(NAs$mask))

#Mdext <- map_dense(input_xyv = Wad_sediment[,c("longitude", "latitude", "mdGrain")],
#                output_xy = NAcoord)

## =============================================================
## find corresponding depths for grain size data
## =============================================================

# ----------------------------------
# use only data from the Waddensea
# ----------------------------------

# bathymetry
batmasked        <- mask_shape(longitude = Wad_depth$longitude,
                               latitude  = Wad_depth$latitude,
                               shape     = Wad_shape)  
batmasked$depth  <- Wad_depth$depth * batmasked$mask 

# grain size
grainmasked <- na.omit(
              mask_shape(coordinates = Wad_sediment[ ,c("longitude", "latitude")],
                         shape       = Wad_shape))

system.time(
  Grain.ext <- map_dense(
                   input_x   = batmasked$longitude, 
                   input_y   = batmasked$latitude, 
                   input_2D  = batmasked$depth,
                   output_xy = grainmasked)
)

Grain.ext <- merge(Grain.ext, Wad_sediment)
with(Grain.ext, plot(v, mdGrain))  # v = value in Grain.ext

# ----------------------------------
# full dataset
# ----------------------------------

system.time(
  Grain.all <- map_dense(
                  input_x   = Wad_depth$longitude, 
                  input_y   = Wad_depth$latitude, 
                  input_2D  = Wad_depth$depth,
                  output_xy = Wad_sediment[, c("longitude", "latitude")])
)

GG <- mask_shape(coordinates = Grain.all[ ,c("longitude", "latitude")],
                 shape       = Wad_shape)
              
GG$mask[is.na(GG$mask)] <- 0

Grain.all <- merge(Grain.all, Wad_sediment)
Grain.all <- merge(Grain.all, GG)

with(Grain.all, plot(v, mdGrain, col=mask+1, pch=".", cex=2,
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

