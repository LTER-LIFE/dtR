\name{Sediment data}
\alias{Sediment}
\title{
  Sediment grain size data for a part of the Waddensea
}

\usage{
data(Sediment)
}


\value{
  a data.frame with the following columns: 
\itemize{
  \item \code{station}, a unique identifier for a station
  \item \code{longitude}, [dgE]
  \item \code{latitude}, [dgN]
  \item \code{mdGrain}, the median grain size, [micrometer]  
 }
}

\description{
 Dataset obtained from rijkswaterstaat, as transformed to netcdf in
 
 https://opendap.deltares.nl/thredds/catalog/opendap/ 
 
 rijkswaterstaat/sedimentatlas_waddenzee/catalog.html.  
 
From the cumulative data (cumphi) versus grain size, the mdGrain was estimated as the median grain size, i.e. where half is larger and smaller, in terms of weight.

A subset of the data, W of 5.01 dg E and south of 53.1 dg N was taken to make this data set small.
}

\note{
a suitable aspect ratio can be estimated as:
1/cos((mean(latitude) * pi)/180)
}
\examples{

## =============================================================
## Plot the median grain size
## =============================================================

asp <- 1/cos((mean(Sediment$latitude)*pi)/180)

with(Sediment, 
  points2D(x      = longitude, 
           y      = latitude, 
           colvar = mdGrain, 
           log    = "c", 
           clab   = "um", 
           main   = "Median grain size",
           pch = ".", las = 1, cex = 2, asp = asp))         

# ---------------------------------------------------
# Plot the bathymetry with Grainsize points on top
# ---------------------------------------------------

plotBathymetry(Marsdiep, 
               pts    = Sediment[,c("longitude", "latitude", "mdGrain")], 
               ptlist = list(cex=1.5), 
               log    = "c", clab = "um",
               main   = "Grain size")

## =============================================================
## find corresponding depths for grain size data
## =============================================================

# Find the depths for the positions of the grain size data
system.time(
  Grain.all <- map_xy(
                  input_x   = Marsdiep$longitude, 
                  input_y   = Marsdiep$latitude, 
                  input_2D  = Marsdiep$depth,
                  output_xy = Sediment[, c("longitude", "latitude")])
)

Grain.all <- merge(Grain.all, Sediment)

# Grain.all has elements x, y, and v (for variable)

with(Grain.all, 
  plot(x    = v, 
       y    = mdGrain, pch = ".", cex = 2,
       main = "sediment-depth relationships"))
}
\seealso{
TO BE FILLED
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

