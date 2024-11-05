\name{interpolate_xyt}
\alias{interpolate_xyt}
\alias{interpolate_xy}
\alias{match_timeseries}
\title{
  Interpolates spatial (timeseries) data -  to be used for spatially sparse datasets.
}

\description{
  \code{interpolate_xy} interpolates spatially-distinct values to other spatial points. 

  \code{interpolate_xyt} interpolates spatially-distinct timeseries to other spatial points. 
  
  For each output position, it finds the three spatially closest data sets in the input data. 
  For the selected spatial input points, the time series values are interpolated to the required time output points. 
  These three data sets are then averaged, using inverse distance weighing.

  \code{match_timeseries} matches a timeseries to a dataset that has an "ID" in its metadata. 
}

\usage{
interpolate_xy (input_xyv, output_xy, barycentric = FALSE,
                nmean = 3, geographic = TRUE, 
                rule = 2, all.output=FALSE)

interpolate_xyt (input_xytv, output_xy, output_t, barycentric = FALSE,
                nmean = 3, ID=NULL, geographic = TRUE, 
                rule = 2, all.output=FALSE)

match_timeseries(input_xytv, data)   

}

\arguments{
  \item{input_xyv }{Data.frame with input data, 3 columns: longitude (or x-values), 
  latitude (or y-values), and value of time series variable. 
  }
  \item{input_xytv }{Data.frame with time series input data, 4 columns: longitude (or x-values), 
  latitude (or y-values), time, value of time series variable. 
  \code{time} is either numeric, or of class Date, posIXct. 
  }
 \item{output_xy }{Two-columned matrix with the latitude (or x-values), and 
  longitude (or y-values) of the output data.
  }
 \item{output_t }{Vector with the times to which the time series should be interpolated.
   Should be of the same type as the input times (3rd column of \code{input})  
   }
 \item{barycentric }{If TRUE, then the interpolation will be done by barocentric coefficients. The default is to use inverse distance approximation.
  }
  \item{nmean }{The number of points over which the new value will be averaged. 
  }
  \item{geographic }{If TRUE then it is assumed that the spatial points are geographic coordinates (latitude, longitude) and an aspect ratio for distances will be used.
  }
 \item{all.output }{If TRUE then information on how the interpolation is done will be written to the attributes.
  }
 \item{ID }{Vector with unique identifiers for the output; should be of same length as nrows of output_xy. Used for labelling the output; If absent output will be labeled "st1", "st2", ...
   Should be of the same type as the input times (3rd column of \code{input})  
   }
  \item{data }{dataset that has an ID in its attributes, from which the output_xy, and ID can be extracted.
  }
   \item{rule }{how to use approx for generating the time series data. An integer (of length 1 or 2) describing how interpolation is to take place outside the interval [min(x), max(x)]. If rule is 1 then NAs are returned for such points and if it is 2, the value at the closest data extreme is used. Use, e.g., rule = 2:1, if the left and right side extrapolation should differ.
  }
 
}
\value{
    Function \code{interpolate_xyt} returns a \code{data.frame} with:
    
    \itemize{
     \item the first column, contains the times, as in \code{output_t},
     \item other columns contain the interpolated values for each of the stations in \code{output_xy}. These columns are either labeled "st1",... or are give the names as provided by argument \code{ID}.
    }    
    The positions of these stations are stored in attribute "stations" of the returned object.
  
}
\note{
For this to work it is necessary that the input consists of spatial points for which a time series of values is defined.

Completely disjoint (x, y, t) points will not be properly interpolated!
}

\examples{

## =============================================================
## 2D-time interpolation of small data set
## =============================================================

# a small data set

times <- seq(1, 10, length.out=1000)

d1   <- cbind(x     = 1, 
              y     = 2,   
              times = times, 
              value = sin(times))
d2   <- cbind(x     = 1, 
              y     = 3,   
              times = times, 
              value = sin(times+0.5))
d3   <- cbind(x     = 2, 
              y     = 2,   
              times = times, 
              value = sin(times-0.5))
d4   <- cbind(x     = 2, 
              y     = 3,   
              times = times, 
              value = sin(times-1))
d5   <- cbind(x     = 2.5, 
              y     = 2.5, 
              times = times, 
              value = sin(times+1))

# combine all to one input data set
input_xytv <- rbind(d1, d2, d3, d4, d5)

# 4 points and 10 time values to interpolate to

output_xy <- cbind(x = runif(n=4, min=1, max=2), 
                   y = runif(n=4, min=2, max=3))
output_t  <- seq(4, 8, length.out=10)

# Interpolate to 2D grid

output_v <- interpolate_xyt(input_xytv = input_xytv, 
                            output_xy  = output_xy, 
                            output_t   = output_t  )
# show results

head(input_xytv)

head(output_xy)
head(output_t)

head(output_v)

attributes(output_v)$stations

## =============================================================
## 2D-time interpolation of sine waves 
## =============================================================

pm <- par(mfrow=c(1,2))
times <- seq(from = 1, to = 10, length.out = 100)

# input: 100 stations with data
Input <- NULL

for (i in 1: 100)
  Input <- rbind(Input, 
                 cbind(x     = runif(1, 0, 10), 
                       y     = runif(1, 0, 10), 
                       times = times,
                       value = sin(times+runif(1, -1, 1)))
                )

# output psotions: 100 stations to be mapped
no        <- 100 
Output    <- cbind(x = runif(no, 1,  9), 
                   y = runif(no, 0, 10))

# output time: long timeseries
output_t  <- seq(from = 4, to = 8, length.out = 200)

system.time(
  output_v <- interpolate_xyt(input_xytv = Input, 
                              output_xy  = Output, 
                              output_t   = output_t)
)

# show results
with(as.data.frame(Input), 
   plot(x    = times, 
        y    = value, 
        main = "input data", 
        pch = ".", cex = 3))

matplot(x    = output_v[,1], 
        y    = output_v[,-1], 
        main = "intepolated", 
        type = "l", lty = 1)

par(mfrow=pm)

## =============================================================
## 2D interpolation of weather data
## =============================================================

Weather.stations <- attributes(Wad_weather)$stations
head(Wad_weather)

# Create daily averages of the Windspeed as input

Wad_weather$Day <- as.Date(Wad_weather$datetime)

mf <- par(mfrow=c(3,1))

WS <- aggregate(x   = Wad_weather$windspeed, 
                by  = list(Wad_weather$Day, Wad_weather$station), 
                FUN = mean)
colnames(WS) <- c("date", "station", "windspeed")               

Input <- merge(
         Weather.stations[ , c("station", "longitude", "latitude")],
         WS              [ , c("station", "date", "windspeed")]
              )
              
# remove first column                 
Input     <- Input[, -1] 

head(Input)   # 4 columns: x, y, t, value

# output: 400 stations, first day of the month

# x-and y values to output
nx   <- 50
ny   <- 50

outx <- seq(4.8,   5.5, length.out=nx)
outy <- seq(52.9, 53.4, length.out=ny)

output_xy <- expand.grid(x= outx,
                         y= outy)
head(output_xy)

# output time value
output_t <- as.Date(c("2021-01-01", "2021-02-01", 
                      "2021-03-01", "2021-04-01", 
                      "2021-05-01", "2021-06-01", 
                      "2021-07-01", "2021-08-01", 
                      "2021-09-01", "2021-10-01", 
                      "2021-11-01", "2021-12-01"))

# interpolate - inverse distances, based on 3 points
system.time(
output_wind <- interpolate_xyt(input_xytv  = Input, 
                               output_xy   = output_xy, 
                               barycentric = FALSE,
                               output_t    = output_t)
)

# inverse distances, based on 6 points
system.time(
output_wind.1 <- interpolate_xyt(input_xytv  = Input, 
                               output_xy   = output_xy, 
                               barycentric = FALSE,
                               nmean       = 6,   # mean of all points
                               output_t    = output_t)
)

# interpolate using barycentric interpolation
system.time(
output_wind.2 <- interpolate_xyt(input_xytv  = Input, 
                               output_xy   = output_xy, 
                               barycentric = TRUE,
                               output_t    = output_t)
)

# interpolate using barycentric interpolation - different triangels
system.time(
output_wind.3 <- interpolate_xyt(input_xytv  = Input, 
                               output_xy   = output_xy, 
                               barycentric = 2,
                               output_t    = output_t)
)

output_stations <- attributes(output_wind)$stations

par(mfrow=c(2, 3))

DD <- subset(Input, subset = date == output_wind[6, 1])

image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = as.numeric(output_wind[6,-1])),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "inverse distance - 3 points")
image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = as.numeric(output_wind.1[6,-1])),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "inverse distance - 6 points")

image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = as.numeric(output_wind.2[6,-1])),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "barycentric interpolation - fixed triangle")
image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = as.numeric(output_wind.3[6,-1])),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "barycentric interpolation - closest triangle")

points2D(x      = DD$longitude, 
         y      = DD$latitude, 
         colvar = DD$windspeed, 
         pch=18, cex=2)

matplot(output_wind[,-1], type="l", lty=1)
# matplot(output_wind.2[,-1], type="l", lty=1)
}
\references{
 Jassby A.D. and T. Platt, 1976. Mathematical formulation of the relationship 
 between photosynthesis and light for phytoplankton. Limnol. Oceanogr. 21:540-547.

 Platt T. and A.D. Jasby, 1976. The relationship between photosynthesis and 
 light for natural assemblages of coastal marine phytoplankton. J. Phycol. 12:421-430.

 Platt T., Gallegos C.L. and W.G. Harrison, 1980. Photoinhibition of photosynthesis 
 in natural assemblages of marine phytoplankton. J. Mar. Res., 38, 687-701
 
 Eilers P.H.C. and J.C.H. Peeters, 1988. A model for the relationship between 
 light intensity and the rate of photosynthesis in phytoplankton. Ecol. Modell., 42, 199-215.

 Webb W.L, Newton M. and D. Starr, 1974. Carbon dioxide exchange of Alnus rubra.
 A mathematical model. Oecologia, Berlin 17, 281-291.
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
