\name{interpolate_xyt}
\alias{interpolate_xyt}
\title{
  Interpolates spatial timeseries data to spatial points
}

\description{

  \code{interpolate_xyt} interpolates spatially-distinct timeseries to other spatial points. 
  The input dataset is typically more sparse than the output. 
  
  For each output position, it finds the three spatially closest data sets in the input data. 
  For the selected spatial input points, the time series values are interpolated to the required time output points. 
  These three data sets are then averaged, using inverse distance weighing.
  
  
  \code{interpolate_dense} interpolates spatially-distinct timeseries to other spatial points, for very dense input datasets. 
  
  The input points are mapped in the grid of the output points, and then for each of the grid cells, the average value is estimated. This requires the output grid to be orthogonal. 

}

\usage{
interpolate_xyt (input_xytv, 
             output_xy, output_x, output_y,
             output_t, 
             nmean = 3, asp = c("geographic", "mean", "none"), 
             rule = 2, ID = NULL)
}

\arguments{
  \item{input_xytv }{Data.frame with input data, 4 columns: longitude (or x-values), 
  latitude (or y-values), time, and value of time series variable. 
  }
  \item{output_xy }{Two-columned matrix with the latitude (or x-values), and 
  longitude (or y-values) of the output data.
  }
  \item{output_x, output_y  }{vectors with \code{x} and \code{y} values, determining the x and y values to which the input data will be mapped. 
  Ignored if \code{output_xy} is given a value.
  }
  \item{output_t }{Vector with the output time-values. 
  }
  \item{nmean }{The number of points over which the new value will be averaged. 
  }
 \item{asp }{The y/x aspect ratio for weighing y versus x. If "geographic" then it is assumed that the spatial points are geographic (spherical) coordinates (latitude, longitude) 
 and an aspect ratio for distances will be used (based on \link{aspect_coord}). 
 If "mean" then the distances will be weighed according to the mean values of x and y; the aspect ratio will be mean(x)/mean(y).  If "none", there will be no rescaling of y- and x- distances (and asp will be set = 1). Also allowed is to give asp a number.
  }
 \item{ID }{Vector with unique identifiers for the output; should be of same length as nrows of output_xy. Used for labelling the output; If absent output will be labeled "st1", "st2", ...
   Should be of the same type as the input times (3rd column of \code{input})  
   }
   \item{rule }{how to use approx for generating the time series data. 
  }
 
}
\value{
    Function \code{interpolate_xy} returns a \code{data.frame} of class \code{dtLife} with:
    
    \itemize{
     \item the first two columns contain the longitude and latitude of the spatial points as input by argument \code{output_xy}.
     \item other columns contain the interpolated values for each of the time points in \code{output_t}. These columns are labeled according to the datetime, cast as a string.  The actual time value, as input by  \code{output_t} is stored in attribute "datetime" of the returned object.
    }    
    
    The positions of the output stations are stored in attribute "stations" of the returned object; the time values are stored in attribute "datetime".

}
\note{
For this to work it is necessary that the input consists of spatial points for which a time series of values is defined.
Completely disjoint (x, y, t) points will not be properly interpolated!

}

\examples{

## =============================================================
## testing on the volcano data set
## =============================================================

# Volcano is a matrix - cast it in x,y,v format

xy          <- expand.grid(x = 1:nrow(volcano), 
                           y = 1:ncol(volcano))

volcano.xyz <- data.frame(xy, 
                          v = as.vector(volcano))

# create a 'time series data set'
volcano.xytz <- rbind(
    data.frame(volcano.xyz[,1:2], time = 1, d = volcano.xyz[,3]    *1), 
    data.frame(volcano.xyz[,1:2], time = 2, d = volcano.xyz[,3]  *0.5), 
    data.frame(volcano.xyz[,1:2], time = 3, d = volcano.xyz[,3]*0.125))

head(volcano.xytz)
par(mfrow = c(2,2))

# interpolate in time, using output_x and output_y
v1 <- interpolate_xyt(input_xytv = volcano.xytz,
                      output_x   = seq(1, nrow(volcano), by = 0.5),
                      output_y   = seq(1, ncol(volcano), by = 0.5),
                      output_t   = seq(1, 3, by = 0.5))

points2D(v1[,1], v1[,2], colvar = v1[,3], pch = 18)

points2D(v1[,1], v1[,2], colvar = v1[,4], pch = 18)

# interpolate in time, using output_xy
v2 <- interpolate_xyt(input_xytv = volcano.xytz,
                      output_xy  = xy,
                      output_t   = 1:3)

points2D(v2[,1], v2[,2], colvar = v2[,3], pch = 18)
points2D(v2[,1], v2[,2], colvar = v2[,4], pch = 18)

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

# 4 spatial points and 10 time values to interpolate to
output_xy <- cbind(x = runif(n = 4, min = 1, max = 2), 
                   y = runif(n = 4, min = 2, max = 3))
output_t  <- seq(4, 8, length.out = 10)

# show inputs
head(input_xytv)
head(output_xy)
head(output_t)

# Interpolate to 2D grid
output_v <- interpolate_xyt(
                        input_xytv = input_xytv, 
                        output_xy  = output_xy, 
                        output_t   = output_t  )
# show results
head(output_v)

attributes(output_v)$stations
attributes(output_v)$datetime

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

# output positions: 100 stations to be mapped
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

times <- attributes(output_v)$datetime
data  <- t(output_v[,-(1:2)])   # remove first two columns, and take transpose
matplot(x    = times, 
        y    = data,
        main = "interpolated", 
        type = "l", lty = 1)

par(mfrow=pm)

## =============================================================
## 2D-time interpolation of weather data (timeseries)
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
output_wind <- interpolate_xyt(
                      input_xytv  = Input, 
                      output_xy   = output_xy, 
                      output_t    = output_t)
)

### NOTE: The output is a matrix; first two columns contain the coordinates

# extract the actual dates first
date <- meta(output_wind)$datetime

# inverse distances, based on 6 points
system.time(
output_wind.1 <- interpolate_xyt(
                      input_xytv  = Input, 
                      output_xy   = output_xy, 
                      nmean       = 6,   # mean of all points
                      output_t    = output_t)
)

DD <- subset(Input, subset = date == "2021-06-01")

par(mfrow=c(2, 2))

image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = output_wind[,8]),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "inverse distance - 3 points")
          
image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = output_wind.1[,8]),
        xlim = range(DD$longitude),
        ylim = range(DD$latitude),
        clim = range(DD$windspeed),
        main = "inverse distance - 6 points")

points2D(x      = DD$longitude, 
         y      = DD$latitude, 
         colvar = DD$windspeed, 
         pch=18, cex=2)

matplot(t(output_wind[,-(1:2)]), type="l", lty=1)

timeseries <- data.frame(datetime = attributes(output_wind)$datetime, 
                         t(output_wind[, -(1:2)]))
                    
meanTS <- average_timeseries(timeseries, 
                             datetime = "datetime", 
                             avgOver  = "year", 
                             value    = colnames(timeseries)[-1])
summary(unlist(meanTS[, -1]))                             


## =============================================================
## 2D interpolation of spatial data
## =============================================================

outx <- seq(  4.57,    7.4, by = 0.005)
outy <- seq(52.86, 53.67, by = 0.005)

outxy <- expand.grid( outx, 
                      outy)
                     
output_xy <- na.omit(mask_shape(coordinates = outxy, 
                                shape       = Shape))

system.time(Sed2D <- 
     interpolate_xy(input_xyv = Sediment, 
                    output_xy = output_xy))

asp <- aspect_coord(Sed2D[,2])

pm <- par(mfrow = c(2, 2))
with(Sediment, points2D(longitude, latitude, 
                        colvar = mdGrain, 
                        main   = "original data",
                        pch = ".", cex = 3,
                        asp = asp))
plot(Shape, add=TRUE)

with(Sed2D, points2D(longitude, latitude, 
                     colvar = mdGrain, 
                     main   = "interpolated data",
                     pch = ".", cex = 3,
                     asp = asp))
plot(Shape, add=TRUE)

with(Sed2D, points2D(longitude, latitude, 
                     colvar = mdGrain, 
                     main   = "interpolated data_2",
                     pch = ".", cex = 3,
                     asp = asp))

plot(Shape, add=TRUE)

## =============================================================
## Water temperature
## =============================================================

head(Wad_watertemp_LR, n = c(3, 3))

WT <- dt_tolong(Wad_watertemp_LR)

Stations <- meta(Wad_watertemp_LR)$stations
Stations

WT <- merge(WT, Stations, by = "station")
WT <- WT[,c("longitude", "latitude", "datetime", "value")]

xyout <- expand.grid(x = seq(4.75,  7.2, length.out = 20),
                     y = seq(53.0, 53.5, length.out = 20))

tout <- seq(from = as.POSIXct("2021-02-01 00:00:00 CET"), 
            to   = as.POSIXct("2021-03-01 00:00:00 CET"), 
            by   = 86400*7)
            
WTint <- interpolate_xyt(input_xytv = WT,
                         output_xy  = xyout,
                         output_t   = tout)

points2D(WTint[,1], WTint[,2], colvar= WTint[,3], pch = 18)

par(mfrow=pm)
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
  \code{\link{interpolate_xy}} for a comparable method as \code{interpolate_xyt}.
  \code{\link{Sediment}} data set with sediment composition.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
