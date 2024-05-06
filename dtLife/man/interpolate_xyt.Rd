\name{interpolate_xyt}
\alias{interpolate_xyt}
\alias{match_timeseries}
\title{
  Interpolates spatial timeseries data -  to be used for spatially sparse datasets.
}

\description{
  \code{interpolate_xyt} interpolates spatially-distinct timeseries to other spatial points. 
  
  For each output position, it finds the three spatially closest data sets in the input data. 
  For the selected spatial input points, the time series values are interpolated to the required time output points. 
  These three data sets are then averaged, using inverse distance weighing.

  \code{match_timeseries} matches a timeseries to a dataset that has an "ID" in its metadata. 
}

\usage{
interpolate_xyt (input.xytv, output.xy, output.t, ID=NULL, wide=TRUE)
match_timeseries(input.xytv, data)   

}

\arguments{
  \item{input.xytv }{Data.frame with time series input data, 4 columns: longitude (or x-values), 
  lagitude (or y-values), time, value of time series varaiable. 
  \code{time} is either numeric, or of class Date, posIXct. 
  }
 \item{output.xy }{Two-columned matrix with the latitude (or x-values), and 
  longitude (or y-values) of the output data.
  }
 \item{output.t }{Vector with the times to which the time series should be interpolated.
   Should be of the same type as the input times (3rd column of \code{input})  
   }
 \item{ID }{Vector with unique identifiers for the output; should be of same length as nrows of output.xy. Used for labelling the output; If absent output will be labeled "st1", "st2", ...
   Should be of the same type as the input times (3rd column of \code{input})  
   }
  \item{data }{dataset that has an ID in its attributes, from which the output.xy, and ID can be extracted.
  }
 \item{wide }{If TRUE, then the output will be in wide format, i.e. one column per point interpolated, if FALSE then the output will be in long format, i.e. 4 columns with (x, y, time, value).
  }
}
\value{
  if argument \code{wide} is \code{TRUE}, then:
  
    Function \code{interpolate_xyt} returns a \code{data.frame} with:
    
    \itemize{
     \item the first column, contains the times, as in \code{output.t},
     \item other columns contain the interpolated values for each of the stations in \code{output.xy}. These columns are either labeled "st1",... or are give the names as provided by argument \code{ID}.
    }    
    The positions of these stations are stored in attribute "stations" of the returned object.
  
  if argument \code{wide} is \code{FALSE}, then:
  
  Function \code{interpolate_xyt} returns a \code{data.frame} with 5 columns: 
  
  \itemize{
     \item \code{ID}, the unique identifier or station name
     \item \code{output.t}, the times, as in argument \code{output.t},
     \item the \code{x}- and \code{y}-coordinates for each of the stations, as input by argument \code{output.xy},
     \item \code{res}, the interpolated value.
  }
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
Input.xytv <- rbind(d1, d2, d3, d4, d5)

# 4 points and 10 time values to interpolate to

Output.xy <- cbind(x = runif(n=4, min=1, max=2), 
                   y = runif(n=4, min=2, max=3))
Output.t  <- seq(4, 8, length.out=10)

# Interpolate to 2D grid

Output.v <- interpolate_xyt(input.xytv = Input.xytv, 
                            output.xy  = Output.xy, 
                            output.t   = Output.t  )
# show results

head(Input.xytv)

head(Output.xy)
head(Output.t)

head(Output.v)

attributes(Output.v)$stations

# same but output now in long format:

Output.v.long <- interpolate_xyt(input.xytv = Input.xytv, 
                                 output.xy  = Output.xy , 
                                 output.t   = Output.t  , 
                                 wide       = FALSE)

head(Output.v.long)

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
Output.t  <- seq(from = 4, to = 8, length.out = 200)

system.time(
  Output.v <- interpolate_xyt(input.xytv = Input, 
                              output.xy  = Output, 
                              output.t   = Output.t)
)

# show results
with(as.data.frame(Input), 
   plot(x    = times, 
        y    = value, 
        main = "input data", 
        pch = ".", cex = 3))

matplot(x    = Output.v[,1], 
        y    = Output.v[,-1], 
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
         Weather.stations[ , c("station", "latitude", "longitude")],
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

Output.xy <- expand.grid(x= outx,
                         y= outy)
head(Output.xy)

# output time value
Output.t <- as.Date(c("2021-01-01", "2021-02-01", 
                      "2021-03-01", "2021-04-01", 
                      "2021-05-01", "2021-06-01", 
                      "2021-07-01", "2021-08-01", 
                      "2021-09-01", "2021-10-01", 
                      "2021-11-01", "2021-12-01"))

# interpolate
Output.wind <- interpolate_xyt(input.xytv  = Input, 
                               output.xy   = Output.xy, 
                               output.t    = Output.t)
                               
Output.stations <- attributes(Output.wind)$stations

par(mfrow=c(3,3))
image2D(x = outx, y = outy, 
        z = matrix(nrow = nx, ncol = ny, 
                   data = as.numeric(Output.wind[6,-1]), 
                   byrow = TRUE))

DD <- subset(Input, subset = date == Output.wind[6,1])
points2D(x      = DD$latitude, 
         y      = DD$longitude, 
         colvar = DD$windspeed, 
         cex = 4, pch = 18)

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
