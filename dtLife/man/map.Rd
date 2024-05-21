\name{map_xy}
\alias{map_xy}
\alias{map_tx}

\title{
  \code{map_xy} haps 2D spatial data to xy or 2D format - to be used for spatially dense datasets.

  \code{map_tx} Maps 2D (space x time) data to xt or 2D format - to be used for dense datasets.
}

\description{
  \code{map_xy} maps 2D data to other spatial points. 
  
  \code{map_ts} maps (time, space) data to other (time, space) points. 
  
  
}
\usage{
map_xy(input.xyv,
       input.x, input.y, input.2D,
       output.xy, 
       output.x, output.y)
       
map_tx(input.txv,
       input.t, input.x, input.2D,
       output.tx,  
       output.t, output.x)
}

\arguments{
  \item{input.xyv }{Data.frame with  input data, 3 columns: latitude (or x-values), 
  longitude (or y-values), and variable value. 
  }
  \item{input.txv }{Data.frame with time series input data, 3 columns: time, position (or x-values), 
  and variable value. 
  }
  \item{input.x, input.y , input.t}{vectors with \code{x} and \code{y} input values (map_xy), or \code{t} and \code{x} input accompanying the 2D grid data. 
  }
  \item{input.2D }{2-dimensional data to interpolate from, a matrix of dimension \code{(length(x), length(y))} or \code{(length(t), length(x))}. 
  }
 \item{output.xy , output.tx }{Two-columned matrix with the latitude (or x-values), and 
  longitude (or y-values) of the output data (\code{map_xy}) or with time and x of the output data.
  }
  \item{output.x, output.y, output.t  }{vectors with \code{x} and \code{y} or \code{t} values, determining the 2D grid data that will be output. 
  }
}

\value{
  Functions \code{map_xy} and \code{map_tx} return either a \code{data.frame} or a \code{list}, depending on the output arguments.
  
  When the output data are specified as \code{output.xy} or \code{output.xt}, then the output will be a \code{data.frame} with 3 columns: 
  
  \itemize{
     \item first and second columns contain the x- and y-, or t - and x coordinates as defined in \code{output.xy} (\code{output.tx}),
     \item third column contains the interpolated values.
  }

  When the output data are specified as two vectors, \code{output.x} and  \code{output.y} (or \code{output.t}), then the output will be a \code{list} with 3 items: 
  
  \itemize{
     \item \code{x} and \code{y} or \code{t} and \code{x} are the two vectors defining the x- and y- (or t-) coordinates as input in \code{output.x} and \code{output.y} (or \code{output.t}),
     \item \code{v} will contain the 2-D matrix with interpolated values.
  }

}

\examples{

## =============================================================
## interpolation from xyv to xyv - small data set
## =============================================================

d1 <- cbind(x = 1,   y = 2,   value = 1)
d2 <- cbind(x = 1,   y = 3,   value = 2)
d3 <- cbind(x = 2,   y = 2,   value = 3)
d4 <- cbind(x = 2,   y = 3,   value = 4)
d5 <- cbind(x = 2.5, y = 2.5, value = 5)

Input.xyv <- rbind(d1, d2, d3, d4, d5)

Output.xy <- cbind(x = runif(n=4, min=1, max=2), 
                   y = runif(n=4, min=2, max=3))

Output.v <- map_xy(input.xyv = Input.xyv, 
                   output.xy = Output.xy)
Input.xyv
Output.xy
Output.v

## =============================================================
## interpolation from xyv to xyv - volcano
## =============================================================

pm <- par(mfrow=c(1,2))
# Volcano is a matrix - cast it in x,y,v format

xy          <- expand.grid(x = 1:nrow(volcano), 
                           y = 1:ncol(volcano))

volcano.xyz <- data.frame(xy, 
                          v = as.vector(volcano))

with(volcano.xyz, points2D(x, y, colvar = v, pch = 18, 
                           main = "volcano, original data"))

head(volcano.xyz)

# select 2000 random points to map

tomap <- data.frame(x = runif(min=1, max=nrow(volcano), n=2000),
                    y = runif(min=1, max=ncol(volcano), n=2000))

# interpolate
system.time(
 volcano.pts <- map_xy(input.xyv   = volcano.xyz, 
                       output.xy   = tomap)
)
head(volcano.pts)
attributes(volcano.pts)$processing

with(volcano.pts, points2D(x      = x, 
                           y      = y, 
                           colvar = v, pch = 18,
                           main   = "volcano, mapped data"))

par(mfrow=pm)

# -----------------------------------------------------------------------
# a test: interpolating points on themselves should be = original points
# -----------------------------------------------------------------------

# this is a bit slow
\dontrun{
system.time(
 volcano.tst <- map_xy(input.xyv   = volcano.xyz, 
                       output.xy   = volcano.xyz[,1:2])
)
summary(volcano.tst[,3]-volcano.xyz[,3])  # all 0

attributes(volcano.pts)$processing
}

## =============================================================
## from 2D to xyv interpolation - volcano
## =============================================================

system.time(
 volcano.pts <- map_xy(input.x   = 1:nrow(volcano),
                       input.y   = 1:ncol(volcano),
                       input.2D  = volcano,
                       output.xy = tomap)
)

head(volcano.pts)
attributes(volcano.pts)$processing

with(volcano.pts, points2D(x      = x, 
                           y      = y, 
                           colvar = v, pch = 18,
                           main   = "volcano, mapped data"))

par(mfrow=pm)

## =============================================================
## from 2D input to xyv output - increase resolution
## =============================================================

# --------------------------------------------------------------
# test: map volcano on itself - should give same as original
# --------------------------------------------------------------

system.time(
 volcano.tst <- map_xy(input.x   = 1:nrow(volcano),
                       input.y   = 1:ncol(volcano),
                       input.2D  = volcano,
                       output.xy = expand.grid(1:nrow(volcano),
                                               1:ncol(volcano))
                               )
)
summary(volcano.tst[, 3] - as.vector(volcano))  # all 0 !

# --------------------------------------------------------------
# double the resolution
# --------------------------------------------------------------

xdbl <- seq(from=1, to=nrow(volcano), by=0.5)
ydbl <- seq(from=1, to=ncol(volcano), by=0.5)

# needs to be done in 3 steps:
xy   <- expand.grid(x = xdbl, 
                    y = ydbl)

system.time(
 volcano.dbl <- map_xy(input.x   = 1:nrow(volcano),
                       input.y   = 1:ncol(volcano),
                       input.2D  = volcano,
                       output.xy = xy)
)

# recast output to matrix
volcdbl <- matrix(nrow = length(xdbl), 
                  ncol = length(ydbl), 
                  data = volcano.dbl$v)
                  
image2D(x=xdbl, y=ydbl, z=volcdbl)

# same resolution as volcano
VV    <- volcdbl[seq(1, nrow(volcdbl), by=2), 
                 seq(1, ncol(volcdbl), by=2)]
                 
image2D(VV-volcano) # SHOULD BE WHITE NOISE or 0

## =============================================================
## from 2D input to 2D output - increase resolution
## =============================================================

# SAME but in one step

system.time(
 volcano.dbl2 <- map_xy(input.x  = 1:nrow(volcano),
                        input.y  = 1:ncol(volcano),
                        input.2D = volcano,
                        output.x = xdbl,
                        output.y = ydbl)
)

VV    <- volcano.dbl2$v[seq(1, nrow(volcano.dbl2$v), by = 2), 
                        seq(1, ncol(volcano.dbl2$v), by = 2)]
image2D(VV-volcano) # SHOULD BE 0

## =============================================================
## find corresponding depths for sediment data
## =============================================================

# ----------------------------------
# use only data from the Waddensea
# ----------------------------------

# bathymetry
batmasked        <- mask_bat(Marsdiep, shape=Shape)  
batmasked$depth  <- Marsdiep$depth * batmasked$mask 

# grain size
Sedmasked <- na.omit(mask_xy(Sediment[ ,c("longitude", "latitude")],
                             shape=Shape))

system.time(
  Sed.ext <- map_xy(input.x   = batmasked$longitude, 
                    input.y   = batmasked$latitude, 
                    input.2D  = batmasked$depth,
                    output.xy = Sedmasked)
)

Sed.ext <- merge(Sed.ext, Sediment)
with(Sed.ext, plot(v, mdGrain))


## =============================================================
## space-time interpolation
## =============================================================

time <- seq(from = as.Date("2001-01-01"), 
            to   = as.Date("2002-01-01"), by = 14)

depth <- seq(0, 1, length.out=10)

# an interesting function
twoD <- outer(time, depth, FUN = function(t, x)
              1*exp(-0.03*x*abs(as.numeric(t-mean(time)))))

dim(twoD)
image2D(x = time, y = depth, z=twoD)

output.x <- seq(0, 1, length.out = 50)
output.t <- seq(from = as.Date("2001-01-01"), 
                to   = as.Date("2001-12-31"), by = 1)
                
twodHR <- map_tx(input.x = depth, input.t = time, 
                 input.2D = twoD,
                 output.x = output.x,
                 output.t = output.t)

image2D(twodHR$v)                                 
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
