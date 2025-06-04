\name{interpolate_xy}
\alias{interpolate_xy}

\title{
 Spatial interpolation in 2 dimensions.
}

\description{
  \code{interpolate_xy} interpolates spatial (xy) values to other spatial points. 
  Uses inverse distance weighing. 
}
\usage{
interpolate_xy (input_xyv,
                input_x, input_y, input_2D,
                output_xy, 
                output_x, output_y,
                nmean = 3, asp = c("geographic", "mean", "none"))
}

\arguments{
  \item{input_xyv }{Data.frame with  input data, 3 columns: latitude (or x-values), 
  longitude (or y-values), and variable value. 
  }
  
  \item{input_x, input_y }{vectors with \code{x} and \code{y} input values, representing rows and columns of the 2D input data. 
  }
  
  \item{input_2D }{2-dimensional data to interpolate from, a matrix of dimension \code{(length(x), length(y))} or \code{(length(t), length(x))}. 
  }
 \item{output_xy }{Two-columned matrix with the latitude (or x-values), and 
  longitude (or y-values) of the output data.
  }
  \item{output_x, output_y  }{vectors with \code{x} and \code{y} values, determining the x and y values to which the input data will be mapped. 
  Ignored if \code{output_xy} is given a value.
  }
 \item{asp }{The y/x aspect ratio for weighing y versus x. If "geographic" then it is assumed that the spatial points are geographic (spherical) coordinates (latitude, longitude) 
 and an aspect ratio for distances will be used (based on \link{aspect_coord}). 
 If "mean" then the distances will be weighed according to the mean values of x and y; the aspect ratio will be mean(x)/mean(y). If "none", there will be no rescaling of y- and x- distances (and asp will be set = 1). Also allowed is to give asp a number.
  }
  \item{nmean }{The number of points over which the new value will be averaged. 
  }
}

\value{
  Function \code{interpolate_xy} returns either a \code{data.frame} or a \code{list}, depending on the output arguments.
  
  When the output data are specified as \code{output_xy}, then the output will be a \code{data.frame} with 3 columns: 
  
  \itemize{
     \item first and second columns contain the x- and y- coordinates as defined in \code{output_xy},
     \item third column contains the interpolated values.
  }

  When the output data are specified as two vectors, \code{output_x} and  \code{output_y}, then the output will be a \code{list} with 3 items: 
  
  \itemize{
     \item \code{x} and \code{y} are the two vectors defining the x- and y- coordinates as input in \code{output_x} and \code{output_y},
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

input_xyv <- rbind(d1, d2, d3, d4, d5)

n <- 4
output_xy <- cbind(x = runif(n = n, min = 1, max = 2), 
                   y = runif(n = n, min = 2, max = 3))

output_v <- interpolate_xy(
                      input_xyv = input_xyv, 
                      output_xy = output_xy)
input_xyv
output_xy
output_v

interpolate_xy(input_xyv = input_xyv, 
               output_xy = input_xyv)
           
## =============================================================
## interpolation from xyv to xyv - volcano
## =============================================================

pm <- par(mfrow = c(1, 2))

# Volcano is a matrix - cast it in x,y,v format

xy          <- expand.grid(x = 1:nrow(volcano), 
                           y = 1:ncol(volcano))

volcano.xyz <- data.frame(xy, 
                          v = as.vector(volcano))

with(volcano.xyz, 
     points2D(x, y, colvar = v, pch = 18, 
              main = "volcano, original data"))

head(volcano.xyz)

# map volcano onto itself - should return the same
# --------------------------------------------------------------

system.time(
 volcano.tst <- interpolate_xy(
                       input_xyv   = volcano.xyz, 
                       output_xy   = volcano.xyz[,1:2])
)

summary(volcano.tst[,3]-volcano.xyz[,3])  # all 0


# select 2000 random points to map
# --------------------------------------------------------------

n <- 2000
tomap <- data.frame(x = runif(min = 1, 
                              max = nrow(volcano), n = n),
                    y = runif(min = 1, 
                              max = ncol(volcano), n = n))

# interpolate
system.time(
 volcano.pts <- interpolate_xy(
                       input_xyv = volcano.xyz, 
                       output_xy = tomap)
)

head(volcano.pts)
attributes(volcano.pts)$processing

with(volcano.pts, points2D(x      = x, 
                           y      = y, 
                           colvar = v, pch = 18,
                           main   = "volcano, mapped data"))

par(mfrow=pm)

## =============================================================
## from 2D to xyv interpolation - volcano
## =============================================================

# map volcano onto itself - should return the same
# --------------------------------------------------------------

system.time(
 volcano.tst2 <- interpolate_xy(
                       input_x   = 1:nrow(volcano),
                       input_y   = 1:ncol(volcano),
                       input_2D  = volcano,
                       output_xy = volcano.xyz[,1:2])
)

summary(volcano.tst2[,3]-volcano.xyz[,3])  # all 0

# select 2000 random points to map
# --------------------------------------------------------------

n <- 2000
tomap <- data.frame(x = runif(min = 1, 
                              max = nrow(volcano), n = n),
                    y = runif(min = 1, 
                              max = ncol(volcano), n = n))
system.time(
 volcano.pts2 <- interpolate_xy(
                          input_x   = 1:nrow(volcano),
                          input_y   = 1:ncol(volcano),
                          input_2D  = volcano,
                          output_xy = tomap)
)

head(volcano.pts2)
attributes(volcano.pts2)$processing

with(volcano.pts2, points2D(x     = x, 
                            y      = y, 
                            colvar = v, pch = 18,
                            main   = "volcano, mapped data"))

par(mfrow=pm)

# double the resolution
# --------------------------------------------------------------

xdbl <- seq(from = 1, to = nrow(volcano), by = 0.5)
ydbl <- seq(from = 1, to = ncol(volcano), by = 0.5)

# interpolation in 3 steps:
xy_dbl   <- expand.grid(x = xdbl, 
                        y = ydbl)

system.time(
 volcano.dbl <- interpolate_xy(
                          input_x   = 1:nrow(volcano),
                          input_y   = 1:ncol(volcano),
                          input_2D  = volcano,
                          output_xy = xy_dbl)
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
## from xyv input to 2D output - increase resolution
## =============================================================

# map volcano onto itself - should return the same
# --------------------------------------------------------------
xy          <- expand.grid(x = 1:nrow(volcano), 
                           y = 1:ncol(volcano))

volcano.xyv <- data.frame(xy, 
                          v = as.vector(volcano))

output_x <- 1:nrow(volcano)
output_y <- 1:ncol(volcano)

system.time(
 volcano.grd <- interpolate_xy(
                         input_xyv  = volcano.xyv ,
                         output_x = output_x, 
                         output_y = output_y)
)

summary(as.vector(volcano.grd$v-volcano))  # all 0

## =============================================================
## from 2D input to 2D output - increase resolution
## =============================================================

# map volcano onto itself - should return the same
# --------------------------------------------------------------

# SAME but in one step

system.time(
 volcano.dbl2 <- interpolate_xy(
                      input_x  = 1:nrow(volcano),
                      input_y  = 1:ncol(volcano),
                      input_2D = volcano,
                      output_x = xdbl,
                      output_y = ydbl)
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
batmasked        <- mask_shape(longitude = Marsdiep$longitude, 
                               latitude  = Marsdiep$latitude, 
                               shape     = Shape)  
batmasked$depth  <- Marsdiep$depth * batmasked$mask 

# grain size
Sedmasked <- na.omit(
      mask_shape(coordinates = Sediment[ ,c("longitude", "latitude")],
                 shape       = Shape))

system.time(
  Sed.ext <- interpolate_xy(
                    input_x   = batmasked$longitude, 
                    input_y   = batmasked$latitude, 
                    input_2D  = batmasked$depth,
                    output_xy = Sedmasked)
)

Sed.ext <- merge(Sed.ext, Sediment)
with(Sed.ext, plot(v, mdGrain))

}

\seealso{
  \code{\link{Sediment}} data set with sediment composition.
  \code{\link{Wad_watertemp_LR}} data set with (low resolution) water temperature from the Waddensea.
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
