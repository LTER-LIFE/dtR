\name{average_xy}
\alias{average_xy}

\title{
 Spatial averaging in 2 dimensions.
}

\description{
  \code{average_xy} averages spatial (xy) values to a spatial grid. 
  To be used when the input grid is denser than the output grid.
}
\usage{
average_xy (input_xyv,
            input_x, input_y, input_2D,
            output_x, output_y)
}

\arguments{
  \item{input_xyv }{Data.frame with  input data, 3 columns: latitude (or x-values), 
  longitude (or y-values), and variable value. 
  }
  
  \item{input_x, input_y }{vectors with \code{x} and \code{y} input values, representing rows and columns of the 2D input data. 
  }
  
  \item{input_2D }{2-dimensional data to interpolate from, a matrix of dimension \code{(length(x), length(y))} or \code{(length(t), length(x))}. 
  }
  \item{output_x, output_y  }{vectors with \code{x} and \code{y} values, determining the output grid to which the input values will be averaged. 
  }
}

\value{
  Function \code{average_xy} returns a \code{list} with 3 items: 
  
  \itemize{
     \item \code{x} and \code{y} are the two vectors defining the x- and y- coordinates as input in \code{output_x} and \code{output_y},
     \item \code{v} will contain the 2-D matrix with averaged values.
  }

}

\examples{

## =============================================================
## from xyv input to 2D output - reduce resolution
## =============================================================

# map volcano onto itself - should return the same
# --------------------------------------------------------------
xy          <- expand.grid(x = 1:nrow(volcano), 
                           y = 1:ncol(volcano))

volcano.xyv <- data.frame(xy, 
                          v = as.vector(volcano))

output_x <- seq(1, nrow(volcano), by = 2)
output_y <- seq(1, ncol(volcano), by = 2)

system.time(
 volcano.grd <- average_xy(
                         input_xyv  = volcano.xyv ,
                         output_x = output_x, 
                         output_y = output_y)
)

image2D(volcano.grd$v)

## =============================================================
## from 2D input to 2D output - increase resolution
## =============================================================

# map volcano onto itself - should return the same
# --------------------------------------------------------------

 volcano2 <- average_xy(
                         input_x  = 1:nrow(volcano) ,
                         input_y  = 1:ncol(volcano) ,
                         input_2D = volcano,
                         output_x = 1:nrow(volcano), 
                         output_y = 1:ncol(volcano))

summary(as.vector(volcano-volcano2$v)) # all 0

# reduce grid resolution

system.time(
 volcano.grd2 <- average_xy(
                         input_x  = 1:nrow(volcano) ,
                         input_y  = 1:ncol(volcano) ,
                         input_2D = volcano,
                         output_x = output_x, 
                         output_y = output_y)
)

## =============================================================
## Reduce bathymetric resolution
## =============================================================

pm <- par(mfrow = c(2,2))
plotBathymetry(Marsdiep)

# mapping on itself
M1 <- average_xy(input_x  = Marsdiep$longitude ,
                 input_y  = Marsdiep$latitude ,
                 input_2D = Marsdiep$depth,
                 output_x = Marsdiep$longitude, 
                 output_y = Marsdiep$latitude)
image2D(x = M1$x, y = M1$y, 
        z = M1$v, asp = aspect_coord(M1$y))      
        
summary(as.vector(M1$v -Marsdiep$depth))       # 0 or NA 

# reducing the grid
output_x <- seq(from = 4.6, to = 5.1, by =  0.01)
output_y <- seq(from = 52.85, to = 53.1, by =  0.01)
M2 <- average_xy(
                 input_x  = Marsdiep$longitude ,
                 input_y  = Marsdiep$latitude ,
                 input_2D = Marsdiep$depth,
                 output_x = output_x, 
                 output_y = output_y)
image2D(x = M2$x, y = M2$y, 
        z = M2$v, asp = aspect_coord(M2$y))                         

# reducing the grid; input as xyv
input_xyv <- cbind(expand.grid(x = Marsdiep$longitude ,
                               y = Marsdiep$latitude),
                   v = as.vector(Marsdiep$depth))

M3 <- average_xy(
                 input_xyv = input_xyv,
                 output_x = output_x, 
                 output_y = output_y)
image2D(x = M3$x, y = M3$y, 
        z = M3$v, asp = aspect_coord(M3$y))                         

par(mfrow = pm)
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
