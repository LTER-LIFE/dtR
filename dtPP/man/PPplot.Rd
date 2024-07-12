\name{Photosynthesis plotting}
\alias{plot.integratedPP}
\alias{image2D.integratedPP}

\title{
  Plots results of the depth-integrated photosynthesis calculation
}

\description{
  \code{plot} and \code{image2D} applied on an integratedPP object plots timeseries output and vertical light profile

  \code{image2D} applied on an integratedPP object plots an image of the time,depth distribution of light and volumetric photosynthesis.
}

\usage{
\method{plot}{integratedPP}(x, ..., type="l", las=1, which = "ts")
\method{image2D}{integratedPP}(z,  las=1, ...)
    
}

\arguments{
 \item{x, z }{An object of class \code{integratedPP} whose contents need to be plotted.
 }
  \item{type }{The type of the plots - the default draws lines. 
  }
  \item{las }{The label alignment; default is to have labels always horizontal.
  } 
  \item{which }{what to plot, one of "ts" (timeseries), "profile" (vertical profile), or both.
  }
  \item{... }{Any arguments passed to the plotting function. 
  } 
 
  
}
\value{
Does not return anything
}

\examples{

## =============================================================
## Changing water height
## =============================================================

depth <- 5              # Mean height of water column

# run model for 10 days; use date format
time <- seq(from = as.POSIXct("2021-08-01 01:00:00"), 
            to   = as.POSIXct("2021-08-11 01:00:00"), 
            by   = 3600)  # output every hour


# dataset for light intensity

st235   <- subset(Wad_weather, 
                  subset= station==235)

It.data <- st235[, c("datetime", "radiation")]

# water height varies in time (~tide)

t.period    <- 12.4*3600   # tidal period in seconds
t.amplitude <- 1           # tidal amplitude in meters

wH      <- t.amplitude * sin(2*pi*as.double(time)/t.period)

Ht.data <- data.frame(t      = time, 
                      height = wH)

PP <- integratedPP(zmax    = depth,   
                   time    = time, 
                   It.data = It.data, 
                   kz      = 0.2,
                   Ht.data = Ht.data)

plot(PP)

image2D(PP, 
        colkey=list(length=0.5, width=0.5))
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
  \code{\link{fitPI}} for fitting functions.

  \code{\link{FRRFdata}}, \code{\link{LABSTAFdata}} for data sets.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
