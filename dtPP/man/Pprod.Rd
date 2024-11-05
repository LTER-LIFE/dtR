\name{Photosynthesis}
\alias{Photosynthesis}
\alias{integratedPP}

\title{
  Estimates depth-integrated photosynthesis
}

\description{
  \code{integratedPP} uses a combination of P-I (photosynthesis-irradiance) parameters,
  a time series of light intensities on the surface of the water or sediment, and 
  light extinction in the water or sediment, to estimate depth-integrated photosynthesis rates.
  
  It can be used on one site (time series), on a spatial grid, and 
  for rivers, estuaries and lakes; in the latter case also the water body's morphology needs
  to be input.
}

\usage{
integratedPP (
    zmin = 0, zmax = 100, nbox = 100, zi = NULL, 
    times = 0:10, convFac = 1.,
    PI.par = c(alpha=0.05, eopt=1000, ps=50), 
    PI.fun = NULL, 
    It.data = data.frame(time=0:10, I=(0:10)*10), 
    It.fun = NULL,
    kz = 0.1, Iz.data = NULL, Iz.fun = NULL, 
    Mz.data = 1, Mz.fun = NULL,
    Ht.data = NULL, Ht.fun = NULL,
    avgOver = NULL, avgTime = 1,
    verbose = TRUE, 
    unit = list(mass = "mgC", length = "m", time = "h", light = "uEinst/m2/s"))   
}

\arguments{
 \item{zmin, zmax }{The minimal and maximum depth over which the photosynthesis needs to be summed, typical unit: [m]. If \code{zi} is not \code{NULL}, then \code{zmin} and \code{zmax} will be estimated as the minimum and maximum of \code{zi} respectively. In the other case, \code{zi} will be estimated based on  \code{zmin}, \code{zmax} and \code{nbox}.
 }
 \item{nbox }{Number of boxes for which the photosynthesis will be calculated and summed. If \code{zi} is not \code{NULL} \code{nbox} will be estimated as the length of \code{zi} - 1.
 }
 \item{zi }{Depth at the box \code{i}nterfaces. If not \code{NULL} then \code{zi} will be used to estimate nbox (length(zi)-1)), \code{zmin} (min(zi)) and \code{zmax} (max(zi))). If NULL, \code{zi} will be calculated using an exponentially increasing grid size, where the size of the first box equals \code{(zmax-zmin)/nbox/10}.
 }
 \item{times }{Time for which the integrated photosynthesis needs to be calculated. The unit of times should be consistent with the units of the time series. Can be of type \code{POSIXct}. If \code{avgOver} is NULL, then \code{times} will also be the output time.
 }
 \item{convFac }{Factor to convert integrated photosynthesis to other units.
 }
 \item{PI.par }{Photosynthesis-irradiance parameters, a \code{vector}, or a \code{data.frame} or \code{matrix}. If a \code{data.frame} or matrix, then there are two options: PI.par varying in depth or PI.par varying over time. In the first case, the name of the first column should be \code{depth} or \code{z}, and the values of PI.par will be interpolated over depth. In the latter case, the name of the first column should be \code{time} (or something else, this is not checked). The parameters will then be interpolated over all output times. Note that when \code{PI.par} is not a vector, then \code{PI.fun} should be accepting a data.frame of parameter values. By default \code{PI.par} contains parameters from the Eilers-Peeters model, i.e. containing the elements called \code{alpha}, \code{eopt}, \code{ps} (in that order).
 }
 \item{PI.fun }{Photosynthesis-irradiance function.  Will be called as PI.fun(I, par), where \code{I} is the light intensity vector or matrix, and \code{par} are the parameters (either a vector or a list with parameters having the same dimensions as I). Function \code{PI.fun} should return the photosynthesis for these light intensities, i.e. the return element of the function should have the same dimensions or length as \code{I}. Defaults to the Eilers-Peeters model (\link{fEP}).
 }
 \item{It.data, It.fun }{Light intensity; either  one value, or a dataseries of light as a function of time (t); typical units are uEinst/m2/s. If not one value, then data set \code{It.data} should be a two-column data.frame or matrix containing (time, It). \code{It.data} will be ignored if \code{It.fun} is present. In that case, function \code{It.fun}, will be called as \code{It.fun(times)}, where times is the time at which the integrated production needs to be estimated. If \code{It.fun} is \code{NULL} (the default), then it will be created based on \code{It.data} as: \code{approxfun(x=It.data, rule=2)}.
 }
 \item{kz }{Light extinction coefficient parameter(s); typical units of [m^{-1}]. 
  \code{kz} will only be used if \code{Iz.data} and \code{Iz.fun} are \code{NULL}. Typically \code{kz} is one value or three values, depending whether one class or two classes of light extinction is assumed. It is also allowed to define a data set of one or three kd values as a function of time. In the latter case, the input will be a matrix or data.frame where the first column is time, and the other columns define the extinction parameters. Based on \code{kz}, an exponential decay of light is estimated as \code{exp(-kz*z)} or \code{kz[1]*exp(-kz[2]*z)+(1-kz[1])*exp(-kz[3]*z)}. See details.
 }
 \item{Iz.data, Iz.fun }{Data set or function that defines light extinction as a function of (water or sediment) depth (z), a dimensionless quantity (fraction) that is typically 1 (or less) near the interface (depth=0) and declines exponentially with depth (z). Data set \code{Iz.data} should be a two-column data.frame or matrix containing (z, Iz). \code{Iz.data} will be ignored and if \code{Iz.fun} is present. In that case, \code{Iz.fun} will be called as \code{Iz.fun(z)}, where z is the depth at the center of the boxes. If \code{Iz.fun} is \code{NULL} (the default), then it will be created as: \code{approxfun(x=Iz.data, rule=2)}. 
 }
 \item{Mz.data, Mz.fun }{Morphology (surface area) as a function of (water or sediment) depth (z), typical units of [m2].  Data set \code{Mz.data} should be a two-column data.frame or matrix containing (z, Mz). \code{Mz.data} will be ignored and if \code{Mz.fun} is present. In that case, function \code{Mz.fun} will be called as \code{Mz.fun(z)}, where z is the depth at the center of the boxes. If \code{Mz.fun} is \code{NULL} (the default), then it will be created as: \code{approxfun(x=Mz.data, rule=2)}. 
 }
 \item{Ht.data, Ht.fun }{Time series or function of (water) elevation; typically in units of [m]. Based on this, the actual depth of the water column will be estimated as zn + Ht. Specification of water elevation makes sense only for shallow areas, where a significant fraction of the light extends till the bottom (at znmax). Note that this will only clip the depth to the calculated depth and set the production below this depth equal to 0. If not \code{NULL}, then data set \code{Ht.data} should be a two-column data.frame or matrix containing (time, Ht). \code{Ht.data} will be ignored and if \code{Ht.fun} is present. In that case, function \code{Ht.fun} will be called as \code{Ht.fun(times)}, where times is the time at which the integrated production needs to be estimated. If \code{Mt.fun} is \code{NULL} (the default), then it will be created as: \code{approxfun(x=Ht.data, rule=2)}. 
 }
 \item{avgOver }{Only when \code{times} is of type POSIXct or POSIXlt: the time instance over which has to be averaged, i.e. one of "sec", "min", "hour", "mon", "year", "mday" (day of the month), "wday" (day of the week), "yday" (day of the year). See \link[dtTools]{average_timeseries}.
 }
 \item{avgTime }{The length of time units over which need to be averaged.
 }
 \item{verbose }{If TRUE, will write warnings to the screen.
 }
 \item{unit }{A vector or list with units of \code{mass}, \code{length}, \code{time}, \code{light}. 
 Names of mass (e.g. ug C), length (e.g. m), time (second), and light intensity. This will be used to generate units of the output variabes, and to label the plots. (note: units generated do not take into account convfac )
 }
}

\details{
Input of time-varying or depth-varying parameters can be either as a data set, or as a function, where the latter has precedence over the former. 

Note that in all cases, the interpolation will use rule = 2, i.e. it will extrapolate outside of the data, where the value at the closest data extreme will be used.

The time-varying parameters are:

 \itemize{
  \item light intensity, defined as \code{It.data} or \code{It.fun}
  \item surface elevation, defined as \code{Ht.data} or \code{Ht.fun}
  \item PI parameters, when \code{PI.par} is a matrix or data.frame with the first column is "time".
  \item kz can also be variable in time, if it is a matrix, with the first column the "times".
 }

The depth-varying parameters are:

 \itemize{
  \item light intensity, defined as \code{Iz.data} or \code{Iz.fun}
  \item morphology, defined as \code{Mt.data} or \code{Mt.fun}
  \item PI parameters, when \code{PI.par} is a matrix or data.frame with the first column is "depth".
 }
}
  
\value{
a \code{list} with several items:
 \itemize{
  \item \code{ts}; the time series: a \code{data.frame} with several columns : 
   \itemize{
    \item times, the output times (argument times), if \code{avgOver} is not NULL, this will be averaged over time. 
    \item PP, the depth-integrated photosynthesis for each time, 
    \item TotalDepth, the sum of \code{zmax} and \code{H.t} (mean depth + surface elevation), 
          i.e. the total depth over which PP was estimated,
    \item Isurf, the light intensity at the surface,
    \item Ibot, the remaining light at the bottom of the grid,
    \item kz, the extinction coefficient (or the mean extinction coefficient if 3 values are input).
    \item If the photosynthetic parameters change over time, their value at each time point.
   }
   
  \item profile, a data.frame with the following columns:
    \itemize{
     \item z, the depth at which photosynthesis has been estimated,
     \item Iz_I0, the relative light profile, Iz/I0; in case light extinction changes over time this will be the mean profile over the time series,
     \item surface, the surface area at each depth (morphology); will be =1 if \code{Mz.data} and \code{Mz.fun} were NULL; 
     \item If the photosynthetic parameters change over depth, their value at each z point
    }
  
  \item light, the 2D matrix of light intensities, of dimension (length(times), length(z)); in case times is one number: a matrix with one row. In case \code{avgOver} is not NULL, this will be averaged over time.
  
  \item prod, the 2D matrix of volumetric production rates, of dimension (length(times), length(z)); in case times is one number: a matrix of one row. In case \code{avgOver} is not NULL, this will be averaged over time.
 }
}


\examples{

## =============================================================================
## Simple photosynthesis calculation
## =============================================================================

# Uses the default settings
PP <- integratedPP()
names(PP)

meta(PP)  # show metadata

PP$ts

attributes(PP)$Description
plot(PP)

## =============================================================================
## changing light intensity, imposing a light depth profile
## =============================================================================

# --------------------------------------
# light profile: two extinction classes
# --------------------------------------
z       <- seq(from=0, to=20, by=0.1)      # depth in m
Iz.data <- data.frame(z  = z, 
                      Iz = 0.5*exp(-0.05*z) + 0.5*exp(-0.15*z))
plot(Iz.data)

# --------------------------------------
# light intensity versus time
# --------------------------------------
st235   <- subset(Wad_weather, 
                  subset= station == 235            &
                          datetime >= "2021-07-31"  &
                          datetime <= "2021-08-06"   )
                    
                    
It.data <- st235[, c("datetime", "radiation")]
plot(It.data)

# --------------------------------------
# height of water is constant, 5 m
# --------------------------------------

depth <- 5

# --------------------------------------
# run model for 4 days; use date format
# --------------------------------------

time <- seq(from = as.POSIXct("2021-08-01 01:00:00"), 
            to   = as.POSIXct("2021-08-05 01:00:00"), 
            by   = 60)  # output every minute


# Photosynthesis with changing irradiance
P1 <- integratedPP(zmax    = depth,      # water depth
                   times   = time,  
                   It.data = It.data,    # irradiance 
                   Iz.data = Iz.data)

# Plot results
par(mfrow=c(2,2))

plot(P1)   # plots P1$ts

plot(P1, which="profile")

with(P1$profile, 
  plot(Iz_I0, z, 
       xlab = "-", ylab = "water depth, m", 
       main = "Light penetration", 
       ylim = rev(range(z)), 
       type = "l")
  )

## =============================================================================
## Changing water height
## =============================================================================

# water height varies in time (~tide)
depth <- 5              # Mean height of water

# run model for 4 days; use date format
time <- seq(from = as.POSIXct("2021-08-01 01:00:00"), 
            to   = as.POSIXct("2021-08-05 01:00:00"), 
            by   = 60)  # output every minute
            
t.period    <- 12.4*3600   # tidal period in seconds
t.amplitude <- 1           # tidal amplitude in meters

wH          <- t.amplitude * sin(2*pi*as.double(time)/t.period)
Ht.data <- data.frame(t      = time, 
                      height = wH)

PP <- integratedPP(zmax    = depth,   
                   times   = time, 
                   It.data = It.data, 
                   Iz.data = Iz.data,
                   Ht.data = Ht.data)

# plot of light, depth, integrated production
plot(PP)

# image2D(PP$light, ylim = c(6,0))

## =============================================================================
## Compare both runs
## =============================================================================

par(mfrow=c(2,2))

# integrated production
plot(PP$ts$times, PP$ts$PP,
       xlab = "times", ylab = "umolC/m2/s", 
       main = "Integrated production", 
       type = "l")
       
lines(P1$ts$times, P1$ts$PP, col="red")

# time-depth  variables

times <- PP$ts$times
depth <- PP$profile$z

image2D(x = times, y = depth, z = PP$light,
        main = "Light - variable water height", 
        ylim = c(6,0))

image2D(x = times, y = depth, z = P1$light,
        main = "Light - constant water height", 
        ylim = c(6,0))

## =============================================================================
## A vertical profile of PI-parameters
## =============================================================================

PP <- integratedPP(
           times    = 0, 
           It.data  = 100,
           PI.par   = data.frame(
                       depth = c(0,   50,  100), 
                       alpha = c(0.5, 0.5, 0.1), 
                       eopt = 180, 
                       ps    = c( 1,   10,   5))
                  )

PP$ts

par(mfrow=c(2,2))
with(PP$profile, {
  
  plot(Iz_I0, z,
     xlab = "-", ylab = "water depth, m", 
     main = "Light penetration", 
     ylim = rev(range(z)), type = "l")
       
  plot(alpha, z,
     xlab = "umolC/m/uEinst", ylab = "water depth, m", 
     main = "alpha", 
     ylim = rev(range(z)), type = "l")

  plot(eopt, z,
     xlab = "uEinst/m2/s", ylab = "water depth, m", 
     main = "eopt, optimal light intensity", 
     ylim = rev(range(z)), type = "l")

  plot(ps, z,
     xlab = "umolC/m3/s", ylab = "water depth, m", 
     main = "ps, maximal photosynthesis rate", 
     ylim = rev(range(z)), type = "l")
})

# only the profiles
plot(PP, which="profile")

## =============================================================================
## A time series of PI-parameters
## =============================================================================

times <- 0:200

PP <- integratedPP(
    times   = times, 
    
    It.data = data.frame(time=times, 
                         I   =pmax(0, sin(2*pi*times/24))*100),
    
    PI.par  = data.frame(
        time = c(0, 50, 200),    # change over time
        alpha = c(0.5, 0.5, 0.1), 
        eopt = 180, 
        ps=c( 1, 10, 5)))
        
head(PP$ts)

par(mfrow=c(2,2))
with(PP$ts, {
  plot(times, Isurf, 
       xlab = "times", ylab = "uEinst/m2/s", 
       main = "Irradiance", type = "l")
  
  plot(times, PP, 
       xlab = "times", ylab = "umolC/m2/s", 
       main = "Integrated production", type = "l")

  plot(times, alpha, 
       xlab = "times", ylab = "/(uEinst/m2/s", 
       main = "alpha", type = "l")
 
  plot(times, ps, 
       xlab = "times", ylab = "ugC/m2/s", 
       main = "ps", type = "l")

})

# similar, but all timeseries variables
plot(PP)

## =============================================================================
## Two extinction classes
## =============================================================================

integratedPP(time    = 0, 
             It.data = 10,  # ct light intensity
             kz      = 0.05)$ts

integratedPP(time    = 0, 
             It.data = 10,
             kz      = c(0.9, 0.05, 0.4))$ts

integratedPP(time    = 0, 
             It.data = 100,
             kz      = c(0.9, 0.05, 0.4))$ts

## =============================================================================
## Extinction parameters vary over time
## =============================================================================

PP2 <- integratedPP(time    = 0:20, 
                    It.data = 10,
                    kz      = data.frame(
                               time = c(0,    20), 
                               kz   = c(0.4, 0.9)))

PP2 <- integratedPP(time    = 0:20, 
                    It.data = 10,
                    kz      = data.frame(
                               time = c(0,    20), 
                               p    = c(0.9, 0.95), 
                               k1   = c(0.05, 0.1), 
                               k2   = c(0.4, 0.9)))

with(PP$profile, 
  plot(Iz_I0, z,
       xlab = "-", ylab = "water depth, m", 
       main = "Average Light penetration", 
       ylim = rev(range(z)), type = "l")
     )
with(PP2$profile, 
      lines(Iz_I0, z, 
            col = "red")
     )

## =============================================================================
## The NIOZ Jetty
## =============================================================================

# --------------------------------
# Chlorophyll
# --------------------------------
# Chlorophyll conc for 2021 do not exist, so we use the 2020 data instead

ChlJetty <- subset(Wad_biogeo_jetty, 
                   subset = datetime >= as.POSIXct("2019-12-15") &
                            datetime <= as.POSIXct("2021-01-01"))

ChlJetty <- ChlJetty [, c("datetime", "Temperature", "Chl", "Secchi")]

# pretend these data are from 2021 (add one year to the datetime)

ChlJetty$datetime <- ChlJetty$datetime+86400*365
                                  
# --------------------------------
# Irradiance
# --------------------------------
# use the weatherdata from station 235 for irradiance
Irr2021 <- subset(Wad_weather, 
                  subset = station==235)[,c("datetime", "radiation")]

# --------------------------------
# extinction coefficient 
# --------------------------------
# kz varies with time; assume relationship with Secchi:
ChlJetty$kz <- 7/ChlJetty$Secchi  # in /m

# --------------------------------
# PI parameters
# --------------------------------
# values for alpha, eopt and ps - 
3
# NOTE: THESE ARE NOT VALUES FOR THE WADDENSEA but for the Ems-Dollard:
# Brinkman & Jacobs, 2023. Gross pelagic primary production in the Ems-Dollard estuary,
# Journal of Sea Research 192 (2023) 102362
#
# alpha and ps are Chl-specific values -> they need to be multiplied with Chl
# at each time step
# eopt and ps depend on temperature


# steepness of the light curve, mgC/mgChl/(uEinst/m2/s)
alpha <- 0.049  

# optimal light intensity (in uEinst/m2/s - original values in W/m2)
eopt  <- 4*(150 + 15*ChlJetty$Temperature)

# assumed temperature dependence for maximal gross production
ps    <- 13*1.06^(ChlJetty$Temperature-20)    # mgC/mgChl/h

# time-variable PI parameters
PI.par <- data.frame(time  = ChlJetty$datetime,
                     alpha = alpha*ChlJetty$Chl,
                     eopt  = eopt,
                     ps    = ps*ChlJetty$Chl)
                     
# --------------------------------
# Integrated production
# --------------------------------

times <- seq(from = as.POSIXct("2021-01-01"), 
             to   = as.POSIXct("2021-12-31"), 
             by   = 3600)

# integrated production, averaged over a day

Pprod_day <- integratedPP(
                 zmax    = 10,            # water depth
                 times   = times,
                 PI.par  = PI.par,
                 It.data = Irr2021,
                 kz      = ChlJetty[, c("datetime", "kz")],
                 avgOver = "day", 
                 avgTime = 1
)                     

plot(Pprod_day)

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
