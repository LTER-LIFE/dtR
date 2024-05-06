\name{Ems-Dollard PI parameters and biogeochemistry}
\alias{Ems_ps}
\alias{Ems_biogeo}
\title{
 Water-column PI parameters and biogeochemical measurements taken in the Ems-Dollard (2012-2013). WMR data
}

\usage{
data(Ems_ps)
data(Ems_biogeo)
}

\value{
Biogeochemical data obtained from the Ems-Dollard at quasi monthly intervals, from februari 2012 till december 2013.

\code{Ems_ps} is  a list with the following items: 
\itemize{
  \item \code{station}, A unique name for the sampling station
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{a_fit, b_fit, c_fit}, phototynthesis-irradiance parameters
  \item \code{a_std_err  b_std_err c_std_err}, phototynthesis-irradiance parameter standard error
  \item \code{Imax, Pmax, alfa}, derived phototynthesis-irradiance parameters (Eilers-Peeters model)
  \item \code{Monod}, phototynthesis-irradiance half-saturation parameter
}

\code{Ems_biogeo} is  a list with the following items: 
\itemize{
  \item \code{station}, A unique name for the sampling station
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{PO4}, Phosphate concentration (mmol/m3)
  \item \code{NO3}, Nitrate concentration (mmol/m3)
  \item \code{NO2}, Nitrite concentration (mmol/m3)
  \item \code{NH4}, Ammonium concentration (mmol/m3)
  \item \code{DON}, Dissolved organic nitrogen (mmol/m3)
  \item \code{DOP}, Dissolved organic phosphorus (mmol/m3)
  \item \code{Si}, Dissolved silicate concentration (mmol/m3)
  \item \code{DIC}, Dissolved inorganic carbon concentration (mmol/m3)
  \item \code{Temperature} Water temperature (dgC)
  \item \code{Salinity}, Water salinity (-)
  \item \code{O2}, water oxygen concentration
  \item \code{Chl}, Chlorophyll concentration (mg/m3), from CTD (Cyclops)
  \item \code{turbidity}, Total turbidity, NTU
  \item \code{Cyanophyta_pgC.l, Diatoms_pgC.l, Dinophyceae_pgC.l, Flagellates_Other_pgC.l}, phytoplankton groups, carbon concentrations
  \item \code{Summ_Phyto_pgC.l, Summ_ug.Chla.l} summed phytoplankton carbon and chlorophyll concentration
    }
}

\description{
Pelagic photosynthesis parameters, and biotic and abiotic conditions in the Ems-Dollard.
}

\details{

Primary production was measured with the 14C technique, from radioactive bicarbonate uptake during 2 hour incubations in a series of light intensities, ranging from 0 to 1400 uE/m2/s. 

These data were fitted with the Eilers-Peeters model. 

Two sets of photosynthetic parameters are provided:

\itemize{

  \item \code{a_fit}, \code{b_fit} and \code{c_fit} are Eilers-Peeters raw parameters, where the photosynthesis (PP) as a function of light intensity (I) is defined as:
  
  \deqn{PP = \frac{I}{a \cdot I^2 + b \cdot I + c}}{PP=I/(a I^2 + b I + c)}
  
  where PP is in \eqn{mgC~(\mu~g chla )^{-1} h^{-1}}{(mg C /(ug chla)/ h}), and I is in \eqn{\mu E~m^{-2} s^{-1}}{uEinst/m2/s}
}

From these parameters, the derived parameters are estimated:
 
\itemize{
 
 \item \code{Pmax}, in \eqn{mgC~(mg~Chl)^{-1}h^{-1}}{mgC/mg Chl/h}, the maximum photosynthesis; \eqn{Pmax=\frac{1}{c \sqrt{a \cdot c} + b}}{Pmax=1/(c sqrt(a c) + b)}
 \item \code{alfa}, in \eqn{mgC~(mg~Chl)^{-1}h^{-1}(\mu E~m^{-2} s^{-1})^{-1}}{mgC/mg chl/h/(uE/m2/s)}, the slope of the PI curve \eqn{alpha=\frac{1}{c}}{alpha=1/c}
 \item \code{Imax}, the optimal light intensity, \eqn{\mu E~m^{-2} s^{-1}}{uE/m2/s}, \eqn{Imax = \sqrt{c/a}}{Imax = sqrt(c/a)}
 }


Extinction coefficients were estimated by fitting depth profiles of light, using the Lambert-Beer equation.

CTD data are from continuous measurements with a PocketBox, from which the following data were extracted:

temperature, 
salinity (derived from conductivity measurements)

oxygen (AADI 3835 optode sensor)

Chromophoric (or coloured) Dissolved Organic Materials (CDOM ), (CYCLOPS‐7 Submersible Fluorimeter)

Yellow substances include humic and fulvic acid -like components. A Yellow Substance sensor is part of the AOA‐equipment (= “Algae Online Analyser")

turbidity, measured by a Seapoint sensor

Algal composition, measured with algae online analyser; distinguishes four algae types: green algae, cyanobacteria (blue‐greens), diatoms and Cryptophycaea. All contents are expressed as micro g chla l‐1; the sum is stored as total chla (μg chla/ 1). Conversion to carbon was done based on average biovolumes.

See Brinkman et al., 2015.

Chlorophyll was also estimated from water samples and sediment samples, using PAM fluorometry (water column) or spectrophotomety (sediment samples).
Chlorophyll-a concentrations were quantified from water samples filtered over a GF/F filter and followed by acetone extraction. 

Phytoplankton cell counts were performed by flow cytometry; cell counts were converted to carbon biomass assuming spherical shapes; a fixed chlorophyll:C ratio is used to covert to chlorophyll.
}

\examples{

# ---------------------------------------------------
# The position of the data sets
# ---------------------------------------------------

stations <- meta(Ems_ps_parameters)$stations

stations

plotBathymetry(Wad_depth, type = "image", col=jet2.col(100),
              main = "Position of Wadden PS data (Ems-Dollard) ", 
              xlim = c(6, 7.2), ylim=c(53.2, 53.8))

points(stations$longitude, stations$latitude, 
       pch=18, cex=3, col="white")

text(x      = stations$longitude, 
     y      = stations$latitude, 
     labels = stations$station,
     cex    = 0.75,
     adj    = 1)


# ---------------------------------------------------
# Plot the data
# ---------------------------------------------------

mf <- par(mfrow=c(2,2))

# A function that first casts the data in wide format and then
# uses matplot to plot them

plotData <- function(var, main=var, pos="bottomright", ...){

    VV <- dttowide(Ems_ps_parameters[, c("datetime", "station", var)], 
                   swap = "station")
    nstat <- ncol(VV)-1 
    
    matplot(VV$datetime, VV[,-1], col = 1:nstat,
        main = main, xlab = "time", 
        las = 1, type = "l", lty=1, ...)
    
    if (! is.na(pos)) 
      legend(x=pos, lty=1, col=1:nstat, cex=0.7,
             legend=colnames(VV)[-1])    
    
}

plotData("Temperature", ylab="dgC")
plotData("Salinity",    ylab="",     pos=NA)
plotData("O2",          ylab="mg/l", pos=NA)
plotData("Turbidity",   ylab="NTU",  pos=NA)

plotData("PO4_P",       ylab="mg/l", pos=NA)
plotData("NH4_N",       ylab="mg/l", pos=NA)
plotData("NO3_N",       ylab="mg/l", pos=NA)
plotData("NO2_N",       ylab="mg/l", pos="topleft")

plotData("Cyclops-Chl-a", ylab="ug/l",                   pos="topleft")
plotData("alfa",          ylab="mgC/mg Chl/h/(uE/m2/s)", pos=NA)
plotData("Pmax",          ylab="mgC/mg Chl/h",           pos=NA)
plotData("Imax",          ylab="uE/m2/s",                pos=NA)

par(mfrow=mf)   

# ---------------------------------------------------
# Estimate station-wide summaries
# ---------------------------------------------------

# The Eilers-Peeters PI parameters

aggregate(x   = Ems_ps_parameters[,c("alfa", "Pmax", "Imax")], 
          by  = list(station = Ems_ps_parameters$station), 
          FUN = mean, na.rm = TRUE)

# The nutrient concentrations
aggregate(x   = Ems_ps_parameters[,c("NH4_N", "NO2_N", "NO3_N", "PO4_P", "Si")], 
          by  = list(station = Ems_ps_parameters$station), 
          FUN = mean, na.rm = TRUE)

# The abiotic conditions, oxygen and chlorophyll
aggregate(x   = Ems_ps_parameters[,c("Temperature", "Salinity", "O2", "Cyclops-Chl-a")], 
          by  = list(station = Ems_ps_parameters$station), 
          FUN = mean, na.rm = TRUE)

# The algal groups
aggregate(x   = Ems_ps_parameters[,c("Cyanophyta_pgC.l", "Diatoms_pgC.l",          
                                     "Dinophyceae_pgC.l", 
                                     "Flagellates_Other_pgC.l")], 
          by  = list(station = Ems_ps_parameters$station), 
          FUN = mean, na.rm = TRUE)
}
\seealso{
  \code{\link{Wad_chlor_RS}} a monthly, spatially-explicit chlorophyll data set.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}

\references{
Brinkman, A.G., Riegman, R., Jacobs, P., Kuhn, S., Meijboom, A. 2015.
Ems-Dollard primary production research: Full data report
Report / IMARES Wageningen UR : C160/14, 629, IMARES Ecosystemen, 

The Water Framework Directive (WFD) requires EU member states to achieve good ecological and chemical status of all designated water bodies (rivers, lakes, transitional and coastal waters) by 2015. Therefore Rijkswaterstaat Waterdienst has initiated the project &#8216;Research mud dynamics Ems Estuary&#8217; . The aim of this project, carried out by Deltares and IMARES, is to (1) improve our knowledge on the mud dynamics in the Ems Estuary, (2) to identify the reasons for the increase in turbidity and (3) to quantify measures to improve the ecological status of the estuary.

IJmuiden, IMARES, https://edepot.wur.nl/350873, ID  - 489709
}
\keyword{utilities}
