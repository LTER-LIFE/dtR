\name{NIOZ Jetty data}
\alias{Wad_biogeo_jetty}
\title{
 Water-column biogeochemical measurements taken from the NIOZ jetty (lon=4.789dgE lat=53.002dgN).
}

\usage{
data(Wad_biogeo_jetty)
}

\value{
Biogeochemical data obtained from the NIOZ jetty, at quasi-weekly intervals.

  a list with the following items: 
\itemize{
  \item \code{sampelID}, A unique identifier
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{Temperature}, Water temperature (dgC)
  \item \code{Salinity}, Water salinity (-)
  \item \code{Secchidepth}, Max depth at which the secchi disk is still visible (m)
  \item \code{Ntot}, Total nitrogen concentration (mmol/m3)
  \item \code{Ptot}, Total phosphorus concentration (mmol/m3)
  \item \code{PO4}, Phosphate concentration (mmol/m3)
  \item \code{NO3}, Nitrate concentration (mmol/m3)
  \item \code{NO2}, Nitrite concentration (mmol/m3)
  \item \code{NH4}, Ammonium concentration (mmol/m3)
  \item \code{DON}, Dissolved organic nitrogen (mmol/m3)
  \item \code{DOP}, Dissolved organic phosphorus (mmol/m3)
  \item \code{Si}, Dissolved silicate concentration (mmol/m3)
  \item \code{DIC}, Dissolved inorganic carbon concentration (mmol/m3)
  \item \code{SPM}, Total suspended particulate matter (mg/m3)
  \item \code{Chl}, Chlorophyll concentration (mg/m3)
    }
}

\description{
 Dataset obtained from the NIOZ monitoring platform. 
 
The NIOZ monitoring platform (jetty, 53 00'06''  N, 4 47'21'' E) ) extends 45 m into the Marsdiep, the westernmost tidal inlet of the Wadden Sea. 

The pier is equipped with automated sensors, which measure a suite of abiotic and biotic parameters such as water temperature, salinity, water level, hyper-spectral light conditions, chlorophyll-a and primary productivity. 

From 1974 onwards, sampling frequency was 40-60 times per year, varying from once or twice a month in winter up to
twice a week during phytoplankton spring blooms.
Sampling was always performed at high tide (+/- 10 min), i.e. when the water was still coming in from the North Sea. 

In the current dataset, only the biotic parameters are selected.

}

\details{

The sampling and measurement details:

\itemize{
\item Water temperature was measured using an oceanographic bucket thermometer. 
\item Salinity was measured by means of Autolab salinometer (Type 601MKIII)
until the mid-1980s and an ENDECO handheld salinometer (Type 2) hereafter.
\item Suspended particulate matter concentrations  were determined by 
filtration over pre-weighted filters. Filters were rinsed with 25 ml
demineralized water and dried overnight at 60 dgC. Several types of filters
were used; cellulose acetate filters with a pore
size of 0.45 microm (1974-1976; 1978-1985) and polycarbonate filters
with a pore size of 0.80 microm (2006-2010) were the most commonly
used ones.
\item Chlorophyll-a concentrations were quantified from water samples filtered 
over MgCO3 coated filters and based on non-acidified spectrophotometric measurements. 
Since 2009, chlorophyll-a was measured using a HPLC. The HPLC-derived values were converted to spectrophotometric values using the empirical equation CHLaSPEC=CHLaHPLC/0.82. 
\item The observed Secchi-disk depth is the depth at which a lowered white disk just
disappears from the observer's sight. This procedure was repeated two to three times, with an estimated accuracy of 5 cm. When the Secchi-disk depth was larger than the actual depth at
the sampling station, measurements were made by holding the stick
at a 45dg angle or performed at an alternative location, i.e. the NIOZ
harbor instead of the NIOZ jetty.
}
}

\examples{

# ---------------------------------------------------
# Plot the position of the data
# ---------------------------------------------------

plotBathymetry(Marsdiep, type = "image",
              pts = c(x=4.789, y=53.002), 
              ptlist = list(pch=18, cex=3),
              main = "Position of NIOZ jetty")

# ---------------------------------------------------
# Plot the high-resolution data
# ---------------------------------------------------

mf <- par(mfrow=c(2,2))

with(Wad_biogeo_jetty, 
   plot(x =datetime, y=Temperature, 
        main="temperature", ylab="dgC", xlab="time",
        las=1, type= "l"))

with(Wad_biogeo_jetty, 
   plot(x = datetime, y = Salinity, 
        main = "salinity", ylab = "", xlab = "time",
        las = 1, type= "l"))
        
with(Wad_biogeo_jetty, 
   plot(x = datetime, y = Secchi, 
        main = "secchi depth", ylab = "m", xlab = "time",
        las = 1, type = "l"))
        
with(Wad_biogeo_jetty, 
   plot(x = datetime, y = SPM, 
        main = "SPM", ylab = "mg/m3", xlab = "time",
        las = 1, type =  "l"))

par(mfrow=c(2,2))        
xlim <- as.POSIXct(c("1990-01-01", "2021-01-01"))

with(Wad_biogeo_jetty, {
   plot(x = datetime, y = PO4, 
        main = "PO4", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type= "l")
   plot(x = datetime, y = NH4, 
        main = "NH4", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type= "l")
   plot(x = datetime, y = NO3, 
        main = "NO3", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type= "l")
   plot(x = datetime, y = NO2, 
        main = "NO2", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type= "l")
})

xlim <- as.POSIXct(c("2000-01-01", "2021-01-01"))
with(Wad_biogeo_jetty, {
   plot(x = datetime, y = DON, 
        main = "DON", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type = "l")
   plot(x = datetime, y = DOP, 
        main = "DOP", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type = "l")
   plot(x = datetime, y = Ntot, 
        main = "Ntot", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type = "l")
   plot(x = datetime, y = Ptot, 
        main = "Ptot", ylab = "mmol/m3", xlab = "time",
        xlim = xlim, las = 1, type = "l")
})

par(mfrow=c(1,1))
with(Wad_biogeo_jetty, 
   plot(x = datetime, y = Chl, 
        main = "chlorophyll", ylab = "mg/m3", xlab = "time",
        las = 1, type = "l"))

pairs(Wad_biogeo_jetty[,8:16], pch = ".", cex = 3)

par(mfrow=mf)   
}
\seealso{
  \code{\link{Wad_chlor_RS}} a monthly, spatially-explicit chlorophyll data set.
}
\author{
  Data provider: Katja Philippart (NIOZ)
  
  Data extraction and preparation: Karline Soetaert <karline.soetaert@nioz.nl>
}

\references{
Cadee GC and Hegeman J 2002. Phytoplankton in the Marsdiep at the end of the 20th
century; 30 years monitoring biomass, primary production, and Phaeocystis
blooms. Journal of Sea Research 48, 97-110.

Philippart CJM, Suhyb Salama Mhd, Kromkamp JC, van der Woerd, HJ,
Zuur, AF and Cadee GC, 2013. Four decades of variability in turbidity in the western Wadden Sea
as derived from corrected Secchi disk readings. Journal of Sea Research 82, 67-79.
}
\keyword{utilities}
