\name{standardizeFRRF}
\alias{standardizeFRRF}

\title{
 Standardize FRRF data and estimate photosynthesis.
}

\description{
 \code{standardizeFRRF} standardizes the FRRF data with the blanc, and calculates the volumetric e-flux \code{JVPII}, according to the absorption method. 
  
}

\usage{
standardizeFRRF(frrf, Fblanc=0, Ka=11800, convJVPII=3600/1000, aLHII_0=NA,
                na.rm=TRUE, verbose=FALSE)
}

\arguments{
  \item{frrf }{An object of type \code{frrf}, as read with FRRFread. Only works for data from the act2 apparatus.
  }
  \item{Fblanc }{The  background fluorescence of the water, i.e. the Fo value of a filtered water sample (without chlorophyll).
  }
  \item{Ka }{Instrument-specific constant, units of m^{-1}. See the manual.
  }
  \item{convJVPII }{Conversion factor for the calculation of the volumetric e-flux JVPII. The default 3.6, converts from umol/m3/s to mmol/m3/hour.
  }
  \item{aLHII_0 }{The cross-sectional surface of the PSII system in the dark (estimated at E=0).
  }
  \item{na.rm }{when \code{TRUE}, will remove the failed data points from the file.
  }
  \item{verbose }{when \code{TRUE}, will print warnings to the screen.
  }
  
}

\note{
Data generated with the  *act2 FRRF* are unstandardized.
Before using these data, it is advised to standardize the data first with the correct blanc, 
and to reestimate the photosynthesis rate JVPII, which may then also be cast in correct units.

If for the LabStaf the blancs are correctly inputted, then the standardisation is already done ant the 
JVPII is already well estimated Uusing the absorption method. 
However, if the blan inputted was = 0, standardization can be done afterwards. 
Unless \code{verbose=FALSE}, a warning will be issued that standardization is attempted on standardized data.
}
\value{
  Returns an updated frrf data.frame where the columns Fo and Fm are replaced with:

\itemize{
 \item Fo=Fo-Fblanc
 \item Fm=Fm-Fblanc
}

The new data.frame is extended with columns labeled:
\itemize{
 \item Fo_uc Fm_uc JVPII_uc, the uncorrected (original) values for Fo, Fm, JVPII. 
 \item Fq = Fm - Fo
 \item a_LHII =Fo*Fm/Fq*(Ka/10^6)), the cross-sectional surface of the PSII system [nm^2]
 \item FqFm = Fq/Fm
 \item JVPII = Fq/ Fm* a_LHII_0 * E * convJVPII, the volumetric photosynthesis
}

Here a_LHII_0 is the cross-sectional surface of the PSII system in the dark. 
If a_LHII is not known at E=0 then it is derived by regressing a_LHII versus E for E < 100.

The value of a_LHII_0 is stored in the attributes of the data.frame.

The objects processing attribute will also be amended with the standardization info 

The attributes of the object will contain the attributes of the original object, amended with a check data.frame that records the "suspicious" value, according to their 
}

\details{
The method for estimating JVPII is the absorption method.

If not specified, then aLHII_0 is estimated at E=0; in case this produces NA, then a linear regression of aLHII0 versus E for E < 100 is performed.
}


\examples{

## =============================================================
## Standardize and estimate volumetric photosynthesis
## act2 data
## =============================================================

# path of the external data of the package
fpath <- system.file("extdata", package="dtPP")

# Read the FRRF file
F1 <- readFRRF(file = "49O_40m_FRRF.csv", 
               dir  = fpath)

# Standardize the FRRF file
F1s  <- standardizeFRRF(F1, 
                        Fblanc = 0.1)

F1s2 <- standardizeFRRF(F1, 
                        Fblanc = 0.15)

head(F1s)

attributes(F1s)$check
attributes(F1s)$processing

with(F1s,  
     plot(E, JVPII, xlim=c(0,200), ylim=c(0,1)))
with(F1s2, 
     lines(E, JVPII, col="red"))

## =============================================================
## Standardize and reestimate volumetric photosynthesis
## Labstaf data
## =============================================================

# read the LabSTAF data:
fpath <- system.file("extdata", package="dtPP")
F2    <- readFRRF(file = "49O_1m_LS.txt", 
                  dir  = fpath, 
                  txt  = "delim")

# show settings of the LabSTAF:
attributes(F2)$settings

# recalculate the JVPII, use same blanc (=0), and alHII estimated by labstaf
F2sb  <- standardizeFRRF(F2, 
                         Fblanc    = 0.0, 
                         convJVPII = 1, 
                         aLHII_0   = 0.02227)

# recalculate the JVPII, use same blanc, but different aLHII_0
F2s   <- standardizeFRRF(F2, 
                         Fblanc    = 0.0, 
                         convJVPII = 1)

# recalculate the JVPII, use different blanc
F2sc  <- standardizeFRRF(F2, 
                         Fblanc    = 0.1, 
                         convJVPII = 1)

with(F2,  
     plot(E, JVPII, type="l", ylim=c(0, 3), las=1, lwd=2))
with(F2s,  
     lines(E, JVPII, col="red", type="b"))
with(F2sb, 
     lines(E, JVPII, col="blue", type="b"))
with(F2sc, 
     lines(E, JVPII, col="green", type="b"))

legend("bottomright", 
       col=c("black", "red", "blue", "green"),
       pch=c(NA, 1,1,1), lwd=c(2,1,1,1), 
       legend=c("labstaf original", "reestimated",
                "reestimated with Labstaf parameters", 
                "different blanc"))

# deviations
# compare with different aLHII_0
with(F2s, summary(JVPII-JVPII_uc))

# compare with labstaf-estimated aLHII_0
with(F2sb, summary(JVPII-JVPII_uc))

# compare with different aLHII_0 and different Fblanc
with(F2sc, summary(JVPII-JVPII_uc))
}

\seealso{
  \code{\link{fitPI}} for fitting.

  \code{\link{fWebb}} for the various PI functions.

  \code{\link{LABSTAFdata}} for a data set generated with the labSTAF.
  
  \code{\link{FRRFdata}} for fitting with normalized and not normalized data.
  
  \code{\link{modFit}} for options of fitting algorithm
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
