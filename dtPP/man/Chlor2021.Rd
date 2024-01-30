\name{Chlor2021}
\alias{Chlor2021}
\title{
  Timeseries of Chl data from remote sensing.
}

\usage{
data(Chlor2021)
}


\value{
  a list with the following items: 
\itemize{
  \item \code{longitude}, longitude at the middle of the grid in x [dgE], length=60
  \item \code{latitude}, latitude at the middle of the grid in y [dgN], length=24
  \item \code{times}, a data.frame with the initial and final times of the data that are averaged, 12 rows
  \item \code{chlorophyll}, the chlorophyll concentration, [mgChl/m3], dim = (60, 24, 12)
  \item \code{source}, reference to the website from which the data were extracted
  \item \code{asp}, the aspect ratio that can be used for plotting
  }
}

\description{
 Dataset obtained from the ocean productivity website from oregon state university
 http://sites.science.oregonstate.edu/ocean.productivity.  
 
 The chlorophyll data were 1080 by 2160 monthly HDF files from 2021, from VIIRS Data, as from 
 http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.chl.viirs.php
 
 A subset of the monthly data (longitude inbetween 2 and 12 dgE, latitude inbetween 52 & 56dgN) was taken.
 
}

\examples{

# ---------------------------------------------------
# Plot the monthly data from februari till november
# ---------------------------------------------------

mf <- par(mfrow=c(3,3))

image2D(x = Chlor2021$longitude, 
        y = Chlor2021$latitude, 
        z = Chlor2021$chlorophyll[,,2:10], 
        clab = "Chl, mg/m3", 
        xlab = "dgE", ylab = "dgN",
        main = month.abb[2:10],
        log = "c", clim = c(0.1, 50),
        las = 1, asp = Chlor2021$asp)
         
par(mfrow=mf)         
}
\seealso{
  \code{\link{fitPI}} for fitting.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}

