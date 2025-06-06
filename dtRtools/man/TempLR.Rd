\name{Water temperature in the wadden sea}
\alias{Wad_watertemp_LR}
\title{
  Low resolution water temperature data from RWS.
}

\usage{
data(Wad_watertemp_LR)
}

\value{
\code{Wad_watertemp_LR} is a data.frame in \code{wide} format of class \code{dtLife}, with the following items: 

\itemize{
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{BOCHTVWTM, BOOMKDP, ...}, the temperature (dgC) at each of the stations; the column names refer to the stations, as described in \code{attributes(object)$stations}.
  }

The object has the following \code{attributes}:
\itemize{
  \item \code{stations}, description of the stations, with their position
  \item \code{variables}, description of the measured variables, and their units
  \item \code{datasource}, data originator
  \item \code{file}, name of the files from which the data came
  \item \code{processing}, track of manipulations performed on the data
  \item \code{fun}, name of the function used to generate the data
  \item \code{format}, wide format
  }  
}

\description{
\code{TempLR} is a low resolution time series for temperature for all stations in the Waddensea, for 2021, at 0.5 to 1 month intervals.
  
Dataset obtained from the Rijkswaterstaat, downloaded from https://waterinfo.rws.nl/

}

\examples{

# ---------------------------------------------------
# Plot the data
# ---------------------------------------------------

plot(Wad_watertemp_LR, mfrow=c(3,4))

# ---------------------------------------------------
# Show the metadata
# ---------------------------------------------------

meta(Wad_watertemp_LR)

}
\seealso{
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

