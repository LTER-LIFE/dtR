\name{Water2021}
\alias{Water2021}
\alias{Wad_watertemp_HR}
\alias{Wad_waterheight_HR}

\title{
  High- resolution water height and water temperature data from RWS.
}

\usage{
data(Wad_watertemp_HR)
}

\note{

The original reported values were averaged values, 5 minutes before till 5 minutes after the recorded time. These timeseries were in 10 minute intervals; to reduce their size, they were averaged over 20 minute intervals.

The values for water height, reported (as integers) in cm were divided by 100 to convert them to m.
}

\value{

\code{Wad_waterheight_HR}, and \code{Wad_watertemp_HR} are \code{data.frames} in \code{wide} format of class \code{dtLife}, with the following items: 
\itemize{
  \item \code{datetime}, Date and time of the measurement, POSIXct format
  \item \code{...}, the temperature (dgC) or water height (m) at each of the stations; the column names refer to the stations, as described in \code{attributes(object)$stations}.
  }

All objects have the following attributes:
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
 \code{Wad_waterheight_HR}, \code{Wad_watertemp_HR}: high resolution water heights and temperature for all stations in the Waddensea,  and for the first 3 months of 2021, at 0.5 hour intervals.
  
 Water height is versus NAP (Normaal Amsterdams Peil) .
 
 Dataset obtained from the Rijkswaterstaat, downloaded from https://waterinfo.rws.nl/

}

\examples{

# ---------------------------------------------------
# Plot the position of the data
# ---------------------------------------------------

with(Wad_depth, 
   image2D(x = longitude, 
           y = latitude, 
           z = depth, clab = "m", 
           main = "Position of high resolutions temperature stationS",
           asp  = asp, las = 1, 
           colkey = list(length=0.5, width=0.5))
           )
           
stats <- attributes(Wad_watertemp_HR)$stations
points(x = stats$longitude, 
       y = stats$latitude, 
       pch = 18, cex = 3, col = "grey")

# ---------------------------------------------------
# Plot the high-resolution data
# ---------------------------------------------------

mf <- par(mfrow=c(2,1))

HT2021 <- Wad_watertemp_HR[,1:2]
colnames(HT2021)[2] <- "T"
HT2021 <- merge(HT2021, Wad_waterheight_HR[,c("datetime", "DENHDVSGR")])
colnames(HT2021)[3] <- "Height"

with(HT2021, 
   plot(x    = datetime, 
        y    = T, 
        main = "temperature 2021", 
        ylab = "dgC", xlab = "time",
        las = 1, type= "l"))

with(HT2021, 
   plot(x    = datetime, 
        y    = Height, 
        main = "water height 2021", ylab = "m", xlab = "time",
        las = 1, type= "l"))

# subset Januari

Tjan <- subset(HT2021, subset = datetime>= "2021-01-01" &
                                datetime < "2021-02-01")

with(Tjan, 
   plot(x    = datetime, 
        y    = T, 
        main = "temp, Januari 2021", 
        xlab = "time", ylab = "dgC", 
        las = 1, type= "l"))
        
Hjan <- subset(HT2021, subset=datetime>= "2021-01-01" &
                              datetime <"2021-02-01")
with(Hjan, 
   plot(x    = datetime, 
        y    = Height, 
        main = "height, Januari 2021", 
        ylab = "m", xlab = "time",
        las = 1, type= "l"))

# ---------------------------------------------------
# Plot Daily-averaged temperature
# ---------------------------------------------------

# average data over one day
Temp2021_day <- average_timeseries(HT2021, 
                                   avgOver = "day", 
                                   value   = "T")

plot(x = Temp2021_day$datetime, 
     y = Temp2021_day$T, 
     main = "daily-averaged temperature", 
     ylab = "dgC", xlab="time",
     las = 1, type = "l")
        
# ---------------------------------------------------
# Plot Hourly-averaged water height
# ---------------------------------------------------

Height2021_h <- average_timeseries(HT2021, 
                                   avgOver = "hour",
                                   avgTime = 1,
                                   value = "Height")

plot(x    = Height2021_h$datetime, 
     y    = Height2021_h$Height, 
     main = "hourly-averaged water height", 
     ylab = "m", xlab = "time",
     las = 1, type = "l", 
     xlim = as.POSIXct(c("2021-01-01","2021-02-01")))

attributes(Height2021_h)$processing
         
par(mfrow=mf) 

# ---------------------------------------------------
# Plot high-resolution data
# ---------------------------------------------------

plot(Wad_waterheight_HR)

# ---------------------------------------------------
# Taking time-averages
# ---------------------------------------------------

# 30 minute resolution
HH <- average_timeseries(
        input   = Wad_weather, 
        avgOver = "min", 
        avgTime = 30,
        datetime = "datetime",
        value   = "temperature")
head(HH)

# day of the month
H2 <- average_timeseries(
        HT2021[, c("datetime", "Height")], 
        avgOver = "mday", 
        avgTime = 1,
        value   = "Height")
       
head(H2)

# monthly averages:
average_timeseries(
        HT2021[, c("datetime", "Height")], 
        avgOver = "mon", 
        avgTime = 1,
        value   = "Height")

# every 1.5 hours
H15 <- average_timeseries(
         input   = HT2021[, c("datetime", "Height")], 
         avgOver = "hour", 
         avgTime = 1.5,
         value   = "Height")
head(H15)
dim(H15)

# same
H15b <- average_timeseries(
           input   = HT2021[, c("datetime", "Height")], 
           avgOver = "sec",  
           avgTime = 3600*1.5,
           value   = "Height")
head(H15b)
dim(H15b)

}
\seealso{
  \code{\link{Bathymetry}} bathymetric data.
  \code{\link{Wad_sediment}} data set with sediment composition.
  \code{\link{Wad_weather}} data set with weather information.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{dataset}

