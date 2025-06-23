\name{Veluwe temperature data}
\alias{Veluwe_mean_tree_temperature}

\title{
  Mean Temperature data accompaning the bud stage data of trees in the Hoge Veluwe, from KNMI
}

\usage{
data(Veluwe_mean_tree_temperature)
}


\value{
  a \code{data.frame} with the following columns: 
  \itemize{
  \item \code{datetime}, date of measurement
  \item \code{temperature}, mean temperature experienced by trees as in \code{Veluwe_tree_budstage}
  }
}

\description{
Mean temperature experienced by trees from the Veluwe bud stage data. 

Based on measured temperatures from 6 stations of the Dutch Meteorological Society (KNMI), the temperature experienced by each tree was first estimated by spatial interpolation. 

Then the average of these temperatures was estimated.
}

\examples{

## =============================================================
## look at the data
## =============================================================

head   (Veluwe_mean_tree_temperature)

summary(Veluwe_mean_tree_temperature)

meta   (Veluwe_mean_tree_temperature)$variables

## =============================================================
## plot the data
## =============================================================

plot(Veluwe_mean_tree_temperature, type = "l")

## =============================================================
## estimate the first day the temperature is above 10 dg
## =============================================================

# remove the incomplete year 1987.
Temp      <- subset(Veluwe_mean_tree_temperature,
                    subset = datetime > "1987-12-31")

Temp$year <- as.POSIXlt(Temp$datetime)$year + 1900

# for each year, find first day where temp = 15 dg (target)
firstWarm <- first_in_period(Temp,
                             datetime = datetime,
                             period = year,
                             value  = temperature,  
                             target = 15)

head(firstWarm)

# day of the year                             
firstWarm$doy <- as.POSIXlt(firstWarm$first)$yday              

with(firstWarm, plot(doy ~ year, 
                     main = "day of first warmth (15dg)", 
                     las  = 1))

abline(with(firstWarm, lm(doy ~ year)),
       lty = 2)
}

\references{
The data were obtained from the KNMI
}

\author{
 \code{authors of the dataset:}
 
 Cherine C. Jantzen https://orcid.org/0009-0006-0723-2682
 
 Stefan J.G. Vriend https://orcid.org/0000-0002-9006-5988
 
 Marcel E. Visser https://orcid.org/0000-0002-1456-1939
 
\code{implementation in R:}
 
 Karline Soetaert https://orcid.org/0000-0003-4603-7100

}
\keyword{data}
\concept{temperature}
\concept{trees}
\concept{ecology}
\concept{plant phenology}
