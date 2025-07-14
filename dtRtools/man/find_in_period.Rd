\name{find_in_period}
\alias{first_in_period}
\alias{last_in_period} 
\alias{count_in_period} 
\alias{all_in_period} 

\title{
  Find time of first, last or all occurrences of a value in a period
}

\description{
    \code{first_in_period} finds the datetime of first occurrence when a value equals a target value within a certain period.
  
    \code{last_in_period} finds the datetime of last occurrence when a value equals a target value within a certain period.
    
    \code{count_in_period} finds the number of times a value equals a target value within a certain period.

    \code{all_in_period} finds the datetime of all occurrence when a value equals a target value.
}

\usage{
first_in_period(x, datetime, value, period, target = 0, by = NULL, 
      type = c("both", "ascending", "descending"))

last_in_period(x, datetime, value, period, target = 0, by = NULL, 
      type = c("both", "ascending", "descending"))

count_in_period(x, datetime, value, period, target = 0, by = NULL, 
      type = c("both", "ascending", "descending"))

all_in_period(x, datetime, value, period, target = 0, by = NULL, 
      type = c("both", "ascending", "descending"))
}

\arguments{
  \item{x }{A data.frame or data.frame extension, then the number of rows should be equal to the length of \code{I}.
  }
  \item{datetime }{the column that has the datetime value. This can be a number, of type \code{posIXct}, or a character string. In the latter case, the datetimes will be converted to  \code{posIXct} using \code{as.Date}.
  }
  \item{period }{the column that defines the period over which the \code{value}s will be matched against \code{target}. 
  }
  \item{value }{the column that has the value to match against \code{target}. 
  }
  \item{target }{the value that needs to be matched. 
  }
  \item{by }{optionally, a column to group by for this operation. 
  }
  \item{type }{allows to select values when the data are in ascending, or descending mode, or both. 
  }
}
\value{
 a \code{data.frame} with several columns : 
 \itemize{
  \item any values as given through the \code{by} argument.
  \item the \code{period} over which the value has been extracted.
  \item a column called \code{ascending} which gives the direction of the data.
  \item a column called \code{first} which gives the datetime of first occurrence of the \code{target} in the period.
  \item the \code{target} value.
 }
 Attributes of the input object are preserved; the data only contain instances where the target has been reached in a period.
}

\details{
 to be added

}

\examples{

# =======================================
# time is expressed in fractional year
# =======================================

# create the timeseries
GNP <- data.frame(time  = seq(from = 1954, to = 1965, , length.out=90), 
                  value = sin(seq(from=1, to=20, length.out=90)) + 
                          round(rnorm(45), 1)
                  )

# year to select from
GNP$year <- as.integer(GNP$time)

# find first occurrence of 0
F0 <- first_in_period(GNP, 
                      datetime = time,     # find first date
                      period   = year,     # in year 
                      value    = value,    # where value
                      target   = 0)        # is equal to 0

# Show results
F0

# find number of occurrences of 0
Fn <- count_in_period(GNP, 
                      datetime = time, 
                      period   = year, 
                      value    = value, 
                      target   = 0)

Fn

# find all occurrences, where trend was descending
F0d <- all_in_period(GNP, 
                     datetime = time, 
                     period   = year, 
                     value    = value, type="descending",
                     target   = 0)

# plot results
plot(GNP$time, GNP$value, type="l")
abline(h = 0,          col = "red")    # target value
abline(v = 1954:1964,  lty = 2)        # beginning of each year
with (F0,  
   points(first, target, col="red"))
with (F0d, 
   points(match, target, col="blue", pch="+"))

# =======================================
# timeseries with actual datetime values
# =======================================

# generate data (high density data)
time  <- seq(from = as.POSIXct("1995-01-01"), 
             to   = as.POSIXct("1998-01-01"), 
             by   = 3600)                      # hourly values
             
value <- sin((1:length(time))/1000)

x     <- data.frame(time  = time, 
                    value = value)
                    
x$year <- as.POSIXlt(x$time)$year + 1900
x$month <- as.POSIXlt(x$time)$mon + 1 

# find first occurrences of value 0.5 *per year*
Fyear <- first_in_period(x, 
                         datetime = time, 
                         value    = value, 
                         period   = year, 
                         target   = 0.5)

# yearly values, ascending
Fa <- first_in_period(x, 
                      datetime = time, 
                      value    = value, 
                      period   = year, 
                      target   = 0.5, 
                      type     = "asc")
                      
plot(time, value, pch = ".")
abline(h = 0.5)   # horizontal line at the target
abline(v = Fyear$first, lwd = 2)
abline(v = Fa$first, col = "red", lty = 2, lwd = 2)

# find first occurrences per *month for each year*

Fmonth <- first_in_period(x, 
                          datetime = time, 
                          value    = value, 
                          period   = month, 
                          by       = year, 
                          target   = 0.5)
# day of the month
Fmonth$dom <- as.POSIXlt(Fmonth$first)$mday # day in month of first occurrence

}

\references{
to be added}
\seealso{
  \code{\link{Veluwe_tree_budstage}} for a data set.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
