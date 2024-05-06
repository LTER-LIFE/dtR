\name{wgs84_to_rd}
\alias{wgs84_to_rd}
\alias{rd_to_wgs84}
\alias{wgs84_to_xy}
\alias{xy_to_wgs84}

\title{
  Converts coordinates from WGS84 to RD- or XY format and vice versa.
}

\description{
  \code{wgs84_to_xy} converts coordinates in UTM-WGS84 (zone31) to XY format (EPSG code should be known). 

  \code{xy_to_wgs84} converts coordinates in XY format (known EPSG) to  UTM-WGS84 (zone31)format. 

  \code{wgs84_to_rd} converts coordinates in UTM-WGS84 (zone31) to Rijksdriehoek (RD) format. 

  \code{rd_to_wgs84} converts coordinates in Rijksdriehoek (RD) to  UTM-WGS84 (zone31)format. 
}

\usage{
wgs84_to_rd (longitude, latitude)
rd_to_wgs84 (X, Y)
wgs84_to_xy (longitude, latitude, EPSG=25831)
xy_to_wgs84 (X, Y, EPSG=25831)
}

\arguments{
  \item{longitude, latitude }{latitude (dgN), and longitude (dgE)  
  }
  \item{X, Y }{positions in RD or XY format (in m). Center =  SE Paris (47dg 58' 29'' N, 3dg 18' 49'' E).
  }
  \item{EPSG }{EPSG number of the XY format. The default is the format as used by RWS.
  }
}
\value{

for \code{rd_to_wgs84} and \code{xy_to_wgs84}:
\itemize{
  \item{longitude, latitude }{latitude (dgN), and longitude (dgE)  
  }
}'

for \code{wgs84_to_rd} and \code{wgs84_to_xy}: 

\itemize{ 
\item{X, Y }{positions in RD or XY format (in m). Center =  SE Paris (47dg 58' 29'' N, 3dg 18' 49'' E).
  }
}

}

\note{
Coordinates in the Netherlands are often expressed in RD (rijksdriehoeksmetingen), 
These are cartesian coordinates, with units of m. The center of RD (where X=0, Y=0) is at 
47dg 58` 29'' N, 3dg 18' 49'' E, located to the SE of Paris.

Rijkswaterstaat often 
}
\examples{

# --------------------------------------
# transformation of one coordinate 
# --------------------------------------

wgs84_to_rd(  6.282173,   58.04366)
rd_to_wgs84(120700.723, 487525.501)

xy_to_wgs84(617198.4, 5869731)
wgs84_to_xy(4.74499, 52.96436)

# --------------------------------------
# transformation of several coordinates
# --------------------------------------

rd_to_wgs84(c(120700, 120800), c(487525, 487525))


#

(RES <- xy_to_wgs84(X=seq(617198.4, 637198.4, length.out=10),
                    Y=5869731))
wgs84_to_xy(RES[,"longitude"], RES[, "latitude"])



}
\references{
F.H. Schreutelkamp, & ir. G.L. Strang van Hees,
 Benaderingsformules voor de transformatie tussen RD- en WGS84-kaartcoordinaten.
 (in Dutch).
}
\seealso{
  \code{\link{Marsdiep}} bathymetric data.
  \code{\link{Wad_weather}} data set with weather information from Waddensea stations.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\keyword{utilities}
