\name{readBathymetry}
\alias{readBathymetry}

\title{
  Reads bathymetric data from the EMODnet bathymetry project.
}

\description{
  \code{readBat} reads netCDF files with bathymetric data, 
  as downloaded from the EMODnet bathymetry project.
}

\usage{
readBathymetry(file, dir="", attr=NULL, lonlim=NULL, latlim=NULL, 
               levels=c(-5, 0, 5, 10, 20), by=1, ...)
}

\arguments{
  \item{file }{The name of a the file to read. This name will be kept as an attribute of the output. 
  }
  \item{dir }{The name of the directory holding the file. 
  }
  \item{lonlim, latlim }{The ranges of longitude, latitude, two valued vector (min, max).
  }
  \item{attr }{A list with information to be added to the attributes of the output.
  }
  \item{levels }{The contour levels.
  }
  \item{by }{To reduce the resolution, set by > 1; e.g. by=2 will reduce the resolution by half.
  }
  \item{... }{Any other arguments passed to the reading function , i.e. dec, quote, fill, etc.
  }
}
\note{
A small example file for bathymetry is in the package's "extdata" directory: 

\itemize{
 \item "BathymetrySmall.nc"
}

The package directory can be grasped with function \code{system.file}. See examples.

}
\note{
In the license text of EMODnet bathymetric data it says:

The data may be used and redistributed for free but is not intended
for legal use, since it may contain inaccuracies. Neither the data
Contributor, ERD, NOAA, nor the United States Government, nor any
of their employees or contractors, makes any warranty, express or
implied, including warranties of merchantability and fitness for a
particular purpose, or assumes any legal liability for the accuracy,
completeness, or usefulness, of this information.
}

\value{
  Function \code{readBathymetry} returns a list with the following items:

\itemize{
 \item latitude, longitude: the coordinates of the depth points, in [dgN], [dgE]; mapping in EPSG:4326
 \item depth: the water depth, [m], estimated as -elevation (Elevation relative to sea level)
 \item source: the name of the originating organization, 
 \item file: the name of the file
 \item variables: the names and units of the variables
 \item asp: aspect ratio for plotting
 \item contours: simple contour lines for plotting
 \item EPSG: the coordinate system (EPSG nr) of the positions,
 \item processing: when the dataset was generated.
}

The object also contains all information passed with argument "attr".
}

\examples{

# path of the external data of the package
fpath <- system.file("extdata", package="dtRtools")

dir(fpath)  # shows all files 

## =============================================================
## Read the bathymetry file
## =============================================================

Bat <- readBathymetry(
            file = "BathymetrySmall.nc", 
            dir  = fpath,
            attr = list(author = "karline soetaert", project = "LTER-life")
               )

names(Bat)

Bat2 <- readBathymetry(
            file   = "BathymetrySmall.nc", 
            dir    = fpath,
            lonlim = c(7.2, 7.4), latlim=c(53.22, 53.30))


## =============================================================
## Plot the bathymetry as an image
## =============================================================

# Use asp to have the correct aspect ratio while plotting:
with(Bat, 
  image2D(x   = longitude, 
          y   = latitude, 
          z   = depth, 
          asp = asp)
    )

# add simple contours
invisible(lapply(Bat$contours, 
                 FUN = lines))


## =============================================================
## Plot the bathymetry as a contours
## =============================================================

with(Bat, 
  contour(x   = longitude, 
          y   = latitude, 
          z   = depth, 
          asp = asp)
    )

}

\references{
"http://www.emodnet-bathymetry.eu/"
}
\seealso{
\link{readKNMI}, \link{readRWS} for other reading functions.
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\details{
The relevant columns of the dutch files mean (taken from the first part of the text files):

          
}
\keyword{utilities}
