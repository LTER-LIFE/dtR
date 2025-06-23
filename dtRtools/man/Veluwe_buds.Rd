\name{Bud burst data}
\alias{Veluwe_tree_budstage}

\title{
  Bud stage data of trees in the Hoge Veluwe, from the department of Animal Ecology/NIOO-KNAW
}

\usage{
data(Veluwe_tree_budstage)
}


\value{
  a \code{data.frame} with the following columns: 
  \itemize{
  \item \code{ID},             Unique ID of tree (integer)
  \item \code{scientificName}, Vernacular name of species of tree
  \item \code{elevation},      Elevation of tree location in metres above standard mean sea level
  \item \code{longitude},      Longitude (east-west position) of tree location in decimal degrees[dgE]
  \item \code{latitude},       Latitude (north-south position) of tree location in decimal degrees[dgN]
  \item \code{datetime},       Date of measurement
  \item \code{value},          Score indicating the bud burst stage of the tree crown
  }
}

\examples{

## =============================================================
## look at the data
## =============================================================

head(Veluwe_tree_budstage)

summary(Veluwe_tree_budstage)

meta(Veluwe_tree_budstage)$variables
meta(Veluwe_tree_budstage)$taxonomy

table(Veluwe_tree_budstage$scientificName)

## =============================================================
## plot the tree positions
## =============================================================

# this is similar to meta(Veluwe_tree_budstage)$ID

# select the trees, and their position
treeposition <- unique(Veluwe_tree_budstage[, 
                                           c("ID", "longitude", "latitude")]) 
# add the number of measurements 
nData        <- data.frame(table(Veluwe_tree_budstage$ID)) 
nData        <- merge(treeposition, nData, by = 1)  

# plot it (will use Freq as the color variable)
plot_dtmap(nData[, -1], pch = 18,
           clab = "# data", 
           main = "Tree position in the Veluwe_tree_budstage data set")

## =============================================================
## estimate the bud burst date for each tree
## =============================================================

Trees      <- Veluwe_tree_budstage

Trees$year <- as.POSIXlt(Trees$datetime)$year + 1900

# for each year (period) and each tree (by), 
# find the first date (datetime) where the bud stage reaches 1 (target)

firstBuds <- first_in_period(Trees,
                             datetime = datetime,
                             period = year,
                             value  = value,        # bud stage value
                             target = 1,            # 1 = bud burst,
                             type   = "ascending",
                             by     = ID)           # calculate for each tree
head(firstBuds)

# day of the year                             
firstBuds$doy <- as.POSIXlt(firstBuds$first)$yday               

with(firstBuds, boxplot(doy ~ year, 
                        main = "day of budburst", 
                        las  = 1))

# mean doy for each tree
# -----------------------

meandoy <- aggregate(firstBuds$doy, 
                     by  = list(ID = firstBuds$ID), 
                     FUN = mean)
colnames(meandoy)[2] <- "doy"                     

# combine with metadata to have the coordinates
meandoy <- merge(meandoy, 
                 meta(Veluwe_tree_budstage)$ID, 
                 by = "ID")

plot_dtmap(meandoy[, c("longitude", "latitude", "doy")],
           pch = 16, 
           clab   = "doy",
           main   = "Average bud burst day in Veluwe")
       
}

\references{
creator = Department of Animal Ecology, Netherlands Institute of Ecology (NIOO-KNAW), 
Wageningen, Netherlands. 

electronicMailAddress = "AnE_Database@nioo.knaw.nl",

\url{https://nioo.knaw.nl/en/facilities/hole-breeding-passerines-monitoring-hoge-veluwe }
}

\description{
The stage of bud burst has been measured from 1988 onwards on several tree species in different areas across the Netherlands to record the phenology of the trees.

Bud stage was visually assessed twice a week from the beginning of April onwards.

In the years 1989, 1990, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, bud stage was scored on a scale from 0 to 3 with increments of 0.25. 
In all other years, only increments of 0.5 were allowed. 

On this scale:

\itemize{

\item 0 means that buds are not yet green, 
\item 0.5 means buds are swollen, 
\item 1 means buds are green and swollen, 
\item 2 means small leaves are sticking out and 
\item 3 means that leaves are bended at the leafbase. 
}

This score was assigned to both the lower part of the tree and the top of the tree, but for the current dataset only the data from the top of the tree were selected.

publication_date of the full dataset is "2023-12-19".

The current dataset only considers: 

the bud burst stage (PO:0025532) of the tree crown of the two Quercus species, and for individuals located in the Hoge Veluwe.

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
\concept{bud burst}
\concept{trees}
\concept{ecology}
\concept{plant phenology}
