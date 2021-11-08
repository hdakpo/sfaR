\name{ricephil}
\alias{ricephil}
\docType{data}
\title{Data on rice production in the Philippines}
\description{This dataset contains annual data collected from 43 smallholder rice producers 
in the Tarlac region of the Philippines between 1990 and 1997.}

\usage{ricephil}
\format{
  A data frame with 344 observations on the following 17 variables.
  \describe{
      \item{YEARDUM}{Time period (1= 1990, ..., 8 = 1997).}
      \item{FARMERCODE}{Farmer code (1, ..., 43).}
      \item{PROD}{Output (tonnes of freshly threshed rice).}
      \item{AREA}{Area planted (hectares).}
      \item{LABOR}{Labor used (man-days of family and hired labor).}
      \item{NPK}{Fertiliser used (kg of active ingredients).}
      \item{OTHER}{Other inputs used (Laspeyres index = 100 for Farm 17 in 1991).}
      \item{PRICE}{Output price (pesos per kg).}
      \item{AREAP}{Rental price of land (pesos per hectare).}
      \item{LABORP}{Labor price (pesos per hired man-day).}
      \item{NPKP}{Fertiliser price (pesos per kg of active ingredient).}
      \item{OTHERP}{Price of other inputs (implicit price index).}
      \item{AGE}{Age of the household head (years).}
      \item{EDYRS}{Education of the household head (years).}
      \item{HHSIZE}{Household size.}
      \item{NADULT}{Number of adults in the household.}
      \item{BANRAT}{Percentage of area classified as bantog (upland) fields.}
  }
}
\details{This dataset is published as supplement to Coelli \emph{et al.} (2005).
While most variables of this dataset were supplied by the International Rice 
Research Institute (IRRI), some were calculated by Coelli \emph{et al.} (2005, see p. 325--326).
The survey is described in Pandey \emph{et al.} (1999).
}

\source{
Supplementary files for Coelli \emph{et al.} (2005),
\url{http://www.uq.edu.au/economics/cepa/crob2005/software/CROB2005.zip}. 
}
\references{
Coelli, T. J., Rao, D. S. P., O'Donnell, C. J., and Battese, G. E. 2005. \emph{An Introduction to Efficiency and Productivity Analysis},
Springer, New York.

Pandey, S., Masciat, P., Velasco, L, and Villano, R. 1999. Risk analysis of a rainfed rice production system system in Tarlac,
Central Luzon, Philippines. \emph{Experimental Agriculture}, \bold{35}:225--237.
}
\examples{
str(ricephil)
summary(ricephil)
}

\keyword{datasets}
