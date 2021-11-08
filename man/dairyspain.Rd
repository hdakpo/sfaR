\name{dairyspain}
\alias{dairyspain}
\docType{data}
\title{Data on Spanish Dairy Farms}

\description{This dataset contains six years of observations on 247 dairy farms in northern Spain, 
drawn from 1993-1998. The original data consist in the farm and year identifications, plus measurements 
on one output (i.e. milk), and four inputs (i.e. cows, land, labor and feed).
}
\usage{dairyspain}
\format{
  A data frame with 1,482 observations on the following 29 variables.
  \describe{
    \item{FARM}{Farm identification.}
    \item{AGEF}{Age of the farmer.}
    \item{YEAR}{Year identification.}
    \item{COWS}{Number of milking cows.}
    \item{LAND}{Agricultural area.}
    \item{MILK}{Milk production.}
    \item{LABOR}{Labor.}
    \item{FEED}{Feed.}
    \item{YIT}{Log of \code{MILK}.}
    \item{X1}{Log of \code{COWS}.}
    \item{X2}{Log of \code{LAND}.}
    \item{X3}{Log of \code{LABOR}.}
    \item{X4}{Log of \code{FEED}.}
    \item{X11}{1/2 * \code{X1}^2.}
    \item{X22}{1/2 * \code{X2}^2.}
    \item{X33}{1/2 * \code{X3}^2.}
    \item{X44}{1/2 * \code{X4}^2.}
    \item{X12}{\code{X1} * \code{X2}.}
    \item{X13}{\code{X1} * \code{X3}.}
    \item{X14}{\code{X1} * \code{X4}.}
    \item{X23}{\code{X2} * \code{X3}.}
    \item{X24}{\code{X2} * \code{X4}.}
    \item{X34}{\code{X3} * \code{X4}.}
    \item{YEAR93}{Dummy for \code{YEAR = 1993}.}
    \item{YEAR94}{Dummy for \code{YEAR = 1994}.}
    \item{YEAR95}{Dummy for \code{YEAR = 1995}.}
    \item{YEAR96}{Dummy for \code{YEAR = 1996}.}
    \item{YEAR97}{Dummy for \code{YEAR = 1997}.}
    \item{YEAR98}{Dummy for \code{YEAR = 1998}.}
  }
}
\details{
This dataset has been used in Alvarez \emph{et al.} (2004). 
The data have been normalized so that the logs of the inputs sum to zero over the 1,482 observations.
}
\source{\url{http://pages.stern.nyu.edu/~wgreene/Econometrics/oldPanelDataSets.htm}}

\references{
Alvarez, A., C. Arias, and W. Greene. 2004. Accounting for unobservables in production models: 
management and inefficiency. \emph{Econometric Society}, \bold{341}:1--20.
}
\examples{
str(dairyspain)
summary(dairyspain)
}

\keyword{datasets}
