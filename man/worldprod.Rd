\name{worldprod}
\alias{worldprod}
\docType{data}
\title{Data on world production}

\description{This dataset provides information on production related variables for eighty-two countries 
over the period 1960–1987.}

\usage{worldprod}

\format{
  A data frame with 2,296 observations on the following 12 variables.
  \describe{
    \item{country}{Country name.}
    \item{code}{Country identification.}
    \item{yr}{Year identification.}
    \item{y}{GDP in 1987 U.S. dollars.}
    \item{k}{Physical capital stock in 1987 U.S. dollars.}
    \item{l}{Labor (number of individuals in the workforce between the ages of 15 and 64).}
    \item{h}{Human capital-adjusted labor.}
    \item{ly}{Log of \code{y}.}
    \item{lk}{Log of \code{k}.}
    \item{ll}{Log of \code{l}.}
    \item{lh}{Log of \code{h}.}
    \item{initStat}{Log of the initial capital to labor ratio of each country, 
       \code{lk} - \code{ll}, measured at the beginning of the sample period.}
  }
}

\details{The dataset is from the World Bank STARS database and has been used in Kumbhakar \emph{et al.} (2014).}

\source{\url{https://sites.google.com/site/sfbook2014/home/for-stata-v12-v13-v14}}

\references{
Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. \emph{A Practitioner's Guide to Stochastic Frontier Analysis Using Stata}. Cambridge University Press.
}

\examples{
str(worldprod)
summary(worldprod)
}

\keyword{datasets}
