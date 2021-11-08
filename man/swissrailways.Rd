\name{swissrailways}
\alias{swissrailways}
\docType{data}
\title{Data on Swiss railways
}
\description{This dataset is an unbalanced panel of 50 Swiss railway companies over the period 1985-1997.}

\usage{swissrailways}
\format{
  A data frame with 605 observations on the following 42 variables.
  \describe{
    \item{ID}{Firm identification.}
    \item{YEAR}{Year identification.}
    \item{NI}{Number of years observed.}
    \item{STOPS}{Number of stops in network.}
    \item{NETWORK}{Network length (in meters).}
    \item{NARROW_T}{Dummy variable for railroads with narrow track.}
    \item{RACK}{Dummy variable for ‘rack rail’ in network.}
    \item{TUNNEL}{Dummy variable for network with tunnels over 300 meters on average.}
    \item{T}{Time indicator, first year = 0.}
    \item{Q2}{Passenger output – passenger km.}
    \item{Q3}{Freight output – ton km.}
    \item{CT}{Total cost (1,000 Swiss franc).}
    \item{PL}{Labor price.}
    \item{PE}{Electricity price.}
    \item{PK}{Capital price.}
    \item{VIRAGE}{1 for railroads with curvy tracks.}
    \item{LNCT}{Log of \code{CT}/\code{PE}.}
    \item{LNQ2}{Log of \code{Q2}.}
    \item{LNQ3}{Log of \code{Q3}.}
    \item{LNNET}{Log of \code{NETWORK}/1000.}
    \item{LNPL}{Log of \code{PL}/\code{PE}.}
    \item{LNPE}{Log of \code{PE}.}
    \item{LNPK}{Log of \code{PK}/\code{PE}.}
    \item{LNSTOP}{Log of \code{STOPS}.}
    \item{MLNQ2}{Mean of \code{LNQ2}.}
    \item{MLNQ3}{Mean of \code{LNQ3}.}
    \item{MLNNET}{Mean of \code{LNNET}.}
    \item{MLNPL}{Mean of \code{LNPL}.}
    \item{MLNPK}{Mean of \code{LNPK}.}
    \item{MLNSTOP}{Mean of \code{LNSTOP}.}
  }
}

\details{
The dataset is extracted from the annual reports of the Swiss Federal Office of Statistics 
on public transport companies and has been used in Farsi \emph{et al.} (2005).
}
\source{
\url{http://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}

\url{http://people.stern.nyu.edu/wgreene/Microeconometrics.htm}
}
\references{
Farsi, M., M. Filippini, and W. Greene. 2005. Efficiency measurement in network industries: Application to the Swiss railway companies. \emph{Journal of Regulatory Economics}, \bold{28}:69--90.
}
\examples{
str(swissrailways)
}

\keyword{datasets}
