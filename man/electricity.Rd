\name{electricity}
\alias{electricity}
\docType{data}
\title{Data on U.S. electric power generation}

\description{This dataset is on electric power generation in the United States.}

\usage{electricity}
\format{
  A data frame with 123 observations on the following 9 variables.
  \describe{
    \item{firm}{Firm identification.}
    \item{cost}{Total cost in 1970, MM USD.}
    \item{output}{Output in million KwH.}
    \item{lprice}{Labor price.}
    \item{lshare}{Labor's cost share.}
    \item{cprice}{Capital price.}
    \item{cshare}{Capital's cost share.}
    \item{fprice}{Fuel price.}
    \item{fshare}{Fuel's cost share.}
  }
}
\details{The dataset is from Christensen and Greene (1976) and has also been used in Greene (1990).}

\source{\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}}

\references{
Christensen, L.R., and W.H. Greene. 1976. Economies of scale in US electric power generation. \emph{The Journal of Political Economy}, \bold{84}:655--676.

Greene, W.H. 1990. A Gamma-distributed stochastic frontier model. \emph{Journal of Econometrics}, \bold{46}:141--163.
}

\examples{
str(electricity)
summary(electricity)
}

\keyword{datasets}
