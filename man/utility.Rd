\name{utility}
\alias{utility}
\docType{data}
\title{Data on U.S. electricity generating plants}

\description{This dataset contains data on fossil fuel fired steam electric
power generation plants in the United States between 1986 and 1996.
}

\usage{utility}
\format{
  A data frame with 791 observations on the following 11 variables.
  \describe{
    \item{firm}{Plant identification.}
    \item{year}{Year identification.}
    \item{y}{Net-steam electric power generation in megawatt-hours.}
    \item{regu}{Dummy variable which takes a value equal to 1 if the power plant is in a 
      state which enacted legislation or issued a regulatory order to implement retail access 
      during the sample period, and 0 otherwise.}
    \item{k}{Capital stock.}
    \item{labor}{Labor and maintenance.}
    \item{fuel}{Fuel.}
    \item{wl}{Labor price.}
    \item{wf}{Fuel price.}
    \item{wk}{Capital price.}
    \item{tc}{Total cost.}
  }
}
\details{The dataset has been used in Kumbhakar \emph{et al.} (2014).}

\source{\url{https://sites.google.com/site/sfbook2014/home/for-stata-v12-v13-v14}}

\references{
Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. \emph{A Practitioner's Guide to 
Stochastic Frontier Analysis Using Stata}. Cambridge University Press.
}

\examples{
str(utility)
summary(utility)
}

\keyword{datasets}
