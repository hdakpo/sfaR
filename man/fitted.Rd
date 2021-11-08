\name{fitted}

\alias{fitted}
\alias{fitted.sfacross}
\alias{fitted.lcmcross}
\title{Extract fitted frontier values of classic or latent class stochastic models}

\description{
\code{\link{fitted}} returns the fitted frontier values from classic or latent class 
stochastic frontier models estimated with \code{\link{sfacross}} or \code{\link{lcmcross}}.
}

\usage{
\method{fitted}{sfacross}(object, ...)

\method{fitted}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{...}{Currently ignored.}
}

\value{
In the case of an object of class \code{'sfacross'}, a vector of fitted values is returned.

In the case of an object of class \code{'lcmcross'}, a data frame containing the 
fitted values for each class is returned where each variable terminates with 
\code{"_c#"}, \code{"#"} being the class number.
}

\note{
The fitted values are ordered in the same way as the corresponding observations
in the dataset used for the estimation.
}

\examples{
## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', data = worldprod)
  fit.cb_2c_h <- fitted(cb_2c_h)
  head(fit.cb_2c_h)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
