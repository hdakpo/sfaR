\name{vcov}

\alias{vcov}
\alias{vcov.sfacross}
\alias{vcov.lcmcross}

\title{Compute variance-covariance matrix of classic or latent class stochastic models}

\description{
\code{\link{vcov}} computes the variance-covariance matrix of the maximum likelihood (ML) 
coefficients of classic or latent class stochastic frontier models estimated by 
\code{\link{sfacross}} or \code{\link{lcmcross}}.
}

\usage{
\method{vcov}{sfacross}(object, extraPar = FALSE, ...)

\method{vcov}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model 
returned by \code{\link{sfacross}} or \code{\link{lcmcross}}}

\item{extraPar}{Logical. Only available for non heteroscedastic models 
returned by \code{\link{sfacross}}. Default = \code{FALSE}. If \code{TRUE}, 
variances and covariances of additional parameters are returned:

\code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}

\code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}

\code{sigmauSq} = \eqn{\exp{(Wu)}} =  \eqn{\exp{(\delta Z_u)}} 

\code{sigmavSq} = \eqn{\exp{(Wv)}} = \eqn{\exp{(\phi Z_v)}}

\code{sigma} = \code{sigmaSq}^0.5

\code{lambda} = \code{lambdaSq}^0.5

\code{sigmau} = \code{sigmauSq}^0.5

\code{sigmav} = \code{sigmavSq}^0.5

\code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})}

\item{...}{Currently ignored}
}

\value{
The variance-covariance matrix of the maximum likelihood coefficients is returned.
}

\details{
The variance-covariance matrix is obtained by the inversion of the
negative Hessian matrix. Depending on the distribution and the
\code{"hessianType"} option, the analytical/numeric Hessian or the bhhh
Hessian or the robust Hessian matrix is evaluated.

The argument \code{extraPar}, is currently available for objects of class \code{'sfacross'}. 
When \code{"extraPar = TRUE"}, the variance-covariance of the additional 
parameters is obtained using the delta method.
}

\examples{
## Using data on Spanish dairy farms
# Cobb Douglas (production function) half normal distribution
cb_s_h <- sfacross(formula = YIT ~ X1 + X2 + X3 + X4, udist = "hnormal", 
    data = dairyspain, S = 1, method = "bfgs")
  vcov(cb_s_h)
  vcov(cb_s_h, extraPar = TRUE)


## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
    data = worldprod, uhet = ~ initStat, S = 1)
  vcov(cb_2c_h)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
