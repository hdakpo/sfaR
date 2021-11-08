\name{coef}

\alias{coef}
\alias{coef.summary.sfacross}
\alias{coef.summary.lcmcross}
\alias{coef.sfacross}
\alias{coef.lcmcross}

\title{Extract coefficients of classic or latent class stochastic models}

\description{
From an object of class \code{'summary.sfacross'} or \code{'summary.lcmcross'}, 
\code{\link{coef}} extracts the coefficients, their standard errors, z-values, and
(asymptotic) P-values.

From on object of class \code{'sfacross'} or \code{'lcmcross'}, it extracts only the 
estimated coefficients.
}
\usage{
\method{coef}{summary.sfacross}(object, ...)

\method{coef}{summary.lcmcross}(object, ...)

\method{coef}{sfacross}(object, extraPar = FALSE, ...)

\method{coef}{lcmcross}(object, extraPar = FALSE, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}, or an object of class 
\code{'summary.sfacross'} or \code{'summary.lcmcross'}.}

\item{extraPar}{Logical (default = \code{FALSE}). Only applies to objects of 
class \code{'sfacross'} or \code{'lcmcross'}. If \code{TRUE}, additional parameters 
are returned:

\code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}

\code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}

\code{sigmauSq} = \eqn{\exp{(Wu)}} = \eqn{\exp{(\delta Z_u)}}

\code{sigmavSq} = \eqn{\exp{(Wv)}} = \eqn{\exp{(\phi Z_v)}}

\code{sigma} = \code{sigmaSq}^0.5

\code{lambda} = \code{lambdaSq}^0.5

\code{sigmau} = \code{sigmauSq}^0.5

\code{sigmav} = \code{sigmavSq}^0.5

\code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})
}

\item{...}{Currently ignored.}
}

\value{
For objects of class \code{'summary.sfacross'} or \code{'summary.lcmcross'}, 
\code{\link{coef}} returns a matrix with four columns. Namely, the estimated
coefficients, their standard errors, z-values, and (asymptotic) P-values.

For objects of class \code{'sfacross'} or \code{'lcmcross'}, \code{\link{coef}} 
returns a numeric vector of the estimated coefficients. If \code{extraPar = TRUE},
additional parameters, detailed in the section \sQuote{Arguments}, are also returned.
In the case of object of class \code{'lcmcross'}, each additional parameter terminates 
with \code{"#"} that represents the class number.
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  coef(tl_u_ts, extraPar = TRUE)
  coef(summary(tl_u_ts))
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
