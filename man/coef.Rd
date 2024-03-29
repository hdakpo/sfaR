% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef}
\alias{coef}
\alias{coef.sfacross}
\alias{coef.summary.sfacross}
\alias{coef.sfalcmcross}
\alias{coef.summary.sfalcmcross}
\alias{coef.sfaselectioncross}
\alias{coef.summary.sfaselectioncross}
\title{Extract coefficients of stochastic frontier models}
\usage{
\method{coef}{sfacross}(object, extraPar = FALSE, ...)

\method{coef}{summary.sfacross}(object, ...)

\method{coef}{sfalcmcross}(object, extraPar = FALSE, ...)

\method{coef}{summary.sfalcmcross}(object, ...)

\method{coef}{sfaselectioncross}(object, extraPar = FALSE, ...)

\method{coef}{summary.sfaselectioncross}(object, ...)
}
\arguments{
\item{object}{A stochastic frontier model returned by \code{\link{sfacross}},
\code{\link{sfalcmcross}}, or \code{\link{sfaselectioncross}}, or an object
of class \code{'summary.sfacross'}, \code{'summary.sfalcmcross'}, or\cr
\code{'summary.sfaselectioncross'}.}

\item{extraPar}{Logical (default = \code{FALSE}). If \code{TRUE}, additional
parameters are returned:

\code{sigmaSq} = \code{sigmauSq} + \code{sigmavSq}

\code{lambdaSq} = \code{sigmauSq}/\code{sigmavSq}

\code{sigmauSq} = \eqn{\exp{(Wu)}} = \eqn{\exp{(\bm{\delta}' \mathbf{Z}_u)}}

\code{sigmavSq} = \eqn{\exp{(Wv)}} = \eqn{\exp{(\bm{\phi}' \mathbf{Z}_v)}}

\code{sigma} = \code{sigmaSq}^0.5

\code{lambda} = \code{lambdaSq}^0.5

\code{sigmau} = \code{sigmauSq}^0.5

\code{sigmav} = \code{sigmavSq}^0.5

\code{gamma} = \code{sigmauSq}/(\code{sigmauSq} + \code{sigmavSq})}

\item{...}{Currently ignored.}
}
\value{
For objects of class \code{'summary.sfacross'},
\code{'summary.sfalcmcross'}, or \code{'summary.sfaselectioncross'},
\code{\link{coef}} returns a matrix with four columns. Namely, the
estimated coefficients, their standard errors, z-values,
and (asymptotic) P-values.

For objects of class \code{'sfacross'}, \code{'sfalcmcross'}, or
\code{'sfaselectioncross'}, \code{\link{coef}} returns a numeric vector of
the estimated coefficients. If \code{extraPar = TRUE}, additional parameters,
detailed in the section \sQuote{Arguments}, are also returned. In the case
of object of class \code{'sfalcmcross'}, each additional
parameter ends with \code{'#'} that represents the class number.
}
\description{
From an object of class \code{'summary.sfacross'},
\code{'summary.sfalcmcross'}, or \code{'summary.sfaselectioncross'},
\code{\link{coef}} extracts the coefficients,
their standard errors, z-values, and (asymptotic) P-values.

From on object of class \code{'sfacross'}, \code{'sfalcmcross'}, or
\code{'sfaselectioncross'}, it extracts only the estimated coefficients.
}
\examples{

\dontrun{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
scaling = TRUE, method = 'mla')
coef(tl_u_ts, extraPar = TRUE)
coef(summary(tl_u_ts))
}

}
\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model
fitting function using cross-sectional or pooled data.

\code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
model fitting function using cross-sectional or pooled data.

\code{\link{sfaselectioncross}} for sample selection in stochastic frontier
model fitting function using cross-sectional or pooled data.
}
\keyword{coefficients}
\keyword{methods}
