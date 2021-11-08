\name{summary}

\alias{summary}
\alias{summary.sfacross}
\alias{print.summary.sfacross}
\alias{summary.lcmcross}
\alias{print.summary.lcmcross}

\title{Summary of results for classic or latent class stochastic models}

\description{Create and print summary results for classic or latent class stochastic 
models returned by \code{\link{sfacross}} and \code{\link{lcmcross}}.
}

\usage{
\method{summary}{sfacross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.sfacross}(x, digits = max(3, getOption("digits") - 2), ...)

\method{summary}{lcmcross}(object, grad = FALSE, ci = FALSE, ...)

\method{print}{summary.lcmcross}(x, digits = max(3, getOption("digits") - 2), ...)
}

\arguments{
\item{object}{An object of either class \code{'sfacross'}, returned by the
function \code{\link{sfacross}}, or class \code{'lcmcross'}, returned by
the function \code{\link{lcmcross}}.}

\item{grad}{Logical. Default = \code{FALSE}. If \code{TRUE}, the gradient 
for the maximum likelihood (ML) estimates of the different parameters is returned.}

\item{ci}{Logical. Default = \code{FALSE}. If \code{TRUE}, the 95\% confidence 
interval for the different parameters (OLS and ML estimates) is returned.}

\item{...}{Currently ignored.}

\item{x}{An object of either class \code{'summary.sfacross'} or
\code{'summary.lcmcross'}.}

\item{digits}{Numeric. Number of digits displayed in values.}
}

\value{
The \code{\link{summary}} method returns a list of class \code{'summary.sfacross'}
or \code{'summary.lcmcross'} that contains the same elements as an object returned by
\code{\link{sfacross}} or \code{\link{lcmcross}} with the following additional elements:

\item{AIC}{Akaike information criterion.}

\item{BIC}{Bayesian information criterion.}

\item{HQIC}{Hannan-Quinn information criterion.}

\item{sigmavSq}{For \code{object} of class \code{'sfacross'}. Variance of the two-sided error
term (\eqn{\sigma_v^2}).}

\item{sigmauSq}{For \code{object} of class \code{'sfacross'}. Parametrization of the variance 
of the one-sided error term (\eqn{\sigma_u^2}).}

\item{Varu}{For \code{object} of class \code{'sfacross'}. Variance of the one-sided error term.}

\item{THETA}{For \code{object} of class \code{'sfacross'} with \code{"udist = uniform"}. 
\eqn{\Theta} value in the case the uniform distribution is defined as: \eqn{u_i \in [0, \Theta]}.}

\item{Eu}{For \code{object} of class \code{'sfacross'}. Expected unconditional inefficiency.}

\item{Expu}{For \code{object} of class \code{'sfacross'}. Expected unconditional efficiency.}

\item{olsRes}{For \code{object} of class \code{'sfacross'}. Matrix of OLS estimates, their
standard errors, t-values, P-values, and when \code{ci = TRUE} their
confidence intervals.}

\item{mlRes}{Matrix of ML estimates, their standard errors, z-values,
asymptotic P-values, and when \code{grad = TRUE} their
gradient, \code{ci = TRUE} their confidence intervals.}

\item{chisq}{For \code{object} of class \code{'sfacross'}. Chi-square statistics of the
difference between the stochastic frontier and the OLS.}

\item{df}{Degree of freedom for the inefficiency model.}
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  summary(tl_u_ts, grad = TRUE, ci = TRUE)
}

\seealso{
  \code{\link{sfacross}}, for the stochastic frontier analysis model fitting function.

  \code{\link{lcmcross}}, for the latent class stochastic frontier analysis model fitting function.

  \code{\link[=coef.sfacross]{coef}} for extracting coefficients of the estimation.

  \code{\link[=efficiencies.sfacross]{efficiencies}} for computing (in-)efficiency estimates.

  \code{\link[=fitted.sfacross]{fitted}} for extracting the fitted frontier values.

  \code{\link[=ic.sfacross]{ic}} for extracting information criteria.

  \code{\link[=logLik.sfacross]{logLik}} for extracting log-likelihood value(s) of the estimation.

  \code{\link[=marginal.sfacross]{marginal}} for computing marginal effects of inefficiency drivers.

  \code{\link[=residuals.sfacross]{residuals}} for extracting residuals of the estimation.

  \code{\link[=vcov.sfacross]{vcov}} for computing the variance-covariance matrix of the coefficients.
  
  \code{\link{skewnessTest}} for implementing skewness test.
}

\keyword{methods}
