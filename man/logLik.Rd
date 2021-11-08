\name{logLik}

\alias{logLik}
\alias{logLik.sfacross}
\alias{logLik.lcmcross}

\title{Extract log-likelihood value of classic or latent class stochastic models}

\description{
\code{\link{logLik}} extracts the log-likelihood value(s) from classic or latent class 
stochastic frontier models estimated with \code{\link{sfacross}} or \code{\link{lcmcross}}.
}

\usage{
\method{logLik}{sfacross}(object, individual = FALSE, ...)

\method{logLik}{lcmcross}(object, individual = FALSE, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{individual}{Logical. If \code{FALSE} (default), the sum of all observations' log-likelihood 
values is returned. If \code{TRUE}, a vector of each observation's log-likelihood value is returned.}

\item{...}{Currently ignored.}
}

\value{
\code{\link{logLik}} returns an object of class \code{'logLik'}, which is either a
numeric matrix with the log-likelihood value (\code{logLik}), the total number of 
observations (\code{Nobs}) and the number of free parameters (\code{df}), when 
\code{individual = FALSE}, 

or a list of elements, containing the log-likelihood of each observation (\code{logLik}), the 
total number of observations (\code{Nobs}) and the number of free parameters (\code{df}), 
when \code{individual = TRUE}. 
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  logLik(tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal", data = worldprod, S = 1)
  logLik(cb_2c_h, individual = TRUE)
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
