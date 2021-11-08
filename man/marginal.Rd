\name{marginal}

\alias{marginal}
\alias{marginal.sfacross}
\alias{marginal.lcmcross}

\title{Compute marginal effects of the inefficiency drivers in classic or latent class stochastic models}

\description{
This function returns marginal effects of the inefficiency drivers 
from classic or latent class stochastic frontier models estimated with
\code{\link{sfacross}} or \code{\link{lcmcross}.}
}

\usage{
\method{marginal}{sfacross}(object, ...)

\method{marginal}{lcmcross}(object, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned 
by \code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{...}{Currently ignored.}
}

\details{
\code{\link{marginal}} operates in the presence of exogenous variables that 
explain inefficiency, namely the inefficiency drivers (\eqn{uhet = ~ Z_u} or \eqn{muhet = ~
 Z_{mu}}).

Two components are computed for each variable: the marginal effects on the
expected inefficiency (\eqn{\frac{\partial E[u]}{\partial Z_{(m)u}}}) and the
marginal effects on the variance of inefficiency (\eqn{\frac{\partial
  V[u]}{\partial Z_{(m)u}}}).

The model also allows the Wang (2002) parametrization of
\eqn{\mu} and \eqn{\sigma_u^2} by the same vector of exogenous variables.
This double parameterization accounts for non-monotonic relationships
between the inefficiency and its drivers.
}

\value{
\code{\link{marginal}} returns a data frame containing the marginal effects of the \eqn{Z_u} variables
on the expected inefficiency (each variable has the prefix \code{"Eu_"})
and on the variance of the inefficiency (each variable has the prefix
\code{"Vu_"}) is returned.

In the case of the latent class model (LCM), each variable terminates with 
\code{"_c#"} where \code{"#"} is the class number.
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu + wl, uhet = ~ regu + wl, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  marg.tl_u_ts <- marginal(tl_u_ts)
  summary(marg.tl_u_ts)

## Using data on eighty-two countries production (DGP)
# LCM Cobb Douglas (production function) half normal distribution
cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = "hnormal",
    data = worldprod, uhet = ~ initStat + h, S = 1, method = "mla")
  marg.cb_2c_h <- marginal(cb_2c_h)
  summary(marg.cb_2c_h)
}

\references{
Wang, H.J. 2002. Heteroscedasticity and non-monotonic efficiency effects of a
stochastic frontier model. \emph{Journal of Productivity Analysis},
\bold{18}:241--253.
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
