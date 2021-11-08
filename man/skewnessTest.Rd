\name{skewnessTest}
\alias{skewnessTest}

\title{Skewness test for stochastic frontier models}

\description{\code{\link{skewnessTest}} computes skewness test for stochastic 
frontier models (i.e. objects of class \code{'sfacross'}).}

\usage{
skewnessTest(object, test = "agostino")
}

\arguments{
\item{object}{An object of class \code{'sfacross'}, returned by \code{\link{sfacross}}.}

\item{test}{A character string specifying the test to implement. If \code{"agostino"} 
(default), D'Agostino skewness test is implemented (D'Agostino and Pearson, 1973). 
If \code{"coelli"}, Coelli skewness test is implemented (Coelli, 1995).}
}

\value{
\code{skewnessTest} returns the results of either the D'Agostino's or the 
Coelli's skewness test.
}

\note{
\code{\link{skewnessTest}} is currently only available for object of class \code{'sfacross'}.
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = "tnormal", muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = "mla")
  skewnessTest(tl_u_ts)
  skewnessTest(tl_u_ts, test = "coelli")
}

\references{
Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
\bold{6}:247--268.

D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
\emph{Biometrika}, \bold{60}:613--622.
}

\keyword{methods}
