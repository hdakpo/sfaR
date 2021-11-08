\name{sfacross}

\alias{sfacross}
\alias{print.sfacross}

\title{Stochastic frontier estimation using cross-section data}

\description{
\code{\link{sfacross}} is a symbolic formula-based function 
for the estimation of stochastic frontier models in the case of cross-sectional 
or pooled cross-section data, using maximum (simulated) likelihood - M(S)L.

The function accounts for heteroscedasticity in both one-sided and
two-sided error terms as in Reifschneider and Stevenson (1991), Caudill and
Ford (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
\emph{et al.} (1991), Huang and Liu (1994) and Battese and Coelli (1995).

Ten distributions are possible for the one-sided error term and nine
optimization algorithms are available.

The truncated normal - normal distribution with scaling property as
in Wang and Schmidt (2002) is also implemented.
}

\usage{
sfacross(formula, muhet, uhet, vhet, logDepVar = TRUE, data, subset, S = 1, 
  udist = "hnormal", scaling = FALSE, start = NULL, method = "bfgs", hessianType = 1,
  simType = "halton", Nsim = 100, prime = 2, burn = 10, antithetics = FALSE,
  seed = 12345, itermax = 2000, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06,
  stepmax = 0.1, qac = "marquardt")
}

\arguments{
\item{formula}{A symbolic description of the model to be estimated based on
the generic function \code{formula} (see section \sQuote{Details}).}

\item{muhet}{A one-part formula to consider heterogeneity in the mean of the
pre-truncated distribution (see section \sQuote{Details}).}

\item{uhet}{A one-part formula to consider heteroscedasticity in the
one-sided error variance (see section \sQuote{Details}).}

\item{vhet}{A one-part formula to consider heteroscedasticity in the
two-sided error variance (see section \sQuote{Details}).}

\item{logDepVar}{Logical. Informs whether the dependent variable is logged 
(\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.}

\item{data}{The data frame containing the data.}

\item{subset}{An optional vector specifying a subset of observations to be
used in the optimization process.}

\item{S}{If \code{S = 1} (default), a production (profit) frontier is
estimated: \eqn{\epsilon_i = v_i-u_i}. If \code{S = -1}, a cost frontier is
estimated: \eqn{\epsilon_i = v_i+u_i}.}

\item{udist}{Character string. Default = \code{"hnormal"}. Distribution specification 
for the one-sided error term. 10 different distributions are available:
\itemize{
  \item \code{"hnormal"}, for the half normal distribution (Aigner \emph{et al.}
1977, Meeusen and Vandenbroeck 1977)
  \item \code{"exponential"}, for the exponential distribution
  \item \code{"tnormal"} for the truncated normal distribution (Stevenson 1980)
  \item \code{"rayleigh"}, for the Rayleigh distribution (Hajargasht 2015)
  \item \code{"uniform"}, for the uniform distribution (Li 1996, Nguyen 2010)
  \item \code{"gamma"}, for the Gamma distribution (Greene 2003)
  \item \code{"lognormal"}, for the log normal distribution (Migon and Medici
2001, Wang and Ye 2020)
  \item \code{"weibull"}, for the Weibull distribution (Tsionas 2007)
  \item \code{"genexponential"}, for the generalized exponential distribution
(Papadopoulos 2020)
  \item \code{"tslaplace"}, for the truncated skewed Laplace distribution (Wang
2012).
}}

\item{scaling}{Logical. Only when \code{udist = "tnormal"} and 
\code{scaling = TRUE}, the scaling property model (Wang and Schmidt 2002) is
estimated. Default = \code{FALSE}. (see section \sQuote{Details}).}

\item{start}{Numeric vector. Optional starting values for the maximum
likelihood (ML) estimation.}

\item{method}{Optimization algorithm used for the estimation. 
Default = \code{"bfgs"}. 9 algorithms are available:
\itemize{
  \item \code{"bfgs"}, for Broyden-Fletcher-Goldfarb-Shanno
(see \code{\link[maxLik:maxBFGS]{maxBFGS}})
  \item \code{"bhhh"}, for Berndt-Hall-Hall-Hausman
(see \code{\link[maxLik:maxNR]{maxBHHH}})
  \item \code{"nr"}, for Newton-Raphson (see \code{\link[maxLik]{maxNR}})
  \item \code{"nm"}, for Nelder-Mead (see \code{\link[maxLik:maxBFGS]{maxNM}})
  \item \code{"ucminf"}, implements a quasi-Newton type with BFGS updating of the
inverse Hessian and soft line search with a trust region type monitoring of
the input to the line search algorithm (see \code{\link[ucminf]{ucminf}})
  \item \code{"mla"}, for general-purpose optimization based on
Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:marqLevAlg]{mla}})
  \item \code{"sr1"}, for Symmetric Rank 1 (see \code{\link[trustOptim]{trust.optim}})
  \item \code{"sparse"}, for trust regions and sparse Hessian
(see \code{\link[trustOptim]{trust.optim}})
  \item \code{"nlminb"}, for optimization using PORT routines
(see \code{\link[stats]{nlminb}})
}}

\item{hessianType}{Integer. If \code{1} (Default), analytic Hessian is
returned for all the distributions except \code{"gamma"},
\code{"lognormal"} and \code{"weibull"} for which the numeric Hessian is
returned. If \code{2}, bhhh Hessian is estimated (\eqn{g'g}). If \code{3}, 
robust Hessian is computed (\eqn{H^{-1}GH^{-1}}).}

\item{simType}{Character string. If \code{simType = "halton"} (Default), Halton 
draws are used for maximum simulated likelihood (MSL). If \code{simType = "ghalton"}, 
Generalized-Halton draws are used for MSL. If \code{simType = "sobol"}, Sobol
draws are used for MSL. If \code{simType = "uniform"}, uniform draws are used
for MSL. (see section \sQuote{Details}).}

\item{Nsim}{Number of draws for MSL.}

\item{prime}{Prime number considered for Halton and
Generalized-Halton draws. Default = \code{2}.}

\item{burn}{Number of the first observations discarded in the case
of Halton draws. Default = \code{10}.}

\item{antithetics}{Logical. Default = \code{FALSE}. If \code{TRUE}, antithetics counterpart
of the uniform draws is computed. (see section \sQuote{Details}).}

\item{seed}{Numeric. Seed for the random draws.}

\item{itermax}{Maximum number of iterations allowed for
optimization. Default = \code{2000}.}

\item{printInfo}{Logical. Print information during optimization. Default =
\code{FALSE}.}

\item{tol}{Numeric. Convergence tolerance. Default = \code{1e-12}.}

\item{gradtol}{Numeric. Convergence tolerance for gradient. Default = \code{1e-06}.}

\item{stepmax}{Numeric. Step max for \code{ucminf} algorithm. Default = \code{0.1}.}

\item{qac}{Character. Quadratic Approximation Correction for \code{"bhhh"}
and \code{"nr"} algorithms. If \code{"stephalving"}, the step length is decreased but
the direction is kept. If \code{"marquardt"} (default), the step length is decreased
while also moving closer to the pure gradient direction. See \code{\link[maxLik:maxNR]{maxBHHH}} and
\code{\link[maxLik]{maxNR}}.}
}

\details{
The stochastic frontier model is defined as:
\deqn{y_i = \alpha + \mathbf{x}'_i\beta + v_i - Su_i}
\deqn{\epsilon_i = v_i -Su_i}

where \eqn{i} is the observation, \eqn{j} is the class, \eqn{y} is the output 
(cost, revenue, profit), \eqn{x} is the vector of main explanatory variables 
(inputs and other control variables), \eqn{u} is the one-sided error term with 
variance \eqn{\sigma_{u}^2}, and \eqn{v} is the two-sided error term with 
variance \eqn{\sigma_{v}^2}.

\code{S = 1} in the case of production (profit) frontier function and 
\code{S = -1} in the case of cost frontier function.

The model is estimated using maximum likelihood (ML) for most distributions
except the Gamma, Weibull and log-normal distributions for which maximum
simulated likelihood (MSL) is used. For this latter, several draws can be
implemented namely Halton, Generalized Halton, Sobol and uniform. In the case
of uniform draws, antithetics can also be computed: first \code{Nsim/2} draws
are obtained, then the \code{Nsim/2} other draws are obtained as
counterpart of one (\code{1-draw}).

To account for heteroscedasticity in the variance parameters of the error
terms, a single part (right) formula can also be specified. To impose the
positivity to these parameters, the variances are modelled as:
\eqn{\sigma^2_u = \exp{(\delta'Z_u)}} or \eqn{\sigma^2_v = \exp{(\phi'Z_v)}}, 
where \eqn{Z_u} and \eqn{Z_v} are the heteroscedasticity variables (inefficiency 
drivers in the case of \eqn{Z_u}) and \eqn{\delta} and \eqn{\phi} the coefficients. 
In the case of heterogeneity in the truncated mean \eqn{\mu}, it is modelled
as \eqn{\mu=\omega'Z_{\mu}}. The scaling property can be applied for the
truncated normal distribution: \eqn{u \sim h(Z_u, \delta)u} where \eqn{u}
follows a truncated normal distribution \eqn{N^+(\tau, \exp{(cu)})}.

In the case of the truncated normal distribution, the convolution of \eqn{u_i}
and \eqn{v_i} is:

\deqn{f(\epsilon_i)=\frac{1}{\sqrt{\sigma_u^2 + \sigma_v^2}}
  \phi\left(\frac{S\epsilon_i + \mu}{\sqrt{\sigma_u^2 + \sigma_v^2}}\right)
  \Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)/\Phi\left(\frac{\mu}{\sigma_u}\right)}

where

\deqn{\mu_{i*}=\frac{\mu\sigma_v^2 - S\epsilon_i\sigma_u^2}{\sigma_u^2 +
  \sigma_v^2}}

and

\deqn{\sigma_*^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}

In the case of the half normal distribution the convolution is obtained by
setting \eqn{\mu=0}.
}

\value{
\code{\link{sfacross}} returns a list of class \code{'sfacross'} containing
the following elements:

\item{call}{The matched call.}

\item{formula}{The estimated model.}

\item{S}{The argument \code{'S'}. See the section \sQuote{Arguments}.}

\item{typeSfa}{Character string. "Stochastic Production/Profit
Frontier, e = v - u" when \code{S = 1} and "Stochastic Cost
Frontier, e = v + u" when \code{S = -1}.}

\item{Nobs}{Number of observations used for optimization.}

\item{nXvar}{Number of explanatory variables in the production or cost frontier.}

\item{nmuZUvar}{Number of variables explaining heterogeneity in the truncated
mean, only if \code{udist = "tnormal"} or \code{"lognormal"}.}

\item{scaling}{The argument \code{'scaling'}. See the section \sQuote{Arguments}.}

\item{logDepVar}{The argument \code{'logDepVar'}. See the section \sQuote{Arguments}.}

\item{nuZUvar}{Number of variables explaining heteroscedasticity in
the one-sided error term.}

\item{nvZVvar}{Number of variables explaining heteroscedasticity in
the two-sided error term.} 

\item{nParm}{Total number of parameters estimated.}

\item{udist}{The argument \code{'udist'}. See the section \sQuote{Arguments}.}

\item{startVal}{Numeric vector. Starting value for M(S)L estimation.}

\item{dataTable}{A data frame (tibble format) containing information on data 
used for optimization along with residuals and fitted values of the OLS and 
M(S)L estimations, and the individual observation log-likelihood.}

\item{olsParam}{Numeric vector. OLS estimates.}

\item{olsStder}{Numeric vector. Standard errors of OLS estimates.}

\item{olsSigmasq}{Numeric. Estimated variance of OLS random error.}

\item{olsLoglik}{Numeric. Log-likelihood value of OLS estimation.}

\item{olsSkew}{Numeric. Skewness of the residuals of the OLS estimation.}

\item{olsM3Okay}{Logical. Indicating whether the residuals of the OLS estimation
have the expected skewness.}

\item{CoelliM3Test}{Coelli's test for OLS residuals skewness. (See Coelli, 1995).}

\item{AgostinoTest}{D'Agostino's test for OLS residuals skewness. (See
D'Agostino and Pearson, 1973).}

\item{optType}{Optimization algorithm used.}

\item{nIter}{Number of iterations of the ML estimation.}

\item{optStatus}{Optimization algorithm termination message.}

\item{startLoglik}{Log-likelihood at the starting values.}

\item{mlLoglik}{Log-likelihood value of the M(S)L estimation.}

\item{mlParam}{Parameters obtained from M(S)L estimation.}

\item{gradient}{Each variable gradient of the M(S)L estimation.}

\item{gradL_OBS}{Matrix. Each variable individual observation gradient of
the M(S)L estimation.}

\item{gradientNorm}{Gradient norm of the M(S)L estimation.}

\item{invHessian}{Covariance matrix of the parameters obtained from
the M(S)L estimation.}

\item{hessianType}{The argument \code{'hessianType'}. See the section \sQuote{Arguments}.}

\item{mlDate}{Date and time of the estimated model.}

\item{simDist}{The argument \code{'simDist'}, only if \code{udist = "gamma"}, 
\code{"lognormal"} or , \code{"weibull"}. See the section \sQuote{Arguments}.}

\item{Nsim}{The argument \code{'Nsim'}, only if \code{udist = "gamma"}, 
\code{"lognormal"} or , \code{"weibull"}. See the section \sQuote{Arguments}.}

\item{FiMat}{Matrix of random draws used for MSL, only if \code{udist = "gamma"}, 
\code{"lognormal"} or , \code{"weibull"}.}
}

\note{For the Halton draws, the code is adapted from the \pkg{mlogit} package.}


\references{
Aigner, D., Lovell, C. A. K., and Schmidt, P. 1977. Formulation and
estimation of stochastic frontier production function models. \emph{Journal
of Econometrics}, \bold{6}(1), 21--37.

Battese, G. E., and Coelli, T. J. 1995. A model for technical inefficiency
effects in a stochastic frontier production function for panel data.
\emph{Empirical Economics}, \bold{20}(2), 325--332.

Caudill, S. B., and Ford, J. M. 1993. Biases in frontier estimation due to
heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.

Caudill, S. B., Ford, J. M., and Gropper, D. M. 1995. Frontier estimation and
firm-specific inefficiency measures in the presence of heteroscedasticity.
\emph{Journal of Business & Economic Statistics}, \bold{13}(1), 105--111.

Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
\bold{6}:247--268.

D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
\emph{Biometrika}, \bold{60}:613--622.

Greene, W. H. 2003. Simulated likelihood estimation of the normal-Gamma
stochastic frontier function. \emph{Journal of Productivity Analysis},
\bold{19}(2-3), 179--190.

Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
cost function. \emph{Journal of Business & Economic Statistics},
\bold{17}(3), 359--363.

Hajargasht, G. 2015. Stochastic frontiers with a Rayleigh distribution.
\emph{Journal of Productivity Analysis}, \bold{44}(2), 199--208.

Huang, C. J., and Liu, J.-T. 1994. Estimation of a non-neutral stochastic
frontier production function. \emph{Journal of Productivity Analysis},
\bold{5}(2), 171--180.

Kumbhakar, S. C., Ghosh, S., and McGuckin, J. T. 1991) A generalized
production frontier approach for estimating determinants of inefficiency in
U.S. dairy farms. \emph{Journal of Business & Economic Statistics},
\bold{9}(3), 279--286.

Li, Q. 1996. Estimating a stochastic production frontier when the adjusted
error is symmetric. \emph{Economics Letters}, \bold{52}(3), 221--228.

Meeusen, W., and Vandenbroeck, J. 1977. Efficiency estimation from
Cobb-Douglas production functions with composed error. \emph{International
Economic Review}, \bold{18}(2), 435--445.

Migon, H. S., and Medici, E. V. 2001. Bayesian hierarchical models for
stochastic production frontier. Lacea, Montevideo, Uruguay.

Nguyen, N. B. 2010. Estimation of technical efficiency in stochastic
frontier analysis. PhD dissertation, Bowling Green State University, August.

Papadopoulos, A. 2021. Stochastic frontier models using the generalized
exponential distribution. \emph{Journal of Productivity Analysis}, \bold{55}:15--29.

Reifschneider, D., and Stevenson, R. 1991. Systematic departures from the
frontier: A framework for the analysis of firm inefficiency.
\emph{International Economic Review}, \bold{32}(3), 715--723.

Stevenson, R. E. 1980. Likelihood Functions for Generalized Stochastic
Frontier Estimation. \emph{Journal of Econometrics}, \bold{13}(1), 57--66.

Tsionas, E. G. 2007. Efficiency measurement with the Weibull stochastic
frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5),
693--706.

Wang, K., and Ye, X. 2020. Development of alternative stochastic frontier
models for estimating time-space prism vertices. \emph{Transportation}.

Wang, H.J., and Schmidt, P. 2002. One-step and two-step estimation of the
effects of exogenous variables on technical efficiency levels. \emph{Journal
of Productivity Analysis}, \bold{18}:129--144.

Wang, J. 2012. A normal truncated skewed-Laplace model in stochastic
frontier analysis. Master thesis, Western Kentucky University, May.
}

\seealso{
  \code{\link[=summary.sfacross]{summary}} for creating and printing summary results.

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

\examples{
## Using data on fossil fuel fired steam electric power generation plants in U.S.
# Translog (cost function) half normal with heteroscedasticity
tl_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
  summary(tl_u_h)

# Translog (cost function) truncated normal with heteroscedasticity
tl_u_t <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'bhhh')
  summary(tl_u_t)

# Translog (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
    log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
    I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
    udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
    scaling = TRUE, method = 'mla')
  summary(tl_u_ts)

## Using data on Philippine rice producers
# Cobb Douglas (production function) generalized exponential, and Weibull distributions
    cb_p_ge <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) + log(NPK) +
    log(OTHER), udist = 'genexponential', data = ricephil, S = 1, method = 'bfgs')
  summary(cb_p_ge)

## Using data on U.S. electric utility industry
# Cost frontier Gamma distribution
tl_u_g <- sfacross(formula = log(cost/fprice) ~ log(output) + I(log(output)^2) + 
    I(log(lprice/fprice)) + I(log(cprice/fprice)), udist = "gamma", uhet = ~ 1, 
    data = electricity, S = -1, method = "bfgs", simType = "halton", Nsim = 200, 
    hessianType = 2) 
  summary(tl_u_g)
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{models}
