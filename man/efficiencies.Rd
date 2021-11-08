\name{efficiencies}

\alias{efficiencies}
\alias{efficiencies.sfacross}
\alias{efficiencies.lcmcross}

\title{Compute conditional (in-)efficiency estimates of classic or latent class stochastic models}

\description{
\code{\link{efficiencies}} returns (in-)efficiency estimates from classic or latent class stochastic 
frontier models estimated with \code{\link{sfacross}} or \code{\link{lcmcross}.}
}

\usage{
\method{efficiencies}{sfacross}(object, level = 0.95, ...)

\method{efficiencies}{lcmcross}(object, level = 0.95, ...)
}

\arguments{
\item{object}{A classic or latent class stochastic frontier model returned by
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{level}{A number between between 0 and 0.9999 used for the computation
of (in-)efficiency confidence intervals (defaut = \code{0.95}). Only used when \code{udist} =
\code{"hnormal"}, \code{"exponential"}, \code{"tnormal"} or \code{"uniform"} in 
\code{\link{sfacross}} or \code{\link{lcmcross}}.}

\item{...}{Currently ignored.}
}

\value{
A data frame that contains individual (in-)efficiency estimates. These are ordered in 
the same way as the corresponding observations in the dataset used for the estimation.

\bold{- For object of class \code{'sfacross'} the following elements are
returned:}

\item{u}{Conditional inefficiency. In the case argument \code{udist} of
\link{sfacross} is set to \code{"uniform"}, two conditional inefficiency
estimates are returned: \code{u1} for the classic conditional inefficiency
following Jondrow \emph{et al.} (1982), and \code{u2} which is obtained when
\eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen,
2010).}

\item{uLB}{Lower bound for conditional inefficiency. Only when the argument
\code{udist} of \link{sfacross} is set to \code{"hnormal"},
\code{"exponential"}, \code{"tnormal"} or \code{"uniform"}.}

\item{uUB}{Upper bound for conditional inefficiency. Only when the argument
\code{udist} of \link{sfacross} is set to \code{"hnormal"},
\code{"exponential"}, \code{"tnormal"} or \code{"uniform"}.}

\item{teJLMS}{\eqn{\exp{(-u)}}. When the argument \code{udist} of
\link{sfacross} is set to \code{"uniform"}, \code{teJLMS1} = \eqn{\exp{(-u1)}}
and \code{teJLMS2} = \eqn{\exp{(-u2)}}. Only when \code{logDepVar = TRUE}.}

\item{m}{Conditional model. Only when the argument \code{udist} of
\link{sfacross} is set to \code{"hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"rayleigh"}.}

\item{teMO}{\eqn{\exp{(-m)}}. Only when, in the function \link{sfacross},
\code{logDepVar = TRUE} and \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, \code{"uniform"}, or \code{"rayleigh"}.}

\item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
the function \link{sfacross}, \code{logDepVar = TRUE} and \code{udist =
"hnormal"}, \code{"exponential"}, \code{"tnormal"}, \code{"genexponential"},
\code{"rayleigh"}, or \code{"tslaplace"}. In the case \code{udist =
"uniform"}, two conditional efficiency estimates are returned: \code{teBC1}
which is the classic conditional efficiency following Battese and Coelli
(1988) and \code{teBC2} when \eqn{\theta/\sigma_v \longrightarrow \infty}
(see Nguyen, 2010).}

\item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
TRUE} and \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"uniform"}.}

\item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
TRUE} and \code{udist = "hnormal"}, \code{"exponential"},
\code{"tnormal"}, or \code{"uniform"}.}

\bold{- For object of class \code{'lcmcross'} the following elements are
returned:}

\item{Group_c}{Most probable class of each observation.}

\item{PosteriorProb_c}{Highest posterior probability.}

\item{PosteriorProb_c#}{Posterior probability of class #.}

\item{PriorProb_c#}{Prior probability of class #.}

\item{u_c}{Conditional inefficiency of the most probable class given the
posterior probability.}

\item{teJLMS_c}{\eqn{\exp{(-u_c)}}. Only when, in the function \link{lcmcross}, \code{logDepVar = TRUE}.}

\item{u_c#}{Conditional inefficiency associated to class #, regardless of \code{Group_c}.}

\item{ineff_c#}{Conditional inefficiency (\code{u_c}) for observations in class # only.}

}

\details{
The conditional inefficiency is obtained following Jondrow \emph{et
al.} (1982) and the conditional efficiency is computed following Battese
and Coelli (1988). In some cases the conditional mode is also returned
(Jondrow \emph{et al.} 1982). The confidence interval is computed following 
Horrace and Schmidt (1996), Hjalmarsson \emph{et al.} (1996), or Berra and 
Sharma (1999) (see \sQuote{Value} section).

In the case of the half normal distribution for the one-sided error term, the formulae 
are as follows (for notations, see the \sQuote{Details} section of 
\code{\link{sfacross}} or \code{\link{lcmcross}}):

\itemize{
\item The conditional inefficiency is
}
\deqn{E\left[u_i|\epsilon_i\right]=\mu_{i*} +
\sigma_*\frac{\phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}{\Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}}

where

\deqn{\mu_{i*}=\frac{-S\epsilon_i\sigma_u^2}{\sigma_u^2 + \sigma_v^2}}

and

\deqn{\sigma_*^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}
\itemize{
\item The Battese and Coelli (1988) conditional efficiency is obtained by:
}
\deqn{E\left[\exp{\left(-u_i\right)}|\epsilon_i\right] = \exp{\left(-\mu_{i*}+\frac{1}{2}\sigma_*^2\right)}
 \frac{\Phi\left(\frac{\mu_{i*}}{\sigma_*}-\sigma_*\right)}{\Phi\left(\frac{\mu_{i*}}{\sigma_*}\right)}}

\itemize{
\item The conditional mode is computed using:
}

\deqn{M\left[u_i|\epsilon_i\right]= \mu_{i*} \quad For \quad \mu_{i*} > 0}

and

\deqn{M\left[u_i|\epsilon_i\right]= 0 \quad For \quad \mu_{i*} \leq 0}
\itemize{
\item The confidence intervals are obtained with:
}

\deqn{\mu_{i*} + I_L\sigma_* \leq E\left[u_i|\epsilon_i\right] \leq \mu_{i*} + I_U\sigma_*}

with \eqn{LB_i = \mu_{i*} + I_L\sigma_*} and \eqn{UB_i = \mu_{i*} + I_U\sigma_*}

and

\deqn{I_L = \Phi^{-1}\left\{1 - \left(1-\frac{\alpha}{2}\right)\left[1-\Phi\left(-\frac{\mu_{i*}}{\sigma_*}\right)\right]\right\}}

and

\deqn{I_U = \Phi^{-1}\left\{1-\frac{\alpha}{2}\left[1-\Phi\left(-\frac{\mu_{i*}}{\sigma_*}\right)\right]\right\}}

Thus

\deqn{\exp{\left(-UB_i\right)} \leq E\left[\exp{\left(-u_i\right)}|\epsilon_i\right] \leq \exp{\left(-LB_i\right)}}
}

\examples{
## Using data on fossil fuel fired steam electric power generation plants in the U.S.
# Translog SFA (cost function) truncated normal with scaling property
tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) + log(wl/wf) + 
    log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) + I(log(wl/wf) * 
    log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)), udist = "tnormal", 
    muhet = ~ regu, uhet = ~ regu, data = utility, S = -1, scaling = TRUE, method = "mla")
  eff.tl_u_ts <- efficiencies(tl_u_ts)
  head(eff.tl_u_ts)
  summary(eff.tl_u_ts)
}

\references{
Battese, G.E., and T.J. Coelli. 1988. Prediction of firm-level technical
efficiencies with a generalized frontier production function and panel data.
\emph{Journal of Econometrics}, \bold{38}:387--399.

Bera, A.K., and S.C. Sharma. 1999. Estimating production uncertainty in
stochastic frontier production function models. \emph{Journal of Productivity
Analysis}, \bold{12}:187-210.

Hjalmarsson, L., S.C. Kumbhakar, and A. Heshmati. 1996. DEA, DFA and SFA: A
comparison. \emph{Journal of Productivity Analysis}, \bold{7}:303-327.

Horrace, W.C., and P. Schmidt. 1996. Confidence statements for efficiency
estimates from stochastic frontier models. \emph{Journal of Productivity
Analysis}, \bold{7}:257-282.

Jondrow, J., C.A.K. Lovell, I.S. Materov, and P. Schmidt. 1982. On the
estimation of technical inefficiency in the stochastic frontier production
function model. \emph{Journal of Econometrics}, \bold{19}:233--238.

Nguyen, N.B. 2010. Estimation of technical efficiency in stochastic frontier
analysis. PhD Dissertation, Bowling Green State University, August.
}

\seealso{
\code{\link{sfacross}}, for the stochastic frontier analysis model fitting
function.

\code{\link{lcmcross}}, for the latent class stochastic frontier analysis
model fitting function.
}

\author{K Hervé Dakpo, Yann Desjeux and Laure Latruffe}

\keyword{methods}
