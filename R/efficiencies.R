################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Efficiency/Inefficiency estimation                                           #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Generalized Zero Inefficiency Stochastic Frontier Analysis        #
#           -Zero inefficiency Stochastic Frontier                             #
#           -Contaminated noise Stochastic Frontier                            #
#           -Multi-Modal Inefficiency Stochastic Frontier Analysis             #
#           -Stochastic/Deterministic Metafrontier Analysis                    #
#           -Sample selection correction for Stochastic Frontier Model         #
#         + Panel data                                                         #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
# Data: Cross sectional data & Pooled data & Panel data                        #
#------------------------------------------------------------------------------#

#' Compute conditional (in-)efficiency estimates of stochastic frontier models
#' 
#' \code{\link{efficiencies}} returns (in-)efficiency estimates of models 
#' estimated with \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' \code{\link{sfagzisfcross}}, \code{\link{sfacnsfcross}}, 
#' \code{\link{sfamisfcross}}, \code{\link{sfazisfcross}}, 
#' \code{\link{sfametacross}}, \code{\link{sfaselectioncross}}, 
#' \code{\link{sfapanel1}}, or \code{\link{sfalcmpanel}}.
#' 
#' @name efficiencies
#' 
#' @aliases efficiencies.sfacross efficiencies.sfalcmcross 
#' efficiencies.sfagzisfcross efficiencies.sfacnsfcross 
#' efficiencies.sfamisfcross efficiencies.sfazisfcross 
#' efficiencies.sfaselectioncross efficiencies.sfametacross
#' efficiencies.sfapanel1 efficiencies.sfalcmpanel
#' 
#' @details In general, the conditional inefficiency is obtained following 
#' Jondrow \emph{et al.} (1982) and the conditional efficiency is computed 
#' following Battese and Coelli (1988). In some cases the conditional mode is 
#' also returned (Jondrow \emph{et al.} 1982). The confidence interval is 
#' computed following Horrace and Schmidt (1996), Hjalmarsson \emph{et al.} 
#' (1996), or Berra and Sharma (1999) (see \sQuote{Value} section).
#'
#' In the case of the half normal distribution for the one-sided error term,
#' the formulae are as follows (for notations, see the \sQuote{Details} section
#' of \code{\link{sfacross}} or \code{\link{sfalcmcross}}):
#'
#' \itemize{ \item The conditional inefficiency is: }
#' 
#' \deqn{E\left\lbrack u_i|\epsilon_i\right
#' \rbrack=\mu_{i\ast} + \sigma_\ast\frac{\phi
#' \left(\frac{\mu_{i\ast}}{\sigma_\ast}\right)}{
#' \Phi\left(\frac{\mu_{i\ast}}{\sigma_\ast}\right)}}
#'
#' where
#' 
#' \deqn{\mu_{i\ast}=\frac{-S\epsilon_i\sigma_u^2}{ \sigma_u^2 + \sigma_v^2}}
#'
#' and
#' 
#' \deqn{\sigma_\ast^2 = \frac{\sigma_u^2 \sigma_v^2}{\sigma_u^2 + \sigma_v^2}}
#' 
#' \itemize{ \item The Battese and Coelli (1988) conditional efficiency is
#' obtained with: } 
#' 
#' \deqn{E\left\lbrack\exp{\left(-u_i\right)}
#' |\epsilon_i\right\rbrack = \exp{\left(-\mu_{i\ast}+
#' \frac{1}{2}\sigma_\ast^2\right)}\frac{\Phi\left(
#' \frac{\mu_{i\ast}}{\sigma_\ast}-\sigma_\ast\right)}{
#' \Phi\left(\frac{\mu_{i\ast}}{\sigma_\ast}\right)}}
#' 
#' \itemize{ \item The reciprocal of the Battese and Coelli (1988) conditional 
#' efficiency is obtained with: } 
#' 
#' \deqn{E\left\lbrack\exp{\left(u_i\right)}
#' |\epsilon_i\right\rbrack = \exp{\left(\mu_{i\ast}+
#' \frac{1}{2}\sigma_\ast^2\right)} \frac{\Phi\left(
#' \frac{\mu_{i\ast}}{\sigma_\ast}+\sigma_\ast\right)}{
#' \Phi\left(\frac{\mu_{i\ast}}{\sigma_\ast}\right)}}
#' 
#' \itemize{ \item The conditional mode is computed using: }
#' 
#' \deqn{M\left\lbrack u_i|\epsilon_i\right
#' \rbrack= \mu_{i\ast} \quad \hbox{For} \quad 
#' \mu_{i\ast} > 0}
#' 
#' and
#' 
#' \deqn{M\left\lbrack u_i|\epsilon_i\right
#' \rbrack= 0 \quad \hbox{For} \quad \mu_{i\ast} \leq 0}
#' 
#' \itemize{ \item The confidence intervals are obtained with: }
#' 
#' \deqn{\mu_{i\ast} + I_L\sigma_\ast \leq 
#' E\left\lbrack u_i|\epsilon_i\right\rbrack \leq 
#' \mu_{i\ast} + I_U\sigma_\ast }
#'
#' with \eqn{LB_i = \mu_{i*} + I_L\sigma_*} and 
#' \eqn{UB_i = \mu_{i*} + I_U\sigma_*}
#'
#' and
#' 
#' \deqn{I_L = \Phi^{-1}\left\lbrace 1 -
#' \left(1-\frac{\alpha}{2}\right)\left\lbrack 1-
#' \Phi\left(-\frac{\mu_{i\ast}}{\sigma_\ast}\right)
#' \right\rbrack\right\rbrace }
#'
#' and
#' 
#' \deqn{I_U = \Phi^{-1}\left\lbrace 1-
#' \frac{\alpha}{2}\left\lbrack 1-\Phi
#' \left(-\frac{\mu_{i\ast}}{\sigma_\ast}\right)
#' \right\rbrack\right\rbrace}
#'
#' Thus
#' 
#' \deqn{\exp{\left(-UB_i\right)} \leq E\left
#' \lbrack\exp{\left(-u_i\right)}|\epsilon_i\right\rbrack 
#' \leq\exp{\left(-LB_i\right)}}
#' 
#' In the case of the sample selection, as underlined in Greene (2010), the 
#' conditional inefficiency could be computed using Jondrow \emph{et al.} 
#' (1982). However, here the conditional (in)efficiency is obtained using the 
#' properties of the closed skew-normal (CSN) distribution (Lai, 2015). The 
#' conditional efficiency can be obtained using the moment generating functions 
#' of a CSN distribution (see Gonzalez-Farias \emph{et al.} (2004)). We have:
#' 
#' \deqn{E\left\lbrack\exp{\left(tu_i\right)}
#' |\epsilon_i\right\rbrack = M_{u|\epsilon}(t)=\frac{\Phi_2\left(
#' \tilde{\mathbf{D}} \tilde{\bm{\Sigma}}t; \tilde{\bm{\kappa}}, 
#' \tilde{\bm{\Delta}} + \tilde{\mathbf{D}}\tilde{\bm{\Sigma}}\tilde{
#' \mathbf{D}}' \right)}{\Phi_2\left(\mathbf{0}; \tilde{\bm{\kappa}}, 
#' \tilde{\bm{\Delta}} + \tilde{\mathbf{D}}\tilde{\bm{\Sigma}}\tilde{
#' \mathbf{D}}'\right)}\exp{\left(t\tilde{\bm{\pi}} + \frac{1}{2}t^2\tilde{
#' \bm{\Sigma}}\right)}}
#' 
#' where \eqn{\tilde{\bm{\pi}} = \frac{-S\epsilon_i\sigma_u^2}{\sigma_v^2 + 
#' \sigma_u^2}}, \eqn{\tilde{\bm{\Sigma}} = \frac{\sigma_v^2\sigma_u^2}{
#' \sigma_v^2 + \sigma_u^2}}, \eqn{\tilde{\mathbf{D}} = \begin{pmatrix} 
#' \frac{S\rho}{\sigma_v} \\ 1 \end{pmatrix}}, \eqn{\tilde{\bm{\kappa}} = 
#' \begin{pmatrix} - \mathbf{Z}'_{si}\bm{\gamma} - 
#' \frac{\rho\sigma_v\epsilon_i}{\sigma_v^2 + \sigma_u^2}\\ 
#' \frac{S\sigma_u^2\epsilon_i}{\sigma_v^2 + \sigma_u^2} \end{pmatrix}}, 
#' \eqn{\tilde{\bm{\Delta}} = \begin{pmatrix}1-\rho^2 & 0 \\ 0 & 0
#' \end{pmatrix}}.
#' 
#' The derivation of the efficiency and the reciprocal efficiency is obtained by 
#' replacing \eqn{t = -1} and \eqn{t =1}, respectively. To obtain the 
#' inefficiency as \eqn{E\left[u_i|\epsilon_i\right]} is more complicated as it 
#' requires the derivation of a multivariate normal cdf. We have:
#' 
#' \deqn{E\left[u_i|\epsilon_i\right] = \left. \frac{\partial M_{u|
#' \epsilon}(t)}{\partial t}\right\rvert_{t = 0}}
#' 
#' Then
#' 
#' \deqn{E\left[u_i|\epsilon_i\right] = \tilde{\bm{\pi}} + 
#' \left(\tilde{\mathbf{D}}\tilde{\bm{\Sigma}}\right)'\frac{\Phi_2^*
#' \left(\mathbf{0}; \tilde{\bm{\kappa}}, \ddot{\bm{\Delta}}\right)}{
#' \Phi_2\left(\mathbf{0}; \tilde{\bm{\kappa}}, \ddot{\bm{\Delta}}\right)}}
#' 
#' where \eqn{\Phi_2^* \left(\mathbf{s}; \tilde{\bm{\kappa}}, \ddot{\bm{\Delta}}
#' \right)= \frac{\partial \Phi_2\left(\mathbf{s}; \tilde{\bm{\kappa}}, 
#' \ddot{\bm{\Delta}} \right)}{\partial \mathbf{s}}}
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#' @param level A number between between 0 and 0.9999 used for the computation
#' of (in-)efficiency confidence intervals (defaut = \code{0.95}). Only used
#' when \code{udist} = \code{'hnormal'}, \code{'exponential'}, \code{'tnormal'}
#' or \code{'uniform'} in \code{\link{sfacross}}, \code{\link{sfapanel1}}, and 
#' \code{\link{sfametacross}}.
#' @param newData Optional data frame that is used to calculate the efficiency 
#' estimates. If NULL (the default), the efficiency estimates are calculated 
#' for the observations that were used in the estimation. In the case of object 
#' of class \code{sfaselectioncross} 
#' @param ... Currently ignored.
#'
#' @return A data frame that contains individual (in-)efficiency estimates.
#' These are ordered in the same way as the corresponding observations in the
#' dataset used for the estimation.
#' 
#' \bold{- For object of class \code{'sfacross'} and \code{'sfapanel1'} the 
#' following elements are returned:}
#'
#' \item{u}{Conditional inefficiency. In the case argument \code{udist} of
#' \link{sfacross} is set to \code{'uniform'}, two conditional inefficiency
#' estimates are returned: \code{u1} for the classic conditional inefficiency
#' following Jondrow \emph{et al.} (1982), and \code{u2} which is obtained when
#' \eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen, 2010).}
#'
#' \item{uLB}{Lower bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross}, or \link{sfapanel1} is set to 
#' \code{'hnormal'}, \code{'exponential'}, \code{'tnormal'} or 
#' \code{'uniform'}.}
#'
#' \item{uUB}{Upper bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross}, or \link{sfapanel1} is set to 
#' \code{'hnormal'}, \code{'exponential'}, \code{'tnormal'} or 
#' \code{'uniform'}.}
#'
#' \item{teJLMS}{\eqn{\exp{(-E[u|\epsilon])}}. When the argument \code{udist} of
#' \link{sfacross}, or \link{sfapanel1} is set to \code{'uniform'}, 
#' \code{teJLMS1} = \eqn{\exp{(-E[u_1|\epsilon])}} and \code{teJLMS2} = 
#' \eqn{\exp{(-E[u_2|\epsilon])}}. Only when \code{logDepVar = TRUE}.}
#'
#' \item{m}{Conditional model. Only when the argument \code{udist} of
#' \link{sfacross}, or \link{sfapanel1} is set to \code{'hnormal'}, 
#' \code{'exponential'}, \code{'tnormal'}, or \code{'rayleigh'}.}
#'
#' \item{teMO}{\eqn{\exp{(-m)}}. Only when, in the function \link{sfacross}, or
#' \link{sfapanel1} \code{logDepVar = TRUE} and \code{udist = 'hnormal'}, 
#' \code{'exponential'}, \code{'tnormal'}, \code{'uniform'}, or 
#' \code{'rayleigh'}.}
#'
#' \item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
#' the function \link{sfacross}, or \link{sfapanel1} \code{logDepVar = TRUE}. 
#' In the case \code{udist = 'uniform'}, two conditional efficiency estimates 
#' are returned: \code{teBC1} which is the classic conditional efficiency 
#' following Battese and Coelli (1988) and \code{teBC2} when 
#' \eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen, 2010).}
#' 
#' \item{teBC_reciprocal}{Reciprocal of Battese and Coelli (1988) conditional 
#' efficiency. Similar to \code{teBC} except that it is computed as 
#' \eqn{E\left[\exp{(u)}|\epsilon\right]}.}
#'
#' \item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, or \link{sfapanel1} 
#' \code{logDepVar = TRUE} and \code{udist = 'hnormal'}, \code{'exponential'}, 
#' \code{'tnormal'}, or \code{'uniform'}.}
#'
#' \item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, or \link{sfapanel1}
#' \code{logDepVar = TRUE} and \code{udist = 'hnormal'}, \code{'exponential'}, 
#' \code{'tnormal'}, or \code{'uniform'}.}
#' 
#' \item{theta}{In the case \code{udist = 'uniform'}. \eqn{u \in [0, \theta]}.}
#' 
#' \bold{- For object of class \code{'sfalcmcross'}, \code{'sfagzisfcross'}, 
#' \code{'sfacnsfcross'}, \code{'sfamisfcross'}, \code{'sfazisfcross'}, or 
#' \code{'sfalcmpanel'} the following elements are returned:}
#'
#' \item{Group_c}{Most probable class for each observation.}
#'
#' \item{PosteriorProb_c}{Highest posterior probability.}
#' 
#' \item{u_c}{Conditional inefficiency of the most probable class given the
#' posterior probability.}
#' 
#' \item{teJLMS_c}{\eqn{\exp{(-E[u_c|\epsilon_c])}}. Only when 
#' \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_c}{\eqn{E\left[\exp{(-u_c)}|\epsilon_c\right]}. Only when 
#' \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c}{\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}. Only 
#' when \code{logDepVar = TRUE}.}
#'
#' \item{PosteriorProb_c#}{Posterior probability of class #.}
#'
#' \item{PriorProb_c#}{Prior probability of class #.}
#'
#' \item{u_c#}{Conditional inefficiency associated to class #, regardless of
#' \code{Group_c}.}
#' 
#' \item{teBC_c#}{Conditional efficiency 
#' (\eqn{E\left[\exp{(-u_c)}|\epsilon_c\right]}) associated to class #, 
#' regardless of \code{Group_c}. Only when \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c#}{Reciprocal conditional efficiency 
#' (\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}) associated to class #, 
#' regardless of \code{Group_c}. Only when \code{logDepVar = TRUE}.}
#'
#' \item{ineff_c#}{Conditional inefficiency (\code{u_c}) for observations in
#' class # only.}
#' 
#' \item{effBC_c#}{Conditional efficiency (\code{teBC_c}) for observations in
#' class # only.}
#' 
#' \item{ReffBC_c#}{Reciprocal conditional efficiency (\code{teBC_reciprocal_c}) 
#' for observations in class # only.}
#' 
#' \item{theta_c#}{In the case \code{udist = 'uniform'}. 
#' \eqn{u \in [0, \theta_{c\#}]}.}
#' 
#' \bold{- For object of class \code{'sfaselectioncross'} the following elements 
#' are returned:}
#' 
#' \item{u}{Conditional inefficiency.}
#'
#' \item{teJLMS}{\eqn{\exp{(-E[u|\epsilon])}}. Only when 
#' \code{logDepVar = TRUE}.}
#'
#' \item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
#' the function \link{sfaselectioncross}, 
#' \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal}{Reciprocal of Battese and Coelli (1988) conditional 
#' efficiency. Similar to \code{teBC} except that it is computed as 
#' \eqn{E\left[\exp{(u)}|\epsilon\right]}.}
#' 
#' \bold{- For object of class \code{'sfametacross'}, the same elements as for 
#' object of class \code{'sfacross'} are returned in addition to the following 
#' additional elements are returned:}
#' 
#' \item{MD}{Metafrontier distance, which is the distance in terms of 
#' inefficiency between the group frontier and the metafrontier.}
#' 
#' \item{MTI}{Metafrontier total inefficiency, which equals the metafrontier
#' distance plus the group inefficiency with its own frontier.}
#' 
#' When \code{logDepVar = TRUE}, additional elements are also returned
#' 
#' \item{TGR_teJLMS}{Technology gap ratio computed using the group efficiency
#' scores computed as \eqn{\exp{\left(E[-u|\epsilon]\right)}}, when 
#' \code{modelType} is any of "hhl14", "bpo04a", "bpo04b", "aos17b", or 
#' "aos17d".}
#' 
#' \item{TGR_teBC}{Technology gap ratio computed using the group efficiency
#' scores computed as \eqn{E\left[\exp{\left(-u\right)}|\epsilon\right]}, when 
#' \code{modelType} is any of "hhl14", "bpo04a", "bpo04b", "aos17b", or 
#' "aos17d".}
#' 
#' \item{MTE_teJLMS}{Metatechnology efficiency computed using the group 
#' efficiency scores computed as \eqn{\exp{\left(E[-u|\epsilon]\right)}}, when 
#' \code{modelType} is any of "hhl14", "bpo04a", "bpo04b", "aos17b", or 
#' "aos17d".}
#' 
#' \item{MTE_teBC}{Metatechnology efficiency computed using the group 
#' efficiency scores computed as \eqn{\exp{\left(E[-u|\epsilon]\right)}}, when 
#' \code{modelType} is any of "hhl14", "bpo04a", "bpo04b", "aos17b", or 
#' "aos17d".}
#' 
#' \item{TGR}{Technology gap ratio, when \code{modelType} is any of "aos17a", 
#' or "aos17c". These two models are based on the use of the unconditional 
#' efficiency.}
#' 
#' \item{MTE}{Metatechnology efficiency, when \code{modelType} is any of 
#' "aos17a", or "aos17c". These two models are based on the use of the 
#' unconditional efficiency.}
#' 
#' 
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfagzisfcross}}, for the generalized zero inefficiency
#'  stochastic frontier analysis model fitting function using cross-sectional or 
#'  pooled data.
#'  
#' \code{\link{sfacnsfcross}}, for the contaminated noise stochastic 
#' frontier analysis model fitting function using cross-sectional data.
#' 
#' \code{\link{sfamisfcross}}, for the multi-modal inefficiency stochastic 
#' frontier analysis model fitting function using cross-sectional data.
#' 
#' \code{\link{sfazisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function using cross-sectional data.
#' 
#' \code{\link{sfametacross}}, for fitting different metafrontier models
#' using cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfapanel1}}, for the first generation stochastic frontier 
#' analysis model fitting function using panel data.
#' 
#' \code{\link{sfalcmpanel}}, for the latent class stochastic frontier analysis
#' model fitting function using panel data.
#'
#' @references Battese, G.E., and T.J. Coelli. 1988. Prediction of firm-level
#' technical efficiencies with a generalized frontier production function and
#' panel data. \emph{Journal of Econometrics}, \bold{38}:387--399.
#'
#' Bera, A.K., and S.C. Sharma. 1999. Estimating production uncertainty in
#' stochastic frontier production function models. \emph{Journal of
#' Productivity Analysis}, \bold{12}:187-210.
#' 
#' Gonzalez-Farias, G., Dominguez-Molina, A., Gupta, A. K., 2004. Additive 
#' properties of skew normal random vectors. 
#' \emph{Journal of Statistical Planning and Inference}. \bold{126}: 521-534.
#' 
#' Greene, W., 2010. A stochastic frontier model with correction 
#' for sample selection. \emph{Journal of Productivity Analysis}. \bold{34}, 
#' 15--24.
#'
#' Hjalmarsson, L., S.C. Kumbhakar, and A. Heshmati. 1996. DEA, DFA and SFA: A
#' comparison. \emph{Journal of Productivity Analysis}, \bold{7}:303-327.
#'
#' Horrace, W.C., and P. Schmidt. 1996. Confidence statements for efficiency
#' estimates from stochastic frontier models. \emph{Journal of Productivity
#' Analysis}, \bold{7}:257-282.
#'
#' Jondrow, J., C.A.K. Lovell, I.S. Materov, and P. Schmidt. 1982. On the
#' estimation of technical inefficiency in the stochastic frontier production
#' function model. \emph{Journal of Econometrics}, \bold{19}:233--238.
#' 
#' Lai, H. P., 2015. Maximum likelihood estimation of the stochastic frontier 
#' model with endogenous switching or sample selection. 
#' \emph{Journal of Productivity Analysis}, \bold{43}: 105-117.
#'
#' Nguyen, N.B. 2010. Estimation of technical efficiency in stochastic frontier
#' analysis. PhD Dissertation, Bowling Green State University, August.
#'
#' @examples
#' 
#' \dontrun{
#' # Using data on fossil fuel fired steam electric power generation plants in 
#' # the U.S.
#' ## Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) + 
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) + 
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)), 
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1, 
#' scaling = TRUE, method = 'mla')
#' eff.tl_u_ts <- efficiencies(tl_u_ts)
#' head(eff.tl_u_ts)
#' summary(eff.tl_u_ts)
#' }
#' 
#' @export
#' @export efficiencies
# @exportS3Method efficiencies sfacross conditional efficiencies sfacross
# ----------
efficiencies.sfacross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$udist == "hnormal") {
    EffRes <- chalfnormeff(object = object, level = level)
  } else {
    if (object$udist == "exponential") {
      EffRes <- cexponormeff(object = object, level = level)
    } else {
      if (object$udist == "gamma") {
        EffRes <- cgammanormeff(object = object, level = level)
      } else {
        if (object$udist == "rayleigh") {
          EffRes <- craynormeff(object = object, level = level)
        } else {
          if (object$udist == "uniform") {
          EffRes <- cuninormeff(object = object, level = level)
          } else {
          if (object$udist == "tnormal") {
            if (object$scaling) {
            EffRes <- ctruncnormscaleff(object = object, level = level)
            } else {
            EffRes <- ctruncnormeff(object = object, level = level)
            }
          } else {
            if (object$udist == "lognormal") {
            EffRes <- clognormeff(object = object, level = level)
            } else {
            if (object$udist == "genexponential") {
              EffRes <- cgenexponormeff(object = object, level = level)
            } else {
              if (object$udist == "tslaplace") {
              EffRes <- ctslnormeff(object = object, level = level)
              } else {
              if (object$udist == "weibull") {
                EffRes <- cweibullnormeff(object = object, level = level)
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffRes)
}

# conditional efficiencies sfalcmcross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfalcmcross
efficiencies.sfalcmcross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the LCM model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  EffRes <- eval(parse(text = paste0("cLCM", object$nClasses, "Chalfnormeff(object = object, level = level)")))
  return(EffRes)
}

# conditional efficiencies sfagzisfcross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfagzisfcross
efficiencies.sfagzisfcross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the GZISF model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  EffRes <- eval(parse(text = paste0("cGZISF", object$nClasses, "Chalfnormeff(object = object, level = level)")))
  return(EffRes)
}

# conditional efficiencies sfacnsfcross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfacnsfcross
efficiencies.sfacnsfcross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the CNSF model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  if (object$sigmauType == "common") {
    if (object$udist == "hnormal") {
      EffRes <- eval(parse(text = paste0("ccnsfhalfnormeff_", object$linkF,
        "(object = object, level = level)")))
    } else {
      if (object$udist == "exponential") {
        eval(parse(text = paste0("ccnsfexponormeff_", object$linkF, "(object = object, level = level)")))
      } else {
        if (object$udist == "gamma") {
          eval(parse(text = paste0("ccnsfgammanormeff_", object$linkF, "(object = object, level = level)")))
        } else {
          if (object$udist == "rayleigh") {
          eval(parse(text = paste0("ccnsfraynormeff_", object$linkF, "(object = object, level = level)")))
          } else {
          if (object$udist == "uniform") {
            eval(parse(text = paste0("ccnsfuninormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
            if (object$udist == "tnormal") {
            eval(parse(text = paste0("ccnsftruncnormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
            if (object$udist == "lognormal") {
              eval(parse(text = paste0("ccnsflognormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "genexponential") {
              eval(parse(text = paste0("ccnsfgenexponormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "tslaplace") {
                eval(parse(text = paste0("ccnsftslnormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
                if (object$udist == "weibull") {
                eval(parse(text = paste0("ccnsfweibullnormeff_",
                  object$linkF, "(object = object, level = level)")))
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  } else {
    if (object$sigmauType == "different") {
      if (object$udist == "hnormal") {
        EffRes <- eval(parse(text = paste0("cmcesfhalfnormeff_", object$linkF,
          "(object = object, level = level)")))
      } else {
        if (object$udist == "exponential") {
          eval(parse(text = paste0("cmcesfexponormeff_", object$linkF, "(object = object, level = level)")))
        } else {
          if (object$udist == "gamma") {
          eval(parse(text = paste0("cmcesfgammanormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
          if (object$udist == "rayleigh") {
            eval(parse(text = paste0("cmcesfraynormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
            if (object$udist == "uniform") {
            eval(parse(text = paste0("cmcesfuninormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
            if (object$udist == "tnormal") {
              eval(parse(text = paste0("cmcesftruncnormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "lognormal") {
              eval(parse(text = paste0("cmcesflognormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "genexponential") {
                eval(parse(text = paste0("cmcesfgenexponormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
                if (object$udist == "tslaplace") {
                eval(parse(text = paste0("cmcesftslnormeff_", object$linkF,
                  "(object = object, level = level)")))
                } else {
                if (object$udist == "weibull") {
                  eval(parse(text = paste0("cmcesfweibullnormeff_",
                  object$linkF, "(object = object, level = level)")))
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffRes)
}

# conditional efficiencies sfamisfcross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfamisfcross
efficiencies.sfamisfcross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the MISF model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  if (object$udist == "hnormal") {
    EffRes <- eval(parse(text = paste0("cmisfhalfnormeff_", object$linkF, "(object = object, level = level)")))
  } else {
    if (object$udist == "exponential") {
      eval(parse(text = paste0("cmisfexponormeff_", object$linkF, "(object = object, level = level)")))
    } else {
      if (object$udist == "gamma") {
        eval(parse(text = paste0("cmisfgammanormeff_", object$linkF, "(object = object, level = level)")))
      } else {
        if (object$udist == "rayleigh") {
          eval(parse(text = paste0("cmisfraynormeff_", object$linkF, "(object = object, level = level)")))
        } else {
          if (object$udist == "uniform") {
          eval(parse(text = paste0("cmisfuninormeff_", object$linkF, "(object = object, level = level)")))
          } else {
          if (object$udist == "tnormal") {
            eval(parse(text = paste0("cmisftruncnormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
            if (object$udist == "lognormal") {
            eval(parse(text = paste0("cmisflognormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
            if (object$udist == "genexponential") {
              eval(parse(text = paste0("cmisfgenexponormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "tslaplace") {
              eval(parse(text = paste0("cmisftslnormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "weibull") {
                eval(parse(text = paste0("cmisfweibullnormeff_", object$linkF,
                "(object = object, level = level)")))
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffRes)
}

# conditional efficiencies sfazisfcross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfazisfcross
efficiencies.sfazisfcross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the ZISF model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  if (object$sigmavType == "common") {
    if (object$udist == "hnormal") {
      EffRes <- eval(parse(text = paste0("czisfhalfnormeff_", object$linkF,
        "(object = object, level = level)")))
    } else {
      if (object$udist == "exponential") {
        eval(parse(text = paste0("czisfexponormeff_", object$linkF, "(object = object, level = level)")))
      } else {
        if (object$udist == "gamma") {
          eval(parse(text = paste0("czisfgammanormeff_", object$linkF, "(object = object, level = level)")))
        } else {
          if (object$udist == "rayleigh") {
          eval(parse(text = paste0("czisfraynormeff_", object$linkF, "(object = object, level = level)")))
          } else {
          if (object$udist == "uniform") {
            eval(parse(text = paste0("czisfuninormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
            if (object$udist == "tnormal") {
            eval(parse(text = paste0("czisftruncnormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
            if (object$udist == "lognormal") {
              eval(parse(text = paste0("czisflognormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "genexponential") {
              eval(parse(text = paste0("czisfgenexponormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "tslaplace") {
                eval(parse(text = paste0("czisftslnormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
                if (object$udist == "weibull") {
                eval(parse(text = paste0("czisfweibullnormeff_",
                  object$linkF, "(object = object, level = level)")))
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  } else {
    if (object$sigmavType == "different") {
      if (object$udist == "hnormal") {
        EffRes <- eval(parse(text = paste0("cmnsfhalfnormeff_", object$linkF,
          "(object = object, level = level)")))
      } else {
        if (object$udist == "exponential") {
          eval(parse(text = paste0("cmnsfexponormeff_", object$linkF, "(object = object, level = level)")))
        } else {
          if (object$udist == "gamma") {
          eval(parse(text = paste0("cmnsfgammanormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
          if (object$udist == "rayleigh") {
            eval(parse(text = paste0("cmnsfraynormeff_", object$linkF,
            "(object = object, level = level)")))
          } else {
            if (object$udist == "uniform") {
            eval(parse(text = paste0("cmnsfuninormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
            if (object$udist == "tnormal") {
              eval(parse(text = paste0("cmnsftruncnormeff_", object$linkF,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "lognormal") {
              eval(parse(text = paste0("cmnsflognormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "genexponential") {
                eval(parse(text = paste0("cmnsfgenexponormeff_", object$linkF,
                "(object = object, level = level)")))
              } else {
                if (object$udist == "tslaplace") {
                eval(parse(text = paste0("cmnsftslnormeff_", object$linkF,
                  "(object = object, level = level)")))
                } else {
                if (object$udist == "weibull") {
                  eval(parse(text = paste0("cmnsfweibullnormeff_",
                  object$linkF, "(object = object, level = level)")))
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffRes)
}

# conditional efficiencies sfaselectioncross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfaselectioncross
efficiencies.sfaselectioncross <- function(object, level = 0.95, newData = NULL,
  ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the selection model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  EffRes <- chalfnormeff_ss(object = object, level = level)
  return(EffRes)
}

# conditional efficiencies sfametacross ----------
#' @rdname efficiencies
#' @export
# @exportS3Method efficiencies sfametacross newData option is not used here
efficiencies.sfametacross <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the metafrontier model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  group_var <- object$dataTable[[object$Ngroup + 1]][object$name_meta_var][, 1]
  group_var_list <- sort(unique(group_var))
  MeffList <- list()
  for (g in group_var_list) {
    object_modified <- object
    object_modified$mlParam <- object$mlParam[, which(group_var_list == g)]
    object_modified$dataTable <- object$dataTable[[which(group_var_list == g)]]
    if (object$udist == "hnormal") {
      EffRes <- chalfnormeff(object = object_modified, level = level)
    } else {
      if (object$udist == "exponential") {
        EffRes <- cexponormeff(object = object_modified, level = level)
      } else {
        if (object$udist == "gamma") {
          EffRes <- cgammanormeff(object = object_modified, level = level)
        } else {
          if (object$udist == "rayleigh") {
          EffRes <- craynormeff(object = object_modified, level = level)
          } else {
          if (object$udist == "uniform") {
            EffRes <- cuninormeff(object = object_modified, level = level)
          } else {
            if (object$udist == "tnormal") {
            if (object$scaling) {
              EffRes <- ctruncnormscaleff(object = object_modified, level = level)
            } else {
              EffRes <- ctruncnormeff(object = object_modified, level = level)
            }
            } else {
            if (object$udist == "lognormal") {
              EffRes <- clognormeff(object = object_modified, level = level)
            } else {
              if (object$udist == "genexponential") {
              EffRes <- cgenexponormeff(object = object_modified, level = level)
              } else {
              if (object$udist == "tslaplace") {
                EffRes <- ctslnormeff(object = object_modified, level = level)
              } else {
                if (object$udist == "weibull") {
                EffRes <- cweibullnormeff(object = object_modified,
                  level = level)
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
    rm(object_modified)
    MeffList[[which(group_var_list == g)]] <- EffRes
  }
  MeffDat <- as.data.frame(matrix(nrow = object$Nobs[object$Ngroup + 1], ncol = ncol(MeffList[[1]])))
  for (g in group_var_list) {
    MeffDat[group_var == g, ] <- MeffList[[which(group_var_list == g)]]
  }
  MeffDat <- cbind(group_var, MeffDat)
  names(MeffDat) <- c(object$name_meta_var, paste0(names(MeffList[[1]]), "_g"))
  # metafrontier efficiencies
  if (object$modelType %in% c("hhl14")) {
    object_modified <- object
    object_modified$mlParam <- object$mlParam[, object$Ngroup + 1]
    object_modified$dataTable <- object$dataTable[[object$Ngroup + 1]]
    # replace reponse variable so ineff is properly computed
    object_modified$dataTable[all.vars(object$formula)[1]] <- object$Yvarm
    if (object$udist == "hnormal") {
      EffRes <- chalfnormeff(object = object_modified, level = level)
    } else {
      if (object$udist == "exponential") {
        EffRes <- cexponormeff(object = object_modified, level = level)
      } else {
        if (object$udist == "gamma") {
          EffRes <- cgammanormeff(object = object_modified, level = level)
        } else {
          if (object$udist == "rayleigh") {
          EffRes <- craynormeff(object = object_modified, level = level)
          } else {
          if (object$udist == "uniform") {
            EffRes <- cuninormeff(object = object_modified, level = level)
          } else {
            if (object$udist == "tnormal") {
            if (object$scaling) {
              EffRes <- ctruncnormscaleff(object = object_modified, level = level)
            } else {
              EffRes <- ctruncnormeff(object = object_modified, level = level)
            }
            } else {
            if (object$udist == "lognormal") {
              EffRes <- clognormeff(object = object_modified, level = level)
            } else {
              if (object$udist == "genexponential") {
              EffRes <- cgenexponormeff(object = object_modified, level = level)
              } else {
              if (object$udist == "tslaplace") {
                EffRes <- ctslnormeff(object = object_modified, level = level)
              } else {
                if (object$udist == "weibull") {
                EffRes <- cweibullnormeff(object = object_modified,
                  level = level)
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
    rm(object_modified)
    if (object$logDepVar == TRUE) {
      MeffDat$TGR_teJLMS <- EffRes$teJLMS
      MeffDat$TGR_teBC <- EffRes$teBC
      MeffDat$MTE_teJLMS <- MeffDat$teJLMS_g * MeffDat$TGR_teJLMS
      MeffDat$MTE_teBC <- MeffDat$teBC_g * MeffDat$TGR_teBC
    } else {
      MeffDat$MD <- EffRes$u
      MeffDat$MTI <- MeffDat$u_g + MeffDat$MD
    }
  } else {
    if (object$modelType %in% c("bpo04a", "bpo04b")) {
      if (object$logDepVar == TRUE) {
        MeffDat$TGR <- exp(object$dataTable[[object$Ngroup + 1]]$mlResiduals)
        MeffDat$MTE_teJLMS <- MeffDat$teJLMS_g * MeffDat$TGR
        MeffDat$MTE_teBC <- MeffDat$teBC_g * MeffDat$TGR
      } else {
        MeffDat$MD <- abs(object$dataTable[[object$Ngroup + 1]]$mlResiduals)
        MeffDat$MTI <- MeffDat$u_g + MeffDat$MD
      }
    } else {
      if (object$modelType %in% c("aos17a", "aos17c")) {
        MeffDat$Mud <- object$Mud
        if (object$logDepVar == TRUE) {
          MeffDat$TGR <- exp(-apply(object$MdMat, 1, mean))
          MeffDat$MTE <- exp(-MeffDat$Mud) * MeffDat$TGR
        } else {
          MeffDat$MD <- apply(object$MdMat, 1, mean)
          MeffDat$MTI <- MeffDat$Mud + MeffDat$MD
        }
      } else {
        if (object$modelType %in% c("aos17b", "aos17d")) {
          if (object$logDepVar == TRUE) {
          MeffDat$TGR <- exp(-apply(object$MdMat, 1, mean))
          MeffDat$MTE_teJLMS <- MeffDat$teJLMS_g * MeffDat$TGR
          MeffDat$MTE_teBC <- MeffDat$teBC_g * MeffDat$TGR
          } else {
          MeffDat$MD <- apply(object$MdMat, 1, mean)
          MeffDat$MTI <- MeffDat$u_g + MeffDat$MD
          }
        }
      }
    }
  }
  return(MeffDat)
}

# conditional efficiencies sfapanel1 ----------
#' @rdname efficiencies
#' @aliases efficiencies.sfapanel1
#' @export
efficiencies.sfapanel1 <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the panel model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  toPaste <- if (object$modelType == "pl81") {
    "pl81"
  } else {
    if (object$modelType == "mbc92") {
      "mbc92"
    } else {
      if (object$modelType == "k90") {
        "k90"
      } else {
        if (object$modelType == "mols93") {
          "mols93"
        } else {
          if (object$modelType %in% c("bc92a", "bc92b", "bc92c", "kw05",
          "c00")) {
          "gzit"
          }
        }
      }
    }
  }
  if (object$udist == "hnormal") {
    EffRes <- eval(parse(text = paste0("phalfnormeff_", toPaste, "(object = object, level = level)")))
  } else {
    if (object$udist == "exponential") {
      EffRes <- eval(parse(text = paste0("pexponormeff_", toPaste, "(object = object, level = level)")))
    } else {
      if (object$udist == "gamma") {
        EffRes <- eval(parse(text = paste0("pgammanormeff_", toPaste, "(object = object, level = level)")))
      } else {
        if (object$udist == "rayleigh") {
          EffRes <- eval(parse(text = paste0("praynormeff_", toPaste, "(object = object, level = level)")))
        } else {
          if (object$udist == "uniform") {
          EffRes <- eval(parse(text = paste0("puninormeff_", toPaste, "(object = object, level = level)")))
          } else {
          if (object$udist == "tnormal") {
            EffRes <- eval(parse(text = paste0("ptruncnormeff_", toPaste, "(object = object, level = level)")))
          } else {
            if (object$udist == "lognormal") {
            EffRes <- eval(parse(text = paste0("plognormeff_", toPaste, "(object = object, level = level)")))
            } else {
            if (object$udist == "genexponential") {
              EffRes <- eval(parse(text = paste0("pgenexponormeff_", toPaste,
              "(object = object, level = level)")))
            } else {
              if (object$udist == "tslaplace") {
              EffRes <- eval(parse(text = paste0("ptslnormeff_", toPaste,
                "(object = object, level = level)")))
              } else {
              if (object$udist == "weibull") {
                EffRes <- eval(parse(text = paste0("pweibullnormeff_", toPaste,
                "(object = object, level = level)")))
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(EffRes)
}

# conditional efficiencies sfalcmpanel ----------
#' @rdname efficiencies
#' @aliases efficiencies.sfalcmpanel
#' @export
efficiencies.sfalcmpanel <- function(object, level = 0.95, newData = NULL, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame", call. = FALSE)
    }
    NobsOldData <- nrow(object$dataTable)
    if (nrow(newData) != NobsOldData) {
      stop("New Data must have the same structure as old data for the LCM model",
        call. = FALSE)
    }
    object$dataTable <- newData
  }
  toPaste <- if (object$modelType == "pl81") {
    "pl81"
  } else {
    if (object$modelType == "mbc92") {
      "mbc92"
    } else {
      if (object$modelType == "k90") {
        "k90"
      } else {
        if (object$modelType == "mols93") {
          "mols93"
        } else {
          if (object$modelType %in% c("bc92a", "bc92b", "bc92c", "kw05",
          "c00")) {
          "gzit"
          }
        }
      }
    }
  }
  EffRes <- eval(parse(text = paste0("pLCM", object$nClasses, "Chalfnormeff_", toPaste,
    "(object = object, level = level)")))
  return(EffRes)
}
