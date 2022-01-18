################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Efficiency/Inefficiency estimation                                           #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Compute conditional (in-)efficiency estimates of stochastic frontier models
#'
#' \code{\link{efficiencies}} returns (in-)efficiency estimates of models estimated with
#' \code{\link{sfacross}}, \code{\link{lcmcross}} or \code{\link{selectioncross}}.
#' 
#' @name efficiencies
#'
#' @details The conditional inefficiency is obtained following Jondrow \emph{et al.}
#' (1982) and the conditional efficiency is computed following Battese and
#' Coelli (1988). In some cases the conditional mode is also returned (Jondrow
#' \emph{et al.} 1982). The confidence interval is computed following Horrace
#' and Schmidt (1996), Hjalmarsson \emph{et al.} (1996), or Berra and Sharma
#' (1999) (see \sQuote{Value} section).
#'
#' In the case of the half normal distribution for the one-sided error term,
#' the formulae are as follows (for notations, see the \sQuote{Details} section
#' of \code{\link{sfacross}}, \code{\link{lcmcross}} or \code{\link{selectioncross}}):
#'
#' \itemize{ \item The conditional inefficiency is }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("E\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right\\\\rbrack=\\\mu_{i\\\ast} + \\\sigma_\\\ast\\\frac{\\\phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}{\\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}") 
#' }
#'
#' where
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("\\\mu_{i\\\ast}=\\\frac{-S\\\epsilon_i\\\sigma_u^2}{\\\sigma_u^2 + \\\sigma_v^2}") 
#' }
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("\\\sigma_\\\ast^2 = \\\frac{\\\sigma_u^2 \\\sigma_v^2}{\\\sigma_u^2 + \\\sigma_v^2}") 
#' }
#' 
#' \itemize{ \item The Battese and Coelli (1988) conditional efficiency is
#' obtained by: } 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("E\\\\left\\\\lbrack\\\exp{\\\\left(-u_i\\\\right)}|\\\\epsilon_i\\\\right\\\\rbrack = \\\exp{\\\\left(-\\\mu_{i\\\ast}+\\\frac{1}{2}\\\sigma_\\\ast^2\\\\right)}\\\frac{\\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}-\\\sigma_\\\ast\\\\right)}{\\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\right)}") 
#' }
#' 
#' \itemize{ \item The reciprocal of the Battese and Coelli (1988) conditional efficiency is
#' obtained by: } 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("E\\\\left\\\\lbrack\\\exp{\\\\left(u_i\\\\right)}|\\\epsilon_i\\\\right\\\\rbrack = \\\exp{\\\\left(\\\mu_{i\\\ast}+\\\frac{1}{2}\\\sigma_\\\ast^2\\\\right)} \\\frac{\\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}+\\\sigma_\\\ast\\\\right)}{\\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}") 
#' }
#' 
#' \itemize{ \item The conditional mode is computed using: }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("M\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right\\\\rbrack= \\\mu_{i\\\ast} \\\quad \\\hbox{For} \\\quad \\\mu_{i\\\ast} > 0") 
#' }
#' 
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("M\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right\\\\rbrack= 0 \\\quad \\\hbox{For} \\\quad \\\mu_{i\\\ast} \\\\leq 0") 
#' } 
#' 
#' \itemize{ \item The confidence intervals are obtained with: }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("\\\mu_{i\\\ast} + I_L\\\sigma_\\\ast \\\\leq E\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right\\\\rbrack \\\\leq \\\mu_{i\\\ast} + I_U\\\sigma_\\\ast") 
#' }
#'
#' with \eqn{LB_i = \mu_{i*} + I_L\sigma_*} and \eqn{UB_i = \mu_{i*} + I_U\sigma_*}
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("I_L = \\\Phi^{-1}\\\\left\\\\lbrace 1 -\\\\left(1-\\\frac{\\\alpha}{2}\\\\right)\\\\left\\\\lbrack 1-\\\Phi\\\\left(-\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)\\\\right\\\\rbrack\\\\right\\\\rbrace") 
#' }
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("I_U = \\\Phi^{-1}\\\\left\\\\lbrace 1-\\\frac{\\\alpha}{2}\\\\left\\\\lbrack 1-\\\Phi\\\\left(-\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)\\\\right\\\\rbrack\\\\right\\\\rbrace") 
#' }
#'
#' Thus
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd("\\\exp{\\\\left(-UB_i\\\\right)} \\\\leq E\\\\left\\\\lbrack\\\exp{\\\\left(-u_i\\\\right)}|\\\epsilon_i\\\\right\\\\rbrack \\\\leq\\\exp{\\\\left(-LB_i\\\\right)}") 
#' }
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{lcmcross}} or \code{\link{selectioncross}}.
#' @param level A number between between 0 and 0.9999 used for the computation
#' of (in-)efficiency confidence intervals (defaut = \code{0.95}). Only used
#' when \code{udist} = \code{'hnormal'}, \code{'exponential'}, \code{'tnormal'}
#' or \code{'uniform'} in \code{\link{sfacross}}. The option is also available for
#' \code{\link{lcmcross}}, or \code{\link{selectioncross}}.
#' @param newData Optional data frame that is used to calculate the efficiency 
#' estimates. If NULL (the default), the efficiency estimates are calculated 
#' for the observations that were used in the estimation.
#' @param ... Currently ignored.
#'
#' @return A data frame that contains individual (in-)efficiency estimates.
#' These are ordered in the same way as the corresponding observations in the
#' dataset used for the estimation.
#'
#' \bold{- For object of class \code{'sfacross'} the following elements are
#' returned:}
#'
#' \item{u}{Conditional inefficiency. In the case argument \code{udist} of
#' \link{sfacross} is set to \code{'uniform'}, two conditional inefficiency
#' estimates are returned: \code{u1} for the classic conditional inefficiency
#' following Jondrow \emph{et al.} (1982), and \code{u2} which is obtained when
#' \eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen, 2010).}
#'
#' \item{uLB}{Lower bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross} is set to \code{'hnormal'},
#' \code{'exponential'}, \code{'tnormal'} or \code{'uniform'}.}
#'
#' \item{uUB}{Upper bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross} is set to \code{'hnormal'},
#' \code{'exponential'}, \code{'tnormal'} or \code{'uniform'}.}
#'
#' \item{teJLMS}{\eqn{\exp{(-E[u|\epsilon])}}. When the argument \code{udist} of
#' \link{sfacross} is set to \code{'uniform'}, \code{teJLMS1} =
#' \eqn{\exp{(-E[u_1|\epsilon])}} and \code{teJLMS2} = \eqn{\exp{(-E[u_2|\epsilon])}}. Only when
#' \code{logDepVar = TRUE}.}
#'
#' \item{m}{Conditional model. Only when the argument \code{udist} of
#' \link{sfacross} is set to \code{'hnormal'}, \code{'exponential'},
#' \code{'tnormal'}, or \code{'rayleigh'}.}
#'
#' \item{teMO}{\eqn{\exp{(-m)}}. Only when, in the function \link{sfacross},
#' \code{logDepVar = TRUE} and \code{udist = 'hnormal'}, \code{'exponential'},
#' \code{'tnormal'}, \code{'uniform'}, or \code{'rayleigh'}.}
#'
#' \item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
#' the function \link{sfacross}, \code{logDepVar = TRUE}. In the case \code{udist =
#' 'uniform'}, two conditional efficiency estimates are returned: \code{teBC1}
#' which is the classic conditional efficiency following Battese and Coelli
#' (1988) and \code{teBC2} when \eqn{\theta/\sigma_v \longrightarrow \infty}
#' (see Nguyen, 2010).}
#' 
#' \item{teBC_reciprocal}{Reciprocal of Battese and Coelli (1988) conditional efficiency. 
#' Similar to \code{teBC} except that it is computed as \eqn{E\left[\exp{(u)}|\epsilon\right]}.}
#'
#' \item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
#' TRUE} and \code{udist = 'hnormal'}, \code{'exponential'}, \code{'tnormal'},
#' or \code{'uniform'}.}
#'
#' \item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, \code{logDepVar =
#' TRUE} and \code{udist = 'hnormal'}, \code{'exponential'}, \code{'tnormal'},
#' or \code{'uniform'}.}
#'
#' \bold{- For object of class \code{'lcmcross'} the following elements are
#' returned:}
#'
#' \item{Group_c}{Most probable class for each observation.}
#'
#' \item{PosteriorProb_c}{Highest posterior probability.}
#' 
#' \item{u_c}{Conditional inefficiency of the most probable class given the
#' posterior probability.}
#' 
#' \item{teJLMS_c}{\eqn{\exp{(-E[u_c|\epsilon_c])}}. Only when, in the function
#' \link{lcmcross}, \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_c}{\eqn{E\left[\exp{(-u_c)}|\epsilon_c\right]}. Only when, in the function
#' \link{lcmcross}, \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c}{\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}. Only when, in the function
#' \link{lcmcross}, \code{logDepVar = TRUE}.}
#'
#' \item{PosteriorProb_c#}{Posterior probability of class #.}
#'
#' \item{PriorProb_c#}{Prior probability of class #.}
#'
#' \item{u_c#}{Conditional inefficiency associated to class #, regardless of
#' \code{Group_c}.}
#' 
#' \item{teBC_c#}{Conditional efficiency (\eqn{E\left[\exp{(-u_c)}|\epsilon_c\right]}) associated to class #, regardless of
#' \code{Group_c}. Only when, in the function \link{lcmcross}, \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c#}{Reciprocal conditional efficiency (\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}) 
#' associated to class #, regardless of \code{Group_c}. 
#' Only when, in the function \link{lcmcross}, \code{logDepVar = TRUE}.}
#'
#' \item{ineff_c#}{Conditional inefficiency (\code{u_c}) for observations in
#' class # only.}
#' 
#' \item{effBC_c#}{Conditional efficiency (\code{teBC_c}) for observations in
#' class # only.}
#' 
#' \item{ReffBC_c#}{Reciprocal conditional efficiency (\code{teBC_reciprocal_c}) for observations in
#' class # only.}
#' 
#' \bold{- For object of class \code{'selectioncross'} the following elements are
#' returned:}
#'
#' \item{u}{Conditional inefficiency.} 
#' \item{uLB}{Lower bound for conditional inefficiency.} 
#' \item{uUB}{Upper bound for conditional inefficiency.} 
#' \item{teJLMS}{\eqn{\exp{(-u)}}. Only when \code{logDepVar = TRUE}.}
#' \item{m}{Conditional mode.} 
#' \item{teMO}{\eqn{\exp{(-m)}}. Only when \code{logDepVar = TRUE}.}
#' \item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when 
#' \code{logDepVar = TRUE}.}
#' \item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
#' efficiency. Only when \code{logDepVar = TRUE}.}
#' \item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
#' efficiency. Only when \code{logDepVar = TRUE}.}
#'
#' @author K Hervé Dakpo, Yann Desjeux, and Laure Latruffe
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{selectioncross}} for sample selection in stochastic frontier model
#' fitting function.
#'
#' @references Battese, G.E., and T.J. Coelli. 1988. Prediction of firm-level
#' technical efficiencies with a generalized frontier production function and
#' panel data. \emph{Journal of Econometrics}, \bold{38}:387--399.
#'
#' Bera, A.K., and S.C. Sharma. 1999. Estimating production uncertainty in
#' stochastic frontier production function models. \emph{Journal of
#' Productivity Analysis}, \bold{12}:187-210.
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
#' Nguyen, N.B. 2010. Estimation of technical efficiency in stochastic frontier
#' analysis. PhD Dissertation, Bowling Green State University, August.
#'
#' @keywords methods efficiencies
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) + log(wl/wf) +
#'     log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) + I(log(wl/wf) *
#'     log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)), udist = 'tnormal',
#'     muhet = ~ regu, uhet = ~ regu, data = utility, S = -1, scaling = TRUE, method = 'mla')
#'   eff.tl_u_ts <- efficiencies(tl_u_ts)
#'   head(eff.tl_u_ts)
#'   summary(eff.tl_u_ts)
#'
#' @aliases efficiencies.sfacross
#' @export
#' @export efficiencies
# conditional efficiencies sfacross ----------
efficiencies.sfacross <- function(object, level = 0.95, newData = NULL,
  ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
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
            EffRes <- ctruncnormscaleff(object = object,
              level = level)
            } else {
            EffRes <- ctruncnormeff(object = object,
              level = level)
            }
          } else {
            if (object$udist == "lognormal") {
            EffRes <- clognormeff(object = object,
              level = level)
            } else {
            if (object$udist == "genexponential") {
              EffRes <- cgenexponormeff(object = object,
              level = level)
            } else {
              if (object$udist == "tslaplace") {
              EffRes <- ctslnormeff(object = object,
                level = level)
              } else {
              if (object$udist == "weibull") {
                EffRes <- cweibullnormeff(object = object,
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
  return(data.frame(EffRes))
}

# conditional efficiencies lcmcross ----------
#' @rdname efficiencies
#' @aliases efficiencies.lcmcross
#' @export
efficiencies.lcmcross <- function(object, level = 0.95, newData = NULL,
  ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$nClasses == 2) {
    EffRes <- cLCM2Chalfnormeff(object = object, level = level)
  } else {
    if (object$nClasses == 3) {
      EffRes <- cLCM3Chalfnormeff(object = object, level = level)
    } else {
      if (object$nClasses == 4) {
        EffRes <- cLCM4Chalfnormeff(object = object,
          level = level)
      } else {
        if (object$nClasses == 5) {
          EffRes <- cLCM5Chalfnormeff(object = object,
          level = level)
        }
      }
    }
  }
  return(data.frame(EffRes))
}

# conditional efficiencies selectioncross ----------
#' @rdname efficiencies
#' @aliases efficiencies.selectioncross
#' @export
efficiencies.selectioncross <- function(object, level = 0.95, newData = NULL,
                                        ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  EffRes <- chalfnormeff_ss(object = object, level = level)
  return(data.frame(EffRes))
}

