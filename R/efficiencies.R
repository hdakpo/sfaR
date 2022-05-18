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
#         -Zero inefficiency stochastic frontier                               #
#         -Contaminated noise stochastic frontier                              #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Compute conditional (in-)efficiency estimates of stochastic frontier models
#'
#' \code{\link{efficiencies}} returns (in-)efficiency estimates of models 
#' estimated with \code{\link{cnsfcross}}, \code{\link{lcmcross}}, 
#' \code{\link{sfacross}}, \code{\link{sfaselectioncross}}, or 
#' \code{\link{zisfcross}}.
#' 
#' @name efficiencies
#'
#' @details The conditional inefficiency is obtained following 
#' Jondrow \emph{et al.} (1982) and the conditional efficiency is computed 
#' following Battese and Coelli (1988). In some cases the conditional mode is 
#' also returned (Jondrow \emph{et al.} 1982). The confidence interval is 
#' computed following Horrace and Schmidt (1996), Hjalmarsson \emph{et al.} 
#' (1996), or Berra and Sharma (1999) (see \sQuote{Value} section).
#'
#' In the case of the half normal distribution for the one-sided error term,
#' the formulae are as follows (for notations, see the \sQuote{Details} section
#' of \code{\link{sfacross}}, \code{\link{lcmcross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}):
#'
#' \itemize{ \item The conditional inefficiency is }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('E\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right
#' \\\\rbrack=\\\mu_{i\\\ast} + \\\sigma_\\\ast\\\frac{\\\phi
#' \\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}{
#' \\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}') 
#' }
#'
#' where
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\mu_{i\\\ast}=\\\frac{-S\\\epsilon_i\\\sigma_u^2}{
#' \\\sigma_u^2 + \\\sigma_v^2}') 
#' }
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\sigma_\\\ast^2 = \\\frac{\\\sigma_u^2 \\\sigma_v^2}{
#' \\\sigma_u^2 + \\\sigma_v^2}') 
#' }
#' 
#' \itemize{ \item The Battese and Coelli (1988) conditional efficiency is
#' obtained by: } 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('E\\\\left\\\\lbrack\\\exp{\\\\left(-u_i\\\\right)}
#' |\\\\epsilon_i\\\\right\\\\rbrack = \\\exp{\\\\left(-\\\mu_{i\\\ast}+
#' \\\frac{1}{2}\\\sigma_\\\ast^2\\\\right)}\\\frac{\\\Phi\\\\left(
#' \\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}-\\\sigma_\\\ast\\\\right)}{
#' \\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\right)}') 
#' }
#' 
#' \itemize{ \item The reciprocal of the Battese and Coelli (1988) conditional 
#' efficiency is obtained by: } 
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('E\\\\left\\\\lbrack\\\exp{\\\\left(u_i\\\\right)}
#' |\\\epsilon_i\\\\right\\\\rbrack = \\\exp{\\\\left(\\\mu_{i\\\ast}+
#' \\\frac{1}{2}\\\sigma_\\\ast^2\\\\right)} \\\frac{\\\Phi\\\\left(
#' \\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}+\\\sigma_\\\ast\\\\right)}{
#' \\\Phi\\\\left(\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)}') 
#' }
#' 
#' \itemize{ \item The conditional mode is computed using: }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('M\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right
#' \\\\rbrack= \\\mu_{i\\\ast} \\\quad \\\hbox{For} \\\quad 
#' \\\mu_{i\\\ast} > 0') 
#' }
#' 
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('M\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right
#' \\\\rbrack= 0 \\\quad \\\hbox{For} \\\quad \\\mu_{i\\\ast} \\\\leq 0') 
#' } 
#' 
#' \itemize{ \item The confidence intervals are obtained with: }
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\mu_{i\\\ast} + I_L\\\sigma_\\\ast \\\\leq 
#' E\\\\left\\\\lbrack u_i|\\\epsilon_i\\\\right\\\\rbrack \\\\leq 
#' \\\mu_{i\\\ast} + I_U\\\sigma_\\\ast') 
#' }
#'
#' with \eqn{LB_i = \mu_{i*} + I_L\sigma_*} and 
#' \eqn{UB_i = \mu_{i*} + I_U\sigma_*}
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('I_L = \\\Phi^{-1}\\\\left\\\\lbrace 1 -
#' \\\\left(1-\\\frac{\\\alpha}{2}\\\\right)\\\\left\\\\lbrack 1-
#' \\\Phi\\\\left(-\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)
#' \\\\right\\\\rbrack\\\\right\\\\rbrace') 
#' }
#'
#' and
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('I_U = \\\Phi^{-1}\\\\left\\\\lbrace 1-
#' \\\frac{\\\alpha}{2}\\\\left\\\\lbrack 1-\\\Phi
#' \\\\left(-\\\frac{\\\mu_{i\\\ast}}{\\\sigma_\\\ast}\\\\right)
#' \\\\right\\\\rbrack\\\\right\\\\rbrace') 
#' }
#'
#' Thus
#' 
#' \Sexpr[results=rd, stage=build]{
#' katex::math_to_rd('\\\exp{\\\\left(-UB_i\\\\right)} \\\\leq E\\\\left
#' \\\\lbrack\\\exp{\\\\left(-u_i\\\\right)}|\\\epsilon_i\\\\right\\\\rbrack 
#' \\\\leq\\\exp{\\\\left(-LB_i\\\\right)}') 
#' }
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{cnsfcross}}, \code{\link{lcmcross}},
#'  \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}.
#' @param level A number between between 0 and 0.9999 used for the computation
#' of (in-)efficiency confidence intervals (defaut = \code{0.95}). Only used
#' when \code{udist} = \code{'hnormal'}, \code{'exponential'}, \code{'tnormal'}
#' or \code{'uniform'} in \code{\link{sfacross}}. The option is also available 
#' for \code{\link{sfaselectioncross}}.
#' @param newData Optional data frame that is used to calculate the efficiency 
#' estimates. If NULL (the default), the efficiency estimates are calculated 
#' for the observations that were used in the estimation.
#' @param ... Currently ignored.
#'
#' @return A data frame that contains individual (in-)efficiency estimates.
#' These are ordered in the same way as the corresponding observations in the
#' dataset used for the estimation.
#' 
#' \bold{- For object of class \code{'cnsfcross'}, or \code{'lcmcross'}, or 
#' \code{'zisfcross'} the following elements are
#' returned:}
#'
#' \item{Group_c}{Most probable class for each observation.}
#'
#' \item{PosteriorProb_c}{Highest posterior probability.}
#' 
#' \item{odRatio}{Posterior odds ratio 
#' \eqn{R_i = \frac{Post. Prob. Class 2}{Post. Prob. Class 1}}. Only for object 
#' of class \code{'zisfcross'}}.
#' 
#' \item{u_c}{Conditional inefficiency of the most probable class given the
#' posterior probability.}
#' 
#' \item{teJLMS_c}{\eqn{\exp{(-E[u_c|\epsilon_c])}}. Only when, in the function
#' \link{cnsfcross}, or \link{lcmcross}, or \link{zisfcross}, 
#' \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_c}{\eqn{E\left[\exp{(-u_c)}|\epsilon_c\right]}. Only when, in the
#'  function \link{cnsfcross}, or \link{lcmcross}, or \link{zisfcross}, 
#'  \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c}{\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}. Only 
#' when, in the function \link{cnsfcross}, or \link{lcmcross}, or 
#' \link{zisfcross}, \code{logDepVar = TRUE}.}
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
#' regardless of \code{Group_c}. Only when, in the function \link{cnsfcross}, 
#' or \link{lcmcross}, or \link{zisfcross}, \code{logDepVar = TRUE}.}
#' 
#' \item{teBC_reciprocal_c#}{Reciprocal conditional efficiency 
#' (\eqn{E\left[\exp{(u_c)}|\epsilon_c\right]}) associated to class #, 
#' regardless of \code{Group_c}. Only when, in the function \link{cnsfcross}, 
#' or \link{lcmcross}, or \link{zisfcross}, \code{logDepVar = TRUE}.}
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
#' \bold{- For object of class \code{'sfacross'} or \code{'sfaselectioncross'} 
#' the following elements are returned:}
#'
#' \item{u}{Conditional inefficiency. In the case argument \code{udist} of
#' \link{sfacross} is set to \code{'uniform'}, two conditional inefficiency
#' estimates are returned: \code{u1} for the classic conditional inefficiency
#' following Jondrow \emph{et al.} (1982), and \code{u2} which is obtained when
#' \eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen, 2010).}
#'
#' \item{uLB}{Lower bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross} is set to \code{'hnormal'},
#' \code{'exponential'}, \code{'tnormal'} or \code{'uniform'}. For object of 
#' class \code{'sfaselectioncross'}, only the \code{'hnormal'} distribution is 
#' available.}
#'
#' \item{uUB}{Upper bound for conditional inefficiency. Only when the argument
#' \code{udist} of \link{sfacross} is set to \code{'hnormal'},
#' \code{'exponential'}, \code{'tnormal'} or \code{'uniform'}. For object of 
#' class \code{'sfaselectioncross'}, only the \code{'hnormal'} distribution is 
#' available.}
#'
#' \item{teJLMS}{\eqn{\exp{(-E[u|\epsilon])}}. When the argument \code{udist} of
#' \link{sfacross} is set to \code{'uniform'}, \code{teJLMS1} =
#' \eqn{\exp{(-E[u_1|\epsilon])}} and \code{teJLMS2} = 
#' \eqn{\exp{(-E[u_2|\epsilon])}}. Only when \code{logDepVar = TRUE}.}
#'
#' \item{m}{Conditional model. Only when the argument \code{udist} of
#' \link{sfacross} is set to \code{'hnormal'}, \code{'exponential'},
#' \code{'tnormal'}, or \code{'rayleigh'}. For object of class 
#' \code{'sfaselectioncross'}, only the \code{'hnormal'} distribution is 
#' available.}
#'
#' \item{teMO}{\eqn{\exp{(-m)}}. Only when, in the function \link{sfacross},
#' \code{logDepVar = TRUE} and \code{udist = 'hnormal'}, \code{'exponential'},
#' \code{'tnormal'}, \code{'uniform'}, or \code{'rayleigh'}. For object of 
#' class \code{'sfaselectioncross'}, only the \code{'hnormal'} distribution is 
#' available.}
#'
#' \item{teBC}{Battese and Coelli (1988) conditional efficiency. Only when, in
#' the function \link{sfacross} or \code{'sfaselectioncross'}, 
#' \code{logDepVar = TRUE}. In the case \code{udist = 'uniform'}, two 
#' conditional efficiency estimates are returned:
#' \code{teBC1} which is the classic conditional efficiency following 
#' Battese and Coelli (1988) and \code{teBC2} when 
#' \eqn{\theta/\sigma_v \longrightarrow \infty} (see Nguyen, 2010).}
#' 
#' \item{teBC_reciprocal}{Reciprocal of Battese and Coelli (1988) conditional 
#' efficiency. Similar to \code{teBC} except that it is computed as 
#' \eqn{E\left[\exp{(u)}|\epsilon\right]}.}
#'
#' \item{teBCLB}{Lower bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, or 
#' \code{'sfaselectioncross'}, \code{logDepVar = TRUE} and 
#' \code{udist = 'hnormal'}, \code{'exponential'}, \code{'tnormal'},
#' or \code{'uniform'}. For object of class \code{'sfaselectioncross'}, only 
#' the \code{'hnormal'} distribution is available.}
#'
#' \item{teBCUB}{Upper bound for Battese and Coelli (1988) conditional
#' efficiency. Only when, in the function \link{sfacross}, or 
#' \code{'sfaselectioncross'}, \code{logDepVar = TRUE} and 
#' \code{udist = 'hnormal'}, \code{'exponential'}, \code{'tnormal'},
#' or \code{'uniform'}. For object of class \code{'sfaselectioncross'}, only 
#' the \code{'hnormal'} distribution is available.}
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{cnsfcross}}, for the contaminated noise stochastic 
#' frontier analysis model fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function.
#' 
#' \code{\link{zisfcross}} for zero inefficiency in stochastic frontier model
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
#' log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) + I(log(wl/wf) *
#' log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)), udist = 'tnormal',
#' muhet = ~ regu, uhet = ~ regu, data = utility, S = -1, scaling = TRUE, method = 'mla')
#' eff.tl_u_ts <- efficiencies(tl_u_ts)
#' head(eff.tl_u_ts)
#' summary(eff.tl_u_ts)
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

# conditional efficiencies sfaselectioncross ----------
#' @rdname efficiencies
#' @aliases efficiencies.sfaselectioncross
#' @export
efficiencies.sfaselectioncross <- function(object, level = 0.95,
  newData = NULL, ...) {
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

# conditional efficiencies zisfcross ----------
#' @rdname efficiencies
#' @aliases efficiencies.zisfcross
#' @export
efficiencies.zisfcross <- function(object, level = 0.95, newData = NULL,
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
  if (object$linkF == "logit") {
    if (object$sigmavType == "common") {
      if (object$udist == "hnormal") {
        EffRes <- czisfhalfnormeff_logit(object = object,
          level = level)
      } else {
        if (object$udist == "exponential") {
          EffRes <- czisfexponormeff_logit(object = object,
          level = level)
        } else {
          if (object$udist == "gamma") {
          EffRes <- czisfgammanormeff_logit(object = object,
            level = level)
          } else {
          if (object$udist == "rayleigh") {
            EffRes <- czisfraynormeff_logit(object = object,
            level = level)
          } else {
            if (object$udist == "uniform") {
            EffRes <- czisfuninormeff_logit(object = object,
              level = level)
            } else {
            if (object$udist == "tnormal") {
              EffRes <- czisftruncnormeff_logit(object = object,
              level = level)
            } else {
              if (object$udist == "lognormal") {
              EffRes <- czisflognormeff_logit(object = object,
                level = level)
              } else {
              if (object$udist == "genexponential") {
                EffRes <- czisfgenexponormeff_logit(object = object,
                level = level)
              } else {
                if (object$udist == "tslaplace") {
                EffRes <- czisftslnormeff_logit(object = object,
                  level = level)
                } else {
                if (object$udist == "weibull") {
                  EffRes <- czisfweibullnormeff_logit(object = object,
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
    } else {
      if (object$sigmavType == "different") {
        if (object$udist == "hnormal") {
          EffRes <- cmnsfhalfnormeff_logit(object = object,
          level = level)
        } else {
          if (object$udist == "exponential") {
          EffRes <- cmnsfexponormeff_logit(object = object,
            level = level)
          } else {
          if (object$udist == "gamma") {
            EffRes <- cmnsfgammanormeff_logit(object = object,
            level = level)
          } else {
            if (object$udist == "rayleigh") {
            EffRes <- cmnsfraynormeff_logit(object = object,
              level = level)
            } else {
            if (object$udist == "uniform") {
              EffRes <- cmnsfuninormeff_logit(object = object,
              level = level)
            } else {
              if (object$udist == "tnormal") {
              EffRes <- cmnsftruncnormeff_logit(object = object,
                level = level)
              } else {
              if (object$udist == "lognormal") {
                EffRes <- cmnsflognormeff_logit(object = object,
                level = level)
              } else {
                if (object$udist == "genexponential") {
                EffRes <- cmnsfgenexponormeff_logit(object = object,
                  level = level)
                } else {
                if (object$udist == "tslaplace") {
                  EffRes <- cmnsftslnormeff_logit(object = object,
                  level = level)
                } else {
                  if (object$udist == "weibull") {
                  EffRes <- cmnsfweibullnormeff_logit(object = object,
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
      }
    }
  } else {
    if (object$linkF == "cauchit") {
      if (object$sigmavType == "common") {
        if (object$udist == "hnormal") {
          EffRes <- czisfhalfnormeff_cauchit(object = object,
          level = level)
        } else {
          if (object$udist == "exponential") {
          EffRes <- czisfexponormeff_cauchit(object = object,
            level = level)
          } else {
          if (object$udist == "gamma") {
            EffRes <- czisfgammanormeff_cauchit(object = object,
            level = level)
          } else {
            if (object$udist == "rayleigh") {
            EffRes <- czisfraynormeff_cauchit(object = object,
              level = level)
            } else {
            if (object$udist == "uniform") {
              EffRes <- czisfuninormeff_cauchit(object = object,
              level = level)
            } else {
              if (object$udist == "tnormal") {
              EffRes <- czisftruncnormeff_cauchit(object = object,
                level = level)
              } else {
              if (object$udist == "lognormal") {
                EffRes <- czisflognormeff_cauchit(object = object,
                level = level)
              } else {
                if (object$udist == "genexponential") {
                EffRes <- czisfgenexponormeff_cauchit(object = object,
                  level = level)
                } else {
                if (object$udist == "tslaplace") {
                  EffRes <- czisftslnormeff_cauchit(object = object,
                  level = level)
                } else {
                  if (object$udist == "weibull") {
                  EffRes <- czisfweibullnormeff_cauchit(object = object,
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
      } else {
        if (object$sigmavType == "different") {
          if (object$udist == "hnormal") {
          EffRes <- cmnsfhalfnormeff_cauchit(object = object,
            level = level)
          } else {
          if (object$udist == "exponential") {
            EffRes <- cmnsfexponormeff_cauchit(object = object,
            level = level)
          } else {
            if (object$udist == "gamma") {
            EffRes <- cmnsfgammanormeff_cauchit(object = object,
              level = level)
            } else {
            if (object$udist == "rayleigh") {
              EffRes <- cmnsfraynormeff_cauchit(object = object,
              level = level)
            } else {
              if (object$udist == "uniform") {
              EffRes <- cmnsfuninormeff_cauchit(object = object,
                level = level)
              } else {
              if (object$udist == "tnormal") {
                EffRes <- cmnsftruncnormeff_cauchit(object = object,
                level = level)
              } else {
                if (object$udist == "lognormal") {
                EffRes <- cmnsflognormeff_cauchit(object = object,
                  level = level)
                } else {
                if (object$udist == "genexponential") {
                  EffRes <- cmnsfgenexponormeff_cauchit(object = object,
                  level = level)
                } else {
                  if (object$udist == "tslaplace") {
                  EffRes <- cmnsftslnormeff_cauchit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "weibull") {
                    EffRes <- cmnsfweibullnormeff_cauchit(object = object,
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
        }
      }
    } else {
      if (object$linkF == "probit") {
        if (object$sigmavType == "common") {
          if (object$udist == "hnormal") {
          EffRes <- czisfhalfnormeff_probit(object = object,
            level = level)
          } else {
          if (object$udist == "exponential") {
            EffRes <- czisfexponormeff_probit(object = object,
            level = level)
          } else {
            if (object$udist == "gamma") {
            EffRes <- czisfgammanormeff_probit(object = object,
              level = level)
            } else {
            if (object$udist == "rayleigh") {
              EffRes <- czisfraynormeff_probit(object = object,
              level = level)
            } else {
              if (object$udist == "uniform") {
              EffRes <- czisfuninormeff_probit(object = object,
                level = level)
              } else {
              if (object$udist == "tnormal") {
                EffRes <- czisftruncnormeff_probit(object = object,
                level = level)
              } else {
                if (object$udist == "lognormal") {
                EffRes <- czisflognormeff_probit(object = object,
                  level = level)
                } else {
                if (object$udist == "genexponential") {
                  EffRes <- czisfgenexponormeff_probit(object = object,
                  level = level)
                } else {
                  if (object$udist == "tslaplace") {
                  EffRes <- czisftslnormeff_probit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "weibull") {
                    EffRes <- czisfweibullnormeff_probit(object = object,
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
        } else {
          if (object$sigmavType == "different") {
          if (object$udist == "hnormal") {
            EffRes <- cmnsfhalfnormeff_probit(object = object,
            level = level)
          } else {
            if (object$udist == "exponential") {
            EffRes <- cmnsfexponormeff_probit(object = object,
              level = level)
            } else {
            if (object$udist == "gamma") {
              EffRes <- cmnsfgammanormeff_probit(object = object,
              level = level)
            } else {
              if (object$udist == "rayleigh") {
              EffRes <- cmnsfraynormeff_probit(object = object,
                level = level)
              } else {
              if (object$udist == "uniform") {
                EffRes <- cmnsfuninormeff_probit(object = object,
                level = level)
              } else {
                if (object$udist == "tnormal") {
                EffRes <- cmnsftruncnormeff_probit(object = object,
                  level = level)
                } else {
                if (object$udist == "lognormal") {
                  EffRes <- cmnsflognormeff_probit(object = object,
                  level = level)
                } else {
                  if (object$udist == "genexponential") {
                  EffRes <- cmnsfgenexponormeff_probit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "tslaplace") {
                    EffRes <- cmnsftslnormeff_probit(object = object,
                    level = level)
                  } else {
                    if (object$udist == "weibull") {
                    EffRes <- cmnsfweibullnormeff_probit(object = object,
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
          }
        }
      } else {
        if (object$linkF == "cloglog") {
          if (object$sigmavType == "common") {
          if (object$udist == "hnormal") {
            EffRes <- czisfhalfnormeff_cloglog(object = object,
            level = level)
          } else {
            if (object$udist == "exponential") {
            EffRes <- czisfexponormeff_cloglog(object = object,
              level = level)
            } else {
            if (object$udist == "gamma") {
              EffRes <- czisfgammanormeff_cloglog(object = object,
              level = level)
            } else {
              if (object$udist == "rayleigh") {
              EffRes <- czisfraynormeff_cloglog(object = object,
                level = level)
              } else {
              if (object$udist == "uniform") {
                EffRes <- czisfuninormeff_cloglog(object = object,
                level = level)
              } else {
                if (object$udist == "tnormal") {
                EffRes <- czisftruncnormeff_cloglog(object = object,
                  level = level)
                } else {
                if (object$udist == "lognormal") {
                  EffRes <- czisflognormeff_cloglog(object = object,
                  level = level)
                } else {
                  if (object$udist == "genexponential") {
                  EffRes <- czisfgenexponormeff_cloglog(object = object,
                    level = level)
                  } else {
                  if (object$udist == "tslaplace") {
                    EffRes <- czisftslnormeff_cloglog(object = object,
                    level = level)
                  } else {
                    if (object$udist == "weibull") {
                    EffRes <- czisfweibullnormeff_cloglog(object = object,
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
          } else {
          if (object$sigmavType == "different") {
            if (object$udist == "hnormal") {
            EffRes <- cmnsfhalfnormeff_cloglog(object = object,
              level = level)
            } else {
            if (object$udist == "exponential") {
              EffRes <- cmnsfexponormeff_cloglog(object = object,
              level = level)
            } else {
              if (object$udist == "gamma") {
              EffRes <- cmnsfgammanormeff_cloglog(object = object,
                level = level)
              } else {
              if (object$udist == "rayleigh") {
                EffRes <- cmnsfraynormeff_cloglog(object = object,
                level = level)
              } else {
                if (object$udist == "uniform") {
                EffRes <- cmnsfuninormeff_cloglog(object = object,
                  level = level)
                } else {
                if (object$udist == "tnormal") {
                  EffRes <- cmnsftruncnormeff_cloglog(object = object,
                  level = level)
                } else {
                  if (object$udist == "lognormal") {
                  EffRes <- cmnsflognormeff_cloglog(object = object,
                    level = level)
                  } else {
                  if (object$udist == "genexponential") {
                    EffRes <- cmnsfgenexponormeff_cloglog(object = object,
                    level = level)
                  } else {
                    if (object$udist == "tslaplace") {
                    EffRes <- cmnsftslnormeff_cloglog(object = object,
                      level = level)
                    } else {
                    if (object$udist == "weibull") {
                      EffRes <- cmnsfweibullnormeff_cloglog(object = object,
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
          }
          }
        }
      }
    }
  }
  return(data.frame(EffRes))
}

# conditional efficiencies cnsfcross ----------
#' @rdname efficiencies
#' @aliases efficiencies.cnsfcross
#' @export
efficiencies.cnsfcross <- function(object, level = 0.95, newData = NULL,
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
  if (object$linkF == "logit") {
    if (object$sigmauType == "common") {
      if (object$udist == "hnormal") {
        EffRes <- ccnsfhalfnormeff_logit(object = object,
          level = level)
      } else {
        if (object$udist == "exponential") {
          EffRes <- ccnsfexponormeff_logit(object = object,
          level = level)
        } else {
          if (object$udist == "gamma") {
          EffRes <- ccnsfgammanormeff_logit(object = object,
            level = level)
          } else {
          if (object$udist == "rayleigh") {
            EffRes <- ccnsfraynormeff_logit(object = object,
            level = level)
          } else {
            if (object$udist == "uniform") {
            EffRes <- ccnsfuninormeff_logit(object = object,
              level = level)
            } else {
            if (object$udist == "tnormal") {
              EffRes <- ccnsftruncnormeff_logit(object = object,
              level = level)
            } else {
              if (object$udist == "lognormal") {
              EffRes <- ccnsflognormeff_logit(object = object,
                level = level)
              } else {
              if (object$udist == "genexponential") {
                EffRes <- ccnsfgenexponormeff_logit(object = object,
                level = level)
              } else {
                if (object$udist == "tslaplace") {
                EffRes <- ccnsftslnormeff_logit(object = object,
                  level = level)
                } else {
                if (object$udist == "weibull") {
                  EffRes <- ccnsfweibullnormeff_logit(object = object,
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
    } else {
      if (object$sigmauType == "different") {
        if (object$udist == "hnormal") {
          EffRes <- cmcesfhalfnormeff_logit(object = object,
          level = level)
        } else {
          if (object$udist == "exponential") {
          EffRes <- cmcesfexponormeff_logit(object = object,
            level = level)
          } else {
          if (object$udist == "gamma") {
            EffRes <- cmcesfgammanormeff_logit(object = object,
            level = level)
          } else {
            if (object$udist == "rayleigh") {
            EffRes <- cmcesfraynormeff_logit(object = object,
              level = level)
            } else {
            if (object$udist == "uniform") {
              EffRes <- cmcesfuninormeff_logit(object = object,
              level = level)
            } else {
              if (object$udist == "tnormal") {
              EffRes <- cmcesftruncnormeff_logit(object = object,
                level = level)
              } else {
              if (object$udist == "lognormal") {
                EffRes <- cmcesflognormeff_logit(object = object,
                level = level)
              } else {
                if (object$udist == "genexponential") {
                EffRes <- cmcesfgenexponormeff_logit(object = object,
                  level = level)
                } else {
                if (object$udist == "tslaplace") {
                  EffRes <- cmcesftslnormeff_logit(object = object,
                  level = level)
                } else {
                  if (object$udist == "weibull") {
                  EffRes <- cmcesfweibullnormeff_logit(object = object,
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
      }
    }
  } else {
    if (object$linkF == "cauchit") {
      if (object$sigmauType == "common") {
        if (object$udist == "hnormal") {
          EffRes <- ccnsfhalfnormeff_cauchit(object = object,
          level = level)
        } else {
          if (object$udist == "exponential") {
          EffRes <- ccnsfexponormeff_cauchit(object = object,
            level = level)
          } else {
          if (object$udist == "gamma") {
            EffRes <- ccnsfgammanormeff_cauchit(object = object,
            level = level)
          } else {
            if (object$udist == "rayleigh") {
            EffRes <- ccnsfraynormeff_cauchit(object = object,
              level = level)
            } else {
            if (object$udist == "uniform") {
              EffRes <- ccnsfuninormeff_cauchit(object = object,
              level = level)
            } else {
              if (object$udist == "tnormal") {
              EffRes <- ccnsftruncnormeff_cauchit(object = object,
                level = level)
              } else {
              if (object$udist == "lognormal") {
                EffRes <- ccnsflognormeff_cauchit(object = object,
                level = level)
              } else {
                if (object$udist == "genexponential") {
                EffRes <- ccnsfgenexponormeff_cauchit(object = object,
                  level = level)
                } else {
                if (object$udist == "tslaplace") {
                  EffRes <- ccnsftslnormeff_cauchit(object = object,
                  level = level)
                } else {
                  if (object$udist == "weibull") {
                  EffRes <- ccnsfweibullnormeff_cauchit(object = object,
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
      } else {
        if (object$sigmauType == "different") {
          if (object$udist == "hnormal") {
          EffRes <- cmcesfhalfnormeff_cauchit(object = object,
            level = level)
          } else {
          if (object$udist == "exponential") {
            EffRes <- cmcesfexponormeff_cauchit(object = object,
            level = level)
          } else {
            if (object$udist == "gamma") {
            EffRes <- cmcesfgammanormeff_cauchit(object = object,
              level = level)
            } else {
            if (object$udist == "rayleigh") {
              EffRes <- cmcesfraynormeff_cauchit(object = object,
              level = level)
            } else {
              if (object$udist == "uniform") {
              EffRes <- cmcesfuninormeff_cauchit(object = object,
                level = level)
              } else {
              if (object$udist == "tnormal") {
                EffRes <- cmcesftruncnormeff_cauchit(object = object,
                level = level)
              } else {
                if (object$udist == "lognormal") {
                EffRes <- cmcesflognormeff_cauchit(object = object,
                  level = level)
                } else {
                if (object$udist == "genexponential") {
                  EffRes <- cmcesfgenexponormeff_cauchit(object = object,
                  level = level)
                } else {
                  if (object$udist == "tslaplace") {
                  EffRes <- cmcesftslnormeff_cauchit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "weibull") {
                    EffRes <- cmcesfweibullnormeff_cauchit(object = object,
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
        }
      }
    } else {
      if (object$linkF == "probit") {
        if (object$sigmauType == "common") {
          if (object$udist == "hnormal") {
          EffRes <- ccnsfhalfnormeff_probit(object = object,
            level = level)
          } else {
          if (object$udist == "exponential") {
            EffRes <- ccnsfexponormeff_probit(object = object,
            level = level)
          } else {
            if (object$udist == "gamma") {
            EffRes <- ccnsfgammanormeff_probit(object = object,
              level = level)
            } else {
            if (object$udist == "rayleigh") {
              EffRes <- ccnsfraynormeff_probit(object = object,
              level = level)
            } else {
              if (object$udist == "uniform") {
              EffRes <- ccnsfuninormeff_probit(object = object,
                level = level)
              } else {
              if (object$udist == "tnormal") {
                EffRes <- ccnsftruncnormeff_probit(object = object,
                level = level)
              } else {
                if (object$udist == "lognormal") {
                EffRes <- ccnsflognormeff_probit(object = object,
                  level = level)
                } else {
                if (object$udist == "genexponential") {
                  EffRes <- ccnsfgenexponormeff_probit(object = object,
                  level = level)
                } else {
                  if (object$udist == "tslaplace") {
                  EffRes <- ccnsftslnormeff_probit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "weibull") {
                    EffRes <- ccnsfweibullnormeff_probit(object = object,
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
        } else {
          if (object$sigmauType == "different") {
          if (object$udist == "hnormal") {
            EffRes <- cmcesfhalfnormeff_probit(object = object,
            level = level)
          } else {
            if (object$udist == "exponential") {
            EffRes <- cmcesfexponormeff_probit(object = object,
              level = level)
            } else {
            if (object$udist == "gamma") {
              EffRes <- cmcesfgammanormeff_probit(object = object,
              level = level)
            } else {
              if (object$udist == "rayleigh") {
              EffRes <- cmcesfraynormeff_probit(object = object,
                level = level)
              } else {
              if (object$udist == "uniform") {
                EffRes <- cmcesfuninormeff_probit(object = object,
                level = level)
              } else {
                if (object$udist == "tnormal") {
                EffRes <- cmcesftruncnormeff_probit(object = object,
                  level = level)
                } else {
                if (object$udist == "lognormal") {
                  EffRes <- cmcesflognormeff_probit(object = object,
                  level = level)
                } else {
                  if (object$udist == "genexponential") {
                  EffRes <- cmcesfgenexponormeff_probit(object = object,
                    level = level)
                  } else {
                  if (object$udist == "tslaplace") {
                    EffRes <- cmcesftslnormeff_probit(object = object,
                    level = level)
                  } else {
                    if (object$udist == "weibull") {
                    EffRes <- cmcesfweibullnormeff_probit(object = object,
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
          }
        }
      } else {
        if (object$linkF == "cloglog") {
          if (object$sigmauType == "common") {
          if (object$udist == "hnormal") {
            EffRes <- ccnsfhalfnormeff_cloglog(object = object,
            level = level)
          } else {
            if (object$udist == "exponential") {
            EffRes <- ccnsfexponormeff_cloglog(object = object,
              level = level)
            } else {
            if (object$udist == "gamma") {
              EffRes <- ccnsfgammanormeff_cloglog(object = object,
              level = level)
            } else {
              if (object$udist == "rayleigh") {
              EffRes <- ccnsfraynormeff_cloglog(object = object,
                level = level)
              } else {
              if (object$udist == "uniform") {
                EffRes <- ccnsfuninormeff_cloglog(object = object,
                level = level)
              } else {
                if (object$udist == "tnormal") {
                EffRes <- ccnsftruncnormeff_cloglog(object = object,
                  level = level)
                } else {
                if (object$udist == "lognormal") {
                  EffRes <- ccnsflognormeff_cloglog(object = object,
                  level = level)
                } else {
                  if (object$udist == "genexponential") {
                  EffRes <- ccnsfgenexponormeff_cloglog(object = object,
                    level = level)
                  } else {
                  if (object$udist == "tslaplace") {
                    EffRes <- ccnsftslnormeff_cloglog(object = object,
                    level = level)
                  } else {
                    if (object$udist == "weibull") {
                    EffRes <- ccnsfweibullnormeff_cloglog(object = object,
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
          } else {
          if (object$sigmauType == "different") {
            if (object$udist == "hnormal") {
            EffRes <- cmcesfhalfnormeff_cloglog(object = object,
              level = level)
            } else {
            if (object$udist == "exponential") {
              EffRes <- cmcesfexponormeff_cloglog(object = object,
              level = level)
            } else {
              if (object$udist == "gamma") {
              EffRes <- cmcesfgammanormeff_cloglog(object = object,
                level = level)
              } else {
              if (object$udist == "rayleigh") {
                EffRes <- cmcesfraynormeff_cloglog(object = object,
                level = level)
              } else {
                if (object$udist == "uniform") {
                EffRes <- cmcesfuninormeff_cloglog(object = object,
                  level = level)
                } else {
                if (object$udist == "tnormal") {
                  EffRes <- cmcesftruncnormeff_cloglog(object = object,
                  level = level)
                } else {
                  if (object$udist == "lognormal") {
                  EffRes <- cmcesflognormeff_cloglog(object = object,
                    level = level)
                  } else {
                  if (object$udist == "genexponential") {
                    EffRes <- cmcesfgenexponormeff_cloglog(object = object,
                    level = level)
                  } else {
                    if (object$udist == "tslaplace") {
                    EffRes <- cmcesftslnormeff_cloglog(object = object,
                      level = level)
                    } else {
                    if (object$udist == "weibull") {
                      EffRes <- cmcesfweibullnormeff_cloglog(object = object,
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
          }
          }
        }
      }
    }
  }
  return(data.frame(EffRes))
}
