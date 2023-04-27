################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Residuals of model (v - S * u)                                               #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract residuals of stochastic frontier models
#'
#' This function returns the residuals' values from stochastic frontier models 
#' estimated with \code{\link{lcmcross}}, \code{\link{sfacross}}, or
#'  \code{\link{sfaselectioncross}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{lcmcross}}, \code{\link{sfacross}}, or
#' \code{\link{sfaselectioncross}}.
#' @param \dots Currently ignored.
#'
#' @name residuals
#'
#' @return When the \code{object} is of \code{'lcmcross'}, 
#' \code{\link{residuals}} returns a data frame containing the residuals values 
#' for each latent class, where each variable terminates with \code{'_c#'}, 
#' \code{'#'} being the class number.
#' 
#' When the \code{object} is of class \code{'sfacross'}, or 
#' \code{'sfaselectioncross'}, \code{\link{residuals}} returns a vector of 
#' residuals values.
#'
#' @note The residuals values are ordered in the same way as the corresponding
#' observations in the dataset used for the estimation.
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function using cross-sectional or pooled data.
#'
#' @keywords methods residuals
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' resid.tl_u_ts <- residuals(tl_u_ts)
#' head(resid.tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod, S = 1)
#' resid.cb_2c_h <- residuals(cb_2c_h)
#' head(resid.cb_2c_h)
#'
#' @aliases residuals.sfacross
#' @export
# residuals from sfacross ----------
residuals.sfacross <- function(object, ...) {
  object$dataTable$mlResiduals
}

# residuals from lcmcross ----------
#' @rdname residuals
#' @aliases residuals.lcmcross
#' @export
residuals.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    object$dataTable[c("mlResiduals_c1",
      "mlResiduals_c2")]
  } else {
    if (object$nClasses == 3) {
      object$dataTable[c("mlResiduals_c1",
        "mlResiduals_c2", "mlResiduals_c3")]
    } else {
      if (object$nClasses == 4) {
        object$dataTable[c("mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3", "mlResiduals_c4")]
      } else {
        if (object$nClasses == 5) {
          object$dataTable[c("mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3", "mlResiduals_c4",
          "mlResiduals_c5")]
        }
      }
    }
  }
}

# residuals from sfaselectioncross ----------
#' @rdname residuals
#' @aliases residuals.sfaselectioncross
#' @export
residuals.sfaselectioncross <- function(object, ...) {
  object$dataTable$mlResiduals
}
