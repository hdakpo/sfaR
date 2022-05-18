################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Residuals of model (v - S * u)                                               #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
#         -Contaminated noise stochastic frontier                              #
#         -Multi-Modal Inefficiency Stochastic Frontier Analysis               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract residuals of stochastic frontier models
#'
#' This function returns the residuals' values from classic or latent class
#' stochastic frontier models estimated with \code{\link{cnsfcross}},
#' \code{\link{lcmcross}}, \code{\link{misfcross}}, \code{\link{sfacross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}.
#'
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{cnsfcross}}, \code{\link{lcmcross}}, \code{\link{cnsfcross}}, 
#' \code{\link{sfacross}}, \code{\link{sfaselectioncross}} or 
#' \code{\link{zisfcross}}.
#' @param ... Currently ignored.
#'
#' @name residuals
#'
#' @return When the \code{object} is of class \code{'lcmcross'},
#' \code{\link{residuals}} returns a data frame containing the residuals values
#' for each latent class, where each variable terminates with \code{'_c#'},
#' \code{'#'} being the class number.
#' 
#' When the \code{object} is of class \code{'cnsfcross'}, 
#' \code{'misfcross'}, \code{'sfacross'}, \code{'sfaselectioncross'} or 
#' \code{'zisfcross'}, \code{\link{residuals}} returns a vector of residuals 
#' values.
#'
#' @note The residuals values are ordered in the same way as the corresponding
#' observations in the dataset used for the estimation.
#'
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{cnsfcross}}, for the contaminated noise stochastic 
#' frontier analysis model fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{misfcross}}, for the multi-modal inefficiency stochastic frontier 
#' analysis model fitting function.
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
    data.frame(select(object$dataTable, "mlResiduals_c1",
      "mlResiduals_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mlResiduals_c1",
        "mlResiduals_c2", "mlResiduals_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3", "mlResiduals_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3", "mlResiduals_c4",
          "mlResiduals_c5"))
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

# residuals from zisfcross ----------
#' @rdname residuals
#' @aliases residuals.zisfcross
#' @export
residuals.zisfcross <- function(object, ...) {
  object$dataTable$mlResiduals
}

# residuals from cnsfcross ----------
#' @rdname residuals
#' @aliases residuals.cnsfcross
#' @export
residuals.cnsfcross <- function(object, ...) {
  object$dataTable$mlResiduals
}

# residuals from misfcross ----------
#' @rdname residuals
#' @aliases residuals.misfcross
#' @export
residuals.misfcross <- function(object, ...) {
  object$dataTable$mlResiduals
}
