################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Fitted values of models                                                      #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
#         -Contaminated noise stochastic frontier                              #
#         -Multi-Modal Inefficiency Stochastic Frontier Analysis               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract fitted values of stochastic frontier models
#'
#' \code{\link{fitted}} returns the fitted frontier values from classic or
#' latent class stochastic frontier models estimated with
#' \code{\link{cnsfcross}}, \code{\link{lcmcross}}, \code{\link{misfcross}}, 
#' \code{\link{sfacross}}, \code{\link{sfaselectioncross}} or 
#' \code{\link{zisfcross}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{cnsfcross}}, \code{\link{lcmcross}}, \code{\link{misfcross}}, 
#' \code{\link{sfacross}}, \code{\link{sfaselectioncross}} or 
#' \code{\link{zisfcross}}.
#' @param ... Currently ignored.
#'
#' @name fitted
#'
#' @return In the case of an object of class \code{'lcmcross'}, a data frame 
#' containing the fitted values for each class is returned where each variable 
#' terminates with \code{'_c#'}, \code{'#'} being the class number.
#' 
#' In the case of an object of class \code{'cnsfcross'}, \code{'misfcross'},
#' \code{'sfacross'}, \code{'sfaselectioncross'} or \code{'zisfcross'}, 
#' a vector of fitted values is returned.
#'
#' @note The fitted values are ordered in the same way as the corresponding
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
#' @keywords methods fitted
#'
#' @examples
#'
#' ## Using data on eighty-two countries production (DGP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod)
#' fit.cb_2c_h <- fitted(cb_2c_h)
#' head(fit.cb_2c_h)
#'
#' @aliases fitted.sfacross
#' @export
# fitted values for sfacross ----------
fitted.sfacross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for lcmcross ----------
#' @rdname fitted
#' @aliases fitted.lcmcross
#' @export
fitted.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(dplyr::select(object$dataTable, "mlFitted_c1",
      "mlFitted_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(dplyr::select(object$dataTable, "mlFitted_c1",
        "mlFitted_c2", "mlFitted_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(dplyr::select(object$dataTable, "mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3", "mlFitted_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(dplyr::select(object$dataTable,
          "mlFitted_c1", "mlFitted_c2", "mlFitted_c3",
          "mlFitted_c4", "mlFitted_c5"))
        }
      }
    }
  }
}

# fitted values for sfaselectioncross ----------
#' @rdname fitted
#' @aliases fitted.sfaselectioncross
#' @export
fitted.sfaselectioncross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for zisfcross ----------
#' @rdname fitted
#' @aliases fitted.zisfcross
#' @export
fitted.zisfcross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for cnsfcross ----------
#' @rdname fitted
#' @aliases fitted.cnsfcross
#' @export
fitted.cnsfcross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for misfcross ----------
#' @rdname fitted
#' @aliases fitted.misfcross
#' @export
fitted.misfcross <- function(object, ...) {
  object$dataTable$mlFitted
}
