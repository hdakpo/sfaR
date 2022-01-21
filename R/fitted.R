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
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract fitted values of stochastic frontier models
#'
#' \code{\link{fitted}} returns the fitted frontier values from classic or
#' latent class stochastic frontier models estimated with
#' \code{\link{sfacross}}, \code{\link{lcmcross}}, \code{\link{selectioncross}} or \code{\link{zisfcross}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{lcmcross}}, \code{\link{selectioncross}} or \code{\link{zisfcross}}.
#' @param ... Currently ignored.
#'
#' @name fitted
#'
#' @return In the case of an object of class \code{'sfacross'}, \code{'selectioncross'} or 
#' \code{'selectioncross'}, a vector of
#' fitted values is returned.
#'
#' In the case of an object of class \code{'lcmcross'}, a data frame containing
#' the fitted values for each class is returned where each variable terminates
#' with \code{'_c#'}, \code{'#'} being the class number.
#'
#' @note The fitted values are ordered in the same way as the corresponding
#' observations in the dataset used for the estimation.
#' 
# @author K Hervé Dakpo, Yann Desjeux, and Laure Latruffe
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
#' \code{\link{zisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function.
#'
#' @keywords methods fitted
#'
#' @examples
#'
#' ## Using data on eighty-two countries production (DGP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', data = worldprod)
#'   fit.cb_2c_h <- fitted(cb_2c_h)
#'   head(fit.cb_2c_h)
#'
#' @aliases fitted.sfacross
#' @export
# fitted values for sfacross ----------
fitted.sfacross <- function(object, ...) {
  object$dataTable$mleFitted
}

# fitted values for lcmcross ----------
#' @rdname fitted
#' @aliases fitted.lcmcross
#' @export
fitted.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(select(object$dataTable, "mlFitted_c1", "mlFitted_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mlFitted_c1",
        "mlFitted_c2", "mlFitted_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3", "mlFitted_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3", "mlFitted_c4",
          "mlFitted_c5"))
        }
      }
    }
  }
}

# fitted values for selectioncross ----------
#' @rdname fitted
#' @aliases fitted.selectioncross
#' @export
fitted.selectioncross <- function(object, ...) {
  object$dataTable$mleFitted
}

# fitted values for zisfcross ----------
#' @rdname fitted
#' @aliases fitted.zisfcross
#' @export
fitted.zisfcross <- function(object, ...) {
  object$dataTable$mleFitted
}
