################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Fitted values of models                                                      #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract fitted values of stochastic frontier models
#'
#' \code{\link{fitted}} returns the fitted frontier values from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' or \code{\link{sfaselectioncross}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, or 
#' \code{\link{sfaselectioncross}}.
#' @param ... Currently ignored.
#'
#' @name fitted
#'
#' @return In the case of an object of class \code{'sfacross'}, or 
#' \code{'sfaselectioncross'}, a vector of fitted values is returned.
#' 
#' In the case of an object of class \code{'sfalcmcross'}, a data frame 
#' containing the fitted values for each class is returned where each variable 
#' ends with \code{'_c#'}, \code{'#'} being the class number.
#'
#' @note The fitted values are ordered in the same way as the corresponding
#' observations in the dataset used for the estimation.
#' 
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfalcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function using cross-sectional or pooled data.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function using cross-sectional or pooled data.
#'
#' @keywords methods fitted
#'
#' @examples
#' 
#' \dontrun{
#' ## Using data on eighty-two countries production (GDP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod)
#' fit.cb_2c_h <- fitted(cb_2c_h)
#' head(fit.cb_2c_h)
#' }
#'
#' @aliases fitted.sfacross
#' @export
# fitted values for sfacross ----------
fitted.sfacross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for sfalcmcross ----------
#' @rdname fitted
#' @aliases fitted.sfalcmcross
#' @export
fitted.sfalcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(object$dataTable[, c("mlFitted_c1", "mlFitted_c2")])
  } else {
    if (object$nClasses == 3) {
      data.frame(object$dataTable[, c("mlFitted_c1", "mlFitted_c2",
        "mlFitted_c3")])
    } else {
      if (object$nClasses == 4) {
        data.frame(object$dataTable[, c("mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3", "mlFitted_c4")])
      } else {
        if (object$nClasses == 5) {
          data.frame(object$dataTable[, c("mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3", "mlFitted_c4",
          "mlFitted_c5")])
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
