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

#' Extract fitted values of stochastic frontier models
#'
#' \code{\link{fitted}} returns the fitted frontier values from stochastic 
#' frontier models estimated with \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' \code{\link{sfagzisfcross}}, \code{\link{sfacnsfcross}}, 
#' \code{\link{sfamisfcross}}, \code{\link{sfazisfcross}}, 
#' \code{\link{sfametacross}}, \code{\link{sfaselectioncross}}, 
#' \code{\link{sfapanel1}}, or \code{\link{sfalcmpanel}}.
#' @param ... Currently ignored.
#'
#' @name fitted
#' 
#' @aliases fitted.sfacross fitted.sfalcmcross fitted.sfagzisfcross
#' fitted.sfacnsfcross fitted.sfamisfcross fitted.sfazisfcross
#' fitted.sfaselectioncross fitted.sfametacross fitted.sfapanel1 
#' fitted.sfalcmpanel
#' 
#' @details Let's consider the standard stochastic frontier model
#' 
#' \deqn{y_i = \alpha + \mathbf{x_i^{\prime}}\bm{\beta} + v_i - Su_i}
#' 
#' The fitted frontier is obtained as
#' 
#' \deqn{\hat{f}_i=\mathbf{x_i^{\prime}}\bm{\hat{\beta}}}
#'
#' @return In the case of an object of class \code{'sfacross'}, 
#' code{'sfacnsfcross'}, \code{'sfamisfcross'}, \code{'sfazisfcross'}, 
#' \code{'sfaselectioncross'}, or \code{'sfapanel1'}, a vector of fitted values 
#' is returned.
#' 
#' In the case of an object of class \code{'sfalcmcross'} 
#' \code{'sfagzisfcross'}, or \code{'sfalcmpanel'}, a data frame containing the 
#' fitted values for each class is returned where each variable ends with 
#' \code{'_c#'}, \code{'#'} being the class number.
#' 
#' For object of class \code{'sfametacross'}, a data frame containing the fitted 
#' values for each each group and the metafrontier is returned. Each variable 
#' name ends with \code{'_g'} where \code{'g'} is the group name. For the 
#' metafrontier \code{'g'} = 'metafrontier'. For \code{'modelType'} = 'aos17a' 
#' or \code{'modelType'} = 'aos17b', the fitted metafrontier is simply the mean
#' values of all the simulated metafrontiers.
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
#' @keywords methods fitted
#'
#' @examples
#' 
#' \dontrun{
#' # Using data on eighty-two countries production (GDP)
#' ## LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod)
#' fit.cb_2c_h <- fitted(cb_2c_h)
#' head(fit.cb_2c_h)
#' }
#' 
#' @export
# @exportS3Method fitted sfacross
# fitted values for sfacross ----------
fitted.sfacross <- function(object, ...) {
  object$dataTable$mlFitted
}

# fitted values for sfalcmcross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfalcmcross
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

# fitted values for sfagzisfcross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfagzisfcross
fitted.sfagzisfcross <- function(...) {
  fitted.sfalcmcross(...)
}

# fitted values for sfacnsfcross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfacnsfcross
fitted.sfacnsfcross <- function(...) {
  fitted.sfacross(...)
}

# fitted values for sfamisfcross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfamisfcross
fitted.sfamisfcross <- function(...) {
  fitted.sfacross(...)
}

# fitted values for sfazisfcross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfazisfcross
fitted.sfazisfcross <- function(...) {
  fitted.sfacross(...)
}

# fitted values for sfaselectioncross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfaselectioncross
fitted.sfaselectioncross <- function(...) {
  fitted.sfacross(...)
}

# fitted values for sfametacross ----------
#' @rdname fitted
#' @export
# @exportS3Method fitted sfametacross
fitted.sfametacross <- function(object, ...) {
  fitMat <- matrix(nrow = object$Nobs[object$Ngroup + 1], ncol = object$Ngroup +
    1)
  namelist <- names(object$dataTable)
  group_var <- object$dataTable[[object$Ngroup + 1]][object$name_meta_var][,
    1]
  group_var_list <- sort(unique(group_var))
  for (g in group_var_list) {
    fitMat[group_var == g, which(group_var_list == g)] <- object$dataTable[[which(group_var_list ==
      g)]][, "mlFitted"]
  }
  fitMat[, object$Ngroup + 1] <- object$dataTable[[object$Ngroup +
    1]][, "mlFitted"]
  fitMat <- as.data.frame(fitMat)
  names(fitMat) <- paste0("mlFitted_", namelist)
  fitMat
}

# fitted values for sfapanel1 ----------
#' @rdname fitted
#' @aliases fitted.sfapanel1
#' @export
fitted.sfapanel1 <- function(...) {
  fitted.sfacross(...)
}

# fitted values for sfalcmpanel ----------
#' @rdname fitted
#' @aliases fitted.sfalcmpanel
#' @export
fitted.sfalcmpanel <- function(...) {
  fitted.sfalcmcross(...)
}
