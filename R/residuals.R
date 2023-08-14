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

#' Extract residuals of stochastic frontier models
#'
#' This function returns the residuals' values from stochastic frontier models 
#' estimated with \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' \code{\link{sfagzisfcross}}, \code{\link{sfacnsfcross}}, 
#' \code{\link{sfamisfcross}}, \code{\link{sfazisfcross}}, 
#' \code{\link{sfametacross}}, \code{\link{sfaselectioncross}}, 
#' \code{\link{sfapanel1}}, or \code{\link{sfalcmpanel}}.
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, 
#' \code{\link{sfaselectioncross}}, \code{\link{sfapanel1}}, or 
#' \code{\link{sfalcmpanel}}.
#' @param \dots Currently ignored.
#'
#' @name residuals
#' 
#' @aliases residuals.sfacross residuals.sfalcmcross residuals.sfagzisfcross 
#' residuals.sfacnsfcross residuals.sfamisfcross residuals.sfazisfcross 
#' residuals.sfaselectioncross residuals.sfametacross residuals.sfapanel1
#' residuals.sfalcmpanel
#'
#' @details Let's consider the standard stochastic frontier model
#'
#' \deqn{y_i = \alpha + \mathbf{x_i^{\prime}}\bm{\beta} + v_i - Su_i}
#' 
#' The fitted frontier is obtained as
#' 
#' \deqn{\hat{f}_i=\mathbf{x_i^{\prime}}\bm{\hat{\beta}}}
#' 
#' then the residuals is obtained as
#' 
#' \deqn{\hat{\epsilon}_i=y_i - \hat{f}_i}
#'
#' @return When the \code{object} is of class \code{'sfacross'}, 
#' code{'sfacnsfcross'}, \code{'sfamisfcross'}, \code{'sfazisfcross'}, 
#' \code{'sfaselectioncross'}, or \code{'sfapanel1'}, \code{\link{residuals}} 
#' returns a vector of residuals values.
#' 
#' When the \code{object} is of class \code{'sfalcmcross'}, 
#' \code{'sfagzisfcross'}, or \code{'sfalcmpanel'}, \code{\link{residuals}} 
#' returns a data frame containing the residuals values for each latent class, 
#' where each variable ends with \code{'_c#'}, \code{'#'} being the class number.
#' 
#' For an object of class \code{'sfametacross'}, a data frame containing the 
#' residuals values for each each group and the metafrontier is returned. Each 
#' variable name ends with \code{'_g'} where \code{'g'} is the group name. For 
#' the metafrontier \code{'g'} = 'metafrontier'. 
#' For \code{'modelType'} = 'aos17a' or \code{'modelType'} = 'aos17b', the 
#' residuals for the metafrontier are simply the mean values of all the 
#' residuals from the simulated metafrontiers.
#'
#' @note The residuals values are ordered in the same way as the corresponding
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
#' @keywords methods residuals
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
#' resid.tl_u_ts <- residuals(tl_u_ts)
#' head(resid.tl_u_ts)
#'
#' # Using data on eighty-two countries production (GDP)
#' ## LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal', 
#' data = worldprod, S = 1)
#' resid.cb_2c_h <- residuals(cb_2c_h)
#' head(resid.cb_2c_h)
#' }
#'
#' @export
# @exportS3Method residuals sfacross
# residuals from sfacross ----------
residuals.sfacross <- function(object, ...) {
  object$dataTable$mlResiduals
}

# residuals from sfalcmcross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfalcmcross
residuals.sfalcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    object$dataTable[c("mlResiduals_c1", "mlResiduals_c2")]
  } else {
    if (object$nClasses == 3) {
      object$dataTable[c("mlResiduals_c1", "mlResiduals_c2",
        "mlResiduals_c3")]
    } else {
      if (object$nClasses == 4) {
        object$dataTable[c("mlResiduals_c1", "mlResiduals_c2",
          "mlResiduals_c3", "mlResiduals_c4")]
      } else {
        if (object$nClasses == 5) {
          object$dataTable[c("mlResiduals_c1", "mlResiduals_c2",
          "mlResiduals_c3", "mlResiduals_c4", "mlResiduals_c5")]
        }
      }
    }
  }
}

# residuals from sfagzisfcross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfagzisfcross
residuals.sfagzisfcross <- function(...) {
  residuals.sfalcmcross(...)
}

# residuals from sfacnsfcross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfacnsfcross
residuals.sfacnsfcross <- function(...) {
  residuals.sfacross(...)
}

# residuals from sfamisfcross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfamisfcross
residuals.sfamisfcross <- function(...) {
  residuals.sfacross(...)
}

# residuals from sfazisfcross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfazisfcross
residuals.sfazisfcross <- function(...) {
  residuals.sfacross(...)
}

# residuals from sfaselectioncross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfaselectioncross
residuals.sfaselectioncross <- function(...) {
  residuals.sfacross(...)
}

# residuals from sfametacross ----------
#' @rdname residuals
#' @export
# @exportS3Method residuals sfametacross
residuals.sfametacross <- function(object, ...) {
  resiMat <- matrix(nrow = object$Nobs[object$Ngroup + 1],
    ncol = object$Ngroup + 1)
  namelist <- names(object$dataTable)
  group_var <- object$dataTable[[object$Ngroup + 1]][object$name_meta_var][,
    1]
  group_var_list <- sort(unique(group_var))
  for (g in group_var_list) {
    resiMat[group_var == g, which(group_var_list == g)] <- object$dataTable[[which(group_var_list ==
      g)]][, "mlResiduals"]
  }
  resiMat[, object$Ngroup + 1] <- object$dataTable[[object$Ngroup +
    1]][, "mlResiduals"]
  resiMat <- as.data.frame(resiMat)
  names(resiMat) <- paste0("mlResiduals_", namelist)
  resiMat
}

# residuals from sfapanel1 ----------
#' @rdname residuals
#' @aliases residuals.sfapanel1
#' @export
residuals.sfapanel1 <- function(...) {
  residuals.sfacross(...)
}

# residuals from sfalcmpanel ----------
#' @rdname residuals
#' @aliases residuals.sfalcmpanel
#' @export
residuals.sfalcmpanel <- function(...) {
  residuals.sfalcmcross(...)
}
