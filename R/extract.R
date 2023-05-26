################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Extract information to be used by texreg R package                           #
# Models: + Cross sectional & Pooled data                                      #
#           -Stochastic Frontier Analysis                                      #
#           -Latent Class Stochastic Frontier Analysis                         #
#           -Sample selection correction for Stochastic Frontier Model         #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Extract frontier information to be used with \bold{texreg} package
#' 
#' Extract coefficients and additional information for stochastic frontier models 
#' returned by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, or 
#' \code{\link{sfaselectioncross}}.
#' 
#' @name extract
#'
#' @param model objects of class \code{'sfacross'}, \code{'sfalcmcross'}, or
#'  \code{'sfaselectioncross'}
#' @param ... Currently ignored
#' 
#' @return A texreg object representing the statistical model.
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
#' model fitting function using cross-sectional data.
#' 
#' @keywords methods extract
#'
#' @examples
#' 
#' hlf <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) + 
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'hnormal', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#
#' trnorm <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) + 
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#' 
#' tscal <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, 
#' S = -1, method = 'bfgs', scaling = TRUE)
#' 
#' expo <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'exponential', uhet = ~ regu, data = utility, S = -1, method = 'bfgs')
#'
#' texreg::screenreg(list(hlf, trnorm, tscal, expo))
#' 
#' @aliases extract.sfacross
#' @export
# results extraction for texreg ------
extract.sfacross <- function(model, ...) {
  co <- coef.sfacross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik +
    log(model$Nobs) * model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co,
    se = se, pvalues = pval, gof.names = gof.names, gof = gof,
    gof.decimal = c(TRUE, TRUE, TRUE, FALSE), model.name = if (model$udist ==
      "tnormal" && model$scaling == TRUE)
      "trunc. scal." else model$udist)
  return(tr)
}

setGeneric("extract", function(model, ...) standardGeneric("extract"),
           package = "texreg")

setMethod("extract", signature = className("sfacross", "sfaR"),
  definition = extract.sfacross)

#' @rdname extract
#' @aliases extract.sfalcmcross
#' @export
extract.sfalcmcross <- function(model, ...) {
  co <- coef.sfalcmcross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik +
    log(model$Nobs) * model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co,
    se = se, pvalues = pval, gof.names = gof.names, gof = gof,
    gof.decimal = c(TRUE, TRUE, TRUE, FALSE), model.name = paste0(model$nClasses,
      " Classes"))
  return(tr)
}

setMethod("extract", signature = className("sfalcmcross", "sfaR"),
  definition = extract.sfalcmcross)

#' @rdname extract
#' @aliases extract.sfaselectioncross
#' @export
extract.sfaselectioncross <- function(model, ...) {
  co <- coef.sfaselectioncross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik +
    log(model$Nobs) * model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co,
    se = se, pvalues = pval, gof.names = gof.names, gof = gof,
    gof.decimal = c(TRUE, TRUE, TRUE, FALSE), model.name = model$lType)
  return(tr)
}

setMethod("extract", signature = className("sfaselectioncross",
  "sfaR"), definition = extract.sfaselectioncross)
