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

#' Extract frontier information to be used with \bold{texreg} package
#' 
#' Extract coefficients and additional information for stochastic frontier 
#' models returned by \code{\link{sfacross}}, \code{\link{sfalcmcross}}, 
#' \code{\link{sfagzisfcross}}, \code{\link{sfacnsfcross}}, 
#' \code{\link{sfamisfcross}}, \code{\link{sfazisfcross}}, 
#' \code{\link{sfametacross}}, \code{\link{sfaselectioncross}}, 
#' \code{\link{sfapanel1}}, or \code{\link{sfalcmpanel}}.
#' 
#' @name extract
#' 
#' @aliases extract.sfacross extract.sfalcmcross extract.sfaselectioncross
#' extract.sfacnsfcross extract.sfalcmcross extract.sfamisfcross
#' extract.sfazisfcross extract.sfametacross extract.sfapanel1 
#' extract.sfalcmpanel
#'
#' @param model objects of class \code{'sfacross'}, \code{'sfalcmcross'}, 
#' \code{'sfagzisfcross'}, \code{'sfacnsfcross'}, \code{'sfamisfcross'}, 
#' \code{'sfazisfcross'}, \code{'sfametacross'}, \code{'sfaselectioncross'}, 
#' \code{'sfapanel1'}, or \code{'sfalcmpanel'}
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
#' model fitting function using cross-sectional data.
#' 
#' \code{\link{sfapanel1}}, for the first generation stochastic frontier 
#' analysis model fitting function using panel data.
#' 
#' \code{\link{sfalcmpanel}}, for the latent class stochastic frontier analysis
#' model fitting function using panel data.
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
#' udist = 'exponential', uhet = ~ regu, data = utility, S = -1, 
#' method = 'bfgs')
#'
#' texreg::screenreg(list(hlf, trnorm, tscal, expo))
#' 
# @export @rawNamespace if (requireNamespace('texreg')) importFrom(texreg,
# extract)
#' @exportS3Method texreg::extract sfacross
# sfacross extraction for texreg ------
extract.sfacross <- function(model, ...) {
  co <- coef.sfacross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = if (model$udist == "tnormal" && model$scaling == TRUE)
      "trunc. scal." else model$udist)
  return(tr)
}

setGeneric("extract", function(model, ...) standardGeneric("extract"), package = "texreg")

setMethod("extract", signature = className("sfacross", "sfaR"), definition = extract.sfacross)

# sfalcmcross extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfalcmcross
extract.sfalcmcross <- function(model, ...) {
  co <- coef.sfalcmcross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("LCM", model$nClasses, " Classes"))
  return(tr)
}

setMethod("extract", signature = className("sfalcmcross", "sfaR"), definition = extract.sfalcmcross)

# sfaselectioncross extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfaselectioncross
extract.sfaselectioncross <- function(model, ...) {
  co <- coef.sfaselectioncross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = model$lType)
  return(tr)
}

setMethod("extract", signature = className("sfaselectioncross", "sfaR"), definition = extract.sfaselectioncross)

# sfacnsfcross extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfacnsfcross
extract.sfacnsfcross <- function(model, ...) {
  co <- coef.sfacnsfcross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("cnsf"))
  return(tr)
}

setMethod("extract", signature = className("sfacnsfcross", "sfaR"), definition = extract.sfacnsfcross)

# sfamisfcross extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfamisfcross
extract.sfamisfcross <- function(model, ...) {
  co <- coef.sfamisfcross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("misf"))
  return(tr)
}

setMethod("extract", signature = className("sfamisfcross", "sfaR"), definition = extract.sfamisfcross)

# sfazisfcross extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfazisfcross
extract.sfazisfcross <- function(model, ...) {
  co <- coef.sfazisfcross(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("zisf"))
  return(tr)
}

setMethod("extract", signature = className("sfazisfcross", "sfaR"), definition = extract.sfazisfcross)

# sfametacross extraction for texreg ------
extract.sfametacross <- function(model, ...) {
  group_var <- model$dataTable[[model$Ngroup + 1]][model$name_meta_var][, 1]
  group_var_list <- sort(unique(group_var))
  frChoice <- displayMenu(c(group_var_list, "metafrontier"), title = "Which frontier do you want? ")
  if ((frChoice == model$Ngroup + 1) && model$modelType %in% c("aos17a", "aos17b",
    "aos17c", "aos17d")) {
    stop("No parameters are estimated for the metafontier for models 'aos17a', 'aos17b', 'aos17c', 'aos17d' \n", 
         call. = FALSE)
  }
  co <- coef.sfametacross(model)[, frChoice]
  names <- names(co)
  se <- sqrt(diag(model$invHessian[[frChoice]]))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = "meta")
  return(tr)
}

setMethod("extract", signature = className("sfametacross", "sfaR"), definition = extract.sfametacross)

# sfapanel1 extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfapanel1
extract.sfapanel1 <- function(model, ...) {
  co <- coef.sfapanel1(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("panel:", model$udist))
  return(tr)
}

setMethod("extract", signature = className("sfapanel1", "sfaR"), definition = extract.sfapanel1)

# sfalcmpanel extraction for texreg ------
#' @rdname extract
#' @exportS3Method texreg::extract sfalcmpanel
extract.sfalcmpanel <- function(model, ...) {
  co <- coef.sfalcmpanel(model)
  names <- names(co)
  se <- sqrt(diag(model$invHessian))
  pval <- 2 * pnorm(-abs(co/se))
  gof <- c(-2 * model$mlLoglik + 2 * model$nParm, -2 * model$mlLoglik + log(model$Nobs) *
    model$nParm, model$mlLoglik, model$Nobs)
  gof.names <- c("AIC", "BIC", "log-likelihood", "Num. obs.")
  tr <- texreg::createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
    gof.names = gof.names, gof = gof, gof.decimal = c(TRUE, TRUE, TRUE, FALSE),
    model.name = paste0("LCM", model$nClasses, " Classes"))
  return(tr)
}

setMethod("extract", signature = className("sfalcmpanel", "sfaR"), definition = extract.sfalcmpanel)
