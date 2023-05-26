################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Skewness tests                                                               #
# Models: -Standard Stochastic Frontier Analysis                               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Skewness test for stochastic frontier models
#' 
#' \code{\link{skewnessTest}} computes skewness test for stochastic frontier
#' models (i.e. objects of class \code{'sfacross'}).
#'
#' @param object An object of class \code{'sfacross'}, returned by
#' \code{\link{sfacross}}.
#' @param test A character string specifying the test to implement. If
#' \code{'agostino'} (default), D'Agostino skewness test is implemented
#' (D'Agostino and Pearson, 1973).  If \code{'coelli'}, Coelli skewness test is
#' implemented (Coelli, 1995).
#'
#' @return \code{skewnessTest} returns the results of either the D'Agostino's
#' or the Coelli's skewness test.
#'
#' @note \code{\link{skewnessTest}} is currently only available for object of
#' class \code{'sfacross'}.
#' 
# @author K Hervé Dakpo, Yann Desjeux, and Laure Latruffe
#'
#' @references Coelli, T. 1995. Estimators and hypothesis tests for a
#' stochastic frontier function - a Monte-Carlo analysis. \emph{Journal of
#' Productivity Analysis}, \bold{6}:247--268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613--622.
#'
#' @keywords methods
#'
#' @examples
#' 
#' \dontrun{
#' ## Using data on fossil fuel fired steam electric power generation plants in the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu, uhet = ~ regu, data = utility, S = -1,
#' scaling = TRUE, method = 'mla')
#' skewnessTest(tl_u_ts)
#' skewnessTest(tl_u_ts, test = 'coelli')
#' }
#'
#' @export
skewnessTest <- function(object, test = "agostino") {
  if (!inherits(object, "sfacross")) {
    stop("Argument 'object' must be of class 'sfacross'")
  }
  if (test == "agostino") {
    object$AgostinoTest
  } else {
    if (test == "coelli") {
      tt <- list(data.name = deparse(substitute(object)),
        statistic = object$CoelliM3Test["z"], p.value = object$CoelliM3Test["p.value"],
        method = "## Coelli's test ##")
      class(tt) <- "htest"
      return(tt)
    } else {
      stop("argument 'test' must be either 'agostino', or 'coelli'",
        call. = FALSE)
    }
  }
}
