################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Marginal impact of Z variables on inefficiency                               #
# Models: -Standard Stochastic Frontier Analysis                               #
#         -Latent Class Stochastic Frontier Analysis                           #
#         -Sample selection correction                                         #
#         -Zero inefficiency stochastic frontier                               #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Marginal effects of the inefficiency drivers in stochastic frontier models
#'
#' This function returns marginal effects of the inefficiency drivers from
#' classic or latent class stochastic frontier models estimated with
#' \code{\link{sfacross}}, \code{\link{lcmcross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}.
#'
#' @details \code{\link{marginal}} operates in the presence of exogenous 
#' variables that explain inefficiency, namely the inefficiency drivers
#'  (\eqn{uhet = ~ Z_u} or \eqn{muhet = ~ Z_{mu}}).
#'
#' Two components are computed for each variable: the marginal effects on the
#' expected inefficiency (\eqn{\frac{\partial E[u]}{\partial Z_{(m)u}}}) and
#' the marginal effects on the variance of inefficiency (\eqn{\frac{\partial
#' V[u]}{\partial Z_{(m)u}}}).
#'
#' The model also allows the Wang (2002) parametrization of \eqn{\mu} and
#' \eqn{\sigma_u^2} by the same vector of exogenous variables. This double
#' parameterization accounts for non-monotonic relationships between the
#' inefficiency and its drivers.
#'
#' @param object A stochastic frontier model returned
#' by \code{\link{sfacross}}, \code{\link{lcmcross}}, 
#' \code{\link{sfaselectioncross}} or \code{\link{zisfcross}}.
#' @param newData Optional data frame that is used to calculate the marginal 
#' effect of \eqn{Z} variables on inefficiency. If NULL (the default), the 
#' marginal estimates are calculated for the observations that were used in the 
#' estimation.
#' @param ... Currently ignored.
#'
#' @name marginal
#'
#' @return \code{\link{marginal}} returns a data frame containing the marginal
#' effects of the \eqn{Z_u} variables on the expected inefficiency (each
#' variable has the prefix \code{'Eu_'}) and on the variance of the
#' inefficiency (each variable has the prefix \code{'Vu_'}) is returned.
#'
#' In the case of the latent class frontier model (LCM) or the zero inefficiency
#' frontier model, each variable terminates with
#' \code{'_c#'} where \code{'#'} is the class number
#' .
# @author K Hervé Dakpo
#'
#' @seealso \code{\link{sfacross}}, for the stochastic frontier analysis model
#' fitting function.
#'
#' \code{\link{lcmcross}}, for the latent class stochastic frontier analysis
#' model fitting function.
#' 
#' \code{\link{sfaselectioncross}} for sample selection in stochastic frontier 
#' model fitting function.
#' 
#' \code{\link{zisfcross}} for zero inefficiency in stochastic frontier model
#' fitting function.
#'
#' @references Wang, H.J. 2002. Heteroscedasticity and non-monotonic efficiency
#' effects of a stochastic frontier model. \emph{Journal of Productivity
#' Analysis}, \bold{18}:241--253.
#'
#' @keywords methods marginal
#'
#' @examples
#'
#' ## Using data on fossil fuel fired steam electric power generation plants in
#' the U.S.
#' # Translog SFA (cost function) truncated normal with scaling property
#' tl_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + I(1/2 * (log(y))^2) +
#' log(wl/wf) + log(wk/wf) + I(1/2 * (log(wl/wf))^2) + I(1/2 * (log(wk/wf))^2) +
#' I(log(wl/wf) * log(wk/wf)) + I(log(y) * log(wl/wf)) + I(log(y) * log(wk/wf)),
#' udist = 'tnormal', muhet = ~ regu + wl, uhet = ~ regu + wl, data = utility, 
#' S = -1, scaling = TRUE, method = 'mla')
#' marg.tl_u_ts <- marginal(tl_u_ts)
#' summary(marg.tl_u_ts)
#'
#' ## Using data on eighty-two countries production (DGP)
#' # LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- lcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#'     data = worldprod, uhet = ~ initStat + h, S = 1, method = 'mla')
#'   marg.cb_2c_h <- marginal(cb_2c_h)
#'   summary(marg.cb_2c_h)
#'
#' @aliases marginal.sfacross
#' @export
#' @export marginal
# marginal effects computation sfacross ----------
marginal.sfacross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$udist == "hnormal") {
    if (object$nuZUvar == 1) {
      stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
        call. = FALSE)
    } else {
      EffMarg <- bind_cols(as_tibble(cmarghalfnorm_Eu(object = object)),
        as_tibble(cmarghalfnorm_Vu(object = object)))
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- bind_cols(as_tibble(cmargexponorm_Eu(object = object)),
          as_tibble(cmargexponorm_Vu(object = object)))
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- bind_cols(as_tibble(cmarggammanorm_Eu(object = object)),
          as_tibble(cmarggammanorm_Vu(object = object)))
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- bind_cols(as_tibble(cmargraynorm_Eu(object = object)),
            as_tibble(cmargraynorm_Vu(object = object)))
          }
        } else {
          if (object$udist == "uniform") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(cmarguninorm_Eu(object = object)),
            as_tibble(cmarguninorm_Vu(object = object)))
          }
          } else {
          if (object$udist == "tnormal") {
            if (object$scaling) {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmargtruncnormscal_Eu(object = object)),
              as_tibble(cmargtruncnormscal_Vu(object = object)))
            }
            } else {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmargtruncnorm_Eu(object = object)),
              as_tibble(cmargtruncnorm_Vu(object = object)))
            }
            }
          } else {
            if (object$udist == "lognormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmarglognorm_Eu(object = object)),
              as_tibble(cmarglognorm_Vu(object = object)))
            }
            } else {
            if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(cmarggenexponorm_Eu(object = object)),
                as_tibble(cmarguninorm_Vu(object = object)))
              }
            } else {
              if (object$udist == "tslaplace") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(cmargtslnorm_Eu(object = object)),
                as_tibble(cmargtslnorm_Vu(object = object)))
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(data.frame(EffMarg))
}

# marginal effects computation lcmcross ----------
#' @rdname marginal
#' @aliases marginal.lcmcross
#' @export
marginal.lcmcross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$nuZUvar == 1) {
    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
      call. = FALSE)
  } else {
    if (object$nClasses == 2) {
      EffMarg <- bind_cols(as_tibble(cmargLCM2Chalfnorm_Eu(object = object)),
        as_tibble(cmargLCM2Chalfnorm_Vu(object = object)))
    } else {
      if (object$nClasses == 3) {
        EffMarg <- bind_cols(as_tibble(cmargLCM3Chalfnorm_Eu(object = object)),
          as_tibble(cmargLCM3Chalfnorm_Vu(object = object)))
      } else {
        if (object$nClasses == 4) {
          EffMarg <- bind_cols(as_tibble(cmargLCM4Chalfnorm_Eu(object = object)),
          as_tibble(cmargLCM4Chalfnorm_Vu(object = object)))
        } else {
          if (object$nClasses == 5) {
          EffMarg <- bind_cols(as_tibble(cmargLCM5Chalfnorm_Eu(object = object)),
            as_tibble(cmargLCM5Chalfnorm_Vu(object = object)))
          }
        }
      }
    }
  }
  return(data.frame(EffMarg))
}

# marginal effects computation sfaselectioncross ----------
#' @rdname marginal
#' @aliases marginal.sfaselectioncross
#' @export
marginal.sfaselectioncross <- function(object, newData = NULL,
  ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$nuZUvar == 1) {
    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
      call. = FALSE)
  } else {
    EffMarg <- bind_cols(as_tibble(cmarghalfnorm_Eu_ss(object = object)),
      as_tibble(cmarghalfnorm_Vu_ss(object = object)))
  }
  return(data.frame(EffMarg))
}

# marginal effects computation zisfcross ----------
#' @rdname marginal
#' @aliases marginal.zisfcross
#' @export
marginal.zisfcross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$linkF == "logit") {
    if (object$sigmavType == "common") {
      if (object$udist == "hnormal") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- bind_cols(as_tibble(czisfmarghalfnorm_Eu_logit(object = object)),
          as_tibble(czisfmarghalfnorm_Vu_logit(object = object)))
        }
      } else {
        if (object$udist == "exponential") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- bind_cols(as_tibble(czisfmargexponorm_Eu_logit(object = object)),
            as_tibble(czisfmargexponorm_Vu_logit(object = object)))
          }
        } else {
          if (object$udist == "gamma") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(czisfmarggammanorm_Eu_logit(object = object)),
            as_tibble(czisfmarggammanorm_Vu_logit(object = object)))
          }
          } else {
          if (object$udist == "rayleigh") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(czisfmargraynorm_Eu_logit(object = object)),
              as_tibble(czisfmargraynorm_Vu_logit(object = object)))
            }
          } else {
            if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(czisfmarguninorm_Eu_logit(object = object)),
              as_tibble(czisfmarguninorm_Vu_logit(object = object)))
            }
            } else {
            if (object$udist == "tnormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(czisfmargtruncnorm_Eu_logit(object = object)),
                as_tibble(czisfmargtruncnorm_Vu_logit(object = object)))
              }
            } else {
              if (object$udist == "lognormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(czisfmarglognorm_Eu_logit(object = object)),
                as_tibble(czisfmarglognorm_Vu_logit(object = object)))
              }
              } else {
              if (object$udist == "genexponential") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(czisfmarggenexponorm_Eu_logit(object = object)),
                  as_tibble(czisfmarggenexponorm_Vu_logit(object = object)))
                }
              } else {
                if (object$udist == "tslaplace") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(czisfmargtslnorm_Eu_logit(object = object)),
                  as_tibble(czisfmargtslnorm_Vu_logit(object = object)))
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    } else {
      if (object$sigmavType == "different") {
        if (object$udist == "hnormal") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- bind_cols(as_tibble(cmnsfmarghalfnorm_Eu_logit(object = object)),
            as_tibble(cmnsfmarghalfnorm_Vu_logit(object = object)))
          }
        } else {
          if (object$udist == "exponential") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(cmnsfmargexponorm_Eu_logit(object = object)),
            as_tibble(cmnsfmargexponorm_Vu_logit(object = object)))
          }
          } else {
          if (object$udist == "gamma") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(cmnsfmarggammanorm_Eu_logit(object = object)),
              as_tibble(cmnsfmarggammanorm_Vu_logit(object = object)))
            }
          } else {
            if (object$udist == "rayleigh") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmargraynorm_Eu_logit(object = object)),
              as_tibble(cmnsfmargraynorm_Vu_logit(object = object)))
            }
            } else {
            if (object$udist == "uniform") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmarguninorm_Eu_logit(object = object)),
                as_tibble(cmnsfmarguninorm_Vu_logit(object = object)))
              }
            } else {
              if (object$udist == "tnormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmargtruncnorm_Eu_logit(object = object)),
                as_tibble(cmnsfmargtruncnorm_Vu_logit(object = object)))
              }
              } else {
              if (object$udist == "lognormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmarglognorm_Eu_logit(object = object)),
                  as_tibble(cmnsfmarglognorm_Vu_logit(object = object)))
                }
              } else {
                if (object$udist == "genexponential") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmarggenexponorm_Eu_logit(object = object)),
                  as_tibble(cmnsfmarggenexponorm_Vu_logit(object = object)))
                }
                } else {
                if (object$udist == "tslaplace") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmargtslnorm_Eu_logit(object = object)),
                    as_tibble(cmnsfmargtslnorm_Vu_logit(object = object)))
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  } else {
    if (object$linkF == "cauchit") {
      if (object$sigmavType == "common") {
        if (object$udist == "hnormal") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- bind_cols(as_tibble(czisfmarghalfnorm_Eu_cauchit(object = object)),
            as_tibble(czisfmarghalfnorm_Vu_cauchit(object = object)))
          }
        } else {
          if (object$udist == "exponential") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(czisfmargexponorm_Eu_cauchit(object = object)),
            as_tibble(czisfmargexponorm_Vu_cauchit(object = object)))
          }
          } else {
          if (object$udist == "gamma") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(czisfmarggammanorm_Eu_cauchit(object = object)),
              as_tibble(czisfmarggammanorm_Vu_cauchit(object = object)))
            }
          } else {
            if (object$udist == "rayleigh") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(czisfmargraynorm_Eu_cauchit(object = object)),
              as_tibble(czisfmargraynorm_Vu_cauchit(object = object)))
            }
            } else {
            if (object$udist == "uniform") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(czisfmarguninorm_Eu_cauchit(object = object)),
                as_tibble(czisfmarguninorm_Vu_cauchit(object = object)))
              }
            } else {
              if (object$udist == "tnormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(czisfmargtruncnorm_Eu_cauchit(object = object)),
                as_tibble(czisfmargtruncnorm_Vu_cauchit(object = object)))
              }
              } else {
              if (object$udist == "lognormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(czisfmarglognorm_Eu_cauchit(object = object)),
                  as_tibble(czisfmarglognorm_Vu_cauchit(object = object)))
                }
              } else {
                if (object$udist == "genexponential") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(czisfmarggenexponorm_Eu_cauchit(object = object)),
                  as_tibble(czisfmarggenexponorm_Vu_cauchit(object = object)))
                }
                } else {
                if (object$udist == "tslaplace") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(czisfmargtslnorm_Eu_cauchit(object = object)),
                    as_tibble(czisfmargtslnorm_Vu_cauchit(object = object)))
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      } else {
        if (object$sigmavType == "different") {
          if (object$udist == "hnormal") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(cmnsfmarghalfnorm_Eu_cauchit(object = object)),
            as_tibble(cmnsfmarghalfnorm_Vu_cauchit(object = object)))
          }
          } else {
          if (object$udist == "exponential") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(cmnsfmargexponorm_Eu_cauchit(object = object)),
              as_tibble(cmnsfmargexponorm_Vu_cauchit(object = object)))
            }
          } else {
            if (object$udist == "gamma") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmarggammanorm_Eu_cauchit(object = object)),
              as_tibble(cmnsfmarggammanorm_Vu_cauchit(object = object)))
            }
            } else {
            if (object$udist == "rayleigh") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmargraynorm_Eu_cauchit(object = object)),
                as_tibble(cmnsfmargraynorm_Vu_cauchit(object = object)))
              }
            } else {
              if (object$udist == "uniform") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmarguninorm_Eu_cauchit(object = object)),
                as_tibble(cmnsfmarguninorm_Vu_cauchit(object = object)))
              }
              } else {
              if (object$udist == "tnormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmargtruncnorm_Eu_cauchit(object = object)),
                  as_tibble(cmnsfmargtruncnorm_Vu_cauchit(object = object)))
                }
              } else {
                if (object$udist == "lognormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmarglognorm_Eu_cauchit(object = object)),
                  as_tibble(cmnsfmarglognorm_Vu_cauchit(object = object)))
                }
                } else {
                if (object$udist == "genexponential") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmarggenexponorm_Eu_cauchit(object = object)),
                    as_tibble(cmnsfmarggenexponorm_Vu_cauchit(object = object)))
                  }
                } else {
                  if (object$udist == "tslaplace") {
                  if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                    EffMarg <- bind_cols(as_tibble(cmnsfmargtslnorm_Eu_cauchit(object = object)),
                    as_tibble(cmnsfmargtslnorm_Vu_cauchit(object = object)))
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    } else {
      if (object$linkF == "probit") {
        if (object$sigmavType == "common") {
          if (object$udist == "hnormal") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- bind_cols(as_tibble(czisfmarghalfnorm_Eu_probit(object = object)),
            as_tibble(czisfmarghalfnorm_Vu_probit(object = object)))
          }
          } else {
          if (object$udist == "exponential") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(czisfmargexponorm_Eu_probit(object = object)),
              as_tibble(czisfmargexponorm_Vu_probit(object = object)))
            }
          } else {
            if (object$udist == "gamma") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(czisfmarggammanorm_Eu_probit(object = object)),
              as_tibble(czisfmarggammanorm_Vu_probit(object = object)))
            }
            } else {
            if (object$udist == "rayleigh") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(czisfmargraynorm_Eu_probit(object = object)),
                as_tibble(czisfmargraynorm_Vu_probit(object = object)))
              }
            } else {
              if (object$udist == "uniform") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(czisfmarguninorm_Eu_probit(object = object)),
                as_tibble(czisfmarguninorm_Vu_probit(object = object)))
              }
              } else {
              if (object$udist == "tnormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(czisfmargtruncnorm_Eu_probit(object = object)),
                  as_tibble(czisfmargtruncnorm_Vu_probit(object = object)))
                }
              } else {
                if (object$udist == "lognormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(czisfmarglognorm_Eu_probit(object = object)),
                  as_tibble(czisfmarglognorm_Vu_probit(object = object)))
                }
                } else {
                if (object$udist == "genexponential") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(czisfmarggenexponorm_Eu_probit(object = object)),
                    as_tibble(czisfmarggenexponorm_Vu_probit(object = object)))
                  }
                } else {
                  if (object$udist == "tslaplace") {
                  if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                    EffMarg <- bind_cols(as_tibble(czisfmargtslnorm_Eu_probit(object = object)),
                    as_tibble(czisfmargtslnorm_Vu_probit(object = object)))
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        } else {
          if (object$sigmavType == "different") {
          if (object$udist == "hnormal") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(cmnsfmarghalfnorm_Eu_probit(object = object)),
              as_tibble(cmnsfmarghalfnorm_Vu_probit(object = object)))
            }
          } else {
            if (object$udist == "exponential") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmargexponorm_Eu_probit(object = object)),
              as_tibble(cmnsfmargexponorm_Vu_probit(object = object)))
            }
            } else {
            if (object$udist == "gamma") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmarggammanorm_Eu_probit(object = object)),
                as_tibble(cmnsfmarggammanorm_Vu_probit(object = object)))
              }
            } else {
              if (object$udist == "rayleigh") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmargraynorm_Eu_probit(object = object)),
                as_tibble(cmnsfmargraynorm_Vu_probit(object = object)))
              }
              } else {
              if (object$udist == "uniform") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmarguninorm_Eu_probit(object = object)),
                  as_tibble(cmnsfmarguninorm_Vu_probit(object = object)))
                }
              } else {
                if (object$udist == "tnormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmargtruncnorm_Eu_probit(object = object)),
                  as_tibble(cmnsfmargtruncnorm_Vu_probit(object = object)))
                }
                } else {
                if (object$udist == "lognormal") {
                  if (object$nmuZUvar == 1 &
                  object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmarglognorm_Eu_probit(object = object)),
                    as_tibble(cmnsfmarglognorm_Vu_probit(object = object)))
                  }
                } else {
                  if (object$udist == "genexponential") {
                  if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                    EffMarg <- bind_cols(as_tibble(cmnsfmarggenexponorm_Eu_probit(object = object)),
                    as_tibble(cmnsfmarggenexponorm_Vu_probit(object = object)))
                  }
                  } else {
                  if (object$udist == "tslaplace") {
                    if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                      call. = FALSE)
                    } else {
                    EffMarg <- bind_cols(as_tibble(cmnsfmargtslnorm_Eu_probit(object = object)),
                      as_tibble(cmnsfmargtslnorm_Vu_probit(object = object)))
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      } else {
        if (object$linkF == "cloglog") {
          if (object$sigmavType == "common") {
          if (object$udist == "hnormal") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- bind_cols(as_tibble(czisfmarghalfnorm_Eu_cloglog(object = object)),
              as_tibble(czisfmarghalfnorm_Vu_cloglog(object = object)))
            }
          } else {
            if (object$udist == "exponential") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(czisfmargexponorm_Eu_cloglog(object = object)),
              as_tibble(czisfmargexponorm_Vu_cloglog(object = object)))
            }
            } else {
            if (object$udist == "gamma") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(czisfmarggammanorm_Eu_cloglog(object = object)),
                as_tibble(czisfmarggammanorm_Vu_cloglog(object = object)))
              }
            } else {
              if (object$udist == "rayleigh") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(czisfmargraynorm_Eu_cloglog(object = object)),
                as_tibble(czisfmargraynorm_Vu_cloglog(object = object)))
              }
              } else {
              if (object$udist == "uniform") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(czisfmarguninorm_Eu_cloglog(object = object)),
                  as_tibble(czisfmarguninorm_Vu_cloglog(object = object)))
                }
              } else {
                if (object$udist == "tnormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(czisfmargtruncnorm_Eu_cloglog(object = object)),
                  as_tibble(czisfmargtruncnorm_Vu_cloglog(object = object)))
                }
                } else {
                if (object$udist == "lognormal") {
                  if (object$nmuZUvar == 1 &
                  object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(czisfmarglognorm_Eu_cloglog(object = object)),
                    as_tibble(czisfmarglognorm_Vu_cloglog(object = object)))
                  }
                } else {
                  if (object$udist == "genexponential") {
                  if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                    EffMarg <- bind_cols(as_tibble(czisfmarggenexponorm_Eu_cloglog(object = object)),
                    as_tibble(czisfmarggenexponorm_Vu_cloglog(object = object)))
                  }
                  } else {
                  if (object$udist == "tslaplace") {
                    if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                      call. = FALSE)
                    } else {
                    EffMarg <- bind_cols(as_tibble(czisfmargtslnorm_Eu_cloglog(object = object)),
                      as_tibble(czisfmargtslnorm_Vu_cloglog(object = object)))
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          } else {
          if (object$sigmavType == "different") {
            if (object$udist == "hnormal") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmarghalfnorm_Eu_cloglog(object = object)),
              as_tibble(cmnsfmarghalfnorm_Vu_cloglog(object = object)))
            }
            } else {
            if (object$udist == "exponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- bind_cols(as_tibble(cmnsfmargexponorm_Eu_cloglog(object = object)),
                as_tibble(cmnsfmargexponorm_Vu_cloglog(object = object)))
              }
            } else {
              if (object$udist == "gamma") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmarggammanorm_Eu_cloglog(object = object)),
                as_tibble(cmnsfmarggammanorm_Vu_cloglog(object = object)))
              }
              } else {
              if (object$udist == "rayleigh") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- bind_cols(as_tibble(cmnsfmargraynorm_Eu_cloglog(object = object)),
                  as_tibble(cmnsfmargraynorm_Vu_cloglog(object = object)))
                }
              } else {
                if (object$udist == "uniform") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmarguninorm_Eu_cloglog(object = object)),
                  as_tibble(cmnsfmarguninorm_Vu_cloglog(object = object)))
                }
                } else {
                if (object$udist == "tnormal") {
                  if (object$nmuZUvar == 1 &
                  object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- bind_cols(as_tibble(cmnsfmargtruncnorm_Eu_cloglog(object = object)),
                    as_tibble(cmnsfmargtruncnorm_Vu_cloglog(object = object)))
                  }
                } else {
                  if (object$udist == "lognormal") {
                  if (object$nmuZUvar == 1 &
                    object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                    EffMarg <- bind_cols(as_tibble(cmnsfmarglognorm_Eu_cloglog(object = object)),
                    as_tibble(cmnsfmarglognorm_Vu_cloglog(object = object)))
                  }
                  } else {
                  if (object$udist == "genexponential") {
                    if (object$nuZUvar == 1) {
                    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                      call. = FALSE)
                    } else {
                    EffMarg <- bind_cols(as_tibble(cmnsfmarggenexponorm_Eu_cloglog(object = object)),
                      as_tibble(cmnsfmarggenexponorm_Vu_cloglog(object = object)))
                    }
                  } else {
                    if (object$udist == "tslaplace") {
                    if (object$nuZUvar ==
                      1) {
                      stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                      call. = FALSE)
                    } else {
                      EffMarg <- bind_cols(as_tibble(cmnsfmargtslnorm_Eu_cloglog(object = object)),
                      as_tibble(cmnsfmargtslnorm_Vu_cloglog(object = object)))
                    }
                    }
                  }
                  }
                }
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  return(data.frame(EffMarg))
}
