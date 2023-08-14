################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Marginal impact of Z variables on inefficiency                               #
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

#' Marginal effects of the inefficiency drivers in stochastic frontier models
#'
#' This function returns marginal effects of the inefficiency drivers from 
#' stochastic frontier models estimated with \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, or 
#' \code{\link{sfaselectioncross}}.
#'
#' @details \code{\link{marginal}} operates in the presence of exogenous 
#' variables that explain inefficiency, namely the inefficiency drivers
#'  (\eqn{uhet = ~ Z_u} or \eqn{muhet = ~ Z_{mu}}).
#'
#' Two components are computed for each variable: the marginal effects on the
#' expected inefficiency (\eqn{\frac{\partial E[u]}{\partial Z_{mu}}}) and
#' the marginal effects on the variance of inefficiency (\eqn{\frac{\partial
#' V[u]}{\partial Z_{mu}}}).
#'
#' The model also allows the Wang (2002) parametrization of \eqn{\mu} and
#' \eqn{\sigma_u^2} by the same vector of exogenous variables. This double
#' parameterization accounts for non-monotonic relationships between the
#' inefficiency and its drivers.
#'
#' @param object A stochastic frontier model returned by \code{\link{sfacross}}, 
#' \code{\link{sfalcmcross}}, \code{\link{sfagzisfcross}}, 
#' \code{\link{sfacnsfcross}}, \code{\link{sfamisfcross}}, 
#' \code{\link{sfazisfcross}}, \code{\link{sfametacross}}, or 
#' \code{\link{sfaselectioncross}}.
#' @param newData Optional data frame that is used to calculate the marginal 
#' effect of \eqn{Z} variables on inefficiency. If NULL (the default), the 
#' marginal estimates are calculated for the observations that were used in the 
#' estimation.
#' @param ... Currently ignored.
#'
#' @name marginal
#'
#' @aliases marginal.sfacross marginal.sfalcmcross marginal.sfagzisfcross 
#' marginal.sfacnsfcross marginal.sfamisfcross marginal.sfazisfcross 
#' marginal.sfaselectioncross marginal.sfametacross
#'
#' @return \code{\link{marginal}} returns a data frame containing the marginal
#' effects of the \eqn{Z_u} variables on the expected inefficiency (each
#' variable has the prefix \code{'Eu_'}) and on the variance of the
#' inefficiency (each variable has the prefix \code{'Vu_'}).
#'
#' In the case of the contaminated noise stochastic frontier (CNSF), or 
#' generalized zero inefficiency stochastic frontier (GZISF), or latent class 
#' stochastic frontier (LCM), or the multi-modal inefficiency stochastic 
#' frontier (MISF) or the zero inefficiency frontier (ZISF), 
#' each variable terminates with \code{'_c#'} where \code{'#'} is the class 
#' number.
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
#' @references Wang, H.J. 2002. Heteroscedasticity and non-monotonic efficiency
#' effects of a stochastic frontier model. \emph{Journal of Productivity
#' Analysis}, \bold{18}:241--253.
#'
#' @keywords methods marginal
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
#' udist = 'tnormal', muhet = ~ regu + wl, uhet = ~ regu + wl, data = utility, 
#' S = -1, scaling = TRUE, method = 'mla')
#' marg.tl_u_ts <- marginal(tl_u_ts)
#' summary(marg.tl_u_ts)
#'
#' # Using data on eighty-two countries production (GDP)
#' ## LCM Cobb Douglas (production function) half normal distribution
#' cb_2c_h <- sfalcmcross(formula = ly ~ lk + ll + yr, udist = 'hnormal',
#'     data = worldprod, uhet = ~ initStat + h, S = 1, method = 'mla')
#'   marg.cb_2c_h <- marginal(cb_2c_h)
#'   summary(marg.cb_2c_h)
#'   }
#' 
#' @export
#' @export marginal
# @exportS3Method marginal sfacross
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
      EffMarg <- data.frame(cbind(cmarghalfnorm_Eu(object = object),
        cmarghalfnorm_Vu(object = object)))
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- data.frame(cbind(cmargexponorm_Eu(object = object),
          cmargexponorm_Vu(object = object)))
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- data.frame(cbind(cmarggammanorm_Eu(object = object),
          cmarggammanorm_Vu(object = object)))
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- data.frame(cbind(cmargraynorm_Eu(object = object),
            cmargraynorm_Vu(object = object)))
          }
        } else {
          if (object$udist == "uniform") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- data.frame(cbind(cmarguninorm_Eu(object = object),
            cmarguninorm_Vu(object = object)))
          }
          } else {
          if (object$udist == "tnormal") {
            if (object$scaling) {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnormscal_Eu(object = object),
              cmargtruncnormscal_Vu(object = object)))
            }
            } else {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnorm_Eu(object = object),
              cmargtruncnorm_Vu(object = object)))
            }
            }
          } else {
            if (object$udist == "lognormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmarglognorm_Eu(object = object),
              cmarglognorm_Vu(object = object)))
            }
            } else {
            if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- data.frame(cbind(cmarggenexponorm_Eu(object = object),
                cmarguninorm_Vu(object = object)))
              }
            } else {
              if (object$udist == "tslaplace") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- data.frame(cbind(cmargtslnorm_Eu(object = object),
                cmargtslnorm_Vu(object = object)))
              }
              } else {
              if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- data.frame(cbind(cmargweibullnorm_Eu(object = object),
                  cmargweibullnorm_Vu(object = object)))
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
  return(EffMarg)
}

# marginal effects computation sfalcmcross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfalcmcross
marginal.sfalcmcross <- function(object, newData = NULL, ...) {
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
      EffMarg <- data.frame(cbind(cmargLCM2Chalfnorm_Eu(object = object),
        cmargLCM2Chalfnorm_Vu(object = object)))
    } else {
      if (object$nClasses == 3) {
        EffMarg <- data.frame(cbind(cmargLCM3Chalfnorm_Eu(object = object),
          cmargLCM3Chalfnorm_Vu(object = object)))
      } else {
        if (object$nClasses == 4) {
          EffMarg <- data.frame(cbind(cmargLCM4Chalfnorm_Eu(object = object),
          cmargLCM4Chalfnorm_Vu(object = object)))
        } else {
          if (object$nClasses == 5) {
          EffMarg <- data.frame(cbind(cmargLCM5Chalfnorm_Eu(object = object),
            cmargLCM5Chalfnorm_Vu(object = object)))
          }
        }
      }
    }
  }
  return(EffMarg)
}

# marginal effects computation sfagzisfcross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfagzisfcross
marginal.sfagzisfcross <- function(...) {
  marginal.sfalcmcross(...)
}

# marginal effects computation sfacnsfcross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfacnsfcross
marginal.sfacnsfcross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$sigmauType == "common") {
    if (object$udist == "hnormal") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmarghalfnorm_Eu_",
          object$linkF, "(object = object), ccnsfmarghalfnorm_Vu_",
          object$linkF, "(object = object)))")))
      }
    } else {
      if (object$udist == "exponential") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmargexponorm_Eu_",
          object$linkF, "(object = object), ccnsfmargexponorm_Vu_",
          object$linkF, "(object = object)))")))
        }
      } else {
        if (object$udist == "gamma") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmarggammanorm_Eu_",
            object$linkF, "(object = object), ccnsfmarggammanorm_Vu_",
            object$linkF, "(object = object)))")))
          }
        } else {
          if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmargraynorm_Eu_",
            object$linkF, "(object = object), ccnsfmargraynorm_Vu_",
            object$linkF, "(object = object)))")))
          }
          } else {
          if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmarguninorm_Eu_",
              object$linkF, "(object = object), ccnsfmarguninorm_Vu_",
              object$linkF, "(object = object)))")))
            }
          } else {
            if (object$udist == "tnormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmargtruncnorm_Eu_",
              object$linkF, "(object = object), ccnsfmargtruncnorm_Vu_",
              object$linkF, "(object = object)))")))
            }
            } else {
            if (object$udist == "lognormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmarglognorm_Eu_",
                object$linkF, "(object = object), ccnsfmarglognorm_Vu_",
                object$linkF, "(object = object)))")))
              }
            } else {
              if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmarggenexponorm_Eu_",
                object$linkF, "(object = object), ccnsfmarggenexponorm_Vu_",
                object$linkF, "(object = object)))")))
              }
              } else {
              if (object$udist == "tslaplace") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmargtslnorm_Eu_",
                  object$linkF, "(object = object), ccnsfmargtslnorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
              } else {
                if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(ccnsfmargweibullnorm_Eu_",
                  object$linkF, "(object = object), ccnsfmargweibullnorm_Vu_",
                  object$linkF, "(object = object)))")))
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
    if (object$sigmauType == "different") {
      if (object$udist == "hnormal") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmarghalfnorm_Eu_",
          object$linkF, "(object = object), cmcesfmarghalfnorm_Vu_",
          object$linkF, "(object = object)))")))
        }
      } else {
        if (object$udist == "exponential") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmargexponorm_Eu_",
            object$linkF, "(object = object), cmcesfmargexponorm_Vu_",
            object$linkF, "(object = object)))")))
          }
        } else {
          if (object$udist == "gamma") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmarggammanorm_Eu_",
            object$linkF, "(object = object), cmcesfmarggammanorm_Vu_",
            object$linkF, "(object = object)))")))
          }
          } else {
          if (object$udist == "rayleigh") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmargraynorm_Eu_",
              object$linkF, "(object = object), cmcesfmargraynorm_Vu_",
              object$linkF, "(object = object)))")))
            }
          } else {
            if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmarguninorm_Eu_",
              object$linkF, "(object = object), cmcesfmarguninorm_Vu_",
              object$linkF, "(object = object)))")))
            }
            } else {
            if (object$udist == "tnormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmargtruncnorm_Eu_",
                object$linkF, "(object = object), cmcesfmargtruncnorm_Vu_",
                object$linkF, "(object = object)))")))
              }
            } else {
              if (object$udist == "lognormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmarglognorm_Eu_",
                object$linkF, "(object = object), cmcesfmarglognorm_Vu_",
                object$linkF, "(object = object)))")))
              }
              } else {
              if (object$udist == "genexponential") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmarggenexponorm_Eu_",
                  object$linkF, "(object = object), cmcesfmarggenexponorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
              } else {
                if (object$udist == "tslaplace") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmargtslnorm_Eu_",
                  object$linkF, "(object = object), cmcesfmargtslnorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
                } else {
                if (object$udist == "weibull") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmcesfmargweibullnorm_Eu_",
                    object$linkF, "(object = object), cmcesfmargweibullnorm_Vu_",
                    object$linkF, "(object = object)))")))
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

# marginal effects computation sfamisfcross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfamisfcross
marginal.sfamisfcross <- function(object, newData = NULL, ...) {
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
      EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmarghalfnorm_Eu_",
        object$linkF, "(object = object), cmisfmarghalfnorm_Vu_",
        object$linkF, "(object = object)))")))
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmargexponorm_Eu_",
          object$linkF, "(object = object), cmisfmargexponorm_Vu_",
          object$linkF, "(object = object)))")))
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmarggammanorm_Eu_",
          object$linkF, "(object = object), cmisfmarggammanorm_Vu_",
          object$linkF, "(object = object)))")))
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmargraynorm_Eu_",
            object$linkF, "(object = object), cmisfmargraynorm_Vu_",
            object$linkF, "(object = object)))")))
          }
        } else {
          if (object$udist == "uniform") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmarguninorm_Eu_",
            object$linkF, "(object = object), cmisfmarguninorm_Vu_",
            object$linkF, "(object = object)))")))
          }
          } else {
          if (object$udist == "tnormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
            1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmargtruncnorm_Eu_",
              object$linkF, "(object = object), cmisfmargtruncnorm_Vu_",
              object$linkF, "(object = object)))")))
            }
          } else {
            if (object$udist == "lognormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmarglognorm_Eu_",
              object$linkF, "(object = object), cmisfmarglognorm_Vu_",
              object$linkF, "(object = object)))")))
            }
            } else {
            if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmarggenexponorm_Eu_",
                object$linkF, "(object = object), cmisfmarggenexponorm_Vu_",
                object$linkF, "(object = object)))")))
              }
            } else {
              if (object$udist == "tslaplace") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmargtslnorm_Eu_",
                object$linkF, "(object = object), cmisfmargtslnorm_Vu_",
                object$linkF, "(object = object)))")))
              }
              } else {
              if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmisfmargweibullnorm_Eu_",
                  object$linkF, "(object = object), cmisfmargweibullnorm_Vu_",
                  object$linkF, "(object = object)))")))
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

# marginal effects computation sfazisfcross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfazisfcross
marginal.sfazisfcross <- function(object, newData = NULL, ...) {
  if (!is.null(newData)) {
    if (!is.data.frame(newData)) {
      stop("argument 'newData' must be of class data.frame")
    }
    object$dataTable <- newData
    object$Nobs <- dim(newData)[1]
  }
  if (object$sigmavType == "common") {
    if (object$udist == "hnormal") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmarghalfnorm_Eu_",
          object$linkF, "(object = object), czisfmarghalfnorm_Vu_",
          object$linkF, "(object = object)))")))
      }
    } else {
      if (object$udist == "exponential") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmargexponorm_Eu_",
          object$linkF, "(object = object), czisfmargexponorm_Vu_",
          object$linkF, "(object = object)))")))
        }
      } else {
        if (object$udist == "gamma") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmarggammanorm_Eu_",
            object$linkF, "(object = object), czisfmarggammanorm_Vu_",
            object$linkF, "(object = object)))")))
          }
        } else {
          if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmargraynorm_Eu_",
            object$linkF, "(object = object), czisfmargraynorm_Vu_",
            object$linkF, "(object = object)))")))
          }
          } else {
          if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmarguninorm_Eu_",
              object$linkF, "(object = object), czisfmarguninorm_Vu_",
              object$linkF, "(object = object)))")))
            }
          } else {
            if (object$udist == "tnormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmargtruncnorm_Eu_",
              object$linkF, "(object = object), czisfmargtruncnorm_Vu_",
              object$linkF, "(object = object)))")))
            }
            } else {
            if (object$udist == "lognormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmarglognorm_Eu_",
                object$linkF, "(object = object), czisfmarglognorm_Vu_",
                object$linkF, "(object = object)))")))
              }
            } else {
              if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmarggenexponorm_Eu_",
                object$linkF, "(object = object), czisfmarggenexponorm_Vu_",
                object$linkF, "(object = object)))")))
              }
              } else {
              if (object$udist == "tslaplace") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmargtslnorm_Eu_",
                  object$linkF, "(object = object), czisfmargtslnorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
              } else {
                if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(czisfmargweibullnorm_Eu_",
                  object$linkF, "(object = object), czisfmargweibullnorm_Vu_",
                  object$linkF, "(object = object)))")))
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
    if (object$sigmavType == "different") {
      if (object$udist == "hnormal") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmarghalfnorm_Eu_",
          object$linkF, "(object = object), cmnsfmarghalfnorm_Vu_",
          object$linkF, "(object = object)))")))
        }
      } else {
        if (object$udist == "exponential") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmargexponorm_Eu_",
            object$linkF, "(object = object), cmnsfmargexponorm_Vu_",
            object$linkF, "(object = object)))")))
          }
        } else {
          if (object$udist == "gamma") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmarggammanorm_Eu_",
            object$linkF, "(object = object), cmnsfmarggammanorm_Vu_",
            object$linkF, "(object = object)))")))
          }
          } else {
          if (object$udist == "rayleigh") {
            if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
            EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmargraynorm_Eu_",
              object$linkF, "(object = object), cmnsfmargraynorm_Vu_",
              object$linkF, "(object = object)))")))
            }
          } else {
            if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmarguninorm_Eu_",
              object$linkF, "(object = object), cmnsfmarguninorm_Vu_",
              object$linkF, "(object = object)))")))
            }
            } else {
            if (object$udist == "tnormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmargtruncnorm_Eu_",
                object$linkF, "(object = object), cmnsfmargtruncnorm_Vu_",
                object$linkF, "(object = object)))")))
              }
            } else {
              if (object$udist == "lognormal") {
              if (object$nmuZUvar == 1 & object$nuZUvar ==
                1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmarglognorm_Eu_",
                object$linkF, "(object = object), cmnsfmarglognorm_Vu_",
                object$linkF, "(object = object)))")))
              }
              } else {
              if (object$udist == "genexponential") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmarggenexponorm_Eu_",
                  object$linkF, "(object = object), cmnsfmarggenexponorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
              } else {
                if (object$udist == "tslaplace") {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmargtslnorm_Eu_",
                  object$linkF, "(object = object), cmnsfmargtslnorm_Vu_",
                  object$linkF, "(object = object)))")))
                }
                } else {
                if (object$udist == "weibull") {
                  if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE)
                  } else {
                  EffMarg <- eval(parse(text = paste0("data.frame(cbind(cmnsfmargweibullnorm_Eu_",
                    object$linkF, "(object = object), cmnsfmargweibullnorm_Vu_",
                    object$linkF, "(object = object)))")))
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

# marginal effects computation sfaselectioncross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfaselectioncross
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
    EffMarg <- data.frame(cbind(cmarghalfnorm_Eu_ss(object = object),
      cmarghalfnorm_Vu_ss(object = object)))
  }
  return(EffMarg)
}

# marginal effects computation sfametacross ----------
#' @rdname marginal
#' @export
# @exportS3Method marginal sfametacross
marginal.sfametacross <- function(object, newData = NULL, ...) {
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
      EffMarg <- data.frame(cbind(cmarghalfnorm_Eu(object = object),
        cmarghalfnorm_Vu(object = object)))
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
      } else {
        EffMarg <- data.frame(cbind(cmargexponorm_Eu(object = object),
          cmargexponorm_Vu(object = object)))
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE)
        } else {
          EffMarg <- data.frame(cbind(cmarggammanorm_Eu(object = object),
          cmarggammanorm_Vu(object = object)))
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
          EffMarg <- data.frame(cbind(cmargraynorm_Eu(object = object),
            cmargraynorm_Vu(object = object)))
          }
        } else {
          if (object$udist == "uniform") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE)
          } else {
            EffMarg <- data.frame(cbind(cmarguninorm_Eu(object = object),
            cmarguninorm_Vu(object = object)))
          }
          } else {
          if (object$udist == "tnormal") {
            if (object$scaling) {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnormscal_Eu(object = object),
              cmargtruncnormscal_Vu(object = object)))
            }
            } else {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmargtruncnorm_Eu(object = object),
              cmargtruncnorm_Vu(object = object)))
            }
            }
          } else {
            if (object$udist == "lognormal") {
            if (object$nmuZUvar == 1 & object$nuZUvar ==
              1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE)
            } else {
              EffMarg <- data.frame(cbind(cmarglognorm_Eu(object = object),
              cmarglognorm_Vu(object = object)))
            }
            } else {
            if (object$udist == "genexponential") {
              if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
              EffMarg <- data.frame(cbind(cmarggenexponorm_Eu(object = object),
                cmarguninorm_Vu(object = object)))
              }
            } else {
              if (object$udist == "tslaplace") {
              if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE)
              } else {
                EffMarg <- data.frame(cbind(cmargtslnorm_Eu(object = object),
                cmargtslnorm_Vu(object = object)))
              }
              } else {
              if (object$udist == "weibull") {
                if (object$nuZUvar == 1) {
                stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                  call. = FALSE)
                } else {
                EffMarg <- data.frame(cbind(cmargweibullnorm_Eu(object = object),
                  cmargweibullnorm_Vu(object = object)))
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
  return(EffMarg)
}
