# marginal effects computation sfacross ----------

marginal.sfacross <- function(object, ...) {
  if (object$udist == "hnormal") {
    if (object$nuZUvar == 1) {
      stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
        call. = FALSE
      )
    } else {
      EffMarg <- bind_cols(
        as_tibble(cmarghalfnorm_Eu(object = object)),
        as_tibble(cmarghalfnorm_Vu(object = object))
      )
    }
  } else {
    if (object$udist == "exponential") {
      if (object$nuZUvar == 1) {
        stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
          call. = FALSE
        )
      } else {
        EffMarg <- bind_cols(
          as_tibble(cmargexponorm_Eu(object = object)),
          as_tibble(cmargexponorm_Vu(object = object))
        )
      }
    } else {
      if (object$udist == "gamma") {
        if (object$nuZUvar == 1) {
          stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
            call. = FALSE
          )
        } else {
          EffMarg <- bind_cols(
            as_tibble(cmarggammanorm_Eu(object = object)),
            as_tibble(cmarggammanorm_Vu(object = object))
          )
        }
      } else {
        if (object$udist == "rayleigh") {
          if (object$nuZUvar == 1) {
            stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
              call. = FALSE
            )
          } else {
            EffMarg <- bind_cols(
              as_tibble(cmargraynorm_Eu(object = object)),
              as_tibble(cmargraynorm_Vu(object = object))
            )
          }
        } else {
          if (object$udist == "uniform") {
            if (object$nuZUvar == 1) {
              stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                call. = FALSE
              )
            } else {
              EffMarg <- bind_cols(
                as_tibble(cmarguninorm_Eu(object = object)),
                as_tibble(cmarguninorm_Vu(object = object))
              )
            }
          } else {
            if (object$udist == "tnormal") {
              if (object$scaling) {
                if (object$nuZUvar == 1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmargtruncnormscal_Eu(object = object)),
                    as_tibble(cmargtruncnormscal_Vu(object = object))
                  )
                }
              } else {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmargtruncnorm_Eu(object = object)),
                    as_tibble(cmargtruncnorm_Vu(object = object))
                  )
                }
              }
            } else {
              if (object$udist == "lognormal") {
                if (object$nmuZUvar == 1 & object$nuZUvar ==
                  1) {
                  stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
                    call. = FALSE
                  )
                } else {
                  EffMarg <- bind_cols(
                    as_tibble(cmarglognorm_Eu(object = object)),
                    as_tibble(cmarglognorm_Vu(object = object))
                  )
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

marginal.lcmcross <- function(object, ...) {
  if (object$nuZUvar == 1) {
    stop("Marginal effects can only be computed from models with exogenous variables that explain inefficiency",
      call. = FALSE
    )
  } else {
    if (object$nClasses == 2) {
      EffMarg <- bind_cols(
        as_tibble(cmargLCM2Chalfnorm_Eu(object = object)),
        as_tibble(cmargLCM2Chalfnorm_Vu(object = object))
      )
    } else {
      if (object$nClasses == 3) {
        EffMarg <- bind_cols(
          as_tibble(cmargLCM3Chalfnorm_Eu(object = object)),
          as_tibble(cmargLCM3Chalfnorm_Vu(object = object))
        )
      } else {
        if (object$nClasses == 4) {
          EffMarg <- bind_cols(
            as_tibble(cmargLCM4Chalfnorm_Eu(object = object)),
            as_tibble(cmargLCM4Chalfnorm_Vu(object = object))
          )
        } else {
          if (object$nClasses == 5) {
            EffMarg <- bind_cols(
              as_tibble(cmargLCM5Chalfnorm_Eu(object = object)),
              as_tibble(cmargLCM5Chalfnorm_Vu(object = object))
            )
          }
        }
      }
    }
  }
  return(data.frame(EffMarg))
}
