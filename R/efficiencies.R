# conditional efficiencies sfacross ----------

efficiencies.sfacross <- function(object, level = 0.95, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (object$udist == "hnormal") {
    EffRes <- chalfnormeff(object = object, level = level)
  } else {
    if (object$udist == "exponential") {
      EffRes <- cexponormeff(object = object, level = level)
    } else {
      if (object$udist == "gamma") {
        EffRes <- cgammanormeff(object = object, level = level)
      } else {
        if (object$udist == "rayleigh") {
          EffRes <- craynormeff(object = object, level = level)
        } else {
          if (object$udist == "uniform") {
            EffRes <- cuninormeff(object = object, level = level)
          } else {
            if (object$udist == "tnormal") {
              if (object$scaling) {
                EffRes <- ctruncnormscaleff(
                  object = object,
                  level = level
                )
              } else {
                EffRes <- ctruncnormeff(
                  object = object,
                  level = level
                )
              }
            } else {
              if (object$udist == "lognormal") {
                EffRes <- clognormeff(
                  object = object,
                  level = level
                )
              } else {
                if (object$udist == "genexponential") {
                  EffRes <- cgenexponormeff(
                    object = object,
                    level = level
                  )
                } else {
                  if (object$udist == "tslaplace") {
                    EffRes <- ctslnormeff(
                      object = object,
                      level = level
                    )
                  } else {
                    if (object$udist == "weibull") {
                      EffRes <- cweibullnormeff(
                        object = object,
                        level = level
                      )
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
return(data.frame(EffRes))
}

# conditional efficiencies lcmcross ----------

efficiencies.lcmcross <- function(object, level = 0.95, ...) {
  if (level < 0 || level > 0.9999) {
    stop("'level' must be between 0 and 0.9999", call. = FALSE)
  }
  if (object$nClasses == 2) {
    EffRes <- cLCM2Chalfnormeff(object = object, level = level)
  } else {
    if (object$nClasses == 3) {
      EffRes <- cLCM3Chalfnormeff(object = object, level = level)
    } else {
      if (object$nClasses == 4) {
        EffRes <- cLCM4Chalfnormeff(
          object = object,
          level = level
        )
      } else {
        if (object$nClasses == 5) {
          EffRes <- cLCM5Chalfnormeff(
            object = object,
            level = level
          )
        }
      }
    }
  }
return(data.frame(EffRes))
}
