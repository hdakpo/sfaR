# information criteria for sfacross ----------

ic.sfacross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
    } else {
      if (IC == "HQIC") {
       obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) * object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep="")
}

# information criteria for lcmcross ----------

ic.lcmcross <- function(object, IC = "AIC", ...) {
  if (!(IC %in% c("AIC", "BIC", "HQIC"))) {
    stop("Unknown information criteria: ", paste(IC), call. = FALSE)
  }
  if (IC == "AIC") {
    obj <- -2 * object$mlLoglik + 2 * object$nParm
  } else {
    if (IC == "BIC") {
      obj <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
    } else {
      if (IC == "HQIC") {
        obj <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) * object$nParm
      }
    }
  }
  message(IC, ": ", prettyNum(obj), sep="")
}
