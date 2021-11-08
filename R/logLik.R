# log likelihood extraction for sfacross ----------

logLik.sfacross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value", call. = FALSE)
  if (individual) {
    LL <- list()
    LL[["logLik"]] <- object$dataTable$logL_OBS
    LL[["Nobs"]] <- object$Nobs
    LL[["df"]] <- object$nParm
  } else {
    LL <- rbind("logLik: " = object$mlLoglik, "Nobs: " = object$Nobs, "df: " = object$nParm)
  }
  return(LL)
}

# log likelihood extraction for lcmcross ----------

logLik.lcmcross <- function(object, individual = FALSE, ...) {
  if (length(individual) != 1 || !is.logical(individual[1]))
    stop("argument 'individual' must be a single logical value", call. = FALSE)
  if (individual) {
    LL <- list()
    LL[["logLik"]] <- object$dataTable$logL_OBS
    LL[["Nobs"]] <- object$Nobs
    LL[["df"]] <- object$nParm
  } else {
    LL <- rbind("logLik: " = object$mlLoglik, "Nobs: " = object$Nobs, "df: " = object$nParm)
  }
  return(LL)
}
