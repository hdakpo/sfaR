# variance covariance matrix for sfacross ----------

vcov.sfacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1])) {
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE
    )
  }
  resCov <- object$invHessian
  if (extraPar) {
    if (object$udist %in% c("tnormal", "lognormal")) {
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    } else {
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nuZUvar +
        1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
    }
    Wu <- mean(as.numeric(crossprod(matrix(delta), t(uHvar))))
    Wv <- mean(as.numeric(crossprod(matrix(phi), t(vHvar))))
    if (object$udist %in% c("tnormal", "lognormal")) {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar > 1) {
        stop("argument 'extraPar' is not available for heteroscedasctic models", call. = FALSE)
      }
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1) {
       stop("argument 'extraPar' is not available for heteroscedasctic models", call. = FALSE)
      }
    }
    jac <- diag(nrow(resCov))
    jac <- rbind(jac, matrix(0, nrow = 9, ncol = ncol(resCov)))
    rownames(jac) <- c(
      rownames(resCov), "sigmaSq", "lambdaSq",
      "sigmauSq", "sigmavSq", "sigma", "lambda", "sigmau",
      "sigmav", "gamma"
    )
    colnames(jac) <- colnames(resCov)
    jac["sigmaSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmaSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["lambdaSq", "Zu_(Intercept)"] <- exp(Wu) / exp(Wv)
    jac["lambdaSq", "Zv_(Intercept)"] <- -exp(Wu + Wv) / exp(2 *
      Wv)
    jac["sigmauSq", "Zu_(Intercept)"] <- exp(Wu)
    jac["sigmavSq", "Zv_(Intercept)"] <- exp(Wv)
    jac["sigma", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu) * (exp(Wu) +
      exp(Wv))^(-1 / 2)
    jac["sigma", "Zv_(Intercept)"] <- 1 / 2 * exp(Wv) * (exp(Wu) +
      exp(Wv))^(-1 / 2)
    jac["lambda", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu / 2) / exp(Wv / 2)
    jac["lambda", "Zv_(Intercept)"] <- -1 / 2 * exp(Wu / 2 +
      Wv / 2) / exp(Wv)
    jac["sigmau", "Zu_(Intercept)"] <- 1 / 2 * exp(Wu / 2)
    jac["sigmav", "Zv_(Intercept)"] <- 1 / 2 * exp(Wv / 2)
    jac["gamma", "Zu_(Intercept)"] <- (exp(Wu) * (exp(Wu) +
      exp(Wv)) - exp(2 * Wu)) / (exp(Wu) + exp(Wv))^2
    jac["gamma", "Zv_(Intercept)"] <- -exp(Wu + Wv) / (exp(Wu) +
      exp(Wv))^2
    resCov <- jac %*% resCov %*% t(jac)
  }
  return(resCov)
}


# variance covariance matrix for lcmcross ----------

vcov.lcmcross <- function(object, ...) {
  resCov <- object$invHessian
  return(resCov)
}
