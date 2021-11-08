# coefficients from sfacross ----

coef.sfacross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value",
         call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    if (object$udist == "tnormal") {
      if (object$scaling) {
        beta <- object$mlParam[1:(object$nXvar)]
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
          (object$nuZUvar - 1))]
        tau <- object$mlParam[object$nXvar + (object$nuZUvar -
          1) + 1]
        cu <- object$mlParam[object$nXvar + (object$nuZUvar -
          1) + 2]
        phi <- object$mlParam[(object$nXvar + (object$nuZUvar -
          1) + 2 + 1):(object$nXvar + (object$nuZUvar -
          1) + 2 + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
        Wu <- cu + 2 * as.numeric(crossprod(matrix(delta),
          t(uHvar[, -1])))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
          object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
          object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      }
    } else {
      if (object$udist == "lognormal") {
        delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
          1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
          object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
          object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
        vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 4)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      } else {
        delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
          object$nuZUvar)]
        phi <- object$mlParam[(object$nXvar + object$nuZUvar +
          1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
        uHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 2)
        vHvar <- model.matrix(object$formula, data = object$dataTable,
          rhs = 3)
        Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
        Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
      }
    }
    if (object$udist == "lognormal" || (object$udist == "tnormal" &
      object$scaling == FALSE)) {
      if (object$nuZUvar > 1 || object$nvZVvar > 1 || object$nmuZUvar >
        1)
        cat("Variances averaged over observations     \n\n")
    } else {
      if (object$nuZUvar > 1 || object$nvZVvar > 1)
        cat("Variances averaged over observations     \n\n")
    }
    cRes <- c(cRes, sigmaSq = mean(exp(Wu)) + mean(exp(Wv)),
      lambdaSq = mean(exp(Wu))/mean(exp(Wv)), sigmauSq = mean(exp(Wu)),
      sigmavSq = mean(exp(Wv)), sigma = sqrt(mean(exp(Wu)) +
        mean(exp(Wv))), lambda = sqrt(mean(exp(Wu))/mean(exp(Wv))),
      sigmau = sqrt(mean(exp(Wu))), sigmav = sqrt(mean(exp(Wv))),
      gamma = mean(exp(Wu))/(mean(exp(Wu)) + mean(exp(Wv))))
  }
  return(cRes)
}

# coefficients from lcmcross ----

coef.lcmcross <- function(object, extraPar = FALSE, ...) {
  if (length(extraPar) != 1 || !is.logical(extraPar[1]))
    stop("argument 'extraPar' must be a single logical value",
      call. = FALSE)
  cRes <- object$mlParam
  if (extraPar) {
    uHvar <- model.matrix(object$formula, data = object$dataTable,
      rhs = 2)
    vHvar <- model.matrix(object$formula, data = object$dataTable,
      rhs = 3)
    delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
      object$nuZUvar)]
    phi1 <- object$mlParam[(object$nXvar + object$nuZUvar +
      1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
    delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
      object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
      object$nvZVvar)]
    phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar +
      object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar +
      2 * object$nvZVvar)]
    Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
    Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
    Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
    Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
    if (object$nClasses == 3) {
      delta3 <- object$mlParam[(3 * object$nXvar + 2 *
        object$nuZUvar + 2 * object$nvZVvar + 1):(3 * object$nXvar +
        3 * object$nuZUvar + 2 * object$nvZVvar)]
      phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar +
        2 * object$nvZVvar + 1):(3 * object$nXvar + 3 *
        object$nuZUvar + 3 * object$nvZVvar)]
      Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
      Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
    } else {
      if (object$nClasses == 4) {
        delta4 <- object$mlParam[(4 * object$nXvar +
          3 * object$nuZUvar + 3 * object$nvZVvar + 1):(4 *
          object$nXvar + 4 * object$nuZUvar + 3 * object$nvZVvar)]
        phi4 <- object$mlParam[(4 * object$nXvar + 4 *
          object$nuZUvar + 3 * object$nvZVvar + 1):(4 *
          object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar)]
        Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
        Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
      } else {
        if (object$nClasses == 5) {
          delta5 <- object$mlParam[(5 * object$nXvar +
          4 * object$nuZUvar + 4 * object$nvZVvar + 1):(5 *
          object$nXvar + 5 * object$nuZUvar + 4 * object$nvZVvar)]
          phi5 <- object$mlParam[(5 * object$nXvar +
          5 * object$nuZUvar + 4 * object$nvZVvar + 1):(5 *
          object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar)]
          Wu5 <- as.numeric(crossprod(matrix(delta5),
          t(uHvar)))
          Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
        }
      }
    }
    if (object$nuZUvar > 1 || object$nvZVvar > 1)
      cat("Variances averaged over observations     \n\n")
    cRes <- c(cRes, sigmaSq1 = mean(exp(Wu1)) + mean(exp(Wv1)),
      lambdaSq1 = mean(exp(Wu1))/mean(exp(Wv1)), sigmauSq1 = mean(exp(Wu1)),
      sigmavSq1 = mean(exp(Wv1)), sigma1 = sqrt(mean(exp(Wu1)) +
        mean(exp(Wv1))), lambda1 = sqrt(mean(exp(Wu1))/mean(exp(Wv1))),
      sigmau1 = sqrt(mean(exp(Wu1))), sigmav1 = sqrt(mean(exp(Wv1))),
      gamma1 = mean(exp(Wu1))/(mean(exp(Wu1)) + mean(exp(Wv1))),
      sigmaSq2 = mean(exp(Wu2)) + mean(exp(Wv2)), lambdaSq2 = mean(exp(Wu2))/mean(exp(Wv2)),
      sigmauSq2 = mean(exp(Wu2)), sigmavSq2 = mean(exp(Wv2)),
      sigma2 = sqrt(mean(exp(Wu2)) + mean(exp(Wv2))), lambda2 = sqrt(mean(exp(Wu2))/mean(exp(Wv2))),
      sigmau2 = sqrt(mean(exp(Wu2))), sigmav2 = sqrt(mean(exp(Wv2))),
      gamma2 = mean(exp(Wu2))/(mean(exp(Wu2)) + mean(exp(Wv2))))
    if (object$nClasses == 3) {
      cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)),
        lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)), sigmauSq3 = mean(exp(Wu3)),
        sigmavSq3 = mean(exp(Wv3)), sigma3 = sqrt(mean(exp(Wu3)) +
          mean(exp(Wv3))), lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
        sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
        gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))))
    } else {
      if (object$nClasses == 4) {
        cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) + mean(exp(Wv3)),
          lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)),
          sigma3 = sqrt(mean(exp(Wu3)) + mean(exp(Wv3))),
          lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) + mean(exp(Wv3))),
          sigmaSq4 = mean(exp(Wu4)) + mean(exp(Wv4)),
          lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)),
          sigmauSq4 = mean(exp(Wu4)), sigmavSq4 = mean(exp(Wv4)),
          sigma4 = sqrt(mean(exp(Wu4)) + mean(exp(Wv4))),
          lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))),
          sigmau4 = sqrt(mean(exp(Wu4))), sigmav4 = sqrt(mean(exp(Wv4))),
          gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) + mean(exp(Wv4))))
      } else {
        if (object$nClasses == 5) {
          cRes <- c(cRes, sigmaSq3 = mean(exp(Wu3)) +
          mean(exp(Wv3)), lambdaSq3 = mean(exp(Wu3))/mean(exp(Wv3)),
          sigmauSq3 = mean(exp(Wu3)), sigmavSq3 = mean(exp(Wv3)),
          sigma3 = sqrt(mean(exp(Wu3)) + mean(exp(Wv3))),
          lambda3 = sqrt(mean(exp(Wu3))/mean(exp(Wv3))),
          sigmau3 = sqrt(mean(exp(Wu3))), sigmav3 = sqrt(mean(exp(Wv3))),
          gamma3 = mean(exp(Wu3))/(mean(exp(Wu3)) +
            mean(exp(Wv3))), sigmaSq4 = mean(exp(Wu4)) +
            mean(exp(Wv4)), lambdaSq4 = mean(exp(Wu4))/mean(exp(Wv4)),
          sigmauSq4 = mean(exp(Wu4)), sigmavSq4 = mean(exp(Wv4)),
          sigma4 = sqrt(mean(exp(Wu4)) + mean(exp(Wv4))),
          lambda4 = sqrt(mean(exp(Wu4))/mean(exp(Wv4))),
          sigmau4 = sqrt(mean(exp(Wu4))), sigmav4 = sqrt(mean(exp(Wv4))),
          gamma4 = mean(exp(Wu4))/(mean(exp(Wu4)) +
            mean(exp(Wv4))), sigmaSq5 = mean(exp(Wu5)) +
            mean(exp(Wv5)), lambdaSq5 = mean(exp(Wu5))/mean(exp(Wv5)),
          sigmauSq5 = mean(exp(Wu5)), sigmavSq5 = mean(exp(Wv5)),
          sigma5 = sqrt(mean(exp(Wu5)) + mean(exp(Wv5))),
          lambda5 = sqrt(mean(exp(Wu5))/mean(exp(Wv5))),
          sigmau5 = sqrt(mean(exp(Wu5))), sigmav5 = sqrt(mean(exp(Wv5))),
          gamma5 = mean(exp(Wu5))/(mean(exp(Wu5)) +
            mean(exp(Wv5))))
        }
      }
    }
  }
  return(cRes)
}
