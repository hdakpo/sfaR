################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gamma-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in seq_len(N)) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P <= 0) {
    return(NA)
  }
  ll <- -1/2 * P * Wu - log(gamma(P)) + exp(Wv)/(2 * exp(Wu)) +
    S * epsilon/sqrt(exp(Wu)) + pnorm(-S * epsilon/sqrt(exp(Wv)) -
    sqrt(exp(Wv)/exp(Wu)), log.p = TRUE) + log(Hi)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for gamma-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
cstgammanorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/2)))^(2/3)
  } else {
    varu <- (-S * m3/2)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1/2 * log((epsiRes^2 - varv)^2)
  dep_v <- 1/2 * log((epsiRes^2 - varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetu$coefficients))) {
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  }
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetv$coefficients))) {
    stop("at least one of the OLS coefficients of 'vhet' is NA: ",
      paste(colnames(vHvar)[is.na(reg_hetv$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  }
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi, P = 1))
}

# Gradient of the likelihood function ----------
#' gradient for gamma-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cgradgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wuwv <- sqrt(exp(Wv)/exp(Wu))
  dwuwv <- dnorm(-(S * (epsilon)/exp(Wv/2) + wuwv))
  pwuwv <- pnorm(-(S * (epsilon)/exp(Wv/2) + wuwv))
  depsi <- dnorm((exp(Wv)/exp(Wu/2) + S * (epsilon))/exp(Wv/2))
  pepsi <- pnorm((exp(Wv)/exp(Wu/2) + S * (epsilon))/exp(Wv/2))
  sigx1 <- (exp(Wv)/exp(Wu/2) + S * (epsilon))
  sigx2 <- (dwuwv/(exp(Wv/2) * pwuwv) - 1/exp(Wu/2))
  sigx3 <- (dwuwv/(exp(Wu) * pwuwv * wuwv))
  sigx4 <- (0.5 * sigx3 - 2 * (exp(Wu)/(2 * exp(Wu))^2)) *
    exp(Wv)
  sigx5 <- (sigx4 - (0.5 * (S * (epsilon)/exp(Wu/2)) + 0.5 *
    P))
  sigx6 <- (exp(Wv)/(2 * exp(Wu)) - (0.5 * (exp(Wv)/(exp(Wu) *
    wuwv)) - 0.5 * (S * (epsilon)/exp(Wv/2))) * dwuwv/pwuwv)
  F1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pepsi, FUN = "*") +
    FiMat
  dqF1 <- dnorm(qnorm(F1))
  F2 <- sweep(qnorm(F1), MARGIN = 1, STATS = exp(Wv/2), FUN = "*")
  F3 <- sweep(F2, MARGIN = 1, STATS = sigx1, FUN = "-")
  sumF3 <- apply(F3^(P - 1), 1, sum)
  F4 <- sweep((1 - FiMat)/dqF1, MARGIN = 1, STATS = depsi,
    FUN = "*")
  F5 <- sweep((0.5 - 0.5 * (F4)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = exp(Wv)/exp(Wu/2), FUN = "*")
  F6 <- sweep(F4, MARGIN = 1, STATS = exp(Wv)/exp(Wu/2) - 0.5 *
    sigx1, FUN = "*") + 0.5 * (F2)
  F7 <- sweep(F6, MARGIN = 1, STATS = exp(Wv)/exp(Wu/2), FUN = "-") *
    F3^(P - 2) * (P - 1)
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    gx[, k] <- apply(sweep(S * (1 - F4) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)/sumF3
  }
  gx <- sweep(Xvar, MARGIN = 1, STATS = S * sigx2, FUN = "*") +
    gx
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    gu[, k] <- apply(sweep(F5, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF3
  }
  gu <- sweep(uHvar, MARGIN = 1, STATS = sigx5, FUN = "*") +
    gu
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    gv[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF3
  }
  gv <- sweep(vHvar, MARGIN = 1, STATS = sigx6, FUN = "*") +
    gv
  gradll <- cbind(gx, gu, gv, (apply((F3)^(P - 1) * log(F3),
    1, sum)/sumF3 - (0.5 * (Wu) + digamma(P))))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for gamma-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
chessgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  ewv_h <- exp(Wv/2)
  ewu_h <- exp(Wu/2)
  wuwv <- sqrt(ewv/ewu)
  dwuwv <- dnorm(-(S * (epsilon)/ewv_h + wuwv))
  pwuwv <- pnorm(-(S * (epsilon)/ewv_h + wuwv))
  depsi <- dnorm((ewv/ewu_h + S * (epsilon))/ewv_h)
  pepsi <- pnorm((ewv/ewu_h + S * (epsilon))/ewv_h)
  sigx1 <- (ewv/ewu_h + S * (epsilon))
  Q <- dim(FiMat)[2]
  F1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pepsi, FUN = "*") +
    FiMat
  qF1 <- qnorm(F1)
  dqF1 <- dnorm(qF1)
  F2 <- sweep(qF1, MARGIN = 1, STATS = ewv_h, FUN = "*")
  F3 <- sweep(F2, MARGIN = 1, STATS = sigx1, FUN = "-")
  sumF3 <- apply(F3^(P - 1), 1, sum)
  F4 <- sweep((1 - FiMat)/dqF1, MARGIN = 1, STATS = depsi,
    FUN = "*")
  F5 <- sweep((0.5 - 0.5 * (F4)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = ewv/ewu_h, FUN = "*")
  F6 <- sweep(F4, MARGIN = 1, STATS = ewv/ewu_h - 0.5 * sigx1,
    FUN = "*") + 0.5 * (F2)
  F7 <- sweep(F6, MARGIN = 1, STATS = ewv/ewu_h, FUN = "-") *
    F3^(P - 2) * (P - 1)
  F8 <- sweep((1 - FiMat) * dqF1 * qF1/dqF1^2, MARGIN = 1,
    STATS = depsi, FUN = "*")
  F9 <- sweep((-F8), MARGIN = 1, STATS = sigx1/ewv_h, FUN = "+")
  F10 <- sweep((1 - FiMat) * (F3)^(P - 2)/(dqF1), MARGIN = 1,
    STATS = depsi/ewv_h, FUN = "*")
  F11 <- ((F3)^(P - 2) + (F3)^(P - 2) * log(F3) * (P - 1))
  F12 <- sweep((F6), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  F13 <- sweep(F9 * (1 - FiMat) * (F3)^(P - 2)/(dqF1), MARGIN = 1,
    STATS = depsi * (ewv/ewu_h - 0.5 * sigx1)/ewv_h, FUN = "*")
  F14 <- sweep((-0.5 * (F8)), MARGIN = 1, STATS = 0.5 * (sigx1/ewv_h),
    FUN = "+")
  F15 <- sweep(((0.5 - 0.5 * F4)^2 * (F3)^(P - 3) * (P - 2) -
    0.5 * (F14 * F10)), MARGIN = 1, STATS = ewv/ewu_h, FUN = "*")
  F16 <- sweep(F14, MARGIN = 1, STATS = (ewv/ewu_h - 0.5 *
    sigx1)/ewv_h, FUN = "*")
  F17 <- sweep((F8), MARGIN = 1, STATS = sigx1/ewv_h, FUN = "-")
  F18 <- sweep(F17, MARGIN = 1, STATS = (ewv/ewu_h - 0.5 *
    sigx1)^2/ewv_h, FUN = "*")
  F19 <- sweep((F18), MARGIN = 1, STATS = 0.5 * (ewv/ewu_h),
    FUN = "+")
  F20 <- sweep(0.5 * (F6), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  sigx7 <- (S * (epsilon)/ewv_h + wuwv)
  sumF3H <- apply((F3)^(P - 1) * log(F3), 1, sum)
  X1 <- matrix(nrow = N, ncol = nXvar)
  X2 <- matrix(nrow = N, ncol = nXvar)
  X3 <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    X1[, k] <- apply(sweep(S * (1 - F4) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)/sumF3
    X2[, k] <- apply(sweep(S * (1 - F4) * (F3)^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum)/sumF3
    X3[, k] <- apply(sweep(S * F11 * (1 - F4), MARGIN = 1,
      STATS = Xvar[, k], FUN = "*"), 1, sum)/sumF3
  }
  ZU1 <- matrix(nrow = N, ncol = nuZUvar)
  ZU2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    ZU1[, k] <- apply(sweep(F5, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF3
    ZU2[, k] <- apply(sweep(F11 * (0.5 - 0.5 * F4), MARGIN = 1,
      STATS = uHvar[, k] * ewv/ewu_h, FUN = "*"), 1, sum)/sumF3
  }
  ZV1 <- matrix(nrow = N, ncol = nvZVvar)
  ZV2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    ZV1[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF3
    ZV2[, k] <- apply(sweep(F12 * F11, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF3
  }
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in seq_len(Q)) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((1 - F4)^2 * (F3)^(P - 3) * (P - 2) - F9 * F10)[,
        r] * (P - 1)/sumF3, FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = S *
      wHvar * ((0.5 - 0.5 * F4) * (1 - F4) * (F3)^(P -
      3) * (P - 2) - 0.5 * (F9 * F10))[, r] * ewv * (P -
      1)/ewu_h/sumF3, FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = S *
      wHvar * (F12 * (1 - F4) * (F3)^(P - 3) * (P - 2) +
      F13)[, r] * (P - 1)/sumF3, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (F15 - 0.5 * ((0.5 - 0.5 * F4) * (F3)^(P - 2)))[,
        r] * ewv * (P - 1)/ewu_h/sumF3, FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (((F16 - 0.5) * F4 + 0.5) * (F3)^(P - 2) + F12 *
        (0.5 - 0.5 * F4) * (F3)^(P - 3) * (P - 2))[,
        r] * ewv * (P - 1)/ewu_h/sumF3, FUN = "*"), vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      ((F19 * F4 + F20) * (F3)^(P - 2) + F12^2 * (F3)^(P -
        3) * (P - 2))[, r] * (P - 1)/sumF3, FUN = "*"),
      vHvar)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) - crossprod(sweep(X2,
    MARGIN = 1, STATS = wHvar/sumF3^2, FUN = "*"), X2) +
    crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * (sigx7/(ewv_h^2 *
      pwuwv) - dwuwv/(ewv_h * pwuwv)^2) * dwuwv, FUN = "*"),
      Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) - crossprod(sweep(X2, MARGIN = 1, STATS = wHvar/sumF3,
    FUN = "*"), ZU1) + crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S * wHvar * (0.5 * ((sigx7/(ewu * pwuwv * wuwv) -
      dwuwv * ewu * wuwv/(ewu * pwuwv * wuwv)^2) * dwuwv *
      ewv/ewv_h) + 0.5/ewu_h), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) - crossprod(sweep(X2,
    MARGIN = 1, STATS = wHvar/sumF3, FUN = "*"), ZV1) - crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * wHvar * ((0.5 * (ewv/(ewu * wuwv)) -
      0.5 * (S * (epsilon)/ewv_h)) * (S * (epsilon)/ewv_h +
      wuwv - dwuwv/pwuwv) + 0.5) * dwuwv/(ewv_h * pwuwv),
    FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(X3,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(X1, MARGIN = 1,
    STATS = wHvar * sumF3H/sumF3, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) - crossprod(sweep(ZU1,
    MARGIN = 1, STATS = wHvar/sumF3, FUN = "*"), ZU1) + crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * ((0.5 * (sigx7/(ewu *
      pwuwv)) - ((0.5 * (dwuwv * ewv/wuwv) + ewu * pwuwv) *
      wuwv - 0.5 * (ewv * pwuwv/wuwv))/(ewu * pwuwv * wuwv)^2) *
      dwuwv) - 2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu/(2 *
      ewu)^2)) * ewv + 0.25 * (S * (epsilon)/ewu_h)), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) -
    crossprod(sweep(ZU1, MARGIN = 1, STATS = wHvar/sumF3,
      FUN = "*"), ZV1) - crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (((0.5 * (ewv/(ewu * wuwv)) - 0.5 * (S *
      (epsilon)/ewv_h)) * (0.5 * sigx7 - 0.5 * (dwuwv/pwuwv))/(ewu *
      wuwv) - 0.5 * ((ewu * wuwv - 0.5 * (ewv/wuwv))/(ewu *
      wuwv)^2)) * dwuwv/pwuwv + 2 * (ewu/(2 * ewu)^2)) *
      ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(ZU2, MARGIN = 1, STATS = wHvar, FUN = "*") -
    sweep(ZU1, MARGIN = 1, STATS = wHvar * sumF3H/sumF3,
      FUN = "*") - 0.5 * sweep(uHvar, MARGIN = 1, STATS = wHvar,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) - crossprod(sweep(ZV1, MARGIN = 1, STATS = wHvar/sumF3,
    FUN = "*"), ZV1) + crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (ewv/(2 * ewu) - ((0.5 * (ewv/(ewu *
      wuwv)) - 0.5 * (S * (epsilon)/ewv_h))^2 * (dwuwv/pwuwv -
      sigx7) + 0.25 * (S * (epsilon)/ewv_h) + 0.5 * ((1/ewu -
      0.5 * (ewv/(ewu * wuwv)^2)) * ewv/wuwv)) * dwuwv/pwuwv),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(ZV2,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(ZV1, MARGIN = 1,
    STATS = wHvar * sumF3H/sumF3, FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(((apply((F3)^(P - 1) * log(F3)^2, 1, sum) -
    sumF3H^2/sumF3)/sumF3 - trigamma(P)) * wHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for gamma-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
gammanormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar,
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType,
  qac) {
  startVal <- if (!is.null(start))
    start else cstgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cgammanormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar)), hessian = 0,
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cgammanormlike, grad = cgradgammanormlike,
      hess = chessgammanormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      hs = function(parm) as(-chessgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      hess = function(parm) -chessgammanormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar)), gradient = function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar)), hessian = function(parm) -chessgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar), control = list(iter.max = itermax,
      trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradgammanormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessgammanormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)
    if (method == "sr1")
      mleObj$hessian <- chessgammanormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)
  }
  mleObj$logL_OBS <- cgammanormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar)
  mleObj$gradL_OBS <- cgradgammanormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for gamma-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
cgammanormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in seq_len(object$Nobs)) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  u <- Hi1/Hi2
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in seq_len(object$Nobs)) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
        1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
        ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
        1))
    }
    teBC <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
      exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
      epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_reciprocal <- exp(-exp(Wv)/exp(Wu/2) - object$S *
      epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
      object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
      Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for gamma-normal distribution
#' @param object object of class sfacross
#' @noRd
cmarggammanorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmarggammanorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
