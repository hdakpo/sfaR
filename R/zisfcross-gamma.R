################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf gamma-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
czisfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P < 0)
    return(NA)
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf gamma-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
cstzisfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + gamma - normal distributions...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgammanormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = printInfo,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat)
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "P",
    paste0("SF_", colnames(Zvar)))
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Gradient of the likelihood function ----------
#' gradient for zisf gamma-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cgradzisfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewv_h <- exp(Wv/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewu_p <- exp(-(Wu * P/2))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2 <- (ewv/ewu_h + S * (epsilon))
  pmusig <- pnorm(-(ewv_h/ewu_h + S * (epsilon)/ewv_h))
  dmusig <- dnorm(sigx2/ewv_h, 0, 1)
  pepsi <- pnorm(sigx2/ewv_h)
  depsi <- dnorm(-(ewv_h/ewu_h + S * (epsilon)/ewv_h), 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewv_h)
  sigx3 <- (depsi/ewv_h - pmusig/ewu_h)
  sigx4 <- (Q * wzdeno * gamma(P))
  sigx5 <- (depsi * ewv_h/ewu_h)
  sigx6 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewv_h, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx12 <- (prC * depsiv/ewv_h + ewu_p * sigx1 * ewz * pmusig *
    sDiv/sigx4)
  sigx8 <- (Q * sigx12 * wzdeno * gamma(P))
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewv_h/ewu_h) -
    0.5 * (S * (epsilon)/ewv_h)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewv_h^2) - 0.5 *
    depsiv)
  sigx11 <- (1/sigx4 - Q * ewz * gamma(P)/sigx4^2)
  sigx13 <- (ewv/ewu_h - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * wzdeno * digamma(P) * gamma(P) *
    sDiv/sigx4^2)
  sigx16 <- (sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
    depsiv/(wzdeno * ewv_h))
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewv_h),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsiv *
    (epsilon)/ewv_h^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx12, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  gu <- sweep(UF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx8,
    FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prC/ewv_h/sigx12, FUN = "*")
  gradll <- cbind(gx, gu, gv, sigx15 * ewu_p * sigx1 * ewz *
    pmusig/sigx12, sweep(Zvar, MARGIN = 1, STATS = sigx16 *
    ewz/sigx12, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf gamma-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
chesszisfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Q <- dim(FiMat)[2]
  ewv <- exp(Wv)
  ewv_h <- exp(Wv/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewu_p <- exp(-(Wu * P/2))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- exp(ewv/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2 <- (ewv/ewu_h + S * (epsilon))
  musig <- (ewv_h/ewu_h + S * (epsilon)/ewv_h)
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(sigx2/ewv_h, 0, 1)
  pepsi <- pnorm(sigx2/ewv_h)
  depsi <- dnorm(-musig, 0, 1)
  depsiv <- dnorm(S * (epsilon)/ewv_h)
  sigx3 <- (depsi/ewv_h - pmusig/ewu_h)
  sigx4 <- (Q * wzdeno * gamma(P))
  sigx5 <- (depsi * ewv_h/ewu_h)
  sigx6 <- (0.5 * (S * (epsilon)/ewu_h) + 0.5 * P + 2 * (ewu *
    ewv/(2 * ewu)^2))
  sigx7 <- (0.5 * sigx5 - sigx6 * pmusig)
  qFi <- qnorm(sweep((1 - FiMat), MARGIN = 1, STATS = pepsi,
    FUN = "*") + FiMat)
  F1 <- dnorm(qFi)
  F2 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig, FUN = "*")
  F3 <- sweep(sweep((qFi), MARGIN = 1, STATS = ewv_h, FUN = "*"),
    MARGIN = 1, STATS = sigx2, FUN = "-")
  sDiv <- apply(F3^(P - 1), 1, sum)
  sigx12 <- (prC * depsiv/ewv_h + ewu_p * sigx1 * ewz * pmusig *
    sDiv/sigx4)
  sigx8 <- (Q * sigx12 * wzdeno * gamma(P))
  sigx9 <- (ewv * pmusig/(2 * ewu) - (0.5 * (ewv_h/ewu_h) -
    0.5 * (S * (epsilon)/ewv_h)) * depsi)
  sigx10 <- (0.5 * (S^2 * depsiv * (epsilon)^2/ewv_h^2) - 0.5 *
    depsiv)
  sigx11 <- (1/sigx4 - Q * ewz * gamma(P)/sigx4^2)
  sigx13 <- (ewv/ewu_h - 0.5 * sigx2)
  sF3 <- apply(F3^(P - 1) * log(F3), 1, sum)
  sigx14 <- (sF3 - 0.5 * (Wu * sDiv))
  sigx15 <- (sigx14/sigx4 - Q * wzdeno * digamma(P) * gamma(P) *
    sDiv/sigx4^2)
  sigx16 <- (sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
    depsiv/(wzdeno * ewv_h))
  F4 <- sweep((1 - FiMat)/F1, MARGIN = 1, STATS = dmusig *
    sigx13, FUN = "*")
  F5 <- sweep(sweep(qFi, MARGIN = 1, STATS = 0.5 * (ewv_h),
    FUN = "*"), MARGIN = 1, STATS = ewv/ewu_h, FUN = "-")
  XF1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF1[, k] <- apply(sweep(S * (1 - F2) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)
  }
  XF2 <- sweep(XF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = S * sigx3 * sDiv, FUN = "*")
  XF3 <- sweep(XF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4,
    FUN = "*")
  XF4 <- sweep(Xvar, MARGIN = 1, STATS = S^2 * prC * depsiv *
    (epsilon)/ewv_h^3, FUN = "*")
  gx <- sweep((XF3 + XF4), MARGIN = 1, STATS = 1/sigx12, FUN = "*")
  UF1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF1[, k] <- apply(sweep((0.5 - 0.5 * (F2)) * F3^(P -
      2) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h * uHvar[,
      k], FUN = "*"), 1, sum)
  }
  UF2 <- sweep(UF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = sigx7 * sDiv, FUN = "*")
  VF1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF1[, k] <- apply(sweep((F4 + F5) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1,
      sum)
  }
  VF2 <- sweep(VF1, MARGIN = 1, STATS = pmusig, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = sigx9 * sDiv, FUN = "*")
  gv <- sweep(VF2, MARGIN = 1, STATS = ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = sigx10 *
    prC/ewv_h/sigx12, FUN = "*")
  F6 <- sweep(-sweep((1 - FiMat) * F1 * qFi/F1^2, MARGIN = 1,
    STATS = dmusig, FUN = "*"), MARGIN = 1, STATS = sigx2/ewv_h,
    FUN = "+")
  F7 <- sweep(F6 * (1 - FiMat) * F3^(P - 2)/F1, MARGIN = 1,
    STATS = dmusig/ewv_h, FUN = "*")
  XF5 <- sweep(XF1, MARGIN = 1, STATS = sigx3, FUN = "*")
  XF6 <- sweep(Xvar, MARGIN = 1, STATS = S * (musig/ewv_h -
    1/ewu_h) * sDiv, FUN = "*")
  XF7 <- sweep((XF6 + XF1), MARGIN = 1, STATS = depsi/ewv_h,
    FUN = "*")
  XF8 <- sweep(XF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  F8 <- sweep(((0.5 - 0.5 * (F2)) * (1 - F2) * F3^(P - 3) *
    (P - 2) - 0.5 * (F7)) * (P - 1), MARGIN = 1, STATS = ewv/ewu_h,
    FUN = "*")
  UF3 <- sweep(UF1, MARGIN = 1, STATS = depsi/ewv_h, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = 1/ewu_h, FUN = "*")
  F9 <- sweep(F6 * (1 - FiMat) * F3^(P - 2)/F1, MARGIN = 1,
    STATS = dmusig * sigx13/ewv_h, FUN = "*")
  XF9 <- sweep(Xvar, MARGIN = 1, STATS = S * depsi/ewv_h, FUN = "*")
  XF10 <- sweep((XF3 + XF4), MARGIN = 1, STATS = pmusig/sigx12,
    FUN = "*")
  XF11 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    XF11[, k] <- apply(sweep(S * (F3^(P - 2) + F3^(P - 2) *
      log(F3) * (P - 1)) * (1 - F2), MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
  }
  XF12 <- sweep(XF1, MARGIN = 1, STATS = (Wu), FUN = "*")
  F10 <- sweep(((1 - FiMat) * F1 * qFi/F1^2), MARGIN = 1, STATS = dmusig,
    FUN = "*")
  F11 <- sweep((-0.5 * F10), MARGIN = 1, STATS = 0.5 * (sigx2/ewv_h),
    FUN = "+")
  F12 <- sweep((F11 * (1 - FiMat) * F3^(P - 2)/F1), MARGIN = 1,
    STATS = dmusig/ewv_h, FUN = "*")
  F13 <- sweep(((0.5 - 0.5 * (F2))^2 * F3^(P - 3) * (P - 2) -
    0.5 * F12), MARGIN = 1, STATS = ewv/ewu_h, FUN = "*")
  F14 <- sweep((F13 - 0.5 * ((0.5 - 0.5 * (F2)) * F3^(P - 2))),
    MARGIN = 1, STATS = ewv * (P - 1)/ewu_h, FUN = "*")
  F15<-sweep(F10,MARGIN=1,STATS=sigx2/ewv_h,FUN="-")
  F16 <- sweep(sweep(F15, MARGIN = 1, STATS = sigx13^2/ewv_h,
    FUN = "*"), MARGIN = 1, STATS = 0.5 * (ewv/ewu_h), FUN = "+")
  F17 <- 0.5 * (sweep(0.5 * (qFi), MARGIN = 1, STATS = ewv_h,
    FUN = "*") + F4)
  F18 <- sweep((F16 * F2 + F17), MARGIN = 1, STATS = ewv/ewu_h,
    FUN = "-")
  F19 <- sweep(F11, MARGIN = 1, STATS = sigx13/ewv_h, FUN = "*")
  F20 <- sweep((((F19 - 0.5) * F2 + 0.5) * (F3)^(P - 2) + (F4 +
    F5) * (0.5 - 0.5 * F2) * (F3)^(P - 3) * (P - 2)), MARGIN = 1,
    STATS = ewv * (P - 1)/ewu_h, FUN = "*")
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S^2 * ((1 - F2)^2 * F3^(P - 3) * (P - 2) - F7) *
        (P - 1))[, r] * pmusig * ewu_p * sigx1 * ewz/sigx4/sigx12,
      FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S * F8)[, r] * pmusig/sigx8 * ewu_p * sigx1 * ewz,
      FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      (S * ((F4 + F5) * (1 - F2) * F3^(P - 3) * (P - 2) +
        F9) * (P - 1))[, r] * pmusig * ewu_p * sigx1 *
      ewz/sigx4/sigx12, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (F14)[, r] * pmusig/sigx8 * ewu_p * sigx1 * ewz,
      FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      F20[, r] * pmusig * ewu_p * sigx1 * ewz/sigx8, FUN = "*"),
      vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      ((F18 * (F3)^(P - 2) + (F4 + F5)^2 * (F3)^(P - 3) *
        (P - 2)) * (P - 1))[, r] * pmusig * ewu_p * sigx1 *
      ewz/sigx4/sigx12, FUN = "*"), vHvar)
  }
  UF4 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    UF4[, k] <- apply(sweep(((F3)^(P - 2) + (F3)^(P - 2) *
      log(F3) * (P - 1)) * (0.5 - 0.5 * F2), MARGIN = 1,
      STATS = uHvar[, k] * ewv/ewu_h, FUN = "*"), 1, sum)
  }
  VF3 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    VF3[, k] <- apply(sweep((F4 + F5) * ((F3)^(P - 2) + (F3)^(P -
      2) * log(F3) * (P - 1)), MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*"), S * (XF5 + XF7 - XF8)) + crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S^2 * prC * depsiv * (S^2 *
      (epsilon)^2/ewv_h^2 - 1)/ewv_h^3/sigx12, FUN = "*"),
    Xvar) - crossprod(sweep(XF3 + XF4, MARGIN = 1, STATS = wHvar/sigx12^2,
    FUN = "*"), XF3 + XF4)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) + crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ewu_p * sigx1 * ewz/sigx8, FUN = "*"), S * UF3) + crossprod(sweep(XF1,
    MARGIN = 1, STATS = wHvar * sigx7/sigx8 * ewu_p * sigx1 *
      ewz, FUN = "*") + sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (musig/ewu_h) - sigx6/ewv_h) * depsi + 0.5 *
    (pmusig/ewu_h)) * sDiv/sigx8 * ewu_p * sigx1 * ewz, FUN = "*"),
    uHvar) - crossprod(sweep((XF3 + XF4), MARGIN = 1, STATS = wHvar *
    Q * wzdeno * gamma(P)/sigx8^2 * ewu_p * sigx1 * ewz,
    FUN = "*"), UF2)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) + crossprod(S * Xvar,
    sweep(VF1, MARGIN = 1, STATS = wHvar * depsi/ewv_h *
      ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") - sweep(VF2,
      MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/ewu_h/sigx12,
      FUN = "*")) + crossprod(sweep(XF1, MARGIN = 1, STATS = wHvar *
    sigx9 * ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = wHvar * S * depsi * (ewv/(2 *
      ewu) - ((0.5 * (ewv_h/ewu_h) - 0.5 * (S * (epsilon)/ewv_h)) *
      musig + 0.5)) * sDiv/ewv_h * ewu_p * sigx1 * ewz/sigx4/sigx12,
      FUN = "*"), vHvar) + crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (0.5 * (S^2 * (epsilon)^2/ewv_h^2 -
      2) - 0.5) * prC * depsiv * (epsilon)/ewv_h^3/sigx12,
    FUN = "*"), vHvar) - crossprod(sweep((XF3 + XF4), MARGIN = 1,
    STATS = wHvar/sigx12, FUN = "*"), gv)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep((XF9 -
    XF10), MARGIN = 1, STATS = wHvar * sigx15 * ewu_p * sigx1 *
    ewz/sigx12, FUN = "*") + (sweep((XF11 - 0.5 * XF12),
    MARGIN = 1, STATS = wHvar * pmusig/sigx4 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") - (sweep(XF1, MARGIN = 1,
    STATS = wHvar * Q * wzdeno * digamma(P) * gamma(P)/sigx4^2 *
      pmusig * ewu_p * sigx1 * ewz/sigx12, FUN = "*") +
    sweep(Xvar, MARGIN = 1, STATS = wHvar * S * sigx15/ewu_h *
      pmusig * ewu_p * sigx1 * ewz/sigx12, FUN = "*"))))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(XF2,
    MARGIN = 1, STATS = wHvar * sigx11 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") - (sweep(gx, MARGIN = 1, STATS = wHvar *
    sigx16 * ewz/sigx12, FUN = "*") + sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * prC * depsiv * (epsilon)/(wzdeno *
      ewv_h^3) * ewz/sigx12, FUN = "*")), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) + crossprod(uHvar, sweep(UF1,
    MARGIN = 1, STATS = wHvar * (depsi * ewv_h/ewu_h - sigx6 *
      pmusig)/sigx8 * ewu_p * sigx1 * ewz, FUN = "*") +
    sweep(uHvar, MARGIN = 1, STATS = wHvar * ((0.5 * (0.5 *
      (ewv_h * musig/ewu_h) - 0.5) - 0.5 * sigx6) * depsi *
      ewv_h/ewu_h - (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) *
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewu_h)) *
      pmusig) * sDiv/sigx8 * ewu_p * sigx1 * ewz, FUN = "*") -
    sweep(UF2, MARGIN = 1, STATS = wHvar * sigx6/sigx8 *
      ewu_p * sigx1 * ewz, FUN = "*")) - crossprod(sweep(UF2,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx8^2 *
      ewu_p * sigx1 * ewz, FUN = "*"), UF2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) +
    crossprod(uHvar, sweep(VF1, MARGIN = 1, STATS = wHvar *
      0.5 * (depsi * ewv_h/ewu_h) * ewu_p * sigx1 * ewz/sigx8,
      FUN = "*") - sweep(VF2, MARGIN = 1, STATS = wHvar *
      sigx6 * ewu_p * sigx1 * ewz/sigx8, FUN = "*")) +
    crossprod(sweep(UF1, MARGIN = 1, STATS = wHvar * sigx9 *
      ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
      MARGIN = 1, STATS = wHvar * ((depsi * ewv_h/(4 *
        (ewu * ewu_h)) - 2 * (ewu * pmusig/(2 * ewu)^2)) *
        ewv - (0.5 * ((0.5 * (ewv_h/ewu_h) - 0.5 * (S *
        (epsilon)/ewv_h)) * musig) - 0.25) * depsi *
        ewv_h/ewu_h) * sDiv * ewu_p * sigx1 * ewz/sigx8,
      FUN = "*"), vHvar) - crossprod(UF2, sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12 * ewu_p *
      sigx1 * ewz/sigx8, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = wHvar * sigx10 * prC/ewv_h/sigx12 * ewu_p * sigx1 *
      ewz/sigx8, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(UF4, MARGIN = 1, STATS = wHvar *
    pmusig * ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sigx7 * sF3 * ewu_p * sigx1 *
      ewz/sigx8, FUN = "*") - (sweep(UF2, MARGIN = 1, STATS = wHvar *
    0.5 * Wu * ewu_p * sigx1 * ewz/sigx8, FUN = "*") + sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * (pmusig * sDiv) * ewu_p *
      sigx1 * ewz/sigx8, FUN = "*")) - sweep(UF2, MARGIN = 1,
    STATS = wHvar * Q * (sigx12 * digamma(P) + sigx15 * ewu_p *
      sigx1 * ewz * pmusig) * wzdeno * gamma(P)/sigx8^2 *
      ewu_p * sigx1 * ewz, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(UF2,
    sweep(Zvar, MARGIN = 1, STATS = wHvar * (1/sigx4 - (sigx16/sigx8 +
      Q * gamma(P)/sigx4^2) * ewz) * ewu_p * sigx1 * ewz/sigx12,
      FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) + crossprod(vHvar, sweep(VF1, MARGIN = 1, STATS = wHvar *
    (ewv * pmusig/(2 * ewu) - 2 * ((0.5 * (ewv_h/ewu_h) -
      0.5 * (S * (epsilon)/ewv_h)) * depsi)) * ewu_p *
    sigx1 * ewz/sigx4/sigx12, FUN = "*") + sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewv/(2 * ewu) * ewu_p * sigx1 * ewz/sigx4/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (ewv * (pmusig - (0.5 * (ewv_h/ewu_h) - 0.5 * (S * (epsilon)/ewv_h)) *
      depsi)/(2 * ewu) - (0.25 * (ewv_h/ewu_h) + 0.25 *
      (S * (epsilon)/ewv_h) - (0.5 * (ewv_h/ewu_h) - 0.5 *
      (S * (epsilon)/ewv_h))^2 * musig) * depsi) * sDiv *
    ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*")) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewv_h^2) - 1) - 0.25) * depsiv *
      (epsilon)^2/ewv_h^2 - 0.5 * sigx10)/ewv_h/sigx12,
    FUN = "*"), vHvar) - crossprod(sweep(VF2, MARGIN = 1,
    STATS = wHvar * ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") +
    sweep(vHvar, MARGIN = 1, STATS = wHvar * sigx10 * prC/ewv_h/sigx12,
      FUN = "*"), sweep(VF2, MARGIN = 1, STATS = wHvar *
    ewu_p * sigx1 * ewz/sigx4/sigx12, FUN = "*") + sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx10 * prC/ewv_h/sigx12,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(VF3,
    MARGIN = 1, STATS = wHvar * pmusig/sigx4 * ewu_p * sigx1 *
      ewz/sigx12, FUN = "*") + sweep(vHvar, MARGIN = 1,
    STATS = wHvar * sigx9 * sF3/sigx4 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") - sweep(0.5 * (VF2), MARGIN = 1, STATS = wHvar *
    Wu/sigx4 * ewu_p * sigx1 * ewz/sigx12, FUN = "*") - (sweep(VF2,
    MARGIN = 1, STATS = wHvar * ewu_p * sigx1 * ewz/sigx4 *
      sigx15 * pmusig/sigx12 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx10 * prC/ewv_h * sigx15 * pmusig/sigx12 * ewu_p *
    sigx1 * ewz/sigx12, FUN = "*") + sweep(VF2, MARGIN = 1,
    STATS = wHvar * Q * wzdeno * digamma(P) * gamma(P)/sigx4^2 *
      ewu_p * sigx1 * ewz/sigx12, FUN = "*")))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
      nZHvar + 1)] <- crossprod(sweep(VF2, MARGIN = 1,
    STATS = wHvar * sigx11 * ewu_p * sigx1 * ewz/sigx12,
    FUN = "*") - (sweep(gv, MARGIN = 1, STATS = wHvar * sigx16 *
    ewz/sigx12, FUN = "*") + sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (0.5 * (S^2 * depsiv * (epsilon)^2/(wzdeno * ewv_h^3)) -
      0.5 * (wzdeno * depsiv * ewv_h/(wzdeno * ewv_h)^2)) *
    prC * ewz/sigx12, FUN = "*")), Zvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(((apply((F3)^(P - 1) * log(F3)^2, 1, sum) -
    0.5 * (Wu * sF3)) * ewu_p * sigx1 * ewz * pmusig/sigx12/sigx4 -
    ((sigx15 * ewu_p * sigx1 * ewz * pmusig/sigx12 + 0.5 *
      (Wu)) * sigx15 + Q * ((2 * sF3 - (0.5 * (Wu) + 2 *
      (Q^2 * wzdeno^2 * digamma(P) * gamma(P)^2/sigx4^2)) *
      sDiv) * digamma(P) + (digamma(P)^2 + trigamma(P)) *
      sDiv) * wzdeno * gamma(P)/sigx4^2) * ewu_p * sigx1 *
      ewz * pmusig/sigx12) * wHvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- colSums(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (sigx11 * sF3 - (sigx16 *
      sigx15 * ewz/sigx12 + (0.5 * (Wu * sigx11) + Q *
      ((2 - 2 * (Q^2 * wzdeno^2 * gamma(P)^2/sigx4^2)) *
        ewz + 1) * digamma(P) * gamma(P)/sigx4^2) * sDiv)) *
      ewu_p * sigx1 * ewz * pmusig/sigx12, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar +
    nvZVvar + nZHvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv_h) +
      ewv_h/(wzdeno * ewv_h)^2) * depsiv - (sigx16^2/sigx12 +
      Q * (2 - 2 * (Q^2 * wzdeno * ewz * gamma(P)^2/sigx4^2)) *
        ewu_p * sigx1 * gamma(P) * pmusig * sDiv/sigx4^2)) *
      ewz + sigx11 * ewu_p * sigx1 * pmusig * sDiv - prC *
      depsiv/(wzdeno * ewv_h)) * ewz/sigx12, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf gamma-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
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
zisfgammanormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar,
  nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    itermax = itermax, printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgammanormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfgammanormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradzisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgammanormlike, grad = cgradzisfgammanormlike,
      hess = chesszisfgammanormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2, iterlim = itermax,
        reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(-chesszisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(czisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hess = function(parm) -chesszisfgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisfgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), hessian = function(parm) -chesszisfgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgammanormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
  }
  mlParam <- if (method %in% c("ucminf", "nlminb")) {
    mleObj$par
  } else {
    if (method == "maxLikAlgo") {
      mleObj$estimate
    } else {
      if (method %in% c("sr1", "sparse")) {
        names(mleObj$solution) <- names(startVal)
        mleObj$solution
      } else {
        if (method == "mla") {
          mleObj$b
        }
      }
    }
  }
  if (hessianType != 2) {
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- chesszisfgammanormlike(mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgammanormlike(mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgammanormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgammanormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

# Conditional efficiencies estimation ----------

czisfgammanormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar +
                        object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
    ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
    ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
                                                                       1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1/Hi2
  u_c2 <- rep(0, object$Nobs) 
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) -
      exp(Wv)
    mui_Ki <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu)) +
      exp(Wv)
    Gi <- numeric(object$Nobs)
    Ki <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi[i] <- mean((mui_Gi[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Gi[i]/sqrt(exp(Wv[i])))))^(P -
                                                                            1))
      Ki[i] <- mean((mui_Ki[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui_Ki[i]/sqrt(exp(Wv[i])))))^(P -
                                                                            1))
    }
    teBC_c1 <- exp(exp(Wv)/exp(Wu/2) + object$S * epsilon +
                     exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) - object$S *
                                          epsilon/exp(Wv/2) - exp(Wv/2)) * Gi/(pnorm(-exp(Wv/2 -
                                                                                            Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv)/exp(Wu/2) - object$S *
                                epsilon + exp(Wv)/2) * pnorm(-exp(Wv/2 - Wu/2) -
                                                               object$S * epsilon/exp(Wv/2) + exp(Wv/2)) * Ki/(pnorm(-exp(Wv/2 -
                                                                                                                            Wu/2) - object$S * epsilon/exp(Wv/2)) * Hi2)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
                                teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
                        NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
                        NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
                     odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
                     teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
                     PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
                     u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
                     PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
                     u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
                     ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
                     effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
                     odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
                     u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
                     u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

czisfmarggammanorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar +
                        object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
    ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
                                                                       1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
                        matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggammanorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar +
                        object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
    ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
                                                                       1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
                        matrix(P * exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
