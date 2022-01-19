################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Convolution: truncated normal - normal                                       #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
czisftruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, 
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon)/(exp(Wu) + 
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + 
    exp(Wv))) * pnorm(mustar/sigmastar)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar*log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf truncatednormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
cstzisftruncnorm <- function(olsObj, epsiRes, nXvar, nmuZUvar, 
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, 
  nZHvar, itermax, printInfo, tol) {
  cat("Initialization: SFA + truncated normal - normal distributions...\n")
  initTrunc <- maxLik(logLik = ctruncnormlike, start = csttruncnorm(olsObj = olsObj, 
    epsiRes = epsiRes, S = S, uHvar = uHvar[, 1, 
      drop = FALSE], nuZUvar = 1, vHvar = vHvar[, 1, drop = FALSE], 
    nvZVvar = 1, nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE]), 
    grad = cgradtruncnormlike, method = "BFGS", control = list(iterlim = itermax, 
      printLevel = printInfo, reltol = tol), nXvar = nXvar, 
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]), 
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar, 
    S = S, wHvar = wHvar, nmuZUvar = 1, muHvar = as.matrix(muHvar[, 1]))
  Esti <- initTrunc$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], rep(0, nmuZUvar - 1), 
    Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar - 1), 
    Esti[nXvar + 3], if (nvZVvar > 1) rep(0, nvZVvar - 1), 
    rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_", 
    colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_", 
    colnames(vHvar)), paste0("SF_", colnames(Zvar)))
names(initTrunc$estimate) <- c(
    names(Esti)[1:nXvar], paste0(
      "Zmu_",
      colnames(muHvar)[1]
    ), paste0("Zu_", colnames(uHvar)[1]),
    paste0("Zv_", colnames(vHvar)[1])
  )
  return(list(StartVal = StartVal, initTrunc = initTrunc))
}

# Gradient of the likelihood function ----------
#' gradient for zisf truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
cgradzisftruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, 
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  dmusig <- dnorm((mu + S * (epsilon))/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/exp(Wu/2))
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  dwu <- dnorm(mu/exp(Wu/2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  muepsi <- (mu + S * (epsilon))
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  wzsq <- (wzdeno * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  muwusq <- (mu * dwu * sqrt(sigma_sq)/exp(Wu/2))
  wzdpsq <- (wzdeno * dmusig * pmusig * pwu/wzsq^2)
  sigx1 <- (mu/ssq - s2sig * (mu * ewv - S * ewu * (epsilon))/ssq^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi * 
    pmusig)
  sigx3 <- (wzdeno * (sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (sigx2 * ewz/sigx3 + S * prC * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (prC * depsi/ewvsr + dmusig * ewz * pmusig/wzsq)
  sigx6 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi * 
    pmusig)
  sigx7 <- (wzsq^2 * exp(Wu/2))
  sigx8 <- (sigx6/sigx3 - wzdeno * dmusig * dwu * pmusig * 
    sqrt(sigma_sq)/sigx7)
  sigx9 <- ((1 - ewu/(sigma_sq)) * ewv/sigmastar)
  sigx10 <- ((0.5 * sigx9 + sigmastar) * (mu * ewv - S * ewu * 
    (epsilon))/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * dpmu - sigx10 * dmu * dmusig)
  sigx12 <- (sigx11 * ewu/wzsq - (0.5 * ewusq - 0.5 * muwusq) * 
    wzdeno * dmusig * pmusig/wzsq^2)
  sigx13 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/(wzdeno * 
    pwu) - 0.5 * wzdpsq)
  sigx14 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 * 
    depsi)
  sigx15 <- (sigx13 * ewv * ewz/sqrt(sigma_sq) + sigx14 * prC/ewvsr)
  sigx16 <- (1/wzsq - ewz * pwu * sqrt(sigma_sq)/wzsq^2)
  sigx17 <- (sigx16 * dmusig * pmusig - prC * depsi/(wzdeno * 
    ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx4/sigx5, 
    FUN = "*"), sweep(muHvar, MARGIN = 1, STATS = sigx8 * 
    ewz/sigx5, FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx12 * 
    ewz/sigx5, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx15/sigx5, 
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17 * 
    ewz/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
chesszisftruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, 
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  sigma_sq <- ewu + ewv
  sigmastar <- sqrt(ewu * ewv/(sigma_sq))
  dmusig <- dnorm((mu + S * (epsilon))/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  depsi <- dnorm(S * (epsilon)/ewvsr)
  pwu <- pnorm(mu/exp(Wu/2))
  dmu <- dnorm((mu * ewv - S * ewu * (epsilon))/ssq)
  dwu <- dnorm(mu/exp(Wu/2))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  prV <- (1 - ewv/(sigma_sq))
  s2sig <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  muepsi <- (mu + S * (epsilon))
  dpmu <- (dmusig * muepsi^2 * pmusig/(sigma_sq)^2)
  wzsq <- (wzdeno * pwu * sqrt(sigma_sq))
  ewusq <- (ewu * pwu/sqrt(sigma_sq))
  muwusq <- (mu * dwu * sqrt(sigma_sq)/exp(Wu/2))
  wzdpsq <- (wzdeno * dmusig * pmusig * pwu/wzsq^2)
  sigx1 <- (mu/ssq - s2sig * (mu * ewv - S * ewu * (epsilon))/ssq^2)
  sigx2 <- (dmu * dmusig * ewu/sigmastar + dmusig * muepsi * 
    pmusig)
  sigx3 <- (wzdeno * (sigma_sq) * pwu * sqrt(sigma_sq))
  sigx4 <- (sigx2 * ewz/sigx3 + S * prC * depsi * (epsilon)/ewvsr^3)
  sigx5 <- (prC * depsi/ewvsr + dmusig * ewz * pmusig/wzsq)
  sigx6 <- (dmu * dmusig * ewv/sigmastar - dmusig * muepsi * 
    pmusig)
  sigx7 <- (wzsq^2 * exp(Wu/2))
  sigx8 <- (sigx6/sigx3 - wzdeno * dmusig * dwu * pmusig * 
    sqrt(sigma_sq)/sigx7)
  sigx9 <- ((1 - ewu/(sigma_sq)) * ewv/sigmastar)
  sigx10 <- ((0.5 * sigx9 + sigmastar) * (mu * ewv - S * ewu * 
    (epsilon))/ssq^2 + S * (epsilon)/ssq)
  sigx11 <- (0.5 * dpmu - sigx10 * dmu * dmusig)
  sigx12 <- (sigx11 * ewu/wzsq - (0.5 * ewusq - 0.5 * muwusq) * 
    wzdeno * dmusig * pmusig/wzsq^2)
  sigx13 <- ((0.5 * dpmu + dmu * dmusig * sigx1)/(wzdeno * 
    pwu) - 0.5 * wzdpsq)
  sigx14 <- (0.5 * (S^2 * depsi * (epsilon)^2/ewvsr^2) - 0.5 * 
    depsi)
  sigx15 <- (sigx13 * ewv * ewz/sqrt(sigma_sq) + sigx14 * prC/ewvsr)
  sigx16 <- (1/wzsq - ewz * pwu * sqrt(sigma_sq)/wzsq^2)
  sigx17 <- (sigx16 * dmusig * pmusig - prC * depsi/(wzdeno * 
    ewvsr))
  sigx18 <- (0.5 * (muepsi^2/(sigma_sq)) - 2) * pmusig/(sigma_sq)
  sigx19 <- (0.5 * dpmu + dmu * dmusig * sigx1) * dwu/((wzdeno * 
    pwu)^2 * exp(Wu/2))
  prU <- (1 - ewu/(sigma_sq))
  muepsiuv <- (mu * ewv - S * ewu * (epsilon))
  sigx20 <- (((2 - muepsi^2/(sigma_sq)) * pmusig + dmu * ewv * 
    muepsi/ssq) * dmusig * muepsi/(sigma_sq)^2)
  ewusr <- exp(Wu/2)
  muqwuq <- (0.5 * ewusq - 0.5 * muwusq)
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar + 
    nZHvar, ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar + 
    nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = wHvar*S^2 * ((((muepsi^2/(sigma_sq) - 1) * pmusig + 
      dmu * ewu * muepsi/ssq) * dmusig + dmu * (dmusig * 
      muepsi - dmusig * muepsiuv/ewv) * ewu/ssq) * ewz/sigx3 + 
      prC * depsi * (S^2 * (epsilon)^2/ewvsr^2 - 1)/ewvsr^3 - 
      sigx4^2/sigx5)/sigx5, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * (((dmu * (dmusig * muepsi - dmusig * 
      muepsiuv/ewv) * ewv/ssq - ((muepsi^2/(sigma_sq) - 
      1) * pmusig + dmu * ewu * muepsi/ssq) * dmusig)/wzsq - 
      wzdeno * sigx2 * dwu * sqrt(sigma_sq)/sigx7)/(sigma_sq) - 
      sigx4 * sigx8/sigx5) * ewz/sigx5, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar*S * 
    ((0.5 * (((muepsi^2/(sigma_sq) - 2) * pmusig + dmu * 
      ewu * muepsi/ssq) * dmusig * muepsi/(sigma_sq)^2) - 
      (sigx10 * dmusig * muepsi/(sigma_sq) + ((0.5 * sigx9 + 
        sigmastar) * ewu/ssq^2 - (sigx10 * muepsiuv/ewv + 
        1/sigmastar)/(sigma_sq)) * dmusig) * dmu) * ewu/wzsq - 
      (sigx12 * sigx4/sigx5 + muqwuq * wzdeno * sigx2/(wzsq^2 * 
        (sigma_sq)))) * ewz/sigx5, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * (((((dmusig * muepsi - dmusig * 
      muepsiuv/ewv) * sigx1/(sigma_sq) - s2sig * dmusig * 
      ewu/ssq^2) * dmu + 0.5 * (((muepsi^2/(sigma_sq) - 
      2) * pmusig + dmu * ewu * muepsi/ssq) * dmusig * 
      muepsi/(sigma_sq)^2))/(wzdeno * pwu) - 0.5 * (wzdeno * 
      sigx2 * pwu/(wzsq^2 * (sigma_sq)))) * ewv * ewz/sqrt(sigma_sq) + 
      S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * 
        prC * depsi * (epsilon)/ewvsr^3 - sigx15 * sigx4/sigx5)/sigx5, 
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * (sigx16 * sigx2/(sigma_sq) - 
      (sigx17 * sigx4/sigx5 + S * prC * depsi * (epsilon)/(wzdeno * 
        ewvsr^3))) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + 
    nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar*(((((1 - 
    muepsi^2/(sigma_sq)) * pmusig + dmu * ewv * muepsi/ssq) * 
    dmusig + dmu * (dmusig * muepsiuv/ewu + dmusig * muepsi) * 
    ewv/ssq)/sigx3 + ((sigx6/(sigx3^2 * ewusr) - 2 * (wzdeno^2 * 
    dmusig * dwu * pmusig * pwu/sigx7^2)) * (sigma_sq) + 
    (dmu * dmusig * ewv/ssq - (dmusig * muepsi/(sigma_sq) + 
      mu * dmusig/ewusr^2) * pmusig)/sigx7) * wzdeno * 
    dwu * sqrt(sigma_sq) + sigx8^2 * ewz/sigx5) * ewz/sigx5), 
    FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, 
    MARGIN = 1, STATS = wHvar*(((0.5 * sigx20 - (((0.5 * sigx9 + 
      sigmastar) * ewv/ssq^2 - sigx10 * muepsiuv/(ewu * 
      (sigma_sq))) * dmusig - sigx10 * dmusig * muepsi/(sigma_sq)) * 
      dmu)/wzsq - sigx11 * wzdeno * dwu * sqrt(sigma_sq)/sigx7) * 
      ewu - ((((0.5 * (ewu/sqrt(sigma_sq)) - 0.5 * ((1 - 
      mu^2/ewusr^2) * sqrt(sigma_sq))) * dmusig * dwu/ewusr - 
      muqwuq * dmusig * muepsi/(sigma_sq)) * pmusig + muqwuq * 
      (dmu * ewv/ssq - 2 * (wzdeno^2 * dwu * (sigma_sq) * 
        pmusig * pwu/sigx7)) * dmusig) * wzdeno/wzsq^2 + 
      sigx12 * sigx8 * ewz/sigx5)) * ewz/sigx5, FUN = "*"), 
    uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, 
    MARGIN = 1, STATS = wHvar*(((((1/ssq - s2sig * ewv/ssq^2) * 
      dmusig - (dmusig * muepsiuv/ewu + dmusig * muepsi) * 
      sigx1/(sigma_sq)) * dmu + 0.5 * sigx20)/(wzdeno * 
      pwu) - (sigx19 + 0.5 * (((1 - 2 * (wzdeno^2 * (sigma_sq) * 
      pwu^2/wzsq^2)) * dmusig * dwu * pmusig/ewusr + sigx6 * 
      pwu/(sigma_sq))/wzsq^2)) * wzdeno) * ewv/sqrt(sigma_sq) - 
      sigx15 * sigx8/sigx5) * ewz/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 
    nvZVvar + nZHvar)] <- crossprod(sweep(muHvar, MARGIN = 1, 
    STATS = wHvar*(sigx16 * dmu * dmusig * ewv/ssq - ((((2 - 2 * 
      (wzdeno^2 * (sigma_sq) * pwu^2/wzsq^2)) * ewz + 1) * 
      dmusig * dwu * sqrt(sigma_sq)/sigx7 + sigx16 * dmusig * 
      muepsi/(sigma_sq)) * pmusig + sigx17 * sigx8 * ewz/sigx5)) * 
      ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), 
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = wHvar*((((0.5 * ((sigx18 - sigx10 * dmu) * 
      ewu) + 0.5 * pmusig) * dmusig * muepsi^2/(sigma_sq)^2 - 
      ((sigx10^2 * ewu * muepsiuv/ssq + ((0.5 * (ewu/(sigma_sq)) - 
        0.5 * (0.5 * prU + ewu/(sigma_sq))) * prU * ewv * 
        muepsiuv/sigmastar - (0.5 * sigx9 + sigmastar) * 
        (2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq) * 
          muepsiuv * sigmastar/ssq^2) + 2 * (S * (epsilon))) * 
        ewu)/ssq^2) * dmusig + sigx10 * (0.5 * (dmusig * 
        ewu * muepsi^2/(sigma_sq)^2) + dmusig)) * dmu)/wzsq - 
      sigx11 * muqwuq * wzdeno/wzsq^2) * ewu - ((((0.5 * 
      (((1 - 0.5 * (ewu/(sigma_sq))) * pwu - 0.5 * (mu * 
        dwu/ewusr)) * ewu/sqrt(sigma_sq)) - 0.5 * (mu * 
      ((0.5 * (mu^2/ewusr^2) - 0.5) * sqrt(sigma_sq) + 
        0.5 * (ewu/sqrt(sigma_sq))) * dwu/ewusr)) * dmusig + 
      0.5 * (muqwuq * dmusig * ewu * muepsi^2/(sigma_sq)^2)) * 
      pmusig - (sigx10 * dmu * ewu + 2 * (muqwuq * wzdeno^2 * 
      pmusig * pwu * sqrt(sigma_sq)/wzsq^2)) * muqwuq * 
      dmusig) * wzdeno/wzsq^2 + sigx12^2 * ewz/sigx5)) * 
      ewz/sigx5, FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), 
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
      nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, 
    STATS = wHvar*((((((sigx10 * dmusig * muepsiuv/sigmastar + 
      0.5 * (dmusig * muepsi^2/(sigma_sq))) * sigx1/(sigma_sq) - 
      ((0.5 * (prU * ewv/(sigma_sq)) + 0.5 * ((ewu/(sigma_sq) - 
        1) * ewv/(sigma_sq) + 1 - 0.5 * (prU * prV))) * 
        muepsiuv/sigmastar + mu * (0.5 * sigx9 + sigmastar) - 
        s2sig * (2 * ((0.5 * sigx9 + sigmastar) * (sigma_sq) * 
          muepsiuv * sigmastar/ssq^2) + S * (epsilon))) * 
        dmusig/ssq^2) * dmu + 0.5 * ((sigx18 - sigx10 * 
      dmu) * dmusig * muepsi^2/(sigma_sq)^2))/(wzdeno * 
      pwu) - 0.5 * (sigx13/(sigma_sq))) * ewu + (0.5 * 
      (mu * sigx19) - 0.5 * ((sigx11 * ewu * pwu - (0.5 * 
      (mu * dwu/ewusr) + 2 * (muqwuq * wzdeno^2 * pwu^2 * 
      sqrt(sigma_sq)/wzsq^2)) * dmusig * pmusig)/wzsq^2)) * 
      wzdeno) * ewv/sqrt(sigma_sq) - sigx15 * sigx12/sigx5) * 
      ewz/sigx5, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), 
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + 
      nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = wHvar*((0.5 * (sigx16 * dmusig * ewu * 
      muepsi^2/(sigma_sq)^2) - ((2 - 2 * (wzdeno^2 * (sigma_sq) * 
      pwu^2/wzsq^2)) * ewz + 1) * muqwuq * dmusig/wzsq^2) * 
      pmusig - (sigx10 * sigx16 * dmu * dmusig * ewu + 
      sigx12 * sigx17 * ewz/sigx5)) * ewz/sigx5, FUN = "*"), 
    Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = wHvar*(((((0.5 * (dmusig * muepsi^2/(sigma_sq)) - 
      dmusig * muepsiuv * sigx1/sigmastar) * ewv * sigx1/(sigma_sq) + 
      dmusig * (mu/ssq - (((3 * (mu) - 2 * (s2sig * (sigma_sq) * 
        muepsiuv * sigmastar/ssq^2)) * ewv - S * ewu * 
        (epsilon)) * s2sig + (0.5 * (ewv/(sigma_sq)) - 
        0.5 * (0.5 * prV + ewv/(sigma_sq))) * prV * ewu * 
        muepsiuv/sigmastar)/ssq^2)) * dmu + (0.5 * ((sigx18 + 
      dmu * sigx1) * ewv) + 0.5 * pmusig) * dmusig * muepsi^2/(sigma_sq)^2)/(wzdeno * 
      pwu) - ((0.5 * (((dmu * sigx1 - wzdeno^2 * pmusig * 
      pwu^2/wzsq^2) * dmusig + 0.5 * dpmu) * ewv) + 0.5 * 
      (dmusig * pmusig)) * wzdeno * pwu/wzsq^2 + 0.5 * 
      (sigx13 * ewv/(sigma_sq)))) * ewv * ewz/sqrt(sigma_sq) + 
      prC * (S^2 * (0.5 * (0.5 * (S^2 * (epsilon)^2/ewvsr^2) - 
        1) - 0.25) * depsi * (epsilon)^2/ewvsr^2 - 0.5 * 
        sigx14)/ewvsr - sigx15^2/sigx5)/sigx5, FUN = "*"), 
    vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + 
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar + 
    1):(nXvar + nmuZUvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = wHvar*(((0.5 * (sigx16 * dmusig * muepsi^2/(sigma_sq)^2) - 
      ((0.5/sqrt(sigma_sq) - wzdeno^2 * pwu^2 * sqrt(sigma_sq)/wzsq^2) * 
        ewz + 0.5 * (wzdeno/sqrt(sigma_sq))) * dmusig * 
        pwu/wzsq^2) * pmusig + sigx16 * dmu * dmusig * 
      sigx1) * ewv - (sigx15 * sigx17/sigx5 + (0.5 * (S^2 * 
      depsi * (epsilon)^2/(wzdeno * ewvsr^3)) - 0.5 * (wzdeno * 
      depsi * ewvsr/(wzdeno * ewvsr)^2)) * prC)) * ewz/sigx5, 
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + 
    nmuZUvar + nuZUvar + nvZVvar + nZHvar), (nXvar + nmuZUvar + 
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 
    nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, 
    STATS = wHvar*((prC * (1/(wzdeno^2 * ewvsr) + ewvsr/(wzdeno * 
      ewvsr)^2) * depsi - (sigx17^2/sigx5 + (2 - 2 * (wzdeno * 
      (sigma_sq) * ewz * pwu^2/wzsq^2)) * dmusig * pmusig * 
      pwu * sqrt(sigma_sq)/wzsq^2)) * ewz + sigx16 * dmusig * 
      pmusig - prC * depsi/(wzdeno * ewvsr)) * ewz/sigx5, 
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf truncatednormal-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
zisftruncnormAlgOpt <- function(start, olsParam, dataTable, S, wHvar, 
  nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, 
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax, 
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start)) 
    start else cstzisftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], 
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, 
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, 
    itermax = itermax, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  startLoglik <- sum(czisftruncnormlike(startVal, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, 
    nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
                         bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
                         nm = function(...) maxNM(...), cg = function(...) maxCG(...),
                         sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal, 
    fn = function(parm) -sum(czisftruncnormlike(parm, nXvar = nXvar, 
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, 
      nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)), 
    gr = function(parm) -colSums(cgradzisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo, 
      maxeval = itermax, stepmax = stepmax, xtol = tol, 
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = czisftruncnormlike, 
    grad = cgradzisftruncnormlike, hess = chesszisftruncnormlike, 
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE, 
    control = list(printLevel = if (printInfo) 2, iterlim = itermax, 
      reltol = tol, tol = tol, qac = qac), nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, 
    nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trust.optim(x = startVal, 
    fn = function(parm) -sum(czisftruncnormlike(parm, nXvar = nXvar, 
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, 
      nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)), 
    gr = function(parm) -colSums(cgradzisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax, 
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, 
      report.level = if (printInfo) 4L else 0, report.precision = 1L)), 
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(czisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), hs = function(parm) as(-chesszisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", 
      control = list(maxit = itermax, cgtol = gradtol, 
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0, 
        report.precision = 1L, preconditioner = 1L)), 
    mla = mla(b = startVal, fn = function(parm) -sum(czisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar)), hess = function(parm) -chesszisftruncnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
      nZHvar = nZHvar), print.info = printInfo, maxiter = itermax, 
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, 
      objective = function(parm) -sum(czisftruncnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, 
        Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisftruncnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, 
        Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesszisftruncnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, 
        Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, 
        trace = printInfo, eval.max = itermax, rel.tol = tol, 
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftruncnormlike(mleObj$par, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
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
    if (method %in% c("ucminf", "nlminb")) 
      mleObj$hessian <- chesszisftruncnormlike(parm = mleObj$par, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, 
        Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1") 
      mleObj$hessian <- chesszisftruncnormlike(parm = mleObj$solution, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, 
        Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftruncnormlike(parm = mlParam, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, 
    nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftruncnormlike(parm = mlParam, 
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, 
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, 
    mleObj = mleObj, mlParam = mlParam, initTrunc = initTrunc))
}

# Conditional efficiencies estimation ----------

czisftruncnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
                                                                  object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar + 1):(object$nXvar +
                                                                                object$nmuZUvar+object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
                         rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon)/(exp(Wu) + 
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) + 
    exp(Wv))) * pnorm(mustar/sigmastar)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
                                                          sigmastar)/pnorm(mustar/sigmastar)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
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
}

# Marginal effects on inefficiencies ----------

czisfmargtruncnorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
                                                                  object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar + 1):(object$nXvar +
                                                                                object$nmuZUvar+object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
                         rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  Lambda <- mu/exp(Wu/2)
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon)/(exp(Wu) + 
                                                             exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) + 
                                                                            exp(Wv))) * pnorm(mustar/sigmastar)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
                      matrix(1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2,
                             ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
                      matrix(exp(Wu/2)/2 * ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
                                              Lambda * (dnorm(Lambda)/pnorm(Lambda))^2), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
                                                             4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
                                                             5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
                   mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                        colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                        colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                         colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargtruncnorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
                                                                  object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar + 1):(object$nXvar +
                                                                                object$nmuZUvar+object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                             object$nvZVvar + 1):(object$nXvar +  object$nmuZUvar+object$nuZUvar +
                                                    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
                         rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
                        rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
                       rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) -
                     (dnorm(Lambda)/pnorm(Lambda))^2)
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon)/(exp(Wu) + 
                                                             exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) + 
                                                                            exp(Wv))) * pnorm(mustar/sigmastar)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
                      matrix(1/exp(Wu/2) * dnorm(Lambda)/pnorm(Lambda) * (m1^2 -
                                                                            m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
                      matrix(exp(Wu) * (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) *
                                          (Lambda + Lambda^3 + (2 + 3 * Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
                                             2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
                             ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
                                                             4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
                                                             5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
                   mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                        colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                        colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
                                         colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(bind_cols(margEff1, margEff2, margEff_c))
}
