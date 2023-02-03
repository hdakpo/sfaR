################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Generalized Zero Inefficiency Stochastic Frontier Analysis            #
# Number of Classes: 2L                                                        #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for gzisf 2 classes halfnormal-normal distribution
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
#' @noRd
cGZISFhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon2/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for gzisf 2 classes halfnormal-normal distribution
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
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
csGZISFfhalfnorm2C <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csthalfnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initHalf <- NULL
  } else {
    cat("Initialization: SFA + halfnormal - normal distributions...\n")
    initHalf <- maxLik::maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradhalfnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = printInfo, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initHalf$estimate
  }
  StartVal <- c(Esti[1:(nXvar + 1)], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0,
    nvZVvar - 1), 0.98 * Esti[1:nXvar], Esti[nXvar + 2],
    if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar],
    paste0("Zv_", colnames(vHvar)), paste0("Cl1_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for gzisf 2 classes halfnormal-normal distribution
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
#' @noRd
cgradGZISFhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewv1 <- exp(Wv1)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  pmusig1 <- pnorm(-musig1)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  depsi1 <- dnorm(epsi1)
  epsi2 <- S * (epsilon2)/ewv2_h
  depsi2 <- dnorm(epsi2)
  sigx1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon1))
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx3 <- (prC * depsi2/ewv2_h + 2 * sigx2)
  sigx4 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prV1 <- (1 - ewv1/(sigma_sq1))
  sigx5 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx6 <- (1/ssq1 - sigx5 * ewu1/ssq1^2)
  sigx7 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx8 <- (0.5 * sigx7 - sigx6 * dmusig1 * depsi1)
  sigx9 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx10 <- (S * sigx8 * (epsilon1)/wzdeno - 0.5 * sigx9)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  sigx11 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  sigx12 <- (sigx11 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx7)
  sigx13 <- (S * sigx12 * (epsilon1)/wzdeno - 0.5 * sigx9)
  sigx14 <- (sigx3 * ewv2_h^3)
  sigx15 <- (S^2 * depsi2 * (epsilon2)^2/ewv2_h^2)
  sigx16 <- (0.5 * sigx15 - 0.5 * depsi2)
  sigx17 <- (sigx3 * ewv2_h)
  sigx18 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx19 <- (sigx18 * depsi1 * pmusig1)
  sigx20 <- (2 * sigx19 - prC * depsi2/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = 2 * (S *
    sigx1 * ewz/sigx4), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = 2 * (ewu1 * ewz * sigx10/s3sq1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv1 * ewz * sigx13/s3sq1),
      FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = S^2 *
      prC * depsi2 * (epsilon2)/sigx14, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx16 * prC/sigx17, FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = sigx20 * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for gzisf 2 classes halfnormal-normal distribution
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
#' @noRd
chessGZISFhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewv1 <- exp(Wv1)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigma_sq1 <- ewu1 + ewv1
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  musig1 <- (S * ewu1 * (epsilon1)/ssq1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  pmusig1 <- pnorm(-musig1)
  epsi1 <- S * (epsilon1)/sqrt(sigma_sq1)
  depsi1 <- dnorm(epsi1)
  epsi2 <- S * (epsilon2)/ewv2_h
  depsi2 <- dnorm(epsi2)
  sigx1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + S * depsi1 *
    pmusig1 * (epsilon1))
  wzsq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (depsi1 * ewz * pmusig1/wzsq1)
  sigx3 <- (prC * depsi2/ewv2_h + 2 * sigx2)
  sigx4 <- (sigx3 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  prU1 <- (1 - ewu1/(sigma_sq1))
  prV1 <- (1 - ewv1/(sigma_sq1))
  sigx5 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  sigx6 <- (1/ssq1 - sigx5 * ewu1/ssq1^2)
  sigx7 <- (S * depsi1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx8 <- (0.5 * sigx7 - sigx6 * dmusig1 * depsi1)
  sigx9 <- (wzdeno * depsi1 * pmusig1/wzsq1^2)
  sigx10 <- (S * sigx8 * (epsilon1)/wzdeno - 0.5 * sigx9)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  sigx11 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  sigx12 <- (sigx11 * dmusig1 * depsi1 * ewu1/ssq1^2 + 0.5 *
    sigx7)
  sigx13 <- (S * sigx12 * (epsilon1)/wzdeno - 0.5 * sigx9)
  sigx14 <- (sigx3 * ewv2_h^3)
  sigx15 <- (S^2 * depsi2 * (epsilon2)^2/ewv2_h^2)
  sigx16 <- (0.5 * sigx15 - 0.5 * depsi2)
  sigx17 <- (sigx3 * ewv2_h)
  sigx18 <- (1/wzsq1 - ewz * sqrt(sigma_sq1)/wzsq1^2)
  sigx19 <- (sigx18 * depsi1 * pmusig1)
  sigx20 <- (2 * sigx19 - prC * depsi2/(wzdeno * ewv2_h))
  sigx21 <- (S * (dmusig1 * ewu1/sigmastar1 + S * pmusig1 *
    (epsilon1)) * (epsilon1)/(sigma_sq1) - pmusig1)
  sigx22 <- (depsi1 * ewu1/ewv1 + depsi1)
  sigx23 <- (wzdeno * sigx1/(wzsq1^2 * (sigma_sq1)))
  sigx24 <- (s3sq1^2 * wzdeno * (sigma_sq1))
  sigx25 <- (S * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx26 <- (0.5 * sigx25 - sigx6 * dmusig1)
  sigx27 <- (S * depsi1 * (S * sigx26 * (epsilon1) - 2 * (pmusig1/(sigma_sq1))) *
    (epsilon1)/(sigma_sq1)^2)
  sigx28 <- (S^2 * sigx18 * depsi1 * (epsilon1)^2/(sigma_sq1)^2)
  sigx29 <- (0.5 * sigx28 - ((0.5/sqrt(sigma_sq1) - wzdeno^2 *
    sqrt(sigma_sq1)/wzsq1^2) * ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) *
    depsi1/wzsq1^2)
  hessll <- matrix(nrow = 2 * nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = 2 * nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 2 * (S^2 * ((depsi1 * sigx21 + S * dmusig1 *
      sigx22 * ewu1 * (epsilon1)/ssq1)/sigx4 - 2 * (sigx1^2 *
      ewz/sigx4^2)) * ewz), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((sigx6 * dmusig1 *
      depsi1 + S * ((0.5 * sigx21 - 0.5 * pmusig1) * depsi1/(sigma_sq1) -
      S * sigx6 * dmusig1 * sigx22 * (epsilon1)) * (epsilon1)/(sigma_sq1))/wzdeno -
      0.5 * sigx23)/s3sq1 - 2 * (sigx1 * ewz * sigx10/sigx24)) *
      ewu1 * ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((S * ((0.5 * sigx21 - 0.5 * pmusig1) * depsi1/(sigma_sq1) +
    S * sigx11 * dmusig1 * sigx22 * ewu1 * (epsilon1)/ssq1^2) *
    (epsilon1)/(sigma_sq1) - sigx11 * dmusig1 * depsi1 *
    ewu1/ssq1^2)/wzdeno - 0.5 * sigx23)/s3sq1 - 2 * (sigx1 *
    ewz * sigx13/sigx24)) * ewv1 * ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = -wHvar * (2 * (S^3 * prC * sigx1 * depsi2 * ewv2_h^3 *
      ewz * (epsilon2)/(sigx14^2 * wzdeno * (sigma_sq1) *
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S * sigx16 * prC *
      sigx1 * ewv2_h * ewz/(sigx17^2 * wzdeno * (sigma_sq1) *
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * sigx18 - 2 * (sigx20 *
      ewz/(sigx3 * wzdeno * sqrt(sigma_sq1)))) * sigx1 *
      ewz/(sigx3 * (sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu1 * (S * (0.5 * sigx27 - (0.5 * (S^2 * sigx6 *
    depsi1 * (epsilon1)^2/(sigma_sq1)^2) - (((0.5 * (ewu1/(sigma_sq1)) +
    1 - 0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1/sigmastar1 +
    (2 - 2 * (sigx5^2 * ewu1 * (sigma_sq1)/ssq1^2)) * sigmastar1)/ssq1^2 +
    S^2 * sigx6^2 * ewu1 * (epsilon1)^2/ssq1) * depsi1) *
    dmusig1) * (epsilon1)/wzdeno - 0.5 * (wzdeno * (S * sigx8 *
    (epsilon1) - wzdeno^2 * depsi1 * pmusig1/wzsq1^2)/wzsq1^2)) +
    S * sigx8 * (epsilon1)/wzdeno - 0.5 * sigx9)/s3sq1 -
    (0.5 * (sigx3/sqrt(sigma_sq1)) + 2 * (ewz * sigx10)) *
      ewu1 * sigx10/s3sq1^2) * ewu1 * ewz), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * (((((0.5 * (prU1 *
      ewv1) - S^2 * sigx11 * sigx6 * ewu1 * (epsilon1)^2)/(sigma_sq1) +
      0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
        1 - 0.5 * (prU1 * prV1))) * depsi1/sigmastar1 +
      0.5 * (S^2 * sigx11 * depsi1 * (epsilon1)^2/(sigma_sq1)^2)) *
      ewu1 + sigx11 * (1 - 2 * (sigx5 * ewu1 * (sigma_sq1) *
      sigmastar1/ssq1^2)) * depsi1) * dmusig1/ssq1^2 +
      0.5 * sigx27) * (epsilon1)/wzdeno - 0.5 * (wzdeno *
      (S * sigx8 * (epsilon1) - wzdeno^2 * depsi1 * pmusig1/wzsq1^2)/wzsq1^2))/s3sq1 -
      (0.5 * (sigx3/sqrt(sigma_sq1)) + 2 * (ewz * sigx10)) *
        sigx13/s3sq1^2) * ewu1 * ewv1 * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (S^2 * prC * depsi2 *
      ewu1 * ewv2_h^3 * ewz * sigx10 * (epsilon2)/(sigx14^2 *
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar +
    nvZVvar + 1):(2 * nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx16 * prC * ewu1 *
      ewv2_h * ewz * sigx10/(sigx17^2 * sqrt(sigma_sq1)))),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar +
    2 * nvZVvar + 1):(2 * nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (2 * (sigx29 * pmusig1 - S * sigx18 * sigx6 * dmusig1 *
      depsi1 * (epsilon1)) - 2 * (sigx20 * ewz * sigx10/s3sq1)) *
    ewu1 * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) -
      0.5 * (0.5 * prV1 + ewv1/(sigma_sq1))) * prV1 + S^2 *
      sigx11^2 * ewu1 * ewv1 * (epsilon1)^2/(ssq1^2 * (sigma_sq1))) *
      depsi1 * ewu1/sigmastar1 + ((0.5 * (S^2 * depsi1 *
      (epsilon1)^2/(sigma_sq1)^2) - 2 * (sigx11 * depsi1 *
      (sigma_sq1) * sigmastar1/ssq1^2)) * ewv1 + depsi1) *
      sigx11) * dmusig1 * ewu1/ssq1^2 + S * (0.5 * (ewv1 *
      (S * (sigx11 * dmusig1 * ewu1/ssq1^2 + 0.5 * sigx25) *
        (epsilon1) - 2 * (pmusig1/(sigma_sq1)))) + 0.5 *
      pmusig1) * depsi1 * (epsilon1)/(sigma_sq1)^2) * (epsilon1)/wzdeno -
      (0.5 * (depsi1 * pmusig1) + 0.5 * (ewv1 * (S * sigx12 *
        (epsilon1) - wzdeno^2 * depsi1 * pmusig1/wzsq1^2))) *
        wzdeno/wzsq1^2)/s3sq1 - (0.5 * (sigx3/sqrt(sigma_sq1)) +
      2 * (ewz * sigx13)) * ewv1 * sigx13/s3sq1^2) * ewv1 *
      ewz), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (S^2 * prC * depsi2 * ewv1 * ewv2_h^3 * ewz * sigx13 *
      (epsilon2)/(sigx14^2 * sqrt(sigma_sq1)))), FUN = "*"),
    Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (2 * (sigx16 * prC * ewv1 * ewv2_h *
      ewz * sigx13/(sigx17^2 * sqrt(sigma_sq1)))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 * nXvar +
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx29 * pmusig1 +
      S * sigx11 * sigx18 * dmusig1 * depsi1 * ewu1 * (epsilon1)/ssq1^2) -
      2 * (sigx20 * ewz * sigx13/s3sq1)) * ewv1 * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * ((S^2 * (epsilon2)^2/ewv2_h^2 -
      1)/sigx14 - S^2 * prC * depsi2 * (epsilon2)^2/sigx14^2) *
      prC * depsi2, FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * ((0.5 * (S^2 * (epsilon2)^2/ewv2_h^2 -
      2) - 0.5)/sigx17 - sigx16 * prC/sigx17^2) * prC *
      depsi2 * (epsilon2)/ewv2_h^2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar), (2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S^2 * (sigx20/sigx3 + 1/wzdeno) *
      prC * depsi2 * ewz * (epsilon2)/sigx14), FUN = "*"),
    Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    2 * nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon2)^2/ewv2_h^2) - 1) - 0.25) * depsi2 *
      (epsilon2)^2/sigx14 - (sigx16 * prC + 0.5 * sigx17) *
      sigx16/sigx17^2), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar +
    2 * nvZVvar), (2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 *
    nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx16 * sigx20/sigx3 +
      0.5 * (S^2 * depsi2 * (epsilon2)^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * depsi2 * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx3), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + 2 * nvZVvar + 1):(2 * nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar), (2 * nXvar + nuZUvar +
    2 * nvZVvar + 1):(2 * nXvar + nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((prC * (1/(wzdeno^2 * ewv2_h) + ewv2_h/(wzdeno * ewv2_h)^2) *
      depsi2 - (sigx20^2/sigx3 + 2 * ((2 - 2 * (wzdeno *
      (sigma_sq1) * ewz/wzsq1^2)) * depsi1 * pmusig1 *
      sqrt(sigma_sq1)/wzsq1^2))) * ewz + 2 * sigx19 - prC *
      depsi2/(wzdeno * ewv2_h)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for gzisf 2 classes halfnormal-normal distribution
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
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
GZISF2ChnormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csGZISFfhalfnorm2C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(cGZISFhalfnormlike2C(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("GZISF 2 Classes Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cGZISFhalfnormlike2C(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradGZISFhalfnormlike2C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cGZISFhalfnormlike2C,
    grad = cgradGZISFhalfnormlike2C, hess = chessGZISFhalfnormlike2C,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cGZISFhalfnormlike2C(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessGZISFhalfnormlike2C(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradGZISFhalfnormlike2C(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessGZISFhalfnormlike2C(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessGZISFhalfnormlike2C(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cGZISFhalfnormlike2C(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradGZISFhalfnormlike2C(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) initHalf = initHalf))
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for gzisf 2 classes halfnormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
cGZISF2Chalfnormeff <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon2/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) *
      pnorm(mustar1/sigmastar1 + sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal effects for for gzisf 2 classes halfnormal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargGZISF2Chalfnorm_Eu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon2/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c <- ifelse(Group_c == 1, margEff_c1, margEff_c2)
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2)
  return(margEff)
}

cmargGZISF2Chalfnorm_Vu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon2/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar],
    nrow = 1), matrix(exp(Wu1) * (1 - (dnorm(0)/pnorm(0))^2),
    ncol = 1))
  colnames(margEff_c1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  margEff_c2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  colnames(margEff_c2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  margEff_c <- ifelse(Group_c == 1, margEff_c1, margEff_c2)
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  margEff <- data.frame(margEff_c, margEff_c1, margEff_c2)
  return(margEff)
}
