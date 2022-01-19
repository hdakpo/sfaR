################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf halfnormal-normal distribution
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
czisfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf halfnormal-normal distribution
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
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
cstzisfhalfnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA halfnormal - normal distribution...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradhalfnormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = printInfo,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("SF_",
    colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------
#' gradient for zisf halfnormal-normal distribution
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
cgradzisfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
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
  dmusig <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  depsi <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (dmusig * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx2 <- (prC * dwsr/ewvsr + 2 * sigx1)
  sigx3 <- (S * dmusig * pmusig * (epsilon)/(sigma_sq)^2)
  sigx4 <- (wzdeno * dmusig * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (2 *
    ((depsi * dmusig * ewu/sigmastar + S * dmusig * pmusig *
      (epsilon)) * ewz/(wzdeno * (sigma_sq)^(3/2))) + S *
    prC * dwsr * (epsilon)/ewvsr^3)/sigx2, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (ewu * ewz * (S * (0.5 * sigx3 -
      (1/(ssq) - (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
        sigmastar) * ewu/(ssq)^2) * depsi * dmusig) *
      (epsilon)/wzdeno - 0.5 * sigx4)/(sigx2 * sqrt(sigma_sq))),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ((0.5 *
    (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 * dwsr) * prC/ewvsr +
    2 * (ewv * ewz * (S * ((0.5 * ((1 - ewv/(sigma_sq)) *
      ewu/sigmastar) + sigmastar) * depsi * dmusig * ewu/(ssq)^2 +
      0.5 * sigx3) * (epsilon)/wzdeno - 0.5 * sigx4)/sqrt(sigma_sq)))/sigx2,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 * ((1/(wzdeno *
    sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno * sqrt(sigma_sq))^2) *
    dmusig * pmusig) - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx2,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf halfnormal-normal distribution
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
chesszisfhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
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
  dmusig <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  ssq <- (sigma_sq) * sigmastar
  pmusig <- pnorm(-(S * ewu * (epsilon)/(ssq)))
  depsi <- dnorm(-(S * ewu * (epsilon)/(ssq)))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigx1 <- (dmusig * ewz * pmusig/(wzdeno * sqrt(sigma_sq)))
  sigx2 <- (prC * dwsr/ewvsr + 2 * sigx1)
  sigx3 <- (S * dmusig * pmusig * (epsilon)/(sigma_sq)^2)
  sigx4 <- (wzdeno * dmusig * pmusig/(wzdeno * sqrt(sigma_sq))^2)
  sigx5 <- (depsi * dmusig * ewu/sigmastar + S * dmusig * pmusig *
    (epsilon))
  sigx6 <- (sigx5 * ewz/(wzdeno * (sigma_sq)^(3/2)))
  sigx7 <- (2 * sigx6 + S * prC * dwsr * (epsilon)/ewvsr^3)
  sigx8 <- (0.5 * ((1 - ewv/(sigma_sq)) * ewu/sigmastar) +
    sigmastar)
  sigx9 <- (wzdeno * sigx5/((wzdeno * sqrt(sigma_sq))^2 * (sigma_sq)))
  sigx10 <- sigx8 * depsi * dmusig * ewu/(ssq)^2
  sigx11 <- (depsi * ewu/sigmastar + S * pmusig * (epsilon))
  sigx12 <- (S * sigx11 * (epsilon)/(sigma_sq) - pmusig)
  sigx13 <- (0.5 * ((1 - ewu/(sigma_sq)) * ewv/sigmastar) +
    sigmastar)
  sigx14 <- (1/(ssq) - sigx13 * ewu/(ssq)^2)
  sigx15 <- (dmusig * ewu/ewv + dmusig)
  ddmu <- (0.5 * sigx3 - sigx14 * depsi * dmusig)
  sigx16 <- (S * ddmu * (epsilon)/wzdeno - 0.5 * sigx4)
  sigx17 <- (S * pmusig * (epsilon)/(sigma_sq)^2)
  sigx18 <- (S * dmusig * (S * (0.5 * sigx17 - sigx14 * depsi) *
    (epsilon) - 2 * (pmusig/(sigma_sq))) * (epsilon)/(sigma_sq)^2)
  sigx19 <- ((1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2) * dmusig * pmusig)
  sigx20 <- (2 * sigx19 - prC * dwsr/(wzdeno * ewvsr))
  sigx21 <- (ewv * ewz * (S * (sigx10 + 0.5 * sigx3) * (epsilon)/wzdeno -
    0.5 * sigx4)/sqrt(sigma_sq))
  sigx22 <- ((0.5 * (S^2 * dwsr * (epsilon)^2/ewvsr^2) - 0.5 *
    dwsr) * prC/ewvsr + 2 * sigx21)
  sigx23 <- (0.5/sqrt(sigma_sq) - wzdeno^2 * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2)
  sigx24 <- (sigx23 * ewz + 0.5 * (wzdeno/sqrt(sigma_sq)))
  sigx25 <- (S^2 * (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/(wzdeno *
    sqrt(sigma_sq))^2) * dmusig * (epsilon)^2/(sigma_sq)^2)
  sigx26 <- (0.5 * sigx25 - sigx24 * dmusig/(wzdeno * sqrt(sigma_sq))^2)
  wzsq <- (wzdeno * sqrt(sigma_sq))^2
  sigx27 <- (wzdeno * (S * ddmu * (epsilon) - wzdeno^2 * dmusig *
    pmusig/wzsq)/wzsq)
  sigx28 <- (1/(wzdeno * sqrt(sigma_sq)) - ewz * sqrt(sigma_sq)/wzsq)
  ewvsq <- (1 - ewv/(sigma_sq))
  sigx29 <- (S * (sigx10 + 0.5 * sigx3) * (epsilon)/wzdeno -
    0.5 * sigx4)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 -
      1)/ewvsr^3 + 2 * ((dmusig * sigx12 + S * depsi *
      sigx15 * ewu * (epsilon)/(ssq)) * ewz/(wzdeno * (sigma_sq)^(3/2))) -
      sigx7^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((sigx14 * depsi *
      dmusig + S * ((0.5 * sigx12 - 0.5 * pmusig) * dmusig/(sigma_sq) -
      S * sigx14 * depsi * sigx15 * (epsilon)) * (epsilon)/(sigma_sq))/wzdeno -
      0.5 * sigx9)/(sigx2 * sqrt(sigma_sq)) - sigx7 * sigx16 *
      sqrt(sigma_sq)/(sigx2 * sqrt(sigma_sq))^2) * ewu *
      ewz), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (2 * (((S * ((0.5 * sigx12 - 0.5 * pmusig) * dmusig/(sigma_sq) +
    S * sigx8 * depsi * sigx15 * ewu * (epsilon)/(ssq)^2) *
    (epsilon)/(sigma_sq) - sigx10)/wzdeno - 0.5 * sigx9) *
    ewv * ewz/sqrt(sigma_sq)) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 -
    2) - 0.5) * prC * dwsr * (epsilon)/ewvsr^3 - sigx22 *
    sigx7/sigx2)/sigx2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx28 * sigx5/(sigma_sq)) -
      (sigx20 * sigx7/sigx2 + S * prC * dwsr * (epsilon)/(wzdeno *
        ewvsr^3))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((ewu * (S * (0.5 * sigx18 - (0.5 * (S^2 * sigx14 *
    dmusig * (epsilon)^2/(sigma_sq)^2) - (((0.5 * (ewu/(sigma_sq)) +
    1 - 0.5 * (0.5 * (1 - ewu/(sigma_sq)) + ewu/(sigma_sq))) *
    (1 - ewu/(sigma_sq)) * ewv/sigmastar + (2 - 2 * (sigx13^2 *
    ewu * (sigma_sq)/(ssq)^2)) * sigmastar)/(ssq)^2 + S^2 *
    sigx14^2 * ewu * (epsilon)^2/(ssq)) * dmusig) * depsi) *
    (epsilon)/wzdeno - 0.5 * sigx27) + S * ddmu * (epsilon)/wzdeno -
    0.5 * sigx4)/(sigx2 * sqrt(sigma_sq)) - (0.5 * (sigx2/sqrt(sigma_sq)) +
    2 * (ewz * sigx16)) * ewu * sigx16/(sigx2 * sqrt(sigma_sq))^2) *
    ewu * ewz), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (ewv * (S * (((((0.5 *
      ((1 - ewu/(sigma_sq)) * ewv) - S^2 * sigx8 * sigx14 *
      ewu * (epsilon)^2)/(sigma_sq) + 0.5 * ((ewu/(sigma_sq) -
      1) * ewv/(sigma_sq) + 1 - 0.5 * ((1 - ewu/(sigma_sq)) *
      ewvsq))) * dmusig/sigmastar + 0.5 * (S^2 * sigx8 *
      dmusig * (epsilon)^2/(sigma_sq)^2)) * ewu + sigx8 *
      (1 - 2 * (sigx13 * ewu * ssq/(ssq)^2)) * dmusig) *
      depsi/(ssq)^2 + 0.5 * sigx18) * (epsilon)/wzdeno -
      (0.5 * sigx27 + 0.5 * (sigx29/(sigma_sq))))) - 2 *
      (sigx22 * sigx16/sigx2)) * ewu * ewz/(sigx2 * sqrt(sigma_sq)),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx26 * pmusig - S *
      sigx28 * sigx14 * depsi * dmusig * (epsilon)) - 2 *
      (sigx20 * ewz * sigx16/(sigx2 * sqrt(sigma_sq)))) *
      ewu * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * (S^2 * dwsr *
      (epsilon)^2/ewvsr^2) - 0.5 * dwsr))/ewvsr + 2 * (ewv *
      ewz * (S * ((((0.5 * (ewv/(sigma_sq)) - 0.5 * (0.5 *
      ewvsq + ewv/(sigma_sq))) * ewvsq + S^2 * sigx8^2 *
      ewu * ewv * (epsilon)^2/((ssq)^2 * (sigma_sq))) *
      dmusig * ewu/sigmastar + ((0.5 * (S^2 * dmusig *
      (epsilon)^2/(sigma_sq)^2) - 2 * (sigx8 * dmusig *
      ssq/(ssq)^2)) * ewv + dmusig) * sigx8) * depsi *
      ewu/(ssq)^2 + S * (0.5 * (ewv * (S * (sigx8 * depsi *
      ewu/(ssq)^2 + 0.5 * sigx17) * (epsilon) - 2 * (pmusig/(sigma_sq)))) +
      0.5 * pmusig) * dmusig * (epsilon)/(sigma_sq)^2) *
      (epsilon)/wzdeno - ((0.5 * (dmusig * pmusig) + 0.5 *
      (ewv * (S * (sigx10 + 0.5 * sigx3) * (epsilon) -
        wzdeno^2 * dmusig * pmusig/wzsq))) * wzdeno/wzsq +
      0.5 * (ewv * sigx29/(sigma_sq))))/sqrt(sigma_sq)) -
      sigx22^2/sigx2)/sigx2, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((sigx26 * pmusig + S * sigx8 * sigx28 * depsi *
      dmusig * ewu * (epsilon)/(ssq)^2) * ewv) - (sigx22 *
      sigx20/sigx2 + (0.5 * (S^2 * dwsr * (epsilon)^2/(wzdeno *
      ewvsr^3)) - 0.5 * (wzdeno * dwsr * ewvsr/(wzdeno *
      ewvsr)^2)) * prC)) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) +
      ewvsr/(wzdeno * ewvsr)^2) * dwsr - (sigx20^2/sigx2 +
      2 * ((2 - 2 * (wzdeno * (sigma_sq) * ewz/wzsq)) *
        dmusig * pmusig * sqrt(sigma_sq)/wzsq))) * ewz +
      2 * sigx19 - prC * dwsr/(wzdeno * ewvsr)) * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf halfnormal-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
zisfhalfnormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfhalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initHalf <- start_st$initHalf
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfhalfnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfhalfnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfhalfnormlike, grad = cgradzisfhalfnormlike,
      hess = chesszisfhalfnormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfhalfnormlike(mleObj$par,
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
      mleObj$hessian <- chesszisfhalfnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfhalfnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfhalfnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfhalfnormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initHalf = initHalf))
}

# Conditional efficiencies estimation ----------

czisfhalfnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
    exp(Wv))) * pnorm(mustar/sigmastar)
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
  return(res)
}

# Marginal effects on inefficiencies ----------

czisfmarghalfnorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
                                                                     exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarghalfnorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
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
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  Pi1 <- 2/sqrt(exp(Wu) + exp(Wv)) * dnorm(object$S * epsilon/sqrt(exp(Wu) +
                                                                     exp(Wv))) * pnorm(mustar/sigmastar)
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
