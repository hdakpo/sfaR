################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 2L                                                        #
# Type:      - Pitt and Lee (1981) specification - PL81                        #
#            - Time Invariant Inefficiency                                     #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for halfnormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
pLCMhalfnormlike2C_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, Zvar, nZHvar, wHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * S * epsilon_i1/(exp(Wv1) + TT * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + TT * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * S * epsilon_i2/(exp(Wv2) + TT * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + TT * exp(Wu2)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(-.Machine$double.xmax), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for halfnormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar_c matrix of Zu variables for pooled data
#' @param vHvar_c matrix of Zv variables for pooled data
#' @param uHvar_p matrix of Zu variables for cross-section
#' @param vHvar_p matrix of Zv variables for cross-section
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar_c vector of weights (weighted likelihood) for pooled data
#' @param wHvar_p vector of weights (weighted likelihood) for cross section
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
psLCMfhalfnorm2C_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
  vHvar_c, uHvar_p, vHvar_p, Yvar, Xvar, S, wHvar_c, wHvar_p, Zvar, nZHvar, pindex,
  TT, whichStart, initIter, initAlg, printInfo, tol) {
  initHalf <- psthalfnorm_pl81(olsObj = olsObj, epsiRes = epsiRes, nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_c[, 1, drop = FALSE], vHvar = vHvar_c[,
      1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, initIter = initIter,
    whichStart = whichStart, initAlg = initAlg, tol = tol, printInfo = printInfo)
  if (whichStart == 1L) {
    Esti <- initHalf$StartVal
    initHalfPanel <- NULL
  } else {
    cat("Initialization: SFA Panel PL81-type + halfnormal-normal distribution...\n")
    initHalfPanel <- maxLik::maxLik(logLik = phalfnormlike_pl81, start = initHalf$StartVal,
      grad = pgradhalfnormlike_pl81, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_p[, 1, drop = FALSE], vHvar = vHvar_p[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
    Esti <- initHalfPanel$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)],
    Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
      1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), paste0("Cl1_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalfPanel = initHalfPanel))
}

# Gradient of the likelihood function ----------
#' gradient for halfnormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
pgradLCMhalfnormlike2C_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  Xepsi_it1 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it1, FUN = "*")
  Xepsi_i1 <- apply(Xepsi_it1, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it2 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it2, FUN = "*")
  Xepsi_i2 <- apply(Xepsi_it2, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  ewu1 <- exp(Wu1)
  ewv1 <- exp(Wv1)
  ewu2 <- exp(Wu2)
  ewv2 <- exp(Wv2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1_h <- exp(Wv1 * TT/2)
  ewv2_h <- exp(Wv2 * TT/2)
  ewz <- exp(Wz)
  sigmasq1 <- sqrt(ewu1 * ewv1/(ewv1 + TT * ewu1))
  sigmasq2 <- sqrt(ewu2 * ewv2/(ewv2 + TT * ewu2))
  wzdeno <- (1 - ewz/(1 + ewz))
  euv1 <- (1 - ewv1/(ewv1 + TT * ewu1))
  euv2 <- (1 - ewv2/(ewv2 + TT * ewu2))
  sigx1_1 <- ((ewv1 + TT * ewu1) * sigmasq1)
  sigx1_2 <- ((ewv2 + TT * ewu2) * sigmasq2)
  sigx2_1 <- (S * ewu1 * epsilon_i1/sigx1_1)
  sigx2_2 <- (S * ewu2 * epsilon_i2/sigx1_2)
  dmu1 <- dnorm(-sigx2_1, 0, 1)
  dmu2 <- dnorm(-sigx2_2, 0, 1)
  pmu1 <- pnorm(-sigx2_1)
  pmu2 <- pnorm(-sigx2_2)
  sigx3_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-sigx2_1)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-sigx2_2)^2)))
  sigx4_1 <- ((1 + ewz) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h)
  sigx4_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx5_1 <- (sigx3_1 * ewz * pmu1 * sigmasq1/sigx4_1)
  sigx5_2 <- (wzdeno * sigx3_2 * pmu2 * sigmasq2/sigx4_2)
  sigx6 <- ((1 + ewz) * (2 * sigx5_2 + 2 * sigx5_1) * (2 * pi)^(TT/2) * ewu1_h *
    ewv1_h)
  sigx7_1 <- ((1 - TT * ewu1/(ewv1 + TT * ewu1)) * ewv1 * pmu1/sigx1_1)
  sigx7_2 <- ((1 - TT * ewu2/(ewv2 + TT * ewu2)) * ewv2 * pmu2/sigx1_2)
  sigx8_1 <- (1/sigx1_1 - (0.5 * ((1 - TT * ewu1/(ewv1 + TT * ewu1)) * ewv1/sigmasq1) +
    TT * sigmasq1) * ewu1/sigx1_1^2)
  sigx8_2 <- (1/sigx1_2 - (0.5 * ((1 - TT * ewu2/(ewv2 + TT * ewu2)) * ewv2/sigmasq2) +
    TT * sigmasq2) * ewu2/sigx1_2^2)
  sigx9_1 <- (0.5 * sigx7_1 + S * sigx8_1 * (S * ewu1 * pmu1 * epsilon_i1/sigx1_1 -
    dmu1) * sigmasq1 * epsilon_i1)
  sigx9_2 <- (0.5 * sigx7_2 + S * sigx8_2 * (S * ewu2 * pmu2 * epsilon_i2/sigx1_2 -
    dmu2) * sigmasq2 * epsilon_i2)
  sigx10_1 <- (sigx9_1 * ewu1/sigx4_1 - 0.5 * ((1 + ewz) * (2 * pi)^(TT/2) * ewu1_h *
    ewv1_h * pmu1 * sigmasq1/sigx4_1^2))
  sigx10_2 <- (sigx9_2 * ewu2/sigx4_2 - 0.5 * ((2 * pi)^(TT/2) * ewu2_h * ewv2_h *
    pmu2 * sigmasq2/sigx4_2^2))
  sigx11_1 <- (0.5 * (euv1 * ewu1/sigmasq1) + sigmasq1)
  sigx11_2 <- (0.5 * (euv2 * ewu2/sigmasq2) + sigmasq2)
  sigx12_1 <- (sigx1_1^2 * (ewv1 + TT * ewu1) * sigmasq1)
  sigx12_2 <- (sigx1_2^2 * (ewv2 + TT * ewu2) * sigmasq2)
  sigx13_1 <- (S^2 * sigx11_1 * ewu1^2 * ewv1 * epsilon_i1^2/sigx12_1)
  sigx13_2 <- (S^2 * sigx11_2 * ewu2^2 * ewv2 * epsilon_i2^2/sigx12_2)
  sigx14_1 <- (euv1 * ewu1 * ewv1 * pmu1/sigx1_1)
  sigx14_2 <- (euv2 * ewu2 * ewv2 * pmu2/sigx1_2)
  sigx15_1 <- (S * sigx11_1 * dmu1 * ewu1 * ewv1 * epsilon_i1/sigx1_1^2 - 0.5 *
    ((2 * sigx13_1 - epsilon_isq1/ewv1) * pmu1))
  sigx15_2 <- (S * sigx11_2 * dmu2 * ewu2 * ewv2 * epsilon_i2/sigx1_2^2 - 0.5 *
    ((2 * sigx13_2 - epsilon_isq2/ewv2) * pmu2))
  sigx16_1 <- (TT * (1 + ewz) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmu1 * sigmasq1/sigx4_1^2)
  sigx16_2 <- (TT * (2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmu2 * sigmasq2/sigx4_2^2)
  sigx17_1 <- ((sigx15_1 * sigmasq1 + 0.5 * sigx14_1)/sigx4_1 - 0.5 * sigx16_1)
  sigx17_2 <- ((sigx15_2 * sigmasq2 + 0.5 * sigx14_2)/sigx4_2 - 0.5 * sigx16_2)
  sigx18 <- ((2 * sigx5_2 + 2 * sigx5_1) * (2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx19 <- (wzdeno * sigx3_2 * pmu2 * sigmasq2/((1 + ewz) * (2 * pi)^(TT/2) *
    ewu2_h * ewv2_h))
  sigx20 <- (2 * ((1/sigx4_1 - (2 * pi)^(TT/2) * ewu1_h * ewv1_h * ewz/sigx4_1^2) *
    sigx3_1 * pmu1 * sigmasq1) - 2 * sigx19)
  Xsig1 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu1 * epsilon_i1/(ewv1 + TT *
    ewu1)), FUN = "*")
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu1 * ewu1/sigx1_1, FUN = "*")
  Xsig3 <- sweep((Xepsi_i1 - 2 * Xsig1), MARGIN = 1, STATS = (pmu1/ewv1), FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu2 * epsilon_i2/(ewv2 + TT *
    ewu2)), FUN = "*")
  Xsig5 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu2 * ewu2/sigx1_2, FUN = "*")
  Xsig6 <- sweep((Xepsi_i2 - 2 * Xsig4), MARGIN = 1, STATS = (pmu2/ewv2), FUN = "*")
  gradll <- cbind(sweep((0.5 * Xsig3 + Xsig2), MARGIN = 1, STATS = -(2 * (sigx3_1 *
    ewz * sigmasq1/sigx6)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx10_1 * sigx3_1 * ewz/(2 * sigx5_2 + 2 * sigx5_1)), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx17_1 * sigx3_1 * ewz/(2 * sigx5_2 + 2 * sigx5_1)),
    FUN = "*"), sweep((0.5 * Xsig6 + Xsig5), MARGIN = 1, STATS = -(2 * (wzdeno *
    sigx3_2 * sigmasq2/sigx18)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx10_2 * wzdeno * sigx3_2/(2 * sigx5_2 + 2 * sigx5_1)), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx17_2 * wzdeno * sigx3_2/(2 * sigx5_2 + 2 * sigx5_1)),
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx20 * ewz/(2 * sigx5_2 + 2 *
    sigx5_1), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for halfnormal-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar_c matrix of Zu variables for pooled data
#' @param vHvar_c matrix of Zv variables for pooled data
#' @param uHvar_p matrix of Zu variables for cross-section
#' @param vHvar_p matrix of Zv variables for cross-section
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar_c vector of weights (weighted likelihood) pooled data
#' @param wHvar_p vector of weights (weighted likelihood) cross-section
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
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
LCM2ChnormAlgOpt_pl81 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax, stepmax,
  tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm2C_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c,
      vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, S = S, Zvar = Zvar, nZHvar = nZHvar, pindex = pindex, TT = TT,
      wHvar_c = wHvar_c, wHvar_p = wHvar_p, initIter = initIter, initAlg = initAlg,
      whichStart = whichStart, tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pLCMhalfnormlike2C_pl81(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("LCM Panel PL81-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike2C_pl81,
    grad = pgradLCMhalfnormlike2C_pl81, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, Zvar = Zvar, nZHvar = nZHvar,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike2C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, Zvar = Zvar, nZHvar = nZHvar,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike2C_pl81(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
      Zvar = Zvar, nZHvar = nZHvar))
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
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike2C_pl81(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike2C_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
    nZHvar = nZHvar)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitHalf = initHalf))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmpanel 2 classes halfnormal-normal distribution
#' @param object object of class lcmpanel
#' @param level level for confidence interval
#' @noRd
pLCM2Chalfnormeff_pl81 <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  pindex <- object$dataTable[, 1:2]
  invariance <- object$invariance
  if (invariance == 1) {
    uHvar_p <- apply(uHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    vHvar_p <- apply(vHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
  } else {
    if (invariance == 2) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
    } else {
      if (invariance == 3) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
      }
    }
  }
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- apply(Zvar, 2, function(x) tapply(x, pindex[, 1], mean))
  TT <- as.numeric(table(pindex[, 1]))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar_p)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar_p)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon_it1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * object$S * epsilon_i1/(exp(Wv1) + TT * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + TT * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar_p)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar_p)))
  epsilon_it2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * object$S * epsilon_i2/(exp(Wv2) + TT * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + TT * exp(Wu2)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  res <- data.frame(levels(pindex[, 1]), res)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  return(res)
}
