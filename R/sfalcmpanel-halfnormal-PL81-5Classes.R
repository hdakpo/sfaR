################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 4L                                                        #
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
pLCMhalfnormlike5C_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, Zvar, nZHvar, wHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)]
  beta5 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)]
  delta5 <- parm[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar)]
  phi5 <- parm[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar)]
  theta1 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar)]
  theta2 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)]
  theta4 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * S * epsilon_i1/(exp(Wv1) + TT * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + TT * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * S * epsilon_i2/(exp(Wv2) + TT * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + TT * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_i3 <- as.numeric(tapply(epsilon_it3, pindex[, 1], sum))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * S * epsilon_i3/(exp(Wv3) + TT * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + TT * exp(Wu3)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon_it4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_i4 <- as.numeric(tapply(epsilon_it4, pindex[, 1], sum))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  mustar4 <- -exp(Wu4) * S * epsilon_i4/(exp(Wv4) + TT * exp(Wu4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wv4) + TT * exp(Wu4)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  epsilon_it5 <- Yvar - as.numeric(crossprod(matrix(beta5), t(Xvar)))
  epsilon_i5 <- as.numeric(tapply(epsilon_it5, pindex[, 1], sum))
  epsilon_isq5 <- as.numeric(tapply(epsilon_it5^2, pindex[, 1], sum))
  mustar5 <- -exp(Wu5) * S * epsilon_i5/(exp(Wv5) + TT * exp(Wu5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wv5) + TT * exp(Wu5)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Pi4 <- 2 * sigmastar4 * exp(-1/2 * (epsilon_isq4/exp(Wv4) - (mustar4/sigmastar4)^2)) *
    pnorm(mustar4/sigmastar4)/((2 * pi)^(TT/2) * exp(Wv4/2 * TT) * exp(Wu4/2))
  Pi5 <- 2 * sigmastar5 * exp(-1/2 * (epsilon_isq5/exp(Wv5) - (mustar5/sigmastar5)^2)) *
    pnorm(mustar5/sigmastar5)/((2 * pi)^(TT/2) * exp(Wv5/2 * TT) * exp(Wu5/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  L <- Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * Pi4 + Probc5 * Pi5
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
psLCMfhalfnorm5C_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
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
    0.98 * Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], 0.98 *
      Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], 0.98 *
      Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], 0.98 *
      Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar +
      2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, 4 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), paste0("Cl1_", colnames(Zvar)), paste0("Cl2_",
      colnames(Zvar)), paste0("Cl3_", colnames(Zvar)), paste0("Cl4_", colnames(Zvar)))
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
pgradLCMhalfnormlike5C_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar +
    3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)]
  beta5 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 4 * nuZUvar +
    4 * nvZVvar)]
  delta5 <- parm[(5 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    4 * nvZVvar)]
  phi5 <- parm[(5 * nXvar + 5 * nuZUvar + 4 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar)]
  theta1 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 1):(5 * nXvar + 5 * nuZUvar +
    5 * nvZVvar + nZHvar)]
  theta2 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + nZHvar + 1):(5 * nXvar +
    5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 2 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar)]
  theta4 <- parm[(5 * nXvar + 5 * nuZUvar + 5 * nvZVvar + 3 * nZHvar + 1):(5 *
    nXvar + 5 * nuZUvar + 5 * nvZVvar + 4 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_i3 <- as.numeric(tapply(epsilon_it3, pindex[, 1], sum))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  epsilon_it4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_i4 <- as.numeric(tapply(epsilon_it4, pindex[, 1], sum))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  epsilon_it5 <- Yvar - as.numeric(crossprod(matrix(beta5), t(Xvar)))
  epsilon_i5 <- as.numeric(tapply(epsilon_it5, pindex[, 1], sum))
  epsilon_isq5 <- as.numeric(tapply(epsilon_it5^2, pindex[, 1], sum))
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it1 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it1, FUN = "*")
  Xepsi_i1 <- apply(Xepsi_it1, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it2 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it2, FUN = "*")
  Xepsi_i2 <- apply(Xepsi_it2, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it3 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it3, FUN = "*")
  Xepsi_i3 <- apply(Xepsi_it3, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it4 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it4, FUN = "*")
  Xepsi_i4 <- apply(Xepsi_it4, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xepsi_it5 <- sweep(-Xvar, MARGIN = 1, STATS = 2 * epsilon_it5, FUN = "*")
  Xepsi_i5 <- apply(Xepsi_it5, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewu3 <- exp(Wu3)
  ewu4 <- exp(Wu4)
  ewu5 <- exp(Wu5)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3 <- exp(Wv3)
  ewv4 <- exp(Wv4)
  ewv5 <- exp(Wv5)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewu3_h <- exp(Wu3/2)
  ewu4_h <- exp(Wu4/2)
  ewu5_h <- exp(Wu5/2)
  ewv1_h <- exp(Wv1 * TT/2)
  ewv2_h <- exp(Wv2 * TT/2)
  ewv3_h <- exp(Wv3 * TT/2)
  ewv4_h <- exp(Wv4 * TT/2)
  ewv5_h <- exp(Wv5 * TT/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  ewz4 <- exp(Wz4)
  sigmasq1 <- sqrt(ewu1 * ewv1/(ewv1 + TT * ewu1))
  sigmasq2 <- sqrt(ewu2 * ewv2/(ewv2 + TT * ewu2))
  sigmasq3 <- sqrt(ewu3 * ewv3/(ewv3 + TT * ewu3))
  sigmasq4 <- sqrt(ewu4 * ewv4/(ewv4 + TT * ewu4))
  sigmasq5 <- sqrt(ewu5 * ewv5/(ewv5 + TT * ewu5))
  wzdeno <- (1 - (ewz1 + ewz2 + ewz3 + ewz4)/(1 + ewz1 + ewz2 + ewz3 + ewz4))
  euv1 <- (1 - ewv1/(ewv1 + TT * ewu1))
  euv2 <- (1 - ewv2/(ewv2 + TT * ewu2))
  euv3 <- (1 - ewv3/(ewv3 + TT * ewu3))
  euv4 <- (1 - ewv4/(ewv4 + TT * ewu4))
  euv5 <- (1 - ewv5/(ewv5 + TT * ewu5))
  sigx1_1 <- ((ewv1 + TT * ewu1) * sigmasq1)
  sigx1_2 <- ((ewv2 + TT * ewu2) * sigmasq2)
  sigx1_3 <- ((ewv3 + TT * ewu3) * sigmasq3)
  sigx1_4 <- ((ewv4 + TT * ewu4) * sigmasq4)
  sigx1_5 <- ((ewv5 + TT * ewu5) * sigmasq5)
  sigx2_1 <- (S * ewu1 * epsilon_i1/sigx1_1)
  sigx2_2 <- (S * ewu2 * epsilon_i2/sigx1_2)
  sigx2_3 <- (S * ewu3 * epsilon_i3/sigx1_3)
  sigx2_4 <- (S * ewu4 * epsilon_i4/sigx1_4)
  sigx2_5 <- (S * ewu5 * epsilon_i5/sigx1_5)
  dmu1 <- dnorm(-sigx2_1, 0, 1)
  dmu2 <- dnorm(-sigx2_2, 0, 1)
  dmu3 <- dnorm(-sigx2_3, 0, 1)
  dmu4 <- dnorm(-sigx2_4, 0, 1)
  dmu5 <- dnorm(-sigx2_5, 0, 1)
  pmu1 <- pnorm(-sigx2_1)
  pmu2 <- pnorm(-sigx2_2)
  pmu3 <- pnorm(-sigx2_3)
  pmu4 <- pnorm(-sigx2_4)
  pmu5 <- pnorm(-sigx2_5)
  piuv1 <- ((2 * pi)^(TT/2) * ewu1_h * ewv1_h)
  piuv2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  piuv3 <- ((2 * pi)^(TT/2) * ewu3_h * ewv3_h)
  piuv4 <- ((2 * pi)^(TT/2) * ewu4_h * ewv4_h)
  piuv5 <- ((2 * pi)^(TT/2) * ewu5_h * ewv5_h)
  sigx3_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-sigx2_1)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-sigx2_2)^2)))
  sigx3_3 <- exp(-(0.5 * (epsilon_isq3/ewv3 - (-sigx2_3)^2)))
  sigx3_4 <- exp(-(0.5 * (epsilon_isq4/ewv4 - (-sigx2_4)^2)))
  sigx3_5 <- exp(-(0.5 * (epsilon_isq5/ewv5 - (-sigx2_5)^2)))
  sigx4_1 <- (sigx3_1 * ewz1 * pmu1 * sigmasq1/piuv1)
  sigx4_2 <- (sigx3_2 * ewz2 * pmu2 * sigmasq2/piuv2)
  sigx4_3 <- (sigx3_3 * ewz3 * pmu3 * sigmasq3/piuv3)
  sigx4_4 <- (sigx3_4 * ewz4 * pmu4 * sigmasq4/piuv4)
  sigx4_5 <- (wzdeno * sigx3_5 * pmu5 * sigmasq5/piuv5)
  sigx5 <- ((2 * sigx4_1 + 2 * sigx4_2 + 2 * sigx4_3 + 2 * sigx4_4)/(1 + ewz1 +
    ewz2 + ewz3 + ewz4) + 2 * sigx4_5)
  sigx6_1 <- (1/sigx1_1 - (0.5 * ((1 - TT * ewu1/(ewv1 + TT * ewu1)) * ewv1/sigmasq1) +
    TT * sigmasq1) * ewu1/sigx1_1^2)
  sigx6_2 <- (1/sigx1_2 - (0.5 * ((1 - TT * ewu2/(ewv2 + TT * ewu2)) * ewv2/sigmasq2) +
    TT * sigmasq2) * ewu2/sigx1_2^2)
  sigx6_3 <- (1/sigx1_3 - (0.5 * ((1 - TT * ewu3/(ewv3 + TT * ewu3)) * ewv3/sigmasq3) +
    TT * sigmasq3) * ewu3/sigx1_3^2)
  sigx6_4 <- (1/sigx1_4 - (0.5 * ((1 - TT * ewu4/(ewv4 + TT * ewu4)) * ewv4/sigmasq4) +
    TT * sigmasq4) * ewu4/sigx1_4^2)
  sigx6_5 <- (1/sigx1_5 - (0.5 * ((1 - TT * ewu5/(ewv5 + TT * ewu5)) * ewv5/sigmasq5) +
    TT * sigmasq5) * ewu5/sigx1_5^2)
  sigx7_1 <- (S * ewu1 * pmu1 * epsilon_i1/sigx1_1 - dmu1)
  sigx7_2 <- (S * ewu2 * pmu2 * epsilon_i2/sigx1_2 - dmu2)
  sigx7_3 <- (S * ewu3 * pmu3 * epsilon_i3/sigx1_3 - dmu3)
  sigx7_4 <- (S * ewu4 * pmu4 * epsilon_i4/sigx1_4 - dmu4)
  sigx7_5 <- (S * ewu5 * pmu5 * epsilon_i5/sigx1_5 - dmu5)
  sigx8_1 <- ((2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmu1 * sigmasq1/piuv1^2)
  sigx8_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmu2 * sigmasq2/piuv2^2)
  sigx8_3 <- ((2 * pi)^(TT/2) * ewu3_h * ewv3_h * pmu3 * sigmasq3/piuv3^2)
  sigx8_4 <- ((2 * pi)^(TT/2) * ewu4_h * ewv4_h * pmu4 * sigmasq4/piuv4^2)
  sigx8_5 <- ((2 * pi)^(TT/2) * ewu5_h * ewv5_h * pmu5 * sigmasq5/piuv5^2)
  sigx9_1 <- ((1 - TT * ewu1/(ewv1 + TT * ewu1)) * ewv1 * pmu1/sigx1_1)
  sigx9_2 <- ((1 - TT * ewu2/(ewv2 + TT * ewu2)) * ewv2 * pmu2/sigx1_2)
  sigx9_3 <- ((1 - TT * ewu3/(ewv3 + TT * ewu3)) * ewv3 * pmu3/sigx1_3)
  sigx9_4 <- ((1 - TT * ewu4/(ewv4 + TT * ewu4)) * ewv4 * pmu4/sigx1_4)
  sigx9_5 <- ((1 - TT * ewu5/(ewv5 + TT * ewu5)) * ewv5 * pmu5/sigx1_5)
  sigx10_1 <- (0.5 * sigx9_1 + S * sigx6_1 * sigx7_1 * sigmasq1 * epsilon_i1)
  sigx10_2 <- (0.5 * sigx9_2 + S * sigx6_2 * sigx7_2 * sigmasq2 * epsilon_i2)
  sigx10_3 <- (0.5 * sigx9_3 + S * sigx6_3 * sigx7_3 * sigmasq3 * epsilon_i3)
  sigx10_4 <- (0.5 * sigx9_4 + S * sigx6_4 * sigx7_4 * sigmasq4 * epsilon_i4)
  sigx10_5 <- (0.5 * sigx9_5 + S * sigx6_5 * sigx7_5 * sigmasq5 * epsilon_i5)
  sigx11_1 <- (sigx10_1 * ewu1/piuv1 - 0.5 * sigx8_1)
  sigx11_2 <- (sigx10_2 * ewu2/piuv2 - 0.5 * sigx8_2)
  sigx11_3 <- (sigx10_3 * ewu3/piuv3 - 0.5 * sigx8_3)
  sigx11_4 <- (sigx10_4 * ewu4/piuv4 - 0.5 * sigx8_4)
  sigx11_5 <- (sigx10_5 * ewu5/piuv5 - 0.5 * sigx8_5)
  sigx12_1 <- (0.5 * (euv1 * ewu1/sigmasq1) + sigmasq1)
  sigx12_2 <- (0.5 * (euv2 * ewu2/sigmasq2) + sigmasq2)
  sigx12_3 <- (0.5 * (euv3 * ewu3/sigmasq3) + sigmasq3)
  sigx12_4 <- (0.5 * (euv4 * ewu4/sigmasq4) + sigmasq4)
  sigx12_5 <- (0.5 * (euv5 * ewu5/sigmasq5) + sigmasq5)
  sigx13_1 <- (S^2 * sigx12_1 * ewu1^2 * ewv1 * epsilon_i1^2/(sigx1_1^2 * (ewv1 +
    TT * ewu1) * sigmasq1))
  sigx13_2 <- (S^2 * sigx12_2 * ewu2^2 * ewv2 * epsilon_i2^2/(sigx1_2^2 * (ewv2 +
    TT * ewu2) * sigmasq2))
  sigx13_3 <- (S^2 * sigx12_3 * ewu3^2 * ewv3 * epsilon_i3^2/(sigx1_3^2 * (ewv3 +
    TT * ewu3) * sigmasq3))
  sigx13_4 <- (S^2 * sigx12_4 * ewu4^2 * ewv4 * epsilon_i4^2/(sigx1_4^2 * (ewv4 +
    TT * ewu4) * sigmasq4))
  sigx13_5 <- (S^2 * sigx12_5 * ewu5^2 * ewv5 * epsilon_i5^2/(sigx1_5^2 * (ewv5 +
    TT * ewu5) * sigmasq5))
  sigx14_1 <- (euv1 * ewu1 * ewv1 * pmu1/sigx1_1)
  sigx14_2 <- (euv2 * ewu2 * ewv2 * pmu2/sigx1_2)
  sigx14_3 <- (euv3 * ewu3 * ewv3 * pmu3/sigx1_3)
  sigx14_4 <- (euv4 * ewu4 * ewv4 * pmu4/sigx1_4)
  sigx14_5 <- (euv5 * ewu5 * ewv5 * pmu5/sigx1_5)
  sigx15_1 <- (TT * (2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmu1 * sigmasq1/piuv1^2)
  sigx15_2 <- (TT * (2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmu2 * sigmasq2/piuv2^2)
  sigx15_3 <- (TT * (2 * pi)^(TT/2) * ewu3_h * ewv3_h * pmu3 * sigmasq3/piuv3^2)
  sigx15_4 <- (TT * (2 * pi)^(TT/2) * ewu4_h * ewv4_h * pmu4 * sigmasq4/piuv4^2)
  sigx15_5 <- (TT * (2 * pi)^(TT/2) * ewu5_h * ewv5_h * pmu5 * sigmasq5/piuv5^2)
  sigx16_1 <- ((2 * sigx13_1 - epsilon_isq1/ewv1) * pmu1)
  sigx16_2 <- ((2 * sigx13_2 - epsilon_isq2/ewv2) * pmu2)
  sigx16_3 <- ((2 * sigx13_3 - epsilon_isq3/ewv3) * pmu3)
  sigx16_4 <- ((2 * sigx13_4 - epsilon_isq4/ewv4) * pmu4)
  sigx16_5 <- ((2 * sigx13_5 - epsilon_isq5/ewv5) * pmu5)
  sigx17_1 <- (S * sigx12_1 * dmu1 * ewu1 * ewv1 * epsilon_i1/sigx1_1^2 - 0.5 *
    sigx16_1)
  sigx17_2 <- (S * sigx12_2 * dmu2 * ewu2 * ewv2 * epsilon_i2/sigx1_2^2 - 0.5 *
    sigx16_2)
  sigx17_3 <- (S * sigx12_3 * dmu3 * ewu3 * ewv3 * epsilon_i3/sigx1_3^2 - 0.5 *
    sigx16_3)
  sigx17_4 <- (S * sigx12_4 * dmu4 * ewu4 * ewv4 * epsilon_i4/sigx1_4^2 - 0.5 *
    sigx16_4)
  sigx17_5 <- (S * sigx12_5 * dmu5 * ewu5 * ewv5 * epsilon_i5/sigx1_5^2 - 0.5 *
    sigx16_5)
  sigx18_1 <- ((sigx17_1 * sigmasq1 + 0.5 * sigx14_1)/piuv1 - 0.5 * sigx15_1)
  sigx18_2 <- ((sigx17_2 * sigmasq2 + 0.5 * sigx14_2)/piuv2 - 0.5 * sigx15_2)
  sigx18_3 <- ((sigx17_3 * sigmasq3 + 0.5 * sigx14_3)/piuv3 - 0.5 * sigx15_3)
  sigx18_4 <- ((sigx17_4 * sigmasq4 + 0.5 * sigx14_4)/piuv4 - 0.5 * sigx15_4)
  sigx18_5 <- ((sigx17_5 * sigmasq5 + 0.5 * sigx14_5)/piuv5 - 0.5 * sigx15_5)
  Xsig1 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu1 * epsilon_i1/(ewv1 + TT *
    ewu1)), FUN = "*")
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu1 * ewu1/sigx1_1, FUN = "*")
  Xsig3 <- sweep((Xepsi_i1 - 2 * Xsig1), MARGIN = 1, STATS = (pmu1/ewv1), FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu2 * epsilon_i2/(ewv2 + TT *
    ewu2)), FUN = "*")
  Xsig5 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu2 * ewu2/sigx1_2, FUN = "*")
  Xsig6 <- sweep((Xepsi_i2 - 2 * Xsig4), MARGIN = 1, STATS = (pmu2/ewv2), FUN = "*")
  Xsig7 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu3 * epsilon_i3/(ewv3 + TT *
    ewu3)), FUN = "*")
  Xsig8 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu3 * ewu3/sigx1_3, FUN = "*")
  Xsig9 <- sweep((Xepsi_i3 - 2 * Xsig7), MARGIN = 1, STATS = (pmu3/ewv3), FUN = "*")
  Xsig10 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu4 * epsilon_i4/(ewv4 + TT *
    ewu4)), FUN = "*")
  Xsig11 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu4 * ewu4/sigx1_4, FUN = "*")
  Xsig12 <- sweep((Xepsi_i4 - 2 * Xsig10), MARGIN = 1, STATS = (pmu4/ewv4), FUN = "*")
  Xsig13 <- sweep(X_iM, MARGIN = 1, STATS = (S^2 * ewu5 * epsilon_i5/(ewv5 + TT *
    ewu5)), FUN = "*")
  Xsig14 <- sweep(X_iM, MARGIN = 1, STATS = S * dmu5 * ewu5/sigx1_5, FUN = "*")
  Xsig15 <- sweep((Xepsi_i5 - 2 * Xsig13), MARGIN = 1, STATS = (pmu5/ewv5), FUN = "*")
  gradll <- cbind(sweep((0.5 * Xsig3 + Xsig2), MARGIN = 1, STATS = -(2 * (sigx3_1 *
    ewz1 * sigmasq1/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4) * (2 * pi)^(TT/2) *
    ewu1_h * ewv1_h))), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx11_1 *
    sigx3_1 * ewz1/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4))), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx18_1 * sigx3_1 * ewz1/(sigx5 * (1 + ewz1 + ewz2 +
      ewz3 + ewz4))), FUN = "*"), sweep((0.5 * Xsig6 + Xsig5), MARGIN = 1,
    STATS = -(2 * (sigx3_2 * ewz2 * sigmasq2/(sigx5 * (1 + ewz1 + ewz2 + ewz3 +
      ewz4) * (2 * pi)^(TT/2) * ewu2_h * ewv2_h))), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (sigx11_2 * sigx3_2 * ewz2/(sigx5 * (1 + ewz1 + ewz2 +
      ewz3 + ewz4))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx18_2 *
    sigx3_2 * ewz2/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4))), FUN = "*"), sweep((0.5 *
    Xsig9 + Xsig8), MARGIN = 1, STATS = -(2 * (sigx3_3 * ewz3 * sigmasq3/(sigx5 *
    (1 + ewz1 + ewz2 + ewz3 + ewz4) * (2 * pi)^(TT/2) * ewu3_h * ewv3_h))), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx11_3 * sigx3_3 * ewz3/(sigx5 *
      (1 + ewz1 + ewz2 + ewz3 + ewz4))), FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = 2 * (sigx18_3 * sigx3_3 * ewz3/(sigx5 * (1 + ewz1 + ewz2 + ewz3 +
        ewz4))), FUN = "*"), sweep((0.5 * Xsig12 + Xsig11), MARGIN = 1, STATS = -(2 *
      (sigx3_4 * ewz4 * sigmasq4/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4) *
        (2 * pi)^(TT/2) * ewu4_h * ewv4_h))), FUN = "*"), sweep(uHvar, MARGIN = 1,
      STATS = 2 * (sigx11_4 * sigx3_4 * ewz4/(sigx5 * (1 + ewz1 + ewz2 + ewz3 +
        ewz4))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx18_4 *
      sigx3_4 * ewz4/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4))), FUN = "*"),
    sweep((0.5 * Xsig15 + Xsig14), MARGIN = 1, STATS = -(2 * (wzdeno * sigx3_5 *
      sigmasq5/(sigx5 * (2 * pi)^(TT/2) * ewu5_h * ewv5_h))), FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = 2 * (sigx11_5 * wzdeno * sigx3_5/sigx5), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx18_5 * wzdeno * sigx3_5/sigx5),
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 * (sigx3_1 * pmu1 * sigmasq1/piuv1) -
      sigx5) * ewz1/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4)), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = (2 * (sigx3_2 * pmu2 * sigmasq2/piuv2) -
      sigx5) * ewz2/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4)), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = (2 * (sigx3_3 * pmu3 * sigmasq3/piuv3) -
      sigx5) * ewz3/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4)), FUN = "*"),
    sweep(Zvar, MARGIN = 1, STATS = (2 * (sigx3_4 * pmu4 * sigmasq4/piuv4) -
      sigx5) * ewz4/(sigx5 * (1 + ewz1 + ewz2 + ewz3 + ewz4)), FUN = "*"))
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
LCM5ChnormAlgOpt_pl81 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, ngZGvar,
  Zvar, nZHvar, Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm5C_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
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
  startLoglik <- sum(pLCMhalfnormlike5C_pl81(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(pLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike5C_pl81,
    grad = pgradLCMhalfnormlike5C_pl81, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, Zvar = Zvar, nZHvar = nZHvar,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike5C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, Zvar = Zvar, nZHvar = nZHvar,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike5C_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike5C_pl81(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike5C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike5C_pl81(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike5C_pl81(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike5C_pl81(parm = mlParam, nXvar = nXvar,
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
#' post. prob. and efficiencies for lcmpanel 5 classes halfnormal-normal distribution
#' @param object object of class lcmpanel
#' @param level level for confidence interval
#' @noRd
pLCM5Chalfnormeff_pl81 <- function(object, level) {
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
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    1):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 4 * object$nuZUvar + 3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 3 * object$nvZVvar +
    1):(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar)]
  beta5 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    1):(5 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar)]
  delta5 <- object$mlParam[(5 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    1):(5 * object$nXvar + 5 * object$nuZUvar + 4 * object$nvZVvar)]
  phi5 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar + 4 * object$nvZVvar +
    1):(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar)]
  theta1 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    1):(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    object$nZHvar + 1):(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    2 * object$nZHvar)]
  theta3 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    2 * object$nZHvar + 1):(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    3 * object$nZHvar)]
  theta4 <- object$mlParam[(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    3 * object$nZHvar + 1):(5 * object$nXvar + 5 * object$nuZUvar + 5 * object$nvZVvar +
    4 * object$nZHvar)]
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
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  epsilon_it1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_i1 <- as.numeric(tapply(epsilon_it1, pindex[, 1], sum))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * object$S * epsilon_i1/(exp(Wv1) + TT * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + TT * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar_p)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar_p)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon_it2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_i2 <- as.numeric(tapply(epsilon_it2, pindex[, 1], sum))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * object$S * epsilon_i2/(exp(Wv2) + TT * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + TT * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon_it3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_i3 <- as.numeric(tapply(epsilon_it3, pindex[, 1], sum))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * object$S * epsilon_i3/(exp(Wv3) + TT * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + TT * exp(Wu3)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar_p)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar_p)))
  Wz4 <- as.numeric(crossprod(matrix(theta4), t(Zvar)))
  epsilon_it4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_i4 <- as.numeric(tapply(epsilon_it4, pindex[, 1], sum))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  mustar4 <- -exp(Wu4) * object$S * epsilon_i4/(exp(Wv4) + TT * exp(Wu4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wv4) + TT * exp(Wu4)))
  Wu5 <- as.numeric(crossprod(matrix(delta5), t(uHvar_p)))
  Wv5 <- as.numeric(crossprod(matrix(phi5), t(vHvar_p)))
  epsilon_it5 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta5), t(Xvar)))
  epsilon_i5 <- as.numeric(tapply(epsilon_it5, pindex[, 1], sum))
  epsilon_isq5 <- as.numeric(tapply(epsilon_it5^2, pindex[, 1], sum))
  mustar5 <- -exp(Wu5) * object$S * epsilon_i5/(exp(Wv5) + TT * exp(Wu5))
  sigmastar5 <- sqrt(exp(Wu5) * exp(Wv5)/(exp(Wv5) + TT * exp(Wu5)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Pi4 <- 2 * sigmastar4 * exp(-1/2 * (epsilon_isq4/exp(Wv4) - (mustar4/sigmastar4)^2)) *
    pnorm(mustar4/sigmastar4)/((2 * pi)^(TT/2) * exp(Wv4/2 * TT) * exp(Wu4/2))
  Pi5 <- 2 * sigmastar5 * exp(-1/2 * (epsilon_isq5/exp(Wv5) - (mustar5/sigmastar5)^2)) *
    pnorm(mustar5/sigmastar5)/((2 * pi)^(TT/2) * exp(Wv5/2 * TT) * exp(Wu5/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc4 <- exp(Wz4)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3) + exp(Wz4))
  Probc5 <- 1 - Probc1 - Probc2 - Probc3 - Probc4
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4 + Pi5 * Probc5)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4 + Pi5 * Probc5)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4 + Pi5 * Probc5)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4 + Pi5 * Probc5)
  Pcond_c5 <- Pi5 * Probc5/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4 + Pi5 * Probc5)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4, Pcond_c5), 1,
    which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c == 2, Pcond_c2, ifelse(Group_c ==
    3, Pcond_c3, ifelse(Group_c == 4, Pcond_c4, Pcond_c5))))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- mustar4 + sigmastar4 * dnorm(mustar4/sigmastar4)/pnorm(mustar4/sigmastar4)
  u_c5 <- mustar5 + sigmastar5 * dnorm(mustar5/sigmastar5)/pnorm(mustar5/sigmastar5)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2, ifelse(Group_c ==
    3, u_c3, ifelse(Group_c == 4, u_c4, u_c5))))
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  ineff_c3 <- ifelse(Group_c == 3, u_c3, NA)
  ineff_c4 <- ifelse(Group_c == 4, u_c4, NA)
  ineff_c5 <- ifelse(Group_c == 5, u_c5, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c3 <- exp(-mustar3 + 1/2 * sigmastar3^2) * pnorm(mustar3/sigmastar3 -
      sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_c4 <- exp(-mustar4 + 1/2 * sigmastar4^2) * pnorm(mustar4/sigmastar4 -
      sigmastar4)/pnorm(mustar4/sigmastar4)
    teBC_c5 <- exp(-mustar5 + 1/2 * sigmastar5^2) * pnorm(mustar5/sigmastar5 -
      sigmastar5)/pnorm(mustar5/sigmastar5)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, ifelse(Group_c == 2, teBC_c2, ifelse(Group_c ==
      3, teBC_c3, ifelse(Group_c == 4, teBC_c4, teBC_c5))))
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    effBC_c3 <- ifelse(Group_c == 3, teBC_c3, NA)
    effBC_c4 <- ifelse(Group_c == 4, teBC_c4, NA)
    effBC_c5 <- ifelse(Group_c == 5, teBC_c5, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c3 <- exp(mustar3 + 1/2 * sigmastar3^2) * pnorm(mustar3/sigmastar3 +
      sigmastar3)/pnorm(mustar3/sigmastar3)
    teBC_reciprocal_c4 <- exp(mustar4 + 1/2 * sigmastar4^2) * pnorm(mustar4/sigmastar4 +
      sigmastar4)/pnorm(mustar4/sigmastar4)
    teBC_reciprocal_c5 <- exp(mustar5 + 1/2 * sigmastar5^2) * pnorm(mustar5/sigmastar5 +
      sigmastar5)/pnorm(mustar5/sigmastar5)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, ifelse(Group_c ==
      2, teBC_reciprocal_c2, ifelse(Group_c == 3, teBC_reciprocal_c3, ifelse(Group_c ==
      4, teBC_reciprocal_c4, teBC_reciprocal_c5))))
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    ReffBC_c3 <- ifelse(Group_c == 3, teBC_reciprocal_c3, NA)
    ReffBC_c4 <- ifelse(Group_c == 4, teBC_reciprocal_c4, NA)
    ReffBC_c5 <- ifelse(Group_c == 5, teBC_reciprocal_c5, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3, u_c3 = u_c3, teBC_c3 = teBC_c3,
      teBC_reciprocal_c3 = teBC_reciprocal_c3, PosteriorProb_c4 = Pcond_c4,
      PriorProb_c4 = Probc4, u_c4 = u_c4, teBC_c4 = teBC_c4, teBC_reciprocal_c4 = teBC_reciprocal_c4,
      PosteriorProb_c5 = Pcond_c5, PriorProb_c5 = Probc5, u_c5 = u_c5, teBC_c5 = teBC_c5,
      teBC_reciprocal_c5 = teBC_reciprocal_c5, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2,
      ineff_c3 = ineff_c3, ineff_c4 = ineff_c4, ineff_c5 = ineff_c5, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, effBC_c3 = effBC_c3, effBC_c4 = effBC_c4, effBC_c5 = effBC_c5,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2, ReffBC_c3 = ReffBC_c3,
      ReffBC_c4 = ReffBC_c4, ReffBC_c5 = ReffBC_c5)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
      u_c3 = u_c3, PosteriorProb_c4 = Pcond_c4, PriorProb_c4 = Probc4, u_c4 = u_c4,
      PosteriorProb_c5 = Pcond_c5, PriorProb_c5 = Probc5, u_c5 = u_c5, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4, ineff_c5 = ineff_c5)
  }
  res <- data.frame(levels(pindex[, 1]), res)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  return(res)
}
