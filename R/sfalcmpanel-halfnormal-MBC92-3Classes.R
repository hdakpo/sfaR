################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 3L                                                        #
# Type:      - Battese and Coelli 1992 - alternative (mbc92)                   #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = 1 + eta1 * (t - T) + eta2 * (t - T)^2                  #
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
pLCMhalfnormlike3C_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, Zvar, nZHvar, wHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1_1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1_2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 3):(2 * nXvar + nuZUvar + nvZVvar +
    2)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 3):(2 * nXvar + 2 * nuZUvar +
    nvZVvar + 2)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 3):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2)]
  eta2_1 <- parm[2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 3]
  eta2_2 <- parm[2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 4]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 4)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 4)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 4)]
  eta3_1 <- parm[3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 5]
  eta3_2 <- parm[3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 6]
  theta1 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 7):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar + 6)]
  theta2 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 7):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar + 6)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- unlist(lapply(TT, FUN = function(x) 1 + eta1_1 * rev(-(0:(x - 1))) +
    eta1_2 * (rev(-(0:(x - 1))))^2))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * S * giepsi1/(exp(Wv1) + gisq1 * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + gisq1 * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- unlist(lapply(TT, FUN = function(x) 1 + eta2_1 * rev(-(0:(x - 1))) +
    eta2_2 * (rev(-(0:(x - 1))))^2))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * S * giepsi2/(exp(Wv2) + gisq2 * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + gisq2 * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- unlist(lapply(TT, FUN = function(x) 1 + eta3_1 * rev(-(0:(x - 1))) +
    eta3_2 * (rev(-(0:(x - 1))))^2))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * S * giepsi3/(exp(Wv3) + gisq3 * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + gisq3 * exp(Wu3)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  L <- Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3
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
psLCMfhalfnorm3C_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
  vHvar_c, uHvar_p, vHvar_p, Yvar, Xvar, S, wHvar_c, wHvar_p, Zvar, nZHvar, pindex,
  TT, whichStart, initIter, initAlg, printInfo, tol) {
  initHalf <- psthalfnorm_mbc92(olsObj = olsObj, epsiRes = epsiRes, nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_c[, 1, drop = FALSE], vHvar = vHvar_c[,
      1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, initIter = initIter,
    whichStart = whichStart, initAlg = initAlg, tol = tol, printInfo = printInfo)
  if (whichStart == 1L) {
    Esti <- initHalf$StartVal
    initHalfPanel <- NULL
  } else {
    cat("Initialization: SFA Panel MBC92-type + halfnormal-normal distribution...\n")
    initHalfPanel <- maxLik::maxLik(logLik = phalfnormlike_mbc92, start = initHalf$StartVal,
      grad = pgradhalfnormlike_mbc92, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_p[, 1, drop = FALSE], vHvar = vHvar_p[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
    Esti <- initHalfPanel$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[(nXvar +
    3):(nXvar + 4)], 0.98 * Esti[1:(nXvar)], 0.98 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar -
    1), Esti[(nXvar + 3):(nXvar + 4)], 0.98 * Esti[1:(nXvar)], 0.98 * Esti[nXvar +
    1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.98 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[(nXvar + 3):(nXvar + 4)], rep(0, 2 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), "eta1", "eta2", names(Esti)[1:nXvar], paste0("Zu_",
      colnames(uHvar_p)), paste0("Zv_", colnames(vHvar_p)), "eta1", "eta2",
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)), paste0("Zv_", colnames(vHvar_p)),
    "eta1", "eta2", paste0("Cl1_", colnames(Zvar)), paste0("Cl2_", colnames(Zvar)))
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
pgradLCMhalfnormlike3C_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1_1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1_2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 3):(2 * nXvar + nuZUvar + nvZVvar +
    2)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 3):(2 * nXvar + 2 * nuZUvar +
    nvZVvar + 2)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 3):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2)]
  eta2_1 <- parm[2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 3]
  eta2_2 <- parm[2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 4]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + 4)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 3 * nuZUvar +
    2 * nvZVvar + 4)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 5):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + 4)]
  eta3_1 <- parm[3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 5]
  eta3_2 <- parm[3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 6]
  theta1 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 7):(3 * nXvar + 3 * nuZUvar +
    3 * nvZVvar + nZHvar + 6)]
  theta2 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + nZHvar + 7):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * nZHvar + 6)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- unlist(lapply(TT, FUN = function(x) 1 + eta1_1 * rev(-(0:(x - 1))) +
    eta1_2 * (rev(-(0:(x - 1))))^2))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- unlist(lapply(TT, FUN = function(x) 1 + eta2_1 * rev(-(0:(x - 1))) +
    eta2_2 * (rev(-(0:(x - 1))))^2))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- unlist(lapply(TT, FUN = function(x) 1 + eta3_1 * rev(-(0:(x - 1))) +
    eta3_2 * (rev(-(0:(x - 1))))^2))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  Xepsi_it1 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it1, FUN = "*")
  Xepsi_i1 <- apply(Xepsi_it1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it2 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it2, FUN = "*")
  Xepsi_i2 <- apply(Xepsi_it2, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it3 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it3, FUN = "*")
  Xepsi_i3 <- apply(Xepsi_it3, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit1 <- sweep(-Xvar, MARGIN = 1, STATS = git1, FUN = "*")
  Xgi1 <- apply(Xgit1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit2 <- sweep(-Xvar, MARGIN = 1, STATS = git2, FUN = "*")
  Xgi2 <- apply(Xgit2, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit3 <- sweep(-Xvar, MARGIN = 1, STATS = git3, FUN = "*")
  Xgi3 <- apply(Xgit3, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit1 <- 2 * (Zit * (1 + Zit * (eta1_1 + eta1_2 * Zit)))
  giZi1 <- as.numeric(tapply(gitZit1, pindex[, 1], sum))
  gitZit2 <- 2 * (Zit * (1 + Zit * (eta2_1 + eta2_2 * Zit)))
  giZi2 <- as.numeric(tapply(gitZit2, pindex[, 1], sum))
  gitZit3 <- 2 * (Zit * (1 + Zit * (eta3_1 + eta3_2 * Zit)))
  giZi3 <- as.numeric(tapply(gitZit3, pindex[, 1], sum))
  gitZitsq1 <- 2 * (Zit^2 * (1 + Zit * (eta1_1 + eta1_2 * Zit)))
  giZisq1 <- as.numeric(tapply(gitZitsq1, pindex[, 1], sum))
  gitZitsq2 <- 2 * (Zit^2 * (1 + Zit * (eta2_1 + eta2_2 * Zit)))
  giZisq2 <- as.numeric(tapply(gitZitsq2, pindex[, 1], sum))
  gitZitsq3 <- 2 * (Zit^2 * (1 + Zit * (eta3_1 + eta3_2 * Zit)))
  giZisq3 <- as.numeric(tapply(gitZitsq3, pindex[, 1], sum))
  Zit_epsit1 <- Zit * epsilon_it1
  Zi_epsi1 <- as.numeric(tapply(Zit_epsit1, pindex[, 1], sum))
  Zit_epsit2 <- Zit * epsilon_it2
  Zi_epsi2 <- as.numeric(tapply(Zit_epsit2, pindex[, 1], sum))
  Zit_epsit3 <- Zit * epsilon_it3
  Zi_epsi3 <- as.numeric(tapply(Zit_epsit3, pindex[, 1], sum))
  Zitsq_epsit1 <- Zit^2 * epsilon_it1
  Zisq_epsi1 <- as.numeric(tapply(Zitsq_epsit1, pindex[, 1], sum))
  Zitsq_epsit2 <- Zit^2 * epsilon_it2
  Zisq_epsi2 <- as.numeric(tapply(Zitsq_epsit2, pindex[, 1], sum))
  Zitsq_epsit3 <- Zit^2 * epsilon_it3
  Zisq_epsi3 <- as.numeric(tapply(Zitsq_epsit3, pindex[, 1], sum))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewu3 <- exp(Wu3)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3 <- exp(Wv3)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewu3_h <- exp(Wu3/2)
  ewv1_h <- exp(Wv1 * TT/2)
  ewv2_h <- exp(Wv2 * TT/2)
  ewv3_h <- exp(Wv3 * TT/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  sigmasq1 <- (ewu1 * gisq1 + ewv1)
  sigmasq2 <- (ewu2 * gisq2 + ewv2)
  sigmasq3 <- (ewu3 * gisq3 + ewv3)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigmasq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigmasq2)
  sigmastar3 <- sqrt(ewu3 * ewv3/sigmasq3)
  ssq1 <- (sigmasq1 * sigmastar1)
  ssq2 <- (sigmasq2 * sigmastar2)
  ssq3 <- (sigmasq3 * sigmastar3)
  musig1 <- (S * ewu1 * giepsi1/ssq1)
  musig2 <- (S * ewu2 * giepsi2/ssq2)
  musig3 <- (S * ewu3 * giepsi3/ssq3)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  pmusig3 <- pnorm(-musig3)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  dmusig3 <- dnorm(-musig3, 0, 1)
  wugsq1 <- (1 - ewu1 * gisq1/sigmasq1)
  wugsq2 <- (1 - ewu2 * gisq2/sigmasq2)
  wugsq3 <- (1 - ewu3 * gisq3/sigmasq3)
  ewvsq1 <- (1 - ewv1/sigmasq1)
  ewvsq2 <- (1 - ewv2/sigmasq2)
  ewvsq3 <- (1 - ewv3/sigmasq3)
  wzdeno <- (1 - (ewz1 + ewz2)/(1 + ewz1 + ewz2))
  sigx1_1 <- (S * ewu1 * pmusig1 * giepsi1/ssq1 - dmusig1)
  sigx1_2 <- (S * ewu2 * pmusig2 * giepsi2/ssq2 - dmusig2)
  sigx1_3 <- (S * ewu3 * pmusig3 * giepsi3/ssq3 - dmusig3)
  sigx2_1 <- (ewvsq1 * ewu1/sigmastar1)
  sigx2_2 <- (ewvsq2 * ewu2/sigmastar2)
  sigx2_3 <- (ewvsq3 * ewu3/sigmastar3)
  sigx3_1 <- (wugsq1 * ewv1 * pmusig1/ssq1)
  sigx3_2 <- (wugsq2 * ewv2 * pmusig2/ssq2)
  sigx3_3 <- (wugsq3 * ewv3 * pmusig3/ssq3)
  sigx4_1 <- (0.5 * (wugsq1 * ewv1/sigmastar1) + sigmastar1 * gisq1)
  sigx4_2 <- (0.5 * (wugsq2 * ewv2/sigmastar2) + sigmastar2 * gisq2)
  sigx4_3 <- (0.5 * (wugsq3 * ewv3/sigmastar3) + sigmastar3 * gisq3)
  sigx5_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-musig1)^2)))
  sigx5_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-musig2)^2)))
  sigx5_3 <- exp(-(0.5 * (epsilon_isq3/ewv3 - (-musig3)^2)))
  sigx6_1 <- ((2 * pi)^(TT/2) * ewu1_h * ewv1_h)
  sigx6_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx6_3 <- ((2 * pi)^(TT/2) * ewu3_h * ewv3_h)
  sigx7_1 <- (sigx5_1 * ewz1 * pmusig1 * sigmastar1/sigx6_1)
  sigx7_2 <- (sigx5_2 * ewz2 * pmusig2 * sigmastar2/sigx6_2)
  sigx7_3 <- (wzdeno * sigx5_3 * pmusig3 * sigmastar3/sigx6_3)
  sigx7bis <- ((2 * sigx7_1 + 2 * sigx7_2)/(1 + ewz1 + ewz2) + 2 * sigx7_3)
  sigx8 <- (sigx7bis * (1 + ewz1 + ewz2) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h)
  sigx9_1 <- (0.5 * sigx3_1 + S * (1/ssq1 - sigx4_1 * ewu1/ssq1^2) * sigx1_1 *
    sigmastar1 * giepsi1)
  sigx9_2 <- (0.5 * sigx3_2 + S * (1/ssq2 - sigx4_2 * ewu2/ssq2^2) * sigx1_2 *
    sigmastar2 * giepsi2)
  sigx9_3 <- (0.5 * sigx3_3 + S * (1/ssq3 - sigx4_3 * ewu3/ssq3^2) * sigx1_3 *
    sigmastar3 * giepsi3)
  sigx10_1 <- ((2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmusig1 * sigmastar1/sigx6_1^2)
  sigx10_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmusig2 * sigmastar2/sigx6_2^2)
  sigx10_3 <- ((2 * pi)^(TT/2) * ewu3_h * ewv3_h * pmusig3 * sigmastar3/sigx6_3^2)
  sigx11_1 <- (sigx9_1 * ewu1/sigx6_1 - 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * ewu2/sigx6_2 - 0.5 * sigx10_2)
  sigx11_3 <- (sigx9_3 * ewu3/sigx6_3 - 0.5 * sigx10_3)
  sigx12_1 <- (S^2 * (0.5 * sigx2_1 + sigmastar1) * ewu1^2 * ewv1 * giepsi1^2/(ssq1^2 *
    sigmasq1 * sigmastar1))
  sigx12_2 <- (S^2 * (0.5 * sigx2_2 + sigmastar2) * ewu2^2 * ewv2 * giepsi2^2/(ssq2^2 *
    sigmasq2 * sigmastar2))
  sigx12_3 <- (S^2 * (0.5 * sigx2_3 + sigmastar3) * ewu3^2 * ewv3 * giepsi3^2/(ssq3^2 *
    sigmasq3 * sigmastar3))
  sigx13_1 <- ((2 * sigx12_1 - epsilon_isq1/ewv1) * pmusig1)
  sigx13_2 <- ((2 * sigx12_2 - epsilon_isq2/ewv2) * pmusig2)
  sigx13_3 <- ((2 * sigx12_3 - epsilon_isq3/ewv3) * pmusig3)
  sigx14_1 <- (ewvsq1 * ewu1 * ewv1 * pmusig1/ssq1)
  sigx14_2 <- (ewvsq2 * ewu2 * ewv2 * pmusig2/ssq2)
  sigx14_3 <- (ewvsq3 * ewu3 * ewv3 * pmusig3/ssq3)
  sigx15_1 <- (TT * (2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmusig1 * sigmastar1/sigx6_1^2)
  sigx15_2 <- (TT * (2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmusig2 * sigmastar2/sigx6_2^2)
  sigx15_3 <- (TT * (2 * pi)^(TT/2) * ewu3_h * ewv3_h * pmusig3 * sigmastar3/sigx6_3^2)
  sigx16_1 <- (S * (0.5 * sigx2_1 + sigmastar1) * dmusig1 * ewu1 * ewv1 * giepsi1/ssq1^2 -
    0.5 * sigx13_1)
  sigx16_2 <- (S * (0.5 * sigx2_2 + sigmastar2) * dmusig2 * ewu2 * ewv2 * giepsi2/ssq2^2 -
    0.5 * sigx13_2)
  sigx16_3 <- (S * (0.5 * sigx2_3 + sigmastar3) * dmusig3 * ewu3 * ewv3 * giepsi3/ssq3^2 -
    0.5 * sigx13_3)
  sigx17_1 <- ((sigx16_1 * sigmastar1 + 0.5 * sigx14_1)/sigx6_1 - 0.5 * sigx15_1)
  sigx17_2 <- ((sigx16_2 * sigmastar2 + 0.5 * sigx14_2)/sigx6_2 - 0.5 * sigx15_2)
  sigx17_3 <- ((sigx16_3 * sigmastar3 + 0.5 * sigx14_3)/sigx6_3 - 0.5 * sigx15_3)
  sigx18_1 <- (sigmastar1 - 0.5 * (ewu1 * ewv1/ssq1))
  sigx18_2 <- (sigmastar2 - 0.5 * (ewu2 * ewv2/ssq2))
  sigx18_3 <- (sigmastar3 - 0.5 * (ewu3 * ewv3/ssq3))
  sigx19_1 <- (Zisq_epsi1/ssq1 - ewu1 * sigx18_1 * giepsi1 * giZisq1/ssq1^2)
  sigx19_2 <- (Zisq_epsi2/ssq2 - ewu2 * sigx18_2 * giepsi2 * giZisq2/ssq2^2)
  sigx19_3 <- (Zisq_epsi3/ssq3 - ewu3 * sigx18_3 * giepsi3 * giZisq3/ssq3^2)
  sigx20_1 <- (ewu1 * ewv1 * pmusig1 * giZisq1/(sigmasq1^2 * sigmastar1))
  sigx20_2 <- (ewu2 * ewv2 * pmusig2 * giZisq2/(sigmasq2^2 * sigmastar2))
  sigx20_3 <- (ewu3 * ewv3 * pmusig3 * giZisq3/(sigmasq3^2 * sigmastar3))
  sigx21_1 <- (S * sigx1_1 * sigmastar1 * sigx19_1 - 0.5 * sigx20_1)
  sigx21_2 <- (S * sigx1_2 * sigmastar2 * sigx19_2 - 0.5 * sigx20_2)
  sigx21_3 <- (S * sigx1_3 * sigmastar3 * sigx19_3 - 0.5 * sigx20_3)
  sigx22_1 <- (Zi_epsi1/ssq1 - ewu1 * sigx18_1 * giepsi1 * giZi1/ssq1^2)
  sigx22_2 <- (Zi_epsi2/ssq2 - ewu2 * sigx18_2 * giepsi2 * giZi2/ssq2^2)
  sigx22_3 <- (Zi_epsi3/ssq3 - ewu3 * sigx18_3 * giepsi3 * giZi3/ssq3^2)
  sigx23_1 <- (ewu1 * ewv1 * pmusig1 * giZi1/(sigmasq1^2 * sigmastar1))
  sigx23_2 <- (ewu2 * ewv2 * pmusig2 * giZi2/(sigmasq2^2 * sigmastar2))
  sigx23_3 <- (ewu3 * ewv3 * pmusig3 * giZi3/(sigmasq3^2 * sigmastar3))
  sigx24_1 <- (S * sigx1_1 * sigmastar1 * sigx22_1 - 0.5 * sigx23_1)
  sigx24_2 <- (S * sigx1_2 * sigmastar2 * sigx22_2 - 0.5 * sigx23_2)
  sigx24_3 <- (S * sigx1_3 * sigmastar3 * sigx22_3 - 0.5 * sigx23_3)
  sigx25 <- (sigx7bis * (1 + ewz1 + ewz2) * (2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx26 <- (sigx7bis * (2 * pi)^(TT/2) * ewu3_h * ewv3_h)
  sigx27 <- (sigx5_1 * pmusig1 * sigmastar1/sigx6_1)
  sigx28 <- (sigx5_2 * pmusig2 * sigmastar2/sigx6_2)
  Xsig1 <- sweep(Xgi1, MARGIN = 1, STATS = (S^2 * ewu1 * giepsi1/sigmasq1), FUN = "*")
  Xsig2 <- sweep(Xgi2, MARGIN = 1, STATS = (S^2 * ewu2 * giepsi2/sigmasq2), FUN = "*")
  Xsig3 <- sweep(Xgi1, MARGIN = 1, STATS = S * dmusig1 * ewu1/ssq1, FUN = "*")
  Xsig4 <- sweep(Xgi2, MARGIN = 1, STATS = S * dmusig2 * ewu2/ssq2, FUN = "*")
  Xsig5 <- sweep((Xepsi_i1 - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv1), FUN = "*")
  Xsig6 <- sweep((Xepsi_i2 - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv2), FUN = "*")
  Xsig7 <- sweep(Xgi3, MARGIN = 1, STATS = (S^2 * ewu3 * giepsi3/sigmasq3), FUN = "*")
  Xsig8 <- sweep(Xgi3, MARGIN = 1, STATS = S * dmusig3 * ewu3/ssq3, FUN = "*")
  Xsig9 <- sweep((Xepsi_i3 - 2 * Xsig7), MARGIN = 1, STATS = (pmusig3/ewv3), FUN = "*")
  gradll <- cbind(sweep((0.5 * Xsig5 + Xsig3), MARGIN = 1, STATS = -(2 * (sigx5_1 *
    ewz1 * sigmastar1/sigx8)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx11_1 * sigx5_1 * ewz1/(sigx7bis * (1 + ewz1 + ewz2))), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx17_1 * sigx5_1 * ewz1/(sigx7bis * (1 + ewz1 +
      ewz2))), FUN = "*"), 2 * (sigx5_1 * ewu1 * ewz1 * sigx24_1/sigx8), 2 *
    (sigx5_1 * ewu1 * ewz1 * sigx21_1/sigx8), sweep((0.5 * Xsig6 + Xsig4), MARGIN = 1,
    STATS = -(2 * (sigx5_2 * ewz2 * sigmastar2/sigx25)), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = 2 * (sigx11_2 * sigx5_2 * ewz2/(sigx7bis * (1 + ewz1 +
      ewz2))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx17_2 *
    sigx5_2 * ewz2/(sigx7bis * (1 + ewz1 + ewz2))), FUN = "*"), 2 * (sigx5_2 *
    ewu2 * ewz2 * sigx24_2/sigx25), 2 * (sigx5_2 * ewu2 * ewz2 * sigx21_2/sigx25),
    sweep((0.5 * Xsig9 + Xsig8), MARGIN = 1, STATS = -(2 * (wzdeno * sigx5_3 *
      sigmastar3/sigx26)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
      (sigx11_3 * wzdeno * sigx5_3/sigx7bis), FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = 2 * (sigx17_3 * wzdeno * sigx5_3/sigx7bis), FUN = "*"), 2 * (wzdeno *
      sigx5_3 * ewu3 * sigx24_3/sigx26), 2 * (wzdeno * sigx5_3 * ewu3 * sigx21_3/sigx26),
    sweep(Zvar, MARGIN = 1, STATS = (2 * sigx27 - sigx7bis) * ewz1/(sigx7bis *
      (1 + ewz1 + ewz2)), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 *
      sigx28 - sigx7bis) * ewz2/(sigx7bis * (1 + ewz1 + ewz2)), FUN = "*"))
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
#' @param Zvar separating variables
#' @param nZHvar number of separating variables
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
LCM3ChnormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax, stepmax,
  tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm3C_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
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
  startLoglik <- sum(pLCMhalfnormlike3C_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("LCM Panel MBC92-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike3C_mbc92,
    grad = pgradLCMhalfnormlike3C_mbc92, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, Zvar = Zvar, nZHvar = nZHvar,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike3C_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, Zvar = Zvar, nZHvar = nZHvar,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike3C_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike3C_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike3C_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike3C_mbc92(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike3C_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike3C_mbc92(parm = mlParam, nXvar = nXvar,
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
#' post. prob. and efficiencies for lcmpanel 3 classes halfnormal-normal distribution
#' @param object object of class lcmpanel
#' @param level level for confidence interval
#' @noRd
pLCM3Chalfnormeff_mbc92 <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta1_1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta1_2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 2]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 3):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar + 2)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    3):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar + 2)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    3):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + 2)]
  eta2_1 <- object$mlParam[2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    3]
  eta2_2 <- object$mlParam[2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    4]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    5):(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + 4)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    5):(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar + 4)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    5):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar + 4)]
  eta3_1 <- object$mlParam[3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    5]
  eta3_2 <- object$mlParam[3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    6]
  theta1 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    7):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar + object$nZHvar +
    6)]
  theta2 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    object$nZHvar + 7):(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    2 * object$nZHvar + 6)]
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
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- unlist(lapply(TT, FUN = function(x) 1 + eta1_1 * rev(-(0:(x - 1))) +
    eta1_2 * (rev(-(0:(x - 1))))^2))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * object$S * giepsi1/(exp(Wv1) + gisq1 * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + gisq1 * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar_p)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar_p)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  epsilon_it2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- unlist(lapply(TT, FUN = function(x) 1 + eta2_1 * rev(-(0:(x - 1))) +
    eta2_2 * (rev(-(0:(x - 1))))^2))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * object$S * giepsi2/(exp(Wv2) + gisq2 * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + gisq2 * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
  epsilon_it3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- unlist(lapply(TT, FUN = function(x) 1 + eta3_1 * rev(-(0:(x - 1))) +
    eta3_2 * (rev(-(0:(x - 1))))^2))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * object$S * giepsi3/(exp(Wv2) + gisq3 * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + gisq3 * exp(Wu3)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2))
  Probc3 <- 1 - Probc1 - Probc2
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3), 1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c == 2, Pcond_c2, Pcond_c3))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  res <- data.frame(levels(pindex[, 1]), Group_c = Group_c, PosteriorProb_c = P_cond_c,
    PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
    PriorProb_c2 = Probc2, u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
    u_c3 = u_c3)
  if (object$logDepVar == TRUE) {
    res <- data.frame(res, mustar1 = mustar1, mustar2 = mustar2, mustar3 = mustar3,
      sigmastar1 = sigmastar1, sigmastar2 = sigmastar2, sigmastar3 = sigmastar3)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u_c1 <- res$u_c1 * git1
  res$u_c2 <- res$u_c2 * git2
  res$u_c3 <- res$u_c3 * git3
  res$u_c <- ifelse(res$Group_c == 1, res$u_c1, ifelse(res$Group_c == 2, res$u_c2,
    res$u_c3))
  res$ineff_c1 <- ifelse(res$Group_c == 1, res$u_c1, NA)
  res$ineff_c2 <- ifelse(res$Group_c == 2, res$u_c2, NA)
  res$ineff_c3 <- ifelse(res$Group_c == 2, res$u_c3, NA)
  if (object$logDepVar == TRUE) {
    res$teJLMS_c <- exp(-res$u_c)
    res$teBC_c1 <- exp(1/2 * res$sigmastar1^2 * git1^2 - res$mustar1 * git1) *
      pnorm(res$mustar1/res$sigmastar1 - res$sigmastar1 * git1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_c2 <- exp(1/2 * res$sigmastar2^2 * git2^2 - res$mustar2 * git2) *
      pnorm(res$mustar2/res$sigmastar2 - res$sigmastar2 * git2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_c3 <- exp(1/2 * res$sigmastar3^2 * git3^2 - res$mustar3 * git3) *
      pnorm(res$mustar3/res$sigmastar3 - res$sigmastar3 * git3)/pnorm(res$mustar3/res$sigmastar3)
    res$teBC_c <- ifelse(res$Group_c == 1, res$teBC_c1, ifelse(res$Group_c ==
      2, res$teBC_c2, res$teBC_c3))
    res$effBC_c1 <- ifelse(res$Group_c == 1, res$teBC_c1, NA)
    res$effBC_c2 <- ifelse(res$Group_c == 2, res$teBC_c2, NA)
    res$effBC_c3 <- ifelse(res$Group_c == 2, res$teBC_c3, NA)
    res$teBC_reciprocal_c1 <- exp(res$mustar1 + 1/2 * res$sigmastar1^2) * pnorm(res$mustar1/res$sigmastar1 +
      res$sigmastar1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_reciprocal_c2 <- exp(res$mustar2 + 1/2 * res$sigmastar2^2) * pnorm(res$mustar2/res$sigmastar2 +
      res$sigmastar2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_reciprocal_c3 <- exp(res$mustar3 + 1/2 * res$sigmastar3^2) * pnorm(res$mustar3/res$sigmastar3 +
      res$sigmastar3)/pnorm(res$mustar3/res$sigmastar3)
    res$teBC_reciprocal_c <- ifelse(res$Group_c == 1, res$teBC_reciprocal_c1,
      ifelse(res$Group_c == 2, res$teBC_reciprocal_c2, res$teBC_reciprocal_c3))
    res$ReffBC_c1 <- ifelse(res$Group_c == 1, res$teBC_reciprocal_c1, NA)
    res$ReffBC_c2 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c2, NA)
    res$ReffBC_c3 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c3, NA)
    res$mustar1 <- NULL
    res$sigmastar1 <- NULL
    res$mustar2 <- NULL
    res$sigmastar2 <- NULL
    res$mustar3 <- NULL
    res$sigmastar3 <- NULL
  }
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
