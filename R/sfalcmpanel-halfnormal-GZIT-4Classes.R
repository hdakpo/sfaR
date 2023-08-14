################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 4L                                                        #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Battese and Coelli 1992 specification:               #
#                          - g(zit) = exp(-eta * (t - T))                      #
#                         Cuesta and Orea (2002), Feng and Serletis (2009)     #
#                          - g(zit) = exp(-eta1 * (t - T) - eta2 * (t - T)^2)  #
#                         Alvarez, Amsler, Orea, Schmidt (2006)                #
#                          - g(zit) = exp(eta * gHvar)                         #
#                         Kumbhakar and Wang 2005 specification:               #
#                          - g(zit) = exp(eta * (t - t1))                      #
#                         Cuesta 2000 specification:                           #
#                          - g(zit) = exp(-eta_i * (t - T))                    #
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
pLCMhalfnormlike4C_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, Zvar, nZHvar, ngZGvar, gHvar, wHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar + ngZGvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar + ngZGvar)]
  eta2 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + ngZGvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * ngZGvar)]
  eta3 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * ngZGvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 3 * ngZGvar)]
  eta4 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * ngZGvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 2 * nZHvar +
    1):((4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 3 * nZHvar))]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- exp(as.numeric(crossprod(matrix(eta1), t(gHvar))))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * S * giepsi1/(exp(Wv1) + gisq1 * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + gisq1 * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- exp(as.numeric(crossprod(matrix(eta2), t(gHvar))))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * S * giepsi2/(exp(Wv2) + gisq2 * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + gisq2 * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- exp(as.numeric(crossprod(matrix(eta3), t(gHvar))))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * S * giepsi3/(exp(Wv3) + gisq3 * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + gisq3 * exp(Wu3)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  epsilon_it4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  git4 <- exp(as.numeric(crossprod(matrix(eta4), t(gHvar))))
  git_epsit4 <- epsilon_it4 * git4
  giepsi4 <- as.numeric(tapply(git_epsit4, pindex[, 1], sum))
  gisq4 <- as.numeric(tapply(git4^2, pindex[, 1], sum))
  mustar4 <- -exp(Wu4) * S * giepsi4/(exp(Wv4) + gisq4 * exp(Wu4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wv4) + gisq4 * exp(Wu4)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Pi4 <- 2 * sigmastar4 * exp(-1/2 * (epsilon_isq4/exp(Wv4) - (mustar4/sigmastar4)^2)) *
    pnorm(mustar4/sigmastar4)/((2 * pi)^(TT/2) * exp(Wv4/2 * TT) * exp(Wu4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  L <- Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * Pi4
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
#' @param modelType specification of inefficiency model G(t)u_i
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar_c vector of weights (weighted likelihood) for pooled data
#' @param wHvar_p vector of weights (weighted likelihood) for cross section
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
psLCMfhalfnorm4C_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
  vHvar_c, uHvar_p, vHvar_p, Yvar, Xvar, ngZGvar, gHvar, S, wHvar_c, wHvar_p, Zvar,
  nZHvar, modelType, pindex, TT, whichStart, initIter, initAlg, printInfo, tol) {
  initHalf <- psthalfnorm_gzit(olsObj = olsObj, epsiRes = epsiRes, nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, ngZGvar = ngZGvar, gHvar = gHvar, uHvar = uHvar_c[,
      1, drop = FALSE], vHvar = vHvar_c[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar_c, modelType = modelType, initIter = initIter, whichStart = whichStart,
    initAlg = initAlg, tol = tol, printInfo = printInfo)
  if (whichStart == 1L) {
    Esti <- initHalf$StartVal
    initHalfPanel <- NULL
  } else {
    cat("Initialization: SFA Panel BC92-type + halfnormal-normal distribution...\n")
    initHalfPanel <- maxLik::maxLik(logLik = phalfnormlike_gzit, start = initHalf$StartVal,
      grad = pgradhalfnormlike_gzit, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_p[, 1, drop = FALSE], vHvar = vHvar_p[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar)
    Esti <- initHalfPanel$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[(nXvar +
    3):(nXvar + ngZGvar + 2)], 0.98 * Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar -
    1), Esti[(nXvar + 3):(nXvar + ngZGvar + 2)], 0.98 * Esti[1:(nXvar)], Esti[nXvar +
    1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[(nXvar + 3):(nXvar + ngZGvar + 2)], 0.98 * Esti[1:(nXvar)],
    Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
      1) rep(0, nvZVvar - 1), Esti[(nXvar + 3):(nXvar + ngZGvar + 2)], rep(0,
      3 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), if (modelType %in% c("bc92a", "kw05")) {
      "eta"
    } else {
      if (modelType == "bc92b") {
        c("eta1", "eta2")
      } else {
        if (modelType == "bc92c") {
          paste0("Zg_", colnames(gHvar))
        } else {
          if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
          }
        }
      }
    }, names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)), paste0("Zv_",
      colnames(vHvar_p)), if (modelType %in% c("bc92a", "kw05")) {
      "eta"
    } else {
      if (modelType == "bc92b") {
        c("eta1", "eta2")
      } else {
        if (modelType == "bc92c") {
          paste0("Zg_", colnames(gHvar))
        } else {
          if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
          }
        }
      }
    }, names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)), paste0("Zv_",
      colnames(vHvar_p)), if (modelType %in% c("bc92a", "kw05")) {
      "eta"
    } else {
      if (modelType == "bc92b") {
        c("eta1", "eta2")
      } else {
        if (modelType == "bc92c") {
          paste0("Zg_", colnames(gHvar))
        } else {
          if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
          }
        }
      }
    }, names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)), paste0("Zv_",
      colnames(vHvar_p)), if (modelType %in% c("bc92a", "kw05")) {
      "eta"
    } else {
      if (modelType == "bc92b") {
        c("eta1", "eta2")
      } else {
        if (modelType == "bc92c") {
          paste0("Zg_", colnames(gHvar))
        } else {
          if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
          }
        }
      }
    }, paste0("Cl1_", colnames(Zvar)), paste0("Cl2_", colnames(Zvar)), paste0("Cl3_",
      colnames(Zvar)))
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
pgradLCMhalfnormlike4C_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, Zvar, nZHvar, ngZGvar, gHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + nuZUvar +
    nvZVvar + ngZGvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + 2 *
    nuZUvar + nvZVvar + ngZGvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + ngZGvar + 1):(2 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar + ngZGvar)]
  eta2 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + ngZGvar + 1):(2 * nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 *
    nXvar + 3 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 2 * ngZGvar)]
  eta3 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 2 * ngZGvar + 1):(3 * nXvar +
    3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 *
    nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 *
    nXvar + 4 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 3 * ngZGvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 3 * ngZGvar)]
  eta4 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * ngZGvar + 1):(4 * nXvar +
    4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 1):(4 *
    nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + nZHvar +
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 2 * nZHvar +
    1):((4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 4 * ngZGvar + 3 * nZHvar))]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- exp(as.numeric(crossprod(matrix(eta1), t(gHvar))))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- exp(as.numeric(crossprod(matrix(eta2), t(gHvar))))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  epsilon_it3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- exp(as.numeric(crossprod(matrix(eta3), t(gHvar))))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  epsilon_it4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  git4 <- exp(as.numeric(crossprod(matrix(eta4), t(gHvar))))
  git_epsit4 <- epsilon_it4 * git4
  giepsi4 <- as.numeric(tapply(git_epsit4, pindex[, 1], sum))
  gisq4 <- as.numeric(tapply(git4^2, pindex[, 1], sum))
  Xepsi_it1 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it1, FUN = "*")
  Xepsi_i1 <- apply(Xepsi_it1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit1 <- sweep(-Xvar, MARGIN = 1, STATS = git1, FUN = "*")
  Xgi1 <- apply(Xgit1, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit1 <- sweep(gHvar, MARGIN = 1, STATS = git1 * epsilon_it1, FUN = "*")
  Zigiepsi1 <- apply(Zitgitepsit1, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq1 <- sweep(gHvar, MARGIN = 1, STATS = 2 * git1^2, FUN = "*")
  Zigisq1 <- apply(Zitgitsq1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it2 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it2, FUN = "*")
  Xepsi_i2 <- apply(Xepsi_it2, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit2 <- sweep(-Xvar, MARGIN = 1, STATS = git2, FUN = "*")
  Xgi2 <- apply(Xgit2, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit2 <- sweep(gHvar, MARGIN = 1, STATS = git2 * epsilon_it2, FUN = "*")
  Zigiepsi2 <- apply(Zitgitepsit2, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq2 <- sweep(gHvar, MARGIN = 1, STATS = 2 * git2^2, FUN = "*")
  Zigisq2 <- apply(Zitgitsq2, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it3 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it3, FUN = "*")
  Xepsi_i3 <- apply(Xepsi_it3, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit3 <- sweep(-Xvar, MARGIN = 1, STATS = git3, FUN = "*")
  Xgi3 <- apply(Xgit3, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit3 <- sweep(gHvar, MARGIN = 1, STATS = git3 * epsilon_it3, FUN = "*")
  Zigiepsi3 <- apply(Zitgitepsit3, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq3 <- sweep(gHvar, MARGIN = 1, STATS = 2 * git3^2, FUN = "*")
  Zigisq3 <- apply(Zitgitsq3, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it4 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it4, FUN = "*")
  Xepsi_i4 <- apply(Xepsi_it4, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit4 <- sweep(-Xvar, MARGIN = 1, STATS = git4, FUN = "*")
  Xgi4 <- apply(Xgit4, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit4 <- sweep(gHvar, MARGIN = 1, STATS = git4 * epsilon_it4, FUN = "*")
  Zigiepsi4 <- apply(Zitgitepsit4, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq4 <- sweep(gHvar, MARGIN = 1, STATS = 2 * git4^2, FUN = "*")
  Zigisq4 <- apply(Zitgitsq4, 2, function(x) tapply(x, pindex[, 1], sum))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewu3 <- exp(Wu3)
  ewu4 <- exp(Wu4)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv3 <- exp(Wv3)
  ewv4 <- exp(Wv4)
  ewvTT1 <- exp(Wv1 * TT/2)
  ewvTT2 <- exp(Wv2 * TT/2)
  ewvTT3 <- exp(Wv3 * TT/2)
  ewvTT4 <- exp(Wv4 * TT/2)
  ewusqrt1 <- exp(Wu1/2)
  ewusqrt2 <- exp(Wu2/2)
  ewusqrt3 <- exp(Wu3/2)
  ewusqrt4 <- exp(Wu4/2)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  sigmasq1 <- (ewu1 * gisq1 + ewv1)
  sigmasq2 <- (ewu2 * gisq2 + ewv2)
  sigmasq3 <- (ewu3 * gisq3 + ewv3)
  sigmasq4 <- (ewu4 * gisq4 + ewv4)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigmasq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigmasq2)
  sigmastar3 <- sqrt(ewu3 * ewv3/sigmasq3)
  sigmastar4 <- sqrt(ewu4 * ewv4/sigmasq4)
  ssqx2_1 <- sigmasq1 * sigmastar1
  ssqx2_2 <- sigmasq2 * sigmastar2
  ssqx2_3 <- sigmasq3 * sigmastar3
  ssqx2_4 <- sigmasq4 * sigmastar4
  musig1 <- (S * ewu1 * giepsi1/(ssqx2_1))
  musig2 <- (S * ewu2 * giepsi2/(ssqx2_2))
  musig3 <- (S * ewu3 * giepsi3/(ssqx2_3))
  musig4 <- (S * ewu4 * giepsi4/(ssqx2_4))
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  pmusig3 <- pnorm(-musig3)
  pmusig4 <- pnorm(-musig4)
  dmusig1 <- dnorm(-musig1)
  dmusig2 <- dnorm(-musig2)
  dmusig3 <- dnorm(-musig3)
  dmusig4 <- dnorm(-musig4)
  wzdeno <- (1 + ewz1 + ewz2 + ewz3)
  Tpi <- (2 * pi)^(TT/2)
  wvsq1 <- (1 - ewv1/sigmasq1)
  wvsq2 <- (1 - ewv2/sigmasq2)
  wvsq3 <- (1 - ewv3/sigmasq3)
  wvsq4 <- (1 - ewv4/sigmasq4)
  wzT1 <- (Tpi * ewusqrt1 * ewvTT1)
  wzT2 <- (Tpi * ewusqrt2 * ewvTT2)
  wzT3 <- (Tpi * ewusqrt3 * ewvTT3)
  wzT4 <- (Tpi * ewusqrt4 * ewvTT4)
  wusq1 <- (1 - ewu1 * gisq1/sigmasq1)
  wusq2 <- (1 - ewu2 * gisq2/sigmasq2)
  wusq3 <- (1 - ewu3 * gisq3/sigmasq3)
  wusq4 <- (1 - ewu4 * gisq4/sigmasq4)
  wzlogit <- (1 - (ewz1 + ewz2 + ewz3)/wzdeno)
  sigx1_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-musig1)^2)))
  sigx1_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-musig2)^2)))
  sigx1_3 <- exp(-(0.5 * (epsilon_isq3/ewv3 - (-musig3)^2)))
  sigx1_4 <- exp(-(0.5 * (epsilon_isq4/ewv4 - (-musig4)^2)))
  sigx2_1 <- (sigx1_1 * ewz1 * pmusig1 * sigmastar1/wzT1)
  sigx2_2 <- (sigx1_2 * ewz2 * pmusig2 * sigmastar2/wzT2)
  sigx2_3 <- (sigx1_3 * ewz3 * pmusig3 * sigmastar3/wzT3)
  sigx2_4 <- (wzlogit * sigx1_4 * pmusig4 * sigmastar4/wzT4)
  sigx3_1 <- (TT * Tpi * ewusqrt1 * ewvTT1 * pmusig1 * sigmastar1/wzT1^2)
  sigx3_2 <- (TT * Tpi * ewusqrt2 * ewvTT2 * pmusig2 * sigmastar2/wzT2^2)
  sigx3_3 <- (TT * Tpi * ewusqrt3 * ewvTT3 * pmusig3 * sigmastar3/wzT3^2)
  sigx3_4 <- (TT * Tpi * ewusqrt4 * ewvTT4 * pmusig4 * sigmastar4/wzT4^2)
  sigx4 <- ((2 * sigx2_1 + 2 * sigx2_2 + 2 * sigx2_3)/wzdeno + 2 * sigx2_4)
  sigx5_1 <- (sigx4 * wzdeno * Tpi * ewusqrt1 * ewvTT1)
  sigx5_2 <- (sigx4 * wzdeno * Tpi * ewusqrt2 * ewvTT2)
  sigx5_3 <- (sigx4 * wzdeno * Tpi * ewusqrt3 * ewvTT3)
  sigx5_4 <- (sigx4 * Tpi * ewusqrt4 * ewvTT4)
  sigx6_1 <- (S * ewu1 * pmusig1 * giepsi1/ssqx2_1 - dmusig1)
  sigx6_2 <- (S * ewu2 * pmusig2 * giepsi2/ssqx2_2 - dmusig2)
  sigx6_3 <- (S * ewu3 * pmusig3 * giepsi3/ssqx2_3 - dmusig3)
  sigx6_4 <- (S * ewu4 * pmusig4 * giepsi4/ssqx2_4 - dmusig4)
  sigx7_1 <- (0.5 * (wvsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (wvsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx7_3 <- (0.5 * (wvsq3 * ewu3/sigmastar3) + sigmastar3)
  sigx7_4 <- (0.5 * (wvsq4 * ewu4/sigmastar4) + sigmastar4)
  Xgsq1 <- Xepsi_i1 - sweep(Xgi1, MARGIN = 1, STATS = 2 * (S^2 * ewu1 * giepsi1/sigmasq1),
    FUN = "*")
  Xgsq2 <- Xepsi_i2 - sweep(Xgi2, MARGIN = 1, STATS = 2 * (S^2 * ewu2 * giepsi2/sigmasq2),
    FUN = "*")
  Xgsq3 <- Xepsi_i3 - sweep(Xgi3, MARGIN = 1, STATS = 2 * (S^2 * ewu3 * giepsi3/sigmasq3),
    FUN = "*")
  Xgsq4 <- Xepsi_i4 - sweep(Xgi4, MARGIN = 1, STATS = 2 * (S^2 * ewu4 * giepsi4/sigmasq4),
    FUN = "*")
  Xsig1 <- sweep(Xgsq1, MARGIN = 1, STATS = 0.5 * (pmusig1/ewv1), FUN = "*") +
    sweep(Xgi1, MARGIN = 1, STATS = S * dmusig1 * ewu1/ssqx2_1, FUN = "*")
  Xsig2 <- sweep(Xgsq2, MARGIN = 1, STATS = 0.5 * (pmusig2/ewv2), FUN = "*") +
    sweep(Xgi2, MARGIN = 1, STATS = S * dmusig2 * ewu2/ssqx2_2, FUN = "*")
  Xsig3 <- sweep(Xgsq3, MARGIN = 1, STATS = 0.5 * (pmusig3/ewv3), FUN = "*") +
    sweep(Xgi3, MARGIN = 1, STATS = S * dmusig3 * ewu3/ssqx2_3, FUN = "*")
  Xsig4 <- sweep(Xgsq4, MARGIN = 1, STATS = 0.5 * (pmusig4/ewv4), FUN = "*") +
    sweep(Xgi4, MARGIN = 1, STATS = S * dmusig4 * ewu4/ssqx2_4, FUN = "*")
  Zsig1 <- sweep(Zigiepsi1, MARGIN = 1, STATS = S * sigx6_1 * sigmastar1/ssqx2_1,
    FUN = "*") - sweep(Zigisq1, MARGIN = 1, STATS = S * sigx6_1 * sigmastar1 *
    ewu1 * (sigmastar1 - 0.5 * (ewu1 * ewv1/ssqx2_1)) * giepsi1/ssqx2_1^2, FUN = "*") -
    sweep(Zigisq1, MARGIN = 1, STATS = 0.5 * (ewu1 * ewv1 * pmusig1/(sigmasq1^2 *
      sigmastar1)), FUN = "*")
  Zsig2 <- sweep(Zigiepsi2, MARGIN = 1, STATS = S * sigx6_2 * sigmastar2/ssqx2_2,
    FUN = "*") - sweep(Zigisq2, MARGIN = 1, STATS = S * sigx6_2 * sigmastar2 *
    ewu2 * (sigmastar2 - 0.5 * (ewu2 * ewv2/ssqx2_2)) * giepsi2/ssqx2_2^2, FUN = "*") -
    sweep(Zigisq2, MARGIN = 1, STATS = 0.5 * (ewu2 * ewv2 * pmusig2/(sigmasq2^2 *
      sigmastar2)), FUN = "*")
  Zsig3 <- sweep(Zigiepsi3, MARGIN = 1, STATS = S * sigx6_3 * sigmastar3/ssqx2_3,
    FUN = "*") - sweep(Zigisq3, MARGIN = 1, STATS = S * sigx6_3 * sigmastar3 *
    ewu3 * (sigmastar3 - 0.5 * (ewu3 * ewv3/ssqx2_3)) * giepsi3/ssqx2_3^2, FUN = "*") -
    sweep(Zigisq3, MARGIN = 1, STATS = 0.5 * (ewu3 * ewv3 * pmusig3/(sigmasq3^2 *
      sigmastar3)), FUN = "*")
  Zsig4 <- sweep(Zigiepsi4, MARGIN = 1, STATS = S * sigx6_4 * sigmastar4/ssqx2_4,
    FUN = "*") - sweep(Zigisq4, MARGIN = 1, STATS = S * sigx6_4 * sigmastar4 *
    ewu4 * (sigmastar4 - 0.5 * (ewu4 * ewv4/ssqx2_4)) * giepsi4/ssqx2_4^2, FUN = "*") -
    sweep(Zigisq4, MARGIN = 1, STATS = 0.5 * (ewu4 * ewv4 * pmusig4/(sigmasq4^2 *
      sigmastar4)), FUN = "*")
  gradll <- cbind(sweep(Xsig1, MARGIN = 1, STATS = -(2 * (sigx1_1 * ewz1 * sigmastar1/sigx5_1)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (((0.5 * (wusq1 * ewv1 *
    pmusig1/ssqx2_1) + S * (1/ssqx2_1 - (0.5 * (wusq1 * ewv1/sigmastar1) + sigmastar1 *
    gisq1) * ewu1/ssqx2_1^2) * sigx6_1 * sigmastar1 * giepsi1) * ewu1/wzT1 -
    0.5 * (Tpi * ewusqrt1 * ewvTT1 * pmusig1 * sigmastar1/wzT1^2)) * sigx1_1 *
    ewz1/(sigx4 * wzdeno)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    ((((S * sigx7_1 * dmusig1 * ewu1 * ewv1 * giepsi1/ssqx2_1^2 - 0.5 * ((2 *
      (S^2 * sigx7_1 * ewu1^2 * ewv1 * giepsi1^2/(ssqx2_1^2 * sigmasq1 * sigmastar1)) -
      epsilon_isq1/ewv1) * pmusig1)) * sigmastar1 + 0.5 * (wvsq1 * ewu1 * ewv1 *
      pmusig1/ssqx2_1))/wzT1 - 0.5 * sigx3_1) * sigx1_1 * ewz1/(sigx4 * wzdeno)),
    FUN = "*"), sweep(Zsig1, MARGIN = 1, STATS = 2 * (sigx1_1 * ewu1 * ewz1/sigx5_1),
    FUN = "*"), sweep(Xsig2, MARGIN = 1, STATS = -(2 * (sigx1_2 * ewz2 * sigmastar2/sigx5_2)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (((0.5 * (wusq2 * ewv2 *
    pmusig2/ssqx2_2) + S * (1/ssqx2_2 - (0.5 * (wusq2 * ewv2/sigmastar2) + sigmastar2 *
    gisq2) * ewu2/ssqx2_2^2) * sigx6_2 * sigmastar2 * giepsi2) * ewu2/wzT2 -
    0.5 * (Tpi * ewusqrt2 * ewvTT2 * pmusig2 * sigmastar2/wzT2^2)) * sigx1_2 *
    ewz2/(sigx4 * wzdeno)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    ((((S * sigx7_2 * dmusig2 * ewu2 * ewv2 * giepsi2/ssqx2_2^2 - 0.5 * ((2 *
      (S^2 * sigx7_2 * ewu2^2 * ewv2 * giepsi2^2/(ssqx2_2^2 * sigmasq2 * sigmastar2)) -
      epsilon_isq2/ewv2) * pmusig2)) * sigmastar2 + 0.5 * (wvsq2 * ewu2 * ewv2 *
      pmusig2/ssqx2_2))/wzT2 - 0.5 * sigx3_2) * sigx1_2 * ewz2/(sigx4 * wzdeno)),
    FUN = "*"), sweep(Zsig2, MARGIN = 1, STATS = 2 * (sigx1_2 * ewu2 * ewz2/sigx5_2),
    FUN = "*"), sweep(Xsig3, MARGIN = 1, STATS = -(2 * (sigx1_3 * ewz3 * sigmastar3/sigx5_3)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (((0.5 * (wusq3 * ewv3 *
    pmusig3/ssqx2_3) + S * (1/ssqx2_3 - (0.5 * (wusq3 * ewv3/sigmastar3) + sigmastar3 *
    gisq3) * ewu3/ssqx2_3^2) * sigx6_3 * sigmastar3 * giepsi3) * ewu3/wzT3 -
    0.5 * (Tpi * ewusqrt3 * ewvTT3 * pmusig3 * sigmastar3/wzT3^2)) * sigx1_3 *
    ewz3/(sigx4 * wzdeno)), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    ((((S * sigx7_3 * dmusig3 * ewu3 * ewv3 * giepsi3/ssqx2_3^2 - 0.5 * ((2 *
      (S^2 * sigx7_3 * ewu3^2 * ewv3 * giepsi3^2/(ssqx2_3^2 * sigmasq3 * sigmastar3)) -
      epsilon_isq3/ewv3) * pmusig3)) * sigmastar3 + 0.5 * (wvsq3 * ewu3 * ewv3 *
      pmusig3/ssqx2_3))/wzT3 - 0.5 * sigx3_3) * sigx1_3 * ewz3/(sigx4 * wzdeno)),
    FUN = "*"), sweep(Zsig3, MARGIN = 1, STATS = 2 * (sigx1_3 * ewu3 * ewz3/sigx5_3),
    FUN = "*"), sweep(Xsig4, MARGIN = 1, STATS = -(2 * (wzlogit * sigx1_4 * sigmastar4/sigx5_4)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (((0.5 * (wusq4 * ewv4 *
    pmusig4/ssqx2_4) + S * (1/ssqx2_4 - (0.5 * (wusq4 * ewv4/sigmastar4) + sigmastar4 *
    gisq4) * ewu4/ssqx2_4^2) * sigx6_4 * sigmastar4 * giepsi4) * ewu4/wzT4 -
    0.5 * (Tpi * ewusqrt4 * ewvTT4 * pmusig4 * sigmastar4/wzT4^2)) * wzlogit *
    sigx1_4/sigx4), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * ((((S *
    sigx7_4 * dmusig4 * ewu4 * ewv4 * giepsi4/ssqx2_4^2 - 0.5 * ((2 * (S^2 *
    sigx7_4 * ewu4^2 * ewv4 * giepsi4^2/(ssqx2_4^2 * sigmasq4 * sigmastar4)) -
    epsilon_isq4/ewv4) * pmusig4)) * sigmastar4 + 0.5 * (wvsq4 * ewu4 * ewv4 *
    pmusig4/ssqx2_4))/wzT4 - 0.5 * sigx3_4) * wzlogit * sigx1_4/sigx4), FUN = "*"),
    sweep(Zsig4, MARGIN = 1, STATS = 2 * (wzlogit * sigx1_4 * ewu4/sigx5_4),
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 * (sigx1_1 * pmusig1 *
      sigmastar1/wzT1) - sigx4) * ewz1/(sigx4 * wzdeno), FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = (2 * (sigx1_2 * pmusig2 * sigmastar2/wzT2) - sigx4) *
        ewz2/(sigx4 * wzdeno), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 *
      (sigx1_3 * pmusig3 * sigmastar3/wzT3) - sigx4) * ewz3/(sigx4 * wzdeno),
      FUN = "*"))
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param Zvar separating variables
#' @param nZHvar number of separating variables
#' @param wHvar_c vector of weights (weighted likelihood) pooled data
#' @param wHvar_p vector of weights (weighted likelihood) cross-section
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @param modelType specification of inefficiency model G(t)u_i
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
LCM4ChnormAlgOpt_gzit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, ngZGvar,
  Zvar, nZHvar, Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, modelType, method, printInfo,
  itermax, stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm4C_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c,
      vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, S = S, Zvar = Zvar, nZHvar = nZHvar, pindex = pindex, TT = TT,
      wHvar_c = wHvar_c, wHvar_p = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar,
      modelType = modelType, initIter = initIter, initAlg = initAlg, whichStart = whichStart,
      tol = tol, printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pLCMhalfnormlike4C_gzit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("LCM Panel BC92-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
      nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike4C_gzit,
    grad = pgradLCMhalfnormlike4C_gzit, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
      nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      Zvar = Zvar, nZHvar = nZHvar, S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
    prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike4C_gzit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(parm)), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, Zvar = Zvar, nZHvar = nZHvar, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike4C_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike4C_gzit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike4C_gzit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike4C_gzit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike4C_gzit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike4C_gzit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitHalf = initHalf))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Posterior probabilities and efficiencies ----------
#' post. prob. and efficiencies for lcmpanel 4 classes halfnormal-normal distribution
#' @param object object of class lcmpanel
#' @param level level for confidence interval
#' @noRd
pLCM4Chalfnormeff_gzit <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta1 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + ncol(object$gHvar))]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + ncol(object$gHvar) +
    1):(2 * object$nXvar + object$nuZUvar + object$nvZVvar + ncol(object$gHvar))]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar) + 1):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar))]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar) + 1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    ncol(object$gHvar))]
  eta2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    ncol(object$gHvar) + 1):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar))]
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar) + 1):(3 * object$nXvar + 2 * object$nuZUvar + 2 *
    object$nvZVvar + 2 * ncol(object$gHvar))]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar) + 1):(3 * object$nXvar + 3 * object$nuZUvar + 2 *
    object$nvZVvar + 2 * ncol(object$gHvar))]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar) + 1):(3 * object$nXvar + 3 * object$nuZUvar + 3 *
    object$nvZVvar + 2 * ncol(object$gHvar))]
  eta3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    2 * ncol(object$gHvar) + 1):(3 * object$nXvar + 3 * object$nuZUvar + 3 *
    object$nvZVvar + 3 * ncol(object$gHvar))]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    3 * ncol(object$gHvar) + 1):(4 * object$nXvar + 3 * object$nuZUvar + 3 *
    object$nvZVvar + 3 * ncol(object$gHvar))]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 3 * object$nvZVvar +
    3 * ncol(object$gHvar) + 1):(4 * object$nXvar + 4 * object$nuZUvar + 3 *
    object$nvZVvar + 3 * ncol(object$gHvar))]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 3 * object$nvZVvar +
    3 * ncol(object$gHvar) + 1):(4 * object$nXvar + 4 * object$nuZUvar + 4 *
    object$nvZVvar + 3 * ncol(object$gHvar))]
  eta4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    3 * ncol(object$gHvar) + 1):(4 * object$nXvar + 4 * object$nuZUvar + 4 *
    object$nvZVvar + 4 * ncol(object$gHvar))]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    4 * ncol(object$gHvar) + 1):(4 * object$nXvar + 4 * object$nuZUvar + 4 *
    object$nvZVvar + 4 * ncol(object$gHvar) + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    4 * ncol(object$gHvar) + object$nZHvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 4 * ncol(object$gHvar) + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 4 * object$nvZVvar +
    4 * ncol(object$gHvar) + 2 * object$nZHvar + 1):(4 * object$nXvar + 4 * object$nuZUvar +
    4 * object$nvZVvar + 4 * ncol(object$gHvar) + 3 * object$nZHvar)]
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
  if (object$modelType == "bc92c") {
    Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  } else {
    Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  }
  Zvar <- apply(Zvar, 2, function(x) tapply(x, pindex[, 1], mean))
  TT <- as.numeric(table(pindex[, 1]))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar_p)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar_p)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  epsilon_it1 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- exp(as.numeric(crossprod(matrix(eta1), t(object$gHvar))))
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
  git2 <- exp(as.numeric(crossprod(matrix(eta2), t(object$gHvar))))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * object$S * giepsi2/(exp(Wv2) + gisq2 * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + gisq2 * exp(Wu2)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar_p)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar_p)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon_it3 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon_isq3 <- as.numeric(tapply(epsilon_it3^2, pindex[, 1], sum))
  git3 <- exp(as.numeric(crossprod(matrix(eta3), t(object$gHvar))))
  git_epsit3 <- epsilon_it3 * git3
  giepsi3 <- as.numeric(tapply(git_epsit3, pindex[, 1], sum))
  gisq3 <- as.numeric(tapply(git3^2, pindex[, 1], sum))
  mustar3 <- -exp(Wu3) * object$S * giepsi3/(exp(Wv2) + gisq3 * exp(Wu3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wv3) + gisq3 * exp(Wu3)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar_p)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar_p)))
  epsilon_it4 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  epsilon_isq4 <- as.numeric(tapply(epsilon_it4^2, pindex[, 1], sum))
  git4 <- exp(as.numeric(crossprod(matrix(eta4), t(object$gHvar))))
  git_epsit4 <- epsilon_it4 * git4
  giepsi4 <- as.numeric(tapply(git_epsit4, pindex[, 1], sum))
  gisq4 <- as.numeric(tapply(git4^2, pindex[, 1], sum))
  mustar4 <- -exp(Wu4) * object$S * giepsi4/(exp(Wv2) + gisq4 * exp(Wu4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wv4) + gisq4 * exp(Wu4)))
  Pi1 <- 2 * sigmastar1 * exp(-1/2 * (epsilon_isq1/exp(Wv1) - (mustar1/sigmastar1)^2)) *
    pnorm(mustar1/sigmastar1)/((2 * pi)^(TT/2) * exp(Wv1/2 * TT) * exp(Wu1/2))
  Pi2 <- 2 * sigmastar2 * exp(-1/2 * (epsilon_isq2/exp(Wv2) - (mustar2/sigmastar2)^2)) *
    pnorm(mustar2/sigmastar2)/((2 * pi)^(TT/2) * exp(Wv2/2 * TT) * exp(Wu2/2))
  Pi3 <- 2 * sigmastar3 * exp(-1/2 * (epsilon_isq3/exp(Wv3) - (mustar3/sigmastar3)^2)) *
    pnorm(mustar3/sigmastar3)/((2 * pi)^(TT/2) * exp(Wv3/2 * TT) * exp(Wu3/2))
  Pi4 <- 2 * sigmastar4 * exp(-1/2 * (epsilon_isq4/exp(Wv4) - (mustar4/sigmastar4)^2)) *
    pnorm(mustar4/sigmastar4)/((2 * pi)^(TT/2) * exp(Wv4/2 * TT) * exp(Wu4/2))
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * Probc3 + Pi4 *
    Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c == 2, Pcond_c2, ifelse(Group_c ==
    3, Pcond_c3, Pcond_c4)))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- mustar4 + sigmastar4 * dnorm(mustar4/sigmastar4)/pnorm(mustar4/sigmastar4)
  res <- data.frame(levels(pindex[, 1]), Group_c = Group_c, PosteriorProb_c = P_cond_c,
    PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
    PriorProb_c2 = Probc2, u_c2 = u_c2, PosteriorProb_c3 = Pcond_c3, PriorProb_c3 = Probc3,
    u_c3 = u_c3, PosteriorProb_c4 = Pcond_c4, PriorProb_c4 = Probc4, u_c4 = u_c4)
  if (object$logDepVar == TRUE) {
    res <- data.frame(res, mustar1 = mustar1, mustar2 = mustar2, mustar3 = mustar3,
      mustar4 = mustar4, sigmastar1 = sigmastar1, sigmastar2 = sigmastar2,
      sigmastar3 = sigmastar3, sigmastar4 = sigmastar4)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u_c1 <- res$u_c1 * git1
  res$u_c2 <- res$u_c2 * git2
  res$u_c3 <- res$u_c3 * git3
  res$u_c4 <- res$u_c4 * git4
  res$u_c <- ifelse(res$Group_c == 1, res$u_c1, ifelse(res$Group_c == 2, res$u_c2,
    ifelse(res$Group_c == 3, res$u_c3, res$u_c4)))
  res$ineff_c1 <- ifelse(res$Group_c == 1, res$u_c1, NA)
  res$ineff_c2 <- ifelse(res$Group_c == 2, res$u_c2, NA)
  res$ineff_c3 <- ifelse(res$Group_c == 2, res$u_c3, NA)
  res$ineff_c4 <- ifelse(res$Group_c == 2, res$u_c4, NA)
  if (object$logDepVar == TRUE) {
    res$teJLMS_c <- exp(-res$u_c)
    res$teBC_c1 <- exp(1/2 * res$sigmastar1^2 * git1^2 - res$mustar1 * git1) *
      pnorm(res$mustar1/res$sigmastar1 - res$sigmastar1 * git1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_c2 <- exp(1/2 * res$sigmastar2^2 * git2^2 - res$mustar2 * git2) *
      pnorm(res$mustar2/res$sigmastar2 - res$sigmastar2 * git2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_c3 <- exp(1/2 * res$sigmastar3^2 * git3^2 - res$mustar3 * git3) *
      pnorm(res$mustar3/res$sigmastar3 - res$sigmastar3 * git3)/pnorm(res$mustar3/res$sigmastar3)
    res$teBC_c4 <- exp(1/2 * res$sigmastar4^2 * git4^2 - res$mustar4 * git4) *
      pnorm(res$mustar4/res$sigmastar4 - res$sigmastar4 * git4)/pnorm(res$mustar4/res$sigmastar4)
    res$teBC_c <- ifelse(res$Group_c == 1, res$teBC_c1, ifelse(res$Group_c ==
      2, res$teBC_c2, res$teBC_c3))
    res$effBC_c1 <- ifelse(res$Group_c == 1, res$teBC_c1, NA)
    res$effBC_c2 <- ifelse(res$Group_c == 2, res$teBC_c2, NA)
    res$effBC_c3 <- ifelse(res$Group_c == 2, res$teBC_c3, NA)
    res$effBC_c4 <- ifelse(res$Group_c == 2, res$teBC_c4, NA)
    res$teBC_reciprocal_c1 <- exp(res$mustar1 + 1/2 * res$sigmastar1^2) * pnorm(res$mustar1/res$sigmastar1 +
      res$sigmastar1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_reciprocal_c2 <- exp(res$mustar2 + 1/2 * res$sigmastar2^2) * pnorm(res$mustar2/res$sigmastar2 +
      res$sigmastar2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_reciprocal_c3 <- exp(res$mustar3 + 1/2 * res$sigmastar3^2) * pnorm(res$mustar3/res$sigmastar3 +
      res$sigmastar3)/pnorm(res$mustar3/res$sigmastar3)
    res$teBC_reciprocal_c4 <- exp(res$mustar4 + 1/2 * res$sigmastar4^2) * pnorm(res$mustar4/res$sigmastar4 +
      res$sigmastar4)/pnorm(res$mustar4/res$sigmastar4)
    res$teBC_reciprocal_c <- res$teBC_reciprocal_c <- ifelse(res$Group_c == 1,
      res$teBC_reciprocal_c1, ifelse(res$Group_c == 2, res$teBC_reciprocal_c2,
        ifelse(res$Group_c == 3, res$teBC_reciprocal_c3, res$teBC_reciprocal_c4)))
    res$ReffBC_c1 <- ifelse(res$Group_c == 1, res$teBC_reciprocal_c1, NA)
    res$ReffBC_c2 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c2, NA)
    res$ReffBC_c3 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c3, NA)
    res$ReffBC_c4 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c4, NA)
    res$mustar1 <- NULL
    res$sigmastar1 <- NULL
    res$mustar2 <- NULL
    res$sigmastar2 <- NULL
    res$mustar3 <- NULL
    res$sigmastar3 <- NULL
    res$mustar4 <- NULL
    res$sigmastar4 <- NULL
  }
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
