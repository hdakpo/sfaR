################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 2L                                                        #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Modified Lee and Schmidt 1993                        #
#                          - g(zit) = exp(-eta_t * (t - T)): g(zit) = 1 for T  #
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
pLCMhalfnormlike2C_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, Zvar, nZHvar, ngZGvar, gHvar, wHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar -
    1)]
  eta1[ngZGvar] <- 1
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + nuZUvar + nvZVvar +
    ngZGvar - 1)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + 2 * nuZUvar +
    nvZVvar + ngZGvar - 1)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + ngZGvar - 1)]
  eta2 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + ngZGvar):(2 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 2 * ngZGvar - 2)]
  eta2[ngZGvar] <- 1
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar - 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar + 2 * ngZGvar - 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
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
psLCMfhalfnorm2C_mols93 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
  vHvar_c, uHvar_p, vHvar_p, Yvar, Xvar, ngZGvar, gHvar, S, wHvar_c, wHvar_p, Zvar,
  nZHvar, pindex, TT, whichStart, initIter, initAlg, printInfo, tol) {
  initHalf <- psthalfnorm_mols93(olsObj = olsObj, epsiRes = epsiRes, nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, ngZGvar = ngZGvar, gHvar = gHvar, uHvar = uHvar_c[,
      1, drop = FALSE], vHvar = vHvar_c[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar_c, initIter = initIter, whichStart = whichStart, initAlg = initAlg,
    tol = tol, printInfo = printInfo)
  if (whichStart == 1L) {
    Esti <- initHalf$StartVal
    initHalfPanel <- NULL
  } else {
    cat("Initialization: SFA Panel MOLS93-type + halfnormal-normal distribution...\n")
    initHalfPanel <- maxLik::maxLik(logLik = phalfnormlike_mols93, start = initHalf$StartVal,
      grad = pgradhalfnormlike_mols93, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_p[, 1, drop = FALSE], vHvar = vHvar_p[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar)
    Esti <- initHalfPanel$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[(nXvar +
    3):(nXvar + ngZGvar + 1)], 0.98 * Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar -
    1), Esti[(nXvar + 3):(nXvar + ngZGvar + 1)], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), paste0("eta_", colnames(gHvar[, -ngZGvar])),
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)), paste0("Zv_", colnames(vHvar_p)),
    paste0("eta_", colnames(gHvar[, -ngZGvar])), paste0("Cl1_", colnames(Zvar)))
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
pgradLCMhalfnormlike2C_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, Zvar, nZHvar, ngZGvar, gHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar -
    1)]
  eta1[ngZGvar] <- 1
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + nuZUvar + nvZVvar +
    ngZGvar - 1)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + 2 * nuZUvar +
    nvZVvar + ngZGvar - 1)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + ngZGvar):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + ngZGvar - 1)]
  eta2 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + ngZGvar):(2 * nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 2 * ngZGvar - 2)]
  eta2[ngZGvar] <- 1
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 2 * ngZGvar - 1):(2 *
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar + 2 * ngZGvar - 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
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
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvTT1 <- exp(Wv1 * TT/2)
  ewvTT2 <- exp(Wv2 * TT/2)
  ewusqrt1 <- exp(Wu1/2)
  ewusqrt2 <- exp(Wu2/2)
  ewz <- exp(Wz)
  sigmasq1 <- (ewu1 * gisq1 + ewv1)
  sigmasq2 <- (ewu2 * gisq2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigmasq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigmasq2)
  ssqx2_1 <- sigmasq1 * sigmastar1
  ssqx2_2 <- sigmasq2 * sigmastar2
  musig1 <- (S * ewu1 * giepsi1/(ssqx2_1))
  musig2 <- (S * ewu2 * giepsi2/(ssqx2_2))
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1)
  dmusig2 <- dnorm(-musig2)
  wzdeno <- (1 + ewz)
  Tpi <- (2 * pi)^(TT/2)
  wvsq1 <- (1 - ewv1/sigmasq1)
  wvsq2 <- (1 - ewv2/sigmasq2)
  wzT1 <- (wzdeno * Tpi * ewusqrt1 * ewvTT1)
  wzT2 <- (Tpi * ewusqrt2 * ewvTT2)
  wusq1 <- (1 - ewu1 * gisq1/sigmasq1)
  wusq2 <- (1 - ewu2 * gisq2/sigmasq2)
  wzlogit <- (1 - ewz/wzdeno)
  sigx1_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-musig1)^2)))
  sigx1_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-musig2)^2)))
  sigx2_1 <- (sigx1_1 * ewz * pmusig1 * sigmastar1/wzT1)
  sigx2_2 <- (wzlogit * sigx1_2 * pmusig2 * sigmastar2/wzT2)
  sigx3_1 <- (TT * wzdeno * Tpi * ewusqrt1 * ewvTT1 * pmusig1 * sigmastar1/wzT1^2)
  sigx3_2 <- (Tpi * ewusqrt2 * ewvTT2 * pmusig2 * sigmastar2/wzT2^2)
  wvsqwu1 <- wvsq1 * ewu1
  wusqwv1 <- wusq1 * ewv1
  wvsqwu2 <- wvsq2 * ewu2
  wusqwv2 <- wusq2 * ewv2
  wuv1 <- ewu1 * ewv1
  wuv2 <- ewu2 * ewv2
  wup1 <- wuv1 * pmusig1
  wup2 <- wuv2 * pmusig2
  sigx4_1 <- (wvsqwu1/sigmastar1)
  sigx4_2 <- (wvsqwu2/sigmastar2)
  sigx5_1 <- (wusqwv1/sigmastar1)
  sigx5_2 <- (wusqwv2/sigmastar2)
  sigZ1_1 <- sweep(Zigisq1, MARGIN = 1, STATS = (wup1/(sigmasq1^2 * sigmastar1)),
    FUN = "*")
  sigZ1_2 <- sweep(Zigisq2, MARGIN = 1, STATS = (wup2/(sigmasq2^2 * sigmastar2)),
    FUN = "*")
  Xgsq1 <- Xepsi_i1 - sweep(Xgi1, MARGIN = 1, STATS = 2 * ewu1 * giepsi1/sigmasq1,
    FUN = "*")
  Xgsq2 <- Xepsi_i2 - sweep(Xgi2, MARGIN = 1, STATS = 2 * ewu2 * giepsi2/sigmasq2,
    FUN = "*")
  sigx7_1 <- sweep(Xgsq1, MARGIN = 1, STATS = pmusig1/ewv1, FUN = "*")
  sigx7_2 <- sweep(Xgsq2, MARGIN = 1, STATS = pmusig2/ewv2, FUN = "*")
  Xsig7_1 <- 0.5 * sigx7_1 + sweep(Xgi1, MARGIN = 1, STATS = S * dmusig1 * ewu1/(ssqx2_1),
    FUN = "*")
  Xsig7_2 <- 0.5 * sigx7_2 + sweep(Xgi2, MARGIN = 1, STATS = S * dmusig2 * ewu2/(ssqx2_2),
    FUN = "*")
  ssig4_1 <- (0.5 * sigx4_1 + sigmastar1)
  ssig4_2 <- (0.5 * sigx4_2 + sigmastar2)
  sigx8_1 <- (S^2 * ssig4_1 * ewu1^2 * ewv1 * giepsi1^2/((ssqx2_1)^2 * ssqx2_1))
  sigx8_2 <- (S^2 * ssig4_2 * ewu2^2 * ewv2 * giepsi2^2/((ssqx2_2)^2 * ssqx2_2))
  sigx9_1 <- ((2 * sigx8_1 - epsilon_isq1/ewv1) * pmusig1)
  sigx9_2 <- ((2 * sigx8_2 - epsilon_isq2/ewv2) * pmusig2)
  sigx10_1 <- (1/wzT1 - Tpi * ewusqrt1 * ewvTT1 * ewz/wzT1^2)
  sigx10_2 <- (wzdeno * Tpi * ewusqrt2 * ewvTT2)
  wzsig1_2 <- wzlogit * sigx1_2
  sigx11_1 <- (sigx10_1 * sigx1_1 * pmusig1 * sigmastar1)
  sigx11_2 <- (wzsig1_2 * pmusig2 * sigmastar2/sigx10_2)
  sigx12_1 <- (wzdeno * Tpi * ewusqrt1 * ewvTT1 * pmusig1 * sigmastar1/wzT1^2)
  sigx12_2 <- (TT * Tpi * ewusqrt2 * ewvTT2 * pmusig2 * sigmastar2/wzT2^2)
  siggsq1 <- (0.5 * sigx5_1 + sigmastar1 * gisq1)
  siggsq2 <- (0.5 * sigx5_2 + sigmastar2 * gisq2)
  sigx13_1 <- (1/(ssqx2_1) - siggsq1 * ewu1/(ssqx2_1)^2)
  sigx13_2 <- (1/(ssqx2_2) - siggsq2 * ewu2/(ssqx2_2)^2)
  sigx14_1 <- (wvsqwu1 * ewv1 * pmusig1/(ssqx2_1))
  sigx14_2 <- (wvsqwu2 * ewv2 * pmusig2/(ssqx2_2))
  sig2sq <- (2 * sigx2_2 + 2 * sigx2_1)
  sigx15 <- (wzdeno * sig2sq * Tpi * ewusqrt1 * ewvTT1)
  sigx16_1 <- (wusqwv1 * pmusig1/(ssqx2_1))
  sigx16_2 <- (wusqwv2 * pmusig2/(ssqx2_2))
  sigx17_1 <- (S * ewu1 * pmusig1 * giepsi1/(ssqx2_1) - dmusig1)
  sigx17_2 <- (S * ewu2 * pmusig2 * giepsi2/(ssqx2_2) - dmusig2)
  sigmix_1 <- (0.5 * sigx16_1 + S * sigx13_1 * sigx17_1 * sigmastar1 * giepsi1)
  sigmix_2 <- (0.5 * sigx16_2 + S * sigx13_2 * sigx17_2 * sigmastar2 * giepsi2)
  sigx18_1 <- (sigmix_1 * ewu1/wzT1 - 0.5 * sigx12_1)
  sigx18_2 <- (sigmix_2 * ewu2/wzT2 - 0.5 * sigx3_2)
  dmussqx2_1 <- dmusig1 * wuv1 * giepsi1/(ssqx2_1)^2
  dmussqx2_2 <- dmusig2 * wuv2 * giepsi2/(ssqx2_2)^2
  sigx19_1 <- (S * ssig4_1 * dmussqx2_1 - 0.5 * sigx9_1)
  sigx19_2 <- (S * ssig4_2 * dmussqx2_2 - 0.5 * sigx9_2)
  sigx20_1 <- sigx18_1 * sigx1_1 * ewz/sig2sq
  sigx20_2 <- sigx18_2 * wzsig1_2/sig2sq
  sigx21_1 <- ((sigx19_1 * sigmastar1 + 0.5 * sigx14_1)/wzT1 - 0.5 * sigx3_1)
  sigx21_2 <- ((sigx19_2 * sigmastar2 + 0.5 * sigx14_2)/wzT2 - 0.5 * sigx12_2)
  sigx22_1 <- sigx21_1 * sigx1_1 * ewz/sig2sq
  sigx22_2 <- sigx21_2 * wzsig1_2/sig2sq
  sigx23_1 <- (sigmastar1 - 0.5 * (wuv1/(ssqx2_1)))
  sigx23_2 <- (sigmastar2 - 0.5 * (wuv2/(ssqx2_2)))
  sigZ2_1 <- sweep(Zigiepsi1, MARGIN = 1, STATS = 1/(ssqx2_1), FUN = "*") - sweep(Zigisq1,
    MARGIN = 1, STATS = ewu1 * sigx23_1 * giepsi1/(ssqx2_1)^2, FUN = "*")
  sigZ2_2 <- sweep(Zigiepsi2, MARGIN = 1, STATS = 1/(ssqx2_2), FUN = "*") - sweep(Zigisq2,
    MARGIN = 1, STATS = ewu2 * sigx23_2 * giepsi2/(ssqx2_2)^2, FUN = "*")
  sigZ3_1 <- sweep(sigZ2_1, MARGIN = 1, STATS = S * sigx17_1 * sigmastar1, FUN = "*") -
    0.5 * sigZ1_1
  sigZ3_2 <- sweep(sigZ2_2, MARGIN = 1, STATS = S * sigx17_2 * sigmastar2, FUN = "*") -
    0.5 * sigZ1_2
  sigx26 <- (sig2sq * Tpi * ewusqrt2 * ewvTT2)
  sigx27 <- (2 * sigx11_1 - 2 * sigx11_2)
  sigx28 <- sigx27 * ewz/sig2sq
  gradll <- cbind(sweep(Xsig7_1, MARGIN = 1, STATS = -(2 * (sigx1_1 * ewz * sigmastar1/sigx15)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx20_1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx22_1), FUN = "*"), sweep(sigZ3_1,
      MARGIN = 1, STATS = 2 * (sigx1_1 * ewu1 * ewz/sigx15), FUN = "*"), sweep(Xsig7_2,
      MARGIN = 1, STATS = -(2 * (wzsig1_2 * sigmastar2/sigx26)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx20_2), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = 2 * (sigx22_2), FUN = "*"), sweep(sigZ3_2, MARGIN = 1,
      STATS = 2 * (wzsig1_2 * ewu2/sigx26), FUN = "*"), sweep(Zvar, MARGIN = 1,
      STATS = sigx28, FUN = "*"))
  return(sweep(gradll[, -c((nXvar + nuZUvar + nvZVvar + ngZGvar), 2 * (nXvar +
    nuZUvar + nvZVvar + ngZGvar))], MARGIN = 1, STATS = wHvar, FUN = "*"))
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
LCM2ChnormAlgOpt_mols93 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, ngZGvar,
  Zvar, nZHvar, Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm2C_mols93(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar_c = uHvar_c,
      vHvar_c = vHvar_c, uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, S = S, Zvar = Zvar, nZHvar = nZHvar, pindex = pindex, TT = TT,
      wHvar_c = wHvar_c, wHvar_p = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar,
      initIter = initIter, initAlg = initAlg, whichStart = whichStart, tol = tol,
      printInfo = printInfo)
    initHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pLCMhalfnormlike2C_mols93(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("LCM Panel MOLS93-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
      nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike2C_mols93,
    grad = pgradLCMhalfnormlike2C_mols93, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
      nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      Zvar = Zvar, nZHvar = nZHvar, S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
    prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike2C_mols93(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(parm)), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, Zvar = Zvar, nZHvar = nZHvar, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike2C_mols93(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_mols93(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_mols93(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)),
        unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike2C_mols93(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike2C_mols93(parm = mlParam, nXvar = nXvar,
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
#' post. prob. and efficiencies for lcmpanel 2 classes halfnormal-normal distribution
#' @param object object of class lcmpanel
#' @param level level for confidence interval
#' @noRd
pLCM2Chalfnormeff_mols93 <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta1 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + ncol(object$gHvar) - 1)]
  eta1[ncol(object$gHvar)] <- 1
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + ncol(object$gHvar)):(2 *
    object$nXvar + object$nuZUvar + object$nvZVvar + ncol(object$gHvar) - 1)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar)):(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar) - 1)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + object$nvZVvar +
    ncol(object$gHvar)):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    ncol(object$gHvar) - 1)]
  eta2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    ncol(object$gHvar)):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar) - 2)]
  eta2[ncol(object$gHvar)] <- 1
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    2 * ncol(object$gHvar) - 1):(2 * object$nXvar + 2 * object$nuZUvar + 2 *
    object$nvZVvar + 2 * ncol(object$gHvar) + object$nZHvar - 2)]
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
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- exp(as.numeric(crossprod(matrix(eta1), t(object$gHvar))))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * object$S * giepsi1/(exp(Wv1) + gisq1 * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + gisq1 * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar_p)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar_p)))
  epsilon_it2 <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- exp(as.numeric(crossprod(matrix(eta2), t(object$gHvar))))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  mustar2 <- -exp(Wu2) * object$S * giepsi2/(exp(Wv2) + gisq2 * exp(Wu2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wv2) + gisq2 * exp(Wu2)))
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
  res <- data.frame(levels(pindex[, 1]), Group_c = Group_c, PosteriorProb_c = P_cond_c,
    PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
    PriorProb_c2 = Probc2, u_c2 = u_c2)
  if (object$logDepVar == TRUE) {
    res <- data.frame(res, mustar1 = mustar1, mustar2 = mustar2, sigmastar1 = sigmastar1,
      sigmastar2 = sigmastar2)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u_c1 <- res$u_c1 * git1
  res$u_c2 <- res$u_c2 * git2
  res$u_c <- ifelse(res$Group_c == 1, res$u_c1, res$u_c2)
  res$ineff_c1 <- ifelse(res$Group_c == 1, res$u_c1, NA)
  res$ineff_c2 <- ifelse(res$Group_c == 2, res$u_c2, NA)
  if (object$logDepVar == TRUE) {
    res$teJLMS_c <- exp(-res$u_c)
    res$teBC_c1 <- exp(1/2 * res$sigmastar1^2 * git1^2 - res$mustar1 * git1) *
      pnorm(res$mustar1/res$sigmastar1 - res$sigmastar1 * git1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_c2 <- exp(1/2 * res$sigmastar2^2 * git2^2 - res$mustar2 * git2) *
      pnorm(res$mustar2/res$sigmastar2 - res$sigmastar2 * git2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_c <- ifelse(res$Group_c == 1, res$teBC_c1, res$teBC_c2)
    res$effBC_c1 <- ifelse(res$Group_c == 1, res$teBC_c1, NA)
    res$effBC_c2 <- ifelse(res$Group_c == 2, res$teBC_c2, NA)
    res$teBC_reciprocal_c1 <- exp(res$mustar1 + 1/2 * res$sigmastar1^2) * pnorm(res$mustar1/res$sigmastar1 +
      res$sigmastar1)/pnorm(res$mustar1/res$sigmastar1)
    res$teBC_reciprocal_c2 <- exp(res$mustar2 + 1/2 * res$sigmastar2^2) * pnorm(res$mustar2/res$sigmastar2 +
      res$sigmastar2)/pnorm(res$mustar2/res$sigmastar2)
    res$teBC_reciprocal_c <- ifelse(res$Group_c == 1, res$teBC_reciprocal_c1,
      res$teBC_reciprocal_c2)
    res$ReffBC_c1 <- ifelse(res$Group_c == 1, res$teBC_reciprocal_c1, NA)
    res$ReffBC_c2 <- ifelse(res$Group_c == 2, res$teBC_reciprocal_c2, NA)
    res$mustar1 <- NULL
    res$sigmastar1 <- NULL
    res$mustar2 <- NULL
    res$sigmastar2 <- NULL
  }
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
