################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Latent Class Stochastic Frontier Model                          #
# Number of Classes: 2L                                                        #
# Type:      - Kumbakhar 1990 (K90)                                            #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = (1 + exp(eta1 * t + eta2 * t^2))^(-1)                  #
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
pLCMhalfnormlike2C_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, Zvar, nZHvar, wHvar) {
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
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 4)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1_1 * seq(1:x) + eta1_2 *
    (seq(1:x))^2)))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  mustar1 <- -exp(Wu1) * S * giepsi1/(exp(Wv1) + gisq1 * exp(Wu1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wv1) + gisq1 * exp(Wu1)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta2_1 * seq(1:x) + eta2_2 *
    (seq(1:x))^2)))
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
#' @param wHvar_c vector of weights (weighted likelihood) for pooled data
#' @param wHvar_p vector of weights (weighted likelihood) for cross section
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
psLCMfhalfnorm2C_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar_c,
  vHvar_c, uHvar_p, vHvar_p, Yvar, Xvar, S, wHvar_c, wHvar_p, Zvar, nZHvar, pindex,
  TT, whichStart, initIter, initAlg, printInfo, tol) {
  initHalf <- psthalfnorm_k90(olsObj = olsObj, epsiRes = epsiRes, nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_c[, 1, drop = FALSE], vHvar = vHvar_c[,
      1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, initIter = initIter,
    whichStart = whichStart, initAlg = initAlg, tol = tol, printInfo = printInfo)
  if (whichStart == 1L) {
    Esti <- initHalf$StartVal
    initHalfPanel <- NULL
  } else {
    cat("Initialization: SFA Panel K90-type + halfnormal-normal distribution...\n")
    initHalfPanel <- maxLik::maxLik(logLik = phalfnormlike_k90, start = initHalf$StartVal,
      grad = pgradhalfnormlike_k90, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, nvZVvar = 1, uHvar = uHvar_p[, 1, drop = FALSE], vHvar = vHvar_p[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
    Esti <- initHalfPanel$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar -
    1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[(nXvar +
    3):(nXvar + 4)], 0.98 * Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) rep(0,
    nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[(nXvar +
    3):(nXvar + 4)], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar_p)),
    paste0("Zv_", colnames(vHvar_p)), "eta1", "eta2", names(Esti)[1:nXvar], paste0("Zu_",
      colnames(uHvar_p)), paste0("Zv_", colnames(vHvar_p)), "eta1", "eta2",
    paste0("Cl1_", colnames(Zvar)))
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
pgradLCMhalfnormlike2C_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 5):(2 * nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 4)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon_it1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon_isq1 <- as.numeric(tapply(epsilon_it1^2, pindex[, 1], sum))
  git1 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1_1 * seq(1:x) + eta1_2 *
    (seq(1:x))^2)))
  git_epsit1 <- epsilon_it1 * git1
  giepsi1 <- as.numeric(tapply(git_epsit1, pindex[, 1], sum))
  gisq1 <- as.numeric(tapply(git1^2, pindex[, 1], sum))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  epsilon_it2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon_isq2 <- as.numeric(tapply(epsilon_it2^2, pindex[, 1], sum))
  git2 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta2_1 * seq(1:x) + eta2_2 *
    (seq(1:x))^2)))
  git_epsit2 <- epsilon_it2 * git2
  giepsi2 <- as.numeric(tapply(git_epsit2, pindex[, 1], sum))
  gisq2 <- as.numeric(tapply(git2^2, pindex[, 1], sum))
  Xepsi_it1 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it1, FUN = "*")
  Xepsi_i1 <- apply(Xepsi_it1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it2 <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it2, FUN = "*")
  Xepsi_i2 <- apply(Xepsi_it2, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit1 <- sweep(-Xvar, MARGIN = 1, STATS = git1, FUN = "*")
  Xgi1 <- apply(Xgit1, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit2 <- sweep(-Xvar, MARGIN = 1, STATS = git2, FUN = "*")
  Xgi2 <- apply(Xgit2, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1_1 <- Zit * exp(Zit * (eta1_1 + eta1_2 * Zit))
  gitd1_epsit_gitsq1 <- -gitd1_1 * epsilon_it1 * git1^2
  gid1_epsi_gisq1 <- as.numeric(tapply(gitd1_epsit_gitsq1, pindex[, 1], sum))
  gitd1_2 <- Zit * exp(Zit * (eta2_1 + eta2_2 * Zit))
  gitd1_epsit_gitsq2 <- -gitd1_2 * epsilon_it2 * git2^2
  gid1_epsi_gisq2 <- as.numeric(tapply(gitd1_epsit_gitsq2, pindex[, 1], sum))
  gitd2_1 <- 2 * Zit * exp(Zit * (eta1_1 + eta1_2 * Zit))
  gitd2_gitcub1 <- -gitd2_1 * git1^3
  gid2_gicub1 <- as.numeric(tapply(gitd2_gitcub1, pindex[, 1], sum))
  gitd2_2 <- 2 * Zit * exp(Zit * (eta2_1 + eta2_2 * Zit))
  gitd2_gitcub2 <- -gitd2_2 * git2^3
  gid2_gicub2 <- as.numeric(tapply(gitd2_gitcub2, pindex[, 1], sum))
  gitd1sq_epsit_gitsq1 <- -Zit^2 * exp(Zit * (eta1_1 + eta1_2 * Zit)) * epsilon_it1 *
    git1^2
  gid1sq_epsi_gisq1 <- as.numeric(tapply(gitd1sq_epsit_gitsq1, pindex[, 1], sum))
  gitd1sq_epsit_gitsq2 <- -Zit^2 * exp(Zit * (eta2_1 + eta2_2 * Zit)) * epsilon_it2 *
    git2^2
  gid1sq_epsi_gisq2 <- as.numeric(tapply(gitd1sq_epsit_gitsq2, pindex[, 1], sum))
  gitd2sq_gitcub1 <- -2 * Zit^2 * exp(Zit * (eta1_1 + eta1_2 * Zit)) * git1^3
  gid2sq_gicub1 <- as.numeric(tapply(gitd2sq_gitcub1, pindex[, 1], sum))
  gitd2sq_gitcub2 <- -2 * Zit^2 * exp(Zit * (eta2_1 + eta2_2 * Zit)) * git2^3
  gid2sq_gicub2 <- as.numeric(tapply(gitd2sq_gitcub2, pindex[, 1], sum))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1_h <- exp(Wv1 * TT/2)
  ewv2_h <- exp(Wv2 * TT/2)
  ewz <- exp(Wz)
  sigmasq1 <- (ewu1 * gisq1 + ewv1)
  sigmasq2 <- (ewu2 * gisq2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigmasq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigmasq2)
  ssq1 <- (sigmasq1 * sigmastar1)
  ssq2 <- (sigmasq2 * sigmastar2)
  musig1 <- (S * ewu1 * giepsi1/ssq1)
  musig2 <- (S * ewu2 * giepsi2/ssq2)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  wugsq1 <- (1 - ewu1 * gisq1/sigmasq1)
  wugsq2 <- (1 - ewu2 * gisq2/sigmasq2)
  ewvsq1 <- (1 - ewv1/sigmasq1)
  ewvsq2 <- (1 - ewv2/sigmasq2)
  wzdeno <- (1 - ewz/(1 + ewz))
  sigx1_1 <- (S * ewu1 * pmusig1 * giepsi1/ssq1 - dmusig1)
  sigx1_2 <- (S * ewu2 * pmusig2 * giepsi2/ssq2 - dmusig2)
  sigx2_1 <- (ewvsq1 * ewu1/sigmastar1)
  sigx2_2 <- (ewvsq2 * ewu2/sigmastar2)
  sigx3_1 <- (wugsq1 * ewv1 * pmusig1/ssq1)
  sigx3_2 <- (wugsq2 * ewv2 * pmusig2/ssq2)
  sigx4_1 <- (0.5 * (wugsq1 * ewv1/sigmastar1) + sigmastar1 * gisq1)
  sigx4_2 <- (0.5 * (wugsq2 * ewv2/sigmastar2) + sigmastar2 * gisq2)
  sigx5_1 <- exp(-(0.5 * (epsilon_isq1/ewv1 - (-musig1)^2)))
  sigx5_2 <- exp(-(0.5 * (epsilon_isq2/ewv2 - (-musig2)^2)))
  sigx6_1 <- ((1 + ewz) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h)
  sigx6_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx7_1 <- (sigx5_1 * ewz * pmusig1 * sigmastar1/sigx6_1)
  sigx7_2 <- (wzdeno * sigx5_2 * pmusig2 * sigmastar2/sigx6_2)
  sigx8 <- ((1 + ewz) * (2 * sigx7_2 + 2 * sigx7_1) * (2 * pi)^(TT/2) * ewu1_h *
    ewv1_h)
  sigx9_1 <- (0.5 * sigx3_1 + S * (1/ssq1 - sigx4_1 * ewu1/ssq1^2) * sigx1_1 *
    sigmastar1 * giepsi1)
  sigx9_2 <- (0.5 * sigx3_2 + S * (1/ssq2 - sigx4_2 * ewu2/ssq2^2) * sigx1_2 *
    sigmastar2 * giepsi2)
  sigx10_1 <- ((1 + ewz) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmusig1 * sigmastar1/sigx6_1^2)
  sigx10_2 <- ((2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmusig2 * sigmastar2/sigx6_2^2)
  sigx11_1 <- (sigx9_1 * ewu1/sigx6_1 - 0.5 * sigx10_1)
  sigx11_2 <- (sigx9_2 * ewu2/sigx6_2 - 0.5 * sigx10_2)
  sigx12_1 <- (S^2 * (0.5 * sigx2_1 + sigmastar1) * ewu1^2 * ewv1 * giepsi1^2/(ssq1^2 *
    sigmasq1 * sigmastar1))
  sigx12_2 <- (S^2 * (0.5 * sigx2_2 + sigmastar2) * ewu2^2 * ewv2 * giepsi2^2/(ssq2^2 *
    sigmasq2 * sigmastar2))
  sigx13_1 <- ((2 * sigx12_1 - epsilon_isq1/ewv1) * pmusig1)
  sigx13_2 <- ((2 * sigx12_2 - epsilon_isq2/ewv2) * pmusig2)
  sigx14_1 <- (ewvsq1 * ewu1 * ewv1 * pmusig1/ssq1)
  sigx14_2 <- (ewvsq2 * ewu2 * ewv2 * pmusig2/ssq2)
  sigx15_1 <- (TT * (1 + ewz) * (2 * pi)^(TT/2) * ewu1_h * ewv1_h * pmusig1 * sigmastar1/sigx6_1^2)
  sigx15_2 <- (TT * (2 * pi)^(TT/2) * ewu2_h * ewv2_h * pmusig2 * sigmastar2/sigx6_2^2)
  sigx16_1 <- (S * (0.5 * sigx2_1 + sigmastar1) * dmusig1 * ewu1 * ewv1 * giepsi1/ssq1^2 -
    0.5 * sigx13_1)
  sigx16_2 <- (S * (0.5 * sigx2_2 + sigmastar2) * dmusig2 * ewu2 * ewv2 * giepsi2/ssq2^2 -
    0.5 * sigx13_2)
  sigx17_1 <- ((sigx16_1 * sigmastar1 + 0.5 * sigx14_1)/sigx6_1 - 0.5 * sigx15_1)
  sigx17_2 <- ((sigx16_2 * sigmastar2 + 0.5 * sigx14_2)/sigx6_2 - 0.5 * sigx15_2)
  sigx18_1 <- (sigmastar1 - 0.5 * (ewu1 * ewv1/ssq1))
  sigx18_2 <- (sigmastar2 - 0.5 * (ewu2 * ewv2/ssq2))
  sigx19_1 <- (gid1_epsi_gisq1/ssq1 - ewu1 * sigx18_1 * gid2_gicub1 * giepsi1/ssq1^2)
  sigx19_2 <- (gid1_epsi_gisq2/ssq2 - ewu2 * sigx18_2 * gid2_gicub2 * giepsi2/ssq2^2)
  sigx20_1 <- (ewu1 * ewv1 * pmusig1 * gid2_gicub1/(sigmasq1^2 * sigmastar1))
  sigx20_2 <- (ewu2 * ewv2 * pmusig2 * gid2_gicub2/(sigmasq2^2 * sigmastar2))
  sigx21_1 <- (S * sigx1_1 * sigmastar1 * sigx19_1 - 0.5 * sigx20_1)
  sigx21_2 <- (S * sigx1_2 * sigmastar2 * sigx19_2 - 0.5 * sigx20_2)
  sigx22_1 <- (gid1sq_epsi_gisq1/ssq1 - ewu1 * sigx18_1 * gid2sq_gicub1 * giepsi1/ssq1^2)
  sigx22_2 <- (gid1sq_epsi_gisq2/ssq2 - ewu2 * sigx18_2 * gid2sq_gicub2 * giepsi2/ssq2^2)
  sigx23_1 <- (ewu1 * ewv1 * pmusig1 * gid2sq_gicub1/(sigmasq1^2 * sigmastar1))
  sigx23_2 <- (ewu2 * ewv2 * pmusig2 * gid2sq_gicub2/(sigmasq2^2 * sigmastar2))
  sigx24_1 <- (S * sigx1_1 * sigmastar1 * sigx22_1 - 0.5 * sigx23_1)
  sigx24_2 <- (S * sigx1_2 * sigmastar2 * sigx22_2 - 0.5 * sigx23_2)
  sigx25 <- ((2 * sigx7_2 + 2 * sigx7_1) * (2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx26 <- ((1 + ewz) * (2 * pi)^(TT/2) * ewu2_h * ewv2_h)
  sigx27 <- (wzdeno * sigx5_2 * pmusig2 * sigmastar2/sigx26)
  sigx28 <- (1/sigx6_1 - (2 * pi)^(TT/2) * ewu1_h * ewv1_h * ewz/sigx6_1^2)
  sigx29 <- (sigx28 * sigx5_1 * pmusig1 * sigmastar1)
  Xsig1 <- sweep(Xgi1, MARGIN = 1, STATS = (S^2 * ewu1 * giepsi1/sigmasq1), FUN = "*")
  Xsig2 <- sweep(Xgi2, MARGIN = 1, STATS = (S^2 * ewu2 * giepsi2/sigmasq2), FUN = "*")
  Xsig3 <- sweep(Xgi1, MARGIN = 1, STATS = S * dmusig1 * ewu1/ssq1, FUN = "*")
  Xsig4 <- sweep(Xgi2, MARGIN = 1, STATS = S * dmusig2 * ewu2/ssq2, FUN = "*")
  Xsig5 <- sweep((Xepsi_i1 - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv1), FUN = "*")
  Xsig6 <- sweep((Xepsi_i2 - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv2), FUN = "*")
  gradll <- cbind(sweep((0.5 * Xsig5 + Xsig3), MARGIN = 1, STATS = -(2 * (sigx5_1 *
    ewz * sigmastar1/sigx8)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 *
    (sigx11_1 * sigx5_1 * ewz/(2 * sigx7_2 + 2 * sigx7_1)), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = 2 * (sigx17_1 * sigx5_1 * ewz/(2 * sigx7_2 + 2 * sigx7_1)),
    FUN = "*"), 2 * (sigx5_1 * ewu1 * ewz * sigx21_1/sigx8), 2 * (sigx5_1 * ewu1 *
    ewz * sigx24_1/sigx8), sweep((0.5 * Xsig6 + Xsig4), MARGIN = 1, STATS = -(2 *
    (wzdeno * sigx5_2 * sigmastar2/sigx25)), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = 2 * (sigx11_2 * wzdeno * sigx5_2/(2 * sigx7_2 + 2 * sigx7_1)), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = 2 * (sigx17_2 * wzdeno * sigx5_2/(2 * sigx7_2 +
      2 * sigx7_1)), FUN = "*"), 2 * (wzdeno * sigx5_2 * ewu2 * sigx21_2/sigx25),
    2 * (wzdeno * sigx5_2 * ewu2 * sigx24_2/sigx25), sweep(Zvar, MARGIN = 1,
      STATS = (2 * sigx29 - 2 * sigx27) * ewz/(2 * sigx7_2 + 2 * sigx7_1),
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
LCM2ChnormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax, stepmax,
  tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psLCMfhalfnorm2C_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
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
  startLoglik <- sum(pLCMhalfnormlike2C_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("LCM Panel K90-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(pLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, gr = function(parm) {
    -colSums(pgradLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pLCMhalfnormlike2C_k90,
    grad = pgradLCMhalfnormlike2C_k90, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, Zvar = Zvar, nZHvar = nZHvar,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, hs = function(parm) {
      as(calculus::jacobian(function(parm) -colSums(pgradLCMhalfnormlike2C_k90(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(parm)),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, Zvar = Zvar, nZHvar = nZHvar, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar))
    }, gradient = function(parm) {
      -colSums(pgradLCMhalfnormlike2C_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar,
        nZHvar = nZHvar))
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradLCMhalfnormlike2C_k90(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_k90(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$par))
    }
    if (method == "sr1") {
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(pgradLCMhalfnormlike2C_k90(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)), unname(mleObj$solution))
    }
  }
  mleObj$logL_OBS <- pLCMhalfnormlike2C_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- pgradLCMhalfnormlike2C_k90(parm = mlParam, nXvar = nXvar,
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
pLCM2Chalfnormeff_k90 <- function(object, level) {
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
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar +
    5):(2 * object$nXvar + 2 * object$nuZUvar + 2 * object$nvZVvar + object$nZHvar +
    4)]
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
  git1 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1_1 * seq(1:x) + eta1_2 *
    (seq(1:x))^2)))
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
  git2 <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta2_1 * seq(1:x) + eta2_2 *
    (seq(1:x))^2)))
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
