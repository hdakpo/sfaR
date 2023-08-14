################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Type:      - Kumbakhar 1990 (K90)                                            #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = (1 + exp(eta1 * t + eta2 * t^2))^(-1)                  #
# Convolution: Truncated-skewed-laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for tsl-normal distribution
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
#' @noRd
ptslnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar, Xvar,
  pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 3]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  mustar2 <- -((1 + lambda) * exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  ll <- log(1 + lambda) - 1/2 * Wu - (TT - 1)/2 * Wv - 1/2 * log(gisq) - (TT -
    1)/2 * log(2 * pi) - log(2 * lambda + 1) + log(2 * exp(-1/2 * (epsilon_isq/exp(Wv) -
    (mustar1/sigmastar)^2)) * pnorm(mustar1/sigmastar) - exp(-1/2 * (epsilon_isq/exp(Wv) -
    (mustar2/sigmastar)^2)) * pnorm(mustar2/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for tsl-normal distribution
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
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
psttslnorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- csttslnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTsl <- NULL
  } else {
    cat("Initialization: SFA + truncated-skewed-laplace-normal distribution...\n")
    initTsl <- maxLik::maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradtslnormlike,
      hess = chesstslnormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initTsl$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3], eta1 = 0.001, eta2 = 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "lambda", "eta1", "eta2")
  return(list(StartVal = StartVal, initTsl = initTsl))
}

# Gradient of the likelihood function ----------
#' gradient for tsl-normal distribution
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
#' @noRd
pgradtslnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 3]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1 <- Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd1_epsit_gitsq <- -gitd1 * epsilon_it * git^2
  gid1_epsi_gisq <- as.numeric(tapply(gitd1_epsit_gitsq, pindex[, 1], sum))
  gitd2 <- 2 * Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd2_gitcub <- -gitd2 * git^3
  gid2_gicub <- as.numeric(tapply(gitd2_gitcub, pindex[, 1], sum))
  gitd1sq_epsit_gitsq <- -Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it *
    git^2
  gid1sq_epsi_gisq <- as.numeric(tapply(gitd1sq_epsit_gitsq, pindex[, 1], sum))
  gitd2sq_gitcub <- -2 * Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gid2sq_gicub <- as.numeric(tapply(gitd2sq_gitcub, pindex[, 1], sum))
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sgsq <- (sqrt(ewv/gisq) * gisq)
  llepsi <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- (((1 + lambda) * ewv/ewu_h + S * giepsi)/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (2 * (sigx3 * pmusig1) - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * ((ewv/ewu_h + S * giepsi) *
    pmusig1))
  sigx6 <- (0.5 * (dmusig2 * ewv/sqrt(ewv/gisq)) - 0.5 * (((1 + lambda) * ewv/ewu_h +
    S * giepsi) * pmusig2))
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx12 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (((1 + lambda) * ewv/ewu_h + S *
    giepsi)/sgsq^2))
  sigx10 <- (llepsi * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx13 <- (2 * (sigx9 * sigx3) - (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) *
    sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * gid1_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx16 <- (((1 + lambda) * ewv/ewu_h + S * giepsi) * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * gid1_epsi_gisq/sgsq - ((1 + lambda) * ewv/ewu_h + S * giepsi) *
    (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx18 <- (2 * (sigx14 * sigx3 * sigx15) - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * gid1sq_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx20 <- (S * gid1sq_epsi_gisq/sgsq - ((1 + lambda) * ewv/ewu_h + S * giepsi) *
    (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx21 <- (2 * (sigx14 * sigx3 * sigx19) - sigx16 * sigx2 * sigx20)
  llepsi2 <- (llepsi * pmusig2 - dmusig2 * ewv/sqrt(ewv/gisq))
  Xsig0 <- sweep(Xgi, MARGIN = 1, STATS = (S * llepsi/gisq), FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * (ewv/ewu_h + S * giepsi)/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sgsq, FUN = "*")
  Xsig4 <- sweep(Xepsi_i - 2 * Xsig0, MARGIN = 1, STATS = (pmusig2/ewv), FUN = "*")
  Xsig5 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sgsq, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig4 + Xsig2), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1, STATS = sigx3, FUN = "*")
  gradll <- cbind(sweep(Xsig7 - 2 * Xsig8, MARGIN = 1, STATS = 1/sigx4, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ((2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) *
      sigx2)/(sigx4 * ewu_h * gisq) - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = (sigx13/sigx4 - 0.5 * (TT - 1)), FUN = "*"), 1/(1 + lambda) -
      (llepsi2 * sigx2/(sigx4 * ewu_h * gisq) + 2/(1 + 2 * lambda)), sigx18/sigx4 -
      0.5 * (gid2_gicub/gisq), sigx21/sigx4 - 0.5 * (gid2sq_gicub/gisq))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for tsl-normal distribution
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
#' @noRd
phesstslnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 3]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1 <- Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd1_epsit_gitsq <- -gitd1 * epsilon_it * git^2
  gid1_epsi_gisq <- as.numeric(tapply(gitd1_epsit_gitsq, pindex[, 1], sum))
  gitd2 <- 2 * Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd2_gitcub <- -gitd2 * git^3
  gid2_gicub <- as.numeric(tapply(gitd2_gitcub, pindex[, 1], sum))
  gitd1sq_epsit_gitsq <- -Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it *
    git^2
  gid1sq_epsi_gisq <- as.numeric(tapply(gitd1sq_epsit_gitsq, pindex[, 1], sum))
  gitd2sq_gitcub <- -2 * Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gid2sq_gicub <- as.numeric(tapply(gitd2sq_gitcub, pindex[, 1], sum))
  gd1Zsqt_gd1t_epsit_gitsq <- -Zit^2 * (1 - 2 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it * git^2
  gd1Zsq_gd1_epsi_gisq <- as.numeric(tapply(gd1Zsqt_gd1t_epsit_gitsq, pindex[,
    1], sum))
  gd1Zcubt_gd1t_epsit_gitsq <- -Zit^3 * (1 - 2 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it * git^2
  gd1Zcub_gd1_epsi_gisq <- as.numeric(tapply(gd1Zcubt_gd1t_epsit_gitsq, pindex[,
    1], sum))
  gd1Zfourt_gd1t_epsit_gitsq <- -Zit^4 * (1 - 2 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * epsilon_it * git^2
  gd1Zfour_gd1_epsi_gisq <- as.numeric(tapply(gd1Zfourt_gd1t_epsit_gitsq, pindex[,
    1], sum))
  gd2Zsqt_gd1t_gitcub <- -2 * Zit^2 * (1 - 3 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gd2Zsq_gd1_gicub <- as.numeric(tapply(gd2Zsqt_gd1t_gitcub, pindex[, 1], sum))
  gd2Zcubt_gd1t_gitcub <- -2 * Zit^3 * (1 - 3 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gd2Zcub_gd1_gicub <- as.numeric(tapply(gd2Zcubt_gd1t_gitcub, pindex[, 1], sum))
  gd2Zfourt_gd1t_gitcub <- -2 * Zit^4 * (1 - 3 * (exp(Zit * (eta1 + eta2 * Zit)) *
    git)) * exp(Zit * (eta1 + eta2 * Zit)) * git^3
  gd2Zfour_gd1_gicub <- as.numeric(tapply(gd2Zfourt_gd1t_gitcub, pindex[, 1], sum))
  Xit_gitd1_gitsq <- sweep(Xvar, MARGIN = 1, STATS = gitd1 * git^2, FUN = "*")
  Xi_gid1_gisq <- apply(Xit_gitd1_gitsq, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xit_gitd1sq_gitsq <- sweep(Xvar, MARGIN = 1, STATS = Zit^2 * exp(Zit * (eta1 +
    eta2 * Zit)) * git^2, FUN = "*")
  Xi_gid1sq_gisq <- apply(Xit_gitd1sq_gitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sqgiv <- sqrt(ewv/gisq)
  sgsq <- (sqgiv * gisq)
  llepsi <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- (((1 + lambda) * ewv/ewu_h + S * giepsi)/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (2 * (sigx3 * pmusig1) - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqgiv) - 0.5 * ((ewv/ewu_h + S * giepsi) * pmusig1))
  sigx6 <- (0.5 * (dmusig2 * ewv/sqgiv) - 0.5 * (((1 + lambda) * ewv/ewu_h + S *
    giepsi) * pmusig2))
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqgiv)
  sigx12 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (((1 + lambda) * ewv/ewu_h + S *
    giepsi)/sgsq^2))
  sigx10 <- (llepsi * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx27 <- (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqgiv)
  sigx13 <- (2 * (sigx9 * sigx3) - sigx27 * sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * gid1_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqgiv - 0.5 *
    (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx16 <- (((1 + lambda) * ewv/ewu_h + S * giepsi) * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * gid1_epsi_gisq/sgsq - ((1 + lambda) * ewv/ewu_h + S * giepsi) *
    (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx18 <- (2 * (sigx14 * sigx3 * sigx15) - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * gid1sq_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqgiv - 0.5 *
    (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx20 <- (S * gid1sq_epsi_gisq/sgsq - ((1 + lambda) * ewv/ewu_h + S * giepsi) *
    (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx21 <- (2 * (sigx14 * sigx3 * sigx19) - sigx16 * sigx2 * sigx20)
  llepsi2 <- (llepsi * pmusig2 - dmusig2 * ewv/sqgiv)
  pdmusig2 <- (pmusig2 - llepsi * dmusig2/sgsq)
  sigx22 <- (llepsi^2 * pmusig2/gisq + ewv * pdmusig2)
  sigx23 <- (0.5 * (llepsi * dmusig2/sgsq) + 0.5 * pdmusig2)
  sigx24 <- (pmusig1 - dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq)
  sigx25 <- (sigx24/sqgiv + dmusig1 * (ewv/ewu_h + S * giepsi)/ewv)
  sigx26 <- (llepsi * dmusig2/ewv + pdmusig2/sqgiv)
  sigx28 <- (sigx7 * (ewv/ewu_h + S * giepsi)/(ewv * gisq) + 0.5/sgsq^2)
  sigx29 <- (0.5 * (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq) + 0.5 * sigx24)
  Xsig0 <- sweep(Xgi, MARGIN = 1, STATS = (S * llepsi/gisq), FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * (ewv/ewu_h + S * giepsi)/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sgsq, FUN = "*")
  Xsig4 <- sweep(Xepsi_i - 2 * Xsig0, MARGIN = 1, STATS = (pmusig2/ewv), FUN = "*")
  Xsig5 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sgsq, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig4 + Xsig2), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1, STATS = sigx3, FUN = "*")
  Xsig9 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx25/gisq, FUN = "*")
  Xsig10 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx14/ewv), FUN = "*")
  Xsig11 <- sweep(Xgi, MARGIN = 1, STATS = (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2,
    FUN = "*")
  Xsig12 <- sweep(Xi_gid1_gisq, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig13 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx26/gisq, FUN = "*")
  Xsig14 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx16/ewv), FUN = "*")
  Xsig15 <- sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig16 <- sweep(Xgi, MARGIN = 1, STATS = (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2,
    FUN = "*")
  Xsig17 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx19, FUN = "*")
  Xsig18 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx15, FUN = "*")
  Xsig19 <- sweep((Xsig12 - Xsig11), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig20 <- sweep((Xsig15 - Xsig16), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig21 <- sweep((Xsig18 + Xsig19), MARGIN = 1, STATS = (sigx3), FUN = "*")
  Xsig22 <- sweep((Xsig17 + Xsig20), MARGIN = 1, STATS = (sigx3), FUN = "*")
  Xsig23 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx18/sigx4, FUN = "*")
  Xsig24 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx21/sigx4, FUN = "*")
  Xsig25 <- sweep((Xsig13 - 0.5 * Xsig14), MARGIN = 1, STATS = sigx17, FUN = "*")
  Xsig26 <- sweep((Xsig13 - 0.5 * Xsig14), MARGIN = 1, STATS = sigx20, FUN = "*")
  Xsig27 <- sweep((Xsig12 - Xsig11), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig28 <- sweep((Xsig15 - Xsig16), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig29 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = 0.5 * (llepsi2/ewv),
    FUN = "*")
  Xsig30 <- sweep(Xgi, MARGIN = 1, STATS = S * pmusig2, FUN = "*")
  Xsig31 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = llepsi2 * ewu_h * gisq/(sigx4 *
    ewu_h * gisq)^2, FUN = "*")
  Xsig32 <- sweep((Xsig30 - Xsig29), MARGIN = 1, STATS = 1/(sigx4 * ewu_h * gisq),
    FUN = "*")
  Xsig33 <- sweep(Xgi, MARGIN = 1, STATS = S * (2 * sigx10 + epsilon_isq/ewv) *
    dmusig2/sgsq, FUN = "*")
  Xsig34 <- sweep(Xgi, MARGIN = 1, STATS = (S * ((1 + lambda)/(ewu_h * gisq) -
    llepsi/sgsq^2)), FUN = "*")
  Xsig35 <- sweep(Xgi, MARGIN = 1, STATS = (S * (1/(ewu_h * gisq) - (ewv/ewu_h +
    S * giepsi)/sgsq^2)), FUN = "*")
  Xsig36 <- sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*")
  Xsig37 <- sweep((2 * Xsig35 + Xsig36), MARGIN = 1, STATS = pmusig1, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = S * sigx8 * dmusig1/sgsq, FUN = "*")
  Xsig38 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx28 * dmusig1 * ewv/sqgiv, FUN = "*")
  Xsig39 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx9/ewv), FUN = "*")
  Xsig40 <- sweep((2 * Xsig34 + Xsig36), MARGIN = 1, STATS = pmusig2, FUN = "*")
  Xsig41 <- sweep(Xgi, MARGIN = 1, STATS = S * (llepsi * sigx12/(ewv * gisq) +
    0.5/sgsq^2) * dmusig2 * ewv/sqgiv, FUN = "*")
  Xsig42 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx27/ewv), FUN = "*")
  Xsig43 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx13/sigx4, FUN = "*")
  Xsig44 <- sweep((0.5 * (Xsig40 - Xsig33) + Xsig41 - 0.5 * Xsig42), MARGIN = 1,
    STATS = sigx2, FUN = "*")
  Xsig45 <- sweep((0.5 * Xsig37 + Xsig38 - 0.5 * Xsig39), MARGIN = 1, STATS = (sigx3),
    FUN = "*")
  Xsig46 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = 0.5 * (sigx6/ewv),
    FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = S * sigx23, FUN = "*")
  Xsig47 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = 0.5 * (sigx5/ewv),
    FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = S * sigx29, FUN = "*")
  Xsig48 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = (2 * (sigx5 * sigx3) -
    sigx6 * (1 + lambda) * sigx2) * ewu_h * gisq/(sigx4 * ewu_h * gisq)^2, FUN = "*")
  Xsig49 <- sweep(Xsig46, MARGIN = 1, STATS = (1 + lambda) * sigx2, FUN = "*") -
    sweep(Xsig47, MARGIN = 1, STATS = 2 * (sigx3), FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 3, ncol = nXvar + nuZUvar +
    nvZVvar + 3)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar * pmusig2 * sigx2/ewv/sigx4))) - crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * 2 * (S^2/gisq) * pmusig2 * sigx2/ewv/sigx4, FUN = "*"), Xgi) -
    crossprod(sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = wHvar * S * dmusig2/sgsq *
      sigx2/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 * Xsig4 +
    Xsig2), MARGIN = 1, STATS = wHvar * sigx2/ewv/sigx4, FUN = "*"), (Xepsi_i -
    2 * Xsig0)) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S^2 * llepsi *
    dmusig2/(sqgiv * gisq^2) * sigx2/ewv/sigx4, FUN = "*"), Xgi)) - 2 * (0.5 *
    (sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar * pmusig1 *
      sigx3/ewv/sigx4))) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
      2 * (S^2/gisq) * pmusig1 * sigx3/ewv/sigx4, FUN = "*"), Xgi) - crossprod(sweep((Xepsi_i -
      2 * Xsig1), MARGIN = 1, STATS = wHvar * S * dmusig1/sgsq * sigx3/ewv/sigx4,
      FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1,
    STATS = wHvar * sigx3/ewv/sigx4, FUN = "*"), (Xepsi_i - 2 * Xsig1)) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S^2 * dmusig1 * (ewv/ewu_h + S * giepsi)/(sqgiv *
      gisq^2) * sigx3/ewv/sigx4, FUN = "*"), Xgi))) - crossprod(sweep((Xsig7 -
    2 * Xsig8), MARGIN = 1, STATS = wHvar/sigx4^2, FUN = "*"), (Xsig7 - 2 * Xsig8))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xsig49, MARGIN = 1,
    STATS = wHvar/(sigx4 * ewu_h * gisq), FUN = "*") - sweep(Xsig48, MARGIN = 1,
    STATS = wHvar, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep((2 *
    Xsig45 - (Xsig43 + Xsig44)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep((Xsig32 - Xsig31),
    MARGIN = 1, STATS = -wHvar * (sigx2), FUN = "*"))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 2] <- colSums(sweep((2 * Xsig21 -
    (Xsig23 + (Xsig25 + Xsig27) * sigx2)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 3] <- colSums(sweep((2 * Xsig22 -
    (Xsig24 + (Xsig26 + Xsig28) * sigx2)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.25 * (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq) -
      0.5 * (0.5 * (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq) - 0.5 * pmusig1)) *
      ewv - 0.5 * (sigx5 * (ewv/ewu_h + S * giepsi)/gisq)) * sigx3) - ((0.25 *
      (llepsi * dmusig2/sgsq) - 0.5 * (0.5 * (llepsi * dmusig2/sgsq) - 0.5 *
      pmusig2)) * ewv - 0.5 * (llepsi * sigx6/gisq)) * (1 + lambda)^2 * sigx2)/(sigx4 *
      ewu_h^2 * gisq) - ((2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) * sigx2)/gisq +
      0.5 * (sigx4 * ewu_h)) * (2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) *
      sigx2) * gisq/(sigx4 * ewu_h * gisq)^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (2 * ((0.5 *
    (0.5 * (sigx8 * dmusig1 * ewv/(ewu_h * sqgiv * gisq)) + 2 * (((0.25 * (ewv/(ewu_h *
      sgsq^2)) - 0.5 * (ewu_h * gisq/(ewu_h * gisq)^2)) * (ewv/ewu_h + S *
      giepsi) - 0.5 * (sigx7 * ewv/ewu_h)) * pmusig1)) - (((0.25 * (ewv/sgsq^2) +
    0.5 * (sigx7 * (ewv/ewu_h + S * giepsi)/gisq))/ewu_h - 0.5 * (ewu_h * gisq/(ewu_h *
    gisq)^2)) * dmusig1 * ewv/sqgiv + 0.5 * (sigx9 * (ewv/ewu_h + S * giepsi)/(ewu_h *
    gisq)))) * sigx3) - ((0.5 * (0.5 * ((2 * sigx10 + epsilon_isq/ewv) * dmusig2 *
    ewv/(ewu_h * sqgiv * gisq)) + 2 * ((llepsi * (0.25 * (ewv/(ewu_h * sgsq^2)) -
    0.5 * (ewu_h * gisq/(ewu_h * gisq)^2)) - 0.5 * (sigx12 * ewv/ewu_h)) * pmusig2)) -
    (((0.25 * (ewv/sgsq^2) + 0.5 * (llepsi * sigx12/gisq))/ewu_h - 0.5 * (ewu_h *
      gisq/(ewu_h * gisq)^2)) * dmusig2 * ewv/sqgiv + 0.5 * (llepsi * sigx27/(ewu_h *
      gisq)))) * (1 + lambda) * sigx2 + sigx13 * (2 * (sigx5 * sigx3) - sigx6 *
    (1 + lambda) * sigx2)/(sigx4 * ewu_h * gisq)))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (llepsi2 * ((2 * (sigx5 * sigx3) - sigx6 * (1 +
      lambda) * sigx2)/gisq + 0.5 * (sigx4 * ewu_h)) * gisq/(sigx4 * ewu_h *
      gisq)^2 + (0.5 * (llepsi2 * llepsi/gisq) + 0.5 * (ewv * pmusig2)) * (1 +
      lambda)/(sigx4 * ewu_h^2 * gisq)) * sigx2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 2] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((((0.5 * (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq) -
      0.5 * pmusig1) * ewv/sqgiv - (0.5 * sigx14 + 0.5 * dmusig1) * (ewv/ewu_h +
      S * giepsi)) * sigx15/gisq + 0.5 * (sigx14 * ewv * (sqgiv - 0.5 * (ewv/sgsq)) *
      gid2_gicub/sgsq^2)) * sigx3) - ((((0.5 * (llepsi * dmusig2/sgsq) - 0.5 *
      pmusig2) * ewv/sqgiv - llepsi * (0.5 * sigx16 + 0.5 * dmusig2)) * sigx17/gisq +
      0.5 * (sigx16 * ewv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)) *
      (1 + lambda) * sigx2 + sigx18 * (2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) *
      sigx2)/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 3] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((((0.5 * (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq) -
      0.5 * pmusig1) * ewv/sqgiv - (0.5 * sigx14 + 0.5 * dmusig1) * (ewv/ewu_h +
      S * giepsi)) * sigx19/gisq + 0.5 * (sigx14 * ewv * (sqgiv - 0.5 * (ewv/sgsq)) *
      gid2sq_gicub/sgsq^2)) * sigx3) - ((((0.5 * (llepsi * dmusig2/sgsq) -
      0.5 * pmusig2) * ewv/sqgiv - llepsi * (0.5 * sigx16 + 0.5 * dmusig2)) *
      sigx20/gisq + 0.5 * (sigx16 * ewv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)) *
      (1 + lambda) * sigx2 + sigx21 * (2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) *
      sigx2)/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((0.5 * (sigx9 * sigx8) + 0.5 * ((2 * ((sigx7/ewu_h - 0.5 * ((1/ewu_h -
      (ewv/ewu_h + S * giepsi) * gisq/sgsq^2) * (ewv/ewu_h + S * giepsi)/sgsq^2)) *
      ewv) - epsilon_isq/ewv) * pmusig1 - sigx7 * sigx8 * dmusig1 * ewv/sqgiv) -
      (1/(ewu_h * gisq) - ((sigx7^2 + 0.5/sgsq^2) * (ewv/ewu_h + S * giepsi) +
        0.5 * ((1/ewu_h - (ewv/ewu_h + S * giepsi) * gisq/sgsq^2) * ewv/sgsq^2) +
        0.5 * sigx7)) * dmusig1 * ewv/sqgiv) * sigx3) - ((0.5 * (sigx27 *
      (2 * sigx10 + epsilon_isq/ewv)) + 0.5 * ((2 * ((sigx12 * (1 + lambda)/ewu_h -
      0.5 * (llepsi * ((1 + lambda)/ewu_h - llepsi * gisq/sgsq^2)/sgsq^2)) *
      ewv) - epsilon_isq/ewv) * pmusig2 - sigx12 * (2 * sigx10 + epsilon_isq/ewv) *
      dmusig2 * ewv/sqgiv) - ((1 + lambda)/(ewu_h * gisq) - ((sigx12^2 + 0.5/sgsq^2) *
      llepsi + 0.5 * (((1 + lambda)/ewu_h - llepsi * gisq/sgsq^2) * ewv/sgsq^2) +
      0.5 * sigx12)) * dmusig2 * ewv/sqgiv) * sigx2 + sigx13^2/sigx4))/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (((((1 + lambda) *
    pmusig2/ewu_h - 0.5 * (dmusig2/sqgiv)) * ewv + 0.5 * (llepsi2 * (2 * sigx10 +
    epsilon_isq/ewv)))/(sigx4 * ewu_h * gisq) - llepsi2 * sigx13 * ewu_h * gisq/(sigx4 *
    ewu_h * gisq)^2) * sigx2), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    2] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (2 * (((((pmusig1/ewu_h -
    sigx7 * dmusig1 * (ewv/ewu_h + S * giepsi)/sqgiv)/gisq - 0.5 * ((ewv/ewu_h +
    S * giepsi) * pmusig1/sgsq^2)) * ewv/sqgiv + sigx7 * dmusig1 * (ewv/ewu_h +
    S * giepsi) + 0.5 * (sigx14 * sigx8)) * sigx15 - ((((0.5/gisq - 0.5 * (1/gisq -
    0.5 * (ewv/sgsq^2)))/sqgiv - (sqgiv - 0.5 * (ewv/sgsq)) * gisq/sgsq^2) *
    (ewv/ewu_h + S * giepsi) + (sqgiv - 0.5 * (ewv/sgsq))/ewu_h) * gid2_gicub +
    0.5 * (S * gid1_epsi_gisq/sqgiv)) * sigx14 * ewv/sgsq^2) * sigx3) - ((((((1 +
    lambda) * pmusig2/ewu_h - llepsi * sigx12 * dmusig2/sqgiv)/gisq - 0.5 * (llepsi *
    pmusig2/sgsq^2)) * ewv/sqgiv + llepsi * sigx12 * dmusig2 + 0.5 * (sigx16 *
    (2 * sigx10 + epsilon_isq/ewv))) * sigx17 - ((((0.5/gisq - 0.5 * (1/gisq -
    0.5 * (ewv/sgsq^2)))/sqgiv - (sqgiv - 0.5 * (ewv/sgsq)) * gisq/sgsq^2) *
    llepsi + (1 + lambda) * (sqgiv - 0.5 * (ewv/sgsq))/ewu_h) * gid2_gicub +
    0.5 * (S * gid1_epsi_gisq/sqgiv)) * sigx16 * ewv/sgsq^2) * sigx2 + sigx18 *
    sigx13/sigx4))/sigx4, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    3] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (2 * (((((pmusig1/ewu_h -
    sigx7 * dmusig1 * (ewv/ewu_h + S * giepsi)/sqgiv)/gisq - 0.5 * ((ewv/ewu_h +
    S * giepsi) * pmusig1/sgsq^2)) * ewv/sqgiv + sigx7 * dmusig1 * (ewv/ewu_h +
    S * giepsi) + 0.5 * (sigx14 * sigx8)) * sigx19 - ((((0.5/gisq - 0.5 * (1/gisq -
    0.5 * (ewv/sgsq^2)))/sqgiv - (sqgiv - 0.5 * (ewv/sgsq)) * gisq/sgsq^2) *
    (ewv/ewu_h + S * giepsi) + (sqgiv - 0.5 * (ewv/sgsq))/ewu_h) * gid2sq_gicub +
    0.5 * (S * gid1sq_epsi_gisq/sqgiv)) * sigx14 * ewv/sgsq^2) * sigx3) - ((((((1 +
    lambda) * pmusig2/ewu_h - llepsi * sigx12 * dmusig2/sqgiv)/gisq - 0.5 * (llepsi *
    pmusig2/sgsq^2)) * ewv/sqgiv + llepsi * sigx12 * dmusig2 + 0.5 * (sigx16 *
    (2 * sigx10 + epsilon_isq/ewv))) * sigx20 - ((((0.5/gisq - 0.5 * (1/gisq -
    0.5 * (ewv/sgsq^2)))/sqgiv - (sqgiv - 0.5 * (ewv/sgsq)) * gisq/sgsq^2) *
    llepsi + (1 + lambda) * (sqgiv - 0.5 * (ewv/sgsq))/ewu_h) * gid2sq_gicub +
    0.5 * (S * gid1sq_epsi_gisq/sqgiv)) * sigx16 * ewv/sgsq^2) * sigx2 + sigx21 *
    sigx13/sigx4))/sigx4, FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 1] <- sum(wHvar *
    (4/(1 + 2 * lambda)^2 - (((llepsi2 * llepsi/gisq + ewv * pmusig2)/(sigx4 *
      ewu_h^2 * gisq) + llepsi2^2 * sigx2/(sigx4 * ewu_h * gisq)^2) * sigx2 +
      1/(1 + lambda)^2)))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 2] <- sum(wHvar *
    (((llepsi2 * sigx18/sigx4 - sigx22 * sigx17/sqgiv)/gisq + sigx16 * ewv *
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2) * sigx2/(sigx4 * ewu_h)))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 3] <- sum(wHvar *
    (((llepsi2 * sigx21/sigx4 - sigx22 * sigx20/sqgiv)/gisq + sigx16 * ewv *
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2) * sigx2/(sigx4 * ewu_h)))
  hessll[nXvar + nuZUvar + nvZVvar + 2, nXvar + nuZUvar + nvZVvar + 2] <- sum(wHvar *
    ((2 * (((((ewv/ewu_h + S * giepsi) * sigx15/(ewv * gisq) - (sqgiv - 0.5 *
      (ewv/sgsq)) * gid2_gicub/sgsq^2) * (ewv/ewu_h + S * giepsi) * pmusig1 +
      (S * pmusig1 * gid1_epsi_gisq - dmusig1 * (ewv/ewu_h + S * giepsi) *
        sigx15)/sgsq) * sigx15 + sigx14 * (S * (gd1Zsq_gd1_epsi_gisq/sgsq -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub * gid1_epsi_gisq/sgsq^2) - (((0.5 *
      ((sqgiv - 0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv * gisq^2)) * ewv * (ewv/ewu_h +
      S * giepsi) * gid2_gicub + S * (sqgiv - 0.5 * (ewv/sgsq)) * gid1_epsi_gisq) *
      gid2_gicub + (ewv/ewu_h + S * giepsi) * (sqgiv - 0.5 * (ewv/sgsq)) *
      (gd2Zsq_gd1_gicub - 2 * (sqgiv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub^2 *
        gisq/sgsq^2)))/sgsq^2)) * sigx3) - ((((llepsi * sigx17/(ewv * gisq) -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2) * llepsi * pmusig2 +
      (S * pmusig2 * gid1_epsi_gisq - llepsi * dmusig2 * sigx17)/sgsq) * sigx17 +
      sigx16 * (S * (gd1Zsq_gd1_epsi_gisq/sgsq - (sqgiv - 0.5 * (ewv/sgsq)) *
        gid2_gicub * gid1_epsi_gisq/sgsq^2) - ((llepsi * (0.5 * ((sqgiv -
        0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv * gisq^2)) * ewv * gid2_gicub +
        S * (sqgiv - 0.5 * (ewv/sgsq)) * gid1_epsi_gisq) * gid2_gicub + llepsi *
        (sqgiv - 0.5 * (ewv/sgsq)) * (gd2Zsq_gd1_gicub - 2 * (sqgiv * (sqgiv -
        0.5 * (ewv/sgsq)) * gid2_gicub^2 * gisq/sgsq^2)))/sgsq^2)) * sigx2 +
      sigx18^2/sigx4))/sigx4 - 0.5 * ((gd2Zsq_gd1_gicub - gid2_gicub^2/gisq)/gisq)))
  hessll[nXvar + nuZUvar + nvZVvar + 2, nXvar + nuZUvar + nvZVvar + 3] <- sum(wHvar *
    ((2 * (((((ewv/ewu_h + S * giepsi) * sigx15/(ewv * gisq) - (sqgiv - 0.5 *
      (ewv/sgsq)) * gid2_gicub/sgsq^2) * (ewv/ewu_h + S * giepsi) * pmusig1 +
      (S * pmusig1 * gid1_epsi_gisq - dmusig1 * (ewv/ewu_h + S * giepsi) *
        sigx15)/sgsq) * sigx19 + sigx14 * (S * (gd1Zcub_gd1_epsi_gisq/sgsq -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub * gid1sq_epsi_gisq/sgsq^2) -
      (((0.5 * ((sqgiv - 0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv * gisq^2)) *
        ewv * (ewv/ewu_h + S * giepsi) * gid2_gicub + S * (sqgiv - 0.5 *
        (ewv/sgsq)) * gid1_epsi_gisq) * gid2sq_gicub + (ewv/ewu_h + S * giepsi) *
        (sqgiv - 0.5 * (ewv/sgsq)) * (gd2Zcub_gd1_gicub - 2 * (sqgiv * (sqgiv -
        0.5 * (ewv/sgsq)) * gid2_gicub * gid2sq_gicub * gisq/sgsq^2)))/sgsq^2)) *
      sigx3) - ((((llepsi * sigx17/(ewv * gisq) - (sqgiv - 0.5 * (ewv/sgsq)) *
      gid2_gicub/sgsq^2) * llepsi * pmusig2 + (S * pmusig2 * gid1_epsi_gisq -
      llepsi * dmusig2 * sigx17)/sgsq) * sigx20 + sigx16 * (S * (gd1Zcub_gd1_epsi_gisq/sgsq -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub * gid1sq_epsi_gisq/sgsq^2) -
      ((llepsi * (0.5 * ((sqgiv - 0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv *
        gisq^2)) * ewv * gid2_gicub + S * (sqgiv - 0.5 * (ewv/sgsq)) * gid1_epsi_gisq) *
        gid2sq_gicub + llepsi * (sqgiv - 0.5 * (ewv/sgsq)) * (gd2Zcub_gd1_gicub -
        2 * (sqgiv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2_gicub * gid2sq_gicub *
          gisq/sgsq^2)))/sgsq^2)) * sigx2 + sigx18 * sigx21/sigx4))/sigx4 -
      0.5 * ((gd2Zcub_gd1_gicub - gid2_gicub * gid2sq_gicub/gisq)/gisq)))
  hessll[nXvar + nuZUvar + nvZVvar + 3, nXvar + nuZUvar + nvZVvar + 3] <- sum(wHvar *
    ((2 * (((((ewv/ewu_h + S * giepsi) * sigx19/(ewv * gisq) - (sqgiv - 0.5 *
      (ewv/sgsq)) * gid2sq_gicub/sgsq^2) * (ewv/ewu_h + S * giepsi) * pmusig1 +
      (S * pmusig1 * gid1sq_epsi_gisq - dmusig1 * (ewv/ewu_h + S * giepsi) *
        sigx19)/sgsq) * sigx19 + sigx14 * (S * (gd1Zfour_gd1_epsi_gisq/sgsq -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub * gid1sq_epsi_gisq/sgsq^2) -
      (((0.5 * ((sqgiv - 0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv * gisq^2)) *
        ewv * (ewv/ewu_h + S * giepsi) * gid2sq_gicub + S * (sqgiv - 0.5 *
        (ewv/sgsq)) * gid1sq_epsi_gisq) * gid2sq_gicub + (ewv/ewu_h + S *
        giepsi) * (sqgiv - 0.5 * (ewv/sgsq)) * (gd2Zfour_gd1_gicub - 2 *
        (sqgiv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub^2 * gisq/sgsq^2)))/sgsq^2)) *
      sigx3) - ((((llepsi * sigx20/(ewv * gisq) - (sqgiv - 0.5 * (ewv/sgsq)) *
      gid2sq_gicub/sgsq^2) * llepsi * pmusig2 + (S * pmusig2 * gid1sq_epsi_gisq -
      llepsi * dmusig2 * sigx20)/sgsq) * sigx20 + sigx16 * (S * (gd1Zfour_gd1_epsi_gisq/sgsq -
      (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub * gid1sq_epsi_gisq/sgsq^2) -
      ((llepsi * (0.5 * ((sqgiv - 0.5 * (ewv/sgsq))/sgsq^2) - 0.5/(sqgiv *
        gisq^2)) * ewv * gid2sq_gicub + S * (sqgiv - 0.5 * (ewv/sgsq)) *
        gid1sq_epsi_gisq) * gid2sq_gicub + llepsi * (sqgiv - 0.5 * (ewv/sgsq)) *
        (gd2Zfour_gd1_gicub - 2 * (sqgiv * (sqgiv - 0.5 * (ewv/sgsq)) * gid2sq_gicub^2 *
          gisq/sgsq^2)))/sgsq^2)) * sigx2 + sigx21^2/sigx4))/sigx4 - 0.5 *
      ((gd2Zfour_gd1_gicub - gid2sq_gicub^2/gisq)/gisq)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for tsl-normal distribution
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
tslnormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttslnorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    initTsl <- start_st$initTsl
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ptslnormlike_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel K90 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(ptslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptslnormlike_k90,
    grad = pgradtslnormlike_k90, hess = phesstslnormlike_k90, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(ptslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(ptslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(ptslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradtslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesstslnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtslnormlike_k90(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phesstslnormlike_k90(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstslnormlike_k90(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptslnormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtslnormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, initTsl = initTsl))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for tsl-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
ptslnormeff_k90 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 2]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 3]
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
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar_p)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar_p)))
  TT <- as.numeric(table(pindex[, 1]))
  epsilon_it <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  mustar2 <- -((1 + lambda) * exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  a <- mustar1/sigmastar
  b <- mustar2/sigmastar
  A <- -1/2 * (epsilon_isq/exp(Wv) - a^2)
  B <- -1/2 * (epsilon_isq/exp(Wv) - b^2)
  u <- (exp(A) * (dnorm(a) * sigmastar + mustar1 * pnorm(a)) - exp(B) * (dnorm(b) *
    sigmastar + mustar2 * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) * pnorm(b))
  res <- data.frame(levels(pindex[, 1]), u = u, mustar1 = mustar1, mustar2 = mustar2,
    sigmastar = sigmastar, A = A, B = B, a = a, b = b)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teBC <- (exp(res$A) * exp(1/2 * res$sigmastar^2 * git^2 - res$mustar1 *
      git) * pnorm(res$a - res$sigmastar * git) - exp(res$B) * exp(1/2 * res$sigmastar^2 *
      git^2 - res$mustar2 * git) * pnorm(res$b - res$sigmastar * git))/(exp(res$A) *
      pnorm(res$a) - exp(res$B) * pnorm(res$b))
    res$teBC_reciprocal <- (exp(res$A) * exp(1/2 * res$sigmastar^2 * git^2 +
      res$mustar1 * git) * pnorm(res$a + res$sigmastar * git) - exp(res$B) *
      exp(1/2 * res$sigmastar^2 * git^2 + res$mustar2 * git) * pnorm(res$b +
      res$sigmastar * git))/(exp(res$A) * pnorm(res$a) - exp(res$B) * pnorm(res$b))
  }
  res$mustar1 <- NULL
  res$mustar2 <- NULL
  res$A <- NULL
  res$B <- NULL
  res$a <- NULL
  res$b <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
