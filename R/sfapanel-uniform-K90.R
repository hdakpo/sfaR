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
# Convolution: uniform - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for uniform-normal distribution
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
puninormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar, Xvar,
  pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -S * giepsi/gisq
  sigmastar <- sqrt(exp(Wv)/gisq)
  theta <- sqrt(12) * exp(Wu/2)
  ll <- log(sigmastar) - log(theta) - TT/2 * Wv - (TT - 1)/2 * log(2 * pi) - 1/2 *
    (epsilon_isq/exp(Wv) - mustar^2/sigmastar^2) + log(pnorm((theta - mustar)/sigmastar) -
    pnorm(-mustar/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for uniform-normal distribution
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
pstuninorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, wHvar, S, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstuninorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initUni <- NULL
  } else {
    cat("Initialization: SFA + uniform-normal distribution...\n")
    initUni <- maxLik::maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgraduninormlike,
      hess = chessuninormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar,
      nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initUni$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "eta1", "eta2")
  return(list(StartVal = StartVal, initUni = initUni))
}

# Gradient of the likelihood function ----------
#' gradient for uniform-normal distribution
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
pgraduninormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
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
  sigx1 <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(sigx1/sqrt(ewv/gisq), 0, 1)
  pmusig1 <- pnorm(sigx1/sqrt(ewv/gisq))
  dmusig2 <- dnorm(S * giepsi/(sqrt(ewv/gisq) * gisq), 0, 1)
  pmusig2 <- pnorm(S * giepsi/(sqrt(ewv/gisq) * gisq))
  sigx2 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq))
  sigx3 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq) * gisq)
  sigx4 <- (0.5 * (S * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.5 *
    (sigx1 * dmusig1))
  sigx5 <- (0.5 * (sigx1 * gid2_gicub) + S * (gid1_epsi_gisq - gid2_gicub * giepsi/gisq))
  sigx6 <- (sqrt(ewv/gisq) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)))
  sigx7 <- (gid1_epsi_gisq/(sqrt(ewv/gisq) * gisq) - sigx6 * gid2_gicub * giepsi/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx8 <- (sigx5 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx7)
  sigx9 <- ((-(S * giepsi/gisq))^2 * gid2_gicub + 2 * (S^2 * (gid1_epsi_gisq -
    gid2_gicub * giepsi/gisq) * giepsi/gisq))
  sigx10 <- (gid1sq_epsi_gisq - gid2sq_gicub * giepsi/gisq)
  sigx11 <- (gid1sq_epsi_gisq/(sqrt(ewv/gisq) * gisq) - sigx6 * gid2sq_gicub *
    giepsi/(sqrt(ewv/gisq) * gisq)^2)
  sigx12 <- ((0.5 * (sigx1 * gid2sq_gicub) + S * sigx10) * dmusig1/(sqrt(ewv/gisq) *
    gisq) - S * dmusig2 * sigx11)
  sigx13 <- ((-(S * giepsi/gisq))^2 * gid2sq_gicub + 2 * (S^2 * sigx10 * giepsi/gisq))
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * (dmusig1 - dmusig2)/sigx3,
    FUN = "*") - 0.5 * (sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*")),
    sweep(uHvar, MARGIN = 1, STATS = (sqrt(12)/2 * (dmusig1 * ewu_h/sigx2) -
      0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx4/sigx2 + 0.5 +
      0.5 * ((epsilon_isq - (-(S * giepsi/gisq))^2 * gisq)/ewv) - 0.5 * TT),
      FUN = "*"), sigx8/(pmusig1 - pmusig2) + 0.5 * (sigx9/ewv) - 0.5 * (gid2_gicub/gisq),
    sigx12/(pmusig1 - pmusig2) + 0.5 * (sigx13/ewv) - 0.5 * (gid2sq_gicub/gisq))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for uniform-normal distribution
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
phessuninormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
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
  Xit_gitd1_gitsq <- sweep(Xvar, MARGIN = 1, STATS = gitd1 * git^2, FUN = "*")
  Xi_gid1_gisq <- apply(Xit_gitd1_gitsq, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xit_gitd1sq_gitsq <- sweep(Xvar, MARGIN = 1, STATS = Zit^2 * exp(Zit * (eta1 +
    eta2 * Zit)) * git^2, FUN = "*")
  Xi_gid1sq_gisq <- apply(Xit_gitd1sq_gitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
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
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sigx1 <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(sigx1/sqrt(ewv/gisq), 0, 1)
  pmusig1 <- pnorm(sigx1/sqrt(ewv/gisq))
  dmusig2 <- dnorm(S * giepsi/(sqrt(ewv/gisq) * gisq), 0, 1)
  pmusig2 <- pnorm(S * giepsi/(sqrt(ewv/gisq) * gisq))
  sigx2 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq))
  sigx3 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq) * gisq)
  sigx4 <- (0.5 * (S * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.5 *
    (sigx1 * dmusig1))
  sigx5 <- (0.5 * (sigx1 * gid2_gicub) + S * (gid1_epsi_gisq - gid2_gicub * giepsi/gisq))
  sigx6 <- (sqrt(ewv/gisq) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)))
  sigx7 <- (gid1_epsi_gisq/(sqrt(ewv/gisq) * gisq) - sigx6 * gid2_gicub * giepsi/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx8 <- (sigx5 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx7)
  sigx9 <- ((-(S * giepsi/gisq))^2 * gid2_gicub + 2 * (S^2 * (gid1_epsi_gisq -
    gid2_gicub * giepsi/gisq) * giepsi/gisq))
  sigx10 <- (gid1sq_epsi_gisq - gid2sq_gicub * giepsi/gisq)
  sigx11 <- (gid1sq_epsi_gisq/(sqrt(ewv/gisq) * gisq) - sigx6 * gid2sq_gicub *
    giepsi/(sqrt(ewv/gisq) * gisq)^2)
  sigx12 <- ((0.5 * (sigx1 * gid2sq_gicub) + S * sigx10) * dmusig1/(sqrt(ewv/gisq) *
    gisq) - S * dmusig2 * sigx11)
  sigx13 <- ((-(S * giepsi/gisq))^2 * gid2sq_gicub + 2 * (S^2 * sigx10 * giepsi/gisq))
  sigx14 <- (gid2sq_gicub * sigx10 + gd2Zfour_gd1_gicub * giepsi)
  sigx15 <- (gd1Zfour_gd1_epsi_gisq - sigx14/gisq)
  sigx16 <- (0.5 * (sigx1 * gid2sq_gicub) + S * sigx10)
  sigx17 <- (sigx6/(sqrt(ewv/gisq) * gisq)^2)
  sigx18 <- (0.5 * sigx17 - 0.5/(sqrt(ewv/gisq) * gisq^2))
  sigx19 <- (sigx18 * ewv * gid2_gicub * gid2sq_gicub + sigx6 * gd2Zcub_gd1_gicub)
  sigx20 <- (gid1_epsi_gisq - gid2_gicub * giepsi/gisq)
  sigx21 <- (gid2_gicub * sigx20 + gd2Zsq_gd1_gicub * giepsi)
  sigx22 <- (gd1Zsq_gd1_epsi_gisq - sigx21/gisq)
  sigx23 <- (gd1Zcub_gd1_epsi_gisq - (gid2sq_gicub * sigx20 + gd2Zcub_gd1_gicub *
    giepsi)/gisq)
  sigx24 <- (0.5 * (sigx1^2/ewv) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)^2))
  sigx25 <- (0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)^2)))
  sigx26 <- (sigx25/sqrt(ewv/gisq) - sigx6 * gisq/(sqrt(ewv/gisq) * gisq)^2)
  sigx27 <- ((pmusig1 - pmusig2)/(sqrt(ewv/gisq) * gisq))
  sigx28 <- ((S * dmusig2 * giepsi/gisq - sigx1 * dmusig1)/(ewv * (pmusig1 - pmusig2) *
    sqrt(ewv/gisq) * gisq) - (dmusig1 - dmusig2)^2/sigx3^2)
  sigx29 <- (sigx1/(ewv * (pmusig1 - pmusig2) * sqrt(ewv/gisq)) + (dmusig1 - dmusig2)/(sigx2^2 *
    gisq))
  sigx30 <- ((0.5 * (dmusig2 * (ewv - S^2 * giepsi^2/gisq)/(sqrt(ewv/gisq) * gisq)^2) -
    0.5 * ((1/gisq - sigx1^2/ewv) * dmusig1))/sigx2 - sigx4 * (dmusig1 - dmusig2)/(sigx2^2 *
    gisq))
  sigx31 <- (sigx6 * gid2_gicub/(sqrt(ewv/gisq) * gisq)^2 + S^2 * sigx7 * giepsi/(ewv *
    gisq))
  sigx32 <- (sigx6 * gid2sq_gicub/(sqrt(ewv/gisq) * gisq)^2 + S^2 * sigx11 * giepsi/(ewv *
    gisq))
  sigx33 <- (sigx5 * sigx1/ewv + 0.5 * (gid2_gicub/gisq))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar + nuZUvar +
    nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * sigx28, FUN = "*"), Xgi) - 0.5 * (sapply(1:nXvar, function(x) {
    crossprod(Xsq[[x]], as.matrix(wHvar/ewv))
  }) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 * (S^2/gisq)/ewv, FUN = "*"),
    Xgi))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = -wHvar * (sqrt(12)/2 * (S * sigx29 * dmusig1 * ewu_h)), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(0.5 *
    (sweep(Xepsi_i, MARGIN = 1, STATS = wHvar/ewv, FUN = "*") - sweep(Xgi, MARGIN = 1,
      STATS = wHvar * 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*")) + sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx30, FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(S * ((((sweep(Xi_gid1_gisq,
    MARGIN = 1, STATS = wHvar * dmusig1/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx33 * dmusig1/(sqrt(ewv/gisq) *
    gisq)/(pmusig1 - pmusig2), FUN = "*")) - sweep(Xgi, MARGIN = 1, STATS = wHvar *
    sigx8 * (dmusig1 - dmusig2)/(pmusig1 - pmusig2)/(sqrt(ewv/gisq) * gisq)/(pmusig1 -
    pmusig2), FUN = "*")) - (sweep(Xi_gid1_gisq, MARGIN = 1, STATS = wHvar/(sqrt(ewv/gisq) *
    gisq) * dmusig2/(pmusig1 - pmusig2), FUN = "*") - sweep(Xgi, MARGIN = 1,
    STATS = wHvar * sigx31 * dmusig2/(pmusig1 - pmusig2), FUN = "*"))) + 0.5 *
    (S * (2 * sweep(Xgi, MARGIN = 1, STATS = wHvar * (gid2_gicub * giepsi/gisq)/(ewv *
      gisq), FUN = "*") + 2 * (sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx20/(ewv *
      gisq), FUN = "*") + (sweep(Xi_gid1_gisq, MARGIN = 1, STATS = wHvar *
      giepsi/(ewv * gisq), FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar *
      gid2_gicub/gisq * giepsi/(ewv * gisq), FUN = "*")))))))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 2] <- colSums(S * (((dmusig1 * (sweep(Xi_gid1sq_gisq,
    MARGIN = 1, STATS = wHvar/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2), FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = wHvar * (sigx16 * sigx1/ewv + 0.5 * (gid2sq_gicub/gisq))/(sqrt(ewv/gisq) *
      gisq)/(pmusig1 - pmusig2), FUN = "*")) - sweep(Xgi, MARGIN = 1, STATS = wHvar *
    sigx12 * (dmusig1 - dmusig2)/(pmusig1 - pmusig2)/(sqrt(ewv/gisq) * gisq)/(pmusig1 -
    pmusig2), FUN = "*")) - (sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = wHvar/(sqrt(ewv/gisq) *
    gisq) * dmusig2/(pmusig1 - pmusig2), FUN = "*") - sweep(Xgi, MARGIN = 1,
    STATS = wHvar * sigx32 * dmusig2/(pmusig1 - pmusig2), FUN = "*"))) + 0.5 *
    (S * (sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 * (gid2sq_gicub * giepsi/gisq)/(ewv *
      gisq), FUN = "*") + 2 * (sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx10/(ewv *
      gisq), FUN = "*") + (sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = wHvar *
      giepsi/(ewv * gisq), FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar *
      gid2sq_gicub/gisq * giepsi/(ewv * gisq), FUN = "*")))))))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sqrt(12)/2 * (((0.5 - sqrt(12)/2 * (sigx1 * ewu_h *
      gisq/ewv))/sigx2 - sqrt(12)/2 * (dmusig1 * ewu_h/sigx2^2)) * dmusig1 *
      ewu_h), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((0.5 *
    ((sqrt(12)/2 - sqrt(12)/2 * (sigx1^2 * gisq/ewv))/sigx2) + sqrt(12)/2 * (sigx4/sigx2^2)) *
    dmusig1 * ewu_h), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/4 * gid2_gicub - sqrt(12)/2 * (sigx5 *
      sigx1 * gisq/ewv))/gisq - sqrt(12)/2 * (sigx8/(pmusig1 - pmusig2))) *
      dmusig1 * ewu_h/sigx2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 2] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/4 * gid2sq_gicub - sqrt(12)/2 * (sigx16 *
      sigx1 * gisq/ewv))/gisq - sqrt(12)/2 * (sigx12/(pmusig1 - pmusig2))) *
      dmusig1 * ewu_h/sigx2, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (S * ((0.5 * (S^2 * giepsi^2) - ewv * gisq)/(sqrt(ewv/gisq) * gisq)^2 +
      1) * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.25 * (sigx1^3 *
      dmusig1 * gisq/ewv))/sigx2 - (((0.5 * sigx27 + 0.5 * (S * dmusig2 * giepsi/(sqrt(ewv/gisq) *
      gisq)^2)) * ewv - 0.5 * (sigx1 * dmusig1)) * sigx4/sigx2^2 + 0.5 * ((epsilon_isq -
      (-(S * giepsi/gisq))^2 * gisq)/ewv))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (((sigx5 * sigx24 *
    dmusig1 - sigx8 * sigx4/(pmusig1 - pmusig2))/sqrt(ewv/gisq) - S * (0.5 *
    (S^2 * sigx7 * giepsi^2) - (sigx26 * gid2_gicub * giepsi + 0.5 * (gid1_epsi_gisq/sqrt(ewv/gisq))) *
    ewv) * dmusig2/(sqrt(ewv/gisq) * gisq)^2)/(pmusig1 - pmusig2) - 0.5 * (sigx9/ewv)),
    FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    2] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (((sigx16 * sigx24 *
    dmusig1 - sigx12 * sigx4/(pmusig1 - pmusig2))/sqrt(ewv/gisq) - S * (0.5 *
    (S^2 * sigx11 * giepsi^2) - (sigx26 * gid2sq_gicub * giepsi + 0.5 * (gid1sq_epsi_gisq/sqrt(ewv/gisq))) *
    ewv) * dmusig2/(sqrt(ewv/gisq) * gisq)^2)/(pmusig1 - pmusig2) - 0.5 * (sigx13/ewv)),
    FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 1] <- sum(wHvar *
    ((((0.5 * (sigx1 * gd2Zsq_gd1_gicub + S * gid2_gicub * sigx20/gisq) + S *
      sigx22 - sigx5^2 * sigx1/ewv)/(sqrt(ewv/gisq) * gisq) - sigx5 * sigx6 *
      gid2_gicub/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 - (sigx8^2/(pmusig1 -
      pmusig2) + S * ((gd1Zsq_gd1_epsi_gisq - S^2 * sigx7^2 * giepsi)/(sqrt(ewv/gisq) *
      gisq) - ((sigx18 * ewv * gid2_gicub^2 + sigx6 * gd2Zsq_gd1_gicub) * giepsi +
      (2 * gid1_epsi_gisq - 2 * (sqrt(ewv/gisq) * sigx6 * gid2_gicub * gisq *
        giepsi/(sqrt(ewv/gisq) * gisq)^2)) * sigx6 * gid2_gicub)/(sqrt(ewv/gisq) *
      gisq)^2) * dmusig2))/(pmusig1 - pmusig2) + 0.5 * (((-(S * giepsi/gisq))^2 *
      gd2Zsq_gd1_gicub + S^2 * (2 * (sigx20^2 + sigx22 * giepsi) + 2 * (gid2_gicub *
      sigx20 * giepsi/gisq))/gisq)/ewv) - 0.5 * ((gd2Zsq_gd1_gicub - gid2_gicub^2/gisq)/gisq)))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 2] <- sum(wHvar *
    ((((0.5 * (sigx1 * gd2Zcub_gd1_gicub + S * gid2sq_gicub * sigx20/gisq) +
      S * sigx23 - sigx5 * sigx16 * sigx1/ewv)/(sqrt(ewv/gisq) * gisq) - sigx16 *
      sigx6 * gid2_gicub/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 - (sigx8 * sigx12/(pmusig1 -
      pmusig2) + S * ((gd1Zcub_gd1_epsi_gisq - S^2 * sigx7 * sigx11 * giepsi)/(sqrt(ewv/gisq) *
      gisq) - (sigx19 * giepsi + sigx6 * (gid2_gicub * gid1sq_epsi_gisq + gid2sq_gicub *
      (gid1_epsi_gisq - 2 * (sqrt(ewv/gisq) * sigx6 * gid2_gicub * gisq * giepsi/(sqrt(ewv/gisq) *
        gisq)^2))))/(sqrt(ewv/gisq) * gisq)^2) * dmusig2))/(pmusig1 - pmusig2) +
      0.5 * (((-(S * giepsi/gisq))^2 * gd2Zcub_gd1_gicub + S^2 * (2 * (sigx20 *
        sigx10 + sigx23 * giepsi) + 2 * (gid2sq_gicub * sigx20 * giepsi/gisq))/gisq)/ewv) -
      0.5 * ((gd2Zcub_gd1_gicub - gid2_gicub * gid2sq_gicub/gisq)/gisq)))
  hessll[nXvar + nuZUvar + nvZVvar + 2, nXvar + nuZUvar + nvZVvar + 2] <- sum(wHvar *
    ((((0.5 * (sigx1 * gd2Zfour_gd1_gicub + S * gid2sq_gicub * sigx10/gisq) +
      S * sigx15 - sigx16^2 * sigx1/ewv)/(sqrt(ewv/gisq) * gisq) - sigx16 *
      sigx6 * gid2sq_gicub/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 - (sigx12^2/(pmusig1 -
      pmusig2) + S * ((gd1Zfour_gd1_epsi_gisq - S^2 * sigx11^2 * giepsi)/(sqrt(ewv/gisq) *
      gisq) - ((sigx18 * ewv * gid2sq_gicub^2 + sigx6 * gd2Zfour_gd1_gicub) *
      giepsi + (2 * gid1sq_epsi_gisq - 2 * (sqrt(ewv/gisq) * sigx6 * gid2sq_gicub *
      gisq * giepsi/(sqrt(ewv/gisq) * gisq)^2)) * sigx6 * gid2sq_gicub)/(sqrt(ewv/gisq) *
      gisq)^2) * dmusig2))/(pmusig1 - pmusig2) + 0.5 * (((-(S * giepsi/gisq))^2 *
      gd2Zfour_gd1_gicub + S^2 * (2 * (sigx10^2 + sigx15 * giepsi) + 2 * (gid2sq_gicub *
      sigx10 * giepsi/gisq))/gisq)/ewv) - 0.5 * ((gd2Zfour_gd1_gicub - gid2sq_gicub^2/gisq)/gisq)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for uniform-normal distribution
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
uninormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstuninorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    Inituni <- start_st$inituni
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(puninormlike_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(puninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgraduninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = puninormlike_k90,
    grad = pgraduninormlike_k90, hess = phessuninormlike_k90, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(puninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(puninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessuninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(puninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessuninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(puninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgraduninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessuninormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgraduninormlike_k90(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessuninormlike_k90(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessuninormlike_k90(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- puninormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgraduninormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, Inituni = Inituni))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for uniform-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
puninormeff_k90 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 2]
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
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  theta <- sqrt(12 * exp(Wu))
  mustar <- -object$S * giepsi/gisq
  sigmastar <- sqrt(exp(Wv)/gisq)
  u1 <- -sigmastar * ((dnorm((theta - mustar)/sigmastar) - dnorm(-mustar/sigmastar))/(pnorm((theta -
    mustar)/sigmastar) - pnorm(-mustar/sigmastar))) + mustar
  u2 <- -sigmastar * (dnorm(-mustar/sigmastar)/(1 - pnorm(-mustar/sigmastar)) +
    mustar/sigmastar)  # when theta/sigmav ---> Infty
  uLB <- sigmastar * qnorm((1 - level)/2 * pnorm((theta - mustar)/sigmastar) +
    (1 - (1 - level)/2) * pnorm(-mustar/sigmastar)) + mustar
  uUB <- sigmastar * qnorm((1 - (1 - level)/2) * pnorm((theta - mustar)/sigmastar) +
    (1 - level)/2 * pnorm(-mustar/sigmastar)) + mustar
  m <- ifelse(-theta < mustar & mustar < 0, -mustar, ifelse(mustar >= 0, 0, theta))
  res <- data.frame(levels(pindex[, 1]), u1 = u1, u2 = u2, uLB = uLB, uUB = uUB,
    m = m, mustar = mustar, sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u1 <- res$u1 * git
  res$u2 <- res$u2 * git
  res$m <- res$m * git
  res$uLB <- res$uLB * git
  res$uUB <- res$uUB * git
  if (object$logDepVar == TRUE) {
    res$teJLMS1 <- exp(-res$u1)
    res$teJLMS2 <- exp(-res$u2)
    res$teMO <- exp(-res$m)
    res$teBC1 <- exp(-mustar * git + sigmastar^2 * git^2/2) * (pnorm((-mustar +
      theta)/sigmastar + sigmastar * git) - pnorm(-mustar/sigmastar + sigmastar *
      git))/(pnorm((theta - mustar)/sigmastar) - pnorm(-mustar/sigmastar))
    res$teBC2 <- exp(-mustar * git + sigmastar^2 * git^2/2) * (1 - pnorm(-mustar/sigmastar +
      sigmastar * git))/(1 - pnorm(-mustar/sigmastar))
    resteBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    res$teBC1_reciprocal <- exp(-mustar * git + sigmastar^2 * git^2/2) * (pnorm((-mustar +
      theta)/sigmastar - sigmastar * git) - pnorm(-mustar/sigmastar - sigmastar *
      git))/(pnorm((theta - mustar)/sigmastar) - pnorm(-mustar/sigmastar))
    res$teBC2_reciprocal <- exp(-mustar * git + sigmastar^2 * git^2/2) * (1 -
      pnorm(-mustar/sigmastar - sigmastar * git))/(1 - pnorm(-mustar/sigmastar))
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
