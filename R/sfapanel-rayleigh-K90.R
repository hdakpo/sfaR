################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Two types: - Kumbakhar 1990 (K90)                                            #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = (1 + exp(eta1 * t + eta2 * t^2))^(-1)                  #
# Convolution: rayleigh - normal                                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for rayleigh-normal distribution
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
praynormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 *
    seq(1:x) + eta2 * (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * S * giepsi/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- 1/2 * (mustar/sigmastar)^2 - epsilon_isq/(2 * exp(Wv)) -
    (TT - 1)/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 * Wu -
    1/2 * log(exp(Wv) + gisq * exp(Wu)) + log(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for rayleigh-normal distribution
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
pstraynorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, printInfo, tol, whichStart,
  initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstraynorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initRay <- NULL
  } else {
    cat("Initialization: SFA + rayleigh-normal distribution...\n")
    initRay <- maxLik::maxLik(logLik = craynormlike, start = cstraynorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradraynormlike, method = initAlg,
      control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar)
    Esti <- initRay$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "eta1",
    "eta2")
  return(list(StartVal = StartVal, initRay = initRay))
}

# Gradient of the likelihood function ----------
#' gradient for rayleigh-normal distribution
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
pgradraynormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 *
    seq(1:x) + eta2 * (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1 <- Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd1_epsit_gitsq <- -gitd1 * epsilon_it * git^2
  gid1_epsi_gisq <- as.numeric(tapply(gitd1_epsit_gitsq, pindex[,
    1], sum))
  gitd2 <- 2 * Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd2_gitcub <- -gitd2 * git^3
  gid2_gicub <- as.numeric(tapply(gitd2_gitcub, pindex[, 1],
    sum))
  gitd1sq_epsit_gitsq <- -Zit^2 * exp(Zit * (eta1 + eta2 *
    Zit)) * epsilon_it * git^2
  gid1sq_epsi_gisq <- as.numeric(tapply(gitd1sq_epsit_gitsq,
    pindex[, 1], sum))
  gitd2sq_gitcub <- -2 * Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) *
    git^3
  gid2sq_gicub <- as.numeric(tapply(gitd2sq_gitcub, pindex[,
    1], sum))
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (S * ewu * giepsi/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig, 0, 1)
  sigx1 <- (sigmastar/ewv - ewu/(ssqx2))
  sigx2 <- (pmusig + S * dmusig * sigx1 * giepsi)
  sigx3 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/sigmasq)
  sigx4 <- (S * giepsi/ewv - sigx2/sigx3)
  sigx5 <- (0.5 * (dmusig * ewv/sigmastar) - S * pmusig * giepsi)
  sigx6 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx7 <- (1/(ssqx2) - (0.5 * sigx6 + sigmastar * gisq) *
    ewu/ssq)
  sigx8 <- (sigx5 * (1 - ewu * gisq/sigmasq)/sigx3 + S^2 *
    sigx7 * ewu * giepsi^2/sigmastar - 0.5 * gisq)
  sigx9 <- ((1 - ewv/sigmasq) * dmusig/sigmastar)
  sigx10 <- (0.5 * sigx9 + S * pmusig * giepsi/sigmasq)
  sigx11 <- (0.5 * ((1 - ewv/sigmasq) * ewu/sigmastar) + sigmastar)
  sigx12 <- ((sigx10/sigx3 - S^2 * sigx11 * ewu * giepsi^2/(ssq *
    sigmastar)) * ewu - 0.5)
  sigx13 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  sigx14 <- (gid1_epsi_gisq/(ssqx2) - ewu * sigx13 * gid2_gicub *
    giepsi/ssq)
  sigx15 <- (0.5 * (dmusig * ewv * gid2_gicub/(ssqx2)) + S^2 *
    dmusig * sigx14 * giepsi)
  sigx16 <- (pmusig * gid1_epsi_gisq - ewu * (pmusig * gid2_gicub/sigmasq +
    S * dmusig * sigx14) * giepsi)
  sigx17 <- ((sigx15 * ewu + S * sigx16)/sigx3 + 0.5 * gid2_gicub)
  sigx18 <- (gid1sq_epsi_gisq/(ssqx2) - ewu * sigx13 * gid2sq_gicub *
    giepsi/ssq)
  sigx19 <- (0.5 * (dmusig * ewv * gid2sq_gicub/(ssqx2)) +
    S^2 * dmusig * sigx18 * giepsi)
  sigx20 <- (sigx19 * ewu + S * (pmusig * gid1sq_epsi_gisq -
    ewu * (pmusig * gid2sq_gicub/sigmasq + S * dmusig * sigx18) *
      giepsi))
  sigx21 <- (sigx20/sigx3 + 0.5 * gid2sq_gicub)
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * ewu *
    sigx4/sigmasq, FUN = "*") - sweep(Xepsi_i, MARGIN = 1,
    STATS = 1/(2 * ewv), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = (sigx8 * ewu/sigmasq - 0.5), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = ((sigx12/sigmasq + 2 * (epsilon_isq/(2 *
      ewv)^2)) * ewv - 0.5 * (TT - 1)), FUN = "*"), ewu *
    (S^2 * ewu * sigx14 * giepsi/sigmastar - sigx17)/sigmasq,
    ewu * (S^2 * ewu * sigx18 * giepsi/sigmastar - sigx21)/sigmasq)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for rayleigh-normal distribution
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
phessraynormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nuZUvar + nvZVvar + 2]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 *
    seq(1:x) + eta2 * (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zit <- unlist(lapply(TT, FUN = function(x) seq(1:x)))
  gitd1 <- Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd1_epsit_gitsq <- -gitd1 * epsilon_it * git^2
  gid1_epsi_gisq <- as.numeric(tapply(gitd1_epsit_gitsq, pindex[,
    1], sum))
  gitd2 <- 2 * Zit * exp(Zit * (eta1 + eta2 * Zit))
  gitd2_gitcub <- -gitd2 * git^3
  gid2_gicub <- as.numeric(tapply(gitd2_gitcub, pindex[, 1],
    sum))
  gitd1sq_epsit_gitsq <- -Zit^2 * exp(Zit * (eta1 + eta2 *
    Zit)) * epsilon_it * git^2
  gid1sq_epsi_gisq <- as.numeric(tapply(gitd1sq_epsit_gitsq,
    pindex[, 1], sum))
  gitd2sq_gitcub <- -2 * Zit^2 * exp(Zit * (eta1 + eta2 * Zit)) *
    git^3
  gid2sq_gicub <- as.numeric(tapply(gitd2sq_gitcub, pindex[,
    1], sum))
  gd1Zsqt_gd1t_epsit_gitsq <- -Zit^2 * (1 - 2 * (exp(Zit *
    (eta1 + eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 *
    Zit)) * epsilon_it * git^2
  gd1Zsq_gd1_epsi_gisq <- as.numeric(tapply(gd1Zsqt_gd1t_epsit_gitsq,
    pindex[, 1], sum))
  gd1Zcubt_gd1t_epsit_gitsq <- -Zit^3 * (1 - 2 * (exp(Zit *
    (eta1 + eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 *
    Zit)) * epsilon_it * git^2
  gd1Zcub_gd1_epsi_gisq <- as.numeric(tapply(gd1Zcubt_gd1t_epsit_gitsq,
    pindex[, 1], sum))
  gd1Zfourt_gd1t_epsit_gitsq <- -Zit^4 * (1 - 2 * (exp(Zit *
    (eta1 + eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 *
    Zit)) * epsilon_it * git^2
  gd1Zfour_gd1_epsi_gisq <- as.numeric(tapply(gd1Zfourt_gd1t_epsit_gitsq,
    pindex[, 1], sum))
  gd2Zsqt_gd1t_gitcub <- -2 * Zit^2 * (1 - 3 * (exp(Zit * (eta1 +
    eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 * Zit)) *
    git^3
  gd2Zsq_gd1_gicub <- as.numeric(tapply(gd2Zsqt_gd1t_gitcub,
    pindex[, 1], sum))
  gd2Zcubt_gd1t_gitcub <- -2 * Zit^3 * (1 - 3 * (exp(Zit *
    (eta1 + eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 *
    Zit)) * git^3
  gd2Zcub_gd1_gicub <- as.numeric(tapply(gd2Zcubt_gd1t_gitcub,
    pindex[, 1], sum))
  gd2Zfourt_gd1t_gitcub <- -2 * Zit^4 * (1 - 3 * (exp(Zit *
    (eta1 + eta2 * Zit)) * git)) * exp(Zit * (eta1 + eta2 *
    Zit)) * git^3
  gd2Zfour_gd1_gicub <- as.numeric(tapply(gd2Zfourt_gd1t_gitcub,
    pindex[, 1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  Xit_gitd1_gitsq <- sweep(Xvar, MARGIN = 1, STATS = gitd1 *
    git^2, FUN = "*")
  Xi_gid1_gisq <- apply(Xit_gitd1_gitsq, 2, function(x) tapply(x,
    pindex[, 1], sum))
  Xit_gitd1sq_gitsq <- sweep(Xvar, MARGIN = 1, STATS = Zit^2 *
    exp(Zit * (eta1 + eta2 * Zit)) * git^2, FUN = "*")
  Xi_gid1sq_gisq <- apply(Xit_gitd1sq_gitsq, 2, function(x) tapply(x,
    pindex[, 1], sum))
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (S * ewu * giepsi/(ssqx2))
  pmusig <- pnorm(-musig)
  dmusig <- dnorm(-musig, 0, 1)
  sigx1 <- (sigmastar/ewv - ewu/(ssqx2))
  sigx2 <- (pmusig + S * dmusig * sigx1 * giepsi)
  sigx3 <- (dmusig * sigmastar - S * ewu * pmusig * giepsi/sigmasq)
  sigx4 <- (S * giepsi/ewv - sigx2/sigx3)
  sigx5 <- (0.5 * (dmusig * ewv/sigmastar) - S * pmusig * giepsi)
  sigx6 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx7 <- (1/(ssqx2) - (0.5 * sigx6 + sigmastar * gisq) *
    ewu/ssq)
  sigx8 <- (sigx5 * (1 - ewu * gisq/sigmasq)/sigx3 + S^2 *
    sigx7 * ewu * giepsi^2/sigmastar - 0.5 * gisq)
  sigx9 <- ((1 - ewv/sigmasq) * dmusig/sigmastar)
  sigx10 <- (0.5 * sigx9 + S * pmusig * giepsi/sigmasq)
  sigx11 <- (0.5 * ((1 - ewv/sigmasq) * ewu/sigmastar) + sigmastar)
  sigx12 <- ((sigx10/sigx3 - S^2 * sigx11 * ewu * giepsi^2/(ssq *
    sigmastar)) * ewu - 0.5)
  sigx13 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  sigx14 <- (gid1_epsi_gisq/(ssqx2) - ewu * sigx13 * gid2_gicub *
    giepsi/ssq)
  sigx15 <- (0.5 * (dmusig * ewv * gid2_gicub/(ssqx2)) + S^2 *
    dmusig * sigx14 * giepsi)
  sigx16 <- (pmusig * gid1_epsi_gisq - ewu * (pmusig * gid2_gicub/sigmasq +
    S * dmusig * sigx14) * giepsi)
  sigx17 <- ((sigx15 * ewu + S * sigx16)/sigx3 + 0.5 * gid2_gicub)
  sigx18 <- (gid1sq_epsi_gisq/(ssqx2) - ewu * sigx13 * gid2sq_gicub *
    giepsi/ssq)
  sigx19 <- (0.5 * (dmusig * ewv * gid2sq_gicub/(ssqx2)) +
    S^2 * dmusig * sigx18 * giepsi)
  sigx20 <- (sigx19 * ewu + S * (pmusig * gid1sq_epsi_gisq -
    ewu * (pmusig * gid2sq_gicub/sigmasq + S * dmusig * sigx18) *
      giepsi))
  sigx21 <- (sigx20/sigx3 + 0.5 * gid2sq_gicub)
  sigx22 <- (1 - ewu * gisq/sigmasq)
  sigx23 <- ((0.5 * sigx22 - 1.5)/(ssqx2) - 0.5 * sigx7)
  ssg6 <- (0.5 * sigx6 + sigmastar * gisq)
  sigx24 <- (sigx23 * ewv - 2 * (ssg6 * ssqx2 * sigx13/ssq))
  sigx25 <- (sigx24 * ewu + 3 * sigmastar)
  sigx26 <- (ssg6 * dmusig * ewv/ssq + S^2 * sigx7 * dmusig *
    giepsi^2/sigmasq)
  sigx27 <- ((((0.5 * sigx22 - 0.5)/(ssqx2) - 0.5 * sigx7) *
    ewv - 2 * (ssg6 * ssqx2 * sigx13/ssq)) * ewu + sigmastar)
  sigx28 <- (S^2 * ewu * sigx18 * giepsi/sigmastar - sigx21)
  sigx29 <- (S^2 * ewu * sigx14 * giepsi/sigmastar - sigx17)
  sigx30 <- ((dmusig + S^2 * sigx11 * dmusig * ewu^2 * ewv *
    giepsi^2/(ssq * ssqx2))/(ssqx2) - sigx11 * dmusig * ewv/ssq)
  sigx31 <- (S * sigx11 * dmusig * ewu * giepsi/ssq - pmusig/sigmasq)
  sigx32 <- ((0.5 * ((1 - ewv/sigmasq)/(ssqx2)) - 0.5 * (1/(ssqx2) -
    sigx11 * ewv/ssq)) * ewu - 2 * (sigx11 * ssqx2 * sigx13/ssq))
  sigx33 <- (0.5 * (sigx13/ssq) - 0.5/(sigmasq^2 * sigmastar))
  sigx34 <- (sigx33 * ewu^2 * ewv * gid2_gicub^2 + sigx13 *
    gd2Zsq_gd1_gicub)
  sigx35 <- (sigx34 * giepsi + (2 * gid1_epsi_gisq - 2 * (ewu *
    ssqx2 * sigx13 * gid2_gicub * giepsi/ssq)) * sigx13 *
    gid2_gicub)
  sigx36 <- (sigx35/ssq + S^2 * ewu * sigx14^2 * giepsi/(ssqx2))
  sigx37 <- (sigx33 * ewu^2 * ewv * gid2sq_gicub^2 + sigx13 *
    gd2Zfour_gd1_gicub)
  sigx38 <- (sigx37 * giepsi + (2 * gid1sq_epsi_gisq - 2 *
    (ewu * ssqx2 * sigx13 * gid2sq_gicub * giepsi/ssq)) *
    sigx13 * gid2sq_gicub)
  sigx39 <- (gd1Zfour_gd1_epsi_gisq/(ssqx2) - (sigx38/ssq +
    S^2 * ewu * sigx18^2 * giepsi/(ssqx2)) * ewu)
  sigx40 <- (sigx33 * ewu^2 * ewv * gid2_gicub * gid2sq_gicub +
    sigx13 * gd2Zcub_gd1_gicub)
  sigx41 <- (sigx40 * giepsi + sigx13 * (gid2_gicub * gid1sq_epsi_gisq +
    gid2sq_gicub * (gid1_epsi_gisq - 2 * (ewu * ssqx2 * sigx13 *
      gid2_gicub * giepsi/ssq))))
  sigx42 <- (gd1Zcub_gd1_epsi_gisq/(ssqx2) - (sigx41/ssq +
    S^2 * ewu * sigx14 * sigx18 * giepsi/(ssqx2)) * ewu)
  sigx43 <- (gd1Zsq_gd1_epsi_gisq/(ssqx2) - sigx36 * ewu)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar +
    nuZUvar + nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S^2 * (1/ewv - (((1 - S^2 * ewu * giepsi^2/(sigmasq *
      ewv)) * sigx1 - ewu/(ssqx2)) * dmusig + ewu * sigx2^2/(sigx3 *
      sigmasq))/sigx3) * ewu/sigmasq, FUN = "*"), Xgi) -
    sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar/(2 *
      ewv))))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (((sigx5 * sigx2/sigx3 +
      0.5 * (S * dmusig * giepsi/sigmastar)) * ewu/sigmasq -
      pmusig) * sigx22/sigx3 + 2 * (S * sigx7 * ewu * giepsi/sigmastar)) *
      ewu/sigmasq, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = wHvar *
    2/(2 * ewv)^2 * ewv, FUN = "*") + sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * (((sigx10 * sigx2/sigx3 - S * (0.5 *
      ((1 - ewv/sigmasq)/ewv) + 1/sigmasq) * dmusig * giepsi/sigmastar) *
      ewu + pmusig)/(sigx3 * sigmasq) - 2 * (S * sigx11 *
      ewu * giepsi/(ssq * sigmastar))) * ewu/sigmasq *
      ewv, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx14 * S * ewu/sigmastar *
      S * ewu/sigmasq, FUN = "*") + sweep(Xi_gid1_gisq,
    MARGIN = 1, STATS = wHvar * giepsi/(ssqx2) * S * ewu/sigmastar *
      S * ewu/sigmasq, FUN = "*") - sweep(Xgi, MARGIN = 1,
    STATS = wHvar * ewu * sigx13 * gid2_gicub/ssq * giepsi *
      S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = wHvar * ((sigx15 * ewu +
      S * sigx16) * sigx2/sigx3 + S * (0.5 * (ewu * gid2_gicub *
      giepsi/sigmasq) - gid1_epsi_gisq) * dmusig/sigmastar -
      pmusig * gid2_gicub) * ewu/sigmasq/sigx3 * S * ewu/sigmasq,
      FUN = "*") - sweep(Xi_gid1_gisq, MARGIN = 1, STATS = wHvar *
    pmusig/sigx3 * S * ewu/sigmasq, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx18 * S * ewu/sigmastar *
      S * ewu/sigmasq, FUN = "*") + sweep(Xi_gid1sq_gisq,
    MARGIN = 1, STATS = wHvar * giepsi/(ssqx2) * S * ewu/sigmastar *
      S * ewu/sigmasq, FUN = "*") - sweep(Xgi, MARGIN = 1,
    STATS = wHvar * ewu * sigx13 * gid2sq_gicub/ssq * giepsi *
      S * ewu/sigmastar * S * ewu/sigmasq, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = wHvar * (sigx20 * sigx2/sigx3 +
      S * (0.5 * (ewu * gid2sq_gicub * giepsi/sigmasq) -
        gid1sq_epsi_gisq) * dmusig/sigmastar - pmusig *
      gid2sq_gicub) * ewu/sigmasq/sigx3 * S * ewu/sigmasq,
      FUN = "*") - sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = wHvar *
    pmusig/sigx3 * S * ewu/sigmasq, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (dmusig * ewv/sigmastar) + ewu * (S^2 * sigx7 *
      dmusig * giepsi^2 - (sigx5 * sigx22/sigx3 + gisq) *
      sigx5/sigmasq) - (0.5 * (0.5 * (sigx22 * dmusig *
      ewv/sigmastar) + S^2 * sigx7 * dmusig * ewu * giepsi^2) +
      S * pmusig * giepsi)) * sigx22/sigx3 + ewu * (S^2 *
      (2/(ssqx2) - (((0.5 * (ewu * gisq/sigmasq) + 2 -
        0.5 * (0.5 * sigx22 + ewu * gisq/sigmasq)) *
        sigx22 * ewv/sigmastar + (4 * gisq - 2 * (ssg6^2 *
        ewu * sigmasq/ssq)) * sigmastar) * ewu/ssq +
        0.5 * (sigx22 * sigx7))) * giepsi^2/sigmastar -
      sigx8 * gisq/sigmasq) - 0.5 * gisq) * ewu/sigmasq,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * ((1 - ewv/sigmasq) *
      dmusig) + 0.5 * ((dmusig * ewv * gisq/sigmasq - S^2 *
      (1 - ewv/sigmasq) * sigx7 * dmusig * ewu * giepsi^2/sigmastar) *
      ewu/sigmasq - 0.5 * (sigx22 * (1 - ewv/sigmasq) *
      dmusig)))/sigmastar + (S * pmusig * giepsi - (sigx10 *
      sigx5 * sigx22/sigx3 + S * (pmusig * gisq/sigmasq +
      S * sigx7 * dmusig * giepsi) * giepsi) * ewu)/sigmasq)/sigx3 -
      (sigx12 * gisq/sigmasq + S^2 * (((0.5 * (sigx22 *
        ewv/sigmasq) + 0.5 * ((ewu * gisq/sigmasq - 1) *
        ewv/sigmasq + 1 - 0.5 * (sigx22 * (1 - ewv/sigmasq))) +
        0.5 * (1 - ewv/sigmasq)) * ewu/sigmastar + sigmastar)/(ssq *
        sigmastar) + sigx11 * (1/(ssq * sigmastar) -
        (0.5 * (ssq * sigx22/(ssqx2)) + 2 * (ssg6 * ewu)) *
          ewu * ewv/(ssq * sigmastar)^2)) * ewu * giepsi^2)) *
      ewu * ewv/sigmasq, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ewu * (ewu * (S^2 * (2 * (gid1_epsi_gisq/(ssqx2)) - ((sigx25 *
    gid2_gicub * giepsi + ssg6 * gid1_epsi_gisq) * ewu/ssq +
    0.5 * (sigx22 * sigx14))) * giepsi/sigmastar - (((0.5 *
    (dmusig * ewv * gid2_gicub/sigmastar) - (sigx15 * ewu +
    S * sigx16) * sigx5 * sigx22/sigx3)/sigmasq + S * (((pmusig *
    gisq/sigmasq + S * sigx7 * dmusig * giepsi) * gid2_gicub/sigmasq +
    S * ((sigx27 * gid2_gicub * giepsi + ssg6 * gid1_epsi_gisq)/ssq +
      S^2 * sigx7 * ewu * sigx14 * giepsi^2/(ssqx2)) *
      dmusig) * ewu + S * (ssg6 * ewu/ssq - 1/(ssqx2)) *
    dmusig * gid1_epsi_gisq - pmusig * gid2_gicub/sigmasq) *
    giepsi - (0.5 * (sigx26 * gid2_gicub) + S^2 * ((sigx27 *
    gid2_gicub * giepsi + ssg6 * gid1_epsi_gisq)/ssq + S^2 *
    sigx7 * ewu * sigx14 * giepsi^2/(ssqx2)) * dmusig * giepsi) *
    ewu)/sigx3 + sigx29 * gisq/sigmasq)) - sigx17)/sigmasq,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ewu * (ewu * (S^2 * (2 * (gid1sq_epsi_gisq/(ssqx2)) -
    ((sigx25 * gid2sq_gicub * giepsi + ssg6 * gid1sq_epsi_gisq) *
      ewu/ssq + 0.5 * (sigx22 * sigx18))) * giepsi/sigmastar -
    (((0.5 * (dmusig * ewv * gid2sq_gicub/sigmastar) - sigx20 *
      sigx5 * sigx22/sigx3)/sigmasq + S * (((pmusig * gisq/sigmasq +
      S * sigx7 * dmusig * giepsi) * gid2sq_gicub/sigmasq +
      S * ((sigx27 * gid2sq_gicub * giepsi + ssg6 * gid1sq_epsi_gisq)/ssq +
        S^2 * sigx7 * ewu * sigx18 * giepsi^2/(ssqx2)) *
        dmusig) * ewu + S * (ssg6 * ewu/ssq - 1/(ssqx2)) *
      dmusig * gid1sq_epsi_gisq - pmusig * gid2sq_gicub/sigmasq) *
      giepsi - (0.5 * (sigx26 * gid2sq_gicub) + S^2 * ((sigx27 *
      gid2sq_gicub * giepsi + ssg6 * gid1sq_epsi_gisq)/ssq +
      S^2 * sigx7 * ewu * sigx18 * giepsi^2/(ssqx2)) *
      dmusig * giepsi) * ewu)/sigx3 + sigx28 * gisq/sigmasq)) -
    sigx21)/sigmasq, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (ewv * (S^2 *
      sigx11 * dmusig * ewu^2 * giepsi^2/(ssq * sigmastar) -
      dmusig)/sigmasq - 0.5 * ((1 - ewv/sigmasq) * dmusig)) +
      0.5 * dmusig) * (1 - ewv/sigmasq)/sigmastar + (ewv *
      (S * sigx31 * giepsi - sigx10^2 * ewu/sigx3) + S *
      pmusig * giepsi)/sigmasq)/sigx3 - S^2 * (sigx11 *
      (1/(ssq * sigmastar) - (0.5 * (ssq * (1 - ewv/sigmasq)/(ssqx2)) +
        2 * (sigx11 * ewv)) * ewu * ewv/(ssq * sigmastar)^2) +
      (0.5 * (ewv/sigmasq) - 0.5 * (0.5 * (1 - ewv/sigmasq) +
        ewv/sigmasq)) * (1 - ewv/sigmasq) * sigmasq/(ssq *
        ewv)) * ewu * giepsi^2) * ewu - (sigx12 * ewv/sigmasq +
      0.5))/sigmasq + (2 - 16 * (ewv^2/(2 * ewv)^2)) *
      epsilon_isq/(2 * ewv)^2) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((0.5 * (sigx30 * gid2_gicub) +
      S * (S * sigx11 * dmusig * gid1_epsi_gisq/ssq - sigx31 *
        gid2_gicub/sigmasq) * giepsi - (sigx15 * ewu +
      S * sigx16) * sigx10/(sigx3 * sigmasq)) * ewv/sigx3 +
      S^2 * ((sigx32 * ewu * gid2_gicub * giepsi + sigx11 *
        gid1_epsi_gisq) * ewv/ssq + 0.5 * ((1 - ewv/sigmasq) *
        sigx14)) * giepsi/sigmastar) * ewu + ewv * sigx29/sigmasq) *
      ewu/sigmasq), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((0.5 * (sigx30 * gid2sq_gicub) +
      S * (S * sigx11 * dmusig * gid1sq_epsi_gisq/ssq -
        sigx31 * gid2sq_gicub/sigmasq) * giepsi - sigx20 *
      sigx10/(sigx3 * sigmasq)) * ewv/sigx3 + S^2 * ((sigx32 *
      ewu * gid2sq_gicub * giepsi + sigx11 * gid1sq_epsi_gisq) *
      ewv/ssq + 0.5 * ((1 - ewv/sigmasq) * sigx18)) * giepsi/sigmastar) *
      ewu + ewv * sigx28/sigmasq) * ewu/sigmasq), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar +
    nvZVvar + 1)] <- sum(wHvar * (ewu * (ewu * (S^2 * ((0.5 *
    (ewu * gid2_gicub * giepsi/sigmasq) + gid1_epsi_gisq) *
    sigx14 + (gd1Zsq_gd1_epsi_gisq/(ssqx2) - sigx35 * ewu/ssq) *
    giepsi)/sigmastar - sigx29 * gid2_gicub/sigmasq) - ((((sigx15 *
    ewu + S * sigx16)^2/(sigx3 * sigmasq) + 0.5 * (((dmusig *
    gd2Zsq_gd1_gicub - S^2 * dmusig * ewu^2 * gid2_gicub *
    sigx14 * giepsi/(ssqx2))/(ssqx2) - dmusig * ewu * sigx13 *
    gid2_gicub^2/ssq) * ewv) + S^2 * (sigx43 * giepsi + gid1_epsi_gisq *
    sigx14) * dmusig) * ewu + S * (pmusig * gd1Zsq_gd1_epsi_gisq -
    (((pmusig * gd2Zsq_gd1_gicub - ewu * (pmusig * gid2_gicub/sigmasq +
      S * dmusig * sigx14) * gid2_gicub)/sigmasq + S *
      dmusig * sigx43) * giepsi + (2 * (S * dmusig * sigx14) +
      pmusig * gid2_gicub/sigmasq) * gid1_epsi_gisq) *
      ewu))/sigx3 + 0.5 * gd2Zsq_gd1_gicub))/sigmasq))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (ewu * (ewu * (S^2 * ((0.5 *
    (ewu * gid2_gicub * giepsi/sigmasq) + gid1_epsi_gisq) *
    sigx18 + (gd1Zcub_gd1_epsi_gisq/(ssqx2) - sigx41 * ewu/ssq) *
    giepsi)/sigmastar - sigx28 * gid2_gicub/sigmasq) - ((((sigx15 *
    ewu + S * sigx16) * sigx20/(sigx3 * sigmasq) + 0.5 *
    (((dmusig * gd2Zcub_gd1_gicub - S^2 * dmusig * ewu^2 *
      gid2sq_gicub * sigx14 * giepsi/(ssqx2))/(ssqx2) -
      dmusig * ewu * sigx13 * gid2_gicub * gid2sq_gicub/ssq) *
      ewv) + S^2 * (sigx42 * giepsi + gid1_epsi_gisq *
    sigx18) * dmusig) * ewu + S * (pmusig * gd1Zcub_gd1_epsi_gisq -
    (((pmusig * gd2Zcub_gd1_gicub - ewu * (pmusig * gid2_gicub/sigmasq +
      S * dmusig * sigx14) * gid2sq_gicub)/sigmasq + S *
      dmusig * sigx42) * giepsi + (pmusig * gid2sq_gicub/sigmasq +
      S * dmusig * sigx18) * gid1_epsi_gisq + S * dmusig *
      sigx14 * gid1sq_epsi_gisq) * ewu))/sigx3 + 0.5 *
    gd2Zcub_gd1_gicub))/sigmasq))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (ewu * (ewu * (S^2 * ((0.5 *
    (ewu * gid2sq_gicub * giepsi/sigmasq) + gid1sq_epsi_gisq) *
    sigx18 + (gd1Zfour_gd1_epsi_gisq/(ssqx2) - sigx38 * ewu/ssq) *
    giepsi)/sigmastar - sigx28 * gid2sq_gicub/sigmasq) -
    (((sigx20^2/(sigx3 * sigmasq) + 0.5 * (((dmusig * gd2Zfour_gd1_gicub -
      S^2 * dmusig * ewu^2 * gid2sq_gicub * sigx18 * giepsi/(ssqx2))/(ssqx2) -
      dmusig * ewu * sigx13 * gid2sq_gicub^2/ssq) * ewv) +
      S^2 * (sigx39 * giepsi + gid1sq_epsi_gisq * sigx18) *
        dmusig) * ewu + S * (pmusig * gd1Zfour_gd1_epsi_gisq -
      (((pmusig * gd2Zfour_gd1_gicub - ewu * (pmusig *
        gid2sq_gicub/sigmasq + S * dmusig * sigx18) *
        gid2sq_gicub)/sigmasq + S * dmusig * sigx39) *
        giepsi + (2 * (S * dmusig * sigx18) + pmusig *
        gid2sq_gicub/sigmasq) * gid1sq_epsi_gisq) * ewu))/sigx3 +
      0.5 * gd2Zfour_gd1_gicub))/sigmasq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for rayleigh-normal distribution
#' @param start starting value for optimization
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
raynormAlgOpt_k90 <- function(start, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo,
  whichStart, initIter, initAlg, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstraynorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar_c, tol = tol, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitRay <- start_st$initRay
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(praynormlike_k90(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel K90 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) {
      -sum(praynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_k90(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = praynormlike_k90,
    grad = pgradraynormlike_k90, hess = phessraynormlike_k90,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_k90(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(praynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_k90(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessraynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(praynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradraynormlike_k90(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessraynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(praynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradraynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phessraynormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradraynormlike_k90(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phessraynormlike_k90(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessraynormlike_k90(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- praynormlike_k90(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradraynormlike_k90(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, InitRay = InitRay))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for rayleigh-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
praynormeff_k90 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  eta2 <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    2]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
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
  epsilon_it <- model.response(model.frame(object$formula,
    data = object$dataTable)) - as.numeric(crossprod(matrix(beta),
    t(Xvar)))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 *
    seq(1:x) + eta2 * (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- -exp(Wu) * object$S * giepsi/(exp(Wv) + gisq *
    exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 +
    sigmastar^2) * pnorm(mustar/sigmastar))/(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0,
    mustar/2 + sqrt(sigmastar^2 + mustar^2/4), 0)
  res <- data.frame(levels(pindex[, 1]), u = u, m = m, mustar = mustar,
    sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(-res$mustar * git + res$sigmastar^2 *
      git^2/2) * (res$sigmastar * git * dnorm(res$mustar/res$sigmastar -
      res$sigmastar * git) + (res$mustar * git - res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) +
      res$mustar * git * pnorm(res$mustar/res$sigmastar))
    res$teBC_reciprocal <- exp(res$mustar * git + res$sigmastar^2 *
      git^2/2) * (res$sigmastar * git * dnorm(res$mustar/res$sigmastar +
      res$sigmastar * git) + (res$mustar * git + res$sigmastar^2 *
      git^2) * pnorm(res$mustar/res$sigmastar + res$sigmastar *
      git))/(res$sigmastar * git * dnorm(res$mustar/res$sigmastar) +
      res$mustar * git * pnorm(res$mustar/res$sigmastar))
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
