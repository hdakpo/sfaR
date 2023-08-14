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
# Convolution: truncated normal - normal                                       #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
ptruncnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar, vHvar,
  muHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 2]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- (mu * exp(Wv) - exp(Wu) * S * giepsi)/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- -TT/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 * log(exp(Wv) + gisq * exp(Wu)) +
    pnorm(mustar/sigmastar, log.p = TRUE) - epsilon_isq/(2 * exp(Wv)) + 1/2 *
    (mustar/sigmastar)^2 - 1/2 * (mu^2/exp(Wu)) - pnorm(mu/sqrt(exp(Wu)), log.p = TRUE)
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for truncatednormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
psttruncnorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, muHvar, nmuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- csttruncnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nmuZUvar = 1,
      nuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTrunc <- NULL
  } else {
    cat("Initialization: SFA + truncated-normal distribution...\n")
    initTrunc <- maxLik::maxLik(logLik = ctruncnormlike, start = csttruncnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nmuZUvar = 1, nuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1,
      vHvar = vHvar[, 1, drop = FALSE]), grad = cgradtruncnormlike, hess = chesstruncnormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nmuZUvar = 1, nuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1,
      vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initTrunc$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nmuZUvar > 1) {
    rep(0, nmuZUvar - 1)
  }, Esti[nXvar + 2], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 3], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_", colnames(muHvar)),
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "eta1", "eta2")
  return(list(StartVal = StartVal, initTrunc = initTrunc))
}

# Gradient of the likelihood function ----------
#' gradient for truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgradtruncnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
  vHvar, muHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 2]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  mug <- (mu * ewv - S * ewu * giepsi)
  dmusig <- dnorm(mug/(ssqx2), 0, 1)
  pmusig <- pnorm(mug/(ssqx2))
  sigx1 <- (mug/ewv + dmusig * ewu/(pmusig * sigmastar))
  sigx2 <- (mug/ewu + dmusig * ewv/(pmusig * sigmastar))
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  sigx3 <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (sigx3 * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (mug/(ssqx2) + dmusig/pmusig)
  sigx6 <- (sigx4 * mug/ssq + S * giepsi/(ssqx2))
  prV <- (1 - ewv/sigmasq)
  sigx7 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx8 <- (mu/(ssqx2) - sigx7 * mug/ssq)
  sigx9 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  gradll <- cbind(sweep(Xepsi_i, MARGIN = 1, STATS = -1/(2 * ewv), FUN = "*") +
    sweep(Xgi, MARGIN = 1, STATS = -S * sigx1/sigmasq, FUN = "*"), sweep(muHvar,
    MARGIN = 1, STATS = (sigx2/sigmasq - (dmu/(ewu_h * pmu) + mu/ewu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = -(((sigx6 * sigx5 + 0.5 * (gisq/sigmasq)) *
      ewu - (0.5 * ((mu)^2/ewu) + 0.5 * (mu * dmu/(ewu_h * pmu))))), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ((sigx5 * sigx8 + 2 * (epsilon_isq/(2 *
      ewv)^2) - 0.5/sigmasq) * ewv - 0.5 * (TT - 1)), FUN = "*"), -(((mug *
      sigx9 * gid2_gicub/ssq + S * gid1_epsi_gisq/(ssqx2)) * sigx5 + 0.5 *
      (gid2_gicub/sigmasq)) * ewu), -(((mug * sigx9 * gid2sq_gicub/ssq + S *
      gid1sq_epsi_gisq/(ssqx2)) * sigx5 + 0.5 * (gid2sq_gicub/sigmasq)) * ewu))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for truncatednormal-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phesstruncnormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
  vHvar, muHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  eta1 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 1]
  eta2 <- parm[nXvar + nmuZUvar + nuZUvar + nvZVvar + 2]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
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
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  Xit_gitd1_gitsq <- sweep(Xvar, MARGIN = 1, STATS = gitd1 * git^2, FUN = "*")
  Xi_gid1_gisq <- apply(Xit_gitd1_gitsq, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xit_gitd1sq_gitsq <- sweep(Xvar, MARGIN = 1, STATS = Zit^2 * exp(Zit * (eta1 +
    eta2 * Zit)) * git^2, FUN = "*")
  Xi_gid1sq_gisq <- apply(Xit_gitd1sq_gitsq, 2, function(x) tapply(x, pindex[,
    1], sum))
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  mug <- (mu * ewv - S * ewu * giepsi)
  dmusig <- dnorm(mug/(ssqx2), 0, 1)
  pmusig <- pnorm(mug/(ssqx2))
  sigx1 <- (mug/ewv + dmusig * ewu/(pmusig * sigmastar))
  sigx2 <- (mug/ewu + dmusig * ewv/(pmusig * sigmastar))
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  sigx3 <- (1 - ewu * gisq/sigmasq)
  sigx4 <- (0.5 * (sigx3 * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (mug/(ssqx2) + dmusig/pmusig)
  sigx6 <- (sigx4 * mug/ssq + S * giepsi/(ssqx2))
  prV <- (1 - ewv/sigmasq)
  sigx7 <- (0.5 * (prV * ewu/sigmastar) + sigmastar)
  sigx8 <- (mu/(ssqx2) - sigx7 * mug/ssq)
  sigx9 <- (sigmastar - 0.5 * (ewu * ewv/(ssqx2)))
  sigx10 <- (sigx1 * dmusig/pmusig - ewu/sigmastar)
  sigx11 <- (1/(ssqx2) - sigx4 * ewu/ssq)
  sigx12 <- (mug * sigx9 * gid2_gicub/ssq + S * gid1_epsi_gisq/(ssqx2))
  sigx13 <- (mug * sigx9 * gid2sq_gicub/ssq + S * gid1sq_epsi_gisq/(ssqx2))
  sigx14 <- (ewv/sigmastar - sigx2 * dmusig/pmusig)
  sigx15 <- (1/(ssqx2) - sigx7 * ewv/ssq)
  sigx16 <- (sigx5 * dmusig/pmusig - 1)
  sigx17 <- (sigx4 * sigmasq * mug * sigmastar/ssq)
  sigx18 <- (0.5 * (mu^2/(ewu_h^3 * pmu)) - (0.5 * (ewu_h * pmu) - 0.5 * (mu *
    dmu))/(ewu_h * pmu)^2)
  sigx19 <- (((0.5 * (sigx3/(ssqx2)) - 0.5 * sigx11) * ewv - 2 * (sigx4 * ssqx2 *
    sigx9/ssq)) * mug - S * sigx9 * giepsi)
  sigx20 <- (((0.5 * (prV/(ssqx2)) - 0.5 * sigx15) * ewu - 2 * (sigx7 * ssqx2 *
    sigx9/ssq)) * mug + mu * sigx9)
  sigx21 <- (1 - sigx5 * dmusig/pmusig)
  sigx22 <- ((0.5 * (sigx9/ssq) - 0.5/(sigmasq^2 * sigmastar)) * ewu * ewv - 2 *
    (ssqx2 * sigx9^2/ssq))
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar + 2, ncol = nXvar +
    nmuZUvar + nuZUvar + nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(-wHvar/(2 *
    ewv)))) + crossprod(sweep(Xgi, MARGIN = 1, STATS = -wHvar * S^2 * ((mug/(ewv *
    pmusig * sigmastar) + dmusig * ewu/(pmusig * sigmastar)^2) * dmusig/sigmasq -
    1/ewv) * ewu/sigmasq, FUN = "*"), Xgi)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * ((mug/(pmusig * sigmastar) + dmusig * ewu * ewv/(pmusig *
      sigmastar)^2) * dmusig/sigmasq - 1)/sigmasq, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = -wHvar * (S * (sigx10 * sigx6/sigmasq + sigx5 * sigx11) *
      ewu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = wHvar * 2/(2 *
    ewv)^2 * ewv, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar * S * (sigx10 *
    sigx8/sigmasq + sigx5 * sigx7 * ewu/ssq) * ewv, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1)] <- -colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx10 * sigx12/sigmasq * ewu, FUN = "*") +
    sweep(Xi_gid1_gisq, MARGIN = 1, STATS = wHvar * S * sigx5/(ssqx2) * ewu,
      FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * S * sigx5 * ewu *
    sigx9 * gid2_gicub/ssq * ewu, FUN = "*"))
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 2)] <- -colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * sigx10 * sigx13/sigmasq * ewu, FUN = "*") +
    sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = wHvar * S * sigx5/(ssqx2) * ewu,
      FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * S * sigx5 * ewu *
    sigx9 * gid2sq_gicub/ssq * ewu, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((1/ewu - (mug/(ewu * pmusig * sigmastar) + dmusig *
      ewv/(pmusig * sigmastar)^2) * dmusig/sigmasq) * ewv/sigmasq + dmu * (dmu/(ewu_h *
      pmu)^2 + mu/(ewu_h^3 * pmu)) - 1/ewu), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((sigx6 *
    sigx14/sigmasq + sigx5 * sigx4 * ewv/ssq) * ewu - (0.5 * (((1 - mu^2/ewu_h^2)/(ewu_h *
    pmu) - mu * dmu/(ewu_h * pmu)^2) * dmu) + mu/ewu))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (sigx5 * sigx15 + sigx14 * sigx8/sigmasq) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1)] <- -colSums(sweep(muHvar, MARGIN = 1, STATS = wHvar * (sigx12 * sigx14/sigmasq +
    sigx5 * ewv * sigx9 * gid2_gicub/ssq) * ewu, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    2)] <- -colSums(sweep(muHvar, MARGIN = 1, STATS = wHvar * (sigx13 * sigx14/sigmasq +
    sigx5 * ewv * sigx9 * gid2sq_gicub/ssq) * ewu, FUN = "*"))
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((((sigx16 * sigx6^2 - 0.5 * (gisq^2/sigmasq^2)) * ewu + ((sigx4 * (mu *
      ewv - (2 * sigx17 + 3 * (S * giepsi)) * ewu) + (0.5 * (ewu * gisq/sigmasq) -
      0.5 * (0.5 * sigx3 + ewu * gisq/sigmasq)) * sigx3 * ewv * mug/sigmastar)/ssq +
      S * giepsi/(ssqx2)) * sigx5 + 0.5 * (gisq/sigmasq)) * ewu + 0.5 * ((mu)^2/ewu) -
      0.5 * (mu * sigx18 * dmu))), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx16 * sigx6 * sigx8 + 0.5 * (gisq/sigmasq^2) -
      ((0.5 * (sigx3 * ewv/sigmasq) + 0.5 * ((ewu * gisq/sigmasq - 1) * ewv/sigmasq +
        1 - 0.5 * (sigx3 * prV))) * mug/sigmastar + mu * sigx4 - sigx7 *
        (2 * sigx17 + S * giepsi)) * sigx5/ssq) * ewu * ewv, FUN = "*"),
    vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1)] <- -colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx19 * gid2_gicub - S * sigx4 * gid1_epsi_gisq) * sigx5/ssq + sigx16 *
      sigx6 * sigx12 - 0.5 * (gid2_gicub * gisq/sigmasq^2)) * ewu + sigx12 *
      sigx5 + 0.5 * (gid2_gicub/sigmasq)) * ewu, FUN = "*"))
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 2)] <- -colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx19 * gid2sq_gicub - S * sigx4 * gid1sq_epsi_gisq) * sigx5/ssq + sigx16 *
      sigx6 * sigx13 - 0.5 * (gid2sq_gicub * gisq/sigmasq^2)) * ewu + sigx13 *
      sigx5 + 0.5 * (gid2sq_gicub/sigmasq)) * ewu, FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx21 * sigx8^2 + 0.5/sigmasq^2 - 16 * (ewv *
      epsilon_isq/(2 * ewv)^4)) * ewv + sigx5 * (mu/(ssqx2) - (((3 * (mu) -
      2 * (sigx7 * sigmasq * mug * sigmastar/ssq)) * ewv - S * ewu * giepsi) *
      sigx7 + (0.5 * (ewv/sigmasq) - 0.5 * (0.5 * prV + ewv/sigmasq)) * prV *
      ewu * mug/sigmastar)/ssq) + 2 * (epsilon_isq/(2 * ewv)^2) - 0.5/sigmasq) *
      ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1)] <- -colSums(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx20 * gid2_gicub - S * sigx7 * gid1_epsi_gisq) * sigx5/ssq +
      sigx12 * sigx21 * sigx8 - 0.5 * (gid2_gicub/sigmasq^2)) * ewu * ewv,
    FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 2)] <- -colSums(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ((sigx20 * gid2sq_gicub - S * sigx7 * gid1sq_epsi_gisq) *
      sigx5/ssq + sigx13 * sigx21 * sigx8 - 0.5 * (gid2sq_gicub/sigmasq^2)) *
      ewu * ewv, FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 1)] <- sum(wHvar * (-(((((sigx22 * mug * gid2_gicub - S * sigx9 *
    gid1_epsi_gisq) * ewu * gid2_gicub + mug * sigx9 * gd2Zsq_gd1_gicub)/ssq +
    S * (gd1Zsq_gd1_epsi_gisq/(ssqx2) - ewu * sigx9 * gid2_gicub * gid1_epsi_gisq/ssq)) *
    sigx5 + sigx16 * sigx12^2 * ewu + 0.5 * ((gd2Zsq_gd1_gicub - ewu * gid2_gicub^2/sigmasq)/sigmasq)) *
    ewu)))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (-(((((sigx22 * mug * gid2_gicub - S * sigx9 *
    gid1_epsi_gisq) * ewu * gid2sq_gicub + mug * sigx9 * gd2Zcub_gd1_gicub)/ssq +
    S * (gd1Zcub_gd1_epsi_gisq/(ssqx2) - ewu * sigx9 * gid2_gicub * gid1sq_epsi_gisq/ssq)) *
    sigx5 + sigx16 * sigx12 * sigx13 * ewu + 0.5 * ((gd2Zcub_gd1_gicub - ewu *
    gid2_gicub * gid2sq_gicub/sigmasq)/sigmasq)) * ewu)))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 2), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (-(((((sigx22 * mug * gid2sq_gicub - S * sigx9 *
    gid1sq_epsi_gisq) * ewu * gid2sq_gicub + mug * sigx9 * gd2Zfour_gd1_gicub)/ssq +
    S * (gd1Zfour_gd1_epsi_gisq/(ssqx2) - ewu * sigx9 * gid2sq_gicub * gid1sq_epsi_gisq/ssq)) *
    sigx5 + sigx16 * sigx13^2 * ewu + 0.5 * ((gd2Zfour_gd1_gicub - ewu * gid2sq_gicub^2/sigmasq)/sigmasq)) *
    ewu)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncatednormal-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
truncnormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, muHvar_c, muHvar_p, nmuZUvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p,
  nvZVvar, Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttruncnorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_c, uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      wHvar = wHvar_c, tol = tol, printInfo = printInfo)
    Inittrunc <- start_st$inittrunc
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ptruncnormlike_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel K90 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(ptruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptruncnormlike_k90,
    grad = pgradtruncnormlike_k90, hess = phesstruncnormlike_k90, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(ptruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(ptruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(ptruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradtruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesstruncnormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtruncnormlike_k90(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phesstruncnormlike_k90(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstruncnormlike_k90(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptruncnormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtruncnormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, Inittrunc = Inittrunc))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncatednormal-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
ptruncnormeff_k90 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  eta1 <- object$mlParam[object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1]
  eta2 <- object$mlParam[object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    2]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  pindex <- object$dataTable[, 1:2]
  invariance <- object$invariance
  if (invariance == 1) {
    muHvar_p <- apply(muHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    uHvar_p <- apply(uHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
    vHvar_p <- apply(vHvar_c, 2, function(x) {
      tapply(x, pindex[, 1], function(u) u[1])
    })
  } else {
    if (invariance == 2) {
      muHvar_p <- apply(muHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[length(u)])
      })
    } else {
      if (invariance == 3) {
        muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
        })
      }
    }
  }
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar_p)))
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
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * giepsi)/(exp(Wv) + gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) * sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB, uUB = uUB, m = m, mustar = mustar,
    sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  res$uLB <- res$uLB * git
  res$uUB <- res$uUB * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(1/2 * res$sigmastar^2 * git^2 - res$mustar * git) * pnorm(res$mustar/res$sigmastar -
      res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    res$teBC_reciprocal <- exp(1/2 * res$sigmastar^2 * git^2 + res$mustar * git) *
      pnorm(res$mustar/res$sigmastar + res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
