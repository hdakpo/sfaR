################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Type:      - Battese and Coelli 1992 - alternative (mbc92)                   #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = 1 + eta1 * (t - T) + eta2 * (t - T)^2                  #
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
ptruncnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
psttruncnorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, muHvar, nmuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter,
  initAlg) {
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
pgradtruncnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit <- 2 * (Zit * (1 + Zit * (eta1 + eta2 * Zit)))
  giZi <- as.numeric(tapply(gitZit, pindex[, 1], sum))
  gitZitsq <- 2 * (Zit^2 * (1 + Zit * (eta1 + eta2 * Zit)))
  giZisq <- as.numeric(tapply(gitZitsq, pindex[, 1], sum))
  Zit_epsit <- Zit * epsilon_it
  Zi_epsi <- as.numeric(tapply(Zit_epsit, pindex[, 1], sum))
  Zitsq_epsit <- Zit^2 * epsilon_it
  Zisq_epsi <- as.numeric(tapply(Zitsq_epsit, pindex[, 1], sum))
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  mug <- (mu * ewv - S * ewu * giepsi)
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  dmusig <- dnorm(mug/(sigmasq * sigmastar), 0, 1)
  pmusig <- pnorm(mug/(sigmasq * sigmastar))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx1 <- (mug/ewv + dmusig * ewu/(pmusig * sigmastar))
  sigx2 <- (mug/ewu + dmusig * ewv/(pmusig * sigmastar))
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (sigx4 * mug/ssq + S * giepsi/(sigmasq * sigmastar))
  sigx6 <- (mug/(sigmasq * sigmastar) + dmusig/pmusig)
  sigx7 <- (0.5 * ((1 - ewv/sigmasq) * ewu/sigmastar) + sigmastar)
  sigx8 <- (mu/(sigmasq * sigmastar) - sigx7 * mug/ssq)
  gradll <- cbind(-(sweep(Xgi, MARGIN = 1, STATS = S * sigx1/sigmasq, FUN = "*") +
    sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*")), sweep(muHvar,
    MARGIN = 1, STATS = (sigx2/sigmasq - (dmu/(ewu_h * pmu) + mu/ewu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = -(((sigx5 * sigx6 + 0.5 * (gisq/sigmasq)) *
      ewu - (0.5 * ((mu)^2/ewu) + 0.5 * (mu * dmu/(ewu_h * pmu))))), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ((sigx6 * sigx8 + 2 * (epsilon_isq/(2 *
      ewv)^2) - 0.5/sigmasq) * ewv - 0.5 * (TT - 1)), FUN = "*"), -(((mug *
      sigx3 * giZi/ssq + S * Zi_epsi/(sigmasq * sigmastar)) * sigx6 + 0.5 *
      (giZi/sigmasq)) * ewu), -(((mug * sigx3 * giZisq/ssq + S * Zisq_epsi/(sigmasq *
      sigmastar)) * sigx6 + 0.5 * (giZisq/sigmasq)) * ewu))
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
phesstruncnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, nmuZUvar, uHvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zit <- unlist(lapply(TT, FUN = function(x) rev(-(0:(x - 1)))))
  gitZit <- 2 * (Zit * (1 + Zit * (eta1 + eta2 * Zit)))
  giZi <- as.numeric(tapply(gitZit, pindex[, 1], sum))
  gitZitsq <- 2 * (Zit^2 * (1 + Zit * (eta1 + eta2 * Zit)))
  giZisq <- as.numeric(tapply(gitZitsq, pindex[, 1], sum))
  Zit_epsit <- Zit * epsilon_it
  Zi_epsi <- as.numeric(tapply(Zit_epsit, pindex[, 1], sum))
  Zitsq_epsit <- Zit^2 * epsilon_it
  Zisq_epsi <- as.numeric(tapply(Zitsq_epsit, pindex[, 1], sum))
  Zisq <- as.numeric(tapply(2 * Zit^2, pindex[, 1], sum))
  Zicub <- as.numeric(tapply(2 * Zit^3, pindex[, 1], sum))
  Zifour <- as.numeric(tapply(2 * Zit^4, pindex[, 1], sum))
  XitZit <- sweep(-Xvar, MARGIN = 1, STATS = Zit, FUN = "*")
  XiZi <- apply(XitZit, 2, function(x) tapply(x, pindex[, 1], sum))
  XitZitsq <- sweep(-Xvar, MARGIN = 1, STATS = Zit^2, FUN = "*")
  XiZisq <- apply(XitZitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewv <- exp(Wv)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  mug <- (mu * ewv - S * ewu * giepsi)
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  dmusig <- dnorm(mug/(sigmasq * sigmastar), 0, 1)
  pmusig <- pnorm(mug/(sigmasq * sigmastar))
  wu2sq <- (1 - ewu * gisq/sigmasq)
  sigx1 <- (mug/ewv + dmusig * ewu/(pmusig * sigmastar))
  sigx2 <- (mug/ewu + dmusig * ewv/(pmusig * sigmastar))
  sigx3 <- (sigmastar - 0.5 * (ewu * ewv/(sigmasq * sigmastar)))
  sigx4 <- (0.5 * (wu2sq * ewv/sigmastar) + sigmastar * gisq)
  sigx5 <- (sigx4 * mug/ssq + S * giepsi/(sigmasq * sigmastar))
  sigx6 <- (mug/(sigmasq * sigmastar) + dmusig/pmusig)
  sigx7 <- (0.5 * ((1 - ewv/sigmasq) * ewu/sigmastar) + sigmastar)
  sigx8 <- (mu/(sigmasq * sigmastar) - sigx7 * mug/ssq)
  sigx9 <- ((1 - ewv/sigmasq)/(sigmasq * sigmastar))
  sigx10 <- (sigx1 * dmusig/pmusig - ewu/sigmastar)
  sigx11 <- (1/(sigmasq * sigmastar) - sigx4 * ewu/ssq)
  sigx12 <- (mug * sigx3 * giZi/ssq + S * Zi_epsi/(sigmasq * sigmastar))
  sigx13 <- (mug * sigx3 * giZisq/ssq + S * Zisq_epsi/(sigmasq * sigmastar))
  sigx14 <- (ewv/sigmastar - sigx2 * dmusig/pmusig)
  sigx15 <- (1/(sigmasq * sigmastar) - sigx7 * ewv/ssq)
  sigx16 <- ((0.5 * (sigx3/ssq) - 0.5/(sigmasq^2 * sigmastar)) * ewu * ewv - 2 *
    (sigmasq * sigmastar * sigx3^2/ssq))
  sigx17 <- (((0.5 * sigx9 - 0.5 * sigx15) * ewu - 2 * (sigx7 * sigmasq * sigmastar *
    sigx3/ssq)) * mug + mu * sigx3)
  sigx18 <- (1 - sigx6 * dmusig/pmusig)
  sigx19 <- (sigx6 * dmusig/pmusig - 1)
  sigx20 <- (((0.5 * (wu2sq/(sigmasq * sigmastar)) - 0.5 * sigx11) * ewv - 2 *
    (sigx4 * sigmasq * sigmastar * sigx3/ssq)) * mug - S * sigx3 * giepsi)
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar + 2, ncol = nXvar +
    nmuZUvar + nuZUvar + nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- -(crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * ((mug/(ewv * pmusig * sigmastar) + dmusig * ewu/(pmusig * sigmastar)^2) *
    dmusig/sigmasq - 1/ewv) * ewu/sigmasq, FUN = "*"), Xgi) + sapply(1:nXvar,
    function(x) crossprod(Xsq[[x]], as.matrix(wHvar/(2 * ewv)))))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * ((mug/(pmusig * sigmastar) + dmusig * ewu * ewv/(pmusig *
      sigmastar)^2) * dmusig/sigmasq - 1)/sigmasq, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = -wHvar * (S * (sigx10 * sigx5/sigmasq + sigx6 * sigx11) *
      ewu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xepsi_i, MARGIN = 1, STATS = wHvar * 2/(2 *
    ewv)^2 * ewv, FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = wHvar * S * (sigx10 *
    sigx8/sigmasq + sigx6 * sigx7 * ewu/ssq) * ewv, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1)] <- -colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx10 * sigx12/sigmasq * S * ewu, FUN = "*") +
    sweep(XiZi, MARGIN = 1, STATS = wHvar/(sigmasq * sigmastar) * sigx6 * S *
      ewu, FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * ewu * sigx3 *
    giZi/ssq * sigx6 * S * ewu, FUN = "*"))
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 2)] <- -colSums(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * sigx10 * sigx13/sigmasq * S * ewu, FUN = "*") +
    sweep(XiZisq, MARGIN = 1, STATS = wHvar/(sigmasq * sigmastar) * sigx6 * S *
      ewu, FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * ewu * sigx3 *
    giZisq/ssq * sigx6 * S * ewu, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((1/ewu - (mug/(ewu * pmusig * sigmastar) + dmusig *
      ewv/(pmusig * sigmastar)^2) * dmusig/sigmasq) * ewv/sigmasq + dmu * (dmu/(ewu_h *
      pmu)^2 + mu/(ewu_h^3 * pmu)) - 1/ewu), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((sigx5 *
    sigx14/sigmasq + sigx6 * sigx4 * ewv/ssq) * ewu - (0.5 * (((1 - mu^2/ewu_h^2)/(ewu_h *
    pmu) - mu * dmu/(ewu_h * pmu)^2) * dmu) + mu/ewu))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (sigx6 * sigx15 + sigx14 * sigx8/sigmasq) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1)] <- colSums(sweep(muHvar, MARGIN = 1, STATS = -wHvar * ((sigx12 * sigx14/sigmasq +
    sigx6 * ewv * sigx3 * giZi/ssq) * ewu), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    2)] <- colSums(sweep(muHvar, MARGIN = 1, STATS = -wHvar * ((sigx13 * sigx14/sigmasq +
    sigx6 * ewv * sigx3 * giZisq/ssq) * ewu), FUN = "*"))
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    -((((sigx19 * sigx5^2 - 0.5 * (gisq^2/sigmasq^2)) * ewu + ((sigx4 * (mu *
      ewv - (2 * (sigx4 * sigmasq * mug * sigmastar/ssq) + 3 * (S * giepsi)) *
      ewu) + (0.5 * (ewu * gisq/sigmasq) - 0.5 * (0.5 * wu2sq + ewu * gisq/sigmasq)) *
      wu2sq * ewv * mug/sigmastar)/ssq + S * giepsi/(sigmasq * sigmastar)) *
      sigx6 + 0.5 * (gisq/sigmasq)) * ewu + 0.5 * ((mu)^2/ewu) - 0.5 * (mu *
      (0.5 * (mu^2/(ewu_h^3 * pmu)) - (0.5 * (ewu_h * pmu) - 0.5 * (mu * dmu))/(ewu_h *
        pmu)^2) * dmu))), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx19 * sigx5 * sigx8 + 0.5 * (gisq/sigmasq^2) -
      ((0.5 * (wu2sq * ewv/sigmasq) + 0.5 * ((ewu * gisq/sigmasq - 1) * ewv/sigmasq +
        1 - 0.5 * (wu2sq * (1 - ewv/sigmasq)))) * mug/sigmastar + mu * sigx4 -
        sigx7 * (2 * (sigx4 * sigmasq * mug * sigmastar/ssq) + S * giepsi)) *
        sigx6/ssq) * ewu * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((((sigx20 * giZi - S * sigx4 * Zi_epsi) * sigx6/ssq + sigx19 * sigx5 * sigx12 -
      0.5 * (gisq * giZi/sigmasq^2)) * ewu + sigx12 * sigx6 + 0.5 * (giZi/sigmasq)) *
      ewu), FUN = "*"))
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((((sigx20 * giZisq - S * sigx4 * Zisq_epsi) * sigx6/ssq + sigx19 * sigx5 *
      sigx13 - 0.5 * (gisq * giZisq/sigmasq^2)) * ewu + sigx13 * sigx6 + 0.5 *
      (giZisq/sigmasq)) * ewu), FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx18 * sigx8^2 + 0.5/sigmasq^2 - 16 * (ewv *
      epsilon_isq/(2 * ewv)^4)) * ewv + sigx6 * (mu/(sigmasq * sigmastar) -
      (((3 * (mu) - 2 * (sigx7 * sigmasq * mug * sigmastar/ssq)) * ewv - S *
        ewu * giepsi) * sigx7 + (0.5 * (ewv/sigmasq) - 0.5 * (0.5 * (1 -
        ewv/sigmasq) + ewv/sigmasq)) * (1 - ewv/sigmasq) * ewu * mug/sigmastar)/ssq) +
      2 * (epsilon_isq/(2 * ewv)^2) - 0.5/sigmasq) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (((sigx17 * giZi - S * sigx7 * Zi_epsi) * sigx6/ssq + sigx12 *
      sigx18 * sigx8 - 0.5 * (giZi/sigmasq^2)) * ewu * ewv), FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (((sigx17 * giZisq - S * sigx7 * Zisq_epsi) * sigx6/ssq +
      sigx13 * sigx18 * sigx8 - 0.5 * (giZisq/sigmasq^2)) * ewu * ewv), FUN = "*"))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 1)] <- sum(wHvar * (-((((sigx16 * mug * giZi - 2 * (S * sigx3 *
    Zi_epsi)) * ewu * giZi + mug * sigx3 * Zisq) * sigx6/ssq + sigx19 * sigx12^2 *
    ewu + 0.5 * ((Zisq - ewu * giZi^2/sigmasq)/sigmasq)) * ewu)))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (-(((((sigx16 * mug * giZi - S * sigx3 * Zi_epsi) *
    giZisq - S * sigx3 * giZi * Zisq_epsi) * ewu + mug * sigx3 * Zicub) * sigx6/ssq +
    sigx19 * sigx12 * sigx13 * ewu + 0.5 * ((Zicub - ewu * giZi * giZisq/sigmasq)/sigmasq)) *
    ewu)))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 2), (nXvar + nmuZUvar + nuZUvar +
    nvZVvar + 2)] <- sum(wHvar * (-((((sigx16 * mug * giZisq - 2 * (S * sigx3 *
    Zisq_epsi)) * ewu * giZisq + mug * sigx3 * Zifour) * sigx6/ssq + sigx19 *
    sigx13^2 * ewu + 0.5 * ((Zifour - ewu * giZisq^2/sigmasq)/sigmasq)) * ewu)))
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
truncnormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, muHvar_c, muHvar_p, nmuZUvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p,
  nvZVvar, Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttruncnorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_c, uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, whichStart = whichStart, initIter = initIter, initAlg = initAlg,
      wHvar = wHvar_c, tol = tol, printInfo = printInfo)
    Inittrunc <- start_st$inittrunc
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ptruncnormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel Modified BC92 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(ptruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptruncnormlike_mbc92,
    grad = pgradtruncnormlike_mbc92, hess = phesstruncnormlike_mbc92, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(ptruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(ptruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(ptruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradtruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesstruncnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtruncnormlike_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phesstruncnormlike_mbc92(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstruncnormlike_mbc92(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptruncnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtruncnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
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
ptruncnormeff_mbc92 <- function(object, level) {
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
