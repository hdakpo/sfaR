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
ptslnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
psttslnorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
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
pgradtslnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  ewu_h <- exp(Wu/2)
  sgsq <- (sqrt(ewv/gisq) * gisq)
  llepsi <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- (llepsi/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  llepsi2 <- (llepsi * pmusig2 - dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (2 * (sigx3 * pmusig1) - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * ((ewv/ewu_h + S * giepsi) *
    pmusig1))
  sigx6 <- (0.5 * (dmusig2 * ewv/sqrt(ewv/gisq)) - 0.5 * (llepsi * pmusig2))
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx12 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (llepsi/sgsq^2))
  sigx10 <- (llepsi * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx13 <- (2 * (sigx9 * sigx3) - (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) *
    sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * Zisq_epsi/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * giZisq/sgsq^2)
  sigx16 <- (llepsi * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * Zisq_epsi/sgsq - llepsi * (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) *
    giZisq/sgsq^2)
  sigx18 <- (2 * (sigx14 * sigx3 * sigx15) - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * Zi_epsi/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) - 0.5 *
    (ewv/sgsq)) * giZi/sgsq^2)
  sigx20 <- (S * Zi_epsi/sgsq - llepsi * (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) *
    giZi/sgsq^2)
  sigx21 <- (2 * (sigx14 * sigx3 * sigx19) - sigx16 * sigx2 * sigx20)
  Xsig0 <- sweep(Xgi, MARGIN = 1, STATS = (S * llepsi/gisq), FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * (ewv/ewu_h + S * giepsi)/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sgsq, FUN = "*")
  Xsig4 <- sweep(Xepsi_i - 2 * Xsig0, MARGIN = 1, STATS = (pmusig2/ewv), FUN = "*")
  Xsig5 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sgsq, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig4 + Xsig2), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1, STATS = sigx3, FUN = "*")
  gradll <- cbind(sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = 1/sigx4, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ((2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) *
      sigx2)/(sigx4 * ewu_h * gisq) - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = (sigx13/sigx4 - 0.5 * (TT - 1)), FUN = "*"), 1/(1 + lambda) -
      (llepsi2 * sigx2/(sigx4 * ewu_h * gisq) + 2/(1 + 2 * lambda)), sigx21/sigx4 -
      0.5 * (giZi/gisq), sigx18/sigx4 - 0.5 * (giZisq/gisq))
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
phesstslnormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  ewu_h <- exp(Wu/2)
  sgsq <- (sqrt(ewv/gisq) * gisq)
  llepsi <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- (llepsi/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  llepsi2 <- (llepsi * pmusig2 - dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (2 * (sigx3 * pmusig1) - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * ((ewv/ewu_h + S * giepsi) *
    pmusig1))
  sigx6 <- (0.5 * (dmusig2 * ewv/sqrt(ewv/gisq)) - 0.5 * (llepsi * pmusig2))
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx12 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (llepsi/sgsq^2))
  sigx10 <- (llepsi * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx13 <- (2 * (sigx9 * sigx3) - (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) *
    sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * Zisq_epsi/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * giZisq/sgsq^2)
  sigx16 <- (llepsi * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * Zisq_epsi/sgsq - llepsi * (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) *
    giZisq/sgsq^2)
  sigx18 <- (2 * (sigx14 * sigx3 * sigx15) - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * Zi_epsi/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) - 0.5 *
    (ewv/sgsq)) * giZi/sgsq^2)
  sigx20 <- (S * Zi_epsi/sgsq - llepsi * (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq)) *
    giZi/sgsq^2)
  sigx21 <- (2 * (sigx14 * sigx3 * sigx19) - sigx16 * sigx2 * sigx20)
  pdmusig2 <- (pmusig2 - llepsi * dmusig2/sgsq)
  sgvsq <- (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq))
  sigx22 <- (llepsi^2 * pmusig2/gisq + ewv * pdmusig2)
  sigx23 <- (0.5 * (llepsi * dmusig2/sgsq) + 0.5 * pdmusig2)
  sigx24 <- (pmusig1 - dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq)
  sigx25 <- (sigx24/sqrt(ewv/gisq) + dmusig1 * (ewv/ewu_h + S * giepsi)/ewv)
  sigx26 <- (llepsi * dmusig2/ewv + pdmusig2/sqrt(ewv/gisq))
  sigx27 <- (dmusig1 * (ewv/ewu_h + S * giepsi)/sgsq)
  sigx28 <- ((0.5 * sigx27 - 0.5 * pmusig1) * ewv/sqrt(ewv/gisq) - (0.5 * sigx14 +
    0.5 * dmusig1) * (ewv/ewu_h + S * giepsi))
  sigx29 <- ((0.5 * (llepsi * dmusig2/sgsq) - 0.5 * pmusig2) * ewv/sqrt(ewv/gisq) -
    llepsi * (0.5 * sigx16 + 0.5 * dmusig2))
  sigx30 <- (2 * (sigx5 * sigx3) - sigx6 * (1 + lambda) * sigx2)
  sigx31 <- (ewu_h * sqrt(ewv/gisq) * gisq)
  sigx32 <- (ewu_h * gisq/(ewu_h * gisq)^2)
  sigx33 <- ((pmusig1/ewu_h - sigx7 * dmusig1 * (ewv/ewu_h + S * giepsi)/sqrt(ewv/gisq))/gisq -
    0.5 * ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq^2))
  sigx34 <- (sigx33 * ewv/sqrt(ewv/gisq) + sigx7 * dmusig1 * (ewv/ewu_h + S * giepsi) +
    0.5 * (sigx14 * sigx8))
  sigx35 <- ((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/sgsq^2)))/sqrt(ewv/gisq) -
    sgvsq * gisq/sgsq^2)
  sigx36 <- (sigx35 * (ewv/ewu_h + S * giepsi) + sgvsq/ewu_h)
  sigx37 <- (((1 + lambda) * pmusig2/ewu_h - llepsi * sigx12 * dmusig2/sqrt(ewv/gisq))/gisq -
    0.5 * (llepsi * pmusig2/sgsq^2))
  sigx38 <- (sigx37 * ewv/sqrt(ewv/gisq) + llepsi * sigx12 * dmusig2 + 0.5 * (sigx16 *
    (2 * sigx10 + epsilon_isq/ewv)))
  sigx39 <- (sigx35 * llepsi + (1 + lambda) * sgvsq/ewu_h)
  sigx40 <- ((ewv/ewu_h + S * giepsi) * sigx19/(ewv * gisq) - sgvsq * giZi/sgsq^2)
  sigx41 <- (sigx40 * (ewv/ewu_h + S * giepsi) * pmusig1 + (S * pmusig1 * Zi_epsi -
    dmusig1 * (ewv/ewu_h + S * giepsi) * sigx19)/sgsq)
  sigx42 <- (0.5 * (sgvsq/sgsq^2) - 0.5/(sqrt(ewv/gisq) * gisq^2))
  sigx43 <- (llepsi * sigx42 * ewv * giZi + S * sgvsq * Zi_epsi)
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
  Xsig11 <- sweep(Xgi, MARGIN = 1, STATS = sgvsq * giZi/sgsq^2, FUN = "*")
  Xsig11bis <- sweep(Xgi, MARGIN = 1, STATS = (1/(ewu_h * gisq) - (ewv/ewu_h +
    S * giepsi)/sgsq^2), FUN = "*")
  Xsig12 <- sweep(XiZi, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig13 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx26/gisq, FUN = "*")
  Xsig14 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx16/ewv), FUN = "*")
  Xsig15 <- sweep(XiZisq, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig16 <- sweep(Xgi, MARGIN = 1, STATS = sgvsq * giZisq/sgsq^2, FUN = "*")
  Xsig17 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx19, FUN = "*")
  Xsig18 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx15, FUN = "*")
  Xsig19 <- sweep((Xsig12 - Xsig11), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig20 <- sweep((Xsig15 - Xsig16), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig21 <- sweep((Xsig17 + Xsig19), MARGIN = 1, STATS = (sigx3), FUN = "*")
  Xsig22 <- sweep((Xsig18 + Xsig20), MARGIN = 1, STATS = (sigx3), FUN = "*")
  Xsig23 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx18/sigx4, FUN = "*")
  Xsig24 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx21/sigx4, FUN = "*")
  Xsig25 <- sweep((Xsig13 - 0.5 * Xsig14), MARGIN = 1, STATS = sigx17, FUN = "*")
  Xsig26 <- sweep((Xsig13 - 0.5 * Xsig14), MARGIN = 1, STATS = sigx20, FUN = "*")
  Xsig27 <- sweep((Xsig12 - Xsig11), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig28 <- sweep((Xsig15 - Xsig16), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig29 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = 0.5 * (llepsi2/ewv),
    FUN = "*")
  Xsig30 <- sweep(Xgi, MARGIN = 1, STATS = S * pmusig2, FUN = "*")
  Xsig31 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx30 * ewu_h * gisq/(sigx4 *
    ewu_h * gisq)^2, FUN = "*")
  Xsig32 <- sweep((Xsig30 - Xsig29), MARGIN = 1, STATS = 1/(sigx4 * ewu_h * gisq),
    FUN = "*")
  Xsig33 <- sweep(Xgi, MARGIN = 1, STATS = S * (2 * sigx10 + epsilon_isq/ewv) *
    dmusig2/sgsq, FUN = "*")
  Xsig34 <- sweep(Xgi, MARGIN = 1, STATS = (S * ((1 + lambda)/(ewu_h * gisq) -
    llepsi/sgsq^2)), FUN = "*")
  Xsig35 <- sweep(Xgi, MARGIN = 1, STATS = S * (sigx7 * (ewv/ewu_h + S * giepsi)/(ewv *
    gisq) + 0.5/sgsq^2) * dmusig1 * ewv/sqrt(ewv/gisq), FUN = "*")
  Xsig36 <- sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*")
  Xsig37 <- sweep((2 * (S * Xsig11bis) + Xsig36), MARGIN = 1, STATS = pmusig1,
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = S * sigx8 * dmusig1/sgsq, FUN = "*")
  Xsig38 <- sweep(Xgi, MARGIN = 1, STATS = S * (0.5 * sigx27 + 0.5 * sigx24), FUN = "*")
  Xsig39 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx9/ewv), FUN = "*")
  Xsig40 <- sweep((2 * Xsig34 + Xsig36), MARGIN = 1, STATS = pmusig2, FUN = "*")
  Xsig41 <- sweep(Xgi, MARGIN = 1, STATS = S * (llepsi * sigx12/(ewv * gisq) +
    0.5/sgsq^2) * dmusig2 * ewv/sqrt(ewv/gisq), FUN = "*")
  Xsig42 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = ((0.5 * sigx11 - sigx12 *
    dmusig2 * ewv/sqrt(ewv/gisq))/ewv), FUN = "*")
  Xsig43 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = sigx13/sigx4, FUN = "*")
  Xsig44 <- sweep((0.5 * (Xsig40 - Xsig33) + Xsig41 - 0.5 * Xsig42), MARGIN = 1,
    STATS = sigx2, FUN = "*")
  Xsig45 <- sweep((0.5 * Xsig37 + Xsig35 - 0.5 * Xsig39), MARGIN = 1, STATS = (sigx3),
    FUN = "*")
  Xsig46 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = 0.5 * (sigx6/ewv),
    FUN = "*") + sweep(Xgi, MARGIN = 1, STATS = S * sigx23, FUN = "*")
  Xsig47 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = 0.5 * (sigx5/ewv),
    FUN = "*") + Xsig38
  Xsig48 <- sweep((Xsig7 - 2 * Xsig8), MARGIN = 1, STATS = llepsi2 * ewu_h * gisq/(sigx4 *
    ewu_h * gisq)^2, FUN = "*")
  Xsig49 <- sweep(Xsig46, MARGIN = 1, STATS = (1 + lambda) * sigx2, FUN = "*") -
    sweep(Xsig47, MARGIN = 1, STATS = 2 * (sigx3), FUN = "*")
  Xsig50 <- sweep((Xsig26 + Xsig27), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig51 <- sweep((Xsig25 + Xsig28), MARGIN = 1, STATS = sigx2, FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 3, ncol = nXvar + nuZUvar +
    nvZVvar + 3)
  hessll[1:nXvar, 1:nXvar] <- (((0.5 * ((sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar * pmusig2 * sigx2/ewv/sigx4))) - 2 * crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * (S^2/gisq) * pmusig2 * sigx2/ewv/sigx4, FUN = "*"),
    Xgi)) - crossprod(sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = wHvar *
    S * dmusig2/sgsq * sigx2/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * crossprod(sweep(((0.5 *
    Xsig4 + Xsig2)), MARGIN = 1, STATS = wHvar * sigx2/ewv/sigx4, FUN = "*"),
    (Xepsi_i - 2 * Xsig0)) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * llepsi * dmusig2/(sqrt(ewv/gisq) * gisq^2) * sigx2/ewv/sigx4, FUN = "*"),
    Xgi))) - 2 * ((0.5 * ((sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar *
    pmusig1 * sigx3/ewv/sigx4))) - 2 * crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    (S^2/gisq) * pmusig1 * sigx3/ewv/sigx4, FUN = "*"), Xgi)) - crossprod(sweep((Xepsi_i -
    2 * Xsig1), MARGIN = 1, STATS = wHvar * S * dmusig1/sgsq * sigx3/ewv/sigx4,
    FUN = "*"), Xgi)) - (crossprod(sweep(0.5 * ((0.5 * Xsig5 + Xsig6)), MARGIN = 1,
    STATS = wHvar * sigx3/ewv/sigx4, FUN = "*"), (Xepsi_i - 2 * Xsig1)) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S^2 * dmusig1 * (ewv/ewu_h + S * giepsi)/(sqrt(ewv/gisq) *
      gisq^2) * sigx3/ewv/sigx4, FUN = "*"), Xgi))))) - crossprod(sweep((Xsig7 -
    2 * Xsig8), MARGIN = 1, STATS = wHvar/sigx4/sigx4, FUN = "*"), (Xsig7 - 2 *
    Xsig8)))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xsig49, MARGIN = 1,
    STATS = wHvar/(sigx4 * ewu_h * gisq), FUN = "*") - sweep(Xsig31, MARGIN = 1,
    STATS = wHvar, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep((2 *
    (Xsig45) - (Xsig43 + Xsig44)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(-(Xsig32 -
    Xsig48), MARGIN = 1, STATS = wHvar * sigx2, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep((2 * Xsig21 -
    (Xsig24 + Xsig50)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 3)] <- colSums(sweep((2 * Xsig22 -
    (Xsig23 + Xsig51)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.25 * sigx27 - 0.5 * (0.5 * sigx27 -
      0.5 * pmusig1)) * ewv - 0.5 * (sigx5 * (ewv/ewu_h + S * giepsi)/gisq)) *
      sigx3) - ((0.25 * (llepsi * dmusig2/sgsq) - 0.5 * (0.5 * (llepsi * dmusig2/sgsq) -
      0.5 * pmusig2)) * ewv - 0.5 * (llepsi * sigx6/gisq)) * (1 + lambda)^2 *
      sigx2)/(sigx4 * ewu_h^2 * gisq) - (sigx30/gisq + 0.5 * (sigx4 * ewu_h)) *
      sigx30 * gisq/(sigx4 * ewu_h * gisq)^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (2 * ((0.5 *
    (0.5 * (sigx8 * dmusig1 * ewv/sigx31) + 2 * (((0.25 * (ewv/(ewu_h * sgsq^2)) -
      0.5 * sigx32) * (ewv/ewu_h + S * giepsi) - 0.5 * (sigx7 * ewv/ewu_h)) *
      pmusig1)) - (((0.25 * (ewv/sgsq^2) + 0.5 * (sigx7 * (ewv/ewu_h + S *
    giepsi)/gisq))/ewu_h - 0.5 * sigx32) * dmusig1 * ewv/sqrt(ewv/gisq) + 0.5 *
    (sigx9 * (ewv/ewu_h + S * giepsi)/(ewu_h * gisq)))) * sigx3) - ((0.5 * (0.5 *
    ((2 * sigx10 + epsilon_isq/ewv) * dmusig2 * ewv/sigx31) + 2 * ((llepsi *
    (0.25 * (ewv/(ewu_h * sgsq^2)) - 0.5 * sigx32) - 0.5 * (sigx12 * ewv/ewu_h)) *
    pmusig2)) - (((0.25 * (ewv/sgsq^2) + 0.5 * (llepsi * sigx12/gisq))/ewu_h -
    0.5 * sigx32) * dmusig2 * ewv/sqrt(ewv/gisq) + 0.5 * (llepsi * (0.5 * sigx11 -
    sigx12 * dmusig2 * ewv/sqrt(ewv/gisq))/(ewu_h * gisq)))) * (1 + lambda) *
    sigx2 + sigx13 * sigx30/(sigx4 * ewu_h * gisq)))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (llepsi2 * (sigx30/gisq + 0.5 * (sigx4 * ewu_h)) *
      gisq/(sigx4 * ewu_h * gisq)^2 + (0.5 * (llepsi2 * llepsi/gisq) + 0.5 *
      (ewv * pmusig2)) * (1 + lambda)/(sigx4 * ewu_h^2 * gisq)) * sigx2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((sigx28 * sigx19/gisq + 0.5 * (sigx14 *
      ewv * sgvsq * giZi/sgsq^2)) * sigx3) - ((sigx29 * sigx20/gisq + 0.5 *
      (sigx16 * ewv * sgvsq * giZi/sgsq^2)) * (1 + lambda) * sigx2 + sigx21 *
      sigx30/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 3)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((sigx28 * sigx15/gisq + 0.5 * (sigx14 *
      ewv * sgvsq * giZisq/sgsq^2)) * sigx3) - ((sigx29 * sigx17/gisq + 0.5 *
      (sigx16 * ewv * sgvsq * giZisq/sgsq^2)) * (1 + lambda) * sigx2 + sigx18 *
      sigx30/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((0.5 * (sigx9 * sigx8) + 0.5 * ((2 * ((sigx7/ewu_h - 0.5 * ((1/ewu_h -
      (ewv/ewu_h + S * giepsi) * gisq/sgsq^2) * (ewv/ewu_h + S * giepsi)/sgsq^2)) *
      ewv) - epsilon_isq/ewv) * pmusig1 - sigx7 * sigx8 * dmusig1 * ewv/sqrt(ewv/gisq)) -
      (1/(ewu_h * gisq) - ((sigx7^2 + 0.5/sgsq^2) * (ewv/ewu_h + S * giepsi) +
        0.5 * ((1/ewu_h - (ewv/ewu_h + S * giepsi) * gisq/sgsq^2) * ewv/sgsq^2) +
        0.5 * sigx7)) * dmusig1 * ewv/sqrt(ewv/gisq)) * sigx3) - ((0.5 *
      ((0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) * (2 * sigx10 +
        epsilon_isq/ewv)) + 0.5 * ((2 * ((sigx12 * (1 + lambda)/ewu_h - 0.5 *
      (llepsi * ((1 + lambda)/ewu_h - llepsi * gisq/sgsq^2)/sgsq^2)) * ewv) -
      epsilon_isq/ewv) * pmusig2 - sigx12 * (2 * sigx10 + epsilon_isq/ewv) *
      dmusig2 * ewv/sqrt(ewv/gisq)) - ((1 + lambda)/(ewu_h * gisq) - ((sigx12^2 +
      0.5/sgsq^2) * llepsi + 0.5 * (((1 + lambda)/ewu_h - llepsi * gisq/sgsq^2) *
      ewv/sgsq^2) + 0.5 * sigx12)) * dmusig2 * ewv/sqrt(ewv/gisq)) * sigx2 +
      sigx13^2/sigx4))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = -wHvar * (((((1 +
    lambda) * pmusig2/ewu_h - 0.5 * (dmusig2/sqrt(ewv/gisq))) * ewv + 0.5 * (llepsi2 *
    (2 * sigx10 + epsilon_isq/ewv)))/(sigx4 * ewu_h * gisq) - llepsi2 * sigx13 *
    ewu_h * gisq/(sigx4 * ewu_h * gisq)^2) * sigx2), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (2 * ((sigx34 *
    sigx19 - (sigx36 * giZi + 0.5 * (S * Zi_epsi/sqrt(ewv/gisq))) * sigx14 *
    ewv/sgsq^2) * sigx3) - ((sigx38 * sigx20 - (sigx39 * giZi + 0.5 * (S * Zi_epsi/sqrt(ewv/gisq))) *
    sigx16 * ewv/sgsq^2) * sigx2 + sigx21 * sigx13/sigx4))/sigx4, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 3)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (2 * ((sigx34 *
    sigx15 - (sigx36 * giZisq + 0.5 * (S * Zisq_epsi/sqrt(ewv/gisq))) * sigx14 *
    ewv/sgsq^2) * sigx3) - ((sigx38 * sigx17 - (sigx39 * giZisq + 0.5 * (S *
    Zisq_epsi/sqrt(ewv/gisq))) * sigx16 * ewv/sgsq^2) * sigx2 + sigx18 * sigx13/sigx4))/sigx4,
    FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (4/(1 + 2 * lambda)^2 - (((llepsi2 * llepsi/gisq + ewv * pmusig2)/(sigx4 *
      ewu_h^2 * gisq) + llepsi2^2 * sigx2/(sigx4 * ewu_h * gisq)^2) * sigx2 +
      1/(1 + lambda)^2)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (((llepsi2 * sigx21/sigx4 - sigx22 * sigx20/sqrt(ewv/gisq))/gisq + sigx16 *
      ewv * sgvsq * giZi/sgsq^2) * sigx2/(sigx4 * ewu_h)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 3)] <- sum(wHvar *
    (((llepsi2 * sigx18/sigx4 - sigx22 * sigx17/sqrt(ewv/gisq))/gisq + sigx16 *
      ewv * sgvsq * giZisq/sgsq^2) * sigx2/(sigx4 * ewu_h)))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    ((2 * ((sigx41 * sigx19 - ((sigx42 * ewv * (ewv/ewu_h + S * giepsi) * giZi +
      S * sgvsq * Zi_epsi) * giZi + ((ewv/ewu_h + S * giepsi) * (Zisq - 2 *
      (sqrt(ewv/gisq) * sgvsq * gisq * giZi^2/sgsq^2)) + S * giZi * Zi_epsi) *
      sgvsq) * sigx14/sgsq^2) * sigx3) - ((((llepsi * sigx20/(ewv * gisq) -
      sgvsq * giZi/sgsq^2) * llepsi * pmusig2 + (S * pmusig2 * Zi_epsi - llepsi *
      dmusig2 * sigx20)/sgsq) * sigx20 - (sigx43 * giZi + (llepsi * (Zisq -
      2 * (sqrt(ewv/gisq) * sgvsq * gisq * giZi^2/sgsq^2)) + S * giZi * Zi_epsi) *
      sgvsq) * sigx16/sgsq^2) * sigx2 + sigx21^2/sigx4))/sigx4 - 0.5 * ((Zisq -
      giZi^2/gisq)/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 3)] <- sum(wHvar *
    ((2 * ((sigx41 * sigx15 - ((sigx42 * ewv * (ewv/ewu_h + S * giepsi) * giZi +
      S * sgvsq * Zi_epsi) * giZisq + ((ewv/ewu_h + S * giepsi) * (Zicub -
      2 * (sqrt(ewv/gisq) * sgvsq * gisq * giZi * giZisq/sgsq^2)) + S * giZi *
      Zisq_epsi) * sgvsq) * sigx14/sgsq^2) * sigx3) - ((((llepsi * sigx20/(ewv *
      gisq) - sgvsq * giZi/sgsq^2) * llepsi * pmusig2 + (S * pmusig2 * Zi_epsi -
      llepsi * dmusig2 * sigx20)/sgsq) * sigx17 - (sigx43 * giZisq + (llepsi *
      (Zicub - 2 * (sqrt(ewv/gisq) * sgvsq * gisq * giZi * giZisq/sgsq^2)) +
      S * giZi * Zisq_epsi) * sgvsq) * sigx16/sgsq^2) * sigx2 + sigx21 * sigx18/sigx4))/sigx4 -
      0.5 * ((Zicub - giZi * giZisq/gisq)/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 3), (nXvar + nuZUvar + nvZVvar + 3)] <- sum(wHvar *
    ((2 * (((((ewv/ewu_h + S * giepsi) * sigx15/(ewv * gisq) - sgvsq * giZisq/sgsq^2) *
      (ewv/ewu_h + S * giepsi) * pmusig1 + (S * pmusig1 * Zisq_epsi - dmusig1 *
      (ewv/ewu_h + S * giepsi) * sigx15)/sgsq) * sigx15 - ((sigx42 * ewv *
      (ewv/ewu_h + S * giepsi) * giZisq + S * sgvsq * Zisq_epsi) * giZisq +
      ((ewv/ewu_h + S * giepsi) * (Zifour - 2 * (sqrt(ewv/gisq) * sgvsq * gisq *
        giZisq^2/sgsq^2)) + S * giZisq * Zisq_epsi) * sgvsq) * sigx14/sgsq^2) *
      sigx3) - ((((llepsi * sigx17/(ewv * gisq) - sgvsq * giZisq/sgsq^2) *
      llepsi * pmusig2 + (S * pmusig2 * Zisq_epsi - llepsi * dmusig2 * sigx17)/sgsq) *
      sigx17 - ((llepsi * sigx42 * ewv * giZisq + S * sgvsq * Zisq_epsi) *
      giZisq + (llepsi * (Zifour - 2 * (sqrt(ewv/gisq) * sgvsq * gisq * giZisq^2/sgsq^2)) +
      S * giZisq * Zisq_epsi) * sgvsq) * sigx16/sgsq^2) * sigx2 + sigx18^2/sigx4))/sigx4 -
      0.5 * ((Zifour - giZisq^2/gisq)/gisq)))
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
tslnormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttslnorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    initTsl <- start_st$initTsl
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ptslnormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel Modified BC92 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(ptslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptslnormlike_mbc92,
    grad = pgradtslnormlike_mbc92, hess = phesstslnormlike_mbc92, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(ptslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(ptslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(ptslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradtslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesstslnormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtslnormlike_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phesstslnormlike_mbc92(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstslnormlike_mbc92(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptslnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtslnormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
