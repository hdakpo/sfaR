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
# Convolution: generalized exponential - normal                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for generalized exponential-normal distribution
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
pgenexponormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  mustar2 <- -(2 * exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  ll <- log(2) - 1/2 * Wu - (TT - 1)/2 * Wv - 1/2 * log(gisq) - (TT - 1)/2 * log(2 *
    pi) + log(exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar1/sigmastar)^2)) * pnorm(mustar1/sigmastar) -
    exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar2/sigmastar)^2)) * pnorm(mustar2/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for generalized exponential-normal distribution
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
pstgenexponorm_k90 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, S, wHvar, printInfo, tol, whichStart, initIter, initAlg) {
  if (whichStart == 1L) {
    Esti <- cstgenexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGenExpo <- NULL
  } else {
    cat("Initialization: SFA + generalized-exponential-normal distribution...\n")
    initGenExpo <- maxLik::maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradgenexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initGenExpo$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "eta1", "eta2")
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for generalized exponential-normal distribution
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
pgradgenexponormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
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
  sgsq <- (sqrt(ewv/gisq) * gisq)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- ((2 * (ewv/ewu_h) + S * giepsi)/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (sigx3 * pmusig1 - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * ((ewv/ewu_h + S * giepsi) *
    pmusig1))
  sigx6 <- (dmusig2 * ewv/sqrt(ewv/gisq) - (2 * (ewv/ewu_h) + S * giepsi) * pmusig2)
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx12 <- (2/(ewu_h * gisq) - 0.5 * ((2 * (ewv/ewu_h) + S * giepsi)/sgsq^2))
  sigx10 <- ((2 * (ewv/ewu_h) + S * giepsi) * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx13 <- (sigx9 * sigx3 - (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) *
    sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * gid1_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx16 <- ((2 * (ewv/ewu_h) + S * giepsi) * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * gid1_epsi_gisq/sgsq - (2 * (ewv/ewu_h) + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx18 <- (sigx14 * sigx3 * sigx15 - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * gid1sq_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx20 <- (S * gid1sq_epsi_gisq/sgsq - (2 * (ewv/ewu_h) + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx21 <- (sigx14 * sigx3 * sigx19 - sigx16 * sigx2 * sigx20)
  Xsig0 <- sweep(Xgi, MARGIN = 1, STATS = (S * (2 * (ewv/ewu_h) + S * giepsi)/gisq),
    FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * (ewv/ewu_h + S * giepsi)/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sgsq, FUN = "*")
  Xsig4 <- sweep(Xepsi_i - 2 * Xsig0, MARGIN = 1, STATS = (pmusig2/ewv), FUN = "*")
  Xsig5 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sgsq, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig4 + Xsig2), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1, STATS = sigx3, FUN = "*")
  gradll <- cbind(sweep(Xsig7 - Xsig8, MARGIN = 1, STATS = 1/sigx4, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ((sigx5 * sigx3 - sigx6 * sigx2)/(sigx4 *
      ewu_h * gisq) - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx13/sigx4 -
      0.5 * (TT - 1)), FUN = "*"), sigx18/sigx4 - 0.5 * (gid2_gicub/gisq),
    sigx21/sigx4 - 0.5 * (gid2sq_gicub/gisq))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for generalized exponential-normal distribution
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
phessgenexponormlike_k90 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar) {
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
  sgsq <- (sqrt(ewv/gisq) * gisq)
  musig1 <- ((ewv/ewu_h + S * giepsi)/sgsq)
  musig2 <- ((2 * (ewv/ewu_h) + S * giepsi)/sgsq)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  sigx2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx3 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  sigx4 <- (sigx3 * pmusig1 - sigx2 * pmusig2)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * ((ewv/ewu_h + S * giepsi) *
    pmusig1))
  sigx6 <- (dmusig2 * ewv/sqrt(ewv/gisq) - (2 * (ewv/ewu_h) + S * giepsi) * pmusig2)
  sigx7 <- (1/(ewu_h * gisq) - 0.5 * ((ewv/ewu_h + S * giepsi)/sgsq^2))
  sigx8 <- (2 * (sigx7 * (ewv/ewu_h + S * giepsi)) + epsilon_isq/ewv)
  sigx9 <- (0.5 * (sigx8 * pmusig1) - sigx7 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx12 <- (2/(ewu_h * gisq) - 0.5 * ((2 * (ewv/ewu_h) + S * giepsi)/sgsq^2))
  sigx10 <- ((2 * (ewv/ewu_h) + S * giepsi) * sigx12)
  sigx11 <- ((2 * sigx10 + epsilon_isq/ewv) * pmusig2)
  sigx13 <- (sigx9 * sigx3 - (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) *
    sigx2)
  sigx14 <- ((ewv/ewu_h + S * giepsi) * pmusig1/sgsq - dmusig1)
  sigx15 <- (S * gid1_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx16 <- ((2 * (ewv/ewu_h) + S * giepsi) * pmusig2/sgsq - dmusig2)
  sigx17 <- (S * gid1_epsi_gisq/sgsq - (2 * (ewv/ewu_h) + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2_gicub/sgsq^2)
  sigx18 <- (sigx14 * sigx3 * sigx15 - sigx16 * sigx2 * sigx17)
  sigx19 <- (S * gid1sq_epsi_gisq/sgsq - (ewv/ewu_h + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx20 <- (S * gid1sq_epsi_gisq/sgsq - (2 * (ewv/ewu_h) + S * giepsi) * (sqrt(ewv/gisq) -
    0.5 * (ewv/sgsq)) * gid2sq_gicub/sgsq^2)
  sigx21 <- (sigx14 * sigx3 * sigx19 - sigx16 * sigx2 * sigx20)
  sigx22 <- (ewv/ewu_h + S * giepsi)
  ddgq <- (dmusig1 * sigx22/sgsq)
  sigx23 <- ((pmusig1 - dmusig1 * sigx22/sgsq)/sqrt(ewv/gisq) + dmusig1 * sigx22/ewv)
  sigx24 <- (sqrt(ewv/gisq) - 0.5 * (ewv/sgsq))
  sigx25 <- (2 * sigx10 + epsilon_isq/ewv)
  sigx26 <- (2 * (ewv/ewu_h) + S * giepsi)
  sigx27 <- (sqrt(ewv/gisq) * gisq^2)
  sigx28 <- (sigx26 * dmusig2/ewv + (pmusig2 - sigx26 * dmusig2/sgsq)/sqrt(ewv/gisq))
  sigx29 <- (0.5 * sigx11 - sigx12 * dmusig2 * ewv/sqrt(ewv/gisq))
  sigx30 <- ((0.5 * ddgq - 0.5 * pmusig1) * ewv/sqrt(ewv/gisq) - (0.5 * sigx14 +
    0.5 * dmusig1) * sigx22)
  sigx31 <- ((sigx26 * dmusig2/sgsq - pmusig2) * ewv - sigx26^2 * pmusig2/gisq)
  sigx32 <- (((pmusig1/ewu_h - sigx7 * dmusig1 * sigx22/sqrt(ewv/gisq))/gisq -
    0.5 * (sigx22 * pmusig1/sgsq^2)) * ewv/sqrt(ewv/gisq) + sigx7 * dmusig1 *
    sigx22 + 0.5 * (sigx14 * sigx8))
  sigx33 <- (((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/sgsq^2)))/sqrt(ewv/gisq) -
    sigx24 * gisq/sgsq^2) * sigx22 + sigx24/ewu_h)
  sigx34 <- (((2 * (pmusig2/ewu_h) - sigx26 * sigx12 * dmusig2/sqrt(ewv/gisq))/gisq -
    0.5 * (sigx26 * pmusig2/sgsq^2)) * ewv/sqrt(ewv/gisq) + sigx26 * sigx12 *
    dmusig2 + 0.5 * (sigx16 * sigx25))
  sigx35 <- (((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/sgsq^2)))/sqrt(ewv/gisq) -
    sigx24 * gisq/sgsq^2) * sigx26 + 2 * (sigx24/ewu_h))
  sigx36 <- ((sigx22 * sigx15/(ewv * gisq) - sigx24 * gid2_gicub/sgsq^2) * sigx22 *
    pmusig1 + (S * pmusig1 * gid1_epsi_gisq - dmusig1 * sigx22 * sigx15)/sgsq)
  sigx43 <- (0.5 * (sigx24/sgsq^2) - 0.5/sigx27)
  sigx37 <- (sigx43 * ewv * sigx22 * gid2_gicub + S * sigx24 * gid1_epsi_gisq)
  sigx38 <- (sqrt(ewv/gisq) * sigx24 * gid2_gicub * gid2sq_gicub * gisq/sgsq^2)
  sigx39 <- (sqrt(ewv/gisq) * sigx24 * gid2_gicub^2 * gisq/sgsq^2)
  sigx40 <- ((sigx26 * sigx17/(ewv * gisq) - sigx24 * gid2_gicub/sgsq^2) * sigx26 *
    pmusig2 + (S * pmusig2 * gid1_epsi_gisq - sigx26 * dmusig2 * sigx17)/sgsq)
  sigx41 <- (sigx43 * sigx26 * ewv * gid2_gicub + S * sigx24 * gid1_epsi_gisq)
  sigx42 <- (gd1Zfour_gd1_epsi_gisq/sgsq - sigx24 * gid2sq_gicub * gid1sq_epsi_gisq/sgsq^2)
  sigx44 <- (sqrt(ewv/gisq) * sigx24 * gid2sq_gicub^2 * gisq/sgsq^2)
  Xsig0 <- sweep(Xgi, MARGIN = 1, STATS = (S * (2 * (ewv/ewu_h) + S * giepsi)/gisq),
    FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * (ewv/ewu_h + S * giepsi)/gisq),
    FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sgsq, FUN = "*")
  Xsig4 <- sweep(Xepsi_i - 2 * Xsig0, MARGIN = 1, STATS = (pmusig2/ewv), FUN = "*")
  Xsig5 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sgsq, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig4 + Xsig2), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig5 + Xsig6), MARGIN = 1, STATS = sigx3, FUN = "*")
  Xsig9 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx23/gisq, FUN = "*")
  Xsig10 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx14/ewv), FUN = "*")
  Xsig11 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx15, FUN = "*")
  Xsig12 <- sweep((Xsig9 - 0.5 * Xsig10), MARGIN = 1, STATS = sigx19, FUN = "*")
  Xsig13 <- sweep(Xi_gid1_gisq, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig14 <- sweep(Xgi, MARGIN = 1, STATS = sigx24 * gid2_gicub/sgsq^2, FUN = "*")
  Xsig15 <- sweep((Xsig13 - Xsig14), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig16 <- sweep((Xsig13 - Xsig14), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig17 <- sweep(((Xepsi_i - 2 * Xsig1)), MARGIN = 1, STATS = sigx9/ewv, FUN = "*")
  Xsig18 <- sweep(Xgi, MARGIN = 1, STATS = (sigx7 * sigx22/(ewv * gisq) + 0.5/sgsq^2) *
    dmusig1 * ewv/sqrt(ewv/gisq), FUN = "*")
  Xsig19 <- sweep(Xgi, MARGIN = 1, STATS = (S * (1/(ewu_h * gisq) - sigx22/sgsq^2)),
    FUN = "*")
  Xsig20 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx6/ewv), FUN = "*")
  Xsig21 <- sweep(Xgi, MARGIN = 1, STATS = S * pmusig2, FUN = "*")
  Xsig22 <- sweep(Xgi, MARGIN = 1, STATS = (0.5 * ddgq + 0.5 * (pmusig1 - dmusig1 *
    sigx22/sgsq)), FUN = "*")
  Xsig23 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx28/gisq, FUN = "*")
  Xsig24 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx16/ewv), FUN = "*")
  Xsig25 <- sweep((Xsig23 - 0.5 * Xsig24), MARGIN = 1, STATS = sigx17, FUN = "*")
  Xsig26 <- sweep((Xsig23 - 0.5 * Xsig24), MARGIN = 1, STATS = sigx20, FUN = "*")
  Xsig27 <- sweep((Xsig25 + Xsig16), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig28 <- sweep((Xsig7 - Xsig8), MARGIN = 1, STATS = sigx21/sigx4, FUN = "*")
  Xsig29 <- sweep((Xsig7 - Xsig8), MARGIN = 1, STATS = sigx18/sigx4, FUN = "*")
  Xsig30 <- sweep((Xsig11 + Xsig15), MARGIN = 1, STATS = sigx3, FUN = "*")
  Xsig31 <- sweep(Xi_gid1sq_gisq, MARGIN = 1, STATS = 1/sgsq, FUN = "*")
  Xsig32 <- sweep(Xgi, MARGIN = 1, STATS = sigx24 * gid2sq_gicub/sgsq^2, FUN = "*")
  Xsig33 <- sweep((Xsig31 - Xsig32), MARGIN = 1, STATS = S * sigx14, FUN = "*")
  Xsig34 <- sweep((Xsig31 - Xsig32), MARGIN = 1, STATS = S * sigx16, FUN = "*")
  Xsig35 <- sweep((Xsig26 + Xsig34), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig36 <- sweep((Xsig12 + Xsig33), MARGIN = 1, STATS = sigx3, FUN = "*")
  Xsig37 <- sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*")
  Xsig38 <- sweep((2 * Xsig19 + Xsig37), MARGIN = 1, STATS = pmusig1, FUN = "*")
  Xsig39 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx8 * dmusig1/sgsq, FUN = "*")
  Xsig40 <- sweep((Xsig7 - Xsig8), MARGIN = 1, STATS = sigx13/sigx4, FUN = "*")
  Xsig41 <- sweep(Xgi, MARGIN = 1, STATS = (S * (2/(ewu_h * gisq) - sigx26/sgsq^2)),
    FUN = "*")
  Xsig42 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx25 * dmusig2/sgsq, FUN = "*")
  Xsig43 <- sweep((2 * Xsig41 + Xsig37), MARGIN = 1, STATS = pmusig2, FUN = "*")
  Xsig44 <- sweep(Xgi, MARGIN = 1, STATS = (sigx26 * sigx12/(ewv * gisq) + 0.5/sgsq^2) *
    dmusig2 * ewv/sqrt(ewv/gisq), FUN = "*")
  Xsig45 <- sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = (sigx29/ewv), FUN = "*")
  Xsig46 <- sweep((0.5 * (Xsig43 - Xsig42) + S * Xsig44 - 0.5 * Xsig45), MARGIN = 1,
    STATS = sigx2, FUN = "*")
  Xsig47 <- sweep((0.5 * (Xsig38 - Xsig39) + S * Xsig18 - 0.5 * Xsig17), MARGIN = 1,
    STATS = sigx3, FUN = "*")
  Xsig48 <- sweep((Xsig7 - Xsig8), MARGIN = 1, STATS = (sigx5 * sigx3 - sigx6 *
    sigx2) * ewu_h * gisq/(sigx4 * ewu_h * gisq)^2, FUN = "*")
  Xsig49 <- sweep((Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx5/ewv), FUN = "*")
  Xsig50 <- sweep((0.5 * Xsig20 + Xsig21), MARGIN = 1, STATS = sigx2, FUN = "*")
  Xsig51 <- sweep((0.5 * Xsig49 + S * Xsig22), MARGIN = 1, STATS = sigx3, FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar + nuZUvar +
    nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar * pmusig2 * sigx2/ewv/sigx4))) - crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * pmusig2 * 2 * (S^2/gisq) * sigx2/ewv/sigx4, FUN = "*"), Xgi) -
    crossprod(sweep((Xepsi_i - 2 * Xsig0), MARGIN = 1, STATS = wHvar * S * dmusig2/sgsq *
      sigx2/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * (crossprod(sweep((0.5 * Xsig4 +
    Xsig2), MARGIN = 1, STATS = wHvar * sigx2/ewv/sigx4, FUN = "*"), (Xepsi_i -
    2 * Xsig0))) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S^2 * sigx26 *
    dmusig2/sigx27 * sigx2/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * (sapply(1:nXvar,
    function(x) crossprod(Xsq[[x]], as.matrix(wHvar * pmusig1 * sigx3/ewv/sigx4))) -
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * pmusig1 * 2 * (S^2/gisq) *
      sigx3/ewv/sigx4, FUN = "*"), Xgi) - crossprod(sweep((Xepsi_i - 2 * Xsig1),
    MARGIN = 1, STATS = wHvar * S * dmusig1/sgsq * sigx3/ewv/sigx4, FUN = "*"),
    Xgi)) - (0.5 * crossprod(sweep(((0.5 * Xsig5 + Xsig6)), MARGIN = 1, STATS = wHvar *
    sigx3/ewv/sigx4, FUN = "*"), (Xepsi_i - 2 * Xsig1)) + crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S^2 * dmusig1 * sigx22/sigx27 * sigx3/ewv/sigx4,
    FUN = "*"), Xgi))) - crossprod(sweep((Xsig7 - Xsig8), MARGIN = 1, STATS = wHvar/sigx4^2,
    FUN = "*"), (Xsig7 - Xsig8))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep((Xsig50 - Xsig51),
    MARGIN = 1, STATS = wHvar/(sigx4 * ewu_h * gisq), FUN = "*") - sweep(Xsig48,
    MARGIN = 1, STATS = wHvar, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep((Xsig47 -
    (Xsig40 + Xsig46)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep((Xsig30 - (Xsig29 +
    Xsig27)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 2] <- colSums(sweep((Xsig36 - (Xsig28 +
    Xsig35)), MARGIN = 1, STATS = wHvar/sigx4, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.25 * ddgq - 0.5 * (0.5 * ddgq - 0.5 * pmusig1)) *
      ewv - 0.5 * (sigx5 * sigx22/gisq)) * sigx3 - sigx2 * (ewv * pmusig2 -
      sigx26 * sigx6/gisq))/(sigx4 * ewu_h^2 * gisq) - ((sigx5 * sigx3 - sigx6 *
      sigx2)/gisq + 0.5 * (sigx4 * ewu_h)) * (sigx5 * sigx3 - sigx6 * sigx2) *
      gisq/(sigx4 * ewu_h * gisq)^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * ((0.5 * (0.5 *
    (sigx8 * dmusig1 * ewv/(ewu_h * sqrt(ewv/gisq) * gisq)) + 2 * (((0.25 * (ewv/(ewu_h *
    sgsq^2)) - 0.5 * (ewu_h * gisq/(ewu_h * gisq)^2)) * sigx22 - 0.5 * (sigx7 *
    ewv/ewu_h)) * pmusig1)) - (((0.25 * (ewv/sgsq^2) + 0.5 * (sigx7 * sigx22/gisq))/ewu_h -
    0.5 * (ewu_h * gisq/(ewu_h * gisq)^2)) * dmusig1 * ewv/sqrt(ewv/gisq) + 0.5 *
    (sigx9 * sigx22/(ewu_h * gisq)))) * sigx3 - (sigx13 * (sigx5 * sigx3 - sigx6 *
    sigx2)/(sigx4 * ewu_h * gisq) + (0.5 * (sigx25 * dmusig2 * ewv/(ewu_h * sqrt(ewv/gisq) *
    gisq) + 2 * (((0.5 * (ewv/(ewu_h * sgsq^2)) - ewu_h * gisq/(ewu_h * gisq)^2) *
    sigx26 - sigx12 * ewv/ewu_h) * pmusig2)) - (((sigx26 * sigx12/gisq + 0.5 *
    (ewv/sgsq^2))/ewu_h - ewu_h * gisq/(ewu_h * gisq)^2) * dmusig2 * ewv/sqrt(ewv/gisq) +
    sigx29 * sigx26/(ewu_h * gisq))) * sigx2))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sigx30 * sigx15/gisq + 0.5 * (sigx14 * ewv *
      sigx24 * gid2_gicub/sgsq^2)) * sigx3 - ((sigx31 * sigx17/sgsq + sigx16 *
      ewv * sigx24 * gid2_gicub/sgsq^2) * sigx2 + sigx18 * (sigx5 * sigx3 -
      sigx6 * sigx2)/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sigx30 * sigx19/gisq + 0.5 * (sigx14 * ewv *
      sigx24 * gid2sq_gicub/sgsq^2)) * sigx3 - ((sigx31 * sigx20/sgsq + sigx16 *
      ewv * sigx24 * gid2sq_gicub/sgsq^2) * sigx2 + sigx21 * (sigx5 * sigx3 -
      sigx6 * sigx2)/(sigx4 * gisq)))/(sigx4 * ewu_h), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (sigx9 * sigx8) + 0.5 * ((2 * ((sigx7/ewu_h - 0.5 * ((1/ewu_h - sigx22 *
      gisq/sgsq^2) * sigx22/sgsq^2)) * ewv) - epsilon_isq/ewv) * pmusig1 -
      sigx7 * sigx8 * dmusig1 * ewv/sqrt(ewv/gisq)) - (1/(ewu_h * gisq) - ((sigx7^2 +
      0.5/sgsq^2) * sigx22 + 0.5 * ((1/ewu_h - sigx22 * gisq/sgsq^2) * ewv/sgsq^2) +
      0.5 * sigx7)) * dmusig1 * ewv/sqrt(ewv/gisq)) * sigx3 - (sigx13^2/sigx4 +
      (0.5 * (sigx29 * sigx25) + 0.5 * ((2 * ((2 * (sigx12/ewu_h) - 0.5 * (sigx26 *
        (2/ewu_h - sigx26 * gisq/sgsq^2)/sgsq^2)) * ewv) - epsilon_isq/ewv) *
        pmusig2 - sigx25 * sigx12 * dmusig2 * ewv/sqrt(ewv/gisq)) - (2/(ewu_h *
        gisq) - ((sigx12^2 + 0.5/sgsq^2) * sigx26 + 0.5 * ((2/ewu_h - sigx26 *
        gisq/sgsq^2) * ewv/sgsq^2) + 0.5 * sigx12)) * dmusig2 * ewv/sqrt(ewv/gisq)) *
        sigx2))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((sigx32 *
    sigx15 - (sigx33 * gid2_gicub + 0.5 * (S * gid1_epsi_gisq/sqrt(ewv/gisq))) *
    sigx14 * ewv/sgsq^2) * sigx3 - ((sigx34 * sigx17 - (sigx35 * gid2_gicub +
    0.5 * (S * gid1_epsi_gisq/sqrt(ewv/gisq))) * sigx16 * ewv/sgsq^2) * sigx2 +
    sigx18 * sigx13/sigx4))/sigx4, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * ((sigx32 *
    sigx19 - (sigx33 * gid2sq_gicub + 0.5 * (S * gid1sq_epsi_gisq/sqrt(ewv/gisq))) *
    sigx14 * ewv/sgsq^2) * sigx3 - ((sigx34 * sigx20 - (sigx35 * gid2sq_gicub +
    0.5 * (S * gid1sq_epsi_gisq/sqrt(ewv/gisq))) * sigx16 * ewv/sgsq^2) * sigx2 +
    sigx21 * sigx13/sigx4))/sigx4, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    (((sigx36 * sigx15 + sigx14 * (S * (gd1Zsq_gd1_epsi_gisq/sgsq - sigx24 *
      gid2_gicub * gid1_epsi_gisq/sgsq^2) - (sigx37 * gid2_gicub + sigx22 *
      sigx24 * (gd2Zsq_gd1_gicub - 2 * sigx39))/sgsq^2)) * sigx3 - ((sigx40 *
      sigx17 + sigx16 * (S * (gd1Zsq_gd1_epsi_gisq/sgsq - sigx24 * gid2_gicub *
      gid1_epsi_gisq/sgsq^2) - (sigx41 * gid2_gicub + sigx26 * sigx24 * (gd2Zsq_gd1_gicub -
      2 * sigx39))/sgsq^2)) * sigx2 + sigx18^2/sigx4))/sigx4 - 0.5 * ((gd2Zsq_gd1_gicub -
      gid2_gicub^2/gisq)/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (((sigx36 * sigx19 + sigx14 * (S * (gd1Zcub_gd1_epsi_gisq/sgsq - sigx24 *
      gid2_gicub * gid1sq_epsi_gisq/sgsq^2) - (sigx37 * gid2sq_gicub + sigx22 *
      sigx24 * (gd2Zcub_gd1_gicub - 2 * sigx38))/sgsq^2)) * sigx3 - ((sigx40 *
      sigx20 + sigx16 * (S * (gd1Zcub_gd1_epsi_gisq/sgsq - sigx24 * gid2_gicub *
      gid1sq_epsi_gisq/sgsq^2) - (sigx41 * gid2sq_gicub + sigx26 * sigx24 *
      (gd2Zcub_gd1_gicub - 2 * sigx38))/sgsq^2)) * sigx2 + sigx18 * sigx21/sigx4))/sigx4 -
      0.5 * ((gd2Zcub_gd1_gicub - gid2_gicub * gid2sq_gicub/gisq)/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    (((((sigx22 * sigx19/(ewv * gisq) - sigx24 * gid2sq_gicub/sgsq^2) * sigx22 *
      pmusig1 + (S * pmusig1 * gid1sq_epsi_gisq - dmusig1 * sigx22 * sigx19)/sgsq) *
      sigx19 + sigx14 * (S * sigx42 - ((sigx43 * ewv * sigx22 * gid2sq_gicub +
      S * sigx24 * gid1sq_epsi_gisq) * gid2sq_gicub + sigx22 * sigx24 * (gd2Zfour_gd1_gicub -
      2 * sigx44))/sgsq^2)) * sigx3 - ((((sigx26 * sigx20/(ewv * gisq) - sigx24 *
      gid2sq_gicub/sgsq^2) * sigx26 * pmusig2 + (S * pmusig2 * gid1sq_epsi_gisq -
      sigx26 * dmusig2 * sigx20)/sgsq) * sigx20 + sigx16 * (S * sigx42 - ((sigx43 *
      sigx26 * ewv * gid2sq_gicub + S * sigx24 * gid1sq_epsi_gisq) * gid2sq_gicub +
      sigx26 * sigx24 * (gd2Zfour_gd1_gicub - 2 * sigx44))/sgsq^2)) * sigx2 +
      sigx21^2/sigx4))/sigx4 - 0.5 * ((gd2Zfour_gd1_gicub - gid2sq_gicub^2/gisq)/gisq)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for generalized exponential-normal distribution
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
genexponormAlgOpt_k90 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstgenexponorm_k90(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    InitGenExpo <- start_st$initGenExpo
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pgenexponormlike_k90(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(pgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pgenexponormlike_k90,
    grad = pgradgenexponormlike_k90, hess = phessgenexponormlike_k90, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(pgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessgenexponormlike_k90(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradgenexponormlike_k90(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessgenexponormlike_k90(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessgenexponormlike_k90(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- pgenexponormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradgenexponormlike_k90(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitGenExpo = InitGenExpo))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for generalized exponential-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
pgenexponormeff_k90 <- function(object, level) {
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
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- 1/unlist(lapply(TT, FUN = function(x) 1 + exp(eta1 * seq(1:x) + eta2 *
    (seq(1:x))^2)))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar1 <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  mustar2 <- -(2 * exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
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
