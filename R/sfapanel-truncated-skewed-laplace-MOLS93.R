################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Modified Lee and Schmidt 1993                        #
#                          - g(zit) = exp(-eta_t * (t - T)): g(zit) = 1 for T  #
# Convolution: truncated skewed laplace - normal                               #
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
ptslnormlike_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar, ngZGvar, gHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  eta[ngZGvar] <- 1
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
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
psttslnorm_mols93 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, ngZGvar, gHvar, S, wHvar, initIter, initAlg, whichStart, printInfo,
  tol) {
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
  }, Esti[nXvar + 3], rep(0.001, ngZGvar - 1))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "P", paste0("eta_", colnames(gHvar[, -ngZGvar])))
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
pgradtslnormlike_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar, gHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  eta[ngZGvar] <- 0
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it, FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmastar <- (sqrt(ewv/gisq) * gisq)
  mustar1 <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- (mustar1/sigmastar)
  pmusig1 <- pnorm(-musig1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  mustar2 <- (ewv/ewu_h + S * giepsi)
  musig2 <- (mustar2/sigmastar)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  expo1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  expo2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx1 <- (dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- (0.5 * sigx1 - 0.5 * (mustar2 * pmusig2))
  sigx3 <- (mustar1 * pmusig1 - dmusig1 * ewv/sqrt(ewv/gisq))
  sigx4 <- (2 * (expo2 * pmusig2) - expo1 * pmusig1)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * (mustar1 * pmusig1))
  sigx6 <- (2 * (sigx2 * expo2) - sigx5 * (1 + lambda) * expo1)
  sigx7 <- (sigx4 * ewu_h * gisq)
  sigx8 <- (1/(ewu_h * gisq) - 0.5 * (mustar2/sigmastar^2))
  sigx9 <- (2 * (sigx8 * mustar2) + epsilon_isq/ewv)
  sigx10 <- (0.5 * (sigx9 * pmusig2) - sigx8 * dmusig2 * ewv/sqrt(ewv/gisq))
  sigx13 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (mustar1/sigmastar^2))
  sigx11 <- (mustar1 * sigx13)
  sigx12 <- ((2 * sigx11 + epsilon_isq/ewv) * pmusig1)
  sigx14 <- (0.5 * sigx12 - sigx13 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx15 <- (2 * (sigx10 * expo2) - sigx14 * expo1)
  sigx16 <- (mustar2 * pmusig2/sigmastar - dmusig2)
  sigx17 <- (sqrt(ewv/gisq) - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = mustar2 * sigx17/sigmastar^2, FUN = "*")
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = mustar1 * sigx17/sigmastar^2, FUN = "*")
  sigx19 <- (mustar1 * pmusig1/sigmastar - dmusig1)
  sigZ3 <- sweep(sigZ1, MARGIN = 1, STATS = 2 * (sigx16 * expo2), FUN = "*") -
    sweep(sigZ2, MARGIN = 1, STATS = sigx19 * expo1, FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar1/gisq), FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar2/gisq), FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i - 2 * Xsig1, MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig4 <- sweep((2 * Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv),
    FUN = "*")
  Xsig5 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sigmastar, FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sigmastar, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1, STATS = expo1, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig4 + Xsig6), MARGIN = 1, STATS = expo2, FUN = "*")
  gradll <- cbind(sweep((Xsig7 - 2 * (Xsig8)), MARGIN = 1, STATS = 1/sigx4, FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx6/sigx7 - 0.5), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = (sigx15/sigx4 - 0.5 * (TT - 1)), FUN = "*"), 1/(1 +
      lambda) - (sigx3 * expo1/sigx7 + 2/(1 + 2 * lambda)), sweep(sigZ3, MARGIN = 1,
      STATS = 1/sigx4, FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = 0.5/gisq,
      FUN = "*"))
  return(sweep(gradll[, -(nXvar + nuZUvar + nvZVvar + ngZGvar + 1)], MARGIN = 1,
    STATS = wHvar, FUN = "*"))
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
phesstslnormlike_mols93 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, pindex, TT, S, wHvar, ngZGvar, gHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  eta[ngZGvar] <- 0
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it, FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  Xzigi <- list()
  for (i in 1:ngZGvar) {
    Xzigi[[i]] <- apply(sweep(-Xvar, MARGIN = 1, STATS = gHvar[, i] * git, FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  Zisqgiepsi <- list()
  for (i in 1:ngZGvar) {
    Zisqgiepsi[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = gHvar[, i] * git *
      epsilon_it, FUN = "*"), 2, function(x) tapply(x, pindex[, 1], sum))
  }
  Zisqgisq <- list()
  for (i in 1:ngZGvar) {
    Zisqgisq[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = 4 * gHvar[, i] *
      git^2, FUN = "*"), 2, function(x) tapply(x, pindex[, 1], sum))
  }
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmastar <- (sqrt(ewv/gisq) * gisq)
  mustar1 <- ((1 + lambda) * ewv/ewu_h + S * giepsi)
  musig1 <- (mustar1/sigmastar)
  pmusig1 <- pnorm(-musig1)
  dmusig1 <- dnorm(-musig1, 0, 1)
  mustar2 <- (ewv/ewu_h + S * giepsi)
  musig2 <- (mustar2/sigmastar)
  pmusig2 <- pnorm(-musig2)
  dmusig2 <- dnorm(-musig2, 0, 1)
  expo1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig1)^2)))
  expo2 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig2)^2)))
  sigx1 <- (dmusig2 * ewv/sqrt(ewv/gisq))
  sigx2 <- (0.5 * sigx1 - 0.5 * (mustar2 * pmusig2))
  sigx3 <- (mustar1 * pmusig1 - dmusig1 * ewv/sqrt(ewv/gisq))
  sigx4 <- (2 * (expo2 * pmusig2) - expo1 * pmusig1)
  sigx5 <- (0.5 * (dmusig1 * ewv/sqrt(ewv/gisq)) - 0.5 * (mustar1 * pmusig1))
  sigx6 <- (2 * (sigx2 * expo2) - sigx5 * (1 + lambda) * expo1)
  sigx7 <- (sigx4 * ewu_h * gisq)
  sigx8 <- (1/(ewu_h * gisq) - 0.5 * (mustar2/sigmastar^2))
  sigx9 <- (2 * (sigx8 * mustar2) + epsilon_isq/ewv)
  sigx10 <- (0.5 * (sigx9 * pmusig2) - sigx8 * dmusig2 * ewv/sqrt(ewv/gisq))
  sigx13 <- ((1 + lambda)/(ewu_h * gisq) - 0.5 * (mustar1/sigmastar^2))
  sigx11 <- (mustar1 * sigx13)
  sigx12 <- ((2 * sigx11 + epsilon_isq/ewv) * pmusig1)
  sigx14 <- (0.5 * sigx12 - sigx13 * dmusig1 * ewv/sqrt(ewv/gisq))
  sigx15 <- (2 * (sigx10 * expo2) - sigx14 * expo1)
  sigx16 <- (mustar2 * pmusig2/sigmastar - dmusig2)
  sigx17 <- (sqrt(ewv/gisq) - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = mustar2 * sigx17/sigmastar^2, FUN = "*")
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(Zigisq,
    MARGIN = 1, STATS = mustar1 * sigx17/sigmastar^2, FUN = "*")
  sigx19 <- (mustar1 * pmusig1/sigmastar - dmusig1)
  sigZ3 <- sweep(sigZ1, MARGIN = 1, STATS = 2 * (sigx16 * expo2), FUN = "*") -
    sweep(sigZ2, MARGIN = 1, STATS = sigx19 * expo1, FUN = "*")
  Xsig1 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar1/gisq), FUN = "*")
  Xsig2 <- sweep(Xgi, MARGIN = 1, STATS = (S * mustar2/gisq), FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i - 2 * Xsig1, MARGIN = 1, STATS = (pmusig1/ewv), FUN = "*")
  Xsig4 <- sweep((2 * Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (pmusig2/ewv),
    FUN = "*")
  Xsig5 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig1/sigmastar, FUN = "*")
  Xsig6 <- sweep(Xgi, MARGIN = 1, STATS = S * dmusig2/sigmastar, FUN = "*")
  Xsig7 <- sweep((0.5 * Xsig3 + Xsig5), MARGIN = 1, STATS = expo1, FUN = "*")
  Xsig8 <- sweep((0.5 * Xsig4 + Xsig6), MARGIN = 1, STATS = expo2, FUN = "*")
  wvgi <- sqrt(ewv/gisq)
  wugi <- (ewu_h * gisq)
  sigx22 <- (pmusig2 - dmusig2 * mustar2/sigmastar)
  sigx23 <- (pmusig1 - mustar1 * dmusig1/sigmastar)
  sigx24 <- (mustar1 * dmusig1/sigmastar)
  sigx25 <- (0.5 * sigx24 - 0.5 * pmusig1)
  lpwu <- (1 + lambda) * pmusig1/ewu_h
  sigx26 <- (lpwu - mustar1 * sigx13 * dmusig1/wvgi)
  sigx27 <- (dmusig2 * mustar2/sigmastar)
  sigwu <- (sigx4 * ewu_h)
  sigx28 <- (sigx4 * ewu_h^2 * gisq)
  sigx29 <- (0.5 * sigx27 - 0.5 * pmusig2)
  sigx30 <- (ewu_h * wvgi * gisq)
  sigx31 <- (ewv/(ewu_h * sigmastar^2))
  sigx34 <- (ewu_h * gisq/wugi^2)
  sigx32 <- (0.25 * sigx31 - 0.5 * sigx34)
  sigx33 <- (ewv/sigmastar^2)
  sigx35 <- (2 * sigx11 + epsilon_isq/ewv)
  sigx36 <- (1/ewu_h - mustar2 * gisq/sigmastar^2)
  sigx37 <- ((1 + lambda)/ewu_h - mustar1 * gisq/sigmastar^2)
  sigx38 <- (pmusig2/ewu_h - sigx8 * dmusig2 * mustar2/wvgi)
  sigZ4 <- sweep(Zigiepsi, MARGIN = 1, STATS = (S/wvgi), FUN = "*")
  sigx40 <- (1/gisq - 0.5 * sigx33)
  sigx41 <- (0.5/gisq - 0.5 * sigx40)
  sigx42 <- (sigx41/wvgi - sigx17 * gisq/sigmastar^2)
  sigZ5 <- sweep(Zigiepsi, MARGIN = 1, STATS = S * pmusig2, FUN = "*") - sweep(sigZ1,
    MARGIN = 1, STATS = dmusig2 * mustar2, FUN = "*")
  sigx45 <- (0.5 * (sigx17/sigmastar^2) - 0.5/(wvgi * gisq^2))
  sigZ6 <- sweep(Zigiepsi, MARGIN = 1, STATS = S * pmusig1, FUN = "*") - sweep(sigZ2,
    MARGIN = 1, STATS = mustar1 * dmusig1, FUN = "*")
  Xsig9 <- sweep((2 * Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx3/ewv), FUN = "*")
  Xsig10 <- sweep((Xsig7 - 2 * (Xsig8)), MARGIN = 1, STATS = sigx3 * ewu_h * gisq/sigx7^2,
    FUN = "*")
  Xsig11 <- sweep(Xgi, MARGIN = 1, STATS = S * (sigx22/wvgi + dmusig2 * mustar2/ewv)/gisq,
    FUN = "*")
  Xsig12 <- sweep((2 * Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (sigx16/ewv),
    FUN = "*")
  ZZ1 <- list()
  for (i in 1:ngZGvar) {
    ZZ1[[i]] <- sweep(Xgi, MARGIN = 1, STATS = sigx17 * Zigisq[, i]/sigmastar^2,
      FUN = "*")
  }
  ZZ2 <- list()
  for (i in 1:ngZGvar) {
    ZZ2[[i]] <- sweep(Xzigi[[i]], MARGIN = 1, STATS = 1/sigmastar, FUN = "*")
  }
  Xsig15 <- sweep(Xgi, MARGIN = 1, STATS = S * (mustar1 * dmusig1/ewv + sigx23/wvgi)/gisq,
    FUN = "*")
  Xsig16 <- sweep((2 * Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx19/ewv),
    FUN = "*")
  ZZ3 <- list()
  for (i in 1:ngZGvar) {
    ZZ3[[i]] <- sweep((ZZ2[[i]] - ZZ1[[i]]), MARGIN = 1, STATS = S * sigx19,
      FUN = "*")
  }
  Xsig18 <- sweep(Xgi, MARGIN = 1, STATS = S * pmusig1, FUN = "*")
  Xsig19 <- sweep(Xgi, MARGIN = 1, STATS = (S * (1/wugi - mustar2/sigmastar^2)),
    FUN = "*")
  Xsig20 <- sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*")
  Xsig21 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx9 * dmusig2/sigmastar, FUN = "*")
  Xsig22 <- sweep(Xgi, MARGIN = 1, STATS = S * (sigx8 * mustar2/(ewv * gisq) +
    0.5/sigmastar^2) * dmusig2 * ewv/wvgi, FUN = "*")
  Xsig23 <- sweep((2 * Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (sigx10/ewv),
    FUN = "*")
  Xsig24 <- sweep((Xsig7 - 2 * (Xsig8)), MARGIN = 1, STATS = sigx15/sigx4, FUN = "*")
  Xsig25 <- sweep(Xgi, MARGIN = 1, STATS = (S * ((1 + lambda)/wugi - mustar1/sigmastar^2)),
    FUN = "*")
  Xsig26 <- sweep(Xgi, MARGIN = 1, STATS = S * sigx35 * dmusig1/sigmastar, FUN = "*")
  Xsig27 <- sweep(Xgi, MARGIN = 1, STATS = S * (mustar1 * sigx13/(ewv * gisq) +
    0.5/sigmastar^2) * dmusig1 * ewv/wvgi, FUN = "*")
  Xsig28 <- sweep((2 * Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx14/ewv),
    FUN = "*")
  Xsig29 <- sweep((2 * Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = (sigx5/ewv), FUN = "*")
  Xsig30 <- sweep(Xgi, MARGIN = 1, STATS = S * (0.5 * sigx24 + 0.5 * sigx23), FUN = "*")
  Xsig31 <- sweep((2 * Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = (sigx2/ewv), FUN = "*")
  Xsig32 <- sweep(Xgi, MARGIN = 1, STATS = S * (0.5 * sigx27 + 0.5 * sigx22), FUN = "*")
  Xsig33 <- sweep((Xsig7 - 2 * (Xsig8)), MARGIN = 1, STATS = sigx6 * ewu_h * gisq/sigx7^2,
    FUN = "*")
  Xsig34 <- sweep((0.5 * Xsig29 + Xsig30), MARGIN = 1, STATS = (1 + lambda) * expo1,
    FUN = "*")
  Xsig35 <- sweep((0.5 * Xsig31 + Xsig32), MARGIN = 1, STATS = (expo2), FUN = "*")
  Xsig36 <- sweep((Xsig34 - 2 * Xsig35), MARGIN = 1, STATS = 1/sigx7, FUN = "*")
  Xsig37 <- sweep((2 * Xsig19 + 2 * Xsig20), MARGIN = 1, STATS = pmusig2, FUN = "*")
  Xsig38 <- sweep((0.5 * (Xsig37 - Xsig21) + Xsig22 - 0.5 * Xsig23), MARGIN = 1,
    STATS = 2 * (expo2), FUN = "*")
  Xsig39 <- sweep((2 * Xsig25 + 2 * Xsig20), MARGIN = 1, STATS = pmusig1, FUN = "*")
  Xsig40 <- sweep((0.5 * (Xsig39 - Xsig26) + Xsig27 - 0.5 * Xsig28), MARGIN = 1,
    STATS = expo1, FUN = "*")
  Xsig41 <- sweep((Xsig38 - (Xsig24 + Xsig40)), MARGIN = 1, STATS = 1/sigx4, FUN = "*")
  Xsig42 <- sweep((Xsig18 - 0.5 * Xsig9), MARGIN = 1, STATS = 1/sigx7, FUN = "*")
  Xsig43 <- sweep((Xsig42 - Xsig10), MARGIN = 1, STATS = expo1, FUN = "*")
  ZZ4 <- list()
  for (i in 1:ngZGvar) {
    ZZ4[[i]] <- sweep((Xsig11 - 0.5 * Xsig12), MARGIN = 1, STATS = sigZ1[, i],
      FUN = "*")
  }
  ZZ5 <- list()
  for (i in 1:ngZGvar) {
    ZZ5[[i]] <- sweep((ZZ2[[i]] - ZZ1[[i]]), MARGIN = 1, STATS = S * sigx16,
      FUN = "*")
  }
  ZZ6 <- list()
  for (i in 1:ngZGvar) {
    ZZ6[[i]] <- sweep((ZZ4[[i]] + ZZ5[[i]]), MARGIN = 1, STATS = 2 * (expo2),
      FUN = "*")
  }
  ZZ7 <- list()
  for (i in 1:ngZGvar) {
    ZZ7[[i]] <- sweep((Xsig7 - 2 * (Xsig8)), MARGIN = 1, STATS = sigZ3[, i]/sigx4,
      FUN = "*")
  }
  ZZ8 <- list()
  for (i in 1:ngZGvar) {
    ZZ8[[i]] <- sweep((Xsig15 - 0.5 * Xsig16), MARGIN = 1, STATS = sigZ2[, i],
      FUN = "*")
  }
  ZZ9 <- list()
  for (i in 1:ngZGvar) {
    ZZ9[[i]] <- sweep((ZZ8[[i]] + ZZ3[[i]]), MARGIN = 1, STATS = expo1, FUN = "*")
  }
  ZZ10 <- list()
  for (i in 1:ngZGvar) {
    ZZ10[[i]] <- sweep((ZZ6[[i]] - (ZZ7[[i]] + ZZ9[[i]])), MARGIN = 1, STATS = wHvar/sigx4,
      FUN = "*")
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + ngZGvar + 1)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(2 *
    Xsq[[x]], as.matrix(wHvar * pmusig1 * expo1/ewv/sigx4))) - crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * 2 * (S^2/gisq) * pmusig1 * expo1/ewv/sigx4, FUN = "*"),
    Xgi) - crossprod(sweep((2 * Xepsi_i - 2 * Xsig1), MARGIN = 1, STATS = wHvar *
    S * dmusig1/sigmastar * expo1/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 *
    Xsig3 + Xsig5), MARGIN = 1, STATS = wHvar * expo1/ewv/sigx4, FUN = "*"),
    (2 * Xepsi_i - 2 * Xsig1)) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * mustar1 * dmusig1/(wvgi * gisq^2) * expo1/ewv/sigx4, FUN = "*"), Xgi)) -
    2 * (0.5 * (sapply(1:nXvar, function(x) crossprod(2 * Xsq[[x]], as.matrix(wHvar *
      pmusig2 * expo2/ewv/sigx4))) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
      2 * (S^2/gisq) * pmusig2 * expo2/ewv/sigx4, FUN = "*"), Xgi) - crossprod(sweep((2 *
      Xepsi_i - 2 * Xsig2), MARGIN = 1, STATS = wHvar * S * dmusig2/sigmastar *
      expo2/ewv/sigx4, FUN = "*"), Xgi)) - (0.5 * crossprod(sweep((0.5 * Xsig4 +
      Xsig6), MARGIN = 1, STATS = wHvar * expo2/ewv/sigx4, FUN = "*"), (2 *
      Xepsi_i - 2 * Xsig2)) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
      S^2 * mustar2 * dmusig2/(wvgi * gisq^2) * expo2/ewv/sigx4, FUN = "*"),
      Xgi))) - crossprod(sweep(Xsig7 - 2 * (Xsig8), MARGIN = 1, STATS = wHvar/sigx4^2,
    FUN = "*"), Xsig7 - 2 * (Xsig8))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod((Xsig36 - Xsig33),
    sweep(uHvar, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(Xsig41,
    sweep(vHvar, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(-(sweep(Xsig43, MARGIN = 1,
    STATS = wHvar, FUN = "*")))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar +
    ngZGvar + 1)] <- sapply(ZZ10, colSums)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.25 * sigx27 - 0.5 * sigx29) * ewv -
      0.5 * (sigx2 * mustar2/gisq)) * expo2) - ((0.25 * sigx24 - 0.5 * sigx25) *
      ewv - 0.5 * (mustar1 * sigx5/gisq)) * (1 + lambda)^2 * expo1)/sigx28 -
      (sigx6/gisq + 0.5 * sigwu) * sigx6 * gisq/sigx7^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (2 * ((0.5 *
    (0.5 * (sigx9 * dmusig2 * ewv/sigx30) + 2 * ((sigx32 * mustar2 - 0.5 * (sigx8 *
      ewv/ewu_h)) * pmusig2)) - (((0.25 * sigx33 + 0.5 * (sigx8 * mustar2/gisq))/ewu_h -
    0.5 * sigx34) * dmusig2 * ewv/wvgi + 0.5 * (sigx10 * mustar2/wugi))) * expo2) -
    ((0.5 * (0.5 * (sigx35 * dmusig1 * ewv/sigx30) + 2 * ((mustar1 * sigx32 -
      0.5 * (sigx13 * ewv/ewu_h)) * pmusig1)) - (((0.25 * sigx33 + 0.5 * (mustar1 *
      sigx13/gisq))/ewu_h - 0.5 * sigx34) * dmusig1 * ewv/wvgi + 0.5 * (mustar1 *
      sigx14/wugi))) * (1 + lambda) * expo1 + sigx15 * sigx6/sigx7))/sigx4,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 1] <- (wHvar *
    ((sigx3 * (sigx6/gisq + 0.5 * sigwu) * gisq/sigx7^2 + (0.5 * (sigx3 * mustar1/gisq) +
      0.5 * (ewv * pmusig1)) * (1 + lambda)/sigx28) * expo1)) %*% uHvar
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + ngZGvar + 1)] <- 2 * (crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * (sigx29 * ewv/wvgi - (0.5 * sigx16 + 0.5 * dmusig2) * mustar2)/gisq *
      expo2/sigwu, FUN = "*"), sigZ1) + crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * 0.5 * (sigx16 * ewv * sigx17/sigmastar^2) * expo2/sigwu,
    FUN = "*"), Zigisq)) - (crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (sigx25 * ewv/wvgi - mustar1 * (0.5 * sigx19 + 0.5 * dmusig1))/gisq * (1 +
    lambda) * expo1/sigwu, FUN = "*"), sigZ2) + crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * 0.5 * (sigx19 * ewv * sigx17/sigmastar^2) * (1 + lambda) *
      expo1/sigwu, FUN = "*"), Zigisq) + crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * sigx6/(sigx4 * gisq)/sigwu, FUN = "*"), sigZ3))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((0.5 * (sigx10 * sigx9) + 0.5 * ((2 * ((sigx8/ewu_h - 0.5 * (sigx36 *
      mustar2/sigmastar^2)) * ewv) - epsilon_isq/ewv) * pmusig2 - sigx8 * sigx9 *
      dmusig2 * ewv/wvgi) - (1/wugi - ((sigx8^2 + 0.5/sigmastar^2) * mustar2 +
      0.5 * (sigx36 * ewv/sigmastar^2) + 0.5 * sigx8)) * dmusig2 * ewv/wvgi) *
      expo2) - ((0.5 * (sigx14 * sigx35) + 0.5 * ((2 * ((sigx13 * (1 + lambda)/ewu_h -
      0.5 * (mustar1 * sigx37/sigmastar^2)) * ewv) - epsilon_isq/ewv) * pmusig1 -
      sigx13 * sigx35 * dmusig1 * ewv/wvgi) - ((1 + lambda)/wugi - ((sigx13^2 +
      0.5/sigmastar^2) * mustar1 + 0.5 * (sigx37 * ewv/sigmastar^2) + 0.5 *
      sigx13)) * dmusig1 * ewv/wvgi) * expo1 + sigx15^2/sigx4))/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), nXvar + nuZUvar + nvZVvar +
    1] <- (wHvar * (-((((lpwu - 0.5 * (dmusig1/wvgi)) * ewv + 0.5 * (sigx3 *
    sigx35))/sigx7 - sigx3 * sigx15 * ewu_h * gisq/sigx7^2) * expo1))) %*% vHvar
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar + 1)] <- 2 * (crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx38/gisq - 0.5 * (mustar2 * pmusig2/sigmastar^2)) *
      ewv/wvgi + sigx8 * dmusig2 * mustar2 + 0.5 * (sigx16 * sigx9)) * expo2/sigx4,
    FUN = "*"), sigZ1) - (crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx42 * mustar2 + sigx17/ewu_h) * sigx16 * ewv/sigmastar^2 * expo2/sigx4,
    FUN = "*"), Zigisq) + crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    0.5 * sigx16 * ewv/sigmastar^2 * expo2/sigx4, FUN = "*"), sigZ4))) - (crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx26/gisq - 0.5 * (mustar1 * pmusig1/sigmastar^2)) *
      ewv/wvgi + mustar1 * sigx13 * dmusig1 + 0.5 * (sigx19 * sigx35)) * expo1/sigx4,
    FUN = "*"), sigZ2) - (crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (sigx42 * mustar1 + (1 + lambda) * sigx17/ewu_h) * sigx19 * ewv/sigmastar^2 *
    expo1/sigx4, FUN = "*"), Zigisq) + crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    0.5 * sigx19 * ewv/sigmastar^2 * expo1/sigx4, FUN = "*"), sigZ4)) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx15/sigx4/sigx4, FUN = "*"), sigZ3))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 1] <- sum(wHvar *
    (4/(1 + 2 * lambda)^2 - (((sigx3 * mustar1/gisq + ewv * pmusig1)/sigx28 +
      sigx3^2 * expo1/sigx7^2) * expo1 + 1/(1 + lambda)^2)))
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar + nvZVvar + 2):(nXvar +
    nuZUvar + nvZVvar + ngZGvar + 1)] <- colSums(sweep(sigZ3, MARGIN = 1, STATS = wHvar *
    sigx3/sigx4/gisq * expo1/sigwu, FUN = "*") - sweep(sigZ2, MARGIN = 1, STATS = wHvar *
    (mustar1^2 * pmusig1/gisq + ewv * sigx23)/wvgi/gisq * expo1/sigwu, FUN = "*") +
    sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx19 * ewv * sigx17/sigmastar^2 *
      expo1/sigwu, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar +
    1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + ngZGvar +
    1)] <- 2 * (crossprod(sweep(sigZ1, MARGIN = 1, STATS = wHvar * mustar2/(ewv *
    gisq) * mustar2 * pmusig2 * expo2/sigx4, FUN = "*"), sigZ1) - crossprod(sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * sigx17/sigmastar^2 * mustar2 * pmusig2 * expo2/sigx4,
    FUN = "*"), sigZ1) + crossprod(sweep(sigZ5, MARGIN = 1, STATS = wHvar/sigmastar *
    expo2/sigx4, FUN = "*"), sigZ1) + S * (sapply(1:ngZGvar, function(x) crossprod(Zisqgiepsi[[x]],
    as.matrix(wHvar * sigx16/sigmastar * expo2/sigx4))) - crossprod(sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * sigx16 * sigx17/sigmastar^2 * expo2/sigx4, FUN = "*"),
    Zigiepsi)) - (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx16 *
    sigx45 * ewv * mustar2/sigmastar^2 * expo2/sigx4, FUN = "*"), Zigisq) + crossprod(sweep(Zigiepsi,
    MARGIN = 1, STATS = wHvar * sigx16 * S * sigx17/sigmastar^2 * expo2/sigx4,
    FUN = "*"), Zigisq) + (sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
    as.matrix(wHvar * sigx16 * mustar2 * sigx17/sigmastar^2 * expo2/sigx4))) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx16 * 2 * (wvgi *
      sigx17 * gisq/sigmastar^2) * mustar2 * sigx17/sigmastar^2 * expo2/sigx4,
      FUN = "*"), Zigisq)))) - (crossprod(sweep(sigZ2, MARGIN = 1, STATS = wHvar *
    mustar1/(ewv * gisq) * mustar1 * pmusig1 * expo1/sigx4, FUN = "*"), sigZ2) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx17/sigmastar^2 *
      mustar1 * pmusig1 * expo1/sigx4, FUN = "*"), sigZ2) + crossprod(sweep(sigZ6,
    MARGIN = 1, STATS = wHvar/sigmastar * expo1/sigx4, FUN = "*"), sigZ2) + S *
    (sapply(1:ngZGvar, function(x) crossprod(Zisqgiepsi[[x]], as.matrix(wHvar *
      sigx19 * expo1/sigmastar/sigx4))) - crossprod(sweep(Zigisq, MARGIN = 1,
      STATS = wHvar * sigx19 * expo1 * sigx17/sigmastar^2/sigx4, FUN = "*"),
      Zigiepsi)) - (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx19 *
    expo1 * mustar1 * sigx45 * ewv/sigmastar^2/sigx4, FUN = "*"), Zigisq) + crossprod(sweep(Zigiepsi,
    MARGIN = 1, STATS = wHvar * sigx19 * expo1 * S * sigx17/sigmastar^2/sigx4,
    FUN = "*"), Zigisq) + (sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
    as.matrix(wHvar * sigx19 * expo1 * mustar1 * sigx17/sigmastar^2/sigx4))) -
    crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx19 * expo1 * mustar1 *
      sigx17 * 2 * (wvgi * sigx17 * gisq/sigmastar^2/sigmastar^2/sigx4), FUN = "*"),
      Zigisq))) + crossprod(sweep(sigZ3, MARGIN = 1, STATS = wHvar/sigx4/sigx4,
    FUN = "*"), sigZ3)) - 0.5 * (sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
    as.matrix(wHvar/gisq))) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar/gisq/gisq,
    FUN = "*"), Zigisq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll[-(nXvar + nuZUvar + nvZVvar + ngZGvar + 1), -(nXvar + nuZUvar +
    nvZVvar + ngZGvar + 1)])
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
tslnormAlgOpt_mols93 <- function(start, randStart, sdStart, olsParam, dataTable,
  S, nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, ngZGvar,
  Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax, stepmax,
  tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttslnorm_mols93(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, ngZGvar = ngZGvar,
      initIter = initIter, initAlg = initAlg, gHvar = gHvar, whichStart = whichStart,
      tol = tol, printInfo = printInfo)
    InitTsl <- start_st$initTsl
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ptslnormlike_mols93(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel MOLS93 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(ptslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptslnormlike_mols93,
    grad = pgradtslnormlike_mols93, hess = phesstslnormlike_mols93, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) {
    -sum(ptslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
    prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(ptslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(ptslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(ptslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradtslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phesstslnormlike_mols93(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtslnormlike_mols93(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
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
      mleObj$hessian <- phesstslnormlike_mols93(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstslnormlike_mols93(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptslnormlike_mols93(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtslnormlike_mols93(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitTsl = InitTsl))
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
ptslnormeff_mols93 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar + 1]
  eta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 2):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$ngZGvar)]
  eta[object$ngZGvar] <- 1
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  gHvar <- object$gHvar
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
  epsilon_it <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
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
