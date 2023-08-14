################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Inefficiency structure: u_it = g(zit)u_i                                     #
#                         Battese and Coelli 1992 specification:               #
#                          - g(zit) = exp(-eta * (t - T))                      #
#                         Cuesta and Orea (2002), Feng and Serletis (2009)     #
#                          - g(zit) = exp(-eta1 * (t - T) - eta2 * (t - T)^2)  #
#                         Alvarez, Amsler, Orea, Schmidt (2006)                #
#                          - g(zit) = exp(eta * gHvar)                         #
#                         Kumbhakar and Wang 2005 specification:               #
#                          - g(zit) = exp(eta * (t - t1))                      #
#                         Cuesta 2000 specification:                           #
#                          - g(zit) = exp(-eta_i * (t - T))                    #
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
puninormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, ngZGvar, gHvar, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
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
#' @param modelType specification of inefficiency model G(t)u_i
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
pstuninorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, ngZGvar, gHvar, wHvar, modelType, S, printInfo, tol, whichStart,
  initIter, initAlg) {
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
  }, if (modelType %in% c("bc92a", "kw05")) {
    0.001
  } else {
    if (modelType == "bc92b") {
      c(0.001, 0.001)
    } else {
      if (modelType == "bc92c") {
        rep(0, ngZGvar)
      } else {
        if (modelType == "c00") {
          rep(0.001, ngZGvar)
        }
      }
    }
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), if (modelType %in% c("bc92a", "kw05")) {
    "eta"
  } else {
    if (modelType == "bc92b") {
      c("eta1", "eta2")
    } else {
      if (modelType == "bc92c") {
        paste0("Zg_", colnames(gHvar))
      } else {
        if (modelType == "c00") {
          paste0("eta_", colnames(gHvar))
        }
      }
    }
  })
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgraduninormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, ngZGvar, gHvar, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it, FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sqgv <- sqrt(ewv/gisq)
  ggsq <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(ggsq/sqgv, 0, 1)
  dmusig2 <- dnorm(S * giepsi/(sqgv * gisq), 0, 1)
  pmusig1 <- pnorm(ggsq/sqgv)
  pmusig2 <- pnorm(S * giepsi/(sqgv * gisq))
  sigx1 <- ((pmusig1 - pmusig2) * sqgv)
  sigx2 <- (dmusig1 - dmusig2)
  sigx3 <- ((epsilon_isq - (-(S * giepsi/gisq))^2 * gisq)/ewv)
  sigx4 <- (S * dmusig2 * ewv * giepsi/(sqgv * gisq)^2)
  sigx5 <- (0.5 * sigx4 - 0.5 * (ggsq * dmusig1))
  sigx6 <- (sqrt(12)/2 * (dmusig1 * ewu_h/sigx1) - 0.5)
  sigx7 <- ((pmusig1 - pmusig2) * sqgv * gisq)
  sigx8 <- (sqgv - 0.5 * (ewv/(sqgv * gisq)))
  sigZ1 <- sweep(Zigisq, MARGIN = 1, STATS = giepsi/gisq, FUN = "*")
  sigZ2 <- sweep(Zigisq, MARGIN = 1, STATS = 0.5 * (ggsq), FUN = "*")
  sigZ3 <- sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = dmusig1/(sqgv *
    gisq), FUN = "*")
  sigZ4 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/(sqgv * gisq), FUN = "*")
  sigZ6 <- sweep(Zigisq, MARGIN = 1, STATS = sigx8 * giepsi/(sqgv * gisq)^2, FUN = "*")
  sigZ7 <- sweep((sigZ4 - sigZ6), MARGIN = 1, STATS = S * dmusig2, FUN = "*")
  sigZ8 <- sweep((sigZ3 - sigZ7), MARGIN = 1, STATS = 1/(pmusig1 - pmusig2), FUN = "*")
  sigZ9 <- sweep(Zigisq, MARGIN = 1, STATS = (-(S * giepsi/gisq))^2, FUN = "*")
  sigZ10 <- sweep((Zigiepsi - sigZ1), MARGIN = 1, STATS = (S^2 * giepsi/gisq),
    FUN = "*")
  sigZ11 <- sweep((sigZ9 + 2 * sigZ10), MARGIN = 1, STATS = 1/ewv, FUN = "*")
  sigZ12 <- sweep(Zigisq, MARGIN = 1, STATS = 0.5/gisq, FUN = "*")
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * sigx2/sigx7, FUN = "*") -
    0.5 * (sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*") - sweep(Xgi,
      MARGIN = 1, STATS = 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*")), sweep(uHvar,
    MARGIN = 1, STATS = sigx6, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx5/sigx1 +
    0.5 + 0.5 * sigx3 - 0.5 * TT), FUN = "*"), sigZ8 + 0.5 * sigZ11 - sigZ12)
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
#' @param gHvar matrix of inefficiency determinants
#' @param ngZGvar number of variables explaining inefficiency
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phessuninormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, ngZGvar, gHvar, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1], sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
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
    Xsq[[i]] <- apply(sweep(2 * Xvar, MARGIN = 1, STATS = Xvar[, i], FUN = "*"),
      2, function(x) tapply(x, pindex[, 1], sum))
  }
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sqgv <- sqrt(ewv/gisq)
  ggsq <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(ggsq/sqgv, 0, 1)
  dmusig2 <- dnorm(S * giepsi/(sqgv * gisq), 0, 1)
  pmusig1 <- pnorm(ggsq/sqgv)
  pmusig2 <- pnorm(S * giepsi/(sqgv * gisq))
  sigx1 <- ((pmusig1 - pmusig2) * sqgv)
  sigx2 <- (dmusig1 - dmusig2)
  sigx3 <- ((epsilon_isq - (-(S * giepsi/gisq))^2 * gisq)/ewv)
  sigx4 <- (S * dmusig2 * ewv * giepsi/(sqgv * gisq)^2)
  sigx5 <- (0.5 * sigx4 - 0.5 * (ggsq * dmusig1))
  sigx6 <- (sqrt(12)/2 * (dmusig1 * ewu_h/sigx1) - 0.5)
  sigx7 <- ((pmusig1 - pmusig2) * sqgv * gisq)
  sigx8 <- (sqgv - 0.5 * (ewv/(sqgv * gisq)))
  sigx9 <- (0.5 * (ggsq^2/ewv) - 0.5 * (ewv/(sqgv * gisq)^2))
  sigx10 <- (0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/(sqgv * gisq)^2)))
  sigZ1 <- sweep(Zigisq, MARGIN = 1, STATS = giepsi/gisq, FUN = "*")
  sigZ2 <- sweep(Zigisq, MARGIN = 1, STATS = 0.5 * (ggsq), FUN = "*")
  sigZ3 <- sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = dmusig1/(sqgv *
    gisq), FUN = "*")
  sigZ4 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/(sqgv * gisq), FUN = "*")
  sigZ6 <- sweep(Zigisq, MARGIN = 1, STATS = sigx8 * giepsi/(sqgv * gisq)^2, FUN = "*")
  sigZ7 <- sweep((sigZ4 - sigZ6), MARGIN = 1, STATS = S * dmusig2, FUN = "*")
  sigZ8 <- sweep((sigZ3 - sigZ7), MARGIN = 1, STATS = 1/(pmusig1 - pmusig2), FUN = "*")
  sigZ9 <- sweep(Zigisq, MARGIN = 1, STATS = (-(S * giepsi/gisq))^2, FUN = "*")
  sigZ10 <- sweep((Zigiepsi - sigZ1), MARGIN = 1, STATS = (S^2 * giepsi/gisq),
    FUN = "*")
  sigZ11 <- sweep((sigZ9 + 2 * sigZ10), MARGIN = 1, STATS = 1/ewv, FUN = "*")
  sigZ12 <- sweep(Zigisq, MARGIN = 1, STATS = 0.5/gisq, FUN = "*")
  sigZ13 <- sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = ggsq/ewv,
    FUN = "*")
  sigZ14 <- sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = (ggsq *
    gisq/ewv), FUN = "*")
  sigZ15 <- sweep((sqrt(12)/4 * Zigisq - sqrt(12)/2 * sigZ14), MARGIN = 1, STATS = 1/gisq,
    FUN = "*")
  sigZ16 <- sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = sigx9 *
    dmusig1, FUN = "*")
  sigZ17 <- sweep((sigZ3 - sigZ7), MARGIN = 1, STATS = sigx5/(pmusig1 - pmusig2),
    FUN = "*")
  sigZ18 <- sweep((sigZ4 - sigZ6), MARGIN = 1, STATS = (S^2 * giepsi^2), FUN = "*")
  sigZ19 <- sweep(Zigisq, MARGIN = 1, STATS = (sigx10/sqgv - sigx8 * gisq/(sqgv *
    gisq)^2) * giepsi, FUN = "*")
  sigZ20 <- sweep(Zigiepsi, MARGIN = 1, STATS = 1/sqgv, FUN = "*")
  sigZ21 <- sweep((sigZ19 + 0.5 * sigZ20), MARGIN = 1, STATS = ewv, FUN = "*")
  sigZ22 <- sweep((0.5 * sigZ18 - sigZ21), MARGIN = 1, STATS = S * dmusig2/(sqgv *
    gisq)^2, FUN = "*")
  sigZ23 <- sweep((sigZ16 - sigZ17), MARGIN = 1, STATS = 1/sqgv, FUN = "*")
  sigZ24 <- sweep(Zigisq, MARGIN = 1, STATS = 2 * (sqgv * sigx8 * gisq * giepsi/(sqgv *
    gisq)^2), FUN = "*")
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar, ncol = nXvar + nuZUvar +
    nvZVvar + ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    ((S * dmusig2 * giepsi/gisq - ggsq * dmusig1)/(ewv * (pmusig1 - pmusig2) *
      sqgv * gisq) - sigx2^2/((pmusig1 - pmusig2) * sqgv * gisq)^2), FUN = "*"),
    Xgi) - 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar/ewv))) -
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 * (S^2/gisq)/ewv, FUN = "*"),
      Xgi))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = -wHvar * (sqrt(12)/2 * (S * (ggsq/(ewv * (pmusig1 - pmusig2) * sqgv) +
      sigx2/(sigx1^2 * gisq)) * dmusig1 * ewu_h)), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(0.5 *
    (sweep(Xepsi_i, MARGIN = 1, STATS = wHvar/ewv, FUN = "*") - sweep(Xgi, MARGIN = 1,
      STATS = wHvar * 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*")) + sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (dmusig2 * (ewv - S^2 * giepsi^2/gisq)/(sqgv *
      gisq)^2) - 0.5 * ((1/gisq - ggsq^2/ewv) * dmusig1))/sigx1 - sigx5 * sigx2/(sigx1^2 *
      gisq)), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    ngZGvar)] <- S * ((((sapply(1:ngZGvar, function(x) {
    crossprod(Xzigi[[x]], as.matrix(wHvar * dmusig1/(sqgv * gisq)/(pmusig1 -
      pmusig2)))
  }) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * dmusig1/(sqgv * gisq)/(pmusig1 -
    pmusig2), FUN = "*"), (sigZ13 + sigZ12))) - crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * sigx2/(pmusig1 - pmusig2)/(sqgv * gisq)/(pmusig1 - pmusig2),
    FUN = "*"), (sigZ3 - sigZ7))) - (sapply(1:ngZGvar, function(x) {
    crossprod(Xzigi[[x]], as.matrix(wHvar/(sqgv * gisq)/(pmusig1 - pmusig2) *
      dmusig2))
  }) - crossprod(Xgi, sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx8/(sqgv * gisq)^2/(pmusig1 -
    pmusig2) * dmusig2, FUN = "*") + sweep((sigZ4 - sigZ6), MARGIN = 1, STATS = wHvar *
    S^2 * giepsi/(ewv * gisq)/(pmusig1 - pmusig2) * dmusig2, FUN = "*")))) +
    0.5 * (S * (2 * ((sapply(1:ngZGvar, function(x) {
      crossprod(Xzigi[[x]], as.matrix(wHvar * giepsi/(ewv * gisq)))
    }) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar/gisq * giepsi/(ewv *
      gisq), FUN = "*"), Zigisq)) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar/(ewv *
      gisq), FUN = "*"), (Zigiepsi - sigZ1))) + crossprod(sweep(2 * (Xgi),
      MARGIN = 1, STATS = wHvar/(ewv * gisq), FUN = "*"), sigZ1))))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sqrt(12)/2 * (((0.5 - sqrt(12)/2 * (ggsq * ewu_h *
      gisq/ewv))/sigx1 - sqrt(12)/2 * (dmusig1 * ewu_h/sigx1^2)) * dmusig1 *
      ewu_h), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((0.5 *
    ((sqrt(12)/2 - sqrt(12)/2 * (ggsq^2 * gisq/ewv))/sigx1) + sqrt(12)/2 * (sigx5/sigx1^2)) *
    dmusig1 * ewu_h), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    dmusig1 * ewu_h/sigx1, FUN = "*"), (sigZ15 - sqrt(12)/2 * (sigZ8)))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (S * ((0.5 * (S^2 * giepsi^2) - ewv * gisq)/(sqgv * gisq)^2 + 1) *
      dmusig2 * ewv * giepsi/(sqgv * gisq)^2) - 0.25 * (ggsq^3 * dmusig1 *
      gisq/ewv))/sigx1 - (((0.5 * ((pmusig1 - pmusig2)/(sqgv * gisq)) + 0.5 *
      (S * dmusig2 * giepsi/(sqgv * gisq)^2)) * ewv - 0.5 * (ggsq * dmusig1)) *
      sigx5/sigx1^2 + 0.5 * ((epsilon_isq - (-(S * giepsi/gisq))^2 * gisq)/ewv))),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- crossprod(vHvar, sweep((sigZ23 -
    sigZ22), MARGIN = 1, STATS = wHvar/(pmusig1 - pmusig2), FUN = "*") - sweep(0.5 *
    sigZ11, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- (((0.5 *
    (sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * ggsq/(sqgv * gisq) * dmusig1/(pmusig1 -
        pmusig2)))
    }) + crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * S/gisq/(sqgv * gisq) *
      dmusig1/(pmusig1 - pmusig2), FUN = "*"), (Zigiepsi - sigZ1))) + S * (sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar/(sqgv * gisq) * dmusig1/(pmusig1 -
        pmusig2)))
    }) - (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar/gisq/(sqgv * gisq) *
    dmusig1/(pmusig1 - pmusig2), FUN = "*"), (Zigiepsi - sigZ1)) + sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * giepsi/gisq/(sqgv * gisq) *
        dmusig1/(pmusig1 - pmusig2)))
    }))) - crossprod(sweep((sigZ2 + S * (Zigiepsi - sigZ1)), MARGIN = 1, STATS = wHvar *
    ggsq/ewv/(sqgv * gisq) * dmusig1/(pmusig1 - pmusig2), FUN = "*"), (sigZ2 +
    S * (Zigiepsi - sigZ1)))) - crossprod(sweep((sigZ2 + S * (Zigiepsi - sigZ1)),
    MARGIN = 1, STATS = wHvar * sigx8/(sqgv * gisq)^2 * dmusig1/(pmusig1 - pmusig2),
    FUN = "*"), Zigisq)) - (crossprod(sweep((sigZ3 - sigZ7), MARGIN = 1, STATS = wHvar/(pmusig1 -
    pmusig2)/(pmusig1 - pmusig2), FUN = "*"), (sigZ3 - sigZ7)) + S * ((sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar/(sqgv * gisq) * dmusig2/(pmusig1 -
        pmusig2)))
    }) - crossprod(sweep((sigZ4 - sigZ6), MARGIN = 1, STATS = wHvar * S^2 * giepsi/(sqgv *
    gisq) * dmusig2/(pmusig1 - pmusig2), FUN = "*"), (sigZ4 - sigZ6))) - ((crossprod(sweep(Zigisq,
    MARGIN = 1, STATS = wHvar * (0.5 * (sigx8/(sqgv * gisq)^2) - 0.5/(sqgv *
      gisq^2)) * ewv * giepsi/(sqgv * gisq)^2 * dmusig2/(pmusig1 - pmusig2),
    FUN = "*"), Zigisq) + sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * sigx8 * giepsi/(sqgv * gisq)^2 *
      dmusig2/(pmusig1 - pmusig2)))
  })) + crossprod(sweep((2 * Zigiepsi - sigZ24), MARGIN = 1, STATS = wHvar * sigx8/(sqgv *
    gisq)^2 * dmusig2/(pmusig1 - pmusig2), FUN = "*"), Zigisq))))) + 0.5 * ((sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar * (-(S * giepsi/gisq))^2/ewv))
    }) + S^2 * (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * 2 * (giepsi/gisq)/gisq/ewv,
    FUN = "*"), (Zigiepsi - sigZ1)) + 2 * ((sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgiepsi[[x]], as.matrix(wHvar * giepsi/gisq/ewv))
  }) - (2 * crossprod(sweep((Zigisq), MARGIN = 1, STATS = wHvar/gisq * giepsi/gisq/ewv,
    FUN = "*"), (Zigiepsi - sigZ1)) + sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * giepsi/gisq * giepsi/gisq/ewv))
  }))) + crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar/gisq/ewv, FUN = "*"),
    (Zigiepsi - sigZ1)))))) - 0.5 * ((sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar/gisq))
  }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar/gisq/gisq, FUN = "*"),
    Zigisq)))
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
#' @param modelType specification of inefficiency model G(t)u_i
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
uninormAlgOpt_gzit <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, modelType,
  ngZGvar, Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstuninorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, modelType = modelType,
      ngZGvar = ngZGvar, initIter = initIter, initAlg = initAlg, gHvar = gHvar,
      whichStart = whichStart, tol = tol, printInfo = printInfo)
    Inituni <- start_st$inituni
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(puninormlike_gzit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p, ngZGvar = ngZGvar, gHvar = gHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel BC92-type Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) {
    -sum(puninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgraduninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = puninormlike_gzit,
    grad = pgraduninormlike_gzit, hess = phessuninormlike_gzit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) {
    -sum(puninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgraduninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
    prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(puninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessuninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(puninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessuninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(puninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgraduninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessuninormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgraduninormlike_gzit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessuninormlike_gzit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessuninormlike_gzit(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- puninormlike_gzit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgraduninormlike_gzit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
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
puninormeff_gzit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar)]
  eta <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar + object$ngZGvar)]
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
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
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
