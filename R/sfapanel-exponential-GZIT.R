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
# Convolution: exponential - normal                                            #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for exponential-normal distribution
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
pexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  mustar <- -(exp(Wv)/(gisq * exp(Wu/2)) + S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
  ll <- 1/2 * Wu - (TT - 1)/2 * Wv - (TT - 1)/2 * log(2 * pi) - 1/2 * log(gisq) -
    epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2 + pnorm(mustar/sigmastar,
    log.p = TRUE)
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for exponential-normal distribution
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
pstexponorm_gzit <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, ngZGvar, gHvar, S, wHvar, modelType, whichStart, initIter, initAlg,
  printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstexponorm(olsObj = olsObj, epsiRes = epsiRes, S = S, nuZUvar = 1,
      uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initExpo <- NULL
  } else {
    cat("Initialization: SFA + exponential-normal distribution...\n")
    initExpo <- maxLik::maxLik(logLik = cexponormlike, start = cstexponorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), grad = cgradexponormlike,
      method = initAlg, control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar)
    Esti <- initExpo$estimate
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
  return(list(StartVal = StartVal, initExpo = initExpo))
}

# Gradient of the likelihood function ----------
#' gradient for exponential-normal distribution
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
pgradexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar, ngZGvar, gHvar) {
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
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it, FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1], sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1], sum))
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sqvg <- sqrt(ewv/gisq)
  sigmastar <- (sqvg * gisq)
  mustar <- (ewv/ewu_h + S * giepsi)
  musig <- (mustar/sigmastar)
  dmusig <- dnorm(-musig, 0, 1)
  pmusig <- pnorm(-musig)
  wug <- (ewu_h * gisq)
  musq <- (mustar/sigmastar^2)
  muwusq <- (1/wug - 0.5 * musq)
  vgsig <- (sqvg - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigisq, MARGIN = 1, STATS = vgsig/sigmastar^2, FUN = "*")
  sigx1 <- (pmusig * sqvg)
  sigx2 <- (mustar/ewv - dmusig/sigx1)
  sigx3 <- (0.5 * (dmusig * ewv/sigx1) - 0.5 * mustar)
  sigx4 <- ((1/ewu_h - dmusig/sigx1) * ewv + S * giepsi)
  sigx5 <- (mustar/sigmastar - dmusig/pmusig)
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(sigZ1,
    MARGIN = 1, STATS = mustar, FUN = "*")
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * sigx2/gisq, FUN = "*") - sweep(Xepsi_i,
    MARGIN = 1, STATS = 1/(2 * ewv), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx3/wug +
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx4 * muwusq + 2 *
    (ewv * epsilon_isq/(2 * ewv)^2) - 0.5 * (TT - 1)), FUN = "*"), sweep(sigZ2,
    MARGIN = 1, STATS = sigx5, FUN = "*") - sweep(Zigisq, MARGIN = 1, STATS = 0.5/gisq,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' optimizations solve for exponential-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
phessexponormlike_gzit <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, pindex, TT, S, wHvar, ngZGvar, gHvar) {
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
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1], sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it, FUN = "*")
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
  ewv <- exp(Wv)
  ewu_h <- exp(Wu/2)
  sqvg <- sqrt(ewv/gisq)
  sigmastar <- (sqvg * gisq)
  mustar <- (ewv/ewu_h + S * giepsi)
  musig <- (mustar/sigmastar)
  dmusig <- dnorm(-musig, 0, 1)
  pmusig <- pnorm(-musig)
  wug <- (ewu_h * gisq)
  musq <- (mustar/sigmastar^2)
  muwusq <- (1/wug - 0.5 * musq)
  vgsig <- (sqvg - 0.5 * (ewv/sigmastar))
  sigZ1 <- sweep(Zigisq, MARGIN = 1, STATS = vgsig/sigmastar^2, FUN = "*")
  sigx1 <- (pmusig * sqvg)
  sigx2 <- (mustar/ewv - dmusig/sigx1)
  sigx3 <- (0.5 * (dmusig * ewv/sigx1) - 0.5 * mustar)
  sigx4 <- ((1/ewu_h - dmusig/sigx1) * ewv + S * giepsi)
  sigx5 <- (mustar/sigmastar - dmusig/pmusig)
  sigZ2 <- sweep(Zigiepsi, MARGIN = 1, STATS = S/sigmastar, FUN = "*") - sweep(sigZ1,
    MARGIN = 1, STATS = mustar, FUN = "*")
  sigx7 <- dmusig * ewv/sigx1^2
  sigx8 <- (dmusig/sigx1^2 - mustar/(ewv * pmusig * sqvg))
  sigx9 <- (0.5 * (mustar/sigx1) - 0.5 * (sigx7)) * dmusig/gisq
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + ngZGvar, ncol = nXvar + nuZUvar +
    nvZVvar + ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    S^2 * (1/ewv - dmusig * sigx8/gisq)/gisq, FUN = "*"), Xgi) - sapply(1:nXvar,
    function(x) crossprod(Xsq[[x]], as.matrix(wHvar/ewv)))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S * (0.5 * (dmusig * (sigx7 - mustar/sigx1)/gisq) - 0.5)/wug,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 2 * (ewv/(2 * ewv)^2), FUN = "*") + sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * ((1 - dmusig * sigx8 * ewv/gisq) * muwusq -
      0.5 * (sigx4/sigmastar^2)), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
    ngZGvar)] <- sapply(1:ngZGvar, function(x) {
    crossprod(Xzigi[[x]], as.matrix(wHvar * S * sigx5/sigmastar))
  }) - crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S * sigx5, FUN = "*"),
    sigZ1) + crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * S * (1/sqvg - dmusig *
    (dmusig/sigx1 - mustar/ewv)/pmusig)/gisq, FUN = "*"), sigZ2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 + 0.5 * (sigx9)) * ewv/(ewu_h^2 * gisq) -
      0.5 * (sigx3 * ewu_h * gisq/wug^2)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (sigx4 *
    (0.25 * (ewv/(ewu_h * sigmastar^2)) - 0.5 * (ewu_h * gisq/wug^2)) - (sigx9 +
    0.5) * muwusq * ewv/ewu_h), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + ngZGvar)] <- crossprod(uHvar, sweep(sigZ1, MARGIN = 1,
    STATS = wHvar * 0.5 * (sigx5 * ewv)/ewu_h, FUN = "*") - sweep(sigZ2, MARGIN = 1,
    STATS = wHvar * ((0.5 * mustar - 0.5 * (dmusig * ewv/sigx1)) * dmusig/pmusig +
      0.5 * (ewv/sqvg))/gisq/ewu_h, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((muwusq * mustar - 1)/sigx1 + (0.5 * (pmusig/sigmastar) - muwusq * dmusig) *
      ewv/sigx1^2) * dmusig + 1/ewu_h) * muwusq + 2 * ((1 - 8 * (ewv^2/(2 *
      ewv)^2)) * epsilon_isq/(2 * ewv)^2) - 0.5 * (sigx4 * (1/ewu_h - mustar *
      gisq/sigmastar^2)/sigmastar^2)) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- crossprod(vHvar, sweep(sigZ2,
    MARGIN = 1, STATS = wHvar * muwusq * (ewv/sqvg - ((dmusig/sigx1 - 1/ewu_h) *
      ewv - S * giepsi) * dmusig/pmusig), FUN = "*") - (sweep(Zigisq, MARGIN = 1,
    STATS = wHvar * (((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/sigmastar^2)))/sqvg -
      vgsig * gisq/sigmastar^2) * mustar + vgsig/ewu_h) * sigx5 * ewv/sigmastar^2,
    FUN = "*") + sweep(Zigiepsi, MARGIN = 1, STATS = wHvar * 0.5 * (S/sqvg) *
    sigx5 * ewv/sigmastar^2, FUN = "*")))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + ngZGvar)] <- S *
    (sapply(1:ngZGvar, function(x) {
      crossprod(Zisqgiepsi[[x]], as.matrix(wHvar * sigx5/sigmastar))
    }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx5 * vgsig/sigmastar^2,
      FUN = "*"), Zigiepsi)) - (crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar *
    sigx5 * (0.5 * (vgsig/sigmastar^2) - 0.5/(sqvg * gisq^2)) * ewv * mustar/sigmastar^2,
    FUN = "*"), Zigisq) + crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar *
    sigx5 * S * vgsig/sigmastar^2, FUN = "*"), Zigisq) + sapply(1:ngZGvar, function(x) {
    crossprod(Zisqgisq[[x]], as.matrix(wHvar * sigx5 * mustar * vgsig/sigmastar^2))
  }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar * sigx5 * mustar * vgsig *
    2 * (sqvg * vgsig * gisq/sigmastar^2)/sigmastar^2, FUN = "*"), Zigisq)) +
    crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar * S/sigmastar, FUN = "*"),
      sigZ2) - (crossprod(sweep(sigZ1, MARGIN = 1, STATS = wHvar * mustar,
    FUN = "*"), sigZ2) + crossprod(sweep(sigZ2, MARGIN = 1, STATS = wHvar * dmusig *
    (dmusig/pmusig - mustar/sigmastar)/pmusig, FUN = "*"), sigZ2)) - 0.5 * (sapply(1:ngZGvar,
    function(x) {
      crossprod(Zisqgisq[[x]], as.matrix(wHvar/gisq))
    }) - crossprod(sweep(Zigisq, MARGIN = 1, STATS = wHvar/gisq/gisq, FUN = "*"),
    Zigisq))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for exponential-normal distribution
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
exponormAlgOpt_gzit <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, gHvar, ngZGvar,
  Yvar, Xvar, pindex, TT, wHvar_c, wHvar_p, whichStart, initIter, initAlg, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, modelType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstexponorm_gzit(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, wHvar = wHvar_c, ngZGvar = ngZGvar,
      gHvar = gHvar, modelType = modelType, tol = tol, printInfo = printInfo)
    InitExpo <- start_st$initExpo
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(pexponormlike_gzit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(pexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = pexponormlike_gzit,
    grad = pgradexponormlike_gzit, hess = phessexponormlike_gzit, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) {
    -sum(pexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
    prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(pexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p),
        "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(pexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(pexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgradexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessexponormlike_gzit(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradexponormlike_gzit(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessexponormlike_gzit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessexponormlike_gzit(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- pexponormlike_gzit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradexponormlike_gzit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam, InitExpo = InitExpo))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
      mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for exponential-normal distribution
#' @param object object of class sfapanel1
#' @param level level for confidence interval
#' @noRd
pexponormeff_gzit <- function(object, level) {
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
  mustar <- -(exp(Wv)/(gisq * exp(Wu/2)) + object$S * giepsi/gisq)
  sigmastar <- sqrt(exp(Wv)/gisq)
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
    res$teBC <- exp(-res$mustar * git + 1/2 * res$sigmastar^2 * git^2) * pnorm(res$mustar/res$sigmastar -
      res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    teBC_reciprocal <- exp(res$mustar * git + 1/2 * res$sigmastar^2 * git^2) *
      pnorm(res$mustar/res$sigmastar + res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
