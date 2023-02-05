################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Two types: - Pitt and Lee (1981) specification - PL81                        #
#            - Time Invariant Inefficiency                                     #
# Convolution: truncated normal - normal                                       #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for truncated normal-normal distribution
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
ptruncnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  muHvar, nmuZUvar, uHvar, vHvar, Yvar, Xvar, pindex, TT, S,
  wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon_i)/(exp(Wv) +
    TT * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
  ll <- -TT/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 * log(exp(Wv) +
    TT * exp(Wu)) + pnorm(mustar/sigmastar, log.p = TRUE) -
    epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2 -
    1/2 * (mu^2/exp(Wu)) - pnorm(mu/sqrt(exp(Wu)), log.p = TRUE)
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for truncated normal-normal distribution
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
psttruncnorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  muHvar, nmuZUvar, nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar,
  whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csttruncnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nmuZUvar = 1, nuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTrunc <- NULL
  } else {
    cat("Initialization: SFA + truncated-normal distribution...\n")
    initTrunc <- maxLik::maxLik(logLik = ctruncnormlike,
      start = csttruncnorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nmuZUvar = 1, nuZUvar = 1, muHvar = muHvar[,
          1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradtruncnormlike, hess = chesstruncnormlike,
      method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nmuZUvar = 1, nuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE], uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE], Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)
    Esti <- initTrunc$estimate
  }
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nmuZUvar >
    1) {
    rep(0, nmuZUvar - 1)
  }, Esti[nXvar + 2], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 3], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  })
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)))
  return(list(StartVal = StartVal, initTrunc = initTrunc))
}

# Gradient of the likelihood function ----------
#' gradient for truncated normal-normal distribution
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
pgradtruncnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  muHvar, nmuZUvar, uHvar, vHvar, Yvar, Xvar, pindex, TT, S,
  wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  ewu <- exp(Wu)
  ewuh <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmasq <- (ewv + TT * ewu)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (mu * ewv - S * ewu * epsilon_i)/(ssqx2)
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dpmu <- dmusig/pmusig
  dmu <- dnorm(mu/exp(Wu/2))
  pmu <- pnorm(mu/exp(Wu/2))
  wvsq <- (1 - ewv/sigmasq)
  wu2sq <- (1 - TT * ewu/sigmasq)
  pmustar <- pmusig * sigmastar
  ttsq <- (TT/sigmasq)
  epsivsq <- (epsilon_isq/(2 * ewv)^2)
  ewup <- ewuh * pmu
  epsi_ss <- S * epsilon_i/(ssqx2)
  muwu <- mu/ewu
  musq <- (musig + dpmu)
  muwusq <- ((mu)^2/ewu)
  sigx0 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx1 <- (mu * ewv - S * ewu * epsilon_i)
  sigx2 <- (sigx1/ewv + dmusig * ewu/(pmustar))
  sigx3 <- (sigx0 * sigx1/ssq + epsi_ss)
  sigx4 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx5 <- (mu/(ssqx2) - sigx4 * sigx1/ssq)
  sigx6 <- (sigx1/ewu + dmusig * ewv/(pmustar))
  gradll <- cbind(sweep(X_iM, MARGIN = 1, STATS = -S * sigx2/sigmasq,
    FUN = "*") - sweep(Xepsi_i, MARGIN = 1, STATS = 2/(2 *
    ewv), FUN = "*"), sweep(muHvar, MARGIN = 1, STATS = (sigx6/sigmasq -
    (dmu/(ewup) + muwu)), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = -(((sigx3 * musq + 0.5 * ttsq) * ewu - (0.5 *
      muwusq + 0.5 * (mu * dmu/(ewup))))), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ((musq * sigx5 + 2 *
      epsivsq - 0.5/sigmasq) * ewv - 0.5 * (TT - 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for truncated normal-normal distribution
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
phesstruncnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  muHvar, nmuZUvar, uHvar, vHvar, Yvar, Xvar, pindex, TT, S,
  wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  X_iM <- apply(-Xvar, 2, function(x) {
    tapply(x, pindex[, 1], sum)
  })
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) {
      tapply(x, pindex[, 1], sum)
    })
  }
  ewu <- exp(Wu)
  ewuh <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmasq <- (ewv + TT * ewu)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  musig <- (mu * ewv - S * ewu * epsilon_i)/(ssqx2)
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dpmu <- dmusig/pmusig
  dmu <- dnorm(mu/exp(Wu/2))
  pmu <- pnorm(mu/exp(Wu/2))
  wvsq <- (1 - ewv/sigmasq)
  wu2sq <- (1 - TT * ewu/sigmasq)
  pmustar <- pmusig * sigmastar
  ttsq <- (TT/sigmasq)
  epsivsq <- (epsilon_isq/(2 * ewv)^2)
  ewup <- ewuh * pmu
  epsi_ss <- S * epsilon_i/(ssqx2)
  muwu <- mu/ewu
  musq <- (musig + dpmu)
  muwusq <- ((mu)^2/ewu)
  sigx0 <- (0.5 * (wu2sq * ewv/sigmastar) + TT * sigmastar)
  sigx1 <- (mu * ewv - S * ewu * epsilon_i)
  sigx2 <- (sigx1/ewv + dmusig * ewu/(pmustar))
  sigx3 <- (sigx0 * sigx1/ssq + epsi_ss)
  sigx4 <- (0.5 * (wvsq * ewu/sigmastar) + sigmastar)
  sigx5 <- (mu/(ssqx2) - sigx4 * sigx1/ssq)
  sigx6 <- (sigx1/ewu + dmusig * ewv/(pmustar))
  sigx7 <- (sigx2 * dpmu - ewu/sigmastar)
  sigx8 <- (sigx1/(pmustar) + dmusig * ewu * ewv/(pmustar)^2)
  sigx9 <- (sigx0 * sigmasq * sigx1 * sigmastar/ssq)
  sigx10 <- (sigx4 * sigmasq * sigx1 * sigmastar/ssq)

  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar,
    ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(X_iM, MARGIN = 1,
    STATS = -S^2 * wHvar * ((sigx1/(ewv * pmustar) + dmusig *
      ewu/(pmustar)^2) * dmusig/sigmasq - 1/ewv) * ewu/sigmasq,
    FUN = "*"), X_iM) - sapply(1:nXvar, function(x) {
    crossprod(Xsq[[x]], as.matrix(wHvar * 2/(2 * ewv)))
  })
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(X_iM,
    MARGIN = 1, STATS = wHvar * S * (sigx8 * dmusig/sigmasq -
      1)/sigmasq, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(X_iM, MARGIN = 1, STATS = -wHvar *
    (S * (sigx7 * sigx3/sigmasq + musq * (1/(ssqx2) - sigx0 *
      ewu/ssq)) * ewu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 2 * (ewv * 2/(2 * ewv)^2),
    FUN = "*"), vHvar) + crossprod(sweep(X_iM, MARGIN = 1,
    STATS = wHvar * S * ewv * (sigx7 * sigx5/sigmasq + musq *
      sigx4 * ewu/ssq), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    ((1/ewu - (sigx1/(ewu * pmustar) + dmusig * ewv/(pmustar)^2) *
      dmusig/sigmasq) * ewv/sigmasq + dmu * (dmu/(ewup)^2 +
      mu/(ewuh^3 * pmu)) - 1/ewu), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx3 * (ewv/sigmastar -
      sigx6 * dpmu)/sigmasq + musq * sigx0 * ewv/ssq) *
      ewu - (0.5 * (((1 - mu^2/ewuh^2)/(ewup) - mu * dmu/(ewup)^2) *
      dmu) + muwu))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (musq * (1/(ssqx2) - sigx4 *
      ewv/ssq) + (ewv/sigmastar - sigx6 * dpmu) * sigx5/sigmasq) *
      ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (((((musq * dpmu - 1) *
      sigx3^2 - 0.5 * (TT^2/sigmasq^2)) * ewu + ((sigx0 *
      (mu * ewv - (2 * sigx9 + 3 * (S * epsilon_i)) * ewu) +
      (0.5 * (TT * ewu/sigmasq) - 0.5 * (0.5 * wu2sq +
        TT * ewu/sigmasq)) * wu2sq * ewv * sigx1/sigmastar)/ssq +
      epsi_ss) * musq + 0.5 * ttsq) * ewu + 0.5 * muwusq -
      0.5 * (mu * (0.5 * (mu^2/(ewuh^3 * pmu)) - (0.5 *
        (ewup) - 0.5 * (mu * dmu))/(ewup)^2) * dmu))),
    FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((musq * dpmu - 1) * sigx3 * sigx5 +
      0.5 * (TT/sigmasq^2) - ((0.5 * (wu2sq * ewv/sigmasq) +
      0.5 * (1 + ewv * (TT * ewu/sigmasq - 1)/sigmasq -
        0.5 * (wvsq * wu2sq))) * sigx1/sigmastar + mu *
      sigx0 - sigx4 * (2 * sigx9 + S * epsilon_i)) * musq/ssq) *
      ewu * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((1 - musq * dpmu) * sigx5^2 +
      0.5/sigmasq^2 - 16 * (ewv * epsilon_isq/(2 * ewv)^4)) *
      ewv + musq * (mu/(ssqx2) - (((3 * (mu) - 2 * sigx10) *
      ewv - S * ewu * epsilon_i) * sigx4 + (0.5 * (ewv/sigmasq) -
      0.5 * (0.5 * wvsq + ewv/sigmasq)) * wvsq * ewu *
      sigx1/sigmastar)/ssq) + 2 * epsivsq - 0.5/sigmasq) *
      ewv, FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncated normal-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar_c matrix of Zmu variables for pooled data
#' @param uHvar_c matrix of Zu variables for pooled data
#' @param vHvar_c matrix of Zv variables for pooled data
#' @param muHvar_p matrix of Zmu variables for cross-section
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
truncnormAlgOpt_pl81 <- function(start, olsParam, dataTable,
  S, nXvar, muHvar_c, muHvar_p, nmuZUvar, uHvar_c, uHvar_p,
  nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c = wHvar_c,
  wHvar_p = wHvar_p, pindex, TT, method, printInfo, itermax,
  stepmax, tol, gradtol, whichStart, initIter, initAlg, hessianType,
  qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttruncnorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar_c, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar_c, tol = tol, whichStart = whichStart,
      initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    Inittrunc <- start_st$inittrunc
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(ptruncnormlike_pl81(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("SFA Panel PL81 Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptruncnormlike_pl81,
    grad = pgradtruncnormlike_pl81, hess = phesstruncnormlike_pl81,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstruncnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phesstruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(ptruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradtruncnormlike_pl81(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phesstruncnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
      pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtruncnormlike_pl81(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
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
      mleObj$hessian <- phesstruncnormlike_pl81(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstruncnormlike_pl81(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptruncnormlike_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtruncnormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, Inittrunc = Inittrunc))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncated normal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ptruncnormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
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
  TT <- as.numeric(table(pindex[, 1]))
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar_p)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar_p)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar_p)))
  epsilon_it <- model.response(model.frame(object$formula,
    data = object$dataTable)) - as.numeric(crossprod(matrix(beta),
    t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon_i)/(exp(Wv) +
    TT * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + TT * exp(Wu)))
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    m <- ifelse(mustar > 0, mustar, 0)
    teMO <- exp(-m)
    teBC <- exp(-mustar + 1/2 * sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar)/pnorm(mustar/sigmastar)
    teBCLB <- exp(-uUB)
    teBCUB <- exp(-uLB)
    teBC_reciprocal <- exp(mustar + 1/2 * sigmastar^2) *
      pnorm(mustar/sigmastar + sigmastar)/pnorm(mustar/sigmastar)
    res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
      uUB = uUB, teJLMS = teJLMS, m = m, teMO = teMO, teBC = teBC,
      teBCLB = teBCLB, teBCUB = teBCUB, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
      uUB = uUB, m = m)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
