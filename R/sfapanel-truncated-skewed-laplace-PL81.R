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
# Convolution: truncated skewed laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for halfnormal-normal distribution
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
ptslnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  mustar1 <- -(exp(Wv)/(TT * exp(Wu/2)) + S * epsilon_i/TT)
  mustar2 <- -((1 + lambda) * exp(Wv)/(TT * exp(Wu/2)) + S *
    epsilon_i/TT)
  sigmastar <- sqrt(exp(Wv)/TT)
  ll <- log(1 + lambda) - 1/2 * Wu - (TT - 1)/2 * Wv - 1/2 *
    log(TT) - (TT - 1)/2 * log(2 * pi) - log(2 * lambda +
    1) + log(2 * exp(-1/2 * (epsilon_isq/exp(Wv) - (mustar1/sigmastar)^2)) *
    pnorm(mustar1/sigmastar) - exp(-1/2 * (epsilon_isq/exp(Wv) -
    (mustar2/sigmastar)^2)) * pnorm(mustar2/sigmastar))
  return(wHvar * ll)
}

# starting value for the log-likelihood ----------
#' starting values for halfnormal-normal distribution
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
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
psttslnorm_pl81 <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, itermax, printInfo,
  tol) {
  cat("Initialization: SFA + truncated-skewed-laplace-normal distribution...\n")
  initTsl <- maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradtslnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = if (printInfo) 2 else 0, reltol = tol),
    nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initTsl$estimate
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar >
    1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 2], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, Esti[nXvar + 3])
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "lambda")
  names(initTsl$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]),
    "lambda")
  return(list(StartVal = StartVal, initTsl = initTsl))
}

# Gradient of the likelihood function ----------
#' gradient for halfnormal-normal distribution
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
pgradtslnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  X_iM <- apply(-Xvar, 2, function(x) tapply(x, pindex[, 1],
    sum))
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  wvwu_epsi <- (ewv/ewu_h + S * epsilon_i)
  epsi_vu <- ((1 + lambda) * ewv/ewu_h + S * epsilon_i)
  wvtt <- ewv/TT
  ssq_v <- (TT * sqrt(wvtt))
  wutt <- TT * ewu_h
  musig <- (epsi_vu/ssq_v)
  dmusig <- dnorm(-musig)
  pmusig <- pnorm(-musig)
  epsi_vu2 <- (wvwu_epsi/ssq_v)
  epsi_vu3 <- (wvwu_epsi/ssq_v^2)
  pepsi <- pnorm(-epsi_vu2)
  depsi <- dnorm(-epsi_vu2, 0, 1)
  sigx1_1 <- (epsi_vu/ssq_v^2)
  sigx1_2 <- (depsi * ewv/sqrt(wvtt))
  sigx2_1 <- ((1 + lambda)/(wutt) - 0.5 * sigx1_1)
  sigx2_2 <- (1/(wutt) - 0.5 * epsi_vu3)
  sigx3_1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq/ewv - (-epsi_vu2)^2)))
  sigx4 <- (2 * (sigx3_2 * pepsi) - sigx3_1 * pmusig)
  sigx5 <- ((0.5 * sigx1_2 - 0.5 * (wvwu_epsi * pepsi)) * sigx3_2)
  sigx6 <- (0.5 * (dmusig * ewv/sqrt(wvtt)) - 0.5 * (epsi_vu *
    pmusig))
  sigx7 <- (2 * sigx5 - sigx6 * (1 + lambda) * sigx3_1)
  ttsigu <- (TT * sigx4 * ewu_h)
  sigx8 <- ((2 * (sigx2_2 * wvwu_epsi) + epsilon_isq/ewv) *
    pepsi)
  sigx9 <- (0.5 * sigx8 - sigx2_2 * depsi * ewv/sqrt(wvtt))
  epsi_vusq <- (2 * (epsi_vu * sigx2_1) + epsilon_isq/ewv)
  sigx10 <- (0.5 * (epsi_vusq * pmusig) - sigx2_1 * dmusig *
    ewv/sqrt(wvtt))
  sigx11 <- (2 * (sigx9 * sigx3_2) - sigx10 * sigx3_1)
  sigx12 <- (epsi_vu * pmusig - dmusig * ewv/sqrt(wvtt))
  Xsig1 <- (sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig/ewv,
    FUN = "*") - sweep(X_iM, MARGIN = 1, STATS = 2 * (S *
    epsi_vu/TT) * pmusig/ewv, FUN = "*"))
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig/ssq_v,
    FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pepsi/ewv,
    FUN = "*") - sweep(X_iM, MARGIN = 1, STATS = 2 * (S *
    wvwu_epsi/TT) * pepsi/ewv, FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = S * depsi/ssq_v,
    FUN = "*")
  Xsig5 <- sweep(0.5 * Xsig1 + Xsig2, MARGIN = 1, STATS = sigx3_1,
    FUN = "*")
  Xsig6 <- sweep((0.5 * Xsig3 + Xsig4), MARGIN = 1, STATS = sigx3_2,
    FUN = "*")
  gradll <- cbind(sweep((Xsig5 - 2 * Xsig6), MARGIN = 1, STATS = 1/sigx4,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx7/ttsigu -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx11/sigx4 -
    0.5 * (TT - 1)), FUN = "*"), 1/(1 + lambda) - (sigx12 *
    sigx3_1/ttsigu + 2/(1 + 2 * lambda)))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for halfnormal-normal distribution
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
phesstslnormlike_pl81 <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, pindex, TT, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[,
    1], sum))
  X_iM <- apply(-Xvar, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[,
      1], sum))
  }
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  wvwu_epsi <- (ewv/ewu_h + S * epsilon_i)
  epsi_vu <- ((1 + lambda) * ewv/ewu_h + S * epsilon_i)
  wvtt <- ewv/TT
  ssq_v <- (TT * sqrt(wvtt))
  wutt <- TT * ewu_h
  musig <- (epsi_vu/ssq_v)
  dmusig <- dnorm(-musig)
  pmusig <- pnorm(-musig)
  epsi_vu2 <- (wvwu_epsi/ssq_v)
  epsi_vu3 <- (wvwu_epsi/ssq_v^2)
  pepsi <- pnorm(-epsi_vu2)
  depsi <- dnorm(-epsi_vu2, 0, 1)
  sigx1_1 <- (epsi_vu/ssq_v^2)
  sigx1_2 <- (depsi * ewv/sqrt(wvtt))
  sigx2_1 <- ((1 + lambda)/(wutt) - 0.5 * sigx1_1)
  sigx2_2 <- (1/(wutt) - 0.5 * epsi_vu3)
  sigx3_1 <- exp(-(0.5 * (epsilon_isq/ewv - (-musig)^2)))
  sigx3_2 <- exp(-(0.5 * (epsilon_isq/ewv - (-epsi_vu2)^2)))
  sigx4 <- (2 * (sigx3_2 * pepsi) - sigx3_1 * pmusig)
  sigx5 <- ((0.5 * sigx1_2 - 0.5 * (wvwu_epsi * pepsi)) * sigx3_2)
  sigx6 <- (0.5 * (dmusig * ewv/sqrt(wvtt)) - 0.5 * (epsi_vu *
    pmusig))
  sigx7 <- (2 * sigx5 - sigx6 * (1 + lambda) * sigx3_1)
  ttsigu <- (TT * sigx4 * ewu_h)
  sigx8 <- ((2 * (sigx2_2 * wvwu_epsi) + epsilon_isq/ewv) *
    pepsi)
  sigx9 <- (0.5 * sigx8 - sigx2_2 * depsi * ewv/sqrt(wvtt))
  epsi_vusq <- (2 * (epsi_vu * sigx2_1) + epsilon_isq/ewv)
  sigx10 <- (0.5 * (epsi_vusq * pmusig) - sigx2_1 * dmusig *
    ewv/sqrt(wvtt))
  sigx11 <- (2 * (sigx9 * sigx3_2) - sigx10 * sigx3_1)
  sigx12 <- (epsi_vu * pmusig - dmusig * ewv/sqrt(wvtt))
  Xsig1 <- (sweep(2 * Xepsi_i, MARGIN = 1, STATS = pmusig/ewv,
    FUN = "*") - sweep(X_iM, MARGIN = 1, STATS = 2 * (S *
    epsi_vu/TT) * pmusig/ewv, FUN = "*"))
  Xsig2 <- sweep(X_iM, MARGIN = 1, STATS = S * dmusig/ssq_v,
    FUN = "*")
  Xsig3 <- sweep(2 * Xepsi_i, MARGIN = 1, STATS = pepsi/ewv,
    FUN = "*") - sweep(X_iM, MARGIN = 1, STATS = 2 * (S *
    wvwu_epsi/TT) * pepsi/ewv, FUN = "*")
  Xsig4 <- sweep(X_iM, MARGIN = 1, STATS = S * depsi/ssq_v,
    FUN = "*")
  Xsig5 <- sweep(0.5 * Xsig1 + Xsig2, MARGIN = 1, STATS = sigx3_1,
    FUN = "*")
  Xsig6 <- sweep((0.5 * Xsig3 + Xsig4), MARGIN = 1, STATS = sigx3_2,
    FUN = "*")
  Xsig7 <- sweep((Xsig5 - 2 * Xsig6), MARGIN = 1, STATS = sigx7/sigx4,
    FUN = "*")
  ddsq <- (depsi * wvwu_epsi/ssq_v)
  sigx13 <- (0.5 * sigx1_2 - 0.5 * (wvwu_epsi * pepsi))
  dmusq <- (epsi_vu * dmusig/ssq_v)
  sigx14 <- (0.5 * dmusq + 0.5 * (pmusig - epsi_vu * dmusig/ssq_v))
  ttusq <- (TT * sigx4 * ewu_h^2)
  sigx15 <- (2 * (sigx2_2 * wvwu_epsi) + epsilon_isq/ewv)
  sigx16 <- (0.25 * (ewv/(ewu_h * ssq_v^2)) - 0.5 * (wutt/(wutt)^2))
  Xsig8 <- sweep(X_iM, MARGIN = 1, STATS = (S * wvwu_epsi/TT),
    FUN = "*")
  Xsig9 <- sweep(X_iM, MARGIN = 1, STATS = (S * epsi_vu/TT),
    FUN = "*")
  ttvq <- (TT^2 * sqrt(wvtt))
  sigx17 <- (epsi_vu * sigx2_1/(TT * ewv) + 0.5/ssq_v^2)
  Xsig10 <- sweep((2 * Xepsi_i - 2 * Xsig9), MARGIN = 1, STATS = 0.5 *
    (sigx6/ewv) * (1 + lambda) * sigx3_1, FUN = "*") + sweep(X_iM,
    MARGIN = 1, STATS = S * sigx14 * (1 + lambda) * sigx3_1,
    FUN = "*")
  Xsig11 <- sweep((2 * Xepsi_i - 2 * Xsig8), MARGIN = 1, STATS = 0.5 *
    (sigx13/ewv) * sigx3_2, FUN = "*") + sweep(X_iM, MARGIN = 1,
    STATS = S * (0.5 * ddsq + 0.5 * (pepsi - depsi * wvwu_epsi/ssq_v)) *
      sigx3_2, FUN = "*")
  Xsig12 <- sweep((Xsig10 - 2 * Xsig11), MARGIN = 1, STATS = 1/ttsigu,
    FUN = "*") - sweep((Xsig5 - 2 * Xsig6), MARGIN = 1, STATS = TT *
    sigx7 * ewu_h/ttsigu^2, FUN = "*")
  Xsig13 <- 0.5 * (sweep(X_iM, MARGIN = 1, STATS = 2 * (S *
    (1/(wutt) - wvwu_epsi/ssq_v^2)) * pepsi, FUN = "*") +
    sweep(Xepsi_i, MARGIN = 1, STATS = 2/ewv * pepsi, FUN = "*") -
    sweep(X_iM, MARGIN = 1, STATS = S * sigx15 * depsi/ssq_v,
      FUN = "*"))
  Xsig14 <- sweep(X_iM, MARGIN = 1, STATS = S * (sigx2_2 *
    wvwu_epsi/(TT * ewv) + 0.5/ssq_v^2) * depsi * ewv/sqrt(wvtt),
    FUN = "*")
  Xsig15 <- sweep((2 * Xepsi_i - 2 * Xsig8), MARGIN = 1, STATS = 0.5 *
    (sigx9/ewv), FUN = "*")
  Xsig16 <- sweep((Xsig5 - 2 * Xsig6), MARGIN = 1, STATS = sigx11/sigx4,
    FUN = "*")
  Xsig17 <- sweep(X_iM, MARGIN = 1, STATS = 2 * (S * ((1 +
    lambda)/(wutt) - epsi_vu/ssq_v^2)) * pmusig, FUN = "*") +
    sweep(Xepsi_i, MARGIN = 1, STATS = 2/ewv * pmusig, FUN = "*")
  Xsig18 <- sweep(X_iM, MARGIN = 1, STATS = S * epsi_vusq *
    dmusig/ssq_v, FUN = "*")
  Xsig19 <- sweep(X_iM, MARGIN = 1, STATS = S * sigx17 * dmusig *
    ewv/sqrt(wvtt), FUN = "*") - sweep((2 * Xepsi_i - 2 *
    Xsig9), MARGIN = 1, STATS = 0.5 * (sigx10/ewv), FUN = "*")
  Xsig20 <- sweep(2 * ((Xsig13 + Xsig14 - Xsig15)), MARGIN = 1,
    STATS = sigx3_2/sigx4, FUN = "*") - (sweep(Xsig16, MARGIN = 1,
    STATS = 1/sigx4, FUN = "*") + sweep((0.5 * (Xsig17 -
    Xsig18) + Xsig19), MARGIN = 1, STATS = sigx3_1/sigx4,
    FUN = "*"))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]],
    as.matrix(wHvar * 2 * pmusig * sigx3_1/ewv/sigx4))) -
    2 * crossprod(sweep(X_iM, MARGIN = 1, STATS = (wHvar *
      S^2 * pmusig/TT) * sigx3_1/ewv/sigx4, FUN = "*"),
      X_iM) - crossprod(sweep((2 * Xepsi_i - 2 * Xsig9),
    MARGIN = 1, STATS = wHvar * S * dmusig/ssq_v * sigx3_1/ewv/sigx4,
    FUN = "*"), X_iM)) - (0.5 * crossprod(sweep((0.5 * Xsig1 +
    Xsig2), MARGIN = 1, STATS = wHvar * sigx3_1/ewv/sigx4,
    FUN = "*"), (2 * Xepsi_i - 2 * Xsig9)) + crossprod(sweep(X_iM,
    MARGIN = 1, STATS = wHvar * S^2 * epsi_vu * dmusig/ttvq *
      sigx3_1/ewv/sigx4, FUN = "*"), X_iM)) - 2 * (0.5 *
    (sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar *
      2 * pepsi * sigx3_2/ewv/sigx4))) - 2 * crossprod(sweep(X_iM,
      MARGIN = 1, STATS = wHvar * (S^2 * pepsi/TT) * sigx3_2/ewv/sigx4,
      FUN = "*"), X_iM) - crossprod(sweep((2 * Xepsi_i -
      2 * Xsig8), MARGIN = 1, STATS = wHvar * S * depsi/ssq_v *
      sigx3_2/ewv/sigx4, FUN = "*"), X_iM)) - (0.5 * crossprod(sweep((0.5 *
    Xsig3 + Xsig4), MARGIN = 1, STATS = wHvar * sigx3_2/ewv/sigx4,
    FUN = "*"), (2 * Xepsi_i - 2 * Xsig8)) + crossprod(sweep(X_iM,
    MARGIN = 1, STATS = wHvar * S^2 * depsi * wvwu_epsi/ttvq *
      sigx3_2/ewv/sigx4, FUN = "*"), X_iM))) - crossprod(sweep((Xsig5 -
    2 * Xsig6), MARGIN = 1, STATS = wHvar/sigx4^2, FUN = "*"),
    (Xsig5 - 2 * Xsig6))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(Xsig12,
    sweep(uHvar, MARGIN = 1, STATS = wHvar, FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(Xsig20, sweep(vHvar, MARGIN = 1,
    STATS = wHvar, FUN = "*"))
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(-(sweep(X_iM,
    MARGIN = 1, STATS = wHvar * S * pmusig/ttsigu * sigx3_1,
    FUN = "*") - sweep((2 * Xepsi_i - 2 * Xsig9), MARGIN = 1,
    STATS = wHvar * 0.5 * (sigx12/ewv)/ttsigu * sigx3_1,
    FUN = "*") - sweep((Xsig5 - 2 * Xsig6), MARGIN = 1, STATS = wHvar *
    TT * sigx12 * ewu_h/ttsigu^2 * sigx3_1, FUN = "*")))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((0.25 * ddsq - 0.5 * (0.5 * ddsq - 0.5 * pepsi)) *
      ewv - 0.5 * (sigx13 * wvwu_epsi/TT)) * sigx3_2) -
      ((0.25 * dmusq - 0.5 * (0.5 * dmusq - 0.5 * pmusig)) *
        ewv - 0.5 * (epsi_vu * sigx6/TT)) * (1 + lambda)^2 *
        sigx3_1)/ttusq - (0.5 * ttsigu + 2 * sigx5 -
      sigx6 * (1 + lambda) * sigx3_1) * sigx7/ttsigu^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((0.5 * (0.5 * (sigx15 *
      depsi * ewv/(wutt * sqrt(wvtt))) + 2 * ((sigx16 *
      wvwu_epsi - 0.5 * (sigx2_2 * ewv/ewu_h)) * pepsi)) -
      (((0.25 * (ewv/ssq_v^2) + 0.5 * (sigx2_2 * wvwu_epsi/TT))/ewu_h -
        0.5 * (wutt/(wutt)^2)) * depsi * ewv/sqrt(wvtt) +
        0.5 * (sigx9 * wvwu_epsi/(wutt)))) * sigx3_2) -
      ((0.5 * (0.5 * (epsi_vusq * dmusig * ewv/(wutt *
        sqrt(wvtt))) + 2 * ((epsi_vu * sigx16 - 0.5 *
        (sigx2_1 * ewv/ewu_h)) * pmusig)) - (((0.25 *
        (ewv/ssq_v^2) + 0.5 * (epsi_vu * sigx2_1/TT))/ewu_h -
        0.5 * (wutt/(wutt)^2)) * dmusig * ewv/sqrt(wvtt) +
        0.5 * (epsi_vu * sigx10/(wutt)))) * (1 + lambda) *
        sigx3_1 + sigx11 * sigx7/ttsigu))/sigx4, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- (wHvar * (sigx12 * (0.5 * ttsigu + 2 * sigx5 -
    sigx6 * (1 + lambda) * sigx3_1)/ttsigu^2 + (0.5 * (sigx12 *
    epsi_vu/TT) + 0.5 * (ewv * pmusig)) * (1 + lambda)/ttusq) *
    sigx3_1) %*% uHvar
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((0.5 * (sigx9 * sigx15) +
      0.5 * ((2 * ((sigx2_2/ewu_h - 0.5 * ((1/ewu_h - TT *
        wvwu_epsi/ssq_v^2) * wvwu_epsi/ssq_v^2)) * ewv) -
        epsilon_isq/ewv) * pepsi - sigx2_2 * sigx15 *
        depsi * ewv/sqrt(wvtt)) - (1/(wutt) - ((sigx2_2^2 +
      0.5/ssq_v^2) * wvwu_epsi + 0.5 * ((1/ewu_h - TT *
      wvwu_epsi/ssq_v^2) * ewv/ssq_v^2) + 0.5 * sigx2_2)) *
      depsi * ewv/sqrt(wvtt)) * sigx3_2) - ((0.5 * (sigx10 *
      epsi_vusq) + 0.5 * ((2 * ((sigx2_1 * (1 + lambda)/ewu_h -
      0.5 * (epsi_vu * ((1 + lambda)/ewu_h - TT * epsi_vu/ssq_v^2)/ssq_v^2)) *
      ewv) - epsilon_isq/ewv) * pmusig - sigx2_1 * epsi_vusq *
      dmusig * ewv/sqrt(wvtt)) - ((1 + lambda)/(wutt) -
      ((sigx2_1^2 + 0.5/ssq_v^2) * epsi_vu + 0.5 * (((1 +
        lambda)/ewu_h - TT * epsi_vu/ssq_v^2) * ewv/ssq_v^2) +
        0.5 * sigx2_1)) * dmusig * ewv/sqrt(wvtt)) *
      sigx3_1 + sigx11^2/sigx4))/sigx4, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- (wHvar * (-(((((1 +
    lambda) * pmusig/ewu_h - 0.5 * (dmusig/sqrt(wvtt))) *
    ewv + 0.5 * (sigx12 * epsi_vusq))/ttsigu - TT * sigx12 *
    sigx11 * ewu_h/ttsigu^2) * sigx3_1))) %*% vHvar
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (4/(1 + 2 * lambda)^2 - (((sigx12 *
    epsi_vu/TT + ewv * pmusig)/ttusq + sigx12^2 * sigx3_1/ttsigu^2) *
    sigx3_1 + 1/(1 + lambda)^2)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for tsl-normal distribution
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
tslnormAlgOpt_pl81 <- function(start, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar,
  Yvar, Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttslnorm_pl81(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_c, vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar_c, itermax = itermax, tol = tol,
      printInfo = printInfo)
    Inittsl <- start_st$inittsl
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(ptslnormlike_pl81(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) {
      -sum(ptslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptslnormlike_pl81,
    grad = pgradtslnormlike_pl81, hess = phesstslnormlike_pl81,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtslnormlike_pl81(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) {
    -sum(ptslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hess = function(parm) {
    -phesstslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
  }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(ptslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradtslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phesstslnormlike_pl81(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtslnormlike_pl81(mleObj$par,
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
        names(mleObj$solution) <- names(startVal)
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
      mleObj$hessian <- phesstslnormlike_pl81(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstslnormlike_pl81(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptslnormlike_pl81(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
    TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtslnormlike_pl81(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, Inittsl = Inittsl))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for tsl-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ptslnormeff_pl81 <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    object$nvZVvar + 1]
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
  epsilon_i <- as.numeric(tapply(epsilon_it, pindex[, 1], sum))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  mustar1 <- -(exp(Wv)/(TT * exp(Wu/2)) + object$S * epsilon_i/TT)
  mustar2 <- -((1 + lambda) * exp(Wv)/(TT * exp(Wu/2)) + object$S *
    epsilon_i/TT)
  sigmastar <- sqrt(exp(Wv)/TT)
  a <- mustar1/sigmastar
  b <- mustar2/sigmastar
  A <- -1/2 * (epsilon_isq/exp(Wv) - a^2)
  B <- -1/2 * (epsilon_isq/exp(Wv) - b^2)
  u <- (exp(A) * (dnorm(a) * sigmastar + mustar1 * pnorm(a)) -
    exp(B) * (dnorm(b) * sigmastar + mustar2 * pnorm(b)))/(exp(A) *
    pnorm(a) - exp(B) * pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (exp(A) * exp(1/2 * sigmastar^2 - mustar1) *
      pnorm(a - sigmastar) - exp(B) * exp(1/2 * sigmastar^2 -
      mustar2) * pnorm(b - sigmastar))/(exp(A) * pnorm(a) -
      exp(B) * pnorm(b))
    teBC_reciprocal <- (exp(A) * exp(1/2 * sigmastar^2 +
      mustar1) * pnorm(a + sigmastar) - exp(B) * exp(1/2 *
      sigmastar^2 + mustar2) * pnorm(b + sigmastar))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    res <- data.frame(levels(pindex[, 1]), u = u, teJLMS = teJLMS,
      teBC = teBC, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(levels(pindex[, 1]), u = u)
  }
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}
