################################################################################
#   
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Panel data                                                             #
# Model: Panel Stochastic Frontier Model                                       #
# Two types: - Battese and Coelli 1992 specification  BC92III                  #
#            - u_it = g(zit)u_i                                                #
#            - g(zit) = exp(eta * gHvar)                                       #
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
ptruncnormlike_bc92III <- function(parm, nXvar, nuZUvar, nvZVvar,
  nmuZUvar, uHvar, vHvar, muHvar, Yvar, Xvar, pindex, TT, S,
  ngZGvar, gHvar, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + ngZGvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- (mu * exp(Wv) - exp(Wu) * S * giepsi)/(exp(Wv) +
    gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  ll <- -TT/2 * log(2 * pi) - (TT - 1)/2 * Wv - 1/2 * log(exp(Wv) +
    gisq * exp(Wu)) + pnorm(mustar/sigmastar, log.p = TRUE) -
    epsilon_isq/(2 * exp(Wv)) + 1/2 * (mustar/sigmastar)^2 -
    1/2 * (mu^2/exp(Wu)) - pnorm(mu/sqrt(exp(Wu)), log.p = TRUE)
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
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param nXvar number of main variables (inputs + env. var)
#' @param wHvar vector of weights (weighted likelihood)
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
psttruncnorm_bc92III <- function(olsObj, epsiRes, nXvar, nmuZUvar,
  nuZUvar, nvZVvar, muHvar, uHvar, vHvar, ngZGvar, gHvar, Yvar,
  Xvar, S, wHvar, itermax, printInfo, tol) {
  cat("Initialization: SFA + truncated-normal distribution...\n")
  initTrunc <- maxLik(logLik = ctruncnormlike, start = csttruncnorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nmuZUvar = 1, nuZUvar = 1,
    muHvar = as.matrix(muHvar[, 1]), uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgradtruncnormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = if (printInfo) 2 else 0, reltol = tol),
    nXvar = nXvar, nmuZUvar = 1, nuZUvar = 1, muHvar = as.matrix(muHvar[,
      1]), uHvar = as.matrix(uHvar[, 1]), nvZVvar = 1,
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initTrunc$estimate
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nmuZUvar >
    1) {
    rep(0, nmuZUvar - 1)
  }, Esti[nXvar + 2], if (nuZUvar > 1) {
    rep(0, nuZUvar - 1)
  }, Esti[nXvar + 3], if (nvZVvar > 1) {
    rep(0, nvZVvar - 1)
  }, rep(0, ngZGvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste("Zg_", colnames(gHvar)))
  names(initTrunc$estimate) <- c(names(Esti)[1:nXvar], paste0("Zmu_",
    colnames(muHvar)[1]), paste0("Zu_", colnames(uHvar)[1]),
    paste0("Zv_", colnames(vHvar)[1]))
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
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgradtruncnormlike_bc92III <- function(parm, nXvar, nuZUvar,
  nvZVvar, nmuZUvar, uHvar, vHvar, muHvar, Yvar, Xvar, pindex,
  TT, S, wHvar, ngZGvar, gHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + ngZGvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1],
    sum))
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  muepsi <- (mu * ewv - S * ewu * giepsi)
  musig <- muepsi/(ssqx2)
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dpmu <- dmusig/pmusig
  dmuu <- dnorm(mu/ewu_h, 0, 1)
  pmuu <- pnorm(mu/ewu_h)
  pmustar <- pmusig * sigmastar
  mup <- (musig + dpmu)
  wup <- (ewu_h * pmuu)
  gsig <- (gisq/sigmasq)
  muwu <- ((mu)^2/ewu)
  wuwv <- (ewu * ewv/(ssqx2))
  zsig <- sweep(Zigisq, MARGIN = 1, STATS = 1/sigmasq, FUN = "*")
  rwv <- (1 - ewv/sigmasq)
  ssuv <- (sigmastar - 0.5 * wuwv)
  szsq <- sweep(Zigisq, MARGIN = 1, STATS = ssuv/ssq, FUN = "*")
  sigx1 <- (muepsi/ewv + dmusig * ewu/(pmustar))
  sigx2 <- (muepsi/ewu + dmusig * ewv/(pmustar))
  sigx3 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx4 <- ((0.5 * sigx3 + sigmastar * gisq) * muepsi/ssq +
    S * giepsi/(ssqx2))
  sigx5 <- (0.5 * (rwv * ewu/sigmastar) + sigmastar)
  sigx6 <- (mu/(ssqx2) - sigx5 * muepsi/ssq)
  sigx7 <- (epsilon_isq/(2 * ewv)^2)
  sigZ1 <- sweep(szsq, MARGIN = 1, STATS = muepsi, FUN = "*") +
    sweep(Zigiepsi, MARGIN = 1, STATS = S/(ssqx2), FUN = "*")
  gradll <- cbind((sweep(Xgi, MARGIN = 1, STATS = -S * sigx1/sigmasq,
    FUN = "*") - sweep(Xepsi_i, MARGIN = 1, STATS = 1/(2 *
    ewv), FUN = "*")), sweep(muHvar, MARGIN = 1, STATS = (sigx2/sigmasq -
    (dmuu/wup + mu/ewu)), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = -((sigx4 * mup + 0.5 * gsig) * ewu - (0.5 * muwu +
      0.5 * (mu * dmuu/wup))), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = ((mup * sigx6 + 2 * sigx7 - 0.5/sigmasq) *
      ewv - 0.5 * (TT - 1)), FUN = "*"), -sweep(sweep(sigZ1,
    MARGIN = 1, STATS = mup, FUN = "*") + 0.5 * zsig, MARGIN = 1,
    STATS = ewu, FUN = "*"))
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
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phesstruncnormlike_bc92III <- function(parm, nXvar, nuZUvar,
  nvZVvar, nmuZUvar, uHvar, vHvar, muHvar, Yvar, Xvar, pindex,
  TT, S, wHvar, ngZGvar, gHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  eta <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + ngZGvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon_it <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsilon_isq <- as.numeric(tapply(epsilon_it^2, pindex[, 1],
    sum))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  Xgit <- sweep(-Xvar, MARGIN = 1, STATS = git, FUN = "*")
  Xgi <- apply(Xgit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xepsi_it <- sweep(-2 * Xvar, MARGIN = 1, STATS = epsilon_it,
    FUN = "*")
  Xepsi_i <- apply(Xepsi_it, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zitgitepsit <- sweep(gHvar, MARGIN = 1, STATS = git * epsilon_it,
    FUN = "*")
  Zigiepsi <- apply(Zitgitepsit, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Zitgitsq <- sweep(gHvar, MARGIN = 1, STATS = 2 * git^2, FUN = "*")
  Zigisq <- apply(Zitgitsq, 2, function(x) tapply(x, pindex[, 1],
    sum))
  Xzigi <- list()
  for (i in 1:ngZGvar) {
    Xzigi[[i]] <- apply(sweep(-Xvar, MARGIN = 1, STATS = gHvar[,
      i] * git, FUN = "*"), 2, function(x) tapply(x, pindex[, 1],
      sum))
  }
  Zisqgiepsi <- list()
  for (i in 1:ngZGvar) {
    Zisqgiepsi[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = gHvar[,
      i] * git * epsilon_it, FUN = "*"), 2, function(x) tapply(x,
      pindex[, 1], sum))
  }
  Zisqgisq <- list()
  for (i in 1:ngZGvar) {
    Zisqgisq[[i]] <- apply(sweep(gHvar, MARGIN = 1, STATS = 4 *
      gHvar[, i] * git^2, FUN = "*"), 2, function(x) tapply(x,
      pindex[, 1], sum))
  }
  Xsq <- list()
  for (i in 1:nXvar) {
    Xsq[[i]] <- apply(sweep(Xvar, MARGIN = 1, STATS = Xvar[,
      i], FUN = "*"), 2, function(x) tapply(x, pindex[, 1],
      sum))
  }
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv <- exp(Wv)
  sigmasq <- (ewu * gisq + ewv)
  sigmastar <- sqrt(ewu * ewv/sigmasq)
  ssq <- (sigmasq * sigmastar)^2
  ssqx2 <- sigmasq * sigmastar
  muepsi <- (mu * ewv - S * ewu * giepsi)
  musig <- muepsi/(ssqx2)
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dpmu <- dmusig/pmusig
  dmuu <- dnorm(mu/ewu_h, 0, 1)
  pmuu <- pnorm(mu/ewu_h)
  pmustar <- pmusig * sigmastar
  mup <- (musig + dpmu)
  wup <- (ewu_h * pmuu)
  gsig <- (gisq/sigmasq)
  muwu <- ((mu)^2/ewu)
  wuwv <- (ewu * ewv/(ssqx2))
  zsig <- sweep(Zigisq, MARGIN = 1, STATS = 1/sigmasq, FUN = "*")
  rwv <- (1 - ewv/sigmasq)
  ssuv <- (sigmastar - 0.5 * wuwv)
  szsq <- sweep(Zigisq, MARGIN = 1, STATS = ssuv/ssq, FUN = "*")
  sigx1 <- (muepsi/ewv + dmusig * ewu/(pmustar))
  sigx2 <- (muepsi/ewu + dmusig * ewv/(pmustar))
  sigx3 <- ((1 - ewu * gisq/sigmasq) * ewv/sigmastar)
  sigx4 <- ((0.5 * sigx3 + sigmastar * gisq) * muepsi/ssq +
    S * giepsi/(ssqx2))
  sigx5 <- (0.5 * (rwv * ewu/sigmastar) + sigmastar)
  sigx6 <- (mu/(ssqx2) - sigx5 * muepsi/ssq)
  sigx7 <- (epsilon_isq/(2 * ewv)^2)
  sigZ1 <- sweep(szsq, MARGIN = 1, STATS = muepsi, FUN = "*") +
    sweep(Zigiepsi, MARGIN = 1, STATS = S/(ssqx2), FUN = "*")
  sigx9 <- (muepsi/(pmustar) + dmusig * ewu * ewv/(pmustar)^2)
  sigx10 <- (sigx1 * dpmu - ewu/sigmastar)
  sigx11 <- (ewv/sigmastar - sigx2 * dpmu)
  sigx12 <- (0.5 * sigx3 + sigmastar * gisq)
  sigx13 <- (sigx12 * sigmasq * muepsi * sigmastar/ssq)
  wugsig <- (ewu * gisq/sigmasq)
  wugsig2 <- (1 - ewu * gisq/sigmasq)
  musq <- (mu^2/(ewu_h^3 * pmuu))
  sigx14 <- (0.5 * musq - (0.5 * wup - 0.5 * (mu * dmuu))/wup^2)
  sqsig <- (gisq/sigmasq^2)
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar +
    ngZGvar, ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar +
    ngZGvar)
  hessll[1:nXvar, 1:nXvar] <- -(crossprod(sweep(Xgi, MARGIN = 1,
    STATS = wHvar * S^2 * ((muepsi/(ewv * pmustar) + dmusig *
      ewu/(pmustar)^2) * dmusig/sigmasq - 1/ewv) * ewu/sigmasq,
    FUN = "*"), Xgi) + sapply(1:nXvar, function(x) {
    crossprod(Xsq[[x]], as.matrix(wHvar/ewv))
  }))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * (sigx9 * dmusig/sigmasq -
      1)/sigmasq, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = -wHvar *
    (S * (sigx10 * sigx4/sigmasq + mup * (1/(ssqx2) - sigx12 *
      ewu/ssq)) * ewu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xepsi_i,
    MARGIN = 1, STATS = wHvar * 2 * (ewv/(2 * ewv)^2), FUN = "*") +
    sweep(Xgi, MARGIN = 1, STATS = wHvar * S * (sigx10 *
      sigx6/sigmasq + mup * sigx5 * ewu/ssq) * ewv, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + ngZGvar)] <- -(crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * (sigx10/sigmasq * S * ewu),
    FUN = "*"), sigZ1) + sapply(1:ngZGvar, function(x) crossprod(Xzigi[[x]],
    as.matrix(wHvar * (mup/(ssqx2) * S * ewu)))) - crossprod(sweep(Xgi,
    MARGIN = 1, STATS = wHvar * (mup * ewu * S * ewu), FUN = "*"),
    szsq))
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    ((1/ewu - (muepsi/(ewu * pmustar) + dmusig * ewv/(pmustar)^2) *
      dmusig/sigmasq) * ewv/sigmasq + dmuu * (dmuu/wup^2 +
      mu/(ewu_h^3 * pmuu)) - 1/ewu), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx4 * sigx11/sigmasq +
      mup * sigx12 * ewv/ssq) * ewu - (0.5 * (((1 - mu^2/ewu_h^2)/wup -
      mu * dmuu/wup^2) * dmuu) + mu/ewu))), FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (mup * (1/(ssqx2) - sigx5 *
      ewv/ssq) + sigx11 * sigx6/sigmasq) * ewv, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar + ngZGvar)] <- -(crossprod(sweep(muHvar, MARGIN = 1,
    STATS = wHvar * sigx11/sigmasq * ewu, FUN = "*"), sigZ1) +
    crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar * mup *
      ewv * ewu, FUN = "*"), szsq))
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (((((mup * dpmu - 1) * sigx4^2 -
      0.5 * (gisq^2/sigmasq^2)) * ewu + ((sigx12 * (mu *
      ewv - (2 * sigx13 + 3 * (S * giepsi)) * ewu) + (0.5 *
      wugsig - 0.5 * (0.5 * wugsig2 + ewu * gisq/sigmasq)) *
      wugsig2 * ewv * muepsi/sigmastar)/ssq + S * giepsi/(ssqx2)) *
      mup + 0.5 * gsig) * ewu + 0.5 * muwu - 0.5 * (mu *
      sigx14 * dmuu))), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((mup * dpmu - 1) * sigx4 * sigx6 + 0.5 *
      sqsig - ((0.5 * (wugsig2 * ewv/sigmasq) + 0.5 * ((ewu *
      gisq/sigmasq - 1) * ewv/sigmasq + 1 - 0.5 * (wugsig2 *
      rwv))) * muepsi/sigmastar + mu * sigx12 - sigx5 *
      (2 * sigx13 + S * giepsi)) * mup/ssq) * ewu * ewv,
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar + ngZGvar)] <- -(crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (wugsig2/(ssqx2)) -
      0.5 * (1/(ssqx2) - sigx12 * ewu/ssq)) * ewv - 2 *
      (sigx12 * ssqx2 * ssuv/ssq)) * muepsi - S * ssuv *
      giepsi) * mup/ssq * ewu * ewu, FUN = "*"), Zigisq) -
    crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * S *
      sigx12 * mup/ssq * ewu * ewu, FUN = "*"), Zigiepsi) +
    crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * (mup *
      dpmu - 1) * sigx4 * ewu * ewu, FUN = "*"), sigZ1) -
    crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * 0.5 *
      (gisq/sigmasq^2) * ewu * ewu, FUN = "*"), Zigisq) +
    crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar * mup *
      ewu, FUN = "*"), sigZ1) + crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 0.5 * ewu, FUN = "*"), zsig))
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((1 - mup * dpmu) * sigx6^2 +
      0.5/sigmasq^2 - 16 * (ewv * epsilon_isq/(2 * ewv)^4)) *
      ewv + mup * (mu/(ssqx2) - (((3 * (mu) - 2 * (sigx5 *
      sigmasq * muepsi * sigmastar/ssq)) * ewv - S * ewu *
      giepsi) * sigx5 + (0.5 * (ewv/sigmasq) - 0.5 * (0.5 *
      rwv + ewv/sigmasq)) * rwv * ewu * muepsi/sigmastar)/ssq) +
      2 * sigx7 - 0.5/sigmasq) * ewv, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + nvZVvar + ngZGvar)] <- -(crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (rwv/(ssqx2)) -
      0.5 * (1/(ssqx2) - sigx5 * ewv/ssq)) * ewu - 2 *
      (sigx5 * ssqx2 * ssuv/ssq)) * muepsi + mu * ssuv) *
      mup/ssq * ewu * ewv, FUN = "*"), Zigisq) - crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * S * sigx5 * mup/ssq * ewu *
      ewv, FUN = "*"), Zigiepsi) + crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - mup * dpmu) * sigx6 *
      ewu * ewv, FUN = "*"), sigZ1) - crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 0.5/sigmasq^2 * ewu * ewv,
    FUN = "*"), Zigisq))
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar + ngZGvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar + ngZGvar)] <- -(crossprod(sweep(Zigisq, MARGIN = 1,
    STATS = wHvar * ((0.5 * (ssuv/ssq) - 0.5/(sigmasq^2 *
      sigmastar)) * ewu * ewv - 2 * (ssqx2 * ssuv^2/ssq)) *
      muepsi * ewu/ssq * mup * ewu, FUN = "*"), Zigisq) -
    crossprod(sweep(Zigiepsi, MARGIN = 1, STATS = wHvar *
      S * ssuv * ewu/ssq * mup * ewu, FUN = "*"), Zigisq) +
    sapply(1:ngZGvar, function(x) crossprod(Zisqgisq[[x]],
      as.matrix(wHvar * muepsi * ssuv/ssq * mup * ewu))) +
    S * (sapply(1:ngZGvar, function(x) crossprod(Zisqgiepsi[[x]],
      as.matrix(wHvar * mup/(ssqx2) * ewu))) - crossprod(sweep(Zigisq,
      MARGIN = 1, STATS = wHvar * ewu * ssuv/ssq * mup *
        ewu, FUN = "*"), Zigiepsi)) + crossprod(sweep(sigZ1,
    MARGIN = 1, STATS = wHvar * (mup * dpmu - 1) * ewu *
      ewu, FUN = "*"), sigZ1) + 0.5 * (sapply(1:ngZGvar,
    function(x) crossprod(Zisqgisq[[x]], as.matrix(wHvar *
      ewu/sigmasq))) - crossprod(sweep(Zigisq, MARGIN = 1,
    STATS = wHvar * ewu/sigmasq/sigmasq * ewu, FUN = "*"),
    Zigisq)))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncatednormal-normal distribution
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
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
truncnormAlgOpt_bc92III <- function(start, olsParam, dataTable,
  S, nXvar, muHvar_c, muHvar_p, nmuZUvar, uHvar_c, uHvar_p,
  gHvar, ngZGvar, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar,
  Xvar, wHvar_c, wHvar_p, pindex, TT, method, printInfo, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- psttruncnorm_bc92III(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      gHvar = gHvar, ngZGvar = ngZGvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_c, uHvar = uHvar_c, vHvar = vHvar_c,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c,
      itermax = itermax, tol = tol, printInfo = printInfo)
    Inittrunc <- start_st$inittrunc
    startVal <- start_st$StartVal
  }
  startLoglik <- sum(ptruncnormlike_bc92III(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_bc92III(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_bc92III(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ptruncnormlike_bc92III,
    grad = pgradtruncnormlike_bc92III, hess = phesstruncnormlike_bc92III,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
    Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S,
    wHvar = wHvar_p), sr1 = trust.optim(x = startVal, fn = function(parm) {
    -sum(ptruncnormlike_bc92III(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtruncnormlike_bc92III(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
    stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
    report.precision = 1L)), sparse = trust.optim(x = startVal,
    fn = function(parm) {
      -sum(ptruncnormlike_bc92III(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgradtruncnormlike_bc92III(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phesstruncnormlike_bc92III(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
        ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S,
        wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) {
    -sum(ptruncnormlike_bc92III(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgradtruncnormlike_bc92III(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hess = function(parm) {
    -phesstruncnormlike_bc92III(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p)
  }, print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) {
    -sum(ptruncnormlike_bc92III(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gradient = function(parm) {
    -colSums(pgradtruncnormlike_bc92III(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar_p, uHvar = uHvar_p, gHvar = gHvar,
      ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
      Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = function(parm) {
    -phesstruncnormlike_bc92III(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar_p,
      uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar,
      vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p)
  }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgradtruncnormlike_bc92III(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
      gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
      Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
      S = S, wHvar = wHvar_p))
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
      mleObj$hessian <- phesstruncnormlike_bc92III(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phesstruncnormlike_bc92III(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
        gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT,
        S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- ptruncnormlike_bc92III(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar_p, gHvar = gHvar, ngZGvar = ngZGvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgradtruncnormlike_bc92III(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar_p, uHvar = uHvar_p,
    gHvar = gHvar, ngZGvar = ngZGvar, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  if (is.null(start)) {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam, Inittrunc = Inittrunc))
  } else {
    return(list(startVal = startVal, startLoglik = startLoglik,
      mleObj = mleObj, mlParam = mlParam))
  }
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncatednormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ptruncnormeff_bc92III <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar)]
  eta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar +
    object$nuZUvar + object$nvZVvar + object$ngZGvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar_c <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  gHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 5)
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
  epsilon_it <- model.response(model.frame(object$formula,
    data = object$dataTable)) - as.numeric(crossprod(matrix(beta),
    t(Xvar)))
  git <- exp(as.numeric(crossprod(matrix(eta), t(gHvar))))
  git_epsit <- epsilon_it * git
  giepsi <- as.numeric(tapply(git_epsit, pindex[, 1], sum))
  gisq <- as.numeric(tapply(git^2, pindex[, 1], sum))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * giepsi)/(exp(Wv) +
    gisq * exp(Wu))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wv) + gisq * exp(Wu)))
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  res <- data.frame(levels(pindex[, 1]), u = u, uLB = uLB,
    uUB = uUB, m = m, mustar = mustar, sigmastar = sigmastar)
  names(res)[1] <- names(pindex)[1]
  res <- merge(pindex, res, by = names(pindex)[1])
  res$u <- res$u * git
  res$m <- res$m * git
  res$uLB <- res$uLB * git
  res$uUB <- res$uUB * git
  if (object$logDepVar == TRUE) {
    res$teJLMS <- exp(-res$u)
    res$teMO <- exp(-res$m)
    res$teBC <- exp(1/2 * res$sigmastar^2 * git^2 - res$mustar *
      git) * pnorm(res$mustar/res$sigmastar - res$sigmastar *
      git)/pnorm(res$mustar/res$sigmastar)
    res$teBCLB <- exp(-res$uUB)
    res$teBCUB <- exp(-res$uLB)
    res$teBC_reciprocal <- exp(1/2 * res$sigmastar^2 * git^2 +
      res$mustar * git) * pnorm(res$mustar/res$sigmastar +
      res$sigmastar * git)/pnorm(res$mustar/res$sigmastar)
  }
  res$mustar <- NULL
  res$sigmastar <- NULL
  res <- pdata.frame(res, names(res)[1:2])
  return(res)
}


