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
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
puninormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
#' @param wHvar vector of weights (weighted likelihood)
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
pstuninorm_mbc92 <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, uHvar, vHvar,
  Yvar, Xvar, wHvar, S, printInfo, tol, whichStart, initIter, initAlg) {
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
  }, 0.001, 0.001)
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "eta1", "eta2")
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
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
pgraduninormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  sigx1 <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(sigx1/sqrt(ewv/gisq), 0, 1)
  pmusig1 <- pnorm(sigx1/sqrt(ewv/gisq))
  dmusig2 <- dnorm(S * giepsi/(sqrt(ewv/gisq) * gisq), 0, 1)
  pmusig2 <- pnorm(S * giepsi/(sqrt(ewv/gisq) * gisq))
  sigx2 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq))
  sigx3 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq) * gisq)
  sigx4 <- (0.5 * (S * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.5 *
    (sigx1 * dmusig1))
  sigx5 <- (0.5 * (sigx1 * giZi) + S * (Zi_epsi - giepsi * giZi/gisq))
  sigx6 <- (sqrt(ewv/gisq) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)))
  sigx7 <- (0.5 * (sigx1 * giZisq) + S * (Zisq_epsi - giepsi * giZisq/gisq))
  sigx8 <- (Zi_epsi/(sqrt(ewv/gisq) * gisq) - sigx6 * giepsi * giZi/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx9 <- (sigx5 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx8)
  sigx10 <- ((-(S * giepsi/gisq))^2 * giZi + 2 * (S^2 * giepsi * (Zi_epsi - giepsi *
    giZi/gisq)/gisq))
  sigx11 <- (Zisq_epsi/(sqrt(ewv/gisq) * gisq) - sigx6 * giepsi * giZisq/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx12 <- ((-(S * giepsi/gisq))^2 * giZisq + 2 * (S^2 * giepsi * (Zisq_epsi -
    giepsi * giZisq/gisq)/gisq))
  sigx13 <- (sigx7 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx11)
  gradll <- cbind(sweep(Xgi, MARGIN = 1, STATS = S * (dmusig1 - dmusig2)/sigx3,
    FUN = "*") - 0.5 * (sweep(Xepsi_i, MARGIN = 1, STATS = 1/ewv, FUN = "*") -
    sweep(Xgi, MARGIN = 1, STATS = 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*")),
    sweep(uHvar, MARGIN = 1, STATS = (sqrt(12)/2 * (dmusig1 * ewu_h/sigx2) -
      0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx4/sigx2 + 0.5 +
      0.5 * ((epsilon_isq - (-(S * giepsi/gisq))^2 * gisq)/ewv) - 0.5 * TT),
      FUN = "*"), sigx9/(pmusig1 - pmusig2) + 0.5 * (sigx10/ewv) - 0.5 * (giZi/gisq),
    sigx13/(pmusig1 - pmusig2) + 0.5 * (sigx12/ewv) - 0.5 * (giZisq/gisq))
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
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param pindex panel indices (ID, TIME)
#' @param TT vector of time of presence
#' @noRd
phessuninormlike_mbc92 <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
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
  sigx1 <- (sqrt(12) * ewu_h + S * giepsi/gisq)
  dmusig1 <- dnorm(sigx1/sqrt(ewv/gisq), 0, 1)
  pmusig1 <- pnorm(sigx1/sqrt(ewv/gisq))
  dmusig2 <- dnorm(S * giepsi/(sqrt(ewv/gisq) * gisq), 0, 1)
  pmusig2 <- pnorm(S * giepsi/(sqrt(ewv/gisq) * gisq))
  sigx2 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq))
  sigx3 <- ((pmusig1 - pmusig2) * sqrt(ewv/gisq) * gisq)
  sigx4 <- (0.5 * (S * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.5 *
    (sigx1 * dmusig1))
  sigx5 <- (0.5 * (sigx1 * giZi) + S * (Zi_epsi - giepsi * giZi/gisq))
  sigx6 <- (sqrt(ewv/gisq) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)))
  sigx7 <- (0.5 * (sigx1 * giZisq) + S * (Zisq_epsi - giepsi * giZisq/gisq))
  sigx8 <- (Zi_epsi/(sqrt(ewv/gisq) * gisq) - sigx6 * giepsi * giZi/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx9 <- (sigx5 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx8)
  sigx10 <- ((-(S * giepsi/gisq))^2 * giZi + 2 * (S^2 * giepsi * (Zi_epsi - giepsi *
    giZi/gisq)/gisq))
  sigx11 <- (Zisq_epsi/(sqrt(ewv/gisq) * gisq) - sigx6 * giepsi * giZisq/(sqrt(ewv/gisq) *
    gisq)^2)
  sigx12 <- ((-(S * giepsi/gisq))^2 * giZisq + 2 * (S^2 * giepsi * (Zisq_epsi -
    giepsi * giZisq/gisq)/gisq))
  sigx13 <- (sigx7 * dmusig1/(sqrt(ewv/gisq) * gisq) - S * dmusig2 * sigx11)
  sigx14 <- (Zisq_epsi - giepsi * giZisq/gisq)
  sigx15 <- (Zicub - giZi * giZisq/gisq)
  sigx16 <- (0.5 * (sigx1^2/ewv) - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)^2))
  sigx17 <- ((0.5/gisq - 0.5 * (1/gisq - 0.5 * (ewv/(sqrt(ewv/gisq) * gisq)^2)))/sqrt(ewv/gisq) -
    sigx6 * gisq/(sqrt(ewv/gisq) * gisq)^2)
  sigx18 <- (Zi_epsi - giepsi * giZi/gisq)
  sigx19 <- (Zifour - giZisq^2/gisq)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 2, ncol = nXvar + nuZUvar +
    nvZVvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar *
    ((S * dmusig2 * giepsi/gisq - sigx1 * dmusig1)/(ewv * (pmusig1 - pmusig2) *
      sqrt(ewv/gisq) * gisq) - (dmusig1 - dmusig2)^2/sigx3^2), FUN = "*"),
    Xgi) - 0.5 * (sapply(1:nXvar, function(x) crossprod(Xsq[[x]], as.matrix(wHvar/ewv))) -
    crossprod(sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 * (S^2/gisq)/ewv, FUN = "*"),
      Xgi))
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xgi, MARGIN = 1,
    STATS = -wHvar * (sqrt(12)/2 * (S * (sigx1/(ewv * (pmusig1 - pmusig2) * sqrt(ewv/gisq)) +
      (dmusig1 - dmusig2)/(sigx2^2 * gisq)) * dmusig1 * ewu_h)), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod((0.5 *
    ((sweep(Xepsi_i, MARGIN = 1, STATS = wHvar/ewv, FUN = "*") - sweep(Xgi, MARGIN = 1,
      STATS = wHvar * 2 * (S^2 * giepsi/gisq)/ewv, FUN = "*"))) + sweep(Xgi,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (dmusig2 * (ewv - S^2 * giepsi^2/gisq)/(sqrt(ewv/gisq) *
      gisq)^2) - 0.5 * ((1/gisq - sigx1^2/ewv) * dmusig1))/sigx2 - sigx4 *
      (dmusig1 - dmusig2)/(sigx2^2 * gisq)), FUN = "*")), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- colSums((S * ((((sweep(XiZi,
    MARGIN = 1, STATS = wHvar * dmusig1/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * (sigx5 * sigx1/ewv +
    0.5 * (giZi/gisq)) * dmusig1/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*")) - sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx9 * (dmusig1 - dmusig2)/(pmusig1 -
    pmusig2)/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2), FUN = "*")) - (sweep(XiZi,
    MARGIN = 1, STATS = wHvar/(sqrt(ewv/gisq) * gisq) * dmusig2/(pmusig1 - pmusig2),
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * (sigx6 * giZi/(sqrt(ewv/gisq) *
    gisq)^2 + S^2 * giepsi * sigx8/(ewv * gisq)) * dmusig2/(pmusig1 - pmusig2),
    FUN = "*"))) + 0.5 * (S * (2 * ((sweep(XiZi, MARGIN = 1, STATS = wHvar *
    giepsi/(ewv * gisq), FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar *
    giZi/gisq * giepsi/(ewv * gisq), FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    sigx18/(ewv * gisq), FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar *
    2 * (giepsi * giZi/gisq)/(ewv * gisq), FUN = "*"))))))
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2)] <- colSums((S * ((((sweep(XiZisq,
    MARGIN = 1, STATS = wHvar * dmusig1/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * (sigx7 * sigx1/ewv +
    0.5 * (giZisq/gisq)) * dmusig1/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*")) - sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx13 * (dmusig1 -
    dmusig2)/(pmusig1 - pmusig2)/(sqrt(ewv/gisq) * gisq)/(pmusig1 - pmusig2),
    FUN = "*")) - (sweep(XiZisq, MARGIN = 1, STATS = wHvar/(sqrt(ewv/gisq) *
    gisq) * dmusig2/(pmusig1 - pmusig2), FUN = "*") - sweep(Xgi, MARGIN = 1,
    STATS = wHvar * (sigx6 * giZisq/(sqrt(ewv/gisq) * gisq)^2 + S^2 * giepsi *
      sigx11/(ewv * gisq)) * dmusig2/(pmusig1 - pmusig2), FUN = "*"))) + 0.5 *
    (S * (2 * ((sweep(XiZisq, MARGIN = 1, STATS = wHvar * giepsi/(ewv * gisq),
      FUN = "*") - sweep(Xgi, MARGIN = 1, STATS = wHvar * giZisq/gisq * giepsi/(ewv *
      gisq), FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar * sigx14/(ewv *
      gisq), FUN = "*")) + sweep(Xgi, MARGIN = 1, STATS = wHvar * 2 * (giepsi *
      giZisq/gisq)/(ewv * gisq), FUN = "*"))))))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * sqrt(12)/2 * (((0.5 - sqrt(12)/2 * (sigx1 * ewu_h *
      gisq/ewv))/sigx2 - sqrt(12)/2 * (dmusig1 * ewu_h/sigx2^2)) * dmusig1 *
      ewu_h), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * ((0.5 *
    ((sqrt(12)/2 - sqrt(12)/2 * (sigx1^2 * gisq/ewv))/sigx2) + sqrt(12)/2 * (sigx4/sigx2^2)) *
    dmusig1 * ewu_h), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/4 * giZi - sqrt(12)/2 * (sigx5 * sigx1 *
      gisq/ewv))/gisq - sqrt(12)/2 * (sigx9/(pmusig1 - pmusig2))) * dmusig1 *
      ewu_h/sigx2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sqrt(12)/4 * giZisq - sqrt(12)/2 * (sigx7 *
      sigx1 * gisq/ewv))/gisq - sqrt(12)/2 * (sigx13/(pmusig1 - pmusig2))) *
      dmusig1 * ewu_h/sigx2, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (S * ((0.5 * (S^2 * giepsi^2) - ewv * gisq)/(sqrt(ewv/gisq) * gisq)^2 +
      1) * dmusig2 * ewv * giepsi/(sqrt(ewv/gisq) * gisq)^2) - 0.25 * (sigx1^3 *
      dmusig1 * gisq/ewv))/sigx2 - (((0.5 * ((pmusig1 - pmusig2)/(sqrt(ewv/gisq) *
      gisq)) + 0.5 * (S * dmusig2 * giepsi/(sqrt(ewv/gisq) * gisq)^2)) * ewv -
      0.5 * (sigx1 * dmusig1)) * sigx4/sigx2^2 + 0.5 * ((epsilon_isq - (-(S *
      giepsi/gisq))^2 * gisq)/ewv))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 1)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (((sigx5 *
    sigx16 * dmusig1 - sigx9 * sigx4/(pmusig1 - pmusig2))/sqrt(ewv/gisq) - S *
    (0.5 * (S^2 * giepsi^2 * sigx8) - (sigx17 * giepsi * giZi + 0.5 * (Zi_epsi/sqrt(ewv/gisq))) *
      ewv) * dmusig2/(sqrt(ewv/gisq) * gisq)^2)/(pmusig1 - pmusig2) - 0.5 *
    (sigx10/ewv)), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar + nuZUvar +
    nvZVvar + 2)] <- colSums(sweep(vHvar, MARGIN = 1, STATS = wHvar * (((sigx7 *
    sigx16 * dmusig1 - sigx13 * sigx4/(pmusig1 - pmusig2))/sqrt(ewv/gisq) - S *
    (0.5 * (S^2 * giepsi^2 * sigx11) - (sigx17 * giepsi * giZisq + 0.5 * (Zisq_epsi/sqrt(ewv/gisq))) *
      ewv) * dmusig2/(sqrt(ewv/gisq) * gisq)^2)/(pmusig1 - pmusig2) - 0.5 *
    (sigx12/ewv)), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 1)] <- sum(wHvar *
    ((((0.5 * (sigx1 * Zisq + S * giZi * sigx18/gisq) - (sigx5^2 * sigx1/ewv +
      S * (giepsi * (Zisq - giZi^2/gisq) + giZi * Zi_epsi)/gisq))/(sqrt(ewv/gisq) *
      gisq) - sigx5 * sigx6 * giZi/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 + S *
      ((((0.5 * (sigx6/(sqrt(ewv/gisq) * gisq)^2) - 0.5/(sqrt(ewv/gisq) * gisq^2)) *
        ewv * giepsi * giZi + sigx6 * Zi_epsi) * giZi + sigx6 * (giepsi *
        (Zisq - 2 * (sqrt(ewv/gisq) * sigx6 * gisq * giZi^2/(sqrt(ewv/gisq) *
          gisq)^2)) + giZi * Zi_epsi))/(sqrt(ewv/gisq) * gisq)^2 + S^2 *
        giepsi * sigx8^2/(sqrt(ewv/gisq) * gisq)) * dmusig2 - sigx9^2/(pmusig1 -
      pmusig2))/(pmusig1 - pmusig2) + 0.5 * (((-(S * giepsi/gisq))^2 * Zisq +
      S^2 * (2 * (giepsi * giZi * sigx18/gisq) + 2 * (Zi_epsi * sigx18 - ((2 *
        Zi_epsi - giepsi * giZi/gisq) * giZi + giepsi * (Zisq - giZi^2/gisq)) *
        giepsi/gisq))/gisq)/ewv) - 0.5 * ((Zisq - giZi^2/gisq)/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    ((((0.5 * (sigx1 * Zicub + S * giZisq * sigx18/gisq) - (sigx5 * sigx7 * sigx1/ewv +
      S * (giepsi * sigx15 + giZisq * Zi_epsi)/gisq))/(sqrt(ewv/gisq) * gisq) -
      sigx7 * sigx6 * giZi/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 + S * ((((0.5 *
      (sigx6/(sqrt(ewv/gisq) * gisq)^2) - 0.5/(sqrt(ewv/gisq) * gisq^2)) *
      ewv * giepsi * giZi + sigx6 * Zi_epsi) * giZisq + sigx6 * (giepsi * (Zicub -
      2 * (sqrt(ewv/gisq) * sigx6 * gisq * giZi * giZisq/(sqrt(ewv/gisq) *
        gisq)^2)) + giZi * Zisq_epsi))/(sqrt(ewv/gisq) * gisq)^2 + S^2 *
      giepsi * sigx8 * sigx11/(sqrt(ewv/gisq) * gisq)) * dmusig2 - sigx9 *
      sigx13/(pmusig1 - pmusig2))/(pmusig1 - pmusig2) + 0.5 * (((-(S * giepsi/gisq))^2 *
      Zicub + S^2 * (2 * (giepsi * giZisq * sigx18/gisq) + 2 * (Zi_epsi * sigx14 -
      giepsi * (giepsi * sigx15 + giZi * sigx14 + giZisq * Zi_epsi)/gisq))/gisq)/ewv) -
      0.5 * (sigx15/gisq)))
  hessll[(nXvar + nuZUvar + nvZVvar + 2), (nXvar + nuZUvar + nvZVvar + 2)] <- sum(wHvar *
    ((((0.5 * (sigx1 * Zifour + S * giZisq * sigx14/gisq) - (sigx7^2 * sigx1/ewv +
      S * (giepsi * sigx19 + giZisq * Zisq_epsi)/gisq))/(sqrt(ewv/gisq) * gisq) -
      sigx7 * sigx6 * giZisq/(sqrt(ewv/gisq) * gisq)^2) * dmusig1 + S * ((((0.5 *
      (sigx6/(sqrt(ewv/gisq) * gisq)^2) - 0.5/(sqrt(ewv/gisq) * gisq^2)) *
      ewv * giepsi * giZisq + sigx6 * Zisq_epsi) * giZisq + sigx6 * (giepsi *
      (Zifour - 2 * (sqrt(ewv/gisq) * sigx6 * gisq * giZisq^2/(sqrt(ewv/gisq) *
        gisq)^2)) + giZisq * Zisq_epsi))/(sqrt(ewv/gisq) * gisq)^2 + S^2 *
      giepsi * sigx11^2/(sqrt(ewv/gisq) * gisq)) * dmusig2 - sigx13^2/(pmusig1 -
      pmusig2))/(pmusig1 - pmusig2) + 0.5 * (((-(S * giepsi/gisq))^2 * Zifour +
      S^2 * (2 * (giepsi * giZisq * sigx14/gisq) + 2 * (Zisq_epsi * sigx14 -
        ((2 * Zisq_epsi - giepsi * giZisq/gisq) * giZisq + giepsi * sigx19) *
          giepsi/gisq))/gisq)/ewv) - 0.5 * (sigx19/gisq)))
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
uninormAlgOpt_mbc92 <- function(start, randStart, sdStart, olsParam, dataTable, S,
  nXvar, uHvar_c, uHvar_p, nuZUvar, vHvar_c, vHvar_p, nvZVvar, Yvar, Xvar, wHvar_c,
  wHvar_p, pindex, TT, method, printInfo, whichStart, initIter, initAlg, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- pstuninorm_mbc92(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_c,
      vHvar = vHvar_c, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar_c, tol = tol,
      whichStart = whichStart, initIter = initIter, initAlg = initAlg, printInfo = printInfo)
    Inituni <- start_st$inituni
    startVal <- start_st$StartVal
  }
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(puninormlike_mbc92(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
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
    -sum(puninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, gr = function(parm) {
    -colSums(pgraduninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
      TT = TT, S = S, wHvar = wHvar_p))
  }, hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
    stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = puninormlike_mbc92,
    grad = pgraduninormlike_mbc92, hess = phessuninormlike_mbc92, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
    Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) {
      -sum(puninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, method = "SR1", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) {
      -sum(puninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hs = function(parm) {
      as(-phessuninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p), "dgCMatrix")
    }, method = "Sparse", control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
      prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) {
      -sum(puninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gr = function(parm) {
      -colSums(pgraduninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hess = function(parm) {
      -phessuninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) {
      -sum(puninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p))
    }, gradient = function(parm) {
      -colSums(pgraduninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar,
        Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p))
    }, hessian = function(parm) {
      -phessuninormlike_mbc92(parm, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar, pindex = pindex,
        TT = TT, S = S, wHvar = wHvar_p)
    }, control = list(iter.max = itermax, trace = if (printInfo) 1 else 0, eval.max = itermax,
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(pgraduninormlike_mbc92(mleObj$par, nXvar = nXvar,
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
      mleObj$hessian <- phessuninormlike_mbc92(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
    if (method == "sr1") {
      mleObj$hessian <- phessuninormlike_mbc92(parm = mleObj$solution, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p,
        Yvar = Yvar, Xvar = Xvar, pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
    }
  }
  mleObj$logL_OBS <- puninormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
  mleObj$gradL_OBS <- pgraduninormlike_mbc92(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar_p, vHvar = vHvar_p, Yvar = Yvar, Xvar = Xvar,
    pindex = pindex, TT = TT, S = S, wHvar = wHvar_p)
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
puninormeff_mbc92 <- function(object, level) {
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
  git <- unlist(lapply(TT, FUN = function(x) 1 + eta1 * rev(-(0:(x - 1))) + eta2 *
    (rev(-(0:(x - 1))))^2))
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
