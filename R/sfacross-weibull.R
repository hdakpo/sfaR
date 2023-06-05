################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: weibull - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for weibull-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cweibullnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  if (k < 0)
    return(NA)
  urMat <- sweep((-log(1 - FiMat))^(1/k), MARGIN = 1, STATS = exp(Wu/2),
    FUN = "*")
  uepsi <- sweep(S * urMat, MARGIN = 1, STATS = epsilon, FUN = "+")
  duepsi <- dnorm(sweep(uepsi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/"))
  ll <- log(apply(sweep(duepsi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/"), 1, mean))
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for weibull-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
cstweibullnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/2)))^(2/3)
  } else {
    varu <- (-S * m3/2)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1/2 * log((epsiRes^2 - varv)^2)
  dep_v <- 1/2 * log((epsiRes^2 - varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetu$coefficients)))
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetv$coefficients)))
    stop("at least one of the OLS coefficients of 'vhet' is NA: ",
      paste(colnames(vHvar)[is.na(reg_hetv$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi, k = 1))
}

# Gradient of the likelihood function ----------
#' gradient for weibull-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cgradweibullnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  gradll <- matrix(nrow = N, ncol = nXvar + nuZUvar + nvZVvar +
    1)
  lFimat <- (-log(1 - FiMat))^(1/k)
  lFiu <- sweep(S * lFimat, MARGIN = 1, STATS = exp(Wu/2),
    FUN = "*")
  lFiuepsi <- sweep(lFiu, MARGIN = 1, STATS = epsilon, FUN = "+")
  dFimat <- dnorm(sweep(lFiuepsi, MARGIN = 1, STATS = 1/exp(Wv/2),
    FUN = "*"))
  dFiv <- sweep(dFimat, MARGIN = 1, STATS = 1/exp(Wv/2), FUN = "*")
  lFi1 <- dFimat * lFiuepsi
  lFi2 <- log(-log(1 - FiMat))
  lFi3 <- lFimat * lFi1
  sigx1 <- sweep(lFi1, MARGIN = 1, STATS = 1/exp(Wv/2)^3, FUN = "*")
  sigx2 <- sweep(S * lFi2 * lFi3, MARGIN = 1, STATS = exp(Wu/2)/(k^2 *
    exp(Wv/2)^3), FUN = "*")
  sigx3 <- sweep(lFi3, MARGIN = 1, STATS = exp(Wu/2)/exp(Wv/2)^3,
    FUN = "*")
  sigx4 <- sweep(sweep(0.5 * lFi1 * lFiuepsi, MARGIN = 1, STATS = 1/exp(Wv/2)^2,
    FUN = "*") - 0.5 * dFimat, MARGIN = 1, STATS = 1/exp(Wv/2),
    FUN = "*")
  sdFiv <- apply(dFiv, 1, sum)
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    gx[, k] <- apply(sweep(sigx1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    gu[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = -(0.5 *
      (S * uHvar[, k])), FUN = "*"), 1, sum)/sdFiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    gv[, k] <- apply(sweep(sigx4, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gradll <- cbind(gx, gu, gv, apply(sigx2, 1, sum)/sdFiv)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for weibull-normal distribution
#' @param parm all parameters to be estimated
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
chessweibullnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv_h <- exp(Wv/2)
  Q <- dim(FiMat)[2]
  lFimat <- (-log(1 - FiMat))^(1/k)
  lFiu <- sweep(S * lFimat, MARGIN = 1, STATS = ewu_h, FUN = "*")
  lFiuepsi <- sweep(lFiu, MARGIN = 1, STATS = epsilon, FUN = "+")
  dFimat <- dnorm(sweep(lFiuepsi, MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*"))
  lFi1sq <- dFimat * lFiuepsi^2
  lFi1 <- dFimat * lFiuepsi
  lFi2 <- log(-log(1 - FiMat))
  lFid <- lFimat * dFimat
  sigx1 <- sweep(lFi1, MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  sigx2 <- sweep(S * lFid * lFi2 * lFiuepsi, MARGIN = 1, STATS = ewu_h/(k^2 *
    ewv_h^3), FUN = "*")
  sigx3 <- sweep(sweep((lFi1sq), MARGIN = 1, STATS = 0.5/ewv_h^2,
    FUN = "*") - 0.5 * dFimat, MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  sigx4 <- sweep(lFid * lFiuepsi, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  dFiv <- sweep(dFimat, MARGIN = 1, STATS = 1/ewv_h, FUN = "*")
  sdFiv <- apply(dFiv, 1, sum)
  sigx5 <- sweep((lFiuepsi^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*") - 1
  sigx6 <- sweep(sigx5 * dFimat, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  sigx7 <- sweep((lFimat * lFiuepsi^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*") - lFimat
  sigx8 <- sweep(sigx7 * dFimat * ewu_h, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  sigx9 <- sweep((lFiuepsi^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*") - 2
  sigx10 <- sweep(dFimat * lFi2, MARGIN = 1, STATS = ewu_h/(k^2 *
    ewv_h^3), FUN = "*")
  sigx11 <- sweep(0.5 * (S * (-log(1 - FiMat))^(2/k) * lFiuepsi),
    MARGIN = 1, STATS = ewu_h/ewv_h^2, FUN = "*")
  sigx12 <- sweep(0.5 * (S * (-log(1 - FiMat))^(2/k)), MARGIN = 1,
    STATS = ewu_h, FUN = "*")
  sigx13 <- sweep(dFimat, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  sigx14 <- sweep((lFimat * lFiuepsi^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*")
  sigx15 <- sweep(dFimat * lFiuepsi, MARGIN = 1, STATS = ewu_h/ewv_h^3,
    FUN = "*")
  sigx16 <- sweep((lFiuepsi^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*")
  sigx17 <- sweep(lFi1sq, MARGIN = 1, STATS = 1/ewv_h^2, FUN = "*")
  sigx18 <- sweep(((0.5 * (0.5 * sigx16 - 1) - 0.25) * sigx17 -
    0.5 * (0.5 * (sigx17) - 0.5 * dFimat)), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  sigx19 <- sweep((lFimat * lFiuepsi^2), MARGIN = 1, STATS = 1/(k^2 *
    ewv_h^5), FUN = "*")
  sigx20 <- sweep((k^2 * lFimat), MARGIN = 1, STATS = ewv_h^3/(k^2 *
    ewv_h^3)^2, FUN = "*")
  sigx21 <- sweep((0.5 * sigx19 - 1.5 * sigx20) * dFimat *
    lFi2 * lFiuepsi, MARGIN = 1, STATS = ewu_h, FUN = "*")
  sigx22 <- sweep(S * (-log(1 - FiMat))^(2/k) * lFiuepsi, MARGIN = 1,
    STATS = ewu_h/ewv_h^2, FUN = "*")
  sigx23 <- sweep(S * (-log(1 - FiMat))^(2/k), MARGIN = 1,
    STATS = ewu_h, FUN = "*")
  sigx24 <- sweep((lFiuepsi * (sigx22 - lFimat) - sigx23) *
    lFi2, MARGIN = 1, STATS = 1/(k^4 * ewv_h^3), FUN = "*")
  sigx25 <- sweep((k * lFimat * lFiuepsi), MARGIN = 1, STATS = ewv_h^3/(k^2 *
    ewv_h^3)^2, FUN = "*")
  sigx26 <- sweep(S * (sigx24 - 2 * sigx25) * dFimat * lFi2,
    MARGIN = 1, STATS = ewu_h, FUN = "*")
  X1 <- matrix(nrow = N, ncol = nXvar)
  X7 <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    X1[, k] <- apply(sweep(sigx1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)
    X7[, k] <- apply(sweep(S * sigx7 * sigx10, MARGIN = 1,
      STATS = Xvar[, k]/sdFiv, FUN = "*"), 1, sum)
  }
  Zu4 <- matrix(nrow = N, ncol = nuZUvar)
  Zu10 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    Zu4[, k] <- apply(sweep(-(0.5 * (S * sigx4)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
    Zu10[, k] <- apply(sweep(S * ((0.5 * lFimat - sigx11) *
      lFiuepsi + sigx12) * sigx10, MARGIN = 1, STATS = uHvar[,
      k]/sdFiv, FUN = "*"), 1, sum)
  }
  Zv3 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    Zv3[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  Zv21 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    Zv21[, k] <- apply(sweep(S * sigx21, MARGIN = 1, STATS = vHvar[,
      k]/sdFiv, FUN = "*"), 1, sum)
  }
  X6 <- list()
  Xu8 <- list()
  Xv9 <- list()
  Zu11 <- list()
  Zuv14 <- list()
  Zv18 <- list()
  for (r in seq_len(Q)) {
    X6[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sigx6[, r]/sdFiv, FUN = "*"), Xvar)
    Xu8[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * sigx8))[, r]/sdFiv, FUN = "*"), uHvar)
    Xv9[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      ((0.5 * sigx9 - 0.5) * sigx1)[, r]/sdFiv, FUN = "*"),
      vHvar)
    Zu11[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * ((0.5 * lFimat - sigx11) * lFiuepsi +
        sigx12) * sigx13))[, r]/sdFiv, FUN = "*"), uHvar)
    Zuv14[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      (S * (0.25 * lFimat + 0.5 * (lFimat - 0.5 * sigx14)) *
        sigx15)[, r]/sdFiv, FUN = "*"), vHvar)
    Zv18[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      sigx18[, r]/sdFiv, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", X6) - crossprod(sweep(X1,
    MARGIN = 1, STATS = wHvar/sdFiv^2, FUN = "*"), X1)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    Xu8) - crossprod(sweep(X1, MARGIN = 1, STATS = wHvar/sdFiv^2,
    FUN = "*"), Zu4)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", Xv9) - crossprod(sweep(X1, MARGIN = 1,
    STATS = wHvar/sdFiv^2, FUN = "*"), Zv3)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(X7,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(X1, MARGIN = 1,
    STATS = wHvar * apply(sigx2, 1, sum)/sdFiv^2, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", Zu11) - crossprod(sweep(Zu4,
    MARGIN = 1, STATS = wHvar/sdFiv^2, FUN = "*"), Zu4)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", Zuv14) -
    crossprod(sweep(Zu4, MARGIN = 1, STATS = wHvar/sdFiv^2,
      FUN = "*"), Zv3)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(Zu10, MARGIN = 1, STATS = wHvar,
    FUN = "*") - sweep(Zu4, MARGIN = 1, STATS = wHvar * apply(sigx2,
    1, sum)/sdFiv^2, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    Zv18) - crossprod(sweep(Zv3, MARGIN = 1, STATS = wHvar/sdFiv^2,
    FUN = "*"), Zv3)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(Zv21,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(Zv3, MARGIN = 1,
    STATS = wHvar * apply(sigx2, 1, sum)/sdFiv^2, FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (apply(sigx26, 1, sum) - apply(sigx2,
    1, sum)^2/sdFiv)/sdFiv)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for weibull-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
weibullnormAlgOpt <- function(start, olsParam, dataTable, S,
  nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar,
  wHvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cweibullnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar)), hessian = 0,
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cweibullnormlike, grad = cgradweibullnormlike,
      hess = chessweibullnormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      hs = function(parm) as(-chessweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      hess = function(parm) -chessweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      gradient = function(parm) -colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
      hessian = function(parm) -chessweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradweibullnormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat, wHvar = wHvar))
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessweibullnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)
    if (method == "sr1")
      mleObj$hessian <- chessweibullnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat, wHvar = wHvar)
  }
  mleObj$logL_OBS <- cweibullnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar)
  mleObj$gradL_OBS <- cgradweibullnormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, N = N, FiMat = FiMat, wHvar = wHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# average efficiency (BC style) evaluation ----------
#' function to estimate unconditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigma standard error of the weibull distribution
#' @param k location parameter
#' @noRd
fnExpUWeiNorm <- function(u, sigma, k) {
  exp(-u) * k/sigma * (u/sigma)^(k - 1) * exp(-(u/sigma)^k)
}

# fn conditional inefficiencies ----------
#' function to estimate conditional inefficiency
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondEffWeibull <- function(u, sigmaU, sigmaV, k, epsilon, S) {
  u * k/(sigmaU * sigmaV) * (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) *
    dnorm((epsilon + S * u)/sigmaV)
}

# fn conditional efficiencies ----------
#' function to estimate conditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondBCEffWeibull <- function(u, sigmaU, sigmaV, k, epsilon,
  S) {
  exp(-u) * k/(sigmaU * sigmaV) * (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) *
    dnorm((epsilon + S * u)/sigmaV)
}

# fn reciprocal conditional efficiencies----------
#' function to estimate conditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param k location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondBCreciprocalEffWeibull <- function(u, sigmaU, sigmaV, k,
  epsilon, S) {
  exp(u) * k/(sigmaU * sigmaV) * (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) *
    dnorm((epsilon + S * u)/sigmaV)
}

# Conditional efficiencies estimation ----------
#' efficiencies for weibull-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
cweibullnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  u <- numeric(object$Nobs)
  density_epsilon_vec <- numeric(object$Nobs)
  for (i in seq_len(object$Nobs)) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon_vec[i] <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u[i] <- hcubature(f = fnCondEffWeibull, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
  }
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- numeric(object$Nobs)
    teBC_reciprocal <- numeric(object$Nobs)
    for (i in seq_len(object$Nobs)) {
      teBC[i] <- hcubature(f = fnCondBCEffWeibull, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), k = k, epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
      teBC_reciprocal[i] <- hcubature(f = fnCondBCreciprocalEffWeibull,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        k = k, epsilon = epsilon[i], S = object$S, vectorInterface = FALSE,
        tol = 1e-15)$integral/density_epsilon_vec[i]
    }
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for weibull-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargweibullnorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmargweibullnorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
