################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: lognormal - normal                                              #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for lognormal-normal distribution
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
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
clognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ll <- numeric(N)
  ur <- list()
  for (i in seq_len(N)) {
    ur[[i]] <- exp(mu[i] + exp(Wu[i]/2) * qnorm(FiMat[i,
      ]))
    ll[i] <- log(mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      S * ur[[i]])/exp(Wv[i]/2))))
  }
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for lognormal-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
cstlognorm <- function(olsObj, epsiRes, S, nmuZUvar, nuZUvar,
  uHvar, muHvar, nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  varu <- tryCatch((nleqslv::nleqslv(x = 0.01, fn = function(x) {
    -exp(9 * x^2/2) + 3 * exp(5 * x^2/2) - 2 * exp(3 * x^2/2) -
      S * m3
  }, method = "Newton")$x)^2, error = function(e) e)
  if (inherits(varu, "error")) {
    varu <- 0.01
  }
  varv <- if ((m2 - exp(varu) * (exp(varu) - 1)) < 0) {
    abs(m2 - exp(varu) * (exp(varu) - 1))
  } else {
    (m2 - exp(varu) * (exp(varu) - 1))
  }
  dep_u <- 1/2 * log((log(1/2 + sqrt(4 * ((epsiRes^2 - varv)^2)^(1/2))/2))^2)
  dep_v <- 1/2 * log((epsiRes^2 - exp(varu) * (exp(varu) -
    1))^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetu$coefficients))) {
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  }
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar,
      drop = FALSE]))
  }
  reg_hetmu <- if (nmuZUvar == 1) {
    lm(epsiRes ~ 1)
  } else {
    lm(epsiRes ~ ., data = as.data.frame(muHvar[, 2:nmuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetmu$coefficients))) {
    stop("at least one of the OLS coefficients of 'muhet' is NA: ",
      paste(colnames(muHvar)[is.na(reg_hetmu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  }
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  omega <- coefficients(reg_hetmu)
  names(omega) <- paste0("Zmu_", colnames(muHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- beta <- c(olsObj[1] + S * exp(varu/2), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, omega, delta, phi))
}

# Gradient of the likelihood function ----------
#' gradient for lognormal-normal distribution
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
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
cgradlognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  qFimat <- qnorm(FiMat)
  WuqFi <- sweep(qFimat, MARGIN = 1, STATS = exp(Wu/2), FUN = "*")
  WumuqFi <- sweep(WuqFi, MARGIN = 1, STATS = mu, FUN = "+")
  WumuqFiepsi <- sweep(S * exp(WumuqFi), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  WuWvmuqFiepsi <- sweep(WumuqFiepsi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/")
  dqFi <- dnorm(WuWvmuqFiepsi)
  WvdqFi <- apply(sweep(dqFi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/"), 1, sum)
  sigx1 <- sweep(dqFi * (WumuqFiepsi), MARGIN = 1, STATS = exp(3 *
    Wv/2), FUN = "/")
  sigx2 <- sweep(dqFi * exp(WumuqFi) * (WumuqFiepsi), MARGIN = 1,
    STATS = exp(3 * Wv/2), FUN = "/")
  sigx3 <- sweep(dqFi * exp(WumuqFi) * qFimat * (WumuqFiepsi),
    MARGIN = 1, STATS = exp(Wu/2)/exp(3 * Wv/2), FUN = "*")
  sigx4 <- sweep(0.5 * (dqFi * (WumuqFiepsi)^2), MARGIN = 1,
    STATS = exp(Wv), FUN = "/")
  sigx5 <- sweep((sigx4 - 0.5 * dqFi), MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    gx[, k] <- apply(sweep(sigx1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/WvdqFi
  }
  gmu <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in seq_len(nmuZUvar)) {
    gmu[, k] <- apply(sweep(sigx2, MARGIN = 1, STATS = -(S *
      muHvar[, k]), FUN = "*"), 1, sum)/WvdqFi
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    gu[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = -(0.5 *
      (S * uHvar[, k])), FUN = "*"), 1, sum)/WvdqFi
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    gv[, k] <- apply(sweep(sigx5, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/WvdqFi
  }
  gradll <- cbind(gx, gmu, gu, gv)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for lognormal-normal distribution
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
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
chesslognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S, N, FiMat, wHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv_h <- exp(Wv/2)
  qFimat <- qnorm(FiMat)
  WuqFi <- sweep(qFimat, MARGIN = 1, STATS = ewu_h, FUN = "*")
  WumuqFi <- sweep(WuqFi, MARGIN = 1, STATS = mu, FUN = "+")
  WumuqFiepsi <- sweep(S * exp(WumuqFi), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  WuWvmuqFiepsi_sq <- (WumuqFiepsi)/ewv_h^2
  WuWvmuqFiepsi_cub <- (WumuqFiepsi)/ewv_h^3
  WuWvmuqFiepsi <- sweep(WumuqFiepsi, MARGIN = 1, STATS = ewv_h,
    FUN = "/")
  dqFi <- dnorm(WuWvmuqFiepsi)
  WvdqFi <- apply(sweep(dqFi, MARGIN = 1, STATS = ewv_h, FUN = "/"),
    1, sum)
  sigx1 <- sweep(dqFi * exp(WumuqFi) * qFimat * WuWvmuqFiepsi_cub,
    MARGIN = 1, STATS = ewu_h, FUN = "*")
  sigx2 <- sweep(dqFi * (WumuqFiepsi)^2, MARGIN = 1, STATS = (1/ewv_h^2),
    FUN = "*")
  sigx3 <- sweep((0.5 * sigx2 - 0.5 * dqFi), MARGIN = 1, STATS = 1/ewv_h,
    FUN = "*")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_len(nXvar)) {
    gx[, k] <- apply(sweep(dqFi * WuWvmuqFiepsi_cub, MARGIN = 1,
      STATS = Xvar[, k], FUN = "*"), 1, sum)
  }
  gmu <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in seq_len(nmuZUvar)) {
    gmu[, k] <- apply(sweep(-(S * dqFi * exp(WumuqFi) * WuWvmuqFiepsi_cub),
      MARGIN = 1, STATS = muHvar[, k], FUN = "*"), 1, sum)
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_len(nuZUvar)) {
    gu[, k] <- apply(sweep(-(0.5 * (S * sigx1)), MARGIN = 1,
      STATS = uHvar[, k], FUN = "*"), 1, sum)
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_len(nvZVvar)) {
    gv[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)
  }
  Q <- dim(FiMat)[2]
  sigx4 <- sweep(((WumuqFiepsi)^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*") - 1
  sigx5 <- sweep(sigx4 * dqFi, MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  sigx6 <- sweep(sigx4 * dqFi * exp(WumuqFi), MARGIN = 1, STATS = 1/ewv_h^3,
    FUN = "*")
  sigx7 <- sweep(sigx4 * dqFi * exp(WumuqFi) * qFimat, MARGIN = 1,
    STATS = ewu_h/ewv_h^3, FUN = "*")
  sigx8 <- sweep(((WumuqFiepsi)^2), MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*") - 2
  sigx9 <- (0.5 * sigx8 - 0.5) * dqFi * WuWvmuqFiepsi_cub
  sigx10 <- sweep(((1 - S * exp(WumuqFi) * WuWvmuqFiepsi_sq) *
    (WumuqFiepsi) + S * exp(WumuqFi)) * dqFi * exp(WumuqFi),
    MARGIN = 1, STATS = 1/ewv_h^3, FUN = "*")
  sigx11 <- sweep(((1 - S * exp(WumuqFi) * WuWvmuqFiepsi_sq) *
    (WumuqFiepsi) + S * exp(WumuqFi)) * dqFi * exp(WumuqFi) *
    qFimat, MARGIN = 1, STATS = ewu_h/ewv_h^3, FUN = "*")
  sigx12 <- sweep((WumuqFiepsi)^2, MARGIN = 1, STATS = 1/ewv_h^2,
    FUN = "*")
  sigx13 <- (0.5 + 0.5 * (2 - sigx12)) * dqFi * exp(WumuqFi) *
    WuWvmuqFiepsi_cub
  sigx14 <- sweep(exp(WumuqFi) * qFimat, MARGIN = 1, STATS = 0.5 *
    (S * ewu_h), FUN = "*")
  sigx15 <- sweep((((0.5 - 0.5 * (S * exp(WumuqFi) * WuWvmuqFiepsi_sq)) *
    WuqFi + 0.5) * (WumuqFiepsi) + sigx14) * dqFi * exp(WumuqFi) *
    qFimat, MARGIN = 1, STATS = ewu_h/ewv_h^3, FUN = "*")
  sigx16 <- (0.25 + 0.5 * (1 - 0.5 * (sigx12))) * sigx1
  sigx17 <- sweep(((0.5 * (0.5 * (sigx12) - 1) - 0.25) * dqFi *
    sigx12 - 0.5 * (0.5 * sigx2 - 0.5 * dqFi)), MARGIN = 1,
    STATS = 1/ewv_h, FUN = "*")
  Xsq <- list()
  Xmu <- list()
  Xu <- list()
  Xv <- list()
  Musq <- list()
  Muu <- list()
  Muv <- list()
  Usq <- list()
  Uv <- list()
  Vsq <- list()
  for (r in seq_len(Q)) {
    Xsq[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sigx5[, r]/WvdqFi, FUN = "*"), Xvar)
    Xmu[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (S * sigx6)[, r]/WvdqFi, FUN = "*"), muHvar)
    Xu[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * sigx7))[, r]/WvdqFi, FUN = "*"), uHvar)
    Xv[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sigx9[, r]/WvdqFi, FUN = "*"), vHvar)
    Musq[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (S * sigx10)[, r]/WvdqFi, FUN = "*"), muHvar)
    Muu[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * sigx11))[, r]/WvdqFi, FUN = "*"), uHvar)
    Muv[[r]] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = S *
      wHvar * sigx13[, r]/WvdqFi, FUN = "*"), vHvar)
    Usq[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
      (0.5 * (S * sigx15))[, r]/WvdqFi, FUN = "*"), uHvar)
    Uv[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = S *
      wHvar * sigx16[, r]/WvdqFi, FUN = "*"), vHvar)
    Vsq[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      sigx17[, r]/WvdqFi, FUN = "*"), vHvar)
  }
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar,
    ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", Xsq) - crossprod(sweep(gx,
    MARGIN = 1, STATS = wHvar/WvdqFi^2, FUN = "*"), gx)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- Reduce("+",
    Xmu) - crossprod(sweep(gx, MARGIN = 1, STATS = wHvar/WvdqFi^2,
    FUN = "*"), gmu)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- Reduce("+", Xu) - crossprod(sweep(gx, MARGIN = 1,
    STATS = wHvar/WvdqFi^2, FUN = "*"), gu)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", Xv) - crossprod(sweep(gx,
    MARGIN = 1, STATS = wHvar/WvdqFi^2, FUN = "*"), gv)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- Reduce("+", Musq) - crossprod(sweep(gmu,
    MARGIN = 1, STATS = wHvar/WvdqFi^2, FUN = "*"), gmu)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+", Muu) -
    crossprod(sweep(gmu, MARGIN = 1, STATS = wHvar/WvdqFi^2,
      FUN = "*"), gu)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+",
    Muv) - crossprod(sweep(gmu, MARGIN = 1, STATS = wHvar/WvdqFi^2,
    FUN = "*"), gv)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- Reduce("+",
    Usq) - crossprod(sweep(gu, MARGIN = 1, STATS = wHvar/WvdqFi^2,
    FUN = "*"), gu)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
      nuZUvar + nvZVvar)] <- Reduce("+", Uv) - crossprod(sweep(gu,
    MARGIN = 1, STATS = wHvar/WvdqFi^2, FUN = "*"), gv)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- Reduce("+", Vsq) -
    crossprod(sweep(gv, MARGIN = 1, STATS = wHvar/WvdqFi^2,
      FUN = "*"), gv)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for lognormal-normal distribution
#' @param start starting value for optimization
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables ()
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @param Yvar vector of dependent variable
#' @param Xvar matrix of main variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param N number of observations
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
lognormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  muHvar, nmuZUvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  wHvar, Yvar, Xvar, method, printInfo, itermax, stepmax, tol,
  gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstlognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar)
  startLoglik <- sum(clognormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    gr = function(parm) -colSums(cgradlognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = clognormlike,
    grad = cgradlognormlike, hess = chesslognormlike, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    gr = function(parm) -colSums(cgradlognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    gr = function(parm) -colSums(cgradlognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    hs = function(parm) as(-chesslognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar),
      "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    gr = function(parm) -colSums(cgradlognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)),
    hess = function(parm) -chesslognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(clognormlike(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar)), gradient = function(parm) -colSums(cgradlognormlike(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar)), hessian = function(parm) -chesslognormlike(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat, wHvar = wHvar), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradlognormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat, wHvar = wHvar))
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
      mleObj$hessian <- chesslognormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)
    if (method == "sr1")
      mleObj$hessian <- chesslognormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat, wHvar = wHvar)
  }
  mleObj$logL_OBS <- clognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)
  mleObj$gradL_OBS <- cgradlognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat, wHvar = wHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# average efficiency (BC style) evaluation ----------
#' function to estimate unconditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigma standard error of the weibull distribution
#' @param mu location parameter
#' @noRd
fnExpULogNorm <- function(u, sigma, mu) {
  1/(u * sigma * sqrt(2 * pi)) * exp(-(log(u) - mu)^2/(2 *
    sigma^2) - u)
}

# fn conditional inefficiencies ----------
#' function to estimate conditional inefficiency
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param mu location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondEffLogNorm <- function(u, sigmaU, sigmaV, mu, epsilon,
  S) {
  1/(sigmaU * sigmaV) * dnorm((log(u) - mu)/sigmaU) * dnorm((epsilon +
    S * u)/sigmaV)
}

# fn conditional efficiencies ----------
#' function to estimate conditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param mu location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondBCEffLogNorm <- function(u, sigmaU, sigmaV, mu, epsilon,
  S) {
  exp(-u)/u * 1/(sigmaU * sigmaV) * dnorm((log(u) - mu)/sigmaU) *
    dnorm((epsilon + S * u)/sigmaV)
}

# fn reciproccal conditional efficiencies ----------
#' function to estimate conditional efficiency (Battese and Coelli style)
#' @param u inefficiency variable over which integration will be done
#' @param sigmaU standard error of the weibull distribution
#' @param sigmaV standard error of the two-sided error component
#' @param mu location parameter
#' @param epsilon composite noise
#' @param S integer for cost/prod estimation
#' @noRd
fnCondBCreciprocalEffLogNorm <- function(u, sigmaU, sigmaV, mu,
  epsilon, S) {
  exp(u)/u * 1/(sigmaU * sigmaV) * dnorm((log(u) - mu)/sigmaU) *
    dnorm((epsilon + S * u)/sigmaV)
}

# Conditional efficiencies estimation ----------
#' efficiencies for lognormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
clognormeff <- function(object, level) {
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
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  u <- numeric(object$Nobs)
  density_epsilon_vec <- numeric(object$Nobs)
  for (i in seq_len(object$Nobs)) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon_vec[i] <- (mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2))))
    u[i] <- hcubature(f = fnCondEffLogNorm, lowerLimit = 0,
      upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
      sigmaV = exp(Wv[i]/2), mu = mu[i], epsilon = epsilon[i],
      S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
  }
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- numeric(object$Nobs)
    teBC_reciprocal <- numeric(object$Nobs)
    for (i in seq_len(object$Nobs)) {
      teBC[i] <- hcubature(f = fnCondBCEffLogNorm, lowerLimit = 0,
        upperLimit = Inf, maxEval = 100, fDim = 1, sigmaU = exp(Wu[i]/2),
        sigmaV = exp(Wv[i]/2), mu = mu[i], epsilon = epsilon[i],
        S = object$S, vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
      teBC_reciprocal[i] <- hcubature(f = fnCondBCreciprocalEffLogNorm,
        lowerLimit = 0, upperLimit = Inf, maxEval = 100,
        fDim = 1, sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2),
        mu = mu[i], epsilon = epsilon[i], S = object$S,
        vectorInterface = FALSE, tol = 1e-15)$integral/density_epsilon_vec[i]
    }
    res <- data.frame(u = u, teJLMS = teJLMS, teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for lognormal-normal distribution
#' @param object object of class sfacross
#' @noRd
cmarglognorm_Eu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(data.frame(margEff))
}

cmarglognorm_Vu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(data.frame(margEff))
}
