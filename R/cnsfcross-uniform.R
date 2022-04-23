################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Contaminated Noise Stochastic Frontier Model                          #
# Two types: - Common inefficiency component (sigma_u)                         #
#            - Mixture composed error (mcesf)                                  # 
# Convolution: uniform - normal                                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf uniform-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u
ccnsfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    S * epsilon)/exp(Wv2/2)) - pnorm(S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u
cmcesfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) +
    S * epsilon)/exp(Wv1/2)) - pnorm(S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) +
    S * epsilon)/exp(Wv2/2)) - pnorm(S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf uniform-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param tol parameter tolerance
#' @noRd
# same sigma_u
cstcnsfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax,
  printInfo, tol) {
  cat("Initialization: SFA + uniform - normal distributions...\n")
  initUni <- maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgraduninormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, nvZVvar = 1,
    uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[,
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initUni$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("CN_", colnames(Zvar)))
  names(initUni$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initUni = initUni))
}

# different sigma_u
cstmcesfuninorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + uniform - normal distributions...\n")
  initUni <- maxLik(logLik = cuninormlike, start = cstuninorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])),
    grad = cgraduninormlike, method = "BFGS", control = list(iterlim = itermax,
      printLevel = printInfo, reltol = tol), nXvar = nXvar,
    nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]),
    vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar)
  Esti <- initUni$estimate
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), paste0("MCE_",
    colnames(Zvar)))
  names(initUni$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initUni = initUni))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf uniform-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u
cgradcnsfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (sqrt(12) * ewu_h + S * epsilon)/ewv1_h
  musig2 <- (sqrt(12) * ewu_h + S * epsilon)/ewv2_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi1 <- S * (epsilon)/ewv1_h
  epsi2 <- S * (epsilon)/ewv2_h
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  pepsi1 <- pnorm(epsi1)
  pepsi2 <- pnorm(epsi2)
  sigx1 <- ((depsi1 - dmusig1) * ewz/(sqrt(12) * (wzdeno *
    ewv1_h)) + 1/sqrt(12) * (prC * (depsi2 - dmusig2)/ewv2_h))
  sigx2 <- (prC * (pmusig2 - pepsi2)/sqrt(12) + ewz * (pmusig1 -
    pepsi1)/(sqrt(12) * wzdeno))
  sigx3 <- (0.5 * (prC * dmusig2/ewv2_h) + dmusig1 * ewz/(2 *
    (wzdeno * ewv1_h)))
  sigx4_1 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * ((sqrt(12) *
    ewu_h + S * epsilon) * dmusig1))
  sigx4_2 <- (0.5 * (S * depsi2 * epsilon) - 0.5 * ((sqrt(12) *
    ewu_h + S * epsilon) * dmusig2))
  sigx5 <- ((1/(sqrt(12) * wzdeno) - sqrt(12) * (ewz/(sqrt(12) *
    wzdeno)^2)) * (pmusig1 - pepsi1) - 1/sqrt(12) * (prC *
    (pmusig2 - pepsi2)/wzdeno))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx3 *
    ewu_h/sigx2 - 0.5), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = sigx4_1 * ewz/(sqrt(12) * (sigx2 * wzdeno * ewv1_h)),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 1/sqrt(12) *
    (sigx4_2 * prC/(sigx2 * ewv2_h)), FUN = "*"), sweep(Zvar,
    MARGIN = 1, STATS = sigx5 * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u
cgradmcesfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  musig1 <- (sqrt(12) * ewu1_h + S * epsilon)/ewv1_h
  musig2 <- (sqrt(12) * ewu2_h + S * epsilon)/ewv2_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi1 <- S * (epsilon)/ewv1_h
  epsi2 <- S * (epsilon)/ewv2_h
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  pepsi1 <- pnorm(epsi1)
  pepsi2 <- pnorm(epsi2)
  sigx1 <- (prC * (depsi2 - dmusig2)/(sqrt(12) * (ewu2_h *
    ewv2_h)) + (depsi1 - dmusig1) * ewz/(sqrt(12) * (wzdeno *
    ewu1_h * ewv1_h)))
  sigx2 <- (prC * (pmusig2 - pepsi2)/(sqrt(12) * ewu2_h) +
    ewz * (pmusig1 - pepsi1)/(sqrt(12) * (wzdeno * ewu1_h)))
  sigx3 <- (dmusig1/(2 * (wzdeno * ewv1_h)) - sqrt(12)/2 *
    (wzdeno * ewu1_h * (pmusig1 - pepsi1)/(sqrt(12) * (wzdeno *
      ewu1_h))^2))
  sigx4 <- (dmusig2/(2 * ewv2_h) - sqrt(12)/2 * (ewu2_h * (pmusig2 -
    pepsi2)/(sqrt(12) * ewu2_h)^2))
  sigx5 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * ((sqrt(12) *
    ewu1_h + S * epsilon) * dmusig1))
  sigx6 <- (sqrt(12) * (sigx2 * wzdeno * ewu1_h * ewv1_h))
  sigx7 <- (0.5 * (S * depsi2 * epsilon) - 0.5 * ((sqrt(12) *
    ewu2_h + S * epsilon) * dmusig2))
  sigx8 <- (sqrt(12) * (sigx2 * ewu2_h * ewv2_h))
  sigx9 <- (1/(sqrt(12) * (wzdeno * ewu1_h)) - sqrt(12) * (ewu1_h *
    ewz/(sqrt(12) * (wzdeno * ewu1_h))^2))
  sigx10 <- (sigx9 * (pmusig1 - pepsi1) - prC * (pmusig2 -
    pepsi2)/(sqrt(12) * (wzdeno * ewu2_h)))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx1/sigx2,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx3 *
    ewz/sigx2, FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = prC *
    sigx4/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx5 *
    ewz/sigx6, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx7 *
    prC/sigx8, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx10 *
    ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf uniform-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u
chesscnsfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (sqrt(12) * ewu_h + S * epsilon)/ewv1_h
  musig2 <- (sqrt(12) * ewu_h + S * epsilon)/ewv2_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi1 <- S * (epsilon)/ewv1_h
  epsi2 <- S * (epsilon)/ewv2_h
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  pepsi1 <- pnorm(epsi1)
  pepsi2 <- pnorm(epsi2)
  sigx1 <- ((depsi1 - dmusig1) * ewz/(sqrt(12) * (wzdeno *
    ewv1_h)) + 1/sqrt(12) * (prC * (depsi2 - dmusig2)/ewv2_h))
  sigx2 <- (prC * (pmusig2 - pepsi2)/sqrt(12) + ewz * (pmusig1 -
    pepsi1)/(sqrt(12) * wzdeno))
  sigx3 <- (0.5 * (prC * dmusig2/ewv2_h) + dmusig1 * ewz/(2 *
    (wzdeno * ewv1_h)))
  sigx4_1 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * ((sqrt(12) *
    ewu_h + S * epsilon) * dmusig1))
  sigx4_2 <- (0.5 * (S * depsi2 * epsilon) - 0.5 * ((sqrt(12) *
    ewu_h + S * epsilon) * dmusig2))
  sigx5 <- ((1/(sqrt(12) * wzdeno) - sqrt(12) * (ewz/(sqrt(12) *
    wzdeno)^2)) * (pmusig1 - pepsi1) - 1/sqrt(12) * (prC *
    (pmusig2 - pepsi2)/wzdeno))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (1/sqrt(12) * (prC * (S * depsi2 * epsilon -
      (sqrt(12) * ewu_h + S * epsilon) * dmusig2)/ewv2_h^3) +
      ewz * (S * depsi1 * epsilon - (sqrt(12) * ewu_h +
        S * epsilon) * dmusig1)/(sqrt(12) * (wzdeno *
        ewv1_h^3)) - sigx1^2/sigx2)/sigx2, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (prC * dmusig2/ewv2_h^3) +
      dmusig1 * ewz/(2 * (wzdeno * ewv1_h^3))) * (sqrt(12) *
      ewu_h + S * epsilon) - sigx1 * sigx3/sigx2) * ewu_h/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * epsilon^2/ewv1_h^2 - 1)) -
    0.5 * (((sqrt(12) * ewu_h + S * epsilon)^2/ewv1_h^2 -
      1) * dmusig1))/(sqrt(12) * (sigx2 * wzdeno * ewv1_h)) -
    sqrt(12) * (sigx1 * sigx4_1 * wzdeno * ewv1_h/(sqrt(12) *
      (sigx2 * wzdeno * ewv1_h))^2)) * ewz, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * 1/sqrt(12) * (S * ((0.5 * (depsi2 * (S^2 *
      epsilon^2/ewv2_h^2 - 1)) - 0.5 * (((sqrt(12) * ewu_h +
      S * epsilon)^2/ewv2_h^2 - 1) * dmusig2))/(sigx2 *
      ewv2_h) - sigx1 * sigx4_2 * ewv2_h/(sigx2 * ewv2_h)^2) *
      prC), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1/(sqrt(12) * wzdeno) -
      sqrt(12) * (ewz/(sqrt(12) * wzdeno)^2)) * (depsi1 -
      dmusig1)/ewv1_h - (sigx5 * sigx1/sigx2 + 1/sqrt(12) *
      (prC * (depsi2 - dmusig2)/(wzdeno * ewv2_h)))) *
      ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 - sigx3 * ewu_h/sigx2) * sigx3 - (sqrt(12)/4 *
      (prC * dmusig2/ewv2_h^3) + sqrt(12)/2 * (dmusig1 *
      ewz)/(2 * (wzdeno * ewv1_h^3))) * (sqrt(12) * ewu_h +
      S * epsilon) * ewu_h) * ewu_h/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewu_h + S * epsilon)^2/ewv1_h^2)) *
      dmusig1)/(sqrt(12) * (sigx2 * wzdeno * ewv1_h)) +
      sqrt(12) * (sigx3 * sigx4_1 * wzdeno * ewv1_h/(sqrt(12) *
        (sigx2 * wzdeno * ewv1_h))^2)) * ewu_h * ewz),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (1/sqrt(12) * ((sigx3 *
      sigx4_2 * ewv2_h/(sigx2 * ewv2_h)^2 + 0.5 * ((sqrt(12)/2 -
      sqrt(12)/2 * ((sqrt(12) * ewu_h + S * epsilon)^2/ewv2_h^2)) *
      dmusig2/(sigx2 * ewv2_h))) * prC * ewu_h)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sqrt(12)/2 * ((1/(sqrt(12) *
      wzdeno) - sqrt(12) * (ewz/(sqrt(12) * wzdeno)^2)) *
      dmusig1/ewv1_h) - (sigx5 * sigx3/sigx2 + 0.5 * (prC *
      dmusig2/(wzdeno * ewv2_h)))) * ewu_h * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      epsilon^3) - 0.25 * ((sqrt(12) * ewu_h + S * epsilon)^3 *
      dmusig1))/(sqrt(12) * (sigx2 * wzdeno * ewv1_h^3)) -
      sqrt(12) * ((sigx4_1 * ewz/sqrt(12) + 0.5 * (sigx2 *
        wzdeno * ewv1_h)) * sigx4_1/(sqrt(12) * (sigx2 *
        wzdeno * ewv1_h))^2)) * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (1/sqrt(12) * (sigx4_1 * sigx4_2 * prC * ewv2_h * ewz)/(sqrt(12) *
      ((sigx2 * ewv2_h)^2 * wzdeno * ewv1_h))), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx4_1 * (1/(sqrt(12) *
      wzdeno) - (sigx5/(sqrt(12) * (sigx2 * wzdeno)) +
      sqrt(12)/(sqrt(12) * wzdeno)^2) * ewz) * ewz/(sigx2 *
      ewv1_h), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * 1/sqrt(12) * (((0.25 * (S^3 * depsi2 *
      epsilon^3) - 0.25 * ((sqrt(12) * ewu_h + S * epsilon)^3 *
      dmusig2))/(sigx2 * ewv2_h^3) - (1/sqrt(12) * (sigx4_2 *
      prC) + 0.5 * (sigx2 * ewv2_h)) * sigx4_2/(sigx2 *
      ewv2_h)^2) * prC), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1/sqrt(12) * (sigx5/sigx2) +
      1/sqrt(12)/wzdeno) * sigx4_2 * prC * ewz/(sigx2 *
      ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((0.577350269189626 * (prC *
      (pmusig2 - pepsi2)/wzdeno^2) - (sigx5^2/sigx2 + (sqrt(12) *
      (1 - 24 * (wzdeno * ewz/(sqrt(12) * wzdeno)^2)) +
      sqrt(12)) * (pmusig1 - pepsi1)/(sqrt(12) * wzdeno)^2)) *
      ewz + (1/(sqrt(12) * wzdeno) - sqrt(12) * (ewz/(sqrt(12) *
      wzdeno)^2)) * (pmusig1 - pepsi1) - 1/sqrt(12) * (prC *
      (pmusig2 - pepsi2)/wzdeno)) * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# different sigma_u
chessmcesfuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  musig1 <- (sqrt(12) * ewu1_h + S * epsilon)/ewv1_h
  musig2 <- (sqrt(12) * ewu2_h + S * epsilon)/ewv2_h
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  epsi1 <- S * (epsilon)/ewv1_h
  epsi2 <- S * (epsilon)/ewv2_h
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  pepsi1 <- pnorm(epsi1)
  pepsi2 <- pnorm(epsi2)
  sigx1 <- (prC * (depsi2 - dmusig2)/(sqrt(12) * (ewu2_h *
    ewv2_h)) + (depsi1 - dmusig1) * ewz/(sqrt(12) * (wzdeno *
    ewu1_h * ewv1_h)))
  sigx2 <- (prC * (pmusig2 - pepsi2)/(sqrt(12) * ewu2_h) +
    ewz * (pmusig1 - pepsi1)/(sqrt(12) * (wzdeno * ewu1_h)))
  sigx3 <- (dmusig1/(2 * (wzdeno * ewv1_h)) - sqrt(12)/2 *
    (wzdeno * ewu1_h * (pmusig1 - pepsi1)/(sqrt(12) * (wzdeno *
      ewu1_h))^2))
  sigx4 <- (dmusig2/(2 * ewv2_h) - sqrt(12)/2 * (ewu2_h * (pmusig2 -
    pepsi2)/(sqrt(12) * ewu2_h)^2))
  sigx5 <- (0.5 * (S * depsi1 * epsilon) - 0.5 * ((sqrt(12) *
    ewu1_h + S * epsilon) * dmusig1))
  sigx6 <- (sqrt(12) * (sigx2 * wzdeno * ewu1_h * ewv1_h))
  sigx7 <- (0.5 * (S * depsi2 * epsilon) - 0.5 * ((sqrt(12) *
    ewu2_h + S * epsilon) * dmusig2))
  sigx8 <- (sqrt(12) * (sigx2 * ewu2_h * ewv2_h))
  sigx9 <- (1/(sqrt(12) * (wzdeno * ewu1_h)) - sqrt(12) * (ewu1_h *
    ewz/(sqrt(12) * (wzdeno * ewu1_h))^2))
  sigx10 <- (sigx9 * (pmusig1 - pepsi1) - prC * (pmusig2 -
    pepsi2)/(sqrt(12) * (wzdeno * ewu2_h)))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (prC * (S * depsi2 * epsilon - (sqrt(12) *
      ewu2_h + S * epsilon) * dmusig2)/(sqrt(12) * (ewu2_h *
      ewv2_h^3)) + ewz * (S * depsi1 * epsilon - (sqrt(12) *
      ewu1_h + S * epsilon) * dmusig1)/(sqrt(12) * (wzdeno *
      ewu1_h * ewv1_h^3)) - sigx1^2/sigx2)/sigx2, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((sqrt(12) * ewu1_h +
      S * epsilon) * dmusig1/(2 * (wzdeno * ewv1_h^2)) -
      sqrt(12)/2 * (wzdeno * (depsi1 - dmusig1) * ewu1_h/(sqrt(12) *
        (wzdeno * ewu1_h))^2))/ewv1_h - sigx1 * sigx3/sigx2) *
      ewz/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((sqrt(12) * ewu2_h +
      S * epsilon) * dmusig2/(2 * ewv2_h^2) - sqrt(12)/2 *
      ((depsi2 - dmusig2) * ewu2_h/(sqrt(12) * ewu2_h)^2))/ewv2_h -
      sigx1 * sigx4/sigx2) * prC/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((0.5 * (depsi1 * (S^2 * epsilon^2/ewv1_h^2 - 1)) -
    0.5 * (((sqrt(12) * ewu1_h + S * epsilon)^2/ewv1_h^2 -
      1) * dmusig1))/sigx6 - sqrt(12) * (sigx1 * sigx5 *
    wzdeno * ewu1_h * ewv1_h/sigx6^2)) * ewz, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (depsi2 * (S^2 *
      epsilon^2/ewv2_h^2 - 1)) - 0.5 * (((sqrt(12) * ewu2_h +
      S * epsilon)^2/ewv2_h^2 - 1) * dmusig2))/sigx8 -
      sqrt(12) * (sigx1 * sigx7 * ewu2_h * ewv2_h/sigx8^2)) *
      prC, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx9 * (depsi1 - dmusig1)/ewv1_h -
      (sigx1 * sigx10/sigx2 + prC * (depsi2 - dmusig2)/(sqrt(12) *
        (wzdeno * ewu2_h * ewv2_h)))) * ewz/sigx2, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    (((sqrt(12)/2 * (((sqrt(12)/2 * (dmusig1/ewv1_h) - 12 *
      (wzdeno^2 * ewu1_h * (pmusig1 - pepsi1)/(sqrt(12) *
        (wzdeno * ewu1_h))^2)) * ewu1_h + 0.5 * (pmusig1 -
      pepsi1)) * wzdeno/(sqrt(12) * (wzdeno * ewu1_h))^2) +
      sqrt(12)/2 * ((sqrt(12) * ewu1_h + S * epsilon) *
        dmusig1)/(2 * (wzdeno * ewv1_h^3))) * ewu1_h +
      sigx3^2 * ewz/sigx2) * ewz/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prC * sigx3 * sigx4 * ewz/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 *
      ((sqrt(12) * ewu1_h + S * epsilon)^2/ewv1_h^2)) *
      dmusig1)/(sqrt(12) * (sigx2 * wzdeno * ewv1_h)) +
      sqrt(12) * ((sigx3 * ewz + 0.5 * sigx2) * sigx5 *
        wzdeno * ewu1_h * ewv1_h/sigx6^2)) * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sqrt(12) * (sigx7 * prC *
      sigx3 * ewu2_h * ewv2_h * ewz/sigx8^2)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((sqrt(12)/2 * (sigx9 * dmusig1/ewv1_h) - (sqrt(12)/2 *
      wzdeno + sqrt(12) * ((0.5 - 12 * (wzdeno^2 * ewu1_h^2/(sqrt(12) *
      (wzdeno * ewu1_h))^2)) * ewz)) * (pmusig1 - pepsi1)/(sqrt(12) *
      (wzdeno * ewu1_h))^2) * ewu1_h - sigx10 * sigx3 *
      ewz/sigx2) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((prC * sigx4^2/sigx2 +
      (sqrt(12)/2 * (((sqrt(12)/2 * (dmusig2/ewv2_h) -
        12 * (ewu2_h * (pmusig2 - pepsi2)/(sqrt(12) *
          ewu2_h)^2)) * ewu2_h + 0.5 * (pmusig2 - pepsi2))/(sqrt(12) *
        ewu2_h)^2) + sqrt(12)/2 * ((sqrt(12) * ewu2_h +
        S * epsilon) * dmusig2)/(2 * ewv2_h^3)) * ewu2_h) *
      prC/sigx2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sqrt(12) * (sigx5 * prC *
      wzdeno * sigx4 * ewu1_h * ewv1_h * ewz/sigx6^2)),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar *
    ((0.5 * ((sqrt(12)/2 - sqrt(12)/2 * ((sqrt(12) * ewu2_h +
      S * epsilon)^2/ewv2_h^2)) * dmusig2)/(sqrt(12) *
      (sigx2 * ewv2_h)) + sqrt(12) * ((prC * sigx4 + 0.5 *
      sigx2) * sigx7 * ewu2_h * ewv2_h/sigx8^2)) * prC),
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx10 * sigx4/sigx2 + dmusig2/(2 *
      (wzdeno * ewv2_h)) - sqrt(12)/2 * (wzdeno * ewu2_h *
      (pmusig2 - pepsi2)/(sqrt(12) * (wzdeno * ewu2_h))^2)) *
      prC * ewz/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi1 *
      epsilon^3) - 0.25 * ((sqrt(12) * ewu1_h + S * epsilon)^3 *
      dmusig1))/(sqrt(12) * (sigx2 * wzdeno * ewu1_h *
      ewv1_h^3)) - sqrt(12) * ((sigx5 * ewz/sqrt(12) +
      0.5 * (sigx2 * wzdeno * ewu1_h * ewv1_h)) * sigx5/sigx6^2)) *
      ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx5 * sigx7 * prC * ewu2_h * ewv2_h *
      ewz/(wzdeno * sigx8^2 * ewu1_h * ewv1_h)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx5 * (1/(sqrt(12) * (wzdeno *
      ewu1_h)) - (sigx10/(sqrt(12) * (sigx2 * wzdeno *
      ewu1_h)) + sqrt(12) * (ewu1_h/(sqrt(12) * (wzdeno *
      ewu1_h))^2)) * ewz) * ewz/(sigx2 * ewv1_h), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.25 * (S^3 * depsi2 *
      epsilon^3) - 0.25 * ((sqrt(12) * ewu2_h + S * epsilon)^3 *
      dmusig2))/(sqrt(12) * (sigx2 * ewu2_h * ewv2_h^3)) -
      sqrt(12) * ((sigx7 * prC/sqrt(12) + 0.5 * (sigx2 *
        ewu2_h * ewv2_h)) * sigx7/sigx8^2)) * prC, FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx10/(sqrt(12) * sigx2) +
      1/(sqrt(12) * wzdeno)) * sigx7 * prC * ewz/(sigx2 *
      ewu2_h * ewv2_h)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((prC * (1/(sqrt(12) * (wzdeno^2 * ewu2_h)) + sqrt(12) *
      (ewu2_h/(sqrt(12) * (wzdeno * ewu2_h))^2)) * (pmusig2 -
      pepsi2) - (sigx10^2/sigx2 + (sqrt(12) * (1 - 24 *
      (wzdeno * ewu1_h^2 * ewz/(sqrt(12) * (wzdeno * ewu1_h))^2)) +
      sqrt(12)) * ewu1_h * (pmusig1 - pepsi1)/(sqrt(12) *
      (wzdeno * ewu1_h))^2)) * ewz + sigx9 * (pmusig1 -
      pepsi1) - prC * (pmusig2 - pepsi2)/(sqrt(12) * (wzdeno *
      ewu2_h))) * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf uniform-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param S integer for cost/prod estimation
#' @param wHvar vector of weights (weighted likelihood)
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# same sigma_u
cnsfuninormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfuninormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfuninormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfuninormlike, grad = cgradcnsfuninormlike,
      hess = chesscnsfuninormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfuninormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ccnsfuninormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfuninormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chesscnsfuninormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsfuninormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsfuninormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfuninormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

# different sigma_u
mcesfuninormAlgOpt <- function(start, olsParam, dataTable, S,
  wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initUni <- start_st$initUni
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfuninormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...), cg = function(...) maxCG(...),
      sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfuninormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfuninormlike, grad = cgradmcesfuninormlike,
      hess = chessmcesfuninormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmcesfuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesfuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfuninormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar))
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
    if (method %in% c("ucminf", "nlminb"))
      mleObj$hessian <- chessmcesfuninormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesfuninormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesfuninormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfuninormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initUni = initUni))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf uniform-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u
ccnsfuninormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- -exp(Wv2/2) * ((dnorm((sqrt(12) * exp(Wu/2) + object$S *
    epsilon)/exp(Wv2/2)) - dnorm(object$S * epsilon/exp(Wv2/2)))/(pnorm((sqrt(12) *
    exp(Wu/2) + object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S *
    epsilon/exp(Wv2/2)))) - object$S * epsilon
  u2_c2 <- exp(Wv2/2) * (dnorm(object$S * epsilon/exp(Wv2/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv2/2))) - object$S * epsilon/exp(Wv2/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((sqrt(12) *
      exp(Wu/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv2)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu/2))/exp(Wv2/2) + exp(Wv2/2)) -
      pnorm(object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)))/(pnorm((sqrt(12) *
      exp(Wu/2) + object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S *
      epsilon/exp(Wv2/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv2)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv2/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv2)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu/2))/exp(Wv2/2) -
        exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2) -
        exp(Wv2/2)))/(pnorm((sqrt(12) * exp(Wu/2) + object$S *
      epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv2)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv2/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

# different sigma_u
cmcesfuninormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u1_c1 <- -exp(Wv1/2) * ((dnorm((sqrt(12) * exp(Wu1/2) + object$S *
    epsilon)/exp(Wv1/2)) - dnorm(object$S * epsilon/exp(Wv1/2)))/(pnorm((sqrt(12) *
    exp(Wu1/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
    epsilon/exp(Wv1/2)))) - object$S * epsilon
  u2_c1 <- exp(Wv1/2) * (dnorm(object$S * epsilon/exp(Wv1/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv1/2))) - object$S * epsilon/exp(Wv1/2))
  u1_c2 <- -exp(Wv2/2) * ((dnorm((sqrt(12) * exp(Wu2/2) + object$S *
    epsilon)/exp(Wv2/2)) - dnorm(object$S * epsilon/exp(Wv2/2)))/(pnorm((sqrt(12) *
    exp(Wu2/2) + object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S *
    epsilon/exp(Wv2/2)))) - object$S * epsilon
  u2_c2 <- exp(Wv2/2) * (dnorm(object$S * epsilon/exp(Wv2/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv2/2))) - object$S * epsilon/exp(Wv2/2))
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  u1_c <- ifelse(Group_c == 1, u1_c1, u1_c2)
  u2_c <- ifelse(Group_c == 1, u2_c1, u2_c2)
  ineff1_c1 <- ifelse(Group_c == 1, u1_c1, NA)
  ineff1_c2 <- ifelse(Group_c == 2, u1_c2, NA)
  ineff2_c1 <- ifelse(Group_c == 1, u2_c1, NA)
  ineff2_c2 <- ifelse(Group_c == 2, u2_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS1_c1 <- exp(-u1_c1)
    teJLMS1_c2 <- exp(-u1_c2)
    teJLMS1_c <- ifelse(Group_c == 1, teJLMS1_c1, teJLMS1_c2)
    teJLMS2_c1 <- exp(-u2_c1)
    teJLMS2_c2 <- exp(-u2_c2)
    teJLMS2_c <- ifelse(Group_c == 1, teJLMS2_c1, teJLMS2_c2)
    teBC1_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv1/2) + exp(Wv1/2)) -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(pnorm((sqrt(12) *
      exp(Wu1/2) + object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_c1 <- exp(object$S * epsilon + exp(Wv1)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_c2 <- exp(object$S * epsilon + exp(Wv2)/2) * (pnorm((object$S *
      epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv2/2) + exp(Wv2/2)) -
      pnorm(object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)))/(pnorm((sqrt(12) *
      exp(Wu2/2) + object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S *
      epsilon/exp(Wv2/2)))
    teBC2_c2 <- exp(object$S * epsilon + exp(Wv2)/2) * (1 -
      pnorm(object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv2/2)))
    teBC1_c <- ifelse(Group_c == 1, teBC1_c1, teBC1_c2)
    teBC2_c <- ifelse(Group_c == 1, teBC2_c1, teBC2_c2)
    effBC1_c1 <- ifelse(Group_c == 1, teBC1_c1, NA)
    effBC1_c2 <- ifelse(Group_c == 2, teBC1_c2, NA)
    effBC2_c1 <- ifelse(Group_c == 1, teBC2_c1, NA)
    effBC2_c2 <- ifelse(Group_c == 2, teBC2_c2, NA)
    teBC1_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu1/2))/exp(Wv1/2) -
        exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2) -
        exp(Wv1/2)))/(pnorm((sqrt(12) * exp(Wu1/2) +
      object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S *
      epsilon/exp(Wv1/2)))
    teBC2_reciprocal_c1 <- exp(-object$S * epsilon + exp(Wv1)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv1/2)))
    teBC1_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv2)/2) *
      (pnorm((object$S * epsilon + sqrt(12) * exp(Wu2/2))/exp(Wv2/2) -
        exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2) -
        exp(Wv2/2)))/(pnorm((sqrt(12) * exp(Wu2/2) +
      object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S *
      epsilon/exp(Wv2/2)))
    teBC2_reciprocal_c2 <- exp(-object$S * epsilon + exp(Wv2)/2) *
      (1 - pnorm(object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)))/(1 -
      pnorm(object$S * epsilon/exp(Wv2/2)))
    teBC1_reciprocal_c <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      teBC1_reciprocal_c2)
    teBC2_reciprocal_c <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      teBC2_reciprocal_c2)
    ReffBC1_c1 <- ifelse(Group_c == 1, teBC1_reciprocal_c1,
      NA)
    ReffBC1_c2 <- ifelse(Group_c == 2, teBC1_reciprocal_c2,
      NA)
    ReffBC2_c1 <- ifelse(Group_c == 1, teBC2_reciprocal_c1,
      NA)
    ReffBC2_c2 <- ifelse(Group_c == 2, teBC2_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u1_c = u1_c, u2_c = u2_c, teJLMS1_c = teJLMS1_c,
      teJLMS2_c = teJLMS2_c, teBC1_c = teBC1_c, teBC2_c = teBC2_c,
      teBC1_reciprocal_c = teBC1_reciprocal_c, teBC2_reciprocal_c = teBC2_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u1_c1 = u1_c1, u2_c1 = u2_c1, teBC1_c1 = teBC1_c1,
      teBC2_c1 = teBC2_c1, teBC1_reciprocal_c1 = teBC1_reciprocal_c1,
      teBC2_reciprocal_c1 = teBC2_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u1_c2 = u1_c2, u2_c2 = u2_c,
      teBC1_c2 = teBC1_c2, teBC2_c2 = teBC2_c2, teBC1_reciprocal_c2 = teBC1_reciprocal_c2,
      teBC2_reciprocal_c2 = teBC2_reciprocal_c2, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2,
      effBC1_c1 = effBC1_c1, effBC2_c1 = effBC2_c1, effBC1_c2 = effBC1_c2,
      effBC2_c2 = effBC2_c2, ReffBC1_c1 = ReffBC1_c1, ReffBC2_c1 = ReffBC2_c1,
      ReffBC1_c2 = ReffBC1_c2, ReffBC2_c2 = ReffBC2_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u1_c = u1_c, u2_c = u2_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u1_c1 = u1_c1, u2_c1 = u2_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u1_c2 = u1_c2, u2_c2 = u2_c, ineff1_c1 = ineff1_c1,
      ineff2_c1 = ineff2_c1, ineff1_c2 = ineff1_c2, ineff2_c2 = ineff2_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for cnsf uniform-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u
ccnsfmarguninorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

ccnsfmarguninorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu/2) * sqrt(12)) * (pnorm((exp(Wu/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# different sigma_u
cmcesfmarguninorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmcesfmarguninorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  Pi1 <- 1/(exp(Wu1/2) * sqrt(12)) * (pnorm((exp(Wu1/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv1/2)) - pnorm(object$S * epsilon/exp(Wv1/2)))
  Pi2 <- 1/(exp(Wu2/2) * sqrt(12)) * (pnorm((exp(Wu2/2) * sqrt(12) +
    object$S * epsilon)/exp(Wv2/2)) - pnorm(object$S * epsilon/exp(Wv2/2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}