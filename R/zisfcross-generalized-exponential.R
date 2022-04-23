################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Two types: - Common noise component (sigma_v)                                #
#            - Different noise component (multimodal noise - mnsf)             #
# Convolution: generalized exponential - normal                                #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf generalized_exponential-normal distribution
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
# Same sigma_v
czisfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# Different sigma_v
cmnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  A <- S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf generalized_exponential-normal distribution
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
# Same sigma_v
cstzisfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + generalized exponential - normal distributions...\n")
  initGenExpo <- maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgenexponormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initGenExpo$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("SF_",
    colnames(Zvar)))
  names(initGenExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Different sigma_v
cstmnsfgenexponorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  itermax, printInfo, tol) {
  cat("Initialization: SFA + generalized exponential - normal distributions...\n")
  initGenExpo <- maxLik(logLik = cgenexponormlike, start = cstgenexponorm(olsObj = olsObj,
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1,
      drop = FALSE]), grad = cgradgenexponormlike, method = "BFGS",
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0,
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[,
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]),
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initGenExpo$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("SF_", colnames(Zvar)))
  names(initGenExpo$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initGenExpo = initGenExpo))
}

# Gradient of the likelihood function ----------
#' gradient for zisf generalized_exponential-normal distribution
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
# Same sigma_v
cgradzisfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) *
    eB)
  sigx2 <- (2 * (sigx1 * ewz/(wzdeno * ewusr)) + S * prC *
    dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * ewz/(wzdeno *
    ewusr)))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 *
    ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (wzdeno * ewusr * (eA * pa - eB * pb)/(wzdeno *
    ewusr)^2)
  sigx8 <- (sigx6/(wzdeno * ewusr) - 0.5 * sigx7)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 *
    (S * (epsilon)/ewvsr)))
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prC/ewvsr + 2 *
    ((eA * sigx10 - sigx11 * eB) * ewz/(wzdeno * ewusr)))
  sigx13 <- ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) *
    (eA * pa - eB * pb))
  sigx14 <- (2 * sigx13 - prC * dwsr/(wzdeno * ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx8 *
    ewz/sigx3), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx12/sigx3,
    FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx14 *
    ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Different sigma_v
cgradmnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  mustar1 <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  eusq1 <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * epsilon/ewu_h))
  evsq1 <- (2 * (ewv1_h/ewu_h) + S * epsilon/ewv1_h)
  depsi <- dnorm(-evsq1, 0, 1)
  pepsi <- pnorm(-evsq1)
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewu_h)
  dpepsi <- (depsi/ewv1_h - 2 * (pepsi/ewu_h))
  epepe <- (eusq1 * pmustar1 - eusq2 * pepsi)
  exusq <- (sigx1 * eusq1 - dpepsi * eusq2)
  ewvu <- (ewu * ewv1/(2 * ewu)^2)
  sigx2 <- (exusq * ewz/(wzdeno * ewu_h))
  sigx3 <- (2 * sigx2 + S * prC * dmustar2 * epsilon/ewv2_h^3)
  sigx4 <- (epepe * ewz/(wzdeno * ewu_h))
  sigx5 <- (prC * dmustar2/ewv2_h + 2 * sigx4)
  sigx6 <- (0.5 * (S * epsilon/ewu_h) + 2 * ewvu)
  sigx7 <- (0.5 * (dmustar1 * ewv1_h/ewu_h) - sigx6 * pmustar1)
  sigx8 <- (depsi * ewv1_h/ewu_h - (2 * (ewv1/ewu) + S * epsilon/ewu_h) *
    pepsi)
  sigx9 <- (sigx7 * eusq1 - sigx8 * eusq2)
  sigx10 <- (wzdeno * ewu_h * epepe/(wzdeno * ewu_h)^2)
  sigx11 <- (sigx9/(wzdeno * ewu_h) - 0.5 * sigx10)
  sigx12 <- (ewv1 * pmustar1/(2 * ewu) - (0.5 * (ewv1_h/ewu_h) -
    0.5 * (S * epsilon/ewv1_h)) * dmustar1)
  sigx13 <- (2 * (ewv1 * pepsi/ewu) - depsi * (ewv1_h/ewu_h -
    0.5 * (S * epsilon/ewv1_h)))
  sigx14 <- (eusq1 * sigx12 - sigx13 * eusq2)
  sigx15 <- (0.5 * (S^2 * dmustar2 * epsilon^2/ewv2_h^2) -
    0.5 * dmustar2)
  sigx16 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx17 <- (sigx16 * epepe)
  sigx18 <- (2 * sigx17 - prC * dmustar2/(wzdeno * ewv2_h))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx5,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * (sigx11 *
    ewz/sigx5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 *
    (sigx14 * ewz/(sigx5 * wzdeno * ewu_h)), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = sigx15 * prC/(sigx5 *
      ewv2_h), FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx18 *
      ewz/sigx5, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf generalized_exponential-normal distribution
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
# Same sigma_v
chesszisfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv <- exp(Wv)
  ewvsr <- exp(Wv/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  da <- dnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  pa <- pnorm(-(ewvsr/ewusr + S * (epsilon)/ewvsr))
  db <- dnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  pb <- pnorm(-(2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr))
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(2 * (ewv/ewu) + 2 * (S * (epsilon)/ewusr))
  eC <- (2 * (ewv/ewu) + S * (epsilon)/ewusr)
  sigx1 <- ((da/ewvsr - pa/ewusr) * eA - (db/ewvsr - 2 * (pb/ewusr)) *
    eB)
  sigx2 <- (2 * (sigx1 * ewz/(wzdeno * ewusr)) + S * prC *
    dwsr * (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + 2 * ((eA * pa - eB * pb) * ewz/(wzdeno *
    ewusr)))
  sigx4 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 *
    ewu)^2))
  sigx5 <- (0.5 * (da * ewvsr/ewusr) - sigx4 * pa)
  sigx6 <- (sigx5 * eA - (db * ewvsr/ewusr - eC * pb) * eB)
  sigx7 <- (wzdeno * ewusr * (eA * pa - eB * pb)/(wzdeno *
    ewusr)^2)
  sigx8 <- (sigx6/(wzdeno * ewusr) - 0.5 * sigx7)
  depsisr <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * (ewvsr/ewusr) - 0.5 * (S * (epsilon)/ewvsr))
  sigx10 <- (ewv * pa/(2 * ewu) - sigx9 * da)
  sigx11 <- (2 * (ewv * pb/ewu) - db * (ewvsr/ewusr - 0.5 *
    (S * (epsilon)/ewvsr)))
  sigx12 <- ((0.5 * depsisr - 0.5 * dwsr) * prC/ewvsr + 2 *
    ((eA * sigx10 - sigx11 * eB) * ewz/(wzdeno * ewusr)))
  sigx13 <- ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) *
    (eA * pa - eB * pb))
  sigx14 <- (2 * sigx13 - prC * dwsr/(wzdeno * ewvsr))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 -
      1)/ewvsr^3 + 2 * (((((ewvsr/ewusr + S * (epsilon)/ewvsr)/ewvsr -
      1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) *
      eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewvsr -
      2/ewusr) * db/ewvsr - 2 * ((db/ewvsr - 2 * (pb/ewusr))/ewusr)) *
      eB) * ewz/(wzdeno * ewusr)) - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 *
      (S * (epsilon)/ewusr) + 2 * (ewu * ewv/(2 * ewu)^2)) *
      pa - 0.5 * (da * ewvsr/ewusr))/ewusr + (0.5 * ((ewvsr/ewusr +
      S * (epsilon)/ewvsr)/ewusr) - sigx4/ewvsr) * da) *
      eA - (((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr)/ewusr -
      eC/ewvsr) * db + (pb - 2 * (db * ewvsr/ewusr - eC *
      pb))/ewusr) * eB)/(wzdeno * ewusr) - (sigx8 * sigx2/sigx3 +
      0.5 * (sigx1 * wzdeno * ewusr/(wzdeno * ewusr)^2))) *
      ewz/sigx3), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (2 * (((da * (ewv/(2 * ewu) - (sigx9 * (ewvsr/ewusr +
    S * (epsilon)/ewvsr) + 0.5))/ewvsr - sigx10/ewusr) *
    eA - ((2 * (ewv/ewu) - ((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) *
    (ewvsr/ewusr - 0.5 * (S * (epsilon)/ewvsr)) + 0.5)) *
    db/ewvsr - 2 * (sigx11/ewusr)) * eB) * ewz/(wzdeno *
    ewusr)) + S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) -
    0.5) * prC * dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (sigx1 * (1/(wzdeno *
      ewusr) - ewusr * ewz/(wzdeno * ewusr)^2)) - (sigx2 *
      sigx14/sigx3 + S * prC * dwsr * (epsilon)/(wzdeno *
      ewvsr^3))) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((0.5 * (0.5 * (ewvsr * (ewvsr/ewusr + S * (epsilon)/ewvsr)/ewusr) -
    0.5) - 0.5 * sigx4) * da * ewvsr/ewusr - (sigx5 * sigx4 +
    (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv/(2 *
      ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * pa)) *
    eA - ((((2 * (ewvsr/ewusr) + S * (epsilon)/ewvsr) * ewvsr -
    S * (epsilon))/ewusr - (0.5 + 2 * (ewv/ewu))) * db *
    ewvsr/ewusr + (0.5 * (S * (epsilon)/ewusr) + 2 * (ewv/ewu)) *
    pb - eC * (db * ewvsr/ewusr - eC * pb)) * eB)/(wzdeno *
    ewusr) - ((0.5 * ((0.5 - wzdeno^2 * ewusr^2/(wzdeno *
    ewusr)^2) * (eA * pa - eB * pb) + sigx5 * eA - (db *
    ewvsr/ewusr - eC * pb) * eB) + 0.5 * sigx6) * wzdeno *
    ewusr/(wzdeno * ewusr)^2 + 2 * (sigx8^2 * ewz/sigx3))) *
    ewz/sigx3), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * ((((da * ewvsr/(4 *
      (ewu * ewusr)) - 2 * (ewu * pa/(2 * ewu)^2)) * ewv -
      ((0.5 * (sigx9 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) -
        0.25) * da * ewvsr/ewusr + sigx4 * sigx10)) *
      eA - (2 * ((db * ewvsr/ewusr - pb) * ewv/ewu) - (((2 *
      (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
      0.5 * (S * (epsilon)/ewvsr)) - 0.5) * db * ewvsr/ewusr +
      sigx11 * eC)) * eB)/(wzdeno * ewusr) - 0.5 * (wzdeno *
      ewusr * (eA * sigx10 - sigx11 * eB)/(wzdeno * ewusr)^2)) -
      2 * (sigx8 * sigx12/sigx3)) * ewz/sigx3, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx6 * (1/(wzdeno *
      ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) - ((0.5 -
      wzdeno^2 * ewusr^2/(wzdeno * ewusr)^2) * ewz + 0.5 *
      wzdeno) * ewusr * (eA * pa - eB * pb)/(wzdeno * ewusr)^2) -
      2 * (sigx8 * sigx14 * ewz/sigx3)) * ewz/sigx3, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (prC * (S^2 * (0.5 * (0.5 *
      (S^2 * (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr *
      (epsilon)^2/ewvsr^2 - 0.5 * (0.5 * depsisr - 0.5 *
      dwsr))/ewvsr + 2 * ((((sigx10/2 + (pa - sigx9 * da)/2) *
      ewv/ewu - (0.25 * (ewvsr/ewusr) + 0.25 * (S * (epsilon)/ewvsr) -
      sigx9^2 * (ewvsr/ewusr + S * (epsilon)/ewvsr)) *
      da) * eA - ((2 * sigx11 + 2 * (pb - db * (ewvsr/ewusr -
      0.5 * (S * (epsilon)/ewvsr)))) * ewv/ewu - (0.25 *
      (S * (epsilon)/ewvsr) + 0.5 * (ewvsr/ewusr) - (2 *
      (ewvsr/ewusr) + S * (epsilon)/ewvsr) * (ewvsr/ewusr -
      0.5 * (S * (epsilon)/ewvsr))^2) * db) * eB) * ewz/(wzdeno *
      ewusr)) - sigx12^2/sigx3)/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + nvZVvar +
      nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (2 * ((1/(wzdeno * ewusr) - ewusr * ewz/(wzdeno * ewusr)^2) *
      (eA * sigx10 - sigx11 * eB)) - (sigx12 * sigx14/sigx3 +
      (0.5 * (S^2 * dwsr * (epsilon)^2/(wzdeno * ewvsr^3)) -
        0.5 * (wzdeno * dwsr * ewvsr/(wzdeno * ewvsr)^2)) *
        prC)) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    nvZVvar + nZHvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewvsr) +
      ewvsr/(wzdeno * ewvsr)^2) * dwsr - (sigx14^2/sigx3 +
      2 * ((2 - 2 * (wzdeno * ewusr^2 * ewz/(wzdeno * ewusr)^2)) *
        ewusr * (eA * pa - eB * pb)/(wzdeno * ewusr)^2))) *
      ewz + 2 * sigx13 - prC * dwsr/(wzdeno * ewvsr)) *
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Different sigma_v
chessmnsfgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
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
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  mustar1 <- (ewv1_h/ewu_h + S * epsilon/ewv1_h)
  dmustar1 <- dnorm(-mustar1, 0, 1)
  dmustar2 <- dnorm(S * epsilon/ewv2_h, 0, 1)
  pmustar1 <- pnorm(-mustar1)
  eusq1 <- exp(ewv1/(2 * ewu) + S * epsilon/ewu_h)
  eusq2 <- exp(2 * (ewv1/ewu) + 2 * (S * epsilon/ewu_h))
  evsq1 <- (2 * (ewv1_h/ewu_h) + S * epsilon/ewv1_h)
  depsi <- dnorm(-evsq1, 0, 1)
  pepsi <- pnorm(-evsq1)
  sigx1 <- (dmustar1/ewv1_h - pmustar1/ewu_h)
  dpepsi <- (depsi/ewv1_h - 2 * (pepsi/ewu_h))
  epepe <- (eusq1 * pmustar1 - eusq2 * pepsi)
  exusq <- (sigx1 * eusq1 - dpepsi * eusq2)
  ewvu <- (ewu * ewv1/(2 * ewu)^2)
  sigx2 <- (exusq * ewz/(wzdeno * ewu_h))
  sigx3 <- (2 * sigx2 + S * prC * dmustar2 * epsilon/ewv2_h^3)
  sigx4 <- (epepe * ewz/(wzdeno * ewu_h))
  sigx5 <- (prC * dmustar2/ewv2_h + 2 * sigx4)
  sigx6 <- (0.5 * (S * epsilon/ewu_h) + 2 * ewvu)
  sigx7 <- (0.5 * (dmustar1 * ewv1_h/ewu_h) - sigx6 * pmustar1)
  sigx8 <- (depsi * ewv1_h/ewu_h - (2 * (ewv1/ewu) + S * epsilon/ewu_h) *
    pepsi)
  sigx9 <- (sigx7 * eusq1 - sigx8 * eusq2)
  sigx10 <- (wzdeno * ewu_h * epepe/(wzdeno * ewu_h)^2)
  sigx11 <- (sigx9/(wzdeno * ewu_h) - 0.5 * sigx10)
  sigx12 <- (ewv1 * pmustar1/(2 * ewu) - (0.5 * (ewv1_h/ewu_h) -
    0.5 * (S * epsilon/ewv1_h)) * dmustar1)
  sigx13 <- (2 * (ewv1 * pepsi/ewu) - depsi * (ewv1_h/ewu_h -
    0.5 * (S * epsilon/ewv1_h)))
  sigx14 <- (eusq1 * sigx12 - sigx13 * eusq2)
  sigx15 <- (0.5 * (S^2 * dmustar2 * epsilon^2/ewv2_h^2) -
    0.5 * dmustar2)
  sigx16 <- (1/(wzdeno * ewu_h) - ewu_h * ewz/(wzdeno * ewu_h)^2)
  sigx17 <- (sigx16 * epepe)
  sigx18 <- (2 * sigx17 - prC * dmustar2/(wzdeno * ewv2_h))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * dmustar2 * (S^2 * epsilon^2/ewv2_h^2 -
      1)/ewv2_h^3 + 2 * ((((mustar1/ewv1_h - 1/ewu_h) *
      dmustar1/ewv1_h - sigx1/ewu_h) * eusq1 - ((evsq1/ewv1_h -
      2/ewu_h) * depsi/ewv1_h - 2 * (dpepsi/ewu_h)) * eusq2) *
      ewz/(wzdeno * ewu_h)) - sigx3^2/sigx5)/sigx5, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * 2 * (S * (((((0.5 + 0.5 *
      (S * epsilon/ewu_h) + 2 * ewvu) * pmustar1 - 0.5 *
      (dmustar1 * ewv1_h/ewu_h))/ewu_h + (0.5 * (mustar1/ewu_h) -
      sigx6/ewv1_h) * dmustar1) * eusq1 - ((evsq1/ewu_h -
      (2 * (ewv1/ewu) + S * epsilon/ewu_h)/ewv1_h) * depsi +
      (pepsi - 2 * sigx8)/ewu_h) * eusq2)/(wzdeno * ewu_h) -
      (sigx11 * sigx3/sigx5 + 0.5 * (exusq * wzdeno * ewu_h/(wzdeno *
        ewu_h)^2))) * ewz/sigx5), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    2 * (S * (((dmustar1 * (ewv1/(2 * ewu) - ((0.5 * (ewv1_h/ewu_h) -
    0.5 * (S * epsilon/ewv1_h)) * mustar1 + 0.5))/ewv1_h -
    sigx12/ewu_h) * eusq1 - ((2 * (ewv1/ewu) - (evsq1 * (ewv1_h/ewu_h -
    0.5 * (S * epsilon/ewv1_h)) + 0.5)) * depsi/ewv1_h -
    2 * (sigx13/ewu_h)) * eusq2)/(sigx5 * wzdeno * ewu_h) -
    wzdeno * sigx3 * ewu_h * sigx14/(sigx5 * wzdeno * ewu_h)^2) *
    ewz), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prC * (S * (0.5 * (S^2 * epsilon^2/ewv2_h^2 -
      2) - 0.5) * dmustar2 * epsilon/(sigx5 * ewv2_h^3) -
      sigx15 * sigx3 * ewv2_h/(sigx5 * ewv2_h)^2), FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (2 * (exusq * sigx16) -
      (sigx3 * sigx18/sigx5 + S * prC * dmustar2 * epsilon/(wzdeno *
        ewv2_h^3))) * ewz/sigx5, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    2 * (((((0.5 * (0.5 * (ewv1_h * mustar1/ewu_h) - 0.5) -
    0.5 * sigx6) * dmustar1 * ewv1_h/ewu_h - (sigx7 * sigx6 +
    (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 *
      ewu)^2) - 0.25 * (S * epsilon/ewu_h)) * pmustar1)) *
    eusq1 - (((evsq1 * ewv1_h - S * epsilon)/ewu_h - (0.5 +
    2 * (ewv1/ewu))) * depsi * ewv1_h/ewu_h + (0.5 * (S *
    epsilon/ewu_h) + 2 * (ewv1/ewu)) * pepsi - (2 * (ewv1/ewu) +
    S * epsilon/ewu_h) * sigx8) * eusq2)/(wzdeno * ewu_h) -
    ((0.5 * ((0.5 - wzdeno^2 * ewu_h^2/(wzdeno * ewu_h)^2) *
      epepe + sigx7 * eusq1 - sigx8 * eusq2) + 0.5 * sigx9) *
      wzdeno * ewu_h/(wzdeno * ewu_h)^2 + 2 * (sigx11^2 *
      ewz/sigx5))) * ewz/sigx5), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((dmustar1 * ewv1_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pmustar1/(2 * ewu)^2)) *
      ewv1 - ((0.5 * ((0.5 * (ewv1_h/ewu_h) - 0.5 * (S *
      epsilon/ewv1_h)) * mustar1) - 0.25) * dmustar1 *
      ewv1_h/ewu_h + sigx6 * sigx12)) * eusq1 - (2 * ((depsi *
      ewv1_h/ewu_h - pepsi) * ewv1/ewu) - ((evsq1 * (ewv1_h/ewu_h -
      0.5 * (S * epsilon/ewv1_h)) - 0.5) * depsi * ewv1_h/ewu_h +
      sigx13 * (2 * (ewv1/ewu) + S * epsilon/ewu_h))) *
      eusq2)/(sigx5 * wzdeno * ewu_h) - (0.5 * sigx5 +
      2 * (sigx11 * ewz)) * wzdeno * ewu_h * sigx14/(sigx5 *
      wzdeno * ewu_h)^2) * ewz), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (2 * (sigx11 * sigx15 *
      prC * ewv2_h * ewz/(sigx5 * ewv2_h)^2)), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (2 * (sigx9 * sigx16 - ((0.5 -
      wzdeno^2 * ewu_h^2/(wzdeno * ewu_h)^2) * ewz + 0.5 *
      wzdeno) * ewu_h * epepe/(wzdeno * ewu_h)^2) - 2 *
      (sigx11 * sigx18 * ewz/sigx5)) * ewz/sigx5, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * 2 * (((((sigx12/2 + (pmustar1 -
      (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * epsilon/ewv1_h)) *
        dmustar1)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * epsilon/ewv1_h) - (0.5 * (ewv1_h/ewu_h) -
      0.5 * (S * epsilon/ewv1_h))^2 * mustar1) * dmustar1) *
      eusq1 - ((2 * sigx13 + 2 * (pepsi - depsi * (ewv1_h/ewu_h -
      0.5 * (S * epsilon/ewv1_h)))) * ewv1/ewu - (0.25 *
      (S * epsilon/ewv1_h) + 0.5 * (ewv1_h/ewu_h) - evsq1 *
      (ewv1_h/ewu_h - 0.5 * (S * epsilon/ewv1_h))^2) *
      depsi) * eusq2)/(sigx5 * wzdeno * ewu_h) - 2 * (sigx14^2 *
      ewz/(sigx5 * wzdeno * ewu_h)^2)) * ewz), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (2 * (sigx15 * prC * sigx14 * ewv2_h * ewz/((sigx5 *
      ewv2_h)^2 * wzdeno * ewu_h))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (2 * sigx16 - 2 * (sigx18 *
      ewz/(sigx5 * wzdeno * ewu_h))) * sigx14 * ewz/sigx5,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (S^2 * (0.5 * (0.5 * (S^2 * epsilon^2/ewv2_h^2) -
      1) - 0.25) * dmustar2 * epsilon^2/(sigx5 * ewv2_h^3) -
      (sigx15 * prC + 0.5 * (sigx5 * ewv2_h)) * sigx15/(sigx5 *
        ewv2_h)^2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx15 * sigx18/sigx5 +
      0.5 * (S^2 * dmustar2 * epsilon^2/(wzdeno * ewv2_h^2)))/ewv2_h -
      0.5 * (wzdeno * dmustar2 * ewv2_h/(wzdeno * ewv2_h)^2)) *
      prC * ewz/sigx5), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * ewv2_h) +
      ewv2_h/(wzdeno * ewv2_h)^2) * dmustar2 - (sigx18^2/sigx5 +
      2 * ((2 - 2 * (wzdeno * ewu_h^2 * ewz/(wzdeno * ewu_h)^2)) *
        ewu_h * epepe/(wzdeno * ewu_h)^2))) * ewz + 2 *
      sigx17 - prC * dmustar2/(wzdeno * ewv2_h)) * ewz/sigx5,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf generalized_exponential-normal distribution
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
# Same sigma_v
zisfgenexponormAlgOpt <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstzisfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(czisfgenexponormlike(startVal, nXvar = nXvar,
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
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(czisfgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = czisfgenexponormlike, grad = cgradzisfgenexponormlike,
      hess = chesszisfgenexponormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(czisfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(czisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesszisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(czisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradzisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesszisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(czisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradzisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesszisfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisfgenexponormlike(mleObj$par,
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
      mleObj$hessian <- chesszisfgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesszisfgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisfgenexponormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisfgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGenExpo = initGenExpo))
}

# Different sigma_v
mnsfgenexponormAlgOpt <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmnsfgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax,
    printInfo = printInfo, tol = tol)
  initGenExpo <- start_st$initGenExpo
  startVal <- start_st$StartVal
  startLoglik <- sum(cmnsfgenexponormlike(startVal, nXvar = nXvar,
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
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cmnsfgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmnsfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmnsfgenexponormlike, grad = cgradmnsfgenexponormlike,
      hess = chessmnsfgenexponormlike, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cmnsfgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmnsfgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmnsfgenexponormlike(mleObj$par,
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
      mleObj$hessian <- chessmnsfgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmnsfgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmnsfgenexponormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmnsfgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGenExpo = initGenExpo))
}

# Conditional efficiencies estimation ----------
#' efficiencies for zisf genexpo-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_v
czisfgenexponormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) -
    exp(B) * (dnorm(b) + b * pnorm(b)))/(exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) *
      pnorm(a - exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) +
      exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(exp(A) * pnorm(a) -
      exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv/2) + exp(Wv)/2) *
      pnorm(a + exp(Wv/2)) - exp(B) * exp(b * exp(Wv/2) +
      exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(exp(A) * pnorm(a) -
      exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Different sigma_v
cmnsfgenexponormeff <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv1/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) -
    exp(B) * (dnorm(b) + b * pnorm(b)))/(exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (exp(A) * exp(-a * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a - exp(Wv1/2)) - exp(B) * exp(-b * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b - exp(Wv1/2)))/(exp(A) * pnorm(a) -
      exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (exp(A) * exp(a * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a + exp(Wv1/2)) - exp(B) * exp(b *
      exp(Wv1/2) + exp(Wv1)/2) * pnorm(b + exp(Wv1/2)))/(exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    teBC_reciprocal_c2 <- rep(1, object$Nobs)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for zisf genexpo-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_v
czisfmarggenexponorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmarggenexponorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + object$nuZUvar +
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Zvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

# Different sigma_v
cmnsfmarggenexponorm_Eu <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

cmnsfmarggenexponorm_Vu <- function(object) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv1)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv1/2) - 2 * exp(Wv1/2)/exp(Wu/2)
  Pi1 <- 2/exp(Wu/2) * (exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv2/2) * dnorm(object$S * epsilon/exp(Wv2/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar -
    1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
