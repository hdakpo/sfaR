################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Zero Inefficiency Stochastic Frontier Model                           #
# Convolution: truncated skewed laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for zisf truncated_skewed_laplace-normal distribution
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
czisftslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, 
  vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + 
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon * 
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0) 
    return(NA)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * 
    exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(wHvar * log(L)))
}

# starting value for the log-likelihood ----------
#' starting values for zisf truncated_skewed_laplace-normal distribution
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
cstzisftslnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, itermax, printInfo, 
  tol) {
  cat("Initialization: SFA + truncated skewed laplace - Normal distributions...\n")
  initTSL <- maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj, 
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 
      1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, 
      drop = FALSE]), grad = cgradtslnormlike, method = "BFGS", 
    control = list(iterlim = itermax, printLevel = if (printInfo) 2 else 0, 
      reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = as.matrix(uHvar[, 
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1]), 
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)
  Esti <- initTSL$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), "lambda", 
    paste0("SF_", colnames(Zvar)))
names(initTSL$estimate) <- c(names(Esti)[1:nXvar], paste0(
    "Zu_",
    colnames(uHvar)[1]
  ), paste0("Zv_", colnames(vHvar)[1]), "lambda")
  return(list(StartVal = StartVal, initTSL = initTSL))
}

# Gradient of the likelihood function ----------
#' gradient for zisf truncated_skewed_laplace-normal distribution
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
cgradzisftslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + 
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  wvwu <- ewvsr/ewusr
  lwv <- ((1 + lambda) * wvwu + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (wvwu + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * 
    (1 + lambda))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzl <- (lwu * wzdeno)
  lwzu <- (1/wzl - lwu * ewz/wzl^2)
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * wvwu) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - 
    (1 + lambda) * pd/ewusr) * eB)
  sigx2 <- ((1 + lambda) * sigx1 * ewz/wzl + S * prC * dwsr * 
    (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + (1 + lambda) * eaeb * ewz/wzl)
  sigx4 <- ((0.5 * (da * wvwu) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * 
    ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - (0.5 * (dd * wvwu) - sigx5 * pd) * 
    (1 + lambda) * eB)
  sigx7 <- (sigx6/wzl - lambda * wzdeno * eaeb * ewusr/wzl^2)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * sigx8 - 0.5 * dwsr) * prC/ewvsr
  sigx10 <- (eA * pawv)
  sigx11 <- ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * 
    dd)
  sigx12 <- (sigx9 + (1 + lambda) * (2 * sigx10 - sigx11 * 
    eB) * ewz/wzl)
  sigx13 <- (2 * (eA * pa) - (llwusq + pd) * eB)
  sigx14 <- (sigx13/wzl - 2 * (wzdeno * (1 + lambda) * eaeb * 
    ewusr/wzl^2))
  sigx15 <- ((1 + lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * 
    ewvsr))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, 
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx7 * 
    (1 + lambda) * ewz/sigx3, FUN = "*"), sweep(vHvar, MARGIN = 1, 
    STATS = sigx12/sigx3, FUN = "*"), sigx14 * ewz/sigx3, 
    sweep(Zvar, MARGIN = 1, STATS = sigx15 * ewz/sigx3, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for zisf truncated_skewed_laplace-normal distribution
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
chesszisftslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + 
    nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusr <- exp(Wu/2)
  ewvsr <- exp(Wv/2)
  ewu <- exp(Wu)
  ewv <- exp(Wv)
  ewz <- exp(Wz)
  wvwu <- ewvsr/ewusr
  lwv <- ((1 + lambda) * wvwu + S * (epsilon)/ewvsr)
  lwu <- (1 + 2 * (lambda * ewusr))
  epsiwv <- (wvwu + S * (epsilon)/ewvsr)
  da <- dnorm(-epsiwv)
  pa <- pnorm(-epsiwv)
  pd <- pnorm(-lwv)
  dd <- dnorm(-lwv)
  dwsr <- dnorm(S * (epsilon)/ewvsr)
  eA <- exp(ewv/(2 * ewu) + S * (epsilon)/ewusr)
  eB <- exp(((1 + lambda) * ewv/(2 * ewu) + S * (epsilon)/ewusr) * 
    (1 + lambda))
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wzl <- (lwu * wzdeno)
  lwzu <- (1/wzl - lwu * ewz/wzl^2)
  ewuwv <- (ewu * ewv/(2 * ewu)^2)
  epsiwu <- (0.5 * (S * (epsilon)/ewusr) + 2 * ewuwv)
  epsiwr <- (0.5 * (wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  pawv <- (ewv * pa/(2 * ewu) - epsiwr * da)
  lwusq <- ((1 + lambda) * ewv/ewu + S * (epsilon)/ewusr)
  llwusq <- (lwusq * pd - dd * wvwu) * (1 + lambda)
  llwvsq <- (0.5 * ((1 + lambda) * wvwu) - 0.5 * (S * (epsilon)/ewvsr))
  eaeb <- (2 * (eA * pa) - eB * pd)
  sigx1 <- (2 * ((da/ewvsr - pa/ewusr) * eA) - (dd/ewvsr - 
    (1 + lambda) * pd/ewusr) * eB)
  sigx2 <- ((1 + lambda) * sigx1 * ewz/wzl + S * prC * dwsr * 
    (epsilon)/ewvsr^3)
  sigx3 <- (prC * dwsr/ewvsr + (1 + lambda) * eaeb * ewz/wzl)
  sigx4 <- ((0.5 * (da * wvwu) - epsiwu * pa) * eA)
  sigx5 <- (0.5 * (S * (epsilon)/ewusr) + 2 * ((1 + lambda) * 
    ewu * ewv/(2 * ewu)^2))
  sigx6 <- (2 * sigx4 - (0.5 * (dd * wvwu) - sigx5 * pd) * 
    (1 + lambda) * eB)
  sigx7 <- (sigx6/wzl - lambda * wzdeno * eaeb * ewusr/wzl^2)
  sigx8 <- (S^2 * dwsr * (epsilon)^2/ewvsr^2)
  sigx9 <- (0.5 * sigx8 - 0.5 * dwsr) * prC/ewvsr
  sigx10 <- (eA * pawv)
  sigx11 <- ((1 + lambda)^2 * ewv * pd/(2 * ewu) - llwvsq * 
    dd)
  sigx12 <- (sigx9 + (1 + lambda) * (2 * sigx10 - sigx11 * 
    eB) * ewz/wzl)
  sigx13 <- (2 * (eA * pa) - (llwusq + pd) * eB)
  sigx14 <- (sigx13/wzl - 2 * (wzdeno * (1 + lambda) * eaeb * 
    ewusr/wzl^2))
  sigx15 <- ((1 + lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * 
    ewvsr))
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + nZHvar + 
    1, ncol = nXvar + nuZUvar + nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = wHvar*S^2 * (prC * dwsr * (S^2 * (epsilon)^2/ewvsr^2 - 
      1)/ewvsr^3 + (1 + lambda) * (2 * (((epsiwv/ewvsr - 
      1/ewusr) * da/ewvsr - (da/ewvsr - pa/ewusr)/ewusr) * 
      eA) - ((lwv/ewvsr - (1 + lambda)/ewusr) * dd/ewvsr - 
      (1 + lambda) * (dd/ewvsr - (1 + lambda) * pd/ewusr)/ewusr) * 
      eB) * ewz/wzl - sigx2^2/sigx3)/sigx3, FUN = "*"), 
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * ((2 * ((((0.5 + 0.5 * (S * (epsilon)/ewusr) + 
      2 * ewuwv) * pa - 0.5 * (da * wvwu))/ewusr + (0.5 * 
      (epsiwv/ewusr) - epsiwu/ewvsr) * da) * eA) - ((0.5 * 
      (lwv/ewusr) - sigx5/ewvsr) * dd + (0.5 * pd - (0.5 * 
      (dd * wvwu) - sigx5 * pd) * (1 + lambda))/ewusr) * 
      (1 + lambda) * eB)/wzl - (sigx2 * sigx7/sigx3 + lambda * 
      wzdeno * sigx1 * ewusr/wzl^2)) * (1 + lambda) * ewz/sigx3, 
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + 
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar*S * 
    ((1 + lambda) * (2 * ((da * (ewv/(2 * ewu) - (epsiwr * 
      epsiwv + 0.5))/ewvsr - pawv/ewusr) * eA) - (((1 + 
      lambda)^2 * ewv/(2 * ewu) - (lwv * llwvsq + 0.5)) * 
      dd/ewvsr - sigx11 * (1 + lambda)/ewusr) * eB) * ewz/wzl + 
      S * (0.5 * (S^2 * (epsilon)^2/ewvsr^2 - 2) - 0.5) * 
        prC * dwsr * (epsilon)/ewvsr^3 - sigx12 * sigx2/sigx3)/sigx3, 
    FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- matrix(colSums(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * ((2 * ((da/ewvsr - pa/ewusr) * 
      eA) - (((lwusq/ewvsr - lwv/ewusr) * dd - (llwusq + 
      2 * pd)/ewusr) * (1 + lambda) + dd/ewvsr) * eB)/wzl - 
      (sigx2 * sigx14/sigx3 + 2 * (wzdeno * (1 + lambda) * 
        sigx1 * ewusr/wzl^2))) * ewz/sigx3, FUN = "*")), 
    ncol = 1)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 2):(nXvar + 
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = wHvar*S * ((1 + lambda) * lwzu * sigx1 - 
      (sigx15 * sigx2/sigx3 + S * prC * dwsr * (epsilon)/(wzdeno * 
        ewvsr^3))) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + 
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar*((2 * 
    (((0.5 * (0.5 * (ewvsr * epsiwv/ewusr) - 0.5) - 0.5 * 
      epsiwu) * da * wvwu - ((0.5 * (da * wvwu) - epsiwu * 
      pa) * epsiwu + (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * 
      ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * 
      pa)) * eA) - ((0.5 * (0.5 * (lwv * (1 + lambda) * 
    wvwu) - 0.5) - 0.5 * (sigx5 * (1 + lambda))) * dd * wvwu - 
    ((0.5 * (dd * wvwu) - sigx5 * pd) * sigx5 * (1 + lambda) + 
      (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * 
        ewu * ewv/(2 * ewu)^2) - 0.25 * (S * (epsilon)/ewusr)) * 
        pd)) * (1 + lambda) * eB)/wzl - (sigx7^2 * (1 + 
    lambda) * ewz/sigx3 + lambda * ((0.5 - 2 * (lambda * 
    lwu * wzdeno^2 * ewusr/wzl^2)) * eaeb + 4 * sigx4 - 2 * 
    ((0.5 * (dd * wvwu) - sigx5 * pd) * (1 + lambda) * eB)) * 
    wzdeno * ewusr/wzl^2)) * (1 + lambda) * ewz/sigx3, FUN = "*"), 
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = wHvar*((2 * (((da * ewvsr/(4 * (ewu * ewusr)) - 
      2 * (ewu * pa/(2 * ewu)^2)) * ewv - ((0.5 * (epsiwr * 
      epsiwv) - 0.25) * da * wvwu + epsiwu * pawv)) * eA) - 
      (((1 + lambda) * dd * ewvsr/(4 * (ewu * ewusr)) - 
        2 * (ewu * pd/(2 * ewu)^2)) * (1 + lambda) * 
        ewv - (sigx11 * sigx5 + (0.5 * (lwv * llwvsq) - 
        0.25) * dd * wvwu)) * (1 + lambda) * eB)/wzl - 
      (sigx12 * sigx7/sigx3 + lambda * wzdeno * (2 * sigx10 - 
        sigx11 * eB) * ewusr/wzl^2)) * (1 + lambda) * 
      ewz/sigx3, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar + 
    1] <- matrix(colSums(sweep(uHvar, MARGIN = 1, STATS = wHvar*((2 * 
    sigx4 - (((0.5 * lwusq - 0.5 * (lwv * wvwu)) * (1 + lambda) + 
    1) * dd * wvwu - ((llwusq + pd) * sigx5 + ((1 + lambda) * 
    ewv/ewu + 0.5 * (S * (epsilon)/ewusr)) * pd)) * (1 + 
    lambda) * eB)/wzl - (sigx7 * sigx14 * (1 + lambda) * 
    ewz/sigx3 + wzdeno * (2 * (((0.5 - 2 * (lambda * lwu * 
    wzdeno^2 * ewusr/wzl^2)) * eaeb + 2 * sigx4 - (0.5 * 
    (dd * wvwu) - sigx5 * pd) * (1 + lambda) * eB) * (1 + 
    lambda)) + lambda * sigx13) * ewusr/wzl^2)) * ewz/sigx3, 
    FUN = "*")), ncol = 1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS =wHvar* (lwzu * sigx6 - (sigx15 * sigx7 * 
      ewz/sigx3 + lambda * ((2 - 2 * (lwu^2 * wzdeno^2/wzl^2)) * 
      ewz + 1) * eaeb * ewusr/wzl^2)) * (1 + lambda) * 
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = wHvar*(prC * (S^2 * (0.5 * (0.5 * (S^2 * 
      (epsilon)^2/ewvsr^2) - 1) - 0.25) * dwsr * (epsilon)^2/ewvsr^2 - 
      0.5 * (0.5 * sigx8 - 0.5 * dwsr))/ewvsr + (1 + lambda) * 
      (2 * (((pawv/2 + (pa - epsiwr * da)/2) * ewv/ewu - 
        (0.25 * (wvwu) + 0.25 * (S * (epsilon)/ewvsr) - 
          epsiwr^2 * epsiwv) * da) * eA) - ((sigx11/2 + 
        (pd - llwvsq * dd)/2) * (1 + lambda)^2 * ewv/ewu - 
        (0.25 * ((1 + lambda) * wvwu) + 0.25 * (S * (epsilon)/ewvsr) - 
          lwv * llwvsq^2) * dd) * eB) * ewz/wzl - sigx12^2/sigx3)/sigx3, 
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    nXvar + nuZUvar + nvZVvar + 1] <- matrix(colSums(sweep(vHvar, 
    MARGIN = 1, STATS = wHvar*((2 * sigx10 - ((((llwusq + pd)/2 + 
      pd) * (1 + lambda) * ewv/ewu - (lwusq * llwvsq + 
      (0.5 - lwv * llwvsq) * wvwu) * dd) * (1 + lambda) - 
      llwvsq * dd) * eB)/wzl - (sigx12 * sigx14/sigx3 + 
      2 * (wzdeno * (1 + lambda) * (2 * sigx10 - sigx11 * 
        eB) * ewusr/wzl^2))) * ewz/sigx3, FUN = "*")), 
    ncol = 1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + 
      nZHvar + 1)] <- crossprod(sweep(vHvar, MARGIN = 1, 
    STATS = wHvar*((1 + lambda) * lwzu * (2 * sigx10 - sigx11 * 
      eB) - (sigx12 * sigx15/sigx3 + (0.5 * (S^2 * dwsr * 
      (epsilon)^2/(wzdeno * ewvsr^3)) - 0.5 * (wzdeno * 
      dwsr * ewvsr/(wzdeno * ewvsr)^2)) * prC)) * ewz/sigx3, 
    FUN = "*"), Zvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar + 
    1] <- sum(-((((((lwv * ewvsr - S * (epsilon))/ewusr - 
    (1 + lambda) * ewv/ewu) * dd * wvwu + ewv * pd/ewu) * 
    (1 + lambda) + (llwusq + 2 * pd) * lwusq - 2 * (dd * 
    wvwu)) * eB/wzl + sigx14^2 * ewz/sigx3 + wzdeno * (2 * 
    (2 * (eA * pa) - ((llwusq + pd) * eB + 4 * (lwu * wzdeno^2 * 
      (1 + lambda) * eaeb * ewusr/wzl^2))) + 2 * sigx13) * 
    ewusr/wzl^2) * ewz/sigx3)*wHvar)
  hessll[nXvar + nuZUvar + nvZVvar + 1, (nXvar + nuZUvar + 
    nvZVvar + 2):(nXvar + nuZUvar + nvZVvar + nZHvar + 1)] <- matrix(colSums(sweep(Zvar, 
    MARGIN = 1, STATS = wHvar*((1/wzl - (((2 - 4 * (lwu^2 * wzdeno^2/wzl^2)) * 
      ewz + 2 * wzdeno) * (1 + lambda) * ewusr + lwu * 
      ewz)/wzl^2) * eaeb - (llwusq * lwzu * eB + sigx15 * 
      sigx14 * ewz/sigx3)) * ewz/sigx3, FUN = "*")), ncol = 1)
  hessll[(nXvar + nuZUvar + nvZVvar + 2):(nXvar + nuZUvar + 
    nvZVvar + nZHvar + 1), (nXvar + nuZUvar + nvZVvar + 2):(nXvar + 
    nuZUvar + nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = wHvar*((prC * (1/(wzdeno^2 * ewvsr) + ewvsr/(wzdeno * 
      ewvsr)^2) * dwsr - (sigx15^2/sigx3 + lwu * (1 + lambda) * 
      (2 - 2 * (lwu^2 * wzdeno * ewz/wzl^2)) * eaeb/wzl^2)) * 
      ewz + (1 + lambda) * lwzu * eaeb - prC * dwsr/(wzdeno * 
      ewvsr)) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for zisf truncated_skewed_laplace-normal distribution
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
zisftslnormAlgOpt <- function(start, olsParam, dataTable, S,wHvar,  
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, 
  Xvar, method, printInfo, itermax, stepmax, tol, gradtol, 
  hessianType, qac) {
  start_st <- if (!is.null(start)) 
    start else cstzisftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], 
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, 
    vHvar = vHvar, nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, 
    Yvar = Yvar, itermax = itermax, printInfo = printInfo, 
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(czisftslnormlike(startVal, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, 
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
                         bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
                         nm = function(...) maxNM(...), cg = function(...) maxCG(...),
                         sann = function(...) maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("ZISF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf(par = startVal, 
    fn = function(parm) -sum(czisftslnormlike(parm, nXvar = nXvar, 
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, 
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, 
    control = list(trace = printInfo, maxeval = itermax, 
      stepmax = stepmax, xtol = tol, grtol = gradtol)), 
    maxLikAlgo = maxRoutine(fn = czisftslnormlike, grad = cgradzisftslnormlike, 
      hess = chesszisftslnormlike, start = startVal, finalHessian = if (hessianType == 
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0, 
        iterlim = itermax, reltol = tol, tol = tol, qac = qac), 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trust.optim(x = startVal, 
      fn = function(parm) -sum(czisftslnormlike(parm, nXvar = nXvar, 
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, 
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", 
      control = list(maxit = itermax, cgtol = gradtol, 
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0, 
        report.precision = 1L)), sparse = trust.optim(x = startVal, 
      fn = function(parm) -sum(czisftslnormlike(parm, nXvar = nXvar, 
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, 
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chesszisftslnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), 
      method = "Sparse", control = list(maxit = itermax, 
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, 
        report.level = if (printInfo) 4L else 0, report.precision = 1L, 
        preconditioner = 1L)), mla = mla(b = startVal, 
      fn = function(parm) -sum(czisftslnormlike(parm, nXvar = nXvar, 
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, 
        Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradzisftslnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chesszisftslnormlike(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, 
      maxiter = itermax, epsa = gradtol, epsb = gradtol), 
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(czisftslnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradzisftslnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesszisftslnormlike(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, 
      trace = printInfo, eval.max = itermax, rel.tol = tol, 
      x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradzisftslnormlike(mleObj$par, 
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
      mleObj$hessian <- chesszisftslnormlike(parm = mleObj$par, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1") 
      mleObj$hessian <- chesszisftslnormlike(parm = mleObj$solution, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- czisftslnormlike(parm = mlParam, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, Zvar = Zvar, 
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradzisftslnormlike(parm = mlParam, 
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, 
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

# Conditional efficiencies estimation ----------

czisftslnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
                                                     object$nvZVvar + 1]
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
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * 
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * 
                                                          exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) -
                      exp(B) * (dnorm(b) + b * pnorm(b)))/(2 * exp(A) * pnorm(a) -
                                                             exp(B) * pnorm(b))
  u_c2 <- rep(0, object$Nobs) 
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) *
                          pnorm(a - exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) +
                                                                exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 * exp(A) *
                                                                                                      pnorm(a) - exp(B) * pnorm(b))
    teBC_c2 <- rep(1, object$Nobs)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A) * exp(a * exp(Wv/2) +
                                                              exp(Wv)/2) * pnorm(a + exp(Wv/2)) - exp(B) * exp(b *
                                                                                                                 exp(Wv/2) + exp(Wv)/2) * pnorm(b + exp(Wv/2)))/(2 *
                                                                                                                                                                   exp(A) * pnorm(a) - exp(B) * pnorm(b))
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
                     odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
                     u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
                     u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

czisfmargtslnorm_Eu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
                             object$nvZVvar + 1]
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
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * 
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * 
                                                          exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
                                                           4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
                              nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}

czisfmargtslnorm_Vu <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
                                                object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
                                                               object$nuZUvar + object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
                             object$nvZVvar + 1]
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
  A <- exp(Wv)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + object$S * epsilon * 
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 * 
                                                          exp(A) * pnorm(a) - exp(B) * pnorm(b))
  Pi2 <- 1/exp(Wv/2) * dnorm(object$S * epsilon/exp(Wv/2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
                                                           8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
                                                                                                                          lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
                                                                                                                                                                             ncol = 1))
  margEff2 <- matrix(0, nrow = object$Nobs, ncol = object$nuZUvar - 1)
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c))
}
