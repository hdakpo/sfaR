################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Standard Stochastic Frontier Analysis                                 #
# Convolution: truncated normal (scaling) - normal                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for truncated normal (scaling)-normal distribution
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
#' @noRd
ctruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar - 1)]
  tau <- parm[nXvar + length(delta) + 1]
  cu <- parm[nXvar + length(delta) + 2]
  phi <- parm[(nXvar + length(delta) + 2 + 1):(nXvar + length(delta) + 2 + nvZVvar)]
  usca <- crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]
  musca <- exp(usca) * tau
  Wusca <- cu + 2 * usca
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- (musca * exp(Wvsca) - exp(Wusca) * S * epsilon)/(exp(Wusca) + exp(Wvsca))
  sigmastar <- sqrt(exp(Wusca) * exp(Wvsca)/(exp(Wusca) + exp(Wvsca)))
  ll <- (-1/2 * log(exp(Wusca) + exp(Wvsca)) + dnorm((musca +
    S * epsilon)/sqrt(exp(Wusca) + exp(Wvsca)), log = TRUE) +
    log(pnorm(mustar/sigmastar)) - log(pnorm(musca/sqrt(exp(Wusca)))))
  RTMB::ADREPORT(ll * wHvar)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for truncated normal (scaling)-normal distribution
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
csttruncnormscal <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  if (S * m3 > 0) {
    ## Coeli (1995) suggests 0.05 for gamma
    varu <- (abs(S * m3 * sqrt(pi/2)/(1 - 4/pi)))^(2/3)
  } else {
    varu <- (S * m3 * sqrt(pi/2)/(1 - 4/pi))^(2/3)
  }
  if (m2 < (pi - 2)/pi * varu) {
    varv <- abs(m2 - (1 - 2/pi) * varu)
  } else {
    varv <- m2 - (1 - 2/pi) * varu
  }
  dep_u <- 1/2 * log(((epsiRes^2 - varv) * pi/(pi - 2))^2)
  dep_v <- 1/2 * log((epsiRes^2 - (1 - 2/pi) * varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar,
      drop = FALSE]))
  }
  if (any(is.na(reg_hetu$coefficients)))
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(if (length(grep("Intercept", colnames(uHvar))) ==
        0)
        colnames(uHvar)[is.na(reg_hetu$coefficients)[-1]] else colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicolinearity",
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
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicolinearity",
      call. = FALSE)
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])
  } else {
    beta <- olsObj
  }
  delta <- rep(0, length(coefficients(reg_hetu)) - 1)
  cu <- unname(coefficients(reg_hetu)[1])
  names(delta) <- paste0("Zscale_", colnames(uHvar)[-1])
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  tau <- mean(epsiRes)
  return(c(beta, delta, tau = tau, cu = cu, phi))
}

# Gradient of the likelihood function ----------
#' gradient for truncated normal (scaling)-normal distribution
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
#' @noRd
cgradtruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar - 1)]
  tau <- parm[nXvar + length(delta) + 1]
  cu <- parm[nXvar + length(delta) + 2]
  phi <- parm[(nXvar + length(delta) + 2 + 1):(nXvar + length(delta) + 2 + nvZVvar)]
  usca <- crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]
  Wusca <- cu + 2 * usca
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- exp(usca)
  .e2 <- exp(Wvsca)
  .e3 <- exp(Wusca)
  .e4 <- exp((Wusca)/2)
  .e5 <- sqrt(.e3 + .e2)
  .e6 <- (S * epsilon + tau * .e1)
  .e7 <- (tau * .e1 * .e2 - S * .e3 * epsilon)
  .e8 <- (.e3 + .e2)
  .e9 <- sqrt(.e3 * .e2/.e8)
  .e10 <- (.e8 * .e9)
  .e11 <- dnorm(.e7/.e10)
  .e12 <- dnorm(.e6/.e5)
  .e13 <- pnorm(.e7/.e10)
  .e14 <- dnorm(tau * .e1/.e4)
  .e15 <- pnorm(tau * .e1/.e4)
  .e16 <- (tau * .e1 * .e2 - 2 * (S * .e3 * epsilon))
  .e17 <- (0.5 * ((2 - 2 * (.e3/.e8)) * .e2/.e9) + 2 * .e9)
  .e18 <- (.e16/.e10 - .e17 * .e3 * .e7/.e10^2)
  .e19 <- (tau * .e1 - .e3 * .e6/.e8)
  .e20 <- (.e12 * .e6 * .e19/.e12 + .e3)
  .e21 <- (.e11 * .e2/(.e13 * .e9) - .e12 * .e6/.e12)
  .e22 <- (0.5 * (.e12 * .e6^2/(.e12 * .e8)) - 0.5)
  .e23 <- (0.5 * ((1 - .e3/.e8) * .e2/.e9) + .e9)
  .e24 <- (.e23 * .e7/.e10^2 + S * epsilon/.e10)
  .e25 <- (.e22/.e8 - .e24 * .e11/.e13)
  .e26 <- (0.5 * ((1 - .e2/.e8) * .e3/.e9) + .e9)
  .e27 <- (tau * .e1/.e10 - .e26 * .e7/.e10^2)
  gradll <- cbind(S * Xvar * (.e12 * .e6/.e12 + .e11 * .e3/(.e13 * .e9))/.e8, uHvar[,
    -1, drop = FALSE] * (.e18 * .e11/.e13 - .e20/.e8), (.e21/.e8 - .e14/(.e4 *
    .e15)) * .e1, .e25 * .e3 + 0.5 * (tau * .e14 * .e1/(.e4 * .e15)), vHvar *
    (.e22/.e8 + .e11 * .e27/.e13) * .e2)
  return(gradll * wHvar)
}

# Hessian of the likelihood function ----------
#' hessian for truncated normal (scaling)-normal distribution
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
#' @noRd
chesstruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar, vHvar, Yvar,
  Xvar, S, wHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar - 1)]
  tau <- parm[nXvar + length(delta) + 1]
  cu <- parm[nXvar + length(delta) + 2]
  phi <- parm[(nXvar + length(delta) + 2 + 1):(nXvar + length(delta) + 2 + nvZVvar)]
  usca <- crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]
  Wusca <- cu + 2 * usca
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- exp(usca)
  .e2 <- exp(Wvsca)
  .e3 <- exp(Wusca)
  .e4 <- exp((Wusca)/2)
  .e5 <- sqrt(.e3 + .e2)
  .e6 <- (S * epsilon + tau * .e1)
  .e7 <- (tau * .e1 * .e2 - S * .e3 * epsilon)
  .e8 <- (.e3 + .e2)
  .e9 <- sqrt(.e3 * .e2/.e8)
  .e10 <- (.e8 * .e9)
  .e11 <- dnorm(.e7/.e10)
  .e12 <- dnorm(.e6/.e5)
  .e13 <- pnorm(.e7/.e10)
  .e14 <- dnorm(tau * .e1/.e4)
  .e15 <- pnorm(tau * .e1/.e4)
  .e16 <- (tau * .e1 * .e2 - 2 * (S * .e3 * epsilon))
  .e17 <- (0.5 * ((2 - 2 * (.e3/.e8)) * .e2/.e9) + 2 * .e9)
  .e18 <- (.e16/.e10 - .e17 * .e3 * .e7/.e10^2)
  .e19 <- (tau * .e1 - .e3 * .e6/.e8)
  .e20 <- (.e12 * .e6 * .e19/.e12 + .e3)
  .e21 <- (.e11 * .e2/(.e13 * .e9) - .e12 * .e6/.e12)
  .e22 <- (0.5 * (.e12 * .e6^2/(.e12 * .e8)) - 0.5)
  .e23 <- (0.5 * ((1 - .e3/.e8) * .e2/.e9) + .e9)
  .e24 <- (.e23 * .e7/.e10^2 + S * epsilon/.e10)
  .e25 <- (.e22/.e8 - .e24 * .e11/.e13)
  .e26 <- (0.5 * ((1 - .e2/.e8) * .e3/.e9) + .e9)
  .e27 <- (tau * .e1/.e10 - .e26 * .e7/.e10^2)
  .e28 <- ((1 - .e12/.e12) * .e6^2/.e8 - 1)
  .e29 <- (((.e6^2/.e8 - 2)/(.e12 * .e8) - .e12 * .e6^2/(.e12 * .e8)^2) * .e12 *
    .e6/.e8)
  .e30 <- (2 - 2 * (.e3/.e8))
  .e31 <- ((2 * (tau * .e1) - .e6^2 * .e19/.e8)/(.e12 * .e8) - (2 * (.e12 * .e3) -
    .e12 * .e6 * .e19) * .e6/(.e12 * .e8)^2)
  .e32 <- (0.5 * (.e31 * .e12 * .e6) - 2 * (.e22 * .e3/.e8))
  .e33 <- (.e7/.e10 + .e11/.e13)
  hessll <- matrix(0, nrow = nXvar + (nuZUvar - 1) + 1 + 1 + nvZVvar, ncol = nXvar +
    (nuZUvar - 1) + 1 + 1 + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((.e28 * .e12/.e12 - (.e7/(.e2 *
    .e13 * .e9) + .e11 * .e3/(.e13 * .e9)^2) * .e11 * .e3/.e8)/.e8) * wHvar,
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + (nuZUvar - 1))] <- crossprod(S * Xvar *
    (((2/.e10 - .e17 * .e3/.e10^2) * .e3 - .e18 * (.e7/.e2 + .e11 * .e3/(.e13 *
      .e9))/.e8) * .e11/.e13 - ((.e6^2/.e8 - 1) * .e19 + (.e3 - .e12 * .e6 *
      .e19/.e12) * .e6/.e8) * .e12/(.e12 * .e8)) * wHvar, uHvar[, -1, drop = FALSE])
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 1] <- crossprod(S * Xvar, (-((.e28 *
    .e12/.e12 + (.e7/(.e13 * .e9) + .e11 * .e3 * .e2/(.e13 * .e9)^2) * .e11/.e8) *
    .e1/.e8)) * wHvar)
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 2] <- crossprod(S * Xvar, ((0.5 * .e29 -
    ((.e23/.e10^2 - .e24 * .e11/(.e8 * .e13 * .e9)) * .e3 - (.e24 * .e7/.e2 +
      1/.e9)/.e8) * .e11/.e13) * .e3) * wHvar)
  hessll[1:nXvar, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 +
    nvZVvar)] <- crossprod(S * Xvar * ((0.5 * .e29 - ((.e26/.e10^2 + .e11 * .e27/(.e8 *
    .e13 * .e9)) * .e3 + .e7 * .e27/(.e8 * .e2)) * .e11/.e13) * .e2) * wHvar,
    vHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + 1):(nXvar + (nuZUvar - 1))] <- crossprod(uHvar[,
    -1, drop = FALSE] * (((tau * .e1 * .e2 - (.e18^2 * .e7 + 4 * (S * .e3 * epsilon)))/.e10 -
    (((.e30 * (.e3/.e8 - 0.5 * (0.5 * .e30 + 2 * (.e3/.e8))) * .e2/.e9 + 2 *
      .e17) * .e7 + .e17 * (2 * (tau * .e1 * .e2) - (2 * (.e17 * .e8 * .e9 *
      .e7/.e10^2) + 4 * (S * epsilon)) * .e3)) * .e3/.e10^2 + .e18^2 * .e11/.e13)) *
    .e11/.e13 - ((((.e12 * .e6 * .e19^2/.e12 - (.e30 * .e6 + tau * .e1) * .e3)/.e8 +
    tau * .e1) * .e6 + (tau * .e1 - .e6^2 * .e19/.e8) * .e19) * .e12/.e12 + (2 -
    2 * (.e20/.e8)) * .e3)/.e8) * wHvar, uHvar[, -1, drop = FALSE])
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar - 1) + 1] <- crossprod(uHvar[,
    -1, drop = FALSE], (((((1/(.e13 * .e9) - (.e18 * .e11 * .e9 + 0.5 * (.e30 *
    .e3 * .e2 * .e13/.e10))/(.e13 * .e9)^2) * .e2 - .e18 * .e7/(.e3 * .e13)) *
    .e11 - (((.e12/.e12 - 1) * .e6^2 * .e19/.e8 + 2 * (tau * .e1) + S * epsilon) *
    .e12/.e12 + 2 * (.e21 * .e3/.e8)))/.e8 + .e14 * (.e4 * .e15/(.e4 * .e15)^2 -
    1/(.e4 * .e15))) * .e1) * wHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar - 1) + 2] <- crossprod(uHvar[,
    -1, drop = FALSE], ((.e32/.e8 + 2 * .e25 - ((.e23 * (tau * .e1 * .e2 - (2 *
    (.e17 * .e8 * .e9 * .e7/.e10^2) + 2 * (S * epsilon)) * .e3) + (0.5 * (.e3/.e8) -
    0.5 * (0.5 * (1 - .e3/.e8) + .e3/.e8)) * .e30 * .e2 * .e7/.e9 - S * .e17 *
    .e3 * epsilon)/.e10^2 - .e24 * .e18 * .e33) * .e11/.e13) * .e3 + 0.5 * (tau *
    (1/(.e4 * .e15) - .e4 * .e15/(.e4 * .e15)^2) * .e14 * .e1)) * wHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar +
    (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(uHvar[, -1, drop = FALSE] * ((.e32/.e8 +
    .e11 * (tau * (1/.e10 - .e17 * .e3/.e10^2) * .e1 - (((0.5 * ((1 - .e2/.e8) *
      (2 - 0.5 * .e30) + 2 * (.e3 * .e2/.e8^2)) + 0.5 * (.e30 * .e2/.e8)) *
      .e3 * .e7/.e9 + .e26 * (tau * .e1 * .e2 - (2 * (.e17 * .e8 * .e9 * .e7/.e10^2) +
      2 * (S * epsilon)) * .e3))/.e10^2 + .e18 * .e33 * .e27))/.e13) * .e2) *
    wHvar, vHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 1] <- sum(((.e14 *
    (.e14/(.e4 * .e15)^2 + tau * .e1/(.e4^3 * .e15)) - (((.e12/.e12 - 1) * .e6^2/.e8 +
    1) * .e12/.e12 + (.e7/(.e3 * .e13 * .e9) + .e11 * .e2/(.e13 * .e9)^2) * .e11 *
    .e2/.e8)/.e8) * .e1^2) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 2] <- sum((((0.5 *
    (((2 - .e6^2/.e8)/(.e12 * .e8) + .e12 * .e6^2/(.e12 * .e8)^2) * .e12 * .e6/.e8) -
    (.e23 * .e2/.e10^2 - .e24 * (.e7/.e3 + .e11 * .e2/(.e13 * .e9))/.e8) * .e11/.e13) *
    .e3 + 0.5 * (((1 - tau^2 * .e1^2/.e4^2)/(.e4 * .e15) - tau * .e14 * .e1/(.e4 *
    .e15)^2) * .e14)) * .e1) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)] <- crossprod(vHvar, ((((1/.e9 - .e7 * .e27/.e3)/.e8 -
    (.e26/.e10^2 + .e11 * .e27/(.e8 * .e13 * .e9)) * .e2) * .e11/.e13 + 0.5 *
    (((2 - .e6^2/.e8)/(.e12 * .e8) + .e12 * .e6^2/(.e12 * .e8)^2) * .e12 * .e6/.e8)) *
    .e1 * .e2) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 2, nXvar + (nuZUvar - 1) + 2] <- sum(((((0.5 *
    ((0.5 * (.e6^2/(.e12 * .e8^3)) - (0.5 * (.e12 * .e6^2/.e8) + .e12)/(.e12 *
      .e8)^2) * .e12 * .e6^2) - .e22/.e8) * .e3 + 0.5 * (.e12 * .e6^2/(.e12 *
    .e8)) - 0.5)/.e8 - ((.e23 * (tau * .e1 * .e2 - (2 * (.e23 * .e8 * .e9 * .e7/.e10^2) +
    3 * (S * epsilon)) * .e3) + (0.5 * (.e3/.e8) - 0.5 * (0.5 * (1 - .e3/.e8) +
    .e3/.e8)) * (1 - .e3/.e8) * .e2 * .e7/.e9)/.e10^2 + .e24^2 * .e33 * .e3 +
    S * epsilon/.e10) * .e11/.e13) * .e3 + 0.5 * (tau * (0.5 * (tau^2 * .e1^2/(.e4^3 *
    .e15)) - (0.5 * (.e4 * .e15) - 0.5 * (tau * .e14 * .e1))/(.e4 * .e15)^2) *
    .e14 * .e1)) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 2, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)] <- crossprod(vHvar, (((.e24 * .e33 * .e27 - ((0.5 * ((1 -
    .e3/.e8) * .e2/.e8) + 0.5 * ((.e3/.e8 - 1) * .e2/.e8 + 1 - 0.5 * ((1 - .e3/.e8) *
    (1 - .e2/.e8)))) * .e7/.e9 + tau * .e23 * .e1 - .e26 * (2 * (.e23 * .e8 *
    .e9 * .e7/.e10^2) + S * epsilon))/.e10^2) * .e11/.e13 + (0.5 * ((0.5 * (.e6^2/(.e12 *
    .e8^3)) - (0.5 * (.e12 * .e6^2/.e8) + .e12)/(.e12 * .e8)^2) * .e12 * .e6^2) -
    .e22/.e8)/.e8) * .e3 * .e2) * wHvar)
  hessll[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar),
    (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(vHvar *
    ((((0.5 * ((0.5 * (.e6^2/(.e12 * .e8^3)) - (0.5 * (.e12 * .e6^2/.e8) + .e12)/(.e12 *
      .e8)^2) * .e12 * .e6^2) - .e22/.e8) * .e2 + 0.5 * (.e12 * .e6^2/(.e12 *
      .e8)) - 0.5)/.e8 + .e11 * (tau * .e1/.e10 - ((((3 * (tau * .e1) - 2 *
      (.e26 * .e8 * .e9 * .e7/.e10^2)) * .e2 - S * .e3 * epsilon) * .e26 +
      (0.5 * (.e2/.e8) - 0.5 * (0.5 * (1 - .e2/.e8) + .e2/.e8)) * (1 - .e2/.e8) *
        .e3 * .e7/.e9)/.e10^2 + .e33 * .e2 * .e27^2))/.e13) * .e2) * wHvar,
    vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for truncated normal (scaling)-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
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
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
truncnormscalAlgOpt <- function(start, randStart, sdStart, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, wHvar, method, derivs, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac, accuracy, stepsize) {
  ## starting values and log likelihood ------
  startVal <- if (!is.null(start))
    start else csttruncnormscal(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], S = S,
    uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar)
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ctruncnormscalike(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S))
  ## automatic differentiation ------
  if (derivs == "ad") {
    if (method %in% c("bfgs", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
    cgradbhhhtruncnormscalAD <- RTMB::MakeADFun(function(p) sum(ctruncnormscalike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = TRUE, silent = TRUE)
    cgradhesstruncnormscalAD <- RTMB::MakeADFun(function(p) sum(ctruncnormscalike(p, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar)), startVal, ADreport = FALSE, silent = TRUE)
    fnGradObstruncnormscal <- function(parm) {
      Ta1 <- RTMB::GetTape(cgradbhhhtruncnormscalAD)
      Ta2 <- RTMB::MakeTape(function(weight) {
        WT <- RTMB::MakeTape(function(x) sum(Ta1(x) * weight), parm)
        (WT$jacfun())(RTMB::advector(parm))
      }, rep(1, nrow(Xvar)))
      t((Ta2$jacfun())(RTMB::advector(rep(1, nrow(Xvar)))))
    }
    ### solve for different algorithms ------
    mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(p) -cgradhesstruncnormscalAD$fn(p),
      gr = function(p) -cgradhesstruncnormscalAD$gr(p), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
        maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
      maxLikAlgo = maxRoutine(fn = cgradhesstruncnormscalAD$fn, grad = cgradhesstruncnormscalAD$gr,
        hess = cgradhesstruncnormscalAD$he, start = startVal, finalHessian = TRUE,
        control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax,
          reltol = tol, tol = tol, qac = qac)), bhhh = maxLik::maxBHHH(fn = cgradbhhhtruncnormscalAD$fn,
        grad = fnGradObstruncnormscal, hess = cgradhesstruncnormscalAD$he, start = startVal,
        finalHessian = TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac)), sr1 = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstruncnormscalAD$fn(p), gr = function(p) -cgradhesstruncnormscalAD$gr(p),
        method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
        fn = function(p) -cgradhesstruncnormscalAD$fn(p), gr = function(p) -cgradhesstruncnormscalAD$gr(p),
        hs = function(p) as(-cgradhesstruncnormscalAD$he(p), "dgCMatrix"), method = "Sparse",
        control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(p) -cgradhesstruncnormscalAD$fn(p),
        gr = function(p) -cgradhesstruncnormscalAD$gr(p), hess = function(p) -cgradhesstruncnormscalAD$he(p),
        print.info = printInfo, maxiter = itermax, epsa = gradtol, epsb = gradtol),
      nlminb = nlminb(start = startVal, objective = function(p) -cgradhesstruncnormscalAD$fn(p),
        gradient = function(p) -cgradhesstruncnormscalAD$gr(p), hessian = function(p) -cgradhesstruncnormscalAD$he(p),
        control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
          eval.max = itermax, rel.tol = tol, x.tol = tol)))
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$gradient <- cgradhesstruncnormscalAD$gr(mleObj$par)
    }
    mlParam <- if (method %in% c("ucminf", "nlminb")) {
      mleObj$par
    } else {
      if (method %in% c("maxLikAlgo", "bhhh")) {
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
        mleObj$hessian <- cgradhesstruncnormscalAD$he(mleObj$par)
      if (method == "sr1")
        mleObj$hessian <- cgradhesstruncnormscalAD$he(mleObj$solution)
      if (method == "mla")
        mleObj$hessian <- cgradhesstruncnormscalAD$he(mleObj$b)
    }
    mleObj$logL_OBS <- cgradbhhhtruncnormscalAD$fn(mlParam)
    mleObj$gradL_OBS <- fnGradObstruncnormscal(mlParam)
    rm(cgradbhhhtruncnormscalAD, cgradhesstruncnormscalAD, fnGradObstruncnormscal)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
    }
    ## numerical derivatives ------
    if (derivs == "numerical") {
      cgradtruncnormscalikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::jacobian(ctruncnormscalike, var = unname(parm), params = list(nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar), accuracy = accuracy,
          stepsize = stepsize)
      }
      chesstruncnormscalikeNum <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
        vHvar, Yvar, Xvar, wHvar, S) {
        calculus::hessian(function(parm) sum(ctruncnormscalike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar)), var = unname(parm),
          accuracy = accuracy, stepsize = stepsize)
      }
      ### solve for different algorithms ------
      mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ctruncnormscalike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
        gr = function(parm) -colSums(cgradtruncnormscalikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
          maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)),
        maxLikAlgo = maxRoutine(fn = ctruncnormscalike, grad = cgradtruncnormscalikeNum,
          hess = chesstruncnormscalikeNum, start = startVal, finalHessian = if (hessianType ==
          2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), method = "SR1", control = list(maxit = itermax,
          cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hs = function(parm) as(-chesstruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"), method = "Sparse",
          control = list(maxit = itermax, cgtol = gradtol, stop.trust.radius = tol,
          prec = tol, report.level = if (printInfo) 2 else 0, report.precision = 1L,
          preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ctruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          gr = function(parm) -colSums(cgradtruncnormscalikeNum(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hess = function(parm) -chesstruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo,
          maxiter = itermax, epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradtruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), hessian = function(parm) -chesstruncnormscalikeNum(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
      if (method %in% c("ucminf", "nlminb")) {
        mleObj$gradient <- colSums(cgradtruncnormscalikeNum(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S))
      }
      mlParam <- if (method %in% c("ucminf", "nlminb")) {
        mleObj$par
      } else {
        if (method %in% c("maxLikAlgo", "bhhh")) {
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
          mleObj$hessian <- chesstruncnormscalikeNum(parm = mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        if (method == "sr1")
          mleObj$hessian <- chesstruncnormscalikeNum(parm = mleObj$solution,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
      mleObj$logL_OBS <- ctruncnormscalike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      mleObj$gradL_OBS <- cgradtruncnormscalikeNum(parm = mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
    } else {
      ## analytical derivatives ------
      if (derivs == "analytical") {
        ### solve for different algorithms ------
        mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
          stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ctruncnormscalike,
          grad = cgradtruncnormscalike, hess = chesstruncnormscalike, start = startVal,
          finalHessian = if (hessianType == 2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
          iterlim = itermax, reltol = tol, tol = tol, qac = qac), nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), sr1 = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hs = function(parm) as(-chesstruncnormscalike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), "dgCMatrix"),
          method = "Sparse", control = list(maxit = itermax, cgtol = gradtol,
          stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
          report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
          fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hess = function(parm) -chesstruncnormscalike(parm, nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Xvar = Xvar, wHvar = wHvar, S = S), print.info = printInfo, maxiter = itermax,
          epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
          objective = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)), gradient = function(parm) -colSums(cgradtruncnormscalike(parm,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)),
          hessian = function(parm) -chesstruncnormscalike(parm, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S), control = list(iter.max = itermax,
          trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
          x.tol = tol)))
        if (method %in% c("ucminf", "nlminb")) {
          mleObj$gradient <- colSums(cgradtruncnormscalike(mleObj$par, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S))
        }
        mlParam <- if (method %in% c("ucminf", "nlminb")) {
          mleObj$par
        } else {
          if (method %in% c("maxLikAlgo", "bhhh")) {
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
          mleObj$hessian <- chesstruncnormscalike(parm = mleObj$par, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
          if (method == "sr1")
          mleObj$hessian <- chesstruncnormscalike(parm = mleObj$solution, nXvar = nXvar,
            nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        }
        mleObj$logL_OBS <- ctruncnormscalike(mlParam, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
        mleObj$gradL_OBS <- cgradtruncnormscalike(parm = mlParam, nXvar = nXvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
          Yvar = Yvar, Xvar = Xvar, wHvar = wHvar, S = S)
      }
    }
  }
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for truncated normal (scaling)-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
ctruncnormscaleff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  phi <- object$mlParam[(object$nXvar + (object$nuZUvar - 1) +
    2 + 1):(object$nXvar + (object$nuZUvar - 1) + 2 + object$nvZVvar)]
  musca <- exp(crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]) * tau
  Wusca <- cu + 2 * crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    crossprod(matrix(beta), t(Xvar))[1, ]
  mustar <- (musca * exp(Wvsca) - exp(Wusca) * object$S * epsilon)/(exp(Wusca) +
    exp(Wvsca))
  sigmastar <- sqrt(exp(Wusca) * exp(Wvsca)/(exp(Wusca) + exp(Wvsca)))
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
    res <- data.frame(u = u, uLB = uLB, uUB = uUB, teJLMS = teJLMS,
      m = m, teMO = teMO, teBC = teBC, teBCLB = teBCLB,
      teBCUB = teBCUB, teBC_reciprocal = teBC_reciprocal)
  } else {
    res <- data.frame(u = u, uLB = uLB, uUB = uUB, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for truncated normal (scaling)-normal distribution
#' @param object object of class sfacross
#' @noRd
cmargtruncnormscal_Eu <- function(object) {
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  hi <- exp(crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ])
  Lambda <- tau/exp(cu/2)
  m1 <- exp(cu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  margEff <- kronecker(matrix(delta, nrow = 1), matrix(m1 *
    hi, ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}

cmargtruncnormscal_Vu <- function(object) {
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  hi2 <- exp(2 * crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ])
  Lambda <- tau/exp(cu/2)
  m2 <- exp(cu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) -
    (dnorm(Lambda)/pnorm(Lambda))^2)
  margEff <- kronecker(matrix(2 * delta, nrow = 1), matrix(m2 *
    hi2))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(data.frame(margEff))
}
