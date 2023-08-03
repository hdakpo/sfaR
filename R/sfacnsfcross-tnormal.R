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
# Link functions: - logit exp(theta * Z)/(1 + exp(theta * Z))                  #
#                 - cauchit 1/pi * atan(theta * Z) + 1/2                       #
#                 - probit pnorm(theta * Z)                                    #
#                 - cloglog 1 - exp(-exp(theta * Z))                           #
# Convolution: truncated normal - normal                                       #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf truncatednormal-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u

## logit specification class membership
ccnsftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
ccnsftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
ccnsftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
ccnsftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv1))) *
    pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + S * epsilon)/sqrt(exp(Wu) + exp(Wv2))) *
    pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# different sigma_u

## logit specification class membership
cmcesftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * S * epsilon)/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmcesftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * S * epsilon)/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmcesftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar, muHvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * S * epsilon)/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmcesftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * S * epsilon)/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(-.Machine$double.xmax), return(wHvar * log(Probc1 *
    Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf truncatednormal-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
# same sigma_u
cstcnsftruncnorm <- function(olsObj, epsiRes, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter,
  initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csttruncnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, uHvar = uHvar[,
      1, drop = FALSE], nuZUvar = 1, vHvar = vHvar[, 1, drop = FALSE], nvZVvar = 1,
      nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE])
    initTrunc <- NULL
  } else {
    cat("Initialization: SFA + truncated normal - normal distribution...\n")
    initTrunc <- maxLik::maxLik(logLik = ctruncnormlike, start = csttruncnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, uHvar = uHvar[, 1, drop = FALSE], nuZUvar = 1,
      vHvar = vHvar[, 1, drop = FALSE], nvZVvar = 1, nmuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE]), grad = cgradtruncnormlike, method = initAlg, control = list(iterlim = initIter,
      printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar, nuZUvar = 1,
      nvZVvar = 1, uHvar = uHvar[, 1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE],
      Yvar = Yvar, Xvar = Xvar, S = S, nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE],
      wHvar = wHvar)
    Esti <- initTrunc$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nmuZUvar > 1) rep(0, nmuZUvar -
    1), Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar +
    3], if (nvZVvar > 1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 3], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_", colnames(muHvar)),
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
      colnames(vHvar)), paste0("CNSF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initTrunc = initTrunc))
}

# different sigma_u
cstmcesftruncnorm <- function(olsObj, epsiRes, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart, initIter,
  initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csttruncnorm(olsObj = olsObj, epsiRes = epsiRes, S = S, uHvar = uHvar[,
      1, drop = FALSE], nuZUvar = 1, vHvar = vHvar[, 1, drop = FALSE], nvZVvar = 1,
      nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE])
    initTrunc <- NULL
  } else {
    cat("Initialization: SFA + truncated normal - normal distribution...\n")
    initTrunc <- maxLik::maxLik(logLik = ctruncnormlike, start = csttruncnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, uHvar = uHvar[, 1, drop = FALSE], nuZUvar = 1,
      vHvar = vHvar[, 1, drop = FALSE], nvZVvar = 1, nmuZUvar = 1, muHvar = muHvar[,
        1, drop = FALSE]), grad = cgradtruncnormlike, method = initAlg, control = list(iterlim = initIter,
      printLevel = if (printInfo) 2 else 0, reltol = tol), nXvar = nXvar, nuZUvar = 1,
      nvZVvar = 1, uHvar = uHvar[, 1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE],
      Yvar = Yvar, Xvar = Xvar, S = S, nmuZUvar = 1, muHvar = muHvar[, 1, drop = FALSE],
      wHvar = wHvar)
    Esti <- initTrunc$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nmuZUvar > 1) rep(0,
    nmuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nmuZUvar > 1) rep(0, nmuZUvar -
    1), 0.95 * Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar - 1), 1.05 *
    Esti[nXvar + 2], if (nuZUvar > 1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar +
    3], if (nvZVvar > 1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 3], if (nvZVvar >
    1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zmu_", colnames(muHvar)),
    paste0("Zmu_", colnames(muHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zu_",
      colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
    paste0("CNSF_", colnames(Zvar)))
  return(list(StartVal = StartVal, initTrunc = initTrunc))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf truncatednormal-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u

## logit specification class membership
cgradcnsftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  starsq1 <- ((sigma_sq1) * sigmastar1)
  starsq2 <- ((sigma_sq2) * sigmastar2)
  ssq1 <- starsq1^2
  ssq2 <- starsq2^2
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  musig1 <- musi1/starsq1
  musig2 <- musi2/starsq2
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  epsi1 <- mupsi/sqrt(sigma_sq1)
  epsi2 <- mupsi/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  sqx <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  wzsrq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (prC * sigx1_2/sqx + sigx1_1 * ewz/wzsq1)
  sigx3 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * ewz * pmusig1/wzsrq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (prC * sigx4_2/sqx + sigx4_1 * ewz/wzsq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/(sigma_sq1)^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/(sigma_sq2)^2)
  usq1 <- (1 - ewu/(sigma_sq1))
  usq2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/ssq1 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/ssq2 + S * (epsilon)/starsq2)
  sigx9_1 <- (wzdeno * depsi1 * pmusig1/wzsrq1^2)
  sigx9_2 <- (pmusig2/(sigma_sq2))
  sigx10_1 <- ((0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)/wzdeno - 0.5 * sigx9_1)
  sigx10_2 <- (0.5 * sigx6_2 - (sigx8_2 * dmusig2 + 0.5 * sigx9_2) * depsi2)
  sigx11 <- (sigx10_1 * ewz/sqrt(sigma_sq1) + sigx10_2 * prC/sqrt(sigma_sq2))
  vsq1 <- (1 - ewv1/(sigma_sq1))
  vsq2 <- (1 - ewv2/(sigma_sq2))
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/ssq1)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/ssq2)
  sigx14_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx14_2 <- (dmusig2 * sigx13_2 - 0.5 * sigx9_2)
  sigx15_1 <- (sigx14_1/wzdeno - 0.5 * sigx9_1)
  sigx15_2 <- (sigx14_2 * depsi2 + 0.5 * sigx6_2)
  sigx16 <- (1/wzsrq1 - ewz * sqrt(sigma_sq1)/wzsrq1^2)
  sigx17 <- (sigx16 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/wzsq2)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  s3sq2 <- (sigx3 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = (sigx5/sigx3 - dmu/(ewu_h * pmu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx11 * ewu/sigx3 + 0.5 * (mu * dmu/(ewu_h *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx15_1 * ewv1 *
      ewz/s3sq1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx15_2 * prC *
      ewv2/s3sq2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx17 * ewz/sigx3,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/sqrt(sigma_sq2) + ewz1 * depsi1 * pmusig1/sqrt(sigma_sq1))
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (ewz2 * sigx4_2/ssqq2 + ewz1 * sigx4_1/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (ewz2 * sigx11_2/sqrt(sigma_sq2) + sigx11_1 * ewz1/sqrt(sigma_sq1))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = (sigx5/sigx2 - dmu/(ewusr * pmu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * ewz1 *
      ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      ewz2 * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9/(pi *
      sigx2 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * pmusig1 * pwZ/sqrt(sigma_sq1))
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - pwZ) * sigx4_2/ssqq2 + sigx4_1 * pwZ/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * pwZ/sqrt(sigma_sq1) + sigx11_2 * (1 - pwZ)/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = (sigx5/sigx2 - dmu/(ewusr * pmu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * pwZ *
      ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      (1 - pwZ) * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9 *
      dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradcnsftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/sqrt(sigma_sq1) + depsi2 * prZ * pmusig2/sqrt(sigma_sq2))
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - prZ) * sigx4_1/ssqq1 + sigx4_2 * prZ/ssqq2)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx11_2 * prZ/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = (sigx5/sigx2 - dmu/(ewusr * pmu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = (sigx14 * ewu/sigx2 + 0.5 * (mu * dmu/(ewusr *
      pmu))), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * (1 -
      prZ) * ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      prZ * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9 *
      prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  muvu1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  muvu2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  musig_sq1 <- muvu1/ssq1^2
  musig_sq2 <- muvu2/ssq2^2
  musig1 <- muvu1/ssq1
  musig2 <- muvu2/ssq2
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmu1 <- pnorm(mu1/ewu1_h)
  pmu2 <- pnorm(mu2/ewu2_h)
  dmu1 <- dnorm(mu1/ewu1_h, 0, 1)
  dmu2 <- dnorm(mu2/ewu2_h, 0, 1)
  musi1 <- (mu1 + S * (epsilon))
  musi2 <- (mu2 + S * (epsilon))
  epsi1 <- musi1/sqrt(sigma_sq1)
  epsi2 <- musi2/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * musi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * musi2 * pmusig2)
  spmu1 <- (wzdeno * (sigma_sq1) * pmu1 * sqrt(sigma_sq1))
  spmu2 <- ((sigma_sq2) * pmu2 * sqrt(sigma_sq2))
  sigx2 <- (prC * sigx1_2/spmu2 + sigx1_1 * ewz/spmu1)
  pmusq1 <- (wzdeno * pmu1 * sqrt(sigma_sq1))
  pmusq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx3 <- (prC * depsi2 * pmusig2/pmusq2 + depsi1 * ewz * pmusig1/pmusq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * musi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * musi2 * pmusig2)
  wup1 <- (pmusq1^2 * ewu1_h)
  wup2 <- (ewu2_h * pmusq2^2)
  sigx5_1 <- (sigx4_1/spmu1 - wzdeno * depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/wup1)
  sigx5_2 <- (sigx4_2/spmu2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/wup2)
  sigx6_1 <- (depsi1 * musi1^2 * pmusig1/(sigma_sq1)^2)
  sigx6_2 <- (depsi2 * musi2^2 * pmusig2/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  pvs1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  pvs2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx7_1 <- (pvs1 * musig_sq1 + S * (epsilon)/ssq1)
  sigx7_2 <- (pvs2 * musig_sq2 + S * (epsilon)/ssq2)
  sigx8_1 <- (0.5 * sigx6_1 - sigx7_1 * dmusig1 * depsi1)
  sigx8_2 <- (0.5 * sigx6_2 - sigx7_2 * dmusig2 * depsi2)
  sigx9_1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  sigx9_2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx10_1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewu1_h)
  sigx10_2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewu2_h)
  sigx11_1 <- (sigx8_1 * ewu1/pmusq1 - (0.5 * sigx9_1 - 0.5 * sigx10_1) * wzdeno *
    depsi1 * pmusig1/pmusq1^2)
  sigx11_2 <- (sigx8_2 * ewu2/pmusq2 - (0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2 *
    pmusig2/pmusq2^2)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx12_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/ssq1 - sigx12_1 * musig_sq1)
  sigx13_2 <- (mu2/ssq2 - sigx12_2 * musig_sq2)
  sigx14_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx14_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx15_1 <- (wzdeno * depsi1 * pmusig1 * pmu1/pmusq1^2)
  sigx15_2 <- (depsi2 * pmusig2 * pmu2/pmusq2^2)
  sigx16_1 <- (sigx14_1/(wzdeno * pmu1) - 0.5 * sigx15_1)
  sigx16_2 <- (sigx14_2/pmu2 - 0.5 * sigx15_2)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  s3sq2 <- (sigx3 * sqrt(sigma_sq2))
  sigx17 <- (1/pmusq1 - ewz * pmu1 * sqrt(sigma_sq1)/pmusq1^2)
  sigx18 <- (wzdeno * pmu2 * sqrt(sigma_sq2))
  sigx19 <- (sigx17 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/sigx18)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx2/sigx3, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = sigx5_1 * ewz/sigx3, FUN = "*"), sweep(muHvar,
      MARGIN = 1, STATS = sigx5_2 * prC/sigx3, FUN = "*"), sweep(uHvar, MARGIN = 1,
      STATS = sigx11_1 * ewz/sigx3, FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx11_2 *
      prC/sigx3, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * ewv1 *
      ewz/s3sq1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 * prC *
      ewv2/s3sq2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx19 * ewz/sigx3,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = sigx5_1 * ewz1/sigx2, FUN = "*"), sweep(muHvar,
      MARGIN = 1, STATS = sigx5_2 * ewz2/sigx2, FUN = "*"), sweep(uHvar, MARGIN = 1,
      STATS = sigx11_1 * ewz1/sigx2, FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx11_2 *
      ewz2/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_1 * ewz1 *
      ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = sigx16_2 *
      ewz2 * ewv2/sigx17_2, FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = sigx9/sigx18,
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = sigx5_1 * pwZ/sigx2, FUN = "*"), sweep(muHvar,
      MARGIN = 1, STATS = sigx5_2 * (1 - pwZ)/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_1 * pwZ/sigx2, FUN = "*"), sweep(uHvar, MARGIN = 1,
      STATS = sigx11_2 * (1 - pwZ)/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx16_1 * pwZ * ewv1/sigx17_1, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx16_2 * (1 - pwZ) * ewv2/sigx17_2, FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx9 * dwZ/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx3/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = sigx5_1 * (1 - prZ)/sigx2, FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = sigx5_2 * prZ/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_1 * (1 - prZ)/sigx2, FUN = "*"), sweep(uHvar,
      MARGIN = 1, STATS = sigx11_2 * prZ/sigx2, FUN = "*"), sweep(vHvar, MARGIN = 1,
      STATS = sigx16_1 * (1 - prZ) * ewv1/sigx17_1, FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = sigx16_2 * prZ * ewv2/sigx17_2, FUN = "*"), sweep(Zvar,
      MARGIN = 1, STATS = sigx9 * prZ * ewz/sigx2, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf truncatednormal-normal distribution
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
#' @param Zvar matrix of separating variables
#' @param nZHvar number of separating variables
#' @noRd
# same sigma_u

## logit specification class membership
chesscnsftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  sigma_sq1 <- ewu + ewv1
  sigma_sq2 <- ewu + ewv2
  sigmastar1 <- sqrt(ewu * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu * ewv2/(sigma_sq2))
  starsq1 <- ((sigma_sq1) * sigmastar1)
  starsq2 <- ((sigma_sq2) * sigmastar2)
  ssq1 <- starsq1^2
  ssq2 <- starsq2^2
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  musig1 <- musi1/starsq1
  musig2 <- musi2/starsq2
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  epsi1 <- mupsi/sqrt(sigma_sq1)
  epsi2 <- mupsi/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  dmu <- dnorm(mu/ewu_h, 0, 1)
  pmu <- pnorm(mu/ewu_h)
  sqx <- ((sigma_sq2) * sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  wzsq1 <- (wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  wzsq2 <- (wzdeno * sqrt(sigma_sq2))
  wzsrq1 <- (wzdeno * sqrt(sigma_sq1))
  sigx2 <- (prC * sigx1_2/sqx + sigx1_1 * ewz/wzsq1)
  sigx3 <- (prC * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * ewz * pmusig1/wzsrq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (prC * sigx4_2/sqx + sigx4_1 * ewz/wzsq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/(sigma_sq1)^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/(sigma_sq2)^2)
  usq1 <- (1 - ewu/(sigma_sq1))
  usq2 <- (1 - ewu/(sigma_sq2))
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/ssq1 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/ssq2 + S * (epsilon)/starsq2)
  sigx9_1 <- (wzdeno * depsi1 * pmusig1/wzsrq1^2)
  sigx9_2 <- (pmusig2/(sigma_sq2))
  sigx10_1 <- ((0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)/wzdeno - 0.5 * sigx9_1)
  sigx10_2 <- (0.5 * sigx6_2 - (sigx8_2 * dmusig2 + 0.5 * sigx9_2) * depsi2)
  sigx11 <- (sigx10_1 * ewz/sqrt(sigma_sq1) + sigx10_2 * prC/sqrt(sigma_sq2))
  vsq1 <- (1 - ewv1/(sigma_sq1))
  vsq2 <- (1 - ewv2/(sigma_sq2))
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/ssq1)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/ssq2)
  sigx14_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx14_2 <- (dmusig2 * sigx13_2 - 0.5 * sigx9_2)
  sigx15_1 <- (sigx14_1/wzdeno - 0.5 * sigx9_1)
  sigx15_2 <- (sigx14_2 * depsi2 + 0.5 * sigx6_2)
  sigx16 <- (1/wzsrq1 - ewz * sqrt(sigma_sq1)/wzsrq1^2)
  sigx17 <- (sigx16 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/wzsq2)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  s3sq2 <- (sigx3 * sqrt(sigma_sq2))
  sigx18_1 <- ((2 - mupsi^2/(sigma_sq1)) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1)
  sigx18_2 <- ((2 - mupsi^2/(sigma_sq2)) * pmusig2 + dmusig2 * ewv2 * mupsi/starsq2)
  sigx19_1 <- ((mupsi^2/(sigma_sq1) - 2) * pmusig1 + dmusig1 * ewu * mupsi/starsq1)
  sigx19_2 <- ((mupsi^2/(sigma_sq2) - 2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2)
  sigx20_1 <- ((mupsi^2/(sigma_sq1) - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1)
  sigx20_2 <- ((mupsi^2/(sigma_sq2) - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2)
  sigx21_1 <- (depsi1 * mupsi - depsi1 * musi1/ewv1)
  sigx21_2 <- (depsi2 * mupsi - depsi2 * musi2/ewv2)
  sigx22_1 <- (sigx19_1 * depsi1 * mupsi/(sigma_sq1)^2)
  sigx22_2 <- (0.5 * sigx19_2 - 0.5 * pmusig2)/(sigma_sq2)
  sigx23_1 <- sigx8_1 * depsi1 * mupsi/(sigma_sq1)
  sigx23_2 <- ((sigma_sq2)^2 * sigmastar2)
  sigx24_1 <- (wzsrq1^2 * (sigma_sq1))
  sigx24_2 <- (wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigx25_1 <- ((0.5 * (mupsi^2/(sigma_sq1)) - 2) * pmusig1/(sigma_sq1) - sigx8_1 *
    dmusig1)
  sigx25_2 <- ((0.5 * (mupsi^2/(sigma_sq2)) - 2) * pmusig2/(sigma_sq2) - sigx8_2 *
    dmusig2)
  sigx26_1 <- (depsi1 * musi1/ewu + depsi1 * mupsi)
  sigx26_2 <- (sigx12_2/ssq2 + 0.5/sigx23_2)
  sigx27_1 <- (wzdeno * sigx4_1/sigx24_1)
  sigx27_2 <- (sigx7_2/ssq2 + 0.5/sigx23_2)
  sigx28_1 <- (sigx18_1 * depsi1 * mupsi/(sigma_sq1)^2)
  sigx28_2 <- (0.5 * sigx18_2 + 0.5 * pmusig2)/(sigma_sq2)
  sigx29_1 <- (sigx8_1 * dmusig1 + wzdeno^2 * pmusig1/wzsrq1^2)
  sigx29_2 <- (sigx8_2 * dmusig2 + pmusig2/(sigma_sq2))
  sigx30_1 <- (sigx7_1 * (sigma_sq1) * musi1 * sigmastar1/ssq1)
  sigx30_2 <- (sigx7_2 * (sigma_sq2) * musi2 * sigmastar2/ssq2)
  sigx31_1 <- ((0.5 * sigx6_1 - sigx29_1 * depsi1) * wzdeno/wzsrq1^2)
  sigx31_2 <- (wzdeno * depsi2 * pmusig2/wzsq2^2)
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((sigx20_1 * depsi1 + dmusig1 * sigx21_1 * ewu/starsq1) * ewz/wzsq1 +
    (sigx20_2 * depsi2 + dmusig2 * sigx21_2 * ewu/starsq2) * prC/sqx - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (prC * (dmusig2 * sigx21_2 * ewv2/starsq2 - sigx20_2 *
      depsi2)/sqx + (dmusig1 * sigx21_1 * ewv1/starsq1 - sigx20_1 * depsi1) *
      ewz/wzsq1 - sigx2 * sigx5/sigx3)/sigx3, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((sigx22_2 - sigx8_2 * dmusig2) * depsi2 *
      mupsi/(sigma_sq2) - (sigx27_2 * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/(sigma_sq2)) *
      dmusig2 * depsi2) * prC/sqrt(sigma_sq2) + ((0.5 * sigx22_1 - (sigx23_1 +
      (sigx7_1 * ewu/ssq1 - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/(sigma_sq1)) *
        depsi1) * dmusig1)/wzdeno - 0.5 * (wzdeno * sigx1_1/sigx24_1)) *
      ewz/sqrt(sigma_sq1) - sigx11 * sigx2/sigx3) * ewu/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((((sigx21_1 *
    sigx13_1/(sigma_sq1) - sigx12_1 * depsi1 * ewu/ssq1) * dmusig1 + 0.5 * sigx22_1)/wzdeno -
    0.5 * (wzdeno * sigx1_1/sigx24_1))/s3sq1 - sigx15_1 * sigx2 * sqrt(sigma_sq1)/s3sq1^2) *
    ewv1 * ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((sigx22_2 + dmusig2 * sigx13_2) * depsi2 * mupsi/(sigma_sq2) - (sigx26_2 *
    ewu + musi2 * sigx13_2/((sigma_sq2) * ewv2)) * dmusig2 * depsi2)/s3sq2 -
    sigx2 * sigx15_2 * sqrt(sigma_sq2)/s3sq2^2) * prC * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx16 * sigx1_1/(sigma_sq1) - (sigx2 * sigx17/sigx3 + prC * sigx1_2/sigx24_2)) *
    ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (dmu * (dmu/(ewu_h * pmu)^2 + mu/(ewu_h^3 * pmu)) -
      ((((1 - mupsi^2/(sigma_sq1)) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) *
        depsi1 + dmusig1 * sigx26_1 * ewv1/starsq1) * ewz/wzsq1 + (((1 -
        mupsi^2/(sigma_sq2)) * pmusig2 + dmusig2 * ewv2 * mupsi/starsq2) *
        depsi2 + dmusig2 * (depsi2 * musi2/ewu + depsi2 * mupsi) * ewv2/starsq2) *
        prC/sqx + sigx5^2/sigx3)/sigx3), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_2 *
    dmusig2 + sigx28_2) * depsi2 * mupsi/(sigma_sq2) - (sigx27_2 * ewv2 - sigx8_2 *
    musi2/(ewu * (sigma_sq2))) * dmusig2 * depsi2) * prC/sqrt(sigma_sq2) + ((0.5 *
    sigx28_1 - ((sigx7_1 * ewv1/ssq1 - sigx8_1 * musi1/(ewu * (sigma_sq1))) *
    depsi1 - sigx23_1) * dmusig1)/wzdeno - 0.5 * sigx27_1) * ewz/sqrt(sigma_sq1) -
    sigx11 * sigx5/sigx3) * ewu/sigx3 + 0.5 * (((1 - mu^2/ewu_h^2)/(ewu_h * pmu) -
    mu * dmu/(ewu_h * pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (((((1/starsq1 - sigx12_1 * ewv1/ssq1) * depsi1 - sigx26_1 * sigx13_1/(sigma_sq1)) *
      dmusig1 + 0.5 * sigx28_1)/wzdeno - 0.5 * sigx27_1)/s3sq1 - sigx15_1 *
      sigx5 * sqrt(sigma_sq1)/s3sq1^2) * ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((sigx28_2 - dmusig2 * sigx13_2) * depsi2 *
      mupsi/(sigma_sq2) + ((1/sigmastar2 - musi2 * sigx13_2/ewu)/(sigma_sq2) -
      sigx26_2 * ewv2) * dmusig2 * depsi2)/s3sq2 - sigx5 * sigx15_2 * sqrt(sigma_sq2)/s3sq2^2) *
      prC * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (sigx16 * sigx4_1/(sigma_sq1) - (sigx5 * sigx17/sigx3 +
      prC * sigx4_2/sigx24_2)) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((0.5 * sigx25_2 - 0.5 * (sigx8_2 * dmusig2 + 0.5 * sigx9_2)) * depsi2 *
      mupsi^2/(sigma_sq2) - 0.5 * sigx10_2) * ewu + 0.5 * (depsi2 * mupsi^2 *
      pmusig2/(sigma_sq2)))/(sigma_sq2) - (((sigx8_2^2 * ewu/starsq2 + sigx7_2/ssq2) *
      musi2 + ((0.5 * (ewu/(sigma_sq2)) - 0.5 * (0.5 * usq2 + ewu/(sigma_sq2))) *
      usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * sigx30_2 + 2 * (S * (epsilon))) *
      ewu)/ssq2 + S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 *
      (sigx29_2 * ewu))/(sigma_sq2)) * depsi2) * prC/sqrt(sigma_sq2) + (((0.5 *
      (sigx25_1 * ewu) + 0.5 * pmusig1) * depsi1 * mupsi^2/(sigma_sq1)^2 -
      ((sigx8_1^2 * ewu * musig1 + ((0.5 * (ewu/(sigma_sq1)) - 0.5 * (0.5 *
        usq1 + ewu/(sigma_sq1))) * usq1 * ewv1 * musi1/sigmastar1 - sigx7_1 *
        (2 * sigx30_1 + 2 * (S * (epsilon))) * ewu)/ssq1) * depsi1 + sigx8_1 *
        (0.5 * (depsi1 * ewu * mupsi^2/(sigma_sq1)^2) + depsi1)) * dmusig1)/wzdeno -
      ((0.5 * (sigx10_1/(sigma_sq1)) + 0.5 * sigx31_1) * ewu + 0.5 * sigx9_1)) *
      ewz/sqrt(sigma_sq1) - sigx11^2 * ewu/sigx3) * ewu/sigx3 + 0.5 * (mu *
      (0.5 * (mu^2/(ewu_h^3 * pmu)) - (0.5 * (ewu_h * pmu) - 0.5 * (mu * dmu))/(ewu_h *
        pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 *
      (depsi1 * mupsi^2/(sigma_sq1))) * sigx13_1/(sigma_sq1) - ((0.5 * (usq1 *
      ewv1/(sigma_sq1)) + 0.5 * ((ewu/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
      1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu * sigx7_1 - sigx12_1 *
      (2 * sigx30_1 + S * (epsilon))) * depsi1/ssq1) * dmusig1 + 0.5 * (sigx25_1 *
      depsi1 * mupsi^2/(sigma_sq1)^2))/wzdeno - 0.5 * sigx31_1)/s3sq1 - (sigx11 *
      sqrt(sigma_sq1) + 0.5 * (sigx3/sqrt(sigma_sq1))) * sigx15_1/s3sq1^2) *
      ewu * ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 *
      (usq2 * ewv2/(sigma_sq2)) + 0.5 * ((ewu/(sigma_sq2) - 1) * ewv2/(sigma_sq2) +
      1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 *
      (2 * sigx30_2 + S * (epsilon)))/ssq2) * dmusig2 + 0.5 * (sigx29_2/(sigma_sq2))) *
      depsi2 + (0.5 * sigx25_2 + 0.5 * sigx14_2) * depsi2 * mupsi^2/(sigma_sq2)^2)/s3sq2 -
      (sigx11 * sqrt(sigma_sq2) + 0.5 * (sigx3/sqrt(sigma_sq2))) * sigx15_2/s3sq2^2) *
      prC * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (sigx16 * depsi1 * mupsi^2/(sigma_sq1)^2) -
      ((0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzsrq1^2) * ewz +
        0.5 * (wzdeno/sqrt(sigma_sq1))) * depsi1/wzsrq1^2) * pmusig1 - (sigx11 *
      sigx17/sigx3 + sigx8_1 * sigx16 * dmusig1 * depsi1 + ((0.5 * sigx6_2 -
      sigx8_2 * dmusig2 * depsi2)/wzdeno - 0.5 * sigx31_2) * prC/sqrt(sigma_sq2))) *
      ewu * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi1 * mupsi^2/(sigma_sq1)) - depsi1 *
      musi1 * sigx13_1/sigmastar1) * ewv1 * sigx13_1/(sigma_sq1) + depsi1 *
      (mu/starsq1 - (((3 * (mu) - 2 * (sigx12_1 * (sigma_sq1) * musi1 * sigmastar1/ssq1)) *
        ewv1 - S * ewu * (epsilon)) * sigx12_1 + (0.5 * (ewv1/(sigma_sq1)) -
        0.5 * (0.5 * vsq1 + ewv1/(sigma_sq1))) * vsq1 * ewu * musi1/sigmastar1)/ssq1)) *
      dmusig1 + (0.5 * (((0.5 * (mupsi^2/(sigma_sq1)) - 2) * pmusig1/(sigma_sq1) +
      dmusig1 * sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi^2/(sigma_sq1)^2)/wzdeno -
      (0.5 * (((dmusig1 * sigx13_1 - wzdeno^2 * pmusig1/wzsrq1^2) * depsi1 +
        0.5 * sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * wzdeno/wzsrq1^2)/s3sq1 -
      (sigx15_1 * ewz + 0.5 * (sigx3/sqrt(sigma_sq1))) * sigx15_1 * ewv1/s3sq1^2) *
      ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx15_1 * sigx15_2 * prC * ewv1 * ewv2 * ewz * sqrt(sigma_sq2)/(s3sq2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (sigx16 * depsi1 * mupsi^2/(sigma_sq1)^2) - ((0.5/sqrt(sigma_sq1) -
      wzdeno^2 * sqrt(sigma_sq1)/wzsrq1^2) * ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) *
      depsi1/wzsrq1^2) * pmusig1 + sigx16 * dmusig1 * depsi1 * sigx13_1 - sigx15_1 *
      sigx17 * ewz/s3sq1) * ewv1 * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/(sigma_sq2)) - 2) * pmusig2/(sigma_sq2) + dmusig2 *
      sigx13_2) + 0.5 * sigx14_2) * ewv2 + 0.5 * pmusig2) * depsi2 * mupsi^2/(sigma_sq2)^2 +
      (dmusig2 * (mu/starsq2 - ((sigx12_2 * (2 * (mu) - 2 * (sigx12_2 * (sigma_sq2) *
        musi2 * sigmastar2/ssq2)) * ewv2 + (0.5 * (ewv2/(sigma_sq2)) - 0.5 *
        (0.5 * vsq2 + ewv2/(sigma_sq2))) * vsq2 * ewu * musi2/sigmastar2)/ssq2 +
        (sigx12_2/ssq2 + ewv2 * sigx13_2^2/starsq2) * musi2)) - (0.5 * ((dmusig2 *
        sigx13_2 - pmusig2/(sigma_sq2)) * ewv2) + 0.5 * pmusig2)/(sigma_sq2)) *
        depsi2)/s3sq2 - (sigx15_2 * prC + 0.5 * (sigx3/sqrt(sigma_sq2))) *
      sigx15_2 * ewv2/s3sq2^2) * prC * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx17 * sigx15_2/sigx3 + (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)/wzdeno -
      0.5 * sigx31_2) * prC * ewv2 * ewz/s3sq2), FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * ((prC * (1/(wzdeno^2 * sqrt(sigma_sq2)) + sqrt(sigma_sq2)/wzsq2^2) *
      depsi2 * pmusig2 - (sigx17^2/sigx3 + (2 - 2 * (wzdeno * (sigma_sq1) *
      ewz/wzsrq1^2)) * depsi1 * pmusig1 * sqrt(sigma_sq1)/wzsrq1^2)) * ewz +
      sigx16 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/wzsq2) * ewz/sigx3,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesscnsftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/sqrt(sigma_sq2) + ewz1 * depsi1 * pmusig1/sqrt(sigma_sq1))
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- (ewz2 * sigx4_2/ssqq2 + ewz1 * sigx4_1/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (ewz2 * sigx11_2/sqrt(sigma_sq2) + sigx11_1 * ewz1/sqrt(sigma_sq1))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * ewz1/ssqq1 +
      (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
        depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
        ewz2/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (ewz2 * (dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) *
      ewv2/starsq2 - ((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
      depsi2)/ssqq2 + ewz1 * (dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) *
      depsi1)/ssqq1 - sigx3 * sigx5/sigx2)/sigx2, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 *
      dmusig1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) *
      dmusig1 * depsi1) * ewz1/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 +
      0.5/(sigma_sq2^2 * sigmastar2)) * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      dmusig2 * depsi2) * ewz2/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) * ewu/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((((0.5 *
    ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) - 0.5 *
    pmusig1)/sigma_sq1 + dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) *
    dmusig1 * depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) -
    0.5 * pmusig2)/sigma_sq2 + dmusig2 * sigx13_2) * depsi2 * mupsi/sigma_sq2 -
    ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 *
      ewv2)) * dmusig2 * depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((sigx1_1/ssqq1 - sigx1_2/ssqq2)/(pi * sigx2 * ((Wz)^2 + 1)) - pi * sigx3 *
    ((Wz)^2 + 1) * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (dmu * (dmu/(ewusr * pmu)^2 + mu/(ewusr^3 * pmu)) -
      ((((1 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) *
        depsi1 + dmusig1 * (depsi1 * musi1/ewu + depsi1 * mupsi) * ewv1/starsq1) *
        ewz1/ssqq1 + (((1 - mupsi^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 *
        mupsi/starsq2) * depsi2 + dmusig2 * (depsi2 * musi2/ewu + depsi2 *
        mupsi) * ewv2/starsq2) * ewz2/ssqq2 + sigx5^2/sigx2)/sigx2), FUN = "*"),
    muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    dmusig1 + (0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
    0.5 * pmusig1)/sigma_sq1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewv1 - sigx8_1 * musi1/(ewu * sigma_sq1)) *
    dmusig1 * depsi1) * ewz1/sqrt(sigma_sq1) + ((sigx8_2 * dmusig2 + (0.5 * ((2 -
    mupsi^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi/starsq2) + 0.5 * pmusig2)/sigma_sq2) *
    depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) *
    ewv2 - sigx8_2 * musi2/(ewu * sigma_sq2)) * dmusig2 * depsi2) * ewz2/sqrt(sigma_sq2) -
    sigx14 * sigx5/sigx2) * ewu/sigx2 + 0.5 * (((1 - mu^2/ewusr^2)/(ewusr * pmu) -
    mu * dmu/(ewusr * pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
      0.5 * pmusig1)/sigma_sq1 - dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 +
      ((1/sigmastar1 - musi1 * sigx13_1/ewu)/sigma_sq1 - (sigx12_1/starsq1^2 +
        0.5/(sigma_sq1^2 * sigmastar1)) * ewv1) * dmusig1 * depsi1)/sigx17_1 -
      sigx5 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * ewz1 * ewv1, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((2 - mupsi^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi/starsq2) + 0.5 * pmusig2)/sigma_sq2 - dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 + ((1/sigmastar2 - musi2 * sigx13_2/ewu)/sigma_sq2 -
      (sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewv2) * dmusig2 *
      depsi2)/sigx17_2 - sigx5 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * ewz2 *
      ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((sigx4_1/ssqq1 - sigx4_2/ssqq2)/(pi * sigx2 *
      ((Wz)^2 + 1)) - pi * sigx5 * ((Wz)^2 + 1) * sigx9/(pi * sigx2 * ((Wz)^2 +
      1))^2), FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
      dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 - 0.5 * sigx11_1) *
      ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 - (((sigx8_1^2 *
      ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
      dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
      ewu))/sigma_sq1) * depsi1) * ewz1/sqrt(sigma_sq1) + ((((0.5 * ((0.5 *
      (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) - 0.5 *
      sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu + 0.5 *
      (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 * ewu/starsq2 +
      sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) - 0.5 * (0.5 * usq2 +
      ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 *
      sigma_sq2 * musi2 * sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu)/starsq2^2 +
      S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 *
      dmusig2 + pmusig2/sigma_sq2) * ewu))/sigma_sq2) * depsi2) * ewz2/sqrt(sigma_sq2) -
      sigx14^2 * ewu/sigx2) * ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 *
      pmu)) - (0.5 * (ewusr * pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)),
    FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_1 * musi1 * sigx13_1/starsq1 - ((0.5 *
      (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 - 1) * ewv1/sigma_sq1 +
      1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon)))/starsq1^2) *
      dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
      depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 -
      sigx8_1 * dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 -
      (sigx14 * sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
      ewz1 * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 *
      (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 +
      1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 *
      (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) *
      dmusig2 + 0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) *
      depsi2 + (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
      ewz2 * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((sigx11_1/sqrt(sigma_sq1) - sigx11_2/sqrt(sigma_sq2))/(pi *
      sigx2 * ((Wz)^2 + 1)) - pi * sigx14 * ((Wz)^2 + 1) * sigx9/(pi * sigx2 *
      ((Wz)^2 + 1))^2) * ewu, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 + dmusig1 * sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 *
      pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 + (dmusig1 * (mu/starsq1 - ((sigx12_1 *
      (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
      ewv1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) *
      vsq1 * ewu * musi1/sigmastar1)/starsq1^2 + (sigx12_1/starsq1^2 + ewv1 *
      sigx13_1^2/starsq1) * musi1)) - (0.5 * ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) *
      ewv1) + 0.5 * pmusig1)/sigma_sq1) * depsi1)/sigx17_1 - (sigx16_1 * ewz1 +
      0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) * ewz1 *
      ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx16_2 * ewz2 * ewz1 * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx16_1 * (1/(pi * sigx2 * ((Wz)^2 + 1)) - pi * ((Wz)^2 + 1) * ewz1 * sigx9/(pi *
    sigx2 * ((Wz)^2 + 1))^2) * ewv1/sqrt(sigma_sq1), FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 + dmusig2 *
      sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 * pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 +
      (dmusig2 * (mu/starsq2 - ((sigx12_2 * (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 *
        musi2 * sigmastar2/starsq2^2)) * ewv2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu * musi2/sigmastar2)/starsq2^2 +
        (sigx12_2/starsq2^2 + ewv2 * sigx13_2^2/starsq2) * musi2)) - (0.5 *
        ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) * ewv2) + 0.5 * pmusig2)/sigma_sq2) *
        depsi2)/sigx17_2 - (sigx16_2 * ewz2 + 0.5 * (sigx2/sqrt(sigma_sq2))) *
      sigx16_2 * ewv2/sigx17_2^2) * ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_2 * (1/(pi * sigx2 * ((Wz)^2 + 1)) + pi * ((Wz)^2 + 1) * ewz2 * sigx9/(pi *
      sigx2 * ((Wz)^2 + 1))^2) * ewv2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = -wHvar * ((2 * (pi * Wz * sigx2) + depsi1 * pmusig1/sqrt(sigma_sq1) -
      depsi2 * pmusig2/sqrt(sigma_sq2)) * sigx9/(pi * sigx2 * ((Wz)^2 + 1))^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesscnsftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/sqrt(sigma_sq2) + depsi1 * pmusig1 * pwZ/sqrt(sigma_sq1))
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - pwZ) * sigx4_2/ssqq2 + sigx4_1 * pwZ/ssqq1)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * pwZ/sqrt(sigma_sq1) + sigx11_2 * (1 - pwZ)/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * pwZ/ssqq1 +
      (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
        depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
        (1 - pwZ)/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((1 - pwZ) * (dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) *
      ewv2/starsq2 - ((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
      depsi2)/ssqq2 + pwZ * (dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) *
      depsi1)/ssqq1 - sigx3 * sigx5/sigx2)/sigx2, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 *
      dmusig1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) *
      dmusig1 * depsi1) * pwZ/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 +
      0.5/(sigma_sq2^2 * sigmastar2)) * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      dmusig2 * depsi2) * (1 - pwZ)/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) *
      ewu/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((((0.5 *
    ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) - 0.5 *
    pmusig1)/sigma_sq1 + dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) *
    dmusig1 * depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    pwZ * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) -
    0.5 * pmusig2)/sigma_sq2 + dmusig2 * sigx13_2) * depsi2 * mupsi/sigma_sq2 -
    ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 *
      ewv2)) * dmusig2 * depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * dwZ/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (dmu * (dmu/(ewusr * pmu)^2 + mu/(ewusr^3 * pmu)) -
      ((((1 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) *
        depsi1 + dmusig1 * (depsi1 * musi1/ewu + depsi1 * mupsi) * ewv1/starsq1) *
        pwZ/ssqq1 + (((1 - mupsi^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 *
        mupsi/starsq2) * depsi2 + dmusig2 * (depsi2 * musi2/ewu + depsi2 *
        mupsi) * ewv2/starsq2) * (1 - pwZ)/ssqq2 + sigx5^2/sigx2)/sigx2),
    FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    dmusig1 + (0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
    0.5 * pmusig1)/sigma_sq1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewv1 - sigx8_1 * musi1/(ewu * sigma_sq1)) *
    dmusig1 * depsi1) * pwZ/sqrt(sigma_sq1) + ((sigx8_2 * dmusig2 + (0.5 * ((2 -
    mupsi^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi/starsq2) + 0.5 * pmusig2)/sigma_sq2) *
    depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) *
    ewv2 - sigx8_2 * musi2/(ewu * sigma_sq2)) * dmusig2 * depsi2) * (1 - pwZ)/sqrt(sigma_sq2) -
    sigx14 * sigx5/sigx2) * ewu/sigx2 + 0.5 * (((1 - mu^2/ewusr^2)/(ewusr * pmu) -
    mu * dmu/(ewusr * pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
      0.5 * pmusig1)/sigma_sq1 - dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 +
      ((1/sigmastar1 - musi1 * sigx13_1/ewu)/sigma_sq1 - (sigx12_1/starsq1^2 +
        0.5/(sigma_sq1^2 * sigmastar1)) * ewv1) * dmusig1 * depsi1)/sigx17_1 -
      sigx5 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * pwZ * ewv1, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((2 - mupsi^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi/starsq2) + 0.5 * pmusig2)/sigma_sq2 - dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 + ((1/sigmastar2 - musi2 * sigx13_2/ewu)/sigma_sq2 -
      (sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewv2) * dmusig2 *
      depsi2)/sigx17_2 - sigx5 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * (1 -
      pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (sigx4_1/ssqq1 - (sigx5 * sigx9/sigx2 + sigx4_2/ssqq2)) *
      dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
      dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 - 0.5 * sigx11_1) *
      ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 - (((sigx8_1^2 *
      ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
      dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
      ewu))/sigma_sq1) * depsi1) * pwZ/sqrt(sigma_sq1) + ((((0.5 * ((0.5 *
      (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) - 0.5 *
      sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu + 0.5 *
      (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 * ewu/starsq2 +
      sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) - 0.5 * (0.5 * usq2 +
      ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 *
      sigma_sq2 * musi2 * sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu)/starsq2^2 +
      S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 *
      dmusig2 + pmusig2/sigma_sq2) * ewu))/sigma_sq2) * depsi2) * (1 - pwZ)/sqrt(sigma_sq2) -
      sigx14^2 * ewu/sigx2) * ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 *
      pmu)) - (0.5 * (ewusr * pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)),
    FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_1 * musi1 * sigx13_1/starsq1 - ((0.5 *
      (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 - 1) * ewv1/sigma_sq1 +
      1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon)))/starsq1^2) *
      dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
      depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 -
      sigx8_1 * dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 -
      (sigx14 * sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
      pwZ * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 *
      (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 +
      1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 *
      (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) *
      dmusig2 + 0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) *
      depsi2 + (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
      (1 - pwZ) * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx11_1/sqrt(sigma_sq1) - (sigx14 * sigx9/sigx2 +
      sigx11_2/sqrt(sigma_sq2))) * dwZ * ewu/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 + dmusig1 * sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 *
      pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 + (dmusig1 * (mu/starsq1 - ((sigx12_1 *
      (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
      ewv1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) *
      vsq1 * ewu * musi1/sigmastar1)/starsq1^2 + (sigx12_1/starsq1^2 + ewv1 *
      sigx13_1^2/starsq1) * musi1)) - (0.5 * ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) *
      ewv1) + 0.5 * pmusig1)/sigma_sq1) * depsi1)/sigx17_1 - (sigx16_1 * pwZ +
      0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) * pwZ *
      ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx16_2 * (1 - pwZ) * pwZ * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx16_1 * (1 - sigx9 * pwZ/sigx2) * dwZ * ewv1/(sigx2 * sqrt(sigma_sq1)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 + dmusig2 *
      sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 * pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 +
      (dmusig2 * (mu/starsq2 - ((sigx12_2 * (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 *
        musi2 * sigmastar2/starsq2^2)) * ewv2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu * musi2/sigmastar2)/starsq2^2 +
        (sigx12_2/starsq2^2 + ewv2 * sigx13_2^2/starsq2) * musi2)) - (0.5 *
        ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) * ewv2) + 0.5 * pmusig2)/sigma_sq2) *
        depsi2)/sigx17_2 - (sigx16_2 * (1 - pwZ) + 0.5 * (sigx2/sqrt(sigma_sq2))) *
      sigx16_2 * ewv2/sigx17_2^2) * (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (((1 - pwZ) * sigx9/sigx2 + 1) * sigx16_2 * dwZ * ewv2/(sigx2 * sqrt(sigma_sq2))),
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = -wHvar * ((sigx9 * dwZ/sigx2 + Wz) * sigx9 * dwZ/sigx2), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesscnsftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi1 <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar)]
  theta <- parm[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu + ewv1)
  sigma_sq2 <- (ewu + ewv2)
  sigmastar1 <- sqrt(ewu * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu * ewv1 - S * ewu * (epsilon))
  musi2 <- (mu * ewv2 - S * ewu * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi <- (mu + S * (epsilon))
  depsi1 <- dnorm(mupsi/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi/sqrt(sigma_sq2), 0, 1)
  dmu <- dnorm(mu/ewusr, 0, 1)
  pmu <- pnorm(mu/ewusr)
  sigx1_1 <- (dmusig1 * depsi1 * ewu/sigmastar1 + depsi1 * mupsi * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu/sigmastar2 + depsi2 * mupsi * pmusig2)
  ssqq1 <- (sigma_sq1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/sqrt(sigma_sq1) + depsi2 * prZ * pmusig2/sqrt(sigma_sq2))
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi * pmusig2)
  sigx5 <- ((1 - prZ) * sigx4_1/ssqq1 + sigx4_2 * prZ/ssqq2)
  sigx6_1 <- (depsi1 * mupsi^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu/sigma_sq1)
  usq2 <- (1 - ewu/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/sqrt(sigma_sq1) - depsi2 * pmusig2/sqrt(sigma_sq2))
  sigx10_1 <- (sigx8_1 * dmusig1 + 0.5 * (pmusig1/sigma_sq1))
  sigx10_2 <- (sigx8_2 * dmusig2 + 0.5 * (pmusig2/sigma_sq2))
  sigx11_1 <- (0.5 * sigx6_1 - sigx10_1 * depsi1)
  sigx11_2 <- (0.5 * sigx6_2 - sigx10_2 * depsi2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14 <- (sigx11_1 * (1 - prZ)/sqrt(sigma_sq1) + sigx11_2 * prZ/sqrt(sigma_sq2))
  sigx15_1 <- (dmusig1 * sigx13_1 - 0.5 * (pmusig1/sigma_sq1))
  sigx15_2 <- (dmusig2 * sigx13_2 - 0.5 * (pmusig2/sigma_sq2))
  sigx16_1 <- (sigx15_1 * depsi1 + 0.5 * sigx6_1)
  sigx16_2 <- (sigx15_2 * depsi2 + 0.5 * sigx6_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar, ncol = nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) * depsi1 +
      dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) * ewu/starsq1) * (1 -
      prZ)/ssqq1 + (((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
      depsi2 + dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) * ewu/starsq2) *
      prZ/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (prZ * (dmusig2 * (depsi2 * mupsi - depsi2 * musi2/ewv2) *
      ewv2/starsq2 - ((mupsi^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) *
      depsi2)/ssqq2 + (1 - prZ) * (dmusig1 * (depsi1 * mupsi - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) *
      depsi1)/ssqq1 - sigx3 * sigx5/sigx2)/sigx2, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((((0.5 * ((mupsi^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu * mupsi/starsq1) - 0.5 * pmusig1)/sigma_sq1 - sigx8_1 *
      dmusig1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 + 0.5/(sigma_sq1^2 *
      sigmastar1)) * ewu - (sigx8_1 * musi1/ewv1 + 1/sigmastar1)/sigma_sq1) *
      dmusig1 * depsi1) * (1 - prZ)/sqrt(sigma_sq1) + (((0.5 * ((mupsi^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) - 0.5 * pmusig2)/sigma_sq2 -
      sigx8_2 * dmusig2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 +
      0.5/(sigma_sq2^2 * sigmastar2)) * ewu - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      dmusig2 * depsi2) * prZ/sqrt(sigma_sq2) - sigx14 * sigx3/sigx2) * ewu/sigx2,
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((((0.5 *
    ((mupsi^2/sigma_sq1 - 2) * pmusig1 + dmusig1 * ewu * mupsi/starsq1) - 0.5 *
    pmusig1)/sigma_sq1 + dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 - ((sigx12_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewu + musi1 * sigx13_1/(sigma_sq1 * ewv1)) *
    dmusig1 * depsi1)/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((0.5 * ((mupsi^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu * mupsi/starsq2) -
    0.5 * pmusig2)/sigma_sq2 + dmusig2 * sigx13_2) * depsi2 * mupsi/sigma_sq2 -
    ((sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewu + musi2 * sigx13_2/(sigma_sq2 *
      ewv2)) * dmusig2 * depsi2)/sigx17_2 - sigx3 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) *
    prZ * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) * prZ * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (dmu * (dmu/(ewusr * pmu)^2 + mu/(ewusr^3 * pmu)) -
      ((((1 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) *
        depsi1 + dmusig1 * (depsi1 * musi1/ewu + depsi1 * mupsi) * ewv1/starsq1) *
        (1 - prZ)/ssqq1 + (((1 - mupsi^2/sigma_sq2) * pmusig2 + dmusig2 *
        ewv2 * mupsi/starsq2) * depsi2 + dmusig2 * (depsi2 * musi2/ewu +
        depsi2 * mupsi) * ewv2/starsq2) * prZ/ssqq2 + sigx5^2/sigx2)/sigx2),
    FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar * ((((sigx8_1 *
    dmusig1 + (0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
    0.5 * pmusig1)/sigma_sq1) * depsi1 * mupsi/sigma_sq1 - ((sigx7_1/starsq1^2 +
    0.5/(sigma_sq1^2 * sigmastar1)) * ewv1 - sigx8_1 * musi1/(ewu * sigma_sq1)) *
    dmusig1 * depsi1) * (1 - prZ)/sqrt(sigma_sq1) + ((sigx8_2 * dmusig2 + (0.5 *
    ((2 - mupsi^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi/starsq2) + 0.5 *
    pmusig2)/sigma_sq2) * depsi2 * mupsi/sigma_sq2 - ((sigx7_2/starsq2^2 + 0.5/(sigma_sq2^2 *
    sigmastar2)) * ewv2 - sigx8_2 * musi2/(ewu * sigma_sq2)) * dmusig2 * depsi2) *
    prZ/sqrt(sigma_sq2) - sigx14 * sigx5/sigx2) * ewu/sigx2 + 0.5 * (((1 - mu^2/ewusr^2)/(ewusr *
    pmu) - mu * dmu/(ewusr * pmu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((2 - mupsi^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi/starsq1) +
      0.5 * pmusig1)/sigma_sq1 - dmusig1 * sigx13_1) * depsi1 * mupsi/sigma_sq1 +
      ((1/sigmastar1 - musi1 * sigx13_1/ewu)/sigma_sq1 - (sigx12_1/starsq1^2 +
        0.5/(sigma_sq1^2 * sigmastar1)) * ewv1) * dmusig1 * depsi1)/sigx17_1 -
      sigx5 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) * (1 - prZ) * ewv1, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((2 - mupsi^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi/starsq2) + 0.5 * pmusig2)/sigma_sq2 - dmusig2 *
      sigx13_2) * depsi2 * mupsi/sigma_sq2 + ((1/sigmastar2 - musi2 * sigx13_2/ewu)/sigma_sq2 -
      (sigx12_2/starsq2^2 + 0.5/(sigma_sq2^2 * sigmastar2)) * ewv2) * dmusig2 *
      depsi2)/sigx17_2 - sigx5 * sigx16_2 * sqrt(sigma_sq2)/sigx17_2^2) * prZ *
      ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar +
    1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (sigx4_1/ssqq1 - (sigx5 * sigx9/sigx2 + sigx4_2/ssqq2)) *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 - sigx8_1 *
      dmusig1) - 0.5 * sigx10_1) * depsi1 * mupsi^2/sigma_sq1 - 0.5 * sigx11_1) *
      ewu + 0.5 * (depsi1 * mupsi^2 * pmusig1/sigma_sq1))/sigma_sq1 - (((sigx8_1^2 *
      ewu/starsq1 + sigx7_1/starsq1^2) * musi1 + ((0.5 * (ewu/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu)/starsq1^2 + S * (epsilon)/starsq1) *
      dmusig1 + (0.5 * pmusig1 - 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1) *
      ewu))/sigma_sq1) * depsi1) * (1 - prZ)/sqrt(sigma_sq1) + ((((0.5 * ((0.5 *
      (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) - 0.5 *
      sigx10_2) * depsi2 * mupsi^2/sigma_sq2 - 0.5 * sigx11_2) * ewu + 0.5 *
      (depsi2 * mupsi^2 * pmusig2/sigma_sq2))/sigma_sq2 - (((sigx8_2^2 * ewu/starsq2 +
      sigx7_2/starsq2^2) * musi2 + ((0.5 * (ewu/sigma_sq2) - 0.5 * (0.5 * usq2 +
      ewu/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 - sigx7_2 * (2 * (sigx7_2 *
      sigma_sq2 * musi2 * sigmastar2/starsq2^2) + 2 * (S * (epsilon))) * ewu)/starsq2^2 +
      S * (epsilon)/starsq2) * dmusig2 + (0.5 * pmusig2 - 0.5 * ((sigx8_2 *
      dmusig2 + pmusig2/sigma_sq2) * ewu))/sigma_sq2) * depsi2) * prZ/sqrt(sigma_sq2) -
      sigx14^2 * ewu/sigx2) * ewu/sigx2 + 0.5 * (mu * (0.5 * (mu^2/(ewusr^3 *
      pmu)) - (0.5 * (ewusr * pmu) - 0.5 * (mu * dmu))/(ewusr * pmu)^2) * dmu)),
    FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_1 * musi1 * sigx13_1/starsq1 - ((0.5 *
      (usq1 * ewv1/sigma_sq1) + 0.5 * ((ewu/sigma_sq1 - 1) * ewv1/sigma_sq1 +
      1 - 0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon)))/starsq1^2) *
      dmusig1 + 0.5 * ((sigx8_1 * dmusig1 + pmusig1/sigma_sq1)/sigma_sq1)) *
      depsi1 + (0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 -
      sigx8_1 * dmusig1) + 0.5 * sigx15_1) * depsi1 * mupsi^2/sigma_sq1^2)/sigx17_1 -
      (sigx14 * sqrt(sigma_sq1) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1/sigx17_1^2) *
      (1 - prZ) * ewu * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((sigx8_2 * musi2 * sigx13_2/starsq2 - ((0.5 *
      (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu/sigma_sq2 - 1) * ewv2/sigma_sq2 +
      1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 + mu * sigx7_2 - sigx12_2 *
      (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) + S * (epsilon)))/starsq2^2) *
      dmusig2 + 0.5 * ((sigx8_2 * dmusig2 + pmusig2/sigma_sq2)/sigma_sq2)) *
      depsi2 + (0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 -
      sigx8_2 * dmusig2) + 0.5 * sigx15_2) * depsi2 * mupsi^2/sigma_sq2^2)/sigx17_2 -
      (sigx14 * sqrt(sigma_sq2) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2/sigx17_2^2) *
      prZ * ewu * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (sigx11_1/sqrt(sigma_sq1) - (sigx14 * sigx9/sigx2 +
      sigx11_2/sqrt(sigma_sq2))) * prZ * ewu * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * ((0.5 * (mupsi^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 + dmusig1 * sigx13_1) + 0.5 * sigx15_1) * ewv1 + 0.5 *
      pmusig1) * depsi1 * mupsi^2/sigma_sq1^2 + (dmusig1 * (mu/starsq1 - ((sigx12_1 *
      (2 * (mu) - 2 * (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) *
      ewv1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 + ewv1/sigma_sq1)) *
      vsq1 * ewu * musi1/sigmastar1)/starsq1^2 + (sigx12_1/starsq1^2 + ewv1 *
      sigx13_1^2/starsq1) * musi1)) - (0.5 * ((dmusig1 * sigx13_1 - pmusig1/sigma_sq1) *
      ewv1) + 0.5 * pmusig1)/sigma_sq1) * depsi1)/sigx17_1 - (sigx16_1 * (1 -
      prZ) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
      (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_1 * sigx16_2 * prZ * (1 - prZ) * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar),
    (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
      2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    sigx16_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ * ewv1 * ewz/(sigx2 * sqrt(sigma_sq1)),
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    ((((0.5 * ((0.5 * (mupsi^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 + dmusig2 *
      sigx13_2) + 0.5 * sigx15_2) * ewv2 + 0.5 * pmusig2) * depsi2 * mupsi^2/sigma_sq2^2 +
      (dmusig2 * (mu/starsq2 - ((sigx12_2 * (2 * (mu) - 2 * (sigx12_2 * sigma_sq2 *
        musi2 * sigmastar2/starsq2^2)) * ewv2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu * musi2/sigmastar2)/starsq2^2 +
        (sigx12_2/starsq2^2 + ewv2 * sigx13_2^2/starsq2) * musi2)) - (0.5 *
        ((dmusig2 * sigx13_2 - pmusig2/sigma_sq2) * ewv2) + 0.5 * pmusig2)/sigma_sq2) *
        depsi2)/sigx17_2 - (sigx16_2 * prZ + 0.5 * (sigx2/sqrt(sigma_sq2))) *
      sigx16_2 * ewv2/sigx17_2^2) * prZ * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar +
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (sigx16_2 * (sigx9 * prZ/sigx2 + 1) * prZ * ewv2 * ewz/(sigx2 * sqrt(sigma_sq2))),
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar + nmuZUvar + nuZUvar +
    2 * nvZVvar + nZHvar), (nXvar + nmuZUvar + nuZUvar + 2 * nvZVvar + 1):(nXvar +
    nmuZUvar + nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1 - (sigx9 * prZ/sigx2 + 1) * ewz) * sigx9 * prZ * ewz/sigx2,
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u

## logit specification class membership
chessmcesftruncnormlike_logit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  ewu1_h <- exp(Wu1/2)
  ewu2_h <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  ssq1 <- ((sigma_sq1) * sigmastar1)
  ssq2 <- ((sigma_sq2) * sigmastar2)
  muvu1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  muvu2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  musig_sq1 <- muvu1/ssq1^2
  musig_sq2 <- muvu2/ssq2^2
  musig1 <- muvu1/ssq1
  musig2 <- muvu2/ssq2
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1, 0, 1)
  dmusig2 <- dnorm(musig2, 0, 1)
  pmu1 <- pnorm(mu1/ewu1_h)
  pmu2 <- pnorm(mu2/ewu2_h)
  dmu1 <- dnorm(mu1/ewu1_h, 0, 1)
  dmu2 <- dnorm(mu2/ewu2_h, 0, 1)
  musi1 <- (mu1 + S * (epsilon))
  musi2 <- (mu2 + S * (epsilon))
  epsi1 <- musi1/sqrt(sigma_sq1)
  epsi2 <- musi2/sqrt(sigma_sq2)
  depsi1 <- dnorm(epsi1, 0, 1)
  depsi2 <- dnorm(epsi2, 0, 1)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * musi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * musi2 * pmusig2)
  spmu1 <- (wzdeno * (sigma_sq1) * pmu1 * sqrt(sigma_sq1))
  spmu2 <- ((sigma_sq2) * pmu2 * sqrt(sigma_sq2))
  sigx2 <- (prC * sigx1_2/spmu2 + sigx1_1 * ewz/spmu1)
  pmusq1 <- (wzdeno * pmu1 * sqrt(sigma_sq1))
  pmusq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx3 <- (prC * depsi2 * pmusig2/pmusq2 + depsi1 * ewz * pmusig1/pmusq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * musi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * musi2 * pmusig2)
  wup1 <- (pmusq1^2 * ewu1_h)
  wup2 <- (ewu2_h * pmusq2^2)
  sigx5_1 <- (sigx4_1/spmu1 - wzdeno * depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/wup1)
  sigx5_2 <- (sigx4_2/spmu2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/wup2)
  sigx6_1 <- (depsi1 * musi1^2 * pmusig1/(sigma_sq1)^2)
  sigx6_2 <- (depsi2 * musi2^2 * pmusig2/(sigma_sq2)^2)
  prU1 <- (1 - ewu1/(sigma_sq1))
  prU2 <- (1 - ewu2/(sigma_sq2))
  pvs1 <- (0.5 * (prU1 * ewv1/sigmastar1) + sigmastar1)
  pvs2 <- (0.5 * (prU2 * ewv2/sigmastar2) + sigmastar2)
  sigx7_1 <- (pvs1 * musig_sq1 + S * (epsilon)/ssq1)
  sigx7_2 <- (pvs2 * musig_sq2 + S * (epsilon)/ssq2)
  sigx8_1 <- (0.5 * sigx6_1 - sigx7_1 * dmusig1 * depsi1)
  sigx8_2 <- (0.5 * sigx6_2 - sigx7_2 * dmusig2 * depsi2)
  sigx9_1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  sigx9_2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx10_1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewu1_h)
  sigx10_2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewu2_h)
  sigx11_1 <- (sigx8_1 * ewu1/pmusq1 - (0.5 * sigx9_1 - 0.5 * sigx10_1) * wzdeno *
    depsi1 * pmusig1/pmusq1^2)
  sigx11_2 <- (sigx8_2 * ewu2/pmusq2 - (0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2 *
    pmusig2/pmusq2^2)
  prV1 <- (1 - ewv1/(sigma_sq1))
  prV2 <- (1 - ewv2/(sigma_sq2))
  sigx12_1 <- (0.5 * (prV1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (prV2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/ssq1 - sigx12_1 * musig_sq1)
  sigx13_2 <- (mu2/ssq2 - sigx12_2 * musig_sq2)
  sigx14_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx14_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx15_1 <- (wzdeno * depsi1 * pmusig1 * pmu1/pmusq1^2)
  sigx15_2 <- (depsi2 * pmusig2 * pmu2/pmusq2^2)
  sigx16_1 <- (sigx14_1/(wzdeno * pmu1) - 0.5 * sigx15_1)
  sigx16_2 <- (sigx14_2/pmu2 - 0.5 * sigx15_2)
  s3sq1 <- (sigx3 * sqrt(sigma_sq1))
  s3sq2 <- (sigx3 * sqrt(sigma_sq2))
  sigx17 <- (1/pmusq1 - ewz * pmu1 * sqrt(sigma_sq1)/pmusq1^2)
  sigx18 <- (wzdeno * pmu2 * sqrt(sigma_sq2))
  sigx19 <- (sigx17 * depsi1 * pmusig1 - prC * depsi2 * pmusig2/sigx18)
  sigx20_1 <- ((musi1^2/(sigma_sq1) - 1) * pmusig1 + dmusig1 * ewu1 * musi1/ssq1)
  sigx20_2 <- ((musi2^2/(sigma_sq2) - 1) * pmusig2 + dmusig2 * ewu2 * musi2/ssq2)
  sigx21_1 <- (depsi1 * musi1 - depsi1 * muvu1/ewv1)
  sigx21_2 <- (depsi2 * musi2 - depsi2 * muvu2/ewv2)
  musq1 <- (musi1^2/(sigma_sq1) - 2)
  musq2 <- (musi2^2/(sigma_sq2) - 2)
  sigx22_1 <- ((musq1 * pmusig1 + dmusig1 * ewu1 * musi1/ssq1) * depsi1 * musi1/(sigma_sq1)^2)
  sigx22_2 <- ((musq2 * pmusig2 + dmusig2 * ewu2 * musi2/ssq2) * depsi2 * musi2/(sigma_sq2)^2)
  sigx23_1 <- sigx7_1 * depsi1 * musi1/(sigma_sq1)
  sigx23_2 <- sigx7_2 * depsi2 * musi2/(sigma_sq2)
  sigx24_1 <- (depsi1 * muvu1/ewu1 + depsi1 * musi1)
  sigx24_2 <- (depsi2 * muvu2/ewu2 + depsi2 * musi2)
  sigx25_1 <- (wzdeno^2 * (sigma_sq1) * pmu1^2/pmusq1^2)
  sigx26_1 <- wzdeno^2 * pmu1^2 * sqrt(sigma_sq1)/pmusq1^2
  sigx27_1 <- (((2 - musi1^2/(sigma_sq1)) * pmusig1 + dmusig1 * ewv1 * musi1/ssq1) *
    depsi1 * musi1/(sigma_sq1)^2)
  sigx27_2 <- (((2 - musi2^2/(sigma_sq2)) * pmusig2 + dmusig2 * ewv2 * musi2/ssq2) *
    depsi2 * musi2/(sigma_sq2)^2)
  hessll <- matrix(nrow = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S^2 * ((sigx20_1 * depsi1 + dmusig1 * sigx21_1 * ewu1/ssq1) * ewz/spmu1 +
    (sigx20_2 * depsi2 + dmusig2 * sigx21_2 * ewu2/ssq2) * prC/spmu2 - sigx2^2/sigx3)/sigx3,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((dmusig1 * sigx21_1 * ewv1/ssq1 - sigx20_1 * depsi1)/pmusq1 -
      wzdeno * sigx1_1 * dmu1 * sqrt(sigma_sq1)/wup1)/(sigma_sq1) - sigx2 *
      sigx5_1/sigx3) * ewz/sigx3, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((dmusig2 * sigx21_2 * ewv2/ssq2 - sigx20_2 *
      depsi2)/pmusq2 - sigx1_2 * dmu2 * sqrt(sigma_sq2)/wup2)/(sigma_sq2) -
      sigx2 * sigx5_2/sigx3) * prC/sigx3, FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * sigx22_1 - (sigx23_1 + (pvs1 * ewu1/ssq1^2 -
      (sigx7_1 * muvu1/ewv1 + 1/sigmastar1)/(sigma_sq1)) * depsi1) * dmusig1) *
      ewu1/pmusq1 - (sigx11_1 * sigx2/sigx3 + (0.5 * sigx9_1 - 0.5 * sigx10_1) *
      wzdeno * sigx1_1/(pmusq1^2 * (sigma_sq1)))) * ewz/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((0.5 *
    sigx22_2 - (sigx23_2 + (pvs2 * ewu2/ssq2^2 - (sigx7_2 * muvu2/ewv2 + 1/sigmastar2)/(sigma_sq2)) *
    depsi2) * dmusig2) * ewu2/pmusq2 - (sigx11_2 * sigx2/sigx3 + (0.5 * sigx9_2 -
    0.5 * sigx10_2) * sigx1_2/((sigma_sq2) * pmusq2^2))) * prC/sigx3, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((((sigx21_1 * sigx13_1/(sigma_sq1) - sigx12_1 * depsi1 * ewu1/ssq1^2) *
    dmusig1 + 0.5 * sigx22_1)/(wzdeno * pmu1) - 0.5 * (wzdeno * sigx1_1 * pmu1/(pmusq1^2 *
    (sigma_sq1))))/s3sq1 - sigx16_1 * sigx2 * sqrt(sigma_sq1)/s3sq1^2) * ewv1 *
    ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ((((sigx21_2 * sigx13_2/(sigma_sq2) - sigx12_2 * depsi2 *
      ewu2/ssq2^2) * dmusig2 + 0.5 * sigx22_2)/pmu2 - 0.5 * (sigx1_2 * pmu2/((sigma_sq2) *
      pmusq2^2)))/s3sq2 - sigx16_2 * sigx2 * sqrt(sigma_sq2)/s3sq2^2) * prC *
      ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx17 * sigx1_1/(sigma_sq1) - (sigx2 *
      sigx19/sigx3 + prC * sigx1_2/(wzdeno * (sigma_sq2) * pmu2 * sqrt(sigma_sq2)))) *
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((((1 - musi1^2/(sigma_sq1)) * pmusig1 + dmusig1 *
      ewv1 * musi1/ssq1) * depsi1 + dmusig1 * sigx24_1 * ewv1/ssq1)/spmu1 +
      ((sigx4_1/(spmu1^2 * ewu1_h) - 2 * (wzdeno^2 * depsi1 * dmu1 * pmusig1 *
        pmu1/wup1^2)) * (sigma_sq1) + (dmusig1 * depsi1 * ewv1/ssq1 - (depsi1 *
        musi1/(sigma_sq1) + mu1 * depsi1/ewu1_h^2) * pmusig1)/wup1) * wzdeno *
        dmu1 * sqrt(sigma_sq1) + sigx5_1^2 * ewz/sigx3) * ewz/sigx3), FUN = "*"),
    muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx5_1 * sigx5_2 * prC * ewz/sigx3^2), FUN = "*"),
    muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 *
    nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * sigx27_1 - ((pvs1 * ewv1/ssq1^2 - sigx7_1 * muvu1/(ewu1 * (sigma_sq1))) *
      depsi1 - sigx23_1) * dmusig1)/pmusq1 - sigx8_1 * wzdeno * dmu1 * sqrt(sigma_sq1)/wup1) *
      ewu1 - ((((0.5 * (ewu1/sqrt(sigma_sq1)) - 0.5 * ((1 - mu1^2/ewu1_h^2) *
      sqrt(sigma_sq1))) * depsi1 * dmu1/ewu1_h - (0.5 * sigx9_1 - 0.5 * sigx10_1) *
      depsi1 * musi1/(sigma_sq1)) * pmusig1 + (0.5 * sigx9_1 - 0.5 * sigx10_1) *
      (dmusig1 * ewv1/ssq1 - 2 * (wzdeno^2 * dmu1 * (sigma_sq1) * pmusig1 *
        pmu1/wup1)) * depsi1) * wzdeno/pmusq1^2 + sigx11_1 * sigx5_1 * ewz/sigx3)) *
    ewz/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_2 * sigx5_1 * prC * ewz/sigx3^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/ssq1 - sigx12_1 * ewv1/ssq1^2) * depsi1 -
      sigx24_1 * sigx13_1/(sigma_sq1)) * dmusig1 + 0.5 * sigx27_1)/(wzdeno *
      pmu1) - (sigx14_1 * dmu1/((wzdeno * pmu1)^2 * ewu1_h) + 0.5 * (((1 -
      2 * sigx25_1) * depsi1 * dmu1 * pmusig1/ewu1_h + sigx4_1 * pmu1/(sigma_sq1))/pmusq1^2)) *
      wzdeno)/s3sq1 - sigx16_1 * sigx5_1 * ewz * sqrt(sigma_sq1)/s3sq1^2) *
      ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * sigx5_1 * prC * ewv2 * ewz * sqrt(sigma_sq2)/s3sq2^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (sigx17 * dmusig1 * depsi1 * ewv1/ssq1 - ((((2 -
      2 * sigx25_1) * ewz + 1) * depsi1 * dmu1 * sqrt(sigma_sq1)/wup1 + sigx17 *
      depsi1 * musi1/(sigma_sq1)) * pmusig1 + sigx19 * sigx5_1 * ewz/sigx3)) *
      ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar +
    2 * nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((((1 -
    musi2^2/(sigma_sq2)) * pmusig2 + dmusig2 * ewv2 * musi2/ssq2) * depsi2 +
    dmusig2 * sigx24_2 * ewv2/ssq2)/spmu2 + ((sigx4_2/(spmu2^2 * ewu2_h) - 2 *
    (depsi2 * dmu2 * pmusig2 * pmu2/wup2^2)) * (sigma_sq2) + (dmusig2 * depsi2 *
    ewv2/ssq2 - (depsi2 * musi2/(sigma_sq2) + mu2 * depsi2/ewu2_h^2) * pmusig2)/wup2) *
    dmu2 * sqrt(sigma_sq2) + sigx5_2^2 * prC/sigx3) * prC/sigx3), FUN = "*"),
    muHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_1 * sigx5_2 * prC * ewz/sigx3^2), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * sigx27_2 - ((pvs2 * ewv2/ssq2^2 - sigx7_2 *
      muvu2/(ewu2 * (sigma_sq2))) * depsi2 - sigx23_2) * dmusig2)/pmusq2 -
      sigx8_2 * dmu2 * sqrt(sigma_sq2)/wup2) * ewu2 - ((((0.5 * (ewu2/sqrt(sigma_sq2)) -
      0.5 * ((1 - mu2^2/ewu2_h^2) * sqrt(sigma_sq2))) * depsi2 * dmu2/ewu2_h -
      (0.5 * sigx9_2 - 0.5 * sigx10_2) * depsi2 * musi2/(sigma_sq2)) * pmusig2 +
      (0.5 * sigx9_2 - 0.5 * sigx10_2) * (dmusig2 * ewv2/ssq2 - 2 * (dmu2 *
        (sigma_sq2) * pmusig2 * pmu2/wup2)) * depsi2)/pmusq2^2 + sigx11_2 *
      sigx5_2 * prC/sigx3)) * prC/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx5_2 * prC * ewv1 * ewz * sqrt(sigma_sq1)/s3sq1^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/ssq2 - sigx12_2 * ewv2/ssq2^2) * depsi2 -
      sigx24_2 * sigx13_2/(sigma_sq2)) * dmusig2 + 0.5 * sigx27_2 - sigx14_2 *
      dmu2/(ewu2_h * pmu2))/pmu2 - 0.5 * (((1 - 2 * ((sigma_sq2) * pmu2^2/pmusq2^2)) *
      depsi2 * dmu2 * pmusig2/ewu2_h + sigx4_2 * pmu2/(sigma_sq2))/pmusq2^2))/s3sq2 -
      sigx16_2 * sigx5_2 * prC * sqrt(sigma_sq2)/s3sq2^2) * prC * ewv2, FUN = "*"),
    vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx19 * sigx5_2/sigx3 + sigx4_2/(wzdeno * (sigma_sq2) * pmu2 * sqrt(sigma_sq2)) -
      wzdeno * depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(sigx18^2 * ewu2_h)) *
      prC * ewz/sigx3), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (musi1^2/(sigma_sq1)) - 2) *
      pmusig1/(sigma_sq1) - sigx7_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      musi1^2/(sigma_sq1)^2 - ((sigx7_1^2 * ewu1 * musig1 + ((0.5 * (ewu1/(sigma_sq1)) -
      0.5 * (0.5 * prU1 + ewu1/(sigma_sq1))) * prU1 * ewv1 * muvu1/sigmastar1 -
      pvs1 * (2 * (pvs1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2) + 2 * (S *
        (epsilon))) * ewu1)/ssq1^2) * depsi1 + sigx7_1 * (0.5 * (depsi1 *
      ewu1 * musi1^2/(sigma_sq1)^2) + depsi1)) * dmusig1)/pmusq1 - sigx8_1 *
      (0.5 * sigx9_1 - 0.5 * sigx10_1) * wzdeno/pmusq1^2) * ewu1 - ((((0.5 *
      (((1 - 0.5 * (ewu1/(sigma_sq1))) * pmu1 - 0.5 * (mu1 * dmu1/ewu1_h)) *
        ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 * ((0.5 * (mu1^2/ewu1_h^2) - 0.5) *
      sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) * dmu1/ewu1_h)) * depsi1 +
      0.5 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) * depsi1 * ewu1 * musi1^2/(sigma_sq1)^2)) *
      pmusig1 - (sigx7_1 * dmusig1 * ewu1 + 2 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) *
      wzdeno^2 * pmusig1 * pmu1 * sqrt(sigma_sq1)/pmusq1^2)) * (0.5 * sigx9_1 -
      0.5 * sigx10_1) * depsi1) * wzdeno/pmusq1^2 + sigx11_1^2 * ewz/sigx3)) *
      ewz/sigx3, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * prC * ewz/sigx3^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (((((sigx7_1 * depsi1 * muvu1/sigmastar1 + 0.5 *
      (depsi1 * musi1^2/(sigma_sq1))) * sigx13_1/(sigma_sq1) - ((0.5 * (prU1 *
      ewv1/(sigma_sq1)) + 0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) +
      1 - 0.5 * (prU1 * prV1))) * muvu1/sigmastar1 + mu1 * pvs1 - sigx12_1 *
      (2 * (pvs1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2) + S * (epsilon))) *
      depsi1/ssq1^2) * dmusig1 + 0.5 * (((0.5 * (musi1^2/(sigma_sq1)) - 2) *
      pmusig1/(sigma_sq1) - sigx7_1 * dmusig1) * depsi1 * musi1^2/(sigma_sq1)^2)) *
      ewu1/(wzdeno * pmu1) + (0.5 * (mu1 * sigx14_1 * dmu1/((wzdeno * pmu1)^2 *
      ewu1_h)) - 0.5 * ((sigx8_1 * ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewu1_h) +
      2 * ((0.5 * sigx9_1 - 0.5 * sigx10_1) * sigx26_1)) * depsi1 * pmusig1)/pmusq1^2)) *
      wzdeno)/s3sq1 - (sigx11_1 * ewz * sqrt(sigma_sq1) + 0.5 * (sigx3 * ewu1/sqrt(sigma_sq1))) *
      sigx16_1/s3sq1^2) * ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx11_1 *
    sigx16_2 * prC * ewv2 * ewz * sqrt(sigma_sq2)/s3sq2^2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((0.5 * (sigx17 * depsi1 * ewu1 * musi1^2/(sigma_sq1)^2) - ((2 - 2 * sigx25_1) *
      ewz + 1) * (0.5 * sigx9_1 - 0.5 * sigx10_1) * depsi1/pmusq1^2) * pmusig1 -
      (sigx7_1 * sigx17 * dmusig1 * depsi1 * ewu1 + sigx11_1 * sigx19 * ewz/sigx3)) *
    ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (musi2^2/(sigma_sq2)) - 2) *
      pmusig2/(sigma_sq2) - sigx7_2 * dmusig2) * ewu2) + 0.5 * pmusig2) * depsi2 *
      musi2^2/(sigma_sq2)^2 - ((sigx7_2^2 * ewu2 * musig2 + ((0.5 * (ewu2/(sigma_sq2)) -
      0.5 * (0.5 * prU2 + ewu2/(sigma_sq2))) * prU2 * ewv2 * muvu2/sigmastar2 -
      pvs2 * (2 * (pvs2 * (sigma_sq2) * muvu2 * sigmastar2/ssq2^2) + 2 * (S *
        (epsilon))) * ewu2)/ssq2^2) * depsi2 + sigx7_2 * (0.5 * (depsi2 *
      ewu2 * musi2^2/(sigma_sq2)^2) + depsi2)) * dmusig2)/pmusq2 - sigx8_2 *
      (0.5 * sigx9_2 - 0.5 * sigx10_2)/pmusq2^2) * ewu2 - ((((0.5 * (((1 -
      0.5 * (ewu2/(sigma_sq2))) * pmu2 - 0.5 * (mu2 * dmu2/ewu2_h)) * ewu2/sqrt(sigma_sq2)) -
      0.5 * (mu2 * ((0.5 * (mu2^2/ewu2_h^2) - 0.5) * sqrt(sigma_sq2) + 0.5 *
        (ewu2/sqrt(sigma_sq2))) * dmu2/ewu2_h)) * depsi2 + 0.5 * ((0.5 *
      sigx9_2 - 0.5 * sigx10_2) * depsi2 * ewu2 * musi2^2/(sigma_sq2)^2)) *
      pmusig2 - (sigx7_2 * dmusig2 * ewu2 + 2 * ((0.5 * sigx9_2 - 0.5 * sigx10_2) *
      pmusig2 * pmu2 * sqrt(sigma_sq2)/pmusq2^2)) * (0.5 * sigx9_2 - 0.5 *
      sigx10_2) * depsi2)/pmusq2^2 + sigx11_2^2 * prC/sigx3)) * prC/sigx3,
    FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx16_1 *
    sigx11_2 * prC * ewv1 * ewz * sqrt(sigma_sq1)/s3sq1^2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx7_2 * depsi2 * muvu2/sigmastar2 + 0.5 * (depsi2 * musi2^2/(sigma_sq2))) *
      sigx13_2/(sigma_sq2) - ((0.5 * (prU2 * ewv2/(sigma_sq2)) + 0.5 * ((ewu2/(sigma_sq2) -
      1) * ewv2/(sigma_sq2) + 1 - 0.5 * (prU2 * prV2))) * muvu2/sigmastar2 +
      mu2 * pvs2 - sigx12_2 * (2 * (pvs2 * (sigma_sq2) * muvu2 * sigmastar2/ssq2^2) +
      S * (epsilon))) * depsi2/ssq2^2) * dmusig2 + 0.5 * (((0.5 * (musi2^2/(sigma_sq2)) -
      2) * pmusig2/(sigma_sq2) - sigx7_2 * dmusig2) * depsi2 * musi2^2/(sigma_sq2)^2)) *
      ewu2 + 0.5 * (mu2 * sigx14_2 * dmu2/(ewu2_h * pmu2)))/pmu2 - 0.5 * ((sigx8_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewu2_h) + 2 * ((0.5 * sigx9_2 - 0.5 *
      sigx10_2) * pmu2^2 * sqrt(sigma_sq2)/pmusq2^2)) * depsi2 * pmusig2)/pmusq2^2))/s3sq2 -
      (sigx11_2 * prC * sqrt(sigma_sq2) + 0.5 * (sigx3 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/s3sq2^2) * prC * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((sigx11_2 * sigx19/sigx3 + sigx8_2 * ewu2/sigx18 - (0.5 *
      sigx9_2 - 0.5 * sigx10_2) * wzdeno * depsi2 * pmusig2/sigx18^2) * prC *
      ewz/sigx3), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (depsi1 * musi1^2/(sigma_sq1)) - depsi1 * muvu1 * sigx13_1/sigmastar1) *
      ewv1 * sigx13_1/(sigma_sq1) + depsi1 * (mu1/ssq1 - (((3 * (mu1) - 2 *
      (sigx12_1 * (sigma_sq1) * muvu1 * sigmastar1/ssq1^2)) * ewv1 - S * ewu1 *
      (epsilon)) * sigx12_1 + (0.5 * (ewv1/(sigma_sq1)) - 0.5 * (0.5 * prV1 +
      ewv1/(sigma_sq1))) * prV1 * ewu1 * muvu1/sigmastar1)/ssq1^2)) * dmusig1 +
      (0.5 * (((0.5 * (musi1^2/(sigma_sq1)) - 2) * pmusig1/(sigma_sq1) + dmusig1 *
        sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * musi1^2/(sigma_sq1)^2)/(wzdeno *
      pmu1) - (0.5 * (((dmusig1 * sigx13_1 - wzdeno^2 * pmusig1 * pmu1^2/pmusq1^2) *
      depsi1 + 0.5 * sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * wzdeno *
      pmu1/pmusq1^2)/s3sq1 - (sigx16_1 * ewz + 0.5 * (sigx3/sqrt(sigma_sq1))) *
      sigx16_1 * ewv1/s3sq1^2) * ewv1 * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * prC * ewv1 * ewv2 * ewz * sqrt(sigma_sq2)/(s3sq2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((0.5 * (sigx17 * depsi1 * musi1^2/(sigma_sq1)^2) -
      ((0.5/sqrt(sigma_sq1) - sigx26_1) * ewz + 0.5 * (wzdeno/sqrt(sigma_sq1))) *
        depsi1 * pmu1/pmusq1^2) * pmusig1 + sigx17 * dmusig1 * depsi1 * sigx13_1 -
      sigx16_1 * sigx19 * ewz/s3sq1) * ewv1 * ewz/sigx3, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * musi2^2/(sigma_sq2)) - depsi2 *
      muvu2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/(sigma_sq2) + depsi2 *
      (mu2/ssq2 - (((3 * (mu2) - 2 * (sigx12_2 * (sigma_sq2) * muvu2 * sigmastar2/ssq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/(sigma_sq2)) -
        0.5 * (0.5 * prV2 + ewv2/(sigma_sq2))) * prV2 * ewu2 * muvu2/sigmastar2)/ssq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (musi2^2/(sigma_sq2)) - 2) * pmusig2/(sigma_sq2) +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * musi2^2/(sigma_sq2)^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/pmusq2^2) * depsi2 +
        0.5 * sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/pmusq2^2)/s3sq2 -
      (sigx16_2 * prC + 0.5 * (sigx3/sqrt(sigma_sq2))) * sigx16_2 * ewv2/s3sq2^2) *
      prC * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx16_2 * sigx19/sigx3 + sigx14_2/(wzdeno *
      pmu2) - 0.5 * (wzdeno * depsi2 * pmusig2 * pmu2/sigx18^2)) * prC * ewv2 *
      ewz/s3sq2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * ((prC * (1/(wzdeno^2 * pmu2 * sqrt(sigma_sq2)) +
      pmu2 * sqrt(sigma_sq2)/sigx18^2) * depsi2 * pmusig2 - (sigx19^2/sigx3 +
      (2 - 2 * (wzdeno * (sigma_sq1) * ewz * pmu1^2/pmusq1^2)) * depsi1 * pmusig1 *
        pmu1 * sqrt(sigma_sq1)/pmusq1^2)) * ewz + sigx17 * depsi1 * pmusig1 -
      prC * depsi2 * pmusig2/sigx18) * ewz/sigx3, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmcesftruncnormlike_cauchit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- (ewz2 * depsi2 * pmusig2/psq2 + ewz1 * depsi1 * pmusig1/psq1)
  sigx3 <- (ewz2 * sigx1_2/ssqq2 + ewz1 * sigx1_1/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  sigx18 <- (pi * sigx2 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      ewz1/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) *
      ewu2/starsq2) * ewz2/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 *
      mupsi1/starsq1) * depsi1)/psq1 - sigx1_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 *
      psq1^2))/sigma_sq1 - sigx3 * sigx5_1/sigx2) * ewz1/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((dmusig2 * (depsi2 * mupsi2 - depsi2 *
      musi2/ewv2) * ewv2/starsq2 - ((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 *
      ewu2 * mupsi2/starsq2) * depsi2)/psq2 - sigx1_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2))/sigma_sq2 - sigx3 * sigx5_2/sigx2) * ewz2/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 *
      depsi1 * mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((0.5 *
    (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 * depsi2 * mupsi2/sigma_sq2 +
    (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      depsi2) * dmusig2) * ewu2/psq2 - (sigx3 * sigx11_2/sigx2 + (0.5 * wup2 -
    0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) * ewz2/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) * sigx13_1/sigma_sq1 - sigx12_1 *
    depsi1 * ewu1/starsq1^2) * dmusig1 + 0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
    dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 *
    (sigx1_1 * pmu1/(sigma_sq1 * psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((depsi2 * mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 -
      sigx12_2 * depsi2 * ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 -
      0.5 * (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 *
      sqrt(sigma_sq2)/sigx17_2^2) * ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx1_1/ssqq1 - sigx1_2/ssqq2)/sigx18 -
      pi * sigx3 * ((Wz)^2 + 1) * sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((((1 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 *
      ewv1 * mupsi1/starsq1) * depsi1 + dmusig1 * (depsi1 * musi1/ewu1 + depsi1 *
      mupsi1) * ewv1/starsq1)/ssqq1 + ((sigx4_1/(ssqq1^2 * ewusr1) - 2 * (depsi1 *
      dmu1 * pmusig1 * pmu1/(ewusr1 * psq1^2)^2)) * sigma_sq1 + (dmusig1 *
      depsi1 * ewv1/starsq1 - (depsi1 * mupsi1/sigma_sq1 + mu1 * depsi1/ewusr1^2) *
      pmusig1)/(ewusr1 * psq1^2)) * dmu1 * sqrt(sigma_sq1) + sigx5_1^2 * ewz1/sigx2) *
      ewz1/sigx2), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx5_1 * sigx5_2 * ewz2 * ewz1/sigx2^2), FUN = "*"),
    muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 *
    nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi1/starsq1) *
      depsi1 * mupsi1/sigma_sq1^2) - ((sigx7_1 * ewv1/starsq1^2 - sigx8_1 *
      musi1/(ewu1 * sigma_sq1)) * depsi1 - sigx8_1 * depsi1 * mupsi1/sigma_sq1) *
      dmusig1)/psq1 - sigx10_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 * psq1^2)) *
      ewu1 - ((((0.5 * (ewu1/sqrt(sigma_sq1)) - 0.5 * ((1 - mu1^2/ewusr1^2) *
      sqrt(sigma_sq1))) * depsi1 * dmu1/ewusr1 - (0.5 * wup1 - 0.5 * mud1) *
      depsi1 * mupsi1/sigma_sq1) * pmusig1 + (0.5 * wup1 - 0.5 * mud1) * (dmusig1 *
      ewv1/starsq1 - 2 * (dmu1 * sigma_sq1 * pmusig1 * pmu1/(ewusr1 * psq1^2))) *
      depsi1)/psq1^2 + sigx11_1 * sigx5_1 * ewz1/sigx2)) * ewz1/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_2 * sigx5_1 * ewz2 * ewz1/sigx2^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq1 - sigx12_1 * ewv1/starsq1^2) *
      depsi1 - (depsi1 * musi1/ewu1 + depsi1 * mupsi1) * sigx13_1/sigma_sq1) *
      dmusig1 + 0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 *
      mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - sigx15_1 * dmu1/(ewusr1 *
      pmu1))/pmu1 - 0.5 * (((1 - 2 * (sigma_sq1 * pmu1^2/psq1^2)) * depsi1 *
      dmu1 * pmusig1/ewusr1 + sigx4_1 * pmu1/sigma_sq1)/psq1^2))/sigx17_1 -
      sigx16_1 * sigx5_1 * ewz1 * sqrt(sigma_sq1)/sigx17_1^2) * ewz1 * ewv1,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * sigx5_1 * ewz2 * ewz1 * ewv2 * sqrt(sigma_sq2)/sigx17_2^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * sigx5_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 *
      sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar +
    2 * nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((((1 -
    mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 +
    dmusig2 * (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * ewv2/starsq2)/ssqq2 +
    ((sigx4_2/(ssqq2^2 * ewusr2) - 2 * (depsi2 * dmu2 * pmusig2 * pmu2/(ewusr2 *
      psq2^2)^2)) * sigma_sq2 + (dmusig2 * depsi2 * ewv2/starsq2 - (depsi2 *
      mupsi2/sigma_sq2 + mu2 * depsi2/ewusr2^2) * pmusig2)/(ewusr2 * psq2^2)) *
      dmu2 * sqrt(sigma_sq2) + sigx5_2^2 * ewz2/sigx2) * ewz2/sigx2), FUN = "*"),
    muHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_1 * sigx5_2 * ewz2 * ewz1/sigx2^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - ((sigx7_2 *
      ewv2/starsq2^2 - sigx8_2 * musi2/(ewu2 * sigma_sq2)) * depsi2 - sigx8_2 *
      depsi2 * mupsi2/sigma_sq2) * dmusig2)/psq2 - sigx10_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2)) * ewu2 - ((((0.5 * (ewu2/sqrt(sigma_sq2)) - 0.5 * ((1 - mu2^2/ewusr2^2) *
      sqrt(sigma_sq2))) * depsi2 * dmu2/ewusr2 - (0.5 * wup2 - 0.5 * mud2) *
      depsi2 * mupsi2/sigma_sq2) * pmusig2 + (0.5 * wup2 - 0.5 * mud2) * (dmusig2 *
      ewv2/starsq2 - 2 * (dmu2 * sigma_sq2 * pmusig2 * pmu2/(ewusr2 * psq2^2))) *
      depsi2)/psq2^2 + sigx11_2 * sigx5_2 * ewz2/sigx2)) * ewz2/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx5_2 * ewz2 * ewz1 * ewv1 * sqrt(sigma_sq1)/sigx17_1^2),
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq2 - sigx12_2 * ewv2/starsq2^2) *
      depsi2 - (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * sigx13_2/sigma_sq2) *
      dmusig2 + 0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 *
      mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - sigx15_2 * dmu2/(ewusr2 *
      pmu2))/pmu2 - 0.5 * (((1 - 2 * (sigma_sq2 * pmu2^2/psq2^2)) * depsi2 *
      dmu2 * pmusig2/ewusr2 + sigx4_2 * pmu2/sigma_sq2)/psq2^2))/sigx17_2 -
      sigx16_2 * sigx5_2 * ewz2 * sqrt(sigma_sq2)/sigx17_2^2) * ewz2 * ewv2,
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (sigx5_2 * (1/sigx18 + pi * ((Wz)^2 + 1) * ewz2 * sigx9/sigx18^2)), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * ewz1/sigx2)) *
      ewz1/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * ewz2 * ewz1/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 *
      (depsi1 * mupsi1^2/sigma_sq1)) * sigx13_1/sigma_sq1 - ((0.5 * (usq1 *
      ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 - 1) * ewv1/sigma_sq1 + 1 -
      0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon))) *
      depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * ewz1 * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx11_1 *
    sigx16_2 * ewz2 * ewz1 * ewv2 * sqrt(sigma_sq2)/sigx17_2^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 * sigx9/sigx18^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * ewu2) + 0.5 * pmusig2) * depsi2 *
      mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 * ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) -
      0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 -
      sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
        2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 + sigx8_2 * (0.5 *
      (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 - sigx10_2 *
      (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 * (ewu2/sigma_sq2)) *
      pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) - 0.5 * (mu2 *
      ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
      ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 +
      2 * ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
      (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * ewz2/sigx2)) *
      ewz2/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx16_1 *
    sigx11_2 * ewz2 * ewz1 * ewv1 * sqrt(sigma_sq1)/sigx17_1^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 +
      mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * ewz2 * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * (1/sigx18 + pi * ((Wz)^2 + 1) * ewz2 * sigx9/sigx18^2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 * sigx13_1/sigmastar1) *
      ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 - (((3 * (mu1) - 2 *
      (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) * ewv1 - S * ewu1 *
      (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 +
      ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) * dmusig1 +
      (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
        sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * ewz1 + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
    ewz1 * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * ewz2 * ewz1 * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1/sigx18 - pi * ((Wz)^2 + 1) * ewz1 *
      sigx9/sigx18^2) * ewv1/sqrt(sigma_sq1), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * ewz2 + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      ewz2 * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * (1/sigx18 + pi * ((Wz)^2 + 1) *
      ewz2 * sigx9/sigx18^2) * ewv2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (pi * Wz * sigx2) + depsi1 * pmusig1/psq1 -
      depsi2 * pmusig2/psq2) * sigx9/sigx18^2), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmcesftruncnormlike_probit <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - pwZ) * depsi2 * pmusig2/psq2 + depsi1 * pmusig1 * pwZ/psq1)
  sigx3 <- ((1 - pwZ) * sigx1_2/ssqq2 + sigx1_1 * pwZ/ssqq1)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      pwZ/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) * ewu2/starsq2) *
      (1 - pwZ)/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 *
      mupsi1/starsq1) * depsi1)/psq1 - sigx1_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 *
      psq1^2))/sigma_sq1 - sigx3 * sigx5_1/sigx2) * pwZ/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((dmusig2 * (depsi2 * mupsi2 - depsi2 *
      musi2/ewv2) * ewv2/starsq2 - ((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 *
      ewu2 * mupsi2/starsq2) * depsi2)/psq2 - sigx1_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2))/sigma_sq2 - sigx3 * sigx5_2/sigx2) * (1 - pwZ)/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 *
      depsi1 * mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((0.5 *
    (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 * depsi2 * mupsi2/sigma_sq2 +
    (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      depsi2) * dmusig2) * ewu2/psq2 - (sigx3 * sigx11_2/sigx2 + (0.5 * wup2 -
    0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) * (1 - pwZ)/sigx2, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) * sigx13_1/sigma_sq1 - sigx12_1 *
    depsi1 * ewu1/starsq1^2) * dmusig1 + 0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
    dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 *
    (sigx1_1 * pmu1/(sigma_sq1 * psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    pwZ * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((depsi2 * mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 -
      sigx12_2 * depsi2 * ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 -
      0.5 * (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 *
      sqrt(sigma_sq2)/sigx17_2^2) * (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) *
      dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((((1 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 *
      ewv1 * mupsi1/starsq1) * depsi1 + dmusig1 * (depsi1 * musi1/ewu1 + depsi1 *
      mupsi1) * ewv1/starsq1)/ssqq1 + ((sigx4_1/(ssqq1^2 * ewusr1) - 2 * (depsi1 *
      dmu1 * pmusig1 * pmu1/(ewusr1 * psq1^2)^2)) * sigma_sq1 + (dmusig1 *
      depsi1 * ewv1/starsq1 - (depsi1 * mupsi1/sigma_sq1 + mu1 * depsi1/ewusr1^2) *
      pmusig1)/(ewusr1 * psq1^2)) * dmu1 * sqrt(sigma_sq1) + sigx5_1^2 * pwZ/sigx2) *
      pwZ/sigx2), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx5_1 * sigx5_2 * (1 - pwZ) * pwZ/sigx2^2),
    FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 *
    nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi1/starsq1) *
      depsi1 * mupsi1/sigma_sq1^2) - ((sigx7_1 * ewv1/starsq1^2 - sigx8_1 *
      musi1/(ewu1 * sigma_sq1)) * depsi1 - sigx8_1 * depsi1 * mupsi1/sigma_sq1) *
      dmusig1)/psq1 - sigx10_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 * psq1^2)) *
      ewu1 - ((((0.5 * (ewu1/sqrt(sigma_sq1)) - 0.5 * ((1 - mu1^2/ewusr1^2) *
      sqrt(sigma_sq1))) * depsi1 * dmu1/ewusr1 - (0.5 * wup1 - 0.5 * mud1) *
      depsi1 * mupsi1/sigma_sq1) * pmusig1 + (0.5 * wup1 - 0.5 * mud1) * (dmusig1 *
      ewv1/starsq1 - 2 * (dmu1 * sigma_sq1 * pmusig1 * pmu1/(ewusr1 * psq1^2))) *
      depsi1)/psq1^2 + sigx11_1 * sigx5_1 * pwZ/sigx2)) * pwZ/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_2 * sigx5_1 * (1 - pwZ) * pwZ/sigx2^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq1 - sigx12_1 * ewv1/starsq1^2) *
      depsi1 - (depsi1 * musi1/ewu1 + depsi1 * mupsi1) * sigx13_1/sigma_sq1) *
      dmusig1 + 0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 *
      mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - sigx15_1 * dmu1/(ewusr1 *
      pmu1))/pmu1 - 0.5 * (((1 - 2 * (sigma_sq1 * pmu1^2/psq1^2)) * depsi1 *
      dmu1 * pmusig1/ewusr1 + sigx4_1 * pmu1/sigma_sq1)/psq1^2))/sigx17_1 -
      sigx16_1 * sigx5_1 * pwZ * sqrt(sigma_sq1)/sigx17_1^2) * pwZ * ewv1,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * sigx5_1 * (1 - pwZ) * pwZ * ewv2 *
      sqrt(sigma_sq2)/sigx17_2^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * sigx5_1 * (1 - sigx9 * pwZ/sigx2) * dwZ/sigx2,
    FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar +
    2 * nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((((1 -
    mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 +
    dmusig2 * (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * ewv2/starsq2)/ssqq2 +
    ((sigx4_2/(ssqq2^2 * ewusr2) - 2 * (depsi2 * dmu2 * pmusig2 * pmu2/(ewusr2 *
      psq2^2)^2)) * sigma_sq2 + (dmusig2 * depsi2 * ewv2/starsq2 - (depsi2 *
      mupsi2/sigma_sq2 + mu2 * depsi2/ewusr2^2) * pmusig2)/(ewusr2 * psq2^2)) *
      dmu2 * sqrt(sigma_sq2) + sigx5_2^2 * (1 - pwZ)/sigx2) * (1 - pwZ)/sigx2),
    FUN = "*"), muHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_1 * sigx5_2 * (1 - pwZ) * pwZ/sigx2^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - ((sigx7_2 *
      ewv2/starsq2^2 - sigx8_2 * musi2/(ewu2 * sigma_sq2)) * depsi2 - sigx8_2 *
      depsi2 * mupsi2/sigma_sq2) * dmusig2)/psq2 - sigx10_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2)) * ewu2 - ((((0.5 * (ewu2/sqrt(sigma_sq2)) - 0.5 * ((1 - mu2^2/ewusr2^2) *
      sqrt(sigma_sq2))) * depsi2 * dmu2/ewusr2 - (0.5 * wup2 - 0.5 * mud2) *
      depsi2 * mupsi2/sigma_sq2) * pmusig2 + (0.5 * wup2 - 0.5 * mud2) * (dmusig2 *
      ewv2/starsq2 - 2 * (dmu2 * sigma_sq2 * pmusig2 * pmu2/(ewusr2 * psq2^2))) *
      depsi2)/psq2^2 + sigx11_2 * sigx5_2 * (1 - pwZ)/sigx2)) * (1 - pwZ)/sigx2,
    FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx5_2 * (1 - pwZ) * pwZ * ewv1 *
      sqrt(sigma_sq1)/sigx17_1^2), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq2 - sigx12_2 * ewv2/starsq2^2) *
      depsi2 - (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * sigx13_2/sigma_sq2) *
      dmusig2 + 0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 *
      mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - sigx15_2 * dmu2/(ewusr2 *
      pmu2))/pmu2 - 0.5 * (((1 - 2 * (sigma_sq2 * pmu2^2/psq2^2)) * depsi2 *
      dmu2 * pmusig2/ewusr2 + sigx4_2 * pmu2/sigma_sq2)/psq2^2))/sigx17_2 -
      sigx16_2 * sigx5_2 * (1 - pwZ) * sqrt(sigma_sq2)/sigx17_2^2) * (1 - pwZ) *
      ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (((1 - pwZ) * sigx9/sigx2 + 1) * sigx5_2 * dwZ/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * pwZ/sigx2)) *
      pwZ/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * (1 - pwZ) * pwZ/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 *
      (depsi1 * mupsi1^2/sigma_sq1)) * sigx13_1/sigma_sq1 - ((0.5 * (usq1 *
      ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 - 1) * ewv1/sigma_sq1 + 1 -
      0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon))) *
      depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * pwZ * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * pwZ * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx11_1 *
    sigx16_2 * (1 - pwZ) * pwZ * ewv2 * sqrt(sigma_sq2)/sigx17_2^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1 - sigx9 * pwZ/sigx2) * dwZ/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * ewu2) + 0.5 * pmusig2) * depsi2 *
      mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 * ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) -
      0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 -
      sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
        2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 + sigx8_2 * (0.5 *
      (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 - sigx10_2 *
      (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 * (ewu2/sigma_sq2)) *
      pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) - 0.5 * (mu2 *
      ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
      ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 +
      2 * ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
      (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * (1 - pwZ)/sigx2)) *
      (1 - pwZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx16_1 *
    sigx11_2 * (1 - pwZ) * pwZ * ewv1 * sqrt(sigma_sq1)/sigx17_1^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 +
      mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * (1 - pwZ) * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * ((1 - pwZ) * sigx9/sigx2 + 1) * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 * sigx13_1/sigmastar1) *
      ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 - (((3 * (mu1) - 2 *
      (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) * ewv1 - S * ewu1 *
      (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 +
      ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) * dmusig1 +
      (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
        sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * pwZ + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
    pwZ * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * (1 - pwZ) * pwZ * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - sigx9 * pwZ/sigx2) * dwZ * ewv1/(sigx2 *
      sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * (1 - pwZ) + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      (1 - pwZ) * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * ((1 - pwZ) * sigx9/sigx2 + 1) *
      dwZ * ewv2/(sigx2 * sqrt(sigma_sq2))), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx9 * dwZ/sigx2 + Wz) * sigx9 * dwZ/sigx2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmcesftruncnormlike_cloglog <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  omega1 <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  omega2 <- parm[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)]
  delta1 <- parm[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)]
  delta2 <- parm[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  sigma_sq1 <- (ewu1 + ewv1)
  sigma_sq2 <- (ewu2 + ewv2)
  sigmastar1 <- sqrt(ewu1 * ewv1/sigma_sq1)
  sigmastar2 <- sqrt(ewu2 * ewv2/sigma_sq2)
  starsq1 <- (sigma_sq1 * sigmastar1)
  starsq2 <- (sigma_sq2 * sigmastar2)
  musi1 <- (mu1 * ewv1 - S * ewu1 * (epsilon))
  musi2 <- (mu2 * ewv2 - S * ewu2 * (epsilon))
  pmusig1 <- pnorm(musi1/starsq1)
  pmusig2 <- pnorm(musi2/starsq2)
  dmusig1 <- dnorm(musi1/starsq1, 0, 1)
  dmusig2 <- dnorm(musi2/starsq2, 0, 1)
  mupsi1 <- (mu1 + S * (epsilon))
  mupsi2 <- (mu2 + S * (epsilon))
  depsi1 <- dnorm(mupsi1/sqrt(sigma_sq1), 0, 1)
  depsi2 <- dnorm(mupsi2/sqrt(sigma_sq2), 0, 1)
  dmu1 <- dnorm(mu1/ewusr1, 0, 1)
  dmu2 <- dnorm(mu2/ewusr2, 0, 1)
  pmu1 <- pnorm(mu1/ewusr1)
  pmu2 <- pnorm(mu2/ewusr2)
  sigx1_1 <- (dmusig1 * depsi1 * ewu1/sigmastar1 + depsi1 * mupsi1 * pmusig1)
  sigx1_2 <- (dmusig2 * depsi2 * ewu2/sigmastar2 + depsi2 * mupsi2 * pmusig2)
  ssqq1 <- (sigma_sq1 * pmu1 * sqrt(sigma_sq1))
  ssqq2 <- (sigma_sq2 * pmu2 * sqrt(sigma_sq2))
  psq1 <- (pmu1 * sqrt(sigma_sq1))
  psq2 <- (pmu2 * sqrt(sigma_sq2))
  sigx2 <- ((1 - prZ) * depsi1 * pmusig1/psq1 + depsi2 * prZ * pmusig2/psq2)
  sigx3 <- ((1 - prZ) * sigx1_1/ssqq1 + sigx1_2 * prZ/ssqq2)
  sigx4_1 <- (dmusig1 * depsi1 * ewv1/sigmastar1 - depsi1 * mupsi1 * pmusig1)
  sigx4_2 <- (dmusig2 * depsi2 * ewv2/sigmastar2 - depsi2 * mupsi2 * pmusig2)
  sigx5_1 <- (sigx4_1/ssqq1 - depsi1 * dmu1 * pmusig1 * sqrt(sigma_sq1)/(ewusr1 *
    psq1^2))
  sigx5_2 <- (sigx4_2/ssqq2 - depsi2 * dmu2 * pmusig2 * sqrt(sigma_sq2)/(ewusr2 *
    psq2^2))
  sigx6_1 <- (depsi1 * mupsi1^2 * pmusig1/sigma_sq1^2)
  sigx6_2 <- (depsi2 * mupsi2^2 * pmusig2/sigma_sq2^2)
  usq1 <- (1 - ewu1/sigma_sq1)
  usq2 <- (1 - ewu2/sigma_sq2)
  sigx7_1 <- (0.5 * (usq1 * ewv1/sigmastar1) + sigmastar1)
  sigx7_2 <- (0.5 * (usq2 * ewv2/sigmastar2) + sigmastar2)
  sigx8_1 <- (sigx7_1 * musi1/starsq1^2 + S * (epsilon)/starsq1)
  sigx8_2 <- (sigx7_2 * musi2/starsq2^2 + S * (epsilon)/starsq2)
  sigx9 <- (depsi1 * pmusig1/psq1 - depsi2 * pmusig2/psq2)
  sigx10_1 <- (0.5 * sigx6_1 - sigx8_1 * dmusig1 * depsi1)
  sigx10_2 <- (0.5 * sigx6_2 - sigx8_2 * dmusig2 * depsi2)
  mud1 <- (mu1 * dmu1 * sqrt(sigma_sq1)/ewusr1)
  mud2 <- (mu2 * dmu2 * sqrt(sigma_sq2)/ewusr2)
  wup1 <- (ewu1 * pmu1/sqrt(sigma_sq1))
  wup2 <- (ewu2 * pmu2/sqrt(sigma_sq2))
  sigx11_1 <- (sigx10_1 * ewu1/psq1 - (0.5 * wup1 - 0.5 * mud1) * depsi1 * pmusig1/psq1^2)
  sigx11_2 <- (sigx10_2 * ewu2/psq2 - (0.5 * wup2 - 0.5 * mud2) * depsi2 * pmusig2/psq2^2)
  vsq1 <- (1 - ewv1/sigma_sq1)
  vsq2 <- (1 - ewv2/sigma_sq2)
  sigx12_1 <- (0.5 * (vsq1 * ewu1/sigmastar1) + sigmastar1)
  sigx12_2 <- (0.5 * (vsq2 * ewu2/sigmastar2) + sigmastar2)
  sigx13_1 <- (mu1/starsq1 - sigx12_1 * musi1/starsq1^2)
  sigx13_2 <- (mu2/starsq2 - sigx12_2 * musi2/starsq2^2)
  sigx14_1 <- (depsi1 * pmusig1 * pmu1/psq1^2)
  sigx14_2 <- (depsi2 * pmusig2 * pmu2/psq2^2)
  sigx15_1 <- (0.5 * sigx6_1 + dmusig1 * depsi1 * sigx13_1)
  sigx15_2 <- (0.5 * sigx6_2 + dmusig2 * depsi2 * sigx13_2)
  sigx16_1 <- (sigx15_1/pmu1 - 0.5 * sigx14_1)
  sigx16_2 <- (sigx15_2/pmu2 - 0.5 * sigx14_2)
  sigx17_1 <- (sigx2 * sqrt(sigma_sq1))
  sigx17_2 <- (sigx2 * sqrt(sigma_sq2))
  hessll <- matrix(nrow = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar,
    ncol = nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    ((((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 * mupsi1/starsq1) *
      depsi1 + dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) * ewu1/starsq1) *
      (1 - prZ)/ssqq1 + (((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 * ewu2 *
      mupsi2/starsq2) * depsi2 + dmusig2 * (depsi2 * mupsi2 - depsi2 * musi2/ewv2) *
      ewu2/starsq2) * prZ/ssqq2 - sigx3^2/sigx2)/sigx2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((dmusig1 * (depsi1 * mupsi1 - depsi1 * musi1/ewv1) *
      ewv1/starsq1 - ((mupsi1^2/sigma_sq1 - 1) * pmusig1 + dmusig1 * ewu1 *
      mupsi1/starsq1) * depsi1)/psq1 - sigx1_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 *
      psq1^2))/sigma_sq1 - sigx3 * sigx5_1/sigx2) * (1 - prZ)/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((dmusig2 * (depsi2 * mupsi2 - depsi2 *
      musi2/ewv2) * ewv2/starsq2 - ((mupsi2^2/sigma_sq2 - 1) * pmusig2 + dmusig2 *
      ewu2 * mupsi2/starsq2) * depsi2)/psq2 - sigx1_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2))/sigma_sq2 - sigx3 * sigx5_2/sigx2) * prZ/sigx2, FUN = "*"),
    muHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
      dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - (sigx8_1 *
      depsi1 * mupsi1/sigma_sq1 + (sigx7_1 * ewu1/starsq1^2 - (sigx8_1 * musi1/ewv1 +
      1/sigmastar1)/sigma_sq1) * depsi1) * dmusig1) * ewu1/psq1 - (sigx3 *
      sigx11_1/sigx2 + (0.5 * wup1 - 0.5 * mud1) * sigx1_1/(sigma_sq1 * psq1^2))) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar * S * ((0.5 *
    (((mupsi2^2/sigma_sq2 - 2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) *
      depsi2 * mupsi2/sigma_sq2^2) - (sigx8_2 * depsi2 * mupsi2/sigma_sq2 +
    (sigx7_2 * ewu2/starsq2^2 - (sigx8_2 * musi2/ewv2 + 1/sigmastar2)/sigma_sq2) *
      depsi2) * dmusig2) * ewu2/psq2 - (sigx3 * sigx11_2/sigx2 + (0.5 * wup2 -
    0.5 * mud2) * sigx1_2/(sigma_sq2 * psq2^2))) * prZ/sigx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (((((depsi1 * mupsi1 - depsi1 * musi1/ewv1) * sigx13_1/sigma_sq1 - sigx12_1 *
    depsi1 * ewu1/starsq1^2) * dmusig1 + 0.5 * (((mupsi1^2/sigma_sq1 - 2) * pmusig1 +
    dmusig1 * ewu1 * mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2))/pmu1 - 0.5 *
    (sigx1_1 * pmu1/(sigma_sq1 * psq1^2)))/sigx17_1 - sigx3 * sigx16_1 * sqrt(sigma_sq1)/sigx17_1^2) *
    (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (((((depsi2 * mupsi2 - depsi2 * musi2/ewv2) * sigx13_2/sigma_sq2 -
      sigx12_2 * depsi2 * ewu2/starsq2^2) * dmusig2 + 0.5 * (((mupsi2^2/sigma_sq2 -
      2) * pmusig2 + dmusig2 * ewu2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2))/pmu2 -
      0.5 * (sigx1_2 * pmu2/(sigma_sq2 * psq2^2)))/sigx17_2 - sigx3 * sigx16_2 *
      sqrt(sigma_sq2)/sigx17_2^2) * prZ * ewv2, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (sigx1_1/ssqq1 - (sigx3 * sigx9/sigx2 + sigx1_2/ssqq2)) *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (((((1 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 *
      ewv1 * mupsi1/starsq1) * depsi1 + dmusig1 * (depsi1 * musi1/ewu1 + depsi1 *
      mupsi1) * ewv1/starsq1)/ssqq1 + ((sigx4_1/(ssqq1^2 * ewusr1) - 2 * (depsi1 *
      dmu1 * pmusig1 * pmu1/(ewusr1 * psq1^2)^2)) * sigma_sq1 + (dmusig1 *
      depsi1 * ewv1/starsq1 - (depsi1 * mupsi1/sigma_sq1 + mu1 * depsi1/ewusr1^2) *
      pmusig1)/(ewusr1 * psq1^2)) * dmu1 * sqrt(sigma_sq1) + sigx5_1^2 * (1 -
      prZ)/sigx2) * (1 - prZ)/sigx2), FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx5_1 * sigx5_2 * prZ * (1 - prZ)/sigx2^2),
    FUN = "*"), muHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 1):(nXvar + 2 *
    nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = wHvar *
    (((0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 * mupsi1/starsq1) *
      depsi1 * mupsi1/sigma_sq1^2) - ((sigx7_1 * ewv1/starsq1^2 - sigx8_1 *
      musi1/(ewu1 * sigma_sq1)) * depsi1 - sigx8_1 * depsi1 * mupsi1/sigma_sq1) *
      dmusig1)/psq1 - sigx10_1 * dmu1 * sqrt(sigma_sq1)/(ewusr1 * psq1^2)) *
      ewu1 - ((((0.5 * (ewu1/sqrt(sigma_sq1)) - 0.5 * ((1 - mu1^2/ewusr1^2) *
      sqrt(sigma_sq1))) * depsi1 * dmu1/ewusr1 - (0.5 * wup1 - 0.5 * mud1) *
      depsi1 * mupsi1/sigma_sq1) * pmusig1 + (0.5 * wup1 - 0.5 * mud1) * (dmusig1 *
      ewv1/starsq1 - 2 * (dmu1 * sigma_sq1 * pmusig1 * pmu1/(ewusr1 * psq1^2))) *
      depsi1)/psq1^2 + sigx11_1 * sigx5_1 * (1 - prZ)/sigx2)) * (1 - prZ)/sigx2,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    (sigx11_2 * sigx5_1 * prZ * (1 - prZ)/sigx2^2), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq1 - sigx12_1 * ewv1/starsq1^2) *
      depsi1 - (depsi1 * musi1/ewu1 + depsi1 * mupsi1) * sigx13_1/sigma_sq1) *
      dmusig1 + 0.5 * (((2 - mupsi1^2/sigma_sq1) * pmusig1 + dmusig1 * ewv1 *
      mupsi1/starsq1) * depsi1 * mupsi1/sigma_sq1^2) - sigx15_1 * dmu1/(ewusr1 *
      pmu1))/pmu1 - 0.5 * (((1 - 2 * (sigma_sq1 * pmu1^2/psq1^2)) * depsi1 *
      dmu1 * pmusig1/ewusr1 + sigx4_1 * pmu1/sigma_sq1)/psq1^2))/sigx17_1 -
      sigx16_1 * sigx5_1 * (1 - prZ) * sqrt(sigma_sq1)/sigx17_1^2) * (1 - prZ) *
      ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * sigx5_1 * prZ * (1 - prZ) * ewv2 *
      sqrt(sigma_sq2)/sigx17_2^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * sigx5_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ *
      ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + nmuZUvar + 1):(nXvar +
    2 * nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar * (((((1 -
    mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 +
    dmusig2 * (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * ewv2/starsq2)/ssqq2 +
    ((sigx4_2/(ssqq2^2 * ewusr2) - 2 * (depsi2 * dmu2 * pmusig2 * pmu2/(ewusr2 *
      psq2^2)^2)) * sigma_sq2 + (dmusig2 * depsi2 * ewv2/starsq2 - (depsi2 *
      mupsi2/sigma_sq2 + mu2 * depsi2/ewusr2^2) * pmusig2)/(ewusr2 * psq2^2)) *
      dmu2 * sqrt(sigma_sq2) + sigx5_2^2 * prZ/sigx2) * prZ/sigx2), FUN = "*"),
    muHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_1 * sigx5_2 * prZ * (1 - prZ)/sigx2^2), FUN = "*"),
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 +
      dmusig2 * ewv2 * mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - ((sigx7_2 *
      ewv2/starsq2^2 - sigx8_2 * musi2/(ewu2 * sigma_sq2)) * depsi2 - sigx8_2 *
      depsi2 * mupsi2/sigma_sq2) * dmusig2)/psq2 - sigx10_2 * dmu2 * sqrt(sigma_sq2)/(ewusr2 *
      psq2^2)) * ewu2 - ((((0.5 * (ewu2/sqrt(sigma_sq2)) - 0.5 * ((1 - mu2^2/ewusr2^2) *
      sqrt(sigma_sq2))) * depsi2 * dmu2/ewusr2 - (0.5 * wup2 - 0.5 * mud2) *
      depsi2 * mupsi2/sigma_sq2) * pmusig2 + (0.5 * wup2 - 0.5 * mud2) * (dmusig2 *
      ewv2/starsq2 - 2 * (dmu2 * sigma_sq2 * pmusig2 * pmu2/(ewusr2 * psq2^2))) *
      depsi2)/psq2^2 + sigx11_2 * sigx5_2 * prZ/sigx2)) * prZ/sigx2, FUN = "*"),
    uHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_1 * sigx5_2 * prZ * (1 - prZ) * ewv1 *
      sqrt(sigma_sq1)/sigx17_1^2), FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = wHvar * (((((1/starsq2 - sigx12_2 * ewv2/starsq2^2) *
      depsi2 - (depsi2 * musi2/ewu2 + depsi2 * mupsi2) * sigx13_2/sigma_sq2) *
      dmusig2 + 0.5 * (((2 - mupsi2^2/sigma_sq2) * pmusig2 + dmusig2 * ewv2 *
      mupsi2/starsq2) * depsi2 * mupsi2/sigma_sq2^2) - sigx15_2 * dmu2/(ewusr2 *
      pmu2))/pmu2 - 0.5 * (((1 - 2 * (sigma_sq2 * pmu2^2/psq2^2)) * depsi2 *
      dmu2 * pmusig2/ewusr2 + sigx4_2 * pmu2/sigma_sq2)/psq2^2))/sigx17_2 -
      sigx16_2 * sigx5_2 * prZ * sqrt(sigma_sq2)/sigx17_2^2) * prZ * ewv2,
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + 2 * nmuZUvar), (nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = -wHvar *
    ((sigx9 * prZ/sigx2 + 1) * sigx5_2 * prZ * ewz/sigx2), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * ewu1) + 0.5 * pmusig1) * depsi1 *
      mupsi1^2/sigma_sq1^2 - ((sigx8_1^2 * ewu1 * musi1/starsq1 + ((0.5 * (ewu1/sigma_sq1) -
      0.5 * (0.5 * usq1 + ewu1/sigma_sq1)) * usq1 * ewv1 * musi1/sigmastar1 -
      sigx7_1 * (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) +
        2 * (S * (epsilon))) * ewu1)/starsq1^2) * depsi1 + sigx8_1 * (0.5 *
      (depsi1 * ewu1 * mupsi1^2/sigma_sq1^2) + depsi1)) * dmusig1)/psq1 - sigx10_1 *
      (0.5 * wup1 - 0.5 * mud1)/psq1^2) * ewu1 - ((((0.5 * (((1 - 0.5 * (ewu1/sigma_sq1)) *
      pmu1 - 0.5 * (mu1 * dmu1/ewusr1)) * ewu1/sqrt(sigma_sq1)) - 0.5 * (mu1 *
      ((0.5 * (mu1^2/ewusr1^2) - 0.5) * sqrt(sigma_sq1) + 0.5 * (ewu1/sqrt(sigma_sq1))) *
      dmu1/ewusr1)) * depsi1 + 0.5 * ((0.5 * wup1 - 0.5 * mud1) * depsi1 *
      ewu1 * mupsi1^2/sigma_sq1^2)) * pmusig1 - (sigx8_1 * dmusig1 * ewu1 +
      2 * ((0.5 * wup1 - 0.5 * mud1) * pmusig1 * pmu1 * sqrt(sigma_sq1)/psq1^2)) *
      (0.5 * wup1 - 0.5 * mud1) * depsi1)/psq1^2 + sigx11_1^2 * (1 - prZ)/sigx2)) *
      (1 - prZ)/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx11_1 * sigx11_2 * prZ * (1 - prZ)/sigx2^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((((sigx8_1 * depsi1 * musi1/sigmastar1 + 0.5 *
      (depsi1 * mupsi1^2/sigma_sq1)) * sigx13_1/sigma_sq1 - ((0.5 * (usq1 *
      ewv1/sigma_sq1) + 0.5 * ((ewu1/sigma_sq1 - 1) * ewv1/sigma_sq1 + 1 -
      0.5 * (usq1 * vsq1))) * musi1/sigmastar1 + mu1 * sigx7_1 - sigx12_1 *
      (2 * (sigx7_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2) + S * (epsilon))) *
      depsi1/starsq1^2) * dmusig1 + 0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) *
      pmusig1/sigma_sq1 - sigx8_1 * dmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)) *
      ewu1 + 0.5 * (mu1 * sigx15_1 * dmu1/(ewusr1 * pmu1)))/pmu1 - 0.5 * ((sigx10_1 *
      ewu1 * pmu1 - (0.5 * (mu1 * dmu1/ewusr1) + 2 * ((0.5 * wup1 - 0.5 * mud1) *
      pmu1^2 * sqrt(sigma_sq1)/psq1^2)) * depsi1 * pmusig1)/psq1^2))/sigx17_1 -
      (sigx11_1 * (1 - prZ) * sqrt(sigma_sq1) + 0.5 * (sigx2 * ewu1/sqrt(sigma_sq1))) *
        sigx16_1/sigx17_1^2) * (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx11_1 *
    sigx16_2 * prZ * (1 - prZ) * ewv2 * sqrt(sigma_sq2)/sigx17_2^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 1):(nXvar + 2 * nmuZUvar + nuZUvar), (nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    sigx11_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((((0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) *
      pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * ewu2) + 0.5 * pmusig2) * depsi2 *
      mupsi2^2/sigma_sq2^2 - ((sigx8_2^2 * ewu2 * musi2/starsq2 + ((0.5 * (ewu2/sigma_sq2) -
      0.5 * (0.5 * usq2 + ewu2/sigma_sq2)) * usq2 * ewv2 * musi2/sigmastar2 -
      sigx7_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
        2 * (S * (epsilon))) * ewu2)/starsq2^2) * depsi2 + sigx8_2 * (0.5 *
      (depsi2 * ewu2 * mupsi2^2/sigma_sq2^2) + depsi2)) * dmusig2)/psq2 - sigx10_2 *
      (0.5 * wup2 - 0.5 * mud2)/psq2^2) * ewu2 - ((((0.5 * (((1 - 0.5 * (ewu2/sigma_sq2)) *
      pmu2 - 0.5 * (mu2 * dmu2/ewusr2)) * ewu2/sqrt(sigma_sq2)) - 0.5 * (mu2 *
      ((0.5 * (mu2^2/ewusr2^2) - 0.5) * sqrt(sigma_sq2) + 0.5 * (ewu2/sqrt(sigma_sq2))) *
      dmu2/ewusr2)) * depsi2 + 0.5 * ((0.5 * wup2 - 0.5 * mud2) * depsi2 *
      ewu2 * mupsi2^2/sigma_sq2^2)) * pmusig2 - (sigx8_2 * dmusig2 * ewu2 +
      2 * ((0.5 * wup2 - 0.5 * mud2) * pmusig2 * pmu2 * sqrt(sigma_sq2)/psq2^2)) *
      (0.5 * wup2 - 0.5 * mud2) * depsi2)/psq2^2 + sigx11_2^2 * prZ/sigx2)) *
      prZ/sigx2, FUN = "*"), uHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -wHvar * (sigx16_1 *
    sigx11_2 * prZ * (1 - prZ) * ewv1 * sqrt(sigma_sq1)/sigx17_1^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((((((sigx8_2 * depsi2 * musi2/sigmastar2 + 0.5 * (depsi2 * mupsi2^2/sigma_sq2)) *
      sigx13_2/sigma_sq2 - ((0.5 * (usq2 * ewv2/sigma_sq2) + 0.5 * ((ewu2/sigma_sq2 -
      1) * ewv2/sigma_sq2 + 1 - 0.5 * (usq2 * vsq2))) * musi2/sigmastar2 +
      mu2 * sigx7_2 - sigx12_2 * (2 * (sigx7_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2) +
      S * (epsilon))) * depsi2/starsq2^2) * dmusig2 + 0.5 * (((0.5 * (mupsi2^2/sigma_sq2) -
      2) * pmusig2/sigma_sq2 - sigx8_2 * dmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)) *
      ewu2 + 0.5 * (mu2 * sigx15_2 * dmu2/(ewusr2 * pmu2)))/pmu2 - 0.5 * ((sigx10_2 *
      ewu2 * pmu2 - (0.5 * (mu2 * dmu2/ewusr2) + 2 * ((0.5 * wup2 - 0.5 * mud2) *
      pmu2^2 * sqrt(sigma_sq2)/psq2^2)) * depsi2 * pmusig2)/psq2^2))/sigx17_2 -
      (sigx11_2 * prZ * sqrt(sigma_sq2) + 0.5 * (sigx2 * ewu2/sqrt(sigma_sq2))) *
        sigx16_2/sigx17_2^2) * prZ * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar),
    (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
      2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx11_2 * (sigx9 * prZ/sigx2 + 1) * prZ * ewz/sigx2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
    (((((0.5 * (depsi1 * mupsi1^2/sigma_sq1) - depsi1 * musi1 * sigx13_1/sigmastar1) *
      ewv1 * sigx13_1/sigma_sq1 + depsi1 * (mu1/starsq1 - (((3 * (mu1) - 2 *
      (sigx12_1 * sigma_sq1 * musi1 * sigmastar1/starsq1^2)) * ewv1 - S * ewu1 *
      (epsilon)) * sigx12_1 + (0.5 * (ewv1/sigma_sq1) - 0.5 * (0.5 * vsq1 +
      ewv1/sigma_sq1)) * vsq1 * ewu1 * musi1/sigmastar1)/starsq1^2)) * dmusig1 +
      (0.5 * (((0.5 * (mupsi1^2/sigma_sq1) - 2) * pmusig1/sigma_sq1 + dmusig1 *
        sigx13_1) * ewv1) + 0.5 * pmusig1) * depsi1 * mupsi1^2/sigma_sq1^2)/pmu1 -
      (0.5 * (((dmusig1 * sigx13_1 - pmusig1 * pmu1^2/psq1^2) * depsi1 + 0.5 *
        sigx6_1) * ewv1) + 0.5 * (depsi1 * pmusig1)) * pmu1/psq1^2)/sigx17_1 -
      (sigx16_1 * (1 - prZ) + 0.5 * (sigx2/sqrt(sigma_sq1))) * sigx16_1 * ewv1/sigx17_1^2) *
    (1 - prZ) * ewv1, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 *
    nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (sigx16_1 * sigx16_2 * prZ * (1 - prZ) * ewv1 * ewv2 * sqrt(sigma_sq2)/(sigx17_2^2 *
      sqrt(sigma_sq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar +
    2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * sigx16_1 * (1 - (1 - prZ) * sigx9/sigx2) * prZ *
      ewv1 * ewz/(sigx2 * sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (((((0.5 * (depsi2 * mupsi2^2/sigma_sq2) - depsi2 *
      musi2 * sigx13_2/sigmastar2) * ewv2 * sigx13_2/sigma_sq2 + depsi2 * (mu2/starsq2 -
      (((3 * (mu2) - 2 * (sigx12_2 * sigma_sq2 * musi2 * sigmastar2/starsq2^2)) *
        ewv2 - S * ewu2 * (epsilon)) * sigx12_2 + (0.5 * (ewv2/sigma_sq2) -
        0.5 * (0.5 * vsq2 + ewv2/sigma_sq2)) * vsq2 * ewu2 * musi2/sigmastar2)/starsq2^2)) *
      dmusig2 + (0.5 * (((0.5 * (mupsi2^2/sigma_sq2) - 2) * pmusig2/sigma_sq2 +
      dmusig2 * sigx13_2) * ewv2) + 0.5 * pmusig2) * depsi2 * mupsi2^2/sigma_sq2^2)/pmu2 -
      (0.5 * (((dmusig2 * sigx13_2 - pmusig2 * pmu2^2/psq2^2) * depsi2 + 0.5 *
        sigx6_2) * ewv2) + 0.5 * (depsi2 * pmusig2)) * pmu2/psq2^2)/sigx17_2 -
      (sigx16_2 * prZ + 0.5 * (sigx2/sqrt(sigma_sq2))) * sigx16_2 * ewv2/sigx17_2^2) *
      prZ * ewv2, FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar +
    1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (sigx16_2 * (sigx9 * prZ/sigx2 + 1) * prZ *
      ewv2 * ewz/(sigx2 * sqrt(sigma_sq2))), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar), (nXvar + 2 * nmuZUvar + 2 * nuZUvar +
    2 * nvZVvar + 1):(nXvar + 2 * nmuZUvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx9 * prZ/sigx2 + 1) * ewz) * sigx9 *
      prZ * ewz/sigx2, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf truncatednormal-normal distribution
#' @param start starting value for optimization
#' @param randStart if random starting values should be used
#' @param sdStart std. Error for random draws for starting values
#' @param olsParam OLS coefficients
#' @param dataTable dataframe contains id of observations
#' @param nXvar number of main variables (inputs + env. var)
#' @param nmuZUvar number of Zmu variables
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param muHvar matrix of Zmu variables
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# same sigma_u

## logit specification class membership
cnsftruncnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsftruncnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsftruncnormlike_logit,
    grad = cgradcnsftruncnormlike_logit, hess = chesscnsftruncnormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chesscnsftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chesscnsftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesscnsftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftruncnormlike_logit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chesscnsftruncnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftruncnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftruncnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftruncnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## cauchit specification class membership
cnsftruncnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsftruncnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsftruncnormlike_cauchit,
    grad = cgradcnsftruncnormlike_cauchit, hess = chesscnsftruncnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chesscnsftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chesscnsftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesscnsftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftruncnormlike_cauchit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chesscnsftruncnormlike_cauchit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftruncnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftruncnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftruncnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## probit specification class membership
cnsftruncnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsftruncnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsftruncnormlike_probit,
    grad = cgradcnsftruncnormlike_probit, hess = chesscnsftruncnormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chesscnsftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_probit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chesscnsftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesscnsftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftruncnormlike_probit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chesscnsftruncnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftruncnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftruncnormlike_probit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftruncnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## cloglog specification class membership
cnsftruncnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(ccnsftruncnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(ccnsftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsftruncnormlike_cloglog,
    grad = cgradcnsftruncnormlike_cloglog, hess = chesscnsftruncnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chesscnsftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(ccnsftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chesscnsftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(ccnsftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chesscnsftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftruncnormlike_cloglog(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chesscnsftruncnormlike_cloglog(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftruncnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftruncnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftruncnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

# Different sigma_u

## logit specification class membership
mcesftruncnormAlgOpt_logit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesftruncnormlike_logit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesftruncnormlike_logit,
    grad = cgradmcesftruncnormlike_logit, hess = chessmcesftruncnormlike_logit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessmcesftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_logit(parm, nXvar = nXvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessmcesftruncnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessmcesftruncnormlike_logit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftruncnormlike_logit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessmcesftruncnormlike_logit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftruncnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftruncnormlike_logit(parm = mlParam, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftruncnormlike_logit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## cauchit specification class membership
mcesftruncnormAlgOpt_cauchit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesftruncnormlike_cauchit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesftruncnormlike_cauchit,
    grad = cgradmcesftruncnormlike_cauchit, hess = chessmcesftruncnormlike_cauchit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessmcesftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cauchit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessmcesftruncnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessmcesftruncnormlike_cauchit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftruncnormlike_cauchit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessmcesftruncnormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftruncnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftruncnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftruncnormlike_cauchit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## probit specification class membership
mcesftruncnormAlgOpt_probit <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesftruncnormlike_probit(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesftruncnormlike_probit,
    grad = cgradmcesftruncnormlike_probit, hess = chessmcesftruncnormlike_probit,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessmcesftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_probit(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessmcesftruncnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessmcesftruncnormlike_probit(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftruncnormlike_probit(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessmcesftruncnormlike_probit(parm = mleObj$par, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftruncnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftruncnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftruncnormlike_probit(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

## cloglog specification class membership
mcesftruncnormAlgOpt_cloglog <- function(start, randStart, sdStart, olsParam, dataTable,
  S, wHvar, nXvar, muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, whichStart, initIter, initAlg,
  tol, gradtol, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart, initIter = initIter,
    initAlg = initAlg, printInfo = printInfo, tol = tol)
  initTrunc <- start_st$initTrunc
  startVal <- start_st$StartVal
  if (randStart)
    startVal <- startVal + rnorm(length(startVal), sd = sdStart)
  startLoglik <- sum(cmcesftruncnormlike_cloglog(startVal, nXvar = nXvar, nuZUvar = nuZUvar,
    nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...), bhhh = function(...) maxLik::maxBHHH(...),
      nr = function(...) maxLik::maxNR(...), nm = function(...) maxLik::maxNM(...),
      cg = function(...) maxLik::maxCG(...), sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal, fn = function(parm) -sum(cmcesftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, control = list(trace = printInfo,
    maxeval = itermax, stepmax = stepmax, xtol = tol, grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cmcesftruncnormlike_cloglog,
    grad = cgradmcesftruncnormlike_cloglog, hess = chessmcesftruncnormlike_cloglog,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0, iterlim = itermax, reltol = tol,
      tol = tol, qac = qac), nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    muHvar = muHvar, nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessmcesftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
      report.precision = 1L, preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(cmcesftruncnormlike_cloglog(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessmcesftruncnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(cmcesftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessmcesftruncnormlike_cloglog(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, trace = printInfo,
    eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftruncnormlike_cloglog(mleObj$par, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
      Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessmcesftruncnormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftruncnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar,
        nmuZUvar = nmuZUvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftruncnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftruncnormlike_cloglog(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, wHvar = wHvar, muHvar = muHvar, nmuZUvar = nmuZUvar,
    Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, mleObj = mleObj,
    mlParam = mlParam, initTrunc = initTrunc))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf truncatednormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsftruncnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## cauchit specification class membership
ccnsftruncnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## probit specification class membership
ccnsftruncnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## cloglog specification class membership
ccnsftruncnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

# different sigma_u

## logit specification class membership
cmcesftruncnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## cauchit specification class membership
cmcesftruncnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## probit specification class membership
cmcesftruncnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

## cloglog specification class membership
cmcesftruncnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- exp(-mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 -
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_c2 <- exp(-mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 -
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(mustar1 + 1/2 * sigmastar1^2) * pnorm(mustar1/sigmastar1 +
      sigmastar1)/pnorm(mustar1/sigmastar1)
    teBC_reciprocal_c2 <- exp(mustar2 + 1/2 * sigmastar2^2) * pnorm(mustar2/sigmastar2 +
      sigmastar2)/pnorm(mustar2/sigmastar2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1, teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1, NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2, NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      teJLMS_c = teJLMS_c, teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c, u_c = u_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
}

# Marginal effects on inefficiencies ----------
#' efficiencies for cnsf truncatednormal-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfmargtruncnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1 - Lambda *
    dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2, ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2)/2 *
    ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) + Lambda * (dnorm(Lambda)/pnorm(Lambda))^2),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtruncnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu/2) *
    dnorm(Lambda)/pnorm(Lambda) * (m1^2 - m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) * (Lambda + Lambda^3 + (2 + 3 * Lambda^2) *
      dnorm(Lambda)/pnorm(Lambda) + 2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
ccnsfmargtruncnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1 - Lambda *
    dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2, ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2)/2 *
    ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) + Lambda * (dnorm(Lambda)/pnorm(Lambda))^2),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtruncnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu/2) *
    dnorm(Lambda)/pnorm(Lambda) * (m1^2 - m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) * (Lambda + Lambda^3 + (2 + 3 * Lambda^2) *
      dnorm(Lambda)/pnorm(Lambda) + 2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
ccnsfmargtruncnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1 - Lambda *
    dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2, ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2)/2 *
    ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) + Lambda * (dnorm(Lambda)/pnorm(Lambda))^2),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtruncnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu/2) *
    dnorm(Lambda)/pnorm(Lambda) * (m1^2 - m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) * (Lambda + Lambda^3 + (2 + 3 * Lambda^2) *
      dnorm(Lambda)/pnorm(Lambda) + 2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
ccnsfmargtruncnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1 - Lambda *
    dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2, ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu/2)/2 *
    ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) + Lambda * (dnorm(Lambda)/pnorm(Lambda))^2),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtruncnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 1):(object$nXvar +
    object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + 1):(object$nXvar + object$nmuZUvar + object$nuZUvar + 2 *
    object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu * exp(Wv1) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv1))
  mustar2 <- (mu * exp(Wv2) - exp(Wu) * object$S * epsilon)/(exp(Wu) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu) * exp(Wv1)/(exp(Wu) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu) * exp(Wv2)/(exp(Wu) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu) + exp(Wv1)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu/sqrt(exp(Wu)))
  Pi2 <- 1/sqrt(exp(Wu) + exp(Wv2)) * dnorm((mu + object$S * epsilon)/sqrt(exp(Wu) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu/sqrt(exp(Wu)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu/2) *
    dnorm(Lambda)/pnorm(Lambda) * (m1^2 - m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), matrix(exp(Wu) *
    (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) * (Lambda + Lambda^3 + (2 + 3 * Lambda^2) *
      dnorm(Lambda)/pnorm(Lambda) + 2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
    ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff1 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff2 <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu], mu_mat[, !idTRUE_mu],
    Wu_mat[, !idTRUE_Wu])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c1")
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c2")
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu], colnames(muHvar)[-1][!idTRUE_mu],
    colnames(uHvar)[-1][!idTRUE_Wu]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# different sigma_u

## logit specification class membership
cmcesfmargtruncnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2,
    ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2,
    ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2)/2 *
    ((1 + Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2)/2 *
    ((1 + Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtruncnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  m1_1 <- exp(Wu1/2) * (Lambda1 + dnorm(Lambda1)/pnorm(Lambda1))
  m2_1 <- exp(Wu1) * (1 - Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2)
  m1_2 <- exp(Wu2/2) * (Lambda2 + dnorm(Lambda2)/pnorm(Lambda2))
  m2_2 <- exp(Wu2) * (1 - Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu1/2) *
    dnorm(Lambda1)/pnorm(Lambda1) * (m1_1^2 - m2_1), ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu2/2) *
    dnorm(Lambda2)/pnorm(Lambda2) * (m1_2^2 - m2_2), ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - 1/2 * dnorm(Lambda1)/pnorm(Lambda1) * (Lambda1 + Lambda1^3 + (2 + 3 *
      Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + 2 * Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2)),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - 1/2 * dnorm(Lambda2)/pnorm(Lambda2) * (Lambda2 + Lambda2^3 + (2 + 3 *
      Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + 2 * Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2)),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmcesfmargtruncnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2,
    ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2,
    ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2)/2 *
    ((1 + Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2)/2 *
    ((1 + Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtruncnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  m1_1 <- exp(Wu1/2) * (Lambda1 + dnorm(Lambda1)/pnorm(Lambda1))
  m2_1 <- exp(Wu1) * (1 - Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2)
  m1_2 <- exp(Wu2/2) * (Lambda2 + dnorm(Lambda2)/pnorm(Lambda2))
  m2_2 <- exp(Wu2) * (1 - Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu1/2) *
    dnorm(Lambda1)/pnorm(Lambda1) * (m1_1^2 - m2_1), ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu2/2) *
    dnorm(Lambda2)/pnorm(Lambda2) * (m1_2^2 - m2_2), ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - 1/2 * dnorm(Lambda1)/pnorm(Lambda1) * (Lambda1 + Lambda1^3 + (2 + 3 *
      Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + 2 * Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2)),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - 1/2 * dnorm(Lambda2)/pnorm(Lambda2) * (Lambda2 + Lambda2^3 + (2 + 3 *
      Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + 2 * Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2)),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmcesfmargtruncnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2,
    ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2,
    ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2)/2 *
    ((1 + Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2)/2 *
    ((1 + Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtruncnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  m1_1 <- exp(Wu1/2) * (Lambda1 + dnorm(Lambda1)/pnorm(Lambda1))
  m2_1 <- exp(Wu1) * (1 - Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2)
  m1_2 <- exp(Wu2/2) * (Lambda2 + dnorm(Lambda2)/pnorm(Lambda2))
  m2_2 <- exp(Wu2) * (1 - Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu1/2) *
    dnorm(Lambda1)/pnorm(Lambda1) * (m1_1^2 - m2_1), ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu2/2) *
    dnorm(Lambda2)/pnorm(Lambda2) * (m1_2^2 - m2_2), ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - 1/2 * dnorm(Lambda1)/pnorm(Lambda1) * (Lambda1 + Lambda1^3 + (2 + 3 *
      Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + 2 * Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2)),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - 1/2 * dnorm(Lambda2)/pnorm(Lambda2) * (Lambda2 + Lambda2^3 + (2 + 3 *
      Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + 2 * Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2)),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmcesfmargtruncnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2,
    ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1 -
    Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2,
    ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1/2)/2 *
    ((1 + Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2/2)/2 *
    ((1 + Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c1")
  colnames(margEff2) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]), "_c2")
  colnames(margEff_c) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]), "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtruncnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + object$nmuZUvar)]
  omega2 <- object$mlParam[(object$nXvar + object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar)]
  delta1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 1):(object$nXvar +
    2 * object$nmuZUvar + object$nuZUvar)]
  delta2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar)]
  theta <- object$mlParam[(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1):(object$nXvar + 2 * object$nmuZUvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable, rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable, rhs = 4)
  Zvar <- model.matrix(object$formula, data = object$dataTable, rhs = 5)
  mu1 <- as.numeric(crossprod(matrix(omega1), t(muHvar)))
  mu2 <- as.numeric(crossprod(matrix(omega2), t(muHvar)))
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar1 <- (mu1 * exp(Wv1) - exp(Wu1) * object$S * epsilon)/(exp(Wu1) + exp(Wv1))
  mustar2 <- (mu2 * exp(Wv2) - exp(Wu2) * object$S * epsilon)/(exp(Wu2) + exp(Wv2))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 1/sqrt(exp(Wu1) + exp(Wv1)) * dnorm((mu1 + object$S * epsilon)/sqrt(exp(Wu1) +
    exp(Wv1))) * pnorm(mustar1/sigmastar1)/pnorm(mu1/sqrt(exp(Wu1)))
  Pi2 <- 1/sqrt(exp(Wu2) + exp(Wv2)) * dnorm((mu2 + object$S * epsilon)/sqrt(exp(Wu2) +
    exp(Wv2))) * pnorm(mustar2/sigmastar2)/pnorm(mu2/sqrt(exp(Wu2)))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  Lambda1 <- mu1/exp(Wu1/2)
  Lambda2 <- mu2/exp(Wu2/2)
  m1_1 <- exp(Wu1/2) * (Lambda1 + dnorm(Lambda1)/pnorm(Lambda1))
  m2_1 <- exp(Wu1) * (1 - Lambda1 * dnorm(Lambda1)/pnorm(Lambda1) - (dnorm(Lambda1)/pnorm(Lambda1))^2)
  m1_2 <- exp(Wu2/2) * (Lambda2 + dnorm(Lambda2)/pnorm(Lambda2))
  m2_2 <- exp(Wu2) * (1 - Lambda2 * dnorm(Lambda2)/pnorm(Lambda2) - (dnorm(Lambda2)/pnorm(Lambda2))^2)
  mu_mat1 <- kronecker(matrix(omega1[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu1/2) *
    dnorm(Lambda1)/pnorm(Lambda1) * (m1_1^2 - m2_1), ncol = 1))
  mu_mat2 <- kronecker(matrix(omega2[2:object$nmuZUvar], nrow = 1), matrix(1/exp(Wu2/2) *
    dnorm(Lambda2)/pnorm(Lambda2) * (m1_2^2 - m2_2), ncol = 1))
  Wu_mat1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1), matrix(exp(Wu1) *
    (1 - 1/2 * dnorm(Lambda1)/pnorm(Lambda1) * (Lambda1 + Lambda1^3 + (2 + 3 *
      Lambda1^2) * dnorm(Lambda1)/pnorm(Lambda1) + 2 * Lambda1 * (dnorm(Lambda1)/pnorm(Lambda1))^2)),
    ncol = 1))
  Wu_mat2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1), matrix(exp(Wu2) *
    (1 - 1/2 * dnorm(Lambda2)/pnorm(Lambda2) * (Lambda2 + Lambda2^3 + (2 + 3 *
      Lambda2^2) * dnorm(Lambda2)/pnorm(Lambda2) + 2 * Lambda2 * (dnorm(Lambda2)/pnorm(Lambda2))^2)),
    ncol = 1))
  idTRUE_mu1 <- substring(names(omega1)[-1], 5) %in% substring(names(delta1)[-1],
    4)
  idTRUE_Wu1 <- substring(names(delta1)[-1], 4) %in% substring(names(omega1)[-1],
    5)
  idTRUE_mu2 <- substring(names(omega2)[-1], 5) %in% substring(names(delta2)[-1],
    4)
  idTRUE_Wu2 <- substring(names(delta2)[-1], 4) %in% substring(names(omega2)[-1],
    5)
  margEff1 <- cbind(mu_mat1[, idTRUE_mu1] + Wu_mat1[, idTRUE_Wu1], mu_mat1[, !idTRUE_mu1],
    Wu_mat1[, !idTRUE_Wu1])
  margEff2 <- cbind(mu_mat2[, idTRUE_mu2] + Wu_mat2[, idTRUE_Wu2], mu_mat2[, !idTRUE_mu2],
    Wu_mat2[, !idTRUE_Wu2])
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  colnames(margEff2) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu2], colnames(muHvar)[-1][!idTRUE_mu2],
    colnames(uHvar)[-1][!idTRUE_Wu2]))
  colnames(margEff_c) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu1], colnames(muHvar)[-1][!idTRUE_mu1],
    colnames(uHvar)[-1][!idTRUE_Wu1]))
  return(data.frame(margEff1, margEff2, margEff_c))
}
