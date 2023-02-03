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
# Convolution: truncated skewed laplace - normal                               #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf tsl-normal distribution
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

## logit specification class membership
ccnsftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
ccnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
ccnsftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
ccnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# different sigma_u

## logit specification class membership
cmcesftslnormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + S * epsilon *
    (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + S * epsilon *
    (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cauchit specification class membership
cmcesftslnormlike_cauchit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + S * epsilon *
    (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + S * epsilon *
    (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## probit specification class membership
cmcesftslnormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + S * epsilon *
    (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + S * epsilon *
    (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

## cloglog specification class membership
cmcesftslnormlike_cloglog <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + S * epsilon *
    (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + S * epsilon *
    (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), return(wHvar *
    log(Probc1 * Pi1 + Probc2 * Pi2)))
}

# starting value for the log-likelihood ----------
#' starting values for cnsf tsl-normal distribution
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
# same sigma_u
cstcnsftslnorm <- function(olsObj, epsiRes, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, whichStart,
  initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csttslnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTSL <- NULL
  } else {
    cat("Initialization: SFA + truncated skewed laplace - normal distributions...\n")
    initTSL <- maxLik::maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradtslnormlike, method = initAlg,
      control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar)
    Esti <- initTSL$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), "lambda", paste0("CN_", colnames(Zvar)))
  return(list(StartVal = StartVal, initTSL = initTSL))
}

# different sigma_u
cstmcesftslnorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- csttslnorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initTSL <- NULL
  } else {
    cat("Initialization: SFA + truncated skewed laplace - normal distributions...\n")
    initTSL <- maxLik::maxLik(logLik = ctslnormlike, start = csttslnorm(olsObj = olsObj,
      epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE]), grad = cgradtslnormlike, method = initAlg,
      control = list(iterlim = initIter, printLevel = if (printInfo) 2 else 0,
        reltol = tol), nXvar = nXvar, nuZUvar = 1, uHvar = uHvar[,
        1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[,
        1, drop = FALSE], Yvar = Yvar, Xvar = Xvar, S = S,
      wHvar = wHvar)
    Esti <- initTSL$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 3], 1.05 *
    Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), "lambda",
    "lambda", paste0("CME_", colnames(Zvar)))
  return(list(StartVal = StartVal, initTSL = initTSL))
}


# Gradient of the likelihood function ----------
#' gradient for cnsf tsl-normal distribution
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

## logit specification class membership
cgradcnsftslnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  epsivu2 <- ((1 + lambda) * ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewu_h) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewu_h) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewv1_h - (1 +
    lambda) * pepsi1/ewu_h) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewv2_h - (1 +
    lambda) * pepsi2/ewu_h) * sigx3_2)
  sigx5 <- (prC * sigx4_2 + sigx4_1 * ewz/wzdeno)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (prC * sigx6_2 + sigx6_1 * ewz/wzdeno)
  sigx8_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx8_2 <- (dmusig2 * ewv2_h/ewu_h)
  sigx9 <- (S * (epsilon)/ewu_h)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewv1_h/ewu_h) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewv2_h/ewu_h) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx15_2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewv1_h/ewu_h) - 0.5 *
    (S * (epsilon)/ewv1_h))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewv2_h/ewu_h) - 0.5 *
    (S * (epsilon)/ewv2_h))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  luh <- (1 + 2 * (lambda * ewu_h))
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewu_h) *
    pepsi1 - depsi1 * ewv1_h/ewu_h)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewu_h) *
    pepsi2 - depsi2 * ewv2_h/ewu_h)
  sigx20 <- (sigx19_1 * sigx3_1 * ewz/wzdeno + sigx19_2 * prC *
    sigx3_2)
  sigx21 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((prC *
    sigx14_2 + sigx14_1 * ewz/wzdeno)/sigx7 - lambda * ewu_h/luh),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (2 * (sigx2_1 *
    sigx16_1) - sigx18_1 * sigx3_1) * ewz/(sigx7 * wzdeno),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = prC * (2 *
    (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)/sigx7, FUN = "*"),
    1/(1 + lambda) - (sigx20/sigx7 + 2 * (ewu_h/luh)), sweep(Zvar,
      MARGIN = 1, STATS = prC * sigx21 * ewz/(sigx7 * wzdeno),
      FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- (ewz2 * sigx4_2 + ewz1 * sigx4_1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (ewz2 * sigx6_2 + ewz1 * sigx6_1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- (ewz2 * sigx14_2 + ewz1 * sigx14_1)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * ewz1 * sigx3_1 + sigx19_2 * ewz2 *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  sigx24 <- (pi * sigx7 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx20/sigx7 -
    lambda * ewusr/luh), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = ewz1 * sigx21_1/sigx7, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = ewz2 * sigx21_2/sigx7, FUN = "*"),
    1/(1 + lambda) - (sigx22/sigx7 + 2 * (ewusr/luh)), sweep(Zvar,
      MARGIN = 1, STATS = sigx23/sigx24, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsftslnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- ((1 - pwZ) * sigx4_2 + sigx4_1 * pwZ)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - pwZ) * sigx6_2 + sigx6_1 * pwZ)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- ((1 - pwZ) * sigx14_2 + sigx14_1 * pwZ)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * sigx3_1 * pwZ + sigx19_2 * (1 - pwZ) *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx20/sigx7 -
    lambda * ewusr/luh), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = pwZ * sigx21_1/sigx7, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = (1 - pwZ) * sigx21_2/sigx7, FUN = "*"),
    1/(1 + lambda) - (sigx22/sigx7 + 2 * (ewusr/luh)), sweep(Zvar,
      MARGIN = 1, STATS = sigx23 * dwZ/sigx7, FUN = "*"))

  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}


## cloglog specification class membership
cgradcnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- ((1 - prZ) * sigx4_1 + sigx4_2 * prZ)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - prZ) * sigx6_1 + sigx6_2 * prZ)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- ((1 - prZ) * sigx14_1 + sigx14_2 * prZ)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * (1 - prZ) * sigx3_1 + sigx19_2 * prZ *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (sigx20/sigx7 -
    lambda * ewusr/luh), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = (1 - prZ) * sigx21_1/sigx7, FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = prZ * sigx21_2/sigx7, FUN = "*"),
    1/(1 + lambda) - (sigx22/sigx7 + 2 * (ewusr/luh)), sweep(Zvar,
      MARGIN = 1, STATS = sigx23 * prZ * ewz/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesftslnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- ((1 + lambda1) * ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  epsi2 <- ((1 + lambda2) * ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv1_h - (1 + lambda1) * pepsi1/ewu1_h)
  sigx2_2 <- (depsi2/ewv2_h - (1 + lambda2) * pepsi2/ewu2_h)
  sigx3_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h) *
    (1 + lambda1))
  sigx4_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h) *
    (1 + lambda2))
  sigx5_1 <- (2 * (sigx1_1 * sigx3_1) - sigx2_1 * sigx4_1)
  sigx5_2 <- (2 * (sigx1_2 * sigx3_2) - sigx2_2 * sigx4_2)
  lau1 <- (1 + 2 * (lambda1 * ewu1_h))
  lau2 <- (1 + 2 * (lambda2 * ewu2_h))
  sigx6 <- (prC * (1 + lambda2) * sigx5_2/lau2 + (1 + lambda1) *
    sigx5_1 * ewz/(lau1 * wzdeno))
  sigx7_1 <- (2 * (sigx3_1 * pmusig1) - sigx4_1 * pepsi1)
  sigx7_2 <- (2 * (sigx3_2 * pmusig2) - sigx4_2 * pepsi2)
  sigx8 <- (prC * (1 + lambda2) * sigx7_2/lau2 + (1 + lambda1) *
    sigx7_1 * ewz/(lau1 * wzdeno))
  sigx9_1 <- (dmusig1 * ewv1_h/ewu1_h)
  sigx9_2 <- (dmusig2 * ewv2_h/ewu2_h)
  siu1 <- (S * (epsilon)/ewu1_h)
  siu2 <- (S * (epsilon)/ewu2_h)
  uv1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  uv2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx10_1 <- ((0.5 * sigx9_1 - (0.5 * siu1 + 2 * uv1) * pmusig1) *
    sigx3_1)
  sigx10_2 <- ((0.5 * sigx9_2 - (0.5 * siu2 + 2 * uv2) * pmusig2) *
    sigx3_2)
  siv1 <- (depsi1 * ewv1_h/ewu1_h)
  siv2 <- (depsi2 * ewv2_h/ewu2_h)
  sigx11_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx11_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  s12v1 <- (0.5 * siv1 - (0.5 * siu1 + 2 * sigx11_1) * pepsi1) *
    (1 + lambda1)
  s12v2 <- (0.5 * siv2 - (0.5 * siu2 + 2 * sigx11_2) * pepsi2) *
    (1 + lambda2)
  sigx12_1 <- (2 * sigx10_1 - s12v1 * sigx4_1)
  sigx12_2 <- (2 * sigx10_2 - (s12v2 * sigx4_2 + lambda2 *
    sigx7_2 * ewu2_h/lau2))
  sigx13_1 <- (sigx12_1/(lau1 * wzdeno) - lambda1 * wzdeno *
    sigx7_1 * ewu1_h/(lau1 * wzdeno)^2)
  sigx14_1 <- (0.5 * (ewv1_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx14_2 <- (0.5 * (ewv2_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv2_h))
  s15pd1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx14_1 * dmusig1)
  s15pd2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx14_2 * dmusig2)
  sigx15_1 <- (sigx3_1 * s15pd1)
  sigx15_2 <- (sigx3_2 * s15pd2)
  sigx16_1 <- ((1 + lambda1) * ewv1_h/ewu1_h)
  sigx16_2 <- ((1 + lambda2) * ewv2_h/ewu2_h)
  sigx17_1 <- (0.5 * sigx16_1 - 0.5 * (S * (epsilon)/ewv1_h))
  sigx17_2 <- (0.5 * sigx16_2 - 0.5 * (S * (epsilon)/ewv2_h))
  s18pd1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) - sigx17_1 *
    depsi1)
  s18pd2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) - sigx17_2 *
    depsi2)
  sigx18_1 <- (2 * sigx15_1 - s18pd1 * sigx4_1)
  sigx18_2 <- (2 * sigx15_2 - s18pd2 * sigx4_2)
  sigx19_1 <- ((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewu1_h)
  sigx19_2 <- ((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewu2_h)
  sigx20_1 <- (sigx19_1 * pepsi1 - depsi1 * ewv1_h/ewu1_h)
  sigx20_2 <- (sigx19_2 * pepsi2 - depsi2 * ewv2_h/ewu2_h)
  sigx21_1 <- (2 * (sigx3_1 * pmusig1) - (sigx20_1 * (1 + lambda1) +
    pepsi1) * sigx4_1)
  sigx21_2 <- ((sigx20_2 * (1 + lambda2) + pepsi2) * sigx4_2 +
    2 * ((1 + lambda2) * sigx7_2 * ewu2_h/lau2))
  sigx22_1 <- (sigx21_1/(lau1 * wzdeno) - 2 * (wzdeno * (1 +
    lambda1) * sigx7_1 * ewu1_h/(lau1 * wzdeno)^2))
  sigx22_2 <- (2 * (sigx3_2 * pmusig2) - sigx21_2)
  lauz1 <- (1/(lau1 * wzdeno) - lau1 * ewz/(lau1 * wzdeno)^2)
  sigx23 <- ((1 + lambda1) * lauz1 * sigx7_1 - prC * (1 + lambda2) *
    sigx7_2/(lau2 * wzdeno))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx6/sigx8,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = sigx13_1 *
    (1 + lambda1) * ewz/sigx8, FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = prC * (1 + lambda2) * sigx12_2/(sigx8 * lau2),
    FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (1 + lambda1) *
    sigx18_1 * ewz/(sigx8 * lau1 * wzdeno), FUN = "*"), sweep(vHvar,
    MARGIN = 1, STATS = prC * (1 + lambda2) * sigx18_2/(sigx8 *
      lau2), FUN = "*"), sigx22_1 * ewz/sigx8, prC * sigx22_2/(sigx8 *
    lau2), sweep(Zvar, MARGIN = 1, STATS = sigx23 * ewz/sigx8,
    FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesftslnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- (ewz2 * (1 + lambda2) * sigx4_2/luh2 + ewz1 * (1 +
    lambda1) * sigx4_1/luh1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (ewz2 * (1 + lambda2) * sigx6_2/luh2 + ewz1 * (1 +
    lambda1) * sigx6_1/luh1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  sigx24 <- (pi * sigx7 * ((Wz)^2 + 1))
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ewz1 * (1 +
    lambda1) * sigx20_1/(sigx7 * luh1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = ewz2 * (1 + lambda2) * sigx20_2/(sigx7 *
      luh2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = ewz1 *
    (1 + lambda1) * sigx21_1/(sigx7 * luh1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = ewz2 * (1 + lambda2) *
      sigx21_2/(sigx7 * luh2), FUN = "*"), ewz1 * sigx22_1/(sigx7 *
      luh1), ewz2 * sigx22_2/(sigx7 * luh2), sweep(Zvar,
      MARGIN = 1, STATS = sigx23/sigx24, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesftslnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- ((1 - pwZ) * (1 + lambda2) * sigx4_2/luh2 + (1 +
    lambda1) * sigx4_1 * pwZ/luh1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - pwZ) * (1 + lambda2) * sigx6_2/luh2 + (1 +
    lambda1) * sigx6_1 * pwZ/luh1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = pwZ * (1 +
    lambda1) * sigx20_1/(sigx7 * luh1), FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = (1 - pwZ) * (1 + lambda2) * sigx20_2/(sigx7 *
      luh2), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = pwZ *
    (1 + lambda1) * sigx21_1/(sigx7 * luh1), FUN = "*"),
    sweep(vHvar, MARGIN = 1, STATS = (1 - pwZ) * (1 + lambda2) *
      sigx21_2/(sigx7 * luh2), FUN = "*"), pwZ * sigx22_1/(sigx7 *
      luh1), (1 - pwZ) * sigx22_2/(sigx7 * luh2), sweep(Zvar,
      MARGIN = 1, STATS = sigx23 * dwZ/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesftslnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- ((1 - prZ) * (1 + lambda1) * sigx4_1/luh1 + (1 +
    lambda2) * sigx4_2 * prZ/luh2)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - prZ) * (1 + lambda1) * sigx6_1/luh1 + (1 +
    lambda2) * sigx6_2 * prZ/luh2)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * sigx5/sigx7,
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (1 - prZ) *
    (1 + lambda1) * sigx20_1/(sigx7 * luh1), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = prZ * (1 + lambda2) *
      sigx20_2/(sigx7 * luh2), FUN = "*"), sweep(vHvar,
      MARGIN = 1, STATS = (1 - prZ) * (1 + lambda1) * sigx21_1/(sigx7 *
        luh1), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = prZ *
      (1 + lambda2) * sigx21_2/(sigx7 * luh2), FUN = "*"),
    (1 - prZ) * sigx22_1/(sigx7 * luh1), prZ * sigx22_2/(sigx7 *
      luh2), sweep(Zvar, MARGIN = 1, STATS = sigx23 * prZ *
      ewz/sigx7, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' hessian for cnsf tsl-normal distribution
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

## logit specification class membership
chesscnsftslnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  musig1 <- (ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewv1_h/ewu_h + S * (epsilon)/ewv1_h)
  epsivu2 <- ((1 + lambda) * ewv2_h/ewu_h + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu_h)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewu_h)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewu_h) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewu_h) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewv1_h - (1 +
    lambda) * pepsi1/ewu_h) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewv2_h - (1 +
    lambda) * pepsi2/ewu_h) * sigx3_2)
  sigx5 <- (prC * sigx4_2 + sigx4_1 * ewz/wzdeno)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (prC * sigx6_2 + sigx6_1 * ewz/wzdeno)
  sigx8_1 <- (dmusig1 * ewv1_h/ewu_h)
  sigx8_2 <- (dmusig2 * ewv2_h/ewu_h)
  sigx9 <- (S * (epsilon)/ewu_h)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewv1_h/ewu_h) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewv2_h/ewu_h) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewv1_h/ewu_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx15_2 <- (0.5 * (ewv2_h/ewu_h) - 0.5 * (S * (epsilon)/ewv2_h))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewv1_h/ewu_h) - 0.5 *
    (S * (epsilon)/ewv1_h))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewv2_h/ewu_h) - 0.5 *
    (S * (epsilon)/ewv2_h))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  luh <- (1 + 2 * (lambda * ewu_h))
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewu_h) *
    pepsi1 - depsi1 * ewv1_h/ewu_h)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewu_h) *
    pepsi2 - depsi2 * ewv2_h/ewu_h)
  sigx20 <- (sigx19_1 * sigx3_1 * ewz/wzdeno + sigx19_2 * prC *
    sigx3_2)
  sigx21 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  sigx22_1 <- (depsi1/ewv1_h - (1 + lambda) * pepsi1/ewu_h)
  sigx22_2 <- (depsi2/ewv2_h - (1 + lambda) * pepsi2/ewu_h)
  lusq <- (1 - 8 * (ewu^2/(2 * ewu)^2))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * (2 * (((musig2/ewv2_h -
      1/ewu_h) * dmusig2/ewv2_h - sigx1_2/ewu_h) * sigx2_2) -
      ((epsivu2/ewv2_h - (1 + lambda)/ewu_h) * depsi2/ewv2_h -
        (1 + lambda) * sigx22_2/ewu_h) * sigx3_2) + (2 *
      (((musig1/ewv1_h - 1/ewu_h) * dmusig1/ewv1_h - sigx1_1/ewu_h) *
        sigx2_1) - ((epsivu1/ewv1_h - (1 + lambda)/ewu_h) *
      depsi1/ewv1_h - (1 + lambda) * sigx22_1/ewu_h) *
      sigx3_1) * ewz/wzdeno - sigx5^2/sigx7)/sigx7, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (prC * (2 * ((((0.5 +
      0.5 * sigx9 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewu_h +
      (0.5 * (musig2/ewu_h) - (0.5 * sigx9 + 2 * sigx10_2)/ewv2_h) *
        dmusig2) * sigx2_2) - ((0.5 * (epsivu2/ewu_h) -
      (0.5 * sigx9 + 2 * sigx12_2)/ewv2_h) * depsi2 + (0.5 *
      pepsi2 - sigx13_2 * (1 + lambda))/ewu_h) * (1 + lambda) *
      sigx3_2) + (2 * ((((0.5 + 0.5 * sigx9 + 2 * sigx10_1) *
      pmusig1 - 0.5 * sigx8_1)/ewu_h + (0.5 * (musig1/ewu_h) -
      (0.5 * sigx9 + 2 * sigx10_1)/ewv1_h) * dmusig1) *
      sigx2_1) - ((0.5 * (epsivu1/ewu_h) - (0.5 * sigx9 +
      2 * sigx12_1)/ewv1_h) * depsi1 + (0.5 * pepsi1 -
      sigx13_1 * (1 + lambda))/ewu_h) * (1 + lambda) *
      sigx3_1) * ewz/wzdeno - (prC * sigx14_2 + sigx14_1 *
      ewz/wzdeno) * sigx5/sigx7)/sigx7, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((2 * ((dmusig1 * (ewv1/(2 * ewu) - (sigx15_1 * musig1 +
    0.5))/ewv1_h - sigx16_1/ewu_h) * sigx2_1) - (((1 + lambda)^2 *
    ewv1/(2 * ewu) - (epsivu1 * sigx17_1 + 0.5)) * depsi1/ewv1_h -
    sigx18_1 * (1 + lambda)/ewu_h) * sigx3_1)/(sigx7 * wzdeno) -
    sigx5 * wzdeno * (2 * (sigx2_1 * sigx16_1) - sigx18_1 *
      sigx3_1)/(sigx7 * wzdeno)^2) * ewz, FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prC * (2 * ((dmusig2 * (ewv2/(2 *
      ewu) - (sigx15_2 * musig2 + 0.5))/ewv2_h - sigx16_2/ewu_h) *
      sigx2_2) - ((((1 + lambda)^2 * ewv2/(2 * ewu) - (epsivu2 *
      sigx17_2 + 0.5)) * depsi2/ewv2_h - sigx18_2 * (1 +
      lambda)/ewu_h) * sigx3_2 + sigx5 * (2 * (sigx2_2 *
      sigx16_2) - sigx18_2 * sigx3_2)/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- (wHvar *
    (-(S * (((((1 + lambda) * ewv1/ewu + S * (epsilon)/ewu_h)/ewv1_h -
      epsivu1/ewu_h) * depsi1 - (sigx19_1 * (1 + lambda) +
      pepsi1)/ewu_h) * sigx3_1 * ewz/wzdeno + ((((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewu_h)/ewv2_h - epsivu2/ewu_h) *
      depsi2 - (sigx19_2 * (1 + lambda) + pepsi2)/ewu_h) *
      prC * sigx3_2 - sigx20 * sigx5/sigx7)/sigx7))) %*%
    Xvar
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((sigx22_2 * sigx3_2 +
      2 * (sigx1_1 * sigx2_1) - (sigx22_1 * sigx3_1 + 2 *
      (sigx1_2 * sigx2_2)))/(sigx7 * wzdeno) - sigx5 *
      wzdeno * sigx21/(sigx7 * wzdeno)^2) * prC * ewz,
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((prC * (2 * (((0.5 * (0.5 * (ewv2_h * musig2/ewu_h) -
      0.5) - 0.5 * (0.5 * sigx9 + 2 * sigx10_2)) * dmusig2 *
      ewv2_h/ewu_h - (sigx11_2 * (0.5 * sigx9 + 2 * sigx10_2) +
      (2 * (lusq * ewu * ewv2/(2 * ewu)^2) - 0.25 * sigx9) *
        pmusig2)) * sigx2_2) - ((0.5 * (0.5 * (epsivu2 *
      (1 + lambda) * ewv2_h/ewu_h) - 0.5) - 0.5 * ((0.5 *
      sigx9 + 2 * sigx12_2) * (1 + lambda))) * depsi2 *
      ewv2_h/ewu_h - (sigx13_2 * (0.5 * sigx9 + 2 * sigx12_2) *
      (1 + lambda) + (2 * (lusq * (1 + lambda) * ewu *
      ewv2/(2 * ewu)^2) - 0.25 * sigx9) * pepsi2)) * (1 +
      lambda) * sigx3_2) + (2 * (((0.5 * (0.5 * (ewv1_h *
      musig1/ewu_h) - 0.5) - 0.5 * (0.5 * sigx9 + 2 * sigx10_1)) *
      dmusig1 * ewv1_h/ewu_h - (sigx11_1 * (0.5 * sigx9 +
      2 * sigx10_1) + (2 * (lusq * ewu * ewv1/(2 * ewu)^2) -
      0.25 * sigx9) * pmusig1)) * sigx2_1) - ((0.5 * (0.5 *
      (epsivu1 * (1 + lambda) * ewv1_h/ewu_h) - 0.5) -
      0.5 * ((0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda))) *
      depsi1 * ewv1_h/ewu_h - (sigx13_1 * (0.5 * sigx9 +
      2 * sigx12_1) * (1 + lambda) + (2 * (lusq * (1 +
      lambda) * ewu * ewv1/(2 * ewu)^2) - 0.25 * sigx9) *
      pepsi1)) * (1 + lambda) * sigx3_1) * ewz/wzdeno -
      (prC * sigx14_2 + sigx14_1 * ewz/wzdeno)^2/sigx7)/sigx7 -
      lambda * (0.5 - lambda * ewu_h/luh) * ewu_h/luh),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((dmusig1 * ewv1_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewv1_h/ewu_h + (0.5 * sigx9 + 2 * sigx10_1) * sigx16_1)) *
      sigx2_1) - (((1 + lambda) * depsi1 * ewv1_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pepsi1/(2 * ewu)^2)) *
      (1 + lambda) * ewv1 - (sigx18_1 * (0.5 * sigx9 +
      2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) - 0.25) *
      depsi1 * ewv1_h/ewu_h)) * (1 + lambda) * sigx3_1)/(sigx7 *
      wzdeno) - (prC * sigx14_2 + sigx14_1 * ewz/wzdeno) *
      wzdeno * (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)/(sigx7 *
      wzdeno)^2) * ewz, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * prC * (2 * (((dmusig2 * ewv2_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) *
      ewv2 - ((0.5 * (sigx15_2 * musig2) - 0.25) * dmusig2 *
      ewv2_h/ewu_h + (0.5 * sigx9 + 2 * sigx10_2) * sigx16_2)) *
      sigx2_2) - ((((1 + lambda) * depsi2 * ewv2_h/(4 *
      (ewu * ewu_h)) - 2 * (ewu * pepsi2/(2 * ewu)^2)) *
      (1 + lambda) * ewv2 - (sigx18_2 * (0.5 * sigx9 +
      2 * sigx12_2) + (0.5 * (epsivu2 * sigx17_2) - 0.25) *
      depsi2 * ewv2_h/ewu_h)) * (1 + lambda) * sigx3_2 +
      (prC * sigx14_2 + sigx14_1 * ewz/wzdeno) * (2 * (sigx2_2 *
        sigx16_2) - sigx18_2 * sigx3_2)/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- (wHvar * (-((((((0.5 * ((1 + lambda) *
    ewv1/ewu + S * (epsilon)/ewu_h) - 0.5 * (epsivu1 * ewv1_h/ewu_h)) *
    (1 + lambda) + 0.5) * depsi1 * ewv1_h/ewu_h - (sigx19_1 *
    (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + ((1 + lambda) *
    ewv1/ewu + 0.5 * sigx9) * pepsi1)) * sigx3_1 * ewz/wzdeno +
    (((0.5 * ((1 + lambda) * ewv2/ewu + S * (epsilon)/ewu_h) -
      0.5 * (epsivu2 * ewv2_h/ewu_h)) * (1 + lambda) +
      0.5) * depsi2 * ewv2_h/ewu_h - (sigx19_2 * (0.5 *
      sigx9 + 2 * sigx12_2) * (1 + lambda) + ((1 + lambda) *
      ewv2/ewu + 0.5 * sigx9) * pepsi2)) * prC * sigx3_2 -
    sigx20 * (prC * sigx14_2 + sigx14_1 * ewz/wzdeno)/sigx7)/sigx7 +
    2 * ((0.5 - lambda * ewu_h/luh) * ewu_h/luh))))) %*%
    uHvar
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx13_2 * sigx3_2 - sigx13_1 * sigx3_1) * (1 + lambda) +
      2 * (sigx11_1 * sigx2_1) - 2 * (sigx11_2 * sigx2_2))/(sigx7 *
      wzdeno) - (prC * sigx14_2 + sigx14_1 * ewz/wzdeno) *
      wzdeno * sigx21/(sigx7 * wzdeno)^2) * prC * ewz,
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_1/2 + (pmusig1 -
      sigx15_1 * dmusig1)/2) * ewv1/ewu - (0.25 * (ewv1_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx15_1^2 * musig1) *
      dmusig1) * sigx2_1) - ((sigx18_1/2 + (pepsi1 - sigx17_1 *
      depsi1)/2) * (1 + lambda)^2 * ewv1/ewu - (0.25 *
      ((1 + lambda) * ewv1_h/ewu_h) + 0.25 * (S * (epsilon)/ewv1_h) -
      epsivu1 * sigx17_1^2) * depsi1) * sigx3_1)/(sigx7 *
      wzdeno) - (2 * (sigx2_1 * sigx16_1) - sigx18_1 *
      sigx3_1)^2 * ewz/(sigx7 * wzdeno)^2) * ewz, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (prC * (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1) *
      (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2) *
      ewz/(sigx7^2 * wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- (wHvar * (-((((sigx19_1 *
    (1 + lambda)/2 + pepsi1) * (1 + lambda) * ewv1/ewu -
    (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewu_h) * sigx17_1 +
      (0.5 - epsivu1 * sigx17_1) * ewv1_h/ewu_h) * depsi1) *
    sigx3_1 - sigx20 * (2 * (sigx2_1 * sigx16_1) - sigx18_1 *
    sigx3_1)/sigx7) * ewz/(sigx7 * wzdeno)))) %*% vHvar
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * prC * (1/(sigx7 * wzdeno) -
      sigx21 * ewz/(sigx7 * wzdeno)^2) * (2 * (sigx2_1 *
      sigx16_1) - sigx18_1 * sigx3_1) * ewz, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prC * (2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu - (0.25 * (ewv2_h/ewu_h) +
      0.25 * (S * (epsilon)/ewv2_h) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - (((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda)^2 * ewv2/ewu - (0.25 *
      ((1 + lambda) * ewv2_h/ewu_h) + 0.25 * (S * (epsilon)/ewv2_h) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2 + prC *
      (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)^2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- (wHvar *
    (-((((sigx19_2 * (1 + lambda)/2 + pepsi2) * (1 + lambda) *
      ewv2/ewu - (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewu_h) *
      sigx17_2 + (0.5 - epsivu2 * sigx17_2) * ewv2_h/ewu_h) *
      depsi2) * sigx3_2 - sigx20 * (2 * (sigx2_2 * sigx16_2) -
      sigx18_2 * sigx3_2)/sigx7) * prC/sigx7))) %*% vHvar
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx18_2 * sigx3_2 - 2 *
      (sigx2_2 * sigx16_2))/(sigx7 * wzdeno) - prC * wzdeno *
      sigx21 * (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)/(sigx7 *
      wzdeno)^2) * prC * ewz, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- sum(wHvar * (4 * (ewu_h^2/luh^2) -
    (((((epsivu1 * ewv1_h - S * (epsilon))/ewu_h - (1 + lambda) *
      ewv1/ewu) * depsi1 * ewv1_h/ewu_h + sigx19_1 * ((1 +
      lambda) * ewv1/ewu + S * (epsilon)/ewu_h) + ewv1 *
      pepsi1/ewu) * sigx3_1 * ewz/wzdeno + (((epsivu2 *
      ewv2_h - S * (epsilon))/ewu_h - (1 + lambda) * ewv2/ewu) *
      depsi2 * ewv2_h/ewu_h + sigx19_2 * ((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewu_h) + ewv2 * pepsi2/ewu) *
      prC * sigx3_2 + sigx20^2/sigx7)/sigx7 + 1/(1 + lambda)^2)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- (wHvar * ((sigx20 * wzdeno * sigx21/(sigx7 * wzdeno)^2 +
    (sigx19_2 * sigx3_2 - sigx19_1 * sigx3_1)/(sigx7 * wzdeno)) *
    prC * ewz)) %*% Zvar
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (prC/(sigx7 * wzdeno) - sigx6_1 *
      ewz/(sigx7 * wzdeno)^2) * prC * sigx21 * ewz, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chesscnsftslnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- (ewz2 * sigx4_2 + ewz1 * sigx4_1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (ewz2 * sigx6_2 + ewz1 * sigx6_1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- (ewz2 * sigx14_2 + ewz1 * sigx14_1)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * ewz1 * sigx3_1 + sigx19_2 * ewz2 *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  sigx24 <- (pi * sigx7 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (ewz2 * (2 * (((musig2/ewvsr2 - 1/ewusr) *
      dmusig2/ewvsr2 - sigx1_2/ewusr) * sigx2_2) - ((epsivu2/ewvsr2 -
      (1 + lambda)/ewusr) * depsi2/ewvsr2 - (1 + lambda) *
      (depsi2/ewvsr2 - (1 + lambda) * pepsi2/ewusr)/ewusr) *
      sigx3_2) + ewz1 * (2 * (((musig1/ewvsr1 - 1/ewusr) *
      dmusig1/ewvsr1 - sigx1_1/ewusr) * sigx2_1) - ((epsivu1/ewvsr1 -
      (1 + lambda)/ewusr) * depsi1/ewvsr1 - (1 + lambda) *
      (depsi1/ewvsr1 - (1 + lambda) * pepsi1/ewusr)/ewusr) *
      sigx3_1) - sigx5^2/sigx7)/sigx7, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (ewz2 * (2 * ((((0.5 +
      0.5 * sigx9 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr +
      (0.5 * (musig2/ewusr) - (0.5 * sigx9 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - ((0.5 * (epsivu2/ewusr) -
      (0.5 * sigx9 + 2 * sigx12_2)/ewvsr2) * depsi2 + (0.5 *
      pepsi2 - sigx13_2 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_2) + ewz1 * (2 * ((((0.5 + 0.5 * sigx9 + 2 *
      sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr + (0.5 *
      (musig1/ewusr) - (0.5 * sigx9 + 2 * sigx10_1)/ewvsr1) *
      dmusig1) * sigx2_1) - ((0.5 * (epsivu1/ewusr) - (0.5 *
      sigx9 + 2 * sigx12_1)/ewvsr1) * depsi1 + (0.5 * pepsi1 -
      sigx13_1 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_1) - sigx20 * sigx5/sigx7)/sigx7, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ewz1 * (2 * ((dmusig1 * (ewv1/(2 * ewu) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr) * sigx2_1) -
    ((((1 + lambda)^2 * ewv1/(2 * ewu) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda)/ewusr) *
      sigx3_1 + sigx5 * sigx21_1/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * ewz2 * (2 * ((dmusig2 * (ewv2/(2 *
      ewu) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr) *
      sigx2_2) - ((((1 + lambda)^2 * ewv2/(2 * ewu) - (epsivu2 *
      sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 * (1 +
      lambda)/ewusr) * sigx3_2 + sigx5 * sigx21_2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S * (((((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr)/ewvsr1 - epsivu1/ewusr) *
      depsi1 - (sigx19_1 * (1 + lambda) + pepsi1)/ewusr) *
      ewz1 * sigx3_1 + ((((1 + lambda) * ewv2/ewu + S *
      (epsilon)/ewusr)/ewvsr2 - epsivu2/ewusr) * depsi2 -
      (sigx19_2 * (1 + lambda) + pepsi2)/ewusr) * ewz2 *
      sigx3_2 - sigx22 * sigx5/sigx7)/sigx7), FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((depsi2/ewvsr2 - (1 +
      lambda) * pepsi2/ewusr) * sigx3_2 + 2 * (sigx1_1 *
      sigx2_1) - ((depsi1/ewvsr1 - (1 + lambda) * pepsi1/ewusr) *
      sigx3_1 + 2 * (sigx1_2 * sigx2_2)))/sigx24 - pi *
      sigx5 * ((Wz)^2 + 1) * sigx23/sigx24^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((ewz2 * (2 * (((0.5 * (0.5 * (ewvsr2 * musig2/ewusr) -
      0.5) - 0.5 * (0.5 * sigx9 + 2 * sigx10_2)) * dmusig2 *
      ewvsr2/ewusr - (sigx11_2 * (0.5 * sigx9 + 2 * sigx10_2) +
      (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 *
        ewu)^2) - 0.25 * sigx9) * pmusig2)) * sigx2_2) -
      ((0.5 * (0.5 * (epsivu2 * (1 + lambda) * ewvsr2/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_2) *
        (1 + lambda))) * depsi2 * ewvsr2/ewusr - (sigx13_2 *
        (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) +
        (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
          ewu * ewv2/(2 * ewu)^2) - 0.25 * sigx9) * pepsi2)) *
        (1 + lambda) * sigx3_2) + ewz1 * (2 * (((0.5 *
      (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 * (0.5 *
      sigx9 + 2 * sigx10_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx11_1 * (0.5 * sigx9 + 2 * sigx10_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * sigx9) * pmusig1)) * sigx2_1) - ((0.5 *
      (0.5 * (epsivu1 * (1 + lambda) * ewvsr1/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_1) *
      (1 + lambda))) * depsi1 * ewvsr1/ewusr - (sigx13_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
        ewv1/(2 * ewu)^2) - 0.25 * sigx9) * pepsi1)) *
      (1 + lambda) * sigx3_1) - sigx20^2/sigx7)/sigx7 -
      lambda * (0.5 - lambda * ewusr/luh) * ewusr/luh),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewz1 * (2 * (((dmusig1 *
      ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pmusig1/(2 *
      ewu)^2)) * ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) *
      dmusig1 * ewvsr1/ewusr + (0.5 * sigx9 + 2 * sigx10_1) *
      sigx16_1)) * sigx2_1) - ((((1 + lambda) * depsi1 *
      ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi1/(2 *
      ewu)^2)) * (1 + lambda) * ewv1 - (sigx18_1 * (0.5 *
      sigx9 + 2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) -
      0.25) * depsi1 * ewvsr1/ewusr)) * (1 + lambda) *
      sigx3_1 + sigx20 * sigx21_1/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ewz2 * (2 * (((dmusig2 *
      ewvsr2/(4 * (ewu * ewusr)) - 2 * (ewu * pmusig2/(2 *
      ewu)^2)) * ewv2 - ((0.5 * (sigx15_2 * musig2) - 0.25) *
      dmusig2 * ewvsr2/ewusr + (0.5 * sigx9 + 2 * sigx10_2) *
      sigx16_2)) * sigx2_2) - ((((1 + lambda) * depsi2 *
      ewvsr2/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi2/(2 *
      ewu)^2)) * (1 + lambda) * ewv2 - (sigx18_2 * (0.5 *
      sigx9 + 2 * sigx12_2) + (0.5 * (epsivu2 * sigx17_2) -
      0.25) * depsi2 * ewvsr2/ewusr)) * (1 + lambda) *
      sigx3_2 + sigx20 * sigx21_2/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((((((0.5 * ((1 + lambda) * ewv1/ewu +
      S * (epsilon)/ewusr) - 0.5 * (epsivu1 * ewvsr1/ewusr)) *
      (1 + lambda) + 0.5) * depsi1 * ewvsr1/ewusr - (sigx19_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + ((1 +
      lambda) * ewv1/ewu + 0.5 * sigx9) * pepsi1)) * ewz1 *
      sigx3_1 + (((0.5 * ((1 + lambda) * ewv2/ewu + S *
      (epsilon)/ewusr) - 0.5 * (epsivu2 * ewvsr2/ewusr)) *
      (1 + lambda) + 0.5) * depsi2 * ewvsr2/ewusr - (sigx19_2 *
      (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) + ((1 +
      lambda) * ewv2/ewu + 0.5 * sigx9) * pepsi2)) * ewz2 *
      sigx3_2 - sigx22 * sigx20/sigx7)/sigx7 + 2 * ((0.5 -
      lambda * ewusr/luh) * ewusr/luh))), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((sigx13_2 * sigx3_2 - sigx13_1 * sigx3_1) * (1 + lambda) +
      2 * (sigx11_1 * sigx2_1) - 2 * (sigx11_2 * sigx2_2))/sigx24 -
      pi * sigx20 * ((Wz)^2 + 1) * sigx23/sigx24^2), FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ewz1 * (2 * (((sigx16_1/2 +
      (pmusig1 - sigx15_1 * dmusig1)/2) * ewv1/ewu - (0.25 *
      (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) -
      sigx15_1^2 * musig1) * dmusig1) * sigx2_1) - (((sigx18_1/2 +
      (pepsi1 - sigx17_1 * depsi1)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * ((1 + lambda) * ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsivu1 * sigx17_1^2) *
      depsi1) * sigx3_1 + ewz1 * sigx21_1^2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (ewz2 * ewz1 * sigx21_1 * sigx21_2/sigx7^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_1 * (1 + lambda)/2 +
      pepsi1) * (1 + lambda) * ewv1/ewu - (((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr) * sigx17_1 + (0.5 -
      epsivu1 * sigx17_1) * ewvsr1/ewusr) * depsi1) * sigx3_1 -
      sigx22 * sigx21_1/sigx7) * ewz1/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1/sigx24 - pi * ((Wz)^2 +
      1) * ewz1 * sigx23/sigx24^2) * sigx21_1, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * ewz2 * (2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu - (0.25 * (ewvsr2/ewusr) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - (((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda)^2 * ewv2/ewu - (0.25 *
      ((1 + lambda) * ewvsr2/ewusr) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2 + ewz2 *
      sigx21_2^2/sigx7))/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_2 * (1 + lambda)/2 +
      pepsi2) * (1 + lambda) * ewv2/ewu - (((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) * sigx17_2 + (0.5 -
      epsivu2 * sigx17_2) * ewvsr2/ewusr) * depsi2) * sigx3_2 -
      sigx22 * sigx21_2/sigx7) * ewz2/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((sigx18_2 * sigx3_2 - 2 *
      (sigx2_2 * sigx16_2))/sigx24 - pi * ((Wz)^2 + 1) *
      ewz2 * sigx23 * sigx21_2/sigx24^2), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- sum(wHvar * (4 * (ewusr^2/luh^2) -
    (((((epsivu1 * ewvsr1 - S * (epsilon))/ewusr - (1 + lambda) *
      ewv1/ewu) * depsi1 * ewvsr1/ewusr + sigx19_1 * ((1 +
      lambda) * ewv1/ewu + S * (epsilon)/ewusr) + ewv1 *
      pepsi1/ewu) * ewz1 * sigx3_1 + (((epsivu2 * ewvsr2 -
      S * (epsilon))/ewusr - (1 + lambda) * ewv2/ewu) *
      depsi2 * ewvsr2/ewusr + sigx19_2 * ((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) + ewv2 * pepsi2/ewu) *
      ewz2 * sigx3_2 + sigx22^2/sigx7)/sigx7 + 1/(1 + lambda)^2)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    ((sigx19_2 * sigx3_2 - sigx19_1 * sigx3_1)/sigx24 + pi *
      sigx22 * ((Wz)^2 + 1) * sigx23/sigx24^2), FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((2 * (sigx2_1 * pmusig1) +
      2 * (pi * Wz * sigx7) + sigx3_2 * pepsi2 - (2 * (sigx2_2 *
      pmusig2) + sigx3_1 * pepsi1)) * sigx23/sigx24^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chesscnsftslnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- ((1 - pwZ) * sigx4_2 + sigx4_1 * pwZ)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - pwZ) * sigx6_2 + sigx6_1 * pwZ)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- ((1 - pwZ) * sigx14_2 + sigx14_1 * pwZ)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * sigx3_1 * pwZ + sigx19_2 * (1 - pwZ) *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((1 - pwZ) * (2 * (((musig2/ewvsr2 -
      1/ewusr) * dmusig2/ewvsr2 - sigx1_2/ewusr) * sigx2_2) -
      ((epsivu2/ewvsr2 - (1 + lambda)/ewusr) * depsi2/ewvsr2 -
        (1 + lambda) * (depsi2/ewvsr2 - (1 + lambda) *
          pepsi2/ewusr)/ewusr) * sigx3_2) + pwZ * (2 *
      (((musig1/ewvsr1 - 1/ewusr) * dmusig1/ewvsr1 - sigx1_1/ewusr) *
        sigx2_1) - ((epsivu1/ewvsr1 - (1 + lambda)/ewusr) *
      depsi1/ewvsr1 - (1 + lambda) * (depsi1/ewvsr1 - (1 +
      lambda) * pepsi1/ewusr)/ewusr) * sigx3_1) - sigx5^2/sigx7)/sigx7,
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 - pwZ) * (2 * ((((0.5 +
      0.5 * sigx9 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr +
      (0.5 * (musig2/ewusr) - (0.5 * sigx9 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - ((0.5 * (epsivu2/ewusr) -
      (0.5 * sigx9 + 2 * sigx12_2)/ewvsr2) * depsi2 + (0.5 *
      pepsi2 - sigx13_2 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_2) + pwZ * (2 * ((((0.5 + 0.5 * sigx9 + 2 *
      sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr + (0.5 *
      (musig1/ewusr) - (0.5 * sigx9 + 2 * sigx10_1)/ewvsr1) *
      dmusig1) * sigx2_1) - ((0.5 * (epsivu1/ewusr) - (0.5 *
      sigx9 + 2 * sigx12_1)/ewvsr1) * depsi1 + (0.5 * pepsi1 -
      sigx13_1 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_1) - sigx20 * sigx5/sigx7)/sigx7, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * pwZ * (2 * ((dmusig1 * (ewv1/(2 * ewu) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr) * sigx2_1) -
    ((((1 + lambda)^2 * ewv1/(2 * ewu) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda)/ewusr) *
      sigx3_1 + sigx5 * sigx21_1/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * (1 - pwZ) * (2 * ((dmusig2 * (ewv2/(2 *
      ewu) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr) *
      sigx2_2) - ((((1 + lambda)^2 * ewv2/(2 * ewu) - (epsivu2 *
      sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 * (1 +
      lambda)/ewusr) * sigx3_2 + sigx5 * sigx21_2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S * (((((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr)/ewvsr1 - epsivu1/ewusr) *
      depsi1 - (sigx19_1 * (1 + lambda) + pepsi1)/ewusr) *
      pwZ * sigx3_1 + ((((1 + lambda) * ewv2/ewu + S *
      (epsilon)/ewusr)/ewvsr2 - epsivu2/ewusr) * depsi2 -
      (sigx19_2 * (1 + lambda) + pepsi2)/ewusr) * (1 -
      pwZ) * sigx3_2 - sigx22 * sigx5/sigx7)/sigx7), FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi2/ewvsr2 - (1 +
      lambda) * pepsi2/ewusr) * sigx3_2 + 2 * (sigx1_1 *
      sigx2_1) - (sigx5 * sigx23/sigx7 + (depsi1/ewvsr1 -
      (1 + lambda) * pepsi1/ewusr) * sigx3_1 + 2 * (sigx1_2 *
      sigx2_2))) * dwZ/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (((1 - pwZ) * (2 * (((0.5 * (0.5 * (ewvsr2 * musig2/ewusr) -
      0.5) - 0.5 * (0.5 * sigx9 + 2 * sigx10_2)) * dmusig2 *
      ewvsr2/ewusr - (sigx11_2 * (0.5 * sigx9 + 2 * sigx10_2) +
      (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 *
        ewu)^2) - 0.25 * sigx9) * pmusig2)) * sigx2_2) -
      ((0.5 * (0.5 * (epsivu2 * (1 + lambda) * ewvsr2/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_2) *
        (1 + lambda))) * depsi2 * ewvsr2/ewusr - (sigx13_2 *
        (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) +
        (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
          ewu * ewv2/(2 * ewu)^2) - 0.25 * sigx9) * pepsi2)) *
        (1 + lambda) * sigx3_2) + pwZ * (2 * (((0.5 *
      (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 * (0.5 *
      sigx9 + 2 * sigx10_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx11_1 * (0.5 * sigx9 + 2 * sigx10_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * sigx9) * pmusig1)) * sigx2_1) - ((0.5 *
      (0.5 * (epsivu1 * (1 + lambda) * ewvsr1/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_1) *
      (1 + lambda))) * depsi1 * ewvsr1/ewusr - (sigx13_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
        ewv1/(2 * ewu)^2) - 0.25 * sigx9) * pepsi1)) *
      (1 + lambda) * sigx3_1) - sigx20^2/sigx7)/sigx7 -
      lambda * (0.5 - lambda * ewusr/luh) * ewusr/luh),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * pwZ * (2 * (((dmusig1 * ewvsr1/(4 *
      (ewu * ewusr)) - 2 * (ewu * pmusig1/(2 * ewu)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr + (0.5 * sigx9 + 2 * sigx10_1) * sigx16_1)) *
      sigx2_1) - ((((1 + lambda) * depsi1 * ewvsr1/(4 *
      (ewu * ewusr)) - 2 * (ewu * pepsi1/(2 * ewu)^2)) *
      (1 + lambda) * ewv1 - (sigx18_1 * (0.5 * sigx9 +
      2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) - 0.25) *
      depsi1 * ewvsr1/ewusr)) * (1 + lambda) * sigx3_1 +
      sigx20 * sigx21_1/sigx7))/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - pwZ) * (2 * (((dmusig2 *
      ewvsr2/(4 * (ewu * ewusr)) - 2 * (ewu * pmusig2/(2 *
      ewu)^2)) * ewv2 - ((0.5 * (sigx15_2 * musig2) - 0.25) *
      dmusig2 * ewvsr2/ewusr + (0.5 * sigx9 + 2 * sigx10_2) *
      sigx16_2)) * sigx2_2) - ((((1 + lambda) * depsi2 *
      ewvsr2/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi2/(2 *
      ewu)^2)) * (1 + lambda) * ewv2 - (sigx18_2 * (0.5 *
      sigx9 + 2 * sigx12_2) + (0.5 * (epsivu2 * sigx17_2) -
      0.25) * depsi2 * ewvsr2/ewusr)) * (1 + lambda) *
      sigx3_2 + sigx20 * sigx21_2/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((((((0.5 * ((1 + lambda) * ewv1/ewu +
      S * (epsilon)/ewusr) - 0.5 * (epsivu1 * ewvsr1/ewusr)) *
      (1 + lambda) + 0.5) * depsi1 * ewvsr1/ewusr - (sigx19_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + ((1 +
      lambda) * ewv1/ewu + 0.5 * sigx9) * pepsi1)) * pwZ *
      sigx3_1 + (((0.5 * ((1 + lambda) * ewv2/ewu + S *
      (epsilon)/ewusr) - 0.5 * (epsivu2 * ewvsr2/ewusr)) *
      (1 + lambda) + 0.5) * depsi2 * ewvsr2/ewusr - (sigx19_2 *
      (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) + ((1 +
      lambda) * ewv2/ewu + 0.5 * sigx9) * pepsi2)) * (1 -
      pwZ) * sigx3_2 - sigx22 * sigx20/sigx7)/sigx7 + 2 *
      ((0.5 - lambda * ewusr/luh) * ewusr/luh))), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((sigx13_2 * sigx3_2 - sigx13_1 * sigx3_1) * (1 + lambda) +
      2 * (sigx11_1 * sigx2_1) - (sigx20 * sigx23/sigx7 +
      2 * (sigx11_2 * sigx2_2))) * dwZ/sigx7, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * pwZ * (2 * (((sigx16_1/2 +
      (pmusig1 - sigx15_1 * dmusig1)/2) * ewv1/ewu - (0.25 *
      (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) -
      sigx15_1^2 * musig1) * dmusig1) * sigx2_1) - (((sigx18_1/2 +
      (pepsi1 - sigx17_1 * depsi1)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * ((1 + lambda) * ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsivu1 * sigx17_1^2) *
      depsi1) * sigx3_1 + pwZ * sigx21_1^2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    ((1 - pwZ) * pwZ * sigx21_1 * sigx21_2/sigx7^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_1 * (1 + lambda)/2 +
      pepsi1) * (1 + lambda) * ewv1/ewu - (((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr) * sigx17_1 + (0.5 -
      epsivu1 * sigx17_1) * ewvsr1/ewusr) * depsi1) * sigx3_1 -
      sigx22 * sigx21_1/sigx7) * pwZ/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx23 * pwZ/sigx7) *
      sigx21_1 * dwZ/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * (1 - pwZ) * (2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu - (0.25 * (ewvsr2/ewusr) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - (((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda)^2 * ewv2/ewu - (0.25 *
      ((1 + lambda) * ewvsr2/ewusr) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2 + (1 -
      pwZ) * sigx21_2^2/sigx7))/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_2 * (1 + lambda)/2 +
      pepsi2) * (1 + lambda) * ewv2/ewu - (((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) * sigx17_2 + (0.5 -
      epsivu2 * sigx17_2) * ewvsr2/ewusr) * depsi2) * sigx3_2 -
      sigx22 * sigx21_2/sigx7) * (1 - pwZ)/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx18_2 * sigx3_2 - ((1 -
      pwZ) * sigx23 * sigx21_2/sigx7 + 2 * (sigx2_2 * sigx16_2))) *
      dwZ/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- sum(wHvar * (4 * (ewusr^2/luh^2) -
    (((((epsivu1 * ewvsr1 - S * (epsilon))/ewusr - (1 + lambda) *
      ewv1/ewu) * depsi1 * ewvsr1/ewusr + sigx19_1 * ((1 +
      lambda) * ewv1/ewu + S * (epsilon)/ewusr) + ewv1 *
      pepsi1/ewu) * pwZ * sigx3_1 + (((epsivu2 * ewvsr2 -
      S * (epsilon))/ewusr - (1 + lambda) * ewv2/ewu) *
      depsi2 * ewvsr2/ewusr + sigx19_2 * ((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) + ewv2 * pepsi2/ewu) *
      (1 - pwZ) * sigx3_2 + sigx22^2/sigx7)/sigx7 + 1/(1 +
      lambda)^2)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (sigx22 * sigx23/sigx7 + sigx19_2 * sigx3_2 - sigx19_1 *
      sigx3_1) * dwZ/sigx7, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23 * dwZ/sigx7 + Wz) *
      sigx23 * dwZ/sigx7), FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chesscnsftslnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda) * ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda) * ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh <- (1 + 2 * (lambda * ewusr))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr)
  sigx2_1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  sigx2_2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  sigx3_1 <- exp(((1 + lambda) * ewv1/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx3_2 <- exp(((1 + lambda) * ewv2/(2 * ewu) + S * (epsilon)/ewusr) *
    (1 + lambda))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda) * pepsi1/ewusr) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda) * pepsi2/ewusr) * sigx3_2)
  sigx5 <- ((1 - prZ) * sigx4_1 + sigx4_2 * prZ)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - prZ) * sigx6_1 + sigx6_2 * prZ)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr)
  sigx9 <- (S * (epsilon)/ewusr)
  sigx10_1 <- (ewu * ewv1/(2 * ewu)^2)
  sigx10_2 <- (ewu * ewv2/(2 * ewu)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda) * ewu * ewv1/(2 * ewu)^2)
  sigx12_2 <- ((1 + lambda) * ewu * ewv2/(2 * ewu)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - (0.5 * sigx9 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - (0.5 * sigx9 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda) *
    sigx3_1)
  sigx14_2 <- (2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda) *
    sigx3_2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda) * ewvsr1/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda) * ewvsr2/ewusr) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda)^2 * ewv1 * pepsi1/(2 * ewu) - sigx17_1 *
    depsi1)
  sigx18_2 <- ((1 + lambda)^2 * ewv2 * pepsi2/(2 * ewu) - sigx17_2 *
    depsi2)
  sigx19_1 <- (((1 + lambda) * ewv1/ewu + S * (epsilon)/ewusr) *
    pepsi1 - depsi1 * ewvsr1/ewusr)
  sigx19_2 <- (((1 + lambda) * ewv2/ewu + S * (epsilon)/ewusr) *
    pepsi2 - depsi2 * ewvsr2/ewusr)
  sigx20 <- ((1 - prZ) * sigx14_1 + sigx14_2 * prZ)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22 <- (sigx19_1 * (1 - prZ) * sigx3_1 + sigx19_2 * prZ *
    sigx3_2)
  sigx23 <- (2 * (sigx2_1 * pmusig1) + sigx3_2 * pepsi2 - (2 *
    (sigx2_2 * pmusig2) + sigx3_1 * pepsi1))
  hessll <- matrix(nrow = nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1, ncol = nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (prZ * (2 * (((musig2/ewvsr2 - 1/ewusr) *
      dmusig2/ewvsr2 - sigx1_2/ewusr) * sigx2_2) - ((epsivu2/ewvsr2 -
      (1 + lambda)/ewusr) * depsi2/ewvsr2 - (1 + lambda) *
      (depsi2/ewvsr2 - (1 + lambda) * pepsi2/ewusr)/ewusr) *
      sigx3_2) + (1 - prZ) * (2 * (((musig1/ewvsr1 - 1/ewusr) *
      dmusig1/ewvsr1 - sigx1_1/ewusr) * sigx2_1) - ((epsivu1/ewvsr1 -
      (1 + lambda)/ewusr) * depsi1/ewvsr1 - (1 + lambda) *
      (depsi1/ewvsr1 - (1 + lambda) * pepsi1/ewusr)/ewusr) *
      sigx3_1) - sigx5^2/sigx7)/sigx7, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (prZ * (2 * ((((0.5 +
      0.5 * sigx9 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr +
      (0.5 * (musig2/ewusr) - (0.5 * sigx9 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - ((0.5 * (epsivu2/ewusr) -
      (0.5 * sigx9 + 2 * sigx12_2)/ewvsr2) * depsi2 + (0.5 *
      pepsi2 - sigx13_2 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_2) + (1 - prZ) * (2 * ((((0.5 + 0.5 * sigx9 +
      2 * sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr +
      (0.5 * (musig1/ewusr) - (0.5 * sigx9 + 2 * sigx10_1)/ewvsr1) *
        dmusig1) * sigx2_1) - ((0.5 * (epsivu1/ewusr) -
      (0.5 * sigx9 + 2 * sigx12_1)/ewvsr1) * depsi1 + (0.5 *
      pepsi1 - sigx13_1 * (1 + lambda))/ewusr) * (1 + lambda) *
      sigx3_1) - sigx20 * sigx5/sigx7)/sigx7, FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * (1 - prZ) * (2 * ((dmusig1 * (ewv1/(2 * ewu) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr) * sigx2_1) -
    ((((1 + lambda)^2 * ewv1/(2 * ewu) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda)/ewusr) *
      sigx3_1 + sigx5 * sigx21_1/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S * prZ * (2 * ((dmusig2 * (ewv2/(2 *
      ewu) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr) *
      sigx2_2) - ((((1 + lambda)^2 * ewv2/(2 * ewu) - (epsivu2 *
      sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 * (1 +
      lambda)/ewusr) * sigx3_2 + sigx5 * sigx21_2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = -wHvar * (S * (((((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr)/ewvsr1 - epsivu1/ewusr) *
      depsi1 - (sigx19_1 * (1 + lambda) + pepsi1)/ewusr) *
      (1 - prZ) * sigx3_1 + ((((1 + lambda) * ewv2/ewu +
      S * (epsilon)/ewusr)/ewvsr2 - epsivu2/ewusr) * depsi2 -
      (sigx19_2 * (1 + lambda) + pepsi2)/ewusr) * prZ *
      sigx3_2 - sigx22 * sigx5/sigx7)/sigx7), FUN = "*"))
  hessll[1:nXvar, (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((depsi2/ewvsr2 - (1 +
      lambda) * pepsi2/ewusr) * sigx3_2 + 2 * (sigx1_1 *
      sigx2_1) - (sigx5 * sigx23/sigx7 + (depsi1/ewvsr1 -
      (1 + lambda) * pepsi1/ewusr) * sigx3_1 + 2 * (sigx1_2 *
      sigx2_2))) * prZ * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((prZ * (2 * (((0.5 * (0.5 * (ewvsr2 * musig2/ewusr) -
      0.5) - 0.5 * (0.5 * sigx9 + 2 * sigx10_2)) * dmusig2 *
      ewvsr2/ewusr - (sigx11_2 * (0.5 * sigx9 + 2 * sigx10_2) +
      (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv2/(2 *
        ewu)^2) - 0.25 * sigx9) * pmusig2)) * sigx2_2) -
      ((0.5 * (0.5 * (epsivu2 * (1 + lambda) * ewvsr2/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_2) *
        (1 + lambda))) * depsi2 * ewvsr2/ewusr - (sigx13_2 *
        (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) +
        (2 * ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) *
          ewu * ewv2/(2 * ewu)^2) - 0.25 * sigx9) * pepsi2)) *
        (1 + lambda) * sigx3_2) + (1 - prZ) * (2 * (((0.5 *
      (0.5 * (ewvsr1 * musig1/ewusr) - 0.5) - 0.5 * (0.5 *
      sigx9 + 2 * sigx10_1)) * dmusig1 * ewvsr1/ewusr -
      (sigx11_1 * (0.5 * sigx9 + 2 * sigx10_1) + (2 * ((1 -
        8 * (ewu^2/(2 * ewu)^2)) * ewu * ewv1/(2 * ewu)^2) -
        0.25 * sigx9) * pmusig1)) * sigx2_1) - ((0.5 *
      (0.5 * (epsivu1 * (1 + lambda) * ewvsr1/ewusr) -
        0.5) - 0.5 * ((0.5 * sigx9 + 2 * sigx12_1) *
      (1 + lambda))) * depsi1 * ewvsr1/ewusr - (sigx13_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + (2 *
      ((1 - 8 * (ewu^2/(2 * ewu)^2)) * (1 + lambda) * ewu *
        ewv1/(2 * ewu)^2) - 0.25 * sigx9) * pepsi1)) *
      (1 + lambda) * sigx3_1) - sigx20^2/sigx7)/sigx7 -
      lambda * (0.5 - lambda * ewusr/luh) * ewusr/luh),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (1 - prZ) * (2 * (((dmusig1 *
      ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pmusig1/(2 *
      ewu)^2)) * ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) *
      dmusig1 * ewvsr1/ewusr + (0.5 * sigx9 + 2 * sigx10_1) *
      sigx16_1)) * sigx2_1) - ((((1 + lambda) * depsi1 *
      ewvsr1/(4 * (ewu * ewusr)) - 2 * (ewu * pepsi1/(2 *
      ewu)^2)) * (1 + lambda) * ewv1 - (sigx18_1 * (0.5 *
      sigx9 + 2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) -
      0.25) * depsi1 * ewvsr1/ewusr)) * (1 + lambda) *
      sigx3_1 + sigx20 * sigx21_1/sigx7))/sigx7, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    nvZVvar + 1):(nXvar + nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * prZ * (2 * (((dmusig2 * ewvsr2/(4 *
      (ewu * ewusr)) - 2 * (ewu * pmusig2/(2 * ewu)^2)) *
      ewv2 - ((0.5 * (sigx15_2 * musig2) - 0.25) * dmusig2 *
      ewvsr2/ewusr + (0.5 * sigx9 + 2 * sigx10_2) * sigx16_2)) *
      sigx2_2) - ((((1 + lambda) * depsi2 * ewvsr2/(4 *
      (ewu * ewusr)) - 2 * (ewu * pepsi2/(2 * ewu)^2)) *
      (1 + lambda) * ewv2 - (sigx18_2 * (0.5 * sigx9 +
      2 * sigx12_2) + (0.5 * (epsivu2 * sigx17_2) - 0.25) *
      depsi2 * ewvsr2/ewusr)) * (1 + lambda) * sigx3_2 +
      sigx20 * sigx21_2/sigx7))/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((((((0.5 * ((1 + lambda) * ewv1/ewu +
      S * (epsilon)/ewusr) - 0.5 * (epsivu1 * ewvsr1/ewusr)) *
      (1 + lambda) + 0.5) * depsi1 * ewvsr1/ewusr - (sigx19_1 *
      (0.5 * sigx9 + 2 * sigx12_1) * (1 + lambda) + ((1 +
      lambda) * ewv1/ewu + 0.5 * sigx9) * pepsi1)) * (1 -
      prZ) * sigx3_1 + (((0.5 * ((1 + lambda) * ewv2/ewu +
      S * (epsilon)/ewusr) - 0.5 * (epsivu2 * ewvsr2/ewusr)) *
      (1 + lambda) + 0.5) * depsi2 * ewvsr2/ewusr - (sigx19_2 *
      (0.5 * sigx9 + 2 * sigx12_2) * (1 + lambda) + ((1 +
      lambda) * ewv2/ewu + 0.5 * sigx9) * pepsi2)) * prZ *
      sigx3_2 - sigx22 * sigx20/sigx7)/sigx7 + 2 * ((0.5 -
      lambda * ewusr/luh) * ewusr/luh))), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((sigx13_2 * sigx3_2 - sigx13_1 * sigx3_1) * (1 + lambda) +
      2 * (sigx11_1 * sigx2_1) - (sigx20 * sigx23/sigx7 +
      2 * (sigx11_2 * sigx2_2))) * prZ * ewz/sigx7, FUN = "*"),
    Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - prZ) * (2 * (((sigx16_1/2 +
      (pmusig1 - sigx15_1 * dmusig1)/2) * ewv1/ewu - (0.25 *
      (ewvsr1/ewusr) + 0.25 * (S * (epsilon)/ewvsr1) -
      sigx15_1^2 * musig1) * dmusig1) * sigx2_1) - (((sigx18_1/2 +
      (pepsi1 - sigx17_1 * depsi1)/2) * (1 + lambda)^2 *
      ewv1/ewu - (0.25 * ((1 + lambda) * ewvsr1/ewusr) +
      0.25 * (S * (epsilon)/ewvsr1) - epsivu1 * sigx17_1^2) *
      depsi1) * sigx3_1 + (1 - prZ) * sigx21_1^2/sigx7))/sigx7,
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar + 2 *
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -wHvar *
    (prZ * (1 - prZ) * sigx21_1 * sigx21_2/sigx7^2), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_1 * (1 + lambda)/2 +
      pepsi1) * (1 + lambda) * ewv1/ewu - (((1 + lambda) *
      ewv1/ewu + S * (epsilon)/ewusr) * sigx17_1 + (0.5 -
      epsivu1 * sigx17_1) * ewvsr1/ewusr) * depsi1) * sigx3_1 -
      sigx22 * sigx21_1/sigx7) * (1 - prZ)/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
      2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - (1 - prZ) * sigx23/sigx7) *
      sigx21_1 * prZ * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(nXvar +
    nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = wHvar * prZ * (2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu - (0.25 * (ewvsr2/ewusr) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - (((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda)^2 * ewv2/ewu - (0.25 *
      ((1 + lambda) * ewvsr2/ewusr) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2 + prZ *
      sigx21_2^2/sigx7))/sigx7, FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((((sigx19_2 * (1 + lambda)/2 +
      pepsi2) * (1 + lambda) * ewv2/ewu - (((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) * sigx17_2 + (0.5 -
      epsivu2 * sigx17_2) * ewvsr2/ewusr) * depsi2) * sigx3_2 -
      sigx22 * sigx21_2/sigx7) * prZ/sigx7), FUN = "*"))
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar), (nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (sigx18_2 * sigx3_2 - (sigx23 *
      sigx21_2 * prZ/sigx7 + 2 * (sigx2_2 * sigx16_2))) *
      prZ * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 1)] <- sum(wHvar * (4 * (ewusr^2/luh^2) -
    (((((epsivu1 * ewvsr1 - S * (epsilon))/ewusr - (1 + lambda) *
      ewv1/ewu) * depsi1 * ewvsr1/ewusr + sigx19_1 * ((1 +
      lambda) * ewv1/ewu + S * (epsilon)/ewusr) + ewv1 *
      pepsi1/ewu) * (1 - prZ) * sigx3_1 + (((epsivu2 *
      ewvsr2 - S * (epsilon))/ewusr - (1 + lambda) * ewv2/ewu) *
      depsi2 * ewvsr2/ewusr + sigx19_2 * ((1 + lambda) *
      ewv2/ewu + S * (epsilon)/ewusr) + ewv2 * pepsi2/ewu) *
      prZ * sigx3_2 + sigx22^2/sigx7)/sigx7 + 1/(1 + lambda)^2)))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 1), (nXvar + nuZUvar +
    2 * nvZVvar + 2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar +
    1)] <- colSums(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (sigx22 * sigx23/sigx7 + sigx19_2 * sigx3_2 - sigx19_1 *
      sigx3_1) * prZ * ewz/sigx7, FUN = "*"))
  hessll[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar + nuZUvar +
    2 * nvZVvar + nZHvar + 1), (nXvar + nuZUvar + 2 * nvZVvar +
    2):(nXvar + nuZUvar + 2 * nvZVvar + nZHvar + 1)] <- crossprod(sweep(Zvar,
    MARGIN = 1, STATS = wHvar * (1 - (sigx23 * prZ/sigx7 +
      1) * ewz) * sigx23 * prZ * ewz/sigx7, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# different sigma_u

## logit specification class membership
chessmcesftslnormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  musig1 <- (ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  musig2 <- (ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  epsi1 <- ((1 + lambda1) * ewv1_h/ewu1_h + S * (epsilon)/ewv1_h)
  epsi2 <- ((1 + lambda2) * ewv2_h/ewu2_h + S * (epsilon)/ewv2_h)
  depsi1 <- dnorm(-epsi1, 0, 1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  pepsi2 <- pnorm(-epsi2)
  sigx1_1 <- (dmusig1/ewv1_h - pmusig1/ewu1_h)
  sigx1_2 <- (dmusig2/ewv2_h - pmusig2/ewu2_h)
  sigx2_1 <- (depsi1/ewv1_h - (1 + lambda1) * pepsi1/ewu1_h)
  sigx2_2 <- (depsi2/ewv2_h - (1 + lambda2) * pepsi2/ewu2_h)
  sigx3_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h)
  sigx3_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h)
  sigx4_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewu1_h) *
    (1 + lambda1))
  sigx4_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewu2_h) *
    (1 + lambda2))
  sigx5_1 <- (2 * (sigx1_1 * sigx3_1) - sigx2_1 * sigx4_1)
  sigx5_2 <- (2 * (sigx1_2 * sigx3_2) - sigx2_2 * sigx4_2)
  lau1 <- (1 + 2 * (lambda1 * ewu1_h))
  lau2 <- (1 + 2 * (lambda2 * ewu2_h))
  sigx6 <- (prC * (1 + lambda2) * sigx5_2/lau2 + (1 + lambda1) *
    sigx5_1 * ewz/(lau1 * wzdeno))
  sigx7_1 <- (2 * (sigx3_1 * pmusig1) - sigx4_1 * pepsi1)
  sigx7_2 <- (2 * (sigx3_2 * pmusig2) - sigx4_2 * pepsi2)
  sigx8 <- (prC * (1 + lambda2) * sigx7_2/lau2 + (1 + lambda1) *
    sigx7_1 * ewz/(lau1 * wzdeno))
  sigx9_1 <- (dmusig1 * ewv1_h/ewu1_h)
  sigx9_2 <- (dmusig2 * ewv2_h/ewu2_h)
  siu1 <- (S * (epsilon)/ewu1_h)
  siu2 <- (S * (epsilon)/ewu2_h)
  uv1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  uv2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx10_1 <- ((0.5 * sigx9_1 - (0.5 * siu1 + 2 * uv1) * pmusig1) *
    sigx3_1)
  sigx10_2 <- ((0.5 * sigx9_2 - (0.5 * siu2 + 2 * uv2) * pmusig2) *
    sigx3_2)
  siv1 <- (depsi1 * ewv1_h/ewu1_h)
  siv2 <- (depsi2 * ewv2_h/ewu2_h)
  sigx11_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx11_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  s12v1 <- (0.5 * siv1 - (0.5 * siu1 + 2 * sigx11_1) * pepsi1) *
    (1 + lambda1)
  s12v2 <- (0.5 * siv2 - (0.5 * siu2 + 2 * sigx11_2) * pepsi2) *
    (1 + lambda2)
  sigx12_1 <- (2 * sigx10_1 - s12v1 * sigx4_1)
  sigx12_2 <- (2 * sigx10_2 - (s12v2 * sigx4_2 + lambda2 *
    sigx7_2 * ewu2_h/lau2))
  sigx13_1 <- (sigx12_1/(lau1 * wzdeno) - lambda1 * wzdeno *
    sigx7_1 * ewu1_h/(lau1 * wzdeno)^2)
  sigx14_1 <- (0.5 * (ewv1_h/ewu1_h) - 0.5 * (S * (epsilon)/ewv1_h))
  sigx14_2 <- (0.5 * (ewv2_h/ewu2_h) - 0.5 * (S * (epsilon)/ewv2_h))
  s15pd1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx14_1 * dmusig1)
  s15pd2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx14_2 * dmusig2)
  sigx15_1 <- (sigx3_1 * s15pd1)
  sigx15_2 <- (sigx3_2 * s15pd2)
  sigx16_1 <- ((1 + lambda1) * ewv1_h/ewu1_h)
  sigx16_2 <- ((1 + lambda2) * ewv2_h/ewu2_h)
  sigx17_1 <- (0.5 * sigx16_1 - 0.5 * (S * (epsilon)/ewv1_h))
  sigx17_2 <- (0.5 * sigx16_2 - 0.5 * (S * (epsilon)/ewv2_h))
  s18pd1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) - sigx17_1 *
    depsi1)
  s18pd2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) - sigx17_2 *
    depsi2)
  sigx18_1 <- (2 * sigx15_1 - s18pd1 * sigx4_1)
  sigx18_2 <- (2 * sigx15_2 - s18pd2 * sigx4_2)
  sigx19_1 <- ((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewu1_h)
  sigx19_2 <- ((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewu2_h)
  sigx20_1 <- (sigx19_1 * pepsi1 - depsi1 * ewv1_h/ewu1_h)
  sigx20_2 <- (sigx19_2 * pepsi2 - depsi2 * ewv2_h/ewu2_h)
  sigx21_1 <- (2 * (sigx3_1 * pmusig1) - (sigx20_1 * (1 + lambda1) +
    pepsi1) * sigx4_1)
  sigx21_2 <- ((sigx20_2 * (1 + lambda2) + pepsi2) * sigx4_2 +
    2 * ((1 + lambda2) * sigx7_2 * ewu2_h/lau2))
  sigx22_1 <- (sigx21_1/(lau1 * wzdeno) - 2 * (wzdeno * (1 +
    lambda1) * sigx7_1 * ewu1_h/(lau1 * wzdeno)^2))
  sigx22_2 <- (2 * (sigx3_2 * pmusig2) - sigx21_2)
  lauz1 <- (1/(lau1 * wzdeno) - lau1 * ewz/(lau1 * wzdeno)^2)
  sigx23 <- ((1 + lambda1) * lauz1 * sigx7_1 - prC * (1 + lambda2) *
    sigx7_2/(lau2 * wzdeno))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * S^2 * (prC * (1 + lambda2) * (2 * (((musig2/ewv2_h -
      1/ewu2_h) * dmusig2/ewv2_h - sigx1_2/ewu2_h) * sigx3_2) -
      ((epsi2/ewv2_h - (1 + lambda2)/ewu2_h) * depsi2/ewv2_h -
        (1 + lambda2) * sigx2_2/ewu2_h) * sigx4_2)/lau2 +
      (1 + lambda1) * (2 * (((musig1/ewv1_h - 1/ewu1_h) *
        dmusig1/ewv1_h - sigx1_1/ewu1_h) * sigx3_1) -
        ((epsi1/ewv1_h - (1 + lambda1)/ewu1_h) * depsi1/ewv1_h -
          (1 + lambda1) * sigx2_1/ewu1_h) * sigx4_1) *
        ewz/(lau1 * wzdeno) - sigx6^2/sigx8)/sigx8, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      siu1 + 2 * uv1) * pmusig1 - 0.5 * sigx9_1)/ewu1_h +
      (0.5 * (musig1/ewu1_h) - (0.5 * siu1 + 2 * uv1)/ewv1_h) *
        dmusig1) * sigx3_1) - ((0.5 * (epsi1/ewu1_h) -
      (0.5 * siu1 + 2 * sigx11_1)/ewv1_h) * depsi1 + (0.5 *
      pepsi1 - s12v1)/ewu1_h) * (1 + lambda1) * sigx4_1)/(lau1 *
      wzdeno) - (sigx6 * sigx13_1/sigx8 + lambda1 * wzdeno *
      sigx5_1 * ewu1_h/(lau1 * wzdeno)^2)) * (1 + lambda1) *
      ewz/sigx8, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      siu2 + 2 * uv2) * pmusig2 - 0.5 * sigx9_2)/ewu2_h +
      (0.5 * (musig2/ewu2_h) - (0.5 * siu2 + 2 * uv2)/ewv2_h) *
        dmusig2) * sigx3_2) - (((0.5 * (epsi2/ewu2_h) -
      (0.5 * siu2 + 2 * sigx11_2)/ewv2_h) * depsi2 + (0.5 *
      pepsi2 - s12v2)/ewu2_h) * (1 + lambda2) * sigx4_2 +
      lambda2 * sigx5_2 * ewu2_h/lau2))/(sigx8 * lau2) -
      sigx6 * lau2 * sigx12_2/(sigx8 * lau2)^2) * prC *
      (1 + lambda2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((2 * ((dmusig1 * (ewv1/(2 * ewu1) - (sigx14_1 *
    musig1 + 0.5))/ewv1_h - s15pd1/ewu1_h) * sigx3_1) - (((1 +
    lambda1)^2 * ewv1/(2 * ewu1) - (epsi1 * sigx17_1 + 0.5)) *
    depsi1/ewv1_h - s18pd1 * (1 + lambda1)/ewu1_h) * sigx4_1)/(sigx8 *
    lau1 * wzdeno) - sigx6 * lau1 * wzdeno * sigx18_1/(sigx8 *
    lau1 * wzdeno)^2) * (1 + lambda1) * ewz, FUN = "*"),
    vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((dmusig2 * (ewv2/(2 *
      ewu2) - (sigx14_2 * musig2 + 0.5))/ewv2_h - s15pd2/ewu2_h) *
      sigx3_2) - (((1 + lambda2)^2 * ewv2/(2 * ewu2) -
      (epsi2 * sigx17_2 + 0.5)) * depsi2/ewv2_h - s18pd2 *
      (1 + lambda2)/ewu2_h) * sigx4_2)/(sigx8 * lau2) -
      sigx6 * lau2 * sigx18_2/(sigx8 * lau2)^2) * prC *
      (1 + lambda2), FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + 2 * nuZUvar + 2 * nvZVvar + 1] <- (wHvar *
    (S * ((2 * (sigx1_1 * sigx3_1) - (((sigx19_1/ewv1_h -
      epsi1/ewu1_h) * depsi1 - (sigx20_1 * (1 + lambda1) +
      2 * pepsi1)/ewu1_h) * (1 + lambda1) + depsi1/ewv1_h) *
      sigx4_1)/(lau1 * wzdeno) - (sigx6 * sigx22_1/sigx8 +
      2 * (wzdeno * (1 + lambda1) * sigx5_1 * ewu1_h/(lau1 *
        wzdeno)^2))) * ewz/sigx8)) %*% Xvar
  hessll[1:nXvar, nXvar + 2 * nuZUvar + 2 * nvZVvar + 2] <- (wHvar *
    (S * ((2 * (sigx1_2 * sigx3_2) - ((((sigx19_2/ewv2_h -
      epsi2/ewu2_h) * depsi2 - (sigx20_2 * (1 + lambda2) +
      2 * pepsi2)/ewu2_h) * (1 + lambda2) + depsi2/ewv2_h) *
      sigx4_2 + 2 * ((1 + lambda2) * sigx5_2 * ewu2_h/lau2)))/(sigx8 *
      lau2) - sigx6 * lau2 * sigx22_2/(sigx8 * lau2)^2) *
      prC)) %*% Xvar
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 + lambda1) * lauz1 *
      sigx5_1 - (sigx6 * sigx23/sigx8 + prC * (1 + lambda2) *
      sigx5_2/(lau2 * wzdeno))) * ewz/sigx8, FUN = "*"),
    Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((0.5 * (0.5 * (ewv1_h * musig1/ewu1_h) - 0.5) -
      0.5 * (0.5 * siu1 + 2 * uv1)) * dmusig1 * ewv1_h/ewu1_h -
      ((0.5 * sigx9_1 - (0.5 * siu1 + 2 * uv1) * pmusig1) *
        (0.5 * siu1 + 2 * uv1) + (2 * ((1 - 8 * (ewu1^2/(2 *
        ewu1)^2)) * ewu1 * ewv1/(2 * ewu1)^2) - 0.25 *
        siu1) * pmusig1)) * sigx3_1) - ((0.5 * (0.5 *
      (epsi1 * (1 + lambda1) * ewv1_h/ewu1_h) - 0.5) -
      0.5 * ((0.5 * siu1 + 2 * sigx11_1) * (1 + lambda1))) *
      depsi1 * ewv1_h/ewu1_h - ((0.5 * siv1 - (0.5 * siu1 +
      2 * sigx11_1) * pepsi1) * (0.5 * siu1 + 2 * sigx11_1) *
      (1 + lambda1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
      (1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2) - 0.25 *
      siu1) * pepsi1)) * (1 + lambda1) * sigx4_1)/(lau1 *
      wzdeno) - (sigx13_1^2 * (1 + lambda1) * ewz/sigx8 +
      lambda1 * ((0.5 - 2 * (lambda1 * lau1 * wzdeno^2 *
        ewu1_h/(lau1 * wzdeno)^2)) * sigx7_1 + 4 * sigx10_1 -
        2 * (s12v1 * sigx4_1)) * wzdeno * ewu1_h/(lau1 *
        wzdeno)^2)) * (1 + lambda1) * ewz/sigx8, FUN = "*"),
    uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (sigx13_1 * prC * lau2 * (1 + lambda1) *
      (1 + lambda2) * sigx12_2 * ewz/(sigx8 * lau2)^2),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((dmusig1 * ewv1_h/(4 *
      (ewu1 * ewu1_h)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (sigx14_1 * musig1) - 0.25) * dmusig1 *
      ewv1_h/ewu1_h + (0.5 * siu1 + 2 * uv1) * s15pd1)) *
      sigx3_1) - (((1 + lambda1) * depsi1 * ewv1_h/(4 *
      (ewu1 * ewu1_h)) - 2 * (ewu1 * pepsi1/(2 * ewu1)^2)) *
      (1 + lambda1) * ewv1 - (s18pd1 * (0.5 * siu1 + 2 *
      sigx11_1) + (0.5 * (epsi1 * sigx17_1) - 0.25) * depsi1 *
      ewv1_h/ewu1_h)) * (1 + lambda1) * sigx4_1)/(sigx8 *
      lau1 * wzdeno) - (sigx13_1 * lau1 * (1 + lambda1) *
      ewz + lambda1 * sigx8 * ewu1_h) * wzdeno * sigx18_1/(sigx8 *
      lau1 * wzdeno)^2) * (1 + lambda1) * ewz, FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (sigx13_1 * prC * lau2 *
      (1 + lambda1) * (1 + lambda2) * sigx18_2 * ewz/(sigx8 *
      lau2)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1] <- (wHvar * (((2 * sigx10_1 - (((0.5 *
    sigx19_1 - 0.5 * (epsi1 * ewv1_h/ewu1_h)) * (1 + lambda1) +
    1) * depsi1 * ewv1_h/ewu1_h - ((sigx20_1 * (1 + lambda1) +
    pepsi1) * (0.5 * siu1 + 2 * sigx11_1) + ((1 + lambda1) *
    ewv1/ewu1 + 0.5 * siu1) * pepsi1)) * (1 + lambda1) *
    sigx4_1)/(lau1 * wzdeno) - (sigx13_1 * sigx22_1 * (1 +
    lambda1) * ewz/sigx8 + wzdeno * (2 * (((0.5 - 2 * (lambda1 *
    lau1 * wzdeno^2 * ewu1_h/(lau1 * wzdeno)^2)) * sigx7_1 +
    2 * sigx10_1 - s12v1 * sigx4_1) * (1 + lambda1)) + lambda1 *
    sigx21_1) * ewu1_h/(lau1 * wzdeno)^2)) * ewz/sigx8)) %*%
    uHvar
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2] <- (wHvar * (-(sigx13_1 * prC * lau2 *
    (1 + lambda1) * sigx22_2 * ewz/(sigx8 * lau2)^2))) %*%
    uHvar
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (lauz1 * sigx12_1 - (sigx23 * sigx13_1 * ewz/sigx8 +
      lambda1 * ((2 - 2 * (lau1^2 * wzdeno^2/(lau1 * wzdeno)^2)) *
        ewz + 1) * sigx7_1 * ewu1_h/(lau1 * wzdeno)^2)) *
    (1 + lambda1) * ewz/sigx8, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewv2_h *
      musig2/ewu2_h) - 0.5) - 0.5 * (0.5 * siu2 + 2 * uv2)) *
      dmusig2 * ewv2_h/ewu2_h - ((0.5 * sigx9_2 - (0.5 *
      siu2 + 2 * uv2) * pmusig2) * (0.5 * siu2 + 2 * uv2) +
      (2 * ((1 - 8 * (ewu2^2/(2 * ewu2)^2)) * ewu2 * ewv2/(2 *
        ewu2)^2) - 0.25 * siu2) * pmusig2)) * sigx3_2) -
      (((0.5 * (0.5 * (epsi2 * (1 + lambda2) * ewv2_h/ewu2_h) -
        0.5) - 0.5 * ((0.5 * siu2 + 2 * sigx11_2) * (1 +
        lambda2))) * depsi2 * ewv2_h/ewu2_h - ((0.5 *
        siv2 - (0.5 * siu2 + 2 * sigx11_2) * pepsi2) *
        (0.5 * siu2 + 2 * sigx11_2) * (1 + lambda2) +
        (2 * ((1 - 8 * (ewu2^2/(2 * ewu2)^2)) * (1 +
          lambda2) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 *
          siu2) * pepsi2)) * (1 + lambda2) * sigx4_2 +
        lambda2 * ((0.5 - lambda2 * ewu2_h/lau2) * sigx7_2 +
          2 * sigx10_2 - s12v2 * sigx4_2) * ewu2_h/lau2))/(sigx8 *
      lau2) - (prC * (1 + lambda2) * sigx12_2 + lambda2 *
      sigx8 * ewu2_h) * sigx12_2/(sigx8 * lau2)^2) * prC *
      (1 + lambda2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prC * lau1 * wzdeno * (1 +
      lambda1) * (1 + lambda2) * sigx12_2 * sigx18_1 *
      ewz/((sigx8 * lau1 * wzdeno)^2 * lau2)), FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((dmusig2 * ewv2_h/(4 * (ewu2 * ewu2_h)) - 2 *
      (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (sigx14_2 *
      musig2) - 0.25) * dmusig2 * ewv2_h/ewu2_h + (0.5 *
      siu2 + 2 * uv2) * s15pd2)) * sigx3_2) - (((1 + lambda2) *
      depsi2 * ewv2_h/(4 * (ewu2 * ewu2_h)) - 2 * (ewu2 *
      pepsi2/(2 * ewu2)^2)) * (1 + lambda2) * ewv2 - (s18pd2 *
      (0.5 * siu2 + 2 * sigx11_2) + (0.5 * (epsi2 * sigx17_2) -
      0.25) * depsi2 * ewv2_h/ewu2_h)) * (1 + lambda2) *
      sigx4_2)/(sigx8 * lau2) - (prC * (1 + lambda2) *
      sigx12_2 + lambda2 * sigx8 * ewu2_h) * sigx18_2/(sigx8 *
      lau2)^2) * prC * (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1] <- (wHvar * (-(sigx22_1 *
    prC * (1 + lambda2) * sigx12_2 * ewz/(sigx8^2 * lau2)))) %*%
    uHvar
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2] <- (wHvar * (((2 * sigx10_2 -
    ((((0.5 * sigx19_2 - 0.5 * (epsi2 * ewv2_h/ewu2_h)) *
      (1 + lambda2) + 1) * depsi2 * ewv2_h/ewu2_h - ((sigx20_2 *
      (1 + lambda2) + pepsi2) * (0.5 * siu2 + 2 * sigx11_2) +
      ((1 + lambda2) * ewv2/ewu2 + 0.5 * siu2) * pepsi2)) *
      sigx4_2 + 2 * (((0.5 - lambda2 * ewu2_h/lau2) * sigx7_2 +
      2 * sigx10_2 - s12v2 * sigx4_2) * ewu2_h/lau2)) *
      (1 + lambda2))/(sigx8 * lau2) - (prC * (1 + lambda2) *
    sigx12_2 + lambda2 * sigx8 * ewu2_h) * sigx22_2/(sigx8 *
    lau2)^2) * prC)) %*% uHvar
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (((sigx23 * sigx12_2/sigx8 +
      (2 * sigx10_2 - s12v2 * sigx4_2)/wzdeno)/lau2 - lambda2 *
      wzdeno * sigx7_2 * ewu2_h/(lau2 * wzdeno)^2) * prC *
      (1 + lambda2) * ewz/sigx8), FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((s15pd1/2 + (pmusig1 -
      sigx14_1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewv1_h/ewu1_h) +
      0.25 * (S * (epsilon)/ewv1_h) - sigx14_1^2 * musig1) *
      dmusig1) * sigx3_1) - ((s18pd1/2 + (pepsi1 - sigx17_1 *
      depsi1)/2) * (1 + lambda1)^2 * ewv1/ewu1 - (0.25 *
      sigx16_1 + 0.25 * (S * (epsilon)/ewv1_h) - epsi1 *
      sigx17_1^2) * depsi1) * sigx4_1)/(sigx8 * lau1 *
      wzdeno) - (1 + lambda1) * sigx18_1^2 * ewz/(sigx8 *
      lau1 * wzdeno)^2) * (1 + lambda1) * ewz, FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (prC * lau2 * (1 + lambda1) * (1 + lambda2) *
      sigx18_1 * sigx18_2 * ewz/((sigx8 * lau2)^2 * lau1 *
      wzdeno)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 1] <- (wHvar * (((2 *
    sigx15_1 - ((((sigx20_1 * (1 + lambda1) + pepsi1)/2 +
    pepsi1) * (1 + lambda1) * ewv1/ewu1 - (sigx19_1 * sigx17_1 +
    (0.5 - epsi1 * sigx17_1) * ewv1_h/ewu1_h) * depsi1) *
    (1 + lambda1) - sigx17_1 * depsi1) * sigx4_1)/(lau1 *
    wzdeno) - (sigx22_1 * ewz/(sigx8 * lau1 * wzdeno) + 2 *
    (wzdeno * ewu1_h/(lau1 * wzdeno)^2)) * (1 + lambda1) *
    sigx18_1) * ewz/sigx8)) %*% vHvar
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    nXvar + 2 * nuZUvar + 2 * nvZVvar + 2] <- (wHvar * (-(prC *
    lau2 * (1 + lambda1) * sigx18_1 * sigx22_2 * ewz/((sigx8 *
    lau2)^2 * lau1 * wzdeno)))) %*% vHvar
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 + lambda1) * (1/(lau1 *
      wzdeno) - (sigx23/(sigx8 * lau1 * wzdeno) + lau1/(lau1 *
      wzdeno)^2) * ewz) * sigx18_1 * ewz/sigx8, FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((s15pd2/2 + (pmusig2 -
      sigx14_2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewv2_h/ewu2_h) +
      0.25 * (S * (epsilon)/ewv2_h) - sigx14_2^2 * musig2) *
      dmusig2) * sigx3_2) - ((s18pd2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda2)^2 * ewv2/ewu2 - (0.25 *
      sigx16_2 + 0.25 * (S * (epsilon)/ewv2_h) - epsi2 *
      sigx17_2^2) * depsi2) * sigx4_2)/(sigx8 * lau2) -
      prC * (1 + lambda2) * sigx18_2^2/(sigx8 * lau2)^2) *
      prC * (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), nXvar + 2 * nuZUvar + 2 * nvZVvar + 1] <- (wHvar *
    (-(sigx22_1 * prC * (1 + lambda2) * sigx18_2 * ewz/(sigx8^2 *
      lau2)))) %*% vHvar
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), nXvar + 2 * nuZUvar + 2 * nvZVvar + 2] <- (wHvar *
    (((2 * sigx15_2 - (((((sigx20_2 * (1 + lambda2) + pepsi2)/2 +
      pepsi2) * (1 + lambda2) * ewv2/ewu2 - (sigx19_2 *
      sigx17_2 + (0.5 - epsi2 * sigx17_2) * ewv2_h/ewu2_h) *
      depsi2) * (1 + lambda2) - sigx17_2 * depsi2) * sigx4_2 +
      2 * ((1 + lambda2) * sigx18_2 * ewu2_h/lau2)))/(sigx8 *
      lau2) - prC * (1 + lambda2) * sigx18_2 * sigx22_2/(sigx8 *
      lau2)^2) * prC)) %*% vHvar
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23/sigx8 + 1/wzdeno) *
      prC * (1 + lambda2) * sigx18_2 * ewz/(sigx8 * lau2)),
    FUN = "*"), Zvar)
  hessll[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1, nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 1] <- sum(wHvar * (-((((((epsi1 *
    ewv1_h - S * (epsilon))/ewu1_h - (1 + lambda1) * ewv1/ewu1) *
    depsi1 * ewv1_h/ewu1_h + ewv1 * pepsi1/ewu1) * (1 + lambda1) +
    (sigx20_1 * (1 + lambda1) + 2 * pepsi1) * sigx19_1 -
    2 * siv1) * sigx4_1/(lau1 * wzdeno) + sigx22_1^2 * ewz/sigx8 +
    wzdeno * (2 * (2 * (sigx3_1 * pmusig1) - ((sigx20_1 *
      (1 + lambda1) + pepsi1) * sigx4_1 + 4 * (lau1 * wzdeno^2 *
      (1 + lambda1) * sigx7_1 * ewu1_h/(lau1 * wzdeno)^2))) +
      2 * sigx21_1) * ewu1_h/(lau1 * wzdeno)^2) * ewz/sigx8)))
  hessll[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1, nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 2] <- sum(wHvar * (-(sigx22_1 *
    prC * lau2 * sigx22_2 * ewz/(sigx8 * lau2)^2)))
  hessll[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1, ((nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar + 2))] <- (wHvar * (((1/(lau1 * wzdeno) -
    (((2 - 4 * (lau1^2 * wzdeno^2/(lau1 * wzdeno)^2)) * ewz +
      2 * wzdeno) * (1 + lambda1) * ewu1_h + lau1 * ewz)/(lau1 *
      wzdeno)^2) * sigx7_1 - (sigx20_1 * (1 + lambda1) *
    lauz1 * sigx4_1 + sigx23 * sigx22_1 * ewz/sigx8)) * ewz/sigx8)) %*%
    Zvar
  hessll[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2, nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 2] <- sum(wHvar * (-(((((((epsi2 *
    ewv2_h - S * (epsilon))/ewu2_h - (1 + lambda2) * ewv2/ewu2) *
    depsi2 * ewv2_h/ewu2_h + ewv2 * pepsi2/ewu2) * (1 + lambda2) +
    (sigx20_2 * (1 + lambda2) + 2 * pepsi2) * sigx19_2 -
    2 * siv2) * sigx4_2 + 2 * (sigx22_2 * ewu2_h/lau2))/(sigx8 *
    lau2) + (prC * sigx22_2 + 2 * (sigx8 * ewu2_h)) * sigx22_2/(sigx8 *
    lau2)^2) * prC)))
  hessll[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2, ((nXvar + 2 *
    nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar + nZHvar + 2))] <- (wHvar * (-(((sigx23 * sigx22_2/sigx8 +
    (2 * (sigx3_2 * pmusig2) - (sigx20_2 * (1 + lambda2) +
      pepsi2) * sigx4_2)/wzdeno)/lau2 - 2 * (wzdeno * (1 +
    lambda2) * sigx7_2 * ewu2_h/(lau2 * wzdeno)^2)) * prC *
    ewz/sigx8))) %*% Zvar
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar + 2), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    (((lau2/(lau2 * wzdeno)^2 + 1/(lau2 * wzdeno^2)) * prC *
      (1 + lambda2) * sigx7_2 - (sigx23^2/sigx8 + lau1 *
      (1 + lambda1) * (2 - 2 * (lau1^2 * wzdeno * ewz/(lau1 *
      wzdeno)^2)) * sigx7_1/(lau1 * wzdeno)^2)) * ewz +
      (1 + lambda1) * lauz1 * sigx7_1 - prC * (1 + lambda2) *
      sigx7_2/(lau2 * wzdeno)) * ewz/sigx8, FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cauchit specification class membership
chessmcesftslnormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- (ewz2 * (1 + lambda2) * sigx4_2/luh2 + ewz1 * (1 +
    lambda1) * sigx4_1/luh1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- (ewz2 * (1 + lambda2) * sigx6_2/luh2 + ewz1 * (1 +
    lambda1) * sigx6_1/luh1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  sigx24 <- (pi * sigx7 * ((Wz)^2 + 1))
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (ewz2 * (1 + lambda2) * (2 * (((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * sigx2_2) -
      ((epsivu2/ewvsr2 - (1 + lambda2)/ewusr2) * depsi2/ewvsr2 -
        (1 + lambda2) * (depsi2/ewvsr2 - (1 + lambda2) *
          pepsi2/ewusr2)/ewusr2) * sigx3_2)/luh2 + ewz1 *
      (1 + lambda1) * (2 * (((musig1/ewvsr1 - 1/ewusr1) *
      dmusig1/ewvsr1 - sigx1_1/ewusr1) * sigx2_1) - ((epsivu1/ewvsr1 -
      (1 + lambda1)/ewusr1) * depsi1/ewvsr1 - (1 + lambda1) *
      (depsi1/ewvsr1 - (1 + lambda1) * pepsi1/ewusr1)/ewusr1) *
      sigx3_1)/luh1 - sigx5^2/sigx7)/sigx7, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_1 + 2 * sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 * sigx9_1 + 2 * sigx10_1)/ewvsr1) *
        dmusig1) * sigx2_1) - (((0.5 * (epsivu1/ewusr1) -
      (0.5 * sigx9_1 + 2 * sigx12_1)/ewvsr1) * depsi1 +
      (0.5 * pepsi1 - sigx13_1 * (1 + lambda1))/ewusr1) *
      (1 + lambda1) * sigx3_1 + lambda1 * sigx4_1 * ewusr1/luh1))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx20_1/(sigx7 * luh1)^2) *
      ewz1 * (1 + lambda1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_2 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 * sigx9_2 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - (((0.5 * (epsivu2/ewusr2) -
      (0.5 * sigx9_2 + 2 * sigx12_2)/ewvsr2) * depsi2 +
      (0.5 * pepsi2 - sigx13_2 * (1 + lambda2))/ewusr2) *
      (1 + lambda2) * sigx3_2 + lambda2 * sigx4_2 * ewusr2/luh2))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx20_2/(sigx7 * luh2)^2) *
      ewz2 * (1 + lambda2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((2 * ((dmusig1 * (ewv1/(2 * ewu1) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr1) * sigx2_1) -
    (((1 + lambda1)^2 * ewv1/(2 * ewu1) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda1)/ewusr1) *
      sigx3_1)/(sigx7 * luh1) - sigx5 * luh1 * sigx21_1/(sigx7 *
    luh1)^2) * ewz1 * (1 + lambda1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((dmusig2 * (ewv2/(2 *
      ewu2) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr2) *
      sigx2_2) - (((1 + lambda2)^2 * ewv2/(2 * ewu2) -
      (epsivu2 * sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 *
      (1 + lambda2)/ewusr2) * sigx3_2)/(sigx7 * luh2) -
      sigx5 * luh2 * sigx21_2/(sigx7 * luh2)^2) * ewz2 *
      (1 + lambda2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_1 * sigx2_1) -
      ((((((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1)/ewvsr1 -
        epsivu1/ewusr1) * depsi1 - (sigx19_1 * (1 + lambda1) +
        2 * pepsi1)/ewusr1) * (1 + lambda1) + depsi1/ewvsr1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx4_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx22_1/(sigx7 * luh1)^2) *
      ewz1, FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_2 * sigx2_2) -
      ((((((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2)/ewvsr2 -
        epsivu2/ewusr2) * depsi2 - (sigx19_2 * (1 + lambda2) +
        2 * pepsi2)/ewusr2) * (1 + lambda2) + depsi2/ewvsr2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx4_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx22_2/(sigx7 * luh2)^2) *
      ewz2, FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * (((1 + lambda1) * sigx4_1/luh1 -
      (1 + lambda2) * sigx4_2/luh2)/sigx24 - pi * sigx5 *
      sigx23 * ((Wz)^2 + 1)/sigx24^2), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) -
      0.5 * (0.5 * sigx9_1 + 2 * sigx10_1)) * dmusig1 *
      ewvsr1/ewusr1 - (sigx11_1 * (0.5 * sigx9_1 + 2 *
      sigx10_1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
      ewu1 * ewv1/(2 * ewu1)^2) - 0.25 * sigx9_1) * pmusig1)) *
      sigx2_1) - (((0.5 * (0.5 * (epsivu1 * (1 + lambda1) *
      ewvsr1/ewusr1) - 0.5) - 0.5 * ((0.5 * sigx9_1 + 2 *
      sigx12_1) * (1 + lambda1))) * depsi1 * ewvsr1/ewusr1 -
      (sigx13_1 * (0.5 * sigx9_1 + 2 * sigx12_1) * (1 +
        lambda1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
        (1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2) - 0.25 *
        sigx9_1) * pepsi1)) * (1 + lambda1) * sigx3_1 +
      lambda1 * ((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
        2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
        sigx3_1) * ewusr1/luh1))/(sigx7 * luh1) - (ewz1 *
      (1 + lambda1) * sigx20_1 + lambda1 * sigx7 * ewusr1) *
      sigx20_1/(sigx7 * luh1)^2) * ewz1 * (1 + lambda1),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * ewz1 * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx20_1 * sigx20_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * sigx9_1 + 2 * sigx10_1) *
      sigx16_1)) * sigx2_1) - (((1 + lambda1) * depsi1 *
      ewvsr1/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pepsi1/(2 *
      ewu1)^2)) * (1 + lambda1) * ewv1 - (sigx18_1 * (0.5 *
      sigx9_1 + 2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) -
      0.25) * depsi1 * ewvsr1/ewusr1)) * (1 + lambda1) *
      sigx3_1)/(sigx7 * luh1) - (ewz1 * (1 + lambda1) *
      sigx20_1 + lambda1 * sigx7 * ewusr1) * sigx21_1/(sigx7 *
      luh1)^2) * ewz1 * (1 + lambda1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * luh2 * (1 +
      lambda1) * (1 + lambda2) * sigx20_1 * sigx21_2/((sigx7 *
      luh2)^2 * luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((2 * (sigx11_1 * sigx2_1) - ((((0.5 *
      ((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) -
      0.5 * (epsivu1 * ewvsr1/ewusr1)) * (1 + lambda1) +
      1) * depsi1 * ewvsr1/ewusr1 - ((sigx19_1 * (1 + lambda1) +
      pepsi1) * (0.5 * sigx9_1 + 2 * sigx12_1) + ((1 +
      lambda1) * ewv1/ewu1 + 0.5 * sigx9_1) * pepsi1)) *
      sigx3_1 + 2 * (((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
      2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
      sigx3_1) * ewusr1/luh1)) * (1 + lambda1))/(sigx7 *
      luh1) - (ewz1 * (1 + lambda1) * sigx20_1 + lambda1 *
      sigx7 * ewusr1) * sigx22_1/(sigx7 * luh1)^2) * ewz1,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * ewz1 * luh2 * (1 + lambda1) *
      sigx20_1 * sigx22_2/((sigx7 * luh2)^2 * luh1)), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 + lambda1) * (1/sigx24 - pi * sigx23 * ((Wz)^2 + 1) *
    ewz1/sigx24^2) * sigx20_1/luh1, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 * sigx9_2 + 2 *
      sigx10_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx11_2 *
      (0.5 * sigx9_2 + 2 * sigx10_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * sigx9_2) *
      pmusig2)) * sigx2_2) - (((0.5 * (0.5 * (epsivu2 *
      (1 + lambda2) * ewvsr2/ewusr2) - 0.5) - 0.5 * ((0.5 *
      sigx9_2 + 2 * sigx12_2) * (1 + lambda2))) * depsi2 *
      ewvsr2/ewusr2 - (sigx13_2 * (0.5 * sigx9_2 + 2 *
      sigx12_2) * (1 + lambda2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * (1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2) -
      0.25 * sigx9_2) * pepsi2)) * (1 + lambda2) * sigx3_2 +
      lambda2 * ((0.5 - lambda2 * ewusr2/luh2) * sigx6_2 +
        2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda2) *
        sigx3_2) * ewusr2/luh2))/(sigx7 * luh2) - (ewz2 *
      (1 + lambda2) * sigx20_2 + lambda2 * sigx7 * ewusr2) *
      sigx20_2/(sigx7 * luh2)^2) * ewz2 * (1 + lambda2),
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * luh1 * (1 +
      lambda1) * (1 + lambda2) * sigx20_2 * sigx21_1/((sigx7 *
      luh1)^2 * luh2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 *
      (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (sigx15_2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      sigx9_2 + 2 * sigx10_2) * sigx16_2)) * sigx2_2) -
      (((1 + lambda2) * depsi2 * ewvsr2/(4 * (ewu2 * ewusr2)) -
        2 * (ewu2 * pepsi2/(2 * ewu2)^2)) * (1 + lambda2) *
        ewv2 - (sigx18_2 * (0.5 * sigx9_2 + 2 * sigx12_2) +
        (0.5 * (epsivu2 * sigx17_2) - 0.25) * depsi2 *
          ewvsr2/ewusr2)) * (1 + lambda2) * sigx3_2)/(sigx7 *
      luh2) - (ewz2 * (1 + lambda2) * sigx20_2 + lambda2 *
      sigx7 * ewusr2) * sigx21_2/(sigx7 * luh2)^2) * ewz2 *
    (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * luh1 * (1 +
      lambda2) * sigx20_2 * sigx22_1/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx11_2 * sigx2_2) -
      ((((0.5 * ((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) -
        0.5 * (epsivu2 * ewvsr2/ewusr2)) * (1 + lambda2) +
        1) * depsi2 * ewvsr2/ewusr2 - ((sigx19_2 * (1 +
        lambda2) + pepsi2) * (0.5 * sigx9_2 + 2 * sigx12_2) +
        ((1 + lambda2) * ewv2/ewu2 + 0.5 * sigx9_2) *
          pepsi2)) * sigx3_2 + 2 * (((0.5 - lambda2 *
        ewusr2/luh2) * sigx6_2 + 2 * (sigx11_2 * sigx2_2) -
        sigx13_2 * (1 + lambda2) * sigx3_2) * ewusr2/luh2)) *
        (1 + lambda2))/(sigx7 * luh2) - (ewz2 * (1 +
      lambda2) * sigx20_2 + lambda2 * sigx7 * ewusr2) *
      sigx22_2/(sigx7 * luh2)^2) * ewz2, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 + lambda2) * (1/sigx24 +
      pi * sigx23 * ((Wz)^2 + 1) * ewz2/sigx24^2) * sigx20_2/luh2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_1/2 + (pmusig1 -
      sigx15_1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - sigx15_1^2 * musig1) *
      dmusig1) * sigx2_1) - ((sigx18_1/2 + (pepsi1 - sigx17_1 *
      depsi1)/2) * (1 + lambda1)^2 * ewv1/ewu1 - (0.25 *
      ((1 + lambda1) * ewvsr1/ewusr1) + 0.25 * (S * (epsilon)/ewvsr1) -
      epsivu1 * sigx17_1^2) * depsi1) * sigx3_1)/(sigx7 *
      luh1) - ewz1 * (1 + lambda1) * sigx21_1^2/(sigx7 *
      luh1)^2) * ewz1 * (1 + lambda1), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (ewz2 * ewz1 * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx21_1 * sigx21_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_1 * sigx16_1) -
      (((((sigx19_1 * (1 + lambda1) + pepsi1)/2 + pepsi1) *
        (1 + lambda1) * ewv1/ewu1 - (((1 + lambda1) *
        ewv1/ewu1 + S * (epsilon)/ewusr1) * sigx17_1 +
        (0.5 - epsivu1 * sigx17_1) * ewvsr1/ewusr1) *
        depsi1) * (1 + lambda1) - sigx17_1 * depsi1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx21_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - ewz1 * (1 + lambda1) * sigx21_1 * sigx22_1/(sigx7 *
      luh1)^2) * ewz1, FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * luh2 * (1 +
      lambda1) * sigx21_1 * sigx22_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 + lambda1) * (1/sigx24 -
      pi * sigx23 * ((Wz)^2 + 1) * ewz1/sigx24^2) * sigx21_1/luh1,
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - ((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda2)^2 * ewv2/ewu2 - (0.25 *
      ((1 + lambda2) * ewvsr2/ewusr2) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2)/(sigx7 *
      luh2) - ewz2 * (1 + lambda2) * sigx21_2^2/(sigx7 *
      luh2)^2) * ewz2 * (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (ewz2 * ewz1 * luh1 * (1 +
      lambda2) * sigx22_1 * sigx21_2/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_2 * sigx16_2) -
      (((((sigx19_2 * (1 + lambda2) + pepsi2)/2 + pepsi2) *
        (1 + lambda2) * ewv2/ewu2 - (((1 + lambda2) *
        ewv2/ewu2 + S * (epsilon)/ewusr2) * sigx17_2 +
        (0.5 - epsivu2 * sigx17_2) * ewvsr2/ewusr2) *
        depsi2) * (1 + lambda2) - sigx17_2 * depsi2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx21_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - ewz2 * (1 + lambda2) * sigx21_2 * sigx22_2/(sigx7 *
      luh2)^2) * ewz2, FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1 + lambda2) * (1/sigx24 +
      pi * sigx23 * ((Wz)^2 + 1) * ewz2/sigx24^2) * sigx21_2/luh2),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- sum(wHvar * (-(((((((epsivu1 *
    ewvsr1 - S * (epsilon))/ewusr1 - (1 + lambda1) * ewv1/ewu1) *
    depsi1 * ewvsr1/ewusr1 + ewv1 * pepsi1/ewu1) * (1 + lambda1) +
    (sigx19_1 * (1 + lambda1) + 2 * pepsi1) * ((1 + lambda1) *
      ewv1/ewu1 + S * (epsilon)/ewusr1) - 2 * (depsi1 *
    ewvsr1/ewusr1)) * sigx3_1 + 2 * (sigx22_1 * ewusr1/luh1))/(sigx7 *
    luh1) + (ewz1 * sigx22_1 + 2 * (sigx7 * ewusr1)) * sigx22_1/(sigx7 *
    luh1)^2) * ewz1)))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-(ewz2 *
    ewz1 * luh2 * sigx22_1 * sigx22_2/((sigx7 * luh2)^2 *
    luh1))))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1/sigx24 - pi * sigx23 * ((Wz)^2 + 1) *
      ewz1/sigx24^2) * sigx22_1/luh1, FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-(((((((epsivu2 *
    ewvsr2 - S * (epsilon))/ewusr2 - (1 + lambda2) * ewv2/ewu2) *
    depsi2 * ewvsr2/ewusr2 + ewv2 * pepsi2/ewu2) * (1 + lambda2) +
    (sigx19_2 * (1 + lambda2) + 2 * pepsi2) * ((1 + lambda2) *
      ewv2/ewu2 + S * (epsilon)/ewusr2) - 2 * (depsi2 *
    ewvsr2/ewusr2)) * sigx3_2 + 2 * (sigx22_2 * ewusr2/luh2))/(sigx7 *
    luh2) + (ewz2 * sigx22_2 + 2 * (sigx7 * ewusr2)) * sigx22_2/(sigx7 *
    luh2)^2) * ewz2)))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = -wHvar * ((1/sigx24 + pi * sigx23 * ((Wz)^2 +
      1) * ewz2/sigx24^2) * sigx22_2/luh2), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar + 2), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    (sigx23 * ((1 + lambda1) * sigx6_1/luh1 + 2 * (pi * Wz *
      sigx7) - (1 + lambda2) * sigx6_2/luh2)/sigx24^2),
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## probit specification class membership
chessmcesftslnormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- ((1 - pwZ) * (1 + lambda2) * sigx4_2/luh2 + (1 +
    lambda1) * sigx4_1 * pwZ/luh1)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - pwZ) * (1 + lambda2) * sigx6_2/luh2 + (1 +
    lambda1) * sigx6_1 * pwZ/luh1)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * ((1 - pwZ) * (1 + lambda2) * (2 * (((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * sigx2_2) -
      ((epsivu2/ewvsr2 - (1 + lambda2)/ewusr2) * depsi2/ewvsr2 -
        (1 + lambda2) * (depsi2/ewvsr2 - (1 + lambda2) *
          pepsi2/ewusr2)/ewusr2) * sigx3_2)/luh2 + pwZ *
      (1 + lambda1) * (2 * (((musig1/ewvsr1 - 1/ewusr1) *
      dmusig1/ewvsr1 - sigx1_1/ewusr1) * sigx2_1) - ((epsivu1/ewvsr1 -
      (1 + lambda1)/ewusr1) * depsi1/ewvsr1 - (1 + lambda1) *
      (depsi1/ewvsr1 - (1 + lambda1) * pepsi1/ewusr1)/ewusr1) *
      sigx3_1)/luh1 - sigx5^2/sigx7)/sigx7, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_1 + 2 * sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 * sigx9_1 + 2 * sigx10_1)/ewvsr1) *
        dmusig1) * sigx2_1) - (((0.5 * (epsivu1/ewusr1) -
      (0.5 * sigx9_1 + 2 * sigx12_1)/ewvsr1) * depsi1 +
      (0.5 * pepsi1 - sigx13_1 * (1 + lambda1))/ewusr1) *
      (1 + lambda1) * sigx3_1 + lambda1 * sigx4_1 * ewusr1/luh1))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx20_1/(sigx7 * luh1)^2) *
      pwZ * (1 + lambda1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_2 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 * sigx9_2 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - (((0.5 * (epsivu2/ewusr2) -
      (0.5 * sigx9_2 + 2 * sigx12_2)/ewvsr2) * depsi2 +
      (0.5 * pepsi2 - sigx13_2 * (1 + lambda2))/ewusr2) *
      (1 + lambda2) * sigx3_2 + lambda2 * sigx4_2 * ewusr2/luh2))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx20_2/(sigx7 * luh2)^2) *
      (1 - pwZ) * (1 + lambda2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((2 * ((dmusig1 * (ewv1/(2 * ewu1) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr1) * sigx2_1) -
    (((1 + lambda1)^2 * ewv1/(2 * ewu1) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda1)/ewusr1) *
      sigx3_1)/(sigx7 * luh1) - sigx5 * luh1 * sigx21_1/(sigx7 *
    luh1)^2) * pwZ * (1 + lambda1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((dmusig2 * (ewv2/(2 *
      ewu2) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr2) *
      sigx2_2) - (((1 + lambda2)^2 * ewv2/(2 * ewu2) -
      (epsivu2 * sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 *
      (1 + lambda2)/ewusr2) * sigx3_2)/(sigx7 * luh2) -
      sigx5 * luh2 * sigx21_2/(sigx7 * luh2)^2) * (1 -
      pwZ) * (1 + lambda2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_1 * sigx2_1) -
      ((((((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1)/ewvsr1 -
        epsivu1/ewusr1) * depsi1 - (sigx19_1 * (1 + lambda1) +
        2 * pepsi1)/ewusr1) * (1 + lambda1) + depsi1/ewvsr1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx4_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx22_1/(sigx7 * luh1)^2) *
      pwZ, FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_2 * sigx2_2) -
      ((((((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2)/ewvsr2 -
        epsivu2/ewusr2) * depsi2 - (sigx19_2 * (1 + lambda2) +
        2 * pepsi2)/ewusr2) * (1 + lambda2) + depsi2/ewvsr2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx4_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx22_2/(sigx7 * luh2)^2) *
      (1 - pwZ), FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 + lambda1) * sigx4_1/luh1 -
      (sigx5 * sigx23/sigx7 + (1 + lambda2) * sigx4_2/luh2)) *
      dwZ/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) -
      0.5 * (0.5 * sigx9_1 + 2 * sigx10_1)) * dmusig1 *
      ewvsr1/ewusr1 - (sigx11_1 * (0.5 * sigx9_1 + 2 *
      sigx10_1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
      ewu1 * ewv1/(2 * ewu1)^2) - 0.25 * sigx9_1) * pmusig1)) *
      sigx2_1) - (((0.5 * (0.5 * (epsivu1 * (1 + lambda1) *
      ewvsr1/ewusr1) - 0.5) - 0.5 * ((0.5 * sigx9_1 + 2 *
      sigx12_1) * (1 + lambda1))) * depsi1 * ewvsr1/ewusr1 -
      (sigx13_1 * (0.5 * sigx9_1 + 2 * sigx12_1) * (1 +
        lambda1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
        (1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2) - 0.25 *
        sigx9_1) * pepsi1)) * (1 + lambda1) * sigx3_1 +
      lambda1 * ((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
        2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
        sigx3_1) * ewusr1/luh1))/(sigx7 * luh1) - (pwZ *
      (1 + lambda1) * sigx20_1 + lambda1 * sigx7 * ewusr1) *
      sigx20_1/(sigx7 * luh1)^2) * pwZ * (1 + lambda1),
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((1 - pwZ) * pwZ * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx20_1 * sigx20_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * sigx9_1 + 2 * sigx10_1) *
      sigx16_1)) * sigx2_1) - (((1 + lambda1) * depsi1 *
      ewvsr1/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pepsi1/(2 *
      ewu1)^2)) * (1 + lambda1) * ewv1 - (sigx18_1 * (0.5 *
      sigx9_1 + 2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) -
      0.25) * depsi1 * ewvsr1/ewusr1)) * (1 + lambda1) *
      sigx3_1)/(sigx7 * luh1) - (pwZ * (1 + lambda1) *
      sigx20_1 + lambda1 * sigx7 * ewusr1) * sigx21_1/(sigx7 *
      luh1)^2) * pwZ * (1 + lambda1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * luh2 *
      (1 + lambda1) * (1 + lambda2) * sigx20_1 * sigx21_2/((sigx7 *
      luh2)^2 * luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((2 * (sigx11_1 * sigx2_1) - ((((0.5 *
      ((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) -
      0.5 * (epsivu1 * ewvsr1/ewusr1)) * (1 + lambda1) +
      1) * depsi1 * ewvsr1/ewusr1 - ((sigx19_1 * (1 + lambda1) +
      pepsi1) * (0.5 * sigx9_1 + 2 * sigx12_1) + ((1 +
      lambda1) * ewv1/ewu1 + 0.5 * sigx9_1) * pepsi1)) *
      sigx3_1 + 2 * (((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
      2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
      sigx3_1) * ewusr1/luh1)) * (1 + lambda1))/(sigx7 *
      luh1) - (pwZ * (1 + lambda1) * sigx20_1 + lambda1 *
      sigx7 * ewusr1) * sigx22_1/(sigx7 * luh1)^2) * pwZ,
    FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * ((1 - pwZ) * pwZ * luh2 * (1 + lambda1) *
      sigx20_1 * sigx22_2/((sigx7 * luh2)^2 * luh1)), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx23 * pwZ/sigx7) * (1 + lambda1) * sigx20_1 *
    dwZ/(sigx7 * luh1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 * sigx9_2 + 2 *
      sigx10_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx11_2 *
      (0.5 * sigx9_2 + 2 * sigx10_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * sigx9_2) *
      pmusig2)) * sigx2_2) - (((0.5 * (0.5 * (epsivu2 *
      (1 + lambda2) * ewvsr2/ewusr2) - 0.5) - 0.5 * ((0.5 *
      sigx9_2 + 2 * sigx12_2) * (1 + lambda2))) * depsi2 *
      ewvsr2/ewusr2 - (sigx13_2 * (0.5 * sigx9_2 + 2 *
      sigx12_2) * (1 + lambda2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * (1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2) -
      0.25 * sigx9_2) * pepsi2)) * (1 + lambda2) * sigx3_2 +
      lambda2 * ((0.5 - lambda2 * ewusr2/luh2) * sigx6_2 +
        2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda2) *
        sigx3_2) * ewusr2/luh2))/(sigx7 * luh2) - ((1 -
      pwZ) * (1 + lambda2) * sigx20_2 + lambda2 * sigx7 *
      ewusr2) * sigx20_2/(sigx7 * luh2)^2) * (1 - pwZ) *
      (1 + lambda2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * luh1 *
      (1 + lambda1) * (1 + lambda2) * sigx20_2 * sigx21_1/((sigx7 *
      luh1)^2 * luh2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 *
      (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (sigx15_2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      sigx9_2 + 2 * sigx10_2) * sigx16_2)) * sigx2_2) -
      (((1 + lambda2) * depsi2 * ewvsr2/(4 * (ewu2 * ewusr2)) -
        2 * (ewu2 * pepsi2/(2 * ewu2)^2)) * (1 + lambda2) *
        ewv2 - (sigx18_2 * (0.5 * sigx9_2 + 2 * sigx12_2) +
        (0.5 * (epsivu2 * sigx17_2) - 0.25) * depsi2 *
          ewvsr2/ewusr2)) * (1 + lambda2) * sigx3_2)/(sigx7 *
      luh2) - ((1 - pwZ) * (1 + lambda2) * sigx20_2 + lambda2 *
      sigx7 * ewusr2) * sigx21_2/(sigx7 * luh2)^2) * (1 -
    pwZ) * (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * luh1 *
      (1 + lambda2) * sigx20_2 * sigx22_1/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx11_2 * sigx2_2) -
      ((((0.5 * ((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) -
        0.5 * (epsivu2 * ewvsr2/ewusr2)) * (1 + lambda2) +
        1) * depsi2 * ewvsr2/ewusr2 - ((sigx19_2 * (1 +
        lambda2) + pepsi2) * (0.5 * sigx9_2 + 2 * sigx12_2) +
        ((1 + lambda2) * ewv2/ewu2 + 0.5 * sigx9_2) *
          pepsi2)) * sigx3_2 + 2 * (((0.5 - lambda2 *
        ewusr2/luh2) * sigx6_2 + 2 * (sigx11_2 * sigx2_2) -
        sigx13_2 * (1 + lambda2) * sigx3_2) * ewusr2/luh2)) *
        (1 + lambda2))/(sigx7 * luh2) - ((1 - pwZ) *
      (1 + lambda2) * sigx20_2 + lambda2 * sigx7 * ewusr2) *
      sigx22_2/(sigx7 * luh2)^2) * (1 - pwZ), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23 * (1 - pwZ)/sigx7 +
      1) * (1 + lambda2) * sigx20_2 * dwZ/(sigx7 * luh2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_1/2 + (pmusig1 -
      sigx15_1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - sigx15_1^2 * musig1) *
      dmusig1) * sigx2_1) - ((sigx18_1/2 + (pepsi1 - sigx17_1 *
      depsi1)/2) * (1 + lambda1)^2 * ewv1/ewu1 - (0.25 *
      ((1 + lambda1) * ewvsr1/ewusr1) + 0.25 * (S * (epsilon)/ewvsr1) -
      epsivu1 * sigx17_1^2) * depsi1) * sigx3_1)/(sigx7 *
      luh1) - pwZ * (1 + lambda1) * sigx21_1^2/(sigx7 *
      luh1)^2) * pwZ * (1 + lambda1), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * ((1 - pwZ) * pwZ * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx21_1 * sigx21_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_1 * sigx16_1) -
      (((((sigx19_1 * (1 + lambda1) + pepsi1)/2 + pepsi1) *
        (1 + lambda1) * ewv1/ewu1 - (((1 + lambda1) *
        ewv1/ewu1 + S * (epsilon)/ewusr1) * sigx17_1 +
        (0.5 - epsivu1 * sigx17_1) * ewvsr1/ewusr1) *
        depsi1) * (1 + lambda1) - sigx17_1 * depsi1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx21_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - pwZ * (1 + lambda1) * sigx21_1 * sigx22_1/(sigx7 *
      luh1)^2) * pwZ, FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * luh2 *
      (1 + lambda1) * sigx21_1 * sigx22_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx23 * pwZ/sigx7) *
      (1 + lambda1) * sigx21_1 * dwZ/(sigx7 * luh1), FUN = "*"),
    Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - ((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda2)^2 * ewv2/ewu2 - (0.25 *
      ((1 + lambda2) * ewvsr2/ewusr2) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2)/(sigx7 *
      luh2) - (1 - pwZ) * (1 + lambda2) * sigx21_2^2/(sigx7 *
      luh2)^2) * (1 - pwZ) * (1 + lambda2), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((1 - pwZ) * pwZ * luh1 *
      (1 + lambda2) * sigx22_1 * sigx21_2/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_2 * sigx16_2) -
      (((((sigx19_2 * (1 + lambda2) + pepsi2)/2 + pepsi2) *
        (1 + lambda2) * ewv2/ewu2 - (((1 + lambda2) *
        ewv2/ewu2 + S * (epsilon)/ewusr2) * sigx17_2 +
        (0.5 - epsivu2 * sigx17_2) * ewvsr2/ewusr2) *
        depsi2) * (1 + lambda2) - sigx17_2 * depsi2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx21_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - (1 - pwZ) * (1 + lambda2) * sigx21_2 * sigx22_2/(sigx7 *
      luh2)^2) * (1 - pwZ), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23 * (1 - pwZ)/sigx7 +
      1) * (1 + lambda2) * sigx21_2 * dwZ/(sigx7 * luh2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- sum(wHvar * (-(((((((epsivu1 *
    ewvsr1 - S * (epsilon))/ewusr1 - (1 + lambda1) * ewv1/ewu1) *
    depsi1 * ewvsr1/ewusr1 + ewv1 * pepsi1/ewu1) * (1 + lambda1) +
    (sigx19_1 * (1 + lambda1) + 2 * pepsi1) * ((1 + lambda1) *
      ewv1/ewu1 + S * (epsilon)/ewusr1) - 2 * (depsi1 *
    ewvsr1/ewusr1)) * sigx3_1 + 2 * (sigx22_1 * ewusr1/luh1))/(sigx7 *
    luh1) + (pwZ * sigx22_1 + 2 * (sigx7 * ewusr1)) * sigx22_1/(sigx7 *
    luh1)^2) * pwZ)))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-((1 -
    pwZ) * pwZ * luh2 * sigx22_1 * sigx22_2/((sigx7 * luh2)^2 *
    luh1))))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx23 * pwZ/sigx7) * sigx22_1 *
      dwZ/(sigx7 * luh1), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-(((((((epsivu2 *
    ewvsr2 - S * (epsilon))/ewusr2 - (1 + lambda2) * ewv2/ewu2) *
    depsi2 * ewvsr2/ewusr2 + ewv2 * pepsi2/ewu2) * (1 + lambda2) +
    (sigx19_2 * (1 + lambda2) + 2 * pepsi2) * ((1 + lambda2) *
      ewv2/ewu2 + S * (epsilon)/ewusr2) - 2 * (depsi2 *
    ewvsr2/ewusr2)) * sigx3_2 + 2 * (sigx22_2 * ewusr2/luh2))/(sigx7 *
    luh2) + ((1 - pwZ) * sigx22_2 + 2 * (sigx7 * ewusr2)) *
    sigx22_2/(sigx7 * luh2)^2) * (1 - pwZ))))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = -wHvar * ((sigx23 * (1 - pwZ)/sigx7 + 1) * sigx22_2 *
      dwZ/(sigx7 * luh2)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar + 2), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = -wHvar *
    ((sigx23 * dwZ/sigx7 + Wz) * sigx23 * dwZ/sigx7), FUN = "*"),
    Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

## cloglog specification class membership
chessmcesftslnormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  lambda1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  lambda2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
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
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  musig1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  musig2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  dmusig1 <- dnorm(-musig1, 0, 1)
  dmusig2 <- dnorm(-musig2, 0, 1)
  pmusig1 <- pnorm(-musig1)
  pmusig2 <- pnorm(-musig2)
  epsivu1 <- ((1 + lambda1) * ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsivu2 <- ((1 + lambda2) * ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsivu1, 0, 1)
  depsi2 <- dnorm(-epsivu2, 0, 1)
  pepsi1 <- pnorm(-epsivu1)
  pepsi2 <- pnorm(-epsivu2)
  luh1 <- (1 + 2 * (lambda1 * ewusr1))
  luh2 <- (1 + 2 * (lambda2 * ewusr2))
  sigx1_1 <- (dmusig1/ewvsr1 - pmusig1/ewusr1)
  sigx1_2 <- (dmusig2/ewvsr2 - pmusig2/ewusr2)
  sigx2_1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  sigx2_2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  sigx3_1 <- exp(((1 + lambda1) * ewv1/(2 * ewu1) + S * (epsilon)/ewusr1) *
    (1 + lambda1))
  sigx3_2 <- exp(((1 + lambda2) * ewv2/(2 * ewu2) + S * (epsilon)/ewusr2) *
    (1 + lambda2))
  sigx4_1 <- (2 * (sigx1_1 * sigx2_1) - (depsi1/ewvsr1 - (1 +
    lambda1) * pepsi1/ewusr1) * sigx3_1)
  sigx4_2 <- (2 * (sigx1_2 * sigx2_2) - (depsi2/ewvsr2 - (1 +
    lambda2) * pepsi2/ewusr2) * sigx3_2)
  sigx5 <- ((1 - prZ) * (1 + lambda1) * sigx4_1/luh1 + (1 +
    lambda2) * sigx4_2 * prZ/luh2)
  sigx6_1 <- (2 * (sigx2_1 * pmusig1) - sigx3_1 * pepsi1)
  sigx6_2 <- (2 * (sigx2_2 * pmusig2) - sigx3_2 * pepsi2)
  sigx7 <- ((1 - prZ) * (1 + lambda1) * sigx6_1/luh1 + (1 +
    lambda2) * sigx6_2 * prZ/luh2)
  sigx8_1 <- (dmusig1 * ewvsr1/ewusr1)
  sigx8_2 <- (dmusig2 * ewvsr2/ewusr2)
  sigx9_1 <- (S * (epsilon)/ewusr1)
  sigx9_2 <- (S * (epsilon)/ewusr2)
  sigx10_1 <- (ewu1 * ewv1/(2 * ewu1)^2)
  sigx10_2 <- (ewu2 * ewv2/(2 * ewu2)^2)
  sigx11_1 <- (0.5 * sigx8_1 - (0.5 * sigx9_1 + 2 * sigx10_1) *
    pmusig1)
  sigx11_2 <- (0.5 * sigx8_2 - (0.5 * sigx9_2 + 2 * sigx10_2) *
    pmusig2)
  sigx12_1 <- ((1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2)
  sigx12_2 <- ((1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2)
  sigx13_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - (0.5 * sigx9_1 +
    2 * sigx12_1) * pepsi1)
  sigx13_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - (0.5 * sigx9_2 +
    2 * sigx12_2) * pepsi2)
  sigx14_1 <- (sigx13_1 * (1 + lambda1) * sigx3_1 + lambda1 *
    sigx6_1 * ewusr1/luh1)
  sigx14_2 <- (sigx13_2 * (1 + lambda2) * sigx3_2 + lambda2 *
    sigx6_2 * ewusr2/luh2)
  sigx15_1 <- (0.5 * (ewvsr1/ewusr1) - 0.5 * (S * (epsilon)/ewvsr1))
  sigx15_2 <- (0.5 * (ewvsr2/ewusr2) - 0.5 * (S * (epsilon)/ewvsr2))
  sigx16_1 <- (ewv1 * pmusig1/(2 * ewu1) - sigx15_1 * dmusig1)
  sigx16_2 <- (ewv2 * pmusig2/(2 * ewu2) - sigx15_2 * dmusig2)
  sigx17_1 <- (0.5 * ((1 + lambda1) * ewvsr1/ewusr1) - 0.5 *
    (S * (epsilon)/ewvsr1))
  sigx17_2 <- (0.5 * ((1 + lambda2) * ewvsr2/ewusr2) - 0.5 *
    (S * (epsilon)/ewvsr2))
  sigx18_1 <- ((1 + lambda1)^2 * ewv1 * pepsi1/(2 * ewu1) -
    sigx17_1 * depsi1)
  sigx18_2 <- ((1 + lambda2)^2 * ewv2 * pepsi2/(2 * ewu2) -
    sigx17_2 * depsi2)
  sigx19_1 <- (((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) *
    pepsi1 - depsi1 * ewvsr1/ewusr1)
  sigx19_2 <- (((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) *
    pepsi2 - depsi2 * ewvsr2/ewusr2)
  sigx20_1 <- (2 * (sigx11_1 * sigx2_1) - sigx14_1)
  sigx20_2 <- (2 * (sigx11_2 * sigx2_2) - sigx14_2)
  sigx21_1 <- (2 * (sigx2_1 * sigx16_1) - sigx18_1 * sigx3_1)
  sigx21_2 <- (2 * (sigx2_2 * sigx16_2) - sigx18_2 * sigx3_2)
  sigx22_1 <- (2 * (sigx2_1 * pmusig1) - ((sigx19_1 * (1 +
    lambda1) + pepsi1) * sigx3_1 + 2 * ((1 + lambda1) * sigx6_1 *
    ewusr1/luh1)))
  sigx22_2 <- (2 * (sigx2_2 * pmusig2) - ((sigx19_2 * (1 +
    lambda2) + pepsi2) * sigx3_2 + 2 * ((1 + lambda2) * sigx6_2 *
    ewusr2/luh2)))
  sigx23 <- ((1 + lambda1) * sigx6_1/luh1 - (1 + lambda2) *
    sigx6_2/luh2)
  hessll <- matrix(nrow = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2, ncol = nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (prZ * (1 + lambda2) * (2 * (((musig2/ewvsr2 -
      1/ewusr2) * dmusig2/ewvsr2 - sigx1_2/ewusr2) * sigx2_2) -
      ((epsivu2/ewvsr2 - (1 + lambda2)/ewusr2) * depsi2/ewvsr2 -
        (1 + lambda2) * (depsi2/ewvsr2 - (1 + lambda2) *
          pepsi2/ewusr2)/ewusr2) * sigx3_2)/luh2 + (1 -
      prZ) * (1 + lambda1) * (2 * (((musig1/ewvsr1 - 1/ewusr1) *
      dmusig1/ewvsr1 - sigx1_1/ewusr1) * sigx2_1) - ((epsivu1/ewvsr1 -
      (1 + lambda1)/ewusr1) * depsi1/ewvsr1 - (1 + lambda1) *
      (depsi1/ewvsr1 - (1 + lambda1) * pepsi1/ewusr1)/ewusr1) *
      sigx3_1)/luh1 - sigx5^2/sigx7)/sigx7, FUN = "*"),
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_1 + 2 * sigx10_1) * pmusig1 - 0.5 * sigx8_1)/ewusr1 +
      (0.5 * (musig1/ewusr1) - (0.5 * sigx9_1 + 2 * sigx10_1)/ewvsr1) *
        dmusig1) * sigx2_1) - (((0.5 * (epsivu1/ewusr1) -
      (0.5 * sigx9_1 + 2 * sigx12_1)/ewvsr1) * depsi1 +
      (0.5 * pepsi1 - sigx13_1 * (1 + lambda1))/ewusr1) *
      (1 + lambda1) * sigx3_1 + lambda1 * sigx4_1 * ewusr1/luh1))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx20_1/(sigx7 * luh1)^2) *
      (1 - prZ) * (1 + lambda1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((((0.5 + 0.5 *
      sigx9_2 + 2 * sigx10_2) * pmusig2 - 0.5 * sigx8_2)/ewusr2 +
      (0.5 * (musig2/ewusr2) - (0.5 * sigx9_2 + 2 * sigx10_2)/ewvsr2) *
        dmusig2) * sigx2_2) - (((0.5 * (epsivu2/ewusr2) -
      (0.5 * sigx9_2 + 2 * sigx12_2)/ewvsr2) * depsi2 +
      (0.5 * pepsi2 - sigx13_2 * (1 + lambda2))/ewusr2) *
      (1 + lambda2) * sigx3_2 + lambda2 * sigx4_2 * ewusr2/luh2))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx20_2/(sigx7 * luh2)^2) *
      prZ * (1 + lambda2), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    S * ((2 * ((dmusig1 * (ewv1/(2 * ewu1) - (sigx15_1 *
    musig1 + 0.5))/ewvsr1 - sigx16_1/ewusr1) * sigx2_1) -
    (((1 + lambda1)^2 * ewv1/(2 * ewu1) - (epsivu1 * sigx17_1 +
      0.5)) * depsi1/ewvsr1 - sigx18_1 * (1 + lambda1)/ewusr1) *
      sigx3_1)/(sigx7 * luh1) - sigx5 * luh1 * sigx21_1/(sigx7 *
    luh1)^2) * (1 - prZ) * (1 + lambda1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * ((dmusig2 * (ewv2/(2 *
      ewu2) - (sigx15_2 * musig2 + 0.5))/ewvsr2 - sigx16_2/ewusr2) *
      sigx2_2) - (((1 + lambda2)^2 * ewv2/(2 * ewu2) -
      (epsivu2 * sigx17_2 + 0.5)) * depsi2/ewvsr2 - sigx18_2 *
      (1 + lambda2)/ewusr2) * sigx3_2)/(sigx7 * luh2) -
      sigx5 * luh2 * sigx21_2/(sigx7 * luh2)^2) * prZ *
      (1 + lambda2), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_1 * sigx2_1) -
      ((((((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1)/ewvsr1 -
        epsivu1/ewusr1) * depsi1 - (sigx19_1 * (1 + lambda1) +
        2 * pepsi1)/ewusr1) * (1 + lambda1) + depsi1/ewvsr1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx4_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - sigx5 * luh1 * sigx22_1/(sigx7 * luh1)^2) *
      (1 - prZ), FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((2 * (sigx1_2 * sigx2_2) -
      ((((((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2)/ewvsr2 -
        epsivu2/ewusr2) * depsi2 - (sigx19_2 * (1 + lambda2) +
        2 * pepsi2)/ewusr2) * (1 + lambda2) + depsi2/ewvsr2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx4_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - sigx5 * luh2 * sigx22_2/(sigx7 * luh2)^2) *
      prZ, FUN = "*"))
  hessll[1:nXvar, (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * S * ((1 + lambda1) * sigx4_1/luh1 -
      (sigx5 * sigx23/sigx7 + (1 + lambda2) * sigx4_2/luh2)) *
      prZ * ewz/sigx7, FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((0.5 * (0.5 * (ewvsr1 * musig1/ewusr1) - 0.5) -
      0.5 * (0.5 * sigx9_1 + 2 * sigx10_1)) * dmusig1 *
      ewvsr1/ewusr1 - (sigx11_1 * (0.5 * sigx9_1 + 2 *
      sigx10_1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
      ewu1 * ewv1/(2 * ewu1)^2) - 0.25 * sigx9_1) * pmusig1)) *
      sigx2_1) - (((0.5 * (0.5 * (epsivu1 * (1 + lambda1) *
      ewvsr1/ewusr1) - 0.5) - 0.5 * ((0.5 * sigx9_1 + 2 *
      sigx12_1) * (1 + lambda1))) * depsi1 * ewvsr1/ewusr1 -
      (sigx13_1 * (0.5 * sigx9_1 + 2 * sigx12_1) * (1 +
        lambda1) + (2 * ((1 - 8 * (ewu1^2/(2 * ewu1)^2)) *
        (1 + lambda1) * ewu1 * ewv1/(2 * ewu1)^2) - 0.25 *
        sigx9_1) * pepsi1)) * (1 + lambda1) * sigx3_1 +
      lambda1 * ((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
        2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
        sigx3_1) * ewusr1/luh1))/(sigx7 * luh1) - ((1 -
      prZ) * (1 + lambda1) * sigx20_1 + lambda1 * sigx7 *
      ewusr1) * sigx20_1/(sigx7 * luh1)^2) * (1 - prZ) *
    (1 + lambda1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prZ * (1 - prZ) * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx20_1 * sigx20_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((dmusig1 * ewvsr1/(4 *
      (ewu1 * ewusr1)) - 2 * (ewu1 * pmusig1/(2 * ewu1)^2)) *
      ewv1 - ((0.5 * (sigx15_1 * musig1) - 0.25) * dmusig1 *
      ewvsr1/ewusr1 + (0.5 * sigx9_1 + 2 * sigx10_1) *
      sigx16_1)) * sigx2_1) - (((1 + lambda1) * depsi1 *
      ewvsr1/(4 * (ewu1 * ewusr1)) - 2 * (ewu1 * pepsi1/(2 *
      ewu1)^2)) * (1 + lambda1) * ewv1 - (sigx18_1 * (0.5 *
      sigx9_1 + 2 * sigx12_1) + (0.5 * (epsivu1 * sigx17_1) -
      0.25) * depsi1 * ewvsr1/ewusr1)) * (1 + lambda1) *
      sigx3_1)/(sigx7 * luh1) - ((1 - prZ) * (1 + lambda1) *
      sigx20_1 + lambda1 * sigx7 * ewusr1) * sigx21_1/(sigx7 *
      luh1)^2) * (1 - prZ) * (1 + lambda1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * luh2 *
      (1 + lambda1) * (1 + lambda2) * sigx20_1 * sigx21_2/((sigx7 *
      luh2)^2 * luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 1)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = wHvar * ((2 * (sigx11_1 * sigx2_1) - ((((0.5 *
      ((1 + lambda1) * ewv1/ewu1 + S * (epsilon)/ewusr1) -
      0.5 * (epsivu1 * ewvsr1/ewusr1)) * (1 + lambda1) +
      1) * depsi1 * ewvsr1/ewusr1 - ((sigx19_1 * (1 + lambda1) +
      pepsi1) * (0.5 * sigx9_1 + 2 * sigx12_1) + ((1 +
      lambda1) * ewv1/ewu1 + 0.5 * sigx9_1) * pepsi1)) *
      sigx3_1 + 2 * (((0.5 - lambda1 * ewusr1/luh1) * sigx6_1 +
      2 * (sigx11_1 * sigx2_1) - sigx13_1 * (1 + lambda1) *
      sigx3_1) * ewusr1/luh1)) * (1 + lambda1))/(sigx7 *
      luh1) - ((1 - prZ) * (1 + lambda1) * sigx20_1 + lambda1 *
      sigx7 * ewusr1) * sigx22_1/(sigx7 * luh1)^2) * (1 -
      prZ), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 2)] <- colSums(sweep(uHvar, MARGIN = 1,
    STATS = -wHvar * (prZ * (1 - prZ) * luh2 * (1 + lambda1) *
      sigx20_1 * sigx22_2/((sigx7 * luh2)^2 * luh1)), FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (1 - sigx23 * (1 - prZ)/sigx7) * (1 + lambda1) * sigx20_1 *
    prZ * ewz/(sigx7 * luh1), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    nuZUvar + 1):(nXvar + 2 * nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((0.5 * (0.5 * (ewvsr2 *
      musig2/ewusr2) - 0.5) - 0.5 * (0.5 * sigx9_2 + 2 *
      sigx10_2)) * dmusig2 * ewvsr2/ewusr2 - (sigx11_2 *
      (0.5 * sigx9_2 + 2 * sigx10_2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * ewu2 * ewv2/(2 * ewu2)^2) - 0.25 * sigx9_2) *
      pmusig2)) * sigx2_2) - (((0.5 * (0.5 * (epsivu2 *
      (1 + lambda2) * ewvsr2/ewusr2) - 0.5) - 0.5 * ((0.5 *
      sigx9_2 + 2 * sigx12_2) * (1 + lambda2))) * depsi2 *
      ewvsr2/ewusr2 - (sigx13_2 * (0.5 * sigx9_2 + 2 *
      sigx12_2) * (1 + lambda2) + (2 * ((1 - 8 * (ewu2^2/(2 *
      ewu2)^2)) * (1 + lambda2) * ewu2 * ewv2/(2 * ewu2)^2) -
      0.25 * sigx9_2) * pepsi2)) * (1 + lambda2) * sigx3_2 +
      lambda2 * ((0.5 - lambda2 * ewusr2/luh2) * sigx6_2 +
        2 * (sigx11_2 * sigx2_2) - sigx13_2 * (1 + lambda2) *
        sigx3_2) * ewusr2/luh2))/(sigx7 * luh2) - (prZ *
      (1 + lambda2) * sigx20_2 + lambda2 * sigx7 * ewusr2) *
      sigx20_2/(sigx7 * luh2)^2) * prZ * (1 + lambda2),
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * luh1 *
      (1 + lambda1) * (1 + lambda2) * sigx20_2 * sigx21_1/((sigx7 *
      luh1)^2 * luh2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar + 2 *
    nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    ((2 * (((dmusig2 * ewvsr2/(4 * (ewu2 * ewusr2)) - 2 *
      (ewu2 * pmusig2/(2 * ewu2)^2)) * ewv2 - ((0.5 * (sigx15_2 *
      musig2) - 0.25) * dmusig2 * ewvsr2/ewusr2 + (0.5 *
      sigx9_2 + 2 * sigx10_2) * sigx16_2)) * sigx2_2) -
      (((1 + lambda2) * depsi2 * ewvsr2/(4 * (ewu2 * ewusr2)) -
        2 * (ewu2 * pepsi2/(2 * ewu2)^2)) * (1 + lambda2) *
        ewv2 - (sigx18_2 * (0.5 * sigx9_2 + 2 * sigx12_2) +
        (0.5 * (epsivu2 * sigx17_2) - 0.25) * depsi2 *
          ewvsr2/ewusr2)) * (1 + lambda2) * sigx3_2)/(sigx7 *
      luh2) - (prZ * (1 + lambda2) * sigx20_2 + lambda2 *
      sigx7 * ewusr2) * sigx21_2/(sigx7 * luh2)^2) * prZ *
    (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * luh1 *
      (1 + lambda2) * sigx20_2 * sigx22_1/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx11_2 * sigx2_2) -
      ((((0.5 * ((1 + lambda2) * ewv2/ewu2 + S * (epsilon)/ewusr2) -
        0.5 * (epsivu2 * ewvsr2/ewusr2)) * (1 + lambda2) +
        1) * depsi2 * ewvsr2/ewusr2 - ((sigx19_2 * (1 +
        lambda2) + pepsi2) * (0.5 * sigx9_2 + 2 * sigx12_2) +
        ((1 + lambda2) * ewv2/ewu2 + 0.5 * sigx9_2) *
          pepsi2)) * sigx3_2 + 2 * (((0.5 - lambda2 *
        ewusr2/luh2) * sigx6_2 + 2 * (sigx11_2 * sigx2_2) -
        sigx13_2 * (1 + lambda2) * sigx3_2) * ewusr2/luh2)) *
        (1 + lambda2))/(sigx7 * luh2) - (prZ * (1 + lambda2) *
      sigx20_2 + lambda2 * sigx7 * ewusr2) * sigx22_2/(sigx7 *
      luh2)^2) * prZ, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23 * prZ/sigx7 + 1) *
      (1 + lambda2) * sigx20_2 * prZ * ewz/(sigx7 * luh2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_1/2 + (pmusig1 -
      sigx15_1 * dmusig1)/2) * ewv1/ewu1 - (0.25 * (ewvsr1/ewusr1) +
      0.25 * (S * (epsilon)/ewvsr1) - sigx15_1^2 * musig1) *
      dmusig1) * sigx2_1) - ((sigx18_1/2 + (pepsi1 - sigx17_1 *
      depsi1)/2) * (1 + lambda1)^2 * ewv1/ewu1 - (0.25 *
      ((1 + lambda1) * ewvsr1/ewusr1) + 0.25 * (S * (epsilon)/ewvsr1) -
      epsivu1 * sigx17_1^2) * depsi1) * sigx3_1)/(sigx7 *
      luh1) - (1 - prZ) * (1 + lambda1) * sigx21_1^2/(sigx7 *
      luh1)^2) * (1 - prZ) * (1 + lambda1), FUN = "*"),
    vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
      2 * nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1,
    STATS = -wHvar * (prZ * (1 - prZ) * luh2 * (1 + lambda1) *
      (1 + lambda2) * sigx21_1 * sigx21_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_1 * sigx16_1) -
      (((((sigx19_1 * (1 + lambda1) + pepsi1)/2 + pepsi1) *
        (1 + lambda1) * ewv1/ewu1 - (((1 + lambda1) *
        ewv1/ewu1 + S * (epsilon)/ewusr1) * sigx17_1 +
        (0.5 - epsivu1 * sigx17_1) * ewvsr1/ewusr1) *
        depsi1) * (1 + lambda1) - sigx17_1 * depsi1) *
        sigx3_1 + 2 * ((1 + lambda1) * sigx21_1 * ewusr1/luh1)))/(sigx7 *
      luh1) - (1 - prZ) * (1 + lambda1) * sigx21_1 * sigx22_1/(sigx7 *
      luh1)^2) * (1 - prZ), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * luh2 *
      (1 + lambda1) * sigx21_1 * sigx22_2/((sigx7 * luh2)^2 *
      luh1)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + nvZVvar),
    (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
      nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (1 - sigx23 * (1 - prZ)/sigx7) *
      (1 + lambda1) * sigx21_1 * prZ * ewz/(sigx7 * luh1),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (((sigx16_2/2 + (pmusig2 -
      sigx15_2 * dmusig2)/2) * ewv2/ewu2 - (0.25 * (ewvsr2/ewusr2) +
      0.25 * (S * (epsilon)/ewvsr2) - sigx15_2^2 * musig2) *
      dmusig2) * sigx2_2) - ((sigx18_2/2 + (pepsi2 - sigx17_2 *
      depsi2)/2) * (1 + lambda2)^2 * ewv2/ewu2 - (0.25 *
      ((1 + lambda2) * ewvsr2/ewusr2) + 0.25 * (S * (epsilon)/ewvsr2) -
      epsivu2 * sigx17_2^2) * depsi2) * sigx3_2)/(sigx7 *
      luh2) - prZ * (1 + lambda2) * sigx21_2^2/(sigx7 *
      luh2)^2) * prZ * (1 + lambda2), FUN = "*"), vHvar)
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 1)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * (prZ * (1 - prZ) * luh1 *
      (1 + lambda2) * sigx22_1 * sigx21_2/((sigx7 * luh1)^2 *
      luh2)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 2)] <- colSums(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * ((2 * (sigx2_2 * sigx16_2) -
      (((((sigx19_2 * (1 + lambda2) + pepsi2)/2 + pepsi2) *
        (1 + lambda2) * ewv2/ewu2 - (((1 + lambda2) *
        ewv2/ewu2 + S * (epsilon)/ewusr2) * sigx17_2 +
        (0.5 - epsivu2 * sigx17_2) * ewvsr2/ewusr2) *
        depsi2) * (1 + lambda2) - sigx17_2 * depsi2) *
        sigx3_2 + 2 * ((1 + lambda2) * sigx21_2 * ewusr2/luh2)))/(sigx7 *
      luh2) - prZ * (1 + lambda2) * sigx21_2 * sigx22_2/(sigx7 *
      luh2)^2) * prZ, FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * nuZUvar +
    2 * nvZVvar), (nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = -wHvar * ((sigx23 * prZ/sigx7 + 1) *
      (1 + lambda2) * sigx21_2 * prZ * ewz/(sigx7 * luh2)),
    FUN = "*"), Zvar)
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 1)] <- sum(wHvar * (-(((((((epsivu1 *
    ewvsr1 - S * (epsilon))/ewusr1 - (1 + lambda1) * ewv1/ewu1) *
    depsi1 * ewvsr1/ewusr1 + ewv1 * pepsi1/ewu1) * (1 + lambda1) +
    (sigx19_1 * (1 + lambda1) + 2 * pepsi1) * ((1 + lambda1) *
      ewv1/ewu1 + S * (epsilon)/ewusr1) - 2 * (depsi1 *
    ewvsr1/ewusr1)) * sigx3_1 + 2 * (sigx22_1 * ewusr1/luh1))/(sigx7 *
    luh1) + ((1 - prZ) * sigx22_1 + 2 * (sigx7 * ewusr1)) *
    sigx22_1/(sigx7 * luh1)^2) * (1 - prZ))))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-(prZ *
    (1 - prZ) * luh2 * sigx22_1 * sigx22_2/((sigx7 * luh2)^2 *
    luh1))))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 1), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = wHvar * (1 - sigx23 * (1 - prZ)/sigx7) * sigx22_1 *
      prZ * ewz/(sigx7 * luh1), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 2)] <- sum(wHvar * (-(((((((epsivu2 *
    ewvsr2 - S * (epsilon))/ewusr2 - (1 + lambda2) * ewv2/ewu2) *
    depsi2 * ewvsr2/ewusr2 + ewv2 * pepsi2/ewu2) * (1 + lambda2) +
    (sigx19_2 * (1 + lambda2) + 2 * pepsi2) * ((1 + lambda2) *
      ewv2/ewu2 + S * (epsilon)/ewusr2) - 2 * (depsi2 *
    ewvsr2/ewusr2)) * sigx3_2 + 2 * (sigx22_2 * ewusr2/luh2))/(sigx7 *
    luh2) + (prZ * sigx22_2 + 2 * (sigx7 * ewusr2)) * sigx22_2/(sigx7 *
    luh2)^2) * prZ)))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 2), (nXvar +
    2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 * nuZUvar +
    2 * nvZVvar + nZHvar + 2)] <- colSums(sweep(Zvar, MARGIN = 1,
    STATS = -wHvar * ((sigx23 * prZ/sigx7 + 1) * sigx22_2 *
      prZ * ewz/(sigx7 * luh2)), FUN = "*"))
  hessll[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar + 2 *
    nuZUvar + 2 * nvZVvar + nZHvar + 2), (nXvar + 2 * nuZUvar +
    2 * nvZVvar + 3):(nXvar + 2 * nuZUvar + 2 * nvZVvar +
    nZHvar + 2)] <- crossprod(sweep(Zvar, MARGIN = 1, STATS = wHvar *
    sigx23 * (1 - (sigx23 * prZ/sigx7 + 1) * ewz) * prZ *
    ewz/sigx7, FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf tsl-normal distribution
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
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @noRd
# Same sigma_u

## logit specification class membership
cnsftslnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsftslnormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsftslnormlike_logit,
      grad = cgradcnsftslnormlike_logit, hess = chesscnsftslnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftslnormlike_logit(mleObj$par,
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
      mleObj$hessian <- chesscnsftslnormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftslnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftslnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftslnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## cauchit specification class membership
cnsftslnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsftslnormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsftslnormlike_cauchit,
      grad = cgradcnsftslnormlike_cauchit, hess = chesscnsftslnormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftslnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chesscnsftslnormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftslnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftslnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftslnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## probit specification class membership
cnsftslnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsftslnormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsftslnormlike_probit,
      grad = cgradcnsftslnormlike_probit, hess = chesscnsftslnormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftslnormlike_probit(mleObj$par,
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
      mleObj$hessian <- chesscnsftslnormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftslnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftslnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftslnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## cloglog specification class membership
cnsftslnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstcnsftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsftslnormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradcnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsftslnormlike_cloglog,
      grad = cgradcnsftslnormlike_cloglog, hess = chesscnsftslnormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chesscnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(ccnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradcnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chesscnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradcnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chesscnsftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsftslnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chesscnsftslnormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chesscnsftslnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- ccnsftslnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsftslnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

# different sigma_u

## logit specification class membership
mcesftslnormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesftslnormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesftslnormlike_logit,
      grad = cgradmcesftslnormlike_logit, hess = chessmcesftslnormlike_logit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesftslnormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesftslnormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftslnormlike_logit(mleObj$par,
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
      mleObj$hessian <- chessmcesftslnormlike_logit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftslnormlike_logit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftslnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftslnormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## cauchit specification class membership
mcesftslnormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesftslnormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesftslnormlike_cauchit,
      grad = cgradmcesftslnormlike_cauchit, hess = chessmcesftslnormlike_cauchit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesftslnormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesftslnormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftslnormlike_cauchit(mleObj$par,
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
      mleObj$hessian <- chessmcesftslnormlike_cauchit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftslnormlike_cauchit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftslnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftslnormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## probit specification class membership
mcesftslnormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesftslnormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesftslnormlike_probit,
      grad = cgradmcesftslnormlike_probit, hess = chessmcesftslnormlike_probit,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesftslnormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesftslnormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftslnormlike_probit(mleObj$par,
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
      mleObj$hessian <- chessmcesftslnormlike_probit(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftslnormlike_probit(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftslnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftslnormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

## cloglog specification class membership
mcesftslnormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  whichStart, initIter, initAlg, hessianType, qac) {
  start_st <- if (!is.null(start))
    start else cstmcesftslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar, uHvar = uHvar,
    nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar,
    nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, whichStart = whichStart,
    initIter = initIter, initAlg = initAlg, printInfo = printInfo,
    tol = tol)
  initTSL <- start_st$initTSL
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesftslnormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    gr = function(parm) -colSums(cgradmcesftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
    hessian = 0, control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesftslnormlike_cloglog,
      grad = cgradmcesftslnormlike_cloglog, hess = chessmcesftslnormlike_cloglog,
      start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
    sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesftslnormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hs = function(parm) as(-chessmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L,
        preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
      fn = function(parm) -sum(cmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gr = function(parm) -colSums(cgradmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hess = function(parm) -chessmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      print.info = printInfo, maxiter = itermax, epsa = gradtol,
      epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      gradient = function(parm) -colSums(cgradmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)),
      hessian = function(parm) -chessmcesftslnormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar),
      control = list(iter.max = itermax, trace = printInfo,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesftslnormlike_cloglog(mleObj$par,
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
      mleObj$hessian <- chessmcesftslnormlike_cloglog(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1")
      mleObj$hessian <- chessmcesftslnormlike_cloglog(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cmcesftslnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesftslnormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initTSL = initTSL))
}

# Conditional efficiencies estimation ----------
#' efficiencies for tsl-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# same sigma_u

## logit specification class membership
ccnsftslnormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
ccnsftslnormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
ccnsftslnormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
ccnsftslnormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# different sigma_u

## logit specification class membership
cmcesftslnormeff_logit <- function(object, level) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
cmcesftslnormeff_cauchit <- function(object, level) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
cmcesftslnormeff_probit <- function(object, level) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
cmcesftslnormeff_cloglog <- function(object, level) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- exp(Wv1/2) * (2 * exp(A1) * (dnorm(a1) + a1 * pnorm(a1)) -
    exp(B1) * (dnorm(b1) + b1 * pnorm(b1)))/(2 * exp(A1) *
    pnorm(a1) - exp(B1) * pnorm(b1))
  u_c2 <- exp(Wv2/2) * (2 * exp(A2) * (dnorm(a2) + a2 * pnorm(a2)) -
    exp(B2) * (dnorm(b2) + b2 * pnorm(b2)))/(2 * exp(A2) *
    pnorm(a2) - exp(B2) * pnorm(b2))
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    teBC_c1 <- (2 * exp(A1) * exp(-a1 * exp(Wv1/2) + exp(Wv1)/2) *
      pnorm(a1 - exp(Wv1/2)) - exp(B1) * exp(-b1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(b1 - exp(Wv1/2)))/(2 * exp(A1) *
      pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_c2 <- (2 * exp(A2) * exp(-a2 * exp(Wv2/2) + exp(Wv2)/2) *
      pnorm(a2 - exp(Wv2/2)) - exp(B2) * exp(-b2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(b2 - exp(Wv2/2)))/(2 * exp(A2) *
      pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- (2 * exp(A1) * exp(a1 * exp(Wv1/2) +
      exp(Wv1)/2) * pnorm(a1 + exp(Wv1/2)) - exp(B1) *
      exp(b1 * exp(Wv1/2) + exp(Wv1)/2) * pnorm(b1 + exp(Wv1/2)))/(2 *
      exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
    teBC_reciprocal_c2 <- (2 * exp(A2) * exp(a2 * exp(Wv2/2) +
      exp(Wv2)/2) * pnorm(a2 + exp(Wv2/2)) - exp(B2) *
      exp(b2 * exp(Wv2/2) + exp(Wv2)/2) * pnorm(b2 + exp(Wv2/2)))/(2 *
      exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, teJLMS_c = teJLMS_c, teBC_c = teBC_c,
      teBC_reciprocal_c = teBC_reciprocal_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, teBC_c1 = teBC_c1,
      teBC_reciprocal_c1 = teBC_reciprocal_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, teBC_c2 = teBC_c2,
      teBC_reciprocal_c2 = teBC_reciprocal_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2, effBC_c1 = effBC_c1, effBC_c2 = effBC_c2,
      ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      u_c = u_c, PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for cnsf tsl-normal distribution
#' @param object object of class sfacross
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfmargtslnorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtslnorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
ccnsfmargtslnorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtslnorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
ccnsfmargtslnorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtslnorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
ccnsfmargtslnorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmargtslnorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  lambda <- object$mlParam[object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 1]
  theta <- object$mlParam[(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + 2):(object$nXvar + object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 1)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B1 <- (1 + lambda)^2 * exp(Wv1)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a1 <- -exp(Wv1/2)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu)) + object$S * epsilon/exp(Wu/2)
  B2 <- (1 + lambda)^2 * exp(Wv2)/(2 * exp(Wu)) + object$S *
    epsilon * (1 + lambda)/exp(Wu/2)
  a2 <- -exp(Wv2/2)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda)/(exp(Wu/2) * (2 * lambda) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# different sigma_u

## logit specification class membership
cmcesfmargtslnorm_Eu_logit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    4 * lambda1 + 2 * lambda1^2)/((1 + lambda1) * (1 + 2 *
    lambda1)), nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    4 * lambda2 + 2 * lambda2^2)/((1 + lambda2) * (1 + 2 *
    lambda2)), nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtslnorm_Vu_logit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    8 * lambda1 + 16 * lambda1^2 + 12 * lambda1^3 + 4 * lambda1^4)/((1 +
    lambda1)^2 * (1 + 2 * lambda1)^2), nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    8 * lambda2 + 16 * lambda2^2 + 12 * lambda2^3 + 4 * lambda2^4)/((1 +
    lambda2)^2 * (1 + 2 * lambda2)^2), nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cauchit specification class membership
cmcesfmargtslnorm_Eu_cauchit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    4 * lambda1 + 2 * lambda1^2)/((1 + lambda1) * (1 + 2 *
    lambda1)), nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    4 * lambda2 + 2 * lambda2^2)/((1 + lambda2) * (1 + 2 *
    lambda2)), nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtslnorm_Vu_cauchit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    8 * lambda1 + 16 * lambda1^2 + 12 * lambda1^3 + 4 * lambda1^4)/((1 +
    lambda1)^2 * (1 + 2 * lambda1)^2), nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    8 * lambda2 + 16 * lambda2^2 + 12 * lambda2^3 + 4 * lambda2^4)/((1 +
    lambda2)^2 * (1 + 2 * lambda2)^2), nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## probit specification class membership
cmcesfmargtslnorm_Eu_probit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    4 * lambda1 + 2 * lambda1^2)/((1 + lambda1) * (1 + 2 *
    lambda1)), nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    4 * lambda2 + 2 * lambda2^2)/((1 + lambda2) * (1 + 2 *
    lambda2)), nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtslnorm_Vu_probit <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    8 * lambda1 + 16 * lambda1^2 + 12 * lambda1^3 + 4 * lambda1^4)/((1 +
    lambda1)^2 * (1 + 2 * lambda1)^2), nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    8 * lambda2 + 16 * lambda2^2 + 12 * lambda2^3 + 4 * lambda2^4)/((1 +
    lambda2)^2 * (1 + 2 * lambda2)^2), nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

## cloglog specification class membership
cmcesfmargtslnorm_Eu_cloglog <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    4 * lambda1 + 2 * lambda1^2)/((1 + lambda1) * (1 + 2 *
    lambda1)), nrow = 1), matrix(exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    4 * lambda2 + 2 * lambda2^2)/((1 + lambda2) * (1 + 2 *
    lambda2)), nrow = 1), matrix(exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmargtslnorm_Vu_cloglog <- function(object) {
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
  lambda1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  lambda2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 2]
  theta <- object$mlParam[(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 3):(object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + object$nZHvar + 2)]
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
  A1 <- exp(Wv1)/(2 * exp(Wu1)) + object$S * epsilon/exp(Wu1/2)
  B1 <- (1 + lambda1)^2 * exp(Wv1)/(2 * exp(Wu1)) + object$S *
    epsilon * (1 + lambda1)/exp(Wu1/2)
  a1 <- -exp(Wv1/2)/exp(Wu1/2) - object$S * epsilon/exp(Wv1/2)
  b1 <- -exp(Wv1/2) * (1 + lambda1)/exp(Wu1/2) - object$S *
    epsilon/exp(Wv1/2)
  A2 <- exp(Wv2)/(2 * exp(Wu2)) + object$S * epsilon/exp(Wu2/2)
  B2 <- (1 + lambda2)^2 * exp(Wv2)/(2 * exp(Wu2)) + object$S *
    epsilon * (1 + lambda2)/exp(Wu2/2)
  a2 <- -exp(Wv2/2)/exp(Wu2/2) - object$S * epsilon/exp(Wv2/2)
  b2 <- -exp(Wv2/2) * (1 + lambda2)/exp(Wu2/2) - object$S *
    epsilon/exp(Wv2/2)
  Pi1 <- (1 + lambda1)/(exp(Wu1/2) * (2 * lambda1) + 1) * (2 *
    exp(A1) * pnorm(a1) - exp(B1) * pnorm(b1))
  Pi2 <- (1 + lambda2)/(exp(Wu2/2) * (2 * lambda2) + 1) * (2 *
    exp(A2) * pnorm(a2) - exp(B2) * pnorm(b2))
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar] * (1 +
    8 * lambda1 + 16 * lambda1^2 + 12 * lambda1^3 + 4 * lambda1^4)/((1 +
    lambda1)^2 * (1 + 2 * lambda1)^2), nrow = 1), matrix(exp(Wu1),
    ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar] * (1 +
    8 * lambda2 + 16 * lambda2^2 + 12 * lambda2^3 + 4 * lambda2^4)/((1 +
    lambda2)^2 * (1 + 2 * lambda2)^2), nrow = 1), matrix(exp(Wu2),
    ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
