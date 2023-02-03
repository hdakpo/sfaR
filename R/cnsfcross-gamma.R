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
# Convolution: gamma - normal                                                  #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for cnsf gamma-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
# same sigma_u

## logit specification class membership
ccnsfgammanormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  if (P < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi2
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cauchit specification class membership
ccnsfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  if (P < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi2
    Probc1 <- 1/pi * atan(Wz) + 1/2
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## probit specification class membership
ccnsfgammanormlike_probit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  if (P < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi2
    Probc1 <- pnorm(Wz)
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cloglog specification class membership
ccnsfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  if (P < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
      exp(-P * Wu/2)/gamma(P) * Hi2
    Probc1 <- 1 - exp(-exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}



## logit specification class membership
cmcesfgammanormlike_logit <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  if (P1 < 0 || P2 < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
      exp(-P1 * Wu1/2)/gamma(P1) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
      exp(-P2 * Wu2/2)/gamma(P2) * Hi2
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cauchit specification class membership
cmcesfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  if (P1 < 0 || P2 < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
      exp(-P1 * Wu1/2)/gamma(P1) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
      exp(-P2 * Wu2/2)/gamma(P2) * Hi2
    Probc1 <- 1/pi * atan(Wz) + 1/2
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## probit specification class membership
cmcesfgammanormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  if (P1 < 0 || P2 < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
      exp(-P1 * Wu1/2)/gamma(P1) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
      exp(-P2 * Wu2/2)/gamma(P2) * Hi2
    Probc1 <- pnorm(Wz)
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

## cloglog specification class membership
cmcesfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1 <- numeric(N)
  Hi2 <- numeric(N)
  for (i in 1:N) {
    Hi1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  if (P1 < 0 || P2 < 0) {
    return(NA)
  } else {
    Pi1 <- exp(S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 * exp(Wu1))) *
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
      exp(-P1 * Wu1/2)/gamma(P1) * Hi1
    Pi2 <- exp(S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 * exp(Wu2))) *
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
      exp(-P2 * Wu2/2)/gamma(P2) * Hi2
    Probc1 <- 1 - exp(-exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA),
      return(wHvar * log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# starting value for the log-likelihood ----------
#' starting values for cnsf gamma-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @param printInfo logical print info during optimization
#' @param whichStart strategy to get starting values
#' @param initIter maximum iterations for initialization
#' @param initAlg algorithm for maxLik  
#' @param tol parameter tolerance
#' @noRd
# same sigma_u
cstcnsfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGamma <- NULL
  } else {
    cat("Initialization: SFA + gamma - normal distributions...\n")
    initGamma <- maxLik::maxLik(logLik = cgammanormlike,
      start = cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradgammanormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE],
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
      FiMat = FiMat)
    Esti <- initGamma$estimate
  }
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_",
    colnames(vHvar)), "P", paste0("CN_", colnames(Zvar)))
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# different sigma_u
cstmcesfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat, whichStart, initIter, initAlg, printInfo, tol) {
  if (whichStart == 1L) {
    Esti <- cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
      S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
      nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE])
    initGamma <- NULL
  } else {
    cat("Initialization: SFA + gamma - normal distributions...\n")
    initGamma <- maxLik::maxLik(logLik = cgammanormlike,
      start = cstgammanorm(olsObj = olsObj, epsiRes = epsiRes,
        S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE],
        nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]),
      grad = cgradgammanormlike, method = initAlg, control = list(iterlim = initIter,
        printLevel = if (printInfo) 2 else 0, reltol = tol),
      nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = uHvar[,
        1, drop = FALSE], vHvar = vHvar[, 1, drop = FALSE],
      Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N,
      FiMat = FiMat)
    Esti <- initGamma$estimate
  }
  StartVal <- c(Esti[1:nXvar], 0.95 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar >
    1) rep(0, nuZUvar - 1), 0.95 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], if (nvZVvar >
    1) rep(0, nvZVvar - 1), 0.95 * Esti[nXvar + 3], 1.05 *
    Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_",
    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("Zv_", colnames(vHvar)), "P",
    "P", paste0("MCE_", colnames(Zvar)))
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Gradient of the likelihood function ----------
#' gradient for cnsf gamma-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
#' @noRd
# same sigma_u

## logit specification class membership
cgradcnsfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wusi1 <- (ewv1/ewusr + S * (epsilon))
  wusi2 <- (ewv2/ewusr + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  evur2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P - 1), 1, sum)
  sigx1 <- (prC * evur2 * pepsi2 * sDiv2 + evur1 * ewz * pepsi1 *
    sDiv1/wzdeno)
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 *
    ewu)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv2/(2 *
    ewu)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu) - (0.5 * (ewvsr1/ewusr) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu) - (0.5 * (ewvsr2/ewusr) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (evur1 * pepsi1 * sDiv1 - evur2 * pepsi2 * sDiv2)
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv2/ewusr, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * ewz/(wzdeno * sigx1)
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prC * evur2/sigx1
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * prC * evur2/sigx1, FUN = "*") + gx2 + sweep(Xvar,
    MARGIN = 1, STATS = S * dpurv1 * sDiv1 * evur1 * ewz/(wzdeno *
      sigx1), FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * evur1 * ewz/(sigx1 *
      wzdeno)
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * prC * evur2/sigx1
  }
  gu <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    prC * evur2/sigx1, FUN = "*") + gu2 + sweep(uHvar, MARGIN = 1,
    STATS = sigx3_1 * sDiv1 * evur1 * ewz/(sigx1 * wzdeno),
    FUN = "*") - 0.5 * P * uHvar
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * ewz/(sigx1 * wzdeno)
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prC * evur2/sigx1
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    evur1 * ewz/(sigx1 * wzdeno), FUN = "*") + gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    prC * evur2/sigx1, FUN = "*") + gv2
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2^(P - 1) * log(F6_2),
    1, sum) * prC * evur2 * pepsi2/sigx1 + apply(F6_1^(P -
    1) * log(F6_1), 1, sum) * evur1 * ewz * pepsi1/(sigx1 *
    wzdeno) - (0.5 * (Wu) + digamma(P)), sweep(Zvar, MARGIN = 1,
    STATS = prC * sigx5 * ewz/(sigx1 * wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradcnsfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  wusi1 <- (ewv1/ewusr + S * (epsilon))
  wusi2 <- (ewv2/ewusr + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  evur2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P - 1), 1, sum)
  sigx1 <- (ewz2 * evur2 * pepsi2 * sDiv2 + ewz1 * evur1 *
    pepsi1 * sDiv1)
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 *
    ewu)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv2/(2 *
    ewu)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu) - (0.5 * (ewvsr1/ewusr) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu) - (0.5 * (ewvsr2/ewusr) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (evur1 * pepsi1 * sDiv1 - evur2 * pepsi2 * sDiv2)
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv2/ewusr, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * ewz1/sigx1
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * ewz2 * evur2/sigx1
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * ewz2 * evur2/sigx1, FUN = "*") + gx2 + sweep(Xvar,
    MARGIN = 1, STATS = S * dpurv1 * sDiv1 * evur1 * ewz1/sigx1,
    FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * evur1 * ewz1/sigx1
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * ewz2 * evur2/sigx1
  }
  gu <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    ewz2 * evur2/sigx1, FUN = "*") + gu2 + sweep(uHvar, MARGIN = 1,
    STATS = sigx3_1 * sDiv1 * evur1 * ewz1/sigx1, FUN = "*") -
    0.5 * P * uHvar
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * ewz1/sigx1
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * ewz2 * evur2/sigx1
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    evur1 * ewz1/sigx1, FUN = "*") + gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    ewz2 * evur2/sigx1, FUN = "*") + gv2
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2^(P - 1) * log(F6_2),
    1, sum) * ewz2 * evur2 * pepsi2/sigx1 + apply(F6_1^(P -
    1) * log(F6_1), 1, sum) * evur1 * ewz1 * pepsi1/sigx1 -
    (0.5 * (Wu) + digamma(P)), sweep(Zvar, MARGIN = 1, STATS = sigx5/(pi *
    sigx1 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradcnsfgammanormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  wusi1 <- (ewv1/ewusr + S * (epsilon))
  wusi2 <- (ewv2/ewusr + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  evur2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P - 1), 1, sum)
  sigx1 <- ((1 - pwZ) * evur2 * pepsi2 * sDiv2 + evur1 * pepsi1 *
    pwZ * sDiv1)
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 *
    ewu)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv2/(2 *
    ewu)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu) - (0.5 * (ewvsr1/ewusr) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu) - (0.5 * (ewvsr2/ewusr) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (evur1 * pepsi1 * sDiv1 - evur2 * pepsi2 * sDiv2)
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv2/ewusr, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * pwZ/sigx1
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * (1 - pwZ) * evur2/sigx1
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * (1 - pwZ) * evur2/sigx1, FUN = "*") + gx2 + sweep(Xvar,
    MARGIN = 1, STATS = S * dpurv1 * sDiv1 * evur1 * pwZ/sigx1,
    FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * evur1 * pwZ/sigx1
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * (1 - pwZ) * evur2/sigx1
  }
  gu <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    (1 - pwZ) * evur2/sigx1, FUN = "*") + gu2 + sweep(uHvar,
    MARGIN = 1, STATS = sigx3_1 * sDiv1 * evur1 * pwZ/sigx1,
    FUN = "*") - 0.5 * P * uHvar
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * pwZ/sigx1
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * (1 - pwZ) * evur2/sigx1
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    evur1 * pwZ/sigx1, FUN = "*") + gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    (1 - pwZ) * evur2/sigx1, FUN = "*") + gv2
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2^(P - 1) * log(F6_2),
    1, sum) * (1 - pwZ) * evur2 * pepsi2/sigx1 + apply(F6_1^(P -
    1) * log(F6_1), 1, sum) * evur1 * pwZ * pepsi1/sigx1 -
    (0.5 * (Wu) + digamma(P)), sweep(Zvar, MARGIN = 1, STATS = dwZ *
    sigx5/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradcnsfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  phi2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(nXvar + nuZUvar +
    2 * nvZVvar)]
  P <- parm[nXvar + nuZUvar + 2 * nvZVvar + 1]
  theta <- parm[(nXvar + nuZUvar + 2 * nvZVvar + 2):(nXvar +
    nuZUvar + 2 * nvZVvar + nZHvar + 1)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu <- exp(Wu)
  ewusr <- exp(Wu/2)
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  wusi1 <- (ewv1/ewusr + S * (epsilon))
  wusi2 <- (ewv2/ewusr + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu) + S * (epsilon)/ewusr)
  evur2 <- exp(ewv2/(2 * ewu) + S * (epsilon)/ewusr)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P - 1), 1, sum)
  sigx1 <- ((1 - prZ) * evur1 * pepsi1 * sDiv1 + prZ * evur2 *
    pepsi2 * sDiv2)
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv1/(2 *
    ewu)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr) + 2 * (ewu * ewv2/(2 *
    ewu)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu) - (0.5 * (ewvsr1/ewusr) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu) - (0.5 * (ewvsr2/ewusr) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (evur1 * pepsi1 * sDiv1 - evur2 * pepsi2 * sDiv2)
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv1/ewusr, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P - 2) * (P -
    1), MARGIN = 1, STATS = ewv2/ewusr, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * (1 - prZ)/sigx1
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prZ * evur2/sigx1
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * prZ * evur2/sigx1, FUN = "*") + gx2 + sweep(Xvar,
    MARGIN = 1, STATS = S * dpurv1 * sDiv1 * evur1 * (1 -
      prZ)/sigx1, FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * evur1 * (1 - prZ)/sigx1
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * prZ * evur2/sigx1
  }
  gu <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    prZ * evur2/sigx1, FUN = "*") + gu2 + sweep(uHvar, MARGIN = 1,
    STATS = sigx3_1 * sDiv1 * evur1 * (1 - prZ)/sigx1, FUN = "*") -
    0.5 * P * uHvar
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * evur1 * (1 - prZ)/sigx1
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P - 2) *
      (P - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prZ * evur2/sigx1
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    evur1 * (1 - prZ)/sigx1, FUN = "*") + gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    prZ * evur2/sigx1, FUN = "*") + gv2
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2^(P - 1) * log(F6_2),
    1, sum) * prZ * evur2 * pepsi2/sigx1 + apply(F6_1^(P -
    1) * log(F6_1), 1, sum) * evur1 * (1 - prZ) * pepsi1/sigx1 -
    (0.5 * (Wu) + digamma(P)), sweep(Zvar, MARGIN = 1, STATS = prZ *
    sigx5 * ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# different sigma_u

## logit specification class membership
cgradmcesfgammanormlike_logit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewusrp1 <- exp(-(Wu1 * P1/2))
  ewusrp2 <- exp(-(Wu2 * P2/2))
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1 + ewz)
  prC <- (1 - ewz/wzdeno)
  wusi1 <- (ewv1/ewusr1 + S * (epsilon))
  wusi2 <- (ewv2/ewusr2 + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  evur2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr1)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr2)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr1 - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr2 - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P1 - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P2 - 1), 1, sum)
  sigx1 <- (prC * ewusrp2 * evur2 * pepsi2 * sDiv2/gamma(P2) +
    ewusrp1 * evur1 * ewz * pepsi1 * sDiv1/(wzdeno * gamma(P1)))
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr1) + 0.5 * P1 + 2 *
    (ewu1 * ewv1/(2 * ewu1)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr2) + 0.5 * P2 + 2 *
    (ewu2 * ewv2/(2 * ewu2)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu1) - (0.5 * (ewvsr1/ewusr1) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu2) - (0.5 * (ewvsr2/ewusr2) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- ((1/(wzdeno * gamma(P1)) - ewz * gamma(P1)/(wzdeno *
    gamma(P1))^2) * ewusrp1 * evur1 * pepsi1 * sDiv1 - prC *
    ewusrp2 * evur2 * pepsi2 * sDiv2/(wzdeno * gamma(P2)))
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P1 - 2) * (P1 -
    1), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P2 - 2) * (P2 -
    1), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * ewusrp1 * evur1 * ewz/(wzdeno *
      gamma(P1) * sigx1)
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prC * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * prC * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gx2 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv1 * sDiv1 *
    ewusrp1 * evur1 * ewz/(wzdeno * gamma(P1) * sigx1), FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * ewusrp1 * evur1 *
      ewz/(wzdeno * gamma(P1) * sigx1)
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * prC * ewusrp2 *
      evur2/(gamma(P2) * sigx1)
  }
  gu1 <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_1 * sDiv1 *
    ewusrp1 * evur1 * ewz/(wzdeno * gamma(P1) * sigx1), FUN = "*")
  gu2 <- gu2 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    prC * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*")
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * ewusrp1 * evur1 * ewz/(wzdeno *
      gamma(P1) * sigx1)
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prC * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    ewusrp1 * evur1 * ewz/(wzdeno * gamma(P1) * sigx1), FUN = "*") +
    gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    prC * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gv2
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, ((apply(F6_1^(P1 -
    1) * log(F6_1), 1, sum) - 0.5 * (Wu1 * sDiv1))/(wzdeno *
    gamma(P1)) - wzdeno * digamma(P1) * gamma(P1) * sDiv1/(wzdeno *
    gamma(P1))^2) * ewusrp1 * evur1 * ewz * pepsi1/sigx1,
    (apply(F6_2^(P2 - 1) * log(F6_2), 1, sum) - (0.5 * (Wu2) +
      digamma(P2)) * sDiv2) * prC * ewusrp2 * evur2 * pepsi2/(sigx1 *
      gamma(P2)), sweep(Zvar, MARGIN = 1, STATS = sigx5 *
      ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cauchit specification class membership
cgradmcesfgammanormlike_cauchit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewusrp1 <- exp(-(Wu1 * P1/2))
  ewusrp2 <- exp(-(Wu2 * P2/2))
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz1 <- (0.5 + atan(Wz)/pi)
  ewz2 <- (0.5 - atan(Wz)/pi)
  wusi1 <- (ewv1/ewusr1 + S * (epsilon))
  wusi2 <- (ewv2/ewusr2 + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  evur2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr1)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr2)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr1 - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr2 - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P1 - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P2 - 1), 1, sum)
  sigx1 <- (ewz2 * ewusrp2 * evur2 * pepsi2 * sDiv2/gamma(P2) +
    ewz1 * ewusrp1 * evur1 * pepsi1 * sDiv1/gamma(P1))
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr1) + 0.5 * P1 + 2 *
    (ewu1 * ewv1/(2 * ewu1)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr2) + 0.5 * P2 + 2 *
    (ewu2 * ewv2/(2 * ewu2)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu1) - (0.5 * (ewvsr1/ewusr1) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu2) - (0.5 * (ewvsr2/ewusr2) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (ewusrp1 * evur1 * pepsi1 * sDiv1/gamma(P1) - ewusrp2 *
    evur2 * pepsi2 * sDiv2/gamma(P2))
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P1 - 2) * (P1 -
    1), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P2 - 2) * (P2 -
    1), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * ewz1 * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * ewz2 * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * ewz2 * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gx2 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv1 * sDiv1 *
    ewz1 * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * ewz1 * ewusrp1 *
      evur1/(gamma(P1) * sigx1)
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * ewz2 * ewusrp2 *
      evur2/(gamma(P2) * sigx1)
  }
  gu1 <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_1 * sDiv1 *
    ewz1 * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*")
  gu2 <- gu2 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    ewz2 * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*")
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * ewz1 * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * ewz2 * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    ewz1 * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*") +
    gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    ewz2 * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gv2
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, (apply(F6_1^(P1 -
    1) * log(F6_1), 1, sum) - (0.5 * (Wu1) + digamma(P1)) *
    sDiv1) * ewz1 * ewusrp1 * evur1 * pepsi1/(sigx1 * gamma(P1)),
    (apply(F6_2^(P2 - 1) * log(F6_2), 1, sum) - (0.5 * (Wu2) +
      digamma(P2)) * sDiv2) * ewz2 * ewusrp2 * evur2 *
      pepsi2/(sigx1 * gamma(P2)), sweep(Zvar, MARGIN = 1,
      STATS = sigx5/(pi * sigx1 * ((Wz)^2 + 1)), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## probit specification class membership
cgradmcesfgammanormlike_probit <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewusrp1 <- exp(-(Wu1 * P1/2))
  ewusrp2 <- exp(-(Wu2 * P2/2))
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  pwZ <- pnorm(Wz)
  dwZ <- dnorm(Wz, 0, 1)
  wusi1 <- (ewv1/ewusr1 + S * (epsilon))
  wusi2 <- (ewv2/ewusr2 + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  evur2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr1)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr2)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr1 - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr2 - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P1 - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P2 - 1), 1, sum)
  sigx1 <- ((1 - pwZ) * ewusrp2 * evur2 * pepsi2 * sDiv2/gamma(P2) +
    ewusrp1 * evur1 * pepsi1 * pwZ * sDiv1/gamma(P1))
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr1) + 0.5 * P1 + 2 *
    (ewu1 * ewv1/(2 * ewu1)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr2) + 0.5 * P2 + 2 *
    (ewu2 * ewv2/(2 * ewu2)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu1) - (0.5 * (ewvsr1/ewusr1) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu2) - (0.5 * (ewvsr2/ewusr2) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (ewusrp1 * evur1 * pepsi1 * sDiv1/gamma(P1) - ewusrp2 *
    evur2 * pepsi2 * sDiv2/gamma(P2))
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P1 - 2) * (P1 -
    1), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P2 - 2) * (P2 -
    1), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * pwZ * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * (1 - pwZ) * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * (1 - pwZ) * ewusrp2 * evur2/(gamma(P2) * sigx1),
    FUN = "*") + gx2 + sweep(Xvar, MARGIN = 1, STATS = S *
    dpurv1 * sDiv1 * pwZ * ewusrp1 * evur1/(gamma(P1) * sigx1),
    FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * pwZ * ewusrp1 *
      evur1/(gamma(P1) * sigx1)
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * (1 - pwZ) * ewusrp2 *
      evur2/(gamma(P2) * sigx1)
  }
  gu1 <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_1 * sDiv1 *
    pwZ * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*")
  gu2 <- gu2 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    (1 - pwZ) * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*")
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * pwZ * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * (1 - pwZ) * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    pwZ * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*") +
    gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    (1 - pwZ) * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gv2
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, (apply(F6_1^(P1 -
    1) * log(F6_1), 1, sum) - (0.5 * (Wu1) + digamma(P1)) *
    sDiv1) * pwZ * ewusrp1 * evur1 * pepsi1/(sigx1 * gamma(P1)),
    (apply(F6_2^(P2 - 1) * log(F6_2), 1, sum) - (0.5 * (Wu2) +
      digamma(P2)) * sDiv2) * (1 - pwZ) * ewusrp2 * evur2 *
      pepsi2/(sigx1 * gamma(P2)), sweep(Zvar, MARGIN = 1,
      STATS = dwZ * sigx5/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## cloglog specification class membership
cgradmcesfgammanormlike_cloglog <- function(parm, nXvar, nuZUvar,
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar,
  N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar +
    nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar +
    2 * nuZUvar + 2 * nvZVvar)]
  P1 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 1]
  P2 <- parm[nXvar + 2 * nuZUvar + 2 * nvZVvar + 2]
  theta <- parm[(nXvar + 2 * nuZUvar + 2 * nvZVvar + 3):(nXvar +
    2 * nuZUvar + 2 * nvZVvar + nZHvar + 2)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu1 <- exp(Wu1)
  ewu2 <- exp(Wu2)
  ewusr1 <- exp(Wu1/2)
  ewusr2 <- exp(Wu2/2)
  ewusrp1 <- exp(-(Wu1 * P1/2))
  ewusrp2 <- exp(-(Wu2 * P2/2))
  ewvsr1 <- exp(Wv1/2)
  ewvsr2 <- exp(Wv2/2)
  ewz <- exp(Wz)
  prZ <- exp(-ewz)
  wusi1 <- (ewv1/ewusr1 + S * (epsilon))
  wusi2 <- (ewv2/ewusr2 + S * (epsilon))
  pmusig1 <- pnorm(wusi1/ewvsr1)
  pmusig2 <- pnorm(wusi2/ewvsr2)
  dmusig1 <- dnorm(wusi1/ewvsr1, 0, 1)
  dmusig2 <- dnorm(wusi2/ewvsr2, 0, 1)
  F1_1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") +
    FiMat
  F1_2 <- sweep((1 - FiMat), MARGIN = 1, STATS = pmusig2, FUN = "*") +
    FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewvsr1/ewusr1 + S * (epsilon)/ewvsr1)
  epsi2 <- (ewvsr2/ewusr2 + S * (epsilon)/ewvsr2)
  depsi1 <- dnorm(-epsi1, 0, 1)
  pepsi1 <- pnorm(-epsi1)
  depsi2 <- dnorm(-epsi2, 0, 1)
  pepsi2 <- pnorm(-epsi2)
  evur1 <- exp(ewv1/(2 * ewu1) + S * (epsilon)/ewusr1)
  evur2 <- exp(ewv2/(2 * ewu2) + S * (epsilon)/ewusr2)
  dpurv1 <- (depsi1/ewvsr1 - pepsi1/ewusr1)
  dpurv2 <- (depsi2/ewvsr2 - pepsi2/ewusr2)
  F4_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1,
    FUN = "*")
  F4_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2,
    FUN = "*")
  F5_1 <- sweep((1 - FiMat)/F3_1, MARGIN = 1, STATS = dmusig1 *
    (ewv1/ewusr1 - 0.5 * wusi1), FUN = "*")
  F5_2 <- sweep((1 - FiMat)/F3_2, MARGIN = 1, STATS = dmusig2 *
    (ewv2/ewusr2 - 0.5 * wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewvsr1), FUN = "*"),
    MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewvsr2), FUN = "*"),
    MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P1 - 1), 1, sum)
  sDiv2 <- apply(F6_2^(P2 - 1), 1, sum)
  sigx1 <- ((1 - prZ) * ewusrp1 * evur1 * pepsi1 * sDiv1/gamma(P1) +
    ewusrp2 * prZ * evur2 * pepsi2 * sDiv2/gamma(P2))
  sigx2_1 <- (0.5 * (S * (epsilon)/ewusr1) + 0.5 * P1 + 2 *
    (ewu1 * ewv1/(2 * ewu1)^2))
  sigx2_2 <- (0.5 * (S * (epsilon)/ewusr2) + 0.5 * P2 + 2 *
    (ewu2 * ewv2/(2 * ewu2)^2))
  sigx3_1 <- (0.5 * (depsi1 * ewvsr1/ewusr1) - sigx2_1 * pepsi1)
  sigx3_2 <- (0.5 * (depsi2 * ewvsr2/ewusr2) - sigx2_2 * pepsi2)
  sigx4_1 <- (ewv1 * pepsi1/(2 * ewu1) - (0.5 * (ewvsr1/ewusr1) -
    0.5 * (S * (epsilon)/ewvsr1)) * depsi1)
  sigx4_2 <- (ewv2 * pepsi2/(2 * ewu2) - (0.5 * (ewvsr2/ewusr2) -
    0.5 * (S * (epsilon)/ewvsr2)) * depsi2)
  sigx5 <- (ewusrp1 * evur1 * pepsi1 * sDiv1/gamma(P1) - ewusrp2 *
    evur2 * pepsi2 * sDiv2/gamma(P2))
  F7_1 <- sweep((0.5 - 0.5 * (F4_1)) * F6_1^(P1 - 2) * (P1 -
    1), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "*")
  F7_2 <- sweep((0.5 - 0.5 * (F4_2)) * F6_2^(P2 - 2) * (P2 -
    1), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5 * (ewvsr1),
    FUN = "*"), MARGIN = 1, STATS = ewv1/ewusr1, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5 * (ewvsr2),
    FUN = "*"), MARGIN = 1, STATS = ewv2/ewusr2, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S * (1 - F4_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi1 * (1 - prZ) * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gx2[, k] <- apply(sweep(S * (1 - F4_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prZ * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gx <- gx1 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv2 *
    sDiv2 * prZ * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gx2 + sweep(Xvar, MARGIN = 1, STATS = S * dpurv1 * sDiv1 *
    (1 - prZ) * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*")
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi1 * (1 - prZ) * ewusrp1 *
      evur1/(gamma(P1) * sigx1)
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum) * pepsi2 * prZ * ewusrp2 *
      evur2/(gamma(P2) * sigx1)
  }
  gu1 <- gu1 + sweep(uHvar, MARGIN = 1, STATS = sigx3_1 * sDiv1 *
    (1 - prZ) * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*")
  gu2 <- gu2 + sweep(uHvar, MARGIN = 1, STATS = sigx3_2 * sDiv2 *
    prZ * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*")
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1 + F8_1) * F6_1^(P1 - 2) *
      (P1 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi1 * (1 - prZ) * ewusrp1 * evur1/(gamma(P1) *
      sigx1)
    gv2[, k] <- apply(sweep((F5_2 + F8_2) * F6_2^(P2 - 2) *
      (P2 - 1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"),
      1, sum) * pepsi2 * prZ * ewusrp2 * evur2/(gamma(P2) *
      sigx1)
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_1 * sDiv1 *
    (1 - prZ) * ewusrp1 * evur1/(gamma(P1) * sigx1), FUN = "*") +
    gv1
  gv2 <- sweep(vHvar, MARGIN = 1, STATS = sigx4_2 * sDiv2 *
    prZ * ewusrp2 * evur2/(gamma(P2) * sigx1), FUN = "*") +
    gv2
  gradll <- cbind(gx, gu1, gu2, gv1, gv2, (apply(F6_1^(P1 -
    1) * log(F6_1), 1, sum) - (0.5 * (Wu1) + digamma(P1)) *
    sDiv1) * (1 - prZ) * ewusrp1 * evur1 * pepsi1/(sigx1 *
    gamma(P1)), (apply(F6_2^(P2 - 1) * log(F6_2), 1, sum) -
    (0.5 * (Wu2) + digamma(P2)) * sDiv2) * prZ * ewusrp2 *
    evur2 * pepsi2/(sigx1 * gamma(P2)), sweep(Zvar, MARGIN = 1,
    STATS = sigx5 * prZ * ewz/sigx1, FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Optimization using different algorithms ----------
#' optimizations solve for cnsf gamma-normal distribution
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
#' @param N number of observations
#' @param FiMat matrix of random draws
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
cnsfgammanormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfgammanormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfgammanormlike_logit,
      grad = cgradcnsfgammanormlike_logit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), parm), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ccnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgammanormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$solution)
  }
  mleObj$logL_OBS <- ccnsfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
cnsfgammanormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfgammanormlike_cauchit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfgammanormlike_cauchit,
      grad = cgradcnsfgammanormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ccnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgammanormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## probit specification class membership
cnsfgammanormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfgammanormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfgammanormlike_probit,
      grad = cgradcnsfgammanormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ccnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgammanormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cloglog specification class membership
cnsfgammanormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstcnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(ccnsfgammanormlike_cloglog(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("CNSF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(ccnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ccnsfgammanormlike_cloglog,
      grad = cgradcnsfgammanormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(ccnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(ccnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ccnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradcnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgammanormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradcnsfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- ccnsfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

# different sigma_u

## logit specification class membership
mcesfgammanormAlgOpt_logit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initWeibull <- start_st$initWeibull
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfgammanormlike_logit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfgammanormlike_logit,
      grad = cgradmcesfgammanormlike_logit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), parm), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfgammanormlike_logit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_logit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), mleObj$solution)
  }
  mleObj$logL_OBS <- cmcesfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfgammanormlike_logit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initWeibull = initWeibull))
}

## cauchit specification class membership
mcesfgammanormAlgOpt_cauchit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfgammanormlike_cauchit(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfgammanormlike_cauchit,
      grad = cgradmcesfgammanormlike_cauchit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfgammanormlike_cauchit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_cauchit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfgammanormlike_cauchit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## probit specification class membership
mcesfgammanormAlgOpt_probit <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfgammanormlike_probit(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar,
    N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfgammanormlike_probit,
      grad = cgradmcesfgammanormlike_probit, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfgammanormlike_probit(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_probit(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfgammanormlike_probit(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

## cloglog specification class membership
mcesfgammanormAlgOpt_cloglog <- function(start, olsParam, dataTable,
  S, wHvar, nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Zvar, nZHvar, Yvar, Xvar, method, printInfo, itermax, stepmax,
  whichStart, initIter, initAlg, tol, gradtol, hessianType,
  qac) {
  start_st <- if (!is.null(start))
    start else cstmcesfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nXvar = nXvar, Xvar = Xvar, Yvar = Yvar,
    whichStart = whichStart, initIter = initIter, initAlg = initAlg,
    printInfo = printInfo, tol = tol)
  initGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  startLoglik <- sum(cmcesfgammanormlike_cloglog(startVal,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  cat("MCESF Estimation...\n")
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(cmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hessian = if (hessianType != 2) 1 else 0,
    control = list(trace = printInfo, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cmcesfgammanormlike_cloglog,
      grad = cgradmcesfgammanormlike_cloglog, start = startVal,
      finalHessian = if (hessianType == 2) "bhhh" else TRUE,
      control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar), sr1 = trustOptim::trust.optim(x = startVal,
      fn = function(parm) -sum(cmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), method = "SR1", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 4L else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(cmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), hs = as(calculus::jacobian(function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), unname(parm)), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(cmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), control = list(iter.max = itermax,
        trace = printInfo, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradmcesfgammanormlike_cloglog(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
      nZHvar = nZHvar))
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
    if (method == "ucminf")
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$par))
    if (method == "sr1")
      mleObj$hessian <- calculus::jacobian(function(parm) colSums(cgradmcesfgammanormlike_cloglog(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
        nZHvar = nZHvar)), unname(mleObj$solution))
  }
  mleObj$logL_OBS <- cmcesfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradmcesfgammanormlike_cloglog(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    S = S, wHvar = wHvar, N = N, FiMat = FiMat, Zvar = Zvar,
    nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam, initGamma = initGamma))
}

# Conditional efficiencies estimation ----------
#' efficiencies for cnsf gamma-normal distribution
#' @param object object of class sfacross
#' @param level level for confidence interval
#' @noRd
# Same sigma_u

## logit specification class membership
ccnsfgammanormeff_logit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
ccnsfgammanormeff_cauchit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
ccnsfgammanormeff_probit <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
ccnsfgammanormeff_cloglog <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Different sigma_u

## logit specification class membership
cmcesfgammanormeff_logit <- function(object, level) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu1/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu2/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu1/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu2/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cauchit specification class membership
cmcesfgammanormeff_cauchit <- function(object, level) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu1/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu2/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu1/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu2/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## probit specification class membership
cmcesfgammanormeff_probit <- function(object, level) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu1/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu2/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu1/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu2/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

## cloglog specification class membership
cmcesfgammanormeff_cloglog <- function(object, level) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c1 <- exp(-u_c1)
    teJLMS_c2 <- exp(-u_c2)
    teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
    mui_Gi1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) -
      exp(Wv1)
    mui_Ki1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1)) +
      exp(Wv1)
    mui_Gi2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) -
      exp(Wv2)
    mui_Ki2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2)) +
      exp(Wv2)
    Gi1 <- numeric(object$Nobs)
    Ki1 <- numeric(object$Nobs)
    Gi2 <- numeric(object$Nobs)
    Ki2 <- numeric(object$Nobs)
    for (i in 1:object$Nobs) {
      Gi1[i] <- mean((mui_Gi1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Ki1[i] <- mean((mui_Ki1[i] + sqrt(exp(Wv1[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki1[i]/sqrt(exp(Wv1[i])))))^(P1 -
        1))
      Gi2[i] <- mean((mui_Gi2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Gi2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
      Ki2[i] <- mean((mui_Ki2[i] + sqrt(exp(Wv2[i])) *
        qnorm(object$FiMat[i, ] + (1 - object$FiMat[i,
          ]) * pnorm(-mui_Ki2[i]/sqrt(exp(Wv2[i])))))^(P2 -
        1))
    }
    teBC_c1 <- exp(exp(Wv1)/exp(Wu1/2) + object$S * epsilon +
      exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) - object$S *
      epsilon/exp(Wv1/2) - exp(Wv1/2)) * Gi1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_c2 <- exp(exp(Wv2)/exp(Wu2/2) + object$S * epsilon +
      exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) - object$S *
      epsilon/exp(Wv2/2) - exp(Wv2/2)) * Gi2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_c <- ifelse(Group_c == 1, teBC_c1, teBC_c2)
    effBC_c1 <- ifelse(Group_c == 1, teBC_c1, NA)
    effBC_c2 <- ifelse(Group_c == 2, teBC_c2, NA)
    teBC_reciprocal_c1 <- exp(-exp(Wv1)/exp(Wu1/2) - object$S *
      epsilon + exp(Wv1)/2) * pnorm(-exp(Wv1/2 - Wu1/2) -
      object$S * epsilon/exp(Wv1/2) + exp(Wv1/2)) * Ki1/(pnorm(-exp(Wv1/2 -
      Wu1/2) - object$S * epsilon/exp(Wv1/2)) * Hi2_1)
    teBC_reciprocal_c2 <- exp(-exp(Wv2)/exp(Wu2/2) - object$S *
      epsilon + exp(Wv2)/2) * pnorm(-exp(Wv2/2 - Wu2/2) -
      object$S * epsilon/exp(Wv2/2) + exp(Wv2/2)) * Ki2/(pnorm(-exp(Wv2/2 -
      Wu2/2) - object$S * epsilon/exp(Wv2/2)) * Hi2_2)
    teBC_reciprocal_c <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      teBC_reciprocal_c2)
    ReffBC_c1 <- ifelse(Group_c == 1, teBC_reciprocal_c1,
      NA)
    ReffBC_c2 <- ifelse(Group_c == 2, teBC_reciprocal_c2,
      NA)
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, teJLMS_c = teJLMS_c,
      teBC_c = teBC_c, teBC_reciprocal_c = teBC_reciprocal_c,
      PosteriorProb_c1 = Pcond_c1, PriorProb_c1 = Probc1,
      u_c1 = u_c1, teBC_c1 = teBC_c1, teBC_reciprocal_c1 = teBC_reciprocal_c1,
      PosteriorProb_c2 = Pcond_c2, PriorProb_c2 = Probc2,
      u_c2 = u_c2, teBC_c2 = teBC_c2, teBC_reciprocal_c2 = teBC_reciprocal_c2,
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, effBC_c1 = effBC_c1,
      effBC_c2 = effBC_c2, ReffBC_c1 = ReffBC_c1, ReffBC_c2 = ReffBC_c2)
  } else {
    res <- data.frame(Group_c = Group_c, PosteriorProb_c = P_cond_c,
      odRatio = odRatio, u_c = u_c, PosteriorProb_c1 = Pcond_c1,
      PriorProb_c1 = Probc1, u_c1 = u_c1, PosteriorProb_c2 = Pcond_c2,
      PriorProb_c2 = Probc2, u_c2 = u_c2, ineff_c1 = ineff_c1,
      ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for cnsf gamma-normal distribution
#' @param object object of class sfacross
#' @noRd
# Same sigma_u

## logit specification class membership
ccnsfmarggammanorm_Eu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarggammanorm_Vu_logit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
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
ccnsfmarggammanorm_Eu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarggammanorm_Vu_cauchit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
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
ccnsfmarggammanorm_Eu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarggammanorm_Vu_probit <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
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
ccnsfmarggammanorm_Eu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

ccnsfmarggammanorm_Vu_cloglog <- function(object) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  phi2 <- object$mlParam[(object$nXvar + object$nuZUvar + object$nvZVvar +
    1):(object$nXvar + object$nuZUvar + 2 * object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + 2 * object$nvZVvar +
    1]
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
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) *
    pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) *
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

# Different sigma_u

## logit specification class membership
cmcesfmarggammanorm_Eu_logit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1/2 * exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2/2 * exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarggammanorm_Vu_logit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1 * exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2 * exp(Wu2), ncol = 1))
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
cmcesfmarggammanorm_Eu_cauchit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1/2 * exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2/2 * exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarggammanorm_Vu_cauchit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1/pi * atan(Wz) + 1/2
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1 * exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2 * exp(Wu2), ncol = 1))
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
cmcesfmarggammanorm_Eu_probit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1/2 * exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2/2 * exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarggammanorm_Vu_probit <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- pnorm(Wz)
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1 * exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2 * exp(Wu2), ncol = 1))
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
cmcesfmarggammanorm_Eu_cloglog <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1/2 * exp(Wu1/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2/2 * exp(Wu2/2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}

cmcesfmarggammanorm_Vu_cloglog <- function(object) {
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
  P1 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
    2 * object$nvZVvar + 1]
  P2 <- object$mlParam[object$nXvar + 2 * object$nuZUvar +
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
  mui1 <- -object$S * epsilon - exp(Wv1)/sqrt(exp(Wu1))
  mui2 <- -object$S * epsilon - exp(Wv2)/sqrt(exp(Wu2))
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1))
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P1 -
      1))
    Hi1_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P2 -
      1))
  }
  Pi1 <- exp(object$S * epsilon/exp(Wu1/2) + exp(Wv1)/(2 *
    exp(Wu1))) * pnorm(-object$S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu1/2)) *
    exp(-P1 * Wu1/2)/gamma(P1) * Hi2_1
  Pi2 <- exp(object$S * epsilon/exp(Wu2/2) + exp(Wv2)/(2 *
    exp(Wu2))) * pnorm(-object$S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu2/2)) *
    exp(-P2 * Wu2/2)/gamma(P2) * Hi2_2
  Probc1 <- 1 - exp(-exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta1[2:object$nuZUvar], nrow = 1),
    matrix(P1 * exp(Wu1), ncol = 1))
  margEff2 <- kronecker(matrix(delta2[2:object$nuZUvar], nrow = 1),
    matrix(P2 * exp(Wu2), ncol = 1))
  margEff_c <- ifelse(Group_c == 1, margEff1, margEff2)
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c1")
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c2")
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1],
    "_c")
  return(data.frame(margEff1, margEff2, margEff_c))
}
