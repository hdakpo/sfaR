################################################################################
#                                                                              #
# R internal functions for the sfaR package                                    #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Data: Cross sectional data & Pooled data                                     #
# Model: Sample Selection Stochastic Frontier Analysis                         #
# Convolution: halfnormal - normal                                             #
#------------------------------------------------------------------------------#

# Log-likelihood ----------
#' log-likelihood for halfnormal-normal distribution + selection bias
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
#' @param PREDICTIONS predictions from the first step probit model
#' @param uBound upper bound for inefficiency (default \code{Inf})
#' @param subdivisions number of divisions for GK quadrature
#' @param intol integration tolerance
#' @param gH Gauss-Hermite quadrature rule (list from fastGHQuad)
#' @param N number of observations
#' @param FiMat matrix of Halton draws
#' @noRd
## Gauss-Kronrod quadrature ----------
chalfnormlike_ss_GK <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  if (rho <= -1 || rho >= 1)
    return(NA)
  int_u <- function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      dnorm((epsilon + S * exp(Wu/2) * u)/exp(Wv/2)) *
        pnorm((rho * (epsilon + S * exp(Wu/2) * u)/exp(Wv/2) +
          PREDICTIONS)/sqrt(1 - rho^2)) * dnorm(u)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  }
  V <- Vectorize(int_u)
  ll <- log(2 * exp(-Wv/2) * V(epsilon, Wu, Wv, S, rho, PREDICTIONS))
  return(ll * wHvar)
}

## Hcubature ----------
chalfnormlike_ss_HCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  if (rho <= -1 || rho >= 1)
    return(NA)
  int_u <- function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      dnorm((epsilon + S * exp(Wu/2) * u)/exp(Wv/2)) *
        pnorm((rho * (epsilon + S * exp(Wu/2) * u)/exp(Wv/2) +
          PREDICTIONS)/sqrt(1 - rho^2)) * dnorm(u)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, epsilon = epsilon, Wu = Wu, Wv = Wv, S = S,
      rho = rho, vectorInterface = TRUE, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  }
  V <- Vectorize(int_u)
  ll <- log(2 * exp(-Wv/2) * V(epsilon, Wu, Wv, S, rho, PREDICTIONS))
  return(ll * wHvar)
}

## Pcubature ----------
chalfnormlike_ss_PCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  if (rho <= -1 || rho >= 1)
    return(NA)
  int_u <- function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      dnorm((epsilon + S * exp(Wu/2) * u)/exp(Wv/2)) *
        pnorm((rho * (epsilon + S * exp(Wu/2) * u)/exp(Wv/2) +
          PREDICTIONS)/sqrt(1 - rho^2)) * dnorm(u)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, epsilon = epsilon, Wu = Wu, Wv = Wv, S = S,
      rho = rho, vectorInterface = TRUE, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  }
  V <- Vectorize(int_u)
  ll <- log(2 * exp(-Wv/2) * V(epsilon, Wu, Wv, S, rho, PREDICTIONS))
  return(ll * wHvar)
}

## Gauss-Hermite quadrature ----------

### Obtained by a change of variables: k = u/sqrt(2)
chalfnormlike_ss_GH <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, gH, N) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  if (rho <= -1 || rho >= 1)
    return(NA)
  f_x_U <- function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    dnorm((epsilon + S * exp(Wu/2) * u * sqrt(2))/exp(Wv/2)) *
      pnorm((rho * (epsilon + S * exp(Wu/2) * u * sqrt(2))/exp(Wv/2) +
        PREDICTIONS)/sqrt(1 - rho^2))/exp(Wv/2) * 2/sqrt(pi)
  }
  ll <- numeric(N)
  for (i in seq_along(ll)) {
    ll[i] <- log(sum(gH$w[gH$x > 0] * f_x_U(u = gH$x[gH$x >
      0], epsilon = epsilon[i], Wu = Wu[i], Wv = Wv[i],
      S = S, rho = rho, PREDICTIONS = PREDICTIONS[i])))
  }
  return(ll * wHvar)
}

## Maximum Simulated Likelihood ----------
chalfnormlike_ss_MSL <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  if (rho <= -1 || rho >= 1)
    return(NA)
  Fi <- numeric(N)
  for (i in seq_along(1:N)) {
    Fi[i] <- mean(dnorm((epsilon[i] + S * exp(Wu[i]/2) *
      qnorm(1/2 + 1/2 * FiMat[i, ]))/exp(Wv[i]/2)) * pnorm((rho *
      (epsilon[i] + S * exp(Wu[i]/2) * qnorm(1/2 + 1/2 *
        FiMat[i, ]))/exp(Wv[i]/2) + PREDICTIONS[i])/sqrt(1 -
      rho^2)) * exp(-Wv[i]/2))
  }
  ll <- log(Fi)
  return(ll * wHvar)
}

# starting value for the log-likelihood ----------
#' starting values for halfnormal-normal distribution + selection bias
#' @param olsObj OLS object
#' @param epsiRes residuals from OLS
#' @param S integer for cost/prod estimation
#' @param nuZUvar number of Zu variables
#' @param nvZVvar number of Zv variables
#' @param uHvar matrix of Zu variables
#' @param vHvar matrix of Zv variables
#' @noRd
csthalfnorm_ss <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- sum(epsiRes^2)/length(epsiRes)
  m3 <- sum(epsiRes^3)/length(epsiRes)
  ## Coelli (1995) suggests 0.05 for gamma when wrong sign
  if (S * m3 > 0) {
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
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])[-length(olsObj)]
  } else {
    beta <- olsObj[-length(olsObj)]
  }
  rho <- olsObj[length(olsObj)]/sqrt(sum(epsiRes^2)/(length((epsiRes) -
    length(olsObj))))
  if (abs(rho) >= 1) {
    if (abs(olsObj[length(olsObj)] < 1)) {
      rho <- olsObj[length(olsObj)]
    } else {
      rho <- sign(rho) * 0.99
    }
  }
  names(rho) <- "rho"
  return(c(beta, delta, phi, rho))
}

# Gradient of the likelihood function ----------
#' log-likelihood for halfnormal-normal distribution + selection bias
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
#' @param PREDICTIONS predictions from the first step probit model
#' @param uBound upper bound for inefficiency (default \code{Inf})
#' @param subdivisions number of divisions for GK quadrature
#' @param intol integration tolerance
#' @param gH Gauss-Hermite quadrature rule (list from fastGHQuad)
#' @param N number of observations
#' @param FiMat matrix of Halton draws
#' @noRd
## Gauss-Kronrod quadrature ----------
cgradhalfnormlike_ss_GK <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = int_beta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS)/Fx, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = int_delta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS)/Fx, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx,
    FUN = "*"), int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## Hcubature ----------
cgradhalfnormlike_ss_HCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = int_beta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS)/Fx, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = int_delta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS)/Fx, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx,
    FUN = "*"), int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## Pcubature ----------
cgradhalfnormlike_ss_PCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = int_beta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS)/Fx, FUN = "*"), sweep(uHvar,
    MARGIN = 1, STATS = int_delta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS)/Fx, FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx,
    FUN = "*"), int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS)/Fx))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## Gauss-Hermite quadrature ----------
cgradhalfnormlike_ss_GH <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, gH, N) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  ewusqr <- exp(Wu/2)
  ewvsqr <- exp(Wv/2)
  sqrho <- sqrt(1 - rho^2)
  ut1 <- tcrossprod(gH$x[gH$x > 0], S * ewusqr)
  ut2 <- sweep(sqrt(2) * ut1, MARGIN = 2, STATS = epsilon,
    FUN = "+")
  ut3 <- dnorm(sweep(ut2, MARGIN = 2, STATS = ewvsqr, FUN = "/"))
  ut4 <- sweep(ut2, MARGIN = 2, STATS = rho/ewvsqr, FUN = "*")
  ut5 <- sweep(ut4, MARGIN = 2, STATS = PREDICTIONS, FUN = "+")
  ut6 <- pnorm(ut5/sqrho)
  ut7 <- sweep(ut2 * ut3 * ut6, MARGIN = 2, STATS = ewvsqr,
    FUN = "/")
  ut8 <- dnorm((ut5)/sqrho, 0, 1)
  ut9 <- (ut7 - rho * ut3 * ut8/sqrho)
  ut10 <- sweep(ut9, MARGIN = 2, STATS = (ewvsqr^2 * sqrt(pi)),
    FUN = "/")
  ut11 <- sweep(ut3 * ut6, MARGIN = 2, STATS = (ewvsqr * sqrt(pi)),
    FUN = "/")
  utSum <- apply(sweep(ut11, MARGIN = 1, STATS = gH$w[gH$x >
    0], FUN = "*"), 2, sum)
  ut12 <- (sqrt(2)/2 * (rho * ut3 * ut8/sqrho) - sqrt(2)/2 *
    (ut7))
  ut13 <- sweep(ut12, MARGIN = 2, STATS = ewusqr/(ewvsqr^2 *
    sqrt(pi)), FUN = "*")
  ut14 <- (0.5 * (ut7) - 0.5 * (rho * ut3 * ut8/sqrho))
  ut15 <- sweep(ut14 * ut2, MARGIN = 2, STATS = (ewvsqr^2 *
    sqrt(pi)), FUN = "/")
  ut16 <- sweep(ut3 * ut6, MARGIN = 2, STATS = ewvsqr * sqrt(pi)/(ewvsqr *
    sqrt(pi))^2, FUN = "*")
  ut17 <- sweep(ut2, MARGIN = 2, STATS = ewvsqr, FUN = "/")
  ut18 <- sweep((ut17 + rho * ut5/(1 - rho^2)) * ut3 * ut8,
    MARGIN = 2, STATS = (ewvsqr * sqrho * sqrt(pi)), FUN = "/")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(sweep(ut10, MARGIN = 1, STATS = gH$w[gH$x >
      0], FUN = "*"), MARGIN = 2, STATS = Xvar[, k], FUN = "*"),
      2, sum)/utSum
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sweep(ut13, MARGIN = 1, STATS = S *
      gH$w[gH$x > 0] * gH$x[gH$x > 0], FUN = "*"), MARGIN = 2,
      STATS = uHvar[, k], FUN = "*"), 2, sum)/utSum
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(sweep((ut15 - 0.5 * ut16), MARGIN = 1,
      STATS = gH$w[gH$x > 0], FUN = "*"), MARGIN = 2, STATS = vHvar[,
      k], FUN = "*"), 2, sum)/utSum
  }
  gradll <- cbind(gx, gu, gv, apply(sweep(ut18, MARGIN = 1,
    STATS = gH$w[gH$x > 0], FUN = "*"), 2, sum)/utSum)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

## Maximum Simulated Likelihood ----------
cgradhalfnormlike_ss_MSL <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, FiMat, N) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusqr <- exp(Wu/2)
  ewvsqr <- exp(Wv/2)
  ewvsqr2 <- exp(-(Wv/2))
  sqrho <- sqrt(1 - rho^2)
  sqrho2 <- (1 - rho^2)
  qFi <- qnorm(0.5 + 0.5 * FiMat)
  F1 <- sweep(qFi, MARGIN = 1, STATS = S * ewusqr, FUN = "*")
  F2 <- sweep(F1, MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- dnorm(sweep(F2, MARGIN = 1, STATS = ewvsqr, FUN = "/"))
  F4 <- sweep(sweep(F2, MARGIN = 1, STATS = rho/ewvsqr, FUN = "*"),
    MARGIN = 1, STATS = PREDICTIONS, FUN = "+")
  F5 <- pnorm(F4/sqrho)
  F6 <- dnorm(F4/sqrho)
  F7 <- sweep(F3 * F5 * F2, MARGIN = 1, STATS = ewvsqr, FUN = "/") -
    F6 * F3 * rho/sqrho
  F8 <- sweep(F7, MARGIN = 1, STATS = ewvsqr2/ewvsqr, FUN = "*")
  sumF <- apply(sweep(F3 * F5, MARGIN = 1, STATS = ewvsqr2,
    FUN = "*"), 1, sum)
  F9 <- F6 * F3 * rho/sqrho
  F10 <- sweep(F3 * F5 * F2, MARGIN = 1, STATS = ewvsqr, FUN = "/")
  F11 <- sweep((0.5 * F9 - 0.5 * F10) * qFi, MARGIN = 1, STATS = ewvsqr2 *
    ewusqr/ewvsqr, FUN = "*")
  F12 <- (sweep(0.5 * F3 * F2^2, MARGIN = 1, STATS = ewvsqr^2,
    FUN = "/") - 0.5 * F3) * F5 - sweep(0.5 * F6 * F3 * F2,
    MARGIN = 1, STATS = (rho/(ewvsqr * sqrho)), FUN = "*")
  F13 <- sweep(F12, MARGIN = 1, STATS = ewvsqr2, FUN = "*")
  F14 <- sweep(F2, MARGIN = 1, STATS = ewvsqr, FUN = "/") +
    F4 * rho/sqrho2
  F15 <- sweep(F14 * F6 * F3, MARGIN = 1, STATS = ewvsqr2/sqrho,
    FUN = "*")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_along(1:nXvar)) {
    gx[, k] <- apply(sweep(F8, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_along(1:nuZUvar)) {
    gu[, k] <- apply(sweep(S * F11, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_along(1:nvZVvar)) {
    gv[, k] <- apply(sweep(F13, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  gradll <- cbind(gx, gu, gv, apply(F15, 1, sum)/sumF)
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}

# Hessian of the likelihood function ----------
#' log-likelihood for halfnormal-normal distribution + selection bias
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
#' @param PREDICTIONS predictions from the first step probit model
#' @param uBound upper bound for inefficiency (default \code{Inf})
#' @param subdivisions number of divisions for GK quadrature
#' @param intol integration tolerance
#' @param gH Gauss-Hermite quadrature rule (list from fastGHQuad)
#' @param N number of observations
#' @param FiMat matrix of Halton draws
#' @noRd
## Gauss-Kronrod quadrature ----------
chesshalfnormlike_ss_GK <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ## f'(x) part
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  ## f''(x) part
  int_betabeta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((pnrho * sigx1/ewvsqr - rho * dnrho2/sqrho) *
        sigx1/ewvsqr - pnrho) * dnrho1 - rho * dnrho2 *
        (dnrho1 * sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3)/sqrho) *
        du * ewvsqr2/ewvsqr^2
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_betadelta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (rho * dnrho2 * (dnrho1 * sigx1/ewvsqr +
        rho * dnrho1 * sigx2/sigx3)/sqrho) - 0.5 * (((pnrho *
        sigx1/ewvsqr - rho * dnrho2/sqrho) * sigx1/ewvsqr -
        pnrho) * dnrho1)) * du * ewvsqr2 * ewusqr/ewvsqr^2)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_betaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * (sigx1^2/ewvsqr^2 - 2) - 0.5) * dnrho1 *
        pnrho * sigx1/ewvsqr - rho * (0.5 * ((dnrho1 *
        sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3) *
        sigx1/ewvsqr - dnrho1) + 0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * dnrho2/sqrho) * du * ewvsqr2/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_betarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
        sigx1/ewvsqr + dnrho1 * (rho * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2 - rho)/sigx3 - 1)) *
        dnrho2 * du * ewvsqr2/(ewvsqr * sqrho))
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_deltadelta <- Vectorize(function(epsilon, Wu, Wv, S,
    rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) - S *
        u * (0.5 * (((0.5 * (rho * dnrho2/sqrho) - 0.5 *
        (pnrho * sigx1/ewvsqr)) * sigx1/ewvsqr + 0.5 *
        pnrho) * dnrho1) + 0.5 * (rho * (0.5 * (dnrho1 *
        sigx1/ewvsqr) + 0.5 * (rho * dnrho1 * sigx2/sigx3)) *
        dnrho2/sqrho)) * ewusqr/ewvsqr) * du * ewvsqr2 *
        ewusqr/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_deltaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.25 + 0.5 * (1 - 0.5 * (sigx1^2/ewvsqr^2))) *
        dnrho1 * pnrho * sigx1/ewvsqr + rho * (0.5 *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) -
        0.5 * (0.5 * dnrho1 - (0.5 * (dnrho1 * sigx1/ewvsqr) +
          0.5 * (rho * dnrho1 * sigx2/sigx3)) * sigx1/ewvsqr)) *
        dnrho2/sqrho) * du * ewvsqr2 * ewusqr/ewvsqr)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_deltarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.5 + rho * (0.5 * rho - 0.5 * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2))/sigx3) * dnrho1 -
        0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
          sigx1/ewvsqr)) * dnrho2 * du * ewvsqr2 * ewusqr/(ewvsqr *
        sqrho))
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_phiphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((((0.5 * (0.5 * (sigx1^2/ewvsqr^2) - 1) - 0.25) *
        dnrho1 * pnrho * sigx1/ewvsqr - 0.5 * (rho *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) *
        dnrho2/sqrho))/ewvsqr - 0.5 * (rho * ((0.5 *
        (dnrho1 * sigx1/ewvsqr) + 0.5 * (rho * dnrho1 *
        sigx2/sigx3)) * sigx1/(ewvsqr^2 * sqrho) - 0.5 *
        (dnrho1 * ewvsqr * sqrho/(ewvsqr * sqrho)^2)) *
        dnrho2)) * sigx1 - 0.5 * ((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * pnrho - 0.5 * (rho * dnrho2 *
        dnrho1 * sigx1/(ewvsqr * sqrho)))) * du * ewvsqr2)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_phirho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) *
        dnrho1 * sigx1/ewvsqr) + dnrho1 * (rho * (0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) -
        0.5 * rho)/sigx3 - 0.5)) * sigx1/ewvsqr - 0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1)) *
        dnrho2 * du * ewvsqr2/sqrho)
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_rhorrho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * (rho -
        (sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) +
        PREDICTIONS + rho * (2 * (sigx1/ewvsqr) + 2 *
        (rho * sigx2/sigx3))) * dnrho2 * dnrho1 * du *
        ewvsqr2/(sigx3 * sqrho))
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    integrate(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lower = 0, upper = uBound, subdivisions = subdivisions,
      epsilon = epsilon, Wu = Wu, Wv = Wv, S = S, rho = rho,
      PREDICTIONS = PREDICTIONS, rel.tol = intol, stop.on.error = FALSE)$value
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (int_betabeta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS) * Fx - int_beta(epsilon, Wu, Wv, S,
      rho, PREDICTIONS)^2)/Fx^2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * (int_betadelta(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_beta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (int_betaphi(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2,
    FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- crossprod(Xvar,
    wHvar * (int_betarho(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (int_deltadelta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_delta(epsilon, Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (int_deltaphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_phi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- crossprod(uHvar, wHvar * (int_deltarho(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon, Wu, Wv,
    S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (int_phiphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_phi(epsilon,
      Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- crossprod(vHvar, wHvar *
    (int_phirho(epsilon, Wu, Wv, S, rho, PREDICTIONS) * Fx -
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon,
        Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (int_rhorrho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) * Fx - int_rho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS)^2)/Fx^2)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

## Hcubature ----------
chesshalfnormlike_ss_HCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ## f'(x) part
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  ## f''(x) part
  int_betabeta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((pnrho * sigx1/ewvsqr - rho * dnrho2/sqrho) *
        sigx1/ewvsqr - pnrho) * dnrho1 - rho * dnrho2 *
        (dnrho1 * sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3)/sqrho) *
        du * ewvsqr2/ewvsqr^2
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betadelta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (rho * dnrho2 * (dnrho1 * sigx1/ewvsqr +
        rho * dnrho1 * sigx2/sigx3)/sqrho) - 0.5 * (((pnrho *
        sigx1/ewvsqr - rho * dnrho2/sqrho) * sigx1/ewvsqr -
        pnrho) * dnrho1)) * du * ewvsqr2 * ewusqr/ewvsqr^2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * (sigx1^2/ewvsqr^2 - 2) - 0.5) * dnrho1 *
        pnrho * sigx1/ewvsqr - rho * (0.5 * ((dnrho1 *
        sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3) *
        sigx1/ewvsqr - dnrho1) + 0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * dnrho2/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
        sigx1/ewvsqr + dnrho1 * (rho * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2 - rho)/sigx3 - 1)) *
        dnrho2 * du * ewvsqr2/(ewvsqr * sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltadelta <- Vectorize(function(epsilon, Wu, Wv, S,
    rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) - S *
        u * (0.5 * (((0.5 * (rho * dnrho2/sqrho) - 0.5 *
        (pnrho * sigx1/ewvsqr)) * sigx1/ewvsqr + 0.5 *
        pnrho) * dnrho1) + 0.5 * (rho * (0.5 * (dnrho1 *
        sigx1/ewvsqr) + 0.5 * (rho * dnrho1 * sigx2/sigx3)) *
        dnrho2/sqrho)) * ewusqr/ewvsqr) * du * ewvsqr2 *
        ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.25 + 0.5 * (1 - 0.5 * (sigx1^2/ewvsqr^2))) *
        dnrho1 * pnrho * sigx1/ewvsqr + rho * (0.5 *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) -
        0.5 * (0.5 * dnrho1 - (0.5 * (dnrho1 * sigx1/ewvsqr) +
          0.5 * (rho * dnrho1 * sigx2/sigx3)) * sigx1/ewvsqr)) *
        dnrho2/sqrho) * du * ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.5 + rho * (0.5 * rho - 0.5 * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2))/sigx3) * dnrho1 -
        0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
          sigx1/ewvsqr)) * dnrho2 * du * ewvsqr2 * ewusqr/(ewvsqr *
        sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phiphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((((0.5 * (0.5 * (sigx1^2/ewvsqr^2) - 1) - 0.25) *
        dnrho1 * pnrho * sigx1/ewvsqr - 0.5 * (rho *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) *
        dnrho2/sqrho))/ewvsqr - 0.5 * (rho * ((0.5 *
        (dnrho1 * sigx1/ewvsqr) + 0.5 * (rho * dnrho1 *
        sigx2/sigx3)) * sigx1/(ewvsqr^2 * sqrho) - 0.5 *
        (dnrho1 * ewvsqr * sqrho/(ewvsqr * sqrho)^2)) *
        dnrho2)) * sigx1 - 0.5 * ((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * pnrho - 0.5 * (rho * dnrho2 *
        dnrho1 * sigx1/(ewvsqr * sqrho)))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phirho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) *
        dnrho1 * sigx1/ewvsqr) + dnrho1 * (rho * (0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) -
        0.5 * rho)/sigx3 - 0.5)) * sigx1/ewvsqr - 0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1)) *
        dnrho2 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rhorrho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * (rho -
        (sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) +
        PREDICTIONS + rho * (2 * (sigx1/ewvsqr) + 2 *
        (rho * sigx2/sigx3))) * dnrho2 * dnrho1 * du *
        ewvsqr2/(sigx3 * sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    hcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (int_betabeta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS) * Fx - int_beta(epsilon, Wu, Wv, S,
      rho, PREDICTIONS)^2)/Fx^2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * (int_betadelta(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_beta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (int_betaphi(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2,
    FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- crossprod(Xvar,
    wHvar * (int_betarho(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (int_deltadelta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_delta(epsilon, Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (int_deltaphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_phi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- crossprod(uHvar, wHvar * (int_deltarho(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon, Wu, Wv,
    S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (int_phiphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_phi(epsilon,
      Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- crossprod(vHvar, wHvar *
    (int_phirho(epsilon, Wu, Wv, S, rho, PREDICTIONS) * Fx -
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon,
        Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (int_rhorrho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) * Fx - int_rho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS)^2)/Fx^2)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

## Pcubature ----------
chesshalfnormlike_ss_PCUB <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, uBound,
  subdivisions, intol) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ## f'(x) part
  int_beta <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((dnrho1 * pnrho * sigx1/ewvsqr - rho * dnrho2 *
        dnrho1/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_delta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (S * u * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) * du *
        ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phi <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * (((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 *
        dnrho1) * pnrho - 0.5 * (rho * dnrho2 * dnrho1 *
        sigx1/(ewvsqr * sqrho))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rho <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      2 * ((sigx1/ewvsqr + rho * sigx2/(1 - rho^2)) * dnrho2 *
        dnrho1 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  ## f''(x) part
  int_betabeta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((pnrho * sigx1/ewvsqr - rho * dnrho2/sqrho) *
        sigx1/ewvsqr - pnrho) * dnrho1 - rho * dnrho2 *
        (dnrho1 * sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3)/sqrho) *
        du * ewvsqr2/ewvsqr^2
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betadelta <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (rho * dnrho2 * (dnrho1 * sigx1/ewvsqr +
        rho * dnrho1 * sigx2/sigx3)/sqrho) - 0.5 * (((pnrho *
        sigx1/ewvsqr - rho * dnrho2/sqrho) * sigx1/ewvsqr -
        pnrho) * dnrho1)) * du * ewvsqr2 * ewusqr/ewvsqr^2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * (sigx1^2/ewvsqr^2 - 2) - 0.5) * dnrho1 *
        pnrho * sigx1/ewvsqr - rho * (0.5 * ((dnrho1 *
        sigx1/ewvsqr + rho * dnrho1 * sigx2/sigx3) *
        sigx1/ewvsqr - dnrho1) + 0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * dnrho2/sqrho) * du * ewvsqr2/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_betarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
        sigx1/ewvsqr + dnrho1 * (rho * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2 - rho)/sigx3 - 1)) *
        dnrho2 * du * ewvsqr2/(ewvsqr * sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltadelta <- Vectorize(function(epsilon, Wu, Wv, S,
    rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * (0.5 * (0.5 * (rho * dnrho2 * dnrho1/sqrho) -
        0.5 * (dnrho1 * pnrho * sigx1/ewvsqr)) - S *
        u * (0.5 * (((0.5 * (rho * dnrho2/sqrho) - 0.5 *
        (pnrho * sigx1/ewvsqr)) * sigx1/ewvsqr + 0.5 *
        pnrho) * dnrho1) + 0.5 * (rho * (0.5 * (dnrho1 *
        sigx1/ewvsqr) + 0.5 * (rho * dnrho1 * sigx2/sigx3)) *
        dnrho2/sqrho)) * ewusqr/ewvsqr) * du * ewvsqr2 *
        ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltaphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.25 + 0.5 * (1 - 0.5 * (sigx1^2/ewvsqr^2))) *
        dnrho1 * pnrho * sigx1/ewvsqr + rho * (0.5 *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) -
        0.5 * (0.5 * dnrho1 - (0.5 * (dnrho1 * sigx1/ewvsqr) +
          0.5 * (rho * dnrho1 * sigx2/sigx3)) * sigx1/ewvsqr)) *
        dnrho2/sqrho) * du * ewvsqr2 * ewusqr/ewvsqr)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_deltarho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (S * u * ((0.5 + rho * (0.5 * rho - 0.5 * ((sigx1/ewvsqr +
        rho * sigx2/sigx3) * sigx2))/sigx3) * dnrho1 -
        0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1 *
          sigx1/ewvsqr)) * dnrho2 * du * ewvsqr2 * ewusqr/(ewvsqr *
        sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phiphi <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((((0.5 * (0.5 * (sigx1^2/ewvsqr^2) - 1) - 0.25) *
        dnrho1 * pnrho * sigx1/ewvsqr - 0.5 * (rho *
        (0.5 * (dnrho1 * sigx1^2/ewvsqr^2) - 0.5 * dnrho1) *
        dnrho2/sqrho))/ewvsqr - 0.5 * (rho * ((0.5 *
        (dnrho1 * sigx1/ewvsqr) + 0.5 * (rho * dnrho1 *
        sigx2/sigx3)) * sigx1/(ewvsqr^2 * sqrho) - 0.5 *
        (dnrho1 * ewvsqr * sqrho/(ewvsqr * sqrho)^2)) *
        dnrho2)) * sigx1 - 0.5 * ((0.5 * (dnrho1 * sigx1^2/ewvsqr^2) -
        0.5 * dnrho1) * pnrho - 0.5 * (rho * dnrho2 *
        dnrho1 * sigx1/(ewvsqr * sqrho)))) * du * ewvsqr2)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_phirho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((0.5 * ((sigx1/ewvsqr + rho * sigx2/sigx3) *
        dnrho1 * sigx1/ewvsqr) + dnrho1 * (rho * (0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) -
        0.5 * rho)/sigx3 - 0.5)) * sigx1/ewvsqr - 0.5 *
        ((sigx1/ewvsqr + rho * sigx2/sigx3) * dnrho1)) *
        dnrho2 * du * ewvsqr2/sqrho)
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_rhorrho <- Vectorize(function(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      dnrho2 <- dnorm(sigx2/sqrho, 0, 1)
      du <- dnorm(u)
      sigx3 <- (1 - rho^2)
      2 * (((sigx1/ewvsqr + rho * sigx2/sigx3) * (rho -
        (sigx1/ewvsqr + rho * sigx2/sigx3) * sigx2) +
        PREDICTIONS + rho * (2 * (sigx1/ewvsqr) + 2 *
        (rho * sigx2/sigx3))) * dnrho2 * dnrho1 * du *
        ewvsqr2/(sigx3 * sqrho))
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  int_u <- Vectorize(function(epsilon, Wu, Wv, S, rho, PREDICTIONS) {
    pcubature(f = function(u, epsilon, Wu, Wv, S, rho, PREDICTIONS) {
      ewusqr <- exp(Wu/2)
      ewvsqr <- exp(Wv/2)
      ewvsqr2 <- exp(-(Wv/2))
      sigx1 <- (S * u * ewusqr + epsilon)
      sigx2 <- (PREDICTIONS + rho * sigx1/ewvsqr)
      sqrho <- sqrt(1 - rho^2)
      pnrho <- pnorm(sigx2/sqrho)
      dnrho1 <- dnorm(sigx1/ewvsqr)
      du <- dnorm(u)
      ewvsqr2 * dnrho1 * pnorm((rho * (epsilon + S * ewusqr *
        u)/ewvsqr + PREDICTIONS)/sqrho) * 2 * du
    }, lowerLimit = 0, upperLimit = uBound, maxEval = subdivisions,
      fDim = 1, vectorInterface = TRUE, epsilon = epsilon,
      Wu = Wu, Wv = Wv, S = S, rho = rho, PREDICTIONS = PREDICTIONS,
      tol = intol)$integral
  })
  Fx <- int_u(epsilon, Wu, Wv, S, rho, PREDICTIONS)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = wHvar * (int_betabeta(epsilon, Wu, Wv, S, rho,
      PREDICTIONS) * Fx - int_beta(epsilon, Wu, Wv, S,
      rho, PREDICTIONS)^2)/Fx^2, FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = wHvar * (int_betadelta(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_beta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
    (int_betaphi(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2,
    FUN = "*"), vHvar)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- crossprod(Xvar,
    wHvar * (int_betarho(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_beta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      int_rho(epsilon, Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
    (int_deltadelta(epsilon, Wu, Wv, S, rho, PREDICTIONS) *
      Fx - int_delta(epsilon, Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2,
    FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = wHvar * (int_deltaphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
      Wu, Wv, S, rho, PREDICTIONS) * int_phi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS))/Fx^2, FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- crossprod(uHvar, wHvar * (int_deltarho(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * Fx - int_delta(epsilon,
    Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon, Wu, Wv,
    S, rho, PREDICTIONS))/Fx^2)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = wHvar * (int_phiphi(epsilon, Wu,
      Wv, S, rho, PREDICTIONS) * Fx - int_phi(epsilon,
      Wu, Wv, S, rho, PREDICTIONS)^2)/Fx^2, FUN = "*"),
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- crossprod(vHvar, wHvar *
    (int_phirho(epsilon, Wu, Wv, S, rho, PREDICTIONS) * Fx -
      int_phi(epsilon, Wu, Wv, S, rho, PREDICTIONS) * int_rho(epsilon,
        Wu, Wv, S, rho, PREDICTIONS))/Fx^2)
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (int_rhorrho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS) * Fx - int_rho(epsilon, Wu, Wv, S, rho,
    PREDICTIONS)^2)/Fx^2)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

## Gauss-Hermite quadrature ----------
chesshalfnormlike_ss_GH <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, gH, N) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- (Yvar - as.numeric(crossprod(matrix(beta), t(Xvar))))
  ewusqr <- exp(Wu/2)
  ewvsqr <- exp(Wv/2)
  sqrho <- sqrt(1 - rho^2)
  ut1 <- tcrossprod(gH$x[gH$x > 0], S * ewusqr)
  ut2 <- sweep(sqrt(2) * ut1, MARGIN = 2, STATS = epsilon,
    FUN = "+")
  ut3 <- dnorm(sweep(ut2, MARGIN = 2, STATS = ewvsqr, FUN = "/"))
  ut4 <- sweep(ut2, MARGIN = 2, STATS = rho/ewvsqr, FUN = "*")
  ut5 <- sweep(ut4, MARGIN = 2, STATS = PREDICTIONS, FUN = "+")
  ut6 <- pnorm(ut5/sqrho)
  ut7 <- sweep(ut2 * ut3 * ut6, MARGIN = 2, STATS = ewvsqr,
    FUN = "/")
  ut8 <- dnorm((ut5)/sqrho, 0, 1)
  ut9 <- (ut7 - rho * ut3 * ut8/sqrho)
  ut10 <- sweep(ut9, MARGIN = 2, STATS = (ewvsqr^2 * sqrt(pi)),
    FUN = "/")
  ut11 <- sweep(ut3 * ut6, MARGIN = 2, STATS = (ewvsqr * sqrt(pi)),
    FUN = "/")
  utSum <- apply(sweep(ut11, MARGIN = 1, STATS = gH$w[gH$x >
    0], FUN = "*"), 2, sum)
  ut12 <- (sqrt(2)/2 * (rho * ut3 * ut8/sqrho) - sqrt(2)/2 *
    (ut7))
  ut13 <- sweep(ut12, MARGIN = 2, STATS = ewusqr/(ewvsqr^2 *
    sqrt(pi)), FUN = "*")
  ut14 <- (0.5 * (ut7) - 0.5 * (rho * ut3 * ut8/sqrho))
  ut15 <- sweep(ut14 * ut2, MARGIN = 2, STATS = (ewvsqr^2 *
    sqrt(pi)), FUN = "/")
  ut16 <- sweep(ut3 * ut6, MARGIN = 2, STATS = ewvsqr * sqrt(pi)/(ewvsqr *
    sqrt(pi))^2, FUN = "*")
  ut17 <- sweep(ut2, MARGIN = 2, STATS = ewvsqr, FUN = "/")
  ut18 <- sweep((ut17 + rho * ut5/(1 - rho^2)) * ut3 * ut8,
    MARGIN = 2, STATS = (ewvsqr * sqrho * sqrt(pi)), FUN = "/")
  ut19 <- sweep(ut2^2, MARGIN = 2, STATS = ewvsqr^2, FUN = "/")
  ut20 <- sweep(rho * ut2 * ut8, MARGIN = 2, STATS = (ewvsqr *
    sqrho), FUN = "/")
  ut21 <- sweep(ut2 * ut3, MARGIN = 2, STATS = ewvsqr, FUN = "/")
  ut22 <- sweep((((ut19 - 1) * ut6 - ut20) * ut3 - rho * (ut21 +
    rho * ut3 * ut5/(1 - rho^2)) * ut8/sqrho), MARGIN = 2,
    STATS = (ewvsqr^3 * sqrt(pi)), FUN = "/")
  ut23 <- sweep((sqrt(2)/2 * (rho * (ut21 + rho * ut3 * ut5/(1 -
    rho^2)) * ut8/sqrho) - sqrt(2)/2 * (((ut19 - 1) * ut6 -
    ut20) * ut3)), MARGIN = 2, STATS = ewusqr/(ewvsqr^3 *
    sqrt(pi)), FUN = "*")
  ut24 <- ((0.5 * ((ut19 - 1) * ut6 - ut20) - 0.5 * ut6) *
    ut3 - 0.5 * (rho * (ut21 + rho * ut3 * ut5/(1 - rho^2)) *
    ut8/sqrho))
  ut25 <- sweep((ut24 * ut17 + 0.5 * (rho * ut3 * ut8/sqrho)),
    MARGIN = 2, STATS = (ewvsqr^2 * sqrt(pi)), FUN = "/")
  ut26 <- sweep(ut9 * sqrt(pi), MARGIN = 2, STATS = (ewvsqr *
    sqrt(pi))^2, FUN = "/")
  ut27 <- sweep(((ut21 + rho * ut3 * ut5/(1 - rho^2)) * (ut17 +
    rho * ut5/(1 - rho^2)) - (1 + rho^2/(1 - rho^2)) * ut3) *
    ut8, MARGIN = 2, STATS = (ewvsqr^2 * sqrho * sqrt(pi)),
    FUN = "/")
  ut28 <- (sqrt(2)/2 * (((sqrt(2)/2 - sqrt(2)/2 * (ut19)) *
    ut6 + sqrt(2)/2 * (ut20)) * ut3) + sqrt(2)/2 * (rho *
    (sqrt(2)/2 * (ut21) + sqrt(2)/2 * (rho * ut3 * ut5/(1 -
      rho^2))) * ut8/sqrho))
  ut29 <- sweep(ut28, MARGIN = 2, STATS = ewusqr/ewvsqr, FUN = "*")
  ut30 <- sweep(S * ut29, MARGIN = 1, STATS = gH$x[gH$x > 0],
    FUN = "*")
  ut31 <- sweep((0.5 * ut12 - ut30), MARGIN = 2, STATS = ewusqr/(ewvsqr^2 *
    sqrt(pi)), FUN = "*")
  ut32 <- ((0.5 * (((sqrt(2)/2 - sqrt(2)/2 * (ut19)) * ut6 +
    sqrt(2)/2 * (ut20)) * ut3) + 0.5 * (rho * (sqrt(2)/2 *
    (ut21) + sqrt(2)/2 * (rho * ut3 * ut5/(1 - rho^2))) *
    ut8/sqrho)) * ut17 + sqrt(2)/2 * ut14)
  ut33 <- sweep(ut32, MARGIN = 2, STATS = (ewvsqr^2 * sqrt(pi)),
    FUN = "/")
  ut34 <- sweep(ut12 * sqrt(pi), MARGIN = 2, STATS = (ewvsqr *
    sqrt(pi))^2, FUN = "/")
  ut35 <- sweep((ut33 - 0.5 * ut34), MARGIN = 2, STATS = ewusqr,
    FUN = "*")
  ut36 <- sweep(((sqrt(2)/2 * (rho^2/(1 - rho^2)) + sqrt(2)/2) *
    ut3 - (ut17 + rho * ut5/(1 - rho^2)) * (sqrt(2)/2 * (ut21) +
    sqrt(2)/2 * (rho * ut3 * ut5/(1 - rho^2)))) * ut8, MARGIN = 2,
    STATS = ewusqr/(ewvsqr^2 * sqrho * sqrt(pi)), FUN = "*")
  ut37 <- sweep((ut2 * ut6), MARGIN = 2, STATS = ewvsqr, FUN = "/")
  ut38 <- sweep((0.5 * (((0.5 * ut37 - 0.5 * (rho * ut8/sqrho)) *
    ut17 - 0.5 * ut6) * ut3) - 0.5 * (rho * (0.5 * (ut21) +
    0.5 * (rho * ut3 * ut5/(1 - rho^2))) * ut8/sqrho)) *
    ut2, MARGIN = 2, STATS = (ewvsqr^3 * sqrt(pi)), FUN = "/")
  ut39 <- sweep(ut14, MARGIN = 2, STATS = ewvsqr^2 * sqrt(pi)/(ewvsqr^2 *
    sqrt(pi))^2, FUN = "*")
  ut40 <- sweep((ut2^2 * ut3), MARGIN = 2, STATS = ewvsqr,
    FUN = "/")
  ut41 <- sweep((ut3), MARGIN = 2, STATS = ewvsqr, FUN = "*")
  ut42 <- sweep(ut6, MARGIN = 2, STATS = pi * ewvsqr^3/(ewvsqr *
    sqrt(pi))^2, FUN = "*")
  ut43 <- sweep((((0.5 * ut40 + 0.5 * ut41) * ut6 - (0.5 *
    (rho * ut2 * ut8/sqrho) + ut42) * ut3) * sqrt(pi)), MARGIN = 2,
    STATS = (ewvsqr * sqrt(pi))^2, FUN = "/")
  ut44 <- sweep(((ut17 + rho * ut5/(1 - rho^2)) * (0.5 * (ut21) +
    0.5 * (rho * ut3 * ut5/(1 - rho^2))) - (0.5 + 0.5 * (rho^2/(1 -
    rho^2))) * ut3) * ut2, MARGIN = 2, STATS = (ewvsqr^2 *
    sqrho * sqrt(pi)), FUN = "/")
  ut45 <- sweep(((ut17 + rho * ut5/(1 - rho^2)) * ut3), MARGIN = 2,
    STATS = ewvsqr * sqrho * sqrt(pi)/(ewvsqr * sqrho * sqrt(pi))^2,
    FUN = "*")
  ut46 <- sweep((rho * (2 * (ut17) + 2 * (rho * ut5/(1 - rho^2))) -
    (ut17 + rho * ut5/(1 - rho^2))^2 * ut5), MARGIN = 2,
    STATS = PREDICTIONS, FUN = "+")
  ut47 <- sweep(ut46, MARGIN = 2, STATS = ((1 - rho^2) * ewvsqr *
    sqrt(pi)), FUN = "/")
  ut48 <- sweep(rho * (ut17 + rho * ut5/(1 - rho^2)), MARGIN = 2,
    STATS = ewvsqr * sqrt(pi)/(ewvsqr * sqrho * sqrt(pi))^2,
    FUN = "*")
  Q <- length(gH$w[gH$x > 0])
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in 1:Q) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sweep(ut22, MARGIN = 1, STATS = gH$w[gH$x > 0], FUN = "*")[r,
        ]/utSum, FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sweep(ut23, MARGIN = 1, STATS = S * gH$w[gH$x > 0] *
        gH$x[gH$x > 0], FUN = "*")[r, ]/utSum, FUN = "*"),
      uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      sweep((ut25 - 0.5 * ut26), MARGIN = 1, STATS = gH$w[gH$x >
        0], FUN = "*")[r, ]/utSum, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      sweep(ut31, MARGIN = 1, STATS = S * gH$w[gH$x > 0] *
        gH$x[gH$x > 0], FUN = "*")[r, ]/utSum, FUN = "*"),
      uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      sweep(ut35, MARGIN = 1, STATS = S * gH$w[gH$x > 0] *
        gH$x[gH$x > 0], FUN = "*")[r, ]/utSum, FUN = "*"),
      vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      sweep((((ut38 - ut39) * ut2 - 0.5 * ut43)), MARGIN = 1,
        STATS = gH$w[gH$x > 0], FUN = "*")[r, ]/utSum,
      FUN = "*"), vHvar)
  }
  X1 <- matrix(nrow = N, ncol = nXvar)
  X2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    X1[, k] <- apply(sweep(sweep(ut10, MARGIN = 1, STATS = gH$w[gH$x >
      0], FUN = "*"), MARGIN = 2, STATS = Xvar[, k], FUN = "*"),
      2, sum)/utSum
    X2[, k] <- apply(sweep(sweep(ut27, MARGIN = 1, STATS = gH$w[gH$x >
      0], FUN = "*"), MARGIN = 2, STATS = Xvar[, k], FUN = "*"),
      2, sum)/utSum
  }
  ZU1 <- matrix(nrow = N, ncol = nuZUvar)
  ZU2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    ZU1[, k] <- apply(sweep(sweep(ut13, MARGIN = 1, STATS = S *
      gH$w[gH$x > 0] * gH$x[gH$x > 0], FUN = "*"), MARGIN = 2,
      STATS = uHvar[, k], FUN = "*"), 2, sum)/utSum
    ZU2[, k] <- apply(sweep(sweep(ut36, MARGIN = 1, STATS = S *
      gH$w[gH$x > 0] * gH$x[gH$x > 0], FUN = "*"), MARGIN = 2,
      STATS = uHvar[, k], FUN = "*"), 2, sum)/utSum
  }
  ZV1 <- matrix(nrow = N, ncol = nvZVvar)
  ZV2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    ZV1[, k] <- apply(sweep(sweep((ut15 - 0.5 * ut16), MARGIN = 1,
      STATS = gH$w[gH$x > 0], FUN = "*"), MARGIN = 2, STATS = vHvar[,
      k], FUN = "*"), 2, sum)/utSum
    ZV2[, k] <- apply(sweep(sweep((ut44 - 0.5 * ut45) * ut8,
      MARGIN = 1, STATS = gH$w[gH$x > 0], FUN = "*"), MARGIN = 2,
      STATS = vHvar[, k], FUN = "*"), 2, sum)/utSum
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) - crossprod(sweep(X1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), X1)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) - crossprod(sweep(X1, MARGIN = 1, STATS = wHvar,
    FUN = "*"), ZU1)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) - crossprod(sweep(X1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), ZV1)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(X2,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(X1, MARGIN = 1,
    STATS = wHvar * apply(sweep(ut18, MARGIN = 1, STATS = gH$w[gH$x >
      0], FUN = "*"), 2, sum)/utSum, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) - crossprod(sweep(ZU1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), ZU1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) -
    crossprod(sweep(ZU1, MARGIN = 1, STATS = wHvar, FUN = "*"),
      ZV1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(ZU2, MARGIN = 1, STATS = wHvar, FUN = "*") -
    sweep(ZU1, MARGIN = 1, STATS = wHvar * apply(sweep(ut18,
      MARGIN = 1, STATS = gH$w[gH$x > 0], FUN = "*"), 2,
      sum)/utSum, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) - crossprod(sweep(ZV1, MARGIN = 1, STATS = wHvar,
    FUN = "*"), ZV1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(ZV2,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(ZV1, MARGIN = 1,
    STATS = wHvar * apply(sweep(ut18, MARGIN = 1, STATS = gH$w[gH$x >
      0], FUN = "*"), 2, sum)/utSum, FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (colSums(sweep((ut47 + ut48) * ut3 *
    ut8/sqrho, MARGIN = 1, STATS = gH$w[gH$x > 0], FUN = "*"))/utSum -
    apply(sweep(ut18, MARGIN = 1, STATS = gH$w[gH$x > 0],
      FUN = "*"), 2, sum)^2/utSum^2))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

## Maximum Simulated Likelihood ----------
chesshalfnormlike_ss_MSL <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, wHvar, S, PREDICTIONS, FiMat, N) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  rho <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ewusqr <- exp(Wu/2)
  ewvsqr <- exp(Wv/2)
  ewvsqr2 <- exp(-(Wv/2))
  sqrho <- sqrt(1 - rho^2)
  sqrho2 <- (1 - rho^2)
  qFi <- qnorm(0.5 + 0.5 * FiMat)
  F1 <- sweep(qFi, MARGIN = 1, STATS = S * ewusqr, FUN = "*")
  F2 <- sweep(F1, MARGIN = 1, STATS = epsilon, FUN = "+")
  F3 <- dnorm(sweep(F2, MARGIN = 1, STATS = ewvsqr, FUN = "/"))
  F4 <- sweep(sweep(F2, MARGIN = 1, STATS = rho/ewvsqr, FUN = "*"),
    MARGIN = 1, STATS = PREDICTIONS, FUN = "+")
  F5 <- pnorm(F4/sqrho)
  F6 <- dnorm(F4/sqrho)
  F7 <- sweep(F3 * F5 * F2, MARGIN = 1, STATS = ewvsqr, FUN = "/") -
    F6 * F3 * rho/sqrho
  F8 <- sweep(F7, MARGIN = 1, STATS = ewvsqr2/ewvsqr, FUN = "*")
  sumF <- apply(sweep(F3 * F5, MARGIN = 1, STATS = ewvsqr2,
    FUN = "*"), 1, sum)
  F9 <- F6 * F3 * rho/sqrho
  F10 <- sweep(F3 * F5 * F2, MARGIN = 1, STATS = ewvsqr, FUN = "/")
  F11 <- sweep((0.5 * F9 - 0.5 * F10) * qFi, MARGIN = 1, STATS = ewvsqr2 *
    ewusqr/ewvsqr, FUN = "*")
  F12 <- (sweep(0.5 * F3 * F2^2, MARGIN = 1, STATS = ewvsqr^2,
    FUN = "/") - 0.5 * F3) * F5 - sweep(0.5 * F6 * F3 * F2,
    MARGIN = 1, STATS = (rho/(ewvsqr * sqrho)), FUN = "*")
  F13 <- sweep(F12, MARGIN = 1, STATS = ewvsqr2, FUN = "*")
  F14 <- sweep(F2, MARGIN = 1, STATS = ewvsqr, FUN = "/") +
    F4 * rho/sqrho2
  F15 <- sweep(F14 * F6 * F3, MARGIN = 1, STATS = ewvsqr2/sqrho,
    FUN = "*")
  F16 <- sweep(F5 * F2, MARGIN = 1, STATS = ewvsqr, FUN = "/")
  F17 <- F6 * rho/sqrho
  F18 <- sweep(F2, MARGIN = 1, STATS = ewvsqr, FUN = "/")
  F19 <- (((F16 - F17) * F18 - F5) * F3 - rho * F6 * (F3 *
    F18 + rho * F3 * F4/sqrho2)/sqrho)
  F20 <- sweep(F19, MARGIN = 1, STATS = ewvsqr2/ewvsqr^2, FUN = "*")
  F21 <- (0.5 * (rho * F6 * (F3 * F18 + rho * F3 * F4/sqrho2)/sqrho) -
    0.5 * (((F16 - F17) * F18 - F5) * F3))
  F22 <- sweep(F21 * qFi, MARGIN = 1, STATS = ewvsqr2 * ewusqr/ewvsqr^2,
    FUN = "*")
  F23 <- sweep(F2^2, MARGIN = 1, STATS = ewvsqr^2, FUN = "/") -
    2
  F24 <- sweep(F3 * F2^2, MARGIN = 1, STATS = ewvsqr^2, FUN = "/")
  F25 <- ((0.5 * F23 - 0.5) * F3 * F16 - rho * (0.5 * ((F3 *
    F18 + rho * F3 * F4/sqrho2) * F18 - F3) + 0.5 * F24 -
    0.5 * F3) * F6/sqrho)
  F26 <- sweep(F25, MARGIN = 1, STATS = ewvsqr2/ewvsqr, FUN = "*")
  F27 <- (F14 * F3 * F18 + F3 * (rho * (F14 * F4 - rho)/sqrho2 -
    1))
  F28 <- sweep(F27 * F6, MARGIN = 1, STATS = ewvsqr2/(ewvsqr *
    sqrho), FUN = "*")
  F29 <- sweep(qFi, MARGIN = 1, STATS = ewusqr/ewvsqr, FUN = "*")
  F30 <- (0.5 * (0.5 * F9 - 0.5 * F10) - S * (0.5 * (((0.5 *
    (F17) - 0.5 * (F16)) * F18 + 0.5 * F5) * F3) + 0.5 *
    (rho * (0.5 * (F3 * F18) + 0.5 * (rho * F3 * F4/sqrho2)) *
      F6/sqrho)) * F29)
  F31 <- sweep(F30 * F29, MARGIN = 1, STATS = ewvsqr2, FUN = "*")
  F32 <- sweep(F2^2, MARGIN = 1, STATS = ewvsqr^2, FUN = "/")
  F33 <- ((0.25 + 0.5 * (1 - 0.5 * F32)) * F3 * F16 + rho *
    (0.5 * (0.5 * F24 - 0.5 * F3) - 0.5 * (0.5 * F3 - (0.5 *
      (F3 * F18) + 0.5 * (rho * F3 * F4/sqrho2)) * F18)) *
    F6/sqrho)
  F34 <- sweep(F33 * F29, MARGIN = 1, STATS = ewvsqr2, FUN = "*")
  F35 <- ((0.5 + rho * (0.5 * rho - 0.5 * (F14 * F4))/sqrho2) *
    F3 - 0.5 * (F14 * F3 * F18))
  F36 <- sweep(F35 * F6 * qFi, MARGIN = 1, STATS = ewvsqr2 *
    ewusqr/(ewvsqr * sqrho), FUN = "*")
  F37 <- ((0.5 * (0.5 * F32 - 1) - 0.25) * F3 * F16 - 0.5 *
    (rho * (0.5 * F24 - 0.5 * F3) * F6/sqrho))
  F38 <- sweep(F2, MARGIN = 1, STATS = ewvsqr^2 * sqrho, FUN = "/")
  F39 <- sweep(F3, MARGIN = 1, STATS = (ewvsqr * sqrho/(ewvsqr *
    sqrho)^2), FUN = "*")
  F40 <- sweep(F37, MARGIN = 1, STATS = ewvsqr, FUN = "/")
  F41 <- ((F40 - 0.5 * (rho * ((0.5 * (F3 * F18) + 0.5 * (rho *
    F3 * F4/sqrho2)) * F38 - 0.5 * F39) * F6)) * F2 - 0.5 *
    F12)
  F42 <- sweep(F41, MARGIN = 1, STATS = ewvsqr2, FUN = "*")
  F43 <- ((0.5 * (F14 * F3 * F18) + F3 * (rho * (0.5 * (F14 *
    F4) - 0.5 * rho)/sqrho2 - 0.5)) * F18 - 0.5 * (F14 *
    F3))
  F44 <- sweep(F43 * F6, MARGIN = 1, STATS = ewvsqr2/sqrho,
    FUN = "*")
  F45 <- rho * (2 * (F18) + 2 * (rho * F4/sqrho2))
  F46 <- F14 * (rho - F14 * F4)
  F47 <- sweep(F46 + F45, MARGIN = 1, STATS = PREDICTIONS,
    FUN = "+")
  F48 <- sweep((F47) * F6 * F3, MARGIN = 1, STATS = ewvsqr2/(sqrho2 *
    sqrho), FUN = "*")
  Q <- dim(FiMat)[2]
  HX1 <- list()
  HXU1 <- list()
  HXV1 <- list()
  HU1 <- list()
  HUV1 <- list()
  HV1 <- list()
  for (r in seq_along(1:Q)) {
    HX1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      F20[, r]/sumF, FUN = "*"), Xvar)
    HXU1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      S * F22[, r]/sumF, FUN = "*"), uHvar)
    HXV1[[r]] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = wHvar *
      F26[, r]/sumF, FUN = "*"), vHvar)
    HU1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      S * F31[, r]/sumF, FUN = "*"), uHvar)
    HUV1[[r]] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = wHvar *
      S * F34[, r]/sumF, FUN = "*"), vHvar)
    HV1[[r]] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = wHvar *
      F42[, r]/sumF, FUN = "*"), vHvar)
  }
  X1 <- matrix(nrow = N, ncol = nXvar)
  X2 <- matrix(nrow = N, ncol = nXvar)
  for (k in seq_along(1:nXvar)) {
    X1[, k] <- apply(sweep(F8, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/sumF
    X2[, k] <- apply(sweep(F28, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  ZU1 <- matrix(nrow = N, ncol = nuZUvar)
  ZU2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in seq_along(1:nuZUvar)) {
    ZU1[, k] <- apply(sweep(S * F11, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF
    ZU2[, k] <- apply(sweep(S * F36, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  ZV1 <- matrix(nrow = N, ncol = nvZVvar)
  ZV2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in seq_along(1:nvZVvar)) {
    ZV1[, k] <- apply(sweep(F13, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF
    ZV2[, k] <- apply(sweep(F44, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF
  }
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- Reduce("+", HX1) - crossprod(sweep(X1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), X1)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- Reduce("+",
    HXU1) - crossprod(sweep(X1, MARGIN = 1, STATS = wHvar,
    FUN = "*"), ZU1)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar +
    nvZVvar)] <- Reduce("+", HXV1) - crossprod(sweep(X1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), ZV1)
  hessll[1:nXvar, nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(X2,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(X1, MARGIN = 1,
    STATS = wHvar * apply(F15, 1, sum)/sumF, FUN = "*"))
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- Reduce("+", HU1) - crossprod(sweep(ZU1,
    MARGIN = 1, STATS = wHvar, FUN = "*"), ZU1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar +
    1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+", HUV1) -
    crossprod(sweep(ZU1, MARGIN = 1, STATS = wHvar, FUN = "*"),
      ZV1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), nXvar + nuZUvar + nvZVvar +
    1] <- colSums(sweep(ZU2, MARGIN = 1, STATS = wHvar, FUN = "*") -
    sweep(ZU1, MARGIN = 1, STATS = wHvar * apply(F15, 1,
      sum)/sumF, FUN = "*"))
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- Reduce("+",
    HV1) - crossprod(sweep(ZV1, MARGIN = 1, STATS = wHvar,
    FUN = "*"), ZV1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar),
    nXvar + nuZUvar + nvZVvar + 1] <- colSums(sweep(ZV2,
    MARGIN = 1, STATS = wHvar, FUN = "*") - sweep(ZV1, MARGIN = 1,
    STATS = wHvar * apply(F15, 1, sum)/sumF, FUN = "*"))
  hessll[nXvar + nuZUvar + nvZVvar + 1, nXvar + nuZUvar + nvZVvar +
    1] <- sum(wHvar * (apply(F48, 1, sum) - apply(F15, 1,
    sum)^2/sumF)/sumF)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------
#' optimizations solve for halfnormal-normal distribution + selection bias
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
# @param selectDum selection variable for first step probit
# model
#' @param wHvar vector of weights (weighted likelihood)
#' @param S integer for cost/prod estimation
#' @param method algorithm for solver
#' @param printInfo logical print info during optimization
#' @param itermax maximum iteration
#' @param stepmax stepmax for ucminf
#' @param tol parameter tolerance
#' @param gradtol gradient tolerance
#' @param hessianType how hessian is computed
#' @param qac qac option for maxLik
#' @param PREDICTIONS predictions from the first step probit model
#' @param uBound upper bound for inefficiency (default \code{Inf})
#' @param subdivisions number of divisions for GK quadrature
#' @param intol integration tolerance
#' @param gH Gauss-Hermite quadrature rule (list from fastGHQuad)
#' @param N number of observations
#' @param FiMat matrix of Halton draws
#' @noRd
## Gauss-Kronrod quadrature ----------
halfnormAlgOpt_ss_GK <- function(start, olsParam, dataTable,
  S, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, selectDum,
  wHvar, PREDICTIONS, uBound, subdivisions, intol, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm_ss(olsObj = olsParam, epsiRes = dataTable[["ols2stepResiduals"]][selectDum ==
    1], S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike_ss_GK(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, uBound = uBound, subdivisions = subdivisions,
    intol = intol))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_GK(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
      subdivisions = subdivisions, intol = intol)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike_ss_GK,
    grad = cgradhalfnormlike_ss_GK, hess = chesshalfnormlike_ss_GK,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_GK(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
      subdivisions = subdivisions, intol = intol)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_GK(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
      subdivisions = subdivisions, intol = intol)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hs = function(parm) as(-chesshalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
      "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_GK(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
      subdivisions = subdivisions, intol = intol)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hess = function(parm) -chesshalfnormlike_ss_GK(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(chalfnormlike_ss_GK(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), gradient = function(parm) -colSums(cgradhalfnormlike_ss_GK(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), hessian = function(parm) -chesshalfnormlike_ss_GK(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike_ss_GK(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol))
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
      mleObj$hessian <- chesshalfnormlike_ss_GK(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike_ss_GK(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
  }
  mleObj$logL_OBS <- chalfnormlike_ss_GK(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, uBound = uBound, subdivisions = subdivisions,
    intol = intol)
  mleObj$gradL_OBS <- cgradhalfnormlike_ss_GK(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

## Hcubature ----------
halfnormAlgOpt_ss_HCUB <- function(start, olsParam, dataTable,
  S, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, selectDum,
  wHvar, PREDICTIONS, uBound, subdivisions, intol, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm_ss(olsObj = olsParam, epsiRes = dataTable[["ols2stepResiduals"]][selectDum ==
    1], S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike_ss_HCUB(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, uBound = uBound, subdivisions = subdivisions,
    intol = intol))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf::ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike_ss_HCUB,
    grad = cgradhalfnormlike_ss_HCUB, hess = chesshalfnormlike_ss_HCUB,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hs = function(parm) as(-chesshalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
      "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hess = function(parm) -chesshalfnormlike_ss_HCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(chalfnormlike_ss_HCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), gradient = function(parm) -colSums(cgradhalfnormlike_ss_HCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), hessian = function(parm) -chesshalfnormlike_ss_HCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike_ss_HCUB(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol))
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
      mleObj$hessian <- chesshalfnormlike_ss_HCUB(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike_ss_HCUB(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
  }
  mleObj$logL_OBS <- chalfnormlike_ss_HCUB(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)
  mleObj$gradL_OBS <- cgradhalfnormlike_ss_HCUB(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

## Pcubature ----------
halfnormAlgOpt_ss_PCUB <- function(start, olsParam, dataTable,
  S, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, selectDum,
  wHvar, PREDICTIONS, uBound, subdivisions, intol, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm_ss(olsObj = olsParam, epsiRes = dataTable[["ols2stepResiduals"]][selectDum ==
    1], S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike_ss_PCUB(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, uBound = uBound, subdivisions = subdivisions,
    intol = intol))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike_ss_PCUB,
    grad = cgradhalfnormlike_ss_PCUB, hess = chesshalfnormlike_ss_PCUB,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), sr1 = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
      stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
      report.precision = 1L)), sparse = trustOptim::trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hs = function(parm) as(-chesshalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
      "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = marqLevAlg::mla(b = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol)),
    hess = function(parm) -chesshalfnormlike_ss_PCUB(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol),
    print.info = printInfo, maxiter = itermax, epsa = gradtol,
    epsb = gradtol), nlminb = nlminb(start = startVal, objective = function(parm) -sum(chalfnormlike_ss_PCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), gradient = function(parm) -colSums(cgradhalfnormlike_ss_PCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)), hessian = function(parm) -chesshalfnormlike_ss_PCUB(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol), control = list(iter.max = itermax,
    trace = if (printInfo) 1 else 0, eval.max = itermax,
    rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike_ss_PCUB(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      uBound = uBound, subdivisions = subdivisions, intol = intol))
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
      mleObj$hessian <- chesshalfnormlike_ss_PCUB(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike_ss_PCUB(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        uBound = uBound, subdivisions = subdivisions,
        intol = intol)
  }
  mleObj$logL_OBS <- chalfnormlike_ss_PCUB(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)
  mleObj$gradL_OBS <- cgradhalfnormlike_ss_PCUB(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, uBound = uBound,
    subdivisions = subdivisions, intol = intol)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

## Gauss-Hermite quadrature ----------
halfnormAlgOpt_ss_GH <- function(start, olsParam, dataTable,
  S, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, selectDum,
  wHvar, PREDICTIONS, gH, N, method, printInfo, itermax, stepmax,
  tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm_ss(olsObj = olsParam, epsiRes = dataTable[["ols2stepResiduals"]][selectDum ==
    1], S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike_ss_GH(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, gH = gH, N = N))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_GH(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, gH = gH, N = N)),
    gr = function(parm) -colSums(cgradhalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike_ss_GH,
    grad = cgradhalfnormlike_ss_GH, hess = chesshalfnormlike_ss_GH,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, gH = gH,
    N = N), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(chalfnormlike_ss_GH(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, gH = gH,
    N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GH(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, gH = gH,
    N = N)), method = "SR1", control = list(maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(chalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N)), hs = function(parm) as(-chesshalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(chalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N)), hess = function(parm) -chesshalfnormlike_ss_GH(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(chalfnormlike_ss_GH(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        gH = gH, N = N)), gradient = function(parm) -colSums(cgradhalfnormlike_ss_GH(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        gH = gH, N = N)), hessian = function(parm) -chesshalfnormlike_ss_GH(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        gH = gH, N = N), control = list(iter.max = itermax,
        trace = if (printInfo) 1 else 0, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike_ss_GH(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      gH = gH, N = N))
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
      mleObj$hessian <- chesshalfnormlike_ss_GH(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        gH = gH, N = N)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike_ss_GH(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        gH = gH, N = N)
  }
  mleObj$logL_OBS <- chalfnormlike_ss_GH(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, gH = gH, N = N)
  mleObj$gradL_OBS <- cgradhalfnormlike_ss_GH(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, gH = gH,
    N = N)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

## Maximum Simulated Likelihood ----------
halfnormAlgOpt_ss_MSL <- function(start, olsParam, dataTable,
  S, nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, selectDum,
  wHvar, PREDICTIONS, FiMat, N, method, printInfo, itermax,
  stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm_ss(olsObj = olsParam, epsiRes = dataTable[["ols2stepResiduals"]][selectDum ==
    1], S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike_ss_MSL(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat, N = N))
  if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
    maxRoutine <- switch(method, bfgs = function(...) maxLik::maxBFGS(...),
      bhhh = function(...) maxLik::maxBHHH(...), nr = function(...) maxLik::maxNR(...),
      nm = function(...) maxLik::maxNM(...), cg = function(...) maxLik::maxCG(...),
      sann = function(...) maxLik::maxSANN(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike_ss_MSL(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
      S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat,
      N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike_ss_MSL,
    grad = cgradhalfnormlike_ss_MSL, hess = chesshalfnormlike_ss_MSL,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat,
    N = N), sr1 = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(chalfnormlike_ss_MSL(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat,
    N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_MSL(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat,
    N = N)), method = "SR1", control = list(maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trustOptim::trust.optim(x = startVal, fn = function(parm) -sum(chalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N)), hs = function(parm) as(-chesshalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = marqLevAlg::mla(b = startVal, fn = function(parm) -sum(chalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N)), gr = function(parm) -colSums(cgradhalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N)), hess = function(parm) -chesshalfnormlike_ss_MSL(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(chalfnormlike_ss_MSL(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        FiMat = FiMat, N = N)), gradient = function(parm) -colSums(cgradhalfnormlike_ss_MSL(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        FiMat = FiMat, N = N)), hessian = function(parm) -chesshalfnormlike_ss_MSL(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        FiMat = FiMat, N = N), control = list(iter.max = itermax,
        trace = if (printInfo) 1 else 0, eval.max = itermax,
        rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike_ss_MSL(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
      FiMat = FiMat, N = N))
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
      mleObj$hessian <- chesshalfnormlike_ss_MSL(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        FiMat = FiMat, N = N)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike_ss_MSL(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS,
        FiMat = FiMat, N = N)
  }
  mleObj$logL_OBS <- chalfnormlike_ss_MSL(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = wHvar,
    S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat, N = N)
  mleObj$gradL_OBS <- cgradhalfnormlike_ss_MSL(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
    wHvar = wHvar, S = S, PREDICTIONS = PREDICTIONS, FiMat = FiMat,
    N = N)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------
#' efficiencies for halfnormal-normal distribution + selection bias
#' @param object object of class selectioncross
#' @param level level for confidence interval
#' @noRd
chalfnormeff_ss <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  rho <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  Xvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ])) - as.numeric(crossprod(matrix(beta), t(Xvar)))
  u <- numeric(object$Nobs)
  if (object$logDepVar == TRUE) {
    teBC <- numeric(object$Nobs)
    teBC_reciprocal <- numeric(object$Nobs)
  }
  for (i in 1:object$Nobs) {
    sigmasq <- exp(Wv[i]) + exp(Wu[i])
    sigmaC <- exp(Wu[i] + Wv[i])/sigmasq
    pi_dot <- -object$S * epsilon[i] * exp(Wu[i])/sigmasq
    kappa_dot <- c(-object$dataTable$PROBIT_PREDICTIONS[i] -
      rho * exp(Wv[i]/2) * epsilon[i]/sigmasq, object$S *
      exp(Wu[i]) * epsilon[i]/sigmasq)
    D_c <- c(object$S * rho/exp(Wv[i]/2), 1)
    D_dot <- D_c * sigmaC
    Delta_dot <- matrix(c(1 - rho^2, rep(0, 3)), nrow = 2,
      ncol = 2) + (matrix(D_c, ncol = 1) %*% matrix(D_c,
      ncol = 2)) * sigmaC
    u[i] <- pi_dot + (D_dot %*% t(mnorm::pmnorm(lower = rep(-Inf,
      2), upper = rep(0, 2), mean = kappa_dot, sigma = Delta_dot,
      grad_upper = TRUE)$grad))/mnorm::pmnorm(lower = rep(-Inf,
      2), upper = rep(0, 2), mean = kappa_dot, sigma = Delta_dot)$prob
    if (object$logDepVar == TRUE) {
      teBC[i] <- mnorm::pmnorm(lower = rep(-Inf, 2), upper = -D_dot,
        mean = kappa_dot, sigma = Delta_dot)$prob/mnorm::pmnorm(lower = rep(-Inf,
        2), upper = rep(0, 2), mean = kappa_dot, sigma = Delta_dot)$prob *
        exp(-pi_dot + 1/2 * sigmaC)
      teBC_reciprocal[i] <- mnorm::pmnorm(lower = rep(-Inf,
        2), upper = D_dot, mean = kappa_dot, sigma = Delta_dot)$prob/mnorm::pmnorm(lower = rep(-Inf,
        2), upper = rep(0, 2), mean = kappa_dot, sigma = Delta_dot)$prob *
        exp(pi_dot + 1/2 * sigmaC)
    }
  }
  if (object$logDepVar == TRUE) {
    res <- data.frame(u = rep(NA, object$Ninit), teJLMS = rep(NA,
      object$Ninit), teBC = rep(NA, object$Ninit), teBC_reciprocal = rep(NA,
      object$Ninit))
    res_eff <- data.frame(u = u, teJLMS = exp(-u), teBC = teBC,
      teBC_reciprocal = teBC_reciprocal)
    res[object$dataTable[all.vars(object$selectionF)[1]] ==
      1, ] <- res_eff
  } else {
    res <- data.frame(u = rep(NA, object$Ninit))
    res_eff <- data.frame(u = u)
    res[object$dataTable[all.vars(object$selectionF)[1]] ==
      1, ] <- res_eff
  }
  return(res)
}

# Marginal effects on inefficiencies ----------
#' marginal impact on efficiencies for halfnormal-normal distribution + selection bias
#' @param object object of class selectioncross
#' @noRd
cmarghalfnorm_Eu_ss <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  res <- matrix(nrow = object$Ninit, ncol = object$nuZUvar)
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  res[object$dataTable[all.vars(object$selectionF)[1]] == 1,
    ] <- margEff
  colnames(res) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(res))
}

cmarghalfnorm_Vu_ss <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable[object$dataTable[all.vars(object$selectionF)[1]] ==
    1, ], rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  res <- matrix(nrow = object$Ninit, ncol = object$nuZUvar)
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  res[object$dataTable[all.vars(object$selectionF)[1]] == 1,
    ] <- margEff
  colnames(res) <- paste0("Eu_", colnames(uHvar)[-1])
  return(data.frame(res))
}
