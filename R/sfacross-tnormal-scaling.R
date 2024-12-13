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
  musca <- exp(crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]) * tau
  Wusca <- cu + 2 * crossprod(matrix(delta), t(uHvar[, -1, drop = FALSE]))[1, ]
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
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- usca
  .e3 <- 2 * .e1 + cu
  .e4 <- exp(.e3)
  .e5 <- exp(Wvsca)
  .e6 <- .e4 + .e5
  .e7 <- exp(.e1)
  .e8 <- tau * .e7
  .e10 <- epsilon
  .e13 <- sqrt(.e4 * .e5/.e6)
  .e14 <- S * .e10
  .e15 <- .e6 * .e13
  .e16 <- .e14 + .e8
  .e18 <- S * .e4 * .e10
  .e19 <- .e8 * .e5
  .e20 <- .e19 - .e18
  .e21 <- .e16/sqrt(.e6)
  .e22 <- .e20/.e15
  .e24 <- exp(.e3/2)
  .e25 <- dnorm(.e21)
  .e26 <- dnorm(.e21)
  .e27 <- dnorm(.e22)
  .e28 <- pnorm(.e22)
  .e29 <- .e8/.e24
  .e30 <- .e15^2
  .e31 <- .e26 * .e16
  .e32 <- (0.5 * (.e26 * .e16^2/(.e25 * .e6)) - 0.5)/.e6
  .e33 <- .e31/.e25
  .e34 <- dnorm(.e29)
  .e35 <- .e24 * pnorm(.e29)
  .e36 <- .e4/.e6
  .e37 <- .e28 * .e13
  gradll <- cbind(S * Xvar * ((.e33 + .e27 * .e4/.e37)/.e6), uHvar[, -1, drop = FALSE] *
    (((.e19 - 2 * .e18)/.e15 - (0.5 * ((2 - 2 * .e36) * .e5/.e13) + 2 * .e13) *
      .e4 * .e20/.e30) * .e27/.e28 - (.e31 * (.e8 - .e4 * .e16/.e6)/.e25 +
      .e4)/.e6), ((.e27 * .e5/.e37 - .e33)/.e6 - .e34/.e35) * .e7, (.e32 -
    ((0.5 * ((1 - .e36) * .e5/.e13) + .e13) * .e20/.e30 + .e14/.e15) * .e27/.e28) *
    .e4 + 0.5 * (tau * .e34 * .e7/.e35), vHvar * ((.e32 + .e27 * (.e8/.e15 -
    (0.5 * ((1 - .e5/.e6) * .e4/.e13) + .e13) * .e20/.e30)/.e28) * .e5))
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
  Wvsca <- crossprod(matrix(phi), t(vHvar))[1, ]
  epsilon <- Yvar - crossprod(matrix(beta), t(Xvar))[1, ]
  .e1 <- usca
  .e3 <- 2 * .e1 + cu
  .e4 <- exp(.e3)
  .e5 <- exp(Wvsca)
  .e6 <- .e4 + .e5
  .e7 <- exp(.e1)
  .e8 <- tau * .e7
  .e10 <- epsilon
  .e11 <- .e4 * .e5
  .e13 <- sqrt(.e11/.e6)
  .e14 <- .e6 * .e13
  .e15 <- S * .e10
  .e16 <- .e8 * .e5
  .e18 <- S * .e4 * .e10
  .e19 <- .e15 + .e8
  .e20 <- .e16 - .e18
  .e21 <- .e19/sqrt(.e6)
  .e22 <- .e20/.e14
  .e23 <- .e4/.e6
  .e24 <- dnorm(.e21)
  .e25 <- .e14^2
  .e26 <- pnorm(.e22)
  .e28 <- dnorm(.e22)
  .e29 <- .e19^2
  .e31 <- exp(.e3/2)
  .e32 <- 2 * .e23
  .e33 <- .e5/.e6
  .e34 <- 2 - .e32
  .e35 <- 1 - .e33
  .e36 <- 1 - .e23
  .e37 <- .e35 * .e4
  .e38 <- .e36 * .e5
  .e39 <- .e24 * .e6
  .e40 <- .e34 * .e5
  .e41 <- 2 * .e13
  .e42 <- .e8/.e31
  .e45 <- 0.5 * (.e40/.e13) + .e41
  .e46 <- .e26 * .e13
  .e50 <- 0.5 * (.e38/.e13) + .e13
  .e52 <- 0.5 * (.e37/.e13) + .e13
  .e53 <- .e45 * .e4
  .e54 <- .e24 * .e29
  .e58 <- pnorm(.e42)
  .e59 <- .e8 - .e4 * .e19/.e6
  .e60 <- .e16 - 2 * .e18
  .e61 <- .e15/.e14
  .e62 <- .e8/.e14
  .e65 <- .e50 * .e20/.e25 + .e61
  .e71 <- .e60/.e14 - .e53 * .e20/.e25
  .e72 <- .e62 - .e52 * .e20/.e25
  .e73 <- .e31 * .e58
  .e74 <- .e24/.e24
  .e75 <- .e39^2
  .e76 <- .e46^2
  .e77 <- dnorm(.e42)
  .e78 <- 0.5 * (.e54/.e39)
  .e79 <- .e24 * .e19
  .e81 <- .e78 - 0.5
  .e82 <- .e28 * .e4
  .e83 <- .e73^2
  .e84 <- .e29/.e6
  .e85 <- .e22 + .e28/.e26
  .e86 <- 0.5 - 0.5 * .e74
  .e87 <- .e28 * .e5
  .e88 <- .e81/.e6
  .e89 <- .e79 * .e59
  .e90 <- .e87/.e46
  .e93 <- (.e86 * .e29/.e6 - 1) * .e24 * .e19/.e24
  .e94 <- .e65 * .e28
  .e97 <- (.e74 - 1) * .e29
  .e99 <- 0.5 * ((0.5 * (.e29/(.e24 * .e6^3)) - (0.5 * (.e54/.e6) + .e24)/.e75) *
    .e24 * .e29) - .e88
  .e100 <- 0.5 * .e34
  .e101 <- .e54/.e75
  .e102 <- .e82/.e46
  .e103 <- .e7^2
  .e104 <- .e65 * .e71
  .e105 <- .e71 * .e20
  .e106 <- .e71 * .e28
  .e119 <- .e53/.e25
  .e121 <- .e6 * .e26 * .e13
  .e123 <- .e29 * .e59/.e6
  .e124 <- .e20 * .e72
  .e125 <- 2 * (.e50 * .e6 * .e13 * .e20/.e25)
  .e126 <- 2 * (.e52 * .e6 * .e13 * .e20/.e25)
  .e127 <- 2 * (.e45 * .e6 * .e13 * .e20/.e25)
  .e128 <- 2 * .e15
  .e129 <- 2 * .e8
  .e130 <- .e89/.e24
  .e132 <- .e31^3 * .e58
  .e134 <- tau * .e77 * .e7
  .e136 <- tau^2 * .e103
  .e140 <- (.e106 * .e13 + 0.5 * (.e34 * .e4 * .e5 * .e26/.e14))/.e76
  .e141 <- .e93 - .e102
  .e142 <- .e93 + .e90
  .e144 <- .e65 * .e85 * .e72
  .e145 <- .e65 * .e20
  .e148 <- ((1 - .e74) * .e29/.e6 - 1) * .e24/.e24
  .e155 <- (.e97/.e6 + 1) * .e24/.e24
  .e156 <- .e71^2
  .e158 <- (.e20/.e46 + .e82 * .e5/.e76) * .e28/.e6
  .e161 <- .e86 * .e19 * .e59 + 2 * .e4
  .e162 <- (0.5 * (((.e129 - .e123)/.e39 - (2 * (.e24 * .e4) - .e89) * .e19/.e75) *
    .e24 * .e19) - 2 * (.e81 * .e4/.e6))/.e6
  .e163 <- .e99/.e6
  .e164 <- (0.5 * (.e38 * .e26/.e14) - .e94 * .e13) * .e4
  .e166 <- .e52/.e25 + .e28 * .e72/.e121
  .e171 <- .e37/.e6
  .e174 <- .e97 * .e59/.e6
  .e175 <- .e6^2
  .e177 <- .e20/.e4 + .e90
  .e178 <- 0.5 * (((.e84 - 2)/.e39 - .e101) * .e24 * .e19/.e6)
  .e179 <- 0.5 * (((2 - .e84)/.e39 + .e101) * .e24 * .e19/.e6)
  .e180 <- 0.5 * (.e36 * .e35)
  .e182 <- 0.5 * (.e37 * .e26/.e14) + .e28 * .e13 * .e72
  .e186 <- 0.5 * .e23 - 0.5 * (0.5 * .e36 + .e23)
  .e188 <- 0.5 * (.e136/.e132) - (0.5 * .e73 - 0.5 * .e134)/.e83
  .e190 <- 1/.e14 - .e119
  .e191 <- 1/.e73
  .e192 <- 1/.e13
  .e193 <- .e125 + .e15
  .e194 <- .e79/.e24
  .e195 <- .e73/.e83
  .e196 <- .e4 * .e26
  .e197 <- .e23 - 0.5 * (.e100 + .e32)
  .e200 <- .e8 - .e126
  .e201 <- .e16 - (.e127 + .e128) * .e4
  hessll <- matrix(0, nrow = nXvar + (nuZUvar - 1) + 1 + 1 + nvZVvar, ncol = nXvar +
    (nuZUvar - 1) + 1 + 1 + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(Xvar * ((.e148 - (.e20/(.e5 * .e26 * .e13) +
    .e82/.e76) * .e28 * .e4/.e6)/.e6) * wHvar, Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + (nuZUvar - 1))] <- crossprod(S * Xvar *
    (((2/.e14 - .e119) * .e4 - .e71 * (.e20/.e5 + .e102)/.e6) * .e28/.e26 - ((.e84 -
      1) * .e59 + (.e4 - .e130) * .e19/.e6) * .e24/.e39) * wHvar, uHvar[, -1,
    drop = FALSE])
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 1] <- crossprod(S * Xvar, (-((.e148 +
    .e158) * .e7/.e6)) * wHvar)
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 2] <- crossprod(S * Xvar, ((.e178 - ((.e50/.e25 -
    .e94/.e121) * .e4 - (.e145/.e5 + .e192)/.e6) * .e28/.e26) * .e4) * wHvar)
  hessll[1:nXvar, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 +
    nvZVvar)] <- crossprod(S * Xvar * ((.e178 - (.e166 * .e4 + .e124/(.e6 * .e5)) *
    .e28/.e26) * .e5) * wHvar, vHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + 1):(nXvar + (nuZUvar - 1))] <- crossprod(uHvar[,
    -1, drop = FALSE] * (((.e16 - (.e156 * .e20 + 4 * .e18))/.e14 - (((.e34 *
    .e197 * .e5/.e13 + 2 * .e45) * .e20 + .e45 * (2 * .e16 - (.e127 + 4 * .e15) *
    .e4)) * .e4/.e25 + .e156 * .e28/.e26)) * .e28/.e26 - ((((.e79 * .e59^2/.e24 -
    (.e34 * .e19 + .e8) * .e4)/.e6 + .e8) * .e19 + (.e8 - .e123) * .e59) * .e24/.e24 +
    (2 - 2 * ((.e130 + .e4)/.e6)) * .e4)/.e6) * wHvar, uHvar[, -1, drop = FALSE])
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar - 1) + 1] <- crossprod(uHvar[,
    -1, drop = FALSE], (((((1/.e46 - .e140) * .e5 - .e105/.e196) * .e28 - ((.e174 +
    .e129 + .e15) * .e24/.e24 + 2 * ((.e90 - .e194) * .e4/.e6)))/.e6 + .e77 *
    (.e195 - .e191)) * .e7) * wHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar - 1) + 2] <- crossprod(uHvar[,
    -1, drop = FALSE], ((.e162 + 2 * (.e88 - .e94/.e26) - ((.e50 * .e201 + .e186 *
    .e34 * .e5 * .e20/.e13 - S * .e45 * .e4 * .e10)/.e25 - .e104 * .e85) * .e28/.e26) *
    .e4 + 0.5 * (tau * (.e191 - .e195) * .e77 * .e7)) * wHvar)
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar +
    (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(uHvar[, -1, drop = FALSE] * ((.e162 +
    .e28 * (tau * .e190 * .e7 - (((0.5 * (.e35 * (2 - .e100) + 2 * (.e11/.e175)) +
      0.5 * (.e40/.e6)) * .e4 * .e20/.e13 + .e52 * .e201)/.e25 + .e71 * .e85 *
      .e72))/.e26) * .e5) * wHvar, vHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 1] <- sum(((.e77 *
    (.e77/.e83 + .e8/.e132) - (.e155 + (.e20/(.e196 * .e13) + .e87/.e76) * .e28 *
    .e5/.e6)/.e6) * .e103) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 2] <- sum((((.e179 -
    (.e50 * .e5/.e25 - .e65 * .e177/.e6) * .e28/.e26) * .e4 + 0.5 * (((1 - .e136/.e31^2)/.e73 -
    .e134/.e83) * .e77)) * .e7) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 1, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)] <- crossprod(vHvar, ((((.e192 - .e124/.e4)/.e6 - .e166 *
    .e5) * .e28/.e26 + .e179) * .e7 * .e5) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 2, nXvar + (nuZUvar - 1) + 2] <- sum((((.e99 *
    .e4 + .e78 - 0.5)/.e6 - ((.e50 * (.e16 - (.e125 + 3 * .e15) * .e4) + .e186 *
    .e36 * .e5 * .e20/.e13)/.e25 + .e65^2 * .e85 * .e4 + .e61) * .e28/.e26) *
    .e4 + 0.5 * (tau * .e188 * .e77 * .e7)) * wHvar)
  hessll[nXvar + (nuZUvar - 1) + 2, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)] <- crossprod(vHvar, (((.e144 - ((0.5 * (.e38/.e6) + 0.5 *
    ((.e23 - 1) * .e5/.e6 + 1 - .e180)) * .e20/.e13 + tau * .e50 * .e7 - .e52 *
    .e193)/.e25) * .e28/.e26 + .e163) * .e4 * .e5) * wHvar)
  hessll[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar),
    (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(vHvar *
    (((.e99 * .e5 + .e78 - 0.5)/.e6 + .e28 * (.e62 - ((((3 * .e8 - .e126) * .e5 -
      .e18) * .e52 + (0.5 * .e33 - 0.5 * (0.5 * .e35 + .e33)) * .e35 * .e4 *
      .e20/.e13)/.e25 + .e85 * .e5 * .e72^2))/.e26) * .e5) * wHvar, vHvar)
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
