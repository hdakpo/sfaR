########################################################
#                                                      #
# truncated normal + normal distributions + scaling    #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

ctruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar - 1)]
  tau <- parm[nXvar + length(delta) + 1]
  cu <- parm[nXvar + length(delta) + 2]
  phi <- parm[(nXvar + length(delta) + 2 + 1):(nXvar + length(delta) +
    2 + nvZVvar)]
  musca <- exp(as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))) * tau
  Wusca <- cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))
  Wvsca <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (musca * exp(Wvsca) - exp(Wusca) * S * epsilon)/(exp(Wusca) +
    exp(Wvsca))
  sigmastar <- sqrt(exp(Wusca) * exp(Wvsca)/(exp(Wusca) + exp(Wvsca)))
  ll <- (-1/2 * log(exp(Wusca) + exp(Wvsca)) + dnorm((musca +
    S * epsilon)/sqrt(exp(Wusca) + exp(Wvsca)), log = TRUE) +
    pnorm(mustar/sigmastar, log.p = TRUE) - pnorm(musca/sqrt(exp(Wusca)),
    log.p = TRUE))
  return(ll)
}

# starting value for the log-likelihood ----------

csttruncnormscal <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  if (S * m3 > 0) {
    ## Coelli (1995) suggests 0.05 for gamma
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
  reg_hetu <- lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar]))
  if (any(is.na(reg_hetu$coefficients)))
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar]))
  }
  if (any(is.na(reg_hetv$coefficients)))
    stop("at least one of the OLS coefficients of 'vhet' is NA: ",
      paste(colnames(vHvar)[is.na(reg_hetv$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
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

cgradtruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar - 1)]
  tau <- parm[nXvar + length(delta) + 1]
  cu <- parm[nXvar + length(delta) + 2]
  phi <- parm[(nXvar + length(delta) + 2 + 1):(nXvar + length(delta) +
    2 + nvZVvar)]
  musca <- exp(as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))) * tau
  Wusca <- cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))
  Wvsca <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Hsca <- as.numeric(crossprod(matrix(delta), t(uHvar[, -1])))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  sigma_sq <- exp(Wusca) + exp(Wvsca)
  muwv <- musca * exp(Wvsca)
  mustar <- (muwv - exp(Wusca) * S * epsilon)/(sigma_sq)
  sigmastar <- sqrt(exp(Wusca) * exp(Wvsca)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dmustar2 <- dnorm((S * (epsilon) + musca)/sqrt(sigma_sq))
  dmu <- dnorm(musca/exp((Wusca)/2))
  pmu <- pnorm(musca/exp((Wusca)/2))
  muwvwu <- muwv - S * exp(Wusca) * (epsilon)
  mustar3 <- (muwvwu)/((sigma_sq) * sigmastar)^2
  mustar4 <- (muwv - 2 * (S * exp(Wusca) * (epsilon)))/((sigma_sq) *
    sigmastar)
  wusq <- exp(Wusca)/(sigma_sq)
  wvsq <- exp(Wvsca)/(sigma_sq)
  dpmusig <- dmusig/pmusig
  pmusigx2 <- pmusig * sigmastar
  wupmu <- exp((Wusca)/2) * pmu
  muepsi <- S * (epsilon) + musca
  muepsix2 <- musca - exp(Wusca) * (muepsi)/(sigma_sq)
  dmusigwu <- dmusig * exp(Wusca)
  dmusigwv <- dmusig * exp(Wvsca)
  sigx2 <- (sigma_sq) * sigmastar
  dmuepsi <- dmustar2 * (muepsi)^2
  dmusq <- (dmustar2 * (sigma_sq))
  sigx3 <- (0.5 * (dmuepsi/dmusq) - 0.5)/(sigma_sq)
  sigx4 <- 0.5 * ((2 - 2 * (wusq)) * exp(Wvsca)/sigmastar) +
    2 * sigmastar
  sigx5 <- 0.5 * ((1 - wvsq) * exp(Wusca)/sigmastar) + sigmastar
  sigx6 <- 0.5 * ((1 - wusq) * exp(Wvsca)/sigmastar) + sigmastar
  sigx7 <- (sigx6) * mustar3 + S * (epsilon)/(sigx2)
  sigx8 <- musca/(sigx2) - (sigx5) * mustar3
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = S * ((muepsi) +
    dmusigwu/(pmusigx2))/(sigma_sq), FUN = "*"), sweep(matrix(uHvar[,
    -1], ncol = nuZUvar - 1), MARGIN = 1, STATS = ((mustar4 -
    (sigx4) * exp(Wusca) * mustar3) * dpmusig - ((muepsi) *
    (muepsix2) + exp(Wusca))/(sigma_sq)), FUN = "*"), ((dmusigwv/(pmusigx2) -
    (muepsi))/(sigma_sq) - dmu/(wupmu)) * exp(Hsca), (sigx3 -
    (sigx7) * dpmusig) * exp(Wusca) + 0.5 * (tau * dmu *
    exp(Hsca)/(wupmu)), sweep(vHvar, MARGIN = 1, STATS = (sigx3 +
    dmusig * (sigx8)/pmusig) * exp(Wvsca), FUN = "*")))
  return(gradll)
}

# Hessian of the likelihood function ----------

chesstruncnormscalike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + (nuZUvar - 1))]
  tau <- parm[nXvar + (nuZUvar - 1) + 1]
  cu <- parm[nXvar + (nuZUvar - 1) + 2]
  phi <- parm[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)]
  musca <- exp(as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))) * tau
  Wusca <- cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))
  Wvsca <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  Hsca <- as.numeric(crossprod(matrix(delta), t(uHvar[, -1])))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  sigma_sq <- exp(Wusca) + exp(Wvsca)
  muwv <- musca * exp(Wvsca)
  mustar <- (muwv - exp(Wusca) * S * epsilon)/(sigma_sq)
  sigmastar <- sqrt(exp(Wusca) * exp(Wvsca)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  dmustar2 <- dnorm((S * (epsilon) + musca)/sqrt(sigma_sq))
  dmu <- dnorm(musca/exp((Wusca)/2))
  pmu <- pnorm(musca/exp((Wusca)/2))
  muwvwu <- muwv - S * exp(Wusca) * (epsilon)
  mustar3 <- (muwvwu)/((sigma_sq) * sigmastar)^2
  mustar4 <- (muwv - 2 * (S * exp(Wusca) * (epsilon)))/((sigma_sq) *
    sigmastar)
  wusq <- exp(Wusca)/(sigma_sq)
  wvsq <- exp(Wvsca)/(sigma_sq)
  dpmusig <- dmusig/pmusig
  pmusigx2 <- pmusig * sigmastar
  wupmu <- exp((Wusca)/2) * pmu
  muepsi <- S * (epsilon) + musca
  muepsix2 <- musca - exp(Wusca) * (muepsi)/(sigma_sq)
  dmusigwu <- dmusig * exp(Wusca)
  dmusigwv <- dmusig * exp(Wvsca)
  sigx2 <- (sigma_sq) * sigmastar
  dmuepsi <- dmustar2 * (muepsi)^2
  dmusq <- (dmustar2 * (sigma_sq))
  sigx3 <- (0.5 * (dmuepsi/dmusq) - 0.5)/(sigma_sq)
  sigx4 <- 0.5 * ((2 - 2 * (wusq)) * exp(Wvsca)/sigmastar) +
    2 * sigmastar
  sigx5 <- 0.5 * ((1 - wvsq) * exp(Wusca)/sigmastar) + sigmastar
  sigx6 <- 0.5 * ((1 - wusq) * exp(Wvsca)/sigmastar) + sigmastar
  sigx7 <- (sigx6) * mustar3 + S * (epsilon)/(sigx2)
  sigx8 <- musca/(sigx2) - (sigx5) * mustar3
  dmuepsix2 <- dmustar2 * (muepsi)
  sigx9 <- 0.5 * ((((muepsi)^2/(sigma_sq) - 2)/dmusq - dmuepsi/dmusq^2) *
    dmuepsix2/(sigma_sq))
  sigx10 <- (0.5 * (((2 * (musca) - (muepsi)^2 * (muepsix2)/(sigma_sq))/dmusq -
    (2 * (dmustar2 * exp(Wusca)) - dmustar2 * (muepsi) *
      (muepsix2)) * (muepsi)/dmusq^2) * dmuepsix2) - 2 *
    ((0.5 * (dmuepsi/dmusq) - 0.5) * wusq))/(sigma_sq)
  sigx11 <- 0.5 * (((2 - (muepsi)^2/(sigma_sq))/dmusq + dmuepsi/dmusq^2) *
    dmuepsix2/(sigma_sq))
  sigx12 <- (sigx5)/(sigx2)^2 + dmusig * (sigx8)/((sigma_sq) *
    pmusigx2)
  sigx13 <- (0.5 * ((0.5 * ((muepsi)^2/(dmustar2 * (sigma_sq)^3)) -
    (0.5 * (dmuepsi/(sigma_sq)) + dmustar2)/dmusq^2) * dmuepsi) -
    sigx3)
  sigx14 <- 0.5 * (dmuepsi/dmusq)
  sigx15 <- mustar4 - (sigx4) * exp(Wusca) * mustar3
  sigx16 <- (muwvwu)/(sigx2) + dpmusig
  hessll <- matrix(nrow = nXvar + (nuZUvar - 1) + 1 + 1 + nvZVvar,
    ncol = nXvar + (nuZUvar - 1) + 1 + 1 + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * ((0 - 1) - ((muwvwu)/(exp(Wvsca) * pmusigx2) +
      dmusigwu/(pmusigx2)^2) * dmusig * wusq)/(sigma_sq),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + (nuZUvar - 1))] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (((2/(sigx2) - (sigx4) * exp(Wusca)/(sigx2)^2) *
      exp(Wusca) - (sigx15) * ((muwvwu)/exp(Wvsca) + dmusigwu/(pmusigx2))/(sigma_sq)) *
      dpmusig - (((muepsi)^2/(sigma_sq) - 1) * (muepsix2) +
      (exp(Wusca) - (muepsi) * (muepsix2)) * (muepsi)/(sigma_sq))/((sigma_sq))),
    FUN = "*"), uHvar[, -1])
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), 1:nXvar] <- t(hessll[1:nXvar,
    (nXvar + 1):(nXvar + (nuZUvar - 1))])
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 1] <- hessll[nXvar +
    (nuZUvar - 1) + 1, 1:nXvar] <- (-(S * ((0 - 1) + ((muwvwu)/(pmusigx2) +
    dmusigwu * exp(Wvsca)/(pmusigx2)^2) * dmusig/(sigma_sq)) *
    exp(Hsca)/(sigma_sq))) %*% Xvar
  hessll[1:nXvar, nXvar + (nuZUvar - 1) + 2] <- hessll[nXvar +
    (nuZUvar - 1) + 2, 1:nXvar] <- (S * (sigx9 - (((sigx6)/(sigx2)^2 -
    (sigx7) * dmusig/((sigma_sq) * pmusigx2)) * exp(Wusca) -
    ((sigx7) * (muwvwu)/exp(Wvsca) + 1/sigmastar)/(sigma_sq)) *
    dpmusig) * exp(Wusca)) %*% Xvar
  hessll[1:nXvar, (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S * (sigx9 - ((sigx12) * exp(Wusca) + (muwvwu) *
      (sigx8)/((sigma_sq) * exp(Wvsca))) * dpmusig) * exp(Wvsca),
    FUN = "*"), vHvar)
  hessll[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar), 1:nXvar] <- t(hessll[1:nXvar, (nXvar +
    (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)])
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + 1):(nXvar +
    (nuZUvar - 1))] <- crossprod(sweep(matrix(uHvar[, -1],
    ncol = (nuZUvar - 1)), MARGIN = 1, STATS = (((muwv - ((sigx15)^2 *
    (muwvwu) + 4 * (S * exp(Wusca) * (epsilon))))/(sigx2) -
    ((((2 - 2 * (wusq)) * (wusq - 0.5 * (0.5 * (2 - 2 * (wusq)) +
      2 * (wusq))) * exp(Wvsca)/sigmastar + 2 * (sigx4)) *
      (muwvwu) + (sigx4) * (2 * (muwv) - (2 * ((sigx4) *
      sigx2 * mustar3) + 4 * (S * (epsilon))) * exp(Wusca))) *
      exp(Wusca)/(sigx2)^2 + (sigx15)^2 * dpmusig)) * dpmusig -
    ((((dmustar2 * (muepsi) * (muepsix2)^2/dmustar2 - ((2 -
      2 * (wusq)) * (muepsi) + musca) * exp(Wusca))/(sigma_sq) +
      musca) * (muepsi) + (musca - (muepsi)^2 * (muepsix2)/(sigma_sq)) *
      (muepsix2)) + (2 - 2 * ((dmustar2 * (muepsi) * (muepsix2)/dmustar2 +
      exp(Wusca))/(sigma_sq))) * exp(Wusca))/(sigma_sq)),
    FUN = "*"), uHvar[, -1])
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar -
    1) + 1] <- hessll[nXvar + (nuZUvar - 1) + 1, (nXvar +
    1):(nXvar + (nuZUvar - 1))] <- (((((1/(pmusigx2) - ((sigx15) *
    dmusig * sigmastar + 0.5 * ((2 - 2 * (wusq)) * exp(Wusca) *
    exp(Wvsca) * pmusig/(sigx2)))/(pmusigx2)^2) * exp(Wvsca) -
    (sigx15) * (muwvwu)/(exp(Wusca) * pmusig)) * dmusig -
    ((2 * (musca) + S * (epsilon)) + 2 * ((dmusigwv/(pmusigx2) -
      dmustar2 * (muepsi)/dmustar2) * wusq)))/(sigma_sq) +
    dmu * (wupmu/(wupmu)^2 - 1/(wupmu))) * exp(Hsca)) %*%
    uHvar[, -1]
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), nXvar + (nuZUvar -
    1) + 2] <- hessll[nXvar + (nuZUvar - 1) + 2, (nXvar +
    1):(nXvar + (nuZUvar - 1))] <- (((sigx10 + 2 * (sigx3 -
    (sigx7) * dpmusig) - (((sigx6) * (muwv - (2 * ((sigx4) *
    sigx2 * mustar3) + 2 * (S * (epsilon))) * exp(Wusca)) +
    (0.5 * (wusq) - 0.5 * (0.5 * (1 - wusq) + wusq)) * (2 -
      2 * (wusq)) * exp(Wvsca) * (muwvwu)/sigmastar - S *
    (sigx4) * exp(Wusca) * (epsilon))/(sigx2)^2 - (sigx7) *
    (sigx15) * (sigx16)) * dpmusig) * exp(Wusca) + 0.5 *
    (tau * (1/(wupmu) - wupmu/(wupmu)^2) * dmu * exp(Hsca)))) %*%
    uHvar[, -1]
  hessll[(nXvar + 1):(nXvar + (nuZUvar - 1)), (nXvar + (nuZUvar -
    1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(sweep(matrix(uHvar[,
    -1], ncol = (nuZUvar - 1)), MARGIN = 1, STATS = (sigx10 +
    dmusig * (tau * (1/(sigx2) - (sigx4) * exp(Wusca)/(sigx2)^2) *
      exp(Hsca) - (((0.5 * ((1 - wvsq) * (2 - 0.5 * (2 -
      2 * (wusq))) + 2 * (exp(Wusca) * wvsq/(sigma_sq))) +
      0.5 * ((2 - 2 * (wusq)) * wvsq)) * exp(Wusca) * (muwvwu)/sigmastar +
      (sigx5) * (muwv - (2 * ((sigx4) * sigx2 * mustar3) +
        2 * (S * (epsilon))) * exp(Wusca)))/(sigx2)^2 +
      (sigx15) * (sigx16) * (sigx8)))/pmusig) * exp(Wvsca),
    FUN = "*"), vHvar)
  hessll[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar), (nXvar + 1):(nXvar + (nuZUvar - 1))] <- t(hessll[(nXvar +
    1):(nXvar + (nuZUvar - 1)), (nXvar + (nuZUvar - 1) + 2 +
    1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)])
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 1] <- sum((dmu *
    (dmu/(wupmu)^2 + musca/(exp((Wusca)/2)^3 * pmu)) - (((dmustar2/dmustar2 -
    1) * (muepsi)^2/(sigma_sq) + 1) + ((muwvwu)/(exp(Wusca) *
    pmusigx2) + dmusigwv/(pmusigx2)^2) * dmusig * wvsq)/(sigma_sq)) *
    exp(Hsca)^2)
  hessll[nXvar + (nuZUvar - 1) + 1, nXvar + (nuZUvar - 1) + 2] <- hessll[nXvar +
    (nuZUvar - 1) + 2, nXvar + (nuZUvar - 1) + 1] <- sum(((sigx11 -
    ((sigx6) * exp(Wvsca)/(sigx2)^2 - (sigx7) * ((muwvwu)/exp(Wusca) +
      dmusigwv/(pmusigx2))/(sigma_sq)) * dpmusig) * exp(Wusca) +
    0.5 * (((1 - tau^2 * exp(Hsca)^2/exp((Wusca)/2)^2)/(wupmu) -
      tau * dmu * exp(Hsca)/(wupmu)^2) * dmu)) * exp(Hsca))
  hessll[nXvar + (nuZUvar - 1) + 1, (nXvar + (nuZUvar - 1) +
    2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)] <- hessll[(nXvar +
    (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar),
    nXvar + (nuZUvar - 1) + 1] <- ((((1/sigmastar - (muwvwu) *
    (sigx8)/exp(Wusca))/(sigma_sq) - (sigx12) * exp(Wvsca)) *
    dpmusig + sigx11) * exp(Hsca) * exp(Wvsca)) %*% vHvar
  hessll[nXvar + (nuZUvar - 1) + 2, nXvar + (nuZUvar - 1) + 2] <- sum(((sigx13 *
    exp(Wusca) + sigx14 - 0.5)/(sigma_sq) - (((sigx6) * (muwv -
    (2 * ((sigx6) * sigx2 * mustar3) + 3 * (S * (epsilon))) *
      exp(Wusca)) + (0.5 * (wusq) - 0.5 * (0.5 * (1 - wusq) +
    wusq)) * (1 - wusq) * exp(Wvsca) * (muwvwu)/sigmastar)/(sigx2)^2 +
    (sigx7)^2 * (sigx16) * exp(Wusca) + S * (epsilon)/(sigx2)) *
    dpmusig) * exp(Wusca) + 0.5 * (tau * (0.5 * (tau^2 *
    exp(Hsca)^2/(exp((Wusca)/2)^3 * pmu)) - (0.5 * (wupmu) -
    0.5 * (tau * dmu * exp(Hsca)))/(wupmu)^2) * dmu * exp(Hsca)))
  hessll[nXvar + (nuZUvar - 1) + 2, (nXvar + (nuZUvar - 1) +
    2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar)] <- hessll[(nXvar +
    (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar - 1) + 2 + nvZVvar),
    nXvar + (nuZUvar - 1) + 2] <- ((((sigx7) * (sigx16) *
    (sigx8) - ((0.5 * ((1 - wusq) * wvsq) + 0.5 * ((wusq -
    1) * wvsq + 1 - 0.5 * ((1 - wusq) * (1 - wvsq)))) * (muwvwu)/sigmastar +
    tau * (sigx6) * exp(Hsca) - (sigx5) * (2 * ((sigx6) *
    sigx2 * mustar3) + S * (epsilon)))/(sigx2)^2) * dpmusig +
    sigx13/(sigma_sq)) * exp(Wusca) * exp(Wvsca)) %*% vHvar
  hessll[(nXvar + (nuZUvar - 1) + 2 + 1):(nXvar + (nuZUvar -
    1) + 2 + nvZVvar), (nXvar + (nuZUvar - 1) + 2 + 1):(nXvar +
    (nuZUvar - 1) + 2 + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = ((sigx13 * exp(Wvsca) + sigx14 -
      0.5)/(sigma_sq) + dmusig * (musca/(sigx2) - ((((3 *
      (musca) - 2 * ((sigx5) * sigx2 * mustar3)) * exp(Wvsca) -
      S * exp(Wusca) * (epsilon)) * (sigx5) + (0.5 * (wvsq) -
      0.5 * (0.5 * (1 - wvsq) + wvsq)) * (1 - wvsq) * exp(Wusca) *
      (muwvwu)/sigmastar)/(sigx2)^2 + (sigx16) * exp(Wvsca) *
      (sigx8)^2))/pmusig) * exp(Wvsca), FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

truncnormscalAlgOpt <- function(start, olsParam, dataTable, S,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csttruncnormscal(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(ctruncnormscalike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradtruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ctruncnormscalike,
    grad = cgradtruncnormscalike, hess = chesstruncnormscalike,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S), sr1 = trust.optim(x = startVal,
    fn = function(parm) -sum(ctruncnormscalike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradtruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ctruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chesstruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(ctruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradtruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chesstruncnormscalike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ctruncnormscalike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgradtruncnormscalike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chesstruncnormscalike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradtruncnormscalike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S))
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
      mleObj$hessian <- chesstruncnormscalike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chesstruncnormscalike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- ctruncnormscalike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradtruncnormscalike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

ctruncnormscaleff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  phi <- object$mlParam[(object$nXvar + (object$nuZUvar - 1) +
    2 + 1):(object$nXvar + (object$nuZUvar - 1) + 2 + object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  musca <- exp(as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))) * tau
  Wusca <- cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1])))
  Wvsca <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
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
    res <- bind_cols(u = u, uLB = uLB, uUB = uUB, teJLMS = teJLMS,
      m = m, teMO = teMO, teBC = teBC, teBCLB = teBCLB,
      teBCUB = teBCUB)
  } else {
    res <- bind_cols(u = u, uLB = uLB, uUB = uUB, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargtruncnormscal_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  hi <- exp(as.numeric(crossprod(matrix(delta), t(uHvar[, -1]))))
  Lambda <- tau/exp(cu/2)
  m1 <- exp(cu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  margEff <- kronecker(matrix(delta, nrow = 1), matrix(m1 *
    hi, ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmargtruncnormscal_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    (object$nuZUvar - 1))]
  tau <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    1]
  cu <- object$mlParam[object$nXvar + (object$nuZUvar - 1) +
    2]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  hi2 <- exp(2 * as.numeric(crossprod(matrix(delta), t(uHvar[,
    -1]))))
  Lambda <- tau/exp(cu/2)
  m2 <- exp(cu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) -
    (dnorm(Lambda)/pnorm(Lambda))^2)
  margEff <- kronecker(matrix(2 * delta, nrow = 1), matrix(m2 *
    hi2))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
