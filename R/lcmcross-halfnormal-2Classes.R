########################################################
#                                                      #
# Latent Class Model (LCM 2C) normal - half normal     #
#                                                      #
########################################################

# Log-likelihood ----------

cLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  L <- Probc1 * Pi1 + Probc2 * Pi2
  ifelse(L <= 0, return(NA), return(log(L)))
}

# starting value for the log-likelihood ----------

csLCMfhalfnorm2C <- function(olsObj, epsiRes, nXvar, nuZUvar, 
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, Zvar, nZHvar, itermax, 
  printInfo, tol) {
  cat("Initialization: SFA + halfnormal - normal distribution...\n")
  initHalf <- maxLik(logLik = chalfnormlike, start = csthalfnorm(olsObj = olsObj, 
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[, 
      1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])), 
    grad = cgradhalfnormlike, method = "BFGS", control = list(iterlim = itermax, 
      printLevel = if (printInfo) 2 else 0, reltol = tol), 
    nXvar = nXvar, nuZUvar = 1, nvZVvar = 1, uHvar = as.matrix(uHvar[, 
      1]), vHvar = as.matrix(vHvar[, 1]), Yvar = Yvar, 
    Xvar = Xvar, S = S)
  Esti <- initHalf$estimate
  StartVal <- c(Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 
    1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 
    1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], Esti[nXvar + 
    1], if (nuZUvar > 1) rep(0, nuZUvar - 1), Esti[nXvar + 
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar], 
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), 
    paste0("Cl1_", colnames(Zvar)))
  names(initHalf$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
    colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------

cgradLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  sigma_sq1 <- exp(Wu1) + exp(Wv1)
  sigma_sq2 <- exp(Wu2) + exp(Wv2)
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(sigma_sq1))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(sigma_sq2))
  dmusig1 <- dnorm(-(S * exp(Wu1) * (epsilon1)/((sigma_sq1) * 
    sigmastar1)))
  dmusig2 <- dnorm(-(S * exp(Wu2) * (epsilon2)/((sigma_sq2) * 
    sigmastar2)))
  pmusig1 <- pnorm(-(S * exp(Wu1) * (epsilon1)/((sigma_sq1) * 
    sigmastar1)))
  pmusig2 <- pnorm(-(S * exp(Wu2) * (epsilon2)/((sigma_sq2) * 
    sigmastar2)))
  depsisq1 <- dnorm(S * (epsilon1)/sqrt(sigma_sq1))
  depsisq2 <- dnorm(S * (epsilon2)/sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsisq1 * exp(Wu1)/sigmastar1 + S * 
    depsisq1 * pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsisq2 * exp(Wu2)/sigmastar2 + S * 
    depsisq2 * pmusig2 * (epsilon2))
  sqsq1 <- ((sigma_sq1) * sigmastar1)
  sqsq2 <- ((sigma_sq2) * sigmastar2)
  sigx2_1 <- (0.5 * ((1 - exp(Wu1)/(sigma_sq1)) * exp(Wv1)/sigmastar1) + 
    sigmastar1)
  sigx2_2 <- (0.5 * ((1 - exp(Wu2)/(sigma_sq2)) * exp(Wv2)/sigmastar2) + 
    sigmastar2)
  sigx3_1 <- (0.5 * ((1 - exp(Wv1)/(sigma_sq1)) * exp(Wu1)/sigmastar1) + 
    sigmastar1)
  sigx3_2 <- (0.5 * ((1 - exp(Wv2)/(sigma_sq2)) * exp(Wu2)/sigmastar2) + 
    sigmastar2)
  wzdeno <- (1 + exp(Wz))
  prC <- (1 - exp(Wz)/wzdeno)
  wzdsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzdsq2 <- (wzdeno * sqrt(sigma_sq2))
  wzlogit <- (prC * depsisq2 * pmusig2/sqrt(sigma_sq2))
  sigx4 <- (2 * wzlogit + 2 * (depsisq1 * exp(Wz) * pmusig1/wzdsq1))
  sigsq_1 <- (sigx4 * sqrt(sigma_sq1))
  sigsq_2 <- (sigx4 * sqrt(sigma_sq2))
  wdpdsq <- (wzdeno * depsisq1 * pmusig1/wzdsq1^2)
  dpepsisq <- (S * depsisq1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5 <- (wzdeno * sigx4 * (sigma_sq1) * sqrt(sigma_sq1))
  sigx6 <- (S * (0.5 * dpepsisq - (1/sqsq1 - sigx2_1 * exp(Wu1)/sqsq1^2) * 
    dmusig1 * depsisq1) * (epsilon1)/wzdeno - 0.5 * wdpdsq)
  sigx7 <- (sigx3_1 * dmusig1 * depsisq1 * exp(Wu1)/sqsq1^2 + 
    0.5 * dpepsisq) * (epsilon1)
  s3q <- sigx4 * (sigma_sq2) * sqrt(sigma_sq2)
  sigx8 <- (1/wzdsq1 - exp(Wz) * sqrt(sigma_sq1)/wzdsq1^2)
  sigx9 <- sigx8 * depsisq1 * pmusig1
  sigx10 <- (1/sqsq2 - sigx2_2 * exp(Wu2)/sqsq2^2)
  sigx11 <- (S * depsisq2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  dpsq2 <- depsisq2 * pmusig2/(sigma_sq2)
  sigx12 <- (S * (0.5 * sigx11 - sigx10 * dmusig2 * depsisq2) * 
    (epsilon2) - 0.5 * (dpsq2))
  sigx13 <- (S * (sigx3_2 * dmusig2 * depsisq2 * exp(Wu2)/sqsq2^2 + 
    0.5 * sigx11) * (epsilon2) - 0.5 * (dpsq2))
  sigx14 <- (prC * depsisq2 * pmusig2/wzdsq2)
  sigx15 <- (2 * (sigx9) - 2 * sigx14) * exp(Wz)
  sigx26 <- (S * sigx7/wzdeno - 0.5 * wdpdsq)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = 2 * (S * 
    sigx1_1 * exp(Wz)/sigx5), FUN = "*"), sweep(uHvar, MARGIN = 1, 
    STATS = 2 * (exp(Wu1) * exp(Wz) * sigx6/sigsq_1), FUN = "*"), 
    sweep(vHvar, MARGIN = 1, STATS = 2 * (exp(Wv1) * exp(Wz) * 
      sigx26/sigsq_1), FUN = "*"), sweep(Xvar, MARGIN = 1, 
      STATS = 2 * (S * prC * sigx1_2/(s3q)), FUN = "*"), 
    sweep(uHvar, MARGIN = 1, STATS = 2 * (prC * exp(Wu2) * 
      sigx12/sigsq_2), FUN = "*"), sweep(vHvar, MARGIN = 1, 
      STATS = 2 * (prC * exp(Wv2) * sigx13/sigsq_2), FUN = "*"), 
    sweep(Zvar, MARGIN = 1, STATS = sigx15/sigx4, FUN = "*"))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessLCMhalfnormlike2C <- function(parm, nXvar, nuZUvar, nvZVvar, 
  uHvar, vHvar, Yvar, Xvar, S, Zvar, nZHvar) {
  beta1 <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi1 <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  beta2 <- parm[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)]
  delta2 <- parm[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + nvZVvar)]
  phi2 <- parm[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  theta <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wz <- as.numeric(crossprod(matrix(theta), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  sigma_sq1 <- exp(Wu1) + exp(Wv1)
  sigma_sq2 <- exp(Wu2) + exp(Wv2)
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(sigma_sq1))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(sigma_sq2))
  dmusig1 <- dnorm(-(S * exp(Wu1) * (epsilon1)/((sigma_sq1) * 
    sigmastar1)))
  dmusig2 <- dnorm(-(S * exp(Wu2) * (epsilon2)/((sigma_sq2) * 
    sigmastar2)))
  pmusig1 <- pnorm(-(S * exp(Wu1) * (epsilon1)/((sigma_sq1) * 
    sigmastar1)))
  pmusig2 <- pnorm(-(S * exp(Wu2) * (epsilon2)/((sigma_sq2) * 
    sigmastar2)))
  depsisq1 <- dnorm(S * (epsilon1)/sqrt(sigma_sq1))
  depsisq2 <- dnorm(S * (epsilon2)/sqrt(sigma_sq2))
  sigx1_1 <- (dmusig1 * depsisq1 * exp(Wu1)/sigmastar1 + S * 
    depsisq1 * pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsisq2 * exp(Wu2)/sigmastar2 + S * 
    depsisq2 * pmusig2 * (epsilon2))
  sqsq1 <- ((sigma_sq1) * sigmastar1)
  sqsq2 <- ((sigma_sq2) * sigmastar2)
  sigx2_1 <- (0.5 * ((1 - exp(Wu1)/(sigma_sq1)) * exp(Wv1)/sigmastar1) + 
    sigmastar1)
  sigx2_2 <- (0.5 * ((1 - exp(Wu2)/(sigma_sq2)) * exp(Wv2)/sigmastar2) + 
    sigmastar2)
  sigx3_1 <- (0.5 * ((1 - exp(Wv1)/(sigma_sq1)) * exp(Wu1)/sigmastar1) + 
    sigmastar1)
  sigx3_2 <- (0.5 * ((1 - exp(Wv2)/(sigma_sq2)) * exp(Wu2)/sigmastar2) + 
    sigmastar2)
  wzdeno <- (1 + exp(Wz))
  prC <- (1 - exp(Wz)/wzdeno)
  wzdsq1 <- (wzdeno * sqrt(sigma_sq1))
  wzdsq2 <- (wzdeno * sqrt(sigma_sq2))
  wzlogit <- (prC * depsisq2 * pmusig2/sqrt(sigma_sq2))
  sigx4 <- (2 * wzlogit + 2 * (depsisq1 * exp(Wz) * pmusig1/wzdsq1))
  sigsq_1 <- (sigx4 * sqrt(sigma_sq1))
  sigsq_2 <- (sigx4 * sqrt(sigma_sq2))
  wdpdsq <- (wzdeno * depsisq1 * pmusig1/wzdsq1^2)
  dpepsisq <- (S * depsisq1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  sigx5 <- (wzdeno * sigx4 * (sigma_sq1) * sqrt(sigma_sq1))
  sigx6 <- (S * (0.5 * dpepsisq - (1/sqsq1 - sigx2_1 * exp(Wu1)/sqsq1^2) * 
    dmusig1 * depsisq1) * (epsilon1)/wzdeno - 0.5 * wdpdsq)
  sigx7 <- (sigx3_1 * dmusig1 * depsisq1 * exp(Wu1)/sqsq1^2 + 
    0.5 * dpepsisq) * (epsilon1)
  s3q <- sigx4 * (sigma_sq2) * sqrt(sigma_sq2)
  sigx8 <- (1/wzdsq1 - exp(Wz) * sqrt(sigma_sq1)/wzdsq1^2)
  sigx9 <- sigx8 * depsisq1 * pmusig1
  sigx10 <- (1/sqsq2 - sigx2_2 * exp(Wu2)/sqsq2^2)
  sigx11 <- (S * depsisq2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  dpsq2 <- depsisq2 * pmusig2/(sigma_sq2)
  sigx12 <- (S * (0.5 * sigx11 - sigx10 * dmusig2 * depsisq2) * 
    (epsilon2) - 0.5 * (dpsq2))
  sigx13 <- (S * (sigx3_2 * dmusig2 * depsisq2 * exp(Wu2)/sqsq2^2 + 
    0.5 * sigx11) * (epsilon2) - 0.5 * (dpsq2))
  sigx14 <- (prC * depsisq2 * pmusig2/wzdsq2)
  sigx15 <- (2 * (sigx9) - 2 * sigx14) * exp(Wz)
  sigx16 <- (S * (dmusig1 * exp(Wu1)/sigmastar1 + S * pmusig1 * 
    (epsilon1)) * (epsilon1)/(sigma_sq1) - pmusig1)
  wvsq1 <- exp(Wv1)/(sigma_sq1)
  wusq1 <- exp(Wu1)/(sigma_sq1)
  sigx17 <- (S * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  wvsq2 <- exp(Wv2)/(sigma_sq2)
  sigx18 <- (0.5 * sigx16 - 0.5 * pmusig1) * depsisq1/(sigma_sq1)
  sigx19 <- 0.5 * (wzdeno * sigx1_1/(wzdsq1^2 * (sigma_sq1)))
  sigx20 <- 0.5 * (S * depsisq1 * (S * (0.5 * sigx17 - (1/sqsq1 - 
    sigx2_1 * exp(Wu1)/sqsq1^2) * dmusig1) * (epsilon1) - 
    2 * (pmusig1/(sigma_sq1))) * (epsilon1)/(sigma_sq1)^2)
  sigx21 <- 0.5 * (wzdeno * depsisq2 * pmusig2/wzdsq2^2)
  sigx22 <- (sigx3_2 * dmusig2 * depsisq2 * exp(Wu2)/sqsq2^2 + 
    0.5 * sigx11)
  sigx23 <- (0.5 * sigx11 - sigx10 * dmusig2 * depsisq2)
  sigx24 <- 0.5 * ((S * sigx23 * (epsilon2) - dpsq2)/(sigma_sq2))
  sigx25 <- (S * (dmusig2 * exp(Wu2)/sigmastar2 + S * pmusig2 * 
    (epsilon2)) * (epsilon2)/(sigma_sq2) - pmusig2)
  sigx26 <- (S * sigx7/wzdeno - 0.5 * wdpdsq)
  wusq2 <- exp(Wu2)/(sigma_sq2)
  sigx27 <- depsisq2 * (epsilon2)^2/(sigma_sq2)^2
  sigx28 <- 0.5 * (wzdeno * (S * (0.5 * dpepsisq - (1/sqsq1 - 
    sigx2_1 * exp(Wu1)/sqsq1^2) * dmusig1 * depsisq1) * (epsilon1) - 
    wzdeno^2 * depsisq1 * pmusig1/wzdsq1^2)/wzdsq1^2)
  sigx29 <- dmusig1 * (depsisq1 * exp(Wu1)/exp(Wv1) + depsisq1) * 
    exp(Wu1) * (epsilon1)
  wusqx2 <- exp(Wu1)/sqsq1^2
  sigx30 <- 1/sqsq1 - sigx2_1 * wusqx2
  sigx31 <- (sigx30) * dmusig1 * depsisq1
  sigx32 <- dmusig2 * (depsisq2 * exp(Wu2)/exp(Wv2) + depsisq2) * 
    exp(Wu2) * (epsilon2)
  sigx33 <- (0.5 * (sigx4/sqrt(sigma_sq2)) + 2 * (prC * sigx12))
  sigx34 <- 0.5 * (S * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  sigx35 <- 0.5 * (S * depsisq2 * (S * (sigx34 - sigx10 * dmusig2) * 
    (epsilon2) - 2 * (pmusig2/(sigma_sq2))) * (epsilon2)/(sigma_sq2)^2)
  sigx36 <- (0.5 * (sigx4/sqrt(sigma_sq1)) + 2 * (exp(Wz) * 
    sigx6))
  hessll <- matrix(nrow = 2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    nZHvar, ncol = 2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    nZHvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S^2 * ((depsisq1 * sigx16 + S * sigx29/sqsq1)/sigx5 - 
      2 * (sigx1_1^2 * exp(Wz)/sigx5^2)) * exp(Wz)), FUN = "*"), 
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * (((sigx31 + S * (sigx18 - 
      S * (sigx30) * dmusig1 * (depsisq1 * exp(Wu1)/exp(Wv1) + 
        depsisq1) * (epsilon1)) * (epsilon1)/(sigma_sq1))/wzdeno - 
      sigx19)/sigsq_1 - 2 * (sigx1_1 * exp(Wz) * sigx6/(sigsq_1^2 * 
      wzdeno * (sigma_sq1)))) * exp(Wu1) * exp(Wz)), FUN = "*"), 
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + 
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = 2 * 
    (S * (((S * (sigx18 + S * sigx3_1 * sigx29/sqsq1^2) * 
      (epsilon1)/(sigma_sq1) - sigx3_1 * dmusig1 * depsisq1 * 
      wusqx2)/wzdeno - sigx19)/sigsq_1 - 2 * (sigx1_1 * 
      exp(Wz) * sigx26/(sigsq_1^2 * wzdeno * (sigma_sq1)))) * 
      exp(Wv1) * exp(Wz)), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = -(4 * (S^2 * prC * sigx1_1 * sigx1_2 * (sigma_sq2) * 
      exp(Wz) * sqrt(sigma_sq2)/((s3q)^2 * wzdeno * (sigma_sq1) * 
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_1 * exp(Wu2) * 
      exp(Wz) * sigx12 * sqrt(sigma_sq2)/(sigsq_2^2 * wzdeno * 
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_1 * exp(Wv2) * 
      exp(Wz) * sigx13 * sqrt(sigma_sq2)/(sigsq_2^2 * wzdeno * 
      (sigma_sq1) * sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = S * (2 * sigx8 - 2 * (sigx15/(wzdeno * 
      sigx4 * sqrt(sigma_sq1)))) * sigx1_1 * exp(Wz)/(sigx4 * 
      (sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + 
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = 2 * 
    (((exp(Wu1) * (S * (sigx20 - (0.5 * (S^2 * (sigx30) * 
      depsisq1 * (epsilon1)^2/(sigma_sq1)^2) - (((0.5 * 
      (wusq1) + 1 - 0.5 * (0.5 * (1 - wusq1) + wusq1)) * 
      (1 - wusq1) * exp(Wv1)/sigmastar1 + (2 - 2 * (sigx2_1^2 * 
      exp(Wu1) * (sigma_sq1)/sqsq1^2)) * sigmastar1)/sqsq1^2 + 
      S^2 * (sigx30)^2 * exp(Wu1) * (epsilon1)^2/sqsq1) * 
      depsisq1) * dmusig1) * (epsilon1)/wzdeno - sigx28) + 
      S * (0.5 * dpepsisq - sigx31) * (epsilon1)/wzdeno - 
      0.5 * wdpdsq)/sigsq_1 - sigx36 * exp(Wu1) * sigx6/sigsq_1^2) * 
      exp(Wu1) * exp(Wz)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((1 - wusq1) * 
      exp(Wv1)) - S^2 * sigx3_1 * (sigx30) * exp(Wu1) * 
      (epsilon1)^2)/(sigma_sq1) + 0.5 * ((wusq1 - 1) * 
      wvsq1 + 1 - 0.5 * ((1 - wusq1) * (1 - wvsq1)))) * 
      depsisq1/sigmastar1 + 0.5 * (S^2 * sigx3_1 * depsisq1 * 
      (epsilon1)^2/(sigma_sq1)^2)) * exp(Wu1) + sigx3_1 * 
      (1 - 2 * (sigx2_1 * exp(Wu1) * (sigma_sq1) * sigmastar1/sqsq1^2)) * 
      depsisq1) * dmusig1/sqsq1^2 + sigx20) * (epsilon1)/wzdeno - 
      sigx28)/sigsq_1 - sigx36 * sigx26/sigsq_1^2) * exp(Wu1) * 
      exp(Wv1) * exp(Wz)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_2 * exp(Wu1) * 
      (sigma_sq2) * exp(Wz) * sigx6 * sqrt(sigma_sq2)/((s3q)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * exp(Wu1) * exp(Wu2) * 
      exp(Wz) * sigx6 * sigx12 * sqrt(sigma_sq2)/(sigsq_2^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * exp(Wu1) * exp(Wv2) * 
      exp(Wz) * sigx13 * sigx6 * sqrt(sigma_sq2)/(sigsq_2^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = (2 * 
    ((0.5 * (S^2 * sigx8 * depsisq1 * (epsilon1)^2/(sigma_sq1)^2) - 
      ((0.5/sqrt(sigma_sq1) - wzdeno^2 * sqrt(sigma_sq1)/wzdsq1^2) * 
        exp(Wz) + 0.5 * (wzdeno/sqrt(sigma_sq1))) * depsisq1/wzdsq1^2) * 
      pmusig1 - S * sigx8 * sigx31 * (epsilon1)) - 2 * 
    (sigx15 * sigx6/sigsq_1)) * exp(Wu1) * exp(Wz)/sigx4, 
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (wvsq1) - 0.5 * 
      (0.5 * (1 - wvsq1) + wvsq1)) * (1 - wvsq1) + S^2 * 
      sigx3_1^2 * exp(Wu1) * exp(Wv1) * (epsilon1)^2/(sqsq1^2 * 
      (sigma_sq1))) * depsisq1 * exp(Wu1)/sigmastar1 + 
      ((0.5 * (S^2 * depsisq1 * (epsilon1)^2/(sigma_sq1)^2) - 
        2 * (sigx3_1 * depsisq1 * (sigma_sq1) * sigmastar1/sqsq1^2)) * 
        exp(Wv1) + depsisq1) * sigx3_1) * dmusig1 * wusqx2 + 
      S * (0.5 * (exp(Wv1) * (S * (sigx3_1 * dmusig1 * 
        wusqx2 + 0.5 * sigx17) * (epsilon1) - 2 * (pmusig1/(sigma_sq1)))) + 
        0.5 * pmusig1) * depsisq1 * (epsilon1)/(sigma_sq1)^2) * 
      (epsilon1)/wzdeno - (0.5 * (depsisq1 * pmusig1) + 
      0.5 * (exp(Wv1) * (S * sigx7 - wzdeno^2 * depsisq1 * 
        pmusig1/wzdsq1^2))) * wzdeno/wzdsq1^2)/sigsq_1 - 
      (0.5 * (sigx4/sqrt(sigma_sq1)) + 2 * (exp(Wz) * sigx26)) * 
        exp(Wv1) * sigx26/sigsq_1^2) * exp(Wv1) * exp(Wz)), 
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(4 * 
    (S * prC * sigx1_2 * (sigma_sq2) * exp(Wv1) * exp(Wz) * 
      sigx26 * sqrt(sigma_sq2)/((s3q)^2 * sqrt(sigma_sq1)))), 
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
      nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, 
    STATS = -(4 * (prC * exp(Wu2) * exp(Wv1) * exp(Wz) * 
      sigx26 * sigx12 * sqrt(sigma_sq2)/(sigsq_2^2 * sqrt(sigma_sq1)))), 
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * exp(Wv1) * exp(Wv2) * 
      exp(Wz) * sigx26 * sigx13 * sqrt(sigma_sq2)/(sigsq_2^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = (2 * ((0.5 * (S^2 * sigx8 * depsisq1 * 
      (epsilon1)^2/(sigma_sq1)^2) - ((0.5/sqrt(sigma_sq1) - 
      wzdeno^2 * sqrt(sigma_sq1)/wzdsq1^2) * exp(Wz) + 
      0.5 * (wzdeno/sqrt(sigma_sq1))) * depsisq1/wzdsq1^2) * 
      pmusig1 + S * sigx3_1 * sigx8 * dmusig1 * depsisq1 * 
      exp(Wu1) * (epsilon1)/sqsq1^2) - 2 * (sigx15 * sigx26/sigsq_1)) * 
      exp(Wv1) * exp(Wz)/sigx4, FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S^2 * ((depsisq2 * sigx25 + S * sigx32/sqsq2)/(s3q) - 
      2 * (prC * sigx1_2^2/(s3q)^2)) * prC), FUN = "*"), 
    Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S * ((sigx10 * dmusig2 * depsisq2 + (S * 
      ((0.5 * sigx25 - 0.5 * pmusig2) * depsisq2/(sigma_sq2) - 
        S * sigx10 * dmusig2 * (depsisq2 * exp(Wu2)/exp(Wv2) + 
          depsisq2) * (epsilon2)) * (epsilon2) - 0.5 * 
      (sigx1_2/(sigma_sq2)))/(sigma_sq2))/sigsq_2 - 2 * 
      (prC * sigx1_2 * sigx12/(sigsq_2^2 * (sigma_sq2)))) * 
      prC * exp(Wu2)), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * (((S * ((0.5 * sigx25 - 
      0.5 * pmusig2) * depsisq2/(sigma_sq2) + S * sigx3_2 * 
      sigx32/sqsq2^2) * (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2) - 
      sigx3_2 * dmusig2 * depsisq2 * exp(Wu2)/sqsq2^2)/sigsq_2 - 
      2 * (prC * sigx1_2 * sigx13/(sigsq_2^2 * (sigma_sq2)))) * 
      prC * exp(Wv2)), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(S * prC * (2 * ((2 * (sigx9) - 
      2 * sigx14)/sigx4) + 2/wzdeno) * sigx1_2 * exp(Wz)/(s3q)), 
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((exp(Wu2) * (S * (sigx35 - 
      (0.5 * (S^2 * sigx10 * sigx27) - (((0.5 * (wusq2) + 
        1 - 0.5 * (0.5 * (1 - wusq2) + wusq2)) * (1 - 
        wusq2) * exp(Wv2)/sigmastar2 + (2 - 2 * (sigx2_2^2 * 
        exp(Wu2) * (sigma_sq2)/sqsq2^2)) * sigmastar2)/sqsq2^2 + 
        S^2 * sigx10^2 * exp(Wu2) * (epsilon2)^2/sqsq2) * 
        depsisq2) * dmusig2) * (epsilon2) - sigx24) + 
      S * sigx23 * (epsilon2) - 0.5 * (dpsq2))/sigsq_2 - 
      sigx33 * exp(Wu2) * sigx12/sigsq_2^2) * prC * exp(Wu2)), 
    FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((1 - wusq2) * 
      exp(Wv2)) - S^2 * sigx3_2 * sigx10 * exp(Wu2) * (epsilon2)^2)/(sigma_sq2) + 
      0.5 * ((wusq2 - 1) * wvsq2 + 1 - 0.5 * ((1 - wusq2) * 
        (1 - wvsq2)))) * depsisq2/sigmastar2 + 0.5 * 
      (S^2 * sigx3_2 * sigx27)) * exp(Wu2) + sigx3_2 * 
      (1 - 2 * (sigx2_2 * exp(Wu2) * (sigma_sq2) * sigmastar2/sqsq2^2)) * 
      depsisq2) * dmusig2/sqsq2^2 + sigx35) * (epsilon2) - 
      sigx24)/sigsq_2 - sigx33 * sigx13/sigsq_2^2) * prC * 
      exp(Wu2) * exp(Wv2)), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(prC * (2 * ((2 * (sigx9) - 2 * 
      sigx14) * sigx12/sigx4) + 2 * (S * sigx23 * (epsilon2)/wzdeno - 
      sigx21)) * exp(Wu2) * exp(Wz)/sigsq_2), FUN = "*"), 
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (wvsq2) - 0.5 * 
      (0.5 * (1 - wvsq2) + wvsq2)) * (1 - wvsq2) + S^2 * 
      sigx3_2^2 * exp(Wu2) * exp(Wv2) * (epsilon2)^2/(sqsq2^2 * 
      (sigma_sq2))) * depsisq2 * exp(Wu2)/sigmastar2 + 
      ((0.5 * (S^2 * sigx27) - 2 * (sigx3_2 * depsisq2 * 
        (sigma_sq2) * sigmastar2/sqsq2^2)) * exp(Wv2) + 
        depsisq2) * sigx3_2) * dmusig2 * exp(Wu2)/sqsq2^2 + 
      S * (0.5 * (exp(Wv2) * (S * (sigx3_2 * dmusig2 * 
        exp(Wu2)/sqsq2^2 + sigx34) * (epsilon2) - 2 * 
        (pmusig2/(sigma_sq2)))) + 0.5 * pmusig2) * depsisq2 * 
        (epsilon2)/(sigma_sq2)^2) * (epsilon2) - (0.5 * 
      (depsisq2 * pmusig2) + 0.5 * (exp(Wv2) * (S * sigx22 * 
      (epsilon2) - dpsq2)))/(sigma_sq2))/sigsq_2 - (0.5 * 
      (sigx4/sqrt(sigma_sq2)) + 2 * (prC * sigx13)) * exp(Wv2) * 
      sigx13/sigsq_2^2) * prC * exp(Wv2)), FUN = "*"), 
    vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(prC * 
    (2 * ((2 * (sigx9) - 2 * sigx14) * sigx13/sigx4) + 2 * 
      (S * sigx22 * (epsilon2)/wzdeno - sigx21)) * exp(Wv2) * 
    exp(Wz)/sigsq_2), FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar + nZHvar), (2 * nXvar + 2 * 
    nuZUvar + 2 * nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, 
    STATS = ((2 * (prC * (1/(wzdeno^2 * sqrt(sigma_sq2)) + 
      sqrt(sigma_sq2)/wzdsq2^2) * depsisq2 * pmusig2) - 
      ((2 * (sigx9) - 2 * sigx14)^2/sigx4 + 2 * ((2 - 2 * 
        (wzdeno * (sigma_sq1) * exp(Wz)/wzdsq1^2)) * 
        depsisq1 * pmusig1 * sqrt(sigma_sq1)/wzdsq1^2))) * 
      exp(Wz) + 2 * (sigx9) - 2 * sigx14) * exp(Wz)/sigx4, 
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

LCM2ChnormAlgOpt <- function(start, olsParam, dataTable, S, nXvar, 
  uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, 
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType, 
  qac, initStart, initAlg, initIter, initFactorLB, initFactorUB) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm2C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar, itermax = itermax, 
      tol = tol, printInfo = printInfo)
    InitHalf <- start_st$initHalf
    startVal <- start_st$StartVal
  }
  if (initStart) {
    startMat <- cbind(startVal, initFactorLB * startVal, 
      initFactorUB * startVal)
    startMat <- cbind(startMat, apply(startMat[, 2:3], 1, 
      which.min), apply(startMat[, 2:3], 1, which.max))
    startMat <- cbind(startMat, ifelse(startMat[, 4] == 1, 
      startMat[, 2], startMat[, 3]), ifelse(startMat[, 
      5] == 1, startMat[, 2], startMat[, 3]))
    initModel <- nlminb(start = startVal, objective = function(parm) -sum(cLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), lower = startMat[, 
      6], upper = startMat[, 7], control = list(iter.max = initIter, 
      trace = if (printInfo) 1 else 0, eval.max = initIter, 
      rel.tol = tol, x.tol = tol))
    startVal <- initModel$par
  }
  startLoglik <- sum(cLCMhalfnormlike2C(startVal, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, Zvar = Zvar, 
    nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...), 
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...), 
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal, 
    fn = function(parm) -sum(cLCMhalfnormlike2C(parm, nXvar = nXvar, 
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, Zvar = Zvar, 
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, 
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax, 
      stepmax = stepmax, xtol = tol, grtol = gradtol)), 
    maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike2C, grad = cgradLCMhalfnormlike2C, 
      hess = chessLCMhalfnormlike2C, start = startVal, 
      finalHessian = if (hessianType == 2) "bhhh" else TRUE, 
      control = list(printLevel = if (printInfo) 2 else 0, 
        iterlim = itermax, reltol = tol, tol = tol, qac = qac), 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), sr1 = trust.optim(x = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", 
      control = list(maxit = itermax, cgtol = gradtol, 
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0, 
        report.precision = 1L)), sparse = trust.optim(x = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), 
      method = "Sparse", control = list(maxit = itermax, 
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, 
        report.level = if (printInfo) 2 else 0, report.precision = 1L, 
        preconditioner = 1L)), mla = mla(b = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessLCMhalfnormlike2C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, 
      maxiter = itermax, epsa = gradtol, epsb = gradtol), 
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessLCMhalfnormlike2C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, 
      trace = if (printInfo) 1 else 0, eval.max = itermax, 
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike2C(mleObj$par, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar))
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
      mleObj$hessian <- chessLCMhalfnormlike2C(parm = mleObj$par, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1") 
      mleObj$hessian <- chessLCMhalfnormlike2C(parm = mleObj$solution, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike2C(parm = mlParam, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, Zvar = Zvar, 
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike2C(parm = mlParam, 
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
    S = S, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, 
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) InitHalf = InitHalf))
}

# Posterior probabilities and efficiencies ----------

cLCM2Chalfnormeff <- function(object, level) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + 
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + 
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
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
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, Pcond_c2)
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c, PosteriorProb_c1 = Pcond_c1, 
      PosteriorProb_c2 = Pcond_c2, PriorProb_c1 = Probc1, PriorProb_c2 = Probc2, u_c = u_c, 
      teJLMS_c = teJLMS_c, u_c1 = u_c1, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c, PosteriorProb_c1 = Pcond_c1, 
      PosteriorProb_c2 = Pcond_c2, PriorProb_c1 = Probc1, PriorProb_c2 = Probc2, u_c = u_c, 
      u_c1 = u_c1, u_c2 = u_c2, ineff_c1 = ineff_c1, ineff_c2 = ineff_c2)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargLCM2Chalfnorm_Eu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + 
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + 
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
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
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c2")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in 1:ncol(margEff_c1)) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], 
      margEff_c2[, c])
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c")
  margEff <- bind_cols(as_tibble(margEff_c), as_tibble(margEff_c1), 
    as_tibble(margEff_c2))
  return(margEff)
}

cmargLCM2Chalfnorm_Vu <- function(object) {
  beta1 <- object$mlParam[1:(object$nXvar)]
  delta1 <- object$mlParam[(object$nXvar + 1):(object$nXvar + 
    object$nuZUvar)]
  phi1 <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar + 
    object$nuZUvar + object$nvZVvar)]
  beta2 <- object$mlParam[(object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar)]
  delta2 <- object$mlParam[(2 * object$nXvar + object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar)]
  phi2 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  theta <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(2 * object$nXvar + 2 * object$nuZUvar + 
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
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu1) * (1 - (dnorm(0)/pnorm(0))^2), 
    ncol = 1))
  colnames(margEff_c1) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu2) * (1 - (dnorm(0)/pnorm(0))^2), 
    ncol = 1))
  colnames(margEff_c2) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c2")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in 1:ncol(margEff_c1)) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], 
      margEff_c2[, c])
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c")
  margEff <- bind_cols(as_tibble(margEff_c), as_tibble(margEff_c1), 
    as_tibble(margEff_c2))
  return(margEff)
}
