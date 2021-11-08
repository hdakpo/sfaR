########################################################
#                                                      #
# Latent Class Model (LCM 4C) normal - half normal     #
#                                                      #
########################################################

# Log-likelihood ----------

cLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, 
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
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(S * epsilon3/sqrt(exp(Wu3) + 
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(S * epsilon4/sqrt(exp(Wu4) + 
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  L <- Probc1 * Pi1 + Probc2 * Pi2 + Probc3 * Pi3 + Probc4 * 
    Pi4
  ifelse(L <= 0, return(NA), return(log(L)))
}

# starting value for the log-likelihood ----------

csLCMfhalfnorm4C <- function(olsObj, epsiRes, nXvar, nuZUvar, 
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
    2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.98 * Esti[1:(nXvar)], 
    Esti[nXvar + 1], if (nuZUvar > 1) rep(0, nuZUvar - 1), 
    Esti[nXvar + 2], if (nvZVvar > 1) rep(0, nvZVvar - 1), 
    0.98 * Esti[1:(nXvar)], Esti[nXvar + 1], if (nuZUvar > 
      1) rep(0, nuZUvar - 1), Esti[nXvar + 2], if (nvZVvar > 
      1) rep(0, nvZVvar - 1), rep(0, 3 * nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar], 
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), 
    names(Esti)[1:nXvar], paste0("Zu_", colnames(uHvar)), 
    paste0("Zv_", colnames(vHvar)), names(Esti)[1:nXvar], 
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)), 
    paste0("Cl1_", colnames(Zvar)), paste0("Cl2_", colnames(Zvar)), 
    paste0("Cl3_", colnames(Zvar)))
  return(list(StartVal = StartVal, initHalf = initHalf))
}

# Gradient of the likelihood function ----------

cgradLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, 
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
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewv1 <- exp(Wv1)
  ewu2 <- exp(Wu2)
  ewv2 <- exp(Wv2)
  ewu3 <- exp(Wu3)
  ewv3 <- exp(Wv3)
  ewu4 <- exp(Wu4)
  ewv4 <- exp(Wv4)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigma_sq3 <- ewu3 + ewv3
  sigma_sq4 <- ewu4 + ewv4
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  sigmastar3 <- sqrt(ewu3 * ewv3/(sigma_sq3))
  sigmastar4 <- sqrt(ewu4 * ewv4/(sigma_sq4))
  dmusig1 <- dnorm(-(S * ewu1 * (epsilon1)/((sigma_sq1) * sigmastar1)))
  dmusig2 <- dnorm(-(S * ewu2 * (epsilon2)/((sigma_sq2) * sigmastar2)))
  dmusig3 <- dnorm(-(S * ewu3 * (epsilon3)/((sigma_sq3) * sigmastar3)))
  dmusig4 <- dnorm(-(S * ewu4 * (epsilon4)/((sigma_sq4) * sigmastar4)))
  pmusig1 <- pnorm(-(S * ewu1 * (epsilon1)/((sigma_sq1) * sigmastar1)))
  pmusig2 <- pnorm(-(S * ewu2 * (epsilon2)/((sigma_sq2) * sigmastar2)))
  pmusig3 <- pnorm(-(S * ewu3 * (epsilon3)/((sigma_sq3) * sigmastar3)))
  pmusig4 <- pnorm(-(S * ewu4 * (epsilon4)/((sigma_sq4) * sigmastar4)))
  depsisq1 <- dnorm(S * (epsilon1)/sqrt(sigma_sq1))
  depsisq2 <- dnorm(S * (epsilon2)/sqrt(sigma_sq2))
  depsisq3 <- dnorm(S * (epsilon3)/sqrt(sigma_sq3))
  depsisq4 <- dnorm(S * (epsilon4)/sqrt(sigma_sq4))
  sigx1_1 <- (dmusig1 * depsisq1 * ewu1/sigmastar1 + S * depsisq1 * 
    pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsisq2 * ewu2/sigmastar2 + S * depsisq2 * 
    pmusig2 * (epsilon2))
  sigx1_3 <- (dmusig3 * depsisq3 * ewu3/sigmastar3 + S * depsisq3 * 
    pmusig3 * (epsilon3))
  sigx1_4 <- (dmusig4 * depsisq4 * ewu4/sigmastar4 + S * depsisq4 * 
    pmusig4 * (epsilon4))
  sqsq1 <- ((sigma_sq1) * sigmastar1)
  sqsq2 <- ((sigma_sq2) * sigmastar2)
  sqsq3 <- ((sigma_sq3) * sigmastar3)
  sqsq4 <- ((sigma_sq4) * sigmastar4)
  sigx2_1 <- (0.5 * ((1 - ewu1/(sigma_sq1)) * ewv1/sigmastar1) + 
    sigmastar1)
  sigx2_2 <- (0.5 * ((1 - ewu2/(sigma_sq2)) * ewv2/sigmastar2) + 
    sigmastar2)
  sigx2_3 <- (0.5 * ((1 - ewu3/(sigma_sq3)) * ewv3/sigmastar3) + 
    sigmastar3)
  sigx2_4 <- (0.5 * ((1 - ewu4/(sigma_sq4)) * ewv4/sigmastar4) + 
    sigmastar4)
  sigx3_1 <- (0.5 * ((1 - ewv1/(sigma_sq1)) * ewu1/sigmastar1) + 
    sigmastar1)
  sigx3_2 <- (0.5 * ((1 - ewv2/(sigma_sq2)) * ewu2/sigmastar2) + 
    sigmastar2)
  sigx3_3 <- (0.5 * ((1 - ewv3/(sigma_sq3)) * ewu3/sigmastar3) + 
    sigmastar3)
  sigx3_4 <- (0.5 * ((1 - ewv4/(sigma_sq4)) * ewu4/sigmastar4) + 
    sigmastar4)
  wzdeno <- (1 + ewz1 + ewz2 + ewz3)
  prC <- (1 - (ewz1 + ewz2 + ewz3)/wzdeno)
  wzdsq1 <- wzdeno * sqrt(sigma_sq1)
  wzdsq2 <- wzdeno * sqrt(sigma_sq2)
  wzdsq3 <- wzdeno * sqrt(sigma_sq3)
  wzlogit1 <- (depsisq1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  wzlogit2 <- (depsisq2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  wzlogit3 <- (depsisq3 * ewz3 * pmusig3/sqrt(sigma_sq3))
  wzlogit4 <- (prC * depsisq4 * pmusig4/sqrt(sigma_sq4))
  sigx4 <- ((2 * wzlogit1 + 2 * wzlogit2 + 2 * wzlogit3)/wzdeno + 
    2 * wzlogit4)
  sigsq_1 <- (sigx4 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigsq_2 <- (sigx4 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigsq_3 <- (sigx4 * wzdeno * (sigma_sq3) * sqrt(sigma_sq3))
  sigsq_4 <- (sigx4 * (sigma_sq4) * sqrt(sigma_sq4))
  dpepsisq1 <- 0.5 * (S * depsisq1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  dpepsisq2 <- 0.5 * (S * depsisq2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  dpepsisq3 <- 0.5 * (S * depsisq3 * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  dpepsisq4 <- 0.5 * (S * depsisq4 * pmusig4 * (epsilon4)/(sigma_sq4)^2)
  wzdsig <- (sigx4 * wzdeno)
  dpsq1 <- (depsisq1 * pmusig1/sqrt(sigma_sq1))
  dpsq2 <- (depsisq2 * pmusig2/sqrt(sigma_sq2))
  dpsq3 <- (depsisq3 * pmusig3/sqrt(sigma_sq3))
  dpsq4 <- (depsisq4 * pmusig4/(sigma_sq4))
  sigwz1 <- sigx4 * wzdsq1
  sigwz2 <- sigx4 * wzdsq2
  sigwz3 <- sigx4 * wzdsq3
  sigwz4 <- sigx4 * sqrt(sigma_sq4)
  sigx5_1 <- 0.5 * (depsisq1 * pmusig1/(sigma_sq1))
  sigx5_2 <- 0.5 * (depsisq2 * pmusig2/(sigma_sq2))
  sigx5_3 <- 0.5 * (depsisq3 * pmusig3/(sigma_sq3))
  sqewu1 <- (1/sqsq1 - sigx2_1 * ewu1/sqsq1^2)
  sqewu2 <- (1/sqsq2 - sigx2_2 * ewu2/sqsq2^2)
  sqewu3 <- (1/sqsq3 - sigx2_3 * ewu3/sqsq3^2)
  sqewu4 <- (1/sqsq4 - sigx2_4 * ewu4/sqsq4^2)
  sigx6_1 <- (S * (dpepsisq1 - sqewu1 * dmusig1 * depsisq1) * 
    (epsilon1) - sigx5_1)
  sigx6_2 <- (S * (dpepsisq2 - sqewu2 * dmusig2 * depsisq2) * 
    (epsilon2) - sigx5_2)
  sigx6_3 <- (S * (dpepsisq3 - sqewu3 * dmusig3 * depsisq3) * 
    (epsilon3) - sigx5_3)
  sigx6_4 <- (S * (dpepsisq4 - sqewu4 * dmusig4 * depsisq4) * 
    (epsilon4) - 0.5 * dpsq4)
  sigx7_1 <- (S * (sigx3_1 * dmusig1 * depsisq1 * ewu1/sqsq1^2 + 
    dpepsisq1) * (epsilon1) - sigx5_1)
  sigx7_2 <- (S * (sigx3_2 * dmusig2 * depsisq2 * ewu2/sqsq2^2 + 
    dpepsisq2) * (epsilon2) - sigx5_2)
  sigx7_3 <- (S * (sigx3_3 * dmusig3 * depsisq3 * ewu3/sqsq3^2 + 
    dpepsisq3) * (epsilon3) - sigx5_3)
  sigx7_4 <- (S * (sigx3_4 * dmusig4 * depsisq4 * ewu4/sqsq4^2 + 
    dpepsisq4) * (epsilon4) - 0.5 * dpsq4)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = 2 * (S * 
    sigx1_1 * ewz1/sigsq_1), FUN = "*"), sweep(uHvar, MARGIN = 1, 
    STATS = 2 * (ewu1 * ewz1 * sigx6_1/(sigwz1)), FUN = "*"), 
    sweep(vHvar, MARGIN = 1, STATS = 2 * (ewv1 * ewz1 * sigx7_1/(sigwz1)), 
      FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = 2 * (S * 
      sigx1_2 * ewz2/sigsq_2), FUN = "*"), sweep(uHvar, 
      MARGIN = 1, STATS = 2 * (ewu2 * ewz2 * sigx6_2/(sigwz2)), 
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * 
      (ewv2 * ewz2 * sigx7_2/(sigwz2)), FUN = "*"), sweep(Xvar, 
      MARGIN = 1, STATS = 2 * (S * sigx1_3 * ewz3/sigsq_3), 
      FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = 2 * 
      (ewu3 * ewz3 * sigx6_3/(sigwz3)), FUN = "*"), sweep(vHvar, 
      MARGIN = 1, STATS = 2 * (ewv3 * ewz3 * sigx7_3/(sigwz3)), 
      FUN = "*"), sweep(Xvar, MARGIN = 1, STATS = 2 * (S * 
      prC * sigx1_4/sigsq_4), FUN = "*"), sweep(uHvar, 
      MARGIN = 1, STATS = 2 * (prC * ewu4 * sigx6_4/(sigwz4)), 
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = 2 * 
      (prC * ewv4 * sigx7_4/(sigwz4)), FUN = "*"), sweep(Zvar, 
      MARGIN = 1, STATS = (2 * dpsq1 - sigx4) * ewz1/wzdsig, 
      FUN = "*"), sweep(Zvar, MARGIN = 1, STATS = (2 * 
      dpsq2 - sigx4) * ewz2/wzdsig, FUN = "*"), sweep(Zvar, 
      MARGIN = 1, STATS = (2 * dpsq3 - sigx4) * ewz3/wzdsig, 
      FUN = "*"))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessLCMhalfnormlike4C <- function(parm, nXvar, nuZUvar, nvZVvar, 
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
  beta3 <- parm[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)]
  delta3 <- parm[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 2 * nvZVvar)]
  phi3 <- parm[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  beta4 <- parm[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)]
  delta4 <- parm[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 3 * nvZVvar)]
  phi4 <- parm[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar)]
  theta1 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)]
  theta2 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)]
  theta3 <- parm[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)]
  Wu1 <- as.numeric(crossprod(matrix(delta1), t(uHvar)))
  Wu2 <- as.numeric(crossprod(matrix(delta2), t(uHvar)))
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- Yvar - as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- Yvar - as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- Yvar - as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- Yvar - as.numeric(crossprod(matrix(beta4), t(Xvar)))
  ewu1 <- exp(Wu1)
  ewv1 <- exp(Wv1)
  ewu2 <- exp(Wu2)
  ewv2 <- exp(Wv2)
  ewu3 <- exp(Wu3)
  ewv3 <- exp(Wv3)
  ewu4 <- exp(Wu4)
  ewv4 <- exp(Wv4)
  ewz1 <- exp(Wz1)
  ewz2 <- exp(Wz2)
  ewz3 <- exp(Wz3)
  sigma_sq1 <- ewu1 + ewv1
  sigma_sq2 <- ewu2 + ewv2
  sigma_sq3 <- ewu3 + ewv3
  sigma_sq4 <- ewu4 + ewv4
  sigmastar1 <- sqrt(ewu1 * ewv1/(sigma_sq1))
  sigmastar2 <- sqrt(ewu2 * ewv2/(sigma_sq2))
  sigmastar3 <- sqrt(ewu3 * ewv3/(sigma_sq3))
  sigmastar4 <- sqrt(ewu4 * ewv4/(sigma_sq4))
  dmusig1 <- dnorm(-(S * ewu1 * (epsilon1)/((sigma_sq1) * sigmastar1)))
  dmusig2 <- dnorm(-(S * ewu2 * (epsilon2)/((sigma_sq2) * sigmastar2)))
  dmusig3 <- dnorm(-(S * ewu3 * (epsilon3)/((sigma_sq3) * sigmastar3)))
  dmusig4 <- dnorm(-(S * ewu4 * (epsilon4)/((sigma_sq4) * sigmastar4)))
  pmusig1 <- pnorm(-(S * ewu1 * (epsilon1)/((sigma_sq1) * sigmastar1)))
  pmusig2 <- pnorm(-(S * ewu2 * (epsilon2)/((sigma_sq2) * sigmastar2)))
  pmusig3 <- pnorm(-(S * ewu3 * (epsilon3)/((sigma_sq3) * sigmastar3)))
  pmusig4 <- pnorm(-(S * ewu4 * (epsilon4)/((sigma_sq4) * sigmastar4)))
  depsisq1 <- dnorm(S * (epsilon1)/sqrt(sigma_sq1))
  depsisq2 <- dnorm(S * (epsilon2)/sqrt(sigma_sq2))
  depsisq3 <- dnorm(S * (epsilon3)/sqrt(sigma_sq3))
  depsisq4 <- dnorm(S * (epsilon4)/sqrt(sigma_sq4))
  sigx1_1 <- (dmusig1 * depsisq1 * ewu1/sigmastar1 + S * depsisq1 * 
    pmusig1 * (epsilon1))
  sigx1_2 <- (dmusig2 * depsisq2 * ewu2/sigmastar2 + S * depsisq2 * 
    pmusig2 * (epsilon2))
  sigx1_3 <- (dmusig3 * depsisq3 * ewu3/sigmastar3 + S * depsisq3 * 
    pmusig3 * (epsilon3))
  sigx1_4 <- (dmusig4 * depsisq4 * ewu4/sigmastar4 + S * depsisq4 * 
    pmusig4 * (epsilon4))
  sqsq1 <- ((sigma_sq1) * sigmastar1)
  sqsq2 <- ((sigma_sq2) * sigmastar2)
  sqsq3 <- ((sigma_sq3) * sigmastar3)
  sqsq4 <- ((sigma_sq4) * sigmastar4)
  sigx2_1 <- (0.5 * ((1 - ewu1/(sigma_sq1)) * ewv1/sigmastar1) + 
    sigmastar1)
  sigx2_2 <- (0.5 * ((1 - ewu2/(sigma_sq2)) * ewv2/sigmastar2) + 
    sigmastar2)
  sigx2_3 <- (0.5 * ((1 - ewu3/(sigma_sq3)) * ewv3/sigmastar3) + 
    sigmastar3)
  sigx2_4 <- (0.5 * ((1 - ewu4/(sigma_sq4)) * ewv4/sigmastar4) + 
    sigmastar4)
  sigx3_1 <- (0.5 * ((1 - ewv1/(sigma_sq1)) * ewu1/sigmastar1) + 
    sigmastar1)
  sigx3_2 <- (0.5 * ((1 - ewv2/(sigma_sq2)) * ewu2/sigmastar2) + 
    sigmastar2)
  sigx3_3 <- (0.5 * ((1 - ewv3/(sigma_sq3)) * ewu3/sigmastar3) + 
    sigmastar3)
  sigx3_4 <- (0.5 * ((1 - ewv4/(sigma_sq4)) * ewu4/sigmastar4) + 
    sigmastar4)
  wzdeno <- (1 + ewz1 + ewz2 + ewz3)
  prC <- (1 - (ewz1 + ewz2 + ewz3)/wzdeno)
  wzdsq1 <- wzdeno * sqrt(sigma_sq1)
  wzdsq2 <- wzdeno * sqrt(sigma_sq2)
  wzdsq3 <- wzdeno * sqrt(sigma_sq3)
  wzlogit1 <- (depsisq1 * ewz1 * pmusig1/sqrt(sigma_sq1))
  wzlogit2 <- (depsisq2 * ewz2 * pmusig2/sqrt(sigma_sq2))
  wzlogit3 <- (depsisq3 * ewz3 * pmusig3/sqrt(sigma_sq3))
  wzlogit4 <- (prC * depsisq4 * pmusig4/sqrt(sigma_sq4))
  sigx4 <- ((2 * wzlogit1 + 2 * wzlogit2 + 2 * wzlogit3)/wzdeno + 
    2 * wzlogit4)
  sigsq_1 <- (sigx4 * wzdeno * (sigma_sq1) * sqrt(sigma_sq1))
  sigsq_2 <- (sigx4 * wzdeno * (sigma_sq2) * sqrt(sigma_sq2))
  sigsq_3 <- (sigx4 * wzdeno * (sigma_sq3) * sqrt(sigma_sq3))
  sigsq_4 <- (sigx4 * (sigma_sq4) * sqrt(sigma_sq4))
  dpepsisq1 <- 0.5 * (S * depsisq1 * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  dpepsisq2 <- 0.5 * (S * depsisq2 * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  dpepsisq3 <- 0.5 * (S * depsisq3 * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  dpepsisq4 <- 0.5 * (S * depsisq4 * pmusig4 * (epsilon4)/(sigma_sq4)^2)
  wzdsig <- (sigx4 * wzdeno)
  dpsq1 <- (depsisq1 * pmusig1/sqrt(sigma_sq1))
  dpsq2 <- (depsisq2 * pmusig2/sqrt(sigma_sq2))
  dpsq3 <- (depsisq3 * pmusig3/sqrt(sigma_sq3))
  dpsq4 <- (depsisq4 * pmusig4/(sigma_sq4))
  sigwz1 <- sigx4 * wzdsq1
  sigwz2 <- sigx4 * wzdsq2
  sigwz3 <- sigx4 * wzdsq3
  sigwz4 <- sigx4 * sqrt(sigma_sq4)
  sigx5_1 <- 0.5 * (depsisq1 * pmusig1/(sigma_sq1))
  sigx5_2 <- 0.5 * (depsisq2 * pmusig2/(sigma_sq2))
  sigx5_3 <- 0.5 * (depsisq3 * pmusig3/(sigma_sq3))
  sqewu1 <- (1/sqsq1 - sigx2_1 * ewu1/sqsq1^2)
  sqewu2 <- (1/sqsq2 - sigx2_2 * ewu2/sqsq2^2)
  sqewu3 <- (1/sqsq3 - sigx2_3 * ewu3/sqsq3^2)
  sqewu4 <- (1/sqsq4 - sigx2_4 * ewu4/sqsq4^2)
  sigx6_1 <- (S * (dpepsisq1 - sqewu1 * dmusig1 * depsisq1) * 
    (epsilon1) - sigx5_1)
  sigx6_2 <- (S * (dpepsisq2 - sqewu2 * dmusig2 * depsisq2) * 
    (epsilon2) - sigx5_2)
  sigx6_3 <- (S * (dpepsisq3 - sqewu3 * dmusig3 * depsisq3) * 
    (epsilon3) - sigx5_3)
  sigx6_4 <- (S * (dpepsisq4 - sqewu4 * dmusig4 * depsisq4) * 
    (epsilon4) - 0.5 * dpsq4)
  sigx7_1 <- (S * (sigx3_1 * dmusig1 * depsisq1 * ewu1/sqsq1^2 + 
    dpepsisq1) * (epsilon1) - sigx5_1)
  sigx7_2 <- (S * (sigx3_2 * dmusig2 * depsisq2 * ewu2/sqsq2^2 + 
    dpepsisq2) * (epsilon2) - sigx5_2)
  sigx7_3 <- (S * (sigx3_3 * dmusig3 * depsisq3 * ewu3/sqsq3^2 + 
    dpepsisq3) * (epsilon3) - sigx5_3)
  sigx7_4 <- (S * (sigx3_4 * dmusig4 * depsisq4 * ewu4/sqsq4^2 + 
    dpepsisq4) * (epsilon4) - 0.5 * dpsq4)
  sigx4wzsq <- sigx4 * wzdeno^2
  sigx8_1 <- (sigx4 * wzdeno * (sigma_sq1)/sqrt(sigma_sq1))
  sigx8_2 <- (sigx4 * wzdeno * (sigma_sq2)/sqrt(sigma_sq2))
  sigx8_3 <- (sigx4 * wzdeno * (sigma_sq3)/sqrt(sigma_sq3))
  sigx9_1 <- (sigx4 * wzdeno/sqrt(sigma_sq1))
  sigx9_2 <- (sigx4 * wzdeno/sqrt(sigma_sq2))
  sigx9_3 <- (sigx4 * wzdeno/sqrt(sigma_sq3))
  wusq1 <- 1 - ewu1/(sigma_sq1)
  wvsq1 <- 1 - ewv1/(sigma_sq1)
  wusq2 <- 1 - ewu2/(sigma_sq2)
  wvsq2 <- 1 - ewv2/(sigma_sq2)
  wusq3 <- 1 - ewu3/(sigma_sq3)
  wvsq3 <- 1 - ewv3/(sigma_sq3)
  wusq4 <- 1 - ewu4/(sigma_sq4)
  wvsq4 <- 1 - ewv4/(sigma_sq4)
  dwp1 <- depsisq1 * ewz1 * pmusig1
  dwp2 <- depsisq2 * ewz2 * pmusig2
  dwp3 <- depsisq3 * ewz3 * pmusig3
  wzlx3 <- (2 * wzlogit1 + 2 * wzlogit2 + 2 * wzlogit3)/wzdeno
  duv1 <- depsisq1 * ewu1/ewv1 + depsisq1
  duv2 <- depsisq2 * ewu2/ewv2 + depsisq2
  duv3 <- depsisq3 * ewu3/ewv3 + depsisq3
  duv4 <- depsisq4 * ewu4/ewv4 + depsisq4
  pepsisq1 <- (S * pmusig1 * (epsilon1)/(sigma_sq1)^2)
  pepsisq2 <- (S * pmusig2 * (epsilon2)/(sigma_sq2)^2)
  pepsisq3 <- (S * pmusig3 * (epsilon3)/(sigma_sq3)^2)
  pepsisq4 <- (S * pmusig4 * (epsilon4)/(sigma_sq4)^2)
  dep1sq1 <- depsisq1/(sigma_sq1)
  dep2sq2 <- depsisq2/(sigma_sq2)
  dep3sq3 <- depsisq3/(sigma_sq3)
  dep4sq4 <- depsisq4/(sigma_sq4)
  sigx10_1 <- (0.5 * (S * (dmusig1 * ewu1/sigmastar1 + S * 
    pmusig1 * (epsilon1)) * (epsilon1)/(sigma_sq1) - pmusig1) - 
    0.5 * pmusig1) * dep1sq1
  sigx10_2 <- (0.5 * (S * (dmusig2 * ewu2/sigmastar2 + S * 
    pmusig2 * (epsilon2)) * (epsilon2)/(sigma_sq2) - pmusig2) - 
    0.5 * pmusig2) * dep2sq2
  sigx10_3 <- (0.5 * (S * (dmusig3 * ewu3/sigmastar3 + S * 
    pmusig3 * (epsilon3)) * (epsilon3)/(sigma_sq3) - pmusig3) - 
    0.5 * pmusig3) * dep3sq3
  sigx10_4 <- (0.5 * (S * (dmusig4 * ewu4/sigmastar4 + S * 
    pmusig4 * (epsilon4)) * (epsilon4)/(sigma_sq4) - pmusig4) - 
    0.5 * pmusig4) * dep4sq4
  sigx11_1 <- (dpepsisq1 - sqewu1 * dmusig1 * depsisq1)
  sigx11_2 <- (dpepsisq2 - sqewu2 * dmusig2 * depsisq2)
  sigx11_3 <- (dpepsisq3 - sqewu3 * dmusig3 * depsisq3)
  sigx11_4 <- (dpepsisq4 - sqewu4 * dmusig4 * depsisq4)
  sigx12_1 <- ((S * sigx11_1 * (epsilon1) - depsisq1 * pmusig1/(sigma_sq1))/(sigma_sq1))
  sigx12_2 <- ((S * sigx11_2 * (epsilon2) - depsisq2 * pmusig2/(sigma_sq2))/(sigma_sq2))
  sigx12_3 <- ((S * sigx11_3 * (epsilon3) - depsisq3 * pmusig3/(sigma_sq3))/(sigma_sq3))
  sigx12_4 <- ((S * sigx11_4 * (epsilon4) - depsisq4 * pmusig4/(sigma_sq4))/(sigma_sq4))
  sigx13_1 <- (S * (sigx3_1 * dmusig1 * depsisq1 * ewu1/sqsq1^2 + 
    dpepsisq1) * (epsilon1) - depsisq1 * pmusig1/(sigma_sq1))
  sigx13_2 <- (S * (sigx3_2 * dmusig2 * depsisq2 * ewu2/sqsq2^2 + 
    dpepsisq2) * (epsilon2) - depsisq2 * pmusig2/(sigma_sq2))
  sigx13_3 <- (S * (sigx3_3 * dmusig3 * depsisq3 * ewu3/sqsq3^2 + 
    dpepsisq3) * (epsilon3) - depsisq3 * pmusig3/(sigma_sq3))
  sigx13_4 <- (S * (sigx3_4 * dmusig4 * depsisq4 * ewu4/sqsq4^2 + 
    dpepsisq4) * (epsilon4) - depsisq4 * pmusig4/(sigma_sq4))
  wzsigx4_1 <- ((sigwz4)^2 * wzdeno * (sigma_sq1)^(3/2))
  wzsigx4_2 <- ((sigwz4)^2 * wzdeno * (sigma_sq2)^(3/2))
  wzsigx4_3 <- ((sigwz4)^2 * wzdeno * (sigma_sq3)^(3/2))
  sigx4wz <- (2 * ((2 * dpsq2 - sigx4)/wzdsig^2) + 2/(sigx4wzsq))
  sigx4w2z <- ((2 - 2 * (ewz1/wzdeno))/wzdsig - 2 * ((2 * dpsq1 - 
    sigx4) * ewz1/wzdsig^2))
  sigx4w3z <- (2 * ((2 * dpsq1 - sigx4)/wzdsig^2) + 2/(sigx4wzsq))
  sigx4w4z <- ((2 - 2 * (ewz2/wzdeno))/wzdsig - 2 * ((2 * dpsq2 - 
    sigx4) * ewz2/wzdsig^2))
  sigx4w5z <- (2 * (wzdeno * (2 * dpsq1 - sigx4)/wzdsig^2) + 
    2/wzdsig)
  sigx4w6z <- (2 * (wzdeno * (2 * dpsq2 - sigx4)/wzdsig^2) + 
    2/wzdsig)
  sigx4w7z <- ((2 * dpsq1 - sigx4) * sqrt(sigma_sq4)/(sigwz4)^2 + 
    1/(sigwz4))
  sigx4w8z <- ((2 * dpsq2 - sigx4) * sqrt(sigma_sq4)/(sigwz4)^2 + 
    1/(sigwz4))
  sigx4w9z <- (0.5 * (sigx4/sqrt(sigma_sq4)) + 2 * (prC * sigx6_4))
  sigx4w10z <- (sigx4 * (sigma_sq4)/sqrt(sigma_sq4))
  sigx4w11z <- (0.5 * (sigx4/sqrt(sigma_sq4)) + 2 * (prC * 
    sigx7_4))
  sigx4w12z <- (2 * ((2 * dpsq3 - sigx4)/wzdsig^2) + 2/(sigx4wzsq))
  sigx4w13z <- ((2 - 2 * (ewz3/wzdeno))/wzdsig - 2 * ((2 * 
    dpsq3 - sigx4) * ewz3/wzdsig^2))
  sigx4w14z <- (2 * (wzdeno * (2 * dpsq3 - sigx4)/wzdsig^2) + 
    2/wzdsig)
  sigx14_1 <- (0.5 * pepsisq1 - sqewu1 * dmusig1)
  sigx14_2 <- (0.5 * pepsisq2 - sqewu2 * dmusig2)
  sigx14_3 <- (0.5 * pepsisq3 - sqewu3 * dmusig3)
  sigx14_4 <- (0.5 * pepsisq4 - sqewu4 * dmusig4)
  sigx15_1 <- (S * depsisq1 * (S * sigx14_1 * (epsilon1) - 
    2 * (pmusig1/(sigma_sq1))) * (epsilon1)/(sigma_sq1)^2)
  sigx15_2 <- (S * depsisq2 * (S * sigx14_2 * (epsilon2) - 
    2 * (pmusig2/(sigma_sq2))) * (epsilon2)/(sigma_sq2)^2)
  sigx15_3 <- (S * depsisq3 * (S * sigx14_3 * (epsilon3) - 
    2 * (pmusig3/(sigma_sq3))) * (epsilon3)/(sigma_sq3)^2)
  sigx15_4 <- (S * depsisq4 * (S * sigx14_4 * (epsilon4) - 
    2 * (pmusig4/(sigma_sq4))) * (epsilon4)/(sigma_sq4)^2)
  sigx16_1 <- (sigx3_1 * dmusig1 * ewu1/sqsq1^2 + 0.5 * pepsisq1)
  sigx16_2 <- (sigx3_2 * dmusig2 * ewu2/sqsq2^2 + 0.5 * pepsisq2)
  sigx16_3 <- (sigx3_3 * dmusig3 * ewu3/sqsq3^2 + 0.5 * pepsisq3)
  sigx16_4 <- (sigx3_4 * dmusig4 * ewu4/sqsq4^2 + 0.5 * pepsisq4)
  s3xq1 <- (sigma_sq1) * sigmastar1/sqsq1^2
  s3xq2 <- (sigma_sq2) * sigmastar2/sqsq2^2
  s3xq3 <- (sigma_sq3) * sigmastar3/sqsq3^2
  s3xq4 <- (sigma_sq4) * sigmastar4/sqsq4^2
  sigx17_1 <- (0.5 * sigx9_1 + 2 * (ewz1 * sigx6_1))
  sigx17_2 <- (0.5 * sigx9_2 + 2 * (ewz2 * sigx6_2))
  sigx17_3 <- (0.5 * sigx9_3 + 2 * (ewz3 * sigx6_3))
  sigx18_1 <- (0.5 * sigx9_1 + 2 * (ewz1 * sigx7_1))
  sigx18_2 <- (0.5 * sigx9_2 + 2 * (ewz2 * sigx7_2))
  sigx18_3 <- (0.5 * sigx9_3 + 2 * (ewz3 * sigx7_3))
  sigx19_1 <- (S * sigx16_1 * (epsilon1) - 2 * (pmusig1/(sigma_sq1)))
  sigx19_2 <- (S * sigx16_2 * (epsilon2) - 2 * (pmusig2/(sigma_sq2)))
  sigx19_3 <- (S * sigx16_3 * (epsilon3) - 2 * (pmusig3/(sigma_sq3)))
  sigx19_4 <- (S * sigx16_4 * (epsilon4) - 2 * (pmusig4/(sigma_sq4)))
  hessll <- matrix(nrow = (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar), ncol = (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar))
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S^2 * ((depsisq1 * (S * (dmusig1 * ewu1/sigmastar1 + 
      S * pmusig1 * (epsilon1)) * (epsilon1)/(sigma_sq1) - 
      pmusig1) + S * dmusig1 * (duv1) * ewu1 * (epsilon1)/sqsq1)/sigsq_1 - 
      2 * (sigx1_1^2 * ewz1/sigsq_1^2)) * ewz1), FUN = "*"), 
    Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * ((sqewu1 * dmusig1 * depsisq1 + 
      (S * (sigx10_1 - S * sqewu1 * dmusig1 * (duv1) * 
        (epsilon1)) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1))/(sigwz1) - 
      2 * (sigx1_1 * ewz1 * sigx6_1/((sigwz1)^2 * (sigma_sq1)))) * 
      ewu1 * ewz1), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + 
    nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = 2 * 
    (S * (((S * (sigx10_1 + S * sigx3_1 * dmusig1 * (duv1) * 
      ewu1 * (epsilon1)/sqsq1^2) * (epsilon1) - 0.5 * (sigx1_1/(sigma_sq1)))/(sigma_sq1) - 
      sigx3_1 * dmusig1 * depsisq1 * ewu1/sqsq1^2)/(sigwz1) - 
      2 * (sigx1_1 * ewz1 * sigx7_1/((sigwz1)^2 * (sigma_sq1)))) * 
      ewv1 * ewz1), FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = -(4 * (S^2 * sigx1_1 * sigx1_2 * (sigma_sq2) * 
      ewz1 * ewz2 * sqrt(sigma_sq2)/(sigsq_2^2 * (sigma_sq1)^(3/2)))), 
    FUN = "*"), Xvar)
  hessll[1:nXvar, (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_1 * ewu2 * ewz1 * 
      ewz2 * sigx6_2 * sqrt(sigma_sq2)/((sigwz2)^2 * (sigma_sq1)^(3/2)))), 
    FUN = "*"), uHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_1 * ewv2 * ewz1 * 
      ewz2 * sigx7_2 * sqrt(sigma_sq2)/((sigwz2)^2 * (sigma_sq1)^(3/2)))), 
    FUN = "*"), vHvar)
  hessll[1:nXvar, (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S^2 * sigx1_1 * sigx1_3 * 
      (sigma_sq3) * ewz1 * ewz3 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      (sigma_sq1)^(3/2)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_1 * ewu3 * ewz1 * 
      ewz3 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * (sigma_sq1)^(3/2)))), 
    FUN = "*"), uHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_1 * ewv3 * ewz1 * 
      ewz3 * sigx7_3 * sqrt(sigma_sq3)/((sigwz3)^2 * (sigma_sq1)^(3/2)))), 
    FUN = "*"), vHvar)
  hessll[1:nXvar, (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S^2 * prC * sigx1_1 * sigx1_4 * 
      (sigma_sq4) * ewz1 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdeno * (sigma_sq1)^(3/2)))), FUN = "*"), Xvar)
  hessll[1:nXvar, (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_1 * ewu4 * 
      ewz1 * sigx6_4 * sqrt(sigma_sq4)/wzsigx4_1)), FUN = "*"), 
    uHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_1 * ewv4 * 
      ewz1 * sigx7_4 * sqrt(sigma_sq4)/wzsigx4_1)), FUN = "*"), 
    vHvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = S * sigx4w2z * sigx1_1 * ewz1/((sigma_sq1)^(3/2)), 
    FUN = "*"), Zvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -(S * 
    sigx4wz * sigx1_1 * ewz1 * ewz2/((sigma_sq1)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[1:nXvar, (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -(S * 
    sigx4w12z * sigx1_1 * ewz1 * ewz3/((sigma_sq1)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar + 
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = 2 * 
    (((ewu1 * (S * (0.5 * sigx15_1 - (0.5 * (S^2 * sqewu1 * 
      depsisq1 * (epsilon1)^2/(sigma_sq1)^2) - (((0.5 * 
      (ewu1/(sigma_sq1)) + 1 - 0.5 * (0.5 * (wusq1) + ewu1/(sigma_sq1))) * 
      (wusq1) * ewv1/sigmastar1 + (2 - 2 * (sigx2_1^2 * 
      ewu1 * (sigma_sq1)/sqsq1^2)) * sigmastar1)/sqsq1^2 + 
      S^2 * sqewu1^2 * ewu1 * (epsilon1)^2/sqsq1) * depsisq1) * 
      dmusig1) * (epsilon1) - 0.5 * sigx12_1) + S * sigx11_1 * 
      (epsilon1) - sigx5_1)/(sigwz1) - sigx17_1 * ewu1 * 
      sigx6_1/(sigwz1)^2) * ewu1 * ewz1), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((wusq1) * 
      ewv1) - S^2 * sigx3_1 * sqewu1 * ewu1 * (epsilon1)^2)/(sigma_sq1) + 
      0.5 * ((ewu1/(sigma_sq1) - 1) * ewv1/(sigma_sq1) + 
        1 - 0.5 * ((wusq1) * (wvsq1)))) * depsisq1/sigmastar1 + 
      0.5 * (S^2 * sigx3_1 * depsisq1 * (epsilon1)^2/(sigma_sq1)^2)) * 
      ewu1 + sigx3_1 * (1 - 2 * (sigx2_1 * ewu1 * s3xq1)) * 
      depsisq1) * dmusig1/sqsq1^2 + 0.5 * sigx15_1) * (epsilon1) - 
      0.5 * sigx12_1)/(sigwz1) - sigx17_1 * sigx7_1/(sigwz1)^2) * 
      ewu1 * ewv1 * ewz1), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 
    nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_2 * ewu1 * (sigma_sq2) * 
      ewz1 * ewz2 * sigx6_1 * sqrt(sigma_sq2)/(sigsq_2^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu1 * ewu2 * ewz1 * ewz2 * 
      sigx6_1 * sigx6_2 * sqrt(sigma_sq2)/((sigwz2)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu1 * ewv2 * ewz1 * ewz2 * 
      sigx7_2 * sigx6_1 * sqrt(sigma_sq2)/((sigwz2)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_3 * ewu1 * (sigma_sq3) * 
      ewz1 * ewz3 * sigx6_1 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu1 * ewu3 * ewz1 * ewz3 * 
      sigx6_1 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu1 * ewv3 * ewz1 * ewz3 * 
      sigx7_3 * sigx6_1 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * ewu1 * 
      (sigma_sq4) * ewz1 * sigx6_1 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq1))), FUN = "*"), Xvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu1 * ewu4 * ewz1 * 
      sigx6_1 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq1))), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu1 * ewv4 * ewz1 * 
      sigx7_4 * sigx6_1 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq1))), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = sigx4w2z * 
    ewu1 * ewz1 * sigx6_1/sqrt(sigma_sq1), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(sigx4wz * ewu1 * ewz1 * ewz2 * 
      sigx6_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(sigx4w12z * ewu1 * ewz1 * ewz3 * 
      sigx6_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (ewv1/(sigma_sq1)) - 
      0.5 * (0.5 * (wvsq1) + ewv1/(sigma_sq1))) * (wvsq1) + 
      S^2 * sigx3_1^2 * ewu1 * ewv1 * (epsilon1)^2/(sqsq1^2 * 
        (sigma_sq1))) * depsisq1 * ewu1/sigmastar1 + 
      ((0.5 * (S^2 * depsisq1 * (epsilon1)^2/(sigma_sq1)^2) - 
        2 * (sigx3_1 * depsisq1 * s3xq1)) * ewv1 + depsisq1) * 
        sigx3_1) * dmusig1 * ewu1/sqsq1^2 + S * (0.5 * 
      (ewv1 * sigx19_1) + 0.5 * pmusig1) * depsisq1 * (epsilon1)/(sigma_sq1)^2) * 
      (epsilon1) - (0.5 * (depsisq1 * pmusig1) + 0.5 * 
      (ewv1 * sigx13_1))/(sigma_sq1))/(sigwz1) - sigx18_1 * 
      ewv1 * sigx7_1/(sigwz1)^2) * ewv1 * ewz1), FUN = "*"), 
    vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
      nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(4 * 
    (S * sigx1_2 * (sigma_sq2) * ewv1 * ewz1 * ewz2 * sigx7_1 * 
      sqrt(sigma_sq2)/(sigsq_2^2 * sqrt(sigma_sq1)))), 
    FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
      nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar, MARGIN = 1, 
    STATS = -(4 * (ewu2 * ewv1 * ewz1 * ewz2 * sigx7_1 * 
      sigx6_2 * sqrt(sigma_sq2)/((sigwz2)^2 * sqrt(sigma_sq1)))), 
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (ewv1 * ewv2 * ewz1 * ewz2 * 
      sigx7_1 * sigx7_2 * sqrt(sigma_sq2)/((sigwz2)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_3 * (sigma_sq3) * 
      ewv1 * ewz1 * ewz3 * sigx7_1 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
      3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (ewu3 * ewv1 * ewz1 * ewz3 * 
      sigx7_1 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
      3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (ewv1 * ewv3 * ewz1 * ewz3 * 
      sigx7_1 * sigx7_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq1)))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
      3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * (sigma_sq4) * 
      ewv1 * ewz1 * sigx7_1 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq1))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
      4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu4 * ewv1 * ewz1 * 
      sigx7_1 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq1))), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
      4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewv1 * ewv4 * ewz1 * 
      sigx7_1 * sigx7_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq1))), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 
      4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = sigx4w2z * ewv1 * ewz1 * sigx7_1/sqrt(sigma_sq1), 
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * 
      nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(sigx4wz * ewv1 * ewz1 * ewz2 * 
      sigx7_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 
      1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(sigx4w12z * ewv1 * ewz1 * ewz3 * 
      sigx7_1/sqrt(sigma_sq1)), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S^2 * ((depsisq2 * (S * (dmusig2 * ewu2/sigmastar2 + 
      S * pmusig2 * (epsilon2)) * (epsilon2)/(sigma_sq2) - 
      pmusig2) + S * dmusig2 * (duv2) * ewu2 * (epsilon2)/sqsq2)/sigsq_2 - 
      2 * (sigx1_2^2 * ewz2/sigsq_2^2)) * ewz2), FUN = "*"), 
    Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = 2 * (S * ((sqewu2 * dmusig2 * depsisq2 + (S * 
      (sigx10_2 - S * sqewu2 * dmusig2 * (duv2) * (epsilon2)) * 
      (epsilon2) - 0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2))/(sigwz2) - 
      2 * (sigx1_2 * ewz2 * sigx6_2/((sigwz2)^2 * (sigma_sq2)))) * 
      ewu2 * ewz2), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * (((S * (sigx10_2 + S * sigx3_2 * 
      dmusig2 * (duv2) * ewu2 * (epsilon2)/sqsq2^2) * (epsilon2) - 
      0.5 * (sigx1_2/(sigma_sq2)))/(sigma_sq2) - sigx3_2 * 
      dmusig2 * depsisq2 * ewu2/sqsq2^2)/(sigwz2) - 2 * 
      (sigx1_2 * ewz2 * sigx7_2/((sigwz2)^2 * (sigma_sq2)))) * 
      ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S^2 * sigx1_2 * sigx1_3 * 
      (sigma_sq3) * ewz2 * ewz3 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      (sigma_sq2)^(3/2)))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_2 * ewu3 * ewz2 * 
      ewz3 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * (sigma_sq2)^(3/2)))), 
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_2 * ewv3 * ewz2 * 
      ewz3 * sigx7_3 * sqrt(sigma_sq3)/((sigwz3)^2 * (sigma_sq2)^(3/2)))), 
    FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S^2 * prC * sigx1_2 * sigx1_4 * 
      (sigma_sq4) * ewz2 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdeno * (sigma_sq2)^(3/2)))), FUN = "*"), Xvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_2 * ewu4 * 
      ewz2 * sigx6_4 * sqrt(sigma_sq4)/wzsigx4_2)), FUN = "*"), 
    uHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_2 * ewv4 * 
      ewz2 * sigx7_4 * sqrt(sigma_sq4)/wzsigx4_2)), FUN = "*"), 
    vHvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * 
    nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(S * sigx4w3z * sigx1_2 * ewz1 * 
      ewz2/((sigma_sq2)^(3/2))), FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = S * sigx4w4z * sigx1_2 * ewz2/((sigma_sq2)^(3/2)), 
    FUN = "*"), Zvar)
  hessll[(nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + 
    nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -(S * 
    sigx4w12z * sigx1_2 * ewz2 * ewz3/((sigma_sq2)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + nuZUvar + nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((ewu2 * (S * (0.5 * sigx15_2 - 
      (0.5 * (S^2 * sqewu2 * depsisq2 * (epsilon2)^2/(sigma_sq2)^2) - 
        (((0.5 * (ewu2/(sigma_sq2)) + 1 - 0.5 * (0.5 * 
          (wusq2) + ewu2/(sigma_sq2))) * (wusq2) * ewv2/sigmastar2 + 
          (2 - 2 * (sigx2_2^2 * ewu2 * (sigma_sq2)/sqsq2^2)) * 
          sigmastar2)/sqsq2^2 + S^2 * sqewu2^2 * ewu2 * 
          (epsilon2)^2/sqsq2) * depsisq2) * dmusig2) * 
      (epsilon2) - 0.5 * sigx12_2) + S * sigx11_2 * (epsilon2) - 
      sigx5_2)/(sigwz2) - sigx17_2 * ewu2 * sigx6_2/(sigwz2)^2) * 
      ewu2 * ewz2), FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + nvZVvar + 
    1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((wusq2) * 
      ewv2) - S^2 * sigx3_2 * sqewu2 * ewu2 * (epsilon2)^2)/(sigma_sq2) + 
      0.5 * ((ewu2/(sigma_sq2) - 1) * ewv2/(sigma_sq2) + 
        1 - 0.5 * ((wusq2) * (wvsq2)))) * depsisq2/sigmastar2 + 
      0.5 * (S^2 * sigx3_2 * depsisq2 * (epsilon2)^2/(sigma_sq2)^2)) * 
      ewu2 + sigx3_2 * (1 - 2 * (sigx2_2 * ewu2 * s3xq2)) * 
      depsisq2) * dmusig2/sqsq2^2 + 0.5 * sigx15_2) * (epsilon2) - 
      0.5 * sigx12_2)/(sigwz2) - sigx17_2 * sigx7_2/(sigwz2)^2) * 
      ewu2 * ewv2 * ewz2), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_3 * ewu2 * (sigma_sq3) * 
      ewz2 * ewz3 * sigx6_2 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu2 * ewu3 * ewz2 * ewz3 * 
      sigx6_2 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 
    1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (ewu2 * ewv3 * ewz2 * ewz3 * 
      sigx7_3 * sigx6_2 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * ewu2 * 
      (sigma_sq4) * ewz2 * sigx6_2 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq2))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu2 * ewu4 * ewz2 * 
      sigx6_2 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq2))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu2 * ewv4 * ewz2 * 
      sigx7_4 * sigx6_2 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq2))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(sigx4w3z * ewu2 * ewz1 * ewz2 * 
      sigx6_2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    2 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = sigx4w4z * 
    ewu2 * ewz2 * sigx6_2/sqrt(sigma_sq2), FUN = "*"), Zvar)
  hessll[(2 * nXvar + nuZUvar + nvZVvar + 1):(2 * nXvar + 2 * 
    nuZUvar + nvZVvar), (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    3 * nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -(sigx4w12z * 
    ewu2 * ewz2 * ewz3 * sigx6_2/sqrt(sigma_sq2)), FUN = "*"), 
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 
    nvZVvar + 1):(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (ewv2/(sigma_sq2)) - 
      0.5 * (0.5 * (wvsq2) + ewv2/(sigma_sq2))) * (wvsq2) + 
      S^2 * sigx3_2^2 * ewu2 * ewv2 * (epsilon2)^2/(sqsq2^2 * 
        (sigma_sq2))) * depsisq2 * ewu2/sigmastar2 + 
      ((0.5 * (S^2 * depsisq2 * (epsilon2)^2/(sigma_sq2)^2) - 
        2 * (sigx3_2 * depsisq2 * s3xq2)) * ewv2 + depsisq2) * 
        sigx3_2) * dmusig2 * ewu2/sqsq2^2 + S * (0.5 * 
      (ewv2 * sigx19_2) + 0.5 * pmusig2) * depsisq2 * (epsilon2)/(sigma_sq2)^2) * 
      (epsilon2) - (0.5 * (depsisq2 * pmusig2) + 0.5 * 
      (ewv2 * sigx13_2))/(sigma_sq2))/(sigwz2) - sigx18_2 * 
      ewv2 * sigx7_2/(sigwz2)^2) * ewv2 * ewz2), FUN = "*"), 
    vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (S * sigx1_3 * (sigma_sq3) * 
      ewv2 * ewz2 * ewz3 * sigx7_2 * sqrt(sigma_sq3)/(sigsq_3^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (ewu3 * ewv2 * ewz2 * ewz3 * 
      sigx7_2 * sigx6_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (ewv2 * ewv3 * ewz2 * ewz3 * 
      sigx7_2 * sigx7_3 * sqrt(sigma_sq3)/((sigwz3)^2 * 
      sqrt(sigma_sq2)))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * (sigma_sq4) * 
      ewv2 * ewz2 * sigx7_2 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq2))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu4 * ewv2 * ewz2 * 
      sigx7_2 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq2))), FUN = "*"), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewv2 * ewv4 * ewz2 * 
      sigx7_2 * sigx7_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq2))), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(sigx4w3z * 
    ewv2 * ewz1 * ewz2 * sigx7_2/sqrt(sigma_sq2)), FUN = "*"), 
    Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = sigx4w4z * ewv2 * ewz2 * sigx7_2/sqrt(sigma_sq2), 
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + nvZVvar + 1):(2 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(sigx4w12z * ewv2 * ewz2 * ewz3 * 
      sigx7_2/sqrt(sigma_sq2)), FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (2 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S^2 * ((depsisq3 * (S * (dmusig3 * 
      ewu3/sigmastar3 + S * pmusig3 * (epsilon3)) * (epsilon3)/(sigma_sq3) - 
      pmusig3) + S * dmusig3 * (duv3) * ewu3 * (epsilon3)/sqsq3)/sigsq_3 - 
      2 * (sigx1_3^2 * ewz3/sigsq_3^2)) * ewz3), FUN = "*"), 
    Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * ((sqewu3 * dmusig3 * depsisq3 + 
      (S * (sigx10_3 - S * sqewu3 * dmusig3 * (duv3) * 
        (epsilon3)) * (epsilon3) - 0.5 * (sigx1_3/(sigma_sq3)))/(sigma_sq3))/(sigwz3) - 
      2 * (sigx1_3 * ewz3 * sigx6_3/((sigwz3)^2 * (sigma_sq3)))) * 
      ewu3 * ewz3), FUN = "*"), uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * (((S * (sigx10_3 + S * sigx3_3 * 
      dmusig3 * (duv3) * ewu3 * (epsilon3)/sqsq3^2) * (epsilon3) - 
      0.5 * (sigx1_3/(sigma_sq3)))/(sigma_sq3) - sigx3_3 * 
      dmusig3 * depsisq3 * ewu3/sqsq3^2)/(sigwz3) - 2 * 
      (sigx1_3 * ewz3 * sigx7_3/((sigwz3)^2 * (sigma_sq3)))) * 
      ewv3 * ewz3), FUN = "*"), vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S^2 * prC * sigx1_3 * sigx1_4 * 
      (sigma_sq4) * ewz3 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdeno * (sigma_sq3)^(3/2)))), FUN = "*"), Xvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_3 * ewu4 * 
      ewz3 * sigx6_4 * sqrt(sigma_sq4)/wzsigx4_3)), FUN = "*"), 
    uHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_3 * ewv4 * 
      ewz3 * sigx7_4 * sqrt(sigma_sq4)/wzsigx4_3)), FUN = "*"), 
    vHvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -(S * 
    sigx4w3z * sigx1_3 * ewz1 * ewz3/((sigma_sq3)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = -(S * sigx4wz * sigx1_3 * ewz2 * ewz3/((sigma_sq3)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(2 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    2 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = S * sigx4w13z * sigx1_3 * ewz3/((sigma_sq3)^(3/2)), 
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 2 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((ewu3 * (S * (0.5 * sigx15_3 - 
      (0.5 * (S^2 * sqewu3 * depsisq3 * (epsilon3)^2/(sigma_sq3)^2) - 
        (((0.5 * (ewu3/(sigma_sq3)) + 1 - 0.5 * (0.5 * 
          (wusq3) + ewu3/(sigma_sq3))) * (wusq3) * ewv3/sigmastar3 + 
          (2 - 2 * (sigx2_3^2 * ewu3 * (sigma_sq3)/sqsq3^2)) * 
          sigmastar3)/sqsq3^2 + S^2 * sqewu3^2 * ewu3 * 
          (epsilon3)^2/sqsq3) * depsisq3) * dmusig3) * 
      (epsilon3) - 0.5 * sigx12_3) + S * sigx11_3 * (epsilon3) - 
      sigx5_3)/(sigwz3) - sigx17_3 * ewu3 * sigx6_3/(sigwz3)^2) * 
      ewu3 * ewz3), FUN = "*"), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((wusq3) * 
      ewv3) - S^2 * sigx3_3 * sqewu3 * ewu3 * (epsilon3)^2)/(sigma_sq3) + 
      0.5 * ((ewu3/(sigma_sq3) - 1) * ewv3/(sigma_sq3) + 
        1 - 0.5 * ((wusq3) * (wvsq3)))) * depsisq3/sigmastar3 + 
      0.5 * (S^2 * sigx3_3 * depsisq3 * (epsilon3)^2/(sigma_sq3)^2)) * 
      ewu3 + sigx3_3 * (1 - 2 * (sigx2_3 * ewu3 * s3xq3)) * 
      depsisq3) * dmusig3/sqsq3^2 + 0.5 * sigx15_3) * (epsilon3) - 
      0.5 * sigx12_3)/(sigwz3) - sigx17_3 * sigx7_3/(sigwz3)^2) * 
      ewu3 * ewv3 * ewz3), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * ewu3 * 
      (sigma_sq4) * ewz3 * sigx6_3 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq3))), FUN = "*"), Xvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu3 * ewu4 * ewz3 * 
      sigx6_3 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq3))), FUN = "*"), uHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu3 * ewv4 * ewz3 * 
      sigx7_4 * sigx6_3 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq3))), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -(sigx4w3z * 
    ewu3 * ewz1 * ewz3 * sigx6_3/sqrt(sigma_sq3)), FUN = "*"), 
    Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(sigx4wz * ewu3 * ewz2 * ewz3 * 
      sigx6_3/sqrt(sigma_sq3)), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 2 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = sigx4w13z * ewu3 * ewz3 * sigx6_3/sqrt(sigma_sq3), 
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    2 * nvZVvar + 1):(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (ewv3/(sigma_sq3)) - 
      0.5 * (0.5 * (wvsq3) + ewv3/(sigma_sq3))) * (wvsq3) + 
      S^2 * sigx3_3^2 * ewu3 * ewv3 * (epsilon3)^2/(sqsq3^2 * 
        (sigma_sq3))) * depsisq3 * ewu3/sigmastar3 + 
      ((0.5 * (S^2 * depsisq3 * (epsilon3)^2/(sigma_sq3)^2) - 
        2 * (sigx3_3 * depsisq3 * s3xq3)) * ewv3 + depsisq3) * 
        sigx3_3) * dmusig3 * ewu3/sqsq3^2 + S * (0.5 * 
      (ewv3 * sigx19_3) + 0.5 * pmusig3) * depsisq3 * (epsilon3)/(sigma_sq3)^2) * 
      (epsilon3) - (0.5 * (depsisq3 * pmusig3) + 0.5 * 
      (ewv3 * sigx13_3))/(sigma_sq3))/(sigwz3) - sigx18_3 * 
      ewv3 * sigx7_3/(sigwz3)^2) * ewv3 * ewz3), FUN = "*"), 
    vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (S * prC * sigx1_4 * (sigma_sq4) * 
      ewv3 * ewz3 * sigx7_3 * sqrt(sigma_sq4)/(sigsq_4^2 * 
      wzdsq3))), FUN = "*"), Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewu4 * ewv3 * ewz3 * 
      sigx7_3 * sigx6_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq3))), FUN = "*"), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(4 * (prC * ewv3 * ewv4 * ewz3 * 
      sigx7_3 * sigx7_4 * sqrt(sigma_sq4)/((sigwz4)^2 * 
      wzdsq3))), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(sigx4w3z * 
    ewv3 * ewz1 * ewz3 * sigx7_3/sqrt(sigma_sq3)), FUN = "*"), 
    Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(sigx4wz * ewv3 * ewz2 * ewz3 * 
      sigx7_3/sqrt(sigma_sq3)), FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = sigx4w13z * ewv3 * ewz3 * sigx7_3/sqrt(sigma_sq3), 
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (3 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S^2 * ((depsisq4 * (S * (dmusig4 * 
      ewu4/sigmastar4 + S * pmusig4 * (epsilon4)) * (epsilon4)/(sigma_sq4) - 
      pmusig4) + S * dmusig4 * (duv4) * ewu4 * (epsilon4)/sqsq4)/sigsq_4 - 
      2 * (prC * sigx1_4^2/sigsq_4^2)) * prC), FUN = "*"), 
    Xvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * ((sqewu4 * dmusig4 * depsisq4 + 
      (S * (sigx10_4 - S * sqewu4 * dmusig4 * (duv4) * 
        (epsilon4)) * (epsilon4) - 0.5 * (sigx1_4/(sigma_sq4)))/(sigma_sq4))/(sigwz4) - 
      2 * (prC * sigx1_4 * sigx6_4/((sigwz4)^2 * (sigma_sq4)))) * 
      prC * ewu4), FUN = "*"), uHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(Xvar, 
    MARGIN = 1, STATS = 2 * (S * (((S * (sigx10_4 + S * sigx3_4 * 
      dmusig4 * (duv4) * ewu4 * (epsilon4)/sqsq4^2) * (epsilon4) - 
      0.5 * (sigx1_4/(sigma_sq4)))/(sigma_sq4) - sigx3_4 * 
      dmusig4 * depsisq4 * ewu4/sqsq4^2)/(sigwz4) - 2 * 
      (prC * sigx1_4 * sigx7_4/((sigwz4)^2 * (sigma_sq4)))) * 
      prC * ewv4), FUN = "*"), vHvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = -(S * 
    prC * sigx4w5z * sigx1_4 * ewz1/((sigma_sq4)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = -(S * prC * sigx4w6z * sigx1_4 * ewz2/((sigma_sq4)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(3 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    3 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Xvar, MARGIN = 1, 
    STATS = -(S * prC * sigx4w14z * sigx1_4 * ewz3/((sigma_sq4)^(3/2))), 
    FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 3 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((ewu4 * (S * (0.5 * sigx15_4 - 
      (0.5 * (S^2 * sqewu4 * depsisq4 * (epsilon4)^2/(sigma_sq4)^2) - 
        (((0.5 * (ewu4/(sigma_sq4)) + 1 - 0.5 * (0.5 * 
          (wusq4) + ewu4/(sigma_sq4))) * (wusq4) * ewv4/sigmastar4 + 
          (2 - 2 * (sigx2_4^2 * ewu4 * (sigma_sq4)/sqsq4^2)) * 
          sigmastar4)/sqsq4^2 + S^2 * sqewu4^2 * ewu4 * 
          (epsilon4)^2/sqsq4) * depsisq4) * dmusig4) * 
      (epsilon4) - 0.5 * sigx12_4) + S * sigx11_4 * (epsilon4) - 
      0.5 * dpsq4)/(sigwz4) - sigx4w9z * ewu4 * sigx6_4/(sigwz4)^2) * 
      prC * ewu4), FUN = "*"), uHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = 2 * (((S * (((((0.5 * ((wusq4) * 
      ewv4) - S^2 * sigx3_4 * sqewu4 * ewu4 * (epsilon4)^2)/(sigma_sq4) + 
      0.5 * ((ewu4/(sigma_sq4) - 1) * ewv4/(sigma_sq4) + 
        1 - 0.5 * ((wusq4) * (wvsq4)))) * depsisq4/sigmastar4 + 
      0.5 * (S^2 * sigx3_4 * depsisq4 * (epsilon4)^2/(sigma_sq4)^2)) * 
      ewu4 + sigx3_4 * (1 - 2 * (sigx2_4 * ewu4 * s3xq4)) * 
      depsisq4) * dmusig4/sqsq4^2 + 0.5 * sigx15_4) * (epsilon4) - 
      0.5 * sigx12_4)/(sigwz4) - sigx4w9z * sigx7_4/(sigwz4)^2) * 
      prC * ewu4 * ewv4), FUN = "*"), vHvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = -(prC * 
    sigx4w5z * ewu4 * ewz1 * sigx6_4/sqrt(sigma_sq4)), FUN = "*"), 
    Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(prC * sigx4w6z * ewu4 * ewz2 * 
      sigx6_4/sqrt(sigma_sq4)), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 3 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(uHvar, 
    MARGIN = 1, STATS = -(prC * sigx4w14z * ewu4 * ewz3 * 
      sigx6_4/sqrt(sigma_sq4)), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    3 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = 2 * (((S * ((((0.5 * (ewv4/(sigma_sq4)) - 
      0.5 * (0.5 * (wvsq4) + ewv4/(sigma_sq4))) * (wvsq4) + 
      S^2 * sigx3_4^2 * ewu4 * ewv4 * (epsilon4)^2/(sqsq4^2 * 
        (sigma_sq4))) * depsisq4 * ewu4/sigmastar4 + 
      ((0.5 * (S^2 * depsisq4 * (epsilon4)^2/(sigma_sq4)^2) - 
        2 * (sigx3_4 * depsisq4 * s3xq4)) * ewv4 + depsisq4) * 
        sigx3_4) * dmusig4 * ewu4/sqsq4^2 + S * (0.5 * 
      (ewv4 * sigx19_4) + 0.5 * pmusig4) * depsisq4 * (epsilon4)/(sigma_sq4)^2) * 
      (epsilon4) - (0.5 * (depsisq4 * pmusig4) + 0.5 * 
      (ewv4 * sigx13_4))/(sigma_sq4))/(sigwz4) - sigx4w11z * 
      ewv4 * sigx7_4/(sigwz4)^2) * prC * ewv4), FUN = "*"), 
    vHvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 
    nZHvar)] <- crossprod(sweep(vHvar, MARGIN = 1, STATS = -(prC * 
    sigx4w5z * ewv4 * ewz1 * sigx7_4/sqrt(sigma_sq4)), FUN = "*"), 
    Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(prC * sigx4w6z * ewv4 * ewz2 * 
      sigx7_4/sqrt(sigma_sq4)), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar), (4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(vHvar, 
    MARGIN = 1, STATS = -(prC * sigx4w14z * ewv4 * ewz3 * 
      sigx7_4/sqrt(sigma_sq4)), FUN = "*"), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar + nZHvar), (4 * nXvar + 4 * 
    nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 4 * nuZUvar + 
    4 * nvZVvar + nZHvar)] <- crossprod(sweep(Zvar, MARGIN = 1, 
    STATS = ((1 - ewz1/wzdeno)/wzdsig - 2 * (dwp1/(wzdsig^2 * 
      sqrt(sigma_sq1)))) * (2 * dpsq1 - sigx4) * ewz1, 
    FUN = "*"), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar + nZHvar), (4 * nXvar + 4 * 
    nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * nXvar + 4 * 
    nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = -(((2 * dpsq1 - sigx4)/(sigx4wzsq) + 
      2 * ((2 * dpsq2 - sigx4) * depsisq1 * pmusig1/(wzdsig^2 * 
        sqrt(sigma_sq1)))) * ewz1 * ewz2), FUN = "*"), 
    Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar + nZHvar), (4 * nXvar + 4 * 
    nuZUvar + 4 * nvZVvar + 2 * nZHvar + 1):(4 * nXvar + 
    4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = -(((2 * dpsq1 - sigx4)/(sigx4wzsq) + 
      2 * ((2 * dpsq3 - sigx4) * depsisq1 * pmusig1/(wzdsig^2 * 
        sqrt(sigma_sq1)))) * ewz1 * ewz3), FUN = "*"), 
    Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 1):(4 * 
      nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = ((1 - ewz2/wzdeno)/wzdsig - 2 * (dwp2/(wzdsig^2 * 
      sqrt(sigma_sq2)))) * (2 * dpsq2 - sigx4) * ewz2, 
    FUN = "*"), Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 
      1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = -(((2 * dpsq2 - sigx4)/(sigx4wzsq) + 
      2 * ((2 * dpsq3 - sigx4) * depsisq2 * pmusig2/(wzdsig^2 * 
        sqrt(sigma_sq2)))) * ewz2 * ewz3), FUN = "*"), 
    Zvar)
  hessll[(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 
    1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar), 
    (4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 2 * nZHvar + 
      1):(4 * nXvar + 4 * nuZUvar + 4 * nvZVvar + 3 * nZHvar)] <- crossprod(sweep(Zvar, 
    MARGIN = 1, STATS = ((1 - ewz3/wzdeno)/wzdsig - 2 * (dwp3/(wzdsig^2 * 
      sqrt(sigma_sq3)))) * (2 * dpsq3 - sigx4) * ewz3, 
    FUN = "*"), Zvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

LCM4ChnormAlgOpt <- function(start, olsParam, dataTable, S, nXvar, 
  uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, Yvar, Xvar, 
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType, 
  qac, initStart, initAlg, initIter, initFactorLB, initFactorUB) {
  if (!is.null(start)) {
    startVal <- start
  } else {
    start_st <- csLCMfhalfnorm4C(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], 
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
    initModel <- nlminb(start = startVal, objective = function(parm) -sum(cLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), lower = startMat[, 
      6], upper = startMat[, 7], control = list(iter.max = initIter, 
      trace = if (printInfo) 1 else 0, eval.max = initIter, 
      rel.tol = tol, x.tol = tol))
    startVal <- initModel$par
  }
  startLoglik <- sum(cLCMhalfnormlike4C(startVal, nXvar = nXvar, 
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
    fn = function(parm) -sum(cLCMhalfnormlike4C(parm, nXvar = nXvar, 
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, Zvar = Zvar, 
      nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = 0, 
    control = list(trace = if (printInfo) 1 else 0, maxeval = itermax, 
      stepmax = stepmax, xtol = tol, grtol = gradtol)), 
    maxLikAlgo = maxRoutine(fn = cLCMhalfnormlike4C, grad = cgradLCMhalfnormlike4C, 
      hess = chessLCMhalfnormlike4C, start = startVal, 
      finalHessian = if (hessianType == 2) "bhhh" else TRUE, 
      control = list(printLevel = if (printInfo) 2 else 0, 
        iterlim = itermax, reltol = tol, tol = tol, qac = qac), 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), sr1 = trust.optim(x = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), method = "SR1", 
      control = list(maxit = itermax, cgtol = gradtol, 
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0, 
        report.precision = 1L)), sparse = trust.optim(x = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), hs = function(parm) as(-chessLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar), "dgCMatrix"), 
      method = "Sparse", control = list(maxit = itermax, 
        cgtol = gradtol, stop.trust.radius = tol, prec = tol, 
        report.level = if (printInfo) 2 else 0, report.precision = 1L, 
        preconditioner = 1L)), mla = mla(b = startVal, 
      fn = function(parm) -sum(cLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)), hess = function(parm) -chessLCMhalfnormlike4C(parm, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar), print.info = printInfo, 
      maxiter = itermax, epsa = gradtol, epsb = gradtol), 
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), gradient = function(parm) -colSums(cgradLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar)), hessian = function(parm) -chessLCMhalfnormlike4C(parm, 
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
      S = S, Zvar = Zvar, nZHvar = nZHvar), control = list(iter.max = itermax, 
      trace = if (printInfo) 1 else 0, eval.max = itermax, 
      rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradLCMhalfnormlike4C(mleObj$par, 
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
      mleObj$hessian <- chessLCMhalfnormlike4C(parm = mleObj$par, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)
    if (method == "sr1") 
      mleObj$hessian <- chessLCMhalfnormlike4C(parm = mleObj$solution, 
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
        S = S, Zvar = Zvar, nZHvar = nZHvar)
  }
  mleObj$logL_OBS <- cLCMhalfnormlike4C(parm = mlParam, nXvar = nXvar, 
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, Zvar = Zvar, 
    nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradLCMhalfnormlike4C(parm = mlParam, 
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
    S = S, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, 
    mleObj = mleObj, mlParam = mlParam, if (is.null(start)) InitHalf = InitHalf))
}

# Posterior probabilities and efficiencies ----------

cLCM4Chalfnormeff <- function(object, level) {
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
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) + 
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) + 
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 
    1, which.max)
  P_cond_c <- ifelse(Group_c == 1, Pcond_c1, ifelse(Group_c == 
    2, Pcond_c2, ifelse(Group_c == 3, Pcond_c3, Pcond_c4)))
  u_c1 <- mustar1 + sigmastar1 * dnorm(mustar1/sigmastar1)/pnorm(mustar1/sigmastar1)
  u_c2 <- mustar2 + sigmastar2 * dnorm(mustar2/sigmastar2)/pnorm(mustar2/sigmastar2)
  u_c3 <- mustar3 + sigmastar3 * dnorm(mustar3/sigmastar3)/pnorm(mustar3/sigmastar3)
  u_c4 <- mustar4 + sigmastar4 * dnorm(mustar4/sigmastar4)/pnorm(mustar4/sigmastar4)
  u_c <- ifelse(Group_c == 1, u_c1, ifelse(Group_c == 2, u_c2, 
    ifelse(Group_c == 3, u_c3, u_c4)))
  ineff_c1 <- ifelse(Group_c == 1, u_c1, NA)
  ineff_c2 <- ifelse(Group_c == 2, u_c2, NA)
  ineff_c3 <- ifelse(Group_c == 3, u_c3, NA)
  ineff_c4 <- ifelse(Group_c == 4, u_c4, NA)
  if (object$logDepVar == TRUE) {
    teJLMS_c <- exp(-u_c)
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c, PosteriorProb_c1 = Pcond_c1, 
      PosteriorProb_c2 = Pcond_c2, PosteriorProb_c3 = Pcond_c3, PosteriorProb_c4 = Pcond_c4, 
      PriorProb_c1 = Probc1, PriorProb_c2 = Probc2, PriorProb_c3 = Probc3, PriorProb_c4 = Probc4, 
      u_c = u_c, teJLMS_c = teJLMS_c, u_c1 = u_c1, u_c2 = u_c2, u_c3 = u_c3, u_c4 = u_c4, 
      ineff_c1 = ineff_c1, ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4)
  } else {
    res <- bind_cols(Group_c = Group_c, PosteriorProb_c = P_cond_c, PosteriorProb_c1 = Pcond_c1, 
      PosteriorProb_c2 = Pcond_c2, PosteriorProb_c3 = Pcond_c3, PosteriorProb_c4 = Pcond_c4, 
      PriorProb_c1 = Probc1, PriorProb_c2 = Probc2, PriorProb_c3 = Probc3, PriorProb_c4 = Probc4, 
      u_c = u_c, u_c1 = u_c1, u_c2 = u_c2, u_c3 = u_c3, u_c4 = u_c4, ineff_c1 = ineff_c1, 
      ineff_c2 = ineff_c2, ineff_c3 = ineff_c3, ineff_c4 = ineff_c4)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargLCM4Chalfnorm_Eu <- function(object) {
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
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) + 
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) + 
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 
    1, which.max)
  margEff_c1 <- kronecker(matrix(delta1[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu1/2) * dnorm(0), ncol = 1))
  colnames(margEff_c1) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c1")
  margEff_c2 <- kronecker(matrix(delta2[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu2/2) * dnorm(0), ncol = 1))
  colnames(margEff_c2) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c2")
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu3/2) * dnorm(0), ncol = 1))
  colnames(margEff_c3) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c3")
  margEff_c4 <- kronecker(matrix(delta4[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu4/2) * dnorm(0), ncol = 1))
  colnames(margEff_c4) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in 1:ncol(margEff_c1)) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], 
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c == 
        3, margEff_c3[, c], margEff_c4[, c])))
  }
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1], 
    "_c")
  margEff <- bind_cols(as_tibble(margEff_c), as_tibble(margEff_c1), 
    as_tibble(margEff_c2), as_tibble(margEff_c3), as_tibble(margEff_c4))
  return(margEff)
}

cmargLCM4Chalfnorm_Vu <- function(object) {
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
  beta3 <- object$mlParam[(2 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar)]
  delta3 <- object$mlParam[(3 * object$nXvar + 2 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar)]
  phi3 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    2 * object$nvZVvar + 1):(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  beta4 <- object$mlParam[(3 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar)]
  delta4 <- object$mlParam[(4 * object$nXvar + 3 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar)]
  phi4 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    3 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar)]
  theta1 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 1):(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar)]
  theta2 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 2 * object$nZHvar)]
  theta3 <- object$mlParam[(4 * object$nXvar + 4 * object$nuZUvar + 
    4 * object$nvZVvar + 2 * object$nZHvar + 1):(4 * object$nXvar + 
    4 * object$nuZUvar + 4 * object$nvZVvar + 3 * object$nZHvar)]
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
  Wu3 <- as.numeric(crossprod(matrix(delta3), t(uHvar)))
  Wu4 <- as.numeric(crossprod(matrix(delta4), t(uHvar)))
  Wv1 <- as.numeric(crossprod(matrix(phi1), t(vHvar)))
  Wv2 <- as.numeric(crossprod(matrix(phi2), t(vHvar)))
  Wv3 <- as.numeric(crossprod(matrix(phi3), t(vHvar)))
  Wv4 <- as.numeric(crossprod(matrix(phi4), t(vHvar)))
  Wz1 <- as.numeric(crossprod(matrix(theta1), t(Zvar)))
  Wz2 <- as.numeric(crossprod(matrix(theta2), t(Zvar)))
  Wz3 <- as.numeric(crossprod(matrix(theta3), t(Zvar)))
  epsilon1 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta1), t(Xvar)))
  epsilon2 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta2), t(Xvar)))
  epsilon3 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta3), t(Xvar)))
  epsilon4 <- model.response(model.frame(object$formula, data = object$dataTable)) - 
    as.numeric(crossprod(matrix(beta4), t(Xvar)))
  mustar1 <- -exp(Wu1) * object$S * epsilon1/(exp(Wu1) + exp(Wv1))
  sigmastar1 <- sqrt(exp(Wu1) * exp(Wv1)/(exp(Wu1) + exp(Wv1)))
  mustar2 <- -exp(Wu2) * object$S * epsilon2/(exp(Wu2) + exp(Wv2))
  sigmastar2 <- sqrt(exp(Wu2) * exp(Wv2)/(exp(Wu2) + exp(Wv2)))
  mustar3 <- -exp(Wu3) * object$S * epsilon3/(exp(Wu3) + exp(Wv3))
  sigmastar3 <- sqrt(exp(Wu3) * exp(Wv3)/(exp(Wu3) + exp(Wv3)))
  mustar4 <- -exp(Wu4) * object$S * epsilon4/(exp(Wu4) + exp(Wv4))
  sigmastar4 <- sqrt(exp(Wu4) * exp(Wv4)/(exp(Wu4) + exp(Wv4)))
  Pi1 <- 2/sqrt(exp(Wu1) + exp(Wv1)) * dnorm(object$S * epsilon1/sqrt(exp(Wu1) + 
    exp(Wv1))) * pnorm(mustar1/sigmastar1)
  Pi2 <- 2/sqrt(exp(Wu2) + exp(Wv2)) * dnorm(object$S * epsilon2/sqrt(exp(Wu2) + 
    exp(Wv2))) * pnorm(mustar2/sigmastar2)
  Pi3 <- 2/sqrt(exp(Wu3) + exp(Wv3)) * dnorm(object$S * epsilon3/sqrt(exp(Wu3) + 
    exp(Wv3))) * pnorm(mustar3/sigmastar3)
  Pi4 <- 2/sqrt(exp(Wu4) + exp(Wv4)) * dnorm(object$S * epsilon4/sqrt(exp(Wu4) + 
    exp(Wv4))) * pnorm(mustar4/sigmastar4)
  Probc1 <- exp(Wz1)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc2 <- exp(Wz2)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc3 <- exp(Wz3)/(1 + exp(Wz1) + exp(Wz2) + exp(Wz3))
  Probc4 <- 1 - Probc1 - Probc2 - Probc3
  Pcond_c1 <- Pi1 * Probc1/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c2 <- Pi2 * Probc2/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c3 <- Pi3 * Probc3/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Pcond_c4 <- Pi4 * Probc4/(Pi1 * Probc1 + Pi2 * Probc2 + Pi3 * 
    Probc3 + Pi4 * Probc4)
  Group_c <- apply(cbind(Pcond_c1, Pcond_c2, Pcond_c3, Pcond_c4), 
    1, which.max)
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
  margEff_c3 <- kronecker(matrix(delta3[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu3) * (1 - (dnorm(0)/pnorm(0))^2), 
    ncol = 1))
  colnames(margEff_c3) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c3")
  margEff_c4 <- kronecker(matrix(delta4[2:object$nuZUvar], 
    nrow = 1), matrix(exp(Wu4) * (1 - (dnorm(0)/pnorm(0))^2), 
    ncol = 1))
  colnames(margEff_c4) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c4")
  margEff_c <- matrix(nrow = nrow(margEff_c1), ncol = ncol(margEff_c1))
  for (c in 1:ncol(margEff_c1)) {
    margEff_c[, c] <- ifelse(Group_c == 1, margEff_c1[, c], 
      ifelse(Group_c == 2, margEff_c2[, c], ifelse(Group_c == 
        3, margEff_c3[, c], margEff_c4[, c])))
  }
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1], 
    "_c")
  margEff <- bind_cols(as_tibble(margEff_c), as_tibble(margEff_c1), 
    as_tibble(margEff_c2), as_tibble(margEff_c3), as_tibble(margEff_c4))
  return(margEff)
}
