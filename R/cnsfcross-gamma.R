######################################################## 
#                                                      #
# Contaminated (cnsf) normal - gamma                   #
#                                                      #
######################################################## 

# Log-likelihood ----------

ccnsfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, 
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
      return(wHvar*log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

cmcesfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, 
                                uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta1 <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  delta2 <- parm[(nXvar + nuZUvar + 1):(nXvar + 2 * nuZUvar)]
  phi1 <- parm[(nXvar + 2 * nuZUvar + 1):(nXvar + 2 * nuZUvar + 
                                            nvZVvar)]
  phi2 <- parm[(nXvar + 2 * nuZUvar + nvZVvar + 1):(nXvar + 2 * 
                                                      nuZUvar + 2 * nvZVvar)]
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
    Pi2 <- exp(S * epsilon/exp(Wu2/2) + exp(Wv1)/(2 * exp(Wu2))) * 
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu2/2)) * 
      exp(-P2 * Wu2/2)/gamma(P2) * Hi2
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
    ifelse(Probc1 * Pi1 + Probc2 * Pi2 <= 0, return(NA), 
           return(wHvar*log(Probc1 * Pi1 + Probc2 * Pi2)))
  }
}

# starting value for the log-likelihood ----------

cstcnsfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar, 
  nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + gamma - normal distributions...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj, 
    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = uHvar[, 1, drop = FALSE], nvZVvar = 1, vHvar = vHvar[, 1, drop = FALSE]), 
    grad = cgradgammanormlike, method = "BFGS", control = list(iterlim = itermax, 
      printLevel = printInfo, reltol = tol), nXvar = nXvar, nuZUvar = 1, 
    nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[, 
      1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat)
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:nXvar], Esti[nXvar + 1], if (nuZUvar > 
    1) rep(0, nuZUvar - 1), 0.95*Esti[nXvar + 2], if (nvZVvar > 
    1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], 
    if (nvZVvar > 1) rep(0, nvZVvar - 1), Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
    colnames(uHvar)), paste0("Zv_", colnames(vHvar)), paste0("Zv_", 
    colnames(vHvar)), "P", paste0("CN_", colnames(Zvar)))
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
                                                              colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]), "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

cstmcesfgammanorm <- function(olsObj, epsiRes, nXvar, nuZUvar, 
                              nvZVvar, uHvar, vHvar, Yvar, Xvar, S, wHvar, Zvar, nZHvar, N, FiMat, itermax, printInfo, tol) {
  cat("Initialization: SFA + gamma - normal distributions...\n")
  initGamma <- maxLik(logLik = cgammanormlike, start = cstgammanorm(olsObj = olsObj, 
                                                                    epsiRes = epsiRes, S = S, nuZUvar = 1, uHvar = as.matrix(uHvar[, 
                                                                                                                                   1]), nvZVvar = 1, vHvar = as.matrix(vHvar[, 1])), 
                      grad = cgradgammanormlike, method = "BFGS", control = list(iterlim = 2000, 
                                                                                 printLevel = 2, reltol = 1e-14), nXvar = nXvar, nuZUvar = 1, 
                      nvZVvar = 1, uHvar = as.matrix(uHvar[, 1]), vHvar = as.matrix(vHvar[, 
                                                                                          1]), Yvar = Yvar, Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat)
  Esti <- initGamma$estimate
  StartVal <- c(Esti[1:nXvar], 0.95*Esti[nXvar + 1], if (nuZUvar > 
                                                         1) rep(0, nuZUvar - 1), 1.05 * Esti[nXvar + 1], if (nuZUvar > 
                                                                                                             1) rep(0, nuZUvar - 1), 0.95*Esti[nXvar + 2], if (nvZVvar > 
                                                                                                                                                               1) rep(0, nvZVvar - 1), 1.05 * Esti[nXvar + 2], 
                if (nvZVvar > 1) rep(0, nvZVvar - 1), 0.95*Esti[nXvar + 3], 1.05*Esti[nXvar + 3], rep(0, nZHvar))
  names(StartVal) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
                                                    colnames(uHvar)), paste0("Zu_", colnames(uHvar)), paste0("Zv_", 
                                                                                                             colnames(vHvar)), paste0("Zv_", colnames(vHvar)), "P", 
                       "P", paste0("MCE_", colnames(Zvar)))
  names(initGamma$estimate) <- c(names(Esti)[1:nXvar], paste0("Zu_", 
                                                              colnames(uHvar)[1]), paste0("Zv_", colnames(vHvar)[1]), "P")
  return(list(StartVal = StartVal, initGamma = initGamma))
}

# Gradient of the likelihood function ----------

cgradcnsfgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, 
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
  ewv1 <- exp(Wv1)
  ewv2 <- exp(Wv2)
  ewu <- exp(Wu)
  ewu_h <- exp(Wu/2)
  ewv1_h <- exp(Wv1/2)
  ewv2_h <- exp(Wv2/2)
  ewz <- exp(Wz)
  wzdeno <- (1+ewz)
  prC <- (1-ewz/wzdeno)
  wusi1 <- (ewv1/ewu_h+S*(epsilon))
  wusi2 <- (ewv2/ewu_h+S*(epsilon))
  musig1 <- wusi1/ewv1_h
  musig2 <- wusi2/ewv2_h
  pmusig1 <- pnorm(musig1)
  pmusig2 <- pnorm(musig2)
  dmusig1 <- dnorm(musig1,0,1)
  dmusig2 <- dnorm(musig2,0,1)
  F1_1 <- sweep((1-FiMat), MARGIN = 1, STATS = pmusig1, FUN = "*") + FiMat
  F1_2 <- sweep((1-FiMat), MARGIN = 1, STATS =pmusig2 , FUN = "*") + FiMat
  F2_1 <- qnorm(F1_1)
  F2_2 <- qnorm(F1_2)
  F3_1 <- dnorm(F2_1)
  F3_2 <- dnorm(F2_2)
  epsi1 <- (ewv1_h/ewu_h+S*(epsilon)/ewv1_h)
  depsi1 <- dnorm(-epsi1,0,1)
  pepsi1 <- pnorm(-epsi1)
  F4_1 <- sweep((1-FiMat)/F3_1, MARGIN = 1, STATS = dmusig1, FUN = "*")
  F4_2 <- sweep((1-FiMat)/F3_2, MARGIN = 1, STATS = dmusig2, FUN = "*")
  F5_1 <- sweep((1-FiMat)/F3_1, MARGIN = 1, STATS = dmusig1*(ewv1/ewu_h-0.5*wusi1), FUN = "*")
  F5_2 <- sweep((1-FiMat)/F3_2, MARGIN = 1, STATS = dmusig2*(ewv2/ewu_h-0.5*wusi2), FUN = "*")
  F6_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = (ewv1_h), FUN = "*"), MARGIN = 1, STATS = wusi1, FUN = "-")
  F6_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = (ewv2_h), FUN = "*"), MARGIN = 1, STATS = wusi2, FUN = "-")
  sDiv1 <- apply(F6_1^(P-1), 1, sum)
  sDiv2 <- apply(F6_2^(P-1), 1, sum)
  sigx1 <- (prC*sDiv2+ewz*sDiv1/wzdeno)
  F7_1 <- sweep((0.5-0.5*(F4_1))*F6_1^(P-2)*(P-1), MARGIN = 1, STATS =ewv1/ewu_h , FUN = "*")
  F7_2 <- sweep((0.5-0.5*(F4_2))*F6_2^(P-2)*(P-1), MARGIN = 1, STATS = ewv2/ewu_h, FUN = "*")
  F8_1 <- sweep(sweep(F2_1, MARGIN = 1, STATS = 0.5*(ewv1_h), FUN = "*"), MARGIN = 1, STATS = ewv1/ewu_h, FUN = "-")
  F8_2 <- sweep(sweep(F2_2, MARGIN = 1, STATS = 0.5*(ewv2_h), FUN = "*"), MARGIN = 1, STATS = ewv2/ewu_h, FUN = "-")
  gx1 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx1[, k] <- apply(sweep(S*(1-F4_1)*F6_1^(P-2)*(P-1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1*wzdeno)
  }
  gx2 <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx2[, k] <- apply(sweep(S*(1-F4_2)*F6_2^(P-2)*(P-1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gx <- sweep(Xvar, MARGIN = 1, STATS = S*(depsi1/(ewv1_h*pepsi1)-1/ewu_h), FUN = "*") + gx1 + gx2
  gu1 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu1[, k] <- apply(sweep(F7_1, MARGIN = 1, STATS = uHvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1*wzdeno)
  }
  gu2 <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu2[, k] <- apply(sweep(F7_2, MARGIN = 1, STATS = uHvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gu <- sweep(uHvar, MARGIN = 1, STATS = ((0.5*(depsi1*ewv1_h/pepsi1)-0.5*(S*(epsilon)))/ewu_h-(0.5*P+2*(ewu*ewv1/(2*ewu)^2))), FUN = "*") + gu1 + gu2 
  gv1 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv1[, k] <- apply(sweep((F5_1+F8_1)*F6_1^(P-2)*(P-1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1, sum) * ewz/(sigx1*wzdeno)
  }
  gv1 <- sweep(vHvar, MARGIN = 1, STATS = (ewv1/(2*ewu)-(0.5*(ewv1_h/ewu_h)-0.5*(S*(epsilon)/ewv1_h))*depsi1/pepsi1), FUN = "*") + gv1
  gv2 <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv2[, k] <- apply(sweep((F5_2+F8_2)*F6_2^(P-2)*(P-1), MARGIN = 1, STATS = vHvar[, k], FUN = "*"), 1, sum) * prC/sigx1
  }
  gradll <- cbind(gx, gu, gv1, gv2, apply(F6_2^(P-1)*log(F6_2), 1, sum) * prC/sigx1 + apply(F6_1^(P-1)*log(F6_1), 1, sum) * ewz/(sigx1*wzdeno)-(0.5*(Wu)+digamma(P)), 
                  sweep(Zvar, MARGIN = 1, STATS = prC*ewz*(sDiv1-sDiv2)/(sigx1*wzdeno), FUN = "*"))
  return(sweep(gradll, MARGIN = 1, STATS = wHvar, FUN = "*"))
}



# Optimization using different algorithms ----------

cnsfgammanormAlgOpt <- function(start, olsParam, dataTable, S, 
                                nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Zvar, nZHvar, 
                                Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol, 
                                hessianType, qac) {
  start_st <- if (!is.null(start)) 
    start else cstcnsfgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]], 
                                S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar, 
                                uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar, nvZVvar = nvZVvar, 
                                nXvar = nXvar, Xvar = Xvar, Yvar = Yvar, itermax = itermax, 
                                printInfo = printInfo, tol = tol)
  InitGamma <- start_st$initGamma
  startVal <- start_st$StartVal
  
  startLoglik <- sum(ccnsfgammanormlike(startVal, nXvar = nXvar, 
                                        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
                                        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N, 
                                        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = maxBFGS, bhhh = maxBHHH, 
                         nr = maxNR, nm = maxNM)
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal, 
                                           fn = function(parm) -sum(ccnsfgammanormlike(parm, nXvar = nXvar, 
                                                                                       nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
                                                                                       vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N, 
                                                                                       FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), gr = function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                                                                                                                          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                          S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                           hessian = if (hessianType != 2) 1 else 0, control = list(trace = printInfo, 
                                                                                                    maxeval = itermax, stepmax = stepmax, xtol = tol, 
                                                                                                    grtol = gradtol)), maxLikAlgo = maxRoutine(fn = ccnsfgammanormlike, 
                                                                                                                                               grad = cgradcnsfgammanormlike, start = startVal, finalHessian = if (hessianType == 
                                                                                                                                                                                                                   2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2, 
                                                                                                                                                                                                                                                       iterlim = itermax, reltol = tol, tol = tol, qac = qac), 
                                                                                                                                               nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                               uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                               S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar), 
                   sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ccnsfgammanormlike(parm, 
                                                                                               nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                               uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                               S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                     gr = function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                         nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                         uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                         S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                     method = "SR1", control = list(maxit = itermax, cgtol = gradtol, 
                                                                    stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 4L else 0, 
                                                                    report.precision = 1L)), sparse = trust.optim(x = startVal, 
                                                                                                                  fn = function(parm) -sum(ccnsfgammanormlike(parm, 
                                                                                                                                                              nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                              uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                              S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                  gr = function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                                                                                                      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                      S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                  hs = function(parm) as(calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                                                                                                                                           nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                                           uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                                           S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                                                            var = (unname(parm)), accuracy = 4), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax, 
                                                                                                                                                                                                                                                 cgtol = gradtol, stop.trust.radius = tol, prec = tol, 
                                                                                                                                                                                                                                                 report.level = if (printInfo) 4L else 0, report.precision = 1L, 
                                                                                                                                                                                                                                                 preconditioner = 1L)), mla = mla(b = startVal, 
                                                                                                                                                                                                                                                                                  fn = function(parm) -sum(ccnsfgammanormlike(parm, 
                                                                                                                                                                                                                                                                                                                              nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                                                                                                                                                              uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                                                                                                                                                              S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                                                                                                                                                                                  gr = function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                                                                                                                                                                                                                                                                      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                                                                                                                                                                      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                                                                                                                                                                      S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                                                                                                                                                                                  print.info = printInfo, maxiter = itermax, epsa = gradtol, 
                                                                                                                                                                                                                                                                                  epsb = gradtol), nlminb = nlminb(start = startVal, 
                                                                                                                                                                                                                                                                                                                   objective = function(parm) -sum(ccnsfgammanormlike(parm, 
                                                                                                                                                                                                                                                                                                                                                                      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                                                                                                                                                                                                      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                                                                                                                                                                                                      S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                                                                                                                                                                                                                   gradient = function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                                                                                                                                                                                                                                                                                                             nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                                                                                                                                                                                                                                                                                                             uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                                                                                                                                                                                                                                                                                                             S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                                                                                                                                                                                                                                                                                                   control = list(iter.max = itermax, trace = printInfo, 
                                                                                                                                                                                                                                                                                                                                  eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradcnsfgammanormlike(mleObj$par, 
                                                      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                      S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar))
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
    if (method == "ucminf") 
      mleObj$hessian <- mleObj$hessian
    if (method == "nlminb") 
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                          S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                           var = (unname(mleObj$par)), accuracy = 4)
    if (method == "sr1") 
      mleObj$hessian <- calculus::jacobian(function(parm) -colSums(cgradcnsfgammanormlike(parm, 
                                                                                          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                                                                          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                                                                          S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)), 
                                           var = (unname(mleObj$solution)), accuracy = 4)
  }
  mleObj$logL_OBS <- ccnsfgammanormlike(parm = mlParam, nXvar = nXvar, 
                                        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, 
                                        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N, 
                                        FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  mleObj$gradL_OBS <- cgradcnsfgammanormlike(parm = mlParam, 
                                             nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
                                             uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, 
                                             S = S, N = N, FiMat = FiMat, Zvar = Zvar, nZHvar = nZHvar)
  return(list(startVal = startVal, startLoglik = startLoglik, 
              mleObj = mleObj, mlParam = mlParam, InitGamma = InitGamma))
}

# Conditional efficiencies estimation ----------

ccnsfgammanormeffI <- function(parm) {
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
    Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * 
      pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) * 
      exp(-P * Wu/2)/gamma(P) * Hi2_1
    Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) * 
      pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) * 
      exp(-P * Wu/2)/gamma(P) * Hi2_2
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  # odRatio <- Pcond_c2/(1 - Pcond_c2)
  u_c1 <- Hi1_1/Hi2_1
  u_c2 <- Hi1_2/Hi2_2
  u_c <- ifelse(Group_c == 1, u_c1, u_c2)
  u_Wc <- Pcond_c1 * u_c1 + Pcond_c2 * u_c2
  teJLMS_c1 <- exp(-u_c1)
  teJLMS_c2 <- exp(-u_c2)
  teJLMS_c <- ifelse(Group_c == 1, teJLMS_c1, teJLMS_c2)
  teJLMS_Wc <- exp(-u_Wc)
  res <- bind_cols(Group_c = Group_c, u_c1 = u_c1, u_c2 = u_c2, 
                   u_c = u_c, u_Wc = u_Wc, teJLMS_c1 = teJLMS_c1, teJLMS_c2 = teJLMS_c2, 
                   teJLMS_c = teJLMS_c, teJLMS_Wc = teJLMS_Wc,
                   Pcond_c1 = Pcond_c1, Pcond_c2 = Pcond_c2)
  return(res)
}

# Marginal effects on inefficiencies ----------

ccnsfmarggammanorm_EuI <- function(object) {
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
  Hi2_1 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i, 
                                                                       ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P - 
                                                                                                                                            1))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i, 
                                                                       ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P - 
                                                                                                                                            1))
  }
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * 
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) * 
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) * 
    pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) * 
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <-  kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), 
                         matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), 
                        matrix(P/2 * exp(Wu/2), ncol = 1))
  margEff_c <- sapply(1:N, FUN = function(x) if (Group_c[x] == 
                                                 1) 
    margEff1[x] else margEff2[x])
  margEff_Wc <- sapply(1:N, FUN = function(x) if (Group_c[x] == 
                                                  1) 
    Pcond_c1[x] * margEff1[x] else Pcond_c2[x] * margEff2[x])
  colnames(margEff1) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Eu_", colnames(uHvar)[-1])
  colnames(margEff_Wc) <- paste0("Eu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c, margEff_Wc))
}

ccnsfmarggammanorm_VuI <- function(object) {
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
  Hi1_1 <- numeric(object$Nobs)
  Hi2_1 <- numeric(object$Nobs)
  Hi1_2 <- numeric(object$Nobs)
  Hi2_2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi2_1[i] <- mean((mui1[i] + sqrt(exp(Wv1[i])) * qnorm(object$FiMat[i, 
                                                                       ] + (1 - object$FiMat[i, ]) * pnorm(-mui1[i]/sqrt(exp(Wv1[i])))))^(P - 
                                                                                                                                            1))
    Hi2_2[i] <- mean((mui2[i] + sqrt(exp(Wv2[i])) * qnorm(object$FiMat[i, 
                                                                       ] + (1 - object$FiMat[i, ]) * pnorm(-mui2[i]/sqrt(exp(Wv2[i])))))^(P - 
                                                                                                                                            1))
  }
  Pi1 <- exp(S * epsilon/exp(Wu/2) + exp(Wv1)/(2 * exp(Wu))) * 
    pnorm(-S * epsilon/exp(Wv1/2) - exp(Wv1/2)/exp(Wu/2)) * 
    exp(-P * Wu/2)/gamma(P) * Hi2_1
  Pi2 <- exp(S * epsilon/exp(Wu/2) + exp(Wv2)/(2 * exp(Wu))) * 
    pnorm(-S * epsilon/exp(Wv2/2) - exp(Wv2/2)/exp(Wu/2)) * 
    exp(-P * Wu/2)/gamma(P) * Hi2_2
  Probc1 <- exp(Wz)/(1 + exp(Wz))
  Probc2 <- 1 - Probc1
    Probc1 <- exp(Wz)/(1 + exp(Wz))
    Probc2 <- 1 - Probc1
  Pcond_c1 <- Probc1 * Pi1/(Probc1 * Pi1 + Probc2 * Pi2)
  Pcond_c2 <- Probc2 * Pi2/(Probc1 * Pi1 + Probc2 * Pi2)
  Group_c <- ifelse(Pcond_c1 > Pcond_c2, 1, 2)
  margEff1 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), 
                        matrix(P * exp(Wu), ncol = 1))
  margEff2 <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1), 
                        matrix(P * exp(Wu), ncol = 1))
  margEff_c <- sapply(1:N, FUN = function(x) if (Group_c[x] == 
                                                 1) 
    margEff1[x] else margEff2[x])
  margEff_Wc <- sapply(1:N, FUN = function(x) if (Group_c[x] == 
                                                  1) 
    Pcond_c1[x] * margEff1[x] else Pcond_c2[x] * margEff2[x])
  colnames(margEff1) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff2) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_c) <- paste0("Vu_", colnames(uHvar)[-1])
  colnames(margEff_Wc) <- paste0("Vu_", colnames(uHvar)[-1])
  return(bind_cols(margEff1, margEff2, margEff_c, margEff_Wc))
}