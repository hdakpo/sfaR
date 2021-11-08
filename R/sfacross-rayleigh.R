########################################################
#                                                      #
# Rayleigh + normal distributions                      #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

craynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-Wu - 1/2 * Wv - (S * epsilon)^2/(2 * exp(Wv)) + 1/2 *
    (mustar/sigmastar)^2 + 1/2 * log(sigmastar^2) + log(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar)))
  return(ll)
}

# starting value for the log-likelihood ----------

cstraynorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  if (S * m3 > 0) {
    varu <- exp(1/2 * log(((2 * (S * m3)^2/(pi * (pi - 3)^2))^2)^(2/6)))
  } else {
    varu <- (2 * m3^2/(pi * (3 - pi)^2))^(2/6)
  }
  if (m2 < ((4 - pi)/2) * varu) {
    varv <- abs(2 - ((4 - pi)/2) * varu)
  } else {
    varv <- m2 - ((4 - pi)/2) * varu
  }
  dep_u <- 1/2 * log(((epsiRes^2 - varv) * 2/(4 - pi))^2)
  dep_v <- 1/2 * log((epsiRes^2 - (4 - pi/2) * varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar]))
  }
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
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu * pi/2), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------

cgradraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  pdmusig <- dmusig * sigmastar - S * exp(Wu) * pmusig * (epsilon)/(sigma_sq)
  sigx2 <- (sigma_sq) * sigmastar
  sigx3 <- (1/(sigx2) - (0.5 * ((1 - exp(Wu)/(sigma_sq)) *
    exp(Wv)/sigmastar) + sigmastar) * exp(Wu)/(sigx2)^2)
  sigx4 <- (sigmastar/exp(Wv) - exp(Wu)/(sigx2))
  sigx5 <- ((sigx2)^2 * sigmastar)
  wusq <- exp(Wu)/(sigma_sq)
  pmusigepsi <- S * pmusig * (epsilon)
  pmusigepsisq <- pmusigepsi/(sigma_sq)
  dmusig2 <- dmusig/sigmastar
  wuepsisq <- exp(Wu) * (epsilon)^2
  sigx6 <- ((0.5 * ((1 - exp(Wv)/(sigma_sq)) * dmusig2) + pmusigepsisq)/(pdmusig) -
    S^2 * (0.5 * ((1 - exp(Wv)/(sigma_sq)) * exp(Wu)/sigmastar) +
      sigmastar) * wuepsisq/sigx5)
  wvsig <- exp(Wv)/sigmastar
  dwsig <- dmusig * wvsig
  dmusig4epsi <- S * dmusig * sigx4 * (epsilon)
  pdmusigsq <- (pdmusig) * (sigma_sq)
  wvsigmasq <- (1 - exp(Wv)/(sigma_sq))
  pdmusigepsi <- pmusig + dmusig4epsi
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = S * (((pdmusigepsi)/(pdmusig) -
    S * (epsilon)/exp(Wv)) * wusq + S * (epsilon)/exp(Wv)),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (((0.5 *
    (dwsig) - pmusigepsi) * exp(Wu)/(pdmusigsq) + 0.5) *
    (1 - wusq) + S^2 * (sigx3) * exp(Wu)^2 * (epsilon)^2/(sigx2) -
    1), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (((sigx6) *
    wusq + 2 * ((S * (epsilon))^2/(2 * exp(Wv))^2)) * exp(Wv) +
    0.5 * (wvsigmasq) - 0.5), FUN = "*")))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessraynormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  pdmusig <- dmusig * sigmastar - S * exp(Wu) * pmusig * (epsilon)/(sigma_sq)
  sigx2 <- (sigma_sq) * sigmastar
  sigx3 <- (1/(sigx2) - (0.5 * ((1 - exp(Wu)/(sigma_sq)) *
    exp(Wv)/sigmastar) + sigmastar) * exp(Wu)/(sigx2)^2)
  sigx4 <- (sigmastar/exp(Wv) - exp(Wu)/(sigx2))
  sigx5 <- ((sigx2)^2 * sigmastar)
  wusq <- exp(Wu)/(sigma_sq)
  pmusigepsi <- S * pmusig * (epsilon)
  pmusigepsisq <- pmusigepsi/(sigma_sq)
  dmusig2 <- dmusig/sigmastar
  wuepsisq <- exp(Wu) * (epsilon)^2
  sigx6 <- ((0.5 * ((1 - exp(Wv)/(sigma_sq)) * dmusig2) + pmusigepsisq)/(pdmusig) -
    S^2 * (0.5 * ((1 - exp(Wv)/(sigma_sq)) * exp(Wu)/sigmastar) +
      sigmastar) * wuepsisq/sigx5)
  wvsig <- exp(Wv)/sigmastar
  dwsig <- dmusig * wvsig
  dmusig4epsi <- S * dmusig * sigx4 * (epsilon)
  pdmusigsq <- (pdmusig) * (sigma_sq)
  wvsigmasq <- (1 - exp(Wv)/(sigma_sq))
  wvsigmasq2 <- exp(Wv)/(sigma_sq)
  pdmusigepsi <- pmusig + dmusig4epsi
  wuwvsig <- (1 - wusq) * wvsig
  wuwvsig2 <- 0.5 * (wuwvsig) + sigmastar
  dmusigwu <- dmusig * exp(Wu)
  dmusig3epsisq <- sigx3 * dmusigwu * (epsilon)^2
  wdmu2sig <- (wvsigmasq) * dmusig2
  dmusigepsi <- dmusig * (epsilon)
  wvsigmasq2wu <- 0.5 * (wvsigmasq * exp(Wu)/sigmastar) + sigmastar
  wuepsi <- exp(Wu) * (epsilon)
  dmusigwuepsi <- dmusigwu * (epsilon)
  wdmu1sig <- (wvsigmasq) * dmusig
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * (((((S^2 * wuepsisq/((sigma_sq) * exp(Wv)) -
      1) * (sigx4) + exp(Wu)/(sigx2)) * dmusig - exp(Wu) *
      (pdmusigepsi)^2/(pdmusigsq))/(pdmusig) + 1/exp(Wv)) *
      wusq - 1/exp(Wv)), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (((pmusig - 0.5 * (S * dmusigwuepsi/(sigx2)))/(pdmusigsq) -
      (0.5 * (dwsig) - pmusigepsi) * exp(Wu) * (pdmusigepsi)/(pdmusigsq)^2) *
      (1 - wusq) - 2 * (S * (sigx3) * wuepsi/(sigx2))) *
      exp(Wu), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (((exp(Wu) * (S * (0.5 * ((wvsigmasq)/exp(Wv)) +
      1/(sigma_sq)) * dmusigepsi/sigmastar - (0.5 * (wdmu2sig) +
      pmusigepsisq) * (pdmusigepsi)/(pdmusig)) - pmusig)/(pdmusigsq) +
      2 * (S * (wvsigmasq2wu) * wuepsi/(sigx5))) * wusq -
      4 * (S * (epsilon)/(2 * exp(Wv))^2)) * exp(Wv), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = ((((0.5 *
    (dwsig) + S * (S * (sigx3) * dmusigwuepsi - pmusig) *
    (epsilon) - ((0.5 * (dwsig) - pmusigepsi) * wusq + 0.5 *
    (0.5 * ((1 - wusq) * dwsig) + S^2 * dmusig3epsisq)))/(pdmusig) -
    0.5)/(sigma_sq) - ((0.5 * (dwsig) - pmusigepsi) * (1 -
    wusq) + pdmusig) * (0.5 * (dwsig) - pmusigepsi) * exp(Wu)/(pdmusigsq)^2) *
    (1 - wusq) + S^2 * ((2 * (sigx3) - ((0.5 * (wusq) + 1 -
    0.5 * (0.5 * (1 - wusq) + wusq)) * wuwvsig + (2 - 2 *
    ((wuwvsig2)^2 * exp(Wu) * (sigma_sq)/(sigx2)^2)) * sigmastar) *
    exp(Wu)/(sigx2)^2)/(sigx2) - (wuwvsig2) * (sigx3) * exp(Wu)/(sigx2)^2) *
    wuepsisq) * exp(Wu), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = (((((0.5 * (exp(Wv) * (S^2 * (wvsigmasq2wu) *
      dmusigwu * wuepsisq/(sigx5) - dmusig)/(sigma_sq) -
      0.5 * (wdmu1sig)) + 0.5 * dmusig) * (wvsigmasq)/sigmastar +
      (exp(Wv) * (S * (S * (wvsigmasq2wu) * dmusigwuepsi/(sigx2)^2 -
        pmusig/(sigma_sq)) * (epsilon) - (0.5 * (wdmu2sig) +
        pmusigepsisq)^2 * exp(Wu)/(pdmusig)) + pmusigepsi)/(sigma_sq))/(pdmusig) -
      ((sigx6) * wvsigmasq2 + S^2 * ((wvsigmasq2wu) * (1/(sigx5) -
        (0.5 * ((sigx2)^2 * (wvsigmasq)/(sigx2)) + 2 *
          ((wvsigmasq2wu) * exp(Wv))) * exp(Wu) * exp(Wv)/(sigx5)^2) +
        (0.5 * (wvsigmasq2) - 0.5 * (0.5 * (wvsigmasq) +
          wvsigmasq2)) * (wvsigmasq) * (sigma_sq)/((sigx2)^2 *
          exp(Wv))) * wuepsisq)) * exp(Wu) - 0.5 * (wvsigmasq))/(sigma_sq) +
      (2 - 16 * (exp(Wv)^2/(2 * exp(Wv))^2)) * (S * (epsilon))^2/(2 *
        exp(Wv))^2) * exp(Wv), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = (((0.5 * (wdmu1sig) + 0.5 * ((dmusig * wvsigmasq2 -
      S^2 * (wvsigmasq) * dmusig3epsisq/sigmastar) * wusq -
      0.5 * ((1 - wusq) * wdmu1sig)))/sigmastar + (pmusigepsi -
      ((0.5 * (wdmu2sig) + pmusigepsisq) * (0.5 * (dwsig) -
        pmusigepsi) * (1 - wusq)/(pdmusig) + S * (pmusig/(sigma_sq) +
        S * (sigx3) * dmusigepsi) * (epsilon)) * exp(Wu))/(sigma_sq))/(pdmusig) +
      0.5/(sigma_sq) - ((sigx6)/(sigma_sq) + S^2 * (((0.5 *
      ((1 - wusq) * wvsigmasq2) + 0.5 * ((wusq - 1) * wvsigmasq2 +
      1 - 0.5 * ((1 - wusq) * (wvsigmasq))) + 0.5 * (wvsigmasq)) *
      exp(Wu)/sigmastar + sigmastar)/(sigx5) + (wvsigmasq2wu) *
      (1/(sigx5) - (0.5 * ((sigx2)^2 * (1 - wusq)/(sigx2)) +
        2 * ((wuwvsig2) * exp(Wu))) * exp(Wu) * exp(Wv)/(sigx5)^2)) *
      (epsilon)^2) * exp(Wu)) * exp(Wu) * wvsigmasq2, FUN = "*"),
    vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

raynormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstraynorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(craynormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(craynormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradraynormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = craynormlike, grad = cgradraynormlike,
      hess = chessraynormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(craynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(craynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chessraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(craynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chessraynormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(craynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgradraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chessraynormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradraynormlike(mleObj$par,
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
      mleObj$hessian <- chessraynormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chessraynormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- craynormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradraynormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

craynormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * object$S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  u <- (mustar * sigmastar * dnorm(mustar/sigmastar) + (mustar^2 +
    sigmastar^2) * pnorm(mustar/sigmastar))/(sigmastar *
    dnorm(mustar/sigmastar) + mustar * pnorm(mustar/sigmastar))
  m <- ifelse(mustar/2 + sqrt(sigmastar^2 + mustar^2/4) > 0,
    mustar/2 + sqrt(sigmastar^2 + mustar^2/4), 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- exp(-mustar + sigmastar^2/2) * (sigmastar * dnorm(mustar/sigmastar -
      sigmastar) + (mustar - sigmastar^2) * pnorm(mustar/sigmastar -
      sigmastar))/(sigmastar * dnorm(mustar/sigmastar) +
      mustar * pnorm(mustar/sigmastar))
    teMO <- exp(-m)
    res <- bind_cols(u = u, teJLMS = teJLMS, teBC = teBC,
      m = m, teMO = teMO)
  } else {
    res <- bind_cols(u = u, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargraynorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * 1/2 * sqrt(pi/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmargraynorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (4 - pi)/2, ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
