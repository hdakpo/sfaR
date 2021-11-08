########################################################
#                                                      #
# truncated normal + normal distributions              #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

ctruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon)/(exp(Wu) +
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-1/2 * log(exp(Wu) + exp(Wv)) + dnorm((mu + S * epsilon)/sqrt(exp(Wu) +
    exp(Wv)), log = TRUE) + pnorm(mustar/sigmastar, log.p = TRUE) -
    pnorm(mu/sqrt(exp(Wu)), log.p = TRUE))
  return(ll)
}

# starting value for the log-likelihood ----------

csttruncnorm <- function(olsObj, epsiRes, S, nmuZUvar, nuZUvar,
  uHvar, muHvar, nvZVvar, vHvar) {
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
  reg_hetmu <- if (nmuZUvar == 1) {
    lm(epsiRes ~ 1)
  } else {
    lm(epsiRes ~ ., data = as.data.frame(muHvar[, 2:nmuZUvar]))
  }
  if (any(is.na(reg_hetmu$coefficients)))
    stop("at least one of the OLS coefficients of 'muhet' is NA: ",
      paste(colnames(muHvar)[is.na(reg_hetmu$coefficients)],
        collapse = ", "), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE)
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  omega <- coefficients(reg_hetmu)
  names(omega) <- paste0("Zmu_", colnames(muHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, omega, delta, phi))
}

# Gradient of the likelihood function ----------

cgradtruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon)/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  mustar2 <- (mu + S * (epsilon))/sqrt(sigma_sq)
  dmustar2 <- dnorm(mustar2)
  dmustar2epsi <- dmustar2 * (mu + S * (epsilon))^2
  dmustar2epsix2 <- dmustar2epsi/(dmustar2 * (sigma_sq))
  dmustar2epsix3 <- (0.5 * (dmustar2epsix2) - 0.5)/(sigma_sq)
  sigx2 <- (sigma_sq) * sigmastar
  sigx3 <- 0.5 * ((1 - exp(Wv)/(sigma_sq)) * exp(Wu)/sigmastar) +
    sigmastar
  wusq <- exp(Wu)/(sigma_sq)
  sigx4 <- (mu * exp(Wv) - S * exp(Wu) * (epsilon))/(sigx2)^2
  sigx5 <- (0.5 * ((1 - wusq) * exp(Wv)/sigmastar) + sigmastar) *
    sigx4
  sigx6 <- sigx5 + S * (epsilon)/(sigx2)
  mustar3 <- mu + S * (epsilon)
  dmu <- dnorm(mu/exp(Wu/2))
  pmu <- pnorm(mu/exp(Wu/2))
  sigmu <- mu/(sigx2) - (sigx3) * sigx4
  pmusigx2 <- pmusig * sigmastar
  pmuwu <- exp(Wu/2) * pmu
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = S * (dmusig *
    exp(Wu)/(pmusigx2) + (mustar3))/(sigma_sq), FUN = "*"),
    sweep(muHvar, MARGIN = 1, STATS = ((dmusig * exp(Wv)/(pmusigx2) -
      (mustar3))/(sigma_sq) - dmu/(pmuwu)), FUN = "*"),
    sweep(uHvar, MARGIN = 1, STATS = ((dmustar2epsix3 - (sigx6) *
      dmusig/pmusig) * exp(Wu) + 0.5 * (mu * dmu/(pmuwu))),
      FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (dmustar2epsix3 +
      dmusig * (sigmu)/pmusig) * exp(Wv), FUN = "*")))
  return(gradll)
}

# Hessian of the likelihood function ----------

chesstruncnormlike <- function(parm, nXvar, nmuZUvar, nuZUvar,
  nvZVvar, muHvar, uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  mustar <- (mu * exp(Wv) - exp(Wu) * S * epsilon)/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(sigma_sq))
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  mustar2 <- (mu + S * (epsilon))/sqrt(sigma_sq)
  dmustar2 <- dnorm(mustar2)
  dmustar2epsi <- dmustar2 * (mu + S * (epsilon))^2
  dmustar2epsix2 <- dmustar2epsi/(dmustar2 * (sigma_sq))
  dmustar2epsix3 <- (0.5 * (dmustar2epsix2) - 0.5)/(sigma_sq)
  sigx2 <- (sigma_sq) * sigmastar
  sigx3 <- 0.5 * ((1 - exp(Wv)/(sigma_sq)) * exp(Wu)/sigmastar) +
    sigmastar
  wusq <- exp(Wu)/(sigma_sq)
  sigx4 <- (mu * exp(Wv) - S * exp(Wu) * (epsilon))/(sigx2)^2
  sigx5 <- (0.5 * ((1 - wusq) * exp(Wv)/sigmastar) + sigmastar) *
    sigx4
  sigx6 <- sigx5 + S * (epsilon)/(sigx2)
  mustar3 <- mu + S * (epsilon)
  dmu <- dnorm(mu/exp(Wu/2))
  pmu <- pnorm(mu/exp(Wu/2))
  sigmu <- mu/(sigx2) - (sigx3) * sigx4
  pmusigx2 <- pmusig * sigmastar
  pmuwu <- exp(Wu/2) * pmu
  dmusqx2 <- (dmustar2 * (sigma_sq))^2
  dmustar2epsix2sq <- dmustar2epsi/dmusqx2
  sigx7 <- 0.5 * ((((mustar3)^2/(sigma_sq) - 2)/(dmustar2 *
    (sigma_sq)) - dmustar2epsix2sq) * dmustar2 * (mustar3)/(sigma_sq))
  sigx8 <- 0.5 * (((2 - (mustar3)^2/(sigma_sq))/(dmustar2 *
    (sigma_sq)) + dmustar2epsix2sq) * dmustar2 * (mustar3)/(sigma_sq))
  muwuepsi <- mu * exp(Wv) - S * exp(Wu) * (epsilon)
  sigx9 <- (muwuepsi)/(sigx2) + dmusig/pmusig
  sigx10 <- 0.5 * ((1 - wusq) * exp(Wv)/sigmastar) + sigmastar
  sigx11 <- (sigma_sq) * (muwuepsi) * sigmastar/(sigx2)^2
  sigmustar <- (muwuepsi)/sigmastar
  sigx12 <- 0.5 * (dmustar2epsi/(sigma_sq)) + dmustar2
  wvsq <- exp(Wv)/(sigma_sq)
  sigx13 <- ((0.5 * ((0.5 * ((mustar3)^2/(dmustar2 * (sigma_sq)^3)) -
    (sigx12)/dmusqx2) * dmustar2epsi) - dmustar2epsix3) *
    exp(Wu) + 0.5 * (dmustar2epsix2) - 0.5)/(sigma_sq)
  sigx14 <- ((sigx3)/(sigx2)^2 + dmusig * (sigmu)/((sigma_sq) *
    pmusigx2))
  sigx15 <- ((0.5 * ((0.5 * ((mustar3)^2/(dmustar2 * (sigma_sq)^3)) -
    (sigx12)/dmusqx2) * dmustar2epsi) - dmustar2epsix3) *
    exp(Wv) + 0.5 * (dmustar2epsix2) - 0.5)/(sigma_sq)
  hessll <- matrix(nrow = nXvar + nmuZUvar + nuZUvar + nvZVvar,
    ncol = nXvar + nmuZUvar + nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * ((0 - 1) - ((muwuepsi)/(exp(Wv) * pmusigx2) +
      dmusig * exp(Wu)/(pmusigx2)^2) * dmusig * wusq)/(sigma_sq),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nmuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -(S * ((0 - 1) + ((muwuepsi)/(pmusigx2) +
      dmusig * exp(Wu) * exp(Wv)/(pmusigx2)^2) * dmusig/(sigma_sq))/(sigma_sq)),
    FUN = "*"), muHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar)] <- crossprod(sweep(Xvar, MARGIN = 1, STATS = S *
    (sigx7 - (((sigx10)/(sigx2)^2 - (sigx6) * dmusig/((sigma_sq) *
      pmusigx2)) * exp(Wu) - ((sigx6) * (muwuepsi)/exp(Wv) +
      1/sigmastar)/(sigma_sq)) * dmusig/pmusig) * exp(Wu),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S * (sigx7 - (sigx14 * exp(Wu) + (muwuepsi) *
      (sigmu)/((sigma_sq) * exp(Wv))) * dmusig/pmusig) *
      exp(Wv), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + 1):(nXvar +
    nmuZUvar)] <- crossprod(sweep(muHvar, MARGIN = 1, STATS = (dmu *
    (dmu/(pmuwu)^2 + mu/(exp(Wu/2)^3 * pmu)) - ((0 + 1) +
    ((muwuepsi)/(exp(Wu) * pmusigx2) + dmusig * exp(Wv)/(pmusigx2)^2) *
      dmusig * wvsq)/(sigma_sq)), FUN = "*"), muHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(uHvar,
    MARGIN = 1, STATS = ((sigx13 - (((sigx10) * (mu * exp(Wv) -
      (2 * ((sigx10) * sigx11) + 3 * (S * (epsilon))) *
        exp(Wu)) + (0.5 * (wusq) - 0.5 * (0.5 * (1 -
      wusq) + wusq)) * (1 - wusq) * exp(Wv) * sigmustar)/(sigx2)^2 +
      (sigx6)^2 * (sigx9) * exp(Wu) + S * (epsilon)/(sigx2)) *
      dmusig/pmusig) * exp(Wu) + 0.5 * (mu * (0.5 * (mu^2/(exp(Wu/2)^3 *
      pmu)) - (0.5 * (pmuwu) - 0.5 * (mu * dmu))/(pmuwu)^2) *
      dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar), (nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
    nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = (sigx15 + dmusig * (mu/(sigx2) -
      ((((3 * (mu) - 2 * ((sigx3) * (sigma_sq) * (muwuepsi) *
        sigmastar/(sigx2)^2)) * exp(Wv) - S * exp(Wu) *
        (epsilon)) * (sigx3) + (0.5 * (exp(Wv)/(sigma_sq)) -
        0.5 * (0.5 * (1 - exp(Wv)/(sigma_sq)) + exp(Wv)/(sigma_sq))) *
        (1 - exp(Wv)/(sigma_sq)) * exp(Wu) * (muwuepsi)/sigmastar)/(sigx2)^2 +
        (sigx9) * exp(Wv) * (mu/(sigx2) - (sigx3) * (muwuepsi)/(sigx2)^2)^2))/pmusig) *
      exp(Wv), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    1):(nXvar + nmuZUvar + nuZUvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = ((sigx8 - ((sigx10) * exp(Wv)/(sigx2)^2 -
      (sigx6) * ((muwuepsi)/exp(Wu) + dmusig * exp(Wv)/(pmusigx2))/(sigma_sq)) *
      dmusig/pmusig) * exp(Wu) + 0.5 * (((1 - mu^2/exp(Wu/2)^2)/(pmuwu) -
      mu * dmu/(pmuwu)^2) * dmu)), FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nmuZUvar), (nXvar + nmuZUvar +
    nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar + nvZVvar)] <- crossprod(sweep(muHvar,
    MARGIN = 1, STATS = (((1/sigmastar - (muwuepsi) * (sigmu)/exp(Wu))/(sigma_sq) -
      sigx14 * exp(Wv)) * dmusig/pmusig + sigx8) * exp(Wv),
    FUN = "*"), vHvar)
  hessll[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar),
    (nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar + nuZUvar +
      nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = (((sigx6) *
    (sigx9) * (sigmu) - ((0.5 * ((1 - wusq) * wvsq) + 0.5 *
    ((wusq - 1) * wvsq + 1 - 0.5 * ((1 - wusq) * (1 - wvsq)))) *
    sigmustar + mu * (sigx10) - (sigx3) * (2 * ((sigx10) *
    sigx11) + S * (epsilon)))/(sigx2)^2) * dmusig/pmusig +
    (0.5 * ((0.5 * ((mustar3)^2/(dmustar2 * (sigma_sq)^3)) -
      (sigx12)/dmusqx2) * dmustar2epsi) - dmustar2epsix3)/(sigma_sq)) *
    exp(Wu) * exp(Wv), FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

truncnormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  muHvar, nmuZUvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar,
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType,
  qac) {
  startVal <- if (!is.null(start))
    start else csttruncnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar)
  startLoglik <- sum(ctruncnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S)), gr = function(parm) -colSums(cgradtruncnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ctruncnormlike, grad = cgradtruncnormlike,
      hess = chesstruncnormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S),
    sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ctruncnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      gr = function(parm) -colSums(cgradtruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      method = "SR1", control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L)), sparse = trust.optim(x = startVal,
      fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S)), gr = function(parm) -colSums(cgradtruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      hs = function(parm) as(-chesstruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S),
        "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(ctruncnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S)), gr = function(parm) -colSums(cgradtruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      hess = function(parm) -chesstruncnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, S = S), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(ctruncnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      gradient = function(parm) -colSums(cgradtruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
      hessian = function(parm) -chesstruncnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S),
      control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradtruncnormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S))
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
      mleObj$hessian <- chesstruncnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)
    if (method == "sr1")
      mleObj$hessian <- chesstruncnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)
  }
  mleObj$logL_OBS <- ctruncnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradtruncnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}


# Conditional efficiencies estimation ----------

ctruncnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nmuZUvar + object$nuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar + object$nvZVvar)]
  Xvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 1)
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  vHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 4)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- (mu * exp(Wv) - exp(Wu) * object$S * epsilon)/(exp(Wu) +
    exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
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

cmargtruncnorm_Eu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Lambda <- mu/exp(Wu/2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(1 - Lambda * dnorm(Lambda)/pnorm(Lambda) - (dnorm(Lambda)/pnorm(Lambda))^2,
      ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2)/2 * ((1 + Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
      Lambda * (dnorm(Lambda)/pnorm(Lambda))^2), ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Eu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(margEff)
}

cmargtruncnorm_Vu <- function(object) {
  omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nmuZUvar)]
  delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
    1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
  muHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 3)
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Lambda <- mu/exp(Wu/2)
  m1 <- exp(Wu/2) * (Lambda + dnorm(Lambda)/pnorm(Lambda))
  m2 <- exp(Wu) * (1 - Lambda * dnorm(Lambda)/pnorm(Lambda) -
    (dnorm(Lambda)/pnorm(Lambda))^2)
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(1/exp(Wu/2) * dnorm(Lambda)/pnorm(Lambda) * (m1^2 -
      m2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - 1/2 * dnorm(Lambda)/pnorm(Lambda) *
      (Lambda + Lambda^3 + (2 + 3 * Lambda^2) * dnorm(Lambda)/pnorm(Lambda) +
        2 * Lambda * (dnorm(Lambda)/pnorm(Lambda))^2)),
      ncol = 1))
  idTRUE_mu <- substring(names(omega)[-1], 5) %in% substring(names(delta)[-1],
    4)
  idTRUE_Wu <- substring(names(delta)[-1], 4) %in% substring(names(omega)[-1],
    5)
  margEff <- cbind(mu_mat[, idTRUE_mu] + Wu_mat[, idTRUE_Wu],
    mu_mat[, !idTRUE_mu], Wu_mat[, !idTRUE_Wu])
  colnames(margEff) <- paste0("Vu_", c(colnames(muHvar)[-1][idTRUE_mu],
    colnames(muHvar)[-1][!idTRUE_mu], colnames(uHvar)[-1][!idTRUE_Wu]))
  return(margEff)
}
