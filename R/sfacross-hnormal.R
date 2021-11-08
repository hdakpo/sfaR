########################################################
#                                                      #
# half-normal + normal distributions                   #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

chalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -exp(Wu) * S * epsilon/(exp(Wu) + exp(Wv))
  sigmastar <- sqrt(exp(Wu) * exp(Wv)/(exp(Wu) + exp(Wv)))
  ll <- (-log(1/2) - 1/2 * log(exp(Wu) + exp(Wv)) + dnorm(epsilon/sqrt(exp(Wu) +
    exp(Wv)), log = TRUE) + pnorm(mustar/sigmastar, log.p = TRUE))
  return(ll)
}

# starting value for the log-likelihood ----------

csthalfnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
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
    beta <- c(olsObj[1] + S * sqrt(varu * 2/pi), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------

cgradhalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvsq <- exp(Wv)/(sigma_sq)
  mustar <- -exp(Wu) * S * epsilon/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * wvsq)
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  sigx2 <- (sigma_sq) * sigmastar
  dmusigwu <- dmusig * exp(Wu)
  pmusig2 <- pmusig * sigmastar
  depsi <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  depsisq <- depsi * (sigma_sq)^2
  depsisqx2 <- 0.5 * (S * depsi * (epsilon)/(depsisq))
  sigx3 <- 0.5 * ((1 - wvsq) * exp(Wu)/sigmastar) + sigmastar
  sigx4 <- 0.5 * ((S^2 * (epsilon)^2/(sigma_sq) - 1)/(depsisq) -
    S^2 * depsi * (sigma_sq) * (epsilon)^2/(depsisq)^2) -
    0.5/(depsisq)
  sigx5 <- 0.5 * (S^2 * (epsilon)^2/(depsi * (sigma_sq)^4)) -
    (0.5 * (S^2 * depsi * (epsilon)^2) + 2 * (depsi * (sigma_sq)))/(depsisq)^2
  sigx5epsi <- S * (sigx5) * depsi * (epsilon)
  wusq <- exp(Wu)/(sigma_sq)
  wuwvsq <- (1 - wusq) * exp(Wv)
  sig2wu <- 1/(sigx2) - (0.5 * (wuwvsq/sigmastar) + sigmastar) *
    exp(Wu)/(sigx2)^2
  sigx6 <- (sig2wu) * dmusig/pmusig
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = S * (dmusigwu/(pmusig2) +
    S * (epsilon))/(sigma_sq), FUN = "*"), sweep(uHvar, MARGIN = 1,
    STATS = exp(Wu) * (S * (depsisqx2 - sigx6) * (epsilon) -
      0.5/(sigma_sq)), FUN = "*"), sweep(vHvar, MARGIN = 1,
    STATS = exp(Wv) * (S * ((sigx3) * dmusigwu/((sigx2)^2 *
      pmusig) + depsisqx2) * (epsilon) - 0.5/(sigma_sq)),
    FUN = "*")))
  return(gradll)
}

# Hessian of the likelihood function ----------

chesshalfnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  sigma_sq <- exp(Wu) + exp(Wv)
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvsq <- exp(Wv)/(sigma_sq)
  mustar <- -exp(Wu) * S * epsilon/(sigma_sq)
  sigmastar <- sqrt(exp(Wu) * wvsq)
  musig <- mustar/sigmastar
  pmusig <- pnorm(musig)
  dmusig <- dnorm(musig)
  sigx2 <- (sigma_sq) * sigmastar
  dmusigwu <- dmusig * exp(Wu)
  pmusig2 <- pmusig * sigmastar
  depsi <- dnorm(S * (epsilon)/sqrt(sigma_sq))
  depsisq <- depsi * (sigma_sq)^2
  depsisqx2 <- 0.5 * (S * depsi * (epsilon)/(depsisq))
  sigx3 <- 0.5 * ((1 - wvsq) * exp(Wu)/sigmastar) + sigmastar
  sigx4 <- 0.5 * ((S^2 * (epsilon)^2/(sigma_sq) - 1)/(depsisq) -
    S^2 * depsi * (sigma_sq) * (epsilon)^2/(depsisq)^2) -
    0.5/(depsisq)
  sigx5 <- 0.5 * (S^2 * (epsilon)^2/(depsi * (sigma_sq)^4)) -
    (0.5 * (S^2 * depsi * (epsilon)^2) + 2 * (depsi * (sigma_sq)))/(depsisq)^2
  sigx5epsi <- S * (sigx5) * depsi * (epsilon)
  wusq <- exp(Wu)/(sigma_sq)
  wuwvsq <- (1 - wusq) * exp(Wv)
  sig2wu <- 1/(sigx2) - (0.5 * (wuwvsq/sigmastar) + sigmastar) *
    exp(Wu)/(sigx2)^2
  sigx6 <- (sig2wu) * dmusig/pmusig
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * (dmusigwu * exp(Wu) * (S * (epsilon)/(exp(Wv) *
      pmusig2) - dmusig/(pmusig2)^2)/(sigma_sq) - 1)/(sigma_sq),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (sigx6 + S * ((sigx4) * depsi -
      (sig2wu) * dmusigwu * (S * (epsilon)/exp(Wv) - dmusig/(pmusig2))/((sigma_sq) *
        pmusig)) * (epsilon)) * exp(Wu), FUN = "*"),
    uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * exp(Wv) * (S * ((sigx3) * dmusigwu *
      exp(Wu) * (S * (epsilon)/((sigx2)^2 * exp(Wv) * pmusig) -
      (sigx2)^2 * dmusig/(((sigx2)^2 * pmusig)^2 * sigmastar))/(sigma_sq) +
      (sigx4) * depsi) * (epsilon) - (sigx3) * dmusigwu/((sigx2)^2 *
      pmusig)), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = ((0.5/(sigma_sq)^2 +
    S * (0.5 * (sigx5epsi) - dmusig * (S * (sig2wu)^2 * (dmusig/pmusig -
      S * exp(Wu) * (epsilon)/(sigx2)) * (epsilon) - ((0.5 *
      (wusq) + 1 - 0.5 * (0.5 * (1 - wusq) + wusq)) * wuwvsq/sigmastar +
      (2 - 2 * ((0.5 * (wuwvsq/sigmastar) + sigmastar)^2 *
        exp(Wu) * (sigma_sq)/(sigx2)^2)) * sigmastar)/(sigx2)^2)/pmusig) *
      (epsilon)) * exp(Wu) + S * (depsisqx2 - sigx6) *
    (epsilon) - 0.5/(sigma_sq)) * exp(Wu), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = ((0.5 * (wvsq) - 0.5)/(sigma_sq) +
      S * ((((0.5 * (wvsq) - 0.5 * (0.5 * (1 - wvsq) +
        wvsq)) * (1 - wvsq) + S^2 * (sigx3)^2 * exp(Wu) *
        exp(Wv) * (epsilon)^2/((sigx2)^2 * (sigma_sq))) *
        exp(Wu)/((sigx2)^2 * pmusig2) + (sigx3) * (1/((sigx2)^2 *
        pmusig) - (sigx3) * (2 * ((sigma_sq) * pmusig2) +
        S * dmusigwu * (epsilon)) * exp(Wv)/((sigx2)^2 *
        pmusig)^2)) * dmusigwu + S * (0.5 * ((sigx5) *
        exp(Wv)) + 0.5/(depsisq)) * depsi * (epsilon)) *
        (epsilon)) * exp(Wv), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = (0.5/(sigma_sq)^2 + S * (((((0.5 * (wuwvsq) -
      S^2 * (sigx3) * (sig2wu) * exp(Wu) * (epsilon)^2)/(sigma_sq) +
      0.5 * ((wusq - 1) * wvsq + 1 - 0.5 * ((1 - wusq) *
        (1 - wvsq))) + 0.5 * (1 - wvsq)) * exp(Wu)/sigmastar +
      sigmastar)/((sigx2)^2 * pmusig) - (sigx3) * (2 *
      ((0.5 * (wuwvsq/sigmastar) + sigmastar) * (sigma_sq) *
        pmusig2) - S * (sigx2)^2 * (sig2wu) * dmusig *
      (epsilon)) * exp(Wu)/((sigx2)^2 * pmusig)^2) * dmusig +
      0.5 * (sigx5epsi)) * (epsilon)) * exp(Wu) * exp(Wv),
    FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

halfnormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csthalfnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(chalfnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = chalfnormlike,
    grad = cgradhalfnormlike, hess = chesshalfnormlike, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S), sr1 = trust.optim(x = startVal,
    fn = function(parm) -sum(chalfnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(chalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chesshalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(chalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradhalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chesshalfnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(chalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgradhalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chesshalfnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradhalfnormlike(mleObj$par,
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
      mleObj$hessian <- chesshalfnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chesshalfnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- chalfnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradhalfnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

chalfnormeff <- function(object, level) {
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
  u <- mustar + sigmastar * dnorm(mustar/sigmastar)/pnorm(mustar/sigmastar)
  uLB <- mustar + qnorm(1 - (1 - (1 - level)/2) * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  uUB <- mustar + qnorm(1 - (1 - level)/2 * (1 - pnorm(-mustar/sigmastar))) *
    sigmastar
  m <- ifelse(mustar > 0, mustar, 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)  ## for cost it is Farrell equivalent
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

cmarghalfnorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu/2) * dnorm(0), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmarghalfnorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (1 - (dnorm(0)/pnorm(0))^2), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
