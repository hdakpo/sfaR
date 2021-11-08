########################################################
#                                                      #
# Weibull + normal distributions                       #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

cweibullnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  if (k < 0)
    return(NA)
  ll <- numeric(N)
  for (i in 1:N) {
    ur <- exp(Wu[i]/2) * (-log(1 - FiMat[i, ]))^(1/k)
    ll[i] <- log(mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv[i]/2))))
  }
  return(ll)
}

# starting value for the log-likelihood ----------

cstweibullnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/2)))^(2/3)
  } else {
    varu <- (-S * m3/2)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1/2 * log((epsiRes^2 - varv)^2)
  dep_v <- 1/2 * log((epsiRes^2 - varu)^2)
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
    beta <- c(olsObj[1] + S * sqrt(varu), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi, k = 1))
}

# Gradient of the likelihood function ----------

cgradweibullnormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  k <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  gradll <- matrix(nrow = N, ncol = nXvar + nuZUvar + nvZVvar +
    1)
  lFimat <- (-log(1 - FiMat))^(1/k)
  lFiu <- sweep(S * lFimat, MARGIN = 1, STATS = exp(Wu/2),
    FUN = "*")
  lFiuepsi <- sweep(lFiu, MARGIN = 1, STATS = epsilon, FUN = "+")
  dFimat <- dnorm(sweep(lFiuepsi, MARGIN = 1, STATS = 1/exp(Wv/2),
    FUN = "*"))
  dFiv <- sweep(dFimat, MARGIN = 1, STATS = 1/exp(Wv/2), FUN = "*")
  lFi1 <- dFimat * lFiuepsi
  lFi2 <- log(-log(1 - FiMat))
  lFi3 <- lFimat * lFi1
  sigx1 <- sweep(lFi1, MARGIN = 1, STATS = 1/exp(Wv/2)^3, FUN = "*")
  sigx2 <- sweep(S * lFi2 * lFi3, MARGIN = 1, STATS = exp(Wu/2)/(k^2 *
    exp(Wv/2)^3), FUN = "*")
  sigx3 <- sweep(lFi3, MARGIN = 1, STATS = exp(Wu/2)/exp(Wv/2)^3,
    FUN = "*")
  sigx4 <- sweep(sweep(0.5 * lFi1 * lFiuepsi, MARGIN = 1, STATS = 1/exp(Wv/2)^2,
    FUN = "*") - 0.5 * dFimat, MARGIN = 1, STATS = 1/exp(Wv/2),
    FUN = "*")
  sdFiv <- apply(dFiv, 1, sum)
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(sigx1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = -(0.5 *
      (S * uHvar[, k])), FUN = "*"), 1, sum)/sdFiv
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(sigx4, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sdFiv
  }
  gradll <- cbind(gx, gu, gv, apply(sigx2, 1, sum)/sdFiv)
  return(gradll)
}

# Optimization using different algorithms ----------

weibullnormAlgOpt <- function(start, olsParam, dataTable, S,
  nXvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar,
  method, printInfo, itermax, stepmax, tol, gradtol, hessianType,
  qac) {
  startVal <- if (!is.null(start))
    start else cstweibullnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cweibullnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
                         bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
                         nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cweibullnormlike,
    grad = cgradweibullnormlike, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat), sr1 = trust.optim(x = startVal,
    fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), hs = function(parm) as(jacobian(function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), parm), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cweibullnormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), gradient = function(parm) -colSums(cgradweibullnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), control = list(iter.max = itermax,
      trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
      x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradweibullnormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat))
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradweibullnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), mleObj$solution)
  }
  mleObj$logL_OBS <- cweibullnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat)
  mleObj$gradL_OBS <- cgradweibullnormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# average efficiency (BC style) evaluation ----------

fnExpUWeiNorm <- function(u, sigma, k) {
  exp(-u) * k/sigma * (u/sigma)^(k - 1) * exp(-(u/sigma)^k)
}

# integral to solve for conditional efficiencies ----------

fnCondEffWeibull <- function(u, sigmaU, sigmaV, k, epsilon, S) {
  u * k/(sigmaU * sigmaV) * (u/sigmaU)^(k - 1) * exp(-(u/sigmaU)^k) *
    dnorm((epsilon + S * u)/sigmaV)
}

# Conditional efficiencies estimation ----------

cweibullnormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
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
  u <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(Wu[i]/2) * (-log(1 - object$FiMat[i, ]))^(1/k)
    density_epsilon <- mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2)))
    u[i] <- integrate(f = fnCondEffWeibull, lower = 0, upper = Inf,
      sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2), k = k,
      epsilon = epsilon[i], S = object$S, rel.tol = 1e-10,
      stop.on.error = FALSE)$value/density_epsilon
  }
  if (object$logDepVar == TRUE){
    teJLMS <- exp(-u)
    res <- bind_cols(u = u, teJLMS = teJLMS)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargweibull_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 1/2,
    nrow = 1), matrix(exp(Wu/2) * gamma(1 + 1/k), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmargweibull_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  k <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * (gamma(1 + 2/k) - (gamma(1 + 1/k))^2),
      ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
