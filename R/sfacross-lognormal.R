########################################################
#                                                      #
# Lognormal + normal distributions                     #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

clognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ll <- numeric(N)
  for (i in 1:N) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(FiMat[i, ]))
    ll[i] <- log(mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      S * ur)/exp(Wv[i]/2))))
  }
  return(ll)
}

# starting value for the log-likelihood ----------

cstlognorm <- function(olsObj, epsiRes, S, nmuZUvar, nuZUvar, uHvar,
  muHvar, nvZVvar, vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  varu <- tryCatch((nleqslv(x = 0.01, fn = function(x) -exp(9 *
    x^2/2) + 3 * exp(5 * x^2/2) - 2 * exp(3 * x^2/2) - S *
    m3, method = "Newton")$x)^2, error = function(e) e)
  if (inherits(varu, "error"))
    varu <- 0.01
  varv <- if ((m2 - exp(varu) * (exp(varu) - 1)) < 0) {
    abs(m2 - exp(varu) * (exp(varu) - 1))
  } else {
    (m2 - exp(varu) * (exp(varu) - 1))
  }
  dep_u <- 1/2 * log((log(1/2 + sqrt(4 * ((epsiRes^2 - varv)^2)^(1/2))/2))^2)
  dep_v <- 1/2 * log((epsiRes^2 - exp(varu) * (exp(varu) -
    1))^2)
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
    beta <- beta <- c(olsObj[1] + S * exp(varu/2), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, omega, delta, phi))
}

# Gradient of the likelihood function ----------

cgradlognormlike <- function(parm, nXvar, nmuZUvar, nuZUvar, nvZVvar,
  muHvar, uHvar, vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  omega <- parm[(nXvar + 1):(nXvar + nmuZUvar)]
  delta <- parm[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar + nuZUvar)]
  phi <- parm[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar + nmuZUvar +
    nuZUvar + nvZVvar)]
  mu <- as.numeric(crossprod(matrix(omega), t(muHvar)))
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  qFimat <- qnorm(FiMat)
  WuqFi <- sweep(qFimat, MARGIN = 1, STATS = exp(Wu/2), FUN = "*")
  WumuqFi <- sweep(WuqFi, MARGIN = 1, STATS = mu, FUN = "+")
  WumuqFiepsi <- sweep(S * exp(WumuqFi), MARGIN = 1, STATS = epsilon,
    FUN = "+")
  WuWvmuqFiepsi <- sweep(WumuqFiepsi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/")
  dqFi <- dnorm(WuWvmuqFiepsi)
  WvdqFi <- apply(sweep(dqFi, MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/"), 1, sum)
  sigx1 <- sweep(dqFi * (WumuqFiepsi), MARGIN = 1, STATS = exp(3 *
    Wv/2), FUN = "/")
  sigx2 <- sweep(dqFi * exp(WumuqFi) * (WumuqFiepsi), MARGIN = 1,
    STATS = exp(3 * Wv/2), FUN = "/")
  sigx3 <- sweep(dqFi * exp(WumuqFi) * qFimat * (WumuqFiepsi),
    MARGIN = 1, STATS = exp(Wu/2)/exp(3 * Wv/2), FUN = "*")
  sigx4 <- sweep(0.5 * (dqFi * (WumuqFiepsi)^2), MARGIN = 1,
    STATS = exp(Wv), FUN = "/")
  sigx5 <- sweep((sigx4 - 0.5 * dqFi), MARGIN = 1, STATS = exp(Wv/2),
    FUN = "/")
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(sigx1, MARGIN = 1, STATS = Xvar[,
      k], FUN = "*"), 1, sum)/WvdqFi
  }
  gmu <- matrix(nrow = N, ncol = nmuZUvar)
  for (k in 1:nmuZUvar) {
    gmu[, k] <- apply(sweep(sigx2, MARGIN = 1, STATS = -(S *
      muHvar[, k]), FUN = "*"), 1, sum)/WvdqFi
  }
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(sigx3, MARGIN = 1, STATS = -(0.5 *
      (S * uHvar[, k])), FUN = "*"), 1, sum)/WvdqFi
  }
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(sigx5, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/WvdqFi
  }
  gradll <- cbind(gx, gmu, gu, gv)
  return(gradll)
}

# Optimization using different algorithms ----------

lognormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  muHvar, nmuZUvar, N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar,
  Yvar, Xvar, method, printInfo, itermax, stepmax, tol, gradtol,
  hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstlognorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar, nmuZUvar = nmuZUvar, muHvar = muHvar)
  startLoglik <- sum(clognormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradlognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = clognormlike,
    grad = cgradlognormlike, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat), sr1 = trust.optim(x = startVal,
    fn = function(parm) -sum(clognormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
      Xvar = Xvar, S = S, N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradlognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(clognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), gr = function(parm) -colSums(cgradlognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), hs = function(parm) as(jacobian(function(parm) -colSums(cgradlognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), parm), "dgCMatrix"), method = "Sparse",
      control = list(maxit = itermax, cgtol = gradtol,
        stop.trust.radius = tol, prec = tol, report.level = if (printInfo) 2 else 0,
        report.precision = 1L, preconditioner = 1L)),
    mla = mla(b = startVal, fn = function(parm) -sum(clognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), gr = function(parm) -colSums(cgradlognormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(clognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), gradient = function(parm) -colSums(cgradlognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), control = list(iter.max = itermax,
        trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
        x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradlognormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat))
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradlognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradlognormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        nmuZUvar = nmuZUvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), mleObj$solution)
  }
  mleObj$logL_OBS <- clognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat)
  mleObj$gradL_OBS <- cgradlognormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, nmuZUvar = nmuZUvar,
    muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
    Xvar = Xvar, S = S, N = N, FiMat = FiMat)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# average efficiency (BC style) evaluation ----------

fnExpULogNorm <- function(u, sigma, mu) {
  1/(u * sigma * sqrt(2 * pi)) * exp(-(log(u) - mu)^2/(2 *
    sigma^2) - u)
}

# integral to solve for conditional efficiencies ----------

fnCondEffLogNorm <- function(u, sigmaU, sigmaV, mu, epsilon,
  S) {
  1/(sigmaU * sigmaV) * dnorm((log(u) - mu)/sigmaU) * dnorm((epsilon +
    S * u)/sigmaV)
}

# Conditional efficiencies estimation ----------

clognormeff <- function(object, level) {
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
  u <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    ur <- exp(mu[i] + exp(Wu[i]/2) * qnorm(object$FiMat[i,
      ]))
    density_epsilon <- (mean(1/exp(Wv[i]/2) * dnorm((epsilon[i] +
      object$S * ur)/exp(Wv[i]/2))))
    u[i] <- integrate(f = fnCondEffLogNorm, lower = 0, upper = Inf,
      sigmaU = exp(Wu[i]/2), sigmaV = exp(Wv[i]/2), mu = mu[i],
      epsilon = epsilon[i], S = object$S)$value/density_epsilon
  }
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    res <- bind_cols(u = u, teJLMS = teJLMS)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmarglognorm_Eu <- function(object) {
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(mu + exp(Wu)/2 + Wu)/2, ncol = 1))
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

cmarglognorm_Vu <- function(object) {
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
  mu_mat <- kronecker(matrix(omega[2:object$nmuZUvar], nrow = 1),
    matrix(2 * (exp(Wu) - 1) * exp(2 * mu + exp(Wu)), ncol = 1))
  Wu_mat <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu) * exp(2 * mu + exp(Wu) + Wu), ncol = 1))
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
