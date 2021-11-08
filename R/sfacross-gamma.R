########################################################
#                                                      #
# Gamma + normal distributions                         #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

cgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mui <- -S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi <- numeric(N)
  for (i in 1:N) {
    Hi[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(FiMat[i,
      ] + (1 - FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  if (P <= 0)
    return(NA)
  ll <- -1/2 * P * Wu - log(gamma(P)) + exp(Wv)/(2 * exp(Wu)) +
    S * epsilon/sqrt(exp(Wu)) + pnorm(-S * epsilon/sqrt(exp(Wv)) -
    sqrt(exp(Wv)/exp(Wu)), log.p = TRUE) + log(Hi)
  return(ll)
}

# starting value for the log-likelihood ----------

cstgammanorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
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
  return(c(beta, delta, phi, P = 1))
}

# Gradient of the likelihood function ----------

cgradgammanormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S, N, FiMat) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  P <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wuwv <- sqrt(exp(Wv)/exp(Wu))
  dwuwv <- dnorm(-(S * (epsilon)/exp(Wv/2) + wuwv))
  pwuwv <- pnorm(-(S * (epsilon)/exp(Wv/2) + wuwv))
  depsi <- dnorm((exp(Wv)/exp(Wu/2) + S * (epsilon))/exp(Wv/2))
  pepsi <- pnorm((exp(Wv)/exp(Wu/2) + S * (epsilon))/exp(Wv/2))
  sigx1 <- (exp(Wv)/exp(Wu/2) + S * (epsilon))
  sigx2 <- (dwuwv/(exp(Wv/2) * pwuwv) - 1/exp(Wu/2))
  sigx3 <- (dwuwv/(exp(Wu) * pwuwv * wuwv))
  sigx4 <- (0.5 * sigx3 - 2 * (exp(Wu)/(2 * exp(Wu))^2)) *
    exp(Wv)
  sigx5 <- (sigx4 - (0.5 * (S * (epsilon)/exp(Wu/2)) + 0.5 *
    P))
  sigx6 <- (exp(Wv)/(2 * exp(Wu)) - (0.5 * (exp(Wv)/(exp(Wu) *
    wuwv)) - 0.5 * (S * (epsilon)/exp(Wv/2))) * dwuwv/pwuwv)
  F1 <- sweep((1 - FiMat), MARGIN = 1, STATS = pepsi, FUN = "*") +
    FiMat
  dqF1 <- dnorm(qnorm(F1))
  F2 <- sweep(qnorm(F1), MARGIN = 1, STATS = exp(Wv/2), FUN = "*")
  F3 <- sweep(F2, MARGIN = 1, STATS = sigx1, FUN = "-")
  sumF3 <- apply(F3^(P - 1), 1, sum)
  F4 <- sweep((1 - FiMat)/dqF1, MARGIN = 1, STATS = depsi,
    FUN = "*")
  F5 <- sweep((0.5 - 0.5 * (F4)) * (F3)^(P - 2) * (P - 1),
    MARGIN = 1, STATS = exp(Wv)/exp(Wu/2), FUN = "*")
  F6 <- sweep(F4, MARGIN = 1, STATS = exp(Wv)/exp(Wu/2) - 0.5 *
    sigx1, FUN = "*") + 0.5 * (F2)
  F7 <- sweep(F6, MARGIN = 1, STATS = exp(Wv)/exp(Wu/2), FUN = "-") *
    F3^(P - 2) * (P - 1)
  gx <- matrix(nrow = N, ncol = nXvar)
  for (k in 1:nXvar) {
    gx[, k] <- apply(sweep(S * (1 - F4) * F3^(P - 2) * (P -
      1), MARGIN = 1, STATS = Xvar[, k], FUN = "*"), 1,
      sum)/sumF3
  }
  gx <- sweep(Xvar, MARGIN = 1, STATS = S * sigx2, FUN = "*") +
    gx
  gu <- matrix(nrow = N, ncol = nuZUvar)
  for (k in 1:nuZUvar) {
    gu[, k] <- apply(sweep(F5, MARGIN = 1, STATS = uHvar[,
      k], FUN = "*"), 1, sum)/sumF3
  }
  gu <- sweep(uHvar, MARGIN = 1, STATS = sigx5, FUN = "*") +
    gu
  gv <- matrix(nrow = N, ncol = nvZVvar)
  for (k in 1:nvZVvar) {
    gv[, k] <- apply(sweep(F7, MARGIN = 1, STATS = vHvar[,
      k], FUN = "*"), 1, sum)/sumF3
  }
  gv <- sweep(vHvar, MARGIN = 1, STATS = sigx6, FUN = "*") +
    gv
  gradll <- cbind(gx, gu, gv, apply((F3)^(P - 1) * log(F3),
    1, sum)/sumF3 - (0.5 * (Wu) + digamma(P)))
  return(gradll)
}

# Optimization using different algorithms ----------

gammanormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  N, FiMat, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstgammanorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cgammanormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
      FiMat = FiMat)), gr = function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cgammanormlike,
    grad = cgradgammanormlike, start = startVal, finalHessian = if (hessianType ==
      2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cgammanormlike(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat)), gr = function(parm) -colSums(cgradgammanormlike(parm,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S, N = N,
    FiMat = FiMat)), method = "SR1", control = list(maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), hs = function(parm) as(jacobian(function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), parm), "dgCMatrix"),
      method = "Sparse", control = list(maxit = itermax,
        cgtol = gradtol, stop.trust.radius = tol, prec = tol,
        report.level = if (printInfo) 2 else 0, report.precision = 1L,
        preconditioner = 1L)), mla = mla(b = startVal,
      fn = function(parm) -sum(cgammanormlike(parm, nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S,
        N = N, FiMat = FiMat)), gr = function(parm) -colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), print.info = printInfo,
      maxiter = itermax, epsa = gradtol, epsb = gradtol),
    nlminb = nlminb(start = startVal, objective = function(parm) -sum(cgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), gradient = function(parm) -colSums(cgradgammanormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S, N = N, FiMat = FiMat)), control = list(iter.max = itermax,
      trace = if (printInfo) 1 else 0, eval.max = itermax, rel.tol = tol,
      x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradgammanormlike(mleObj$par,
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
      mleObj$hessian <- jacobian(function(parm) colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), mleObj$par)
    if (method == "sr1")
      mleObj$hessian <- jacobian(function(parm) colSums(cgradgammanormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S, N = N, FiMat = FiMat)), mleObj$solution)
  }
  mleObj$logL_OBS <- cgammanormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat)
  mleObj$gradL_OBS <- cgradgammanormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S, N = N, FiMat = FiMat)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

cgammanormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
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
  mui <- -object$S * epsilon - exp(Wv)/sqrt(exp(Wu))
  Hi1 <- numeric(object$Nobs)
  Hi2 <- numeric(object$Nobs)
  for (i in 1:object$Nobs) {
    Hi1[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P))
    Hi2[i] <- mean((mui[i] + sqrt(exp(Wv[i])) * qnorm(object$FiMat[i,
      ] + (1 - object$FiMat[i, ]) * pnorm(-mui[i]/sqrt(exp(Wv[i])))))^(P -
      1))
  }
  u <- Hi1/Hi2
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    res <- bind_cols(u = u, teJLMS = teJLMS)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmarggammanorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P/2 * exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmarggammanorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  P <- object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
    1]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(P * exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
