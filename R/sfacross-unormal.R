########################################################
#                                                      #
# uniform + normal distributions                       #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

cuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ll <- -Wu + log(pnorm((exp(Wu) + S * epsilon)/exp(Wv/2)) -
    pnorm(S * epsilon/exp(Wv/2)))
  return(ll)
}

# starting value for the log-likelihood ----------

cstuninorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
  vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m4 <- moment(epsiRes, order = 4)
  if ((m2^2 - m4) < 0) {
    theta <- (abs(120 * (3 * m2^2 - m4)))^(1/4)
    varu <- theta^2/12
  } else {
    theta <- (120 * (3 * m2^2 - m4))^(1/4)
    varu <- theta^2/12
  }
  if ((m2 - varu) < 0) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1/2 * log(((epsiRes^2 - varv) * 12)^2) * 1/2
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
    beta <- c(olsObj[1] + S * theta/2, olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------

cgraduninormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsiv <- S * (epsilon)/exp(Wv/2)
  epsiu <- exp(Wu) + S * (epsilon)
  epsiuv <- (epsiu)/exp(Wv/2)
  depsiv <- dnorm(epsiv)
  depsiuv <- dnorm(epsiuv)
  pepsiv <- pnorm(epsiv)
  pepsiuv <- pnorm(epsiuv)
  sigx1 <- 0.5 * (S * depsiv * (epsilon)) - 0.5 * (depsiuv *
    (epsiu))
  sigx2 <- depsiv - depsiuv
  sigx3 <- pepsiuv - pepsiv
  sigx4 <- exp(Wv/2) * (sigx3)
  depsiuvx2 <- depsiuv * exp(Wu)
  gradll <- (cbind(sweep(Xvar, MARGIN = 1, STATS = S * (sigx2)/(sigx4),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = (depsiuvx2/(sigx4) -
    1), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx1)/(sigx4),
    FUN = "*")))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessuninormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  epsiv <- S * (epsilon)/exp(Wv/2)
  epsiu <- exp(Wu) + S * (epsilon)
  epsiuv <- (epsiu)/exp(Wv/2)
  depsiv <- dnorm(epsiv)
  depsiuv <- dnorm(epsiuv)
  pepsiv <- pnorm(epsiv)
  pepsiuv <- pnorm(epsiuv)
  sigx1 <- 0.5 * (S * depsiv * (epsilon)) - 0.5 * (depsiuv *
    (epsiu))
  sigx2 <- depsiv - depsiuv
  sigx3 <- pepsiuv - pepsiv
  sigx4 <- exp(Wv/2) * (sigx3)
  depsiuvx2 <- depsiuv * exp(Wu)
  sigx5 <- exp(Wv/2)^3 * (sigx3)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * ((S * depsiv * (epsilon) - depsiuv * (epsiu))/(sigx5) -
      (sigx2)^2/(sigx4)^2), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * ((epsiu)/(sigx5) - (sigx2)/(sigx4)^2) *
      depsiuvx2, FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * ((0.5 * (depsiv * (S^2 * (epsilon)^2/exp(Wv/2)^2 -
      1)) - 0.5 * (((epsiu)^2/exp(Wv/2)^2 - 1) * depsiuv))/(sigx4) -
      (sigx1) * (sigx2)/(sigx4)^2), FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = ((1 -
    exp(Wu) * epsiuv/exp(Wv/2))/(sigx4) - depsiuvx2/(sigx4)^2) *
    depsiuvx2, FUN = "*"), uHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = -(((sigx1)/(sigx4)^2 + 0.5 * ((1 - epsiuv^2)/(sigx4))) *
      depsiuvx2), FUN = "*"), vHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = ((0.25 * (S^3 * depsiv * (epsilon)^3) -
      0.25 * (depsiuv * (epsiu)^3))/(sigx5) - (0.5 * (sigx4) +
      sigx1) * (sigx1)/(sigx4)^2), FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

uninormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstuninorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cuninormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cuninormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgraduninormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = cuninormlike, grad = cgraduninormlike,
      hess = chessuninormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(cuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgraduninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgraduninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chessuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(cuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgraduninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chessuninormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgraduninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chessuninormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgraduninormlike(mleObj$par,
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
      mleObj$hessian <- chessuninormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chessuninormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- cuninormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgraduninormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

cuninormeff <- function(object, level) {
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
  theta <- exp(Wu)
  u1 <- -exp(Wv/2) * ((dnorm((theta + object$S * epsilon)/exp(Wv/2)) -
    dnorm(object$S * epsilon/exp(Wv/2)))/(pnorm((theta +
    object$S * epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))) -
    object$S * epsilon
  u2 <- exp(Wv/2) * (dnorm(object$S * epsilon/exp(Wv/2))/(1 -
    pnorm(object$S * epsilon/exp(Wv/2))) - object$S * epsilon/exp(Wv/2))  # when theta/sigmav ---> Infty
  uLB <- exp(Wv/2) * qnorm((1 - level)/2 * pnorm((theta + object$S *
    epsilon)/exp(Wv/2)) + (1 - (1 - level)/2) * pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon
  uUB <- exp(Wv/2) * qnorm((1 - (1 - level)/2) * pnorm((theta +
    object$S * epsilon)/exp(Wv/2)) + (1 - level)/2 * pnorm(object$S *
    epsilon/exp(Wv/2))) - object$S * epsilon
  m <- ifelse(-theta < object$S * epsilon & object$S * epsilon <
    0, -object$S * epsilon, ifelse(object$S * epsilon >=
    0, 0, theta))
  if (object$logDepVar == TRUE) {
    teJLMS1 <- exp(-u1)
    teJLMS2 <- exp(-u2)
    teMO <- exp(-m)
    teBC1 <- exp(object$S * epsilon + exp(Wv)/2) * (pnorm((object$S *
      epsilon + theta)/exp(Wv/2) + exp(Wv/2)) - pnorm(object$S *
      epsilon/exp(Wv/2) + exp(Wv/2)))/(pnorm((theta + object$S *
      epsilon)/exp(Wv/2)) - pnorm(object$S * epsilon/exp(Wv/2)))
    teBC2 <- exp(object$S * epsilon + exp(Wv)/2) * (1 - pnorm(object$S *
      epsilon/exp(Wv/2) + exp(Wv/2)))/(1 - pnorm(object$S *
      epsilon/exp(Wv/2)))
    teBCLB <- exp(-uUB)
    teBCUB <- exp(-uLB)
    res <- bind_cols(u1 = u1, u2 = u2, uLB = uLB, uUB = uUB,
      teJLMS1 = teJLMS1, teJLMS2 = teJLMS2, m = m, teMO = teMO,
      teBC1 = teBC1, teBC2 = teBC2, teBCLB = teBCLB, teBCUB = teBCUB)
  } else {
    res <- bind_cols(u1 = u1, u2 = u2, uLB = uLB, uUB = uUB,
      m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmarguninorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(sqrt(3)/2 * exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmarguninorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
