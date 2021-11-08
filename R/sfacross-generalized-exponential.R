########################################################
#                                                      #
# Generalized-exponential + normal distributions       #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

cgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  ll <- log(2) - 1/2 * Wu + log(exp(A) * pnorm(a) - exp(B) *
    pnorm(b))
  return(ll)
}

# starting value for the log-likelihood ----------

cstgenexponorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar,
  nvZVvar, vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3/9)))^(2/3)
  } else {
    varu <- (-S * m3/9)^(2/3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - 5/4 * varu)
  } else {
    varv <- m2 - 5/4 * varu
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
    beta <- c(olsObj[1] + S * sqrt(varu) * 3/2, olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------

cgradgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvwu <- exp(Wv/2)/exp(Wu/2)
  wvwuepsi <- wvwu + S * (epsilon)/exp(Wv/2)
  wvwuepsix2 <- 2 * (wvwu) + S * (epsilon)/exp(Wv/2)
  a <- -(wvwuepsi)
  pa <- pnorm(a)
  da <- dnorm(a)
  b <- -(wvwuepsix2)
  pb <- pnorm(b)
  db <- dnorm(b)
  eA <- exp(exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2))
  eB <- exp(2 * (exp(Wv)/exp(Wu)) + 2 * (S * (epsilon)/exp(Wu/2)))
  eC <- 2 * (exp(Wv)/exp(Wu)) + S * (epsilon)/exp(Wu/2)
  epsiv <- 0.5 * (S * (epsilon)/exp(Wv/2))
  epsiuv <- 0.5 * (wvwu) - epsiv
  pda <- da/exp(Wv/2) - pa/exp(Wu/2)
  sigx1 <- (pda) * eA - (db/exp(Wv/2) - 2 * (pb/exp(Wu/2))) *
    eB
  pab <- eA * pa - eB * pb
  sigx2 <- 0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (exp(Wu) *
    exp(Wv)/(2 * exp(Wu))^2)
  sigx3 <- (0.5 * (da * wvwu) - (sigx2) * pa) * eA - (db *
    wvwu - (eC) * pb) * eB
  sigx5 <- 2 * (exp(Wv) * pb/exp(Wu)) - db * (wvwu - epsiv)
  sigx4 <- eA * (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) * da) -
    (sigx5) * eB
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (sigx1)/(pab),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((sigx3)/(pab) -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx4)/(pab),
    FUN = "*"))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessgenexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar,
  uHvar, vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvwu <- exp(Wv/2)/exp(Wu/2)
  wvwuepsi <- wvwu + S * (epsilon)/exp(Wv/2)
  wvwuepsix2 <- 2 * (wvwu) + S * (epsilon)/exp(Wv/2)
  a <- -(wvwuepsi)
  pa <- pnorm(a)
  da <- dnorm(a)
  b <- -(wvwuepsix2)
  pb <- pnorm(b)
  db <- dnorm(b)
  eA <- exp(exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2))
  eB <- exp(2 * (exp(Wv)/exp(Wu)) + 2 * (S * (epsilon)/exp(Wu/2)))
  eC <- 2 * (exp(Wv)/exp(Wu)) + S * (epsilon)/exp(Wu/2)
  epsiv <- 0.5 * (S * (epsilon)/exp(Wv/2))
  epsiuv <- 0.5 * (wvwu) - epsiv
  pda <- da/exp(Wv/2) - pa/exp(Wu/2)
  sigx1 <- (pda) * eA - (db/exp(Wv/2) - 2 * (pb/exp(Wu/2))) *
    eB
  pab <- eA * pa - eB * pb
  sigx2 <- 0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (exp(Wu) *
    exp(Wv)/(2 * exp(Wu))^2)
  sigx3 <- (0.5 * (da * wvwu) - (sigx2) * pa) * eA - (db *
    wvwu - (eC) * pb) * eB
  sigx5 <- 2 * (exp(Wv) * pb/exp(Wu)) - db * (wvwu - epsiv)
  sigx4 <- eA * (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) * da) -
    (sigx5) * eB
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * ((((wvwuepsi)/exp(Wv/2) - 1/exp(Wu/2)) *
      da/exp(Wv/2) - (pda)/exp(Wu/2)) * eA - ((((wvwuepsix2)/exp(Wv/2) -
      2/exp(Wu/2)) * db/exp(Wv/2) - 2 * ((db/exp(Wv/2) -
      2 * (pb/exp(Wu/2)))/exp(Wu/2))) * eB + (sigx1)^2/(pab)))/(pab),
    FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * ((((0.5 + sigx2) * pa - 0.5 *
      (da * wvwu))/exp(Wu/2) + (0.5 * ((wvwuepsi)/exp(Wu/2)) -
      (sigx2)/exp(Wv/2)) * da) * eA - ((((wvwuepsix2)/exp(Wu/2) -
      (eC)/exp(Wv/2)) * db + (pb - 2 * (db * wvwu - (eC) *
      pb))/exp(Wu/2)) * eB + (sigx3) * (sigx1)/(pab)))/(pab),
    FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * ((da * (exp(Wv)/(2 * exp(Wu)) -
      ((epsiuv) * (wvwuepsi) + 0.5))/exp(Wv/2) - (exp(Wv) *
      pa/(2 * exp(Wu)) - (epsiuv) * da)/exp(Wu/2)) * eA -
      (((2 * (exp(Wv)/exp(Wu)) - ((wvwuepsix2) * (wvwu -
        epsiv) + 0.5)) * db/exp(Wv/2) - 2 * ((sigx5)/exp(Wu/2))) *
        eB + (sigx1) * (sigx4)/(pab)))/(pab), FUN = "*"),
    vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = (((0.5 *
    (0.5 * (exp(Wv/2) * (wvwuepsi)/exp(Wu/2)) - 0.5) - 0.5 *
    (sigx2)) * da * wvwu - ((0.5 * (da * wvwu) - (sigx2) *
    pa) * (sigx2) + (2 * ((1 - 8 * (exp(Wu)^2/(2 * exp(Wu))^2)) *
    exp(Wu) * exp(Wv)/(2 * exp(Wu))^2) - 0.25 * (S * (epsilon)/exp(Wu/2))) *
    pa)) * eA - (((((wvwuepsix2) * exp(Wv/2) - S * (epsilon))/exp(Wu/2) -
    (0.5 + 2 * (exp(Wv)/exp(Wu)))) * db * wvwu + (0.5 * (S *
    (epsilon)/exp(Wu/2)) + 2 * (exp(Wv)/exp(Wu))) * pb -
    (eC) * (db * wvwu - (eC) * pb)) * eB + (sigx3)^2/(pab)))/(pab),
    FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = ((((exp(Wv) * pa/(2 * exp(Wu)) -
      (epsiuv) * da)/2 + (pa - (epsiuv) * da)/2) * exp(Wv)/exp(Wu) -
      (0.25 * (wvwu) + 0.25 * (S * (epsilon)/exp(Wv/2)) -
        (epsiuv)^2 * (wvwuepsi)) * da) * eA - (((2 *
      (sigx5) + 2 * (pb - db * (wvwu - epsiv))) * exp(Wv)/exp(Wu) -
      (0.25 * (S * (epsilon)/exp(Wv/2)) + 0.5 * (wvwu) +
        b * (wvwu - epsiv)^2) * db) * eB + (sigx4)^2/(pab)))/(pab),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = (((da * exp(Wv/2)/(4 * (exp(Wu) * exp(Wu/2))) -
      2 * (exp(Wu) * pa/(2 * exp(Wu))^2)) * exp(Wv) - ((0.5 *
      ((epsiuv) * (wvwuepsi)) - 0.25) * da * wvwu + (sigx2) *
      (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) * da))) *
      eA - ((sigx3) * (sigx4)/(pab) + (2 * ((db * wvwu -
      pb) * exp(Wv)/exp(Wu)) - (((wvwuepsix2) * (wvwu -
      epsiv) - 0.5) * db * wvwu + (sigx5) * (eC))) * eB))/(pab),
    FUN = "*"), vHvar)
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

genexponormAlgOpt <- function(start, olsParam, dataTable, S,
  nXvar, uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method,
  printInfo, itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else cstgenexponorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(cgenexponormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hessian = 0, control = list(trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol)), maxLikAlgo = maxRoutine(fn = cgenexponormlike,
    grad = cgradgenexponormlike, hess = chessgenexponormlike,
    start = startVal, finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S), sr1 = trust.optim(x = startVal,
    fn = function(parm) -sum(cgenexponormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chessgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(cgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chessgenexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(cgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgradgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chessgenexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradgenexponormlike(mleObj$par,
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
      mleObj$hessian <- chessgenexponormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chessgenexponormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- cgenexponormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradgenexponormlike(parm = mlParam,
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

cgenexponormeff <- function(object, level) {
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
  A <- object$S * epsilon/exp(Wu/2) + exp(Wv)/(2 * exp(Wu))
  B <- 2 * object$S * epsilon/exp(Wu/2) + 2 * exp(Wv)/exp(Wu)
  a <- -object$S * epsilon/exp(Wv/2) - exp(Wv/2)/exp(Wu/2)
  b <- -object$S * epsilon/exp(Wv/2) - 2 * exp(Wv/2)/exp(Wu/2)
  u <- exp(Wv/2) * (exp(A) * (dnorm(a) + a * pnorm(a)) - exp(B) *
    (dnorm(b) + b * pnorm(b)))/(exp(A) * pnorm(a) - exp(B) *
    pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) * pnorm(a -
      exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) + exp(Wv)/2) *
      pnorm(b - exp(Wv/2)))/(exp(A) * pnorm(a) - exp(B) *
      pnorm(b))
    res <- bind_cols(u = u, teJLMS = teJLMS, teBC = teBC)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmarggenexponorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 3/4,
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmarggenexponorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 5/4,
    nrow = 1), matrix(exp(Wu), ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
