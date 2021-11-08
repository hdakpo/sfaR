########################################################
#                                                      #
# Truncated skewed laplace + normal distributions      #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

ctslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu)) + S * epsilon/exp(Wu/2)
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) + S * epsilon *
    (1 + lambda)/exp(Wu/2)
  a <- -exp(Wv/2)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - S * epsilon/exp(Wv/2)
  if (lambda < 0) {
    return(NA)
  } else {
    ll <- log(1 + lambda) - log(2 * lambda + 1) - 1/2 * Wu +
      log(2 * exp(A) * pnorm(a) - exp(B) * pnorm(b))
    return(ll)
  }
}

# starting value for the log-likelihood ----------

csttslnorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
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
  return(c(beta, delta, phi, lambda = 0.03))
}

# Gradient of the likelihood function ----------

cgradtslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvwu <- exp(Wv/2)/exp(Wu/2)
  d <- (1 + lambda) * wvwu + S * (epsilon)/exp(Wv/2)
  pa <- pnorm(-(wvwu + S * (epsilon)/exp(Wv/2)))
  da <- dnorm(-(wvwu + S * (epsilon)/exp(Wv/2)))
  pb <- pnorm(-(d))
  db <- dnorm(-(d))
  eB <- exp(((1 + lambda) * exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2)) *
    (1 + lambda))
  eA <- exp(exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2))
  eC <- (1 + lambda) * exp(Wv)/exp(Wu) + S * (epsilon)/exp(Wu/2)
  epsiv <- 0.5 * (S * (epsilon)/exp(Wv/2))
  epsiuv <- 0.5 * (wvwu) - epsiv
  sigx1 <- 0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (exp(Wu) *
    exp(Wv)/(2 * exp(Wu))^2)
  papb <- 2 * (eA * pa) - eB * pb
  epsivl <- 0.5 * ((1 + lambda) * wvwu) - epsiv
  sigx2 <- 2 * (eA * (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) *
    da)) - ((1 + lambda)^2 * exp(Wv) * pb/(2 * exp(Wu)) -
    (epsivl) * db) * eB
  epsiuvl <- (1 + lambda) * exp(Wu) * exp(Wv)/(2 * exp(Wu))^2
  epsiuvlx2 <- (0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (epsiuvl))
  sigx3 <- 2 * ((0.5 * (da * wvwu) - (sigx1) * pa) * eA) -
    (0.5 * (db * wvwu) - epsiuvlx2 * pb) * (1 + lambda) *
      eB
  sigx4 <- db/exp(Wv/2) - (1 + lambda) * pb/exp(Wu/2)
  dpa <- da/exp(Wv/2) - pa/exp(Wu/2)
  sigx5 <- 2 * ((dpa) * eA) - (sigx4) * eB
  pdbuv <- (eC) * pb - db * wvwu
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (sigx5)/(papb),
    FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((sigx3)/(papb) -
    0.5), FUN = "*"), sweep(vHvar, MARGIN = 1, STATS = (sigx2)/(papb),
    FUN = "*"), 1/(1 + lambda) - ((pdbuv) * eB/(papb) + 2/(1 +
    2 * lambda)))
  return(gradll)
}

# Hessian of the likelihood function ----------

chesstslnormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
  vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  lambda <- parm[nXvar + nuZUvar + nvZVvar + 1]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  wvwu <- exp(Wv/2)/exp(Wu/2)
  d <- (1 + lambda) * wvwu + S * (epsilon)/exp(Wv/2)
  pa <- pnorm(-(wvwu + S * (epsilon)/exp(Wv/2)))
  da <- dnorm(-(wvwu + S * (epsilon)/exp(Wv/2)))
  pb <- pnorm(-(d))
  db <- dnorm(-(d))
  eB <- exp(((1 + lambda) * exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2)) *
    (1 + lambda))
  eA <- exp(exp(Wv)/(2 * exp(Wu)) + S * (epsilon)/exp(Wu/2))
  eC <- (1 + lambda) * exp(Wv)/exp(Wu) + S * (epsilon)/exp(Wu/2)
  epsiv <- 0.5 * (S * (epsilon)/exp(Wv/2))
  epsiuv <- 0.5 * (wvwu) - epsiv
  sigx1 <- 0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (exp(Wu) *
    exp(Wv)/(2 * exp(Wu))^2)
  papb <- 2 * (eA * pa) - eB * pb
  epsivl <- 0.5 * ((1 + lambda) * wvwu) - epsiv
  sigx2 <- 2 * (eA * (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) *
    da)) - ((1 + lambda)^2 * exp(Wv) * pb/(2 * exp(Wu)) -
    (epsivl) * db) * eB
  epsiuvl <- (1 + lambda) * exp(Wu) * exp(Wv)/(2 * exp(Wu))^2
  epsiuvlx2 <- (0.5 * (S * (epsilon)/exp(Wu/2)) + 2 * (epsiuvl))
  sigx3 <- 2 * ((0.5 * (da * wvwu) - (sigx1) * pa) * eA) -
    (0.5 * (db * wvwu) - epsiuvlx2 * pb) * (1 + lambda) *
      eB
  sigx4 <- db/exp(Wv/2) - (1 + lambda) * pb/exp(Wu/2)
  dpa <- da/exp(Wv/2) - pa/exp(Wu/2)
  sigx5 <- 2 * ((dpa) * eA) - (sigx4) * eB
  pdbuv <- (eC) * pb - db * wvwu
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar + 1, ncol = nXvar +
    nuZUvar + nvZVvar + 1)
  hessll[1:nXvar, 1:nXvar] <- crossprod(sweep(Xvar, MARGIN = 1,
    STATS = S^2 * (2 * ((((wvwu + S * (epsilon)/exp(Wv/2))/exp(Wv/2) -
      1/exp(Wu/2)) * da/exp(Wv/2) - (dpa)/exp(Wu/2)) *
      eA) - ((((d)/exp(Wv/2) - (1 + lambda)/exp(Wu/2)) *
      db/exp(Wv/2) - (1 + lambda) * (sigx4)/exp(Wu/2)) *
      eB + (sigx5)^2/(papb)))/(papb), FUN = "*"), Xvar)
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (2 * ((((0.5 + sigx1) * pa -
      0.5 * (da * wvwu))/exp(Wu/2) + (0.5 * ((wvwu + S *
      (epsilon)/exp(Wv/2))/exp(Wu/2)) - (sigx1)/exp(Wv/2)) *
      da) * eA) - (((0.5 * ((d)/exp(Wu/2)) - epsiuvlx2/exp(Wv/2)) *
      db + (0.5 * pb - (0.5 * (db * wvwu) - epsiuvlx2 *
      pb) * (1 + lambda))/exp(Wu/2)) * (1 + lambda) * eB +
      (sigx3) * (sigx5)/(papb)))/(papb), FUN = "*"), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (2 * ((da * (exp(Wv)/(2 * exp(Wu)) -
      ((epsiuv) * (wvwu + S * (epsilon)/exp(Wv/2)) + 0.5))/exp(Wv/2) -
      (exp(Wv) * pa/(2 * exp(Wu)) - (epsiuv) * da)/exp(Wu/2)) *
      eA) - ((((1 + lambda)^2 * exp(Wv)/(2 * exp(Wu)) -
      ((d) * (epsivl) + 0.5)) * db/exp(Wv/2) - ((1 + lambda)^2 *
      exp(Wv) * pb/(2 * exp(Wu)) - (epsivl) * db) * (1 +
      lambda)/exp(Wu/2)) * eB + (sigx5) * (sigx2)/(papb)))/(papb),
    FUN = "*"), vHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + nvZVvar + 1)] <- matrix(colSums(sweep(Xvar,
    MARGIN = 1, STATS = -(S * (((eC)/exp(Wv/2) - (d)/exp(Wu/2)) *
      db - (((pdbuv) * (1 + lambda) + pb)/exp(Wu/2) + (pdbuv) *
      (sigx5)/(papb))) * eB/(papb)), FUN = "*")), ncol = 1)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = (2 *
    (((0.5 * (0.5 * (exp(Wv/2) * (wvwu + S * (epsilon)/exp(Wv/2))/exp(Wu/2)) -
      0.5) - 0.5 * (sigx1)) * da * wvwu - ((0.5 * (da *
      wvwu) - (sigx1) * pa) * (sigx1) + (2 * ((1 - 8 *
      (exp(Wu)^2/(2 * exp(Wu))^2)) * exp(Wu) * exp(Wv)/(2 *
      exp(Wu))^2) - 0.25 * (S * (epsilon)/exp(Wu/2))) *
      pa)) * eA) - (((0.5 * (0.5 * ((d) * (1 + lambda) *
    wvwu) - 0.5) - 0.5 * (epsiuvlx2 * (1 + lambda))) * db *
    wvwu - ((0.5 * (db * wvwu) - epsiuvlx2 * pb) * epsiuvlx2 *
    (1 + lambda) + (2 * ((1 - 8 * (exp(Wu)^2/(2 * exp(Wu))^2)) *
    epsiuvl) - 0.25 * (S * (epsilon)/exp(Wu/2))) * pb)) *
    (1 + lambda) * eB + (sigx3)^2/(papb)))/(papb), FUN = "*"),
    uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = (2 * ((((exp(Wv) * pa/(2 * exp(Wu)) -
      (epsiuv) * da)/2 + (pa - (epsiuv) * da)/2) * exp(Wv)/exp(Wu) -
      (0.25 * (wvwu) + 0.25 * (S * (epsilon)/exp(Wv/2)) -
        (epsiuv)^2 * (wvwu + S * (epsilon)/exp(Wv/2))) *
        da) * eA) - (((((1 + lambda)^2 * exp(Wv) * pb/(2 *
      exp(Wu)) - (epsivl) * db)/2 + (pb - (epsivl) * db)/2) *
      (1 + lambda)^2 * exp(Wv)/exp(Wu) - (0.25 * ((1 +
      lambda) * wvwu) + 0.25 * (S * (epsilon)/exp(Wv/2)) -
      (d) * (epsivl)^2) * db) * eB + (sigx2)^2/(papb)))/(papb),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(sweep(uHvar, MARGIN = 1,
    STATS = (2 * ((((0.5 * (da * wvwu) - (sigx1) * pa)/(2 *
      exp(Wu)) - 2 * (exp(Wu) * pa/(2 * exp(Wu))^2)) *
      exp(Wv) + ((epsiuv) * (sigx1) + 0.5 * ((0.5 - (epsiuv) *
      (wvwu + S * (epsilon)/exp(Wv/2))) * wvwu)) * da) *
      eA) - ((((epsivl) * epsiuvlx2 + 0.5 * ((0.5 - (d) *
      (epsivl)) * wvwu)) * db + ((0.5 * (db * wvwu) - epsiuvlx2 *
      pb) * (1 + lambda)/(2 * exp(Wu)) - 2 * (exp(Wu) *
      pb/(2 * exp(Wu))^2)) * (1 + lambda) * exp(Wv)) *
      (1 + lambda) * eB + (sigx3) * (sigx2)/(papb)))/(papb),
    FUN = "*"), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + nvZVvar +
    1)] <- matrix(colSums(sweep(uHvar, MARGIN = 1, STATS = -((((0.5 *
    (eC) - 0.5 * ((d) * wvwu)) * (1 + lambda) + 0.5) * db *
    wvwu - ((pdbuv) * (epsiuvlx2 * (1 + lambda) + (sigx3)/(papb)) +
    ((1 + lambda) * exp(Wv)/exp(Wu) + 0.5 * (S * (epsilon)/exp(Wu/2))) *
      pb)) * eB/(papb)), FUN = "*")), ncol = 1)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + nvZVvar + 1)] <- matrix(colSums(sweep(vHvar,
    MARGIN = 1, STATS = -((((pdbuv) * (1 + lambda)/2 + pb) *
      (1 + lambda) * exp(Wv)/exp(Wu) - (((eC) * (epsivl) +
      (0.5 - (d) * (epsivl)) * wvwu) * db + (pdbuv) * (sigx2)/(papb))) *
      eB/(papb)), FUN = "*")), ncol = 1)
  hessll[(nXvar + nuZUvar + nvZVvar + 1), (nXvar + nuZUvar + nvZVvar +
    1)] <- sum(4/(1 + 2 * lambda)^2 - ((((pdbuv) * eB/(papb) +
    eC) * (pdbuv) + (((d) * exp(Wv/2) - S * (epsilon))/exp(Wu/2) -
    (1 + lambda) * exp(Wv)/exp(Wu)) * db * wvwu + exp(Wv) *
    pb/exp(Wu)) * eB/(papb) + 1/(1 + lambda)^2))
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll<-(hessll+(hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

tslnormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
  uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method, printInfo,
  itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start))
    start else csttslnorm(olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
    S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
    nvZVvar = nvZVvar)
  startLoglik <- sum(ctslnormlike(startVal, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...))
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(par = startVal,
    fn = function(parm) -sum(ctslnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    gr = function(parm) -colSums(cgradtslnormlike(parm, nXvar = nXvar,
      nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S)),
    hessian = 0, control = list(trace = if (printInfo) 1 else 0, maxeval = itermax,
      stepmax = stepmax, xtol = tol, grtol = gradtol)),
    maxLikAlgo = maxRoutine(fn = ctslnormlike, grad = cgradtslnormlike,
      hess = chesstslnormlike, start = startVal, finalHessian = if (hessianType ==
        2) "bhhh" else TRUE, control = list(printLevel = if (printInfo) 2 else 0,
        iterlim = itermax, reltol = tol, tol = tol, qac = qac),
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), sr1 = trust.optim(x = startVal, fn = function(parm) -sum(ctslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), method = "SR1", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L)),
    sparse = trust.optim(x = startVal, fn = function(parm) -sum(ctslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hs = function(parm) as(-chesstslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), "dgCMatrix"), method = "Sparse", control = list(maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L,
      preconditioner = 1L)), mla = mla(b = startVal, fn = function(parm) -sum(ctslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), gr = function(parm) -colSums(cgradtslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S)), hess = function(parm) -chesstslnormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S), print.info = printInfo, maxiter = itermax,
      epsa = gradtol, epsb = gradtol), nlminb = nlminb(start = startVal,
      objective = function(parm) -sum(ctslnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), gradient = function(parm) -colSums(cgradtslnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)), hessian = function(parm) -chesstslnormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S), control = list(iter.max = itermax, trace = if (printInfo) 1 else 0,
        eval.max = itermax, rel.tol = tol, x.tol = tol)))
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradtslnormlike(mleObj$par,
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
      mleObj$hessian <- chesstslnormlike(parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
    if (method == "sr1")
      mleObj$hessian <- chesstslnormlike(parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S)
  }
  mleObj$logL_OBS <- ctslnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  mleObj$gradL_OBS <- cgradtslnormlike(parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S)
  return(list(startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam))
}

# Conditional efficiencies estimation ----------

ctslnormeff <- function(object, level) {
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
  lambda <- Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  A <- exp(Wv)/(2 * exp(Wu) + object$S * epsilon/exp(Wu/2))
  B <- (1 + lambda)^2 * exp(Wv)/(2 * exp(Wu) + object$S * epsilon *
    (1 + lambda)/exp(Wu/2))
  a <- -exp(Wv/2)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  b <- -exp(Wv/2) * (1 + lambda)/exp(Wu/2) - object$S * epsilon/exp(Wv/2)
  u <- exp(Wv/2) * (2 * exp(A) * (dnorm(a) + a * pnorm(a)) -
    exp(B) * (dnorm(b) + b * dnorm(b)))/(2 * exp(A) * pnorm(a) -
    exp(B) * pnorm(b))
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teBC <- (2 * exp(A) * exp(-a * exp(Wv/2) + exp(Wv)/2) *
      pnorm(a - exp(Wv/2)) - exp(B) * exp(-b * exp(Wv/2) +
      exp(Wv)/2) * pnorm(b - exp(Wv/2)))/(2 * exp(A) *
      pnorm(a) - exp(B) * pnorm(b))
    res <- bind_cols(u = u, teJLMS = teJLMS, teBC = teBC)
  } else {
    res <- bind_cols(u = u)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargtslnorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  lambda <- Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    4 * lambda + 2 * lambda^2)/((1 + lambda) * (1 + 2 * lambda)),
    nrow = 1), matrix(exp(Wu/2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmargtslnorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula, data = object$dataTable,
    rhs = 2)
  lambda <- Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * (1 +
    8 * lambda + 16 * lambda^2 + 12 * lambda^3 + 4 * lambda^4)/((1 +
    lambda)^2 * (1 + 2 * lambda)^2), nrow = 1), matrix(exp(Wu),
    ncol = 1))
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
