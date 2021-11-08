########################################################
#                                                      #
# Exponential + normal distributions                   #
#                                                      #
#                                                      #
########################################################

# Log-likelihood ----------

cexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
                          vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  ll <- (-Wu / 2 + log(pnorm(-S * epsilon / sqrt(exp(Wv)) - sqrt(exp(Wv) / exp(Wu)))) +
    S * epsilon / sqrt(exp(Wu)) + exp(Wv) / (2 * exp(Wu)))
  return(ll)
}

# starting value for the log-likelihood ----------

cstexponorm <- function(olsObj, epsiRes, S, nuZUvar, uHvar, nvZVvar,
                        vHvar) {
  m2 <- moment(epsiRes, order = 2)
  m3 <- moment(epsiRes, order = 3)
  if (S * m3 > 0) {
    varu <- (abs((-S * m3 / 2)))^(2 / 3)
  } else {
    varu <- (-S * m3 / 2)^(2 / 3)
  }
  if (m2 < varu) {
    varv <- abs(m2 - varu)
  } else {
    varv <- m2 - varu
  }
  dep_u <- 1 / 2 * log((epsiRes^2 - varv)^2)
  dep_v <- 1 / 2 * log((epsiRes^2 - varu)^2)
  reg_hetu <- if (nuZUvar == 1) {
    lm(log(varu) ~ 1)
  } else {
    lm(dep_u ~ ., data = as.data.frame(uHvar[, 2:nuZUvar]))
  }
  if (any(is.na(reg_hetu$coefficients))) {
    stop("At least one of the OLS coefficients of 'uhet' is NA: ",
      paste(colnames(uHvar)[is.na(reg_hetu$coefficients)],
        collapse = ", "
      ), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE
    )
  }
  reg_hetv <- if (nvZVvar == 1) {
    lm(log(varv) ~ 1)
  } else {
    lm(dep_v ~ ., data = as.data.frame(vHvar[, 2:nvZVvar]))
  }
  if (any(is.na(reg_hetv$coefficients))) {
    stop("at least one of the OLS coefficients of 'vhet' is NA: ",
      paste(colnames(vHvar)[is.na(reg_hetv$coefficients)],
        collapse = ", "
      ), ". This may be due to a singular matrix due to potential perfect multicollinearity",
      call. = FALSE
    )
  }
  delta <- coefficients(reg_hetu)
  names(delta) <- paste0("Zu_", colnames(uHvar))
  phi <- coefficients(reg_hetv)
  names(phi) <- paste0("Zv_", colnames(vHvar))
  if (names(olsObj)[1] == "(Intercept)") {
    beta <- c(olsObj[1] + S * sqrt(varu), olsObj[-1])
  } else {
    beta <- olsObj
  }
  return(c(beta, delta, phi))
}

# Gradient of the likelihood function ----------

cgradexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
                              vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -(S * (epsilon) / exp(Wv / 2) + sqrt(exp(Wv) / exp(Wu)))
  pmustar <- pnorm(mustar)
  dmustar <- dnorm(mustar)
  sigx <- exp(Wu) * pmustar * sqrt(exp(Wv) / exp(Wu))
  sigx2 <- dmustar / (exp(Wv / 2) * pmustar)
  pdmustar <- dmustar / pmustar
  su_sv <- sqrt(exp(Wv) / exp(Wu))
  sv_epsi <- S * (epsilon) / exp(Wv / 2)
  su_epsi <- S * (epsilon) / exp(Wu / 2)
  sigx3 <- 0.5 * (exp(Wv) / (exp(Wu) * su_sv)) - 0.5 * (sv_epsi)
  gradll <- cbind(sweep(Xvar, MARGIN = 1, STATS = S * (sigx2 -
    1 / exp(Wu / 2)), FUN = "*"), sweep(uHvar, MARGIN = 1, STATS = ((0.5 *
    (dmustar / (sigx)) - 1 / (2 * exp(Wu))) * exp(Wv) - (0.5 +
    0.5 * (su_epsi))), FUN = "*"), sweep(vHvar,
    MARGIN = 1,
    STATS = (exp(Wv) / (2 * exp(Wu)) - (sigx3) * pdmustar),
    FUN = "*"
  ))
  return(gradll)
}

# Hessian of the likelihood function ----------

chessexponormlike <- function(parm, nXvar, nuZUvar, nvZVvar, uHvar,
                              vHvar, Yvar, Xvar, S) {
  beta <- parm[1:(nXvar)]
  delta <- parm[(nXvar + 1):(nXvar + nuZUvar)]
  phi <- parm[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)]
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- Yvar - as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -(S * (epsilon) / exp(Wv / 2) + sqrt(exp(Wv) / exp(Wu)))
  pmustar <- pnorm(mustar)
  dmustar <- dnorm(mustar)
  sigx <- exp(Wu) * pmustar * sqrt(exp(Wv) / exp(Wu))
  sigx2 <- dmustar / (exp(Wv / 2) * pmustar)
  pdmustar <- dmustar / pmustar
  su_sv <- sqrt(exp(Wv) / exp(Wu))
  sv_epsi <- S * (epsilon) / exp(Wv / 2)
  su_epsi <- S * (epsilon) / exp(Wu / 2)
  sigx3 <- 0.5 * (exp(Wv) / (exp(Wu) * su_sv)) - 0.5 * (sv_epsi)
  hessll <- matrix(nrow = nXvar + nuZUvar + nvZVvar, ncol = nXvar +
    nuZUvar + nvZVvar)
  hessll[1:nXvar, 1:nXvar] <- crossprod(
    sweep(Xvar,
      MARGIN = 1,
      STATS = S^2 * ((sv_epsi + su_sv) / (exp(Wv / 2)^2 * pmustar) -
        sigx2 / (exp(Wv / 2) * pmustar)) * dmustar, FUN = "*"
    ),
    Xvar
  )
  hessll[1:nXvar, (nXvar + 1):(nXvar + nuZUvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = S * (0.5 * (((sv_epsi + su_sv) / (sigx) -
      dmustar * exp(Wu) * su_sv / (sigx)^2) * dmustar * exp(Wv / 2)) +
      0.5 / exp(Wu / 2)), FUN = "*"
  ), uHvar)
  hessll[1:nXvar, (nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(Xvar,
    MARGIN = 1, STATS = -(S * ((sigx3) * (sv_epsi + su_sv -
      pdmustar) + 0.5) * sigx2), FUN = "*"
  ), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + 1):(nXvar +
    nuZUvar)] <- crossprod(sweep(uHvar, MARGIN = 1, STATS = ((0.5 *
    ((0.5 * ((sv_epsi + su_sv) / (exp(Wu) * pmustar)) - ((0.5 *
      (dmustar * exp(Wv) / su_sv) + exp(Wu) * pmustar) *
      su_sv - 0.5 * (exp(Wv) * pmustar / su_sv)) / (sigx)^2) *
      dmustar) + (1 / (2 * exp(Wu)))) * exp(Wv) + 0.25 *
    (su_epsi)), FUN = "*"), uHvar)
  hessll[(nXvar + nuZUvar + 1):(nXvar + nuZUvar + nvZVvar), (nXvar +
    nuZUvar + 1):(nXvar + nuZUvar + nvZVvar)] <- crossprod(sweep(vHvar,
    MARGIN = 1, STATS = (exp(Wv) / (2 * exp(Wu)) - ((sigx3)^2 *
      (pdmustar + mustar) + 0.25 * (sv_epsi) + 0.5 * ((1 / exp(Wu) -
        0.5 * (exp(Wv) / (exp(Wu) * su_sv)^2)) * exp(Wv) / su_sv)) *
      pdmustar), FUN = "*"
  ), vHvar)
  hessll[(nXvar + 1):(nXvar + nuZUvar), (nXvar + nuZUvar + 1):(nXvar +
    nuZUvar + nvZVvar)] <- crossprod(
    sweep(uHvar,
      MARGIN = 1,
      STATS = -((((sigx3) * (0.5 * (sv_epsi + su_sv) - 0.5 *
        (pdmustar)) / (exp(Wu) * su_sv) - 0.5 * ((exp(Wu) *
        su_sv - 0.5 * (exp(Wv) / su_sv)) / (exp(Wu) * su_sv)^2)) *
        pdmustar + 1 / (2 * exp(Wu))) * exp(Wv)), FUN = "*"
    ),
    vHvar
  )
  hessll[lower.tri(hessll)] <- t(hessll)[lower.tri(hessll)]
  # hessll <- (hessll + (hessll))/2
  return(hessll)
}

# Optimization using different algorithms ----------

exponormAlgOpt <- function(start, olsParam, dataTable, S, nXvar,
                           uHvar, nuZUvar, vHvar, nvZVvar, Yvar, Xvar, method, printInfo,
                           itermax, stepmax, tol, gradtol, hessianType, qac) {
  startVal <- if (!is.null(start)) {
    start
  } else {
    cstexponorm(
      olsObj = olsParam, epsiRes = dataTable[["olsResiduals"]],
      S = S, uHvar = uHvar, nuZUvar = nuZUvar, vHvar = vHvar,
      nvZVvar = nvZVvar
    )
  }
  startLoglik <- sum(cexponormlike(startVal,
    nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S
  ))
  if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
    maxRoutine <- switch(method, bfgs = function(...) maxBFGS(...),
      bhhh = function(...) maxBHHH(...), nr = function(...) maxNR(...),
      nm = function(...) maxNM(...)
    )
    method <- "maxLikAlgo"
  }
  mleObj <- switch(method, ucminf = ucminf(
    par = startVal,
    fn = function(parm) {
      -sum(cexponormlike(parm,
        nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S
      ))
    },
    gr = function(parm) {
      -colSums(cgradexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, hessian = 0, control = list(
      trace = if (printInfo) 1 else 0,
      maxeval = itermax, stepmax = stepmax, xtol = tol,
      grtol = gradtol
    )
  ), maxLikAlgo = maxRoutine(
    fn = cexponormlike,
    grad = cgradexponormlike, hess = chessexponormlike, start = startVal,
    finalHessian = if (hessianType == 2) "bhhh" else TRUE,
    control = list(
      printLevel = if (printInfo) 2 else 0,
      iterlim = itermax, reltol = tol, tol = tol, qac = qac
    ),
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
    vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S
  ), sr1 = trust.optim(
    x = startVal,
    fn = function(parm) {
      -sum(cexponormlike(parm,
        nXvar = nXvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, S = S
      ))
    },
    gr = function(parm) {
      -colSums(cgradexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, method = "SR1", control = list(
      maxit = itermax,
      cgtol = gradtol, stop.trust.radius = tol, prec = tol,
      report.level = if (printInfo) 2 else 0, report.precision = 1L
    )
  ),
  sparse = trust.optim(x = startVal, fn = function(parm) {
    -sum(cexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S
    ))
  }, gr = function(parm) {
    -colSums(cgradexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S
    ))
  }, hs = function(parm) {
    as(-chessexponormlike(parm,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S
    ), "dgCMatrix")
  }, method = "Sparse", control = list(
    maxit = itermax,
    cgtol = gradtol, stop.trust.radius = tol, prec = tol,
    report.level = if (printInfo) 2 else 0, report.precision = 1L,
    preconditioner = 1L
  )), mla = mla(
    b = startVal, fn = function(parm) {
      -sum(cexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, gr = function(parm) {
      -colSums(cgradexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, hess = function(parm) {
      -chessexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      )
    }, print.info = printInfo, maxiter = itermax,
    epsa = gradtol, epsb = gradtol
  ), nlminb = nlminb(
    start = startVal,
    objective = function(parm) {
      -sum(cexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, gradient = function(parm) {
      -colSums(cgradexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      ))
    }, hessian = function(parm) {
      -chessexponormlike(parm,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      )
    }, control = list(
      iter.max = itermax, trace = if (printInfo) 1 else 0,
      eval.max = itermax, rel.tol = tol, x.tol = tol
    )
  )
  )
  if (method %in% c("ucminf", "nlminb")) {
    mleObj$gradient <- colSums(cgradexponormlike(mleObj$par,
      nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
      uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
      S = S
    ))
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
    if (method %in% c("ucminf", "nlminb")) {
      mleObj$hessian <- chessexponormlike(
        parm = mleObj$par,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      )
    }
    if (method == "sr1") {
      mleObj$hessian <- chessexponormlike(
        parm = mleObj$solution,
        nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
        uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
        S = S
      )
    }
  }
  mleObj$logL_OBS <- cexponormlike(
    parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S
  )
  mleObj$gradL_OBS <- cgradexponormlike(
    parm = mlParam, nXvar = nXvar,
    nuZUvar = nuZUvar, nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar,
    Yvar = Yvar, Xvar = Xvar, S = S
  )
  return(list(
    startVal = startVal, startLoglik = startLoglik,
    mleObj = mleObj, mlParam = mlParam
  ))
}

# Conditional efficiencies estimation ----------

cexponormeff <- function(object, level) {
  beta <- object$mlParam[1:(object$nXvar)]
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  phi <- object$mlParam[(object$nXvar + object$nuZUvar + 1):(object$nXvar +
    object$nuZUvar + object$nvZVvar)]
  Xvar <- model.matrix(object$formula,
    data = object$dataTable,
    rhs = 1
  )
  uHvar <- model.matrix(object$formula,
    data = object$dataTable,
    rhs = 2
  )
  vHvar <- model.matrix(object$formula,
    data = object$dataTable,
    rhs = 3
  )
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  epsilon <- model.response(model.frame(object$formula, data = object$dataTable)) -
    as.numeric(crossprod(matrix(beta), t(Xvar)))
  mustar <- -object$S * epsilon - exp(Wv) / sqrt(exp(Wu))
  u <- mustar + sqrt(exp(Wv)) * dnorm(mustar / sqrt(exp(Wv))) / pnorm(mustar / sqrt(exp(Wv)))
  uLB <- mustar + qnorm(1 - (1 - (1 - level) / 2) * (1 - pnorm(-mustar / sqrt(exp(Wv))))) *
    sqrt(exp(Wv))
  uUB <- mustar + qnorm(1 - (1 - level) / 2 * (1 - pnorm(-mustar / sqrt(exp(Wv))))) *
    sqrt(exp(Wv))
  m <- ifelse(mustar > 0, mustar, 0)
  if (object$logDepVar == TRUE) {
    teJLMS <- exp(-u)
    teMO <- exp(-m)
    teBC <- exp(-mustar + 1 / 2 * exp(Wv)) * pnorm(mustar / sqrt(exp(Wv)) -
      sqrt(exp(Wv))) / pnorm(mustar / sqrt(exp(Wv)))
    teBCLB <- exp(-uUB)
    teBCUB <- exp(-uLB)
    res <- bind_cols(
      u = u, uLB = uLB, uUB = uUB, teJLMS = teJLMS,
      m = m, teMO = teMO, teBC = teBC, teBCLB = teBCLB,
      teBCUB = teBCUB
    )
  } else {
    res <- bind_cols(u = u, uLB = uLB, uUB = uUB, m = m)
  }
  return(res)
}

# Marginal effects on inefficiencies ----------

cmargexponorm_Eu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula,
    data = object$dataTable,
    rhs = 2
  )
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(matrix(delta[2:object$nuZUvar] * 1 / 2,
    nrow = 1
  ), matrix(exp(Wu / 2), ncol = 1))
  colnames(margEff) <- paste0("Eu_", colnames(uHvar)[-1])
  return(margEff)
}

cmargexponorm_Vu <- function(object) {
  delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
    object$nuZUvar)]
  uHvar <- model.matrix(object$formula,
    data = object$dataTable,
    rhs = 2
  )
  Wu <- as.numeric(crossprod(matrix(delta), t(uHvar)))
  margEff <- kronecker(
    matrix(delta[2:object$nuZUvar], nrow = 1),
    matrix(exp(Wu), ncol = 1)
  )
  colnames(margEff) <- paste0("Vu_", colnames(uHvar)[-1])
  return(margEff)
}
