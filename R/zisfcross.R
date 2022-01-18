# ZISF estimation for cross sectional (pooled) data
# ----------

zisfcross <- function(formula, muhet, uhet, vhet, thet, logDepVar = TRUE,
  data, subset, weights, wscale = TRUE, S = 1L, udist = "hnormal",
  start = NULL, method = "bfgs", hessianType = 1, simType = "halton",
  Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE, seed = 12345,
  itermax = 2000L, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06,
  stepmax = 0.1, qac = "marquardt") {
  # u distribution check -------
  udist <- tolower(udist)
  if (!(udist %in% c("hnormal", "exponential", "tnormal", "rayleigh",
    "uniform", "gamma", "lognormal", "weibull", "genexponential",
    "tslaplace"))) {
    stop("Unknown inefficiency distribution: ", paste(udist),
      call. = FALSE)
  }
  # Formula manipulation -------
  if (length(Formula(formula))[2] != 1) {
    stop("argument 'formula' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mc), nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$drop.unused.levels <- TRUE
  formula <- interCheckMain(formula = formula)
  if (!missing(muhet)) {
    muhet <- clhsCheck_mu(formula = muhet, scaling = FALSE)
  } else {
    muhet <- ~1
  }
  if (!missing(uhet)) {
    uhet <- clhsCheck_u(formula = uhet, scaling = FALSE)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- clhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  if (!missing(thet)) {
    thet <- clhsCheck_q(formula = thet)
  } else {
    thet <- ~1
  }
  formula <- formDist_zisfcross(udist = udist, formula = formula,
    muhet = muhet, uhet = uhet, vhet = vhet, qhet = thet)
  # Generate required datasets -------
  if (missing(data)) {
    data <- environment(formula)
  }
  mc$formula <- formula
  mc$na.action <- na.pass
  mc[[1L]] <- quote(model.frame)
  mc <- eval(mc, parent.frame())
  validObs <- rowSums(is.na(mc) | is.infinite.data.frame(mc)) ==
    0
  Yvar <- model.response(mc, "numeric")
  Yvar <- Yvar[validObs]
  mtX <- terms(formula, data = data, rhs = 1)
  Xvar <- model.matrix(mtX, mc)
  Xvar <- Xvar[validObs, , drop = FALSE]
  nXvar <- ncol(Xvar)
  N <- nrow(Xvar)
  if (N == 0L) {
    stop("0 (non-NA) cases", call. = FALSE)
  }
  wHvar <- as.vector(model.weights(mc))
  if (length(wscale) != 1 || !is.logical(wscale[1])) {
    stop("argument 'wscale' must be a single logical value",
      call. = FALSE)
  }
  if (!is.null(wHvar)) {
    if (!is.numeric(wHvar)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    } else {
      if (any(wHvar < 0 | is.na(wHvar)))
        stop("missing or negative weights not allowed",
          call. = FALSE)
    }
    if (wscale) {
      wHvar <- wHvar/sum(wHvar) * N
    }
  } else {
    wHvar <- rep(1, N)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (",
      length(Yvar), ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula, data = data,
      rhs = 2))
    muHvar <- model.matrix(mtmuH, mc)
    muHvar <- muHvar[validObs, , drop = FALSE]
    nmuZUvar <- ncol(muHvar)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
    mtZ <- delete.response(terms(formula, data = data, rhs = 5))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
    mtZ <- delete.response(terms(formula, data = data, rhs = 4))
    Zvar <- model.matrix(mtZ, mc)
    Zvar <- Zvar[validObs, , drop = FALSE]
    nZHvar <- ncol(Zvar)
  }
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Zero Inefficiency Production/Profit Frontier, e = v - u"
  } else {
    "Zero Inefficiency Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  # Number of parameters -------
  nParm <- if (udist %in% c("tnormal", "lognormal")) {
    nXvar + nmuZUvar + nuZUvar + nvZVvar + nZHvar
  } else {
    if (udist %in% c("gamma", "weibull", "tslaplace")) {
      nXvar + nuZUvar + nvZVvar + nZHvar + 1
    } else {
      nXvar + nuZUvar + nvZVvar + nZHvar
    }
  }
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ",
        nParm, " parameters", call. = FALSE)
    }
  }
  if (nParm > N) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm",
    "cg", "sann", "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method), call. = FALSE)
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L,
    2L, 3L))) {
    stop("argument 'hessianType' must equal either 1 or 2 or 3",
      call. = FALSE)
  }
  # Draws for SML -------
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    if (!(simType %in% c("halton", "ghalton", "sobol", "uniform"))) {
      stop("Unknown or non-available random draws method",
        call. = FALSE)
    }
    if (!is.numeric(Nsim) || length(Nsim) != 1) {
      stop("argument 'Nsim' must be a single numeric scalar",
        call. = FALSE)
    }
    if (!is.numeric(burn) || length(burn) != 1) {
      stop("argument 'burn' must be a single numeric scalar",
        call. = FALSE)
    }
    if (!is_prime(prime)) {
      stop("argument 'prime' must be a single prime number",
        call. = FALSE)
    }
    if (length(antithetics) != 1 || !is.logical(antithetics[1])) {
      stop("argument 'antithetics' must be a single logical value",
        call. = FALSE)
    }
    if (antithetics && (Nsim%%2) != 0) {
      Nsim <- Nsim + 1
    }
    simDist <- if (simType == "halton") {
      "Halton"
    } else {
      if (simType == "ghalton") {
        "Generalized Halton"
      } else {
        if (simType == "sobol") {
          "Sobol"
        } else {
          if (simType == "uniform") {
          "Uniform"
          }
        }
      }
    }
    cat("Initialization of", Nsim, simDist, "draws per observation ...\n")
    FiMat <- drawMat(N = N, Nsim = Nsim, simType = simType,
      prime = prime, burn = burn + 1, antithetics = antithetics,
      seed = seed)
  }
  # Other optimization options -------
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar",
      call. = FALSE)
  }
  if (itermax != round(itermax)) {
    stop("argument 'itermax' must be an integer", call. = FALSE)
  }
  if (itermax <= 0) {
    stop("argument 'itermax' must be positive", call. = FALSE)
  }
  itermax <- as.integer(itermax)
  if (length(printInfo) != 1 || !is.logical(printInfo[1])) {
    stop("argument 'printInfo' must be a single logical value",
      call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1) {
    stop("argument 'tol' must be numeric", call. = FALSE)
  }
  if (tol < 0) {
    stop("argument 'tol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(gradtol) || length(gradtol) != 1) {
    stop("argument 'gradtol' must be numeric", call. = FALSE)
  }
  if (gradtol < 0) {
    stop("argument 'gradtol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(stepmax) || length(stepmax) != 1) {
    stop("argument 'stepmax' must be numeric", call. = FALSE)
  }
  if (stepmax < 0) {
    stop("argument 'stepmax' must be non-negative", call. = FALSE)
  }
  if (!(qac %in% c("marquardt", "stephalving"))) {
    stop("argument 'qac' must be either 'marquardt' or 'stephalving'",
      call. = FALSE)
  }
  # Step 1: OLS -------
  olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
    if (dim(Xvar)[2] == 1) {
      lm(Yvar ~ 1)
    } else {
      lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]), weights = wHvar)
    }
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar), weights = wHvar)
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ",
      paste(colnames(Xvar)[is.na(olsRes$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  olsParam <- c(olsRes$coefficients)
  if (inherits(data, "plm.dim")) {
    dataTable <- data[validObs, 1:2]
  } else {
    dataTable <- data.frame(IdObs = c(1:sum(validObs)))
  }
  dataTable <- as_tibble(cbind(dataTable, data[, all.vars(terms(formula))],
    weights = wHvar))
  dataTable <- mutate(dataTable, olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes))
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are ", if (S == 1) {
      " right"
    } else {
      "left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
      call. = FALSE)
  }
  FunArgs <- if (udist == "tnormal") {
    list(start = start, olsParam = olsParam, dataTable = dataTable,
      nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
      nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
      vHvar = vHvar, Yvar = Yvar, Zvar = Zvar, nZHvar = nZHvar,
      Xvar = Xvar, S = S, wHvar = wHvar, method = method,
      printInfo = printInfo, itermax = itermax, stepmax = stepmax,
      tol = tol, gradtol = gradtol, hessianType = hessianType,
      qac = qac)
  } else {
    if (udist == "lognormal") {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvar, Zvar = Zvar, nZHvar = nZHvar,
        Xvar = Xvar, S = S, wHvar = wHvar, N = N, FiMat = FiMat,
        method = method, printInfo = printInfo, itermax = itermax,
        stepmax = stepmax, tol = tol, gradtol = gradtol,
        hessianType = hessianType, qac = qac)
    } else {
      if (udist %in% c("gamma", "weibull")) {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Zvar = Zvar, nZHvar = nZHvar, Xvar = Xvar,
          S = S, wHvar = wHvar, N = N, FiMat = FiMat,
          method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
      } else {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar,
          Zvar = Zvar, nZHvar = nZHvar, Xvar = Xvar,
          S = S, wHvar = wHvar, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol,
          gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      }
    }
  }
  ## MLE run -------
  mleList <- tryCatch(switch(udist, hnormal = do.call(zisfhalfnormAlgOpt,
    FunArgs), exponential = do.call(zisfexponormAlgOpt, FunArgs),
    tnormal = do.call(zisftruncnormAlgOpt, FunArgs), rayleigh = do.call(zisfraynormAlgOpt,
      FunArgs), gamma = do.call(zisfgammanormAlgOpt, FunArgs),
    uniform = do.call(zisfuninormAlgOpt, FunArgs), lognormal = do.call(zisflognormAlgOpt,
      FunArgs), weibull = do.call(zisfweibullnormAlgOpt,
      FunArgs), genexponential = do.call(zisfgenexponormAlgOpt,
      FunArgs), tslaplace = do.call(zisftslnormAlgOpt,
      FunArgs)), error = function(e) e)
  if (inherits(mleList, "error")) {
    stop("The current error occurs during optimization:\n",
      mleList$message, call. = FALSE)
  }
  # Inverse Hessian + other -------
  mleList$invHessian <- vcovObj(mleObj = mleList$mleObj, hessianType = hessianType,
    method = method, nParm = nParm)
  mleList <- c(mleList, if (method == "ucminf") {
    list(type = "ucminf maximization", nIter = unname(mleList$mleObj$info["neval"]),
      status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
      gradient = mleList$mleObj$gradient)
  } else {
    if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
      list(type = mleList$mleObj$type, nIter = mleList$mleObj$iterations,
        status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method == "sr1") {
        list(type = "SR1 maximization", nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
          gradient = mleList$mleObj$gradient)
      } else {
        if (method == "mla") {
          list(type = "Levenberg-Marquardt maximization",
          nIter = mleList$mleObj$ni, status = switch(mleList$mleObj$istop,
            `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached",
            `4` = "algorithm encountered a problem in the function computation"),
          mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad)
        } else {
          if (method == "sparse") {
          list(type = "Sparse Hessian maximization",
            nIter = mleList$mleObj$iterations, status = mleList$mleObj$status,
            mleLoglik = -mleList$mleObj$fval, gradient = mleList$mleObj$gradient)
          } else {
          if (method == "nlminb") {
            list(type = "nlminb maximization", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
            gradient = mleList$mleObj$gradient)
          }
          }
        }
      }
    }
  })
  # quick renaming -------
  if (udist %in% c("tnormal", "lognormal")) {
    names(mleList$startVal) <- fName_mu_zisfcross(Xvar = Xvar,
      muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, Zvar = Zvar)
  } else {
    names(mleList$startVal) <- fName_uv_zisfcross(Xvar = Xvar,
      udist = udist, uHvar = uHvar, vHvar = vHvar, Zvar = Zvar)
  }
  names(mleList$mlParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- N
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuZUvar <- nmuZUvar
  }
  returnObj$nZHvar <- nZHvar
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  if (is.null(start)) {
    if (udist == "hnormal") {
      returnObj$InitHalf <- mleList$InitHalf
    } else {
      if (udist == "exponential") {
        returnObj$initExpo <- mleList$initExpo
      } else {
        if (udist == "tnormal") {
          returnObj$initTrunc <- mleList$initTrunc
        } else {
          if (udist == "rayleigh") {
          returnObj$initRay <- mleList$initRay
          } else {
          if (udist == "uniform") {
            returnObj$initUni <- mleList$initUni
          } else {
            if (udist == "gamma") {
            returnObj$initGamma <- mleList$initGamma
            } else {
            if (udist == "lognormal") {
              returnObj$initLog <- mleList$initLog
            } else {
              if (udist == "weibull") {
              returnObj$initWeibull <- mleList$initWeibull
              } else {
              if (udist == "genexponential") {
                returnObj$initGenExpo <- mleList$initGenExpo
              } else {
                if (udist == "tslaplace") {
                returnObj$initTSL <- mleList$initTSL
                }
              }
              }
            }
            }
          }
          }
        }
      }
    }
  }
  returnObj$isWeights <- !all.equal(wHvar, rep(1, N))
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$gradient <- mleList$gradient
  returnObj$gradL_OBS <- mleList$mleObj$gradL_OBS
  returnObj$gradientNorm <- sqrt(sum(mleList$gradient^2))
  returnObj$invHessian <- mleList$invHessian
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic/Numeric Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    } else {
      if (hessianType == 3) {
        "Robust Hessian"
      }
    }
  }
  returnObj$mlDate <- mleDate
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat
  }
  rm(mleList)
  class(returnObj) <- "zisfcross"
  return(returnObj)
}

# print for zisfcross ----------

print.zisfcross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call, width.cutoff = 500))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(zisfdist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}
