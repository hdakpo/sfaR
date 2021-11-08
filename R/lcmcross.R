# LCM estimation for cross sectional data ----------

lcmcross <- function(formula, uhet, vhet, thet, logDepVar = TRUE, data, subset,
                     S = 1L, udist = "hnormal", start = NULL, lcmClasses = 2,
                     method = "bfgs", hessianType = 1, itermax = 2000L,
                     printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
                     qac = "marquardt", initStart = FALSE, initAlg = "nlminb",
                     initIter = 100, initFactorLB = 0.5, initFactorUB = 1.5) {
  # u distribution check -------
  udist <- tolower(udist)
  if (udist != "hnormal") {
    stop("Currently latent class model only handles half-normal distribution ... ", 
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
  if (!missing(uhet)) {
    uhet <- lhsCheck_u(formula = uhet, scaling = FALSE)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- lhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  if (!missing(thet)) {
    thet <- lhsCheck_t(formula = thet)
  } else {
    thet <- ~1
  }
  formula <- formDist_lcmcross(formula = formula, uhet = uhet, 
    vhet = vhet, thet = thet)
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
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (", 
      length(Yvar), ") must be the same to the number of observations of the exogenous variables (", 
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
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
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier", 
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Latent Class Production/Profit Frontier, e = v - u"
  } else {
    "Latent Class Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value", 
      call. = FALSE)
  }
  if (length(initStart) != 1 || !is.logical(initStart[1])) {
    stop("argument 'initStart' must be a single logical value", 
      call. = FALSE)
  }
  if (!(lcmClasses %in% 2:5)) {
    stop("argument 'lcmClasses' must be comprised between 2 and 5", 
      call. = FALSE)
  }
  # Number of parameters -------
  nParm <- lcmClasses * (nXvar + nuZUvar + nvZVvar) + (lcmClasses - 
    1) * nZHvar
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
    "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ", 
      paste(method), call. = FALSE)
  }
  initAlg <- tolower(initAlg)
  if (initAlg != "nlminb") {
    stop("Only 'nlminb' is available as initialization algorithm ")
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L, 
    2L, 3L))) {
    stop("argument 'hessianType' must equal either 1 or 2 or 3", 
      call. = FALSE)
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
  if (!is.numeric(initIter) || length(initIter) != 1) {
    stop("argument 'initIter' must be a single numeric scalar", 
      call. = FALSE)
  }
  if (initIter != round(initIter)) {
    stop("argument 'initIter' must be an integer", call. = FALSE)
  }
  if (initIter <= 0) {
    stop("argument 'initIter' must be positive", call. = FALSE)
  }
  initIter <- as.integer(initIter)
  if (!is.numeric(initFactorLB) || length(initFactorLB) != 
    1) {
    stop("argument 'initFactorLB' must be a single numeric value", 
      call. = FALSE)
  }
  if (!is.numeric(initFactorUB) || length(initFactorUB) != 
    1) {
    stop("argument 'initFactorUB' must be a single numeric value", 
      call. = FALSE)
  }
  # Step 1: OLS -------
  olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
    lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]))
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar))
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ", 
      paste(colnames(Xvar)[is.na(olsRes$coefficients)], 
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity", 
      call. = FALSE)
  }
  olsParam <- c(olsRes$coefficients)
  # olsSigmasq <- summary(olsRes)$sigma^2
  # olsStder <- sqrt(diag(vcov(olsRes)))
  # olsLoglik <- logLik(olsRes)[1]
  if (inherits(data, "plm.dim")) {
    dataTable <- data[validObs, 1:2]
  } else {
    dataTable <- data.frame(IdObs = c(1:sum(validObs)))
  }
  dataTable <- as_tibble(cbind(dataTable, data[validObs, all.vars(terms(formula))]))
  dataTable <- mutate(dataTable, olsResiduals = residuals(olsRes), 
    olsFitted = fitted(olsRes))
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  # olsM3Okay <- if (S * olsSkew < 0) { 'Residuals have the
  # 'right' skeweness' } else { 'Residuals have the 'wrong'
  # skeweness' }
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are ", if (S == 1) {
      " right"
    } else {
      "left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'", 
      call. = FALSE)
  }
  # CoelliM3Test <- c(z = moment(dataTable[['olsResiduals']],
  # order = 3)/sqrt(6 * moment(dataTable[['olsResiduals']],
  # order = 2)^3/N), p.value = 2 *
  # pnorm(-abs(moment(dataTable[['olsResiduals']], order =
  # 3)/sqrt(6 * moment(dataTable[['olsResiduals']], order =
  # 2)^3/N)))) AgostinoTest <-
  # dagoTest(dataTable[['olsResiduals']])@test Step 2: MLE
  # arguments -------
  FunArgs <- list(start = start, olsParam = olsParam, dataTable = dataTable, 
    nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar, 
    uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Zvar = Zvar, 
    nZHvar = nZHvar, Xvar = Xvar, S = S, method = method, 
    printInfo = printInfo, itermax = itermax, stepmax = stepmax, 
    tol = tol, gradtol = gradtol, hessianType = hessianType, 
    qac = qac, initStart = initStart, initAlg = initAlg, 
    initIter = initIter, initFactorLB = initFactorLB, initFactorUB = initFactorUB)
  ## MLE run -------
  mleList <- tryCatch(switch(as.character(lcmClasses), `2` = do.call(LCM2ChnormAlgOpt, 
    FunArgs), `3` = do.call(LCM3ChnormAlgOpt, FunArgs), `4` = do.call(LCM4ChnormAlgOpt, 
    FunArgs), `5` = do.call(LCM5ChnormAlgOpt, FunArgs)), 
    error = function(e) e)
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
    if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
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
  names(mleList$startVal) <- fName_lcmcross(Xvar = Xvar, uHvar = uHvar, 
    vHvar = vHvar, Zvar = Zvar, nZHvar = nZHvar, lcmClasses = lcmClasses)
  names(mleList$mlParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mleDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals_c1 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]), 
    t(Xvar)))
  dataTable$mlFitted_c1 <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]), 
    t(Xvar)))
  dataTable$mlResiduals_c2 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(nXvar + 
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]), 
    t(Xvar)))
  dataTable$mlFitted_c2 <- as.numeric(crossprod(matrix(mleList$mlParam[(nXvar + 
    nuZUvar + nvZVvar + 1):(2 * nXvar + nuZUvar + nvZVvar)]), 
    t(Xvar)))
  if (lcmClasses == 3) {
    dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
    dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
      nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
      2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
  } else {
    if (lcmClasses == 4) {
      dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
        2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
        nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * nXvar + 
        2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
      dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 * 
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
        3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
      dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 * 
        nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * nXvar + 
        3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
    } else {
      if (lcmClasses == 5) {
        dataTable$mlResiduals_c3 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
          nXvar + 2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c3 <- as.numeric(crossprod(matrix(mleList$mlParam[(2 * 
          nXvar + 2 * nuZUvar + 2 * nvZVvar + 1):(3 * 
          nXvar + 2 * nuZUvar + 2 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c4 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(3 * 
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
          nXvar + 3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c4 <- as.numeric(crossprod(matrix(mleList$mlParam[(3 * 
          nXvar + 3 * nuZUvar + 3 * nvZVvar + 1):(4 * 
          nXvar + 3 * nuZUvar + 3 * nvZVvar)]), t(Xvar)))
        dataTable$mlResiduals_c5 <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[(4 * 
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * 
          nXvar + 4 * nuZUvar + 4 * nvZVvar)]), t(Xvar)))
        dataTable$mlFitted_c5 <- as.numeric(crossprod(matrix(mleList$mlParam[(4 * 
          nXvar + 4 * nuZUvar + 4 * nvZVvar + 1):(5 * 
          nXvar + 4 * nuZUvar + 4 * nvZVvar)]), t(Xvar)))
      }
    }
  }
  dataTable$logL_OBS <- mleList$mleObj$logL_OBS
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- N
  returnObj$nXvar <- nXvar
  returnObj$nZHvar <- nZHvar
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  # returnObj$olsParam <- olsParam
  # returnObj$olsStder <- olsStder
  # returnObj$olsSigmasq <- olsSigmasq
  # returnObj$olsLoglik <- olsLoglik
  # returnObj$olsSkew <- olsSkew
  # returnObj$olsM3Okay <- olsM3Okay
  if (is.null(start)) {
    returnObj$InitHalf <- mleList$InitHalf
  }
  # returnObj$CoelliM3Test <- CoelliM3Test
  # returnObj$AgostinoTest <- AgostinoTest
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$initStart <- initStart
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$nClasses <- lcmClasses
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
  rm(mleList)
  class(returnObj) <- "lcmcross"
  # print.lcmcross(returnObj)
  return(returnObj)
}

# print for lcmcross ----------

print.lcmcross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call, width.cutoff = 500))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat("Normal-Half Normal Latent Class Stochastic Frontier Model", 
    "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}
