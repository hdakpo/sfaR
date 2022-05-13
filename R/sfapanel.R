# SFA estimation for panel data ---------- shall I create
# ghet or simply use uhet
# instead of modelType maybe use ineffType
sfapanel <- function(formula, muhet, uhet, vhet, logDepVar = TRUE,
  data, idVar = NULL, timeVar = NULL, subset, weights, wscale = TRUE,
  S = 1L, modelType = "bc92a", udist = "hnormal", start = NULL,
  invariance = 2L, method = "bfgs", hessianType = 2L, simType = "halton",
  Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE, seed = 12345,
  itermax = 2000, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06,
  stepmax = 0.1, qac = "marquardt") {
  # panel model check -------
  modelType <- tolower(modelType)
  if (!(modelType %in% c("pl81", "bc92a", "bc92b", "bc92c",
    "k90", "kw05", "cu00"))) {
    stop("Unknown SFA panel model: ", paste(modelType), call. = FALSE)
  }
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
  m <- match(c("formula", "data", "subset", "weights"), names(mc),
    nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$drop.unused.levels <- TRUE
  formula <- interCheckMain(formula = formula)
  if (modelType %in% c("bc92a", "bc92b", "k90", "kw05", "cu00")) {
    muhet <- ~1
    uhet <- ~1
    if (!missing(vhet)) {
      vhet <- clhsCheck_v(formula = vhet)
    } else {
      vhet <- ~1
    }
  } else {
    if (modelType == "pl81") {
      if (!missing(muhet)) {
        muhet <- plhsCheck_mu(formula = muhet)
      } else {
        muhet <- ~1
      }
      if (!missing(uhet)) {
        uhet <- plhsCheck_u(formula = uhet)
      } else {
        uhet <- ~1
      }
      if (!missing(vhet)) {
        vhet <- clhsCheck_v(formula = vhet)
      } else {
        vhet <- ~1
      }
    } else {
      if (modelType == "bc92c") {
        muhet <- ~1
        uhet <- plhsCheck_u_bc92c(formula = uhet)
      }
    }
  }
  formula <- formDist_sfacross(udist = udist, formula = formula,
    muhet = muhet, uhet = uhet, vhet = vhet)
  # Generate required datasets -------
  if (missing(data))
    data <- environment(formula)
  if (!inherits(data, "pdata.frame")) {
    if (is.null(idVar) & is.null(timeVar)) {
      stop("'data' must be of class 'pdata.frame' or arguments 'idVar' & 'timeVar' must be provided",
        call. = FALSE)
    } else {
      if (is.null(idVar) & !is.null(timeVar)) {
        stop("Argument 'idVar' must be provided", call. = FALSE)
      } else {
        if (!is.null(idVar) & is.null(timeVar)) {
          stop("Arguments 'timeVar' must be provided",
          call. = FALSE)
        } else {
          if (!is.null(idVar) & !is.null(timeVar)) {
          checkNames(c(idVar, timeVar), names(data))
          data <- pdata.frame(data, index = c(idVar,
            timeVar))
          }
        }
      }
    }
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
  NT <- nrow(Xvar)
  pindex <- index(data)
  TT <- as.numeric(table(pindex[, 1]))
  N <- length(TT)
  if (NT == 0L) {
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
      wHvar_c <- wHvar/sum(wHvar) * NT
      varw_p <- if (invariance == 1) {
        as.vector(tapply(wHvar, pindex[, 1], function(u) u[1]))
      } else {
        if (invariance == 2) {
          as.vector(tapply(wHvar, pindex[, 1], function(u) u[length(u)]))
        } else {
          if (invariance == 3) {
          as.vector(tapply(wHvar, pindex[, 1], mean))
          }
        }
      }
      wHvar_p <- N * varw_p/sum(varw_p)
    }
  } else {
    wHvar_c <- rep(1, NT)
    wHvar_p <- rep(1, N)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (",
      length(Yvar), ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
  # check for hetero in U (first, last, mean) -------
  if (!(invariance %in% 1:3))
    stop("'invariance' can either be 1, 2 or 3", call. = FALSE)
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula, data = data,
      rhs = 2))
    muHvar_c <- model.matrix(mtmuH, mc)
    muHvar_c <- muHvar_c[validObs, , drop = FALSE]
    nmuZUvar <- ncol(muHvar_c)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar_c <- model.matrix(mtuH, mc)
    uHvar_c <- uHvar_c[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar_c)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar_c <- model.matrix(mtvH, mc)
    vHvar_c <- vHvar_c[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar_c)
    if (invariance == 1) {
      muHvar_p <- apply(muHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          muHvar_p <- apply(muHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], mean)
          })
        }
      }
    }
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar_c <- model.matrix(mtuH, mc)
    uHvar_c <- uHvar_c[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar_c)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar_c <- model.matrix(mtvH, mc)
    vHvar_c <- vHvar_c[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar_c)
    if (invariance == 1) {
      uHvar_p <- apply(uHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
      vHvar_p <- apply(vHvar_c, 2, function(x) {
        tapply(x, pindex[, 1], function(u) u[1])
      })
    } else {
      if (invariance == 2) {
        uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
        vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) u[length(u)])
        })
      } else {
        if (invariance == 3) {
          uHvar_p <- apply(uHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) mean(u))
          })
          vHvar_p <- apply(vHvar_c, 2, function(x) {
          tapply(x, pindex[, 1], function(u) mean(u))
          })
        }
      }
    }
  }
  if (modelType == "bc92a") {
    gHvar <- cbind(eta = unlist(lapply(TT, FUN = function(x) {
      -rev(-(0:(x - 1)))
    })))
    ngZGvar <- dim(gHvar)[2]
  } else {
    if (modelType == "bc92b") {
      gHvar <- cbind(eta1 = unlist(lapply(TT, FUN = function(x) {
        -rev(-(0:(x - 1)))
      })), eta2 = unlist(lapply(TT, FUN = function(x) {
        -rev(-(0:(x - 1)))^2
      })))
      ngZGvar <- dim(gHvar)[2]
    } else {
      if (modelType == "bc92c") {
        gHvar <- uHvar_c
        ngZGvar <- dim(gHvar)[2]
      } else {
        if (modelType == "kw05") {
          gHvar <- cbind(eta = unlist(lapply(TT, FUN = function(x) (seq(1,
          x) - 1))))
          ngZGvar <- dim(gHvar)[2]
        } else {
          if (modelType == "cu00") {
          iHvar <- model.matrix(~-1 + as.factor(data[,
            names(attr(data, "index"))[1]]), rhs = 1)
          nameID <- levels(data[, names(attr(data,
            "index"))[1]])
          colnames(iHvar) <- nameID
          gHvar <- sweep(iHvar, MARGIN = 1, STATS = unlist(lapply(TT,
            FUN = function(x) {
            -rev(-(0:(x - 1)))
            })), FUN = "*")
          ngZGvar <- dim(gHvar)[2]
          }
        }
      }
    }
  }
  # check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Stochastic Production/Profit Frontier, e = v - u"
  } else {
    "Stochastic Cost Frontier, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value",
      call. = FALSE)
  }
  # Number of parameters -------
  if (modelType == "pl81") {
    nParm <- if (udist %in% c("tnormal", "lognormal")) {
      nXvar + nmuZUvar + nuZUvar + nvZVvar
    } else {
      if (udist %in% c("gamma", "weibull", "tslaplace")) {
        nXvar + nuZUvar + nvZVvar + 1
      } else {
        nXvar + nuZUvar + nvZVvar
      }
    }
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05",
      "cu00")) {
      nParm <- if (udist %in% c("tnormal", "lognormal")) {
        nXvar + nmuZUvar + ngZGvar + nvZVvar + 1  # +1 for nuZUvar
      } else {
        if (udist %in% c("gamma", "weibull", "tslaplace")) {
          nXvar + ngZGvar + nvZVvar + 1 + 1  # +1 for nuZUvar
        } else {
          nXvar + ngZGvar + nvZVvar + 1  # +1 for nuZUvar
        }
      }
    } else {
      if (modelType == "k90") {
        nParm <- if (udist == "tnormal") {
          nXvar + nmuZUvar + nuZUvar + nvZVvar + 2
        } else {
          if (udist == "lognormal") {
          nXvar + nmuZUvar + nuZUvar + nvZVvar + 2
          } else {
          if (udist %in% c("gamma", "weibull", "tslaplace")) {
            nXvar + nuZUvar + nvZVvar + 1 + 2
          } else {
            nXvar + nuZUvar + nvZVvar + 2
          }
          }
        }
      }
    }
  }
  # checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ",
        nParm, " parameters", call. = FALSE)
    }
  }
  if (nParm > NT) {
    stop("Model has more parameters than observations", call. = FALSE)
  }
  # check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm",
    "cg", "sann", "sr1", "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ",
      paste(method), call. = FALSE)
  }
  # check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L,
    2L))) {
    stop("argument 'hessianType' must equal either 1 or 2",
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
    cat("Initialization of", Nsim, simDist, "draws per observation/cross-section ...\n")
    FiMat_N <- drawMat(N = N, Nsim = Nsim, simType = simType,
      prime = prime, burn = burn + 1, antithetics = antithetics,
      seed = seed)
    FiMat_NT <- drawMat(N = NT, Nsim = Nsim, simType = simType,
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
      lm(Yvar ~ ., data = as.data.frame(Xvar[, -1]), weights = wHvar_c)
    }
  } else {
    lm(Yvar ~ -1 + ., data = as.data.frame(Xvar), weights = wHvar_c)
  }
  if (any(is.na(olsRes$coefficients))) {
    stop("at least one of the OLS coefficients is NA: ",
      paste(colnames(Xvar)[is.na(olsRes$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
      call. = FALSE)
  }
  olsParam <- c(olsRes$coefficients)
  dataTable <- data[validObs, 1:2]
  dataTable <- as_tibble(cbind(dataTable, data[, all.vars(terms(formula))],
    weights = wHvar_c))
  dataTable <- mutate(dataTable, olsResiduals = residuals(olsRes),
    olsFitted = fitted(olsRes))
  olsSkew <- skewness(dataTable[["olsResiduals"]])
  if (S * olsSkew > 0) {
    warning("The residuals of the OLS are", if (S == 1) {
      " right"
    } else {
      " left"
    }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
      call. = FALSE)
  }
  # Step 2: MLE arguments -------
  if (modelType %in% c("pl81", "k90")) {
    FunArgs <- if (udist == "tnormal") {
      list(start = start, olsParam = olsParam, dataTable = dataTable,
        nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c,
        vHvar_c = vHvar_c, muHvar_p = muHvar_p, uHvar_p = uHvar_p,
        vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
        S = S, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
        pindex = pindex, TT = TT, method = method, printInfo = printInfo,
        itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType,
        qac = qac)
    } else {
      if (udist == "lognormal") {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c,
          vHvar_c = vHvar_c, muHvar_p = muHvar_p, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
          wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S,
          pindex = pindex, TT = TT, N = N, NT = NT, FiMat_N = FiMat_N,
          FiMat_NT = FiMat_NT, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol,
          gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      } else {
        if (udist %in% c("gamma", "weibull")) {
          list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
          wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S,
          pindex = pindex, TT = TT, N = N, NT = NT,
          FiMat_N = FiMat_N, FiMat_NT = FiMat_NT, method = method,
          printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
        } else {
          list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar_c = uHvar_c, vHvar_c = vHvar_c, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
          wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S,
          pindex = pindex, TT = TT, method = method,
          printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
        }
      }
    }
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05",
      "cu00")) {
      FunArgs <- if (udist == "tnormal") {
        list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c,
          vHvar_c = vHvar_c, muHvar_p = muHvar_p, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
          modelType = modelType, ngZGvar = ngZGvar, gHvar = gHvar,
          S = S, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
          pindex = pindex, TT = TT, method = method,
          printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType,
          qac = qac)
      } else {
        if (udist == "lognormal") {
          list(start = start, olsParam = olsParam, dataTable = dataTable,
          nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar_c = muHvar_c, uHvar_c = uHvar_c,
          vHvar_c = vHvar_c, muHvar_p = muHvar_p, uHvar_p = uHvar_p,
          vHvar_p = vHvar_p, Yvar = Yvar, Xvar = Xvar,
          wHvar_c = wHvar_c, wHvar_p = wHvar_p, S = S,
          modelType = modelType, ngZGvar = ngZGvar,
          gHvar = gHvar, pindex = pindex, TT = TT,
          N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT,
          method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
        } else {
          if (udist %in% c("gamma", "weibull")) {
          list(start = start, olsParam = olsParam,
            dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar,
            nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
            uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
            Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
            S = S, modelType = modelType, ngZGvar = ngZGvar,
            gHvar = gHvar, pindex = pindex, TT = TT,
            N = N, NT = NT, FiMat_N = FiMat_N, FiMat_NT = FiMat_NT,
            method = method, printInfo = printInfo,
            itermax = itermax, stepmax = stepmax, tol = tol,
            gradtol = gradtol, hessianType = hessianType,
            qac = qac)
          } else {
          list(start = start, olsParam = olsParam,
            dataTable = dataTable, nXvar = nXvar, nuZUvar = nuZUvar,
            nvZVvar = nvZVvar, uHvar_c = uHvar_c, vHvar_c = vHvar_c,
            uHvar_p = uHvar_p, vHvar_p = vHvar_p, Yvar = Yvar,
            Xvar = Xvar, wHvar_c = wHvar_c, wHvar_p = wHvar_p,
            S = S, modelType = modelType, ngZGvar = ngZGvar,
            gHvar = gHvar, pindex = pindex, TT = TT,
            method = method, printInfo = printInfo,
            itermax = itermax, stepmax = stepmax, tol = tol,
            gradtol = gradtol, hessianType = hessianType,
            qac = qac)
          }
        }
      }
    }
  }
  ## MLE run -------
  if (modelType == "pl81") {
    mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_pl81,
      FunArgs), tnormal = do.call(truncnormAlgOpt_pl81,
      FunArgs), exponential = do.call(exponormAlgOpt_pl81,
      FunArgs), rayleigh = do.call(raynormAlgOpt_pl81,
      FunArgs), genexponential = do.call(genexponormAlgOpt_pl81,
      FunArgs), tslaplace = do.call(tslnormAlgOpt_pl81,
      FunArgs)), error = function(e) e)
  } else {
    if (modelType %in% c("bc92a", "bc92b", "bc92c", "kw05",
      "cu00")) {
      mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_bc92I,
        FunArgs), tnormal = do.call(truncnormAlgOpt_bc92I,
        FunArgs), exponential = do.call(exponormAlgOpt_bc92I,
        FunArgs), rayleigh = do.call(raynormAlgOpt_bc92I,
        FunArgs), genexponential = do.call(genexponormAlgOpt_bc92I,
        FunArgs), tslaplace = do.call(tslnormAlgOpt_bc92I,
        FunArgs)), error = function(e) e)
    } else {
      if (modelType == "k90") {
        mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt_k90,
          FunArgs), tnormal = do.call(truncnormAlgOpt_k90,
          FunArgs), exponential = do.call(exponormAlgOpt_k90,
          FunArgs), rayleigh = do.call(raynormAlgOpt_k90,
          FunArgs), genexponential = do.call(genexponormAlgOpt_k90,
          FunArgs), tslaplace = do.call(tslnormAlgOpt_k90,
          FunArgs)), error = function(e) e)
      }
    }
  }
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
    names(mleList$startVal) <- fName_mu_sfapanel(Xvar = Xvar,
      udist = udist, muHvar = muHvar_p, uHvar = uHvar_p,
      vHvar = vHvar_p, modelType = modelType)
  } else {
    names(mleList$startVal) <- fName_uv_sfapanel(Xvar = Xvar,
      udist = udist, uHvar = uHvar_p, vHvar = vHvar_p,
      modelType = modelType)
  }
  names(mleList$mlParam) <- names(mleList$startVal)
  rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
  names(mleList$gradient) <- names(mleList$mlParam)
  colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
  # Return object -------
  mlDate <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M")
  dataTable$mlResiduals <- Yvar - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  dataTable$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
    t(Xvar)))
  datlogL <- data.frame(levels(pindex[, 1]), logL_OBS = mleList$mleObj$logL_OBS/TT)
  names(datlogL)[1] <- names(pindex)[1]
  dataTable <- merge(dataTable, datlogL, by = names(pindex)[1])
  dataTable <- pdata.frame(dataTable, names(dataTable)[1:2])
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- NT
  returnObj$Nid <- N
  returnObj$Vtime <- TT  # comment
  returnObj$Ntime <- mean(TT)  # comment
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuZUvar <- nmuZUvar
  }
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  # ngZGvar must also be accounted for
  returnObj$logDepVar <- logDepVar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- mleList$startVal
  returnObj$dataTable <- dataTable
  returnObj$modelType <- modelType
  returnObj$invariance <- invariance
  if (is.null(start)) {
    if (udist == "hnormal") {
      returnObj$initHalf <- mleList$initHalf
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
  returnObj$isWeights <- !all.equal(wHvar_p, rep(1, N))
  returnObj$optType <- mleList$type
  returnObj$nIter <- mleList$nIter
  returnObj$optStatus <- mleList$status
  returnObj$startLoglik <- mleList$startLoglik
  returnObj$mlLoglik <- mleList$mleLoglik
  returnObj$mlParam <- mleList$mlParam
  returnObj$gradient <- mleList$gradient
  datGradL <- mleList$mleObj$gradL_OBS
  datGradL <- cbind(levels(pindex[, 1]), datGradL)
  colnames(datGradL)[1] <- names(pindex)[1]
  returnObj$gradL_OBS <- datGradL
  returnObj$gradientNorm <- sqrt(sum(mleList$gradient^2))
  returnObj$invHessian <- mleList$invHessian
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  returnObj$mlDate <- mlDate
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat_N
  }
  rm(mleList)
  class(returnObj) <- "sfapanel"
  return(returnObj)
}

# print for sfapanel ----------
print.sfapanel <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(sfadist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}
