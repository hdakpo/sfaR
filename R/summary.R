# summary for sfacross ----------

summary.sfacross <- function(object, grad = FALSE, ci = FALSE, ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE
    )
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE
    )
  }
  object$AIC <- -2 * object$mlLoglik + 2 * object$nParm
  object$BIC <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  object$HQIC <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) * object$nParm
  if (object$udist == "tnormal") {
    if (object$scaling) {
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        (object$nuZUvar - 1))]
      tau <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 1]
      cu <- object$mlParam[object$nXvar + (object$nuZUvar -
        1) + 2]
      phi <- object$mlParam[(object$nXvar + (object$nuZUvar -
        1) + 2 + 1):(object$nXvar + (object$nuZUvar -
        1) + 2 + object$nvZVvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    } else {
      omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nmuZUvar)]
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    }
  } else {
    if (object$udist == "lognormal") {
      omega <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nmuZUvar)]
      delta <- object$mlParam[(object$nXvar + object$nmuZUvar +
        1):(object$nXvar + object$nmuZUvar + object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nmuZUvar +
        object$nuZUvar + 1):(object$nXvar + object$nmuZUvar +
        object$nuZUvar + object$nvZVvar)]
      muHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 4
      )
    } else {
      delta <- object$mlParam[(object$nXvar + 1):(object$nXvar +
        object$nuZUvar)]
      phi <- object$mlParam[(object$nXvar + object$nuZUvar +
        1):(object$nXvar + object$nuZUvar + object$nvZVvar)]
      uHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 2
      )
      vHvar <- model.matrix(object$formula,
        data = object$dataTable,
        rhs = 3
      )
    }
  }
  mu <- if (object$udist == "tnormal") {
    if (object$scaling) {
      exp(as.numeric(crossprod(matrix(delta), t(uHvar[
        ,
        -1
      ])))) * tau
    } else {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    }
  } else {
    if (object$udist == "lognormal") {
      as.numeric(crossprod(matrix(omega), t(muHvar)))
    } else {
      NULL
    }
  }
  P <- if (object$udist == "gamma") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  k <- if (object$udist == "weibull") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  lambda <- if (object$udist == "tslaplace") {
    object$mlParam[object$nXvar + object$nuZUvar + object$nvZVvar +
      1]
  } else {
    NULL
  }
  Wu <- if (object$udist == "tnormal" & object$scaling == TRUE) {
    cu + 2 * as.numeric(crossprod(matrix(delta), t(uHvar[
      ,
      -1
    ])))
  } else {
    as.numeric(crossprod(matrix(delta), t(uHvar)))
  }
  Wv <- as.numeric(crossprod(matrix(phi), t(vHvar)))
  object$sigmavSq <- mean(exp(Wv))
  object$sigmauSq <- mean(exp(Wu))
  object$Varu <- varuFun(
    object = object, mu = mu, P = P, k = k,
    lambda = lambda
  )
  if (object$udist == "uniform") {
    object$THETA <- sqrt(12 * object$sigmauSq)
  }
  object$Eu <- euFun(
    object = object, mu = mu, P = P, k = k,
    lambda = lambda
  )
  object$Expu <- eExpuFun(
    object = object, mu = mu, P = P,
    k = k, lambda = lambda
  )
  # OLS estimates and stder, p-values, CI
  dfOLS <- object$Nobs - object$nXvar
  if (ci) {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 6)
    colnames(olsRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "t value", "Pr(>|t|)"
    )
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1] - qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 4] <- olsRes[, 1] + qt(0.975, df = dfOLS) *
      olsRes[, 2]
    olsRes[, 5] <- olsRes[, 1] / olsRes[, 2]
    olsRes[, 6] <- 2 * pt(-abs(olsRes[, 5]), df = dfOLS)
  } else {
    olsRes <- matrix(nrow = object$nXvar + 1, ncol = 4)
    colnames(olsRes) <- c(
      "Coefficient", "Std. Error", "t value",
      "Pr(>|t|)"
    )
    olsRes[, 1] <- c(object$olsParam, object$olsSigmasq)
    olsRes[, 2] <- c(object$olsStder, NA)
    olsRes[, 3] <- olsRes[, 1] / olsRes[, 2]
    olsRes[, 4] <- 2 * pt(-abs(olsRes[, 3]), df = dfOLS)
  }
  row.names(olsRes) <- c(names(object$olsParam), "sigmaSq")
  object$olsRes <- olsRes
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mlRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mlRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)"
    )
    mlRes[, 1] <- object$mlParam
    mlRes[, 2] <- sqrt(diag(object$invHessian))
    mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[
      ,
      2
    ]
    mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[
      ,
      2
    ]
    mlRes[, 5] <- object$gradient
    mlRes[, 6] <- mlRes[, 1] / mlRes[, 2]
    mlRes[, 7] <- 2 * pnorm(-abs(mlRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mlRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mlRes) <- c(
        "Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)"
      )
      mlRes[, 1] <- object$mlParam
      mlRes[, 2] <- sqrt(diag(object$invHessian))
      mlRes[, 3] <- object$gradient
      mlRes[, 4] <- mlRes[, 1] / mlRes[, 2]
      mlRes[, 5] <- 2 * pnorm(-abs(mlRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mlRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mlRes) <- c(
          "Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)"
        )
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[
          ,
          2
        ]
        mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[
          ,
          2
        ]
        mlRes[, 5] <- mlRes[, 1] / mlRes[, 2]
        mlRes[, 6] <- 2 * pnorm(-abs(mlRes[, 5]))
      } else {
        mlRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mlRes) <- c(
          "Coefficient", "Std. Error",
          "z value", "Pr(>|z|)"
        )
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] / mlRes[, 2]
        mlRes[, 4] <- 2 * pnorm(-abs(mlRes[, 3]))
      }
    }
  }
  row.names(mlRes) <- names(object$startVal)
  object$mlRes <- mlRes
  object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nXvar - object$nvZVvar
  class(object) <- "summary.sfacross"
  return(object)
}

# print summary for sfacross ----------

print.summary.sfacross <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  mlRes <- x$mlRes
  if (dim(mlRes)[2] == 4) {
    mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
      digits = digits, format = "f"
    ))
    mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
      digits = digits, format = "f"
    ))
    mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
      digits = digits, format = "f"
    ))
    mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
      digits = digits, format = "e"
    ))
  } else {
    if (dim(mlRes)[2] == 5) {
      mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
        digits = digits, format = "f"
      ))
      mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
        digits = digits, format = "f"
      ))
      mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
        digits = digits, format = "e"
      ))
      mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
        digits = digits, format = "f"
      ))
      mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
        digits = digits, format = "e"
      ))
    } else {
      if (dim(mlRes)[2] == 6) {
        mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1], 
        digits = digits, format = "f"))
        mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2], 
        digits = digits, format = "f"))
        mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3], 
        digits = digits, format = "f"))
        mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4], 
        digits = digits, format = "f"))
        mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5], 
        digits = digits, format = "f"))
        mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6], 
        digits = digits, format = "e"))
      } else {
        if (dim(mlRes)[2] == 7) {
          mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
          digits = digits, format = "f"))
          mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
          digits = digits, format = "f"))
          mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
          digits = digits, format = "f"))
          mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
          digits = digits, format = "f"))
          mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
          digits = digits, format = "e"))
          mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
          digits = digits, format = "f"))
          mlRes[, 7] <- as.numeric(formatC(x$mlRes[, 7],
          digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mlRes) <- formatC(row.names(mlRes),
    width = max(nchar(row.names(mlRes))),
    flag = "-"
  )
  mlRes1 <- mlRes[1:x$nXvar, , drop = FALSE]
  if (x$udist == "tnormal") {
    if (x$scaling) {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + (x$nuZUvar -
        1)), , drop = FALSE]
      mlRes3 <- mlRes[x$nXvar + (x$nuZUvar - 1) + 1, , drop = FALSE]
      mlRes4 <- mlRes[x$nXvar + (x$nuZUvar - 1) + 2, , drop = FALSE]
      mlRes5 <- mlRes[(x$nXvar + (x$nuZUvar - 1) + 2 +
        1):(x$nXvar + (x$nuZUvar - 1) + 2 + x$nvZVvar), , drop = FALSE]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar), , drop = FALSE]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar), , drop = FALSE]
    }
  } else {
    if (x$udist == "lognormal") {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nmuZUvar), ,
        drop = FALSE
      ]
      mlRes3 <- mlRes[(x$nXvar + x$nmuZUvar + 1):(x$nXvar +
        x$nmuZUvar + x$nuZUvar), , drop = FALSE]
      mlRes4 <- mlRes[(x$nXvar + x$nmuZUvar + x$nuZUvar +
        1):(x$nXvar + x$nmuZUvar + x$nuZUvar + x$nvZVvar), ,
      drop = FALSE
      ]
    } else {
      mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), ,
        drop = FALSE
      ]
      mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar +
        x$nuZUvar + x$nvZVvar), , drop = FALSE]
      if (x$udist %in% c("gamma", "weibull", "tslaplace")) {
        mlRes4 <- mlRes[x$nXvar + x$nuZUvar + x$nvZVvar +
          1, , drop = FALSE]
      }
    }
  }
  lengthSum <- nchar(sfadist(x$udist)) + 10
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(sfadist(x$udist), "\n")
  cat(
    "Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
      nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n"
  )
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
  collapse = ""
  ), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"
  ), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
      digits = digits, format = "e"
    ))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"
  ), "\n")
  cat(
    "Estimation based on:", paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
      3), collapse = ""), "N = ", x$Nobs, "and K = ", x$nParm,
    "\n"
  )
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC / x$Nobs,
      digits = 3,
      format = "f"
    )) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"
  ), "AIC/N  = ", formatC(x$AIC / x$Nobs,
    digits = 3, format = "f"
  ), "\n")
  cat(
    paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
      digits = 1, format = "f"
    )) - nchar("BIC/N  = ") - nchar(formatC(x$BIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(
    paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
      digits = 1, format = "f"
    )) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Variances: Sigma-squared(v)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(x$sigmavSq,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$sigmavSq,
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma(v)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(v)   = ") - nchar(formatC(sqrt(x$sigmavSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmavSq),
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma-squared(u)   = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(x$sigmauSq,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$sigmauSq,
    digits = digits, format = "f"
  ), "\n")
  cat("           Sigma(u)           = ", paste0(rep(" ", lengthSum -
    nchar("Variances: Sigma-squared(u)   = ") - nchar(formatC(sqrt(x$sigmauSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmauSq),
    digits = digits, format = "f"
  ), "\n")
  cat(
    "Sigma = Sqrt[(s^2(u)+s^2(v))] = ", paste0(rep(" ", lengthSum -
      nchar("Sigma = Sqrt[(s^2(u)+s^2(v))] = ") - nchar(formatC(sqrt(x$sigmavSq +
        x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(sqrt(x$sigmavSq + x$sigmauSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat(
    "Gamma = sigma(u)^2/sigma^2    = ", paste0(rep(" ", lengthSum -
      nchar("Gamma = sigma(u)^2/sigma^2    = ") - nchar(formatC(x$sigmauSq / (x$sigmavSq +
        x$sigmauSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$sigmauSq / (x$sigmavSq + x$sigmauSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat("Lambda = sigma(u)/sigma(v)    = ", paste0(rep(" ", lengthSum -
    nchar("Lambda = sigma(u)/sigma(v)    = ") - nchar(formatC(sqrt(x$sigmauSq / x$sigmavSq),
      digits = digits, format = "f"
    ))), collapse = ""), formatC(sqrt(x$sigmauSq / x$sigmavSq),
    digits = digits, format = "f"
  ), "\n")
  cat(
    "Var[u]/{Var[u]+Var[v]}        = ", paste0(rep(" ", lengthSum -
      nchar("Var[u]/{Var[u]+Var[v]}        = ") - nchar(formatC(x$Varu / (x$Varu +
        x$sigmavSq), digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu / (x$Varu + x$sigmavSq),
      digits = digits,
      format = "f"
    ), "\n"
  )
  cat(
    "Var[e]                        = ", paste0(rep(" ", lengthSum -
      nchar("Var[e]                        = ") - nchar(formatC(x$Varu +
        x$sigmavSq, digits = digits, format = "f"))), collapse = ""),
    formatC(x$Varu + x$sigmavSq, digits = digits, format = "f"),
    "\n"
  )
  if (x$udist == "uniform") {
    cat("THETA                         = ", paste0(rep(
      " ",
      lengthSum - nchar("THETA                         = ") -
        nchar(formatC(x$theta, digits = digits, format = "f"))
    ),
    collapse = ""
    ), formatC(x$theta,
      digits = digits,
      format = "f"
    ), "\n")
  }
  if (x$nuZUvar > 1 || x$nvZVvar > 1) {
    cat("Variances averaged over observations \n")
  }
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat("Average inefficiency E[u]     = ", paste0(rep(" ", lengthSum -
    nchar("Average inefficiency E[u]     = ") - nchar(formatC(x$Eu,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$Eu,
    digits = digits, format = "f"
  ), "\n")
  cat("Average efficiency E[exp(-u)] = ", paste0(rep(" ", lengthSum -
    nchar("Average efficiency E[exp(-u)] = ") - nchar(formatC(x$Expu,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$Expu,
    digits = digits, format = "f"
  ), "\n")
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("-----[ Tests vs. No Inefficiency ]-----\n")
  cat("Likelihood Ratio Test of Inefficiency\n")
  cat(
    "Deg. freedom for inefficiency model", paste0(rep(
      " ",
      lengthSum - nchar("Deg. freedom for inefficiency model") -
        nchar(formatC(x$df, digits = digits, format = "d"))
    ),
    collapse = ""
    ), formatC(x$df, digits = digits, format = "d"),
    "\n"
  )
  cat("Log Likelihood for OLS Log(H0) = ", paste0(rep(
    " ",
    lengthSum - nchar("Log Likelihood for OLS Log(H0) = ") -
      nchar(formatC(x$olsLoglik, digits = digits, format = "f"))
  ),
  collapse = ""
  ), formatC(x$olsLoglik,
    digits = digits,
    format = "f"
  ), "\n")
  cat("LR statistic: ", "\n")
  cat(
    "Chisq = 2*[LogL(H0)-LogL(H1)]  = ", paste0(rep(
      " ",
      lengthSum - nchar("Chisq = 2*[LogL(H0)-LogL(H1)]  = ") -
        nchar(formatC(x$chisq, digits = digits, format = "f"))
    ),
    collapse = ""
    ), formatC(x$chisq, digits = digits, format = "f"),
    "\n"
  )
  cat("Kodde-Palm C*:       95%:", formatC(qchibarsq(0.95,
    df = x$df
  ), digits = digits, format = "f"), paste0(rep(
    " ",
    lengthSum - nchar("Kodde-Palm C*:       95%:") - nchar(formatC(qchibarsq(0.95,
      df = x$df
    ), digits = digits, format = "f")) - nchar(formatC(qchibarsq(0.99,
      df = x$df
    ), digits = digits, format = "f")) - nchar("99%") -
      3
  ), collapse = ""), "99%:", formatC(qchibarsq(0.99,
    df = x$df
  ), digits = digits, format = "f"), "\n")
  cat("Coelli (1995) skewness test on OLS residuals\n")
  cat("M3T                            = ", paste0(rep(
    " ",
    lengthSum - nchar("M3T                            = ") -
      nchar(formatC(x$CoelliM3Test[1],
        digits = digits,
        format = "f"
      ))
  ), collapse = ""), formatC(x$CoelliM3Test[1],
    digits = digits, format = "f"
  ), "\n")
  cat("final maximum likelihood estimates \n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA", width = lengthSum +
    2 + switch(dimCoefTable, `4` = 18, `5` = 31, `6` = 43,
      `7` = 57
    )), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes1, P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  if (x$udist == "tnormal") {
    if (x$scaling) {
      cat(centerText("Scaling property parameters", width = lengthSum +
        2 + switch(dimCoefTable, `4` = 18, `5` = 31,
          `6` = 43, `7` = 57
        )), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes4,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes5,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes4,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    }
  } else {
    if (x$udist == "lognormal") {
      cat(centerText("Location parameter [offset mu] in u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes2,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes3,
        P.values = TRUE, digits = digits,
        signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (two-sided error)",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes4,
        P.values = TRUE, digits = digits,
        signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      if (x$udist == "gamma") {
        cat(
          centerText("Parameter in variance of u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes2,
          P.values = TRUE, digits = digits,
          signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (two-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes3,
          P.values = TRUE, digits = digits,
          signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Location parameter P in u (one-sided error)",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes4,
          P.values = TRUE, digits = digits,
          signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
      } else {
        if (x$udist == "weibull") {
          cat(
            centerText("Parameter in variance of u (one-sided error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes2,
            P.values = TRUE, digits = digits,
            signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (two-sided error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes3,
            P.values = TRUE, digits = digits,
            signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Shape parameter k in u (one-sided error)",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes4,
            P.values = TRUE, digits = digits,
            signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
        } else {
          if (x$udist == "tslaplace") {
            cat(
              centerText("Parameter in variance of u (one-sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mlRes2,
              P.values = TRUE, digits = digits,
              signif.legend = FALSE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Parameters in variance of v (two-sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mlRes3,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Skewness parameter 'lambda' in u (one-sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            printCoefmat(mlRes4,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
          } else {
            cat(
              centerText("Parameter in variance of u (one-sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mlRes2,
              P.values = TRUE, digits = digits,
              signif.legend = FALSE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            cat(
              centerText("Parameters in variance of v (two-sided error)",
                width = lengthSum + 2 + switch(dimCoefTable,
                  `4` = 18, `5` = 31, `6` = 43, `7` = 57
                )
              ),
              "\n"
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
            printCoefmat(mlRes3,
              P.values = TRUE, digits = digits,
              signif.legend = TRUE
            )
            cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )),
            collapse = ""
            ), "\n")
          }
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  invisible(x)
}

# summary for lcmcross ----------

summary.lcmcross <- function(object, grad = FALSE, ci = FALSE, ...) {
  if (length(grad) != 1 || !is.logical(grad[1])) {
    stop("argument 'grad' must be a single logical value",
      call. = FALSE
    )
  }
  if (length(ci) != 1 || !is.logical(ci[1])) {
    stop("argument 'ci' must be a single logical value",
      call. = FALSE
    )
  }
  object$AIC <- -2 * object$mlLoglik + 2 * object$nParm
  object$BIC <- -2 * object$mlLoglik + log(object$Nobs) * object$nParm
  object$HQIC <- -2 * object$mlLoglik + 2 * log(log(object$Nobs)) * object$nParm
  # MLE estimates and stder, p-values, CI, Gradient
  if (grad && ci) {
    mlRes <- matrix(nrow = object$nParm, ncol = 7)
    colnames(mlRes) <- c(
      "Coefficient", "Std. Error", "binf",
      "bsup", "gradient", "z-value", "Pr(>|z|)"
    )
    mlRes[, 1] <- object$mlParam
    mlRes[, 2] <- sqrt(diag(object$invHessian))
    mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[
      ,
      2
    ]
    mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[
      ,
      2
    ]
    mlRes[, 5] <- object$gradient
    mlRes[, 6] <- mlRes[, 1] / mlRes[, 2]
    mlRes[, 7] <- 2 * pnorm(-abs(mlRes[, 6]))
  } else {
    if (grad == TRUE && ci == FALSE) {
      mlRes <- matrix(nrow = object$nParm, ncol = 5)
      colnames(mlRes) <- c(
        "Coefficient", "Std. Error",
        "gradient", "z-value", "Pr(>|z|)"
      )
      mlRes[, 1] <- object$mlParam
      mlRes[, 2] <- sqrt(diag(object$invHessian))
      mlRes[, 3] <- object$gradient
      mlRes[, 4] <- mlRes[, 1] / mlRes[, 2]
      mlRes[, 5] <- 2 * pnorm(-abs(mlRes[, 4]))
    } else {
      if (grad == FALSE && ci == TRUE) {
        mlRes <- matrix(nrow = object$nParm, ncol = 6)
        colnames(mlRes) <- c(
          "Coefficient", "Std. Error",
          "binf", "bsup", "z-value", "Pr(>|z|)"
        )
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] - qnorm(0.975) * mlRes[
          ,
          2
        ]
        mlRes[, 4] <- mlRes[, 1] + qnorm(0.975) * mlRes[
          ,
          2
        ]
        mlRes[, 5] <- mlRes[, 1] / mlRes[, 2]
        mlRes[, 6] <- 2 * pnorm(-abs(mlRes[, 5]))
      } else {
        mlRes <- matrix(nrow = object$nParm, ncol = 4)
        colnames(mlRes) <- c(
          "Coefficient", "Std. Error",
          "z value", "Pr(>|z|)"
        )
        mlRes[, 1] <- object$mlParam
        mlRes[, 2] <- sqrt(diag(object$invHessian))
        mlRes[, 3] <- mlRes[, 1] / mlRes[, 2]
        mlRes[, 4] <- 2 * pnorm(-abs(mlRes[, 3]))
      }
    }
  }
  row.names(mlRes) <- names(object$startVal)
  object$mlRes <- mlRes
  # object$chisq <- 2 * (object$mlLoglik - object$olsLoglik)
  object$df <- object$nParm - object$nClasses * object$nXvar -
    object$nClasses * object$nvZVvar - object$nZHvar *
    (object$nClasses - 1)
  class(object) <- "summary.lcmcross"
  return(object)
}

# print summary for lcmcross ----------

print.summary.lcmcross <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  mlRes <- x$mlRes
  if (dim(mlRes)[2] == 4) {
    mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
      digits = digits, format = "f"
    ))
    mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
      digits = digits, format = "f"
    ))
    mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
      digits = digits, format = "f"
    ))
    mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
      digits = digits, format = "e"
    ))
  } else {
    if (dim(mlRes)[2] == 5) {
      mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
        digits = digits, format = "f"
      ))
      mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
        digits = digits, format = "f"
      ))
      mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
        digits = digits, format = "e"
      ))
      mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
        digits = digits, format = "f"
      ))
      mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
        digits = digits, format = "e"
      ))
    } else {
      if (dim(mlRes)[2] == 6) {
        mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
        digits = digits, format = "f"))
        mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
        digits = digits, format = "f"))
        mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
        digits = digits, format = "f"))
        mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
        digits = digits, format = "f"))
        mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
        digits = digits, format = "f"))
        mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
        digits = digits, format = "e"))
      } else {
        if (dim(mlRes)[2] == 7) {
          mlRes[, 1] <- as.numeric(formatC(x$mlRes[, 1],
          digits = digits, format = "f"))
          mlRes[, 2] <- as.numeric(formatC(x$mlRes[, 2],
          digits = digits, format = "f"))
          mlRes[, 3] <- as.numeric(formatC(x$mlRes[, 3],
          digits = digits, format = "f"))
          mlRes[, 4] <- as.numeric(formatC(x$mlRes[, 4],
          digits = digits, format = "f"))
          mlRes[, 5] <- as.numeric(formatC(x$mlRes[, 5],
          digits = digits, format = "e"))
          mlRes[, 6] <- as.numeric(formatC(x$mlRes[, 6],
          digits = digits, format = "f"))
          mlRes[, 7] <- as.numeric(formatC(x$mlRes[, 7],
          digits = digits, format = "e"))
        }
      }
    }
  }
  row.names(mlRes) <- formatC(row.names(mlRes),
    width = max(nchar(row.names(mlRes))),
    flag = "-"
  )
  mlRes1 <- mlRes[1:x$nXvar, ]
  mlRes2 <- mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), , drop = FALSE]
  mlRes3 <- mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar + x$nuZUvar +
    x$nvZVvar), , drop = FALSE]
  sfaModel <- "Normal-Half Normal Latent Class Stochastic Frontier Model"
  lengthSum <- nchar(sfaModel) # + 10
  dimCoefTable <- as.character(dim(x$mlRes)[2])
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(sfaModel, "\n")
  cat(
    "Dependent Variable:", paste0(rep(" ", lengthSum - nchar("Dependent Variable:") -
      nchar(paste0(attr(x$formula, "lhs")))), collapse = ""),
    paste0(attr(x$formula, "lhs")), "\n"
  )
  cat("Log likelihood solver:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood solver:") - nchar(x$optType)),
  collapse = ""
  ), x$optType, "\n")
  cat("Log likelihood iter:", paste0(rep(" ", lengthSum - nchar("Log likelihood iter:") -
    nchar(x$nIter)), collapse = ""), x$nIter, "\n")
  cat("Log likelihood value:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood value:") - nchar(formatC(x$mlLoglik,
      digits = digits, format = "f"
    ))), collapse = ""), formatC(x$mlLoglik,
    digits = digits, format = "f"
  ), "\n")
  cat("Log likelihood gradient norm:", paste0(rep(" ", lengthSum -
    nchar("Log likelihood gradient norm:") - nchar(formatC(x$gradientNorm,
      digits = digits, format = "e"
    ))), collapse = ""), formatC(x$gradientNorm,
    digits = digits, format = "e"
  ), "\n")
  cat(
    "Estimation based on:", paste0(rep(" ", lengthSum - nchar("Estimation based on:") -
      nchar(x$Nobs) - nchar(x$nParm) - nchar("N = ") - nchar("and K = ") -
      3), collapse = ""), "N = ", x$Nobs, "and K = ", x$nParm,
    "\n"
  )
  cat("Inf. Cr:", paste0(rep(" ", lengthSum - nchar("Inf. Cr:") -
    nchar("AIC  = ") - nchar(formatC(x$AIC, digits = 1, format = "f")) -
    nchar("AIC/N  = ") - nchar(formatC(x$AIC / x$Nobs,
      digits = 3,
      format = "f"
    )) - 3), collapse = ""), "AIC  = ", formatC(x$AIC,
    digits = 1, format = "f"
  ), "AIC/N  = ", formatC(x$AIC / x$Nobs,
    digits = 3, format = "f"
  ), "\n")
  cat(
    paste0(rep(" ", lengthSum - nchar("BIC  = ") - nchar(formatC(x$BIC,
      digits = 1, format = "f"
    )) - nchar("BIC/N  = ") - nchar(formatC(x$BIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "BIC  = ",
    formatC(x$BIC, digits = 1, format = "f"), "BIC/N  = ",
    formatC(x$BIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(
    paste0(rep(" ", lengthSum - nchar("HQIC = ") - nchar(formatC(x$HQIC,
      digits = 1, format = "f"
    )) - nchar("HQIC/N = ") - nchar(formatC(x$HQIC / x$Nobs,
      digits = 3, format = "f"
    )) - 2), collapse = ""), "HQIC = ",
    formatC(x$HQIC, digits = 1, format = "f"), "HQIC/N = ",
    formatC(x$HQIC / x$Nobs, digits = 3, format = "f"), "\n"
  )
  cat(paste0(rep("-", lengthSum + 2), collapse = ""), "\n")
  cat(x$typeSfa, "\n")
  cat("Latent class model with", x$nClasses, "latent classes \n")
  cat("final maximum likelihood estimates \n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[1:x$nXvar, , drop = FALSE],
    P.values = TRUE,
    digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[(x$nXvar + 1):(x$nXvar + x$nuZUvar), ,
    drop = FALSE
  ], P.values = TRUE, digits = digits, signif.legend = FALSE)
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 1",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + 1):(x$nXvar + x$nuZUvar +
    x$nvZVvar), , drop = FALSE],
  P.values = TRUE, digits = digits,
  signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Deterministic Component of SFA for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[(x$nXvar + x$nuZUvar + x$nvZVvar + 1):(2 *
    x$nXvar + x$nuZUvar + x$nvZVvar), , drop = FALSE],
  P.values = TRUE,
  digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameter in variance of u (one-sided error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[(2 * x$nXvar + x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar), , drop = FALSE],
  P.values = TRUE, digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  cat(centerText("Parameters in variance of v (two-sided error) for latent class 2",
    width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
      `5` = 31, `6` = 43, `7` = 57
    )
  ), "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + x$nvZVvar +
    1):(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar), , drop = FALSE],
  P.values = TRUE, digits = digits, signif.legend = FALSE
  )
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  if (x$nClasses == 2) {
    cat(centerText("Estimated prior probabilities for class membership",
      width = lengthSum + 2 + switch(dimCoefTable, `4` = 18,
        `5` = 31, `6` = 43, `7` = 57
      )
    ), "\n")
    cat(
      paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57
      )), collapse = ""),
      "\n"
    )
    printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar + 2 *
      x$nvZVvar + 1):(2 * x$nXvar + 2 * x$nuZUvar + 2 * x$nvZVvar +
      x$nZHvar), , drop = FALSE],
    P.values = TRUE, digits = digits,
    signif.legend = TRUE
    )
    cat(
      paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
        `4` = 18, `5` = 31, `6` = 43, `7` = 57
      )), collapse = ""),
      "\n"
    )
  } else {
    if (x$nClasses == 3) {
      cat(centerText("Deterministic Component of SFA for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameter in variance of u (one-sided error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Parameters in variance of v (two-sided error) for latent class 3",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar), , drop = FALSE],
      P.values = TRUE,
      digits = digits, signif.legend = FALSE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      cat(centerText("Estimated prior probabilities for class membership",
        width = lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )
      ), "\n")
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
      printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
        3 * x$nvZVvar + 2 * x$nZHvar), , drop = FALSE],
      P.values = TRUE, digits = digits, signif.legend = TRUE
      )
      cat(
        paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
          `4` = 18, `5` = 31, `6` = 43, `7` = 57
        )), collapse = ""),
        "\n"
      )
    } else {
      if (x$nClasses == 4) {
        cat(
          centerText("Deterministic Component of SFA for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameter in variance of u (one-sided error) for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (two-sided error) for latent class 3",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Deterministic Component of SFA for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameter in variance of u (one-sided error) for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Parameters in variance of v (two-sided error) for latent class 4",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar), , drop = FALSE],
        P.values = TRUE,
        digits = digits, signif.legend = FALSE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        cat(
          centerText("Estimated prior probabilities for class membership",
            width = lengthSum + 2 + switch(dimCoefTable,
              `4` = 18, `5` = 31, `6` = 43, `7` = 57
            )
          ),
          "\n"
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
        printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
          4 * x$nvZVvar + 3 * x$nZHvar), , drop = FALSE],
        P.values = TRUE, digits = digits, signif.legend = TRUE
        )
        cat(
          paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )), collapse = ""),
          "\n"
        )
      } else {
        if (x$nClasses == 5) {
          cat(
            centerText("Deterministic Component of SFA for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(2 * x$nXvar + 2 * x$nuZUvar +
            2 * x$nvZVvar + 1):(3 * x$nXvar + 2 * x$nuZUvar +
            2 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one-sided error) for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 2 * x$nuZUvar +
            2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
            2 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (two-sided error) for latent class 3",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
            2 * x$nvZVvar + 1):(3 * x$nXvar + 3 * x$nuZUvar +
            3 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Deterministic Component of SFA for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(3 * x$nXvar + 3 * x$nuZUvar +
            3 * x$nvZVvar + 1):(4 * x$nXvar + 3 * x$nuZUvar +
            3 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one-sided error) for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 3 * x$nuZUvar +
            3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
            3 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (two-sided error) for latent class 4",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
            3 * x$nvZVvar + 1):(4 * x$nXvar + 4 * x$nuZUvar +
            4 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Deterministic Component of SFA for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(4 * x$nXvar + 4 * x$nuZUvar +
            4 * x$nvZVvar + 1):(5 * x$nXvar + 4 * x$nuZUvar +
            4 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameter in variance of u (one-sided error) for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 4 * x$nuZUvar +
            4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
            4 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Parameters in variance of v (two-sided error) for latent class 5",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
            4 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar), , drop = FALSE],
          P.values = TRUE,
          digits = digits, signif.legend = FALSE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          cat(
            centerText("Estimated prior probabilities for class membership",
              width = lengthSum + 2 + switch(dimCoefTable,
                `4` = 18, `5` = 31, `6` = 43, `7` = 57
              )
            ),
            "\n"
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
          printCoefmat(mlRes[(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 1):(5 * x$nXvar + 5 * x$nuZUvar +
            5 * x$nvZVvar + 4 * x$nZHvar), , drop = FALSE],
          P.values = TRUE, digits = digits, signif.legend = TRUE
          )
          cat(paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
            `4` = 18, `5` = 31, `6` = 43, `7` = 57
          )),
          collapse = ""
          ), "\n")
        }
      }
    }
  }
  cat(x$mlDate, "\n")
  cat("Log likelihood status:", x$optStatus, "\n")
  cat(
    paste0(rep("-", lengthSum + 2 + switch(dimCoefTable,
      `4` = 18, `5` = 31, `6` = 43, `7` = 57
    )), collapse = ""),
    "\n"
  )
  invisible(x)
}
