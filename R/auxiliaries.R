# Intercept check in formulas ----------

interCheckMain <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 0)
    stop("'formula' has no left hand side", call. = FALSE)
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required", call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("main model is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# intercept check in heteroscedasticity ----------

lhsCheck_u <- function(formula, scaling) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required for heteroscedasticity in u",
      call. = FALSE)
  } else {
    if (scaling == FALSE && attr(terM, "intercept") == 0)
      warning("heteroscedasticity in u is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

lhsCheck_mu <- function(formula, scaling) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required for heterogeneity in mu",
      call. = FALSE)
  } else {
    if (scaling == FALSE && attr(terM, "intercept") == 0)
      warning("heterogeneity in mu is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

lhsCheck_v <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required for heteroscedasticity in v",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("heteroscedasticity in v is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

lhsCheck_t <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required in the logit form in LCM ",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("logit form in LCM is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Formulas depending on the distribution ----------

formDist_sfacross <- function(udist, formula, muhet, uhet, vhet) {
  if (udist %in% c("tnormal", "lognormal")) {
    formula <- as.Formula(formula, muhet, uhet, vhet)
  } else {
    formula <- as.Formula(formula, uhet, vhet)
  }
  return(formula)
}

formDist_lcmcross <- function(formula, uhet, vhet, thet) {
  formula <- as.Formula(formula, uhet, vhet, thet)
  return(formula)
}

# Check infinite values ----------

is.infinite.data.frame <- function(x) {
  y <- do.call(cbind, lapply(x, is.infinite))
  if (.row_names_info(x) > 0L)
    rownames(y) <- row.names(x)
  y
}

# names for variables ----------

fName_mu_sfacross <- function(Xvar, udist, muHvar, uHvar, vHvar,
  scaling) {
  c(colnames(Xvar), if (udist == "tnormal") {
    if (scaling) {
      if (colnames(muHvar)[1] == "(Intercept)" || colnames(uHvar)[1] ==
        "(Intercept)") {
        c(paste0("Zscale_", colnames(muHvar)[-1]), "tau",
          "cu", paste0("Zv_", colnames(vHvar)))
      } else {
        c(paste0("Zscale_", colnames(muHvar)), "tau",
          "cu", paste0("Zv_", colnames(vHvar)))
      }
    } else {
      c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
        colnames(uHvar)), paste0("Zv_", colnames(vHvar)))
    }
  } else {
    if (udist == "lognormal") {
      c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
        colnames(uHvar)), paste0("Zv_", colnames(vHvar)))
    }
  })
}

fName_uv_sfacross <- function(Xvar, udist, uHvar, vHvar) {
  c(colnames(Xvar), if (udist == "gamma") {
    c(paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
      "P")
  } else {
    if (udist == "weibull") {
      c(paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
        "k")
    } else {
      if (udist == "tslaplace") {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "lambda")
      } else {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)))
      }
    }
  })
}

fName_lcmcross <- function(Xvar, uHvar, vHvar, Zvar, nZHvar, lcmClasses) {
  c(rep(c(colnames(Xvar), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar))), lcmClasses), paste0(rep(paste0("Cl",
    1:(lcmClasses - 1)), each = nZHvar), "_", rep(colnames(Zvar),
    lcmClasses - 1)))
}

# Halton sequence (code from mlogit) ----------

halton <- function(prime, length, drop) {
  halt <- 0
  t <- 0
  while (length(halt) < length + drop) {
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime -
      1, 1)/prime^t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

# matrix of draws ----------

drawMat <- function(N, Nsim, simType, prime, burn, seed, antithetics) {
  if (simType == "halton") {
    matDraw <- matrix(halton(prime = prime, length = (Nsim *
      N), drop = burn), nrow = N, ncol = Nsim, byrow = TRUE)
  } else {
    if (simType == "ghalton") {
      nthPrime <- generate_primes(min = 0, max = prime)
      idPrime <- which(prime == nthPrime)
      set.seed(seed)
      matDraw <- matrix(ghalton(n = Nsim * N, d = idPrime,
        method = "generalized"))
    } else {
      if (simType == "sobol") {
        matDraw <- matrix(sobol(n = Nsim * N, dim = 1,
          scrambling = 3, seed = seed), nrow = N, ncol = Nsim,
          byrow = TRUE)
      } else {
        if (simType == "uniform") {
          set.seed(seed)
          if (antithetics) {
          u1 <- matrix(runif(n = (Nsim * N)/2), nrow = N,
            ncol = Nsim/2, byrow = TRUE)
          u2 <- 1 - u1
          matDraw <- cbind(u1, u2)
          } else {
          matDraw <- matrix(runif(n = Nsim * N), nrow = N,
            ncol = Nsim, byrow = TRUE)
          }
        }
      }
    }
  }
}

# Inverse hessian ----------

invHess_fun <- function(hess) {
  if (!is.null(hess)) {
    hessev <- tryCatch(abs(eigen(hess, symmetric = TRUE,
      only.values = TRUE)$values), error = function(e) e)
    if (inherits(hessev, "error")) {
      warning("cannot invert hessian, using eigenvalues",
        call. = FALSE, immediate. = TRUE)
      nParm <- dim(hess)[1]
      invhess <- matrix(Inf, nParm, nParm)
    } else {
      if (min(hessev) > (1e-12 * max(hessev))) {
        ## is 1e-12 relatively acceptable!!! to what values solve does
        ## not work
        invhess <- solve(-hess)
        invhess <- (invhess + t(invhess))/2
      } else {
        warning("hessian is singular for 'qr.solve' switching to 'ginv'",
          call. = FALSE, immediate. = TRUE)
        invhess <- ginv(-as.matrix(hess))
        invhess <- (invhess + t(invhess))/2
      }
    }
    invhess
  } else return(NULL)
}

vcovObj <- function(mleObj, hessianType, method, nParm) {
  if (hessianType == 1) {
    if (method == "mla") {
      invhess <- matrix(nrow = nParm, ncol = nParm)
      invhess[upper.tri(invhess, diag = TRUE)] <- mleObj$v
      invhess[lower.tri(invhess)] <- t(invhess)[lower.tri(invhess)]
      invhess <- (invhess + t(invhess))/2
    } else {
      if (method == "sparse") {
        invhess <- invHess_fun(hess = -mleObj$hessian)
      } else {
        invhess <- invHess_fun(hess = mleObj$hessian)
      }
    }
  } else {
    if (hessianType == 2) {
      if (method %in% c("bfgs", "bhhh", "nr", "nm")) {
        invhess <- invHess_fun(hess = mleObj$hessian)
      } else {
        hess <- -crossprod(mleObj$gradL_OBS)
        invhess <- invHess_fun(hess = hess)
      }
    } else {
      if (hessianType == 3) {
        if (method == "mla") {
          invhess <- matrix(nrow = nParm, ncol = nParm)
          invhess[upper.tri(invhess, diag = TRUE)] <- mleObj$v
          invhess[lower.tri(invhess)] <- t(invhess)[lower.tri(invhess)]
          invhess <- (invhess + t(invhess))/2
          invhess <- invhess %*% crossprod(mleObj$gradL_OBS) %*%
          invhess  # H^(-1)GH^(-1) G: outer product of gradient
        } else {
          if (method == "sparse") {
          invhess <- invHess_fun(hess = -mleObj$hessian)
          invhess <- invhess %*% crossprod(mleObj$gradL_OBS) %*%
            invhess
          } else {
          invhess <- invHess_fun(hess = mleObj$hessian)
          invhess <- invhess %*% crossprod(mleObj$gradL_OBS) %*%
            invhess
          }
        }
      }
    }
  }
  invhess
}

# SFA + distribution ----------

sfadist <- function(udist) {
  switch(udist, tnormal = "Truncated-Normal Normal Stochastic Frontier Model",
    hnormal = "Normal-Half Normal Stochastic Frontier Model",
    exponential = "Exponential Normal Stochastic Frontier Model",
    rayleigh = "Rayleigh Normal Stochastic Frontier Model",
    uniform = "Uniform Normal Stochastic Frontier Model",
    gamma = "Gamma Normal Stochastic Frontier Model",
    lognormal = "Log-Normal Normal Stochastic Frontier Model",
    weibull = "Weibull Normal Stochastic Frontier Model",
    genexponential = "Generalized-Exponential Normal Stochastic Frontier Model",
    tslaplace = "Truncated Skewed-Laplace Normal Stochastic Frontier Model")
}

# variance of u ----------

varuFun <- function(object, mu, P, lambda, k) {
  if (object$udist == "hnormal") {
    object$sigmauSq * (pi - 2)/pi
  } else {
    if (object$udist == "tnormal") {
      a <- (pnorm(-mean(mu)/sqrt(object$sigmauSq)))^(-1)
      mean(mu)^2 * a/2 * (1 - a/2) + a/2 * (pi - a)/pi *
        object$sigmauSq
    } else {
      if (object$udist == "exponential") {
        object$sigmauSq
      } else {
        if (object$udist == "rayleigh") {
          (4 - pi)/2 * object$sigmauSq
        } else {
          if (object$udist == "gamma") {
          P * object$sigmauSq
          } else {
          if (object$udist == "lognormal") {
            (exp(object$sigmauSq) - 1) * exp(2 * mean(mu) +
            object$sigmauSq)
          } else {
            if (object$udist == "uniform") {
            object$sigmauSq
            } else {
            if (object$udist == "genexponential") {
              5/4 * object$sigmauSq
            } else {
              if (object$udist == "tslaplace") {
              object$sigmauSq * (1 + 8 * lambda +
                16 * lambda^2 + 12 * lambda^3 +
                4 * lambda^4)/((1 + lambda)^2 *
                (1 + 2 * lambda)^2)
              } else {
              if (object$udist == "weibull") {
                object$sigmauSq * (gamma(1 + 2/k) -
                (gamma(1 + 1/k))^2)
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

# expected value of u ----------

euFun <- function(object, mu, P, lambda, k) {
  if (object$udist == "hnormal") {
    sqrt(object$sigmauSq) * sqrt(2/pi)
  } else {
    if (object$udist == "tnormal") {
      -exp(-mean(mu)^2/(2 * object$sigmauSq)) * (sqrt(pi) *
        (sqrt(2) * mean(mu) * erfc(mean(mu)/(sqrt(2 *
          object$sigmauSq))) - 2^(3/2) * mean(mu)) *
        exp(mean(mu)^2/(2 * object$sigmauSq)) - 2 * sqrt(object$sigmauSq))/(sqrt(pi) *
        (sqrt(2) * erf(mean(mu)/sqrt(2 * object$sigmauSq)) +
          sqrt(2)))
    } else {
      if (object$udist == "exponential") {
        sqrt(object$sigmauSq)
      } else {
        if (object$udist == "rayleigh") {
          sqrt(object$sigmauSq) * sqrt(pi/2)
        } else {
          if (object$udist == "gamma") {
          P * sqrt(object$sigmauSq)
          } else {
          if (object$udist == "lognormal") {
            exp(mean(mu) + object$sigmauSq/2)
          } else {
            if (object$udist == "uniform") {
            sqrt(12 * object$sigmauSq)/2
            } else {
            if (object$udist == "genexponential") {
              3/2 * sqrt(object$sigmauSq)
            } else {
              if (object$udist == "tslaplace") {
              sqrt(object$sigmauSq) * (1 + 4 *
                lambda + 2 * lambda^2)/((1 + lambda) *
                (1 + 2 * lambda))
              } else {
              if (object$udist == "weibull") {
                sqrt(object$sigmauSq) * gamma(1 +
                1/k)
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

# expected value of E(exp(-u)) ----------

eExpuFun <- function(object, mu, P, lambda, k) {
  if (object$udist == "hnormal") {
    2 * (1 - pnorm(sqrt(object$sigmauSq))) * exp(object$sigmauSq/2)
  } else {
    if (object$udist == "tnormal") {
      -sqrt(pi) * (exp(object$sigmauSq/2) * erf((sqrt(2) *
        object$sigmauSq - sqrt(2) * mean(mu))/(2 * sqrt(object$sigmauSq))) -
        exp(object$sigmauSq/2))/(sqrt(pi) * (exp(mean(mu)) *
        erf(mean(mu)/sqrt(2 * object$sigmauSq)) + exp(mean(mu))))
    } else {
      if (object$udist == "exponential") {
        1/(1 + sqrt(object$sigmauSq))
      } else {
        if (object$udist == "gamma") {
          ((1/sqrt(object$sigmauSq))/(1 + 1/sqrt(object$sigmauSq)))^P
        } else {
          if (object$udist == "lognormal") {
          integrate(fnExpULogNorm, lower = 0, upper = Inf,
            rel.tol = 1e-10, stop.on.error = FALSE,
            sigma = sqrt(object$sigmauSq), mu = mean(mu))$value
          } else {
          if (object$udist == "uniform") {
            (1 - exp(-object$theta))/object$theta
          } else {
            if (object$udist == "rayleigh") {
            1 - (exp(object$sigmauSq/2) * sqrt(2 *
              pi) * sqrt(object$sigmauSq) * erfc(sqrt(object$sigmauSq)/sqrt(2)))/2
            } else {
            if (object$udist == "genexponential") {
              2/((sqrt(object$sigmauSq) + 1) * (sqrt(object$sigmauSq) +
              2))
            } else {
              if (object$udist == "tslaplace") {
              (1 + lambda)/(2 * lambda + 1) * (2/(sqrt(object$sigmauSq) +
                1) - 1/(1 + sqrt(object$sigmauSq) +
                lambda))
              } else {
              if (object$udist == "weibull") {
                integrate(fnExpUWeiNorm, lower = 0,
                upper = Inf, rel.tol = 1e-10,
                stop.on.error = FALSE, sigma = sqrt(object$sigmauSq),
                k = k)$value
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

# Center text strings (adapted from gdata package) ----------

trimChar <- function(s, recode.factor = TRUE, ...) {
  s <- sub(pattern = "^[[:blank:]]+", replacement = "", x = s)
  s <- sub(pattern = "[[:blank:]]+$", replacement = "", x = s)
  s
}

centerText <- function(x, width) {
  retval <- vector(length = length(x), mode = "character")
  for (i in 1:length(x)) {
    text <- trimChar(x[i])
    textWidth <- nchar(text)
    nspaces <- floor((width - textWidth)/2)
    spaces <- paste(rep(" ", nspaces), sep = "", collapse = "")
    retval[i] <- paste(spaces, text, sep = "", collapse = "\n")
  }
  retval
}


# S3methods

efficiencies <- function(object, ...) {
  UseMethod("efficiencies")
}

marginal <- function(object, ...) {
  UseMethod("marginal")
}

ic <- function(object, ...) {
  UseMethod("ic")
}


setClass("dagoTest",
    representation(
        call = "call",
        data = "list",
        test = "list",
        title = "character")
)

setMethod("show", "dagoTest",
      function(object)
{  # Unlike print the argument for show is 'object'.
    x = object

    # Title:
    cat("## ", "D'Agostino's  Test", " ##\n", sep = "")

    # Test Results:
    test = x@test
    cat("\nTest Results:\n", sep = "")

    # Statistic:
    if (!is.null(test$statistic)) {
        statistic = test$statistic
        Names = names(statistic)
        cat("  STATISTIC:\n")
        for (i in 1:length(Names)) {
            if (!is.na(statistic[i])) {
                cat(paste("    ", Names[i], ": ",
                    round(statistic[i], digits = 4), "\n", sep = "" ) )
            }
        }
    }

    # P-Value:
    if (!is.null(test$p.value)) {
        pval = test$p.value
        Names = names(pval)
        if (Names[1] == "") space = "" else space = ": "
        cat("  P.VALUE:\n")
        for (i in 1:length(Names)) {
            if (!is.na(pval[i])) {
                if (class(version) != "Sversion") {
                    cat(paste("    ", Names[i], space,
                    format.pval(pval[i], digits = 4), " \n", sep = "" ) )
                } else {
                    cat(paste("    ", Names[i], space,
                    round(pval[i], digits = 4), " \n", sep = "" ) )
                }
            }
        }
    }
})
