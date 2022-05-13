################################################################################
#                                                                              #
# R auxiliary functions for the sfaR package                                   #
#                                                                              #
################################################################################

# Intercept check for main formula ----------
#' @param formula main formula of model
#' @noRd
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

# Intercept check for selection formula ----------
#' @param formula main formula of model
#' @noRd
interCheckSelection <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 0)
    stop("'formula' has no left hand side", call. = FALSE)
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required", call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("Selection model is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Intercept hetero. in u (cross-section) ----------
#' @param formula formula for heteroscedasticity in inefficiency term
#' @param scaling logical for scaling property
#' @noRd
clhsCheck_u <- function(formula, scaling) {
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

# Intercept hetero. in mu (cross section) ----------
#' @param formula formula for heterogeneity in inefficiency term
#' @param scaling logical for scaling property
#' @noRd
clhsCheck_mu <- function(formula, scaling) {
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

# Intercept hetero. in v (cross section) ----------
#' @param formula formula for heteroscedasticity in noise component
#' @noRd
clhsCheck_v <- function(formula) {
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

# Intercept separating variables in LCM ----------
#' @param formula formula for logit form in LCM
#' @noRd
clhsCheck_t <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required in the logit form in LCM",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("logit form in LCM is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Intercept separating variables in ZISF ----------
#' @param formula formula for logit form in ZISF
#' @noRd
clhsCheck_q <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required in the logit form in ZISF",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("logit form in ZISF is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Intercept hetero. in u (panel data pl81) ----------
#' @param formula main formula of model
#' @noRd
plhsCheck_u <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required for heteroscedasticity in u",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("heteroscedasticity in u is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Intercept hetero. in mu (panel data pl81) ----------
#' @param formula formula for heterogeneity in inefficiency term
#' @noRd
plhsCheck_mu <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required for heterogeneity in mu",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("heterogeneity in mu is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# remove int. in uhet (panel data bc92c) ----------
#' @param formula formula for group variable
#' @noRd
plhsCheck_u_bc92c <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (attr(terM, "intercept") == 1) {
    formula <- formula(paste0(" ~ 0 + ", attr(terM, "term.labels")))
  }
  return(formula)
}

# Int. sep. var. in ZISF part of NLC ----------
#' @param formula formula for logit form in ZISF part of NLC
#' @noRd
clhsCheck_sf_nlc <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0 && attr(terM,
    "intercept") == 0) {
    stop("at least one exogenous variable is required in the logit form in ZISF part of NLC",
      call. = FALSE)
  } else {
    if (attr(terM, "intercept") == 0)
      warning("logit form in ZISF part of NLC is estimated without intercept",
        call. = FALSE)
  }
  return(formula)
}

# Remove intercept in ghet (metacross) ----------
#' @param formula formula for group variable
#' @noRd
clhsCheck_meta <- function(formula) {
  terM <- terms(formula(formula))
  if (attr(terM, "response") == 1)
    formula[[2]] <- NULL
  if (length(attr(terM, "term.labels")) == 0) {
    stop("a group variable must be specified for metafrontier estimation",
      call. = FALSE)
  }
  if (length(attr(terM, "term.labels")) > 1) {
    stop("only one group variable must be provided for metafrontier estimation",
      call. = FALSE)
  }
  if (attr(terM, "intercept") == 1) {
    formula <- formula(paste0(" ~ 0 + ", attr(terM, "term.labels")))
  }
  return(formula)
}

# Formulas depending on the distribution ----------
#' @param udist inefficiency term distribution
#' @param formula formula for model
#' @param muhet heterogeneity in mu
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @noRd
formDist_sfacross <- function(udist, formula, muhet, uhet, vhet) {
  if (udist %in% c("tnormal", "lognormal")) {
    formula <- as.Formula(formula, muhet, uhet, vhet)
  } else {
    formula <- as.Formula(formula, uhet, vhet)
  }
  return(formula)
}

#' @param formula formula for model
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @param thet separating variables for LCM
#' @noRd
formDist_lcmcross <- function(formula, uhet, vhet, thet) {
  formula <- as.Formula(formula, uhet, vhet, thet)
  return(formula)
}

#' @param udist inefficiency term distribution
#' @param formula formula for model
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @noRd
formDist_selectioncross <- function(udist, formula, uhet, vhet) {
  formula <- as.Formula(formula, uhet, vhet)
  return(formula)
}

#' @param udist inefficiency term distribution
#' @param formula formula for model
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @param qhet separating variables for ZISF
#' @noRd
formDist_zisfcross <- function(udist, formula, muhet, uhet, vhet,
  qhet) {
  if (udist %in% c("tnormal", "lognormal")) {
    formula <- as.Formula(formula, muhet, uhet, vhet, qhet)
  } else {
    formula <- as.Formula(formula, uhet, vhet, qhet)
  }
  return(formula)
}

#' @param formula formula for model
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @param thet separating variables for LCM part in NLC
#' @param sfhet separating variables for ZISF part in NLC
#' 
#' @noRd
formDist_nlccross <- function(formula, uhet, vhet, thet, sfhet) {
  formula <- as.Formula(formula, uhet, vhet, thet, sfhet)
  return(formula)
}

#' @param udist inefficiency term distribution
#' @param formula formula for model
#' @param uhet heteroscedasticity in u
#' @param vhet heteroscedasticity in v
#' @param ghet group variable for metafrontier
#' @noRd
formDist_metacross <- function(udist, formula, muhet, uhet, vhet,
  ghet) {
  if (udist %in% c("tnormal", "lognormal")) {
    formula <- as.Formula(formula, muhet, uhet, vhet, ghet)
  } else {
    formula <- as.Formula(formula, uhet, vhet, ghet)
  }
  return(formula)
}

# Check infinite values ----------
#' @param x data frame
#' @noRd
is.infinite.data.frame <- function(x) {
  y <- do.call(cbind, lapply(x, is.infinite))
  if (.row_names_info(x) > 0L)
    rownames(y) <- row.names(x)
  y
}

# names for variables ----------
#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param udist inefficiency distribution
#' @param muHvar heterogeneity variables in mu
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @param scaling logical for scaling property
#' @noRd
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

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param udist inefficiency distribution
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @noRd
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

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param udist inefficiency distribution
#' @param muHvar heterogeneity variables in mu
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @noRd
fName_mu_sfapanel <- function(Xvar, udist, muHvar, uHvar, vHvar,
  modelType) {
  if (modelType == "pl81") {
    c(colnames(Xvar), if (udist == "gamma") {
      c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
        colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
        "P")
    } else {
      if (udist == "weibull") {
        c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
          "k")
      } else {
        if (udist == "tslaplace") {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
          "lambda")
        } else {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)))
        }
      }
    })
  } else {
    if (modelType == "bc92a") {
      c(colnames(Xvar), if (udist == "gamma") {
        c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
          "P", "eta")
      } else {
        if (udist == "weibull") {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
          "k", "eta")
        } else {
          if (udist == "tslaplace") {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
            colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
            "lambda", "eta")
          } else {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
            colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
            "eta")
          }
        }
      })
    } else {
      if (modelType %in% c("bc92b", "k90")) {
        c(colnames(Xvar), if (udist == "gamma") {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
          colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
          "P", "eta1", "eta2")
        } else {
          if (udist == "weibull") {
          c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
            colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
            "k", "eta1", "eta2")
          } else {
          if (udist == "tslaplace") {
            c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
            colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
            "lambda", "eta1", "eta2")
          } else {
            c(paste0("Zmu_", colnames(muHvar)), paste0("Zu_",
            colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
            "eta1", "eta2")
          }
          }
        })
      }
    }
  }
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param udist inefficiency distribution
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @noRd
fName_uv_sfapanel <- function(Xvar, udist, uHvar, vHvar, modelType) {
  if (modelType == "pl81") {
    c(colnames(Xvar), if (udist == "gamma") {
      c(paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
        "P")
    } else {
      if (udist == "weibull") {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "k")
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
  } else {
    if (modelType == "bc92a") {
      c(colnames(Xvar), if (udist == "gamma") {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "P", "eta")
      } else {
        if (udist == "weibull") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "k", "eta")
        } else {
          if (udist == "tslaplace") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), "lambda", "eta")
          } else {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), "eta")
          }
        }
      })
    } else {
      if (modelType %in% c("bc92b", "k90")) {
        c(colnames(Xvar), if (udist == "gamma") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "P", "eta1", "eta2")
        } else {
          if (udist == "weibull") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), "k", "eta1", "eta2")
          } else {
          if (udist == "tslaplace") {
            c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), "lambda", "eta1", "eta2")
          } else {
            c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), "eta1", "eta2")
          }
          }
        })
      }
    }
  }
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @param Zvar separating variables in LCM
#' @param nZHvar number of separating variables in LCM
#' @param lcmClasses number of classes in LCM
#' @noRd
fName_lcmcross <- function(Xvar, uHvar, vHvar, Zvar, nZHvar,
  lcmClasses) {
  c(rep(c(colnames(Xvar), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar))), lcmClasses), paste0(rep(paste0("Cl",
    1:(lcmClasses - 1)), each = nZHvar), "_", rep(colnames(Zvar),
    lcmClasses - 1)))
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @noRd
fName_uvr_selectioncross <- function(Xvar, uHvar, vHvar) {
  c(colnames(Xvar), c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), "rho"))
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param muHvar heterogeneity variables in mu
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @param Zvar separating variables in ZISF
#' @noRd
fName_mu_zisfcross <- function(Xvar, muHvar, uHvar, vHvar, Zvar,
  sigmavType) {
  if (sigmavType == "common") {
    c(colnames(Xvar), c(paste0("Zmu_", colnames(muHvar)),
      paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
      paste0("SF_", colnames(Zvar))))
  } else {
    if (sigmavType == "different") {
      c(colnames(Xvar), c(paste0("Zmu_", colnames(muHvar)),
        paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
        paste0("ZI_", colnames(Zvar))))
    }
  }
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param udist inefficiency distribution
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @param Zvar separating variables in ZISF
#' @noRd
fName_uv_zisfcross <- function(Xvar, udist, uHvar, vHvar, Zvar,
  sigmavType) {
  if (sigmavType == "common") {
    c(colnames(Xvar), if (udist == "gamma") {
      c(paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
        "P", paste0("ZI_", colnames(Zvar)))
    } else {
      if (udist == "weibull") {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "k", paste0("ZI_", colnames(Zvar)))
      } else {
        if (udist == "tslaplace") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), "lambda", paste0("Cl_",
          colnames(Zvar)))
        } else {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), paste0("ZI_", colnames(Zvar)))
        }
      }
    })
  } else {
    if (sigmavType == "different") {
      c(colnames(Xvar), if (udist == "gamma") {
        c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
          "P", paste0("ZI_", colnames(Zvar)))
      } else {
        if (udist == "weibull") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
          colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
          "k", paste0("ZI_", colnames(Zvar)))
        } else {
          if (udist == "tslaplace") {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
            "lambda", paste0("ZI_", colnames(Zvar)))
          } else {
          c(paste0("Zu_", colnames(uHvar)), paste0("Zv_",
            colnames(vHvar)), paste0("Zv_", colnames(vHvar)),
            paste0("ZI_", colnames(Zvar)))
          }
        }
      })
    }
  }
}

#' @param Xvar variables in main formula (e.g. inputs/outputs)
#' @param uHvar heteroscedasticity variables in u
#' @param vHvar heteroscedasticity variables in v
#' @param Zvar matrix of separating variables for latent classes
#' @param nZHvar number of separating variables for latent classes
#' @param Qvar matrix of separating variables for inefficient class
#' @param nQHvar number of separating variables for inefficient class
#' @param lcmClasses number of classes in LCM
#' @noRd
fName_nlccross <- function(Xvar, uHvar, vHvar, Zvar, nZHvar,
  Qvar, nQHvar, lcmClasses) {
  c(colnames(Xvar), paste0("Zu_", colnames(uHvar)), paste0("Zv_",
    colnames(vHvar)), paste0("SF1_", colnames(Qvar)), colnames(Xvar),
    paste0("Zu_", colnames(uHvar)), paste0("Zv_", colnames(vHvar)),
    paste0("SF2_", colnames(Qvar)), paste0(rep(paste0("Cl",
      1:(lcmClasses - 1)), each = nZHvar), "_", rep(colnames(Zvar),
      lcmClasses - 1)))
}

# Compute skewness ----------
#' @param x vector for which skewness is computed
#' @noRd
skewness <- function(x) {
  n <- length(x)
  (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
}

# Halton sequence (code from mlogit) ----------
#' @param prime prime number
#' @param length length of halton sequence
#' @param drop number of first observations to drop
#' @noRd
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

list_primes <- randtoolbox::get.primes(1e+05)

#' @param prime prime number
#' @noRd
is_prime <- function(prime) {
  prime %in% list_primes
}

#' @param N number of observations
#' @param Nsim number of draw per observation
#' @param prime prime number
#' @param burn number of first observations dropped
#' @param seed seed for the draw
#' @param antithetics logical for antithetics draws
#' @noRd
drawMat <- function(N, Nsim, simType, prime, burn, seed, antithetics) {
  if (simType == "halton") {
    matDraw <- matrix(halton(prime = prime, length = (Nsim *
      N), drop = burn), nrow = N, ncol = Nsim, byrow = TRUE)
  } else {
    if (simType == "ghalton") {
      idPrime <- which(list_primes == prime)
      set.seed(seed)
      matDraw <- matrix(ghalton(n = Nsim * N, d = idPrime,
        method = "generalized"), nrow = N, ncol = Nsim,
        byrow = TRUE)
    } else {
      if (simType == "sobol") {
        matDraw <- matrix(sobol(n = Nsim * N, dim = 1,
          scrambling = 1, seed = seed), nrow = N, ncol = Nsim,
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
  return(matDraw)
}

# Inverse hessian ----------
#' @param X matrix X
#' @noRd
ginvsfaR <- function(X, tol = sqrt(.Machine$double.eps)) {
  if (!is.numeric(X))
    stop("'hessian' must be a numeric matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u)) else if (!any(Positive))
    array(0, dim(X)[2L:1L]) else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
    t(Xsvd$u[, Positive, drop = FALSE]))
}

#' @param hess hessian matrix
#' @noRd
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
        ## is 1e-12 relatively acceptable!!! to what values
        ## solve does not work
        invhess <- solve(-hess)
        invhess <- (invhess + t(invhess))/2
      } else {
        warning("hessian is singular for 'qr.solve' switching to 'ginv'",
          call. = FALSE, immediate. = TRUE)
        invhess <- ginvsfaR(-as.matrix(hess))
        invhess <- (invhess + t(invhess))/2
      }
    }
    invhess
  } else return(NULL)
}

#' @param mleObj optimization object
#' @param hessianType numeric for type of inverse hessian
#' @param method optimization algorithm
#' @param nParm number of parameters in optimization
#' @noRd
vcovObj <- function(mleObj, hessianType, method, nParm) {
  if (hessianType == 1) {
    if (method == "mla") {
      invhess <- matrix(nrow = nParm, ncol = nParm)
      invhess[upper.tri(invhess, diag = TRUE)] <- mleObj$v
      invhess[lower.tri(invhess)] <- t(invhess)[lower.tri(invhess)]
      invhess <- (invhess + t(invhess))/2
    } else {
      if (method %in% c("sparse", "trust")) {
        invhess <- invHess_fun(hess = -mleObj$hessian)
      } else {
        invhess <- invHess_fun(hess = mleObj$hessian)
      }
    }
  } else {
    if (hessianType == 2) {
      if (method %in% c("bfgs", "bhhh", "nr", "cg", "nm")) {
        invhess <- invHess_fun(hess = mleObj$hessian)
      } else {
        hess <- -crossprod(mleObj$gradL_OBS)
        invhess <- invHess_fun(hess = hess)
      }
    }
  }
  invhess
}

# SFA + distribution ----------
#' @param udist inefficiency distribution
#' @noRd
sfadist <- function(udist) {
  switch(udist, tnormal = "Truncated-Normal Normal SF Model",
    hnormal = "Normal-Half Normal SF Model", exponential = "Exponential Normal SF Model",
    rayleigh = "Rayleigh Normal SF Model", uniform = "Uniform Normal SF Model",
    gamma = "Gamma Normal SF Model", lognormal = "Log-Normal Normal SF Model",
    weibull = "Weibull Normal SF Model", genexponential = "Generalized-Exponential Normal SF Model",
    tslaplace = "Truncated Skewed-Laplace Normal SF Model")
}

# ZISF + distribution ----------
#' @param udist inefficiency distribution
#' @noRd
zisfdist <- function(udist) {
  switch(udist, tnormal = "Truncated-Normal Normal ZISF model",
    hnormal = "Normal-Half Normal ZISF Model", exponential = "Exponential Normal ZISF model",
    rayleigh = "Rayleigh Normal ZISF model", uniform = "Uniform Normal ZISF model",
    gamma = "Gamma Normal ZISF model", lognormal = "Log-Normal Normal ZISF model",
    weibull = "Weibull Normal ZISF model", genexponential = "Generalized-Exponential Normal ZISF model",
    tslaplace = "Truncated Skewed-Laplace Normal ZISF model")
}

# variance of u ----------
#' @param object object from sfacross/lcmcross ...
#' @param mu mu in truncated normal and log-normal distribution
#' @param P Shape parameter in gamma distribution
#' @param lambda parameter in truncated skewed laplace distribution
#' @param k parameter in weibull distribution
#' @noRd
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
#' @param x vector for error function
#' @param object object from sfacross/lcmcross ...
#' @param mu mu in truncated normal and log-normal distribution
#' @param P Shape parameter in gamma distribution
#' @param lambda parameter in truncated skewed laplace distribution
#' @param k parameter in weibull distribution
#' @noRd
erf <- function(x) {
  # 2*pnorm(sqrt(2)*x)-1 # or
  pchisq(2 * x^2, 1) * sign(x)
}

erfc <- function(x) {
  # 1 - erf(x)
  2 * pnorm(-sqrt(2) * x)
}

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
#' @param object object from sfacross/lcmcross ...
#' @param mu mu in truncated normal and log-normal distribution
#' @param P Shape parameter in gamma distribution
#' @param lambda parameter in truncated skewed laplace distribution
#' @param k parameter in weibull distribution
#' @noRd
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

# Likelihood for probit model ----------

## Likelihood ----------
probit_likelihood <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  wHvar * (Yvar * log(pnorm(Z)) + (1 - Yvar) * log(pnorm(-Z)))
}

## Gradient ----------
probit_gradient <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  gx <- wHvar * (Yvar * dnorm(Z)/pnorm(Z) - (1 - Yvar) * dnorm(-Z)/pnorm(-Z))
  sweep(Xvar, MARGIN = 1, STATS = gx, FUN = "*")
}

# Likelihood for logit model ----------

## Likelihood ----------
logit_likelihood <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  wHvar * (Yvar * log(exp(Z)/(1 + exp(Z))) + (1 - Yvar) * log(1/(1 +
    exp(Z))))
}

## Gradient ----------
logit_gradient <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  gx <- wHvar * (Yvar * (1 - exp(Z)/(1 + exp(Z))) - (1 - Yvar) *
    exp(Z)/(1 + exp(Z)))
  sweep(Xvar, MARGIN = 1, STATS = gx, FUN = "*")
}

# Likelihood for cauchit model ----------

## Likelihood ----------
cauchit_likelihood <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  wHvar * (Yvar * log(1/pi * atan(Z) + 1/2) + (1 - Yvar) *
    log(1/2 - 1/pi * atan(Z)))
}

## Gradient ----------
cauchit_gradient <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  gx <- wHvar * (Yvar/(0.5 + atan(Z)/pi) - (1 - Yvar)/(0.5 -
    atan(Z)/pi))/(pi * ((Z)^2 + 1))
  sweep(Xvar, MARGIN = 1, STATS = gx, FUN = "*")
}

# Likelihood for cloglog model ----------

## Likelihood ----------
cloglog_likelihood <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  wHvar * (Yvar * log(1 - exp(-exp(Z))) + (1 - Yvar) * log(exp(-exp(Z))))
}

## Gradient ----------
cloglog_gradient <- function(beta, Xvar, Yvar, wHvar) {
  Z <- as.numeric(crossprod(matrix(beta), t(Xvar)))
  gx <- wHvar * exp(Z) * (Yvar * (1 + exp(-exp(Z))/(1 - exp(-exp(Z)))) -
    1)
  sweep(Xvar, MARGIN = 1, STATS = gx, FUN = "*")
}

# Center text strings (gdata package) ----------
#' @param s string
#' @noRd
trimChar <- function(s, recode.factor = TRUE, ...) {
  s <- sub(pattern = "^[[:blank:]]+", replacement = "", x = s)
  s <- sub(pattern = "[[:blank:]]+$", replacement = "", x = s)
  s
}

#' @param x strings
#' @param width width
#' @noRd
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

# mix-chisquare distribution (emdbook) ----------
#' @param p numeric vector of positive values
#' @param q numeric vector of quantiles (0-1)
#' @param df degrees of freedom (positive integer)
#' @param mix mixture parameter: fraction of distribution that is chi-square(n-1) distributed
#' @param lower.tail return lower tail values?
#' @param log.p return log probabilities?
#' @noRd
pchibarsq <- function(p, df = 1, mix = 0.5, lower.tail = TRUE,
  log.p = FALSE) {
  df <- rep(df, length.out = length(p))
  mix <- rep(mix, length.out = length(p))
  c1 <- ifelse(df == 1, if (lower.tail)
    1 else 0, pchisq(p, df - 1, lower.tail = lower.tail))
  c2 <- pchisq(p, df, lower.tail = lower.tail)
  r <- mix * c1 + (1 - mix) * c2
  if (log.p)
    log(r) else r
}

qchibarsq <- function(q, df = 1, mix = 0.5) {
  n <- max(length(q), length(df), length(mix))
  df <- rep(df, length.out = n)
  mix <- rep(mix, length.out = n)
  q <- rep(q, length.out = n)
  tmpf2 <- function(q, df, mix) {
    if (df > 1) {
      tmpf <- function(x) {
        pchibarsq(x, df, mix) - q
      }
      uniroot(tmpf, lower = qchisq(q, df - 1), upper = qchisq(q,
        df))$root
    } else {
      newq <- (q - mix)/(1 - mix)
      ifelse(newq < 0, 0, qchisq(newq, df = 1))
    }
  }
  mapply(tmpf2, q, df, mix)
}

# D'Agostino normality test (fBasics) ----------
#' @param x numeric vector of data values
#' @noRd
.skewness.test <- function(x) {
  n <- length(x)
  if (n < 8)
    stop("Sample size must be at least 8 for skewness test")
  meanX <- mean(x)
  s <- sqrt(mean((x - meanX)^2))
  a3 <- mean((x - meanX)^3)/s^3
  SD3 <- sqrt(6 * (n - 2)/((n + 1) * (n + 3)))
  U3 <- a3/SD3
  b <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3))/((n -
    2) * (n + 5) * (n + 7) * (n + 9))
  W2 <- sqrt(2 * (b - 1)) - 1
  delta <- 1/sqrt(log(sqrt(W2)))
  a <- sqrt(2/(W2 - 1))
  Z3 <- delta * log((U3/a) + sqrt((U3/a)^2 + 1))
  pZ3 <- 2 * (1 - pnorm(abs(Z3), 0, 1))
  names(Z3) <- "Z3"
  RVAL <- list(statistic = Z3, p.value = pZ3)
  return(RVAL)
}

.kurtosis.test <- function(x) {
  n <- length(x)
  if (n < 20)
    stop("Sample size must be at least 20 for kurtosis test")
  meanX <- mean(x)
  s <- sqrt(mean((x - meanX)^2))
  a4 <- mean((x - meanX)^4)/s^4
  SD4 <- sqrt(24 * (n - 2) * (n - 3) * n/((n + 1)^2 * (n +
    3) * (n + 5)))
  U4 <- (a4 - 3 + 6/(n + 1))/SD4
  B <- (6 * (n * n - 5 * n + 2)/((n + 7) * (n + 9))) * sqrt((6 *
    (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
  A <- 6 + (8/B) * ((2/B) + sqrt(1 + 4/(B^2)))
  jm <- sqrt(2/(9 * A))
  pos <- ((1 - 2/A)/(1 + U4 * sqrt(2/(A - 4))))^(1/3)
  Z4 <- (1 - 2/(9 * A) - pos)/jm
  pZ4 <- 2 * (1 - pnorm(abs(Z4), 0, 1))
  names(Z4) <- "Z4"
  RVAL <- list(statistic = Z4, p.value = pZ4)
  return(RVAL)
}

.omnibus.test <- function(x) {
  n <- length(x)
  if (n < 20)
    stop("sample size must be at least 20 for omnibus test")
  meanX <- mean(x)
  s <- sqrt(mean((x - meanX)^2))
  a3 <- mean((x - meanX)^3)/s^3
  a4 <- mean((x - meanX)^4)/s^4
  SD3 <- sqrt(6 * (n - 2)/((n + 1) * (n + 3)))
  SD4 <- sqrt(24 * (n - 2) * (n - 3) * n/((n + 1)^2 * (n +
    3) * (n + 5)))
  U3 <- a3/SD3
  U4 <- (a4 - 3 + 6/(n + 1))/SD4
  b <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3))/((n -
    2) * (n + 5) * (n + 7) * (n + 9))
  W2 <- sqrt(2 * (b - 1)) - 1
  delta <- 1/sqrt(log(sqrt(W2)))
  a <- sqrt(2/(W2 - 1))
  Z3 <- delta * log((U3/a) + sqrt((U3/a)^2 + 1))
  B <- (6 * (n * n - 5 * n + 2)/((n + 7) * (n + 9))) * sqrt((6 *
    (n + 3) * (n + 5))/(n * (n - 2) * (n - 3)))
  A <- 6 + (8/B) * ((2/B) + sqrt(1 + 4/(B^2)))
  jm <- sqrt(2/(9 * A))
  pos <- ((1 - 2/A)/(1 + U4 * sqrt(2/(A - 4))))^(1/3)
  Z4 <- (1 - 2/(9 * A) - pos)/jm
  omni <- Z3^2 + Z4^2
  pomni <- 1 - pchisq(omni, 2)
  RVAL <- list(statistic = omni, p.value = pomni)
  return(RVAL)
}

setClass("fHTEST", representation(call = "call", data = "list",
  test = "list", title = "character", description = "character"))

dagoTest <- function(x) {
  x <- as.vector(x)
  call <- match.call()
  test <- .omnibus.test(x)
  skew <- .skewness.test(x)
  kurt <- .kurtosis.test(x)
  test$data.name <- "ols residuals"
  PVAL <- c(test$p.value, skew$p.value, kurt$p.value)
  names(PVAL) <- c("Omnibus  Test", "Skewness Test", "Kurtosis Test")
  test$p.value <- PVAL
  STATISTIC <- c(test$statistic, skew$statistic, kurt$statistic)
  names(STATISTIC) <- c("Chi2 | Omnibus", "Z3  | Skewness",
    "Z4  | Kurtosis")
  test$statistic <- STATISTIC
  class(test) <- "list"
  title <- "D'Agostino Normality Test"
  new("fHTEST", call = call, data = list(x = x), test = test,
    title = as.character(title))
}

# S3methods ----------
#' @param object sfacross, lcmcross, selectioncross ... objects
#' @noRd
efficiencies <- function(object, ...) {
  UseMethod("efficiencies", object)
}

#' @param object sfacross, lcmcross, selectioncross ... objects
#' @noRd
ic <- function(object, ...) {
  UseMethod("ic", object)
}

#' @param object sfacross, lcmcross, selectioncross ... objects
#' @noRd
marginal <- function(object, ...) {
  UseMethod("marginal", object)
}

# #' @param object sfacross, lcmcross, selectioncross ...
# objects #' @noRd nobs <- function(x, ...) {
# UseMethod('nobs', x) }

setClass("dagoTest", representation(call = "call", data = "list",
  test = "list", title = "character"))

setMethod("show", "dagoTest", function(object) {
  # Unlike print the argument for show is 'object'.
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
        cat(paste("    ", Names[i], ": ", round(statistic[i],
          digits = 4), "\n", sep = ""))
      }
    }
  }
  # P-Value:
  if (!is.null(test$p.value)) {
    pval = test$p.value
    Names = names(pval)
    if (Names[1] == "")
      space = "" else space = ": "
    cat("  P.VALUE:\n")
    for (i in 1:length(Names)) {
      if (!is.na(pval[i])) {
        if (class(version) != "Sversion") {
          cat(paste("    ", Names[i], space, format.pval(pval[i],
          digits = 4), " \n", sep = ""))
        } else {
          cat(paste("    ", Names[i], space, round(pval[i],
          digits = 4), " \n", sep = ""))
        }
      }
    }
  }
})
