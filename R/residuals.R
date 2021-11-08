# residuals from sfacross ----------

residuals.sfacross <- function(object, ...) {
  object$dataTable$mlResiduals
}

# residuals from lcmcross ----------

residuals.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(select(object$dataTable, "mlResiduals_c1",
      "mlResiduals_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mlResiduals_c1",
        "mlResiduals_c2", "mlResiduals_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3",
          "mlResiduals_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mlResiduals_c1",
          "mlResiduals_c2", "mlResiduals_c3",
          "mlResiduals_c4", "mlResiduals_c5"))
        }
      }
    }
  }
}
