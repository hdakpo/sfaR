# fitted values for sfacross ----------

fitted.sfacross <- function(object, ...) {
  object$dataTable$mleFitted
}

# fitted values for lcmcross ----------

fitted.lcmcross <- function(object, ...) {
  if (object$nClasses == 2) {
    data.frame(select(object$dataTable, "mlFitted_c1",
      "mlFitted_c2"))
  } else {
    if (object$nClasses == 3) {
      data.frame(select(object$dataTable, "mlFitted_c1",
        "mlFitted_c2", "mlFitted_c3"))
    } else {
      if (object$nClasses == 4) {
        data.frame(select(object$dataTable, "mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3",
          "mlFitted_c4"))
      } else {
        if (object$nClasses == 5) {
          data.frame(select(object$dataTable, "mlFitted_c1",
          "mlFitted_c2", "mlFitted_c3",
          "mlFitted_c4", "mlFitted_c5"))
        }
      }
    }
  }
}
