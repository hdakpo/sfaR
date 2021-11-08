# skewness test for sfacross ----------

skewnessTest <- function(object, test = "agostino") {
if (! inherits(object,"sfacross")) {
        stop("Argument 'object' must be of class 'sfacross'")
}
  if (test == "agostino") {
    object$AgostinoTest
  } else {
    if (test == "coelli") {
	tt=list("data.name" = deparse(substitute(object)),
	"statistic" = object$CoelliM3Test["z"],
			"p.value" = object$CoelliM3Test["p.value"],
			"method" = "## Coelli's test ##")
	class(tt) <- "htest"
	return(tt)
		} else {
      stop("argument 'test' must be either 'agostino', or 'coelli'",
           call. = FALSE)
    }
  }
}

# # Previous version
# skewnessTest <- function(object, test = "agostino") {
# if (! inherits(object,"sfacross")) {
        # stop("Argument 'object' must be of class 'sfacross'")
# }
  # if (test == "agostino") {
    # object$AgostinoTest
  # } else {
    # if (test == "coelli") {
        # cat("## Coelli's test ##\n", sep = "")
        # cat("\nTest Results:\n", sep = "")
        # cat("  STATISTIC (z):", round(object$CoelliM3Test["z"], digits = 4), "\n")
        # cat("  P.VALUE:", format.pval(object$CoelliM3Test["p.value"], digits = 4), "\n")
    # } else {
      # stop("argument 'test' must be either 'agostino', or 'coelli'",
           # call. = FALSE)
    # }
  # }
# }
