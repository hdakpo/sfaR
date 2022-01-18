# globalVariables(c("EffRes", "LogNorm"))

.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste0("* Please cite the 'sfaR' package as:\n", 
        "  Dakpo KH., Desjeux Y., and Latruffe L. (2022). sfaR: Stochastic Frontier Analysis Routines. R package version 1.0.0.\n", 
        "See also: citation(\"sfaR\")\n\n", 
        "* For any questions, suggestions, or comments on the 'sfaR' package, please make use of Tracker facilities at:\n",
        "  https://github.com/hdakpo/sfaR/issues"))
}
