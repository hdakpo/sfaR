# globalVariables(c("EffRes", "LogNorm"))

.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste0("* Please cite the 'sfaR' package as:\n", 
        "  Dakpo KH., Desjeux Y. and Latruffe L. (2021). sfaR: Stochastic Frontier Analysis using R. R package version 0.1.0.\n", 
        "See also: citation(\"sfaR\")\n\n", 
        "* For any questions, suggestions, or comments on the 'sfaR' package, please make use of Tracker facilities at:\n",
        "  https://r-forge.r-project.org/projects/sfar/"))
}
