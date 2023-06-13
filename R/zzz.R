.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste0("* Please cite the 'sfaR' package as:\n", 
        "  Dakpo KH., Desjeux Y., Henningsen A., and Latruffe L. (2023). sfaR: Stochastic Frontier Analysis Routines. R package version 1.0.0.\n\n", 
        "See also: citation(\"sfaR\")\n\n", 
        "* For any questions, suggestions, or comments on the 'sfaR' package, you can contact directly the authors or visit:\n",
        "  https://github.com/hdakpo/sfaR/issues"))
}
