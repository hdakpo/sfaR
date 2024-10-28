sfaRStartupMessage <- function()
{
  msg <- c(paste0(
    "           ****           *******  
          /**/           /**////** 
  ****** ******  ******  /**   /** 
 **//// ///**/  //////** /*******  
//*****   /**    ******* /**///**  
 /////**  /**   **////** /**  //** 
 ******   /**  //********/**   //**
//////    //    //////// //     //    version 1.0.1"),
    "\n\n* Please cite the 'sfaR' package as:\n", 
    "  Dakpo KH., Desjeux Y., Henningsen A., and Latruffe L. (2024). sfaR: Stochastic Frontier Analysis Using R. R package version 1.0.1.\n\n", 
    "See also: citation(\"sfaR\")\n\n", 
    "* For any questions, suggestions, or comments on the 'sfaR' package, you can contact directly the authors or visit:",
    "  https://github.com/hdakpo/sfaR/issues")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  msg <- sfaRStartupMessage()
  # if (!interactive())
  #   msg[1] <- paste("Package 'sfaR' version", "1.0.1")
  packageStartupMessage(msg)      
  invisible()
}