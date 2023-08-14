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
//////    //    //////// //     //    version 2.0.0.9000"),
"\n\n* Please cite the 'sfaR' package as:\n", 
"  Dakpo KH., Desjeux Y., Henningsen A., and Latruffe L. (2023). sfaR: Stochastic Frontier Analysis Using R. R package version 2.0.0.9000.\n\n", 
"See also: citation(\"sfaR\")\n\n", 
"* For any questions, suggestions, or comments on the 'sfaR' package, you can contact directly the authors or visit:",
"  https://github.com/hdakpo/sfaR/issues")
  return(msg)
}

.onAttach <- function(lib, pkg) {
  msg <- sfaRStartupMessage()
  # if (!interactive())
  #   msg[1] <- paste("Package 'sfaR' version", "2.0.0.9000")
  packageStartupMessage(msg)      
  invisible()
}
