################################################################################
#                                                                              #
# data examples in sfaR package                                                #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# data list                                                                    #
#         -Data on Spanish dairy farms                                         #
#         -Data on U.S. electric power generation                              #
#         -Data on rice production in the Philippines                          #
#         -Data on Swiss railway companies                                     #
#         -Data on U.S. electricity generating plants                          #
#         -Data on world production                                            #
#         -Data on Norwegian dairy farms                                       #
# Data type: Cross sectional data & Panel data                                 #
#------------------------------------------------------------------------------#

#' Data on Spanish dairy farms
#'
#' This dataset contains six years of observations on 247 dairy farms in
#' northern Spain, drawn from 1993-1998. The original data consist in the farm
#' and year identifications, plus measurements on one output (i.e. milk), and
#' four inputs (i.e. cows, land, labor and feed).
#'
#' @details This dataset has been used in Alvarez \emph{et al.} (2004).  The data have
#' been normalized so that the logs of the inputs sum to zero over the 1,482
#' observations.
#'
#' @name dairyspain
#' @docType data
#' @format A data frame with 1,482 observations on the following 29 variables.
#' \describe{ \item{FARM}{Farm identification.} \item{AGEF}{Age of the farmer.}
#' \item{YEAR}{Year identification.} \item{COWS}{Number of milking cows.}
#' \item{LAND}{Agricultural area.} \item{MILK}{Milk production.}
#' \item{LABOR}{Labor.} \item{FEED}{Feed.} \item{YIT}{Log of \code{MILK}.}
#' \item{X1}{Log of \code{COWS}.} \item{X2}{Log of \code{LAND}.} \item{X3}{Log
#' of \code{LABOR}.} \item{X4}{Log of \code{FEED}.} \item{X11}{1/2 *
#' \code{X1}^2.} \item{X22}{1/2 * \code{X2}^2.} \item{X33}{1/2 * \code{X3}^2.}
#' \item{X44}{1/2 * \code{X4}^2.} \item{X12}{\code{X1} * \code{X2}.}
#' \item{X13}{\code{X1} * \code{X3}.} \item{X14}{\code{X1} * \code{X4}.}
#' \item{X23}{\code{X2} * \code{X3}.} \item{X24}{\code{X2} * \code{X4}.}
#' \item{X34}{\code{X3} * \code{X4}.} \item{YEAR93}{Dummy for \code{YEAR =
#' 1993}.} \item{YEAR94}{Dummy for \code{YEAR = 1994}.} \item{YEAR95}{Dummy for
#' \code{YEAR = 1995}.} \item{YEAR96}{Dummy for \code{YEAR = 1996}.}
#' \item{YEAR97}{Dummy for \code{YEAR = 1997}.} \item{YEAR98}{Dummy for
#' \code{YEAR = 1998}.} }
#' @references Alvarez, A., C. Arias, and W. Greene. 2004. Accounting for
#' unobservables in production models: management and inefficiency.
#' \emph{Econometric Society}, \bold{341}:1--20.
#' @source
#' \url{https://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}
#' @keywords datasets
#' @examples
#'
#' str(dairyspain)
#' summary(dairyspain)
NULL

#' Data on U.S. electric power generation
#'
#' This dataset is on electric power generation in the United States.
#'
#' @details The dataset is from Christensen and Greene (1976) and has also been used in
#' Greene (1990).
#'
#' @name electricity
#' @docType data
#' @format A data frame with 123 observations on the following 9 variables.
#' \describe{ \item{firm}{Firm identification.} \item{cost}{Total cost in 1970,
#' MM USD.} \item{output}{Output in million KwH.} \item{lprice}{Labor price.}
#' \item{lshare}{Labor's cost share.} \item{cprice}{Capital price.}
#' \item{cshare}{Capital's cost share.} \item{fprice}{Fuel price.}
#' \item{fshare}{Fuel's cost share.} }
#' @references Christensen, L.R., and W.H. Greene. 1976. Economies of scale in
#' US electric power generation. \emph{The Journal of Political Economy},
#' \bold{84}:655--676.
#'
#' Greene, W.H. 1990. A Gamma-distributed stochastic frontier model.
#' \emph{Journal of Econometrics}, \bold{46}:141--163.
#' @source \url{https://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}
#' @keywords datasets
#' @examples
#'
#' str(electricity)
#' summary(electricity)
NULL

#' Data on rice production in the Philippines
#'
#' This dataset contains annual data collected from 43 smallholder rice
#' producers in the Tarlac region of the Philippines between 1990 and 1997.
#'
#' This dataset is published as supplement to Coelli \emph{et al.} (2005).
#' While most variables of this dataset were supplied by the International Rice
#' Research Institute (IRRI), some were calculated by Coelli \emph{et al.}
#' (2005, see p. 325--326). The survey is described in Pandey \emph{et al.}
#' (1999).
#'
#' @name ricephil
#' @docType data
#' @format A data frame with 344 observations on the following 17 variables.
#' \describe{ \item{YEARDUM}{Time period (1= 1990, ..., 8 = 1997).}
#' \item{FARMERCODE}{Farmer code (1, ..., 43).} \item{PROD}{Output (tonnes of
#' freshly threshed rice).} \item{AREA}{Area planted (hectares).}
#' \item{LABOR}{Labor used (man-days of family and hired labor).}
#' \item{NPK}{Fertiliser used (kg of active ingredients).} \item{OTHER}{Other
#' inputs used (Laspeyres index = 100 for Farm 17 in 1991).}
#' \item{PRICE}{Output price (pesos per kg).} \item{AREAP}{Rental price of land
#' (pesos per hectare).} \item{LABORP}{Labor price (pesos per hired man-day).}
#' \item{NPKP}{Fertiliser price (pesos per kg of active ingredient).}
#' \item{OTHERP}{Price of other inputs (implicit price index).} \item{AGE}{Age
#' of the household head (years).} \item{EDYRS}{Education of the household head
#' (years).} \item{HHSIZE}{Household size.} \item{NADULT}{Number of adults in
#' the household.} \item{BANRAT}{Percentage of area classified as bantog
#' (upland) fields.} }
#' @references Coelli, T. J., Rao, D. S. P., O'Donnell, C. J., and Battese, G.
#' E. 2005. \emph{An Introduction to Efficiency and Productivity Analysis},
#' Springer, New York.
#'
#' Pandey, S., Masciat, P., Velasco, L, and Villano, R. 1999. Risk analysis of
#' a rainfed rice production system system in Tarlac, Central Luzon,
#' Philippines. \emph{Experimental Agriculture}, \bold{35}:225--237.
# @source Supplementary files for Coelli \emph{et al.} (2005),
# \url{http://www.uq.edu.au/economics/cepa/crob2005/software/CROB2005.zip}.
#' @keywords datasets
#' @examples
#'
#' str(ricephil)
#' summary(ricephil)
NULL

#' Data on Swiss railway companies
#'
#' This dataset is an unbalanced panel of 50 Swiss railway companies over the
#' period 1985-1997.
#'
#' The dataset is extracted from the annual reports of the Swiss Federal Office
#' of Statistics on public transport companies and has been used in Farsi
#' \emph{et al.} (2005).
#'
#' @name swissrailways
#' @docType data
#' @format A data frame with 605 observations on the following 42 variables.
#' \describe{ \item{ID}{Firm identification.} \item{YEAR}{Year identification.}
#' \item{NI}{Number of years observed.} \item{STOPS}{Number of stops in
#' network.} \item{NETWORK}{Network length (in meters).} \item{NARROW_T}{Dummy
#' variable for railroads with narrow track.} \item{RACK}{Dummy variable for
#' ‘rack rail’ in network.} \item{TUNNEL}{Dummy variable for network with
#' tunnels over 300 meters on average.} \item{T}{Time indicator, first year =
#' 0.} \item{Q2}{Passenger output – passenger km.} \item{Q3}{Freight output
#' – ton km.} \item{CT}{Total cost (1,000 Swiss franc).} \item{PL}{Labor
#' price.} \item{PE}{Electricity price.} \item{PK}{Capital price.}
#' \item{VIRAGE}{1 for railroads with curvy tracks.} \item{LNCT}{Log of
#' \code{CT}/\code{PE}.} \item{LNQ2}{Log of \code{Q2}.} \item{LNQ3}{Log of
#' \code{Q3}.} \item{LNNET}{Log of \code{NETWORK}/1000.} \item{LNPL}{Log of
#' \code{PL}/\code{PE}.} \item{LNPE}{Log of \code{PE}.} \item{LNPK}{Log of
#' \code{PK}/\code{PE}.} \item{LNSTOP}{Log of \code{STOPS}.} \item{MLNQ2}{Mean
#' of \code{LNQ2}.} \item{MLNQ3}{Mean of \code{LNQ3}.} \item{MLNNET}{Mean of
#' \code{LNNET}.} \item{MLNPL}{Mean of \code{LNPL}.} \item{MLNPK}{Mean of
#' \code{LNPK}.} \item{MLNSTOP}{Mean of \code{LNSTOP}.} }
#' @references Farsi, M., M. Filippini, and W. Greene. 2005. Efficiency
#' measurement in network industries: Application to the Swiss railway
#' companies. \emph{Journal of Regulatory Economics}, \bold{28}:69--90.
#' @source
#' \url{https://pages.stern.nyu.edu/~wgreene/Text/Edition7/tablelist8new.htm}
#'
#\url{http://people.stern.nyu.edu/wgreene/Microeconometrics.htm}
#' @keywords datasets
#' @examples
#'
#' str(swissrailways)
NULL

#' Data on U.S. electricity generating plants
#'
#' This dataset contains data on fossil fuel fired steam electric power
#' generation plants in the United States between 1986 and 1996.
#'
#' The dataset has been used in Kumbhakar \emph{et al.} (2014).
#'
#' @name utility
#' @docType data
#' @format A data frame with 791 observations on the following 11 variables.
#' \describe{ \item{firm}{Plant identification.} \item{year}{Year
#' identification.} \item{y}{Net-steam electric power generation in
#' megawatt-hours.} \item{regu}{Dummy variable which takes a value equal to 1
#' if the power plant is in a state which enacted legislation or issued a
#' regulatory order to implement retail access during the sample period, and 0
#' otherwise.} \item{k}{Capital stock.} \item{labor}{Labor and maintenance.}
#' \item{fuel}{Fuel.} \item{wl}{Labor price.} \item{wf}{Fuel price.}
#' \item{wk}{Capital price.} \item{tc}{Total cost.} }
#' @references Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. \emph{A
#' Practitioner's Guide to Stochastic Frontier Analysis Using Stata}. Cambridge
#' University Press.
#' @source
#' \url{https://sites.google.com/view/sfbook-stata/home}
#' @keywords datasets
#' @examples
#'
#' str(utility)
#' summary(utility)
NULL

#' Data on world production
#'
#' This dataset provides information on production related variables for
#' eighty-two countries over the period 1960–1987.
#'
#' The dataset is from the World Bank STARS database and has been used in
#' Kumbhakar \emph{et al.} (2014).
#'
#' @name worldprod
#' @docType data
#' @format A data frame with 2,296 observations on the following 12 variables.
#' \describe{ \item{country}{Country name.} \item{code}{Country
#' identification.} \item{yr}{Year identification.} \item{y}{GDP in 1987 U.S.
#' dollars.} \item{k}{Physical capital stock in 1987 U.S. dollars.}
#' \item{l}{Labor (number of individuals in the workforce between the age of
#' 15 and 64).} \item{h}{Human capital-adjusted labor.} \item{ly}{Log of
#' \code{y}.} \item{lk}{Log of \code{k}.} \item{ll}{Log of \code{l}.}
#' \item{lh}{Log of \code{h}.} \item{initStat}{Log of the initial capital to
#' labor ratio of each country, \code{lk} - \code{ll}, measured at the
#' beginning of the sample period.} }
#' @references Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. \emph{A
#' Practitioner's Guide to Stochastic Frontier Analysis Using Stata}. Cambridge
#' University Press.
#' @source
#' \url{https://sites.google.com/view/sfbook-stata/home}
#' @keywords datasets
#' @examples
#'
#' str(worldprod)
#' summary(worldprod)
NULL

#' Data on Norwegian dairy farms
#'
#' This dataset contains nine years (1998-2006) of information on Norwegian dairy
#' farms.
#' 
#' @name dairynorway
#' @docType data
#' @format A data frame with 2,727 observations on the following 23 variables.
#' \describe{ \item{farmid}{Farm identification.} \item{year}{Year
#' identification.} \item{y1}{Milk sold (1000 liters).} \item{y2}{Meat (1000 NOK).} 
#' \item{y3}{Support payments (1000 NOK).} \item{y4}{Other outputs (1000 NOK).} 
#' \item{p1}{Milk price (NOK/liter).} \item{p2}{Meat price (cattle index).} 
#' \item{p3}{Support payments price (CP index).} \item{p4}{Other outputs price index.} 
#' \item{x1}{Land (decare (daa) = 0.1 ha).} \item{x2}{Labour (1000 hours).} 
#' \item{x3}{Purchase feed (1000 NOK).} \item{x4}{Other variable costs (1000 NOK).} 
#' \item{x5}{Cattle capital (1000 NOK).} \item{x6}{Other capital (1000 NOK).} 
#' \item{w1}{Land price (NOK/daa).} \item{w2}{Labour price (NOK/hour).} 
#' \item{w3}{Feed price index.} \item{w4}{Other variable cost index.} 
#' \item{w5}{Cattle capital rent.} \item{w6}{Other capital rent and depreciation.} 
#' \item{tc}{Total cost.}}
#' @references Kumbhakar, S.C., H.J. Wang, and A. Horncastle. 2014. \emph{A
#' Practitioner's Guide to Stochastic Frontier Analysis Using Stata}. Cambridge
#' University Press.
#' @source
#' \url{https://sites.google.com/view/sfbook-stata/home}
#' @keywords datasets
#' @examples
#'
#' str(dairynorway)
#' summary(dairynorway)
NULL

