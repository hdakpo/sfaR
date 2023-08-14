################################################################################
#                                                                              #
# R functions for the sfaR package                                             #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# Model: Deterministic/Stochastic Metafrontier Analysis                        #
# Data: Cross sectional data & Pooled data                                     #
#------------------------------------------------------------------------------#

#' Metafrontier estimation using cross-sectional or pooled data
#' 
#' @description
#' \code{\link{sfametacross}} is a symbolic formula-based function for the
#' estimation of stochastic and deterministic metafrontiers models in the case 
#' of cross-sectional or pooled cross-sectional data, using maximum (simulated) 
#' likelihood - M(S)L, linear or quadratic programming, or simulation depending 
#' on the case.
#'
#' The function accounts for heteroscedasticity in both one-sided and two-sided
#' error terms as in Reifschneider and Stevenson (1991), Caudill and Ford
#' (1993), Caudill \emph{et al.} (1995) and Hadri (1999), but also
#' heterogeneity in the mean of the pre-truncated distribution as in Kumbhakar
#' \emph{et al.} (1991), Huang and Liu (1994) and Battese and Coelli (1995).
#'
#' Depending on the model, ten distributions are possible for the one-sided 
#' error term and eleven optimization algorithms are available.
#' 
#' @aliases sfametacross bread.sfametacross estfun.sfametacross 
#' print.sfametacross
#'
#' @param formula A symbolic description of the model to be estimated based on
#' the generic function \code{formula} (see section \sQuote{Details}).
#' @param muhet A one-part formula to consider heterogeneity in the mean of the
#' pre-truncated distribution (see section \sQuote{Details}).
#' @param uhet A one-part formula to consider heteroscedasticity in the
#' one-sided error variance (see section \sQuote{Details}).
#' @param vhet A one-part formula to consider heteroscedasticity in the
#' two-sided error variance (see section \sQuote{Details}).
#' @param ghet A one-part formula specifying the group variable for the 
#' metafrontier estimation
#' @param logDepVar Logical. Informs whether the dependent variable is logged
#' (\code{TRUE}) or not (\code{FALSE}). Default = \code{TRUE}.
#' @param data The data frame containing the data.
#' @param subset An optional vector specifying a subset of observations to be
#' used in the optimization process.
#' @param weights An optional vector of weights to be used for weighted 
#' log-likelihood. Should be \code{NULL} or numeric vector with positive values. 
#' When \code{NULL}, a numeric vector of 1 is used.
#' @param wscale Logical. When \code{weights} is not \code{NULL}, a scaling 
#' transformation is used such that the \code{weights} sum to the sample 
#' size. Default \code{TRUE}. When \code{FALSE} no scaling is used.
#' @param S If \code{S = 1} (default), a production (profit) frontier is
#' estimated: \eqn{\epsilon_i = v_i-u_i}. If \code{S = -1}, a cost frontier is
#' estimated: \eqn{\epsilon_i = v_i+u_i}.
#' @param modelType Metafrontier models that can be estimated. Five models are
#' allowed. \code{'hhl14'} estimates the stochastic metafrontier following
#' Huang et al. (2014) - Default. \code{'bpo04a'} and \code{'bpo04b'} estimate
#' the deterministic metafrontier following Battese et al. (2004), and using
#' linear and quadratic programming, respectively. \code{'aos17a'} and 
#' \code{'aos17b'} estimate the stochastic metafrontier based on 
#' Amsler et al. (2017) using simulation. \code{'aos17a'} computes the meta-
#' distance assuming the unconditional inefficiency term, while \code{'aos17b'}
#' assumes that the inefficiency term is conditional on the overall error
#' term. \code{'aos17a'} and \code{'aos17b'} assumes that the two-sided error
#' terms \eqn{v_1, v_2, \cdots, v_G}  are mutually independent for the \eqn{G}
#' groups. Models \code{'aos17c'} and \code{'aos17d'} reproduces \code{'aos17a'}
#' and \code{'aos17b'}, respectively by assuming an equiproportionnate 
#' correlation. The user is prompted to provide a value for the correlation when
#' any of these two models is chosen.
#' @param udist Character string. Default = \code{'hnormal'}. Distribution
#' specification for the one-sided error term. 10 different distributions are
#' available: \itemize{ \item \code{'hnormal'}, for the half normal
#' distribution (Aigner \emph{et al.} 1977, Meeusen and Vandenbroeck 1977)
#' \item \code{'exponential'}, for the exponential distribution \item
#' \code{'tnormal'} for the truncated normal distribution (Stevenson 1980)
#' \item \code{'rayleigh'}, for the Rayleigh distribution (Hajargasht 2015)
#' \item \code{'uniform'}, for the uniform distribution (Li 1996, Nguyen 2010)
#' \item \code{'gamma'}, for the Gamma distribution (Greene 2003) \item
#' \code{'lognormal'}, for the log normal distribution (Migon and Medici 2001,
#' Wang and Ye 2020) \item \code{'weibull'}, for the Weibull distribution
#' (Tsionas 2007) \item \code{'genexponential'}, for the generalized
#' exponential distribution (Papadopoulos 2020) \item \code{'tslaplace'}, for
#' the truncated skewed Laplace distribution (Wang 2012). } For 
#' \code{'modelType'} = 'aos17b', only \code{'udist'} = 'hnormal' is available.
#' @param start Numeric vector. Optional starting values for the maximum
#' likelihood (ML) estimation.
#' @param randStart Logical. Define if random starting values should be used for
#' M(S)L estimation. New starting values are obtained as old ones + draws from
#' normal distribution with std. deviation of 0.01. \code{'seed'} is not 
#' considered here, then each run will provide different starting values 
#' (unless a seed is set by the user before the run).
#' @param method Optimization algorithm used for the estimation.  Default =
#' \code{'bfgs'}. 11 algorithms are available: \itemize{ \item \code{'bfgs'},
#' for Broyden-Fletcher-Goldfarb-Shanno (see
#' \code{\link[maxLik:maxBFGS]{maxBFGS}}) \item \code{'bhhh'}, for
#' Berndt-Hall-Hall-Hausman (see \code{\link[maxLik:maxBHHH]{maxBHHH}}) \item
#' \code{'nr'}, for Newton-Raphson (see \code{\link[maxLik:maxNR]{maxNR}})
#' \item \code{'nm'}, for Nelder-Mead (see \code{\link[maxLik:maxNM]{maxNM}}) 
#' \item \code{'cg'}, for Conjugate Gradient (see 
#' \code{\link[maxLik:maxCG]{maxCG}}) \item \code{'sann'}, for Simulated 
#' Annealing (see \code{\link[maxLik:maxSANN]{maxSANN}}) \item \code{'ucminf'}, 
#' for a quasi-Newton type optimization with BFGS updating of the inverse 
#' Hessian and soft line search with a trust region type monitoring of
#' the input to the line search algorithm (see 
#' \code{\link[ucminf:ucminf]{ucminf}})
#' \item \code{'mla'}, for general-purpose optimization based on
#' Marquardt-Levenberg algorithm (see \code{\link[marqLevAlg:mla]{mla}})
#' \item \code{'sr1'}, for Symmetric Rank 1 (see
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item \code{'sparse'}, 
#' for trust regions and sparse Hessian (see 
#' \code{\link[trustOptim:trust.optim]{trust.optim}}) \item
#' \code{'nlminb'}, for optimization using PORT routines (see
#' \code{\link[stats:nlminb]{nlminb}})}
#' @param hessianType Integer. If \code{1} (Default), analytic Hessian is
#' returned for all the distributions. If \code{2},
#' bhhh Hessian is estimated (\eqn{g'g}).
#' @param metaSim Integer. Number of simulations conducted for the estimation of
#' the stochastic metafrontier in the case of models \code{'aos17a'} and 
#' \code{'aos17b'}. For models \code{'bpo04a'} and \code{'bpo04b'}, 
#' \code{'metaSim'} is used for the simulation of the parameters standard 
#' errors. Default 5000.
#' @param simType Character string. If \code{simType = 'halton'} (Default),
#' Halton draws are used for maximum simulated likelihood (MSL). If
#' \code{simType = 'ghalton'}, Generalized-Halton  draws are used for MSL. If
#' \code{simType = 'sobol'} or \code{simType = 'rsobol'}, Sobol draws or 
#' randomized Sobol draws are used for MSL, respectively. If 
#' \code{simType = 'richtmyer'}, or \code{simType = 'rrichtmyer'}, 
#' Richtmyer sequence or randomized Richtmyer sequence is used for MSL 
#' estimation, respectively. If \code{simType = 'mlhs'}, modified latin 
#' hypercube sampling is performed. If \code{simType = 'uniform'}, uniform draws 
#' are used for MSL. (see section \sQuote{Details}).
#' @param Nsim Number of draws for MSL. Default 100.
#' @param prime Prime number considered for quasi-random sequences. The first
#' 10,000 primes are available. Default = \code{2}.
#' @param burn Number of the first observations discarded in the case of Halton,
#' Richtmyer, and Sobol draws. For the randomized versions, the burn is not
#' considered. Default = \code{10}.
#' @param antithetics Logical. Default = \code{FALSE}. If \code{TRUE},
#' antithetics counterpart of the uniform draws is computed. (see section
#' \sQuote{Details}).
#' @param seed Numeric. Seed for the random draws.
#' @param itermax Maximum number of iterations allowed for optimization.
#' Default = \code{2000}.
#' @param printInfo Logical. Print information during optimization. Default =
#' \code{FALSE}.
#' @param tol Numeric. Convergence tolerance. Default = \code{1e-12}.
#' @param gradtol Numeric. Convergence tolerance for gradient. Default =
#' \code{1e-06}.
#' @param stepmax Numeric. Step max for \code{ucminf} algorithm. Default =
#' \code{0.1}.
#' @param qac Character. Quadratic Approximation Correction for \code{'bhhh'}
#' and \code{'nr'} algorithms. If \code{'stephalving'}, the step length is
#' decreased but the direction is kept. If \code{'marquardt'} (default), the
#' step length is decreased while also moving closer to the pure gradient
#' direction. See \code{\link[maxLik:maxBHHH]{maxBHHH}} and
#' \code{\link[maxLik:maxNR]{maxNR}}.
#' @param x an object of class sfametacross (returned by the function 
#' \code{\link{sfametacross}}).
#' @param ... additional arguments of frontier are passed to sfametacross; 
#' additional arguments of the print, bread, and estfun methods are currently
#'  ignored.
#'  
#' @details
#' The stochastic frontier model for the cross-sectional data is defined for 
#' unit \eqn{i} belonging to group \eqn{d_i} as
#' 
#' \deqn{y_{i} = \mathbf{x}'_i\bm{\beta}_{d_i} + v_{i,{d_i}} - Su_{i,{d_i}}}
#' 
#' where \eqn{S=1} in the case of a production function, and \eqn{S=-1} in the 
#' case of a cost function. Even though, we only observe \eqn{g = d_i}, the 
#' concept of the metafrontier is based on assessing by how much unit \eqn{i} 
#' could have produced if it had used the technology of a different group. 
#' Conceptually, we have 
#' 
#' \deqn{y_{ig} = \mathbf{x}'_i\bm{\beta}_g + v_{ig} - Su_{ig} \;\; \text{with} 
#' \;\; g= 1, \cdots, G}
#' 
#' The first step is the estimation of each group frontier. The second step, 
#' which is the metafrontier estimation varies depending on the methodology.
#' 
#' Battese et al. (2004) suggested a deterministic way for estimating the 
#' metafrontier. The metafrontier is defined as
#' 
#' \deqn{y_i^* =\mathbf{x}'_i\bm{\beta}^*}
#' 
#' such that 
#' 
#' \deqn{\mathbf{x}'_i\bm{\beta}^* \geq \mathbf{x}'_i\bm{\beta}_g}
#' 
#' The parameters \eqn{\bm{\beta}^*} can be obtained by solving either a linear 
#' program:
#' 
#' \deqn{\begin{matrix} 
#' \min L \equiv \sum_{i = 1} \left(\mathbf{x}'_i\bm{\beta}^* - 
#' \mathbf{x}'_i\bm{\beta}_g\right) \\[1em]
#' \text{s.t.} \; \mathbf{x}'_i\bm{\beta}^* \geq \mathbf{x}'_i\bm{\beta}_g
#'  \end{matrix}}
#'  
#' or a quadratic program:
#' 
#' \deqn{\begin{matrix} 
#' \min L \equiv \sum_{i = 1} \left(\mathbf{x}'_i\bm{\beta}^* - 
#' \mathbf{x}'_i\bm{\beta}_g\right)^2 \\[1em]
#' \text{s.t.} \; \mathbf{x}'_i\bm{\beta}^* \geq \mathbf{x}'_i\bm{\beta}_g
#'  \end{matrix}}
#' 
#' The thechnology gap ratio is obtained as (assuming that \eqn{y_i} is 
#' expressed in logarithm)
#' 
#' \deqn{TGR_i = \frac{\exp{\left(\mathbf{x}'_i\bm{\beta}_g\right)}}{
#' \exp{\left(\mathbf{x}'_i\bm{\beta}^*\right)}}}
#' 
#' Huang et al. (2014) suggested a stochastic version of the metafrontier. Let's
#' assume that \eqn{\ln f_i^g=\mathbf{x}'_i\bm{\beta}_g} and 
#' \eqn{\ln f_i^M = \mathbf{x}'_i\bm{\beta}^*}. Huang et al. (2014) assumed the 
#' following relation
#' 
#' \deqn{\ln f_i^g = \ln f_i^M - u_{i,M}}
#' 
#' where \eqn{\exp{\left(u_{i,M}\right)}} is the TGR.
#' 
#' Since \eqn{\ln f_i^g} is not observed, Huang et al. (2014) suggested the 
#' following estimation
#' 
#' \deqn{\ln \hat{f}_i^g = \ln f_i^M - u_{i,M} + v_{i,M}}
#' 
#' where \eqn{\ln \hat{f}_i^g} is the predicted frontier obtained from the 
#' individual group frontier.
#' 
#' The new quasi maximum likelihood estimator requires a correction of the 
#' standard error in order to account for heteroscedasticity. The robust
#' sandwich-form of the variance-covariance matrix estimator can be used to 
#' correct for the standard errors.
#' 
#' A new stochastic metafrontier has also been suggested in 
#' Amsler et al. (2017). Based on the previous developments, The stochastic 
#' frontiers corresponding to each group are defined as
#' 
#' \deqn{f_{ig} = \mathbf{x}'_i\bm{\beta}_g + v_{ig} \;\; \text{with} 
#' \;\; g= 1, \cdots, G}
#' 
#' Amsler et al. (2017) defined the stochastic metafrontier in the case of a 
#' production function as 
#' 
#' \deqn{f_i = \max\left[f_{i1}, \cdots, f_{iG}\right]}
#' 
#' The following decomposition can then be obtained 
#' 
#' \deqn{\left(f_i - y_i\right)=\left(f_{i,d_i} - y_i\right) + 
#' \left(f_i - f_{i,d_i}\right)}
#' 
#' This equation is also equivalent to
#' 
#' \deqn{U_i=U_{i, d_i} + M_{i, d_i}}
#' 
#' where \eqn{U_{i, d_i}} is the one-sided technical inefficiency for unit 
#' \eqn{i} in the stochastic frontier model for group \eqn{d_i}, and 
#' \eqn{M_{i, d_i}}, which is the metafrontier distance, can be evaluated
#' based on different assumptions in relation to the total error term 
#' (\eqn{\epsilon = v-u}).
#' 
#' For simplicity, let's write
#' 
#' \deqn{U=U_{d} + M_{d}}
#' 
#' and assume that \eqn{U}, \eqn{U_d}, and \eqn{M_d} can be evaluated 
#' unconditionnally to \eqn{\epsilon}.
#' 
#' Therefore, in terms of expected values, we have
#' 
#' \deqn{\mu = \mu_d + \tau_d}
#' 
#' In the case \eqn{u_d} follows a half normal distribution, we have
#' 
#' \deqn{\mu_d = E\left[U_d\right] = \sqrt{\frac{2}{\pi}}\sigma_{u,d}}
#' 
#' On the other hand \eqn{\tau_d = E\left[M_d\right]}. As underlined by 
#' Amsler et al. 2017, treating the frontiers as stochastic given 
#' \eqn{\mathbf{x}} yields 
#' 
#' \deqn{f_s \sim \left[\mathbf{x}'\bm{\beta} + \mathcal{N}
#' \left(0, \sigma_{v,s}^2\right)\right]}
#' 
#' or equivalently
#' 
#' \deqn{f_s \sim \mathcal{N}\left(\mathbf{x}'\bm{\beta}, 
#' \sigma_{v,s}^2\right)}
#' 
#' Then \eqn{f = \max\left(f_1, \cdots, f_S\right)} is the maximum of a set of 
#' normal distributions. The expectation of a set of \eqn{P} normal random 
#' variables is not known. Therefore Amsler et al. 2017 suggested the use of 
#' simulation.
#' 
#' For replication \eqn{r (r = 1, \cdots, R)}, take draws 
#' \eqn{f_1^{(r)}, \cdots, f_S^{(r)}} from the normal distributions of 
#' \eqn{f_1, \cdots, f_S}, and calculate the metafrontier 
#' \eqn{f^{(r)} = \max(f_1^{(r)}, \cdots, f_S^{(r)})}. Then 
#' \eqn{M_d^{(r)} = f^{(r)}-f_d^{(r)}}. Finally 
#' 
#' \deqn{\tau_d = \sum_{r = 1}^{R}M_d^{(r)}}
#' 
#' Conditioning now on \eqn{\epsilon}, we have 
#' \eqn{\mu_d^*=E\left[U|\epsilon\right]}, which can be obtained following 
#' Jondrow et al. (1982). On the other hand, for 
#' \eqn{\tau_d^*=E\left[M_d|\epsilon\right]}, we need to derive the distribution 
#' of \eqn{v|\epsilon}. 
#' 
#' In the case of the half normal distribution, we have
#' 
#' \deqn{f(v|\epsilon)=\frac{1}{\sigma^*\Phi\left(-a^*\right)\sqrt{2\pi}}
#' \exp\left\{-\frac{1}{2}\left[\frac{v-\mu}{\sigma_*}\right]^2\right\}}
#' 
#' where \eqn{-\frac{\epsilon - \mu}{\sigma^*}=-a^*} and 
#' \eqn{\mu = \frac{\sigma_v^2}{\sigma_u^2 + \sigma_v^2}\epsilon}
#' 
#' The previous pdf is equivalent to a normal distribution 
#' \eqn{\mathcal{N}\left(\mu, \sigma_*^2\right)} truncated on the left at 
#' \eqn{\epsilon}.
#' 
#' As previously, the evaluation of \eqn{\tau_d^*} requires simulation. It  is 
#' worth noting that for state \eqn{d}, conditioning on \eqn{\epsilon} matters, 
#' while for other states \eqn{s\neq d}, conditioning on \eqn{\epsilon} does not 
#' matter. Therefore, for each replication \eqn{r}, and for \eqn{s\neq d}, 
#' draw \eqn{f_s^{(r)}} as explained in the previous subsection. However, for 
#' group \eqn{d}, draw \eqn{f_d^{(r)}} equal to 
#' \eqn{\mathbf{x}'\bm{\beta}} plus a draw from the distribution of 
#' \eqn{v|\epsilon}.
#' 
#' Both the approaches suggested by Amsler et al. (2017) can be extended to 
#' assuming that the \eqn{v_{is} are not independent.}
#'
#' @return \code{\link{sfametacross}} returns a list of class 
#' \code{'sfametacross'} containing the following elements:
#'
#' \item{call}{The matched call.}
#'
#' \item{formula}{A multi-part formula of the estimated model.}
#' 
#' \item{S}{The argument \code{'S'}. See the section \sQuote{Arguments}.}
#'
#' \item{typeSfa}{Character string. 'Stochastic Production/Profit Frontier, e =
#' v - u' when \code{S = 1} and 'Stochastic Cost Frontier, e = v + u' when
#' \code{S = -1}.}
#' 
#' \item{Nobs}{Vector containing the number of observations for each group, in 
#' addition to the metafrontier estimation.}
#' 
#' \item{Ngroup}{Number of groups for the metafrontier.}
#' 
#' \item{name_meta_var}{Character sting indicating the name of variable 
#' containing the groups of observations.}
#' 
#' \item{nXvar}{Number of explanatory variables in the production or cost
#' frontier.}
#'
#' \item{nmuZUvar}{Number of variables explaining heterogeneity in the
#' truncated mean, only if \code{udist = 'tnormal'} or \code{'lognormal'}.}
#' 
#' \item{Yvarm}{Numeric vector of group predicted frontiers for estimating
#' metafrontiers.}
#' 
#' \item{logDepVar}{The argument \code{'logDepVar'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{nuZUvar}{Number of variables explaining heteroscedasticity in the
#' one-sided error term.}
#'
#' \item{nvZVvar}{Number of variables explaining heteroscedasticity in the
#' two-sided error term.}
#'
#' \item{nParm}{Total number of parameters estimated.}
#'
#' \item{udist}{The argument \code{'udist'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{startVal}{Numeric matrix containing the starting values used for the 
#' estimation of the group frontiers and also the metafrontier in the case
#' \code{'modelType'} = 'hll14'.}
#' 
#' \item{modelType}{Character string indicating the type of metafrontier that
#' is estimated. See the section
#' \sQuote{Arguments}.}
#' 
#' \item{olsParam}{List of numeric vectors. OLS estimates for group frontiers, 
#' and metafrontier in the case \code{'modelType'} = 'hll14'.}
#'
#' \item{olsStder}{List of numeric vectors. Standard errors of OLS estimates
#' for group frontiers, and metafrontier in the case 
#' \code{'modelType'} = 'hll14'.}
#'
#' \item{olsSigmasq}{Numeric. Estimated variance of OLS random error for group 
#' frontiers, and metafrontier in the case \code{'modelType'} = 'hll14'.}
#'
#' \item{olsLoglik}{List. Log-likelihood value of OLS estimation for group 
#' frontiers, and metafrontier in the case \code{'modelType'} = 'hll14'.}
#'
#' \item{olsSkew}{List. Skewness of the residuals of the OLS estimation.}
#'
#' \item{olsM3Okay}{List of character string. Indicate whether the residuals of 
#' the OLS estimation have the expected skewness for group frontiers, 
#' and metafrontier in the case \code{'modelType'} = 'hll14'.}
#'
#' \item{CoelliM3Test}{List. Coelli's test for OLS residuals skewness for group 
#' frontiers, and metafrontier in the case \code{'modelType'} = 'hll14'. 
#' (See Coelli, 1995).}
#'
#' \item{AgostinoTest}{List. D'Agostino's test for OLS residuals skewness for 
#' group frontiers, and metafrontier in the case \code{'modelType'} = 'hll14'. 
#' (See D'Agostino and Pearson, 1973).}
#' 
#' \item{isWeights}{Logical. If \code{TRUE} weighted log-likelihood is
#' maximized.}
#'
#' \item{optType}{Optimization algorithm(s) used for group frontiers and 
#' metafrontier.}
#'
#' \item{nIter}{Integer vector. Number of iterations of the ML estimation.}
#'
#' \item{optStatus}{String vector. Optimization algorithm termination message.}
#'
#' \item{startLoglik}{Numeric vector. Log-likelihood at the starting values.}
#'
#' \item{mlLoglik}{Numeric vector. Log-likelihood value of the M(S)L 
#' estimation.}
#'
#' \item{mlParam}{Numeric Matrix. Parameters obtained from M(S)L estimation, and 
#' LP or QP estimation.}
#' 
#'  \item{dataTable}{List of data frames. Each data frame contains information 
#'  on data used for optimization (group frontiers and metafrontiers) along with 
#'  residuals and fitted values of the OLS and M(S)L estimations, and the 
#'  individual observation log-likelihood (whenever likelihood is used). 
#'  When \code{weights} is specified an additional variable is also provided in 
#' \code{dataTable}.}
#'
#' \item{gradient}{Numeric matrix. Each variable gradient of the M(S)L 
#' estimation.}
#'
#' \item{gradL_OBS}{List of numeric matrices. Each variable individual 
#' observation gradient of the M(S)L estimation.}
#'
#' \item{gradientNorm}{Numeric vector. Gradient norm of the M(S)L estimation.}
#'
#' \item{invHessian}{List of numeric matrices. Covariance matrix of the 
#' parameters obtained from the M(S)L estimation.}
#' 
#' \item{conditionNums}{List of numeric matrices. Condition number adding 
#' columns one by one for the different groups and in the case 
#' \code{'modelType'} = 'hhl14', for the metafrontier.}
#' 
#' \item{hessianType}{The argument \code{'hessianType'}. See the section
#' \sQuote{Arguments}.}
#' 
#' \item{metaSim}{Number of simulations used in the case \code{'modelType'} = 
#' 'bpo04a', 'bpo04b', 'aos17a', or 'aos17b'.}
#' 
#' \item{lpSimRes}{Matrix of parameters obtained from the LP run in the 
#' case \code{'modelType'} = 'bpo04a'. Used to obtain standard errors of 
#' parameters.}
#' 
#' \item{qpSimRes}{Matrix of parameters obtained from the QP run in the 
#' case \code{'modelType'} = 'bpo04b'. Used to obtain standard errors of 
#' parameters.}
#' 
#' \item{MdMat}{Matrix of metafrontier distance as described in 
#' Amsler et al. (2017).}
#' 
#' \item{Mud}{Numeric vector for the unconditional inefficiency when
#' \code{'modelType'} = 'aos17a'.}
#' 
#' \item{groupFrontierMat}{Matrix of simulated group frontiers in the case
#' \code{'modelType'} = 'aos17a', or 'aos17b'.}
#' 
#' \item{metaFrontierMat}{Matrix of simulated metafrontier in the case
#' \code{'modelType'} = 'aos17a', or 'aos17b'.}
#' 
#' \item{mlDate}{Date and time of the estimated models.}
#'
#' \item{simDist}{The argument \code{'simDist'}, only if \code{udist =
#' 'gamma'}, \code{'lognormal'} or , \code{'weibull'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{Nsim}{The argument \code{'Nsim'}, only if \code{udist = 'gamma'},
#' \code{'lognormal'} or , \code{'weibull'}. See the section
#' \sQuote{Arguments}.}
#'
#' \item{FiMat}{List of matrices of random draws used for MSL, only if 
#' \code{udist = 'gamma'}, \code{'lognormal'} or , \code{'weibull'}.}
#' 
#' @export
#' 
#' @references Aigner, D., Lovell, C. A. K., and Schmidt, P. 1977. Formulation
#' and estimation of stochastic frontier production function models.
#' \emph{Journal of Econometrics}, \bold{6}(1), 21--37.
#' 
#' Amsler, C., O'Donnell, C. J., Schmidt, P., 2017. Stochastic metafrontiers, 
#' \emph{Econometric Reviews}, \bold{36}, 1007-1020.
#'
#' Battese, G. E., and Coelli, T. J. 1995. A model for technical inefficiency
#' effects in a stochastic frontier production function for panel data.
#' \emph{Empirical Economics}, \bold{20}(2), 325--332.
#' 
#' Battese, G. E., Rao, D. S. P., O'Donnell, C. J., 2004. A metafrontier 
#' production function for estimation of technical efficiencies and technology 
#' gaps for firms operating under different technologies, 
#' \emph{Journal of Productivity Analysis}, \bold{21}, 91-103.
#'
#' Caudill, S. B., and Ford, J. M. 1993. Biases in frontier estimation due to
#' heteroscedasticity. \emph{Economics Letters}, \bold{41}(1), 17--20.
#'
#' Caudill, S. B., Ford, J. M., and Gropper, D. M. 1995. Frontier estimation
#' and firm-specific inefficiency measures in the presence of
#' heteroscedasticity. \emph{Journal of Business & Economic Statistics},
#' \bold{13}(1), 105--111.
#'
#' Coelli, T. 1995. Estimators and hypothesis tests for a stochastic frontier
#' function - a Monte-Carlo analysis. \emph{Journal of Productivity Analysis},
#' \bold{6}:247--268.
#'
#' D'Agostino, R., and E.S. Pearson. 1973. Tests for departure from normality.
#' Empirical results for the distributions of \eqn{b_2} and \eqn{\sqrt{b_1}}.
#' \emph{Biometrika}, \bold{60}:613--622.
#'
#' Greene, W. H. 2003. Simulated likelihood estimation of the normal-Gamma
#' stochastic frontier function. \emph{Journal of Productivity Analysis},
#' \bold{19}(2-3), 179--190.
#'
#' Hadri, K. 1999. Estimation of a doubly heteroscedastic stochastic frontier
#' cost function. \emph{Journal of Business & Economic Statistics},
#' \bold{17}(3), 359--363.
#'
#' Hajargasht, G. 2015. Stochastic frontiers with a Rayleigh distribution.
#' \emph{Journal of Productivity Analysis}, \bold{44}(2), 199--208.
#'
#' Huang, C. J., and Liu, J.-T. 1994. Estimation of a non-neutral stochastic
#' frontier production function. \emph{Journal of Productivity Analysis},
#' \bold{5}(2), 171--180.
#' 
#' Huang, C. J., Huang, T. H., Liu, N. H., 2014. A new approach to estimating 
#' the metafrontier production function based on a stochastic frontier 
#' framework, \emph{Journal of Productivity Analysis}, \bold{42}, 241-254.
#' 
#' Jondrow, J., C.A.K. Lovell, I.S. Materov, and P. Schmidt. 1982. On the
#' estimation of technical inefficiency in the stochastic frontier production
#' function model. \emph{Journal of Econometrics}, \bold{19}:233--238.
#'
#' Kumbhakar, S. C., Ghosh, S., and McGuckin, J. T. 1991) A generalized
#' production frontier approach for estimating determinants of inefficiency in
#' U.S. dairy farms. \emph{Journal of Business & Economic Statistics},
#' \bold{9}(3), 279--286.
#'
#' Li, Q. 1996. Estimating a stochastic production frontier when the adjusted
#' error is symmetric. \emph{Economics Letters}, \bold{52}(3), 221--228.
#'
#' Meeusen, W., and Vandenbroeck, J. 1977. Efficiency estimation from
#' Cobb-Douglas production functions with composed error. \emph{International
#' Economic Review}, \bold{18}(2), 435--445.
#'
#' Migon, H. S., and Medici, E. V. 2001. Bayesian hierarchical models for
#' stochastic production frontier. Lacea, Montevideo, Uruguay.
#'
#' Nguyen, N. B. 2010. Estimation of technical efficiency in stochastic
#' frontier analysis. PhD dissertation, Bowling Green State University, August.
#' 
#' Owusu, E. S., Bravo-Ureta, B. E., 2022. Gender and Productivity Differentials 
#' in Smallholder Groundnut Farming in Malawi: Accounting for Technology 
#' Differences, \emph{The Journal of Development Studies} \bold{58}, 989-1013.
#'
#' Papadopoulos, A. 2021. Stochastic frontier models using the generalized
#' exponential distribution. \emph{Journal of Productivity Analysis},
#' \bold{55}:15--29.
#'
#' Reifschneider, D., and Stevenson, R. 1991. Systematic departures from the
#' frontier: A framework for the analysis of firm inefficiency.
#' \emph{International Economic Review}, \bold{32}(3), 715--723.
#'
#' Stevenson, R. E. 1980. Likelihood Functions for Generalized Stochastic
#' Frontier Estimation. \emph{Journal of Econometrics}, \bold{13}(1), 57--66.
#'
#' Tsionas, E. G. 2007. Efficiency measurement with the Weibull stochastic
#' frontier. \emph{Oxford Bulletin of Economics and Statistics}, \bold{69}(5),
#' 693--706.
#'
#' Wang, K., and Ye, X. 2020. Development of alternative stochastic frontier
#' models for estimating time-space prism vertices. \emph{Transportation}.
#'
#' Wang, H.J., and Schmidt, P. 2002. One-step and two-step estimation of the
#' effects of exogenous variables on technical efficiency levels. \emph{Journal
#' of Productivity Analysis}, \bold{18}:129--144.
#'
#' Wang, J. 2012. A normal truncated skewed-Laplace model in stochastic
#' frontier analysis. Master thesis, Western Kentucky University, May.
#'
#' @keywords models optimize cross-section likelihood
#'
#' @examples
#' \dontrun{
#' # Replication of the results in Owusu and Bravo-Ureta (2022)
#' ## download data
#' 
#' data <- read.csv('https://tandf.figshare.com/ndownloader/files/32434483/filename.csv')
#' 
#' formula <- gnut ~ land + labor + seed + agrochem + usehirlab + imprseed + 
#' earlypln + latepln + soil3 + extens + socnet + farmtyp3 + farmtyp2 + dist1 + 
#' dist2
#'   
#' ## Run stochastic metafrontier following Huang et al. (2014)
#' 
#' metaHHL14 <- sfametacross(formula = formula, ghet = ~ male, data = data, 
#' modelType = 'hhl14')
#' 
#' summary(metaHHL14)
#' 
#' ### Type '3' to obtain the summary for the metafrontier
#' 
#' ## Run stochastic metafrontiers following Amsler et al. (2017)
#' ### possibility to setup a progress bar for the simulation
#' handlers(global = TRUE)
#' handlers('progress')
#' 
#' metaAOS17B <- sfametacross(formula = formula, ghet = ~ male, data = data, 
#' modelType = 'aos17b', metaSim = 5000)
#' 
#' summary(metaAOS17B)
#'  
#' }
sfametacross <- function(formula, muhet, uhet, vhet, ghet, logDepVar = TRUE, data,
  subset, weights, wscale = TRUE, S = 1L, modelType = "hhl14", udist = "hnormal",
  start = NULL, randStart = FALSE, method = "bfgs", hessianType = 1, metaSim = 5000,
  simType = "halton", Nsim = 300, prime = 2L, burn = 10, antithetics = FALSE, seed = 12345,
  itermax = 2000L, printInfo = FALSE, tol = 1e-12, gradtol = 1e-06, stepmax = 0.1,
  qac = "marquardt") {
  # metafrontier model check -------
  modelType <- tolower(modelType)
  if (!(modelType %in% c("hhl14", "bpo04a", "bpo04b", "aos17a", "aos17b", "aos17c",
    "aos17d"))) {
    stop("Unknown SFA metafrontier model: ", paste(modelType), call. = FALSE)
  }
  # meta model and udist check -------
  if (modelType %in% c("aos17b", "aos17d") && udist != "hnormal")
    stop("only 'hnormal' distribution is available when 'modelType = aos17b or aos17d'",
      call. = FALSE)
  # get correlation value from user when modelType = aos17c or aos17d -------

  # Formula manipulation -------
  if (length(Formula(formula))[2] != 1) {
    stop("argument 'formula' must have one RHS part", call. = FALSE)
  }
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mc), nomatch = 0L)
  mc <- mc[c(1L, m)]
  mc$drop.unused.levels <- TRUE
  formula <- interCheckMain(formula = formula, data = data)
  if (!missing(muhet)) {
    muhet <- clhsCheck_mu(formula = muhet)
  } else {
    muhet <- ~1
  }
  if (!missing(uhet)) {
    uhet <- clhsCheck_u(formula = uhet)
  } else {
    uhet <- ~1
  }
  if (!missing(vhet)) {
    vhet <- clhsCheck_v(formula = vhet)
  } else {
    vhet <- ~1
  }
  if (!missing(ghet)) {
    ghet <- clhsCheck_meta(formula = ghet)
  } else {
    stop("argument 'ghet' must be provided for metafrontier estimation", call. = FALSE)
  }
  formula <- formDist_sfametacross(udist = udist, formula = formula, muhet = muhet,
    uhet = uhet, vhet = vhet, ghet = ghet)
  # Generate required datasets -------
  if (missing(data)) {
    data <- environment(formula)
  }
  mc$formula <- formula
  mc$na.action <- na.pass
  mc[[1L]] <- quote(model.frame)
  mc <- eval(mc, parent.frame())
  validObs <- rowSums(is.na(mc) | is.infinite.data.frame(mc)) == 0
  Yvar <- model.response(mc, "numeric")
  Yvar <- Yvar[validObs]
  mtX <- terms(formula, data = data, rhs = 1)
  Xvar <- model.matrix(mtX, mc)
  Xvar <- Xvar[validObs, , drop = FALSE]
  nXvar <- ncol(Xvar)
  N <- nrow(Xvar)
  if (N == 0L) {
    stop("0 (non-NA) cases", call. = FALSE)
  }
  if (length(Yvar) != nrow(Xvar)) {
    stop(paste("the number of observations of the dependent variable (", length(Yvar),
      ") must be the same to the number of observations of the exogenous variables (",
      nrow(Xvar), ")", sep = ""), call. = FALSE)
  }
  if (udist %in% c("tnormal", "lognormal")) {
    mtmuH <- delete.response(terms(formula, data = data, rhs = 2))
    muHvar <- model.matrix(mtmuH, mc)
    muHvar <- muHvar[validObs, , drop = FALSE]
    nmuZUvar <- ncol(muHvar)
    mtuH <- delete.response(terms(formula, data = data, rhs = 3))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 4))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
  } else {
    mtuH <- delete.response(terms(formula, data = data, rhs = 2))
    uHvar <- model.matrix(mtuH, mc)
    uHvar <- uHvar[validObs, , drop = FALSE]
    nuZUvar <- ncol(uHvar)
    mtvH <- delete.response(terms(formula, data = data, rhs = 3))
    vHvar <- model.matrix(mtvH, mc)
    vHvar <- vHvar[validObs, , drop = FALSE]
    nvZVvar <- ncol(vHvar)
  }
  ## metafrontier variable -------
  group_var <- model.frame(ghet, data = data, na.action = na.pass)
  group_var <- group_var[validObs, , drop = FALSE]
  name_meta_var <- names(group_var)
  group_var_list <- sort(unique(group_var[, 1]))
  Ngroup <- length(group_var_list)
  if (Ngroup == 1) {
    stop("there must be at leat 2 groups for metafrontier estimation", call. = FALSE)
  }
  ## weights for groups/metafrontier likelihood -------
  wH <- as.vector(model.weights(mc))
  if (length(wscale) != 1 || !is.logical(wscale[1])) {
    stop("argument 'wscale' must be a single logical value", call. = FALSE)
  }
  if (!is.null(wH)) {
    if (!is.numeric(wH)) {
      stop("'weights' must be a numeric vector", call. = FALSE)
    } else {
      if (any(wH < 0 | is.na(wH)))
        stop("missing or negative weights not allowed", call. = FALSE)
    }
    if (wscale) {
      # scaling considers separated group
      wHvar <- numeric(N)
      for (g in group_var_list) {
        wHvar[group_var == g] <- wH[group_var == g]/sum(wH[group_var == g]) *
          length(wH[group_var == g])
      }
      # scaling for metafrontier
      wHvarm <- wH/sum(wH) * N
    }
  } else {
    wHvar <- wHvarm <- rep(1, N)
  }
  # Check other supplied options -------
  if (length(S) != 1 || !(S %in% c(-1L, 1L))) {
    stop("argument 'S' must equal either 1 or -1: 1 for production or profit frontier
   and -1 for cost frontier",
      call. = FALSE)
  }
  typeSfa <- if (S == 1L) {
    "Metafrontier for Production/Profit function, e = v - u"
  } else {
    "Metafrontier for Cost Frontier function, e = v + u"
  }
  if (length(logDepVar) != 1 || !is.logical(logDepVar[1])) {
    stop("argument 'logDepVar' must be a single logical value", call. = FALSE)
  }
  # Number of parameters -------
  nParm <- if (udist %in% c("tnormal", "lognormal")) {
    nXvar + nmuZUvar + nuZUvar + nvZVvar
  } else {
    if (udist %in% c("gamma", "weibull", "tslaplace")) {
      nXvar + nuZUvar + nvZVvar + 1
    } else {
      nXvar + nuZUvar + nvZVvar
    }
  }
  # Checking starting values when provided -------
  if (!is.null(start)) {
    if (length(start) != nParm) {
      stop("Wrong number of initial values: model has ", nParm, " parameters",
        call. = FALSE)
    }
  }
  # Set std. Error when random start is allowed -------
  sdStart <- if (randStart)
    0.01 else NULL
  # Check algorithms -------
  method <- tolower(method)
  if (!(method %in% c("ucminf", "bfgs", "bhhh", "nr", "nm", "cg", "sann", "sr1",
    "mla", "sparse", "nlminb"))) {
    stop("Unknown or non-available optimization algorithm: ", paste(method),
      call. = FALSE)
  }
  # Check hessian type
  if (length(hessianType) != 1 || !(hessianType %in% c(1L, 2L))) {
    stop("argument 'hessianType' must equal either 1 or 2", call. = FALSE)
  }
  if (!is.numeric(metaSim) || length(metaSim) != 1) {
    stop("argument 'metaSim' must be a single numeric scalar", call. = FALSE)
  }
  # SML arguments -------
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    if (!(simType %in% c("halton", "ghalton", "sobol", "rsobol", "richtmyer",
      "rrichtmyer", "uniform", "mlhs"))) {
      stop("Unknown or non-available random draws method", call. = FALSE)
    }
    if (!is.numeric(Nsim) || length(Nsim) != 1) {
      stop("argument 'Nsim' must be a single numeric scalar", call. = FALSE)
    }
    if (!is.numeric(burn) || length(burn) != 1) {
      stop("argument 'burn' must be a single numeric scalar", call. = FALSE)
    }
    if (!is_prime(prime)) {
      stop("argument 'prime' must be a single prime number", call. = FALSE)
    }
    if (length(antithetics) != 1 || !is.logical(antithetics[1])) {
      stop("argument 'antithetics' must be a single logical value", call. = FALSE)
    }
    if (antithetics && (Nsim%%2) != 0) {
      Nsim <- Nsim + 1
    }
    simDist <- if (simType == "halton") {
      "Halton"
    } else {
      if (simType == "ghalton") {
        "Generalized Halton"
      } else {
        if (simType == "sobol") {
          "Sobol"
        } else {
          if (simType == "rsobol") {
          "Randomized Sobol"
          } else {
          if (simType == "richtmyer") {
            "Richtmyer"
          } else {
            if (simType == "rrichtmyer") {
            "Randomized Richtmyer"
            } else {
            if (simType == "uniform") {
              "Uniform"
            } else {
              if (simType == "mlhs") {
              "Modified Latin Hypercube Sampling"
              }
            }
            }
          }
          }
        }
      }
    }
  }
  # Other optimization options -------
  if (!is.numeric(itermax) || length(itermax) != 1) {
    stop("argument 'itermax' must be a single numeric scalar", call. = FALSE)
  }
  if (itermax != round(itermax)) {
    stop("argument 'itermax' must be an integer", call. = FALSE)
  }
  if (itermax <= 0) {
    stop("argument 'itermax' must be positive", call. = FALSE)
  }
  itermax <- as.integer(itermax)
  if (length(printInfo) != 1 || !is.logical(printInfo[1])) {
    stop("argument 'printInfo' must be a single logical value", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1) {
    stop("argument 'tol' must be numeric", call. = FALSE)
  }
  if (tol < 0) {
    stop("argument 'tol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(gradtol) || length(gradtol) != 1) {
    stop("argument 'gradtol' must be numeric", call. = FALSE)
  }
  if (gradtol < 0) {
    stop("argument 'gradtol' must be non-negative", call. = FALSE)
  }
  if (!is.numeric(stepmax) || length(stepmax) != 1) {
    stop("argument 'stepmax' must be numeric", call. = FALSE)
  }
  if (stepmax < 0) {
    stop("argument 'stepmax' must be non-negative", call. = FALSE)
  }
  if (!(qac %in% c("marquardt", "stephalving"))) {
    stop("argument 'qac' must be either 'marquardt' or 'stephalving'", call. = FALSE)
  }
  if (inherits(data, "plm.dim")) {
    dataTablem <- data[validObs, names(index(data))]
  } else {
    dataTablem <- data.frame(IdObs = c(1:sum(validObs)))
  }
  dataTablem <- cbind(dataTablem, data[, all.vars(terms(formula))], weights = wHvar)
  # Looping over each group -------
  startVal <- list()
  olsParam <- list()
  olsStder <- list()
  olsSigmasq <- list()
  olsLoglik <- list()
  olsSkew <- list()
  olsM3Okay <- list()
  CoelliM3Test <- list()
  AgostinoTest <- list()
  dataTable <- list()
  optType <- list()
  nIter <- list()
  optStatus <- list()
  startLoglik <- list()
  mlLoglik <- list()
  mlParam <- list()
  gradient <- list()
  gradL_OBS <- list()
  gradientNorm <- list()
  invHessian <- list()
  conditionNums <- list()
  FiMat <- list()
  mlDate <- list()
  for (g in group_var_list) {
    ## sub-sample the different variables
    Yvarg <- Yvar[group_var[, 1] == g]
    Xvarg <- Xvar[group_var[, 1] == g, , drop = FALSE]
    Ng <- nrow(Xvarg)
    if (udist %in% c("tnormal", "lognormal")) {
      muHvarg <- muHvar[group_var[, 1] == g, , drop = FALSE]
      uHvarg <- uHvar[group_var[, 1] == g, , drop = FALSE]
      vHvarg <- vHvar[group_var[, 1] == g, , drop = FALSE]
    } else {
      uHvarg <- uHvar[group_var[, 1] == g, , drop = FALSE]
      vHvarg <- vHvar[group_var[, 1] == g, , drop = FALSE]
    }
    wHvarg <- wHvar[group_var[, 1] == g]
    # Step 1: OLS -------
    olsRes <- if (colnames(Xvarg)[1] == "(Intercept)") {
      if (dim(Xvarg)[2] == 1) {
        lm(Yvarg ~ 1)
      } else {
        lm(Yvarg ~ ., data = as.data.frame(Xvarg[, -1]), weights = wHvarg)
      }
    } else {
      lm(Yvarg ~ -1 + ., data = as.data.frame(Xvarg), weights = wHvarg)
    }
    if (any(is.na(olsRes$coefficients))) {
      stop("at least one of the OLS coefficients is NA: ", paste(colnames(Xvarg)[is.na(olsRes$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
        call. = FALSE)
    }
    names(olsRes$coefficients) <- colnames(Xvar)
    olsParamg <- c(olsRes$coefficients)
    olsSkewg <- skewness(residuals(olsRes))
    olsSigmasqg <- summary(olsRes)$sigma^2
    olsStderg <- sqrt(diag(vcov(olsRes)))
    olsLoglikg <- logLik(olsRes)[1]
    olsM3Okayg <- if (S * olsSkewg < 0) {
      "Residuals have the expected skeweness"
    } else {
      "Residuals do not have the expected skeweness"
    }
    if (S * olsSkewg > 0) {
      warning(paste0("The residuals of the OLS for group = ", g, " are"), if (S ==
        1) {
        " right"
      } else {
        " left"
      }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
        call. = FALSE)
    }
    m2 <- mean((residuals(olsRes) - mean(residuals(olsRes)))^2)
    m3 <- mean((residuals(olsRes) - mean(residuals(olsRes)))^3)
    CoelliM3Testg <- c(z = m3/sqrt(6 * m2^3/N), p.value = 2 * pnorm(-abs(m3/sqrt(6 *
      m2^3/N))))
    AgostinoTestg <- dagoTest(residuals(olsRes))
    class(AgostinoTestg) <- "dagoTest"
    if (udist %in% c("gamma", "lognormal", "weibull")) {
      # Initialization of draws for each group
      cat("Initialization of", Nsim, simDist, "draws per observation ...\n")
      FiMatg <- drawMatUniDim(N = Ng, Nsim = Nsim, simType = simType, prime = prime,
        burn = burn + 1, antithetics = antithetics, seed = seed)
    }
    dataTableg <- dataTablem[group_var == g, ]
    dataTableg <- cbind(dataTableg, olsResiduals = residuals(olsRes), olsFitted = fitted(olsRes))
    FunArgs <- if (udist == "tnormal") {
      list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParamg,
        dataTable = dataTableg, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
        nvZVvar = nvZVvar, muHvar = muHvarg, uHvar = uHvarg, vHvar = vHvarg,
        Yvar = Yvarg, Xvar = Xvarg, S = S, wHvar = wHvarg, method = method,
        printInfo = printInfo, itermax = itermax, stepmax = stepmax, tol = tol,
        gradtol = gradtol, hessianType = hessianType, qac = qac)
    } else {
      if (udist == "lognormal") {
        list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParamg,
          dataTable = dataTableg, nXvar = nXvar, nmuZUvar = nmuZUvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, muHvar = muHvarg, uHvar = uHvarg, vHvar = vHvarg,
          Yvar = Yvarg, Xvar = Xvarg, S = S, wHvar = wHvarg, N = Ng, FiMat = FiMatg,
          method = method, printInfo = printInfo, itermax = itermax, stepmax = stepmax,
          tol = tol, gradtol = gradtol, hessianType = hessianType, qac = qac)
      } else {
        if (udist %in% c("gamma", "weibull")) {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParamg,
          dataTable = dataTableg, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvarg, vHvar = vHvarg, Yvar = Yvarg, Xvar = Xvarg, S = S,
          wHvar = wHvarg, N = Ng, FiMat = FiMatg, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
        } else {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParamg,
          dataTable = dataTableg, nXvar = nXvar, nuZUvar = nuZUvar, nvZVvar = nvZVvar,
          uHvar = uHvarg, vHvar = vHvarg, Yvar = Yvarg, Xvar = Xvarg, S = S,
          wHvar = wHvarg, method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
          qac = qac)
        }
      }
    }
    ## MLE run -------
    cat("Stochastic Frontier for group = ", g, "\n")
    mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt, FunArgs),
      exponential = do.call(exponormAlgOpt, FunArgs), tnormal = do.call(truncnormAlgOpt,
        FunArgs), rayleigh = do.call(raynormAlgOpt, FunArgs), gamma = do.call(gammanormAlgOpt,
        FunArgs), uniform = do.call(uninormAlgOpt, FunArgs), lognormal = do.call(lognormAlgOpt,
        FunArgs), weibull = do.call(weibullnormAlgOpt, FunArgs), genexponential = do.call(genexponormAlgOpt,
        FunArgs), tslaplace = do.call(tslnormAlgOpt, FunArgs)), error = function(e) e)
    if (inherits(mleList, "error")) {
      stop("The current error occurs during optimization:\n", mleList$message,
        call. = FALSE)
    }
    # Inverse Hessian + other -------
    mleList$invHessian <- vcovObj(mleObj = mleList$mleObj, hessianType = hessianType,
      method = method, nParm = nParm)
    mleList <- c(mleList, if (method == "ucminf") {
      list(type = "ucminf max.", nIter = unname(mleList$mleObj$info["neval"]),
        status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
        list(type = substr(mleList$mleObj$type, 1, 27), nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
          gradient = mleList$mleObj$gradient)
      } else {
        if (method == "sr1") {
          list(type = "SR1 max.", nIter = mleList$mleObj$iterations, status = mleList$mleObj$status,
          mleLoglik = -mleList$mleObj$fval, gradient = mleList$mleObj$gradient)
        } else {
          if (method == "mla") {
          list(type = "Lev. Marquardt max.", nIter = mleList$mleObj$ni,
            status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached", `4` = "algorithm encountered a problem in the function computation"),
            mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad)
          } else {
          if (method == "sparse") {
            list(type = "Sparse Hessian max.", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
            gradient = mleList$mleObj$gradient)
          } else {
            if (method == "nlminb") {
            list(type = "nlminb max.", nIter = mleList$mleObj$iterations,
              status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
              gradient = mleList$mleObj$gradient)
            }
          }
          }
        }
      }
    })
    # quick renaming -------
    if (udist %in% c("tnormal", "lognormal")) {
      names(mleList$startVal) <- fName_mu_sfacross(Xvar = Xvar, udist = udist,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, scaling = FALSE)
    } else {
      names(mleList$startVal) <- fName_uv_sfacross(Xvar = Xvar, udist = udist,
        uHvar = uHvar, vHvar = vHvar)
    }
    names(mleList$mlParam) <- names(mleList$startVal)
    rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
    names(mleList$gradient) <- names(mleList$mlParam)
    colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
    mlDateg <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
    # Return objects by loop -------
    dataTableg$mlResiduals <- Yvarg - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
      t(Xvarg)))
    dataTableg$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
      t(Xvarg)))
    dataTableg$logL_OBS <- mleList$mleObj$logL_OBS
    dataTable[[which(group_var_list == g)]] <- dataTableg
    startVal[[which(group_var_list == g)]] <- mleList$startVal
    olsParam[[which(group_var_list == g)]] <- olsParamg
    olsStder[[which(group_var_list == g)]] <- olsStderg
    olsSigmasq[[which(group_var_list == g)]] <- olsSigmasqg
    olsLoglik[[which(group_var_list == g)]] <- olsLoglikg
    olsSkew[[which(group_var_list == g)]] <- olsSkewg
    olsM3Okay[[which(group_var_list == g)]] <- olsM3Okayg
    CoelliM3Test[[which(group_var_list == g)]] <- CoelliM3Testg
    AgostinoTest[[which(group_var_list == g)]] <- AgostinoTestg
    optType[[which(group_var_list == g)]] <- mleList$type
    nIter[[which(group_var_list == g)]] <- mleList$nIter
    optStatus[[which(group_var_list == g)]] <- mleList$status
    startLoglik[[which(group_var_list == g)]] <- mleList$startLoglik
    mlLoglik[[which(group_var_list == g)]] <- mleList$mleLoglik
    mlParam[[which(group_var_list == g)]] <- mleList$mlParam
    gradient[[which(group_var_list == g)]] <- mleList$gradient
    gradL_OBS[[which(group_var_list == g)]] <- mleList$mleObj$gradL_OBS
    gradientNorm[[which(group_var_list == g)]] <- sqrt(sum(mleList$gradient^2))
    invHessian[[which(group_var_list == g)]] <- mleList$invHessian
    conditionNums[[which(group_var_list == g)]] <- condiNum(mleObj = mleList$mleObj,
      method = method, nParm = nParm)
    if (udist %in% c("gamma", "lognormal", "weibull")) {
      FiMat[[which(group_var_list == g)]] <- FiMatg
    }
    mlDate[[which(group_var_list == g)]] <- mlDateg
    rm(mleList)
  }
  # metafrontier parameters -------
  Yvarm <- numeric(N)
  for (g in group_var_list) {
    Yvarm[group_var == g] <- dataTable[[which(group_var_list == g)]]$mlFitted
  }
  ## hhl14 model -------
  if (modelType == "hhl14") {
    ## Step 1: OLS -------
    olsRes <- if (colnames(Xvar)[1] == "(Intercept)") {
      if (dim(Xvar)[2] == 1) {
        lm(Yvarm ~ 1)
      } else {
        lm(Yvarm ~ ., data = as.data.frame(Xvar[, -1]), weights = wHvarm)
      }
    } else {
      lm(Yvarm ~ -1 + ., data = as.data.frame(Xvar), weights = wHvarm)
    }
    if (any(is.na(olsRes$coefficients))) {
      stop("at least one of the OLS coefficients is NA: ", paste(colnames(Xvar)[is.na(olsRes$coefficients)],
        collapse = ", "), "This may be due to a singular matrix
   due to potential perfect multicollinearity",
        call. = FALSE)
    }
    olsParam[[Ngroup + 1]] <- c(olsRes$coefficients)
    olsSkew[[Ngroup + 1]] <- skewness(residuals(olsRes))
    olsSigmasq[[Ngroup + 1]] <- summary(olsRes)$sigma^2
    olsStder[[Ngroup + 1]] <- sqrt(diag(vcov(olsRes)))
    olsLoglik[[Ngroup + 1]] <- logLik(olsRes)[1]
    olsM3Okay[[Ngroup + 1]] <- if (S * olsSkew[[Ngroup + 1]] < 0) {
      "Residuals have the expected skeweness"
    } else {
      "Residuals do not have the expected skeweness"
    }
    if (S * olsSkew[[Ngroup + 1]] > 0) {
      warning(paste0("The residuals of the OLS for metafrontier ", " are"),
        if (S == 1) {
          " right"
        } else {
          " left"
        }, "-skewed. This may indicate the absence of inefficiency or
  model misspecification or sample 'bad luck'",
        call. = FALSE)
    }
    m2 <- mean((residuals(olsRes) - mean(residuals(olsRes)))^2)
    m3 <- mean((residuals(olsRes) - mean(residuals(olsRes)))^3)
    CoelliM3Test[[Ngroup + 1]] <- c(z = m3/sqrt(6 * m2^3/N), p.value = 2 * pnorm(-abs(m3/sqrt(6 *
      m2^3/N))))
    AgostinoTestm <- dagoTest(residuals(olsRes))
    class(AgostinoTestm) <- "dagoTest"
    AgostinoTest[[Ngroup + 1]] <- AgostinoTestm
    if (udist %in% c("gamma", "lognormal", "weibull")) {
      # Initialization of draws for metafrontier
      cat("Initialization of", Nsim, simDist, "draws per observation ...\n")
      FiMat[[Ngroup + 1]] <- drawMatUniDim(N = N, Nsim = Nsim, simType = simType,
        prime = prime, burn = burn + 1, antithetics = antithetics, seed = seed)
    }
    dataTable[[Ngroup + 1]] <- cbind(dataTablem, olsResiduals = residuals(olsRes),
      olsFitted = fitted(olsRes))
    FunArgs <- if (udist == "tnormal") {
      list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam[[Ngroup +
        1]], dataTable = dataTable[[Ngroup + 1]], nXvar = nXvar, nmuZUvar = nmuZUvar,
        nuZUvar = nuZUvar, nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
        vHvar = vHvar, Yvar = Yvarm, Xvar = Xvar, S = S, wHvar = wHvarm,
        method = method, printInfo = printInfo, itermax = itermax, stepmax = stepmax,
        tol = tol, gradtol = gradtol, hessianType = hessianType, qac = qac)
    } else {
      if (udist == "lognormal") {
        list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam[[Ngroup +
          1]], dataTable = dataTable[[Ngroup + 1]], nXvar = nXvar, nmuZUvar = nmuZUvar,
          nuZUvar = nuZUvar, nvZVvar = nvZVvar, muHvar = muHvar, uHvar = uHvar,
          vHvar = vHvar, Yvar = Yvarm, Xvar = Xvar, S = S, wHvar = wHvarm,
          N = N, FiMat = FiMat[[Ngroup + 1]], method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
      } else {
        if (udist %in% c("gamma", "weibull")) {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam[[Ngroup +
          1]], dataTable = dataTable[[Ngroup + 1]], nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvarm,
          Xvar = Xvar, S = S, wHvar = wHvarm, N = N, FiMat = FiMat[[Ngroup +
            1]], method = method, printInfo = printInfo, itermax = itermax,
          stepmax = stepmax, tol = tol, gradtol = gradtol, hessianType = hessianType,
          qac = qac)
        } else {
          list(start = start, randStart = randStart, sdStart = sdStart, olsParam = olsParam[[Ngroup +
          1]], dataTable = dataTable[[Ngroup + 1]], nXvar = nXvar, nuZUvar = nuZUvar,
          nvZVvar = nvZVvar, uHvar = uHvar, vHvar = vHvar, Yvar = Yvarm,
          Xvar = Xvar, S = S, wHvar = wHvarm, method = method, printInfo = printInfo,
          itermax = itermax, stepmax = stepmax, tol = tol, gradtol = gradtol,
          hessianType = hessianType, qac = qac)
        }
      }
    }
    ## QMLE run -------
    cat("Stochastic Metafrontier: model hhl14")
    mleList <- tryCatch(switch(udist, hnormal = do.call(halfnormAlgOpt, FunArgs),
      exponential = do.call(exponormAlgOpt, FunArgs), tnormal = do.call(truncnormAlgOpt,
        FunArgs), rayleigh = do.call(raynormAlgOpt, FunArgs), gamma = do.call(gammanormAlgOpt,
        FunArgs), uniform = do.call(uninormAlgOpt, FunArgs), lognormal = do.call(lognormAlgOpt,
        FunArgs), weibull = do.call(weibullnormAlgOpt, FunArgs), genexponential = do.call(genexponormAlgOpt,
        FunArgs), tslaplace = do.call(tslnormAlgOpt, FunArgs)), error = function(e) e)
    if (inherits(mleList, "error")) {
      stop("The current error occurs during optimization:\n", mleList$message,
        call. = FALSE)
    }
    ## sandwich vcov Hessian + other -------
    invhess1 <- vcovObj(mleObj = mleList$mleObj, hessianType = 1, method = method,
      nParm = nParm)
    grgr <- crossprod(mleList$mleObj$gradL_OBS)
    mleList$invHessian <- invhess1 %*% grgr %*% invhess1
    mleList <- c(mleList, if (method == "ucminf") {
      list(type = "ucminf max.", nIter = unname(mleList$mleObj$info["neval"]),
        status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$value,
        gradient = mleList$mleObj$gradient)
    } else {
      if (method %in% c("bfgs", "bhhh", "nr", "nm", "cg", "sann")) {
        list(type = substr(mleList$mleObj$type, 1, 27), nIter = mleList$mleObj$iterations,
          status = mleList$mleObj$message, mleLoglik = mleList$mleObj$maximum,
          gradient = mleList$mleObj$gradient)
      } else {
        if (method == "sr1") {
          list(type = "SR1 max.", nIter = mleList$mleObj$iterations, status = mleList$mleObj$status,
          mleLoglik = -mleList$mleObj$fval, gradient = mleList$mleObj$gradient)
        } else {
          if (method == "mla") {
          list(type = "Lev. Marquardt max.", nIter = mleList$mleObj$ni,
            status = switch(mleList$mleObj$istop, `1` = "convergence criteria were satisfied",
            `2` = "maximum number of iterations was reached", `4` = "algorithm encountered a problem in the function computation"),
            mleLoglik = -mleList$mleObj$fn.value, gradient = mleList$mleObj$grad)
          } else {
          if (method == "sparse") {
            list(type = "Sparse Hessian max.", nIter = mleList$mleObj$iterations,
            status = mleList$mleObj$status, mleLoglik = -mleList$mleObj$fval,
            gradient = mleList$mleObj$gradient)
          } else {
            if (method == "nlminb") {
            list(type = "nlminb max.", nIter = mleList$mleObj$iterations,
              status = mleList$mleObj$message, mleLoglik = -mleList$mleObj$objective,
              gradient = mleList$mleObj$gradient)
            }
          }
          }
        }
      }
    })
    ## quick renaming -------
    if (udist %in% c("tnormal", "lognormal")) {
      names(mleList$startVal) <- fName_mu_sfacross(Xvar = Xvar, udist = udist,
        muHvar = muHvar, uHvar = uHvar, vHvar = vHvar, scaling = FALSE)
    } else {
      names(mleList$startVal) <- fName_uv_sfacross(Xvar = Xvar, udist = udist,
        uHvar = uHvar, vHvar = vHvar)
    }
    rownames(mleList$invHessian) <- colnames(mleList$invHessian) <- names(mleList$mlParam)
    names(mleList$mlParam) <- names(mleList$startVal)
    names(mleList$gradient) <- names(mleList$mlParam)
    colnames(mleList$mleObj$gradL_OBS) <- names(mleList$mlParam)
    ## return meta objects -------
    mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
    startVal[[Ngroup + 1]] <- mleList$startVal
    optType[[Ngroup + 1]] <- mleList$type
    nIter[[Ngroup + 1]] <- mleList$nIter
    optStatus[[Ngroup + 1]] <- mleList$status
    startLoglik[[Ngroup + 1]] <- mleList$startLoglik
    mlLoglik[[Ngroup + 1]] <- mleList$mleLoglik
    mlParam[[Ngroup + 1]] <- mleList$mlParam
    gradient[[Ngroup + 1]] <- mleList$gradient
    gradL_OBS[[Ngroup + 1]] <- mleList$mleObj$gradL_OBS
    gradientNorm[[Ngroup + 1]] <- sqrt(sum(mleList$gradient^2))
    invHessian[[Ngroup + 1]] <- mleList$invHessian
    conditionNums[[Ngroup + 1]] <- condiNum(mleObj = mleList$mleObj, method = method,
      nParm = nParm)
    dataTable[[Ngroup + 1]]$mlResiduals <- Yvarm - as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
      t(Xvar)))
    dataTable[[Ngroup + 1]]$mlFitted <- as.numeric(crossprod(matrix(mleList$mlParam[1:nXvar]),
      t(Xvar)))
    dataTable[[Ngroup + 1]]$logL_OBS <- mleList$mleObj$logL_OBS
    rm(mleList)
  } else {
    ## bpo04a model -------
    if (modelType == "bpo04a") {
      cat("Stochastic Metafrontier: model bpo04a \n")
      ### LP -------
      metaLPres <- metaLPfun(N = N, nXvar = nXvar, Xvar = S * Xvar, Yvarm = S *
        Yvarm, varnames = names(mlParam[[1]])[1:nXvar])
      #### simulation for se ----------
      cat("Model bpo04a: simulation for standard errors \n")
      randFitted <- matrix(nrow = N, ncol = metaSim)
      set.seed(seed)
      for (g in group_var_list) {
        randDrawMat <- mnorm::rmnorm(metaSim, mean = mlParam[[which(group_var_list ==
          g)]][1:nXvar], sigma = invHessian[[which(group_var_list == g)]][1:nXvar,
          1:nXvar])
        Xvarg <- Xvar[group_var[, 1] == g, , drop = FALSE]
        randFitted[group_var == g, ] <- tcrossprod(Xvarg, randDrawMat)
      }
      # progressr::handlers('progress')
      p <- progressr::progressor(metaSim)
      metaSimVec <- seq_len(metaSim)
      lpSimRes <- list()
      for (sim in metaSimVec) {
        lpSimRes[[sim]] <- metaLPfun(N = N, nXvar = nXvar, Xvar = S * Xvar,
          Yvarm = S * randFitted[, sim], varnames = names(mlParam[[1]])[1:nXvar])$parRes
        p(sprintf("Simulation %g", metaSimVec[sim]))
      }
      lpSimRes <- do.call(rbind, lpSimRes)
      rm(randFitted)
      ### return meta objects -------
      mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
      mlParam[[Ngroup + 1]] <- c(metaLPres$parRes, rep(NA, length(mlParam[[Ngroup]]) -
        nXvar))
      startVal[[Ngroup + 1]] <- NA
      optType[[Ngroup + 1]] <- "Linear Programming: lpsolve"
      nIter[[Ngroup + 1]] <- metaLPres$iter
      optStatus[[Ngroup + 1]] <- lpStatus(metaLPres$statusCode)
      startLoglik[[Ngroup + 1]] <- NA
      mlLoglik[[Ngroup + 1]] <- NA
      gradient[[Ngroup + 1]] <- NA
      gradL_OBS[[Ngroup + 1]] <- NA
      gradientNorm[[Ngroup + 1]] <- NA
      vcovSim <- var(lpSimRes)
      colnames(vcovSim) <- rownames(vcovSim) <- names(mlParam[[1]])[1:nXvar]
      invHessian[[Ngroup + 1]] <- vcovSim
      conditionNums[[Ngroup + 1]] <- NA
      mlFitted <- as.numeric(crossprod(matrix(metaLPres$parRes), t(Xvar)))
      mlResiduals <- Yvarm - mlFitted
      dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
        mlFitted = mlFitted)
      dataTable[[Ngroup + 1]]$logL_OBS <- NA
    } else {
      ## bpo04b model -------
      if (modelType == "bpo04b") {
        cat("Stochastic Metafrontier: model bpo04b", "\n")
        ### QP ------
        parRes <- lsei::lsi(a = Xvar, b = Yvarm, e = S * Xvar, f = S * Yvarm)
        #### simulation for se ----------
        cat("Model bpo04b: simulation for standard errors \n")
        randFitted <- matrix(nrow = N, ncol = metaSim)
        set.seed(seed)
        for (g in group_var_list) {
          randDrawMat <- mnorm::rmnorm(metaSim, mean = mlParam[[which(group_var_list ==
          g)]][1:nXvar], sigma = invHessian[[which(group_var_list == g)]][1:nXvar,
          1:nXvar])
          Xvarg <- Xvar[group_var[, 1] == g, , drop = FALSE]
          randFitted[group_var == g, ] <- tcrossprod(Xvarg, randDrawMat)
        }
        # progressr::handlers('progress')
        p <- progressr::progressor(metaSim)
        metaSimVec <- seq_len(metaSim)
        qpSimRes <- list()
        for (sim in metaSimVec) {
          qpSimRes[[sim]] <- lsei::lsi(a = Xvar, b = randFitted[, sim], e = Xvar,
          f = randFitted[, sim])
          p(sprintf("Simulation %g", metaSimVec[sim]))
        }
        qpSimRes <- do.call(rbind, qpSimRes)
        colnames(qpSimRes) <- names(mlParam[[1]])[1:nXvar]
        rm(randFitted)
        ### return meta objects -------
        mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
        names(parRes) <- names(mlParam[[1]])[1:nXvar]
        mlParam[[Ngroup + 1]] <- c(parRes, rep(NA, length(mlParam[[Ngroup]]) -
          nXvar))
        startVal[[Ngroup + 1]] <- NA
        optType[[Ngroup + 1]] <- "Quadratic Programming: lsei"
        nIter[[Ngroup + 1]] <- NA
        optStatus[[Ngroup + 1]] <- NA
        startLoglik[[Ngroup + 1]] <- NA
        mlLoglik[[Ngroup + 1]] <- NA
        gradient[[Ngroup + 1]] <- NA
        gradL_OBS[[Ngroup + 1]] <- NA
        gradientNorm[[Ngroup + 1]] <- NA
        vcovSim <- var(qpSimRes)
        colnames(vcovSim) <- rownames(vcovSim) <- names(mlParam[[1]])[1:nXvar]
        invHessian[[Ngroup + 1]] <- vcovSim
        conditionNums[[Ngroup + 1]] <- NA
        mlFitted <- as.numeric(crossprod(matrix(parRes), t(Xvar)))
        mlResiduals <- Yvarm - mlFitted
        dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
          mlFitted = mlFitted)
        dataTable[[Ngroup + 1]]$logL_OBS <- NA
      } else {
        ## aos17a model -------
        if (modelType == "aos17a") {
          cat("Stochastic Metafrontier: model aos17a", "\n")
          ### simulation -------
          betaMat <- do.call(cbind, mlParam)
          fittedMat <- Xvar %*% betaMat[1:nXvar, , drop = FALSE]
          epsilonMat <- sweep(-fittedMat, MARGIN = 1, STATS = Yvar, FUN = "+")
          if (udist %in% c("tnormal", "lognormal")) {
          muMat <- muHvar %*% betaMat[(nXvar + 1):(nXvar + nmuZUvar), ,
            drop = FALSE]
          WuMat <- uHvar %*% betaMat[(nXvar + nmuZUvar + 1):(nXvar + nmuZUvar +
            nuZUvar), , drop = FALSE]
          WvMat <- vHvar %*% betaMat[(nXvar + nmuZUvar + nuZUvar + 1):(nXvar +
            nmuZUvar + nuZUvar + nvZVvar), , drop = FALSE]
          } else {
          WuMat <- uHvar %*% betaMat[(nXvar + 1):(nXvar + nuZUvar), , drop = FALSE]
          WvMat <- vHvar %*% betaMat[(nXvar + nuZUvar + 1):(nXvar + nuZUvar +
            nvZVvar), , drop = FALSE]
          }
          # progressr::handlers('progress')
          set.seed(seed)
          p <- progressr::progressor(N)
          MetaList <- list()
          for (i in seq_len(N)) {
          draws <- mnorm::rmnorm(n = metaSim, mean = fittedMat[i, ], sigma = diag(exp(WvMat[i,
            ])))
          if (S == 1) {
            metaFrontier <- apply(draws, 1, max)
          } else {
            metaFrontier <- apply(draws, 1, min)
          }
          groupFrontier <- draws[, which(group_var_list == group_var[i,
            ])]
          MetaList[[i]] <- list(groupFrontier, metaFrontier)
          p(sprintf("Simulation for observation %g", i))
          }
          groupFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
          FUN = function(z) MetaList[[z]][[1]]))
          metaFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
          FUN = function(z) MetaList[[z]][[2]]))
          MdMat <- if (S == 1) {
          metaFrontierMat - groupFrontierMat
          } else {
          groupFrontierMat - metaFrontierMat
          }
          colnames(groupFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
          colnames(metaFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
          colnames(MdMat) <- paste0("Md_Sim#", 1:metaSim)
          Mud <- numeric(N)
          for (j in seq_len(N)) {
          Mud[j] <- if (udist == "hnormal") {
            sqrt(2/pi) * exp(WuMat[j, which(group_var_list == group_var[j,
            ])]/2)
          } else {
            if (udist == "tnormal") {
            mean(muMat[, which(group_var_list == group_var[j, ])]) +
              exp(WuMat[j, which(group_var_list == group_var[j, ])]/2) *
              dnorm(mean(muMat[, which(group_var_list == group_var[j,
                ])])/exp(WuMat[j, which(group_var_list == group_var[j,
                ])]/2))/pnorm(mean(muMat[, which(group_var_list ==
              group_var[j, ])])/exp(WuMat[j, which(group_var_list ==
              group_var[j, ])]/2))
            } else {
            if (udist == "exponential") {
              exp(WuMat[j, which(group_var_list == group_var[j, ])]/2)
            } else {
              if (udist == "rayleigh") {
              exp(WuMat[j, which(group_var_list == group_var[j, ])]/2) *
                sqrt(pi/2)
              } else {
              if (udist == "gamma") {
                betaMat[nXvar + nuZUvar + nvZVvar + 1, which(group_var_list ==
                group_var[j, ])] * exp(WuMat[j, which(group_var_list ==
                group_var[j, ])]/2)
              } else {
                if (udist == "lognormal") {
                exp(mean(muMat[, which(group_var_list == group_var[j,
                  ])]) + WuMat[j, which(group_var_list == group_var[j,
                  ])]/2)
                } else {
                if (udist == "uniform") {
                  sqrt(12 * exp(WuMat[j, which(group_var_list ==
                  group_var[j, ])]/2))/2
                } else {
                  if (udist == "genexponential") {
                  3/2 * exp(WuMat[j, which(group_var_list == group_var[j,
                    ])]/2)
                  } else {
                  if (udist == "tslaplace") {
                    lambda <- betaMat[nXvar + nuZUvar + nvZVvar +
                    1, which(group_var_list == group_var[j, ])]
                    exp(WuMat[j, which(group_var_list == group_var[j,
                    ])]/2) * (1 + 4 * lambda + 2 * lambda^2)/((1 +
                    lambda) * (1 + 2 * lambda))
                  } else {
                    if (udist == "weibull") {
                    exp(WuMat[j, which(group_var_list == group_var[j,
                      ])]/2) * gamma(1 + 1/betaMat[nXvar + nuZUvar +
                      nvZVvar + 1, which(group_var_list == group_var[j,
                      ])])
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
          ### return meta objects -------
          mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
          mlParam[[Ngroup + 1]] <- NA
          startVal[[Ngroup + 1]] <- NA
          optType[[Ngroup + 1]] <- NA
          nIter[[Ngroup + 1]] <- NA
          optStatus[[Ngroup + 1]] <- NA
          startLoglik[[Ngroup + 1]] <- NA
          mlLoglik[[Ngroup + 1]] <- NA
          gradient[[Ngroup + 1]] <- NA
          gradL_OBS[[Ngroup + 1]] <- NA
          gradientNorm[[Ngroup + 1]] <- NA
          invHessian[[Ngroup + 1]] <- NA
          conditionNums[[Ngroup + 1]] <- NA
          mlFitted <- apply(metaFrontierMat, 1, mean)
          mlResiduals <- Yvarm - mlFitted
          dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
          mlFitted = mlFitted)
          dataTable[[Ngroup + 1]]$logL_OBS <- NA
        } else {
          ## aos17b model -------
          if (modelType == "aos17b") {
          cat("Stochastic Metafrontier: model aos17b", "\n")
          ### simulation -------
          betaMat <- do.call(cbind, mlParam)
          fittedMat <- Xvar %*% betaMat[1:nXvar, , drop = FALSE]
          WuMat <- uHvar %*% betaMat[(nXvar + 1):(nXvar + nuZUvar), , drop = FALSE]
          WvMat <- vHvar %*% betaMat[(nXvar + nuZUvar + 1):(nXvar + nuZUvar +
            nvZVvar), , drop = FALSE]
          epsilonMat <- sweep(-fittedMat, MARGIN = 1, STATS = Yvar, FUN = "+")
          # progressr::handlers('progress')
          set.seed(seed)
          p <- progressr::progressor(N)
          MetaList <- list()
          for (i in seq_len(N)) {
            sigmavsq <- exp(WvMat[i, ])
            sigmausq <- exp(WuMat[i, ])
            sigmasq <- sigmausq + sigmavsq
            if (S == 1) {
            vdraws <- tmvtnorm::rtmvnorm(n = metaSim, mean = sigmavsq *
              epsilonMat[i, ]/sigmasq, sigma = diag(sigmausq * sigmavsq/sigmasq),
              lower = epsilonMat[i, ], algorithm = "gibbs")
            } else {
            vdraws <- tmvtnorm::rtmvnorm(n = metaSim, mean = sigmavsq *
              epsilonMat[i, ]/sigmasq, sigma = diag(sigmausq * sigmavsq/sigmasq),
              upper = epsilonMat[i, ], algorithm = "gibbs")
            }
            draws <- mnorm::rmnorm(metaSim, mean = fittedMat[i, ], sigma = diag(sigmavsq))
            draws[, which(group_var_list == group_var[i, ])] <- fittedMat[i,
            which(group_var_list == group_var[i, ])] + vdraws[, which(group_var_list ==
            group_var[i, ])]
            groupFrontier <- draws[, which(group_var_list == group_var[i,
            ])]
            if (S == 1) {
            metaFrontier <- apply(draws, 1, max)
            } else {
            metaFrontier <- apply(draws, 1, min)
            }
            MetaList[[i]] <- list(groupFrontier, metaFrontier)
            p(sprintf("Simulation for observation %g", i))
          }
          groupFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
            FUN = function(z) MetaList[[z]][[1]]))
          metaFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
            FUN = function(z) MetaList[[z]][[2]]))
          MdMat <- if (S == 1) {
            metaFrontierMat - groupFrontierMat
          } else {
            groupFrontierMat - metaFrontierMat
          }
          colnames(groupFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
          colnames(metaFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
          colnames(MdMat) <- paste0("Md_Sim#", 1:metaSim)
          Mud <- numeric(N)
          for (obs in seq_len(N)) {
            Mud[obs] <- sqrt(2/pi) * exp(WuMat[obs, which(group_var_list ==
            group_var[obs, ])]/2)
          }
          ### return meta objects -------
          mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
          mlParam[[Ngroup + 1]] <- NA
          startVal[[Ngroup + 1]] <- NA
          optType[[Ngroup + 1]] <- NA
          nIter[[Ngroup + 1]] <- NA
          optStatus[[Ngroup + 1]] <- NA
          startLoglik[[Ngroup + 1]] <- NA
          mlLoglik[[Ngroup + 1]] <- NA
          gradient[[Ngroup + 1]] <- NA
          gradL_OBS[[Ngroup + 1]] <- NA
          gradientNorm[[Ngroup + 1]] <- NA
          invHessian[[Ngroup + 1]] <- NA
          conditionNums[[Ngroup + 1]] <- NA
          mlFitted <- apply(metaFrontierMat, 1, mean)
          mlResiduals <- Yvarm - mlFitted
          dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
            mlFitted = mlFitted)
          dataTable[[Ngroup + 1]]$logL_OBS <- NA
          } else {
          ## aos17c model -------
          if (modelType == "aos17c") {
            cat("Stochastic Metafrontier: model aos17c", "\n")
            corval <- inputNumber("Enter a value for the equicorrelation: ")
            if (corval < 0.01 || corval > 0.99)
            stop("Equicorrelation value must be comprise between 0 and 1",
              call. = FALSE)
            ### simulation -------
            betaMat <- do.call(cbind, mlParam)
            fittedMat <- Xvar %*% betaMat[1:nXvar, , drop = FALSE]
            epsilonMat <- sweep(-fittedMat, MARGIN = 1, STATS = Yvar, FUN = "+")
            if (udist %in% c("tnormal", "lognormal")) {
            muMat <- muHvar %*% betaMat[(nXvar + 1):(nXvar + nmuZUvar),
              , drop = FALSE]
            WuMat <- uHvar %*% betaMat[(nXvar + nmuZUvar + 1):(nXvar +
              nmuZUvar + nuZUvar), , drop = FALSE]
            WvMat <- vHvar %*% betaMat[(nXvar + nmuZUvar + nuZUvar +
              1):(nXvar + nmuZUvar + nuZUvar + nvZVvar), , drop = FALSE]
            } else {
            WuMat <- uHvar %*% betaMat[(nXvar + 1):(nXvar + nuZUvar),
              , drop = FALSE]
            WvMat <- vHvar %*% betaMat[(nXvar + nuZUvar + 1):(nXvar +
              nuZUvar + nvZVvar), , drop = FALSE]
            }
            # progressr::handlers('progress')
            set.seed(seed)
            p <- progressr::progressor(N)
            MetaList <- list()
            for (i in seq_len(N)) {
            draws <- mnorm::rmnorm(n = metaSim, mean = fittedMat[i, ],
              sigma = fillCov(stdvec = exp(WvMat[i, ]/2), corval = corval))
            if (S == 1) {
              metaFrontier <- apply(draws, 1, max)
            } else {
              metaFrontier <- apply(draws, 1, min)
            }
            groupFrontier <- draws[, which(group_var_list == group_var[i,
              ])]
            MetaList[[i]] <- list(groupFrontier, metaFrontier)
            p(sprintf("Simulation for observation %g", i))
            }
            groupFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
            FUN = function(z) MetaList[[z]][[1]]))
            metaFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
            FUN = function(z) MetaList[[z]][[2]]))
            MdMat <- if (S == 1) {
            metaFrontierMat - groupFrontierMat
            } else {
            groupFrontierMat - metaFrontierMat
            }
            colnames(groupFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
            colnames(metaFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
            colnames(MdMat) <- paste0("Md_Sim#", 1:metaSim)
            Mud <- numeric(N)
            for (j in seq_len(N)) {
            Mud[j] <- if (udist == "hnormal") {
              sqrt(2/pi) * exp(WuMat[j, which(group_var_list == group_var[j,
              ])]/2)
            } else {
              if (udist == "tnormal") {
              mean(muMat[, which(group_var_list == group_var[j, ])]) +
                exp(WuMat[j, which(group_var_list == group_var[j, ])]/2) *
                dnorm(mean(muMat[, which(group_var_list == group_var[j,
                  ])])/exp(WuMat[j, which(group_var_list == group_var[j,
                  ])]/2))/pnorm(mean(muMat[, which(group_var_list ==
                group_var[j, ])])/exp(WuMat[j, which(group_var_list ==
                group_var[j, ])]/2))
              } else {
              if (udist == "exponential") {
                exp(WuMat[j, which(group_var_list == group_var[j, ])]/2)
              } else {
                if (udist == "rayleigh") {
                exp(WuMat[j, which(group_var_list == group_var[j,
                  ])]/2) * sqrt(pi/2)
                } else {
                if (udist == "gamma") {
                  betaMat[nXvar + nuZUvar + nvZVvar + 1, which(group_var_list ==
                  group_var[j, ])] * exp(WuMat[j, which(group_var_list ==
                  group_var[j, ])]/2)
                } else {
                  if (udist == "lognormal") {
                  exp(mean(muMat[, which(group_var_list == group_var[j,
                    ])]) + WuMat[j, which(group_var_list == group_var[j,
                    ])]/2)
                  } else {
                  if (udist == "uniform") {
                    sqrt(12 * exp(WuMat[j, which(group_var_list ==
                    group_var[j, ])]/2))/2
                  } else {
                    if (udist == "genexponential") {
                    3/2 * exp(WuMat[j, which(group_var_list ==
                      group_var[j, ])]/2)
                    } else {
                    if (udist == "tslaplace") {
                      lambda <- betaMat[nXvar + nuZUvar + nvZVvar +
                      1, which(group_var_list == group_var[j,
                      ])]
                      exp(WuMat[j, which(group_var_list == group_var[j,
                      ])]/2) * (1 + 4 * lambda + 2 * lambda^2)/((1 +
                      lambda) * (1 + 2 * lambda))
                    } else {
                      if (udist == "weibull") {
                      exp(WuMat[j, which(group_var_list ==
                        group_var[j, ])]/2) * gamma(1 + 1/betaMat[nXvar +
                        nuZUvar + nvZVvar + 1, which(group_var_list ==
                        group_var[j, ])])
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
            ### return meta objects -------
            mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
            mlParam[[Ngroup + 1]] <- NA
            startVal[[Ngroup + 1]] <- NA
            optType[[Ngroup + 1]] <- NA
            nIter[[Ngroup + 1]] <- NA
            optStatus[[Ngroup + 1]] <- NA
            startLoglik[[Ngroup + 1]] <- NA
            mlLoglik[[Ngroup + 1]] <- NA
            gradient[[Ngroup + 1]] <- NA
            gradL_OBS[[Ngroup + 1]] <- NA
            gradientNorm[[Ngroup + 1]] <- NA
            invHessian[[Ngroup + 1]] <- NA
            conditionNums[[Ngroup + 1]] <- NA
            mlFitted <- apply(metaFrontierMat, 1, mean)
            mlResiduals <- Yvarm - mlFitted
            dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
            mlFitted = mlFitted)
            dataTable[[Ngroup + 1]]$logL_OBS <- NA
          } else {
            ## aos17d model -------
            if (modelType == "aos17d") {
            corval <- inputNumber("Enter a value for the equicorrelation: ")
            if (corval < 0.01 || corval > 0.99)
              stop("Equicorrelation value must be comprise between 0 and 1",
              call. = FALSE)
            cat("Stochastic Metafrontier: model aos17d", "\n")
            ### simulation -------
            betaMat <- do.call(cbind, mlParam)
            fittedMat <- Xvar %*% betaMat[1:nXvar, , drop = FALSE]
            WuMat <- uHvar %*% betaMat[(nXvar + 1):(nXvar + nuZUvar),
              , drop = FALSE]
            WvMat <- vHvar %*% betaMat[(nXvar + nuZUvar + 1):(nXvar +
              nuZUvar + nvZVvar), , drop = FALSE]
            epsilonMat <- sweep(-fittedMat, MARGIN = 1, STATS = Yvar,
              FUN = "+")
            # progressr::handlers('progress')
            set.seed(seed)
            p <- progressr::progressor(N)
            MetaList <- list()
            for (i in seq_len(N)) {
              sigmavsq <- exp(WvMat[i, ])
              sigmausq <- exp(WuMat[i, ])
              sigmasq <- sigmausq + sigmavsq
              if (S == 1) {
              vdraws <- tmvtnorm::rtmvnorm(n = metaSim, mean = sigmavsq *
                epsilonMat[i, ]/sigmasq, sigma = diag(sigmausq * sigmavsq/sigmasq),
                lower = epsilonMat[i, ], algorithm = "gibbs")
              } else {
              vdraws <- tmvtnorm::rtmvnorm(n = metaSim, mean = sigmavsq *
                epsilonMat[i, ]/sigmasq, sigma = diag(sigmausq * sigmavsq/sigmasq),
                upper = epsilonMat[i, ], algorithm = "gibbs")
              }
              draws <- mnorm::rmnorm(metaSim, mean = fittedMat[i, ],
              sigma = diag(sigmavsq))
              draws[, which(group_var_list == group_var[i, ])] <- fittedMat[i,
              which(group_var_list == group_var[i, ])] + vdraws[, which(group_var_list ==
              group_var[i, ])]
              # re-draw for non 'i' group
              SigMat <- fillCov(stdvec = exp(WvMat[i, ]/2), corval = corval)
              drawList <- lapply(seq_len(metaSim), FUN = function(x) mnorm::rmnorm(1,
              mean = SigMat[which(group_var_list == group_var[i, ]),
                -which(group_var_list == group_var[i, ])]/SigMat[which(group_var_list ==
                group_var[i, ]), which(group_var_list == group_var[i,
                ])] * vdraws[x, which(group_var_list == group_var[i,
                ])], sigma = SigMat[-which(group_var_list == group_var[i,
                ]), -which(group_var_list == group_var[i, ]), drop = FALSE] -
                1/SigMat[which(group_var_list == group_var[i, ]), which(group_var_list ==
                group_var[i, ])] * tcrossprod(SigMat[which(group_var_list ==
                group_var[i, ]), -which(group_var_list == group_var[i,
                ])])))
              draws[, which(group_var_list != group_var[i, ])] <- do.call(rbind,
              drawList)
              groupFrontier <- draws[, which(group_var_list == group_var[i,
              ])]
              if (S == 1) {
              metaFrontier <- apply(draws, 1, max)
              } else {
              metaFrontier <- apply(draws, 1, min)
              }
              MetaList[[i]] <- list(groupFrontier, metaFrontier)
              p(sprintf("Simulation for observation %g", i))
            }
            groupFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
              FUN = function(z) MetaList[[z]][[1]]))
            metaFrontierMat <- do.call(rbind, lapply(seq_len(length(MetaList)),
              FUN = function(z) MetaList[[z]][[2]]))
            MdMat <- if (S == 1) {
              metaFrontierMat - groupFrontierMat
            } else {
              groupFrontierMat - metaFrontierMat
            }
            colnames(groupFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
            colnames(metaFrontierMat) <- paste0("Md_Sim#", 1:metaSim)
            colnames(MdMat) <- paste0("Md_Sim#", 1:metaSim)
            Mud <- numeric(N)
            for (obs in seq_len(N)) {
              Mud[obs] <- sqrt(2/pi) * exp(WuMat[obs, which(group_var_list ==
              group_var[obs, ])]/2)
            }
            ### return meta objects -------
            mlDate[[Ngroup + 1]] <- format(Sys.time(), "Model was estimated on : %b %a %d, %Y at %H:%M:%S")
            mlParam[[Ngroup + 1]] <- NA
            startVal[[Ngroup + 1]] <- NA
            optType[[Ngroup + 1]] <- NA
            nIter[[Ngroup + 1]] <- NA
            optStatus[[Ngroup + 1]] <- NA
            startLoglik[[Ngroup + 1]] <- NA
            mlLoglik[[Ngroup + 1]] <- NA
            gradient[[Ngroup + 1]] <- NA
            gradL_OBS[[Ngroup + 1]] <- NA
            gradientNorm[[Ngroup + 1]] <- NA
            invHessian[[Ngroup + 1]] <- NA
            conditionNums[[Ngroup + 1]] <- NA
            mlFitted <- apply(metaFrontierMat, 1, mean)
            mlResiduals <- Yvarm - mlFitted
            dataTable[[Ngroup + 1]] <- cbind(dataTablem, mlResiduals = mlResiduals,
              mlFitted = mlFitted)
            dataTable[[Ngroup + 1]]$logL_OBS <- NA
            }
          }
          }
        }
      }
    }
  }
  # Return overall object -------
  returnObj <- list()
  returnObj$call <- cl
  returnObj$formula <- formula
  returnObj$S <- S
  returnObj$typeSfa <- typeSfa
  returnObj$Nobs <- c(table(group_var), N)
  names(returnObj$Nobs) <- c(group_var_list, "metafrontier")
  returnObj$Ngroup <- Ngroup
  returnObj$name_meta_var <- name_meta_var
  returnObj$nXvar <- nXvar
  if (udist %in% c("tnormal", "lognormal")) {
    returnObj$nmuZUvar <- nmuZUvar
  }
  returnObj$Yvarm <- Yvarm
  returnObj$logDepVar <- logDepVar
  returnObj$nuZUvar <- nuZUvar
  returnObj$nvZVvar <- nvZVvar
  returnObj$nParm <- nParm
  returnObj$udist <- udist
  returnObj$startVal <- do.call(cbind, startVal)
  colnames(returnObj$startVal) <- c(group_var_list, "metafrontier")
  returnObj$modelType <- modelType
  returnObj$olsParam <- olsParam  # list
  returnObj$olsStder <- olsStder  # list
  returnObj$olsSigmasq <- olsSigmasq  # list
  returnObj$olsLoglik <- olsLoglik  # list
  returnObj$olsSkew <- olsSkew  # list
  returnObj$olsM3Okay <- olsM3Okay  # list
  returnObj$CoelliM3Test <- CoelliM3Test  # list
  returnObj$AgostinoTest <- AgostinoTest  # list
  if (modelType %in% c("bpo04a", "bpo04b", "aos17a", "aos17b", "aos17c", "aos17d")) {
    names(returnObj$olsParam) <- group_var_list
    names(returnObj$olsStder) <- group_var_list
    names(returnObj$olsSigmasq) <- group_var_list
    names(returnObj$olsLoglik) <- group_var_list
    names(returnObj$olsSkew) <- group_var_list
    names(returnObj$olsM3Okay) <- group_var_list
    names(returnObj$CoelliM3Test) <- group_var_list
    names(returnObj$AgostinoTest) <- group_var_list
  } else {
    names(returnObj$olsParam) <- c(group_var_list, "metafrontier")
    names(returnObj$olsStder) <- c(group_var_list, "metafrontier")
    names(returnObj$olsSigmasq) <- c(group_var_list, "metafrontier")
    names(returnObj$olsLoglik) <- c(group_var_list, "metafrontier")
    names(returnObj$olsSkew) <- c(group_var_list, "metafrontier")
    names(returnObj$olsM3Okay) <- c(group_var_list, "metafrontier")
    names(returnObj$CoelliM3Test) <- c(group_var_list, "metafrontier")
    names(returnObj$AgostinoTest) <- c(group_var_list, "metafrontier")
  }
  returnObj$isWeights <- !all.equal(wHvar, rep(1, N))
  returnObj$optType <- unique(do.call(c, optType))
  returnObj$nIter <- do.call(c, nIter)
  names(returnObj$nIter) <- c(group_var_list, "metafrontier")
  returnObj$optStatus <- do.call(c, optStatus)
  names(returnObj$optStatus) <- c(group_var_list, "metafrontier")
  returnObj$startLoglik <- do.call(c, startLoglik)
  names(returnObj$startLoglik) <- c(group_var_list, "metafrontier")
  returnObj$mlLoglik <- do.call(c, mlLoglik)
  names(returnObj$mlLoglik) <- c(group_var_list, "metafrontier")
  returnObj$mlParam <- do.call(cbind, mlParam)
  colnames(returnObj$mlParam) <- c(group_var_list, "metafrontier")
  returnObj$dataTable <- dataTable
  names(returnObj$dataTable) <- c(group_var_list, "metafrontier")
  returnObj$gradient <- do.call(cbind, gradient)
  colnames(returnObj$gradient) <- c(group_var_list, "metafrontier")
  returnObj$gradL_OBS <- gradL_OBS  # list
  names(returnObj$gradL_OBS) <- c(group_var_list, "metafrontier")
  returnObj$gradientNorm <- do.call(c, gradientNorm)
  names(returnObj$gradientNorm) <- c(group_var_list, "metafrontier")
  returnObj$invHessian <- invHessian  # list
  names(returnObj$invHessian) <- c(group_var_list, "metafrontier")
  returnObj$conditionNums <- conditionNums  # list
  names(returnObj$conditionNums) <- c(group_var_list, "metafrontier")
  returnObj$hessianType <- if (hessianType == 1) {
    "Analytic Hessian"
  } else {
    if (hessianType == 2) {
      "BHHH Hessian"
    }
  }
  if (modelType %in% c("bpo04a", "bpo04b", "aos17a", "aos17b", "aos17c", "aos17d")) {
    returnObj$metaSim <- metaSim
  }
  if (modelType == "bpo04a") {
    returnObj$lpSimRes <- lpSimRes
  } else {
    if (modelType == "bpo04b") {
      returnObj$qpSimRes <- qpSimRes
    } else {
      if (modelType %in% c("aos17a", "aos17b", "aos17c", "aos17d")) {
        returnObj$MdMat <- MdMat
        returnObj$Mud <- Mud
        returnObj$groupFrontierMat <- groupFrontierMat
        returnObj$metaFrontierMat <- metaFrontierMat
      }
    }
  }
  returnObj$mlDate <- mlDate  #list
  if (udist %in% c("gamma", "lognormal", "weibull")) {
    returnObj$simDist <- simDist
    returnObj$Nsim <- Nsim
    returnObj$FiMat <- FiMat  # list
    names(returnObj$FiMat) <- c(group_var_list, "metafrontier")
  }
  class(returnObj) <- "sfametacross"
  gc()
  return(returnObj)
}

# print for sfametacross ----------
#' @rdname sfametacross
#' @exportS3Method print sfametacross
print.sfametacross <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call))
  cat("\n\n")
  cat("Likelihood estimates using", x$optType, "\n")
  cat(sfacrossdist(x$udist), "\n")
  cat("Status:", x$optStatus, "\n\n")
  cat(x$typeSfa, "\n")
  print.default(format(x$mlParam), print.gap = 2, quote = FALSE)
  invisible(x)
}

# Bread for Sandwich Estimator ----------
#' @rdname sfametacross
#' @exportS3Method sandwich::bread sfametacross
bread.sfametacross <- function(x, ...) {
  group_var <- x$dataTable[[x$Ngroup + 1]][x$name_meta_var][, 1]
  group_var_list <- sort(unique(group_var))
  frChoice <- displayMenu(c(group_var_list, "metafrontier"), title = "Which frontier do you want? ")
  if ((frChoice == x$Ngroup + 1) && x$modelType %in% c("aos17a", "aos17b", "aos17c",
    "aos17d")) {
    stop("No parameters are estimated for the metafontier for models 'aos17a', 'aos17b', 'aos17c', 'aos17d' \n",
      call. = FALSE)
  }
  ## something is still missing here depending on the type of frontier
  if (x$hessianType == "Analytic Hessian") {
    return(x$invHessian[[frChoice]] * x$Nobs[frChoice])
  } else {
    cat("Computing Analytical Hessian \n")
    Yvar <- model.response(model.frame(x$formula, data = x$dataTable))
    Xvar <- model.matrix(x$formula, rhs = 1, data = x$dataTable)
    if (x$udist %in% c("tnormal", "lognormal")) {
      muHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      uHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 4, data = x$dataTable)
    } else {
      uHvar <- model.matrix(x$formula, rhs = 2, data = x$dataTable)
      vHvar <- model.matrix(x$formula, rhs = 3, data = x$dataTable)
    }
    if (x$udist == "hnormal") {
      hessAnalytical <- chesshalfnormlike(parm = x$mlParam[, frChoice], nXvar = ncol(Xvar),
        nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
        Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights, S = x$S)
    } else {
      if (x$udist == "exponential") {
        hessAnalytical <- chessexponormlike(parm = x$mlParam[, frChoice],
          nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
          uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
          S = x$S)
      } else {
        if (x$udist == "tnormal") {
          if (x$scaling == TRUE) {
          hessAnalytical <- chesstruncnormscalike(parm = x$mlParam[, frChoice],
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
            uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          hessAnalytical <- chesstruncnormlike(parm = x$mlParam[, frChoice],
            nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar), nuZUvar = ncol(uHvar),
            nvZVvar = ncol(vHvar), muHvar = muHvar, uHvar = uHvar, vHvar = vHvar,
            Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights, S = x$S)
          }
        } else {
          if (x$udist == "rayleigh") {
          hessAnalytical <- chessraynormlike(parm = x$mlParam[, frChoice],
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
            uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          } else {
          if (x$udist == "uniform") {
            hessAnalytical <- chessuninormlike(parm = x$mlParam[, frChoice],
            nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
            uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
            S = x$S)
          } else {
            if (x$udist == "genexponential") {
            hessAnalytical <- chessgenexponormlike(parm = x$mlParam[,
              frChoice], nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              wHvar = x$dataTable$weights, S = x$S)
            } else {
            if (x$udist == "tslaplace") {
              hessAnalytical <- chesstslnormlike(parm = x$mlParam[, frChoice],
              nXvar = ncol(Xvar), nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar),
              uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
              wHvar = x$dataTable$weights, S = x$S)
            } else {
              if (x$udist == "gamma") {
              hessAnalytical <- chessgammanormlike(parm = x$mlParam[,
                frChoice], nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
                Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
                S = x$S, N = x$Nobs, FiMat = x$FiMat)
              } else {
              if (x$udist == "weibull") {
                hessAnalytical <- chessweibullnormlike(parm = x$mlParam[,
                frChoice], nXvar = ncol(Xvar), nuZUvar = ncol(uHvar),
                nvZVvar = ncol(vHvar), uHvar = uHvar, vHvar = vHvar,
                Yvar = Yvar, Xvar = Xvar, wHvar = x$dataTable$weights,
                S = x$S, N = x$Nobs, FiMat = x$FiMat)
              } else {
                if (x$udist == "lognormal") {
                hessAnalytical <- chesslognormlike(parm = x$mlParam[,
                  frChoice], nXvar = ncol(Xvar), nmuZUvar = ncol(muHvar),
                  nuZUvar = ncol(uHvar), nvZVvar = ncol(vHvar), muHvar = muHvar,
                  uHvar = uHvar, vHvar = vHvar, Yvar = Yvar, Xvar = Xvar,
                  wHvar = x$dataTable$weights, S = x$S, N = x$Nobs,
                  FiMat = x$FiMat)
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
    invHess <- invHess_fun(hess = hessAnalytical)
    colnames(invHess) <- rownames(invHess) <- names(x$mlParam)
    return(invHess * x$Nobs)
  }
}

# Gradients Evaluated at each Observation ----------
#' @rdname sfametacross
#' @exportS3Method sandwich::estfun sfametacross
estfun.sfametacross <- function(x, ...) {
  group_var <- x$dataTable[[x$Ngroup + 1]][x$name_meta_var][, 1]
  group_var_list <- sort(unique(group_var))
  frChoice <- displayMenu(c(group_var_list, "metafrontier"), title = "Which frontier do you want? ")
  if ((frChoice == x$Ngroup + 1) && x$modelType %in% c("aos17a", "aos17b", "aos17c",
    "aos17d")) {
    stop("No parameters are estimated for the metafontier for models 'aos17a', 'aos17b', 'aos17c', 'aos17d' \n",
      call. = FALSE)
  }
  return(x$gradL_OBS[[frChoice]])
}