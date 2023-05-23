library("sfaR")
data("utility")
data("ricephil")
data("electricity")

## Using data on fossil fuel fired steam electric power
## generation plants in U.S.  Cobb-Douglas (cost function)
## half normal with heteroscedasticity
cd_u_h <- sfacross(formula = log(tc/wf) ~ log(y) + log(wl/wf) +
  log(wk/wf), udist = "hnormal", uhet = ~regu, data = utility,
  S = -1, method = "bfgs")
logLik(cd_u_h)
all.equal(c(logLik(cd_u_h)), sum(logLik(cd_u_h, individual = TRUE)[["logLik"]]))
round(coef(summary(cd_u_h)), 3)
t(sapply(efficiencies(cd_u_h), function(x) round(summary(x),
  3)))

# Cobb-Douglas (cost function) truncated normal with
# heteroscedasticity
cd_u_t <- sfacross(formula = log(tc/wf) ~ log(y) + log(wl/wf) +
  log(wk/wf), udist = "tnormal", muhet = ~regu, data = utility,
  S = -1, method = "bhhh")
logLik(cd_u_t)
all.equal(c(logLik(cd_u_t)), sum(logLik(cd_u_t, individual = TRUE)[["logLik"]]))
round(coef(summary(cd_u_t)), 3)
t(sapply(efficiencies(cd_u_t), function(x) round(summary(x),
  3)))

# Cobb-Douglas (cost function) truncated normal with
# scaling property
cd_u_ts <- sfacross(formula = log(tc/wf) ~ log(y) + log(wl/wf) +
  log(wk/wf), udist = "tnormal", muhet = ~regu, uhet = ~regu,
  data = utility, S = -1, scaling = TRUE, method = "mla")
logLik(cd_u_ts)
all.equal(c(logLik(cd_u_ts)), sum(logLik(cd_u_ts, individual = TRUE)[["logLik"]]))
round(coef(summary(cd_u_ts)), 3)
t(sapply(efficiencies(cd_u_ts), function(x) round(summary(x),
  3)))

## Using data on Philippine rice producers Cobb Douglas
## (production function) generalized exponential, and
## Weibull distributions
cd_p_ge <- sfacross(formula = log(PROD) ~ log(AREA) + log(LABOR) +
  log(NPK) + log(OTHER), udist = "genexponential", data = ricephil,
  S = 1, method = "bfgs")
logLik(cd_p_ge)
all.equal(c(logLik(cd_p_ge)), sum(logLik(cd_p_ge, individual = TRUE)[["logLik"]]))
round(coef(summary(cd_p_ge)), 3)
t(sapply(efficiencies(cd_p_ge), function(x) round(summary(x),
  3)))

## Using data on U.S. electric utility industry Cost
## frontier Gamma distribution
cd_u_g <- sfacross(formula = log(cost/fprice) ~ log(output) +
  log(lprice/fprice) + log(cprice/fprice), udist = "gamma",
  uhet = ~1, data = electricity, S = -1, method = "bfgs", simType = "halton",
  Nsim = 200, hessianType = 2)
logLik(cd_u_g)
all.equal(c(logLik(cd_u_g)), sum(logLik(cd_u_g, individual = TRUE)[["logLik"]]))
round(coef(summary(cd_u_g)), 3)
t(sapply(efficiencies(cd_u_g), function(x) round(summary(x),
  3)))
