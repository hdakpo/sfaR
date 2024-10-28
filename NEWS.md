# sfaR 1.0.1

## BUG FIXES

* In tests `sfacross` returns different output

***
# sfaR 1.0.0
Changes in 'sfaR' version to 1.0.0 (2023-06-13).

## NEW FEATURES

* Sample selection stochastic frontier model is introduced with function
`sfaselectioncross`.

* `texreg` package can be used for results output.

* Robust variance-covariance matrix can be obtained using `lmtest` package.

* Add reciprocal of efficiency.

* Add Battese and Coelli (1988) efficiency in the case of `sfalcmcross`.

## BUG FIXES

* `sfacross` previously returns the wrong sign of the gradient. Now the correct 
sign is returned.

## DEPRECATED & DEFUNCT

* `lcmcross` is now replaced by `sfalcmcross`. Associated methods are 
modified accordingly to the new class `sfalcmcross`.

## OTHER USER-VISIBLE CHANGES

* Package maintainer has changed from Yann Desjeux to K Hervé Dakpo.

* Arne Henningsen is now co-author of the package. Welcome onboard!

* Remove dependencies to packages `dplyr`, `emdbook`, `fBasics`, `gsl`, `MASS`,
`moments`, `numDeriv`, and `primes`.

* In the case of the truncated skewed Laplace distribution, starting value for
`lambda` is changed from 0.03 to 1.
 
* In the case of the uniform distribution, `theta` is reparameterized from
`theta = exp(Wu)` to `theta = sqrt(12) * exp(Wu/2)`.
 
* Remove `initStart` option from `lcmcross`. It is replaced with `whichStart` option.

* Add Roxygen comments.

***
# sfaR 0.1.1
Changes in 'sfaR' version 0.1.1 (2022-05-05).

## NEW FEATURES

* The package is enriched with a `NEWS` file.

## BUG FIXES

* The `fitted()` function now works properly for objects of class `sfacross`.

* When `extraPar = TRUE`, the `coef()` function now works properly on objects of class `lcmcross` for models with a number of classes (argument `lcmClasses`) equal to 4 or more.

## DEPRECATED & DEFUNCT
None.

## OTHER USER-VISIBLE CHANGES

* The `efficiencies()` documentation is enriched with an `lcmcross` example.

* Updated `DESCRIPTION` file.

* Updated startup message.

***
# sfaR 0.1.0
Changes in 'sfaR' version 0.1.0 (released on CRAN 2021-05-04).

## NEW FEATURES

* When the function `efficiencies()` is applied to an object of class `lcmcross` an additional value is returned. This new value, `ineff_c#`, returns the conditional inefficiency for observations classified under class# only. Therefore a value is assigned to `ineff_c#` only if the observations is classified under class# (i.e. `Group_c = #`).

## BUG FIXES

* In `efficiencies()` function applied to an object of class `lcmcross`, the values for `u_c`, `teJLMS_c`, `Group_c`, and `PosteriorProb_c` are now returned correctly. In the previous version (v 0.0.91) all observations were applied the same values when the argument `lcmClasses` was set to 3 and more. This bug is now fixed.

## DEPRECATED & DEFUNCT
None.

## OTHER USER-VISIBLE CHANGES

* Updated `DESCRIPTION` file.

* In the `efficiencies()` function applied to an object of class `lcmcross`, the returned elements have been reordered in order to better reflect the logic of the `lcmcross` model.

* Updated `efficiencies()` documentation to incorporate the new `ineff_c#` elements returned by the function, as well as the reordering of the returned values.

* Updated startup message.

***
# sfaR 0.0.91
First public release (released on CRAN 2021-03-05).
