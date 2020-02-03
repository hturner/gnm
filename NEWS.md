Changes in gnm 1.1-1
====================

Bug fixes
---------

 * `confint.profile.gnm()` now works for a single parameter  ![(#10)](https://github.com/hturner/gnm/issues/10)
 * `gnm()` now works when `eliminate` is specified as `NULL` ![(#14)](https://github.com/hturner/gnm/issues/14)
 * allow `set` argument of `getContrasts()` to be numeric, as documented
 * convert old tests to unit tests.
 
Changes in gnm 1.1-0
====================

Changes in behaviour
--------------------

 * generalize `Diag` and `Symm` to work with factors with levels that are not in alphabetical order. Factor levels are now only sorted by `Diag` if the input factors have different sets of levels.
 * make `se()` generic to allow methods to be added (e.g. as in **logmult** package).

Improvements
------------

 * C routines now registered to avoid accidental clashes with other packages.
   - use of `R_forceSymbols` routine requires R >= 3.0.0.
 * use jss.bst vs chicago.bst in vignette.
 * update imports to include recommended packages.
 * avoid warnings regarding recycling a length 1 array.
 * avoid using `print` or `cat` outside print methods, so print output always optional (e.g. by setting `verbose` argument or using `suppressMessages`).

Bug fixes
---------

 * environment of formula preserved when using `gnm()` with `eliminate` argument.
 * allow factor response in binomial gnms.
 * allow matrix response in quasibinomial gnms.
 * make `x = FALSE` also work when using `gnm()` with `eliminate` argument.
 * allow response in formula to be specified as an expression (e.g. `D/E`) when formula uses `instances`.

Changes in gnm 1.0-8
====================

Improvements
-------------

 * now use lazy data loading.
 * improvements to vignette (with thanks to Michael Friendly).
 * copyright notices added to source files to clarify authorship, with
   appropriate credit given to contributors in Rd and DESCRIPTION files.

Changes in behaviour
--------------------

 * predict.gnm now includes eliminate term in predictions on new data.

Bug fixes
---------

 * `expandCategorical` now works when there are no covariates in the data.
 * better handling of single-column model matrices.
 * `predict.gnm` now works with `se.fit = TRUE` for models with eliminated 
   terms; correctly handles new data without all levels of homogeneous 
   factors present, and respects contrasts settings.
 * environment of formula preserved when using instances.

Changes in gnm 1.0-7
====================

Bug fixes
---------

 * corrected use of `anova.glmlist` in `anova.gnm`.


Changes in gnm 1.0-6
====================

Bug fixes
---------

 * added catch for when deviance becomes `NaN`.

Changes in gnm 1.0-5
====================

Improvements
------------

 * eliminated coefficients now returned as named vector.
 * step-quartering introduced to start-up iterations to avoid increasing 
   deviance.

Changes in behaviour
--------------------

 * `gnm` no longer restarts if algorithm fails - better to provide 
   improved starting values in this case. 

Bug fixes
---------

 * fixed bug in the way `etastart` is used to initialise the linear 
   parameters when `eliminate` is used as well. 


Changes in gnm 1.0-4
====================

Bug fixes
---------

 * restarting mechanism now reinitialises correctly.
 * removed call to external C function that is no longer available.
 * `gnm` now works with `eliminate` argument when remaining linear part of 
   predictor only involves one parameter.

Changes in gnm 1.0-3
====================

Improvements
------------

 * `pickCoef` extended to allow fixed pattern matching and to optionally
   return actual coefficients rather than their indices.
 * `gnm` now looks for exact match in coefficient names when a single 
   character string is passed to the `constrain` argument before treating 
   as regular expression.
 * `hatvalues.gnm` has been reimplemented to work more efficiently for large
   model matrices.
 * `"nonlin"` terms defined for homogeneous factors will now accept factors
   specified as an interaction (using `:`).

Bug fixes
---------

 * results now returned in original order for models fitted with `eliminate`
   argument.
 * bug introduced into `residSVD` reverted so that now correctly aggregates 
   working residuals.
 * `anova.gnm` now works when model is a single `"nonlin"` term.


Changes in gnm 1.0-2
====================

Improvements
------------

 * factors specified as homogenous in nonlin functions can now be 
   specified as interactions of factors.

Bug fixes
---------

 * fixed bug so that variables handled correctly in nonlinTerms.
 * corrected rank calculation for constrained models.
 * removed calls to Internal


Changes in gnm 1.0-1
====================

New Features
------------

 * added `meanResiduals` function
 * `check` argument added to getContrasts.

Improvements
------------

 * added example SVD calculation to `?wheat`; also in vignette.

Bug fixes
---------

 * added `update.gnm` so that nonlinear terms were no ordered as linear, 
   first order terms.


Changes in gnm 1.0-0
====================

Improvements
------------

 * eliminated coefficients now treated entirely separately, in particular 
   the design matrix no longer has columns for these coefficients, making
   the algorithm far more efficient for models with many eliminated
   coefficients.

 * more reliable calculation of rank

Changes in Behaviour
--------------------

 * `ofInterest` and `constrain` now index non-eliminated coefficients only.

 * eliminated coefficients now returned as attribute of returned 
   coefficient vector.

 * `"lsMethod"` argument to `gnm` removed as now the LAPACK routines are  
   always used to determine the least squares solution at the heart of 
   the fitting algorithm. Hence `qrSolve` and `cholInv` deprecated.

 * the `"eliminate"`, `"onlyFirstCol"` and `"onlyNonElim"` arguments to
   `MPinv` have been removed as no longer used.

Bug fixes
---------

 * `etastart` now works for models with no linear parameters.

 * `anova` now ignores terms that are completely constrained.

Changes in gnm 0.10-0
=====================

Improvements
------------

 * `mustart`/`etastart` now used to obtain starting values for linear and 
   nonlinear parameters separately, improving performance.

Changes in Behaviour
--------------------

 * `expandCategorical` now groups together individuals with common covariate
   values, by default. New `group` argument added to switch this behaviour.

Bug fixes
---------

 * `print.profile.gnm` now prints full result.

 * data now read in correctly for Lee-Carter example in vignette.

Changes in gnm 0.9-9
====================

New Features
------------

 * `etastart` and `mustart` arguments added to `gnm`.


Changes in gnm 0.9-8
====================

Improvements
------------

 * `gnm` now returns `data` argument as `glm` does.


Changes in gnm 0.9-7
====================

Bug fixes
---------

 * more minor corrections in documentation.


Changes in gnm 0.9-6
====================

Bug fixes
---------

 * minor corrections in documentation.


Changes in gnm 0.9-5
====================	

Improvements
------------

 * `getContrasts` can now estimate _scaled_ contrasts with more flexibility in
   how the reference level is defined. 
 * changed tolerance level in checkEstimable to `1e6 * .Machine$double.eps` as
   previous tolerance too strict for some examples. 

Changes in Behaviour
--------------------

 * `getContrasts` now only handles one set of parameters at a time.
 * use of `Const` is now restricted to the symbolic predictors of `"nonlin"` 
   functions.
 * `Nonlin` - the wrapper function for plug-in functions - is now defunct. Use
   `"nonlin"` functions to specify custom nonlinear terms.

Bug fixes
---------

 * `plot.gnm` now uses standardised Pearson residuals for plot `which = 5` so 
   that the Cook's distance contours are correct


Changes in gnm 0.9-4
====================

New Features
------------

 * `predict` now implemented for `"gnm"` objects	   

Improvements
------------

 * results formatted as contingency tables where appropriate by extractor 
   functions (`fitted`, etc), rather than `gnm`

Changes in Behaviour
--------------------

 * default for `match` argument of `nonlinTerms` now zero vector (i.e. no
   matching to arguments of `call` by default)

Bug fixes
---------

 * `termPredictors` now works on `"gnm"` objects fitted with `glm.fit`
 * intercept removed when `eliminate` argument of `gnm` is non-`NULL`
 * models with all parameters eliminated now summarised sensibly
 * `Diag` and `Symm` now work for factors of length 1
 * as`gnm` now returns object with `"gnm"`-type terms component
 * print method for `"profile.gnm"` objects now exported


Changes in gnm 0.9-3
====================

New Features
------------

 * added `DrefWeights` for computing the weights in a diagonal reference 
   term and the corresponding standard errors.
    

Improvements
------------

 * `"assign"` attribute now attached to the parameter vector when passed to 
   start functions defined by `"nonlin"` functions, specifying the 
   correspondence between parameters and predictors in the nonlinear term.

Bug fixes
---------

 * start function in `Dref` now identifies weight parameters correctly.

 * can now evaluate term predictors for `"nonlin"` terms that depend on 
   covariates.


Changes in gnm 0.9-2
====================    

Improvements
------------

 * Calls to `"nonlin"` functions now evaluated in the same environment and 
   enclosure as call to create model frame, so `"nonlin"` functions should 
   be able to find variables in gnm calls - potentially useful for setting 
   starting values.

Bug fixes
---------

 * `gnm` algorithm now reinitiates correctly when restarting after 
   non-convergence.

 * `gnm` now works correctly when a model is specified with nonlinear terms 
   in between linear terms.


Changes in gnm 0.9-1
====================

New Features
------------

 * introduction of functions of class `"nonlin"` for the unified specification 
   of nonlinear terms. `Mult`, `Exp`, `Dref` and `MultHomog` have all been
   converted to functions of this class.

 * added `Inv` to specify the reciprocal of a predictor.

 * added `Const` to specify a constant in a predictor.

 * added `instances` to specify multiple instances of a nonlinear term.
    

Improvements
------------

 * nonlinear terms can now be nested.

 * `Exp` can now be used outside of `Mult` or to exponentiate part of a 
   constituent multiplier.

Changes in Behaviour
--------------------

 * to accommodate the increased functionality introduced by `"nonlin"` 
   functions, new labelling conventions have been introduced. In 
   particular, most `"nonlin"` functions use argument-matched parameter 
   labels.

 * in the new implementation of `Dref` the `formula` argument has been
   re-named `delta` to provide more informative parameter labels under the
   new conventions.

Bug fixes
---------

 * specifying `ofInterest = "[?]"` in `gnm` now works as documented.


Changes in gnm 0.8-5
====================

New Features
------------

    none in this release

Improvements
------------

 * added a new `ridge` argument to gnm, to allow some control over
   the Levenberg-Marquardt regularization of the internal least squares
   calculation

 * changed the default ridge constant to 1e-8 (from 1e-5), to increase
   speed of convergence (especially in cases where there are infinite
   parameter estimates)

 * modified the `"qr"` method so that it no longer checks for rank deficiency
   (it was both unreliable, and not necessary since the matrix is regularized 
   prior to solving)

 * substantial speed improvements in model fitting when there are large 
   numbers of eliminated parameters, achieved mainly via a new internal 
   function `cholInv1`.  Corresponding example timings changed in the 
   Overview document (vignette).

 * speed improvements in `vcov.gnm` when there are eliminated parameters; new
   logical argument `use.eliminate` gives control over this

 * in `getContrasts`, added new arguments `dispersion` and `use.eliminate`, 
   both of which are passed on to `vcov`

 * implemented faster alternatives to `ifelse` in `gnmFit`

 * speed gains from use of `tcrossprod`.  Because of this the gnm
   package now requires R 2.3.0 or later.

Changes in Behaviour
--------------------

 * in `gnm`, changed the default value of argument x to `TRUE` (it was
   previously `FALSE`)

 * in `checkEstimable`, changed the name of the first argument from 
   `coefMatrix` to `combMatrix` (to reflect better that it is a matrix of 
   coefficient *combinations*); and changed the default tolerance value to 
   one which should give more reliable results.  Also, more fundamentally, 
   changed the check to be whether combinations are in the column span of 
   `crossprod(X)` instead of the row span of `X`; the results should be the same, 
   but the new version is much faster for large n.

 * `model.matrix.gnm` no longer passes extra arguments to gnm as it's unlikely 
   to be useful/sensible. For the same reasons it will not pass extra 
   arguments to `model.frame`, unlike `model.matrix.lm`

 * `getContrasts` now results in a list only when the `sets` argument itself
   is a list;  otherwise (i.e., normally) the result is a single object 
   (rather than a list of objects) of class `qv`

Bug fixes
---------

 * fixed a bug in internal function `quick.glm.fit`, which greatly improves
   its performance.  Also changed the default value of the `nIter` argument
   from 3 to 2.

 * fixed a small bug in `demo(gnm)`

 * fixed a bug in `vcov.gnm`, which previously gave an error when data were 
   of class `"table"`)

 * fixed `summary.gnm` so that it now takes proper account of the dispersion
   argument

 * in `se`, added new arguments `Vcov` and `dispersion`; the latter fixes a 
   bug, while the former minimizes wasted computation in `summary.gnm`

 * fixed bug in `model.matrix.gnm` so that it can compute the model matrix 
   even when original data is not available - unless model frame has not been
   saved.  Original data still needed to update model frame - this is the 
   same as for glms, etc.

 * fixed bug in `gnm` so that reconstructing `"table"`-class data works for 
   models with weights/offsets 


Changes in gnm 0.8-4
====================

New Features
------------

 * added `"gnm"` methods for `profile` and `confint`. Use of `alpha` argument 
   differs slightly from `"glm"` methods: see help files.

 * `constrain` argument to `gnm` now supplemented by `constrainTo` 
   argument, allowing specification of values to which parameters should be 
   constrained.

 * `gnm` now has `ofInterest` argument to specify a subset of coefficients 
   which are of interest - returned in `ofInterest` component of `"gnm"` 
   object as named numeric vector. `print` summaries of model object/its 
   components extracted by accessor functions only print coefficients of 
   interest and (where appropriate) methods for `"gnm"` objects select 
   coefficients of interest by default.

 * added `ofInterest` and `ofInterest<-` to extract/replace `ofInterest`
   component of `"gnm"` object.

 * added `parameters` which returns coefficient vector with constrained
   parameters replaced by their constrained value.

 * added `pickCoef` function to aid selection of coefficients - returns 
   numeric indices of coefficients selected by Tk dialog or regular 
   expression matching.

Improvements
------------

 * `constrain` argument to `gnm` now accepts a regular expression to match 
   against coefficient names.


Changes in Behaviour
--------------------

 * `constrain` component of `"gnm"` objects is now a numeric, rather than 
   logical, vector of indices.

 * all `"gnm"` methods for which a subset of the coefficients may be 
   specified by numeric indices now interpret those indices as 
   referencing the full coefficient vector (not just non-eliminated
   parameters).

 * `gnm` now preserves order of terms rather than moving all linear 
   terms to the start (this fixes bug in `anova.gnm`).

 * the `"pick"` option for the `constrain` argument to `gnm` and the 
   `estimate` argument to `se` has been replaced by `"[?]"` to avoid
   possible conflict with coefficient names/regular expressions.

Bug Fixes
---------

 * fixed bug in `se` so will now work for single parameter.

 * fixed bug in `summary.gnm` so will now work for models with one 
   parameter.

 * fixed bug in `anova` so that rows of returned table are correct for 
   models with eliminated terms.

 * fixed bug in `eliminate` so that it now accepts interactions.

 * fixed bug in `MPinv` so that it works for models in which all parameters 
   are eliminated.


Changes in gnm 0.8-3
====================

Improvements
------------

 * improved use of functions from other packages

Bug Fixes
---------

 * fixed bug in `asGnm.lm` where object not fully identified

 * corrected maintainer address in DESCRIPTION!


Changes in gnm 0.8-2
====================

New Features
------------

 * added demonstration script to run using `demo`

 * added package help file, opened by package?gnm

Improvements
------------

 * improved existing documentation


Changes in gnm 0.8-1
====================

New Features
------------

 * added the `method` argument to `MPinv`, to allow the method of 
   calculation to be specified.  Permitted values are `"svd"` to
   compute the pseudo-inverse by singular value decomposition, 
   and `"chol"` to use the Cholesky decomposition instead.  The latter
   is valid only for symmetric matrices, but is usually faster
   and more accurate.

 * added the `lsMethod` argument to `gnm`, to allow specification of 
   the numerical method used for least-squares calculations in the
   core of the iterative algorithm.  Permitted options are `"chol"` and 
   `"qr"`.

 * added new function `qrSolve`, which behaves like `base::qr.coef` but 
   in the non-full-rank case gives the minimum-length solution rather
   than an arbitrary solution determined by pivoting.

 * added `.onUnload` so that compiled code is unloaded when namespace of 
   package is unloaded using `unloadNamespace`.

 * added `coef` argument to `model.matrix.gnm` so that the model matrix 
   can be evaluated at any specified value of the parameter vector. 

 * added as`gnm` generic to coerce linear model objects to gnm objects.

 * added `exitInfo` for printing numerical details of last iteration on 
   non-convergence of `gnm`.

 * added new dataset, friend, to illustrate a workaround to fit a
   homogeneous RC(2) using `gnm` - documented in help file for `MultHomog`.

Improvements
------------

 * `gnm` now takes less time per (main) iteration, due to improvements 
   made internally in the iterative algorithm.  These include 
   pre-scaling of the local design matrix, and Levenberg-Marquardt 
   adjustment of the least-squares solvers so that rank determination
   is no longer necessary.  

 * the default convergence tolerance has been tightened (from 1e-4 to 1e-6)

 * modified `model.matrix.gnm` so it can be used when only the namespace 
   of gnm is loaded.

Bug Fixes
---------

 * fixed bug in `gnm` so that `subset` now works with table data.

 * fixed bug in `model.matrix.gnm` so can construct model matrix from `"gnm"` 
   object even when original call not made in `.GlobalEnv`.

 * fixed bug in the examples on help page for `House2001` data.

 * fixed bug so that `formula` in `gnm` now accepts `.` in formulae 
   even when `eliminate = NULL`.

 * fixed bug in `getContrasts`, so that the first two columns
   of the qvframe component of each element of the result list
   are correctly named as "estimate" and "SE", as required for
   objects of class "qv".


Changes in gnm 0.8-0
=======================

New Features
------------

 * added `"model.matrix"` option for `method` argument of `gnm` so that 
   model matrix can be obtained much faster. The new method is used in 
   `model.matrix.gnm` and `vcov.gnm`.

 * added new utility function `residSVD`, to facilitate the calculation
   of good starting values for parameters in certain `Mult` terms.

 * added new dataset `House2001`, to illustrate the use of `gnm` in 
   Rasch-type scaling of legislator votes.

 * added new utility function `expandCategorical` for expanding data frame
   on the basis of a categorical variable.

 * added `formula.gnm` method - returns formula from `"gnm"` object excluding
   the `eliminate`d factor where necessary.

Improvements
------------

 * `gnm` now takes less time to run due to improvements made in internal 
   functions.

 * the fitting algorithm used by `gnm` now copes better with zero-valued 
   residuals.

 * output given by `gnm` when `trace = TRUE` or `verbose = TRUE` is now 
   displayed as it is generated on console-based versions of R. 

 * `plot.gnm` now includes option `which = 5` as in `plot.lm` in R >= 2.2.0. 
   Now has separate help page.	

 * the `constrain` argument to `gnm` now accepts the names of parameters.

 * the `formula` argument to `gnm` now accepts `.` as described in 
   `?terms.formula`, ignoring eliminated factor if in `data`.

 * interface for `se` extended - can now use to find standard errors for 
   all parameters or (a selection of) individual parameters in a gnm model.

 * made it possible to use `gnm` with alternative fitting function.

 * `".Environment"` attribute now attached to `"gnm"` objects so that gnm 
   package loaded when workspace containing `"gnm"` objects is loaded.	

	
Changes in Behaviour
--------------------

 * start-up iterations now only update column of design matrix required in 
   next iteration. Therefore plug-in functions using the default start-up 
   procedure for nonlinear parameters need a `localDesignFunction` with the 
   argument `ind` specifying the column that should be returned.

 * modified output given by `gnm` when `trace = TRUE`: now prints initial
   deviance and the deviance at the end of each iteration.

 * modified updates of linear parameters in starting procedure: now offset 
   contribution of fully specified terms only.

 * results of `summary.gnm`, `vcov.gnm` and `coef.gnm` now include any 
   eliminated parameters. Print methods have been added for `"vcov.gnm"` and 
   `"coef.gnm"` objects so that any eliminated parameters are not shown.

 * `Mult` terms are no longer split into components by `anova.gnm`,
   `termPredictors.gnm`, `labels.gnm` or the `"assign"` attribute of the 
   model matrix - consistent with `terms` output. 

 * the `eliminate` argument to `gnm` must now be an expression that 
   evaluates to a factor - this reverts the extension of 0.7-2.	

 * when using `gnm` with `constrain = "pick"`, the name(s) of the chosen 
   parameter(s) will replace `"pick"` in the returned model call.

 * `getContrasts` now uses first level of a factor as the reference level
   (by default).

 * `gnmControl` replaced by arguments to `gnm`.

 * `gnm` now uses `glm.fit` for linear models (with control parameters at 
   the `gnm` defaults) unless `eliminate` is non-`NULL`.

 * `vcov.gnm` and `summary.gnm` now return variance-covariance matrices 
   including any aliased parameters.

 * `summary.gnm` now returns standard errors with test statistics etc, 
   where estimated parameters are identified.

Bug Fixes
---------

 * fixed bug in `summary.gnm`, `anova.gnm`, `termPredictors.gnm` and 
   `model.matrix.gnm` where search for model variables was incorrect. 

 * fixed bug preventing estimation of weight parameters in `Dref` terms and 
   changed default starting values so that these parameters no longer sum 
   to one or appear to be estimable. 

 * corrected options for `method` argument in `gnm` help file: replaced 
   `method = "coef"` with `method = "coefNames"`.

 * fixed bug in `gnm` so that it can handle tables with missing values when
   formatting components of fit.

 * `hatvalues.gnm` now works for objects produced from table data.

 * `residuals.gnm` now returns table not matrix when `type = "deviance"` 
   for `"gnm"` objects produced from table data.

 * `hatvalues.gnm`, `cooks.distance.gnm` and `plot.gnm` now handle cases 
   which are fitted exactly (giving a hat value of 1).

 * example fitting proportional odds model in `backPain` help file now works.

 * fixed bug in `Mult` terms so that an offset can be added to a constituent
   multiplier without an unspecified intercept being added also.

 * `gnm` argument `constrain = "pick"` now allows selection of more than 
   one constraint and is compatible with use of `eliminate`.

 * `gnm` can now fit models which only have the term specified by `eliminate`.


Changes in gnm 0.7-2
=======================

Improvements
------------

 * Extended use of the `eliminate` argument of `gnm` to allow crossed 
   factors - this also fixes bug which occurred when interactions were 
   eliminated in the presence of lower order terms involving other factors 

Changes in Behaviour
--------------------

 * `vcov` returned by `gnm` now has no rank attribute (as before, the
   rank is returned as the separate component `rank`).

Bug Fixes
---------

 * Changed the calculation of `df.residual` returned by `gnm` to 
   correctly take account of zero-weighted observations (as in `glm`).

 * When `gnm` is called with arguments `x = TRUE` or `VCOV = TRUE`, the 
   returned matrices now include columns of zeros for constrained 
   parameters.

 * Corrected evaluation of model frame in `gnm` so that if data is missing, 
   variables are taken from `environment(formula)`, as documented. Modified 
   evaluation of plug-in functions to be consistent with this, i.e. 
   objects are taken from `environment(formula)` if not in model frame.

 * `MPinv` now checks that the diagonal elements of an `eliminate`d 
   submatrix are all non-zero and reports an error otherwise.



Changes in gnm 0.7-1
=======================

New Features
------------

 * `Topo` introduced for creating topological interaction factors.

 * `anova` implemented for objects of class `c("gnm", "glm")`.


Improvements
------------

 * Diagnostic messages given by `gnm` have been improved.

 * Step-halving introduced in main iterations of `gnm` to ensure deviance 
   is reduced at every iteration.

 * `getContrasts` now (additionally) reports quasi standard errors, when
   available.

 * Calls to `gnm` plug-in functions are now evaluated in the environment 
   of the model frame and the enclosing environment of the parent frame 
   of the call to `gnm`. This means that variables can be found in a 
   more standard fashion.


Changes in Behaviour
--------------------

 * The `data` argument of `Nonlin` is defunct: `Nonlin` now identifies 
   variables to be added to the model frame as those passed to unspecified
   arguments of the plug-in function or those identified by a companion 
   function to the plug-in, which is of a specified format.

 * The (optional) `start` object returned by a plug-in function can no 
   longer be a function, only a vector. However it may now include `NA` 
   values, to indicate parameters which may be treated as linear for the 
   purpose of finding starting values, given the non-`NA` values.


Bug Fixes
---------

 * The `eliminate` argument of `gnm` now handles functions of variables in 
   the given formula e.g. `~ strata(A, B), ~ as.factor(A):as.factor(B)`, etc. 

 * `gnm` was giving an error for models with either no linear parameters, 
   or none specified by the `start` argument, this is now fixed.

 * Long calls to plug-in functions caused problems in parsing the model 
   formula: now fixed.

 * `gnm` now only restarts after failing if there are unspecified nonlinear
   parameters.

 * `gnm` now returns `NULL` if model fails.

 * Bug fixed in calculation of starting values for `gnm` that occurred when some
   parameters were constrained.
	
