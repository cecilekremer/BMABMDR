# EFSA Bayesian Benchmark Dose Response Platform

## BMABMDR 0.1.5

* included option to compare all Quantal models to saturated model

## BMABMDR 0.1.4

* added option for informative priors in covariate analysis (not to be included in WebApp)

## BMABMDR 0.1.3

* fixed covariates likelihood fold change parameter

* fixed null model start values

## BMABMDR 0.1.2

* fixed start value BMD if 0

## BMABMDR 0.1.1

* changed mistake in BMD start value calculation in covariate analysis

## BMABMDR 0.1.0

* fixed error in basic.plot()

## BMABMDR 0.0.0.9087

* fixed error in modelTestQ()

## BMABMDR 0.0.0.9086

* fixed error in covariate analysis when using only lognormal distribution

## BMABMDR 0.0.0.9085

* fixed error in covariate analysis when using only lognormal distribution

## BMABMDR 0.0.0.9084

* updated modelTest w/ possibility to test each model against saturated model

## BMABMDR 0.0.0.9083

* implemented modelTest w/ possibility to test each model against saturated model

## BMABMDR 0.0.0.9083

* allow litters of size 1 for clustered quantal data

## BMABMDR 0.0.0.9082

* fixed dataprep for maxDose<1

## BMABMDR 0.0.0.9081

* fixed anydoseresponseQ()

## BMABMDR 0.0.0.9080

* fixed error in modelTestQ()

## BMABMDR 0.0.0.9079

* implemented anydoseresponseQ without using 'brms' package

## BMABMDR 0.0.0.9078

* fix error in modelTest()

## BMABMDR 0.0.0.9077

* fix error in Laplace function

## BMABMDR 0.0.0.9076

* fixed error in quantile() function

## BMABMDR 0.0.0.9075

* added output in Laplace function (for use in simulations)

* fixed error in solving Hessian

## BMABMDR 0.0.0.9074

* updated modelTest function

## BMABMDR 0.0.0.9073

* provide error message if clustered = TRUE but some litters contain only one observation

## BMABMDR 0.0.0.9072

* fixed error in plot for heavily left-skewed distribution

## BMABMDR 0.0.0.9071

* fixed error in sampling

## BMABMDR 0.0.0.9070

* use default start values for MCMC if optimizing step not working

## BMABMDR 0.0.0.9069

* change in weights calculation to prevent error in solve(), now using Moore-Penrose pseudoinverse of hessian

## BMABMDR 0.0.0.9068

* fixed warnings that trigger an error in R-4.3

## BMABMDR 0.0.0.9067

* added 'Covariate' column to output of full.laplace_MA_Cov in summary table

## BMABMDR 0.0.0.9065

* added option for extended dose range upper value (default = 3)

* fixed error Covariate analysis (BMD prior if maxDose/2 < 0.5)

## BMABMDR 0.0.0.9064

* added model-averaged posterior of bkg and maxy (not relevant for R4EU platform)

* fixed error in plot function

## BMABMDR 0.0.0.9063

* fixed graphical presentation of BMD posterior

## BMABMDR 0.0.0.9062

* fixed code that occasionally triggers error in ModelTest

* plot_prior: truncated prior for d

* implemented user-defined prior for d (in PREP_DATA_*: option prior.d = 'custom')

* added License file

## BMABMDR 0.0.0.9061

* added option to plot function to indicate whether only converged models (conv = TRUE) should be plotted, defaults to FALSE

* set MA_conv equal to MA if all models converged

* updated plot legends

## BMABMDR 0.0.0.9060

* fixed likelihood L4_Q clustered

* fix compatibility with ggplot2 v3.4.2

## BMABMDR 0.0.0.9059

* warning message sampling specifying whether models are excluded for Bridge or Laplace

* exclude models that are not included in MA (due to fitting problems) from plots

* reversed BF in modelTest to correspond to Platform documentation

* fixed error induced by modelTest (likelihood calculation)

## BMABMDR 0.0.0.9058

* reversed BF in modelTest

* updated package documentation

## BMABMDR 0.0.0.9057

* Fixed covariate plot

## BMABMDR 0.0.0.9056

* Updated quantal covariate analysis

## BMABMDR 0.0.0.9055

* Fixed error continuous covariate analysis (to do for quantal)

## BMABMDR 0.0.0.9054

* Updated output of Shapiro Wilks test to be more specific

* Added covariate levels to summary table in full.laplace_MA_Cov()

## BMABMDR 0.0.0.9053

* Updated x-range for plot_priorQ

* Changed covariate plot to log10-scale for dose 

## BMABMDR 0.0.0.9052

* Updated plots to allow for non-zero first dose

* Updated continuous PREP_DATA functions to allow for non-zero first dose

## BMABMDR 0.0.0.9051

* Fixed plot for quantal covariates (accounting for different dose levels between covariate levels)

* Updated implementation for continuous covariates accounting for different dose levels

## BMABMDR 0.0.0.9050

* Fixed BMD shape for informative prior (should not default to 4)

## BMABMDR 0.0.0.9049

* Fixed covariates plot and summary tables to be on original dose range

## BMABMDR 0.0.0.9048

* Fixed dependency issue in anydoseresponse for continuous data

## BMABMDR 0.0.0.9047

* Updated anydoseresponseQ() to use the brms package --> 'use.mcmc' argument is gone now

* Model-averaged posterior is truncated AFTER calculating the credible interval

## BMABMDR 0.0.0.9046

* Updated plotting functions to use correct data (arithmetic vs geometric)

* Model-averaged posterior truncated at maxDose^2; when this happens, a warning is included in the output as 'p.msg'

* Model-specific BMD plot also truncated at maxDose^2

## Update 29/11/2022

* Added 'log' option in plot.BMADR() function (defaults to FALSE, if TRUE the response is plotted on log10 scale)

* Fixed gamlss error for continuous data

## Update 29/09/2022

* Fixed small error in covariate implementation (Laplace function)

## Update 27/09/2022

* LL.R updated to include likelihood of null model for quantal and clustered continuous data

* anydoseresponse.R updated to include option to use Laplace (use.mcmc = FALSE) for clustered quantal and clustered continuous data

## Updates 22/09/2022

* PREP_DATA_N_C and PREP_DATA_LN_C now return 4 additional objects: shapiro.p, shapiro.msg, bartlett.p, bartlett.msg (containing the p-value for shapiro/bartlett test and the corresponding message)

* anydoseresponseQ now works with Laplace for clustered data by specifying 'cluster = TRUE, use.mcmc = FALSE'

## Updates 27/07/2022

* Merged plot_prior.R scripts for clustered and unclustered endpoints

* Added getBMD.R script for covariates

* Added fun_cov_selection.R script for covariates

* Added plot_model_fit.R script for covariates

* Merged output.R scripts for all endpoints

* Merged classes.R scripts for all endpoints

* Merged diagnostics.R scripts for all endpoints

* Merged anydoseresponse.R scripts for all endpoints

* Merged DRM.R scripts for all endpoints

* Merged FUNs.R scripts for all endpoints

* Added NtoLN() and LNtoN() functions to FUNs.R script

* Merged fun_Data.R scripts for all endpoints

* Merged fun_Laplace.R scripts for all endpoints

* Merged fct_optim.R scripts for all endpoints

* Merged LL.R scripts for all endpoints

* Merged fun_modelTest.R scripts for all endpoints

* Merged fct_mcmc.R scripts for all endpoints

* Merged fun_Sampling.R scripts for all endpoints

* Merged plot_prior.R scripts for all endpoints


## General Updates

* Added help files to all the functions 

* Added News.md file

* Added tests and examples

* Added two internal datasets


## Specific Updates

* Added function to get the available models. 

* Updated methods for continuous endpoints in case of negative geometric means (this change affects almost all files)

* Internal datsets are now available for testing the functions in the package

* The available datasets are: immunotoxicity(continuous-increasing), LearningMemory(continuous-decreasing)

* To get the default prior for the maximum response (similar for LN):

    data\_N <- PREP\_DATA\_N(summ.data, sumstats = T, q = q)
    
    min.maxresp <- data\_N\$data\$priorlb[3]*data\_N\$data\$priormu[1]
    
    mode.maxresp <- data\_N\$data\$priormu[3]*data\_N\$data\$priormu[1]
    
    max.maxresp <- data\_N\$data\$priorub[3]*data\_N\$data\$priormu[1]
