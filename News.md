# EFSA Bayesian Benchmark Dose Response Platform

## Update 29/11/2022

* Added 'log' option in plot.BMADR() function (defaults to TRUE to plot response on log10 scale)

* Fixed gamlss error for continuous data

## Update 29/09/2022

* Fixed small error in covariate implementation (Laplace function)

## Update 27/09/2022

* LL.R updated to include likelihood of null model for quantal and clustered continuous data

* anydoseresponse.R updated to include option to use Laplace (use.mcmc = TRUE) for clustered quantal and clustered continuous data

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
