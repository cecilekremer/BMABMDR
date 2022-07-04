# EFSA Bayesian Benchmark Dose Response Platform

## General Updates

* Added help files to all the functions 

* Added News.md file

* Added tests and examples

* Added two internal datasets

## Updates 25/6/2022

* Added functions and stan models for clustered continuous data (modelTest not yet implemented)

* Fixed anydoseresponseQ() function

* Fixed problem with plots

* Added option for clustered data to plot.BMADR()

## Updates 2/6/2022

* Updated Stan models for continuous and quantal endpoints

* FUNs.R : change in flat() function input

* fun_DataN.R : change in PREP_DATA_N() function input + one additional element in output list

* fun_DataLN.R : change in PREP_DATA_LN() function input + one additional element in output list

* fun_DataQ.R : change in PREP_DATA_QA() function input + one additional element in output list

* fun_Laplace.R : only internal changes

* fun_LaplaceQ.R : internal changes + added included models to output

* fct_optim.R : only internal changes

* fct_optimQ.R : only internal changes

* fun_Sampling.R : only internal changes

* fun_SamplingQ.R : only internal changes + added included models to output

* fct_mcmc.R : only internal changes

* fct_mcmcQ.R : only internal changes

* added progress bar increase for Laplace and Sampling method


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
    
