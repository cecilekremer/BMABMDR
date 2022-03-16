# EFSA Bayesian Benchmark Dose Response Platform

## General Updates

* Added help files to all the functions 

* Added News.md file

* Added tests and examples

* Updated methods for continuous endpoints in case of negative geometric means (this change affects almost all files)

* 


## Specific Updates

* To get the default prior for the maximum response (similar for LN):

    data\_N <- PREP\_DATA\_N(summ.data, sumstats = T, q = q)
    
    min.maxresp <- data\_N\$data\$priorlb[3]*data\_N\$data\$priormu[1]
    
    mode.maxresp <- data\_N\$data\$priormu[3]*data\_N\$data\$priormu[1]
    
    max.maxresp <- data\_N\$data\$priorub[3]*data\_N\$data\$priormu[1]

*

*

