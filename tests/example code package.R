rm(list=ls())

# install package from zip file
# install.packages("~/GitHub/BMABMDR_0.0.0.9011.tar.gz", repos = NULL, type = "source")

library(BMABMDR)
library(gamlss)
# library(posterior)
# library(RColorBrewer)
# library(ggpubr)

## available models?
get_models('continuous')
get_models('quantal')

# sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
#
# Matrix products: default
#
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
# [4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252
#
# attached base packages:
#   [1] parallel  splines   stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] gamlss_5.3-4       nlme_3.1-152       gamlss.dist_5.3-2  MASS_7.3-54        gamlss.data_6.0-1  BMABMDR_0.0.0.9011
#
# loaded via a namespace (and not attached):
#   [1] matrixStats_0.61.0   RColorBrewer_1.1-2   rstan_2.21.2         numDeriv_2016.8-1.1  tensorA_0.36.2       tools_4.1.0          backports_1.4.0
# [8] utf8_1.2.2           R6_2.5.1             DBI_1.1.1            colorspace_2.0-2     raster_3.4-13        withr_2.4.3          sp_1.4-5
# [15] tidyselect_1.1.1     gridExtra_2.3        prettyunits_1.1.1    processx_3.5.2       Brobdingnag_1.2-6    curl_4.3.2           compiler_4.1.0
# [22] cli_3.0.1            unmarked_1.1.1       posterior_1.0.1      scales_1.1.1         checkmate_2.0.0      mvtnorm_1.1-2        AICcmodavg_2.3-1
# [29] mc2d_0.1-21          callr_3.7.0          stringr_1.4.0        StanHeaders_2.21.0-7 foreign_0.8-81       rio_0.5.27           pkgconfig_2.0.3
# [36] bbmle_1.0.24         rlang_0.4.12         readxl_1.3.1         VGAM_1.1-5           generics_0.1.0       farver_2.1.0         jsonlite_1.7.2
# [43] dplyr_1.0.7          zip_2.2.0            car_3.0-11           distributional_0.2.2 inline_0.3.19        magrittr_2.0.1       loo_2.4.1
# [50] Matrix_1.3-3         Rcpp_1.0.7           munsell_0.5.0        fansi_0.5.0          abind_1.4-5          lifecycle_1.0.1      stringi_1.7.4
# [57] carData_3.0-4        pkgbuild_1.3.1       plyr_1.8.6           grid_4.1.0           bdsmatrix_1.3-4      forcats_0.5.1        crayon_1.4.2
# [64] lattice_0.20-44      haven_2.4.3          hms_1.1.0            ps_1.6.0             pillar_1.6.4         ggpubr_0.4.0         ggsignif_0.6.2
# [71] codetools_0.2-18     stats4_4.1.0         rstantools_2.1.1     glue_1.4.2           V8_3.6.0             data.table_1.14.0    RcppParallel_5.1.4
# [78] vctrs_0.3.8          cellranger_1.1.0     gtable_0.3.0         purrr_0.3.4          tidyr_1.1.3          assertthat_0.2.1     ggplot2_3.3.5
# [85] openxlsx_4.2.4       xtable_1.8-4         broom_0.7.9          coda_0.19-4          rstatix_0.7.0        survival_3.2-11      tibble_3.1.4
# [92] ellipsis_0.3.2       bridgesampling_1.1-2

#######################
### CONTINUOUS DATA ###

dose = c(0,6.25,12.5,25,50,100)
mean = c(10.87143,10.16669,10.81050,10.41179,12.38305,18.47681)
sd = c(1.804554,1.805939,3.858265,1.626007,2.045695,2.322449)
n = rep(10,6)

summ.data = data.frame(x = dose, y = mean, s = sd, n = n)
plot(summ.data$x, summ.data$y, type = 'l')

# Test for dose-response effect
anydoseresponseN(summ.data$x, summ.data$y, summ.data$s, summ.data$n) # normal distribution
anydoseresponseLN(summ.data$x, summ.data$y, summ.data$s, summ.data$n) # lognormal distribution

### ANALYSIS ###

# sampling specification
ndr=30000
nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

# prior model weights
prior.weights = c(rep(1,4), rep(1,4), rep(0,8))


# bmr
q = 0.1

pvec = c(0.05,0.5,0.95)


## Preparing data

# uninformative
data_N = PREP_DATA_N(summ.data,
                     sumstats = T,
                     q = q)
data_LN = PREP_DATA_LN(summ.data,
                       sumstats = T,
                       q = q)

# informative, for example:
# data_N = PREP_DATA_N(summ.data, sumstats = T, q = q, bkg = c(8,10.58,12), maxy = c(18,20.28,22), prior.BMD = c(20,40,60), shape.BMD = 4)
# data_LN = PREP_DATA_LN(summ.data, sumstats = T, q = q, bkg = c(8,10.58,12), maxy = c(18,20.28,22), prior.BMD = c(20,40,60), shape.BMD = 4)


## Laplace approximation

FLBMD=full.laplace_MA(data_N,
                      data_LN,
                      prior.weights,
                      ndraws=ndr,
                      seed=123,
                      pvec=pvec,
                      plot=F)
# MA estimates
FLBMD$MA
# model weights
round(FLBMD$weights,4)
# model-specific fit
FLBMD$E4_N
# test whether best-fitting model fits wel (BF < 10 means equally well as saturated model; BF > 10 means best fit is better than saturated model)
FLBMD$bf

# output as dataframe/list
BMDWeights(FLBMD, 'continuous')
summary.BMADR(FLBMD)

# plot output
pFLBMD = plot.BMADR(FLBMD, weight_type = "LP", include_data = T, all = F, title = '')
pFLBMD$BMDs
pFLBMD$weights
pFLBMD$model_fit_N
pFLBMD$model_fit_LN
pFLBMD$model_fit
pFLBMD$MA_fit

# plot prior vs posterior
plot_prior(FLBMD, data_N$data, "E4_N", parms = T)
plot_prior(FLBMD, data_N$data, "E4_N", parms = F)
plot_prior(FLBMD, data_N$data, "P4_N", parms = T)
plot_prior(FLBMD, data_LN$data, "L4_LN", parms = T)


## Sampling

SBMD = sampling_MA(data_N, data_LN,
                   prior.weights,
                   ndraws=ndr, nrchains=nrch,
                   nriterations=nriter, warmup=wu, delta=dl,
                   treedepth=trd, seed=sd, pvec=pvec,
                   plot=F)
# MA estimates
SBMD$MA_bridge_sampling
SBMD$MA_laplace
# convergence & divergence
SBMD$convergence
SBMD$divergences*100 # percentage of iterations that were divergent
# model-specific fit
SBMD$E4_N
# test whether best-fitting model fits wel (BF < 10 means equally well as saturated model; BF > 10 means best fit is better than saturated model)
SBMD$bf

# output as dataframe/list
BMDWeights(SBMD, 'continuous')
summary.BMADR(SBMD)

# plot output
pSBMD = plot.BMADR(SBMD, weight_type = "BS", include_data = T, all = F, title = '')
pSBMD$BMDs
pSBMD$weights
pSBMD$model_fit_N
pSBMD$model_fit_LN
pSBMD$model_fit
pSBMD$MA_fit

# plot prior vs posterior
plot_prior(SBMD, data_N$data, "E4_N", parms = T)
plot_prior(SBMD, data_LN$data, "P4_LN", parms = T)



####################
### QUANTAL DATA ###

dose = c(0, 5, 15, 50, 100)
y = c(0, 4, 6, 5, 12)
n = c(20, 20, 20, 20, 20)

summ.data = data.frame(x = dose, y = y, n = n)

# Test for dose-response effect (not optimal yet)
# anydoseresponseQ(summ.data$x, summ.data$y, summ.data$n)

### ANALYSIS ###

# sampling specification
ndr=30000
nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

# prior model weights
prior.weights = rep(1,8)

# bmr
q = 0.1

pvec = c(0.05,0.5,0.95)


## Preparing data

# uninformative
data_Q = PREP_DATA_QA(summ.data,
                      sumstats = T,
                      q = q)

# informative, for example:
# data_Q = PREP_DATA_QA(summ.data, sumstats = T, q = q, prior.BMD = c(10,16,25), shape.BMD = 4)


## Laplace approximation

FLBMD_Q = full.laplaceQ_MA(data_Q,
                           prior.weights,
                           ndraws=ndr,
                           seed=123,
                           pvec=pvec)
# MA estimates
FLBMD_Q$MA
# model weights
round(FLBMD_Q$weights,4)
# model-specific fit
FLBMD_Q$E4_Q
# test whether best-fitting model fits wel (BF < 10 means equally well as saturated model; BF > 10 means best fit is better than saturated model)
# FLBMD_Q$bf

# output as dataframe/list
BMDWeights(FLBMD_Q, 'quantal')
summary.BMADRQ(FLBMD_Q)

# plot output
pFLBMD_Q = plot.BMADRQ(FLBMD_Q, weight_type = "LP", include_data = T, all = F, title = '')
pFLBMD_Q$BMDs
pFLBMD_Q$weights
pFLBMD_Q$model_fit
pFLBMD_Q$MA_fit

# plot prior vs posterior
plot_priorQ(FLBMD_Q, data_Q$data, "E4_Q")
plot_priorQ(FLBMD_Q, data_Q$data, "P4_Q")
plot_priorQ(FLBMD_Q, data_Q$data, "L4_Q")


## Sampling

SBMD_Q = samplingQ_MA(data_Q,
                      prior.weights,
                      ndraws=ndr, nrchains=nrch,
                      nriterations=nriter, warmup=wu, delta=dl,
                      treedepth=trd, seed=sd, pvec=pvec)
# MA estimates
SBMD_Q$MA_bridge_sampling
SBMD_Q$MA_laplace
# convergence & divergence
SBMD_Q$convergence
SBMD_Q$divergences*100 # percentage of iterations that were divergent

# output as dataframe/list
BMDWeights(SBMD_Q, 'quantal')
summary.BMADRQ(SBMD_Q)

# plot output
pSBMD_Q = plot.BMADRQ(SBMD_Q, weight_type = "BS", include_data = T, all = F, title = '')
pSBMD_Q$BMDs
pSBMD_Q$weights
pSBMD_Q$model_fit
pSBMD_Q$MA_fit

# plot prior vs posterior
plot_priorQ(SBMD_Q, data_Q$data, "E4_Q")

# clustered Quantal Data

clusterdata <- data.frame(
  dose = c(rep(c(0, 0.2, 0.6, 6, 60, 120), c(26, 22, 24, 26, 19, 16))),
  y = c(1,1,0,0,0,0,0,2,1,3,0,0,0,0,0,0, 0, 1,0,0,0,0,
        2,0,1,0,3, 0,3,0,1,0,1,1,0,4,0,0,1, 0,  0,  0,  0,  2,  0,  0,  0,
        0,  0,  0,  4,  0,  4,  0,  1,  0,  0,  0,  0,1,  3, 1,  0,  0,  5,
        0,  0,  1,  5,  0,  1,  0,  0,  6,  0,  1,  1,  0, 2,  0,  0,  0,  1,
        11,  0,  5,  5,  0,  3,  0,  0,  3,  0,  0,  0,  2,  5,  3,  6,  1,  1, 10,  3,
        4,  2,  3,  2,  2,  2,  9,  3,  2,  1,  0,  1,  3,  2, 3,  3,  7,  0,
        5,  3,  2,  6,  3,  1,  2,  6,  3,  3, 2,  2),
  n = c(12,12,12,13,1, 13,10,14, 12, 12, 14, 13, 12, 11, 15, 14, 11, 11, 12, 11, 14, 12, 13, 13, 14, 13,
        8, 14, 10, 13, 14, 12, 10, 10, 11, 10, 12, 12, 14, 14,  9,  8, 12,  9, 13,  9, 13, 12, 12, 13, 10,
        13,  9, 11,  9, 10, 12, 14,  5,  9,  9,  8,  7, 13, 14, 12, 13, 15, 11, 11,  7, 14,  8,  9, 12,
        6,  9,  9, 13, 10,  4,  8, 10, 11,  7, 10, 10,  7, 12,  7, 11, 14,  4,  7,  6,  2,  5,  9,  6,  1,
        1, 10,  3,  4,  3,  5,  2,  4,  2,  9,  3,  6,  1,  4,  1, 10,  2, 3,  3,  7,  1,  5,  3,  2,  6,
        8,  1,  2,  6,  3,  3,  2, 2),
  liter = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1)

)

clusterdataQ <- PREP_DATA_QA(data = clusterdata, sumstats = TRUE,
                             q = 0.1, bkg = NULL, shape.a = 4, shape.BMD = 0.0001,
                             cluster = TRUE)
pw <- rep(1, 8)

### Laplace Approximation
testbb_laplace <- full.laplaceQ_MA(data.Q = clusterdataQ, prior.weights = pw)

# MA estimates
testbb_laplace$MA
# model weights
round(testbb_laplace$weights,4)
# model-specific fit
testbb_laplace$E4_Q

# output as dataframe/list
BMDWeights(testbb_laplace, 'quantal')
summary.BMADRQ(testbb_laplace)

# plot output
pSBMD_CQ <- plot(testbb_laplace, weight_type = 'LP', title = 'Quantal (Full Laplace)', all = F)
pSBMD_CQ$BMDs
pSBMD_CQ$weights
pSBMD_CQ$model_fit
pSBMD_CQ$MA_fit

# plot prior vs posterior
plot_priorQ(testbb_laplace, data = clusterdataQ$data,
            model_name = "QE4_Q")

### Sampling

testbb_sampling <- samplingQ_MA(data.Q = clusterdataQ, prior.weights = pw)

testbb_sampling$MA_bridge_sampling
testbb_sampling$MA_laplace
# convergence & divergence
testbb_sampling$convergence
testbb_sampling$divergences*100 # percentage of iterations that were divergent

# output as dataframe/list
BMDWeights(testbb_sampling, 'quantal')
summary.BMADRQ(testbb_sampling)


# plot output
pFBMD_CQ <- plot(testbb_sampling, weight_type = 'LP', title = 'Quantal (Sampling)', all = F)
pFBMD_CQ$BMDs
pFBMD_CQ$weights
pFBMD_CQ$model_fit
pFBMD_CQ$MA_fit




