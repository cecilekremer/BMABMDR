
# Continuous data
dose = c(0,6.25,12.5,25,50,100)
mean = c(10.87143,10.16669,10.81050,10.41179,12.38305,18.47681)
sd = c(1.804554,1.805939,3.858265,1.626007,2.045695,2.322449)
n = rep(10,6)

summ.data = data.frame(x = dose, y = mean, s = sd, n = n)


test_that("Dose-response effect", {
    
    # normal distribution
    expect_is(anydoseresponseN(summ.data$x, summ.data$y, summ.data$s, summ.data$n),
      "list")
    # lognormal distribution
    expect_is(anydoseresponseLN(summ.data$x, summ.data$y, summ.data$s, summ.data$n),
      "list")
    
  })


context("Analysis continuous")


# sampling specification
ndr=30000
nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

# prior model weights
prior.weights = c(rep(1,8), rep(1,8))

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


test_that("Laplace approximation", {
    
    FLBMD=full.laplace_MA(data_N,
      data_LN,
      prior.weights,
      ndraws=ndr,
      seed=123,
      pvec=pvec,
      plot=F)
    
    # MA estimates
    expect_equal(FLBMD$MA, c(BMDL = 20.95212, BMD = 36.93779, BMDU = 54.18734), tolerance = 1E-05)
    # model weights
    expect_is(FLBMD$weights, "numeric")
    # model-specific fit
    expect_is(FLBMD$E4_N, "matrix")
    # test whether best-fitting model fits wel (BF < 10 means equally well as saturated model; BF > 10 means best fit is better than saturated model)
    expect_equal(FLBMD$bf, 41.60845)
    
    # plot output
    pFLBMD = plot.BMADR(FLBMD, weight_type = "LP", include_data = T, all = F, title = '')
    expect_is(pFLBMD, "list")
    expect_equal(names(pFLBMD), c("BMDs", "weights", "model_fit_N", "model_fit_LN", "model_fit", "MA_fit"))

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
    
  })

test_that("Sampling", {
    
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
    
  })
