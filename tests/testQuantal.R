

dataDir <- '~/GitHub/BMABMDR/data/'


test_that("Dose-response effect quantal", {

  # Quantal data
  dose = c(0, 5, 15, 50, 100)
  y = c(0, 4, 6, 5, 12)
  n = c(20, 20, 20, 20, 20)

  summ.data = data.frame(x = dose, y = y, n = n)

  anydoseresponseQ(summ.data$x, summ.data$y, summ.data$n)

})


test_that("Data prep quantal", {

  # Quantal data
  dose = c(0, 5, 15, 50, 100)
  y = c(0, 4, 6, 5, 12)
  n = c(20, 20, 20, 20, 20)

  summ.data = data.frame(x = dose, y = y, n = n)

  # sampling specification
  ndr=30000
  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

  # prior model weights
  prior.weights = rep(1,8)

  # bmr
  q = 0.1

  pvec = c(0.05,0.5,0.95)

  # uninformative
  data_Q = PREP_DATA_QA(summ.data,
                        sumstats = T,
                        q = q)

})


# informative, for example:
# data_Q = PREP_DATA_QA(summ.data, sumstats = T, q = q, prior.BMD = c(10,16,25), shape.BMD = 4)


test_that("Laplace approximation quantal", {

  # Quantal data
  dose = c(0, 5, 15, 50, 100)
  y = c(0, 4, 6, 5, 12)
  n = c(20, 20, 20, 20, 20)

  summ.data = data.frame(x = dose, y = y, n = n)

  # sampling specification
  ndr=30000
  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

  # prior model weights
  prior.weights = rep(1,8)

  # bmr
  q = 0.1

  pvec = c(0.05,0.5,0.95)

  # uninformative
  data_Q = PREP_DATA_QA(summ.data,
                        sumstats = T,
                        q = q,
                        extended = T)

  FLBMD_Q = full.laplaceQ_MA(data_Q,
                             prior.weights,
                             ndraws=ndr,
                             seed=123,
                             pvec=pvec)

  FLBMD_Q$w.msg
  FLBMD_Q$p.msg
  # MA estimates
  FLBMD_Q$MA
  # model weights
  round(FLBMD_Q$weights,4)
  # model-specific fit
  FLBMD_Q$E4_Q
  # test whether best-fitting model fits wel (BF < 10 means equally well as saturated model; BF > 10 means best fit is better than saturated model)
  FLBMD_Q$bf
  FLBMD_Q$gof_check

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

})

test_that("Sampling quantal", {

  # Quantal data
  dose = c(0, 5, 15, 50, 100)
  y = c(0, 4, 6, 5, 12)
  n = c(20, 20, 20, 20, 20)

  summ.data = data.frame(x = dose, y = y, n = n)

  # sampling specification
  ndr=30000
  nrch=3;nriter=3000;wu=1000;dl=0.8;trd=10;sd=123

  # prior model weights
  prior.weights = rep(1,8)

  # bmr
  q = 0.1

  pvec = c(0.05,0.5,0.95)

  # uninformative
  data_Q = PREP_DATA_QA(summ.data,
                        sumstats = T,
                        q = q,
                        extended = T)

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

  # plot output
  pSBMD_Q = plot.BMADRQ(SBMD_Q, weight_type = "LP", include_data = T, all = F, title = '')
  pSBMD_Q$BMDs
  pSBMD_Q$weights
  pSBMD_Q$model_fit
  pSBMD_Q$MA_fit

  # plot prior vs posterior
  plot_priorQ(SBMD_Q, data_Q$data, "E4_Q")

})

