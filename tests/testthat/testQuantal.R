


# test_that("Dose-response effect quantal", {
#
#   # Quantal data
#   dose = c(0, 5, 15, 50, 100)
#   y = c(0, 4, 6, 5, 12)
#   n = c(20, 20, 20, 20, 20)
#
#   summ.data = data.frame(x = dose, y = y, n = n)
#
#   adr <- anydoseresponseQ(summ.data$x, summ.data$y, summ.data$n)
#   expect_type(adr, 'list')
#
# })


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
  expect_type(data_Q, 'list')

})


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

  # MA estimates
  expect_s3_class(FLBMD_Q, 'BMADRQ')
  expect_true(is.numeric(FLBMD_Q$MA))
  # model weights
  expect_true(is.numeric(FLBMD_Q$weights))
  # model-specific fit
  expect_true(is.matrix(FLBMD_Q$E4_Q))

  # plot output
  pFLBMD_Q = plot.BMADRQ(FLBMD_Q, weight_type = "LP", include_data = T, all = F, title = '')
  expect_type(pFLBMD_Q, 'list')
  expect_equal(names(pFLBMD_Q), c("BMDs", "weights", "model_fit", "MA_fit"))

  # plot prior vs posterior
  expect_true(is.ggplot(plot_priorQ(FLBMD_Q, data_Q$data, "E4_Q")))

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

  expect_s3_class(SBMD_Q, 'BMADRQ')

  # MA estimates
  expect_true(is.numeric(SBMD_Q$MA_bridge_sampling))
  # model-specific fit
  expect_true(is.matrix(SBMD_Q$E4_Q))
  # plot output
  pSBMD_Q = plot.BMADRQ(SBMD_Q, weight_type = "LP", include_data = T, all = F, title = '')
  expect_type(pSBMD_Q, 'list')
  expect_equal(names(pSBMD_Q), c("BMDs", "weights", "model_fit", "MA_fit"))

  # plot prior vs posterior
  expect_true(is.ggplot(plot_priorQ(SBMD_Q, data_Q$data, "E4_Q")))

})

