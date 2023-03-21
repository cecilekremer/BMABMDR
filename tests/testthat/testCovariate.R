# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################

test_that("Summary continuous data with covariate", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data.test <- read.csv(file.path(dataDir,'test_data.csv'), sep = ';')

  summ.data <- data.frame(
    x = data.test$Dose,
    y = data.test$Mean,
    s = data.test$SD,
    n = data.test$N,
    cov = data.test$group
  )
  q = 0.2
  prior.weights = rep(1,16)

  # Fit models
  FLBMD <- full.laplace_MA_Cov(summ.data,
                               sumstats = T,
                               sd = T, # option not used for Quantal data
                               q = q,
                               prior.d = 'N11',
                               extended = F,
                               ndraws = 30000,
                               seed = 123,
                               pvec = c(0.05, 0.5, 0.95),
                               prior.weights = prior.weights)

  expect_type(FLBMD, 'list')
  expect_true(is.numeric(FLBMD$MA))
  expect_true(is.data.frame(FLBMD$summary))

  # Plots
  pt <- basic.plot(FLBMD, model_name = 'E4_N', increasing = T)
  expect_true(is.ggplot(pt))


})

test_that("Individual continuous data with covariate", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  load(file.path(dataDir, "das1.rda"))
  data.test <- das1$data

  ind.data <- data.frame(
    x = data.test$Dose,
    y = data.test$LDH,
    #          s = data.test$SD,
    #          n = data.test$N,
    cov = data.test$sex
  )
  ind.data <- ind.data[which(!is.na(ind.data$y)),]

  q = 0.05
  prior.weights = rep(1,16)

  # Fit models TODO: returns an error
  FLBMD <- full.laplace_MA_Cov(
    data = ind.data,
    sumstats = FALSE,
    sd = TRUE, # option not used for Quantal data
    q = q,
    prior.d = 'N11'
  )

  expect_type(FLBMD, 'list')
  expect_true(is.data.frame(FLBMD$summary))

  # Plots
  pt <- basic.plot(FLBMD, model_name = 'E4_N', increasing = T)
  expect_true(is.ggplot(pt))

})


test_that("quantal data with covariate", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.csv(file.path(dataDir, "example_quantal.csv"))

  # Each covariate level should be present in each dose group
  data.input <- data.frame(
    dose = rep(data$dose, 2),
    y = rep(data$response, 2),
    n = rep(data$size, 2),
    covariate = c(rep('Male', 5), rep('Female', 5))
  )
  # anydoseresponseQ(data.input$dose, data.input$y, data.input$n, cluster = FALSE, use.mcmc = FALSE)


  q = 0.1

  # Fit models TODO: returns an error
  modelFit <- full.laplace_MA_Q_Cov(
    data = data.input,
    q = q
  )

  expect_true(is.data.frame(modelFit$summary))

  # Plots
  pt <- basic.plotQ(modelFit, model_name = 'E4_Q')
  expect_true(is.ggplot(pt))

  # Fit without covariate
  data.input.Q <- PREP_DATA_QA(data = data.input, q = q, sumstats = T)
  modelFit <- full.laplaceQ_MA(data.input.Q, prior.weights = rep(1,8), pvec = c(0.05,0.5,0.95))
  expect_s3_class(modelFit, 'BMADRQ')
  # modelFit$MA
  # modelFit$gof_check

})

