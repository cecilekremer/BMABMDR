

test_that("Prepare data (individual continuous)", {


  # Read in data
  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.csv(file.path(dataDir, "EggShellThickness.csv"))

  summ.data <- data.frame(
    x = data[,"Exposure"],
    y = data[, "Eggshell.Thickness"]
  )

  # Prepare data
  q = 0.05
  sumstats = FALSE
  sd = TRUE
  extended = TRUE  #TRUE

  # Uninformative
  data_N = PREP_DATA_N(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    extended = extended
  )

  data_LN = PREP_DATA_LN(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    extended = extended
  )

  expect_type(data_N, 'list')
  expect_type(data_LN, 'list')

})

test_that("Prepare data (summary continuous)", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.table(file.path(dataDir, "test_data_cont.txt"), header = T, sep = "\t", dec = ".")

  summ.data <- data.frame(
    x = data[, "doseN"],
    y = data[, "Response"],
    s = data[, "Stdev"],
    n = data[, "N"]
  )

  ADRN <- anydoseresponseN(summ.data$x, summ.data$y, summ.data$s, summ.data$n)
  ADRLN <- anydoseresponseLN(summ.data$x, summ.data$y, summ.data$s, summ.data$n)

  expect_type(ADRN, 'list')
  expect_type(ADRLN, 'list')

  ## Preparing data
  # bmr
  q = 0.05
  sumstats = TRUE
  sd = TRUE
  extended = TRUE

  # Uninformative
  data_N = PREP_DATA_N(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'N11',
    extended = extended
  )

  data_LN = PREP_DATA_LN(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'N11',
    extended = extended
  )

  # data_N = PREP_DATA_N(summ.data, sumstats = T, q = q, sd = T, prior.d = 'N11', extended = T)#,

  expect_type(data_N, 'list')
  expect_type(data_LN, 'list')

  # informative, for example:
  # data_N = PREP_DATA_N(summ.data, sumstats = T, q = q, bkg = c(8,10.58,12), maxy = c(18,20.28,22), prior.BMD = c(20,40,60), shape.BMD = 4)
  # data_LN = PREP_DATA_LN(summ.data, sumstats = T, q = q, bkg = c(8,10.58,12), maxy = c(18,20.28,22), prior.BMD = c(20,40,60), shape.BMD = 4)

})

test_that("Laplace approximation summary cont data", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.table(file.path(dataDir, "test_data_cont.txt"), header = T, sep = "\t", dec = ".")

  summ.data <- data.frame(
    x = data[, "doseN"],
    y = data[, "Response"],
    s = data[, "Stdev"],
    n = data[, "N"]
  )

  # summ.data$s = min(summ.data$s)

  ## Preparing data
  # bmr
  q = 0.05
  sumstats = TRUE
  sd = TRUE
  extended = TRUE

  # Uninformative
  data_N = PREP_DATA_N(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'EPA',
    maxy = c(12,12.1,12.5),
    extended = extended
  )

  data_LN = PREP_DATA_LN(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'EPA',
    maxy = c(12,12.1,12.5),
    extended = extended
  )

  expect_type(data_N, 'list')
  expect_type(data_LN, 'list')

  # sampling specification
  ndr=30000
  nrch=3
  nriter=3000
  wu=1000
  dl=0.8
  trd=10
  sd=123

  # prior model weights
  prior.weights = c(rep(1,8), rep(1,8))

  pvec = c(0.05,0.5,0.95)

  # Fit models
  FLBMD=full.laplace_MA(
    data_N,
    data_LN,
    prior.weights,
    ndraws=ndr,
    seed=123,
    pvec=pvec,
    plot=F)


  # MA estimates
  expect_s3_class(FLBMD, 'BMADR')
  expect_true(is.numeric(FLBMD$MA))
  # model weights
  expect_true(is.numeric(FLBMD$weights))
  # model-specific fit
  expect_true(is.matrix(FLBMD$E4_N))

  # plot output
  pFLBMD = plot.BMADR(FLBMD, weight_type = "LP", include_data = T, all = F, title = '', log = F)
  expect_type(pFLBMD, "list")
  expect_equal(names(pFLBMD), c("BMDs", "weights", "model_fit_N", "model_fit_LN", "model_fit", "MA_fit"))

  # plot prior vs posterior
  expect_true(is.ggplot(plot_prior(FLBMD, data_N$data, "E4_N", parms = T)))

})


test_that("Dose-response effect continuous", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.table(file.path(dataDir, "test_data_cont.txt"), header = T, sep = "\t", dec = ".")

  summ.data <- data.frame(
    x = data[, "doseN"],
    y = data[, "Response"],
    s = data[, "Stdev"],
    n = data[, "N"]
  )

  adr <- anydoseresponseN(summ.data$x, summ.data$y, summ.data$s, summ.data$n)
  adrln <- anydoseresponseLN(summ.data$x, summ.data$y, summ.data$s, summ.data$n)

  expect_type(adr, 'list')
  expect_type(adrln, 'list')

})


test_that("Sampling continuous data", {

  dataDir <- '~/GitHub/BMABMDR/data/'
  data <- read.table(file.path(dataDir, "test_data_cont.txt"), header = T, sep = "\t", dec = ".")

  summ.data <- data.frame(
    x = data[, "doseN"],
    y = data[, "Response"],
    s = data[, "Stdev"],
    n = data[, "N"]
  )

  summ.data$s <- min(summ.data$s)

  # sampling specification
  ndr=30000
  nrch=3
  nriter=3000
  wu=1000
  dl=0.8
  trd=10
  sd=123

  # prior model weights
  prior.weights = c(rep(1,8), rep(0,8))

  pvec = c(0.05,0.5,0.95)

  ## Preparing data
  # bmr
  q = 0.05
  sumstats = TRUE
  sd = TRUE
  extended = TRUE

  # Uninformative
  data_N = PREP_DATA_N(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'EPA',
    extended = extended
  )

  data_LN = PREP_DATA_LN(
    summ.data,
    sumstats = sumstats,
    q = q,
    sd = sd,
    prior.d = 'EPA',
    extended = extended
  )

  expect_type(data_N, 'list')
  expect_type(data_LN, 'list')

  SBMD = sampling_MA(data_N, data_LN,
                     prior.weights,
                     ndraws=ndr, nrchains=nrch,
                     nriterations=nriter, warmup=wu, delta=dl,
                     treedepth=trd, seed=sd, pvec=pvec)

  expect_s3_class(SBMD, 'BMADR')

  # MA estimates
  expect_true(is.numeric(SBMD$MA_bridge_sampling))
  # convergence & divergence
  # SBMD$convergence
  # SBMD$divergences*100 # percentage of iterations that were divergent
  # model-specific fit
  expect_true(is.matrix(SBMD$E4_N))

  # plot output
  pSBMD = plot.BMADR(SBMD, weight_type = "LP", include_data = T, all = F, title = '', log = F)
  expect_type(pSBMD, 'list')
  expect_equal(names(pSBMD), c("BMDs", "weights", "model_fit_N", "model_fit", "MA_fit"))

  # plot prior vs posterior
  expect_true(is.ggplot(plot_prior(SBMD, data_N$data, "E4_N", parms = T)))

})


