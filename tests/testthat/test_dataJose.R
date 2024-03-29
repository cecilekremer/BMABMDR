# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################


test_that("Logit and probit are fit",{

  dataDir <- '~/GitHub/BMABMDR/data/'

  # Load data
  summ.data <- read.csv(file.path(dataDir,'examplecontinuouswithcovariate.csv'))

  # Subset data
  summ.data <- summ.data[summ.data$Covariate == "Male", c("Dose", "Response3", "SD", "N")]

  # Argument list
  argList <- list(
    data = summ.data,
    sumstats = TRUE,
    sd = TRUE,
    q = 0.1
  )

  # Bartlett test
  ## normal distribution
  b.Normal <- bartlett(sd = summ.data[,3], n = summ.data[,4])

  ## lognormal distribution
  ### Convert arithmetic to geometric summary stats
  summ.data.LN <- data.frame(x = summ.data[,1],
                             y = NtoLN(summ.data[,2], summ.data[,3])[1:length(summ.data[,1])],
                             s = NtoLN(summ.data[,2], summ.data[,3])[(length(summ.data[,1])+1):(2*length(summ.data[,1]))],
                             n = summ.data[,4]
  )
  b.Lognormal <- bartlett(sd = log(summ.data.LN[,3]), n = summ.data.LN[,4])

  ########## Bartlett is rejected for logarithmic scale

  # Data for default analyses
  data_N = do.call("PREP_DATA_N", argList)
  data_LN = do.call("PREP_DATA_LN", argList)

  ## Data for sensitivity analysis
  summ.data.LN.sens1 <- summ.data.LN
  summ.data.LN.sens1$s <- min(summ.data.LN.sens1$s) # minimum geometric sd
  argList$data <- summ.data.LN.sens1
  argList <- append(argList, list(geom.stats = T)) # specify that geometric summary stats are given as input

  data_LN_sens1 <- do.call("PREP_DATA_LN", argList)

  summ.data.LN.sens2 <- summ.data.LN
  summ.data.LN.sens2$s <- max(summ.data.LN.sens2$s) # maximum geometric sd
  argList$data <- summ.data.LN.sens2

  data_LN_sens2 <- do.call("PREP_DATA_LN", argList)

  analysisData <- list(data = list(
    default = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
    defaultAdaptedWeights = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(0,8))),  # set weights to 0 for Lognormal
    minVar = list(data_N = data_N, data_LN = data_LN_sens1, prior.weights = c(rep(1,8), rep(1,8))),
    maxVar = list(data_N = data_N, data_LN = data_LN_sens2, prior.weights = c(rep(1,8), rep(1,8)))

  ), warning = NULL
  )

  # Fit models
  argListFit <- list(
    data.N = analysisData$data$minVar$data_N,
    data.LN = analysisData$data$minVar$data_LN,
    prior.weights = c(rep(0,8), rep(1,8)),
    ndraws = 30000
  )

  fittedModels <- do.call("full.laplace_MA", argListFit)

  modelSummary <- summary.BMADR(fittedModels, type = "continuous")

  # modelSummary$BMDWeights

})



# test_that("Test iterating full.laplace_MA", {
#
#   dataDir <- '~/GitHub/BMABMDR/data/'
#
#   load(file.path(dataDir, "argList_fail.RData"))
#   argListFail <- argList
#   rm(argList)
#   load(file.path(dataDir, "argList_pass.RData"))
#   argListPass <- argList
#
#   identical(argListFail, argListPass)
#   expect_equal(argListFail, argListPass)
#
#   for (i in 1:10) {
#     modelFit <- do.call("full.laplace_MA",
#                         argListPass[!names(argListPass) %in% c("nrchains", "nriterations", "warmup")])
#   }
#
# })


test_that("Bridge Sampling works on Response3", {

  dataDir <- '~/GitHub/BMABMDR/data/'


  # Load data
  summ.data <- read.csv(file.path(dataDir,'examplecontinuouswithcovariate.csv'))

  # Subset data
  summ.data <- summ.data[summ.data$Covariate == "Male", c("Dose", "Response3", "SD", "N")]

  # Argument list
  argList <- list(
    data = summ.data,
    sumstats = TRUE,
    sd = TRUE,
    q = 0.1
  )

  # Bartlett test
  ## normal distribution
  b.Normal <- bartlett(sd = summ.data[,3], n = summ.data[,4])

  ## lognormal distribution
  ### Convert arithmetic to geometric summary stats
  summ.data.LN <- data.frame(x = summ.data[,1],
                             y = NtoLN(summ.data[,2], summ.data[,3])[1:length(summ.data[,1])],
                             s = NtoLN(summ.data[,2], summ.data[,3])[(length(summ.data[,1])+1):(2*length(summ.data[,1]))],
                             n = summ.data[,4]
  )
  b.Lognormal <- bartlett(sd = log(summ.data.LN[,3]), n = summ.data.LN[,4])

  ########## Bartlett is rejected for logarithmic scale

  # Data for default analyses
  data_N = do.call("PREP_DATA_N", argList)
  data_LN = do.call("PREP_DATA_LN", argList)

  ## Data for sensitivity analysis
  summ.data.LN.sens1 <- summ.data.LN
  summ.data.LN.sens1$s <- min(summ.data.LN.sens1$s) # minimum geometric sd
  argList$data <- summ.data.LN.sens1
  argList <- append(argList, list(geom.stats = T)) # specify that geometric summary stats are given as input

  data_LN_sens1 <- do.call("PREP_DATA_LN", argList)

  anydoseresponseLN(summ.data.LN.sens1$x, summ.data.LN.sens1$y, summ.data.LN.sens1$s, summ.data.LN.sens1$n)

  summ.data.LN.sens2 <- summ.data.LN
  summ.data.LN.sens2$s <- max(summ.data.LN.sens2$s) # maximum geometric sd
  argList$data <- summ.data.LN.sens2

  data_LN_sens2 <- do.call("PREP_DATA_LN", argList)

  analysisData <- list(data = list(
    default = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
    defaultAdaptedWeights = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(0,8))),  # set weights to 0 for Lognormal
    minVar = list(data_N = data_N, data_LN = data_LN_sens1, prior.weights = c(rep(1,8), rep(1,8))),
    maxVar = list(data_N = data_N, data_LN = data_LN_sens2, prior.weights = c(rep(1,8), rep(1,8)))

  ), warning = NULL
  )

  selectedDistribution <- "Lognormal"
  # Remove defaultAdaptedWeights scenario if Distribution type has been selected and the barlett test is rejected for that distribution
  if (length(selectedDistribution) < 2){

    if(all(selectedDistribution == "Lognormal" && b.Lognormal[2] < 0.05)){

      analysisData$data$defaultAdaptedWeights <- NULL

    } else if(selectedDistribution == "Normal" & b.Normal[2] < 0.05) {

      analysisData$data$defaultAdaptedWeights <- NULL

    }

  }

  # Fit models using BS
  ## Error thrown when i = 2 or 3, meaning that there is a problem when the minimum and maximum variance is used
  for (i in seq_along(analysisData$data)) {
    argListFit <- list(
      data.N = analysisData$data[[i]]$data_N,
      data.LN = analysisData$data[[i]]$data_LN,
      prior.weights = c(rep(0,8), rep(1,8))
    )

    # fittedModels <- do.call("sampling_MA", argListFit)
    fittedModels <- do.call("full.laplace_MA", argListFit)

  }

})











