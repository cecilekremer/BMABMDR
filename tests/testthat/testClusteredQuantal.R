# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################

test_that("Fit models clustered quantal",{

  dataDir <- '~/GitHub/BMABMDR/data/'

  data <- read.csv(file.path(dataDir, "example_clusteredQuantal.csv"))

  orderedData <- data[,
                      c("dose", "y", "n", "litter")
  ]

  # anydoseresponseQ(orderedData$dose, orderedData$y, orderedData$n, cluster = T)

  argListPrep <- list(
    data = orderedData,
    q = 0.1,
    extended = TRUE,
    cluster = TRUE
  )

  prepData <- do.call("PREP_DATA_QA", argListPrep)

  # Fit models
  argListFit <- list(
    data.Q = prepData,
    prior.weights = rep(1,8)
  )

  modelFit <- do.call("full.laplaceQ_MA", argListFit)
  expect_s3_class(modelFit, 'BMADRQ')

  pFLBMD_Q = plot.BMADRQ(modelFit, weight_type = "LP", include_data = T, all = F, title = '')
  expect_type(pFLBMD_Q, 'list')
  expect_equal(names(pFLBMD_Q), c("BMDs", "weights","model_fit2", "model_fit", "MA_fit2", "MA_fit"))


})
