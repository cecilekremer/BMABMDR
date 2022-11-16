# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################


test_that("Fit models clustered quantal",{

      data <- read.csv(file.path(dataDir, "example_clusteredQuantal.csv"))

      orderedData <- data[,
          c("dose", "y", "n", "litter")
      ]

      anydoseresponseQ(orderedData$dose, orderedData$y, orderedData$n, cluster = T)

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
    modelFit$MA
    modelFit$gof_check

    pFLBMD_Q = plot.BMADRQ(modelFit, weight_type = "LP", include_data = T, all = F, title = '')
    pFLBMD_Q$BMDs
    pFLBMD_Q$weights
    pFLBMD_Q$model_fit
    pFLBMD_Q$MA_fit

    # pFLBMD_Q$model_fit2
    # pFLBMD_Q$MA_fit2

    })
