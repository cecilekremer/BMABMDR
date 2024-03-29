# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################

test_that("LDH response is fitted",{

   dataDir <- '~/GitHub/BMABMDR/data/'
   load(file.path(file.path(dataDir, "das1.rda")))

   indData <- na.omit(das1$data[, c("Dose", "LDH")]) #BW instead of LDH works

   argListPrepData <- list(
      data = indData,
      sumstats = FALSE,
      sd = TRUE,
      q = 0.05,
      extended = TRUE
   )

   data_N <- do.call("PREP_DATA_N", argListPrepData)
   data_LN <-  do.call("PREP_DATA_LN", argListPrepData)

   # Fit models

   argListFit <- list(
      data.N = data_N,
      data.LN = data_LN,
      prior.weights = rep(1,16)
   )

   modelFit <- do.call("full.laplace_MA", argListFit)
   # modelFit <- do.call("sampling_MA", argListFit)
   # modelFit$w.msg

   expect_s3_class(modelFit, 'BMADR')

   # Plots

   modelPlot <- plot.BMADR(
      mod.obj = modelFit,
      weight_type = "LP",
      type = ifelse(argListFit$data.N$data$is_increasing == 1, "increasing", "decreasing"),
      title = '', all = F, log = T
   )

   expect_type(modelPlot, 'list')
   expect_equal(names(modelPlot), c("BMDs", "weights", "model_fit_N", "model_fit_LN", "model_fit", "MA_fit"))

})
