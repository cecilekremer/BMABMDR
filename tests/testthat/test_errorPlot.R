# Project: bmdBayesian
#
# Author: wverlinden
###############################################################################

#' Function to generate the analysis data for summary data taking into account the outcome of the Bartlett test
#' Order of columns needs to be: dose, response/mean, sd/se, n
#'
#' @param argList list with arguments for PREP_DATA-functions
#' @param selectedDistribution character vector holding the selected distributions (either "Normal", "Lognormal" or both)
#' @param pval numeric indicating alpha level
#'
#' @importFrom BMABMDR bartlett NtoLN bartlett NLN_test
#'
#' @return list with prepared data for each analysis
#' @export



test_that("Prepare data",{


  getAnalysisData <- function(argList, selectedDistribution, pval = 0.05){

    summ.data <- argList$data

    # Convert to SD if necessary to avoid errors when using functions requiring SD (such as bartlett, NtoLN)
    if (argList$sd == FALSE){
      summ.data[,3] <- summ.data[,3]*sqrt(summ.data[,4])
      argList$sd <- TRUE
      argList$data <- summ.data
    }

    ## Bartlett test
    # normal distribution
    b.Normal <- bartlett(sd = summ.data[,3], n = summ.data[,4])

    # lognormal distribution
    ## Convert arithmetic to geometric summary stats
    summ.data.LN <- data.frame(x = summ.data[,1],
                               y = NtoLN(summ.data[,2], summ.data[,3])[1:length(summ.data[,1])],
                               s = NtoLN(summ.data[,2], summ.data[,3])[(length(summ.data[,1])+1):(2*length(summ.data[,1]))],
                               n = summ.data[,4]
    )

    b.Lognormal <- bartlett(sd = log(summ.data.LN[,3]), n = summ.data.LN[,4])

    ## Different scenario's depending on outcome of Bartlett test
    analysisData <- if (b.Normal[2] >= pval & b.Lognormal[2] >= pval){

      # If assumption met for both distributions
      data_N <- do.call("PREP_DATA_N", argList)
      data_LN <- do.call("PREP_DATA_LN", argList)

      prior.weights <- c(rep(1,8), rep(1,8)) # Do analysis for both distributions

      list(data = list(
        default = list(data_N = data_N, data_LN = data_LN, prior.weights = prior.weights, warning = NULL)
      ), warning = NULL
      )

    } else if(b.Normal[2] >= pval & b.Lognormal[2] < pval){

      ## If assumption met for Normal only
      # Data for default analyses with Normal only
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

      # Make sublists per scenario

      list(data = list(
        default = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
        defaultAdaptedWeights = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(0,8))),  # set weights to 0 for Lognormal
        minVar = list(data_N = data_N, data_LN = data_LN_sens1, prior.weights = c(rep(1,8), rep(1,8))),
        maxVar = list(data_N = data_N, data_LN = data_LN_sens2, prior.weights = c(rep(1,8), rep(1,8)))

      ), warning = NULL
      )

    } else if(b.Normal[2] < pval & b.Lognormal[2] >= pval){

      ## If assumption met for Lognormal only
      # Data for default analyses with Lognormal only
      data_N = do.call("PREP_DATA_N", argList)
      data_LN = do.call("PREP_DATA_LN", argList)

      ## Data for sensitivity analyses
      summ.data.N.sens1 <- summ.data
      summ.data.N.sens1$s <- min(summ.data$s) # minimum arithmetic sd
      argList$data <- summ.data.N.sens1

      data_N_sens1 = do.call("PREP_DATA_N", argList)

      summ.data.N.sens2 <- summ.data
      summ.data.N.sens2$s <- max(summ.data$s) # maximum arithmetic sd
      argList$data <- summ.data.N.sens2

      data_N_sens2 = do.call("PREP_DATA_N", argList)

      list(data = list(
        default = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
        defaultAdaptedWeights = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(0,8), rep(1,8))),  # set weights to 0 for Normal
        minVar = list(data_N = data_N_sens1, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
        maxVar = list(data_N = data_N_sens2, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8)))
      ), warning = NULL
      )
    } else if(b.Normal[2] < pval & b.Lognormal[2] < pval){

      ## If assumption met for none
      # Data for regular analyses with both distributions, but warning should be given about violation of assumptions
      data_N = do.call("PREP_DATA_N", argList)
      data_LN = do.call("PREP_DATA_LN", argList)

      # Data for analysis with smallest SD
      summ.data.N.sens1 <- summ.data
      summ.data.N.sens1$s <- min(summ.data$s)
      argList$data <- summ.data.N.sens1

      data_N_sens1 = do.call("PREP_DATA_N", argList)

      summ.data.LN.sens1 <- summ.data.LN
      summ.data.LN.sens1$s <- min(summ.data.LN.sens1$s)
      argList$data <- summ.data.LN.sens1
      argList <- append(argList, list(geom.stats = T)) # specify that geometric summary stats are given as input

      data_LN_sens1 <- do.call("PREP_DATA_LN", argList)

      # Data for analysis with largest SD
      summ.data.N.sens2 <- summ.data
      summ.data.N.sens2$s <- max(summ.data.N.sens2$s)
      argList$geom.stats = FALSE # revert back to false for normal data
      argList$data <- summ.data.N.sens2

      data_N_sens2 = do.call("PREP_DATA_N", argList)

      summ.data.LN.sens2 <- summ.data.LN
      summ.data.LN.sens2$s <- max(summ.data.LN.sens2$s)
      argList$data <- summ.data.LN.sens2
      argList$geom.stats = TRUE # specify that geometric summary stats are given as input

      data_LN_sens2 <- do.call("PREP_DATA_LN", argList)

      list(data = list(
        default = list(data_N = data_N, data_LN = data_LN, prior.weights = c(rep(1,8), rep(1,8))),
        minVar = list(data_N = data_N_sens1, data_LN = data_LN_sens1, prior.weights = c(rep(1,8), rep(1,8))),
        maxVar = list(data_N = data_N_sens2, data_LN = data_LN_sens2, prior.weights = c(rep(1,8), rep(1,8)))
      ), warning = TRUE
      )

    }

    # Remove defaultAdaptedWeights scenario if any Distribution type has been unchecked
    if (length(selectedDistribution) < 2){
      analysisData$data$defaultAdaptedWeights <- NULL
    }

    return(analysisData)
  }

  dataDir <- '~/GitHub/BMABMDR/data/'

  data <- read.table(file.path(dataDir, "test_data_cont.txt"), header = T, sep = "\t", dec = ".")

  # data <- read.table(file.path(dataDir, "test new EFSA software.txt"), header = T, sep = "\t", dec = ".")

  summ.data <- data.frame(
    x = data[, "doseN"],
    y = data[, "Response"],
    s = data[, "Stdev"],
    n = data[, "N"]
  )

  ## Preparing data
  # bmr
  argList = list(
    data = summ.data,
    q = 0.05,
    sumstats = TRUE,
    sd = TRUE,
    extended = TRUE,
    prior.d = 'N11'
  )

  # Perform sensitivity analyses
  analysisData <- getAnalysisData(
    argList = argList,
    selectedDistribution = "Lognormal"
  )

  ## Fit models
  # sampling specification
  ndraws = 30000
  seed = 123

  # prior model weights
  prior.weights = c(rep(0,8), rep(1,8))
  pvec = c(0.05,0.5,0.95)

  # Fit models
  FLBMD = full.laplace_MA(
    data.N = analysisData$data$minVar$data_N,
    data.LN = analysisData$data$minVar$data_LN,
    prior.weights = prior.weights,
    ndraws = ndraws,
    seed = seed,
    pvec = pvec
  )

  ## Create plots
  pFLBMD = plot.BMADR(
    mod.obj = FLBMD,
    weight_type = "LP",
    log = FALSE,
    include_data = TRUE,
    all = FALSE,
    title = ''
  )

  expect_type(pFLBMD, 'list')

  SBMD = sampling_MA(
    data.N = analysisData$data$minVar$data_N,
    data.LN = analysisData$data$minVar$data_LN,
    prior.weights = prior.weights,
    ndraws = ndraws,
    seed = seed,
    pvec = pvec
  )

  pSBMD = plot.BMADR(SBMD, 'increasing', weight_type = 'BS', all = F, title = '')
  expect_type(pSBMD, 'list')

})
