#' Function to summarize individual data
#'
#' @param data data in order of: dose, response, covariate OR dose, response, litter
#' @param type either 'continuous' or 'quantal'
#' @param cluster logical indicating whether the data include a litter effect (defaults to FALSE)
#' @param covar logical indicating whether the data include a covariate (defaults to FALSE)
#'
#' @return dataframe of summary data
#'
#' @export summarize.indiv.data
summarize.indiv.data <- function(data,
                                 type = c('continuous','quantal'), cluster = FALSE, covar = FALSE){
  require('dplyr')

  if(type == 'continuous'){

    if(cluster == FALSE && covar == FALSE){

      data = data[order(data[, 1]), ]
      doses = data[, 1]
      maxDose = max(doses)
      dose.a=sort(unique(doses))
      N=length(dose.a)
      mean.a=rep(NA,N)
      sd.a=rep(NA,N)
      n.a=rep(NA,N)
      y = data[, 2]
      for (iu in (1:N)){
        mean.a[iu]=mean(y[doses==dose.a[iu]])
        sd.a[iu]=sd(y[doses==dose.a[iu]])
        n.a[iu]=sum(doses==dose.a[iu])
      }
      dose.a = dose.a/maxDose

      summary.data <- data.frame(Dose = dose.a*maxDose,
                                 Response = mean.a,
                                 SD = sd.a,
                                 N = n.a)


    }else if(cluster == FALSE && covar == TRUE){

      data = data[order(data[, 1]), ]
      data$dose = data[,1]
      data$resp = data[,2]
      data$cov = data[,3]
      indiv.data <- data %>%
        dplyr::group_by(dose, cov) %>%
        dplyr::arrange(by_group = dose) %>%
        dplyr::summarise(mean = mean(resp, na.rm = T), sd = sd(resp, na.rm=T), n = n())
      dose.a = indiv.data$dose
      maxDose = max(dose.a)
      mean.a = indiv.data$mean
      sd.a = indiv.data$sd
      n.a = indiv.data$n
      N = length(dose.a)
      dose.a = dose.a/maxDose
      covar = indiv.data$cov

      summary.data <- data.frame(Dose = dose.a*maxDose,
                                 Response = mean.a,
                                 SD = sd.a,
                                 N = n.a,
                                 Covariate = covar)

    }else if(cluster == TRUE && covar == FALSE){

      library(dplyr)
      indiv.data <- data.frame(dose = data[,1],
                               response = data[,2],
                               litter = data[,3])
      indiv.data <- indiv.data %>%
        # dplyr::group_by(dose, litter) %>%
        dplyr::group_by(dose) %>%
        dplyr::arrange(by_group = dose) %>% # order dose groups
        dplyr::summarise(mean.a = mean(response), sd.a = sd(response), n.a = n())
      # indiv.data <- indiv.data %>%
      #   dplyr::summarise()
      #
      #   dplyr::mutate(cluster = dplyr::cur_group_id(),
      #                 count = n())


      summary.data <- data.frame(Dose = indiv.data$dose,
                                    Response = indiv.data$mean.a,
                                    SD = indiv.data$sd.a,
                                    N = indiv.data$n.a)


    }else if(cluster == TRUE && covar == TRUE){
      stop('Covariates not implemented for clustered data')
    }

  }else if(type == 'quantal'){

    if(cluster == FALSE && covar == FALSE){

      doses = data[, 1]
      dose.a = unique(doses)
      N = length(dose.a)
      maxDose = max(dose.a)
      y.a = rep(NA, length(unique(doses)))
      n.a = rep(NA, length(unique(doses)))
      ybin = data[, 2]
      nbin = data[, 3]
      for(iu in 1:N){
        y.a[iu] = sum(ybin[doses==dose.a[iu]])
        n.a[iu] = sum(nbin[doses==dose.a[iu]])
      }
      dose.a = dose.a/maxDose

      summary.data <- data.frame(Dose = dose.a*maxDose,
                                 n.Events = y.a,
                                 N = n.a)

    }else if(cluster == FALSE && covar == TRUE){

      data = data[order(data[, 1]), ]
      data$dose = data[,1]
      data$ybin = data[,2]
      data$cov = data[,3]
      indiv.data <- data %>%
        dplyr::group_by(dose, cov) %>%
        dplyr::arrange(by_group = dose) %>%
        dplyr::summarise(y.a = sum(ybin, na.rm = T), n.a = n())
      dose.a = indiv.data$dose
      maxDose = max(dose.a)
      y.a = indiv.data$y.a
      n.a = indiv.data$n
      N = length(dose.a)
      dose.a = dose.a/maxDose
      covar = indiv.data$cov
      dose.a = dose.a/maxDose

      summary.data <- data.frame(Dose = dose.a*maxDose,
                                 n.Events = y.a,
                                 N = n.a,
                                 Covariate = covar)

    }else if(cluster == TRUE && covar == FALSE){

      doses = data[, 1]
      dose.a = unique(doses)
      maxDose = max(dose.a)
      N = length(dose.a)
      # litter = data[,3]
      # y.a = rep(NA, length(unique(litter)))
      # n.a = rep(NA, length(unique(litter)))
      # dose.a = rep(NA, length(unique(litter)))
      y.a = rep(NA, length(unique(dose.a)))
      n.a = rep(NA, length(unique(dose.a)))
      ybin = data[, 2]
      nbin = data[, 3]
      # id = 1
      # for(iu in unique(litter)){
      #   y.a[id] = sum(ybin[litter == iu])
      #   n.a[id] = sum(litter == iu)
      #   dose.a[id] = unique(doses[litter == iu])
      #   id = id + 1
      # }
      for(iu in 1:N){
        y.a[iu] = sum(ybin[doses==dose.a[iu]])
        n.a[iu] = sum(nbin[doses==dose.a[iu]])
      }
      dose.a = dose.a/maxDose

      summary.data <- data.frame(Dose = dose.a*maxDose,
                                 n.Events = y.a,
                                 N = n.a)


    }else if(cluster == TRUE && covar == TRUE){
      stop('Covariates not implemented for clustered data')
    }

  }

  return(summary.data)

}
