
#' Basic plotting for the methods using covariates
#'
#' @param x BMDBMA model object
#' @param model_name model for which the fit should be plotted (one of get_models('continuous'))
#' @param increasing logical indicating whether the continuous data are increasing (T) or decreasing (F)
#'
#' @return a basic plot
#'
#' @export basic.plot
#'
basic.plot <- function(x, model_name, increasing){

  covar <- as.factor(unique(x$data[,5]))
  temp <- temp1 <- temp2 <- temp3 <- NULL
  x$data$geom.y <- NtoLN(x$data$y, x$data$s)[1:dim(x$data)[1]]
  x$data$geom.s <- NtoLN(x$data$y, x$data$s)[(dim(x$data)[1]+1):(2*dim(x$data)[1])]
  maxDose <- x$maxDose

  dose = unique(x$data$x)

  x$data$x2 = ifelse(x$data$x == 0, log10(min(x$data$x[x$data$x != 0])/4), log10(x$data$x))

  if(min(dose) == 0){
    plot.labs = dose
    plot.breaks = c(log10(min(dose[dose!=0])/4), log10(dose[2:length(dose)]))
  }else{
    plot.labs = c(0, dose)
    plot.breaks = c(log10(min(dose[dose!=0])/4), log10(dose))
  }

  # if(increasing == TRUE){

  model_pars <- get('x')[[paste0('par',model_name)]]
  if(increasing == TRUE){
    DRM <- get(paste0('DRM.', model_name, 'I'))
  }  else{
    DRM <- get(paste0('DRM.', model_name, 'D'))
  }

  if(!is.na(x$summary$Submodel[x$summary$Model==model_name][1])){

    if(x$summary$Submodel[x$summary$Model==model_name][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  'par3',
                  'par4[1]')

        if(grepl('_LN',model_name)){
          temp <- c(temp, exp(DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift)))
        }else{
          temp <- c(temp, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift))
        }
      }
      if(grepl('_LN',model_name)){
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = geom.y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = geom.y - geom.s, ymax = geom.y + geom.s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }else{
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = y - s, ymax = y + s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }
    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  'par3',
                  paste0('par4[',i,']'))
        if(grepl('_LN',model_name)){
          temp1 <- c(temp1, exp(DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift)))
        }else{
          temp1 <- c(temp1, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift))
        }
      }
      if(grepl('_LN',model_name)){
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp1,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = geom.y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov))          +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = geom.y - geom.s, ymax = geom.y + geom.s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }else{
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp1,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = y - s, ymax = y + s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }

    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'all'){

      for(i in 1:length(covar)){

        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  'par3',
                  paste0('par4[',i,']')
        )
        if(grepl('_LN',model_name)){
          temp2 <- c(temp2, exp(DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift)))
        }else{
          temp2 <- c(temp2, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift))
        }

      }
      if(grepl('_LN',model_name)){
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp2,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = geom.y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = geom.y - geom.s, ymax = geom.y + geom.s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }else{
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                               y = temp2,
                               cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y, colour = cov))+
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = y - s, ymax = y + s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }

    }else{
      if(grepl('_LN',model_name)){
        temp3 <- c(temp3, exp(DRM(par = model_pars[c(1,2,9,4)], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift)))
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)), y = temp3)
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = geom.y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = geom.y - geom.s, ymax = geom.y + geom.s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }else{
        temp3 <- c(temp3, DRM(par = model_pars[c(1,2,9,4)], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q, x$shift))
        dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)), y = temp3)
        p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y, colour = cov)) +
          geom_line(data = dataTemp, aes(x = x, y = y)) +
          geom_errorbar(data = x$data, mapping = aes(x = x2, ymin = y - s, ymax = y + s, colour = cov),
                        size = 1, width = NA, linetype = 'dotted') +
          scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
      }
    }

    subtit <- paste0(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model==model_name][1])," ("),
                     paste0('Weight = ', round(x$weights[model_name], 4)),")")
    pp <- p + labs(title = paste0('Fitted model: ', model_name),
                   subtitle = subtit,
                   x = 'Dose', y = 'Response', colour = 'Covariate')
  }

  return(pp)
}

#' @rdname basic.plot
#' @export
basic.plotQ <- function(x, model_name){

  covar <- unique(x$data[,4])
  model_pars <- get('x')[[paste0('par',model_name)]]
  DRM <- get(paste0('DRM.', model_name))
  temp <- temp1 <- temp2 <- temp3 <- NULL
  maxDose <- x$maxDose

  dose = unique(x$data$x)

  x$data$x2 = ifelse(x$data$x == 0, log10(min(x$data$x[x$data$x != 0])/4), log10(x$data$x))

  if(min(dose) == 0){
    plot.labs = dose
    plot.breaks = c(log10(min(dose[dose!=0])/4), log10(dose[2:length(dose)]))
  }else{
    plot.labs = c(0, dose)
    plot.breaks = c(log10(min(dose[dose!=0])/4), log10(dose))
  }

  if(!is.na(x$summary$Submodel[x$summary$Model==model_name][1])){


    if(x$summary$Submodel[x$summary$Model==model_name][1] == 'background'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',

                  'par3[1]')
        temp <- c(temp, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q))
      }
      dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)),
                             y = temp,
                             cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
      p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y/n, colour = cov)) +
        geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
        scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)

    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),

                  paste0('par3[',i,']'))
        # lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), x$q), col = i, lwd = 2)
        temp1 <- c(temp1, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q))
      }
      dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)), y = temp1,
                             cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
      p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y/n, colour = cov)) +
        geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
        scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)

    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'all'){

      for(i in 1:length(covar)){

        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),

                  paste0('par3[',i,']')
        )

        # lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), x$q), col = i, lwd = 2)
        temp2 <- c(temp2, DRM(par = model_pars[pars], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q))
      }
      dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)), y = temp2,
                             cov = rep(covar, each = length(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)))))
      p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y/n, colour = cov)) +
        geom_line(data = dataTemp, aes(x = x, y = y, colour = cov)) +
        scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)

    }else{
      # lines(seq(0,1,0.01), DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), x$q), col = 1, lwd = 2)
      temp3 <- c(temp3, DRM(par = model_pars[c(1,2,9,4)], (10^(log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01))))/maxDose, x$q))
      dataTemp <- data.frame(x = log10(seq(ifelse(min(x$data$x) == 0, min(x$data$x[x$data$x != 0])/4, min(x$data$x)), maxDose, 0.01)), y = temp3)
      p <- ggplot() + geom_point(data = x$data, aes(x = x2, y = y/n, colour = cov)) +
        geom_line(data = dataTemp, aes(x = x, y = y)) +
        scale_x_continuous(trans = 'log10', labels = plot.labs, breaks = plot.breaks)
    }

    subtit <- paste0(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model==model_name][1])," ("),
                     paste0('Weight = ', round(x$weights[model_name], 4)),")")
    pp <- p + labs(title = paste0('Fitted model: ', model_name),
                   subtitle = subtit,
                   x = 'Dose', y = 'Response', colour = 'Covariate')

  }


  return(pp)
}
