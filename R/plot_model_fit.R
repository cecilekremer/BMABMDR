
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

  covar <- unique(x$data[,5])

  if(increasing == TRUE){

    model_pars <- get('x')[[paste0('par',model_name)]]
    DRM <- get(paste0('DRM.', model_name, 'I'))

    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]),
         x$data$y[x$data$cov == covar[1]], main = paste0('Fitted model: ', model_name),
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = 'Dose', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]),
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    if(!is.na(x$summary$Submodel[x$summary$Model==model_name][1])){

      if(x$summary$Submodel[x$summary$Model==model_name][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift)), col = i, lwd = 2)
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)
          }
        }
      }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift)), col = i, lwd = 2)
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)
          }
        }
      }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'all'){

        for(i in 1:length(covar)){

          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift)), col = i, lwd = 2)
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)
          }

        }

      }else{
        if(grepl('_LN',model_name)){
          lines(seq(0,1,0.01), exp(DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), 0.1, x$shift)), col = 1, lwd = 2)
        }else{
          lines(seq(0,1,0.01), DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), 0.1, x$shift), col = 1, lwd = 2)
        }
      }

      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model==model_name][1])),
                          paste0('Weight ', round(x$weights[model_name], 4))), cex = 1.3, bty = 'n')
    }

  }else{

    model_pars <- get('x')[[paste0('par',model_name)]]
    DRM <- get(paste0('DRM.', model_name, 'D'))

    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]),
         x$data$y[x$data$cov == covar[1]], main =  paste0('Fitted model: ', model_name),
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = 'Dose', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]),
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }

    if(!is.na(x$summary$Submodel[x$summary$Model==model_name][1])){

      if(x$summary$Submodel[x$summary$Model==model_name][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2))
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)

          }
        }
      }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift)), col = i, lwd = 2)
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)
          }
        }
      }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'all'){

        for(i in 1:length(covar)){

          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          if(grepl('_LN',model_name)){
            lines(seq(0,1,0.01), exp(DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift)), col = i, lwd = 2)
          }else{
            lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1, x$shift), col = i, lwd = 2)
          }
        }

      }else{
        if(grepl('_LN',model_name)){
          lines(seq(0,1,0.01), exp(DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), 0.1, x$shift)), col = 1, lwd = 2)
        }else{
          lines(seq(0,1,0.01), DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), 0.1, x$shift), col = 1, lwd = 2)
        }
      }

      legend('topright', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model==model_name][1])),
                           paste0('Weight ', round(x$weights[model_name], 4))), cex = 1.3, bty = 'n')
    }
  }


}

#' @rdname basic.plot
#' @export
basic.plotQ <- function(x, model_name){

  covar <- unique(x$data[,5])
  model_pars <- get('x')[[paste0('par',model_name)]]
  DRM <- get(paste0('DRM.', model_name))

  ## E4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]),
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = paste0('Fitted model: ', model_name),
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]),
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }

  if(!is.na(x$summary$Submodel[x$summary$Model==model_name][1])){

    if(x$summary$Submodel[x$summary$Model==model_name][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',

                  'par3[1]')
        lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),

                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model==model_name][1] == 'all'){

      for(i in 1:length(covar)){

        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),

                  paste0('par3[',i,']')
        )

        lines(seq(0,1,0.01), DRM(par = model_pars[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)

      }

    }else{
      lines(seq(0,1,0.01), DRM(par = model_pars[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }

    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model==model_name][1])),
                        paste0('Weight ', round(x$weights[model_name], 4))), cex = 1.3, bty = 'n')
  }
}
