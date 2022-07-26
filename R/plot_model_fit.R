
#' Basic plotting for the methods using covariates
#'
#' @param x BMDBMA model object
#' @param increasing logical indicating whether the continuous data are increasing (T) or decreasing (F)
#'              
#' @return a basic plot
#' 
#' @export basic.plot
#' 
basic.plot <- function(x, increasing){
  
  par(mfrow = c(4,4))
  
  covar <- rownames(x$MA)
  
  if(increasing == T){
    ### INCREASING = T
    
    ## E4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: E4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='E4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.E4_NI(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.E4_NI(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.E4_NI(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.E4_NI(par = x$parE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='E4_N'][1])), paste0('Weight ', round(x$weights[1], 4))), cex = 1.3, bty = 'n')
    }
    
    ## IE4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: IE4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='IE4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.IE4_NI(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.IE4_NI(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.IE4_NI(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.IE4_NI(par = x$parIE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='IE4_N'][1])), paste0('Weight ', round(x$weights[2], 4))), cex = 1.3, bty = 'n')
    }
    
    ## H4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: H4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='H4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.H4_NI(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.H4_NI(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.H4_NI(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.H4_NI(par = x$parH4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='H4_N'][1])), paste0('Weight ', round(x$weights[3], 4))), cex = 1.3, bty = 'n')
    }
    
    ## LN4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: LN4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='LN4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.LN4_NI(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.LN4_NI(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.LN4_NI(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.LN4_NI(par = x$parLN4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='LN4_N'][1])), paste0('Weight ', round(x$weights[4], 4))), cex = 1.3, bty = 'n')
    }
    
    ## G4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: G4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='G4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.G4_NI(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.G4_NI(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.G4_NI(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.G4_NI(par = x$parG4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='G4_N'][1])), paste0('Weight ', round(x$weights[5], 4))), cex = 1.3, bty = 'n')
    }
    
    ## QE4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: QE4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='QE4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.QE4_NI(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.QE4_NI(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.QE4_NI(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.QE4_NI(par = x$parQE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='QE4_N'][1])), paste0('Weight ', round(x$weights[6], 4))), cex = 1.3, bty = 'n')
    }
    
    ## P4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: P4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='P4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.P4_NI(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.P4_NI(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.P4_NI(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.P4_NI(par = x$parP4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='P4_N'][1])), paste0('Weight ', round(x$weights[7], 4))), cex = 1.3, bty = 'n')
    }
    
    ## L4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: L4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='L4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.L4_NI(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.L4_NI(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.L4_NI(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.L4_NI(par = x$parL4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='L4_N'][1])), paste0('Weight ', round(x$weights[8], 4))), cex = 1.3, bty = 'n')
    }
    
    ## E4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: E4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='E4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.E4_LNI(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.E4_LNI(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.E4_LNI(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.E4_LNI(par = x$parE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='E4_LN'][1])), paste0('Weight ', round(x$weights[9], 4))), cex = 1.3, bty = 'n')
    }
    
    ## IE4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: IE4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='IE4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.IE4_LNI(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.IE4_LNI(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.IE4_LNI(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.IE4_LNI(par = x$parIE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='IE4_LN'][1])), paste0('Weight ', round(x$weights[10], 4))), cex = 1.3, bty = 'n')
    }
    
    ## H4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: H4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='H4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.H4_LNI(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.H4_LNI(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.H4_LNI(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.H4_LNI(par = x$parH4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='H4_LN'][1])), paste0('Weight ', round(x$weights[11], 4))), cex = 1.3, bty = 'n')
    }
    
    ## LN4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: LN4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='LN4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.LN4_LNI(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.LN4_LNI(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.LN4_LNI(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.LN4_LNI(par = x$parLN4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='LN4_LN'][1])), paste0('Weight ', round(x$weights[12], 4))), cex = 1.3, bty = 'n')
    }
    
    ## G4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: G4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='G4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.G4_LNI(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.G4_LNI(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.G4_LNI(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.G4_LNI(par = x$parG4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='G4_LN'][1])), paste0('Weight ', round(x$weights[13], 4))), cex = 1.3, bty = 'n')
    }
    
    ## QE4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: QE4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='QE4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.QE4_LNI(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.QE4_LNI(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.QE4_LNI(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.QE4_LNI(par = x$parQE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='QE4_LN'][1])), paste0('Weight ', round(x$weights[14], 4))), cex = 1.3, bty = 'n')
    }
    
    ## P4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: P4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='P4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.P4_LNI(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.P4_LNI(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.P4_LNI(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.P4_LNI(par = x$parP4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='P4_LN'][1])), paste0('Weight ', round(x$weights[15], 4))), cex = 1.3, bty = 'n')
    }
    
    ## L4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: L4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='L4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.L4_LNI(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.L4_LNI(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.L4_LNI(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.L4_LNI(par = x$parL4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='L4_LN'][1])), paste0('Weight ', round(x$weights[16], 4))), cex = 1.3, bty = 'n')
    }
    
    ### DECREASING = T
  }else{
    ## E4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: E4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='E4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.E4_ND(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.E4_ND(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='E4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.E4_ND(par = x$parE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.E4_ND(par = x$parE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='E4_N'][1])), paste0('Weight ', round(x$weights[1], 4))), cex = 1.3, bty = 'n')
    }
    
    ## IE4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: IE4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='IE4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.IE4_ND(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.IE4_ND(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='IE4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.IE4_ND(par = x$parIE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.IE4_ND(par = x$parIE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='IE4_N'][1])), paste0('Weight ', round(x$weights[2], 4))), cex = 1.3, bty = 'n')
    }
    
    ## H4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: H4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='H4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.H4_ND(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.H4_ND(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='H4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.H4_ND(par = x$parH4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.H4_ND(par = x$parH4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='H4_N'][1])), paste0('Weight ', round(x$weights[3], 4))), cex = 1.3, bty = 'n')
    }
    
    ## LN4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: LN4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='LN4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.LN4_ND(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.LN4_ND(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='LN4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.LN4_ND(par = x$parLN4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.LN4_ND(par = x$parLN4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='LN4_N'][1])), paste0('Weight ', round(x$weights[4], 4))), cex = 1.3, bty = 'n')
    }
    
    ## G4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: G4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='G4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.G4_ND(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.G4_ND(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='G4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.G4_ND(par = x$parG4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.G4_ND(par = x$parG4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='G4_N'][1])), paste0('Weight ', round(x$weights[5], 4))), cex = 1.3, bty = 'n')
    }
    
    ## QE4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: QE4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='QE4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.QE4_ND(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.QE4_ND(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='QE4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.QE4_ND(par = x$parQE4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.QE4_ND(par = x$parQE4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='QE4_N'][1])), paste0('Weight ', round(x$weights[6], 4))), cex = 1.3, bty = 'n')
    }
    
    ## P4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: P4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='P4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.P4_ND(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.P4_ND(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='P4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.P4_ND(par = x$parP4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.P4_ND(par = x$parP4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='P4_N'][1])), paste0('Weight ', round(x$weights[7], 4))), cex = 1.3, bty = 'n')
    }
    
    ## L4_N
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: L4_N', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='L4_N'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), DRM.L4_ND(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), DRM.L4_ND(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='L4_N'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), DRM.L4_ND(par = x$parL4_N[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), DRM.L4_ND(par = x$parL4_N[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='L4_N'][1])), paste0('Weight ', round(x$weights[8], 4))), cex = 1.3, bty = 'n')
    }
    
    ## E4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: E4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='E4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.E4_LND(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.E4_LND(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='E4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.E4_LND(par = x$parE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.E4_LND(par = x$parE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='E4_LN'][1])), paste0('Weight ', round(x$weights[9], 4))), cex = 1.3, bty = 'n')
    }
    
    ## IE4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: IE4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='IE4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.IE4_LND(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.IE4_LND(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='IE4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.IE4_LND(par = x$parIE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.IE4_LND(par = x$parIE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='IE4_LN'][1])), paste0('Weight ', round(x$weights[10], 4))), cex = 1.3, bty = 'n')
    }
    
    ## H4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: H4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='H4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.H4_LND(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.H4_LND(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='H4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.H4_LND(par = x$parH4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.H4_LND(par = x$parH4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='H4_LN'][1])), paste0('Weight ', round(x$weights[11], 4))), cex = 1.3, bty = 'n')
    }
    
    ## LN4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: LN4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='LN4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.LN4_LND(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.LN4_LND(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='LN4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.LN4_LND(par = x$parLN4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.LN4_LND(par = x$parLN4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='LN4_LN'][1])), paste0('Weight ', round(x$weights[12], 4))), cex = 1.3, bty = 'n')
    }
    
    ## G4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: G4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='G4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.G4_LND(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.G4_LND(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='G4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.G4_LND(par = x$parG4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.G4_LND(par = x$parG4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='G4_LN'][1])), paste0('Weight ', round(x$weights[13], 4))), cex = 1.3, bty = 'n')
    }
    
    ## QE4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: QE4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='QE4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.QE4_LND(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.QE4_LND(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='QE4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.QE4_LND(par = x$parQE4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.QE4_LND(par = x$parQE4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='QE4_LN'][1])), paste0('Weight ', round(x$weights[14], 4))), cex = 1.3, bty = 'n')
    }
    
    ## P4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: P4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='P4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.P4_LND(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.P4_LND(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='P4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.P4_LND(par = x$parP4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.P4_LND(par = x$parP4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='P4_LN'][1])), paste0('Weight ', round(x$weights[15], 4))), cex = 1.3, bty = 'n')
    }
    
    ## L4_LN
    plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
         x$data$y[x$data$cov == covar[1]], main = 'Fitted model: L4_LN', 
         ylim = c(min(x$data$y)-0.5, max(x$data$y)+0.5), xlab = '', ylab = 'Response')
    j = 2
    for(i in 2:length(covar)){
      points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
             x$data$y[x$data$cov == covar[i]], pch = j, col = j)
      j = j + 1
    }
    
    if(!is.na(x$summary$Submodel[x$summary$Model=='L4_LN'][1])){
      
      if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'a_sigma2'){
        for(i in 1:length(covar)){
          pars <- c(paste0('par1[',i,']'),
                    'par2[1]',
                    'par3',
                    'par4[1]')
          lines(seq(0,1,0.01), exp(DRM.L4_LND(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        } 
      }else if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'BMD_d'){
        for(i in 1:length(covar)){
          pars <- c('par1[1]',
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']'))
          lines(seq(0,1,0.01), exp(DRM.L4_LND(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
        }
      }else if(x$summary$Submodel[x$summary$Model=='L4_LN'][1] == 'all'){
        
        for(i in 1:length(covar)){
          
          pars <- c(paste0('par1[',i,']'),
                    paste0('par2[',i,']'),
                    'par3',
                    paste0('par4[',i,']')
          )
          
          lines(seq(0,1,0.01), exp(DRM.L4_LND(par = x$parL4_LN[pars], seq(0,1,0.01), 0.1, shift = 0)), col = i, lwd = 2)
          
        }
        
      }else{
        lines(seq(0,1,0.01), exp(DRM.L4_LND(par = x$parL4_LN[c(1,2,9,4)], seq(0,1,0.01), 0.1, shift = 0)), col = 1, lwd = 2)
      }
      
      legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='L4_LN'][1])), paste0('Weight ', round(x$weights[16], 4))), cex = 1.3, bty = 'n')
    }
  }
}

#' @rdname basic.plot
#' @export
basic.plotQ <- function(x){
  
  par(mfrow = c(2,4))
  
  covar <- rownames(x$MA)
  
  ## E4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: E4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='E4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='E4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.E4_Q(par = x$parE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='E4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.E4_Q(par = x$parE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='E4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.E4_Q(par = x$parE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.E4_Q(par = x$parE4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='E4_Q'][1])), paste0('Weight ', round(x$weights[1], 4))), cex = 1.3, bty = 'n')
  }
  
  ## IE4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: IE4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='IE4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='IE4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.IE4_Q(par = x$parIE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='IE4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.IE4_Q(par = x$parIE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='IE4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.IE4_Q(par = x$parIE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.IE4_Q(par = x$parIE4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='IE4_Q'][1])), paste0('Weight ', round(x$weights[2], 4))), cex = 1.3, bty = 'n')
  }
  
  ## H4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: H4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='H4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='H4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.H4_Q(par = x$parH4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='H4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.H4_Q(par = x$parH4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='H4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.H4_Q(par = x$parH4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.H4_Q(par = x$parH4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='H4_Q'][1])), paste0('Weight ', round(x$weights[3], 4))), cex = 1.3, bty = 'n')
  }
  
  ## LN4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: LN4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='LN4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='LN4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.LN4_Q(par = x$parLN4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='LN4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.LN4_Q(par = x$parLN4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='LN4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.LN4_Q(par = x$parLN4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.LN4_Q(par = x$parLN4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='LN4_Q'][1])), paste0('Weight ', round(x$weights[4], 4))), cex = 1.3, bty = 'n')
  }
  
  ## G4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: G4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='G4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='G4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.G4_Q(par = x$parG4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='G4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.G4_Q(par = x$parG4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='G4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.G4_Q(par = x$parG4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.G4_Q(par = x$parG4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='G4_Q'][1])), paste0('Weight ', round(x$weights[5], 4))), cex = 1.3, bty = 'n')
  }
  
  ## QE4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: QE4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='QE4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='QE4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.QE4_Q(par = x$parQE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='QE4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.QE4_Q(par = x$parQE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='QE4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.QE4_Q(par = x$parQE4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.QE4_Q(par = x$parQE4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='QE4_Q'][1])), paste0('Weight ', round(x$weights[6], 4))), cex = 1.3, bty = 'n')
  }
  
  ## P4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: P4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='P4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='P4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.P4_Q(par = x$parP4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='P4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.P4_Q(par = x$parP4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='P4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.P4_Q(par = x$parP4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.P4_Q(par = x$parP4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='P4_Q'][1])), paste0('Weight ', round(x$weights[7], 4))), cex = 1.3, bty = 'n')
  }
  
  ## L4_Q
  plot(x$data$x[x$data$cov == covar[1]]/max(x$data$x[x$data$cov == covar[1]]), 
       x$data$y[x$data$cov == covar[1]]/x$data$n[x$data$cov == covar[1]], main = 'Fitted model: L4_Q', 
       ylim = c(0,1), xlab = '', ylab = 'Response')
  j = 2
  for(i in 2:length(covar)){
    points(x$data$x[x$data$cov == covar[i]]/max(x$data$x[x$data$cov == covar[i]]), 
           x$data$y[x$data$cov == covar[i]]/x$data$n[x$data$cov == covar[1]], pch = j, col = j)
    j = j + 1
  }
  
  if(!is.na(x$summary$Submodel[x$summary$Model=='L4_Q'][1])){
    
    if(x$summary$Submodel[x$summary$Model=='L4_Q'][1] == 'a_sigma2'){
      for(i in 1:length(covar)){
        pars <- c(paste0('par1[',i,']'),
                  'par2[1]',
                  
                  'par3[1]')
        lines(seq(0,1,0.01), DRM.L4_Q(par = x$parL4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      } 
    }else if(x$summary$Submodel[x$summary$Model=='L4_Q'][1] == 'BMD_d'){
      for(i in 1:length(covar)){
        pars <- c('par1[1]',
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']'))
        lines(seq(0,1,0.01), DRM.L4_Q(par = x$parL4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
      }
    }else if(x$summary$Submodel[x$summary$Model=='L4_Q'][1] == 'all'){
      
      for(i in 1:length(covar)){
        
        pars <- c(paste0('par1[',i,']'),
                  paste0('par2[',i,']'),
                  
                  paste0('par3[',i,']')
        )
        
        lines(seq(0,1,0.01), DRM.L4_Q(par = x$parL4_Q[pars], seq(0,1,0.01), 0.1), col = i, lwd = 2)
        
      }
      
    }else{
      lines(seq(0,1,0.01), DRM.L4_Q(par = x$parL4_Q[c(1,2,9,4)], seq(0,1,0.01), 0.1), col = 1, lwd = 2)
    }
    
    legend('topleft', c(paste0('Best submodel: ', paste0(x$summary$Submodel[x$summary$Model=='L4_Q'][1])), paste0('Weight ', round(x$weights[8], 4))), cex = 1.3, bty = 'n')
  }
  
}