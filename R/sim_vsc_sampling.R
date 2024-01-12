
args <- commandArgs(trailingOnly = TRUE)

# out = args[1] #working directory
# cat(",out=",out)
cores = as.numeric(args[1]) #number of cores to run in parallel
cat(",cores=",cores)
# simID = args[3]
# cat(",simID=", simID)

library(BMABMDR)
library(mc2d)
library(bbmle)
library(foreach)
library(doParallel)

## Functions
source('/scratch/leuven/326/vsc32693/testEFSA/FUNS/sim_functions.R')
source('/scratch/leuven/326/vsc32693/testEFSA/FUNS/fun_fit_pert.R')
source('/scratch/leuven/326/vsc32693/testEFSA/FUNS/fun_cum_meta.R')

load('/scratch/leuven/326/vsc32693/testEFSA/DATA/simdata_hill_logN_g5_sim5.RData')

# Sampling specifications
q = 0.22
prior.weights = rep(1, 16)
ndr = 30000
nrch = 3; nriter = 3000; wu = 1000; dl = 0.8; trd = 10; sd = 123
pvec = c(0.05,0.5,0.95)

cl <- makeCluster(cores)
registerDoParallel(cl)

nsim <- 100

foreach(s = 1:nsim, .packages=c("gamlss","mvtnorm","methods","rstan","bbmle","bridgesampling","mc2d","bayesplot","posterior","bbmle","BMABMDR")) %dopoar% {

  datasets <- sim_data[sim_data[,61] == s, ]
  
  step <- 1
  
  ## Create summary datasets
  for(d in 1:3){
    dose.a=sort(unique(doses))
    N=length(dose.a)
    mean.a=rep(NA,N)
    sd.a=rep(NA,N)
    n.a=rep(NA,N)
    y = sim_data[d, 1:60]
    for (iu in (1:N)){
      mean.a[iu]=mean(y[doses==dose.a[iu]])
      sd.a[iu]=sd(y[doses==dose.a[iu]])
      n.a[iu]=sum(doses==dose.a[iu])
    }
    summ.data = data.frame(x = dose.a, y = mean.a, s = sd.a, n = n.a)
    assign(paste0('dataset_', d), summ.data)
  }
  
  ## New data
  dose.a=sort(unique(doses))
  N=length(dose.a)
  mean.a=rep(NA,N)
  sd.a=rep(NA,N)
  n.a=rep(NA,N)
  y = sim_data[4, 1:60]
  for (iu in (1:N)){
    mean.a[iu]=mean(y[doses==dose.a[iu]])
    sd.a[iu]=sd(y[doses==dose.a[iu]])
    n.a[iu]=sum(doses==dose.a[iu])
  }
  summ.data = data.frame(x = dose.a, y = mean.a, s = sd.a, n = n.a)
  assign('new_data', summ.data)
  
  ## Default analysis new dataset
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
  SBMD_default <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                  nrchains = nrch, nriterations = nriter,
                                  warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  # setwd(out)
  # write(c(s, SBMD_default$weights_bridge_sampling, SBMD_default$MA_bridge_sampling, SBMD_default$MA_post_bs, SBMD_default$convergence,
  #         quantile(SBMD_default$bkg_post_bs, seq(0,1,0.005)), quantile(SBMD_default$maxy_post_bs, seq(0,1,0.005)),
  #         as.vector(SBMD_default$E4_N['par.d',]), as.vector(SBMD_default$IE4_N['par.d',]), as.vector(SBMD_default$H4_N['par.d',]), as.vector(SBMD_default$LN4_N['par.d',]),
  #         as.vector(SBMD_default$G4_N['par.d',]), as.vector(SBMD_default$QE4_N['par.d',]), as.vector(SBMD_default$P4_N['par.d',]), as.vector(SBMD_default$L4_N['par.d',]),
  #         as.vector(SBMD_default$E4_LN['par.d',]), as.vector(SBMD_default$IE4_LN['par.d',]), as.vector(SBMD_default$H4_LN['par.d',]), as.vector(SBMD_default$LN4_LN['par.d',]),
  #         as.vector(SBMD_default$G4_LN['par.d',]), as.vector(SBMD_default$QE4_LN['par.d',]), as.vector(SBMD_default$P4_LN['par.d',]), as.vector(SBMD_default$L4_LN['par.d',])
  #         ), file = paste0('SBMD_default_sim', s, '.txt'), ncolumns = 687, append = F)
  
  out_default <- c(SBMD_default$MA_bridge_sampling, quantile(SBMD_default$bkg_post_bs, pvec), quantile(SBMD_default$maxy_post_bs, pvec))
  res.bkg <- out_default[c(4:6)]
  res.maxy <- out_default[c(7:9)]
  res.bmd <- out_default[c(1:3)]
  res.bkg <- as.data.frame(t(res.bkg))
  names(res.bkg) <- c('lower', 'median', 'upper')
  res.bkg$est <- rep('bkg', 1)
  res.maxy <- as.data.frame(t(res.maxy))
  names(res.maxy) <- c('lower', 'median', 'upper')
  res.maxy$est <- rep('maxy', 1)
  res.bmd <- as.data.frame(t(res.bmd))
  names(res.bmd) <- c('lower', 'median', 'upper')
  res.bmd$est <- rep('bmd', 1)
  
  res.all.default <- rbind(res.bkg, res.maxy, res.bmd)
  save(res.all.default, file = paste0('/scratch/leuven/326/vsc32693/testEFSA/OUT/output_default_sim', s, '.RData'))
  
  step <- 2
  
  ## Analyse historical datasets
  out.pert <- matrix(NA, nrow = 3, ncol = 24)
  for(d in 1:3){
    summ.data <- get(paste0('dataset_', d))
    data_N <- PREP_DATA_N(data = summ.data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
    data_LN <- PREP_DATA_LN(data = summ.data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
    SBMD <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                            nrchains = nrch, nriterations = nriter,
                            warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    # write(c(s, d, SBMD$weights_bridge_sampling, SBMD$MA_bridge_sampling, SBMD$MA_post_bs, SBMD$convergence,
    #         quantile(SBMD$bkg_post_bs, seq(0,1,0.005)), quantile(SBMD$maxy_post_bs, seq(0,1,0.005)),
    #         as.vector(SBMD$E4_N['par.d',]), as.vector(SBMD$IE4_N['par.d',]), as.vector(SBMD$H4_N['par.d',]), as.vector(SBMD$LN4_N['par.d',]),
    #         as.vector(SBMD$G4_N['par.d',]), as.vector(SBMD$QE4_N['par.d',]), as.vector(SBMD$P4_N['par.d',]), as.vector(SBMD$L4_N['par.d',]),
    #         as.vector(SBMD$E4_LN['par.d',]), as.vector(SBMD$IE4_LN['par.d',]), as.vector(SBMD$H4_LN['par.d',]), as.vector(SBMD$LN4_LN['par.d',]),
    #         as.vector(SBMD$G4_LN['par.d',]), as.vector(SBMD$QE4_LN['par.d',]), as.vector(SBMD$P4_LN['par.d',]), as.vector(SBMD$L4_LN['par.d',])
    # ), file = paste0('SBMD_dataset',d,'_sim', s, '.txt'), ncolumns = 687+1, append = F)
    
    ### Fit PERT
    out.pert[d, ] <- as.numeric(fit_pert(SBMD))
    
  }
  out.all <- as.data.frame(out.pert)
  out.all$simID <- s
  colnames(out.all) <- c('BMD.min','BMD.mode','BMD.max','BMD.shape',
                         'BKG.min','BKG.mode','BKG.max','BKG.shape',
                         'MAXY.min','MAXY.mode','MAXY.max','MAXY.shape',
                         'pert.bmd.mu', 'pert.bmd.var', 'pert.bmd.var.discounted', 'pert.bmd.s.discounted',
                         'pert.bkg.mu', 'pert.bkg.var', 'pert.bkg.var.discounted', 'pert.bkg.s.discounted',
                         'pert.maxy.mu', 'pert.maxy.var', 'pert.maxy.var.discounted', 'pert.maxy.s.discounted',
                         'sim.ID')
  # write.csv(out.all, file = paste0('parEst_histData_sim', s, '.RData'), append = F, quote = F, col.names = T, row.names = F, sep = ',')
  
  ##--------------------------------------------------------
  ## Analyse new data w/ priors from each historical dataset
  ##--------------------------------------------------------
  
  step <- 3
  
  par_estimates <- out.all
  
  ### FITTED PERT
  
  output <- matrix(nrow = 3, ncol = 45)
  for(d in 1:3){
    ## Informative prior on background
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'BKG.shape'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'BKG.shape'])
    SBMD_bkg <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on max response
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'MAXY.shape'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'MAXY.shape'])
    SBMD_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                 nrchains = nrch, nriterations = nriter,
                                 warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on BMD
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                          shape.BMD = par_estimates[d, 'BMD.shape'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                            shape.BMD = par_estimates[d, 'BMD.shape'])
    SBMD_bmd <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on bkg and maxy
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'BKG.shape'],
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'MAXY.shape'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'BKG.shape'],
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'MAXY.shape'])
    SBMD_bkg_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                     nrchains = nrch, nriterations = nriter,
                                     warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on bkg, maxy and BMD
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'BKG.shape'],
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'MAXY.shape'],
                          prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                          shape.BMD = par_estimates[d, 'BMD.shape'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'BKG.shape'],
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'MAXY.shape'],
                            prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                            shape.BMD = par_estimates[d, 'BMD.shape'])
    SBMD_all <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    # Save credible intervals
    output[d, ] <- c(SBMD_bkg$MA_bridge_sampling, SBMD_maxy$MA_bridge_sampling, SBMD_bmd$MA_bridge_sampling, SBMD_bkg_maxy$MA_bridge_sampling,
                     SBMD_all$MA_bridge_sampling,
                     quantile(SBMD_bkg$bkg_post_bs, pvec), quantile(SBMD_maxy$bkg_post_bs, pvec), quantile(SBMD_bmd$bkg_post_bs, pvec),
                     quantile(SBMD_bkg_maxy$bkg_post_bs, pvec), quantile(SBMD_all$bkg_post_bs, pvec),
                     quantile(SBMD_bkg$maxy_post_bs, pvec), quantile(SBMD_maxy$maxy_post_bs, pvec), quantile(SBMD_bmd$maxy_post_bs, pvec),
                     quantile(SBMD_bkg_maxy$maxy_post_bs, pvec), quantile(SBMD_all$maxy_post_bs, pvec))
    
    ## TO DO: SAVE OUTPUT regarding coverage etc.
    
  }
  # Background
  bkg.result <- output[, c(16:30)]
  bkg_inf_bkg <- bkg.result[, c(1:3)]
  bkg_inf_maxy <- bkg.result[, c(4:6)]
  bkg_inf_bmd <- bkg.result[, c(7:9)]
  bkg_inf_bkgmaxy <- bkg.result[, c(10:12)]
  bkg_inf_all <- bkg.result[, c(13:15)]
  bkg_results <- rbind(bkg_inf_bkg, bkg_inf_maxy, bkg_inf_bmd, bkg_inf_bkgmaxy, bkg_inf_all)
  bkg_results <- as.data.frame(bkg_results)
  names(bkg_results) <- c('lower','median','upper')
  bkg_results$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  bkg_results$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  bkg_results$pert <- 'fitted'
  # Max response
  maxy.result <- output[, c(31:45)]
  maxy_inf_bkg <- maxy.result[, c(1:3)]
  maxy_inf_maxy <- maxy.result[, c(4:6)]
  maxy_inf_bmd <- maxy.result[, c(7:9)]
  maxy_inf_bkgmaxy <- maxy.result[, c(10:12)]
  maxy_inf_all <- maxy.result[, c(13:15)]
  maxy_results <- rbind(maxy_inf_bkg, maxy_inf_maxy, maxy_inf_bmd, maxy_inf_bkgmaxy, maxy_inf_all)
  maxy_results <- as.data.frame(maxy_results)
  names(maxy_results) <- c('lower','median','upper')
  maxy_results$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  maxy_results$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  maxy_results$pert <- 'fitted'
  # BMD
  bmd.result <- output[, c(1:15)]
  bmd_inf_bkg <- bmd.result[, c(1:3)]
  bmd_inf_maxy <- bmd.result[, c(4:6)]
  bmd_inf_bmd <- bmd.result[, c(7:9)]
  bmd_inf_bkgmaxy <- bmd.result[, c(10:12)]
  bmd_inf_all <- bmd.result[, c(13:15)]
  bmd_results <- rbind(bmd_inf_bkg, bmd_inf_maxy, bmd_inf_bmd, bmd_inf_bkgmaxy, bmd_inf_all)
  bmd_results <- as.data.frame(bmd_results)
  names(bmd_results) <- c('lower','median','upper')
  bmd_results$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  bmd_results$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  bmd_results$pert <- 'fitted'
  
  ### DISCOUNTED PERT
  
  step <- 4
  
  output2 <- matrix(nrow = 3, ncol = 45)
  for(d in 1:3){
    ## Informative prior on background
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'pert.bkg.s.discounted'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'pert.bkg.s.discounted'])
    SBMD_bkg <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on max response
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'pert.maxy.s.discounted'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'pert.maxy.s.discounted'])
    SBMD_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                 nrchains = nrch, nriterations = nriter,
                                 warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on BMD
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                          shape.BMD = par_estimates[d, 'pert.bmd.s.discounted'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                            shape.BMD = par_estimates[d, 'pert.bmd.s.discounted'])
    SBMD_bmd <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on bkg and maxy
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'pert.bkg.s.discounted'],
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'pert.maxy.s.discounted'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'pert.bkg.s.discounted'],
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'pert.maxy.s.discounted'])
    SBMD_bkg_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                     nrchains = nrch, nriterations = nriter,
                                     warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    ## Informative prior on bkg, maxy and BMD
    data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                          shape.a = par_estimates[d, 'pert.bkg.s.discounted'],
                          maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                          shape.c = par_estimates[d, 'pert.maxy.s.discounted'],
                          prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                          shape.BMD = par_estimates[d, 'pert.bmd.s.discounted'])
    data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                            bkg = c(par_estimates[d, 'BKG.min'], par_estimates[d, 'BKG.mode'], par_estimates[d, 'BKG.max']),
                            shape.a = par_estimates[d, 'pert.bkg.s.discounted'],
                            maxy = c(par_estimates[d, 'MAXY.min'], par_estimates[d, 'MAXY.mode'], par_estimates[d, 'MAXY.max']),
                            shape.c = par_estimates[d, 'pert.maxy.s.discounted'],
                            prior.BMD = c(par_estimates[d, 'BMD.min'], par_estimates[d, 'BMD.mode'], par_estimates[d, 'BMD.max']),
                            shape.BMD = par_estimates[d, 'pert.bmd.s.discounted'])
    SBMD_all <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                nrchains = nrch, nriterations = nriter,
                                warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
    
    # Save credible intervals
    output2[d, ] <- c(SBMD_bkg$MA_bridge_sampling, SBMD_maxy$MA_bridge_sampling, SBMD_bmd$MA_bridge_sampling, SBMD_bkg_maxy$MA_bridge_sampling,
                      SBMD_all$MA_bridge_sampling,
                      quantile(SBMD_bkg$bkg_post_bs, pvec), quantile(SBMD_maxy$bkg_post_bs, pvec), quantile(SBMD_bmd$bkg_post_bs, pvec),
                      quantile(SBMD_bkg_maxy$bkg_post_bs, pvec), quantile(SBMD_all$bkg_post_bs, pvec),
                      quantile(SBMD_bkg$maxy_post_bs, pvec), quantile(SBMD_maxy$maxy_post_bs, pvec), quantile(SBMD_bmd$maxy_post_bs, pvec),
                      quantile(SBMD_bkg_maxy$maxy_post_bs, pvec), quantile(SBMD_all$maxy_post_bs, pvec))
    
    ## TO DO: SAVE OUTPUT regarding coverage etc.
    
  }
  # Background
  bkg.result2 <- output2[, c(16:30)]
  bkg_inf_bkg <- bkg.result2[, c(1:3)]
  bkg_inf_maxy <- bkg.result2[, c(4:6)]
  bkg_inf_bmd <- bkg.result2[, c(7:9)]
  bkg_inf_bkgmaxy <- bkg.result2[, c(10:12)]
  bkg_inf_all <- bkg.result2[, c(13:15)]
  bkg_results2 <- rbind(bkg_inf_bkg, bkg_inf_maxy, bkg_inf_bmd, bkg_inf_bkgmaxy, bkg_inf_all)
  bkg_results2 <- as.data.frame(bkg_results2)
  names(bkg_results2) <- c('lower','median','upper')
  bkg_results2$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  bkg_results2$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  bkg_results2$pert <- 'discounted'
  # Max response
  maxy.result2 <- output2[, c(31:45)]
  maxy_inf_bkg <- maxy.result2[, c(1:3)]
  maxy_inf_maxy <- maxy.result2[, c(4:6)]
  maxy_inf_bmd <- maxy.result2[, c(7:9)]
  maxy_inf_bkgmaxy <- maxy.result2[, c(10:12)]
  maxy_inf_all <- maxy.result2[, c(13:15)]
  maxy_results2 <- rbind(maxy_inf_bkg, maxy_inf_maxy, maxy_inf_bmd, maxy_inf_bkgmaxy, maxy_inf_all)
  maxy_results2 <- as.data.frame(maxy_results2)
  names(maxy_results2) <- c('lower','median','upper')
  maxy_results2$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  maxy_results2$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  maxy_results2$pert <- 'discounted'
  # BMD
  bmd.result2 <- output2[, c(1:15)]
  bmd_inf_bkg <- bmd.result2[, c(1:3)]
  bmd_inf_maxy <- bmd.result2[, c(4:6)]
  bmd_inf_bmd <- bmd.result2[, c(7:9)]
  bmd_inf_bkgmaxy <- bmd.result2[, c(10:12)]
  bmd_inf_all <- bmd.result2[, c(13:15)]
  bmd_results2 <- rbind(bmd_inf_bkg, bmd_inf_maxy, bmd_inf_bmd, bmd_inf_bkgmaxy, bmd_inf_all)
  bmd_results2 <- as.data.frame(bmd_results2)
  names(bmd_results2) <- c('lower','median','upper')
  bmd_results2$inf <- c(rep('bkg', 3), rep('maxy', 3), rep('bmd', 3), rep('bkg_maxy', 3), rep('all', 3))
  bmd_results2$analysis <- rep(c('Dataset 1','Dataset 2','Dataset 3'), 5)
  bmd_results2$pert <- 'discounted'
  
  ## Combine results of separate analyses
  results_data <- rbind(bkg_results, bkg_results2,
                        maxy_results, maxy_results2,
                        bmd_results, bmd_results2)
  results_data$data <- NA
  save(results_data, file = paste0('/scratch/leuven/326/vsc32693/testEFSA/OUT/output_per_data_sim', s, '.RData'))
  
  ##--------------------------------------------------------
  ## Analyse new data w/ pooled historical datasets
  ##--------------------------------------------------------
  
  step <- 5
  
  ## Combine summary data
  pooled_data <- rbind(dataset_1, dataset_2, dataset_3)
  # save(pooled_data, file = paste0('./pooled_summ_data_sim', s, '.RData'))
  
  ## Analyse 'historical' dataset
  data_N <- PREP_DATA_N(data = pooled_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
  data_LN <- PREP_DATA_LN(data = pooled_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11')
  SBMD_poolHist <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                   nrchains = nrch, nriterations = nriter,
                                   warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## TO DO: SAVE OUTPUT for coverage etc.
  
  ## Estimate shape of PERT for each parameter
  pars.pooled <- fit_pert(SBMD_poolHist)
  # write.table(pars.pooled, file = paste0("./estPars_pooledData_sim", s, ".csv"), append = F, quote = F, col.names = F, row.names = F, sep = ",")
  
  ### Analyse new data: fitted PERT
  
  ## Informative prior on background
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['BKG.shape']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['BKG.shape']))
  SBMD_bkg <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on max response
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['MAXY.shape']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['MAXY.shape']))
  SBMD_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                               nrchains = nrch, nriterations = nriter,
                               warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on BMD
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                        shape.BMD = as.numeric(pars.pooled['BMD.shape']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                          shape.BMD = as.numeric(pars.pooled['BMD.shape']))
  SBMD_bmd <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on bkg and maxy
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['BKG.shape']),
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['MAXY.shape']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['BKG.shape']),
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['MAXY.shape']))
  SBMD_bkg_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                   nrchains = nrch, nriterations = nriter,
                                   warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on bkg, maxy and BMD
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['BKG.shape']),
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['MAXY.shape']),
                        prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                        shape.BMD = as.numeric(pars.pooled['BMD.shape']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['BKG.shape']),
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['MAXY.shape']),
                          prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                          shape.BMD = as.numeric(pars.pooled['BMD.shape']))
  SBMD_all <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  # Save credible intervals
  output.fitted <- c(SBMD_bkg$MA_bridge_sampling, SBMD_maxy$MA_bridge_sampling, SBMD_bmd$MA_bridge_sampling, SBMD_bkg_maxy$MA_bridge_sampling,
                     SBMD_all$MA_bridge_sampling,
                     quantile(SBMD_bkg$bkg_post_bs, pvec), quantile(SBMD_maxy$bkg_post_bs, pvec), quantile(SBMD_bmd$bkg_post_bs, pvec),
                     quantile(SBMD_bkg_maxy$bkg_post_bs, pvec), quantile(SBMD_all$bkg_post_bs, pvec),
                     quantile(SBMD_bkg$maxy_post_bs, pvec), quantile(SBMD_maxy$maxy_post_bs, pvec), quantile(SBMD_bmd$maxy_post_bs, pvec),
                     quantile(SBMD_bkg_maxy$maxy_post_bs, pvec), quantile(SBMD_all$maxy_post_bs, pvec))
  
  ## TO DO: SAVE OUTPUT regarding coverage etc.
  
  ### Analyse new data: discounted PERT
  step <- 6
  ## Informative prior on background
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']))
  SBMD_bkg <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on max response
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']))
  SBMD_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                               nrchains = nrch, nriterations = nriter,
                               warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on BMD
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                        shape.BMD = as.numeric(pars.pooled['pert.bmd.s.discounted']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                          shape.BMD = as.numeric(pars.pooled['pert.bmd.s.discounted']))
  SBMD_bmd <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on bkg and maxy
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']),
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']),
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']))
  SBMD_bkg_maxy <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                   nrchains = nrch, nriterations = nriter,
                                   warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  ## Informative prior on bkg, maxy and BMD
  data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                        bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                        shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']),
                        maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                        shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']),
                        prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                        shape.BMD = as.numeric(pars.pooled['pert.bmd.s.discounted']))
  data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                          bkg = c(as.numeric(pars.pooled['BKG.min']), as.numeric(pars.pooled['BKG.mode']), as.numeric(pars.pooled['BKG.max'])),
                          shape.a = as.numeric(pars.pooled['pert.bkg.s.discounted']),
                          maxy = c(as.numeric(pars.pooled['MAXY.min']), as.numeric(pars.pooled['MAXY.mode']), as.numeric(pars.pooled['MAXY.max'])),
                          shape.c = as.numeric(pars.pooled['pert.maxy.s.discounted']),
                          prior.BMD = c(as.numeric(pars.pooled['BMD.min']), as.numeric(pars.pooled['BMD.mode']), as.numeric(pars.pooled['BMD.max'])),
                          shape.BMD = as.numeric(pars.pooled['pert.bmd.s.discounted']))
  SBMD_all <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                              nrchains = nrch, nriterations = nriter,
                              warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
  
  # Save credible intervals
  output.discounted <- c(SBMD_bkg$MA_bridge_sampling, SBMD_maxy$MA_bridge_sampling, SBMD_bmd$MA_bridge_sampling, SBMD_bkg_maxy$MA_bridge_sampling,
                         SBMD_all$MA_bridge_sampling,
                         quantile(SBMD_bkg$bkg_post_bs, pvec), quantile(SBMD_maxy$bkg_post_bs, pvec), quantile(SBMD_bmd$bkg_post_bs, pvec),
                         quantile(SBMD_bkg_maxy$bkg_post_bs, pvec), quantile(SBMD_all$bkg_post_bs, pvec),
                         quantile(SBMD_bkg$maxy_post_bs, pvec), quantile(SBMD_maxy$maxy_post_bs, pvec), quantile(SBMD_bmd$maxy_post_bs, pvec),
                         quantile(SBMD_bkg_maxy$maxy_post_bs, pvec), quantile(SBMD_all$maxy_post_bs, pvec))
  
  ## TO DO: SAVE OUTPUT regarding coverage etc.
  
  output.pooled <- rbind(output.fitted, output.discounted)
  
  # Background
  bkg.result <- output.pooled[, c(16:30)]
  bkg_inf_bkg <- bkg.result[, c(1:3)]
  bkg_inf_maxy <- bkg.result[, c(4:6)]
  bkg_inf_bmd <- bkg.result[, c(7:9)]
  bkg_inf_bkgmaxy <- bkg.result[, c(10:12)]
  bkg_inf_all <- bkg.result[, c(13:15)]
  bkg_results_pooled <- rbind(bkg_inf_bkg, bkg_inf_maxy, bkg_inf_bmd, bkg_inf_bkgmaxy, bkg_inf_all)
  bkg_results_pooled <- as.data.frame(bkg_results_pooled)
  names(bkg_results_pooled) <- c('lower','median','upper')
  bkg_results_pooled$inf <- c(rep('bkg', 2), rep('maxy', 2), rep('bmd', 2), rep('bkg_maxy', 2), rep('all', 2))
  bkg_results_pooled$pert <- rep(c('fitted','discounted'), 5)
  # Max response
  maxy.result <- output.pooled[, c(31:45)]
  maxy_inf_bkg <- maxy.result[, c(1:3)]
  maxy_inf_maxy <- maxy.result[, c(4:6)]
  maxy_inf_bmd <- maxy.result[, c(7:9)]
  maxy_inf_bkgmaxy <- maxy.result[, c(10:12)]
  maxy_inf_all <- maxy.result[, c(13:15)]
  maxy_results_pooled <- rbind(maxy_inf_bkg, maxy_inf_maxy, maxy_inf_bmd, maxy_inf_bkgmaxy, maxy_inf_all)
  maxy_results_pooled <- as.data.frame(maxy_results_pooled)
  names(maxy_results_pooled) <- c('lower','median','upper')
  maxy_results_pooled$inf <- c(rep('bkg', 2), rep('maxy', 2), rep('bmd', 2), rep('bkg_maxy', 2), rep('all', 2))
  maxy_results_pooled$pert <- rep(c('fitted','discounted'), 5)
  # BMD
  bmd.result <- output.pooled[, c(1:15)]
  bmd_inf_bkg <- bmd.result[, c(1:3)]
  bmd_inf_maxy <- bmd.result[, c(4:6)]
  bmd_inf_bmd <- bmd.result[, c(7:9)]
  bmd_inf_bkgmaxy <- bmd.result[, c(10:12)]
  bmd_inf_all <- bmd.result[, c(13:15)]
  bmd_results_pooled <- rbind(bmd_inf_bkg, bmd_inf_maxy, bmd_inf_bmd, bmd_inf_bkgmaxy, bmd_inf_all)
  bmd_results_pooled <- as.data.frame(bmd_results_pooled)
  names(bmd_results_pooled) <- c('lower','median','upper')
  bmd_results_pooled$inf <- c(rep('bkg', 2), rep('maxy', 2), rep('bmd', 2), rep('bkg_maxy', 2), rep('all', 2))
  bmd_results_pooled$pert <- rep(c('fitted','discounted'), 5)
  
  results_pooled <- rbind(bkg_results_pooled, maxy_results_pooled, bmd_results_pooled)
  results_pooled$analysis <- 'Pooled'
  save(results_pooled, file = paste0('/scratch/leuven/326/vsc32693/testEFSA/OUT/output_pooled_sim', s, '.RData'))
  
  
  ##--------------------------------------------------------
  ## Analyse new data w/ cumulative historical datasets
  ##--------------------------------------------------------
  
  step <- 7
  
  for(f in c("fitted","discounted")){
    for(p in c("bkg","maxy","bmd","bkg_maxy","all")){
      
      # Cumulative prior
      meta <- fun_cum_meta(dataset_1, dataset_2, dataset_3, inf = p, pert = f) 
      prior.meta <- meta$pars
      dat.order <- c('dataset1','dataset2','dataset3','new')
      
      # Analyze 'new' dataset
      
      ## TO DO: SAVE OUTPUT for coverage etc.
      if(f == "fitted"){
        
        if(p == "bkg"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$BKG.shape)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$BKG.shape)
          
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          
          # save(SBMD_data_new, file = './cumulative_out/cum_meta_SBMD_data_new_bkg.RData')
          
          output.data.new.bkg <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bkg, file = './cumulative_out/output_cum_data_new_bkg.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bkg)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bkg.fitted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bkg.fitted, file = './output_cum_meta_fitted_bkg.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "maxy"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$MAXY.shape)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$MAXY.shape)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out/cum_meta_SBMD_data_new_maxy.RData')
          
          output.data.new.maxy <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.maxy, file = './cumulative_out/output_cum_data_new_maxy.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.maxy)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.maxy.fitted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.maxy.fitted, file = './output_cum_meta_fitted_maxy.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "bmd"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                shape.BMD = prior.meta$BMD.shape)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                  shape.BMD = prior.meta$BMD.shape)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out/cum_meta_SBMD_data_new_bmd.RData')
          
          output.data.new.bmd <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bmd, file = './cumulative_out/output_cum_data_new_bmd.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bmd)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bmd.fitted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bmd.fitted, file = './output_cum_meta_fitted_bmd.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "bkg_maxy"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$BKG.shape,
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$MAXY.shape)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$BKG.shape,
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$MAXY.shape)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out/cum_meta_SBMD_data_new_bkgmaxy.RData')
          
          output.data.new.bkgmaxy <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bkgmaxy, file = './cumulative_out/output_cum_data_new_bkgmaxy.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bkgmaxy)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bkgmaxy.fitted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bkgmaxy.fitted, file = './output_cum_meta_fitted_bkgmaxy.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "all"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$BKG.shape,
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$MAXY.shape,
                                prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                shape.BMD = prior.meta$BMD.shape)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$BKG.shape,
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$MAXY.shape,
                                  prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                  shape.BMD = prior.meta$BMD.shape)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out/cum_meta_SBMD_data_new_all.RData')
          
          output.data.new.all <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.all, file = './cumulative_out/output_cum_data_new_all.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.all)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.all.fitted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.all.fitted, file = './output_cum_meta_fitted_all.RData')
          
          rm(SBMD_data_new)
          
        }
        
      }else if(f == "discounted"){
        
        if(p == "bkg"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$pert.bkg.s.discounted)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$pert.bkg.s.discounted)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out_disc/cum_meta_SBMD_data_new_bkg.RData')
          
          output.data.new.bkg <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bkg, file = './cumulative_out_disc/output_cum_data_new_bkg.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bkg)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bkg.discounted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bkg.discounted, file = './output_cum_meta_discounted_bkg.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "maxy"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$pert.maxy.s.discounted)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$pert.maxy.s.discounted)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out_disc/cum_meta_SBMD_data_new_maxy.RData')
          
          output.data.new.maxy <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.maxy, file = './cumulative_out_disc/output_cum_data_new_maxy.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.maxy)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.maxy.discounted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.maxy.discounted, file = './output_cum_meta_discounted_maxy.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "bmd"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                shape.BMD = prior.meta$pert.bmd.s.discounted)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                  shape.BMD = prior.meta$pert.bmd.s.discounted)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out_disc/cum_meta_SBMD_data_new_bmd.RData')
          
          output.data.new.bmd <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bmd, file = './cumulative_out_disc/output_cum_data_new_bmd.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bmd)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bmd.discounted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bmd.discounted, file = './output_cum_meta_discounted_bmd.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "bkg_maxy"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$pert.bkg.s.discounted,
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$pert.maxy.s.discounted)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$pert.bkg.s.discounted,
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$pert.maxy.s.discounted)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out_disc/cum_meta_SBMD_data_new_bkgmaxy.RData')
          
          output.data.new.bkgmaxy <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.bkgmaxy, file = './cumulative_out_disc/output_cum_data_new_bkgmaxy.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.bkgmaxy)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.bkgmaxy.discounted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.bkgmaxy.discounted, file = './output_cum_meta_discounted_bkgmaxy.RData')
          
          rm(SBMD_data_new)
          
        }else if(p == "all"){
          
          data_N <- PREP_DATA_N(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                shape.a = prior.meta$pert.bkg.s.discounted,
                                maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                shape.c = prior.meta$pert.maxy.s.discounted,
                                prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                shape.BMD = prior.meta$pert.bmd.s.discounted)
          data_LN <- PREP_DATA_LN(data = new_data, sumstats = TRUE, sd = TRUE, extended = TRUE, q = q, prior.d = 'N11',
                                  bkg = c(prior.meta$BKG.min, prior.meta$BKG.mode, prior.meta$BKG.max),
                                  shape.a = prior.meta$pert.bkg.s.discounted,
                                  maxy = c(prior.meta$MAXY.min, prior.meta$MAXY.mode, prior.meta$MAXY.max),
                                  shape.c = prior.meta$pert.maxy.s.discounted,
                                  prior.BMD = c(prior.meta$BMD.min, prior.meta$BMD.mode, prior.meta$BMD.max),
                                  shape.BMD = prior.meta$pert.bmd.s.discounted)
          SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr, 
                                           nrchains = nrch, nriterations = nriter,
                                           warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
          n <- 1
          while(class(SBMD_data_new)[1] == 'try-error' && n < 6){
            SBMD_data_new <- try(sampling_MA(data.N = data_N, data.LN = data_LN, prior.weights = prior.weights, ndraws = ndr,
                                             nrchains = nrch, nriterations = nriter,
                                             warmup = wu, delta = dl, treedepth = trd, seed = sd, pvec = pvec))
            n <- n + 1
          }
          # save(SBMD_data_new, file = './cumulative_out_disc/cum_meta_SBMD_data_new_all.RData')
          
          output.data.new.all <- c(SBMD_data_new$MA_bridge_sampling, quantile(SBMD_data_new$bkg_post_bs, pvec), quantile(SBMD_data_new$maxy_post_bs, pvec))
          # save(output.data.new.all, file = './cumulative_out_disc/output_cum_data_new_all.RData')
          
          output.data.all <- rbind(meta$out, output.data.new.all)
          
          res.bkg <- output.data.all[, c(4:6)]
          res.maxy <- output.data.all[, c(7:9)]
          res.bmd <- output.data.all[, c(1:3)]
          res.bkg <- as.data.frame(res.bkg)
          names(res.bkg) <- c('lower', 'median', 'upper')
          res.bkg$data <- dat.order
          res.bkg$est <- rep('bkg', 4)
          res.maxy <- as.data.frame(res.maxy)
          names(res.maxy) <- c('lower', 'median', 'upper')
          res.maxy$data <- dat.order
          res.maxy$est <- rep('maxy', 4)
          res.bmd <- as.data.frame(res.bmd)
          names(res.bmd) <- c('lower', 'median', 'upper')
          res.bmd$data <- dat.order
          res.bmd$est <- rep('bmd', 4)
          
          res.all.all.discounted <- rbind(res.bkg, res.maxy, res.bmd)
          # save(res.all.all.discounted, file = './output_cum_meta_discounted_all.RData')
          
          rm(SBMD_data_new)
          
        }
      }
    }
  }
  
  res.all.bkg.fitted$pert <- 'fitted'
  res.all.bkg.discounted$pert <- 'discounted'
  res.all.bkg <- rbind(res.all.bkg.fitted, res.all.bkg.discounted)
  res.all.maxy.fitted$pert <- 'fitted'
  res.all.maxy.discounted$pert <- 'discounted'
  res.all.maxy <- rbind(res.all.maxy.fitted, res.all.maxy.discounted)
  res.all.bmd.fitted$pert <- 'fitted'
  res.all.bmd.discounted$pert <- 'discounted'
  res.all.bmd <- rbind(res.all.bmd.fitted, res.all.bmd.discounted)
  res.all.bkgmaxy.fitted$pert <- 'fitted'
  res.all.bkgmaxy.discounted$pert <- 'discounted'
  res.all.bkgmaxy <- rbind(res.all.bkgmaxy.fitted, res.all.bkgmaxy.discounted)
  res.all.all.fitted$pert <- 'fitted'
  res.all.all.discounted$pert <- 'discounted'
  res.all.all <- rbind(res.all.all.fitted, res.all.all.discounted)
  res.all.bkg$inf <- 'bkg'
  res.all.maxy$inf <- 'maxy'
  res.all.bmd$inf <- 'bmd'
  res.all.bkgmaxy$inf <- 'bkg_maxy'
  res.all.all$inf <- 'all'
  res.all <- rbind(res.all.bkg, res.all.maxy, res.all.bmd, res.all.bkgmaxy, res.all.all)
  
  res.all$analysis <- 'Cumulative'
  # res.all <- res.all[res.all$data == 'new',]
  # res.all <- res.all[,-4]
  save(res.all, file = paste0('/scratch/leuven/326/vsc32693/testEFSA/OUT/output_cumulative_sim', s, '.RData'))
  
  
  
}

stopCluster(cl)
