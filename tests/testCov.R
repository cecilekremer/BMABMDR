

#####################################
### SIMULATED DATA WITH COVARIATE ###
#####################################

### CONTINUOUS
##############

n=100; dp=2; doses=c(rep(0,n),rep(c(0,dp^(-4),dp^(-3),dp^(-2),dp^(-1),1),n))
simsigma=1; q=0.1
pvec=c(0.05,0.5,0.95)
BMD.true=0.12 # flat
min.resp = 5.3
ctrue = 0.65
max.resp = min.resp*ctrue
par1=c((min.resp),(BMD.true),ctrue,1) # N data
set.seed(2)
y1=rnorm(rep(1,length(doses)),
         mean=DRM.E4_ND(par=par1,
                        x=doses,q=q),
         sd=simsigma)
min.resp = 5.3 # 10.6
max.resp = min.resp*ctrue
BMD.true = 0.08
simsigma = 1 #0.8
par2=c((min.resp),(BMD.true),ctrue,1.5) # N data
y2=rnorm(rep(1,length(doses)),
         mean=DRM.E4_ND(par=par2,
                        x=doses,q=q),
         sd=simsigma)
par(mfrow = c(1,1))
plot(doses, y1, xlab = 'Dose', ylab = 'Response', main = 'True model = E4_N',
     ylim = c(min(y1,y2)-0.5, max(y1,y2)+0.5))
points(doses, y2, pch = 2, col = 2)

# preparing the summary statistics
dose.a1=sort(unique(doses)); N1=length(dose.a1)
mean.a1=rep(NA,N1); sd.a1=rep(NA,N1); n.a1=rep(NA,N1); sex1=rep(NA,N1)
for (iu in (1:N1)){
  mean.a1[iu]=mean(y1[doses==dose.a1[iu]])
  sd.a1[iu]=sd(y1[doses==dose.a1[iu]])
  n.a1[iu]=sum(doses==dose.a1[iu])
  sex1[iu]=1
}
dose.a2=sort(unique(doses)); N2=length(dose.a2)
mean.a2=rep(NA,N2); sd.a2=rep(NA,N2); n.a2=rep(NA,N2); sex2=rep(NA,N2)
for (iu in (1:N2)){
  mean.a2[iu]=mean(y2[doses==dose.a2[iu]])
  sd.a2[iu]=sd(y2[doses==dose.a2[iu]])
  n.a2[iu]=sum(doses==dose.a2[iu])
  sex2[iu]=2
}
summ.data = data.frame(x = c(dose.a1, dose.a2), y = c(mean.a1, mean.a2), s = c(sd.a1, sd.a2),
                       n = c(n.a1, n.a2), sex = c(sex1, sex2))
lines(summ.data$x[summ.data$sex==1], summ.data$y[summ.data$sex==1], lty = 2, lwd = 2)
lines(summ.data$x[summ.data$sex==2], summ.data$y[summ.data$sex==2], col = 2, lwd = 2)
# legend('topright', c('a = 5.3, BMD = 0.12, d = 1, sigma = 1', 'a = 10.6, BMD = 0.08, d = 1.5, sigma = 0.8'), col = c(1,2), pch = c(1,2), bty = 'n')
legend('topright', c('a = 5.3, BMD = 0.12, d = 1, sigma = 1', 'a = 5.3, BMD = 0.08, d = 1.5, sigma = 1'), col = c(1,2), pch = c(1,2), bty = 'n')


prior.weights = rep(1,16)
FLBMD <- full.laplace_MA_Cov(
  data = summ.data,
  sumstats = TRUE,
  sd = TRUE, # option not used for Quantal data
  q = q,
  prior.d = 'N11'
)
FLBMD$MA
FLBMD$summary
# Plots
pt <- list()
k = 1
for (i in seq_along(get_models(type = "continuous"))){
  pt[[k]] <- basic.plot(FLBMD, model_name = get_models("continuous")[i], increasing = T)
  k = k+1
}
pt

### QUANTAL
###########

# a1 = 0.1; a2 = 0.1
# BMD1 = 0.25; BMD2 = 0.16
# d1 = d2 = 1
# q = 0.1
# source('~/BBMD/TRYOUT CODE/QUANTAL DR/simData.R')
# sim_data1 <- sim_DRQ_data(DRM = "Exponential", DRM_par = c(a1, BMD1, d1),
#                           # n = c(20, 20, 20, 20, 20),
#                           n = rep(100, 5),
#                           dose = c(0, 5, 15, 50, 100), q = q)
# sim_data1$sex = rep(1, 5)
# # sim_data2 <- sim_DRQ_data(DRM = "Exponential", DRM_par = c(a2, BMD2, d2),
# #                           # n = c(20, 20, 20, 20, 20),
# #                           n = rep(100, 5),
# #                           dose = c(0, 5, 15, 50, 100), q = q)
# sim_data2 = sim_data1
# sim_data2$sex = rep(2, 5)
# summ.data <- data.frame(dose = c(sim_data1$dose.a, sim_data2$dose.a),
#                         y = c(sim_data1$y.a, sim_data2$y.a),
#                         n = c(sim_data1$n.a, sim_data2$n.a),
#                         sex = c(sim_data1$sex, sim_data2$sex))
# par(mfrow=c(1,1))
# plot(summ.data$dose[summ.data$sex==1], summ.data$y[summ.data$sex==1]/summ.data$n[summ.data$sex==1], type = 'l', ylim = c(0,1))
# lines(summ.data$dose[summ.data$sex==2], summ.data$y[summ.data$sex==2]/summ.data$n[summ.data$sex==2], col = 2)
#
# modelFit <- full.laplace_MA_Q_Cov(
#   data = summ.data,
#   q = q
# )
# modelFit$summary
# modelFit$MA
#
# # Plots --> MODEL FIT?? covariate effects not reflected in fit
# for (i in seq_along(get_models(type = "quantal"))) {
#   basic.plotQ(modelFit, model_name = get_models(type = "quantal")[i])
# }


#####################################
### REAL DATA WITH COVARIATE      ###
#####################################

### CONTINUOUS
##############

dataDir <- '~/GitHub/BMABMDR/tests/'
load(file.path(dataDir, "das1.rda"))
data.test <- das1$data
ind.data <- data.frame(
  x = data.test$Dose,
  y = data.test$LDH,
  #          s = data.test$SD,
  #          n = data.test$N,
  cov = data.test$sex
)
ind.data <- ind.data[which(!is.na(ind.data$y)),]
q = 0.05
prior.weights = rep(1,16)
FLBMD <- full.laplace_MA_Cov(
  data = ind.data,
  sumstats = FALSE,
  sd = TRUE, # option not used for Quantal data
  q = q,
  prior.d = 'N11'
)
FLBMD$MA
FLBMD$summary
# Plots
# for (i in seq_along(get_models(type = "continuous"))){
#   basic.plot(FLBMD, model_name = get_models("continuous")[i], increasing = T)
# }
basic.plot(FLBMD, model_name = 'E4_N', increasing = T)

### QUANTAL
###########

data <- read.csv(file.path(dataDir, "example_quantal.csv"))

# Each covariate level should be present in each dose group
data.input <- data.frame(
  dose = rep(data$dose, 2),
  y = rep(data$response, 2),
  n = rep(data$size, 2),
  covariate = c(rep('Male', 5), rep('Female', 5))
)
# data.input$y[1:5] <- data.input$y[1:5] + 2

q = 0.1
# anydoseresponseQ(data.input$dose, data.input$y, data.input$n)

# Fit without covariate
data.input.Q <- PREP_DATA_QA(data = data.input, q = q, sumstats = T)
modelFit <- full.laplaceQ_MA(data.input.Q, prior.weights = rep(1,8), pvec = c(0.05,0.5,0.95))
modelFit$MA

# Fit with covariate
modelFit <- full.laplace_MA_Q_Cov(
  data = data.input,
  q = q
)
modelFit$summary
modelFit$MA

# Plots --> MODEL FIT?? covariate effects not reflected in fit
# for (i in seq_along(get_models(type = "quantal"))) {
#   basic.plotQ(modelFit, model_name = get_models(type = "quantal")[i])
# }
basic.plotQ(modelFit, 'E4_Q')

