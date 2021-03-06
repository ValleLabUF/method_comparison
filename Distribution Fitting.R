
### Denis example

cutoff=c(-pi/2,0,pi/2)
bins=c(0.7,0.2,0.1,0)
plot(bins,type='h')

obj.function=function(param){
  mu1=param[1]
  sd1=param[2]
  bins.estim=rep(NA,4)
  bins.estim[1]=pnorm(cutoff[1],mean=mu1,sd=sd1)
  bins.estim[2]=pnorm(cutoff[2],mean=mu1,sd=sd1)-pnorm(cutoff[1],mean=mu1,sd=sd1)
  bins.estim[3]=pnorm(cutoff[3],mean=mu1,sd=sd1)-pnorm(cutoff[2],mean=mu1,sd=sd1)
  bins.estim[4]=pnorm(Inf,mean=mu1,sd=sd1)-pnorm(cutoff[3],mean=mu1,sd=sd1)
  mean(abs(bins-bins.estim))  #Mean Absolute Error
}

optim(c(1,.001),obj.function)






### Snail kite example
source('aux functions VM.R')  #source Denis' code that fixes pvonmises() and uses box-constrained                               #optimization

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")

set.seed(1)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(viridis)
library(circular)
library(wesanderson)



source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')





############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Snail Kite Gridded Data_TOHO_behav2.csv', header = T, sep = ',')
dat$date<- dat$date %>% as_datetime()
dat.list<- df.to.list(dat)  #for later behavioral assignment

nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- get.summary.stats_behav(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1 

#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

res=LDA_behavior_gibbs(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn, ndata.types=ndata.types)


#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(res = res, dat_red = dat_red)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors

#Plot histograms of frequency data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(param ~ behav, scales = "free_y")






### SL

#Specify the mass distributions per behavior for SL
behav.res.SL<- behav.res %>% filter(param == "SL")  #select only SL hists

#define bin number and limits for step lengths
max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins


#functions for determing probs and optimization
gamma.function=function(param){
  shape1=param[1]
  rate1=param[2]
  
  bins.estim=rep(NA,(length(dist.bin.lims) - 1))
  
  for (i in 2:length(dist.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)
    } else {
      bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)-
        pgamma(dist.bin.lims[i-1],shape=shape1,rate=rate1)
    }
  }
  
  bins.estim
}

gamma.function.optim=function(param, probs){
  shape1=param[1]
  rate1=param[2]
  
  bins.estim=rep(NA,(length(dist.bin.lims) - 1))
  
  for (i in 2:length(dist.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)
    } else {
      bins.estim[i-1]=pgamma(dist.bin.lims[i],shape=shape1,rate=rate1)-
        pgamma(dist.bin.lims[i-1],shape=shape1,rate=rate1)
    }
  }
  
  mean(abs(probs-bins.estim))  #Mean Absolute Error
}


weibull.function=function(param){
  shape1=param[1]
  scale1=param[2]
  
  bins.estim=rep(NA,(length(dist.bin.lims) - 1))
  
  for (i in 2:length(dist.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pweibull(dist.bin.lims[i],shape=shape1,scale=scale1)
    } else {
      bins.estim[i-1]=pweibull(dist.bin.lims[i],shape=shape1,scale=scale1)-
        pweibull(dist.bin.lims[i-1],shape=shape1,scale=scale1)
    }
  }
  
  bins.estim
}

weibull.function.optim=function(param, probs){
  shape1=param[1]
  scale1=param[2]
  
  bins.estim=rep(NA,(length(dist.bin.lims) - 1))
  
  for (i in 2:length(dist.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pweibull(dist.bin.lims[i],shape=shape1,scale=scale1)
    } else {
      bins.estim[i-1]=pweibull(dist.bin.lims[i],shape=shape1,scale=scale1)-
        pweibull(dist.bin.lims[i-1],shape=shape1,scale=scale1)
    }
  }
  
  mean(abs(probs-bins.estim))  #Mean Absolute Error
}


#Initialize optim functions using 100 different combinations of both params for gamma and weibull

store.SL.probs<- matrix(NA, 3*3*max(behav.res.SL$bin), 4)  #3 behaviors for 3 sources of 5 bins
colnames(store.SL.probs)<- c("prop", "bin", "source", "behav")
store.SL.probs[,2]<- rep(1:max(behav.res.SL$bin), times=3*3)
store.SL.probs[,3]<- rep(rep(c("Discretized","Gamma","Weibull"), each=max(behav.res.SL$bin)), 3)
store.SL.probs[,4]<- rep(1:3, each=3*max(behav.res.SL$bin))

#initial values for gamma distrib
param1<- vector()
param1[1]<- 0.0625
for (i in 2:10) {
  param1[i]<- 2*param1[i-1]
}
param1

param2<- vector()
param2[1]<- 0.0001
for (i in 2:10) {
  param2[i]<- 10*param2[i-1]
}
param2

gamma.params<- expand.grid(param1, param2)


#initial values for weibull distrib
param1<- vector()
param1[1]<- 0.0625
for (i in 2:10) {
  
  param1[i]<- 2*param1[i-1]
}
param1


param2<- vector()
param2[1]<- 25
for (i in 2:10) {
  
  param2[i]<- 2*param2[i-1]
}
param2

weibull.params<- expand.grid(param1, param2)


for (i in 1:length(unique(behav.res$behav))) {
  behav.prop<- behav.res.SL %>% filter(behav == i) %>% dplyr::select(prop)
  
  
  SL.gamma.res<- matrix(NA, nrow(gamma.params), 3)
  colnames(SL.gamma.res)<- c("param1","param2","value")
  for (j in 1:nrow(gamma.params)) {
  SL.gamma.fit<- optim(c(gamma.params[j,1], gamma.params[j,2]), gamma.function.optim,
                       probs = behav.prop$prop, method = "Nelder-Mead")
  SL.gamma.res[j,1:2]<- SL.gamma.fit$par
  SL.gamma.res[j,3]<- SL.gamma.fit$value
  }
  
  
  SL.weibull.res<- matrix(NA, nrow(weibull.params), 3)
  colnames(SL.weibull.res)<- c("param1","param2","value")
  for (k in 1:nrow(weibull.params)) {
  SL.weibull.fit<- optim(c(weibull.params[k,1], weibull.params[k,2]), weibull.function.optim,
                         probs = behav.prop$prop, method = "Nelder-Mead")
  SL.weibull.res[k,1:2]<- SL.weibull.fit$par
  SL.weibull.res[k,3]<- SL.weibull.fit$value
  }
  
  ind<- which(store.SL.probs[,4] == i)
  store.SL.probs[ind,1]<- c(behav.prop$prop,
                            gamma.function(SL.gamma.res[which.min(SL.gamma.res[,"value"]),
                                                        c("param1","param2")]),
                            weibull.function(SL.weibull.res[which.min(SL.weibull.res[,"value"]),
                                                            c("param1","param2")]))
}

store.SL.probs<- as.data.frame(store.SL.probs)
store.SL.probs$prop<- as.numeric(as.character(store.SL.probs$prop))
store.SL.probs$behav<- factor(store.SL.probs$behav, levels = c(1:3))
levels(store.SL.probs$behav)<- c("Encamped","Exploratory","Transit")


#Step Length comparison
ggplot(store.SL.probs, aes(x=bin, y=prop)) +
  geom_area(aes(fill=source, color = source, group = source),
            position = position_dodge(width = 0), stat="identity", alpha = 0.35) +
  labs(x="Bins", y="Proportion") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(guide = F) +
  theme_bw() +
  facet_wrap(~behav) +
  guides(fill = guide_legend(title = "Distribution", override.aes = list(alpha = 1))) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))






### TA

#Specify the mass distributions per behavior for TA
behav.res.TA<- behav.res %>% filter(param == "TA")  #select only TA hists

#define bin number and limits for turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins


# vm.function=function(param){
#   mu1=param[1]
#   kappa1=exp(param[2])
#   
#   bins.estim=rep(NA,(length(angle.bin.lims) - 1))
#   
#   for (i in 2:length(angle.bin.lims)) {
#     if (i-1 == 1) {
#       bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)
#     } else {
#       bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)-
#         pvonmises(angle.bin.lims[i-1],mu=mu1,kappa=kappa1)
#     }
#   }
#   
#   bins.estim
# }
# 
# vm.function.optim=function(param, probs){
#   mu1=param[1]
#   kappa1=exp(param[2])
#   
#   bins.estim=rep(NA,(length(angle.bin.lims) - 1))
#   
#   for (i in 2:length(angle.bin.lims)) {
#     if (i-1 == 1) {
#       bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)
#     } else {
#       bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)-
#         pvonmises(angle.bin.lims[i-1],mu=mu1,kappa=kappa1)
#     }
#   }
#   
#   mean(abs(probs-bins.estim))  #Mean Absolute Error
# }


wc.function=function(param){
  mu1=circular(param[1],type='angles',units='radians')
  rho1=exp(param[2])/(1+exp(param[2]))
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    pwrappedcauchy=integrate(dwrappedcauchy,
                             angle.bin.lims[i-1],
                             angle.bin.lims[i],
                             mu=mu1,
                             rho=rho1)
    bins.estim[i-1]=pwrappedcauchy$value
  }
  
  bins.estim
}

wc.function.optim=function(param, probs){  #for optimization
  bins.estim=wc.function(param)
  mean(abs(probs-bins.estim))  #Mean Absolute Error
}


#Initialize optim functions using 100 different combinations of both params for von mises and wrapped cauchy

store.TA.probs<- matrix(NA, 3*3*max(behav.res.TA$bin), 4)  #3 behaviors for 3 sources of 8 bins
colnames(store.TA.probs)<- c("prop", "bin", "source", "behav")
store.TA.probs[,2]<- rep(1:max(behav.res.TA$bin), times=3*3)
store.TA.probs[,3]<- rep(rep(c("Discretized","von Mises","wrapped Cauchy"),
                             each=max(behav.res.TA$bin)), 3)
store.TA.probs[,4]<- rep(1:3, each=3*max(behav.res.TA$bin))


#initial values for von mises distrib
param1<- seq(-pi, pi, length.out = 10)
param2<- c(seq(0.1, 0.9, length.out = 5),
           seq(1, 100, length.out = 5))

VM.init.params<- expand.grid(param1, param2)

#initial values for wrapped cauchy distrib
param1<- seq(-pi, pi, length.out = 10)
param2<- seq(0.1, 1, length.out = 10)

WC.init.params<- expand.grid(param1, param2)


for (i in 1:length(unique(behav.res$behav))) {
  print(i)
  
  behav.prop<- behav.res.TA %>% filter(behav == i) %>% dplyr::select(prop)
  
  TA.vm.res<- matrix(NA, nrow(VM.init.params), 3)
  colnames(TA.vm.res)<- c("param1","param2","MAE")
  for (j in 1:nrow(VM.init.params)) {
    param=c(VM.init.params[j,1], VM.init.params[j,2])
    TA.vm.fit<- optim(par=param, vm.function.optim,
                      probs = behav.prop$prop, method = "L-BFGS-B",
                      lower=c(-pi,0.001),upper=c(pi,10000)) #box-constrained optimization
    TA.vm.res[j,1:2]<- TA.vm.fit$par
    TA.vm.res[j,3]<- TA.vm.fit$value
  }
  
  TA.wc.res<- matrix(NA, nrow(WC.init.params), 3)
  colnames(TA.wc.res)<- c("param1","param2","MAE")
  for (k in 1:nrow(WC.init.params)) {
    param=c(WC.init.params[k,1], WC.init.params[k,2])
    TA.wc.fit<- optim(par=param, wc.function.optim,
                      probs = behav.prop$prop, method = "Nelder-Mead")
    TA.wc.res[k,1:2]<- TA.wc.fit$par
    TA.wc.res[k,3]<- TA.wc.fit$value
  }
  
  ind<- which(store.TA.probs[,4] == i)
  store.TA.probs[ind,1]<- c(behav.prop$prop,
                            vm.function(TA.vm.res[which.min(TA.vm.res[,"MAE"]),
                                                     c("param1","param2")]),
                            wc.function(TA.wc.res[which.min(TA.wc.res[,"MAE"]),
                                                     c("param1","param2")]))
}

store.TA.probs<- as.data.frame(store.TA.probs)
store.TA.probs$prop<- as.numeric(as.character(store.TA.probs$prop))
store.TA.probs$behav<- factor(store.TA.probs$behav, levels = c(1:3))
levels(store.TA.probs$behav)<- c("Encamped","Exploratory","Transit")


#Check that all distribs sum to 1
store.TA.probs %>% group_by(behav, source) %>% summarise(sum=sum(prop))

#Turning Angle comparison
ggplot(store.TA.probs, aes(x=bin, y=prop)) +
  geom_area(aes(fill=source, color = source, group = source),
            position = position_dodge(width = 0), stat="identity", alpha = 0.35) +
  labs(x="Bins", y="Proportion") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(guide = F) +
  theme_bw() +
  facet_wrap(~behav) +
  guides(fill = guide_legend(title = "Distribution", override.aes = list(alpha = 1))) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))







### Snow leopard example

############################
#### Load and Prep Data ####
############################

set.seed(1)

#get data
dat<- read.csv('Snow Leopard Data_behav.csv', header = T, sep = ',')
dat$date<- dat$date %>% as_datetime()
dat.list<- df.to.list(dat)  #for later behavioral assignment

nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- get.summary.stats_behav(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1 

#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

res=LDA_behavior_gibbs(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn, ndata.types=ndata.types)


#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(res = res, dat_red = dat_red)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors

#Plot histograms of frequency data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(param ~ behav, scales = "free_y")



### SL

#Specify the mass distributions per behavior for SL
behav.res.SL<- behav.res %>% filter(param == "SL")  #select only SL hists

#define bin number and limits for step lengths
max.dist=max(dat[dat$dt == 10800,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 10800,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins



#Initialize optim functions using 100 different combinations of both params for gamma and weibull

store.SL.probs<- matrix(NA, 3*3*max(behav.res.SL$bin), 4)  #3 behaviors for 3 sources of 5 bins
colnames(store.SL.probs)<- c("prop", "bin", "source", "behav")
store.SL.probs[,2]<- rep(1:max(behav.res.SL$bin), times=3*3)
store.SL.probs[,3]<- rep(rep(c("Discretized","Gamma","Weibull"), each=max(behav.res.SL$bin)), 3)
store.SL.probs[,4]<- rep(1:3, each=3*max(behav.res.SL$bin))

#initial values for gamma distrib
param1<- vector()
param1[1]<- 0.0625
for (i in 2:10) {
  param1[i]<- 2*param1[i-1]
}
param1

param2<- vector()
param2[1]<- 0.0001
for (i in 2:10) {
  param2[i]<- 10*param2[i-1]
}
param2

gamma.params<- expand.grid(param1, param2)


#initial values for weibull distrib
param1<- vector()
param1[1]<- 0.25
for (i in 2:10) {
  
  param1[i]<- 2*param1[i-1]
}
param1


param2<- vector()
param2[1]<- 0.0625
for (i in 2:10) {
  
  param2[i]<- 2*param2[i-1]
}
param2

weibull.params<- expand.grid(param1, param2)


for (i in 1:length(unique(behav.res$behav))) {
  behav.prop<- behav.res.SL %>% filter(behav == i) %>% dplyr::select(prop)
  
  
  SL.gamma.res<- matrix(NA, nrow(gamma.params), 3)
  colnames(SL.gamma.res)<- c("param1","param2","value")
  for (j in 1:nrow(gamma.params)) {
    SL.gamma.fit<- optim(c(gamma.params[j,1], gamma.params[j,2]), gamma.function.optim,
                         probs = behav.prop$prop, method = "Nelder-Mead")
    SL.gamma.res[j,1:2]<- SL.gamma.fit$par
    SL.gamma.res[j,3]<- SL.gamma.fit$value
  }
  
  
  SL.weibull.res<- matrix(NA, nrow(weibull.params), 3)
  colnames(SL.weibull.res)<- c("param1","param2","value")
  for (k in 1:nrow(weibull.params)) {
    SL.weibull.fit<- optim(c(weibull.params[k,1], weibull.params[k,2]), weibull.function.optim,
                           probs = behav.prop$prop, method = "Nelder-Mead")
    SL.weibull.res[k,1:2]<- SL.weibull.fit$par
    SL.weibull.res[k,3]<- SL.weibull.fit$value
  }
  
  ind<- which(store.SL.probs[,4] == i)
  store.SL.probs[ind,1]<- c(behav.prop$prop,
                            gamma.function(SL.gamma.res[which.min(SL.gamma.res[,"value"]),
                                                        c("param1","param2")]),
                            weibull.function(SL.weibull.res[which.min(SL.weibull.res[,"value"]),
                                                            c("param1","param2")]))
}

store.SL.probs<- as.data.frame(store.SL.probs)
store.SL.probs$prop<- as.numeric(as.character(store.SL.probs$prop))
store.SL.probs$behav<- factor(store.SL.probs$behav, levels = c(2,3,1))
levels(store.SL.probs$behav)<- c("Resting","Encamped","Exploratory")


#Step Length comparison
ggplot(store.SL.probs, aes(x=bin, y=prop)) +
  geom_area(aes(fill=source, color = source, group = source),
            position = position_dodge(width = 0), stat="identity", alpha = 0.35) +
  labs(x="Bins", y="Proportion") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(guide = F) +
  theme_bw() +
  facet_wrap(~behav) +
  guides(fill = guide_legend(title = "Distribution", override.aes = list(alpha = 1))) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))








### TA

#Specify the mass distributions per behavior for TA
behav.res.TA<- behav.res %>% filter(param == "TA")  #select only TA hists


#Initialize optim functions using 100 different combinations of both params for von mises and wrapped cauchy

store.TA.probs<- matrix(NA, 3*3*max(behav.res.TA$bin), 4)  #3 behaviors for 3 sources of 8 bins
colnames(store.TA.probs)<- c("prop", "bin", "source", "behav")
store.TA.probs[,2]<- rep(1:max(behav.res.TA$bin), times=3*3)
store.TA.probs[,3]<- rep(rep(c("Discretized","von Mises","wrapped Cauchy"),
                             each=max(behav.res.TA$bin)), 3)
store.TA.probs[,4]<- rep(1:3, each=3*max(behav.res.TA$bin))


#initial values for von mises distrib
param1<- seq(-pi, pi, length.out = 10)
param2<- c(seq(0.1, 0.9, length.out = 5),
           seq(1, 100, length.out = 5))

VM.init.params<- expand.grid(param1, param2)

#initial values for wrapped cauchy distrib
param1<- seq(-pi, pi, length.out = 10)
param2<- seq(0.1, 1, length.out = 10)

WC.init.params<- expand.grid(param1, param2)


for (i in 1:length(unique(behav.res$behav))) {
  print(i)
  
  behav.prop<- behav.res.TA %>% filter(behav == i) %>% dplyr::select(prop)
  
  TA.vm.res<- matrix(NA, nrow(VM.init.params), 3)
  colnames(TA.vm.res)<- c("param1","param2","MAE")
  for (j in 1:nrow(VM.init.params)) {
    param=c(VM.init.params[j,1], VM.init.params[j,2])
    TA.vm.fit<- optim(par=param, vm.function.optim,
                      probs = behav.prop$prop, method = "L-BFGS-B",
                      lower=c(-pi,0.001),upper=c(pi,10000)) #box-constrained optimization
    TA.vm.res[j,1:2]<- TA.vm.fit$par
    TA.vm.res[j,3]<- TA.vm.fit$value
  }
  
  TA.wc.res<- matrix(NA, nrow(WC.init.params), 3)
  colnames(TA.wc.res)<- c("param1","param2","MAE")
  for (k in 1:nrow(WC.init.params)) {
    param=c(WC.init.params[k,1], WC.init.params[k,2])
    TA.wc.fit<- optim(par=param, wc.function.optim,
                      probs = behav.prop$prop, method = "Nelder-Mead")
    TA.wc.res[k,1:2]<- TA.wc.fit$par
    TA.wc.res[k,3]<- TA.wc.fit$value
  }
  
  ind<- which(store.TA.probs[,4] == i)
  store.TA.probs[ind,1]<- c(behav.prop$prop,
                            vm.function(TA.vm.res[which.min(TA.vm.res[,"MAE"]),
                                                  c("param1","param2")]),
                            wc.function(TA.wc.res[which.min(TA.wc.res[,"MAE"]),
                                                  c("param1","param2")]))
}

store.TA.probs<- as.data.frame(store.TA.probs)
store.TA.probs$prop<- as.numeric(as.character(store.TA.probs$prop))
store.TA.probs$behav<- factor(store.TA.probs$behav, levels = c(2,3,1))
levels(store.TA.probs$behav)<- c("Resting","Encamped","Exploratory")


#Check that all distribs sum to 1
store.TA.probs %>% group_by(behav, source) %>% summarise(sum=sum(prop))


#Turning Angle comparison
ggplot(store.TA.probs, aes(x=bin, y=prop)) +
  geom_area(aes(fill=source, color = source, group = source),
            position = position_dodge(width = 0), stat="identity", alpha = 0.35) +
  labs(x="Bins", y="Proportion") +
  scale_fill_viridis_d() +
  scale_color_viridis_d(guide = F) +
  theme_bw() +
  facet_wrap(~behav) +
  guides(fill = guide_legend(title = "Distribution", override.aes = list(alpha = 1))) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))
