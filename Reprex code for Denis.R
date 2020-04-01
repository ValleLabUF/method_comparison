set.seed(1)

library(tidyverse)
library(viridis)
library(circular)


#######################
## Load and viz data ##
#######################

behav.res<- read.csv('Denis snow leopard reprex data.csv', header = T, sep = ',')


#Plot histograms of frequency data
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(param ~ behav, scales = "free_y")





###############
## Functions ##
###############

vm.function=function(param){  #for calculating probability masses
  mu1=param[1]
  kappa1=exp(param[2])
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)
    } else {
      bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)-
        pvonmises(angle.bin.lims[i-1],mu=mu1,kappa=kappa1)
    }
  }
  
  bins.estim
}

#----------------------------------
vm.function.optim=function(param, probs){  #for optimization
  mu1=param[1]
  kappa1=exp(param[2])
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    if (i-1 == 1) {
      bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)
    } else {
      bins.estim[i-1]=pvonmises(angle.bin.lims[i],mu=mu1,kappa=kappa1)-
        pvonmises(angle.bin.lims[i-1],mu=mu1,kappa=kappa1)
    }
  }
  
  mean(abs(probs-bins.estim))  #Mean Absolute Error
}



###################################################
## Optimize functions based on von Mises distrib ##
###################################################

#Focus on turning angles
behav.res.TA<- behav.res %>% filter(param == "TA")  #select only TA hists

#define bin number and limits for turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins


#Initialize optim functions using 100 different combinations of both params for von mises and wrapped cauchy
param1<- seq(-pi, pi, length.out = 10)
param2<- seq(-3, 3, length.out = 10)

TA.init.params<- expand.grid(param1, param2)


#Only fit for behavior 3
behav.prop<- behav.res.TA %>% filter(behav == 3) %>% dplyr::select(prop)
  
#Run optimization
TA.vm.res<- matrix(NA, nrow(TA.init.params), 3)
colnames(TA.vm.res)<- c("param1","param2","MAE")
for (j in 1:nrow(TA.init.params)) {
  TA.vm.fit<- optim(c(TA.init.params[j,1], TA.init.params[j,2]), vm.function.optim,
                    probs = behav.prop$prop, method = "Nelder-Mead")
  TA.vm.res[j,1:2]<- TA.vm.fit$par
  TA.vm.res[j,3]<- TA.vm.fit$value
}


store.TA.probs<- matrix(NA, 1*2*max(behav.res.TA$bin), 4)  #1 behavior for 2 sources of 8 bins
colnames(store.TA.probs)<- c("prop", "bin", "source", "behav")
store.TA.probs[,2]<- rep(1:max(behav.res.TA$bin), times=1*2)
store.TA.probs[,3]<- rep(rep(c("Discretized","von Mises"),
                             each=max(behav.res.TA$bin)), 1)
store.TA.probs[,4]<- rep(3, each=2*max(behav.res.TA$bin))

#Select params w/ lowest MAE
store.TA.probs[,1]<- c(behav.prop$prop,
                          vm.function(TA.vm.res[which.min(TA.vm.res[,"MAE"]),
                                                c("param1","param2")]))


store.TA.probs<- as.data.frame(store.TA.probs)
store.TA.probs$prop<- as.numeric(as.character(store.TA.probs$prop))
store.TA.probs$behav<- factor(store.TA.probs$behav)



#################
## Viz results ##
#################

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
