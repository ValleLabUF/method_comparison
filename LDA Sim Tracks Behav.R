#######################
#### Run LDA Model ####
#######################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_LDA_behavior")


library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(ggnewscale)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

#get data
dat<- read.csv("CRW_MM_tsegs.csv", as.is = T)  #mixed-membership sim
dat.list<- df.to.list(dat, ind = "id")  #for later behavioral assignment
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

#run Gibbs sampler
obs.list<- df.to.list(obs, ind = "id")
res<- list()
elapsed.time<- vector()

for (i in 1:length(obs.list)) {
  start.time<- Sys.time()
  res[[i]]=LDA_behavior_gibbs(dat=obs.list[[i]], gamma1=gamma1, alpha=alpha,
                              ngibbs=ngibbs, nmaxclust=nmaxclust,
                              nburn=nburn, ndata.types=ndata.types)
  end.time<- Sys.time()
  elapsed.time[i]<- difftime(end.time, start.time, units = "min")
}

#Check traceplot of log likelihood
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
plot(res[[i]]$loglikel, type='l', main = paste("ID",names(obs.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Extract and plot proportions of behaviors per time segment
theta.post<- map(res, pluck, "theta") %>% 
  map(., function(x) x[(nburn+1):ngibbs,])  #extract samples from posterior
theta.estim<- theta.post %>% 
  map(., colMeans) %>% 
  map(., ~matrix(., ncol = nmaxclust)) #calc mean of posterior

#export theta.estim if need to re-analyze later
# theta.estim_export<- theta.estim %>% map(., as.data.frame) %>% bind_rows(., .id = 'id')
# write.csv(theta.estim_export, "theta_estim.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("theta_estim.csv", as.is=T)
# theta.estim<- theta.estim %>% group_split(., id, keep=F)



#boxplots
par(mfrow=c(2,2), ask = T)
for (i in 1:length(res)) {
  boxplot(theta.estim[[i]], xlab="Behavior", ylab="Proportion of Total Behavior",
          main = paste("ID",names(obs.list)[i]))
}
par(mfrow=c(1,1), ask=F)

#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
purrr::map(theta.estim, function(x) round(colSums(x)/nrow(x), digits = 3)) %>% 
  purrr::map(., ~sum(.[1:3])) #%>% 
  # unlist() %>% 
  # data.frame() %>% 
  # summarise(mean=mean(.), sd=sd(.))
## First 3 behaviors have mean of 97.5% and SD=0.04

## Viz histograms from model
behav.res<- purrr::map(res, get_behav_hist, dat_red = dat_red) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors


#Plot histograms of proportion data; order color scale from slow to fast
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
  ggplot(behav.res[[i]], aes(x = bin, y = prop, fill = as.factor(behav))) +
    geom_bar(stat = 'identity') +
    labs(x = "\nBin", y = "Proportion\n", title = names(obs.list)[i]) +
    theme_bw() +
    theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
          axis.text.x.bottom = element_text(size = 12),
          strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
    scale_fill_viridis_d(guide = F) +
    facet_grid(behav ~ param, scales = "free_x")
  )
}
par(ask=F)



## Viz behavior over time
#Assign behaviors (via theta) to each time segment
theta.estim<- map(theta.estim, . %>% 
      as.data.frame() %>%  #convert to DF
      dplyr::select(., 1:3) %>%   #only select 1st three behaviors (cols)
      mutate(row_sum = rowSums(.)) %>%  #calculate sum of these behavior proportions
      mutate_at(1:3, ~ ./row_sum) %>%  #normalize proportions for 3 behaviors
      dplyr::select(-row_sum) %>%  #remove row_sum col
      mutate(tseg = 1:nrow(.)) %>%  #add col for time segment
      rename('1' = V1, '2' = V2, '3' = V3) %>%  #rename cols as numbers
      dplyr::select(tseg, 1, 2, 3))  #reorder cols
names(theta.estim)<- names(obs.list)


#calc obs per tseg using SL bins (more reliable than TA)
nobs.list<- map(obs.list, . %>% 
                  mutate(n = rowSums(dplyr::select(., str_subset(names(.), pattern = "y1")))) %>%
                  dplyr::select(., -c(str_subset(names(.), pattern = "y")))) 


#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim2<- list()
for (i in 1:length(dat.list)) {
  theta.estim2[[i]]<- dat.list[[i]] %>% 
    drop_na(behav_fine) %>% 
    mutate(date=time1-1) %>% 
    aug_behav_df(dat = .,
                 theta.estim = theta.estim[[i]],
                 nobs = nobs.list[[i]])
}



#Need to manually inspect all histograms and assign proper order (slowest to fastest)
# behav.order<- list(c(3,1,2), c(2,1,3), c(3,2,1), c(2,1,3), c(3,2,1),
#                    c(2,3,1), c(1,3,2), c(2,1,3), c(2,3,1), c(2,3,1),
#                    c(2,3,1), c(2,1,3), c(2,1,3), c(1,2,3), c(1,3,2),
#                    c(3,1,2), c(2,1,3), c(2,3,1), c(3,2,1), c(3,1,2))
behav.order<- list(c(1,3,2), c(1,2,3), c(2,3,1), c(1,3,2), c(2,3,1),
                   c(1,3,2), c(2,3,1), c(2,3,1), c(2,3,1), c(1,2,3),
                   c(2,1,3), c(1,2,3), c(2,1,3), c(2,3,1), c(1,3,2),
                   c(1,3,2), c(1,2,3), c(2,3,1), c(2,3,1), c(2,3,1))

names(behav.order)<- names(theta.estim)


#Change into long format
theta.estim.long<- list()
for (i in 1:length(theta.estim2)) {
  theta.estim.long[[i]]<- theta.estim2[[i]] %>% 
    pivot_longer(cols = c(-tseg, -time1, -date), values_to = "prop", names_to = "behavior") %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}





#generate long form of true behavior for MIXED-MEMBERSHIP SIM
true.behavior.long<- map(dat.list, . %>% 
                           drop_na("behav_fine") %>% 
                           mutate(true.tseg = rep(1:(track_length[1]/100), each = 100)) %>% 
                           mutate_at("behav_fine", as.factor) %>% 
                           group_by(true.tseg, behav_fine) %>% 
                           count(behav_fine, .drop = FALSE) %>% 
                           mutate(prop = n/100) %>% 
                           rename(behavior = behav_fine) %>% 
                           map_df(., rep, 100) %>% 
                           arrange(true.tseg) %>% 
                           mutate(time1 = rep(1:(nrow(.)/3), each = 3)))


## Plot traces of true and modeled behavior proportions across time segments
par(ask=T)
for (i in 1:length(behav.res)) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = theta.estim.long[[i]],
                aes(x=date, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = true.behavior.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n", title = names(obs.list)[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)


#Plot scatterplot and compare to 1:1 line
# scatter.prop.comp<- data.frame(behavior = theta.estim.long$behavior,
#                                prop.estim = theta.estim.long$prop,
#                                prop.true = true.behavior.long$prop)
# 
# ggplot(scatter.prop.comp, aes(prop.estim, prop.true, color = behavior)) +
#   geom_point(size = 3, alpha = 0.5) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_bw() +
#   facet_wrap(~behavior)


#assign behavior from sim to data
dat2<- list()
for (i in 1:length(theta.estim.long)) {  #assign behaviors to all obs
  theta.estim.long[[i]]<- cbind(id = names(dat.list)[i], theta.estim.long[[i]])
  dat2[[i]]<- assign_behav(dat.list = dat.list[i], theta.estim.long = theta.estim.long[[i]],
                           behav.names = c("Encamped","ARS","Transit"))
}

dat2<- map_dfr(dat2, `[`)


## Plot tracks with modeled behaviors

ggplot(data = dat2 %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat2 %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")



## Viz elapsed time

time<- elapsed.time %>% 
  unlist() %>% 
  data.frame(time = .)

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))


ggplot(time, aes(track_length, time)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Elapsed Time (min)") +
  theme_bw()




#export results
# write.csv(dat2, "Modeled MM Sim Tracks w Behav.csv", row.names = F)  #for mixed-membership sim
# write.csv(time, "LDA_elapsed_time.csv", row.names = F)  #units = min