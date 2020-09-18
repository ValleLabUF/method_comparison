#######################
#### Run LDA Model ####
#######################

set.seed(123)

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
dat<- read.csv("CRW_MM_tsegs_weird.csv", as.is = T)  #mixed-membership sim
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
# write.csv(theta.estim_export, "theta_estim_weird.csv", row.names = F)

#read-in data
# theta.estim<- read.csv("theta_estim_weird.csv", as.is=T)
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
behav.res<- purrr::map(res, get_behav_hist, nburn = nburn, ngibbs = ngibbs,
                       nmaxclust = nmaxclust, var.names = c("Step Length", "Turning Angle")) %>% 
  purrr::map(., function(x) x[x$behav <= 3,])  #only select the top 3 behaviors
names(behav.res)<- names(dat.list)
behav.res_exp<- bind_rows(behav.res, .id = "id")

#export behav.res values
# write.csv(behav.res_exp, "CRW MM LDA Phi values_weird.csv", row.names = F)


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
behav.order<- list(c(1,3,2), c(1,3,2), c(2,3,1), c(2,1,3), c(1,3,2),
                   c(2,1,3), c(1,2,3), c(2,1,3), c(1,2,3), c(2,1,3),
                   c(1,2,3), c(2,1,3), c(1,2,3), c(1,3,2), c(1,3,2),
                   c(1,2,3), c(2,3,1), c(1,2,3), c(3,2,1), c(3,2,1))


names(behav.order)<- names(theta.estim)

behav.order_exp<- bind_rows(behav.order) %>% 
  t() %>% 
  data.frame() %>% 
  mutate(id = names(behav.order))
# write.csv(behav.order_exp, "CRW MM LDA behavior order_weird.csv", row.names = F)


#Change into long format
theta.estim.long<- list()
for (i in 1:length(theta.estim2)) {
  theta.estim.long[[i]]<- theta.estim2[[i]] %>% 
    pivot_longer(cols = c(-tseg, -time1, -date), values_to = "prop", names_to = "behavior") %>% 
    mutate_at("behavior", as.factor) %>% 
    mutate_at("behavior", ~recode(., 'X1' = behav.order[[i]][1],
                                  'X2' = behav.order[[i]][2], 'X3' = behav.order[[i]][3]))
}





##True proportions by simulation ID
true.behavior.long<- list()
for (i in 1:length(dat.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(dat.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav.coarse = rep(dat.list[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(dat.list[[i]]$track_length[1]),each = 3))
  
  true.behavior.long[[i]]$prop<- 0.1
  
  cond<- true.behavior.long[[i]][,"behav.coarse"]
  ind<- which(true.behavior.long[[i]][,"behavior"] == cond)
  
  true.behavior.long[[i]][ind,"prop"]<- 0.8
}
names(true.behavior.long)<- names(dat.list)


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
# write.csv(dat2, "Modeled MM Sim Tracks w Behav_weird.csv", row.names = F)  #for mixed-membership sim
# write.csv(time, "LDA_elapsed_time_weird.csv", row.names = F)  #units = min





##################
#### Figure 4 ####
##################

# Load breakpoints
bayes.brkpts<- read.csv("Bayesian_allbreakpts_weird.csv")
true.brkpts_weird<- read.csv("CRW_MM_sim_brkpts_weird.csv")


## Part a: segmentation heatmap
data<- dat.list$`2_3`[,c("id","SL","TA")]  #ID 2_3
nbins<- c(5,8)

## Function to transform bin numbers into binary matrix for heatmap
behav.seg.image=function(dat, nbins) {  #Transform single var vectors into pres/abs matrices for                                         #heatmap; nbins is vector of bins per param in order
  dat<- dat[,-1]  #remove id col
  behav.list<- map2(list(dat), nbins, ~matrix(0, nrow = nrow(.x), ncol = .y))
  for (i in 1:length(behav.list)) {
    for (j in 1:nrow(dat)){
      behav.list[[i]][j,dat[,i][j]]=1
    }
  }
  
  names(behav.list)<- names(dat)
  behav.list
}


behav.heat<- behav.seg.image(data, nbins)


SL<- data.frame(behav.heat$SL)
names(SL)<- 1:nbins[1]
SL<- SL %>% gather(key, value) %>% mutate(time=rep(dat.list[["2_3"]]$time1, times=nbins[1]),
                                          behav=rep("Step Length", nrow(data)*nbins[1]))

TA<- data.frame(behav.heat$TA)
names(TA)<- 1:nbins[2]
TA<- TA %>% gather(key, value) %>% mutate(time=rep(dat.list[["2_3"]]$time1, times=nbins[2]),
                                          behav=rep("Turning Angle", nrow(data)*nbins[2]))

behav.heat_long<- rbind(SL,TA)
behav.heat_long$value<- factor(behav.heat_long$value)
levels(behav.heat_long$value)<- c("Unused","Used")

breakpt<- bayes.brkpts %>% 
  filter(id == unique(data$id) & type == "Model") %>% 
  dplyr::select(brks) %>% 
  rename(breaks = brks)

true.brks<- data.frame(t(true.brkpts_weird[12,-1])) %>%
  drop_na() %>% 
  rename(breaks = X12) %>% 
  pull(breaks)

ticks.bottom<- data.frame(x=true.brks, y=0.5, xend=true.brks, yend=1.5)
ticks.top<- data.frame(x=rep(true.brks, 2), y=NA, xend=rep(true.brks, 2), yend=NA)
ticks.top$y<- rep(7.5, length(true.brks))
ticks.top$yend<- rep(8.5, length(true.brks))
ticks.top$behav<- rep(c("Step Length","Turning Angle"), each = length(true.brks))


p.seg<- ggplot(behav.heat_long, aes(x=time)) +
  facet_wrap(~behav, scales = 'free', nrow = 2) +
  scale_fill_viridis_d('') +
  scale_y_discrete(expand = c(0,0)) +
  geom_vline(data = breakpt, aes(xintercept = breaks - 0.5), color = viridis(n=9)[7],
             size = 1, alpha = 1) +
  geom_segment(data = ticks.top, aes(x = x, xend = xend, y = y, yend = yend),
               color = "black", size = 1, alpha = 1) +
  geom_segment(data = ticks.bottom, aes(x = x, xend = xend, y = y, yend = yend), color = "black",
               size = 1, alpha = 1) +
  labs(x = "\nTime", y = "\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 12, face = 'bold'),
        plot.title = element_text(size = 20, hjust = 0, vjust = -6),
        plot.margin = margin(0, 1, 0.5, 0.5, "cm"),
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.position = "top",
        legend.text = element_text(
          margin = margin(r = 15, unit = "pt")))




## Part b: behavior histogram
library(circular)
source('helper functions.R')

#calculate true proportions of SL and TA by behavior for all bins for ID 2_2
SL.params_weird<- data.frame(par1 = c(0.25, 2, exp(2)), par2 = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0), par2 = c(0.5, pi, 1))

true.b_weird<- extract.behav.props_weird(params = list(SL.params_weird, TA.params),
                                         lims = list(dist.bin.lims, angle.bin.lims),
                                         behav.names = c("Encamped","ARS","Transit"))

bayes.b<- behav.res[[12]]
bayes.b$behav<- bayes.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Encamped") %>%
  str_replace_all(., "3", "Transit") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.hist<- ggplot(true.b_weird, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = bayes.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
             color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")




##Part c: generate time series plots comparing behavior proportions
dat1<- theta.estim.long[[12]]
dat1$behavior<- dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

bayes.list<- df_to_list(dat2, "id")
true.behavior.long<- list()
for (i in 1:length(bayes.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(bayes.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav_coarse = rep(bayes.list[[i]]$behav_coarse,
                                                          each = 3),
                                       behav_fine = rep(bayes.list[[i]]$behav_fine,
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(bayes.list[[i]]$track_length[1]),each = 3))
  
  true.behavior.long[[i]]$prop<- 0.1
  
  cond<- true.behavior.long[[i]][,"behav_coarse"]
  ind<- which(true.behavior.long[[i]][,"behavior"] == cond)
  
  true.behavior.long[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long[[i]][cond2, "prop"]<- true.behavior.long[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long)<- names(bayes.list)

true.dat1<- true.behavior.long[[12]]
true.dat1$behavior<- true.dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

p.prop<- ggplot() +
  geom_path(data = true.dat1,
            aes(x=time1, y=prop, color = behavior),
            size = 1.5, linejoin = "round", lineend = "round") +
  scale_color_manual(values = c(viridis(n=20)[c(1,9)], "gold3"), guide=F) +
  new_scale_color() +
  geom_path(data = dat1,
            aes(x=time1-1, y=prop, color = behavior),
            size = 0.75, linejoin = "round", lineend = "round") +
  scale_color_manual(values = c(viridis(n=20)[c(7,13)], "gold2"), guide=F) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  facet_wrap(~behavior, nrow = 3)




## Make composite
library(gridExtra)

# png("Figure 4 (results from sim).png", width = 14, height = 10, units = "in", res = 330)
grid.arrange(p.seg, p.hist, p.prop, heights = c(0.2, 1, 0.1, 1),
             layout_matrix = rbind(c(NA, NA),
                                   c(1, 2),
                                   c(NA, NA),
                                   c(3, 3)))
# dev.off()
