################################
#### Run Segmentation Model ####
################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_behavior")


library(tidyverse)
library(tictoc)
library(furrr)
library(viridis)
library(lubridate)
library(ggforce)


source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')


setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

#load and manipulate data
dat<- read.csv("CRW_HC_sim.csv", as.is = T)  #for hard-clustering sim
# dat<- read.csv("CRW_MM_sim.csv", as.is = T)  #for mixed-membership sim
dat$dt<- 3600
dat$id<- 1
dat.list<- df.to.list(dat=dat, ind = "id")
names(dat)[3:4]<- c("dist","rel.angle")

#filter data for tstep of interest
behav.list<- behav.prep(dat=dat, tstep = 3600)  #add move params and filter by 3600 s interval

#define bin number and limits for turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

#define bin number and limits for step lengths
max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins


#Viz limits on continuous vars
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df, aes(x=dist)) +
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = dist.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nStep Length", y = "Density\n")

ggplot(behav.df, aes(x=rel.angle)) +
  geom_density(fill = "indianred") +
  geom_vline(xintercept = angle.bin.lims, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = "\nTurning Angle (rad)", y = "Density\n")



#assign bins to obs
behav.list<- map(behav.list, discrete_move_par, lims = list(dist.bin.lims, angle.bin.lims),
                 varIn = c("dist", "rel.angle"), varOut = c("SL", "TA"))
behav.list2<- lapply(behav.list, function(x) subset(x, select = c(id, SL, TA)))  #retain id and parameters on which to segment



#Viz discretization of params
behav.df2<- map_dfr(behav.list2, `[`)
behav.df2<- behav.df2 %>% gather(key, value, -id)

param.prop<- behav.df2 %>%
  group_by(key, value) %>%
  summarise(n=n()) %>%
  mutate(prop=n/nrow(behav.df)) %>%
  ungroup()  #if don't ungroup after grouping, ggforce won't work

param.prop<- param.prop[-c(6,15),]
param.prop[1:5, "value"]<- ((diff(dist.bin.lims)/2) + dist.bin.lims[1:5])
param.prop[6:13, "value"]<- (diff(angle.bin.lims)/2) + angle.bin.lims[1:8]


ggplot(data = param.prop %>% filter(key == "SL"), aes(value, prop)) +
  geom_bar(stat = "identity", width = (diff(dist.bin.lims)-0.025),
           fill = "lightblue", color = "black") +
  facet_zoom(xlim = c(0,5)) +
  labs(x = "Step Length", y = "Proportion") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


ggplot(data = param.prop %>% filter(key == "TA"), aes(value, prop)) +
  geom_bar(stat = "identity", fill = "indianred", color = "black") +
  labs(x = "Turning Angle (radians)", y = "Proportion") +
  theme_bw() +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



## Run RJMCMC

#prior
alpha = 1

ngibbs = 40000

plan(multisession)
dat.res<- behavior_segment(data = behav.list2, ngibbs = ngibbs, nbins = c(5,8), alpha = alpha)
plan(sequential)  #closes background workers
#takes 9 min for 40000 iterations


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(behav.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)

## Heatmaps
plot.heatmap(data = behav.list, nbins = c(5,8), brkpts = brkpts, dat.res = dat.res,
             type = "behav", title = F, legend = T)


## Compare True vs Modeled Breakpoints
model.brkpts<- na.omit(t(brkpts[-1])) %>% as.vector()
# true.brkpts<- which(diff(as.numeric(as.factor(dat$true.behav))) != 0) - 1  #for hard-clustering sim
true.brkpts<- which(diff(as.numeric(as.factor(dat$behav_coarse))) != 0) - 1  #for mixed-membership sim
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))

accuracy<- matrix(NA,length(model.brkpts),1)
for (i in 1:length(model.brkpts)) {  #assign brkpts as accurate or inaccurate
  
  tmp<- c(model.brkpts[i] - (20:0), model.brkpts[i] + (1:20)) %in% true.brkpts %>% sum()
  
  if (tmp == 0) {
    accuracy[i]<- "Inaccurate"
  } else {
    accuracy[i]<- "Accurate"
  }
}

if (sum(abs(diff(model.brkpts)) < 20) >= 0) {  #identify duplicate brkpts
  ind<- which(abs(diff(model.brkpts)) <= 20)
  ind<- sort(c(ind, ind+1))
}

ind.acc<- ind[which(accuracy[ind] == "Accurate")]
ind.inacc<- ind[which(accuracy[ind] == "Inaccurate")]
accuracy[ind.acc]<- "Accurate Duplicate"
accuracy[ind.inacc]<- "Inaccurate Duplicate"
accuracy<- c(rep("True",length(true.brkpts)), accuracy)


#identify missing breakpoints from model
status<- matrix(NA,length(true.brkpts),1)
for (i in 1:length(true.brkpts)) {
  
  tmp<- c(true.brkpts[i] - (50:0), true.brkpts[i] + (1:50)) %in% model.brkpts %>% sum()
  
  if (tmp == 0) {
    status[i]<- "Missing"
  } else {
    status[i]<- "Present"
  }
}

miss.ind<- which(status =="Missing")
status.miss<- data.frame(brks = true.brkpts[miss.ind], type = rep("Model", length(miss.ind)),
                         acc = rep("Missing", length(miss.ind)))


all.brkpts$acc<- accuracy
all.brkpts<- rbind(all.brkpts, status.miss)

#for hard-clustering
# ggplot(all.brkpts, aes(x=brks, y=type, color = acc)) +
#   geom_point(size=3) +
#   theme_bw() +
#   labs(x="Time", y="Type") +
#   theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10), legend.position = "top") +
#   scale_color_manual("Accuracy", values = c("forestgreen","lightgreen","firebrick","black"))

# #for mixed-membership
ggplot(all.brkpts, aes(x=brks, y=type, color = acc, shape = acc)) +
  geom_point(size=3) +
  theme_bw() +
  labs(x="Time", y="Type") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10),
        legend.position = "top") +
  scale_color_manual("Accuracy", values = c("forestgreen","lightgreen","firebrick","black")) +
  scale_shape_manual("Accuracy", values = c(16,16,4,16))


## Calculate percentage of each measure of breakpoint accuracy for summary

all.brkpts %>% filter(type == "Model") %>% group_by(acc) %>% tally() %>%
  mutate(freq = n/sum(n))



dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF


#export results for hard-clustering and mixed-membership simulations
# write.csv(dat_out, "CRW_HC_tsegs.csv", row.names = F)
# write.csv(dat_out, "CRW_MM_tsegs_multinom.csv", row.names = F)