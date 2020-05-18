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

source('helper functions.R')


#load and manipulate data
dat<- read.csv("CRW_MM_sim_multinom.csv", as.is = T)  #for mixed-membership sims
true.brkpts<- read.csv("CRW_MM_sim_brkpts.csv", as.is = T)
dat$dt<- 3600
# dat$id<- 1
dat.list<- df.to.list(dat=dat, ind = "id")
names(dat)[4:5]<- c("dist","rel.angle")

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
# track_length == 1000; takes 3.5 min for 40000 iterations
dat.res_1k<- behavior_segment(data = behav.list2[1:5], ngibbs = ngibbs, nbins = c(5,8),
                           alpha = alpha)
# track_length == 5000; takes 14 min for 40000 iterations
dat.res_5k<- behavior_segment(data = behav.list2[6:10], ngibbs = ngibbs, nbins = c(5,8),
                           alpha = alpha)
# track_length == 10000; takes 25 min for 40000 iterations
dat.res_10k<- behavior_segment(data = behav.list2[11:15], ngibbs = ngibbs, nbins = c(5,8),
                           alpha = alpha)
# track_length == 50000; takes 138 min for 40000 iterations
dat.res_50k<- behavior_segment(data = behav.list2[16:20], ngibbs = ngibbs, nbins = c(5,8),
                           alpha = alpha)
plan(sequential)  #closes background workers


## If sims analyzed separately, merge all runs together in single list
dat.res<- mapply(rbind, dat.res_1k, dat.res_5k, dat.res_10k, dat.res_50k, SIMPLIFY=FALSE)

## Reclassify and restructure data as needed
dat.res$nbrks[,2:ncol(dat.res$nbrks)]<- apply(dat.res$nbrks[,2:ncol(dat.res$nbrks)], 2,
                      function(x) as.numeric(as.character(x)))
dat.res$LML[,2:ncol(dat.res$LML)]<- apply(dat.res$LML[,2:ncol(dat.res$LML)], 2,
                    function(x) as.numeric(as.character(x)))
dat.res$brkpts<- dat.res$brkpts %>% split(., row(.)) %>% flatten()  #deals with dimensional list
names(dat.res$brkpts)<- dat.res$nbrks[,1]


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res$nbrks, type = "nbrks")
traceplot(data = dat.res$LML, type = "LML")


##Determine maximum likelihood (ML) for selecting breakpoints
ML<- getML(dat.res$LML, nburn = 20000)
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML)

## Heatmaps
plot.heatmap(data = behav.list[4], nbins = c(5,8), brkpts = brkpts, type = "behav", title = T,
             legend = T)




# Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:nrow(brkpts)) {
all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = brkpts[i,-1], true.brkpts = true.brkpts[i,-1],
                            acc.tol = 20, dup.tol = 20, miss.tol = 50)
}


## Calculate percentage of each measure of breakpoint accuracy for summary
brkpt.acc<- lapply(all.brkpts, function(x) x %>%
                        filter(type == "Model") %>%
                        group_by(acc) %>%
                        tally() %>%
                        mutate(freq = n/sum(n)))


# Calculate mean proportion (and SD) of different accuracy levels across all 5 sims per duration
tmp<- seq(5, 20, by = 5)
summ.stats<- list()

for (i in 1:length(tmp)) {
  summ.stats[[i]]<- brkpt.acc[(tmp[i]-4):tmp[i]] %>% map_dfr(`[`)
}

# % accurate breakpoints (accurate + duplicate)
map(summ.stats, . %>%
      group_by(acc) %>%
      summarise(mean = mean(freq), sd = sd(freq), n = sum(n)) %>%
      mutate(prop = n/sum(n)) %>% 
      filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
      summarise(prop = sum(prop)))  #and calculate accuracy across all sims combined

# number of breakpoints missed
map(all.brkpts, . %>% 
      filter(acc == "Missing") %>% 
      tally()) %>% 
  unlist() %>% 
  sum()

# total number of true breakpoints 
map(all.brkpts, . %>% 
      filter(type == "True") %>% 
      tally()) %>% 
  unlist() %>% 
  sum()





#Compare elapsed time
time<- dat.res$elapsed.time
time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))
time$time<- time$time %>% 
  as.character() %>% 
  str_remove_all(" mins") %>%
  as.numeric()


ggplot(time, aes(track_length, time)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Elapsed Time (min)") +
  theme_bw()



#assign time seg and make as DF
dat_out<- map(behav.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)

#export breakpoints for easier reference and elapsed time for method comparison
names(all.brkpts)<- names(dat.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')


#export results for hard-clustering and mixed-membership simulations
# write.csv(dat_out, "CRW_MM_tsegs.csv", row.names = F)
# write.csv(all.brkpts, "Bayesian_allbreakpts.csv", row.names = F)
# write.csv(time, "Bayesian_elapsed_time.csv", row.names = F)  #units = min
