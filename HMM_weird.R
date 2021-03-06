###############################################################
###############################################################
#Movement ecology introduction
###############################################################
###############################################################

library(raster)        #for raster manipulation
library(rgdal)         #for shp data, projections; version 1.3-4 used
library(adehabitatLT)  #for interpreting trajectories
library(momentuHMM)
library(tidyverse)

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_behavior")
source('helper functions.R')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

#------------------------------#
#Hidden Markov model
#------------------------------#

d<- read.csv("CRW_MM_sim_weird.csv", as.is = T)
names(d)[1]<- "ID"
d.list<- df.to.list(d, ind = "ID")
d.list<- map(d.list, ~mutate(., time = seq(c(ISOdate(2020,6,12)), by = "hour",
                                           length.out = n())))  #create fake hourly times

#format data
move.d<- map(d.list, ~prepData(.[,1:3], type = "UTM", coordNames = c("x","y")))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(move.d)) {
  plot(y~x, data = move.d[[i]], main = paste("ID",names(d.list)[i]), type="o", col="gray65",
       asp=1, cex=0, pch=19, xlab="Longitude",ylab="Latitude")
}
par(mfrow=c(1,1), ask=F)



### Test with 2-4 behavioral states and then perform order selection via AIC/BIC
source('Iterate HMMs.R')

set.seed(3)
hmm.res_1k<- run.HMMs(move.d, 1:5)
names(hmm.res_1k)<- names(move.d)[1:5]

set.seed(3)
hmm.res_5k<- run.HMMs(move.d, 6:10)
hmm.res_5k<- hmm.res_5k[6:10]
names(hmm.res_5k)<- names(move.d)[6:10]

set.seed(1)
hmm.res_10k<- run.HMMs(move.d, 11:15)
hmm.res_10k<- hmm.res_10k[11:15]
names(hmm.res_10k)<- names(move.d)[11:15]

set.seed(1)
hmm.res_50k<- run.HMMs(move.d, 16:20)
hmm.res_50k<- hmm.res_50k[16:20]
names(hmm.res_50k)<- names(move.d)[16:20]


## Combine all model results
hmm.res<- c(hmm.res_1k, hmm.res_5k, hmm.res_10k, hmm.res_50k)



### Extract param values from fitted HMMs

hmm.params<- map(hmm.res, ~{map(.$models[2], getPar0) %>% 
    map(., ~pluck(., "Par")) %>% 
    flatten(.)
})

#step length params
hmm.step.params<- map(hmm.params, ~pluck(., "step")[1:6]) %>% 
  bind_rows() %>% 
  t() %>% 
  as.data.frame()
names(hmm.step.params)<- c("mean_1","mean_2","mean_3","sd_1","sd_2","sd_3")

#turning angle params
hmm.angle.params<- map(hmm.params, ~pluck(., "angle")) %>% 
  bind_rows() %>% 
  t() %>% 
  as.data.frame()
names(hmm.angle.params)<- c("mean_1","mean_2","mean_3","concentration_1",
                            "concentration_2","concentration_3")

#save model params
# write.csv(hmm.step.params, "HMM result step params_weird.csv")
# write.csv(hmm.angle.params, "HMM result angle params_weird.csv")




#Compare elapsed times
time<- map(hmm.res, . %>% pluck("elapsed.time")) %>% 
  map_dfr(., `[`) %>% 
  t() %>% 
  data.frame(stringsAsFactors=FALSE) %>% 
  rename(time = '.') %>% 
  mutate(., time_num = parse_number(time), time_unit = word(time, 2)) %>% 
  transmute(time = case_when(time_unit == "secs" ~ time_num/60,
                             time_unit == "mins" ~ time_num))

time$track_length<- rep(c('1k','5k','10k','50k'), each = 5) %>% 
  factor(., levels = c('1k','5k','10k','50k'))


ggplot(time, aes(track_length, time)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Elapsed Time (min)") +
  theme_bw()




## View output from each model and inspect pseudo-residuals

for (i in 1:length(hmm.res)) {
  plot(hmm.res[[i]]$models[[1]])
  plot(hmm.res[[i]]$models[[2]])
  plot(hmm.res[[i]]$models[[3]])
}

par(ask=T)
for (i in 1:length(hmm.res)) {
  plotPR(hmm.res[[i]]$models[[1]])
  plotPR(hmm.res[[i]]$models[[2]])
  plotPR(hmm.res[[i]]$models[[3]])
}
par(ask=F)


## Make inference via AIC
for (i in 1:length(hmm.res)) {
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  print(AIC(k.2, k.3, k.4))
  print(AICweights(k.2, k.3, k.4))  #4 states is far and away the best model
}

# Identify K per AIC
k.optim_AIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                3, 3, 3, 3, 4,
                4, 4, 3, 3, 3)

## Make inference via BIC
# BIC = -2*logL + p*log(T)

for (i in 1:length(hmm.res)) {
  
  k.2<- hmm.res[[i]]$models[[1]]  # 2 states
  k.3<- hmm.res[[i]]$models[[2]]  # 3 states  
  k.4<- hmm.res[[i]]$models[[3]]  # 4 states
  
  print(names(hmm.res)[i])
  #K=2
  print(2*k.2$mod$minimum) + (length(k.2$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=3
  print(2*k.3$mod$minimum) + (length(k.3$mod$estimate)*log(nrow(d.list[[i]])))  
  #K=4
  print(2*k.4$mod$minimum) + (length(k.4$mod$estimate)*log(nrow(d.list[[i]])))  
}

# Identify K per BIC
k.optim_BIC<- c(3, 3, 3, 3, 3,
                3, 3, 3, 3, 3,
                4, 3, 3, 3, 4,
                4, 4, 3, 3, 3)

# Identify K per AIC, BIC, and density distributions
k.optim<- c(3, 3, 3, 3, 3,
            3, 3, 3, 3, 3,
            3, 3, 3, 3, 3,
            3, 3, 3, 3, 3)

table(k.optim)/20
# AIC suggested 3 sims w/ 4 states
# BIC suggested 4 sims w/ 4 states
# holistic evaluation resulted in all sims using 3 states


### Sticking w/ K=3 for direct comparison

hmm.states<- hmm.res %>% 
  map(., ~pluck(., 1, 2)) %>% 
  map(viterbi) %>% 
  map2(d.list, ., ~cbind(.x, hmm.state = c(NA, .y[-length(.y)]))) %>% 
  bind_rows()




########################
### Model Validation ###
########################

### Coarse-scale

# Overall
hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == hmm.state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 78.1% to 88.4%

# For 'Resting' behavior
rest.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 1) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 1 & hmm.state == 1) %>% 
  tally() %>% 
  left_join(., rest.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 80.2% to 90.2%

# For 'ARS' behavior
ars.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 2) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 2 & hmm.state == 2) %>% 
  tally() %>% 
  left_join(., ars.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 72.5% to 88.0%

# For 'Transit' behavior
transit.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 3) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_coarse == 3 & hmm.state == 3) %>% 
  tally() %>% 
  left_join(., transit.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 74.0% to 90.0%



### Fine-scale

# Overall
hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == hmm.state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 90.3% to 94.1%

# For 'Resting' behavior
rest.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 1) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 1 & hmm.state == 1) %>% 
  tally() %>% 
  left_join(., rest.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 92.9% to 97.7%

# For 'ARS' behavior
ars.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 2) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 2 & hmm.state == 2) %>% 
  tally() %>% 
  left_join(., ars.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 76.5% to 89.4%

# For 'Transit' behavior
transit.size<- hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 3) %>% 
  tally()

hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == 3 & hmm.state == 3) %>% 
  tally() %>% 
  left_join(., transit.size, by = "ID") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 93.1% to 98.9%


# Export data and results
# write.csv(hmm.states, "HMM results_weird.csv", row.names = F)
# write.csv(time, "HMM_elapsed_time_weird.csv", row.names = F)  #units = min