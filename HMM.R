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

d<- read.csv("CRW_MM_sim_multinom.csv", as.is = T)
names(d)[1]<- "ID"
d.list<- df.to.list(d, ind = "ID")
d.list<- map(d.list, ~mutate(., time = seq(c(ISOdate(2020,5,4)), by = "hour",
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

set.seed(2)


hmm.res<- list()

for (j in 1:length(move.d)) {
  
  start.time<- Sys.time()
  
  # Empty list for order selection
  k.models<- list()
  
  
  ## K = 2
  print(paste("Simulation", j))
  
  stateNames <- c("Encamped","Exploratory")
  whichzero <- which(move.d[[j]]$step == 0)
  propzero <- length(whichzero)/nrow(move.d[[j]])
  zeromass0 <- c(propzero, 0.05)        #for zero distances by state
  
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.1, 5),
                     max = c(2, 20))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.1, 2),
                   max = c(1, 8))
  # Turning angle mean
  angleMean0 <- c(pi, 0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.4, 0.5),
                     max = c(0.99, 0.99))
  # Fit model
  if(propzero > 0) {  #don't include zero mass if no 0s present
    stepPar0 <- c(stepMean0, stepSD0, zeromass0)
  } else {
    stepPar0 <- c(stepMean0, stepSD0)
  }
  anglePar0 <- c(angleMean0,angleCon0)
  k.models[[1]] <- fitHMM(data = move.d[[j]], nbStates = 2, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          retryFits = 30)
  
  
  
  
  ## K = 3
  print(paste("Simulation", j))
  
  stateNames <- c("Resting","ARS","Transit")
  whichzero <- which(move.d[[j]]$step == 0)
  propzero <- length(whichzero)/nrow(move.d[[j]])
  zeromass0 <- c(propzero, 0.01, 0.05)        #for zero distances by state
  
  # Step length mean
  stepMean0 <- runif(3,
                     min = c(0.1, 0.5, 5),
                     max = c(2, 5, 20))
  # Step length standard deviation
  stepSD0 <- runif(3,
                   min = c(0.1, 1, 2),
                   max = c(1, 2, 8))
  # Turning angle mean
  angleMean0 <- c(pi, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(3,
                     min = c(0.4, 0.01, 0.5),
                     max = c(0.99, 0.3, 0.99))
  # Fit model
  if(propzero > 0) {  #don't include zero mass if no 0s present
    stepPar0 <- c(stepMean0, stepSD0, zeromass0)
  } else {
    stepPar0 <- c(stepMean0, stepSD0)
  }
  anglePar0 <- c(angleMean0,angleCon0)
  k.models[[2]] <- fitHMM(data = move.d[[j]], nbStates = 3, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          retryFits = 30)
  
  
  
  
  
  ## K = 4
  print(paste("Simulation", j))
  
  stateNames <- c("Resting","ARS","Exploratory","Transit")
  whichzero <- which(move.d[[j]]$step == 0)
  propzero <- length(whichzero)/nrow(move.d[[j]])
  zeromass0 <- c(propzero, 0.01, 0.03, 0.05)        #for zero distances by state
  
  # Step length mean
  stepMean0 <- runif(4,
                     min = c(0.1, 0.5, 3, 8),
                     max = c(2, 5, 10, 20))
  # Step length standard deviation
  stepSD0 <- runif(4,
                   min = c(0.1, 1, 2, 4),
                   max = c(1, 2, 4, 8))
  # Turning angle mean
  angleMean0 <- c(pi, pi, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(4,
                     min = c(0.4, 0.01, 0.01, 0.5),
                     max = c(0.99, 0.5, 0.5, 0.99))
  # Fit model
  if(propzero > 0) {  #don't include zero mass if no 0s present
    stepPar0 <- c(stepMean0, stepSD0, zeromass0)
  } else {
    stepPar0 <- c(stepMean0, stepSD0)
  }
  anglePar0 <- c(angleMean0,angleCon0)
  k.models[[3]] <- fitHMM(data = move.d[[j]], nbStates = 4, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          retryFits = 30)
  
  
  
  end.time<- Sys.time()
  elapsed.time<- difftime(end.time, start.time, units = "min")
  
  hmm.res[[j]]<- list(models = k.models, elapsed.time = elapsed.time)
}

names(hmm.res)<- names(move.d)

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
k.optim_AIC<- c(3, 4, 2, 3, 3,
                3, 3, 4, 3, 2,
                2, 4, 4, 3, 4,
                3, 4, 3, 3, 3)

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
k.optim_BIC<- c(3, 4, 2, 3, 4,
                3, 3, 4, 3, 2,
                2, 4, 4, 3, 4,
                3, 4, 3, 3, 3)


# Identify K per AIC, BIC, and density distributions
k.optim<- c(3, 4, 2, 3, 3,
            3, 3, 4, 3, 2,
            2, 4, 4, 3, 4,
            3, 4, 3, 3, 3)

table(k.optim)/20
# 55% (n=11) of simulations were estimated to exhibit 3 behaviors
# 15% (n=3) were estimated to exhibit 2 behaviors
# 30% (n=6) were estimated to exhibit 4 behaviors


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
# accuracy ranges from 6.8% to 89.6%

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
# accuracy ranges from 4.3% to 100.0%

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
# accuracy ranges from 8.4% to 100.0%

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
# accuracy ranges from 5.3% to 100.0%



### Fine-scale

# Overall
hmm.states %>% 
  group_by(track_length, ID) %>% 
  filter(., behav_fine == hmm.state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 4.8% to 92.6%

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
# accuracy ranges from 0.04% to 100.0%

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
# accuracy ranges from 1.7% to 100.0%

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
# accuracy ranges from 0.23% to 100.0%


# Export data and results
# write.csv(hmm.states, "HMM results.csv", row.names = F)
# write.csv(time, "HMM_elapsed_time.csv", row.names = F)  #units = min