#########################
### Method Comparison ###
#########################

library(tidyverse)
library(wesanderson)
library(lubridate)
library(cowplot)
library(viridis)
library(ggnewscale)

source('helper functions.R')


# Load elapsed time
seg.time<- read.csv("Bayesian_elapsed_time_weird.csv")
lda.time<- read.csv("LDA_elapsed_time_weird.csv")
bcpa.time<- read.csv("BCPA_elapsed_time_weird.csv")
hmm.time<- read.csv("HMM_elapsed_time_weird.csv")

# Load breakpoints
bayes.brkpts<- read.csv("Bayesian_allbreakpts_weird.csv")
bcpa.brkpts<- read.csv("BCPA_allbrkpts_weird.csv")

# Load results
# bayes.res<- read.csv("Modeled MM Sim Tracks w Behav.csv")
# hmm.res<- read.csv("HMM results.csv")
bayes.res_weird<- read.csv("Modeled MM Sim Tracks w Behav_weird.csv")
hmm.res_weird<- read.csv("HMM results_weird.csv")

# Load true breakpoints
# true.brkpts<- read.csv("CRW_MM_sim_brkpts.csv", as.is = T)
true.brkpts_weird<- read.csv("CRW_MM_sim_brkpts_weird.csv", as.is = T)

# Load helper functions
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_behavior")
source('helper functions.R')
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

############################
### Compare Elapsed Time ###
############################

#Add times together for segmentation and LDA model
bayes.time<- data.frame(time = seg.time$time + lda.time$time,
                        track_length = seg.time$track_length)

time<- rbind(bayes.time, bcpa.time, hmm.time)
time$method<- rep(c("Bayesian", "BCPA", "HMM"), each = 20)
time$track_length<- time$track_length %>% 
  factor(., levels = c('1k','5k','10k','50k'))


p.time<- ggplot(time, aes(track_length, time, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "Elapsed Time (min)\n") +
  scale_x_discrete(labels = c(1000,5000,10000,50000)) +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,3,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,3,5)]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = c(0.15,0.85),
        legend.background = element_blank(),
        legend.text = element_text(size = 12))


###########################
### Compare Breakpoints ###
###########################

bayes.brkpts$method<- rep("Bayesian", nrow(bayes.brkpts))
bcpa.brkpts$method<- rep("BCPA", nrow(bcpa.brkpts))
all.brkpts<- rbind(bayes.brkpts, bcpa.brkpts)


# brkpt.acc<- all.brkpts %>% 
#   group_by(method, id, acc) %>% 
#   filter(type == "Model") %>% 
#   tally() %>% 
#   mutate(freq = n/sum(n)) %>% 
#   mutate(track_length = case_when(str_detect(id, "_1") ~ "1k",
#                                   str_detect(id, "_2") ~ "5k",
#                                   str_detect(id, "_3") ~ "10k",
#                                   str_detect(id, "_4") ~ "50k")) %>% 
#   summarise(mean = mean(freq), sd = sd(freq), n = sum(n)) %>%
#   mutate(prop = n/sum(n)) %>% 
#   filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
#   summarise(prop = sum(prop)) %>%  #and calculate accuracy across all sims combined
#   ungroup()

brkpt.acc<- all.brkpts %>% 
  group_by(method, id, acc) %>% 
  filter(type == "Model") %>% 
  tally() %>% 
  mutate(freq = n/sum(n)) %>% 
  mutate(track_length = case_when(str_detect(id, "_1") ~ "1000",
                                  str_detect(id, "_2") ~ "5000",
                                  str_detect(id, "_3") ~ "10000",
                                  str_detect(id, "_4") ~ "50000")) %>% 
  # summarise(mean = mean(freq), sd = sd(freq), n = sum(n)) %>%
  group_by(method, track_length, id) %>% 
  filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.acc$track_length<- brkpt.acc$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))


#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
p.brk<- ggplot(brkpt.acc, aes(track_length, freq, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "Proportion of Accurate Breakpoints\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,3)], guide = F) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,3)], guide = F) +
  ylim(0,1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.grid = element_blank())


#calc mean accuracy per track_length and method
brkpt.acc %>% 
  group_by(method, track_length) %>% 
  summarise(mean=mean(freq))

# number of breakpoints missed
all.brkpts %>% 
  filter(acc == "Missing") %>% 
  group_by(method) %>% 
  tally()

# total number of true breakpoints 
all.brkpts %>% 
  filter(type == "True") %>% 
  group_by(method) %>% 
  tally()



#####################################################################
### Compare Accuracy of Behavior Estimates (Gamma/Wrapped Cauchy) ###
#####################################################################

# Assign identifiers by method and make consistent behavior colname
bayes.res$method<- rep("Bayesian", nrow(bayes.res))
hmm.res$method<- rep("HMM", nrow(hmm.res))

bayes.res<- bayes.res %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = c("Encamped","ARS","Transit"))) %>%
  mutate_at("state", as.numeric)
hmm.res<- hmm.res %>% rename(state = hmm.state, id = ID)


# Modify hmm.res to be same as bayes.res format
# calc proportions of behaviors by true time segment and then identify dominant behavior
hmm.res2<- hmm.res %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  dplyr::select(-c(time2, n)) %>% 
  pivot_wider(names_from = behavior, values_from = prop) %>% 
  mutate(track_length = max(time1)) %>% 
  ungroup() %>% 
  mutate(state = apply(.[,4:6], 1, which.max)) %>% 
  mutate_at("id", as.character) %>% 
  arrange(track_length, id) %>% 
  dplyr::select(-track_length)

hmm.res3<- hmm.res %>% 
  dplyr::select(-state) %>% 
  df.to.list("id") %>% 
  map2(.,
       df.to.list(hmm.res2, "id") %>% map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y)) %>% 
  bind_rows()
  

# Combine all datasets
res<- rbind(bayes.res[,c("id","behav_fine","behav_coarse","track_length","state","method")],
            hmm.res3[,c("id","behav_fine","behav_coarse","track_length","state","method")])



## Overall

#Coarse-scale behavior
res %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))

summ.stats_coarse<- res %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  ungroup()

summ.stats_coarse$track_length<- summ.stats_coarse$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))

p.coarse<- ggplot(summ.stats_coarse, aes(track_length, acc, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  ylim(0,1) +
  labs(x="", y = "Accuracy of Behavior Estimates\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




#### Compare Accuracy of Bayesian and HMM Proportion Estimates ####

##True proportions by simulation ID
bayes.list<- df.to.list(bayes.res, "id")

true.behavior.long<- list()
for (i in 1:length(bayes.list)) {
  true.behavior.long[[i]]<- data.frame(true.tseg = rep(1:(bayes.list[[i]]$track_length[1]/100),
                                                       each = 300),
                                       behav_coarse = rep(bayes.list[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behav_fine = rep(bayes.list[[i]]$behav_fine[-1],
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

## True proportions for HMMs (from time segments using true breakpoints)
hmm.props<- hmm.res %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  ungroup()
  
par(ask=T)
for (i in 1:length(unique(as.character(hmm.props$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data = hmm.props %>% filter(id == unique(as.character(hmm.res$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = true.behavior.long[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(hmm.res$id))[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)





## True proportions for Bayesian model (from modeled time segments)
bayes.props<- bayes.res %>% 
  rename(Encamped = X1, ARS = X2, Transit = X3) %>% 
  drop_na() %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",
               values_to = "prop") %>% 
  dplyr::select(id, tseg, behavior, prop, time1) %>% 
  mutate_at("behavior", ~recode(., 'Encamped' = 1, 'ARS' = 2, 'Transit' = 3))
bayes.props$time1<- bayes.props$time1 - 1


## Calculate RMSE
true.behavior<- true.behavior.long %>% bind_rows(.id = "id")

hmm.rmse<- vector()
for (i in 1:length(unique(hmm.res$id))) {
  ind<- unique(as.character(hmm.res$id))[i]
  
  hmm.rmse[i]<- sqrt(sum((hmm.props[hmm.props$id == ind, "prop"] - 
    true.behavior[true.behavior$id == ind, "prop"])^2) / nrow(hmm.props[hmm.props$id == ind,]))
}


bayes.rmse<- vector()
for (i in 1:length(unique(bayes.res$id))) {
  ind<- unique(as.character(bayes.res$id))[i]
  
  bayes.rmse[i]<- sqrt(sum((bayes.props[bayes.props$id == ind, "prop"] - 
                      true.behavior[true.behavior$id == ind, "prop"])^2) /
    nrow(bayes.props[bayes.props$id == ind,]))
}


rmse.df<- data.frame(id = rep(unique(as.character(hmm.res$id)), 2),
                 track_length = factor(rep(rep(c(1000,5000,10000,50000), each = 5), 2),
                                       levels = c("1000","5000","10000","50000")),
                 rmse = c(bayes.rmse, hmm.rmse),
                 method = rep(c("Bayesian","HMM"), each = 20))



#Plot results

p.rmse<- ggplot(rmse.df, aes(track_length, rmse, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "RMSE\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  ylim(0, 0.175) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




######################################################################################
### Compare Accuracy of Behavior Estimates (Truncated & Log Normal/Wrapped Cauchy) ###
######################################################################################


# Assign identifiers by method and make consistent behavior colname
bayes.res_weird$method<- rep("Bayesian", nrow(bayes.res_weird))
hmm.res_weird$method<- rep("HMM", nrow(hmm.res_weird))

bayes.res_weird<- bayes.res_weird %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = c("Encamped","ARS","Transit"))) %>%
  mutate_at("state", as.numeric)
hmm.res_weird<- hmm.res_weird %>% rename(state = hmm.state, id = ID)


# Modify hmm.res_weird to be same as bayes.res_weird format
# calc proportions of behaviors by true time segment and then identify dominant behavior
hmm.res_weird2<- hmm.res_weird %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  dplyr::select(-c(time2, n)) %>% 
  pivot_wider(names_from = behavior, values_from = prop) %>% 
  mutate(track_length = max(time1)) %>% 
  ungroup() %>% 
  mutate(state = apply(.[,4:6], 1, which.max)) %>% 
  mutate_at("id", as.character) %>% 
  arrange(track_length, id) %>% 
  dplyr::select(-track_length)

hmm.res_weird3<- hmm.res_weird %>% 
  dplyr::select(-state) %>% 
  df.to.list("id") %>% 
  map2(.,
       df.to.list(hmm.res_weird2, "id") %>% map(., ~rbind(c(unique(.$id), rep(NA, 7)), .)),
       ~cbind(.x, .y[,-1])) %>% 
  bind_rows()


# Combine all datasets
res_weird<- rbind(bayes.res_weird[,c("id","behav_fine","behav_coarse","track_length","state",
                                     "method")],
            hmm.res_weird3[,c("id","behav_fine","behav_coarse","track_length","state","method")])



## Overall

#Coarse-scale behavior
res_weird %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))

summ.stats_coarse_weird<- res_weird %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_coarse == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  ungroup()

summ.stats_coarse_weird$track_length<- summ.stats_coarse_weird$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))

p.coarse_weird<- ggplot(summ.stats_coarse_weird, aes(track_length,acc,fill = method,
                                                     color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  ylim(0,1) +
  labs(x="\nTrack Length (observations)", y = "Accuracy of Behavior Estimates\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




#### Compare Accuracy of Bayesian and HMM Proportion Estimates ####

##True proportions by simulation ID
bayes.list_weird<- df.to.list(bayes.res_weird, "id")

true.behavior.long_weird<- list()
for (i in 1:length(bayes.list_weird)) {
  true.behavior.long_weird[[i]]<- 
    data.frame(true.tseg = rep(1:(bayes.list_weird[[i]]$track_length[1]/100), each = 300),
                                       behav_coarse = rep(bayes.list_weird[[i]]$behav_coarse[-1],
                                                          each = 3),
                                       behav_fine = rep(bayes.list_weird[[i]]$behav_fine[-1],
                                                        each = 3),
                                       behavior = rep(1:3, 1000),
                                       time1 = rep(1:(bayes.list_weird[[i]]$track_length[1]),
                                                   each = 3))
  
  true.behavior.long_weird[[i]]$prop<- 0.1
  
  cond<- true.behavior.long_weird[[i]][,"behav_coarse"]
  ind<- which(true.behavior.long_weird[[i]][,"behavior"] == cond)
  
  true.behavior.long_weird[[i]][ind,"prop"]<- 0.8
  
  #add 0 or 1 for pure segments
  ind1<- true.behavior.long_weird[[i]] %>% 
    drop_na() %>% 
    group_by(true.tseg, behav_fine) %>% 
    tally() %>% 
    mutate(prop.true = n/sum(n))
  
  ind2<- ind1[which(ind1$prop.true == 1),]
  
  for (j in 1:nrow(ind2)) {
    cond2<- which(true.behavior.long_weird[[i]]$true.tseg == as.numeric(ind2[j,"true.tseg"]))
    true.behavior.long_weird[[i]][cond2, "prop"]<- true.behavior.long_weird[[i]][cond2,] %>% 
      mutate_at("prop", ~case_when(behavior == as.numeric(ind2[j,"behav_fine"]) ~ 1,
                                   behavior != as.numeric(ind2[j,"behav_fine"]) ~ 0)) %>% 
      dplyr::pull(prop)
  }
}
names(true.behavior.long_weird)<- names(bayes.list_weird)

## True proportions for HMMs (from time segments using true breakpoints)
hmm.props_weird<- hmm.res_weird %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg_hmm, brkpts = true.brkpts_weird) %>% 
  bind_rows() %>% 
  drop_na() %>%
  mutate_at("state", as.factor) %>%
  group_by(id, tseg, state) %>% 
  count(state, .drop = FALSE) %>% 
  group_by(id, tseg) %>% 
  mutate(prop = n/sum(n)) %>% 
  rename(behavior = state) %>% 
  uncount(sum(n), .id = "time2") %>%
  arrange(id, tseg, time2) %>%
  group_by(id) %>%
  mutate(time1 = rep(1:(n()/3), each = 3)) %>% 
  ungroup()

par(ask=T)
for (i in 1:length(unique(as.character(hmm.props_weird$id)))) {
  print(
    #Plot overlapping traces
    ggplot() +
      geom_line(data=hmm.props_weird %>% filter(id == unique(as.character(hmm.res_weird$id))[i]),
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 1) +
      scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
      new_scale_color() +
      geom_line(data = true.behavior.long_weird[[i]],
                aes(x=time1, y=prop, color = as.character(behavior)),
                size = 0.55) +
      scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
      labs(x = "\nObservation", y = "Proportion of Behavior\n",
           title = unique(as.character(hmm.res_weird$id))[i]) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 12, face = "bold")) +
      scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
      facet_wrap(~behavior, nrow = 3)
  )
}
par(ask=F)





## True proportions for Bayesian model (from modeled time segments)
bayes.props_weird<- bayes.res_weird %>% 
  rename(Encamped = X1, ARS = X2, Transit = X3) %>% 
  drop_na() %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",
               values_to = "prop") %>% 
  dplyr::select(id, tseg, behavior, prop, time1) %>% 
  mutate_at("behavior", ~recode(., 'Encamped' = 1, 'ARS' = 2, 'Transit' = 3))
bayes.props_weird$time1<- bayes.props_weird$time1 - 1


## Calculate RMSE
true.behavior_weird<- true.behavior.long_weird %>% bind_rows(.id = "id")

hmm.rmse_weird<- vector()
for (i in 1:length(unique(hmm.res_weird$id))) {
  ind<- unique(as.character(hmm.res_weird$id))[i]
  
  hmm.rmse_weird[i]<- sqrt(sum((hmm.props_weird[hmm.props_weird$id == ind, "prop"] - 
                           true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) / nrow(hmm.props_weird[hmm.props_weird$id == ind,]))
}


bayes.rmse_weird<- vector()
for (i in 1:length(unique(bayes.res_weird$id))) {
  ind<- unique(as.character(bayes.res_weird$id))[i]
  
  bayes.rmse_weird[i]<- sqrt(sum((bayes.props_weird[bayes.props_weird$id == ind, "prop"] - 
                             true.behavior_weird[true.behavior_weird$id == ind, "prop"])^2) /
                        nrow(bayes.props_weird[bayes.props_weird$id == ind,]))
}


rmse.df_weird<- data.frame(id = rep(unique(as.character(hmm.res_weird$id)), 2),
                    track_length = factor(rep(rep(c(1000,5000,10000,50000), each = 5), 2),
                                          levels = c("1000","5000","10000","50000")),
                    rmse = c(bayes.rmse_weird, hmm.rmse_weird),
                    method = rep(c("Bayesian","HMM"), each = 20))

#summarize results
rmse.df_weird %>% 
  group_by(method, track_length) %>% 
  summarise(mean=mean(rmse))


#Plot results

p.rmse_weird<- ggplot(rmse.df_weird, aes(track_length, rmse, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "RMSE\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  ylim(0, 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 10))




plot_grid(NULL, NULL, NULL,
          p.time, NULL, p.brk,
          # NULL, NULL, NULL,
          # p.coarse, NULL, p.rmse,
          NULL, NULL, NULL,
          p.coarse_weird, NULL, p.rmse_weird,
          align = "hv", nrow = 4, rel_widths = c(1,0.1,1), rel_heights = c(0.2,1,0.1,1))

# ggsave("Figure 3 (method comparison).png", width = 12, height = 12, units = "in", dpi = 330)