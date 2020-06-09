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
seg.time<- read.csv("Bayesian_elapsed_time.csv")
lda.time<- read.csv("LDA_elapsed_time.csv")
bcpa.time<- read.csv("BCPA_elapsed_time.csv")
hmm.time<- read.csv("HMM_elapsed_time.csv")
# hmm2.time<- read.csv("HMM2_elapsed_time.csv")

# Load breakpoints
bayes.brkpts<- read.csv("Bayesian_allbreakpts.csv")
bcpa.brkpts<- read.csv("BCPA_allbrkpts.csv")

# Load results
bayes.res<- read.csv("Modeled MM Sim Tracks w Behav.csv")
# bcpa.res<- read.csv("Clustered BCPA data.csv")
hmm.res<- read.csv("HMM results.csv")
# hmm2.res<- read.csv("HMM2 results.csv")

# Load true breakpoints
true.brkpts<- read.csv("CRW_MM_sim_brkpts.csv", as.is = T)

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
  labs(x="\nTrack Length (observations)", y = "Elapsed Time (min)\n") +
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
  labs(x="\nTrack Length (observations)", y = "Proportion of Accurate Breakpoints\n") +
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



##############################################
### Compare Accuracy of Behavior Estimates ###
##############################################

# Assign identifiers by method and make consistent behavior colname
bayes.res$method<- rep("Bayesian", nrow(bayes.res))
# bcpa.res$method<- rep("BCPA/K-means", nrow(bcpa.res))
hmm.res$method<- rep("HMM", nrow(hmm.res))
# hmm2.res$method<- rep("HMM2", nrow(hmm2.res))

bayes.res<- bayes.res %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = c("Encamped","ARS","Transit"))) %>%
  mutate_at("state", as.numeric)
# bcpa.res<- bcpa.res %>% rename(state = cluster)
hmm.res<- hmm.res %>% rename(state = hmm.state, id = ID)
# hmm2.res<- hmm2.res %>% rename(state = hmm.state, id = ID)


# Modify hmm.res to be same as bayes.res format
# calc proportions of behaviors by true time segment and then identify dominant behavior
hmm.res2<- hmm.res %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg, brkpts = true.brkpts) %>% 
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
            # bcpa.res[,c("id","behav_fine","behav_coarse","track_length","state","method")],
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
summ.stats_coarse<- summ.stats_coarse %>% 
  filter(method == "Bayesian" | method == "HMM")  #don't compare BCPA behavior

p.coarse<- ggplot(summ.stats_coarse, aes(track_length, acc, fill = method, color = method)) +
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




#Fine-scale behavior
# res %>% 
#   drop_na() %>% 
#   group_by(method, track_length, id) %>% 
#   filter(., behav_fine == state) %>% 
#   tally() %>% 
#   mutate(acc = n/track_length) %>% 
#   summarise(min=min(acc), max=max(acc), mean=mean(acc))
# 
# summ.stats_fine<- res %>% 
#   drop_na() %>% 
#   group_by(method, track_length, id) %>% 
#   filter(., behav_fine == state) %>% 
#   tally() %>% 
#   mutate(acc = n/track_length) %>% 
#   ungroup()
# 
# summ.stats_fine$track_length<- summ.stats_fine$track_length %>% 
#   factor(., levels = c('1000','5000','10000','50000'))
# summ.stats_fine<- summ.stats_fine %>% 
#   filter(method == "Bayesian" | method == "HMM" | method == "HMM2")  #don't compare BCPA behavior
# 
# 
# p.fine<- ggplot(summ.stats_fine, aes(track_length, acc, fill = method, color = method)) +
#   geom_boxplot() +
#   stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
#                position = position_dodge(0.75),
#                fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
#   ylim(0,1) +
#   labs(x="\nTrack Length (observations)", y = "Accuracy of Behavior Estimates\n") +
#   scale_fill_manual("", values = wes_palette("Zissou1")[c(1,3,5)]) +
#   scale_color_manual("", values = wes_palette("Zissou1")[c(1,3,5)]) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14),
#         legend.position = "n",
#         legend.text = element_text(size = 10))



plot_grid(NULL, NULL, NULL,
          p.time, NULL, p.brk,
          NULL, NULL, NULL,
          p.coarse, NULL, NULL,
          align = "hv", nrow = 4, rel_widths = c(1,0.1,1), rel_heights = c(0.2,1,0.1,1))

# ggsave("Figure 3 (method comparison).png", width = 12, height = 8, units = "in", dpi = 330)





#### Compare Accuracy of Bayesian and HMM Proportion Estimates ####

##True proportions by simulation ID
bayes.list<- df.to.list(bayes.res, "id")
true.behavior.long<- purrr::map(bayes.list, . %>% 
                                  drop_na("behav_fine") %>% 
                                  mutate(true.tseg = rep(1:(track_length[1]/100),each = 100)) %>%
                                  mutate_at("behav_fine", as.factor) %>% 
                                  group_by(true.tseg, behav_fine) %>% 
                                  count(behav_fine, .drop = FALSE) %>% 
                                  mutate(prop = n/100) %>% 
                                  rename(behavior = behav_fine) %>% 
                                  map_df(., rep, 100) %>%
                                  arrange(true.tseg) %>%
                                  mutate(time1 = rep(1:(nrow(.)/3), each = 3)))

## True proportions for HMMs (from time segments using true breakpoints)
hmm.props<- hmm.res %>% 
  df.to.list(., ind = "id") %>% 
  purrr::map(., ~mutate(., time1 = 1:nrow(.))) %>% 
  purrr::map(., assign.time.seg, brkpts = true.brkpts) %>% 
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
  drop_na() %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",
               values_to = "prop") %>% 
  dplyr::select(id, tseg, behavior, prop, time1) %>% 
  mutate_at("behavior", ~recode(., 'Encamped' = 1, 'ARS' = 2, 'Transit' = 3))
bayes.props$time1<- bayes.props$time1 - 1


## Calculate SSE
true.behavior<- true.behavior.long %>% bind_rows(.id = "id")

hmm.sse<- vector()
for (i in 1:length(unique(as.character(hmm.res$id)))) {
  ind<- unique(as.character(hmm.res$id))[i]
  
  hmm.sse[i]<- sum((hmm.props[hmm.props$id == ind, "prop"] - 
    true.behavior[true.behavior$id == ind, "prop"])^2)
}


bayes.sse<- vector()
for (i in 1:length(unique(as.character(bayes.res$id)))) {
  ind<- unique(as.character(bayes.res$id))[i]
  
  bayes.sse[i]<- sum((bayes.props[bayes.props$id == ind, "prop"] - 
                      true.behavior[true.behavior$id == ind, "prop"])^2)
}


sse.df<- data.frame(id = rep(unique(as.character(hmm.res$id)), 2),
                 track_length = factor(rep(rep(c(1000,5000,10000,50000), each = 5), 2),
                                       levels = c("1000","5000","10000","50000")),
                 sse = c(bayes.sse, hmm.sse),
                 method = rep(c("Bayesian","HMM"), each = 20))



#Plot results

ggplot(sse.df, aes(track_length, sse, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "Sum of Squared Errors\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "top",
        legend.text = element_text(size = 10))
