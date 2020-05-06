#########################
### Method Comparison ###
#########################

library(tidyverse)
library(wesanderson)
library(lubridate)


# Load elapsed time
seg.time<- read.csv("Bayesian_elapsed_time.csv")
lda.time<- read.csv("LDA_elapsed_time.csv")
bcpa.time<- read.csv("BCPA_elapsed_time.csv")
hmm.time<- read.csv("HMM_elapsed_time.csv")

# Load breakpoints
bayes.brkpts<- read.csv("Bayesian_allbreakpts.csv")
bcpa.brkpts<- read.csv("BCPA_allbrkpts.csv")

# Load results
bayes.res<- read.csv("Modeled MM Sim Tracks w Behav.csv")
bcpa.res<- read.csv("Clustered BCPA data.csv")
hmm.res<- read.csv("HMM results.csv")



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


ggplot(time, aes(track_length, time, fill = method)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Elapsed Time (min)") +
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) +
  theme_bw()


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
  mutate(track_length = case_when(str_detect(id, "_1") ~ "1k",
                                  str_detect(id, "_2") ~ "5k",
                                  str_detect(id, "_3") ~ "10k",
                                  str_detect(id, "_4") ~ "50k")) %>% 
  # summarise(mean = mean(freq), sd = sd(freq), n = sum(n)) %>%
  group_by(method, track_length, id) %>% 
  filter(acc == "Accurate" | acc == "Accurate Duplicate") %>% 
  summarise(freq = sum(freq)) %>%  #and calculate accuracy across all sims combined
  ungroup()

brkpt.acc$track_length<- brkpt.acc2$track_length %>% 
  factor(., levels = c('1k','5k','10k','50k'))


#Accuracy (includes 'accurate' and 'accurate duplicate' classifications)
ggplot(brkpt.acc, aes(track_length, freq, fill = method)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Proportion of Accurate Breakpoints") +
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) +
  theme_bw()


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
bcpa.res$method<- rep("BCPA/K-means", nrow(bcpa.res))
hmm.res$method<- rep("HMM", nrow(hmm.res))

bayes.res<- bayes.res %>% rename(state = behav)
bcpa.res<- bcpa.res %>% rename(state = cluster)
hmm.res<- hmm.res %>% rename(state = hmm.state, id = ID)

# Combine all datasets
res<- rbind(bayes.res[,c("id","behav_fine","behav_coarse","track_length","state","method")],
            bcpa.res[,c("id","behav_fine","behav_coarse","track_length","state","method")],
            hmm.res[,c("id","behav_fine","behav_coarse","track_length","state","method")])



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

summ.stats_coarse$track_length<- summ.stats$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))


ggplot(summ.stats_coarse, aes(track_length, acc, fill = method)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Accuracy of Behavior Estimates") +
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) +
  theme_bw()




#Fine-scale behavior
res %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_fine == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))

summ.stats_fine<- res %>% 
  drop_na() %>% 
  group_by(method, track_length, id) %>% 
  filter(., behav_fine == state) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  ungroup()

summ.stats_fine$track_length<- summ.stats$track_length %>% 
  factor(., levels = c('1000','5000','10000','50000'))


ggplot(summ.stats_fine, aes(track_length, acc, fill = method)) +
  geom_boxplot() +
  labs(x="Track Length", y = "Accuracy of Behavior Estimates") +
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) +
  ylim(0,1) +
  theme_bw()
