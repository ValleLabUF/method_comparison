library(bcpa)
library(tidyverse)
library(cluster)
library(factoextra)

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_behavior")
source('helper functions.R')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")
source('helper functions.R')

set.seed(1)



#################
### Prep Data ###
#################

d<- read.csv("CRW_MM_sim_multinom.csv", as.is = T)
true.brkpts<- read.csv("CRW_MM_sim_brkpts.csv", as.is = T)

d.list<- df.to.list(d, ind = "id")
d.list<- map(d.list, ~mutate(., time = seq(c(ISOdate(2020,4,29)), by = "hour",
                                           length.out = n())))  #create fake hourly times 

mytrack<- map(d.list, ~MakeTrack(.$x, .$y, .$time))

#Plot tracks
par(mfrow=c(2,2), ask=T)
for (i in 1:length(mytrack)) {
  plot(mytrack[[i]], main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)



################
### Run BCPA ###
################

Simp1<- map(mytrack, GetVT)

Simp.ws<- list()
elapsed.time<- list()
cp<- list()
cw=30
#clusterwidth= the number of times a changepoint must be identified in the moving window for it to count. 1=all changepoints
#windowsize = # points to consider each time the model is run, k=2: default sensitivity (lower=less sensitive)

for (i in 1:length(Simp1)) {
  start.time<- Sys.time()
  
  Simp.ws[[i]]<- WindowSweep(Simp1[[i]], "V*cos(Theta)", windowsize=80, progress=T, K=2)
  cp[[i]]<-ChangePointSummary(Simp.ws[[i]], clusterwidth=cw)
  
  end.time<- Sys.time()
  elapsed.time[[i]]<- difftime(end.time, start.time, units = "min")
}

names(elapsed.time)<- names(d.list)

#Plot results
par(mfrow=c(2,2), ask=T)
for (i in 1:length(Simp.ws)) {
  plot(Simp.ws[[i]], type="flat", clusterwidth=cw, main = paste("ID",names(d.list)[i]))
}
par(mfrow=c(1,1), ask=F)


#Plot diagnostics
par(ask=T)
for (i in 1:length(Simp.ws)) {
  DiagPlot(Simp.ws[[i]])
}
par(ask=F)


# Store model results in DFs
y<- list()
for (i in 1:length(d.list)) {
  y[[i]]<- cbind(d.list[[i]][3:nrow(d.list[[i]]),],
                 t = Simp.ws[[i]]$t)  #number of lines in Simp.ws
  
  y[[i]]$group<-0
  len<- length(cp[[i]]$phases$t0)
  
  for (j in 1:len){ #at least the number of changepoints
    ind1 = cp[[i]]$phases$t0[j] 
    ind2 = cp[[i]]$phases$t1[j] 
    test1 = (y[[i]]$t>=ind1 & y[[i]]$t<=ind2)
    y[[i]]$group[test1]=j  #assign time segments to obs
  }
}



## Compare True vs Modeled Breakpoints
all.brkpts<- list()
for (i in 1:length(cp)) {
  all.brkpts[[i]]<- brkpt.accuracy(model.brkpts = round(cp[[i]]$breaks$middle, 0),
                                   true.brkpts = true.brkpts[i,-1],
                                   acc.tol = 20, dup.tol = 20, miss.tol = 50)
}


## Calculate percentage of each measure of breakpoint accuracy for summary
brkpt.acc<- map(all.brkpts, . %>%
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

## 1117 of 1976 brkpts missing (56.5% missing)


#Compare elapsed time
time<- map_dfr(elapsed.time, `[`) %>% 
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


#export breakpoints for easier reference and elapsed time for method comparison
names(all.brkpts)<- names(d.list)
all.brkpts<- bind_rows(all.brkpts, .id = 'id')
# write.csv(all.brkpts, "BCPA_allbrkpts.csv", row.names = F)
# write.csv(time, "BCPA_elapsed_time.csv", row.names = F)  #units = min





##########################
### K-means Clustering ###
##########################

# Prep data
d2<- map(cp, . %>% 
           pluck("phases") %>% 
           dplyr::select(mu.hat, s.hat, rho.hat) %>% 
           mutate(group = 1:nrow(.)))

# K-Means Cluster Analysis
par(ask=T)
for (i in seq_along(d2)) {
  print(
    fviz_nbclust(d2[[i]], kmeans, k.max = 5, method = "wss")
  )
}
par(ask=F)

# Identify K per elbow method; only 3 (from shortest tracks) found 3 states to be optimal
k.optim<- c(3, 4, 3, 3, 4,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2,
            2, 2, 2, 2, 2)

# Run k-means using values of k
kmeans.fit<- list()
tictoc::tic()
kmeans.fit<- map(d2, ~kmeans(.x[,1:3], 3))
tictoc::toc()

# Recode behavior numbers to match with true behavior numbers (using cluster means)
# print(kmeans.fit)
behav.order<- list(c(3,1,2), c(3,1,3), c(3,1,2), c(1,2,3), c(1,2,3),
                   c(2,1,3), c(1,2,3), c(1,2,3), c(2,1,3), c(3,1,3),
                   c(3,3,1), c(1,3,2), c(3,1,2), c(3,1,3), c(2,1,3),
                   c(2,1,3), c(1,2,3), c(3,2,1), c(2,1,3), c(3,2,1))


# append cluster assignment
bcpa.clust<- pmap(list(d2, kmeans.fit), ~data.frame(.x, cluster = .y$cluster))
for (i in seq_along(bcpa.clust)) {
  bcpa.clust[[i]]<- bcpa.clust[[i]] %>% 
    mutate_at("cluster", as.factor) %>% 
    mutate_at("cluster", ~recode(., '1' = behav.order[[i]][1],
                                 '2' = behav.order[[i]][2], '3' = behav.order[[i]][3]))
}
bcpa.dat<- pmap(list(y, bcpa.clust), ~merge(.x, .y[,c("group","cluster")], by = "group")) %>% 
  bind_rows()



########################
### Model Validation ###
########################

# Overall
bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == cluster) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 45.9% to 78.9%

# For 'Resting' behavior
rest.size<- bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 1) %>% 
  tally()

bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 1 & cluster == 1) %>% 
  tally() %>% 
  left_join(., rest.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 23.7% to 100.0%

# For 'ARS' behavior
ars.size<- bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 2) %>% 
  tally()

bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 2 & cluster == 2) %>% 
  tally() %>% 
  left_join(., ars.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 31.5% to 95.5%

# For 'Transit' behavior
transit.size<- bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 3) %>% 
  tally()

bcpa.dat %>% 
  group_by(track_length, id) %>% 
  filter(., behav_coarse == 3 & cluster == 3) %>% 
  tally() %>% 
  left_join(., transit.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 35.2% to 100.0%




# Export data
# write.csv(bcpa.dat, "Clustered BCPA data.csv", row.names = F)