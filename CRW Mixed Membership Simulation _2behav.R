#####################################
#####################################
#### Mixed Membership Simulation ####
#####################################
#####################################

## 2 behaviors for HMM

library(tidyverse)
library(circular)

source('Simulation Functions.R')

### Simulate full track ###

#define behaviors and randomly sample 50 (for 50 time segments)
#weight probs so that behavior 1 (Encamped) occurs 80% and behavior 2 (Transit) occurs 20%

set.seed(4)

#create vector of dominant behaviors per time segment (50 segments)
behav<- sample(c(1,2), 50, replace = TRUE, prob = c(0.90, 0.10))
table(behav) #check freq

#randomly choose 3 segments of each behavior to be 'pure' instead of mixed
behav1.pure<- sample(which(behav==1), 3, replace = FALSE)
behav2.pure<- sample(which(behav==2), 3, replace = FALSE)

pure<- c(behav1.pure, behav2.pure) %>% sort()


#create vector of behaviors within each time segment (duration of 100 steps)
behav.full<- vector("list", length(behav))
for (i in 1:length(behav)) {
  if (i %in% pure) {
    behav.full[[i]]<- rep(behav[i], 100)
  }else if (behav[i] == 1) {
    behav.full[[i]]<- sample(c(1,2), 100, replace = TRUE, prob = c(0.8, 0.2))
  } else if (behav[i] == 2) {
    behav.full[[i]]<- sample(c(1,2), 100, replace = TRUE, prob = c(0.2, 0.8))
  } 
}
behav.full<- unlist(behav.full)

SL.params<- data.frame(shape=c(2, 10), scale = c(1, 1))
TA.params<- data.frame(mu=c(pi, 0), rho = c(0.8, 0.8))


#simulate track
#n=duration of each observation (behav.full), behav is a vector of behaviors, SL.params and TA.params are DFs of the necessary params from which to generate distributions for SL and TA from gamma and wrapped cauchy distribs, and Z0 is the initial location

track<- CRW.sim(n=1, behav = behav.full, SL.params = SL.params, TA.params = TA.params, Z0=c(0,0))
names(track)[5]<- "behav_fine"
track$behav_fine<- factor(track$behav)
levels(track$behav_fine)<- c("Encamped","Transit")
track<- track %>% mutate(behav_coarse = c(NA, rep(behav, each=100)) %>% factor())
levels(track$behav_coarse)<- c("Encamped","Transit")

true.brkpts<- which(diff(behav) != 0) * 100


# plot that puppy
ggplot(data = track[-1,], aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))


#compare distributions of SL and TA among fine-scale behaviors
ggplot(track, aes(SL)) +
  geom_density(aes(fill=behav_fine), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()

ggplot(track, aes(TA)) +
  geom_density(aes(fill=behav_fine), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()


#compare distributions of SL and TA among coarse-scale behaviors
ggplot(track, aes(SL)) +
  geom_density(aes(fill=behav_coarse), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()

ggplot(track, aes(TA)) +
  geom_density(aes(fill=behav_coarse), alpha = 0.6, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw()





# write.csv(track, "CRW_MM_sim_2behav.csv", row.names = F)


