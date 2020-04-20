library(tidyverse)

sim_hard<- read.csv("Sim1 hard.csv", as.is = T)
sim_mixed<- read.csv("Sim2 mixed.csv", as.is = T)
snki_res<- read.csv("snki movements with HMM states.csv", as.is = T)
dat_MM<- read.csv("Modeled MM Sim Tracks w Behav.csv", as.is = T)



########################
#### Simulated Data ####
########################

sim_mixed$behav_fine<- factor(sim_mixed$behav_fine, levels = c("Resting", "ARS", "Transit"))
sim_mixed$behav_coarse<- factor(sim_mixed$behav_coarse, levels = c("Resting", "ARS", "Transit"))
sim_mixed$hmm<- factor(sim_mixed$hmm, levels = 1:3)
levels(sim_mixed$hmm)<- c("Resting", "ARS", "Transit")
dat_MM$behav<- factor(dat_MM$behav, levels = c("Resting", "ARS", "Transit"))

sim_mixed<- sim_mixed[!is.na(sim_mixed$behav_fine),]



### Check accuracy of HMM model

#Fine-scale
true.b.fine<- sim_mixed$behav_fine %>%
  as.numeric()
hmm.b<- sim_mixed$hmm %>%
  as.numeric()

(which(true.b.fine == hmm.b) %>% length()) / length(true.b.fine)  #91.4% accurate


#Coarse-scale
true.b.coarse<- sim_mixed$behav_coarse %>%
  as.numeric()

(which(true.b.coarse == hmm.b) %>% length()) / length(true.b.coarse)  #87.5% accurate





### Accuracy of our Bayesian model

#Fine-scale
model.b<- dat_MM$behav[-1] %>%
  as.numeric()

(which(true.b.fine == model.b) %>% length()) / length(true.b.fine)  #82.5% accurate


#Coarse-scale

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)  #98.9% accurate







### Plot maps for results from each model

#HMM

ggplot() +
  geom_path(data = sim_mixed, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = sim_mixed, aes(x, y, fill=hmm), size=2.5, pch=21) +
  scale_fill_viridis_d("Behavior") +
  geom_point(data = sim_mixed[1,], aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = sim_mixed[nrow(sim_mixed),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_equal()




#Bayesian model

ggplot() +
  geom_path(data = dat_MM, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat_MM[-1,], aes(x, y, fill=behav), size=2.5, pch=21) +
  scale_fill_viridis_d("Behavior") +
  geom_point(data = dat_MM[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat_MM[nrow(dat_MM),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_equal()








#############################
#### Snail Kite Analysis ####
#############################

### Plot maps for results from each model

#HMM

ggplot() +
  geom_path(data = snki_res, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = snki_res, aes(x, y, fill=factor(state)), size=2.5,
             pch=21) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_equal()  +
  facet_wrap(~id)
