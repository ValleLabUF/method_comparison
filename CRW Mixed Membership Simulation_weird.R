############################################
#### Simulations w/ Weird Distributions ####
############################################

library(tidyverse)
library(tictoc)

source('Simulation Functions.R')
source('helper functions.R')

### Simulate full track w/ different numbers of time segments ###
## Create tracks w/ 10, 50, 100, and 500 segments while keeping 100 obs per tseg
## Generate 5 different versions of each of these tracks (all else being equal)

#define behaviors and sample them from Categorical distribution
#weight probs so that behavior 1 (Resting) occurs 50%, behavior 2 (Area-restricted search) occurs 30%, and behavior 3 (Transit) occurs 20%

set.seed(2)


#simulate track
ntseg<- c(10, 50, 100, 500)
nstep<- 100
SL.params<- data.frame(lo = c(0, 0, 0),
                       hi = c(Inf, Inf, Inf),
                       mu=c(0.25, 2, exp(2)),
                       sig = c(1, 2, exp(3)))
TA.params<- data.frame(par1 = c(0.5, -pi, 0),
                       par2 = c(0.5, pi, 1))


tic()
track.sim<- CRW.sim2(nsim=5, ntseg = ntseg, nstep = nstep, SL.params = SL.params,
                    TA.params = TA.params, Z0=c(0,0))
toc()
#takes ~36 s to run


#extract tracks
tracks<- track.sim$tracks %>%
  modify_depth(2, ~modify_at(., "id", as.character)) %>%
  modify_depth(1, ~map_dfr(., `[`)) %>%
  map_dfr(`[`)

#extract breakpoints
max.length<- track.sim$brkpts %>%  #to set max number of columns for DF
  modify_depth(2, ~length(.)) %>%
  unlist() %>%
  max()

brkpts<- track.sim$brkpts %>%
  modify_depth(2, function(x) {c(x, rep(NA, max.length-length(x)))}) %>%
  modify_depth(1, ~map_dfr(., `[`)) %>%
  map(t) %>%
  map(as.data.frame) %>% 
  map_dfr(., `[`)
brkpts<- cbind(id = unique(tracks$id), brkpts)
names(brkpts)<- c('id', paste0("Brk_",1:(ncol(brkpts)-1)))


## Plot tracks ##

ggplot(data = tracks %>%
         group_by(id) %>%
         slice(2:n()) %>%
         ungroup(), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = tracks %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  # coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(track_length ~ id, scales = "free")




## Compare distributions of SL and TA among fine-scale behaviors ##

# behav_fine
ggplot(tracks[tracks$track_length == 1000,], aes(SL, color=behav_fine)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(SL, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(tracks[tracks$track_length == 1000,], aes(TA, color=behav_fine)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(TA, color=behav_fine),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


# behav_coarse
ggplot(tracks[tracks$track_length == 1000,], aes(SL, color=behav_coarse)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(SL, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()


ggplot(tracks[tracks$track_length == 1000,], aes(TA, color=behav_coarse)) +
  geom_line(stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data =  tracks[tracks$track_length == 5000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 10000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  geom_line(data = tracks[tracks$track_length == 50000,], aes(TA, color=behav_coarse),
            stat = "density", alpha = 0.5, na.rm = T) +
  scale_color_viridis_d("Behavior") +
  theme_bw()






# write.csv(tracks, "CRW_MM_sim_weird.csv", row.names = F)
# write.csv(brkpts, "CRW_MM_sim_brkpts_weird.csv", row.names = F)

