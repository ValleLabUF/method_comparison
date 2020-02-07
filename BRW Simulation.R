library(tidyverse)
library(circular)

source('Simulation Functions.R')

#######################################
#######################################
#### Biased Random Walk Simulation ####
#######################################
#######################################


#randomly sample coordinates for 50 ACs

#create vector of coordinates
set.seed(2)
AC.x<- sample(30, 10, replace = FALSE)
set.seed(1)
AC.y<- sample(30, 10, replace = FALSE)
AC<- data.frame(x=AC.x, y=AC.y)

ggplot(AC, aes(x, y)) +
  geom_point(size = 2) +
  theme_bw() +
  coord_equal()


#simulate track
#(n=250 per phase for 30 phases; 'a' and 'b' are shape and scale params for gamma distribution, Z.center is a matrix/DF of AC locations, Z0 is the initial location, and rho is the concentration param of a wrapped cauchy distribution that governs how attracted the simulation is to ACs)

set.seed(3)
track<- multiBRW.sim(n=1000, a = 1, b = 1, nphases = 15, Z.center = AC, Z0 = c(3,4), rho = 0.8)
track$time1<- 1:nrow(track)  #add variable for time


ggplot(data = track, aes(x, y)) +
  geom_path(size=0.5, color='gray75') +
  geom_point(size = 2, alpha = 0.3) +
  geom_point(data = AC, aes(x,y), color = "gold", size = 4) +
  geom_point(data = track[1,], aes(x,y), color = "green", pch = 1, size = 3, stroke = 1.5) +
  geom_point(data = track[nrow(track),], aes(x,y), color = "red", pch = 2, size = 3, stroke = 1.5) +
  theme_bw() +
  coord_equal()





#export
write.csv(track, "BRW_sim.csv", row.names = F)

