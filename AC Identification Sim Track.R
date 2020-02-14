###################################
#### Identify Activity Centers ####
###################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/activcenter_subset_locations")

set.seed(10)

#load libraries and read important functions
library('Rcpp')
library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(progress)


sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')
source('helper functions.R')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")


#load data
dat<- read.csv("BRW_tsegs.csv", as.is = TRUE)
obs<- get.summary.stats_obs(dat)  #frequency of visitation in each location (column) for each time segment (row)
obs1<- as.matrix(obs[,-1])  #for proper use by model


#geographical coordinates of locations
grid<- raster(extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)) + 5)
res(grid)<- 5
grid[]<- 0

grid.cell.locs<- coordinates(grid) %>% data.frame()
names(grid.cell.locs)<- c("x", "y")
grid.cell.locs$grid.cell<- 1:length(grid)
grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]


#Define initial activity centers (top 10 by # of obs)
tmp<- colSums(obs[,-1]) %>% data.frame(grid.cell = colnames(obs[,-1]), nobs = .) %>%
  arrange(desc(nobs)) %>% slice(n=1:10) %>% dplyr::select(grid.cell)
tmp<- tmp$grid.cell %>% as.character() %>% as.numeric()
ind<- sample(tmp, size = 10, replace = F)

ac.coord.init<- grid.coord[ind,]

#top 10
ac.coord.init2<- grid.coord[tmp,]

#potential locations for activity centers (AC)
possib.ac=grid.coord #these don't have to be identical (i.e., we can define AC's on a coarser grid)

### Run Gibbs sampler

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=10
gamma1=0.1

#run gibbs sampler
options(warn=2)

res=gibbs.activity.center(dat=obs1,grid.coord=grid.coord[,-3],n.ac=n.ac,
                          ac.coord.init=ac.coord.init[,-3],gamma1=gamma1,
                          possib.ac=possib.ac[,-3])

#plot output and look at frequency of AC visitation
plot(res$logl,type='l')
plot(res$phi,type='l')



##############################################
### Extract AC Coordinates and Assignments ###
##############################################

##use ACs from iteration with max log likelihood (after burn-in)
ML<- res$logl %>% order(decreasing = T)
ML<- ML[ML > 500][1]
ac<- res$z[ML,]
ac.coords<- matrix(NA, length(unique(ac)), 2)
colnames(ac.coords)<- c("x","y")
tmp<- res$coord[ML,]

for (i in 1:length(unique(ac))) {
  ac.coords[i,]<- round(c(tmp[i], tmp[i+length(unique(ac))]), 1)
}

ac.coords<- data.frame(ac.coords, ac = 1:length(unique(ac)))

#rearrange to match order/labelling of true ACs; order subject to change w each model run
ac.coords2<- data.frame(ac.coords[c(5,1,3,4,10,2,8,6,7,9),1:2], ac=1:length(unique(ac)))

table(ac)


############################
### Add ACs to Dataframe ###
############################

tseg.length<- dat %>% group_by(id, tseg) %>% tally()
tseg.length<- tseg.length$n
ac.aug<- rep(ac, times = tseg.length)

dat$model.ac<- ac.aug

#change ID of true.ac to match output from model for comparison
set.seed(2)
AC.x<- sample(30, 10, replace = FALSE)
set.seed(1)
AC.y<- sample(30, 10, replace = FALSE)
AC<- data.frame(x=AC.x, y=AC.y)


AC$true.ac<- 1:10
AC$model.ac<- c(5,1,3,4,10,2,8,6,7,9)  #order is subject to change w each run of model
for (i in 1:nrow(dat)) {
  dat$true.ac_mod[i]<- AC$model.ac[dat$true.ac[i]]
}

#Calculate number of obs per AC
dat %>% group_by(id) %>% dplyr::select(true.ac_mod) %>% table()
dat %>% group_by(id) %>% dplyr::select(model.ac) %>% table()


## Map
grid<- raster(extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)) + 5)
res(grid)<- 5
grid[]<- 0
borders<- rasterToPolygons(grid, dissolve = F)
borders_f<- fortify(borders)

ggplot() +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat, aes(x, y), color="grey45", size=2, alpha = 0.5) +
  geom_point(data = ac.coords, aes(x, y), color = "blue", size = 4, pch = 17, stroke = 1) +
  geom_point(data = AC, aes(x, y), color = "gold", size = 3, pch = 1, stroke = 1) +
  scale_color_viridis_c("Activity Center") +
  labs(x = "X", y = "Y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_equal()







####################################
#### Evaluate Accuracy of Model ####
####################################

#calc distances between true and modeled ACs
dist<- rep(NA, 10)
tmp<- (AC[,1:2] - ac.coords2[, 1:2])^2
dist<- sqrt(tmp[,1] + tmp[,2])

range(dist)  #min = 0.76, max = 2.86
mean(dist)  #mean = 1.84



#assignment of ACs to obs
tmp<- which(dat$true.ac_mod - dat$model.ac == 0) %>% length()
tmp/nrow(dat)  #98.7% accuracy for AC assignment
