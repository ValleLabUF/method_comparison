################################
#### Run Segmentation Model ####
################################



### Data Prep for Segmentation Model ###

library(dplyr)
library(ggplot2)
library(lubridate)
library(sp)
library(raster)
library(rgdal)


dat<- read.csv("BRW_sim.csv", as.is = TRUE)

### Create Grid to Discretize Space

# 5 unit res w 2.5 unit buffer on each side
grid<- raster(extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)) + 5)
res(grid)<- 5
grid[]<- 0
dat$grid.cell<- cellFromXY(grid, dat[,c("x","y")])


### Plot all points over grid

#Create grid cell borders
borders<- rasterToPolygons(grid, dissolve = F)
borders_f<- fortify(borders)

#Calc points per cell
tab<- table(cellFromXY(grid, dat[,c("x","y")]))
grid[as.numeric(names(tab))] <- tab
grid_f<- as.data.frame(grid, xy = TRUE)
names(grid_f)[3]<- "count"


#plot points over grid
ggplot() +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat, aes(x=x, y=y), color = "firebrick", size=0.5, alpha=0.5) +
  labs(x = "X", y = "Y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_equal()

#plot density surface of points in grid
ggplot() +
  geom_tile(data=grid_f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "X", y = "Y") +
  theme_bw() +
  coord_equal()


#heatmap of grid cell occupancy over time
nloc<- ncell(grid)
nobs<- nrow(dat)
obs<- matrix(0, nobs, nloc)

for (i in 1:nrow(dat)) {
  obs[i, dat$grid.cell[i]]<- 1
}

obs<- data.frame(obs)
names(obs)<- 1:nloc
obs.long<- obs %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs.long$key<- as.numeric(obs.long$key)
obs.long$value<- factor(obs.long$value)
levels(obs.long$value)<- c("Absence","Presence")

#View grid cell use over time
ggplot(obs.long, aes(x=time, y=key, fill=value)) +
  geom_tile() +
  scale_fill_viridis_d("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Observations", y = "Grid Cell") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        title = element_text(size = 20))






setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_model")


library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)

source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/method_comparison")

## Prepare data

dat$id<- 1
dat.list<- df.to.list(dat = dat)

#only select necessary cols and re-number grid cell IDs
dat.long<- map_dfr(dat.list, `[`) %>% dplyr::select(id, grid.cell, time1)  #create DF
names(dat.long)[2]<- "loc.id"
dat.long$loc.id<- dat.long$loc.id %>% factor()
levels(dat.long$loc.id)<- 1:length(unique(dat.long$loc.id))  #change from raw to modified cell ID
dat.long$loc.id<- dat.long$loc.id %>% as.character() %>% as.numeric()

#convert back to list
dat.list2<- df.to.list(dat.long)



ngibbs = 10000

#prior
alpha=1

## Run Gibbs sampler
plan(multisession)
dat.res<- space_segment(data = dat.list2, ngibbs = ngibbs, alpha = alpha)
###Takes 3 min to run for 10000 iterations


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- unique(dat.long$id)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)


## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- matrix(c(1,dat.res$brkpts[[1]][[ML]]), 1, length(dat.res$brkpts[[1]][[ML]])+1) %>% data.frame()
names(brkpts)[1]<- "id"


## Heatmaps

#since trouble w/ line of code in plot.heatmap(), including here with modified code
data<- dat.list2[[1]]

#re-define loc.id based only on those visited by this individual
uni.loc=unique(data$loc.id)
aux=data.frame(loc.id=uni.loc,loc.id1=1:length(uni.loc))
dat1=merge(data,aux,all=T)
dat1$loc.id=dat1$loc.id1
dat2=dat1[order(dat1$time1),c('loc.id','time1')]

nloc<- length(uni.loc)
nobs<- nrow(data)
obs<- matrix(0, nobs, nloc)

for (i in 1:nrow(data)) {
  obs[i, dat2$loc.id[i]]<- 1
}

obs<- data.frame(obs)
names(obs)<- 1:nloc
obs.long<- obs %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
obs.long$key<- as.numeric(obs.long$key)
obs.long$value<- factor(obs.long$value)
levels(obs.long$value)<- c("Absence","Presence")

ind=which(unique(data$id) == brkpts$id)
breakpt<- brkpts[ind,-1] %>% t() %>% data.frame()
names(breakpt)<- "breaks"


ggplot(obs.long, aes(x=time, y=key, fill=value)) +
  geom_tile() +
  scale_fill_viridis_d("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7], size = 0.15) +
  labs(x = "Observations", y = "Grid Cell") +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        title = element_text(size = 20))




#Compare true and modeled breakpoints

true.brkpts<- which(diff(dat$true.ac) != 0)
model.brkpts<- t(breakpt) %>% as.vector()
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))

accuracy<- matrix(NA,length(model.brkpts),1)
for (i in 1:length(model.brkpts)) {
  
  tmp<- c(model.brkpts[i] - (10:0), model.brkpts[i] + (1:10)) %in% true.brkpts %>% sum()
  
  if (tmp == 0) {
    accuracy[i]<- "Inaccurate"
  } else {
    accuracy[i]<- "Accurate"
  }
}

if (sum(abs(diff(model.brkpts)) <= 10) > 0) {
  ind<- which(abs(diff(model.brkpts)) <= 10)
  ind<- sort(c(ind, ind+1))
}

ind.acc<- ind[which(accuracy[ind] == "Accurate")]
ind.inacc<- ind[which(accuracy[ind] == "Inaccurate")]
accuracy[ind.acc]<- "Accurate Duplicate"
accuracy[ind.inacc]<- "Inaccurate Duplicate"
accuracy<- c(rep("True",length(true.brkpts)), accuracy)

all.brkpts$acc<- accuracy

ggplot(all.brkpts, aes(x=brks, y=type, color = acc)) +
  geom_point(size=3) +
  theme_bw() +
  labs(x="Time", y="Type") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10), legend.position = "top") +
  scale_color_manual("Accuracy", values = c("forestgreen","lightgreen","firebrick","salmon","black"))



######################################
#### Assign Spatial Time Segments ####
######################################

breakpt<- brkpts[,-1] %>% as.numeric(.[1,])
breakpt1=c(0,breakpt,Inf)
n=length(breakpt1)
res=matrix(NA,nrow(dat.list[[1]]),1)
for (i in 2:n){
  ind=which(dat.list2[[1]]$time1>=breakpt1[i-1] & dat.list2[[1]]$time1<breakpt1[i])
  res[ind,]=i-1
}
dat_out<- dat.list[[1]]
dat_out$tseg<- as.vector(res)





#export
write.csv(dat_out, "BRW_tsegs.csv", row.names = F)
