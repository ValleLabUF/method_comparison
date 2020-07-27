##################
#### Figure 1 ####
##################


#Running this script is dependent on having run all sets of analyses from generation of simulation to clustering via LDA and plotting results


dat2_focal<- dat2 %>% filter(id == '2_2')

## Sim Track Map
ggplot(data = dat2_focal %>% slice(2:n()), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(fill = "grey45", pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2_focal %>% slice(n=1), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat2_focal %>% slice(n()), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_void() +
  theme(axis.title = element_blank(),
        legend.position = c(0.9,0.2)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_cartesian()


ggsave("Figure 1a.png", width = 6, height = 4, units = "in", dpi = 330)




## Discretized Movmt Param Distribs
behav.df<- map_dfr(behav.list, `[`)

ggplot(behav.df %>% filter(id == '2_2'), aes(x=dist)) +
  geom_density(fill = "grey55") +
  geom_segment(data = data.frame(dist.bin.lims), aes(x=dist.bin.lims, xend=dist.bin.lims,
                                                     y=0, yend=0.4),
               linetype = "dashed", lwd = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Step Length", y = "Density") +
  annotate("segment", x = mean(dist.bin.lims[1:2]), xend = -2, y = 0.4, yend = 0.42,
           colour = "black") +
  annotate("segment", x = mean(dist.bin.lims[2:3]), xend = mean(dist.bin.lims[2:3]),
           y = 0.4, yend = 0.42, colour = "black") +
  annotate("text",
           x = c(-3, mean(dist.bin.lims[2:3]), mean(dist.bin.lims[3:4]),
                                  mean(dist.bin.lims[4:5]), mean(dist.bin.lims[5:6])),
           y = c(0.44, 0.44, 0.4, 0.4, 0.4),
           label = paste("Bin", 1:5),
           size = 7)

ggsave("Figure 1b.png", width = 9, height = 6, units = "in", dpi = 330)


ggplot(behav.df %>% filter(id == '2_2'), aes(x=rel.angle)) +
  geom_density(fill = "grey55") +
  geom_segment(data = data.frame(angle.bin.lims), aes(x=angle.bin.lims, xend=angle.bin.lims,
                                                     y=0, yend=0.25),
               linetype = "dashed", lwd = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Turning Angle (rad)", y = "Density") +
  annotate("text",
           x = (pi/8) + angle.bin.lims[1:8],
           y = rep(0.25, 8),
           label = paste("Bin", 1:8),
           size = 7)

ggsave("Figure 1c.png", width = 9, height = 6, units = "in", dpi = 330)






## Segmentation Heatmap
# plot.heatmap(data = list(behav.list$`2_2`[1:5000,]), nbins = c(5,8), brkpts = brkpts[7,],
             # dat.res = dat.res, type = "behav", title = F, legend = T)

nbins = c(5,8)
data<- behav.list$`2_2`
behav.heat<- behav.seg.image(data, nbins)

SL<- data.frame(behav.heat$SL)
names(SL)<- 1:nbins[1]
SL<- SL %>% gather(key, value) %>% mutate(time=rep(data$time1,
                                                   times=nbins[1]),
                                          behav=rep("SL",nrow(data)*nbins[1]))

TA<- data.frame(behav.heat$TA)
names(TA)<- 1:nbins[2]
TA<- TA %>% gather(key, value) %>% mutate(time=rep(data$time1,
                                                   times=nbins[2]),
                                          behav=rep("TA",nrow(data)*nbins[2]))

behav.heat_long<- rbind(SL,TA)
behav.heat_long$value<- factor(behav.heat_long$value)
levels(behav.heat_long$value)<- c("Unoccupied","Occupied")

ind=which(unique(data$id) == brkpts$id)
breakpt<- brkpts[ind,-1] %>% purrr::discard(is.na) %>% t() %>% data.frame()
names(breakpt)<- "breaks"


behav.heat_long$behav<- factor(behav.heat_long$behav)
levels(behav.heat_long$behav)<- c("Step Length", "Turning Angle")

ggplot(behav.heat_long, aes(x=time, y=key, fill=value)) +
    geom_tile() +
    facet_wrap(~behav, scales = 'free', nrow = 2) +
    scale_fill_viridis_d('') +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = breakpt, aes(xintercept = breaks),
               color = viridis(n=9)[7], size = 1) +
    labs(x = "Time", y = "Bin") +
    theme_bw() +
    theme(axis.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 16, face = 'bold'),
          plot.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 16))

ggsave("Figure 1d.png", width = 7, height = 5, units = "in", dpi = 330)
  
  
  
  
## Behavior Histograms
behav.res_focal<- behav.res[[7]]  #for simulation ID 2_2
behav.res_focal$param<- factor(behav.res_focal$param)
levels(behav.res_focal$param)<- c("Step Length", "Turning Angle")
behav.res_focal<- behav.res_focal %>% filter(behav <=3)
behav.res_focal$behav<- factor(behav.res_focal$behav, levels = c(3,1,2))
levels(behav.res_focal$behav)<- c("Resting","ARS","Transit")

ggplot(behav.res_focal, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "Bin", y = "Proportion") +
  scale_x_continuous(breaks = 1:8) +
  theme_bw() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.bottom = element_text(size = 16),
        strip.text = element_text(size = 16),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(behav ~ param, scales = "free_x")

ggsave("Figure 1e.png", width = 6, height = 4, units = "in", dpi = 330)






## Time Series Behavior Proportions
theta.estim.long_focal<- theta.estim.long[[7]]  #for simulation ID 2_2
theta.estim.long_focal$behavior<- factor(theta.estim.long_focal$behavior, levels = c(3,1,2))
levels(theta.estim.long_focal$behavior)<- c("Resting","ARS","Transit")


ggplot(theta.estim.long_focal) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "Time", y = "Proportion of Behavior") +
  scale_fill_viridis_d("") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 16),
        plot.margin = margin(5,20,5,5))

ggsave("Figure 1f.png", width = 6, height = 4, units = "in", dpi = 330)




## Sim Track Annotated Map
ggplot(data = dat2_focal %>% slice(2:n()), aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2_focal %>% slice(n=1), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat2_focal %>% slice(n()), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  scale_fill_viridis_d("", labels = c("Resting","ARS","Transit")) +
  theme_void() +
  theme(axis.title = element_blank(),
        legend.position = c(0.9,0.2)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  coord_cartesian()

ggsave("Figure 1g.png", width = 6, height = 4, units = "in", dpi = 330)








##################
#### Figure 2 ####
##################

library(sf)
library(lwgeom)
library(tidyverse)
library(nhdplusTools)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(raster)


## Load hydro data

setwd("~/Documents/NHD Seamless Geodatabase/NHDPlusNationalData")

download_dir<- download_nhdplushr(nhd_dir = getwd(), 
                                  hu_list = c("0307","0308","0309","0310","0311","0312","0313",
                                              "0314"), download_files = TRUE)
waterbody<- read_sf(file.path(download_dir, "nhd_hr_FL.gpkg"), "NHDWaterbody") %>%
  st_cast("MULTIPOLYGON")
waterbody2<- waterbody %>%
  filter(AreaSqKM > 50) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 250)

#import FL layer
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida") %>% st_transform(CRS("+init=epsg:32617"))

#mask flowline and waterbody to FL layer
waterbody_fl<- st_intersection(st_make_valid(waterbody2), fl) %>% 
  filter(GNIS_Name != "The Everglades")


## Load snail kite data

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/git_segmentation_model")

dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")
dat$date<- as_datetime(dat$date)




ggplot() +
  geom_sf(data = fl, fill = "grey45", color = NA) +
  geom_sf(data = waterbody_fl, fill = "white", color = NA) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  scale_color_viridis_d("", guide=F) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10))

setwd("~/Documents/Manuscripts/Bayesian Behavior Estimation Model/Figures")
# ggsave("Figure 2.png", width = 4, height = 6, units = "in", dpi = 330)





##################
#### Figure 4 ####
##################


## Part a: segmentation heatmap
data<- behav.list2$`2_2`  #ID 2_2
nbins<- c(5,8)
behav.heat<- behav.seg.image(data, nbins)

SL<- data.frame(behav.heat$SL)
names(SL)<- 1:nbins[1]
SL<- SL %>% gather(key, value) %>% mutate(time=rep(behav.list[["2_2"]]$time1, times=nbins[1]),
                                          behav=rep("Step Length", nrow(data)*nbins[1]))

TA<- data.frame(behav.heat$TA)
names(TA)<- 1:nbins[2]
TA<- TA %>% gather(key, value) %>% mutate(time=rep(behav.list[["2_2"]]$time1, times=nbins[2]),
                                          behav=rep("Turning Angle", nrow(data)*nbins[2]))

behav.heat_long<- rbind(SL,TA)
behav.heat_long$value<- factor(behav.heat_long$value)
levels(behav.heat_long$value)<- c("Unoccupied","Occupied")

breakpt<- bayes.brkpts %>% 
  filter(id == unique(data$id) & type == "Model") %>% 
  dplyr::select(brks) %>% 
  rename(breaks = brks)

true.brks<- data.frame(t(true.brkpts[7,-1])) %>%
  drop_na() %>% 
  rename(breaks = X7) %>% 
  pull(breaks)

ticks.bottom<-data.frame(x=true.brks, y=0.5, xend=true.brks, yend=1.5)
ticks.top<- data.frame(x=rep(true.brks, 2), y=NA, xend=rep(true.brks, 2), yend=NA)
# ticks.top$y<- rep(c(4.5,7.5), each = length(true.brks))
# ticks.top$yend<- rep(c(5.5,8.5), each = length(true.brks))
ticks.top$y<- rep(7.5, length(true.brks))
ticks.top$yend<- rep(8.5, length(true.brks))
ticks.top$behav<- rep(c("Step Length","Turning Angle"), each = length(true.brks))


p.seg<- ggplot(behav.heat_long, aes(x=time#, y=key
                                    )) +
  # geom_tile(fill="n") +
  facet_wrap(~behav, scales = 'free', nrow = 2) +
  scale_fill_viridis_d('') +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_vline(data = breakpt, aes(xintercept = breaks - 0.5), color = viridis(n=9)[7],
             size = 1, alpha = 1) +
  geom_segment(data = ticks.top, aes(x = x, xend = xend, y = y, yend = yend#, group = behav
                                     ),
               color = "black", size = 1, alpha = 1) +
  geom_segment(data = ticks.bottom, aes(x = x, xend = xend, y = y, yend = yend), color = "black",
               size = 1, alpha = 1) +
  labs(x = "\nTime", y = "\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 12, face = 'bold'),
        plot.title = element_text(size = 20, hjust = 0, vjust = -6),
        plot.margin = margin(0, 1, 0.5, 0.5, "cm"),
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.position = "top",
        legend.text = element_text(
          margin = margin(r = 15, unit = "pt")))




## Part b: behavior histogram
library(circular)

#calculate true proportions of SL and TA by behavior for all bins for ID 2_2
SL.params<- data.frame(par1 = c(0.25, 2, 10), par2 = c(1, 1, 1))
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))
true.b<- extract.behav.props(params = list(SL.params, TA.params),
                             lims = list(dist.bin.lims,angle.bin.lims),
                             behav.names = c("Encamped","ARS","Transit"))

bayes.b<- behav.res[[7]] %>% rename(., var = param)
bayes.b$behav<- bayes.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Encamped") %>%
  str_replace_all(., "3", "Transit") %>%
  factor(., levels = c("Encamped","ARS","Transit"))
  


p.hist<- ggplot(true.b, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = bayes.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
             color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")




##Part c: generate time series plots comparing behavior proportions
dat1<- theta.estim.long[[7]]
dat1$behavior<- dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

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

true.dat1<- true.behavior.long[[7]]
true.dat1$behavior<- true.dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

p.prop<- ggplot() +
  geom_line(data = dat1,
            aes(x=date, y=prop, color = behavior),
            size = 1.5) +
  scale_color_manual(values = c(viridis(n=20)[c(1,9)], "gold3"), guide=F) +
  new_scale_color() +
  geom_line(data = true.dat1,
            aes(x=time1, y=prop, color = behavior),
            size = 0.75) +
  scale_color_manual(values = c(viridis(n=20)[c(7,13)], "gold2"), guide=F) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) +
  facet_wrap(~behavior, nrow = 3)




## Make composite
library(gridExtra)
setwd("~/Documents/Manuscripts/Bayesian Behavior Estimation Model/Figures")

# png("Figure 4 (results from sim).png", width = 14, height = 10, units = "in", res = 330)
grid.arrange(p.seg, p.hist, p.prop, heights = c(0.2, 1, 0.1, 1),
             layout_matrix = rbind(c(NA, NA),
                                   c(1, 2),
                                   c(NA, NA),
                                   c(3, 3)))
# dev.off()




##################
#### Figure X ####
##################
library(circular)
library(wesanderson)
library(cowplot)
source('helper functions.R')

### Define bin limits
dat<- read.csv("CRW_MM_sim_multinom.csv", as.is = T)
dat$dt<- 3600
names(dat)[4:5]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat[dat$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(dat[dat$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins


### Gamma/Wrapped Cauchy

#Import distribution params for SL and TA from all HMMs
hmm.SL.params<- read.csv("HMM result gamma params.csv", as.is = T)
hmm.TA.params<- read.csv("HMM result wrapped Cauchy params.csv", as.is = T)

#Manipulate param dfs to reformat for extract.behav.props()
hmm.SL.params2<- list()
for (i in 1:nrow(hmm.SL.params)) {
  hmm.SL.params2[[i]]<- data.frame(par1 = as.numeric(hmm.SL.params[i,2:4]),
                                   par2 = as.numeric(hmm.SL.params[i,5:7]))
}
names(hmm.SL.params2)<- hmm.SL.params$X

hmm.TA.params2<- list()
for (i in 1:nrow(hmm.TA.params)) {
  hmm.TA.params2[[i]]<- data.frame(par1 = as.numeric(hmm.TA.params[i,2:4]),
                                   par2 = as.numeric(hmm.TA.params[i,5:7]))
}
names(hmm.TA.params2)<- hmm.TA.params$X



#convert mean and sd to shape and rate params for gamma dist
for (j in 1:length(hmm.SL.params2)) {
  for (i in 1:nrow(hmm.SL.params2[[j]])) {
    shape<- (hmm.SL.params2[[j]][i,1]^2) / (hmm.SL.params2[[j]][i,2]^2)
    rate<- hmm.SL.params2[[j]][i,1] / (hmm.SL.params2[[j]][i,2]^2)
    
    params<- c(shape, rate)
    
    hmm.SL.params2[[j]][i,]<- params
  }
}


SL.params<- data.frame(par1 = c(0.25, 2, 10), par2 = c(1, 1, 1))
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))

true.b<- extract.behav.props(params = list(SL.params, TA.params),
                             lims = list(dist.bin.lims, angle.bin.lims),
                             behav.names = c("Encamped","ARS","Transit"))


hmm.b<- map2(hmm.SL.params2, hmm.TA.params2, ~extract.behav.props(params = list(.x, .y),
                            lims = list(dist.bin.lims, angle.bin.lims),
                            behav.names = c("Encamped","ARS","Transit"))
)



### Bayesian
behav.res<-  read.csv("CRW MM LDA Phi values.csv", as.is = T)
behav.res<- df.to.list(behav.res, "id")

behav.order<- read.csv("CRW MM LDA behavior order.csv", as.is = T)
behav.order<- df.to.list(behav.order, "id")
behav.order<- map(behav.order, as.numeric)


#calculate true proportions of SL and TA by behavior for all bins for ID 2_2 for Bayesian model
bayes.b<- behav.res[[20]] %>% rename(., var = param)
bayes.b$behav<- bayes.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Transit") %>% 
  str_replace_all(., "2", "ARS") %>%
  str_replace_all(., "3", "Encamped") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.gamWC<- ggplot(true.b, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b[[20]], aes(x=bin,y=prop,group = behav), pch=21, size = 2,
             fill="grey45", color="black", stroke=1) +
  geom_point(data = bayes.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2,
             fill="grey85", color="black", stroke=1) +
  labs(x = "", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")



### Calculate root mean square error (RMSE)

#Bayesian
bayes.rmse<- list()
for (i in 1:length(behav.res)) {
  bayes.b<- behav.res[[i]] %>% 
    rename(., var = param)
  
  behavs<- c("Encamped","ARS","Transit")
  
  bayes.b$behav<- bayes.b$behav %>% 
    factor(levels = behav.order[[i]]) %>% 
    str_replace_all(., "1", behavs[behav.order[[i]][1]]) %>% 
    str_replace_all(., "2", behavs[behav.order[[i]][2]]) %>%
    str_replace_all(., "3", behavs[behav.order[[i]][3]]) %>%
    factor(., levels = c("Encamped","ARS","Transit"))
  bayes.b<- bayes.b %>% 
    group_by(var) %>% 
    arrange(behav) %>% 
    ungroup()
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(bayes.b$var))) {
    tmp<- bayes.b %>% 
      filter(var == unique(bayes.b$var)[j])
    true.tmp<- true.b %>% 
      filter(var == unique(bayes.b$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  bayes.rmse[[i]]<- data.frame(vec)
  names(bayes.rmse)[i]<- names(behav.res)[i]
}
bayes.rmse<- bind_rows(bayes.rmse) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




#HMM
hmm.rmse<- list()
for (i in 1:length(behav.res)) {
  hmm.b<- extract.behav.props(params = list(hmm.SL.params2[[i]], hmm.TA.params2[[i]]),
                              lims = list(dist.bin.lims, angle.bin.lims),
                              behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(hmm.b$var))) {
    tmp<- hmm.b %>% 
      filter(var == unique(hmm.b$var)[j])
    true.tmp<- true.b %>% 
      filter(var == unique(hmm.b$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  hmm.rmse[[i]]<- data.frame(vec)
  names(hmm.rmse)[i]<- names(behav.res)[i]
}
hmm.rmse<- bind_rows(hmm.rmse) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))



rmse.df<- data.frame(id = rep(rep(names(behav.res), 2), each = 2),
                    track_length = factor(rep(rep(rep(c(1000,5000,10000,50000),
                                                      each = 5), 2), each = 2),
                                          levels = c("1000","5000","10000","50000")),
                    rmse = rbind(bayes.rmse, hmm.rmse),
                    method = rep(c("Bayesian","HMM"), each = 40))


p.rmse<- ggplot(rmse.df, aes(track_length, rmse.value, fill = method, color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="", y = "RMSE\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "top",
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_wrap(~rmse.var, ncol = 2)





### Truncated_Log Normal/Wrapped Cauchy

### Define bin limits
dat_weird<- read.csv("CRW_MM_sim_weird.csv", as.is = T)
dat_weird$dt<- 3600
names(dat_weird)[4:5]<- c("dist","rel.angle")

angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(dat_weird[dat_weird$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(dat_weird[dat_weird$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins

hmm.SL.params_weird<- read.csv("HMM result step params_weird.csv", as.is = T)  #weird
hmm.TA.params_weird<- read.csv("HMM result angle params_weird.csv", as.is = T)  #weird

#Manipulate param dfs to reformat for extract.behav.props()
hmm.SL.params2_weird<- list()
for (i in 1:nrow(hmm.SL.params_weird)) {
  hmm.SL.params2_weird[[i]]<- data.frame(par1 = as.numeric(hmm.SL.params_weird[i,2:4]),
                                   par2 = as.numeric(hmm.SL.params_weird[i,5:7]))
}
names(hmm.SL.params2_weird)<- hmm.SL.params_weird$X

hmm.TA.params2_weird<- list()
for (i in 1:nrow(hmm.TA.params_weird)) {
  hmm.TA.params2_weird[[i]]<- data.frame(par1 = as.numeric(hmm.TA.params_weird[i,2:4]),
                                   par2 = as.numeric(hmm.TA.params_weird[i,5:7]))
}
names(hmm.TA.params2_weird)<- hmm.TA.params_weird$X



#convert mean and sd to shape and rate params for gamma dist
for (j in 1:length(hmm.SL.params2_weird)) {
  for (i in 1:nrow(hmm.SL.params2_weird[[j]])) {
    shape<- (hmm.SL.params2_weird[[j]][i,1]^2) / (hmm.SL.params2_weird[[j]][i,2]^2)
    rate<- hmm.SL.params2_weird[[j]][i,1] / (hmm.SL.params2_weird[[j]][i,2]^2)
    
    params<- c(shape, rate)
    
    hmm.SL.params2_weird[[j]][i,]<- params
  }
}


SL.params_weird<- data.frame(par1 = c(0.25, 2, exp(2)), par2 = c(1, 2, exp(3)))  #for weird distrib
TA.params<- data.frame(par1 = c(pi, pi, 0), par2 = c(0.8, 0, 0.8))

true.b_weird<- extract.behav.props_weird(params = list(SL.params_weird, TA.params),  #for weird distrib
                             lims = list(dist.bin.lims, angle.bin.lims),
                             behav.names = c("Encamped","ARS","Transit"))


hmm.b_weird<- map2(hmm.SL.params2_weird, hmm.TA.params2_weird,
                   ~extract.behav.props(params = list(.x, .y),
                                        lims = list(dist.bin.lims, angle.bin.lims),
                                        behav.names = c("Encamped","ARS","Transit"))
)



## Bayesian
behav.res_weird<-  read.csv("CRW MM LDA Phi values_weird.csv", as.is = T)
behav.res_weird<- df.to.list(behav.res_weird, "id")

behav.order_weird<- read.csv("CRW MM LDA behavior order_weird.csv", as.is = T)
behav.order_weird<- df.to.list(behav.order_weird, "id")
behav.order_weird<- map(behav.order_weird, as.numeric)


#calculate true proportions of SL and TA by behavior for all bins for ID 2_2 for Bayesian model
bayes.b_weird<- behav.res_weird[[20]] %>% rename(., var = param)
bayes.b_weird$behav<- bayes.b_weird$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Transit") %>% 
  str_replace_all(., "2", "ARS") %>%
  str_replace_all(., "3", "Encamped") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.normWC<- ggplot(true.b_weird, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = hmm.b_weird[[20]], aes(x=bin,y=prop,group=behav), pch=21, size=2,
             fill="grey45", color="black", stroke=1) +
  geom_point(data = bayes.b_weird, aes(x=bin, y=prop, group = behav), pch=21, size = 2,
             fill="grey85", color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")



### Calculate root mean square error (RMSE)

#Bayesian
bayes.rmse_weird<- list()
for (i in 1:length(behav.res_weird)) {
  bayes.b_weird<- behav.res_weird[[i]] %>% 
    rename(., var = param)
  
  behavs<- c("Encamped","ARS","Transit")
  
  bayes.b_weird$behav<- bayes.b_weird$behav %>% 
    factor(levels = behav.order_weird[[i]]) %>% 
    str_replace_all(., "1", behavs[behav.order_weird[[i]][1]]) %>% 
    str_replace_all(., "2", behavs[behav.order_weird[[i]][2]]) %>%
    str_replace_all(., "3", behavs[behav.order_weird[[i]][3]]) %>%
    factor(., levels = c("Encamped","ARS","Transit"))
  bayes.b_weird<- bayes.b_weird %>% 
    group_by(var) %>% 
    arrange(behav) %>% 
    ungroup()
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(bayes.b_weird$var))) {
    tmp<- bayes.b_weird %>% 
      filter(var == unique(bayes.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(bayes.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  bayes.rmse_weird[[i]]<- data.frame(vec)
  names(bayes.rmse_weird)[i]<- names(behav.res)[i]
}
bayes.rmse_weird<- bind_rows(bayes.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))




#HMM
hmm.rmse_weird<- list()
for (i in 1:length(behav.res_weird)) {
  hmm.b_weird<- extract.behav.props(params = list(hmm.SL.params2_weird[[i]],
                                                  hmm.TA.params2_weird[[i]]),
                              lims = list(dist.bin.lims, angle.bin.lims),
                              behav.names = c("Encamped","ARS","Transit"))
  
  #create vector to store RMSE for SL and TA separately
  vec<- vector()
  
  for (j in 1:length(unique(hmm.b_weird$var))) {
    tmp<- hmm.b_weird %>% 
      filter(var == unique(hmm.b_weird$var)[j])
    true.tmp<- true.b_weird %>% 
      filter(var == unique(hmm.b_weird$var)[j])
    
    vec[j]<- sqrt(sum((tmp$prop - true.tmp$prop)^2) / nrow(tmp))  #RMSE
  }
  
  hmm.rmse_weird[[i]]<- data.frame(vec)
  names(hmm.rmse_weird)[i]<- names(behav.res)[i]
}
hmm.rmse_weird<- bind_rows(hmm.rmse_weird) %>% 
  rename(value = vec) %>% 
  mutate(var = rep(c("Step Length", "Turning Angle"), times = 20))



rmse.df_weird<- data.frame(id = rep(rep(names(behav.res), 2), each = 2),
                     track_length = factor(rep(rep(rep(c(1000,5000,10000,50000),
                                                       each = 5), 2), each = 2),
                                           levels = c("1000","5000","10000","50000")),
                     rmse = rbind(bayes.rmse_weird, hmm.rmse_weird),
                     method = rep(c("Bayesian","HMM"), each = 40))


p.rmse_weird<- ggplot(rmse.df_weird, aes(track_length, rmse.value, fill = method,
                                         color = method)) +
  geom_boxplot() +
  stat_summary(geom = "crossbar", width = 0.6, fatten=1.5, color="black",
               position = position_dodge(0.75),
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  labs(x="\nTrack Length (observations)", y = "RMSE\n") +
  scale_fill_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  scale_color_manual("", values = wes_palette("Zissou1")[c(1,5)]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "n",
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")) +
  facet_wrap(~rmse.var, ncol = 2)




## Make composite
library(gridExtra)
setwd("~/Documents/Manuscripts/Bayesian Behavior Estimation Model/Figures")

png("Figure S1 (rmse from sim).png", width = 14.7, height = 10.5, units = "in", res = 330)
grid.arrange(p.gamWC, p.rmse, p.normWC, p.rmse_weird, heights = c(0.2, 1, 0.1, 1),
             widths = c(1, 0.2, 0.2, 0.2, 1, 0.1),
             layout_matrix = rbind(c(NA, NA, NA, NA, NA, NA),
                                   c(1, NA, 2, 2, 2, NA),
                                   c(NA, NA, NA, NA, NA, NA),
                                   c(3, NA, 4, 4, 4, NA)))
dev.off()






### Make plot comparing HMM gamma distribs against true generating trunc/log norm distrib
hmm.SL.df_weird<- hmm.SL.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.SL.params_weird<- data.frame(id = 1:3,
                            SL.params_weird,
                            behavior = c("Encamped","ARS","Transit"),
                            track_length = "True")


plot_data_SL<- 
  pmap_df(hmm.SL.df_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dgamma(x, shape = par1, rate = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_SL$behavior<- factor(plot_data_SL$behavior, levels = c("Encamped","ARS","Transit"))


true_plot_data_SL<- 
  pmap_df(true.SL.params_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(0, 40, by = 0.025),
                   y = dtnorm(x, mean1 = par1, sd1 = par2, lo = 0, hi = Inf),
                   behavior = behavior,
                   track_length = track_length)
          })
true_plot_data_SL$behavior<- factor(true_plot_data_SL$behavior,
                                    levels = c("Encamped","ARS","Transit"))


# Plot overlapping densities per behavior (step lengths)
p.gamma<- ggplot(data = plot_data_SL, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_SL, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nStep Length (units)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  scale_y_continuous(breaks = c(0, 0.5, 1.00)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3)





### Make plot comparing HMM wrapped Cauchy distribs against true generating WC distrib
hmm.TA.df_weird<- hmm.TA.params2_weird %>% 
  bind_rows(., .id = "id") %>% 
  mutate(behavior = rep(c("Encamped","ARS","Transit"), 20),
         track_length = factor(rep(c(1000,5000,10000,50000), each = 15),
                               levels = c("1000","5000","10000","50000")))

true.TA.params_weird<- data.frame(id = 1:3,
                            TA.params,
                            behavior = c("Encamped","ARS","Transit"),
                            track_length = "True")


plot_data_TA<- 
  pmap_df(hmm.TA.df_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (pi/512)),
                   y = dwrappedcauchy(x, mu = circular(par1), rho = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
plot_data_TA$behavior<- factor(plot_data_TA$behavior, levels = c("Encamped","ARS","Transit"))


true_plot_data_TA<- 
  pmap_df(true.TA.params_weird,
          function(id, par1, par2, behavior, track_length) {
            tibble(id = id,
                   x = seq(-pi, pi, by = (pi/512)),
                   y = dwrappedcauchy(x, mu = circular(par1), rho = par2),
                   behavior = behavior,
                   track_length = track_length)
          })
true_plot_data_TA$behavior<- factor(true_plot_data_TA$behavior,
                                    levels = c("Encamped","ARS","Transit"))


# Plot overlapping densities per behavior (step lengths)
p.wc<- ggplot(data = plot_data_TA, aes(color = track_length)) +
  geom_line(aes(group = id, x = x, y = y), alpha = .4) + 
  geom_line(data = true_plot_data_TA, aes(group = id, x = x, y = y), color = "black",
            size = 0.9) +
  labs(x = "\nTurning Angle (radians)", y = "Density\n") +
  scale_color_brewer("Track Length", palette = "Dark2") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.75))) +
  facet_wrap(~behavior, nrow = 3)



## Create composite plot w/ shared legend
library(cowplot)

p.comp<- plot_grid(p.gamma + theme(legend.position="none"),
                   NA,
                   p.wc + theme(legend.position="none"),
                   align = 'v',
                   rel_widths = c(1, 0.1, 1),
                   hjust = -1,
                   nrow = 1)

# extract the legend from one of the plots
legend.comp<- get_legend(p.gamma + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
plot_grid(p.comp, legend.comp, ncol = 1, rel_heights = c(1, 0.1))

# ggsave("Figure 5.png", width = 7, height = 6, units = "in", dpi = 330)













### Extra figures for ESA presentation ###

theta.estim_df<- theta.estim[[7]] %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:7, names_to = "behavior", values_to = "prop") %>% 
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:7

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1), 
              alpha = 0.3) +
  # scale_color_brewer("", palette = "Dark2", guide = F) +
  # scale_fill_brewer("", palette = "Dark2", guide = F) +
  labs(x="\nBehavior", y="Proportion of Time Segment\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# ggsave("Theta boxplot.png", width = 6, height = 4, units = "in", dpi = 330)






qual_pal<- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"),
             brewer.pal(8, "Accent"), brewer.pal(12, "Paired"))
qual_pal<- my.cols(44)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



ggplot(data = dat2_focal %>% slice(2:n()), aes(x,y)) +
  geom_path(aes(color = factor(tseg)), size = 0.75) +
  scale_color_manual(values = sample(qual_pal, 44)) +
  # geom_point(fill = "grey45", pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = dat2_focal %>% slice(n=1), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat2_focal %>% slice(n()), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_void() +
  theme(axis.title = element_blank(),
        legend.position = "") +
  coord_cartesian()

# ggsave("Track segments.png", width = 6, height = 4, units = "in", dpi = 330)







ggplot(data = dat2_focal %>% slice(2:n()), aes(x,y)) +
  geom_path(aes(color = factor(behav_coarse), group = tseg), size = 0.75) +
  scale_color_brewer("", palette = "Dark2", breaks = 1:3,
                     labels = c("Encamped","ARS","Transit")) +
  geom_point(data = dat2_focal %>% slice(n=1), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat2_focal %>% slice(n()), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  theme_void() +
  theme(axis.title = element_blank(),
        legend.position = "") +
  coord_cartesian()

# ggsave("Track clustered.png", width = 6, height = 4, units = "in", dpi = 330)
