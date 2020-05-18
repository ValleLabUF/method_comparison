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




## Discretized Movmt Param Histograms
dat2_focal_long<- dat2_focal %>% 
  dplyr::select(id,SL,TA) %>% 
  gather(key, value, -id) %>% 
  group_by(key, value) %>%
  summarise(n=n()) %>% 
  drop_na() %>% 
  mutate(prop=n/sum(n)) %>%
  ungroup()

names(dat2_focal_long)[1:2]<- c("param","bin")
dat2_focal_long$param<- factor(dat2_focal_long$param)
levels(dat2_focal_long$param)<- c("Step Length", "Turning Angle")

ggplot(data = dat2_focal_long, aes(bin, prop)) +
  geom_bar(stat = "identity", fill = "grey45", color = "black") +
  labs(x = "Bin", y = "Proportion") +
  theme_bw() +
  facet_wrap(~param, nrow = 2, scales = "free_x") +
  scale_x_continuous(breaks = 1:8) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))

ggsave("Figure 1b.png", width = 6, height = 4, units = "in", dpi = 330)




## Segmentation Heatmap
# plot.heatmap(data = list(behav.list$`2_2`[1:5000,]), nbins = c(5,8), brkpts = brkpts[7,],
             dat.res = dat.res, type = "behav", title = F, legend = T)

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

ggsave("Figure 1c.png", width = 7, height = 5, units = "in", dpi = 330)
  
  
  
  
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

ggsave("Figure 1d.png", width = 6, height = 4, units = "in", dpi = 330)






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

ggsave("Figure 1e.png", width = 6, height = 4, units = "in", dpi = 330)




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

ggsave("Figure 1f.png", width = 6, height = 4, units = "in", dpi = 330)








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
data<- behav.list$`2_2`  #ID 2_2
nbins<- c(5,8)
behav.heat<- behav.seg.image(data, nbins)

SL<- data.frame(behav.heat$SL)
names(SL)<- 1:nbins[1]
SL<- SL %>% gather(key, value) %>% mutate(time=rep(data$time1, times=nbins[1]),
                                          behav=rep("Step Length", nrow(data)*nbins[1]))

TA<- data.frame(behav.heat$TA)
names(TA)<- 1:nbins[2]
TA<- TA %>% gather(key, value) %>% mutate(time=rep(data$time1, times=nbins[2]),
                                          behav=rep("Turning Angle", nrow(data)*nbins[2]))

behav.heat_long<- rbind(SL,TA)
behav.heat_long$value<- factor(behav.heat_long$value)
levels(behav.heat_long$value)<- c("Unoccupied","Occupied")

ind=which(unique(data$id) == brkpts$id)
breakpt<- brkpts[ind,-1] %>% purrr::discard(is.na) %>% t() %>% data.frame()
names(breakpt)<- "breaks"



p.seg<- ggplot(behav.heat_long, aes(x=time, y=key, fill=value)) +
  geom_tile(alpha = 0.5) +
  facet_wrap(~behav, scales = 'free', nrow = 2) +
  scale_fill_viridis_d('') +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_vline(data = breakpt, aes(xintercept = breaks - 0.5), color = viridis(n=9)[7],
             size = 1, alpha = 1) +
  geom_vline(data = data.frame(t(true.brkpts[7,-1])) %>% drop_na() %>% rename(breaks = X7),
             aes(xintercept = breaks - 0.5), color = "blue", size = 0.4, alpha = 1) +
  labs(x = "\nTime", y = "Bin\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = 'bold'),
        plot.title = element_text(size = 20, hjust = 0, vjust = -6),
        plot.margin = margin(0, 1, 0.5, 0.5, "cm"),
        panel.grid = element_blank(),
        legend.justification = "right",
        legend.position = "top",
        legend.text = element_text(
          margin = margin(r = 15, unit = "pt")))




## Part b: behavior histogram

#calculate true proportions of SL and TA by behavior for all bins for ID 2_2
tmp1<- dat %>% filter(id == "2_2")
b.list<- list()

for (j in 1:max(tmp1$behav_coarse, na.rm = T)) {
  tmp2<- tmp1 %>% filter(behav_coarse == j)
  
  true.SL.props<- vector()
  true.TA.props<- vector()
  for (i in 2:length(dist.bin.lims)) {
    true.SL.props[i-1]<- length(which(tmp2$dist < dist.bin.lims[i] & 
                                        tmp2$dist > dist.bin.lims[i-1])) / nrow(tmp2)
  }
  for (i in 2:length(angle.bin.lims)) {
    true.TA.props[i-1]<- length(which(tmp2$rel.angle < angle.bin.lims[i] & 
                                        tmp2$rel.angle > angle.bin.lims[i-1])) / nrow(tmp2)
  }
  b.list[[j]]<- data.frame(behav = j, param = rep(c("Step Length","Turning Angle"), c(5,8)),
                           bin = c(1:5,1:8), prop = c(true.SL.props, true.TA.props))
}

true.b<- bind_rows(b.list)
true.b$behav<- true.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

behav.res<- behav.res[[7]]
behav.res$param<- str_replace_all(behav.res$param, "SL", "Step Length") %>% 
  str_replace_all(., "TA", "Turning Angle")
behav.res$behav<- behav.res$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

p.hist<- ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  geom_point(data = true.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
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
  facet_grid(behav ~ param, scales = "free_x")




##Part c: generate time series plots comparing behavior proportions
dat1<- theta.estim.long[[7]]
dat1$behavior<- dat1$behavior %>% 
  as.character() %>% 
  str_replace_all(., "1", "Encamped") %>% 
  str_replace_all(., "2", "ARS") %>% 
  str_replace_all(., "3", "Transit") %>% 
  factor(., levels = c("Encamped","ARS","Transit"))

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
            size = 1) +
  scale_color_manual(values = viridis(n=20)[c(1,9,18)], guide=F) +
  new_scale_color() +
  geom_line(data = true.dat1,
            aes(x=time1, y=prop, color = behavior),
            size = 0.55) +
  scale_color_manual(values = viridis(n=20)[c(7,13,20)], guide=F) +
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
grid.arrange(p.seg, p.hist, p.prop, heights = c(0.2, 1, 0.1, 1), layout_matrix = rbind(c(NA, NA),
                                                                                       c(1, 2),
                                                                                       c(NA, NA),
                                                                                       c(3, 3)))
# dev.off()

# ggsave("Figure 4 (results from sim).png", width = 4, height = 6, units = "in", dpi = 330)