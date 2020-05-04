#Running this script is dependent on having run all sets of analyses from generation of simulation to clustering via LDA and plotting results


## Sim Track Map
ggplot(data = track[-1,], aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(fill = "gray45", pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank())

ggsave("Figure 1a.png", width = 6, height = 4, units = "in", dpi = 330)




## Discretized Movmt Param Histograms
bin<- c(1:5, 1:8)
param.prop$bin<- bin
param.prop$key<- factor(param.prop$key)
levels(param.prop$key)<- c("Step Length", "Turning Angle")

ggplot(data = param.prop, aes(bin, prop)) +
  geom_bar(stat = "identity", fill = "grey45", color = "black") +
  labs(x = "Bin", y = "Proportion") +
  theme_bw() +
  facet_wrap(~key, nrow = 2, scales = "free_x") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 14, face = "bold"))

ggsave("Figure 1b.png", width = 6, height = 4, units = "in", dpi = 330)




## Segmentation Heatmap
plot.heatmap(data = list(behav.list[[1]][1:2500,]), nbins = c(5,8), brkpts = brkpts[,1:18],
             dat.res = dat.res, type = "behav", title = F, legend = T)

nbins = c(5,8)
data<- behav.list[[1]]
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
    labs(x = "Observations", y = "Bin", title = paste("ID",
                                                      unique(data$id))) +
    theme_bw() +
    theme(axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 11),
          strip.text = element_text(size = 14, face = 'bold'),
          plot.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 14),
          plot.margin = margin(5,20,5,5))

ggsave("Figure 1c.png", width = 6, height = 4, units = "in", dpi = 330)
  
  
  
  
## Behavior Histograms
behav.res$param<- factor(behav.res$param)
levels(behav.res$param)<- c("Step Length", "Turning Angle")

ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
    geom_bar(stat = 'identity') +
    labs(x = "Bin", y = "Proportion") +
    theme_bw() +
    theme(axis.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x.bottom = element_text(size = 14),
          strip.text = element_text(size = 14),
          strip.text.x = element_text(face = "bold"),
          panel.grid = element_blank()) +
    scale_fill_manual(values = viridis(n=3), guide = F) +
    facet_grid(param ~ behav, scales = "fixed")

ggsave("Figure 1d.png", width = 6, height = 4, units = "in", dpi = 330)






## Time Series Behavior Proportions
ggplot(theta.estim.long) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "Observation", y = "Proportion of Behavior") +
  scale_fill_viridis_d("") +
  theme_bw() + scale_x_continuous(expand = c(0,0)) +scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.bottom = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14),
        plot.margin = margin(5,20,5,5))

ggsave("Figure 1e.png", width = 6, height = 4, units = "in", dpi = 330)




## Sim Track Annotated Map
ggplot(data = track[-1,], aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 14)))

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

#plot all tracks
ggplot() +
  geom_sf(data = fl, fill = "grey80", color = NA) +
  geom_sf(data = waterbody_fl, fill = "grey45", color = NA) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  scale_color_viridis_d("", guide=F) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey45"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10))

setwd("~/Documents/Manuscripts/Bayesian Behavior Estimation Model/Figures")
ggsave("Figure 2.png", width = 4, height = 6, units = "in", dpi = 330)
