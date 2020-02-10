#Compare KDE against modeled ACs

library(tidyverse)
library(adehabitatHR)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ks)
library(raster)
library(sf)
library(lubridate)


dat<- read.csv("Snail Kite Gridded Data_AC_TOHO.csv", as.is = TRUE)


###################
### Modeled ACs ###
###################

ac.coords<- read.csv("Activity Center Coordinates_TOHO.csv", header = T, sep = ',')

## Map

#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida")
fl<- sf::st_transform(fl, crs = "+init=epsg:32617") #change projection to UTM 17N

# lakes
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
lakes10<- sf::st_transform(lakes10, crs = "+init=epsg:32617") %>%
  sf::st_crop(xmin = min(dat$x-20000), xmax = max(dat$x+20000), ymin = min(dat$y-20000),
              ymax = max(dat$y+20000))

nests<- dat %>% group_by(id) %>% dplyr::select(c(id, x, y)) %>% slice(n=1)


# ACs
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_point(data = dat, aes(x, y, fill = "Raw"), shape = 21, size = 1,
             alpha = 0.2) +
  geom_point(data = nests, aes(x, y, fill = "Nests"), shape = 24, size = 3.5, alpha = 0.7) +
  geom_point(data = ac.coords, aes(x, y, fill = "ACs"), shape = 21, size = 3, alpha = 0.8) +
  labs(x="Longitude", y="Latitude") +
  scale_fill_manual("", values = c(viridis(n=3)[1],"red","grey55")) +
  guides(fill = guide_legend(override.aes = list(shape = c(21,24,21)))) +
  theme_bw()




###########
### KDE ###
###########

## H plug-in
dat.coords<- dat[,c("x","y")]
Hpi.dat<- Hpi(x=dat.coords, nstage = 2)
kde.hpi<- kde(x=dat.coords, H=Hpi.dat, compute.cont = T)
kd_df <- expand.grid(x=kde.hpi$eval.points[[1]], y=kde.hpi$eval.points[[2]]) %>% 
  mutate(z = c(kde.hpi$estimate)/sum(kde.hpi$estimate))

ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "black", alpha = 0.5) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_tile(data=kd_df, aes(x, y, fill=z), alpha=0.75) +
  geom_point(data = ac.coords, aes(x, y), color = viridis(n=9)[7], size = 2, pch=1, stroke=0.5) +
  scale_fill_viridis_c("Proportional \nDensity \nEstimate") +
  theme_bw() +
  labs(x="Longitude", y="Latitude")



############################
### ACs vs KDE over Time ###
############################


### ACs
dat2<- dat
dat2$date<- as.POSIXct(dat2$date, format = "%Y-%m-%d")
dat2<- dat2[order(dat2$date),]


ggplot() +
  geom_point(data = dat2, aes(x=date, y=ac)) +
  labs(x = "\nTime", y = "Activity Center\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = 1:20, trans = "reverse")


## AC Heatmap by month and year

dat2<- dat2 %>% mutate(month = lubridate::month(date), year = lubridate::year(date))
dat.sum<- dat2 %>% group_by(year, month, ac) %>% tally() %>% group_by(year, month) %>%
  mutate(N=sum(n)) %>% mutate(prop = n/N)
dat.sum$date<- as.Date(paste0(dat.sum$year,"-", dat.sum$month,"-01"), "%Y-%m-%d")

foo<- matrix(0, 20*length(unique(dat.sum$date)), 3)
colnames(foo)<- c("date","ac","prop")
foo[,1]<- rep(unique(dat.sum$date), each=20)
foo[,2]<- rep(1:20, length(unique(dat.sum$date)))
for (i in 1:nrow(dat.sum)) {
  ind<- which(dat.sum$date[i] == foo[,1] & dat.sum$ac[i] == foo[,2])
foo[ind,3]<- dat.sum$prop[i]
}
foo<- data.frame(foo)

ggplot(data = foo, aes(x=lubridate::as_date(date), y=ac)) +
  geom_tile(aes(fill=prop), width = 31) +
  scale_fill_viridis_c("Proportion of\nObservations\nper Month") +
  scale_y_continuous(trans = "reverse", breaks = 1:20, expand = c(0,0)) +
  scale_x_date(date_labels = "%b %Y", expand = c(0,0)) +
  labs(x="Time", y="Activity Center") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


### KDE (by season)

#function to calculate seasons
getSeason <- function(DATES) {
  WS <- as.Date("2012-12-21", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-20",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-20",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-22",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}

dat2<- dat2 %>% mutate(season = getSeason(dat2$date), year = year(dat2$date))

## determine grid dims for kde estimates
grid_5<- raster(extent(min(dat$x-5000), max(dat$x+5000),
                       min(dat$y-5000), max(dat$y+5000)))
res(grid_5)<- 5000
proj4string(grid_5)<- CRS("+init=epsg:32617")
grid_5[]<- 0
dim(grid_5)  #96 x 52


#calculate each kde
kde_season<- vector("list", 16)
years<- unique(dat2$year)
seasons<- unique(dat2$season)[c(4,1:3)]
oo=1

for (i in 1:length(years)) {
  for (j in 1:length(seasons)) {
    
    ## H plug-in
    dat.coords<- dat2 %>% filter(year == years[i] & season == seasons[j]) %>% dplyr::select(x,y)
    
    if (nrow(dat.coords) > 0) {
      Hpi.dat<- Hpi(x=dat.coords, nstage = 2)
      kde.hpi<- kde(x=dat.coords, H=Hpi.dat, compute.cont = T, gridsize = c(96,52), xmin = c(min(dat2$x), min(dat2$y)), xmax = c(max(dat2$x), max(dat2$y)))
      kd_df <- expand.grid(x=kde.hpi$eval.points[[1]], y=kde.hpi$eval.points[[2]]) %>% 
        mutate(z = c(kde.hpi$estimate)/sum(kde.hpi$estimate), season = seasons[j], year = years[i])
    } else {
      kd_df<- 0
    }
    
    kde_season[[oo]]<- kd_df
    
    oo=oo+1
  }
}


names(kde_season)<- paste(rep(seasons, length(years)), rep(years, each=length(seasons)))
kde_season<- kde_season[2:16]
kde_season<- map_dfr(kde_season, `[`)
kde_season$season<- factor(kde_season$season, levels = seasons)
kde_season$year<- factor(kde_season$year, levels = years)




#extract ACs w each associated time interval
ac_season<- vector("list", 16)
oo=1

for (i in 1:length(years)) {
  for (j in 1:length(seasons)) {
    
    ind<- dat2 %>% filter(year == years[i] & season == seasons[j]) %>% dplyr::select(ac) %>% unique()
    
    if (length(ind) > 0) {
      tmp<- ac.coords[ac.coords$ac %in% ind$ac,] %>% mutate(season = seasons[j], year = years[i])
    } else {
      tmp<- 0
    }
    
    ac_season[[oo]]<- tmp
    
    oo=oo+1
  }
}


ac_season<- ac_season[2:16]
ac_season<- map_dfr(ac_season, `[`)
ac_season$season<- factor(ac_season$season, levels = seasons)
ac_season$year<- factor(ac_season$year, levels = years)




lakes<- as(lakes10, "Spatial") %>% fortify()  #need to convert from sf so it shows up w facet

ggplot() +
  geom_sf(data = fl) +
  geom_polygon(data = lakes, aes(long, lat, group = group), fill = "black", alpha = 0.5) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_tile(data=kde_season, aes(x=x, y=y, fill=z), alpha=0.75) +
  # geom_point(data = ac_season, aes(x, y), color = viridis(n=9)[7], size = 2, pch=1, stroke=0.25) +
  scale_fill_viridis_c("Proportional \nDensity \nEstimate") +
  theme_bw() +
  theme(axis.text = element_blank(), strip.text = element_text(face = "bold")) +
  labs(x="Longitude", y="Latitude") +
  facet_grid(year~season)


