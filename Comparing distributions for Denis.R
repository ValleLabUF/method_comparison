library(tidyverse)
library(gridExtra)

#Load function to estimate proportions per bin
extract.behav.props=function(dat, id1, behav.col, dist.bin.lims, angle.bin.lims) {
  #only defined for step lengths ('dist') and turning angles ('rel.angle')
  #id must be in quotes
  #behav.col must be col name in quotes
  #dist.bin.lims and angle.bin.lims are numeric vectors
  
  #number of bins for both params are already defined as 5 (SL) and 8 (TA)
  #order of states is set as encamped, ARS, transit
  
  tmp1<- dat %>% filter(id == id1)
  behav.res<- list()
  behav.col<- rlang::sym(behav.col)
  
  for (j in 1:max(tmp1 %>% dplyr::select(!!behav.col), na.rm = T)) {
    tmp2<- tmp1[behav.col == j,]
    tmp2<- tmp1 %>% dplyr::filter(!!behav.col == j)
    
    SL.props<- vector()
    TA.props<- vector()
    for (i in 2:length(dist.bin.lims)) {
      SL.props[i-1]<- length(which(tmp2$dist < dist.bin.lims[i] & 
                                     tmp2$dist > dist.bin.lims[i-1])) / nrow(tmp2)
    }
    for (i in 2:length(angle.bin.lims)) {
      TA.props[i-1]<- length(which(tmp2$rel.angle < angle.bin.lims[i] & 
                                     tmp2$rel.angle > angle.bin.lims[i-1])) / nrow(tmp2)
    }
    behav.res[[j]]<- data.frame(behav = j, param = rep(c("Step Length","Turning Angle"), c(5,8)),
                                bin = c(1:5,1:8), prop = c(SL.props, TA.props))
  }
  
  b<- bind_rows(behav.res)
  b$behav<- b$behav %>% 
    as.character() %>% 
    str_replace_all(., "1", "Encamped") %>% 
    str_replace_all(., "2", "ARS") %>% 
    str_replace_all(., "3", "Transit") %>% 
    factor(., levels = c("Encamped","ARS","Transit"))
  
  b
}





# Load results
bayes.res<- read.csv("Modeled MM Sim Tracks w Behav.csv")
hmm.res<- read.csv("HMM results.csv")

# Assign identifiers by method and make consistent behavior colname
bayes.res$method<- rep("Bayesian", nrow(bayes.res))
hmm.res$method<- rep("HMM", nrow(hmm.res))

# Rename variables
bayes.res<- bayes.res %>% 
  rename(state = behav) %>% 
  mutate_at("state", ~factor(., levels = c("Encamped","ARS","Transit"))) %>%
  mutate_at("state", as.numeric)
hmm.res<- hmm.res %>% 
  rename(state = hmm.state, id = ID)



# Define bin limits for step lengths and turning angles
angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins

max.dist=max(bayes.res[bayes.res$dt == 3600,]$dist, na.rm = T)
dist.bin.lims=quantile(bayes.res[bayes.res$dt == 3600,]$dist, c(0,0.25,0.50,0.75,0.90), na.rm=T)
dist.bin.lims=c(dist.bin.lims, max.dist)  #5 bins




# True bin proportions for track 3_1 of 1000 observations
true.b<- extract.behav.props(dat = dat, id1 = '3_1', behav.col = "behav_coarse",
                             dist.bin.lims = dist.bin.lims,
                             angle.bin.lims = angle.bin.lims)

### HMM (Gamma/Wrapped Cauchy)

names(hmm.res)[4:5]<- c("dist","rel.angle")

hmm.b<- extract.behav.props(dat = hmm.res, id1 = "3_1", behav.col = "state",
                            dist.bin.lims = dist.bin.lims,
                            angle.bin.lims = angle.bin.lims)


p.hmm<- ggplot(hmm.b, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = true.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
             color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n", title = "HMM (Gamma/Wrapped Cauchy)") +
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
  facet_grid(behav ~ param, scales = "free_x")




### Bayesian
#calculate true proportions of SL and TA by behavior for all bins for ID 3_1 for Bayesian model
bayes.b<- structure(list(bin = c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 
                                 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 
                                 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
                         behav = c("1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "3", "3",
                                   "3", "3", "3", "1", "1", "1", "1", "1", "1", "1", "1", "2",
                                   "2", "2", "2", "2", "2", "2", "2", "3", "3", "3", "3", "3",
                                   "3", "3", "3"), 
                         prop = c(0.0244181355537063, 0.237450916437443, 0.567028087396184, 
                                  0.128523433602, 0.042579427010667, 0.560148249880302, 0.373240762652426, 
                                  0.00352772121478304, 0.0203461111900491, 0.0427371550624402, 
                                  0.00162533866151291, 0.00148855013640278, 0.00192196940763541, 
                                  0.568813957669013, 0.426150184125436, 0.144805239025516, 
                                  0.0518179415483269, 0.145203395099022, 0.136083900910635, 
                                  0.127084752481579, 0.137597061525012, 0.137584878570597, 
                                  0.119822830839312, 0.45867037096118, 0.0583719833638732, 
                                  0.00168940006871556, 0.0333529157157388, 0.0348198743526581, 
                                  0.000905693133189075, 0.00773201568027708, 0.404457746724369, 
                                  0.00608267068668819, 0.0215055457943176, 0.0148851264414491, 
                                  0.414490988654745, 0.504125406262491, 0.0234112521143246, 
                                  0.00375547903500086, 0.0117435310109843),
                         param = c("Step Length", "Step Length", "Step Length", "Step Length",
                                   "Step Length", "Step Length", "Step Length", "Step Length",
                                   "Step Length", "Step Length", "Step Length", "Step Length",
                                   "Step Length", "Step Length", "Step Length", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle", "Turning Angle",
                                   "Turning Angle", "Turning Angle")),
                    row.names = c(NA, -39L), class = c("tbl_df", "tbl", "data.frame"))


bayes.b$behav<- bayes.b$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "ARS") %>% 
  str_replace_all(., "2", "Encamped") %>%
  str_replace_all(., "3", "Transit") %>%
  factor(., levels = c("Encamped","ARS","Transit"))



p.bayes<- ggplot(bayes.b, aes(x = bin, y = prop, fill = behav)) +
  geom_bar(stat = 'identity') +
  geom_point(data = true.b, aes(x=bin, y=prop, group = behav), pch=21, size = 2, fill="white",
             color="black", stroke=1) +
  labs(x = "\nBin", y = "Proportion\n", title = "Bayesian") +
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
  facet_grid(behav ~ param, scales = "free_x")



# Plot distributions side by side
grid.arrange(p.bayes, p.hmm, ncol = 2)
