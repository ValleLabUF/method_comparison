library(bcpa)

set.seed(1)



d<-track

d$time <- seq(c(ISOdate(2020,3,20)), by = "hour", length.out = 5001)#create fake hourly times for BCPA

mytrack <- MakeTrack(d$x,d$y,d$time)
plot(mytrack)

Simp1 <- GetVT(mytrack)

tic()
Simp.ws <- WindowSweep(Simp1, "V*cos(Theta)", windowsize=80, progress=T, K=2)
toc()
#windowsize = # points to consider each time the model is run, k=2: default sensitivity (lower=less sensitive)
head(Simp.ws$ws)#show each breakpoint in time, the supported model, mu, s, and rho on each side of the track

cw=30 #clusterwidth= the number of times a changepoint must be identified in the moving window for it to count. 1=all changepoints.
cp<-ChangePointSummary(Simp.ws, clusterwidth=cw)
plot(Simp.ws, type="flat", clusterwidth=cw)
DiagPlot(Simp.ws)#normal?

cp #all the changepoints
y<-cbind(d[3:nrow(d),], Simp.ws$t)#number of lines in Simp.ws

names(y)[8]<-"t" #7 for hard clustering, 8 for mixed

y$group<-0
for (i in 1:1800){ #at least the number of changepoints
  ind1 = cp$phases$t0[i] 
  ind2 = cp$phases$t1[i] 
  test1 = (y$t>=ind1 & y$t<=ind2)
  y$group[test1]=i
}

library(viridis)
plot(Simp.ws, type="flat",xaxt="none",xlab=NA, clusterwidth=cw, ylab="Persistence velocity", main="CRW mixed membership simulation: behavioral changepoint analysis")

ggplot(y, aes(t,SL)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill = factor(group)), pch = 21, size = 2, alpha = 0.5, na.rm = T) +
  scale_fill_viridis("group",discrete=TRUE) +
  theme_bw()+ xlab("Timestep") + ylab("Step length (km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),legend.position = "none")+ 
  geom_vline(xintercept=c(cp$phases$t1), linetype="dashed", color = "red", size=0.7)

ggplot(y, aes(t,SL)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_fine), pch = 21, size = 2, alpha = 0.5, na.rm = T) +
  scale_fill_viridis_d("Behavior") +
  theme_bw() + xlab("Timestep") + ylab("Step length (km)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())



## Compare True vs Modeled Breakpoints
model.brkpts<- round(cp$breaks$middle, 0)
true.brkpts<- which(diff(as.numeric(as.factor(d$behav_coarse))) != 0) - 1  #for mixed-membership sim
all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts), type = rep(c("True","Model"),
                                                                        c(length(true.brkpts),
                                                                          length(model.brkpts))))

accuracy<- matrix(NA,length(model.brkpts),1)
for (i in 1:length(model.brkpts)) {  #assign brkpts as accurate or inaccurate
  
  tmp<- c(model.brkpts[i] - (20:0), model.brkpts[i] + (1:20)) %in% true.brkpts %>% sum()
  
  if (tmp == 0) {
    accuracy[i]<- "Inaccurate"
  } else {
    accuracy[i]<- "Accurate"
  }
}

if (sum(abs(diff(model.brkpts)) < 20) >= 0) {  #identify duplicate brkpts
  ind<- which(abs(diff(model.brkpts)) <= 20)
  ind<- sort(c(ind, ind+1))
}

ind.acc<- ind[which(accuracy[ind] == "Accurate")]
ind.inacc<- ind[which(accuracy[ind] == "Inaccurate")]
accuracy[ind.acc]<- "Accurate Duplicate"
accuracy[ind.inacc]<- "Inaccurate Duplicate"
accuracy<- c(rep("True",length(true.brkpts)), accuracy)


#identify missing breakpoints from model
status<- matrix(NA,length(true.brkpts),1)
for (i in 1:length(true.brkpts)) {
  
  tmp<- c(true.brkpts[i] - (50:0), true.brkpts[i] + (1:50)) %in% model.brkpts %>% sum()
  
  if (tmp == 0) {
    status[i]<- "Missing"
  } else {
    status[i]<- "Present"
  }
}

miss.ind<- which(status =="Missing")
status.miss<- data.frame(brks = true.brkpts[miss.ind], type = rep("Model", length(miss.ind)),
                         acc = rep("Missing", length(miss.ind)))


all.brkpts$acc<- accuracy
all.brkpts<- rbind(all.brkpts, status.miss)

# #for mixed-membership
ggplot(all.brkpts, aes(x=brks, y=type, color = acc, shape = acc)) +
  geom_point(size=3) +
  theme_bw() +
  labs(x="Time", y="Type") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10),
        legend.position = "top") +
  scale_color_manual("Accuracy", values = c("forestgreen","firebrick","firebrick","black")) +
  scale_shape_manual("Accuracy", values = c(16,16,4,16))


## Calculate percentage of each measure of breakpoint accuracy for summary

all.brkpts %>% filter(type == "Model") %>% group_by(acc) %>% tally() %>%
  mutate(freq = n/sum(n))




#clustering
d2<- cp$phases %>% dplyr::select(mu.hat, s.hat, rho.hat) %>% mutate(group = 1:nrow(.))


# K-Means Cluster Analysis
kmeans.res<- matrix(NA, 7, 1)
for (i in 1:7) {
fit<- kmeans(d2[,1:3], i) 
kmeans.res[i,]<- fit$betweenss/fit$totss
}
plot(kmeans.res)
#elbow appears to occur at k = 2

fit<- kmeans(d2[,1:3], 2)
fit ## cluster 1 = 'Resting/ARS', cluster 2 = 'transit'


# append cluster assignment
bcpa.clust<- data.frame(d2, cluster = factor(fit$cluster, levels = 1:2))
levels(bcpa.clust$cluster)<- c("Resting_ARS","Transit")
bcpa.dat<- merge(y, bcpa.clust[,c("group","cluster")], by= "group")
# write.csv(bcpa.dat, "Clustered BCPA data.csv", row.names = F)


## Plot of original sim track
ggplot(data = bcpa.dat, aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=behav_coarse), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))



## Plot of BCPA-Kmeans results
ggplot(data = bcpa.dat, aes(x,y)) +
  geom_path(color = "gray75") +
  geom_point(aes(fill=cluster), pch = 21, size = 2.5, alpha = 0.7) +
  geom_point(data = track[1,], aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = track[nrow(track),], aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  coord_equal() +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16)) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))



## Validation

# Overall
true.b.coarse<- bcpa.dat$behav_coarse %>% as.numeric()
true.b.coarse[true.b.coarse == 2]<- 1
model.b<- bcpa.dat$cluster %>% as.numeric()
model.b[model.b == 2]<- 3

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 96.7%  #for sim track from multinom where 'resting' and 'ARS' merged

# For 'Resting/ARS' behavior
true.b.coarse_rest<- which(true.b.coarse == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.coarse_rest %in% model.b_rest) %>% length()) / length(true.b.coarse_rest)
# 96.3% accuracy for 'Resting/ARS' compared to 'Resting/ARS'

# For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 98.6% accuracy for 'Transit' at coarse scale