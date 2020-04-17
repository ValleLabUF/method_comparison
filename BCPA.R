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
plot(Simp.ws, type="flat", clusterwidth=cw, mu.where = "topright")
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



#clustering
d2<- cp$phases %>% dplyr::select(mu.hat, s.hat, rho.hat) %>% mutate(group = 1:nrow(.))


# K-Means Cluster Analysis
kmeans.res<- matrix(NA, 7, 1)
for (i in 1:7) {
fit<- kmeans(d2[,1:3], i) 
kmeans.res[i,]<- fit$betweenss/fit$totss
}
plot(kmeans.res)
#elbow appears to occur at k =3

fit<- kmeans(d2[,1:3], 3)
fit ## cluster 1 = 'resting', cluster 2 = 'ARS', cluster 3 = 'transit'


# append cluster assignment
bcpa.clust<- data.frame(d2, cluster = factor(fit$cluster, levels = 1:3))
levels(bcpa.clust$cluster)<- c("Resting","ARS","Transit")
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
model.b<- bcpa.dat$cluster %>% as.numeric()

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 59.4% accuracy when including all different behaviors together at coarse scale


# For 'Resting' behavior
true.b.coarse_rest<- which(true.b.coarse == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.coarse_rest %in% model.b_rest) %>% length()) / length(true.b.coarse_rest)
# 36.1% accuracy for 'Resting' at coarse scale

# For 'ARS' behavior
true.b.coarse_exp<- which(true.b.coarse == 2)
model.b_exp<- which(model.b == 2)
(which(true.b.coarse_exp %in% model.b_exp) %>% length()) / length(true.b.coarse_exp)
# 74.3% accuracy for 'ARS' at coarse scale

# For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 97.7% accuracy for 'Transit' at coarse scale