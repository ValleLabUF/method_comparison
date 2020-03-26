library(bcpa)
setwd("C:\\Users\\cpoli\\Dropbox\\Dissertation\\Behavioral partitioning work")

d<-track

d$time <- seq(c(ISOdate(2020,3,20)), by = "hour", length.out = 5001)#create fake hourly times for BCPA

mytrack <- MakeTrack(d$x,d$y,d$time)
plot(mytrack)

Simp1 <- GetVT(mytrack)

Simp.ws <- WindowSweep(Simp1, "V*cos(Theta)", windowsize=80, progress=T, K=2)#windowsize = # points to consider each time the model is run, k=2: default sensitivity (lower=less sensitive)
head(Simp.ws$ws)#show each breakpoint in time, the supported model, mu, s, and rho on each side of the track

cw=1 #clusterwidth= the number of times a changepoint must be identified in the moving window for it to count. 1=all changepoints.
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



#clustering
library(doBy)
SL<-summaryBy(SL ~ group, data = y, FUN = list(mean))#also try FUN= median
TA<-summaryBy(TA ~ group, data = y, FUN = list(mean))
d<-as.data.frame(na.omit(scale(cbind(SL[,2],TA[,2]))))#try scaled and not scaled


# K-Means Cluster Analysis
fit <- kmeans(d, 3) # 3 cluster solution
# get cluster means
aggregate(d,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(d, fit$cluster)

colnames(mydata)<-c("SL","TA","BCPA.fine.cluster")
write.csv(merge(x, y, by= "BCPA.fine.cluster"),"clustered segments kmeans.csv" )



# Ward Hierarchical Clustering
dt <- dist(d, method = "euclidean") # distance matrix
fit <- hclust(dt, method="ward.D")
plot(fit) # display dendogram
groups <- cutree(fit, k=3) # cut tree into 3 clusters
# draw dendogram with red borders around the 3 clusters
rect.hclust(fit, k=3, border="red")
x<-cbind(groups, d[1])
colnames(x)<-c("BCPA.fine.cluster","SL")
write.csv(merge(x, y, by= "BCPA.fine.cluster"),"clustered segments Ward.csv" )

