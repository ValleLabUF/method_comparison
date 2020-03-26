###############################################################
###############################################################
#Movement ecology introduction
###############################################################
###############################################################

library(raster)        #for raster manipulation
library(rgdal)         #for shp data, projections; version 1.3-4 used
library(adehabitatLT)  #for interpreting trajectories
library(momentuHMM)
library(dplyr)

#------------------------------#
#Hidden Markov model
#------------------------------#

da<-track
da$time <- seq(c(ISOdate(2020,3,20)), by = "hour", length.out = 5001)#create fake hourly times for BCPA

da$time<- as.POSIXct(da$time, format = "%Y/%m/%d %H:%M:%S")
utmcoord <- SpatialPoints(da[,1:2])
da$x<-attr(utmcoord,"coords")[,1]
da$y<-attr(utmcoord,"coords")[,2]

#format
move.d <- prepData(data.frame(ID=1, x = da$x, y = da$y), 
                   type="UTM", coordNames = c("x","y"))#assigned an id of 1 because it's one track
#check data
plot(y~x, data=move.d, type="o", col="gray65", asp=1, cex=0, pch=19, xlab="Longitude",ylab="Latitude")

#choose starting values from data
hist(move.d$step, col = 'skyblue4', breaks = 50)
mean(move.d$step, na.rm=T);sd(move.d$step, na.rm=T)
acf(move.d$step[!is.na(move.d$step)])
hist(move.d$angle, col = 'purple', breaks = 50)

#zero distances-needed for 'zermass0' starting values
whichzero <- which(move.d$step == 0)
propzero <- length(whichzero)/nrow(move.d)
stateNames <- c("Resting","Exploratory","Transit")

set.seed(12345)
# Number of tries with different starting values
niter <- 50
# Save list of fitted models
allm <- list()
for(i in 1:niter) {
  # Step length mean
  stepMean0 <- runif(3,
                     min = c(0.1, 0.5, 5),
                     max = c(2, 5, 20))
  # Step length standard deviation
  stepSD0 <- runif(3,
                   min = c(0.1, 1, 2),
                   max = c(1, 2, 8))
  # Turning angle mean
  angleMean0 <- c(pi, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(3,
                     min = c(0.4, 0.01, 0.5),
                     max = c(0.99, 0.3, 0.99))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0,angleCon0)# 
  allm[[i]] <- fitHMM(data = move.d, nbStates = 3, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = list(step = "gamma", angle = "wrpcauchy"),
                      formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                      estAngleMean = list(angle=TRUE),
                      stateNames = stateNames)
}


# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
allnllk

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)

# Best fitting model
mbest <- allm[[whichbest]]
mbest

#examine, plot
print(mbest)
plot(mbest)
plotStates(mbest, ask = FALSE)

ss<-viterbi(mbest)#extract states
fin<-cbind(da,ss)
write.csv(fin,"Sim2 mixed hmm.csv")
