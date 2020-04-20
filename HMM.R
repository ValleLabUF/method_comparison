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
library(tictoc)

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
stateNames <- c("Resting","ARS","Transit")

### Test with 2-4 behavioral states and then perform order selection via AIC/BIC

set.seed(12345)
# Empty list for order selection
k.models<- list()


tic()
## K = 2

# Set number of iterations, save list of fitted models, and name behavioral states
niter <- 50
allm <- list()
stateNames <- c("Encamped","Exploratory")
for(i in 1:niter) {
  print(i)
  
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.1, 5),
                     max = c(2, 20))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.1, 2),
                   max = c(1, 8))
  # Turning angle mean
  angleMean0 <- c(pi, 0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.4, 0.5),
                     max = c(0.99, 0.99))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0,angleCon0)
  allm[[i]] <- fitHMM(data = move.d, nbStates = 2, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = list(step = "gamma", angle = "wrpcauchy"),
                      formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                      estAngleMean = list(angle=TRUE),
                      stateNames = stateNames)
}

# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)

# Best fitting model
k.models[[1]] <- allm[[whichbest]]







## K = 3

# Set number of iterations, save list of fitted models, and name behavioral states
niter <- 50
allm <- list()
stateNames <- c("Resting","ARS","Transit")

for(i in 1:niter) {
  print(i)
  
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
  anglePar0 <- c(angleMean0,angleCon0)
  allm[[i]] <- fitHMM(data = move.d, nbStates = 3, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = list(step = "gamma", angle = "wrpcauchy"),
                      formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                      estAngleMean = list(angle=TRUE),
                      stateNames = stateNames)
}


# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)

# Best fitting model
k.models[[2]] <- allm[[whichbest]]








## K = 4

# Set number of iterations, save list of fitted models, and name behavioral states
niter <- 50
allm <- list()
stateNames <- c("Resting","ARS","Exploratory","Transit")

for(i in 1:niter) {
  print(i)
  
  # Step length mean
  stepMean0 <- runif(4,
                     min = c(0.1, 0.5, 3, 8),
                     max = c(2, 5, 10, 20))
  # Step length standard deviation
  stepSD0 <- runif(4,
                   min = c(0.1, 1, 2, 4),
                   max = c(1, 2, 4, 8))
  # Turning angle mean
  angleMean0 <- c(pi, pi, 0, 0)
  # Turning angle concentration
  angleCon0 <- runif(4,
                     min = c(0.4, 0.01, 0.01, 0.5),
                     max = c(0.99, 0.5, 0.5, 0.99))
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0,angleCon0)
  allm[[i]] <- fitHMM(data = move.d, nbStates = 4, 
                      Par0 = list(step = stepPar0, angle = anglePar0),
                      dist = list(step = "gamma", angle = "wrpcauchy"),
                      formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                      estAngleMean = list(angle=TRUE),
                      stateNames = stateNames)
}


# Extract likelihoods of fitted models
allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))

# Index of best fitting model (smallest negative log-likelihood)
whichbest <- which.min(allnllk)

# Best fitting model
k.models[[3]] <- allm[[whichbest]]

toc()


## View output from each model and inspect pseudo-residuals
plot(k.models[[1]]); plotPR(k.models[[1]])
plot(k.models[[2]]); plotPR(k.models[[2]])
plot(k.models[[3]]); plotPR(k.models[[3]])



## Make inference via AIC
k.2<- k.models[[1]]  # 2 states
k.3<- k.models[[2]]  # 3 states  
k.4<- k.models[[3]]  # 4 states

AIC(k.2, k.3, k.4)
AICweights(k.2, k.3, k.4)  #4 states is far and away the best model


## Make inference via BIC
# BIC = -2*logL + p*log(T)

#K=2
(2*k.2$mod$minimum) + (length(k.2$mod$estimate)*log(nrow(da)))  #25731.71
#K=3
(2*k.3$mod$minimum) + (length(k.3$mod$estimate)*log(nrow(da)))  #23150.92
#K=4
(2*k.4$mod$minimum) + (length(k.4$mod$estimate)*log(nrow(da)))  #22829.74; BEST


### Sticking w/ K=3

#examine, plot
plot(k.3)
plotStates(k.3, ask = FALSE)

ss<-viterbi(k.3)[1:5000]  #extract states
fin<-cbind(da, hmm = c(NA, ss))
write.csv(fin,"Sim2 mixed hmm.csv")
