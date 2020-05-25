
run.HMMs2 = function(list1, elements) {
  hmm.res<- list()
  
  for (j in elements) {
    start.time<- Sys.time()
    
    # Empty list for order selection
    k.models<- list()
    
    
    ## K = 2
    
    allm<- list()
    niter<- 30
    stateNames <- c("Encamped","Exploratory")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=2"))
      print(paste("Iteration", i))
      
      # Step length shape param
      stepShape0 <- runif(2,
                         min = c(0.1, 1),
                         max = c(0.8, 5))
      # Step length scale param
      stepScale0 <- runif(2,
                       min = c(0.1, 0.5),
                       max = c(0.8, 2))
      # Turning angle mean
      angleMean0 <- c(pi, 0)
      # Turning angle concentration
      angleCon0 <- c(5, 5)
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepShape0, stepScale0, zeromass0)
      } else {
        stepPar0 <- c(stepShape0, stepScale0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 2, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "weibull", angle = "vm"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[1]] <- allm[[whichbest]]
    
    
    
    
    ## K = 3
    
    allm<- list()
    niter<- 30
    stateNames <- c("Resting","ARS","Transit")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, 0, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=3"))
      print(paste("Iteration", i))
      
      # Step length shape param
      stepShape0 <- runif(3,
                         min = c(0.1, 1, 2),
                         max = c(0.5, 2, 5))
      # Step length scale param
      stepScale0 <- runif(3,
                       min = c(0.01, 1, 8),
                       max = c(0.2, 5, 13))
      # Turning angle mean
      angleMean0 <- c(pi, 0, 0)
      # Turning angle concentration
      angleCon0 <- runif(3,
                         min = c(1, 0.01, 1),
                         max = c(8, 0.3, 8))
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepShape0, stepScale0, zeromass0)
      } else {
        stepPar0 <- c(stepShape0, stepScale0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 3, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "weibull", angle = "vm"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")  
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[2]] <- allm[[whichbest]]
    
    
    
    
    
    ## K = 4
    
    allm<- list()
    niter<- 30
    stateNames <- c("Resting","ARS","Exploratory","Transit")
    whichzero <- which(list1[[j]]$step == 0)
    propzero <- length(whichzero)/nrow(list1[[j]])
    zeromass0 <- c(propzero, propzero, 0, 0)        #for zero distances by state
    
    for (i in 1:niter) {
      print(paste("Simulation", j))
      print(paste("K=4"))
      print(paste("Iteration", i))
      
      # Step length shape param
      stepShape0 <- runif(4,
                          min = c(0.1, 0.5, 1, 2),
                          max = c(0.5, 1, 2, 5))
      # Step length scale param
      stepScale0 <- runif(4,
                          min = c(0.1, 1, 5, 5),
                          max = c(0.6, 3, 15, 15))
      # Turning angle mean
      angleMean0 <- c(pi, pi, 0, 0)
      # Turning angle concentration
      angleCon0 <- runif(4,
                         min = c(0.5, 0.01, 0.5, 0.5),
                         max = c(0.99, 0.5, 0.9, 0.99))
      # Fit model
      if(propzero > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepShape0, stepScale0, zeromass0)
      } else {
        stepPar0 <- c(stepShape0, stepScale0)
      }
      anglePar0 <- c(angleMean0, angleCon0)
      allm[[i]] <- fitHMM(data = list1[[j]], nbStates = 4, 
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          dist = list(step = "weibull", angle = "vm"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = stateNames,
                          optMethod = "Nelder-Mead")  
    }
    
    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(allm, function(m) m$mod$minimum))
    
    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)
    
    # Best fitting model
    k.models[[3]] <- allm[[whichbest]]
    
    
    
    end.time<- Sys.time()
    elapsed.time<- difftime(end.time, start.time, units = "min")
    
    hmm.res[[j]]<- list(models = k.models, elapsed.time = elapsed.time)
  }
  
  hmm.res
}
