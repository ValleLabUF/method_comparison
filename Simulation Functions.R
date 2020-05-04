### CRW model

CRW.sim=function(nsim, ntseg, nstep, SL.params, TA.params, Z0) {  
  #nsim=number of simulations
  #ntseg=number of time segments
  #nstep=number of steps per time segment
  #nbehav=number of behavioral states ==> derive from nrow SL.params
  #SL.params=df of shape and scale params
  #TA.params=df of mean TA and concen. param
  #Z0=initial location
  
  #uses gamma and wrapped cauchy distribs
  #behaviors params must be in order
  #for simulating w/ 3 behavioral states
  
  track.list<- list()  #store all tracks of all durations
  brkpt.list<- list()  #store all true breakpoints
  
  for (k in 1:length(ntseg)) {
    
    behav<- list()  #to store behaviors for each time segment
    behav.all<- list()  #to store behaviors for each observation
    sim.list<- list()  #to store all sim tracks of same ntseg
    sim.brkpts<- list()  #to store all sim breakpoints of same ntseg
    
    #randomly choose n segments of each behavior to be 'pure' instead of mixed
    npure<- ifelse(ntseg[k] <= 10, 0,
                   ifelse(ntseg[k] > 10 & ntseg[k] <=50, 3,
                          ifelse(ntseg[k] > 50 & ntseg[k] <=100, 6,
                                 9)))
    
    for (j in 1:nsim) {
      
      #create vector of dominant behaviors per time segment (50 segments)
      behav[[j]]<- rmultinom(ntseg[k], 1, c(0.5, 0.30, 0.20)) %>% 
        apply(., 2, function(x) which(x == 1)) #change from matrix to vector
      # table(behav[[k]])/ntseg[k] #check freq
      
      #assign positions for pure time segments (if any)
      behav1.pure<- sample(which(behav[[j]]==1), npure, replace = FALSE)
      behav2.pure<- sample(which(behav[[j]]==2), npure, replace = FALSE)
      behav3.pure<- sample(which(behav[[j]]==3), npure, replace = FALSE)
      
      pure<- c(behav1.pure, behav2.pure, behav3.pure) %>% sort()
      
      
      #create vector of behaviors within each time segment
      behav.full<- list()
      for (i in 1:length(behav[[j]])) {
        if (i %in% pure) {
          behav.full[[i]]<- rep(behav[[j]][i], nstep)
        }else if (behav[[j]][i] == 1) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.8, 0.1, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 2) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.1, 0.8, 0.1))  %>% 
            apply(., 2, function(x) which(x == 1))
        } else if (behav[[j]][i] == 3) {
          behav.full[[i]]<- rmultinom(nstep, 1, c(0.1, 0.1, 0.8))  %>% 
            apply(., 2, function(x) which(x == 1))
        }
      }
      
      behav.all[[j]]<- unlist(behav.full)
      
      
      
      #create vector of step lengths
      SL<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          SL[[i]]<- rgamma(1, shape = SL.params[1,1], scale = SL.params[1,2])  #Rest
        } else if (behav.all[[j]][i] == 2) {
          SL[[i]]<- rgamma(1, shape = SL.params[2,1], scale = SL.params[2,2])  #Exploratory
        } else {
          SL[[i]]<- rgamma(1, shape = SL.params[3,1], scale = SL.params[3,2])  #Transit
        }
      }
      SL<- unlist(SL)
      
      
      #create vector of turning angles
      TA<- list()
      for (i in 1:length(behav.all[[j]])) {
        if (behav.all[[j]][i] == 1) {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[1,1]), rho=TA.params[1,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Rest
        } else if (behav.all[[j]][i] == 2) {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[2,1]), rho=TA.params[2,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Exploratory
        } else {
          TA[[i]]<- rwrappedcauchy(1, mu=circular(TA.params[3,1]), rho=TA.params[3,2]) %>%
            ifelse(. > pi, .-(2*pi), .)  #Transit
        }
      }
      TA<- unlist(TA)
      
      
      # cumulative angle
      Phi <- cumsum(TA)
      
      # step length components
      dX <- SL*cos(Phi)
      dY <- SL*sin(Phi)
      
      # actual X-Y values
      X <- c(Z0[1], Z0[1] + cumsum(dX))
      Y <- c(Z0[2], Z0[2] + cumsum(dY))
      track<- data.frame(id = as.character(paste0(j,"_",k)), x = X, y = Y, SL = c(NA,SL),
                         TA = c(NA, TA),
                         behav_fine = as.factor(c(NA, behav.all[[j]])),
                         behav_coarse = c(NA, rep(behav[[j]], each=nstep)) %>% factor(),
                         track_length = length(behav.all[[j]]))
      
      sim.list[[j]]<- track
      
      sim.brkpts[[j]]<- which(diff(behav[[j]]) != 0) * nstep
    }
    names(sim.list)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    names(sim.brkpts)<- sim.list %>% modify_depth(1, ~unique(.$id)) %>% unlist()
    
    track.list[[k]]<- sim.list
    brkpt.list[[k]]<- sim.brkpts
  }
  names(track.list)<- ntseg
  names(brkpt.list)<- ntseg
  
  list(tracks = track.list, brkpts = brkpt.list)
}

#----------------------------------------------

### BRW model (single phase)

BRW.sim=function (n = 50, a, b, Z.center = c(0,0), Z0 = c(3,4), rho) {
  #based on BCRW model described in Bailey et al. 2018 "Navigational efficiency in a biased and correlated random walk model of individual animal movement". Ecology. 99(1): 217-223.
  
  Z <- matrix(c(Z0, rep(NA, n*2 - 2)), n, 2, byrow = T)  #matrix of simulated locations
  ref.pt<- c(1,0)  #for reference vector
  
  #function to calc angle between 2 vectors
  angle <- function(x,y){
    dot.prod <- x%*%y 
    norm.x <- norm(x,type="2")
    norm.y <- norm(y,type="2")
    theta <- acos(dot.prod / (norm.x * norm.y))
    as.numeric(theta)
  }
  
  for (i in 2:n) {
    u<- as.matrix(c(Z.center-Z[i-1,]))  #vector formed by location and AC
    v<- as.matrix(ref.pt)  #reference vector for positive x-axis (representing 0 rad)
    
    #angle between loc-AC vector and the positive x-axis multiplied by the sign of the vector (x-component)
    if (u[1] < 0 & u[2] < 0) {
      ang.sign<- -1  #1st quadrant
    } else if (u[1] > 0 & u[2] < 0) {
      ang.sign<- -1  #2nd quadrant
    } else if (u[1] > 0 & u[2] > 0) {
      ang.sign<- 1  #3rd quadrant
    } else {
      ang.sign<- 1  #4th quadrant
    }
    
    #angle directed towards AC from location
    omega<- angle(t(u), v) * ang.sign
    if (is.na(omega))
      omega<- 0
    
    #generate navigation error (-pi < theta < pi) so that mean of cos(theta) is in [0,1]
    phi <- rnorm(1, 0, 0.5)
    
    #angle towards AC + error
    mu <- omega + phi
    if (mu < 0)
      mu<- 2 * pi + (omega + phi)
    theta <- rwrappedcauchy(1, mu = circular(mu), rho = rho) %>% as.numeric()
    
    # new step
    dist<- rgamma(1, a, b)  #step length
    dx <- dist * cos(theta)
    dy <- dist * sin(theta)
    
    # actual X-Y values
    x <- Z[i-1, 1] + dx
    y <- Z[i-1, 2] + dy
    Z[i,] <- c(x,y)
  }
  track <- data.frame(x = Z[,1], y = Z[,2])
  
  return(track)
}

#----------------------------------------------

### multiBRW model (multi-phase)

multiBRW.sim=function (Z0, Z.centers, n, nphases, a, b, rho, ...) {
  
  Z.list <- list()
  Z.centers<- as.matrix(Z.centers)
  
  #for visiting each Z.center (AC) once depending on nphases
  if (nphases > nrow(Z.centers) & nphases %% nrow(Z.centers) == 0) {
    
    fold<- ceiling(nphases / nrow(Z.centers))  #number of times to sample
    ind<- c(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = FALSE),
            rep(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = TRUE), fold-1))
    
  } else if (nphases > nrow(Z.centers) & nphases %% nrow(Z.centers) != 0) {
    
    fold<- ceiling(nphases / nrow(Z.centers))  #number of times to sample
    tmp<- nphases %% nrow(Z.centers)  #remainder
    ind<- c(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = FALSE),
            rep(sample(x = nrow(Z.centers), size = nrow(Z.centers), replace = TRUE), fold-2),
            sample(x = nrow(Z.centers), size = tmp, replace = TRUE))
    
  } else {
    
    ind<- sample(x = nrow(Z.centers), size = nphases, replace = FALSE)
    
  }
  
  for (i in 1:nphases) {
    if (i == 1) {
      Z0 = Z0
    } else {
      Z0 = Z.list[[i - 1]][n,] %>% as.numeric()
    }
    
    Z.list[[i]] <- BRW.sim(n = n, a = a, b = b, 
                           Z.center = Z.centers[ind[i],], 
                           Z0 = Z0, rho = rho)
  }
  
  track<- do.call(rbind.data.frame, Z.list)
  track$true.ac<- rep(ind, each = n)
  
  return(track)
}