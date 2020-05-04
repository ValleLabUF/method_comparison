
brkpt.accuracy<- function(model.brkpts, true.brkpts, acc.tol, dup.tol, miss.tol) {
  
  model.brkpts<- as.numeric(model.brkpts) %>% na.omit() %>% as.vector()
  true.brkpts<- as.numeric(true.brkpts) %>% na.omit() %>% as.vector()
  all.brkpts<- data.frame(brks = c(true.brkpts, model.brkpts),
                          type = rep(c("True","Model"),c(length(true.brkpts),
                                                         length(model.brkpts))))
  
  accuracy<- matrix(NA, length(model.brkpts), 1)
  for (i in 1:length(model.brkpts)) {  #assign brkpts as accurate or inaccurate
    
    tmp<- c(model.brkpts[i] - (acc.tol:0), model.brkpts[i] + (1:acc.tol)) %in% true.brkpts %>%
      sum()
    
    if (tmp == 0) {
      accuracy[i]<- "Inaccurate"
    } else {
      accuracy[i]<- "Accurate"
    }
  }
  
  if (sum(abs(diff(model.brkpts)) <= dup.tol) >= 0) {  #identify duplicate brkpts
    ind<- which(abs(diff(model.brkpts)) <= dup.tol)
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
    
    tmp<- c(true.brkpts[i] - (miss.tol:0), true.brkpts[i] + (1:miss.tol)) %in% model.brkpts %>%
      sum()
    
    if (tmp == 0) {
      status[i]<- "Missing"
    } else {
      status[i]<- "Present"
    }
  }
  
  miss.ind<- which(status =="Missing")
   
  
  all.brkpts$acc<- accuracy
  all.brkpts[miss.ind, "acc"]<- "Missing"
  
  all.brkpts
}
