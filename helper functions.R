
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
#--------------------------------------

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
