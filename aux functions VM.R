dVM=function(x,mu1,kappa1){
  dvonmises(x,mu=mu1,kappa=kappa1)
}

#pvonmises gives some weird results. I will numerically integrate myself
integrated.dVM=function(lo,hi,mu1,kappa1){
  integrate(dVM,lower=lo,upper=hi,mu1=mu1,kappa1=kappa1)  
}

vm.function=function(param){  #for calculating probability masses
  mu1=circular(param[1],type='angles',units='radians')
  kappa1=param[2]
  
  bins.estim=rep(NA,(length(angle.bin.lims) - 1))
  
  for (i in 2:length(angle.bin.lims)) {
    pvonmises1 = integrated.dVM(lo=angle.bin.lims[i-1],
                                hi=angle.bin.lims[i],
                                mu1=mu1,
                                kappa1=kappa1)
    bins.estim[i-1] = pvonmises1$value
  }
  
  unlist(bins.estim)
}

vm.function.optim=function(param, probs){  #for optimization
  bins.estim=vm.function(param)
  mean(abs(probs-bins.estim))  #Mean Absolute Error
}


