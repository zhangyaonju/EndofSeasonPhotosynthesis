### retrieve phenology
splinefit <- function(y, tin , w, freedom=7, ...){
  # xpred.out <- spline(t, y, xout = tout)$y %>% zoo(., tout)
  n <- length(y)
  #tin <- t-floor(t)
  ylu<-range(y)
  if (length(y)==92){
    freedom=9
  }else if(length(y)>100){
    freedom=11
  }
  tout <- tin
  fit <- smooth.spline(tin, y, w, df = freedom)
  #z <- predict(fit)$y
  xpred.out <- predict(fit, tout)$y
  
  return(xpred.out)
}

get_MR<-function(y){
  ### get multi year average seasonal cycle
  n_obs<-length(y[[1]])
  if (n_obs>94){
    n_obs<-365
  }
  obs_array<-array(NA,dim=c(n_obs,length(y)))
  for (i in 1:length(y)){
    obs_array[,i]<-y[[i]][1:n_obs]
  }
  m_season<-apply(obs_array,1,mean,na.rm=T)
  change_rate<-diff(m_season)
  #### get threshold 30%
  v_range<-range(m_season)
  v_range[1]<-max(0,v_range[1])
  if (v_range[2]<=0.05){
    return(rep(NA,3))
  }
  thr<-v_range[1]+0.3*(diff(v_range))
  mr<-c(max(m_season,na.rm=T),m_season[which.min(change_rate)[1]+1],thr)
  return(mr)
}

remove_not_consective<-function(vec,n_min){
  # get the summary
  tf_summ<-rle(vec)
  toremove<-which(tf_summ$values==T&tf_summ$lengths<n_min)
  if (length(toremove)==0){
    return(vec)
  }
  end_id<-cumsum(tf_summ$lengths)
  start_id<-c(1,end_id[1:length(end_id)-1]+1)
  for (i in 1:length(toremove)){
    vec[start_id[toremove[i]]:end_id[toremove[i]]]<-F
  }
  return(vec)
}

thresholdextract<-function(fitval,tin,threshval){
  ## peak
  maxdate<-which.max(fitval)[1]
  if(any(max(fitval)<threshval)){
    return(tin[maxdate])
  }
  # green up
  greenup<-diff(fitval)>=0
  greenup<-c(greenup[1],greenup)
  browndown<-diff(fitval)<=0
  browndown<-c(browndown[1],browndown)
  #### revised for fluxdata
  #### for sos , the one near peak is used, for EOS, the one later is used.
  gt_thr<-fitval>=threshval
  ### the consective period is 30 days, for CSIF, it is 8 obs.
  if (length(tin)>95){
    cons<-30
  }else{
    cons<-8
  }
  good_gs<-remove_not_consective(gt_thr,cons) ##
  gs_all<-which(browndown&good_gs)
  if(any(gs_all>maxdate)){
    eos_init<-max(gs_all[gs_all>=maxdate])
    eos<-(get_accurate(fitval[eos_init:(eos_init+1)],threshval)-0.5)/length(tin)+tin[eos_init]
  }else(
    eos<-NA
  )
  return(eos)
}

get_accurate<-function(nearby,thr){
  return((nearby[1]-thr)/(nearby[1]-nearby[2]))
}

pheno_wSpline<-function(y, t, w, mr){
  fitted_y<-splinefit(y,t,w)   #fitted_y<-splinefit(y2011$y,y2011$t,y2011$w)
  ### extract values from threshold
  ### because the 
  phen<-thresholdextract(fitted_y,t,mr[3])
  return(phen)
}


