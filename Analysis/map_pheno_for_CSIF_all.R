args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)
uid <- as.numeric(myargs)

# get phenology for globe -------------------------------------------------
library(ncdf4)
library(abind)
setwd("/rigel/glab/users/zy2309/")
csif_f<-list.files("./DATA/bg_csif/clear_inst_SIF_4day_HD/",full.names = T)
# load the CSIF data ------------------------------------------------------
load_sif<-function(ilat){
  # 
  csif<-array(NA,dim=c(720,10,92*17))
  for (i in 1:length(csif_f)){
    ncin<-nc_open(csif_f[i])
    temp_dat<-ncvar_get(ncin,varid="clear_daily_sif")
    csif[,,(i*92-91):(i*92)]<-temp_dat[,(ilat*10-9):(ilat*10),]
  }
  ## set NA to 0
  csif[is.na(csif)]<-0
  # smooth the csif year by year from 2001 to 2017
  smoothed_csif<-array(NA,dim=c(92*17,720,10))
  for (y in 1:17){
    if (y!=1 & y!=17){
      year_dat<-csif[,,(y*92-183):(y*92+92)]
    }else if (y==1){
      year_dat<-abind(csif[,,1:92],csif[,,1:184],along=3)
    }else if (y==17){
      year_dat<-abind(csif[,,1381:1564],csif[,,1473:1564],along=3)
    }
    smoothed_csif[(y*92-91):(y*92),,]<-apply(year_dat,c(1,2),smooth_sif)
  }
  return(smoothed_csif)
}

# 
get_threshold<-function(dat){

  dim(dat)<-c(92,17,720,10)
  min_year<-apply(dat,c(2,3,4),min)
  max_year<-apply(dat,c(2,3,4),max)
  median_min<-apply(min_year,c(2,3),median)
  median_max<-apply(max_year,c(2,3),median)
  
  threshold<-median_min*0.7+median_max*0.3
  ### if variability is too small, it is considered as evergreen and threshold is set to zero
  threshold[(median_min>median_max*0.3)|(median_max<0.05)]<-NA
  return(threshold)
}


smooth_sif<-function(vec){
  n<-length(vec)
  if (sum(vec==0)==n){
    return(rep(NA,92))
  }
  df<-round(n/9)
  sm_vec<-smooth.spline(x = rep(c((1:91)*4-2,364),3)+rep((0:2)*365,each=92),y=vec,df = df)
  return(sm_vec$y[93:184])
}

get_acc<-function(x1,y1,y2,thr){
  if(is.na(x1)|is.infinite(x1)){
    return(NA)
  }
  x<-(thr-y1)/(y2-y1)
  return(x1+x)
}


### phenology pixel
retrieve_pixel<-function(vec){
  ### find the peak
  if(is.na(vec[185])){
    return(c(NA,NA))
  }
  peak_d<-which.max(vec[45:130])+44
  ## the peak should between the 47:136
  if(length(peak_d)>1){
    peak_d<-median(peak_d)
  }
  if(vec[peak_d]<vec[185]){
    return(c(NA,NA))
  }
  # get the sos between 1:peak_d
  sos<-max(which(vec[1:peak_d]<vec[185]))
  # get the eos between peak:184
  eos<-min(which(vec[peak_d:184]<vec[185]))+peak_d-1
  sos_acc<-get_acc(sos,vec[sos],vec[sos+1],vec[185])
  eos_acc<-get_acc(eos-1,vec[eos-1],vec[eos],vec[185])
  return(c(sos_acc,eos_acc))
}

##
adjust_north<-function(vec){
  vec[is.na(vec)]<- -1000
  res<-vec
  res[vec<45]<-vec[vec<45]*4-183
  res[vec>=45&vec<=47]<-vec[vec>45&vec<=47]*2+88-181
  res[vec>47]<-vec[vec>47]*4-186
  return(res)
}

##
adjust_south<-function(vec){
  vec[is.na(vec)]<- -1000
  res<-vec
  res[vec<91]<-vec[vec<91]*4-367
  res[vec>=91&vec<=93]<-vec[vec>91&vec<=93]*2+180-365
  res[vec>93]<-vec[vec>93]*4-370
  return(res)
}

### retrieve_phenology
retrieve_pheno<-function(smooth_dat,threshold,uid){
  sos_all<-array(NA,dim=c(17,720,10))
  eos_all<-array(NA,dim=c(17,720,10))
  dim(threshold)<-c(1,720,10)
  for (y in 1:17){
    if (y!=1 & y!=17){
      year_dat<-smooth_dat[(y*92-183):(y*92+92),,]
    }else if (y==1){
      year_dat<-abind(smooth_dat[1:92,,],smooth_dat[1:184,,],along=1)
    }else if (y==17){
      year_dat<-abind(smooth_dat[1381:1564,,],smooth_dat[1473:1564,,],along=1)
    }
    
    if (uid>18){
      ## for north use half+full+half
      pheno_dat<-abind(year_dat[47:230,,],threshold,along=1)
      sos_eos_temp<-apply(pheno_dat,c(2,3),retrieve_pixel)
      ### for the current year, if negative values, it represent the sos starts previous year.

      sos_eos<-apply(sos_eos_temp,c(2,3),adjust_north)
      
      sos_all[y,,]<-sos_eos[1,,]
      eos_all[y,,]<-sos_eos[2,,]
    }else{
      ### for south use previous+current
      pheno_dat<-abind(year_dat[1:184,,],threshold,along=1)
      sos_eos_temp<-apply(pheno_dat,c(2,3),retrieve_pixel)
      ### two breaks if(x<91) y=x*4-367
      
      sos_eos<-apply(sos_eos_temp,c(2,3),adjust_south)
      # for the current year, if negative values, it represent the sos starts previous year.
      sos_all[y,,]<-sos_eos[1,,]
      eos_all[y,,]<-sos_eos[2,,]
    }
  }
  sos_all[sos_all< -1000]<-NA
  eos_all[eos_all< -1000]<-NA
  save(sos_all,eos_all,threshold,file = paste0("./PROJECT/EOS_SM/CSIF/pheno/",formatC(uid,width=2,flag="0"),
                                               ".pheno_daily_csif.RData"))
}


#### save the smoothed SIF

smoothed_sif<-load_sif(uid)
save(smoothed_sif,file=paste0("./PROJECT/EOS_SM/CSIF/smoothed_sif/",formatC(uid,width=2,flag="0"),
                              ".smoothed_clear_daily_csif.RData"))

thresh<-get_threshold(smoothed_sif)

## retrieve phenology
### read in three years of data first
retrieve_pheno(smoothed_sif,thresh,uid)


