# args = commandArgs(trailingOnly=F)
# print(args)
# myargs <-sub('-','',args[length(args)])
# print(myargs)
# uid <- as.numeric(myargs)

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
# new method 4/09/2019   #{
# get_threshold<-function(dat){
#   if (sum(is.na(dat))>92){
#     return(rep(NA,3))
#   }
#   dim(dat)<-c(92,17)
#   msc<-apply(dat,1,mean,na.rm=T)
#   min_d<-apply(dat,2,min,na.rm=T)
#   max_d<-apply(dat,2,max,na.rm=T)
#   median_min<-max(median(min_d,na.rm=T),0)
#   median_max<-median(max_d,na.rm=T)
#   threshold<-median_min*0.7+median_max*0.3
#   if(median_min>median_max*0.3|median_max<0.05){
#     return(rep(NA,3))
#   }
#   two_cycle<-c(msc,msc)
#   max_ind<-which.max(msc)[1]
#   one_cycle<-two_cycle[max_ind:(max_ind+92)]
#   min_ind<-which.min(one_cycle)[1]
#   
#   ###### need to work on find the min period
#   below_thr<-which(one_cycle<threshold)
#   if (length(below_thr)==0){
#     return(c(threshold,max_ind,min_ind))
#   }
#   sep_list<-tapply(below_thr, cumsum(c(TRUE, diff(below_thr) != 1)), identity)
#   median_of_min<-as.vector(unlist(lapply(sep_list,function(x)floor(mean(x[min_ind %in% x])))))
#   median_min<-(median_of_min[!is.na(median_of_min)]+max_ind-1)
#   if(max_ind>80){
#     median_min<-median_min-92
#     max_ind<-max_ind-92
#   }
#   return(c(threshold,max_ind,median_min))
# }
#}
###not use 
get_threshold<-function(dat){
  dim(dat)<-c(92,17,720,10)
  min_year<-apply(dat,c(2,3,4),min)
  max_year<-apply(dat,c(2,3,4),max)
  median_min<-apply(min_year,c(2,3),median)
  median_min[median_min<0]<-0
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



### new method 4/09/19
# retrieve_pixel<-function(pheno_vec){
#   if (sum(is.na(pheno_vec))>92){
#     return(c(NA,NA))
#   }
#   min_start<-pheno_vec[279]
#   max_ind<-pheno_vec[278]
#   thre<-pheno_vec[277]
#   new_max<-max_ind+92-(min_start-10)
#   if(is.na(thre)){
#     return(c(NA,NA))
#   }
#   ##select a full year cycle for current year with 10 more obs before and 10 more after
#   one_cycle<-pheno_vec[(min_start-10):(min_start+101)]
#   peak_d<-which.max(one_cycle)
#   if(length(peak_d)>1){
#     peak_d<-median(peak_d)
#   }
#   if(one_cycle[peak_d]<thre){
#     return(c(NA,NA))
#   }
#   ### first guess the multi year values
#   increase<-diff(one_cycle[1:new_max])>0
#   increase<-c(increase[1],increase)
#   decrease<-diff(one_cycle[new_max:112])<0
#   decrease<-c(decrease[1],decrease)
#   gt1<-one_cycle[1:new_max]>thre
#   sos<-which(c(0,diff(gt1))&increase==1)
#   gt2<-one_cycle[new_max:112]>thre
#   eos<-which(c(0,diff(gt2))&decrease==1)+(new_max-1)
#   ## second guess
#   if (length(sos)==0){
#     increase<-diff(one_cycle[1:peak_d])>0
#     increase<-c(increase[1],increase)
#     gt1<-one_cycle[1:peak_d]>thre
#     sos<-which(c(0,diff(gt1))&increase==1)
#   }
#   if (length(eos)==0){
#     decrease<-diff(one_cycle[peak_d:112])<0
#     decrease<-c(decrease[1],decrease)
#     gt2<-one_cycle[peak_d:112]>thre
#     eos<-which(c(0,diff(gt2))&decrease==1)+(new_max-1)
#   }
#   sos<-min(sos)
#   eos<-max(eos)
#   sos_acc<-get_acc(sos-1,one_cycle[sos-1],one_cycle[sos],thre)+min_start-10-92 -1
#   eos_acc<-get_acc(eos,one_cycle[eos],one_cycle[eos+1],thre)+min_start-10-92-1
#   return(c(sos_acc,eos_acc))
# }

### phenology pixel   not use
retrieve_pixel<-function(vec){
  ### find the peak the new length of vec is 92+20 obs
  if(is.na(vec[113])){
    return(c(NA,NA))
  }
  peak_d<-which.max(vec[11:102])+10
  ## the peak should between the 47:136
  if(length(peak_d)>1){
    peak_d<-median(peak_d)
  }
  if(vec[peak_d]<vec[113]){
    return(c(NA,NA))
  }
  # get the sos between 1:peak_d
  ### new method on 4/10/19
  increase<-diff(vec[1:peak_d])>0
  increase<-c(increase[1],increase)
  decrease<-diff(vec[peak_d:112])<0
  decrease<-c(decrease[1],decrease)
  gt1<-vec[1:peak_d]>vec[113]
  sos<-which(c(0,diff(gt1))&increase==1)
  gt2<-vec[peak_d:112]>vec[113]
  eos<-which(c(0,diff(gt2))&decrease==1)+(peak_d-1)
  sos<-min(sos)
  eos<-max(eos)
  sos_acc<-get_acc(sos-1,vec[sos-1],vec[sos],vec[113])
  eos_acc<-get_acc(eos,vec[eos],vec[eos+1],vec[113])
  # sos<-max(which(vec[1:peak_d]<vec[185]))
  # get the eos between peak:184
  # eos<-min(which(vec[peak_d:184]<vec[185]))+peak_d-1
  # sos_acc<-get_acc(sos,vec[sos],vec[sos+1],vec[185])
  # eos_acc<-get_acc(eos-1,vec[eos-1],vec[eos],vec[185])
  return(c(sos_acc,eos_acc))
}

##
adjust_north<-function(vec){
  vec[is.na(vec)]<- -1000
  res<-vec
  ## 83:194   break point at 92 i.e.9, start doy is 82*4
  res[vec<9]<-vec[vec<9]*4+82*4-2-365
  res[vec>=9&vec<=11]<-vec[vec>=9&vec<=11]*2.5+339.5-365
  res[vec>11]<-vec[vec>11]*4-5+82*4-365
  return(res)
}

##
adjust_south<-function(vec){
  vec[is.na(vec)]<- -1000
  res<-vec
  ## 37:148   break point at 92 i.e. 55, start doy is 36*4
  res[vec<55]<-vec[vec<55]*4-2+36*4-365
  res[vec>=55&vec<=57]<-vec[vec>=55&vec<=57]*2.5+80.5+36*4-365
  res[vec>57]<-vec[vec>57]*4-5+36*4-365
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
      pheno_dat<-abind(year_dat[83:194,,],threshold,along=1)
      sos_eos_temp<-apply(pheno_dat,c(2,3),retrieve_pixel)
      ### for the current year, if negative values, it represent the sos starts previous year.

      sos_eos<-apply(sos_eos_temp,c(2,3),adjust_north)

      sos_all[y,,]<-sos_eos[1,,]
      eos_all[y,,]<-sos_eos[2,,]
    }else{
      ### for south use previous+current
      pheno_dat<-abind(year_dat[37:148,,],threshold,along=1)
      sos_eos_temp<-apply(pheno_dat,c(2,3),retrieve_pixel)
      ### two breaks if(x<91) y=x*4-367

      sos_eos<-apply(sos_eos_temp,c(2,3),adjust_south)
      # for the current year, if negative values, it represent the sos starts previous year.
      sos_all[y,,]<-sos_eos[1,,]
      eos_all[y,,]<-sos_eos[2,,]
    }
    
    # ## modified 4/09/19
    # pheno_dat<-abind(year_dat,threshold,along=1)
    # sos_eos_temp<-apply(pheno_dat,c(2,3),retrieve_pixel)
    # 
    # sos_all[y,,]<-sos_eos_temp[1,,]*4
    # eos_all[y,,]<-sos_eos_temp[2,,]*4
  }
  sos_all[sos_all< -1000]<-NA
  eos_all[eos_all< -1000]<-NA
  save(sos_all,eos_all,threshold,file = paste0("./PROJECT/EOS_SM/CSIF/pheno/",formatC(uid,width=2,flag="0"),
                                               ".pheno_daily_csif.RData"))
}


#### save the smoothed SIF

# smoothed_sif<-load_sif(uid)
# save(smoothed_sif,file=paste0("./PROJECT/EOS_SM/CSIF/smoothed_sif/",formatC(uid,width=2,flag="0"),
#                               ".smoothed_clear_daily_csif.RData"))


for (uid  in 1:36){
load(paste0("./PROJECT/EOS_SM/CSIF/smoothed_sif/",formatC(uid,width=2,flag="0"),
                                           ".smoothed_clear_daily_csif.RData"))

#thresh<-apply(smoothed_sif,c(2,3),get_threshold)
  thresh<-get_threshold(smoothed_sif)
# for (i in 1:720){
#   for (j in 1:10){
#     get_threshold(smoothed_sif[,i,j])
#   }
# }
## retrieve phenology
### read in three years of data first

retrieve_pheno(smoothed_sif,thresh,uid)
}

