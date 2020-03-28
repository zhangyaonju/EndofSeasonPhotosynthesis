### get pre eos climate data from remote sensing datasets
####   calcualte the 4day Ta Rn
args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ### year starts from 2003 to 2016 total of 14 years

library(ncdf4)
library(abind)
library(lubridate)
setwd("/rigel/glab/users/zy2309/")


## read two years of data 
get_RS_data<-function(files,k){
  ncin<-nc_open(files[1])
  dat<-ncvar_get(ncin,varid=ncin$var[[k]])
  al<-which(dim(dat)%%180!=0)
  if (length(al)==0){
    al=3
  }
  for (i in 2:length(files)){
    ncin<-nc_open(files[i])
    temp1<-ncvar_get(ncin,varid=ncin$var[[k]])
    dat<-abind(dat,temp1,along=al)
  }
  return(dat)
}


load("./PROJECT/EOS_SM/CSIF/mean_pheno.RData")
## read the PAR data
par_files<-list.files("./DATA/BESS_daily_HD/",full.names = T)
lst_files<-list.files("./DATA/MYD11C2_HD/",full.names = T)
###prec_files<-list.files("./DATA/MSWEP/",full.names = T)
prec_files<-list.files("./DATA/GPCP_daily/",full.names = T)
## using the GPCP dataset instead of 
year<-uid+2002

### the temp is 8day, PAR and prec are daily

par_f<-par_files[substr(basename(par_files),1,4)==year|substr(basename(par_files),1,4)==year-1]
par<-get_RS_data(par_f,1)
#prec_f<-prec_files[substr(basename(prec_files),1,4)==year|substr(basename(prec_files),1,4)==year-1]
## for the gpcp data, the long range is 0-360 deg and spatial resolutio is 1 deg
prec_f<-prec_files[substr(basename(prec_files),24,27)==year|substr(basename(prec_files),24,27)==year-1]
prec<-get_RS_data(prec_f,6) ### for the GPCP data, the dimension are also vars and prec is the 6th variable
prec_temp<-abind(prec[181:360,180:1,],prec[1:180,180:1,],along=1)


lst_f<-lst_files[substr(basename(lst_files),1,4)==year|substr(basename(lst_files),1,4)==year-1]
lst_day<-get_RS_data(lst_f,1)
lst_night<-get_RS_data(lst_f,2)

#### get the preseason data

get_var_period<-function(vec){
  n<-length(vec)
  if(is.na(vec[n])){
    return(rep(NA,4))
  }
  end<-min(730,vec[n]+365)
  t0<-mean(vec[(end-15):end],na.rm=T)
  t1<-mean(vec[(end-30):end],na.rm=T)
  t2<-mean(vec[(end-60):end],na.rm=T)
  t3<-mean(vec[(end-90):end],na.rm=T)
  return(c(t0,t1,t2,t3))
}

### get mean for the 8day 
get_mean<-function(dat,s1,e1){
  s1<-max(1,s1)
  sint<-ceiling(s1-0.5)
  eint<-floor(e1-0.5)
  sum_dat<-sum(dat[(sint+2):eint])+
    ((0.5+sint-s1)*dat[sint]+(1.5-sint+s1)*dat[sint+1])*(0.5+sint-s1)/2+
    0.5*dat[sint+1]+
    ((e1-0.5-eint)*dat[eint+2]+(2.5+eint-e1)*dat[eint+1])*(e1-0.5-eint)/2+
    0.5*dat[eint+1]
  return(sum_dat/(e1-s1))
}

get_var_period_8day<-function(vec){
  if(is.na(vec[93])){
    return(rep(NA,4))
  }
  end<-min(92,vec[93]/8+46)
  t0<-get_mean(vec,end-15/8,end)
  t1<-get_mean(vec,end-30/8,end)
  t2<-get_mean(vec,end-60/8,end)
  t3<-get_mean(vec,end-90/8,end)
  return(c(t0,t1,t2,t3))
}

###for LST
al<-which(dim(lst_day)%%360!=0)
adat<-abind(lst_day,mean_eos,along=al)
lst_day_eos<-apply(adat,setdiff(c(1,2,3),al),get_var_period_8day)
adat<-abind(lst_night,mean_eos,along=al)
lst_night_eos<-apply(adat,setdiff(c(1,2,3),al),get_var_period_8day)

### for prec
#$%^ this is for the MSWEP data 
# al<-which(dim(prec)%%360!=0)
# n<-dim(prec)[al]
# adat<-abind(prec[,360:1,(1+leap_year(year-1)):(n-leap_year(year))],mean_eos,along=al)
# tp_eos<-apply(adat,setdiff(c(1,2,3),al),get_var_period)

#$%^ this is for the GPCP v1.3 data
al<-3
n<-dim(prec)[al]
precd<-array(NA,dim=c(720,360,n))
for (i in 0:3){
  precd[1:360*2-i%%2,1:180*2-floor(i/2),]<-prec_temp
}
adat<-abind(precd[,,(1+leap_year(year-1)):(n-leap_year(year))],mean_eos,along=al)
tp_eos<-apply(adat,setdiff(c(1,2,3),al),get_var_period)

### for par
al<-which(dim(par)%%360!=0)
n<-dim(par)[al]
adat<-abind(par[,,(1+leap_year(year-1)):(n-leap_year(year))],mean_eos,along=al)
par_eos<-apply(adat,setdiff(c(1,2,3),al),get_var_period)

save(lst_day_eos,lst_night_eos,tp_eos,par_eos,
     file=paste0("./PROJECT/EOS_SM/Climate/eos_RS2/",year,".climate_dat.RData"))


