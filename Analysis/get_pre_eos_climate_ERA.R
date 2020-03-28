####   calcualte the 4day Ta Rn
args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)

# get the climate variables for the pre eos -------------------------------
library(lubridate)
library(ncdf4)
library(abind)
setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
t2m_f<-list.files("./Climate/ERA/t2m/",full.names = T)
tp_f<-list.files("./Climate/ERA/tp/",full.names = T)
par_f<-list.files("./Climate/ERA/par/",full.names = T)

### eos_data
load("./CSIF/mean_pheno.RData")
#load("./PROJECT/EOS_SM/CSIF/pheno_results.RData")

get_era<-function(files){
  ncin<-nc_open(files[1])
  temp1<-ncvar_get(ncin,varid=ncin$var[[1]])
  ncin<-nc_open(files[2])
  temp2<-ncvar_get(ncin,varid=ncin$var[[1]])
  dat<-abind(temp1,temp2,along=3)
  #dat2<-dat[,360:1,]
  return(dat)
}

get_var_period<-function(vec){
  if(is.na(vec[731])){
    return(rep(NA,4))
  }
  end<-min(730,round(vec[731]+365))
  t0<-mean(vec[(end-15):end])
  t1<-mean(vec[(end-30):end])
  t2<-mean(vec[(end-60):end])
  t3<-mean(vec[(end-90):end])
  return(c(t0,t1,t2,t3))
}

y=uid+2000
#get_pre_eos<-function(y){
## read two years of data from each variable
## one for previous year and one for current year
t2m<-get_era(t2m_f[1:2+y-2001])
tp<-get_era(tp_f[1:2+y-2001])
par<-get_era(par_f[1:2+y-2001])
n<-dim(t2m)[3]
t2m_eos_dat<-abind(t2m[,,(1+leap_year(y-1)):(n-leap_year(y))],mean_eos,along=3)
#temp<-t2m[,,365+240]
#save(temp,file=paste0('./',y,"240.temp.RData"))
#print(dim(t2m_eos_dat))
t2m_eos<-apply(t2m_eos_dat,c(1,2),get_var_period)
tp_eos_dat<-abind(tp[,,(1+leap_year(y-1)):(n-leap_year(y))],mean_eos,along=3)
tp_eos<-apply(tp_eos_dat,c(1,2),get_var_period)
par_eos_dat<-abind(par[,,(1+leap_year(y-1)):(n-leap_year(y))],mean_eos,along=3)
par_eos<-apply(par_eos_dat,c(1,2),get_var_period)

save(t2m_eos,tp_eos,par_eos,file=paste0("./Climate/eos_ERA/",y,".climate_dat.RData"))
#}

#get_pre_eos(uid+2000)

########  calculate the mean T for preEOS
if (F){
  setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
  f<-list.files("./Climate/eos_ERA/",full.names = T)
  pre_eos_max_temp<-array(NA,dim=c(4,720,360))
  pre_eos_mean_temp<-array(NA,dim=c(4,720,360))
  for (i in 1:4){
    m_year_dat<-array(NA,dim=c(length(f),720,360))
    for (j in 1:length(f)){
      load(f[j])
      m_year_dat[j,,]<-t2m_eos[i,,]
    }
    pre_eos_max_temp[i,,]<-apply(m_year_dat,c(2,3),max,na.rm=T)
    pre_eos_mean_temp[i,,]<-apply(m_year_dat,c(2,3),mean,na.rm=T)
  }
  save(pre_eos_max_temp,pre_eos_mean_temp,file="./analysis/pre_eos_mean_max_temp.RData")
}
  
  
  
