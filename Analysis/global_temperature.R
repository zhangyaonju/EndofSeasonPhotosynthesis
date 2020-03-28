### calculate mean annual temperature
### calculate the culmulative tempearture >0 and >5 deg

args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ## uid for the four different lags

setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
library(ncdf4)
library(parallel)
era_t_f<-list.files("./Climate/ERA/t2m/",full.names = T)
year=1999+uid

cl<-makeCluster(getOption("cl.cores",8))

ncin<-nc_open(era_t_f[uid])
temp<-ncvar_get(ncin,varid=ncin$var[[1]])

mean_t<-parApply(cl,temp,c(1,2),mean,na.rm=T)

### temp>0
temp0<-temp
temp0[temp0<0]<-NA
sum_t0<-parApply(cl,temp0,c(1,2),sum,na.rm=T)

temp5<-temp
temp5[temp5<5]<-NA

sum_t5<-parApply(cl,temp5,c(1,2),sum,na.rm=T)

save(mean_t,sum_t0,sum_t5,file=paste0("./Climate/mean_annual_T/",year,".temperature.RData"))



