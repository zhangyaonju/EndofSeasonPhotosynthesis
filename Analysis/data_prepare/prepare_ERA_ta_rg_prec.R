####   calcualte the 4day Ta Rn
args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)

library(ncdf4)
library(abind)
library(lubridate)
setwd("/rigel/glab/users/zy2309/")


exportnc<-function(dat,var,unit,outfile,year){
  dimx<-ncdim_def(name = "longitude", unit="",vals=seq(-179.75,179.75,0.5))
  dimy<-ncdim_def(name = "latitude", unit="",vals=seq(-89.75,89.75,0.5))
  
  nday<-leap_year(year)
  dimt<-ncdim_def(name = "day",units = "",vals = rep(year,nday+365)*1000+1:(365+nday))
  var_tre<-ncvar_def(name = var,units = unit,dim = list(dimx,dimy,dimt),
                     missval = -9999,longname = var,compression=9)
  dat[is.na(dat)]<- -9999
  ncout<-nc_create(filename = outfile,vars = list(var_tre))
  ncvar_put(nc = ncout,varid = var_tre,vals = dat)
  nc_close(ncout)
}

era_f1<-list.files("./DATA/ERA/",pattern="^2000",full.names = T)
era_f2<-list.files("./DATA/ERA/",pattern="^2015",full.names = T)
vars<-c("par","tp","t2m")
ydays1<-c(0.5, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365, 365, 366, 365, 365)*2
ydays2<-c(0.5, 365, 366, 365)*2
#get_climate_date<-function(file1,file2,var){
file1<-era_f1[uid]
file2<-era_f2[uid]
var<-vars[uid]
  for (y in 2000:2017){
    if (y==2000){
      ncin<-nc_open(file1)
      var_dat<-ncvar_get(ncin,varid=ncin$var[[1]])
    }else if(y==2015){
      ncin<-nc_open(file2)
      var_dat<-ncvar_get(ncin,varid=ncin$var[[1]])
    }
    if (y< 2015){
      yeardat<-var_dat[,,sum(ydays1[1:(y-1999)]):(sum(ydays1[1:(y-1998)])-1)]
    }else{
      yeardat<-var_dat[,,sum(ydays2[1:(y-2014)]):(sum(ydays2[1:(y-2013)])-1)]
    }
    n<-dim(yeardat)[3]
    #### convert half day to day
    dim(yeardat)<-c(720,361,2,n/2)
    yearday<-apply(yeardat,c(1,2,4),mean)
    
    firsthalf<-(yearday[1:360,1:360,]+yearday[1:360,2:361,]+
                  yearday[2:361,1:360,]+yearday[2:361,2:361,])/4
    secondhalf<-(yearday[361:720,1:360,]+yearday[361:720,2:361,]+
                   yearday[c(362:720,1),1:360,]+yearday[c(362:720,1),2:361,])/4
    var_daily_dat<-abind(secondhalf[,360:1,],firsthalf[,360:1,],along=1)
    if (var=="par"){
      var_daily_dat<-var_daily_dat/12/3600  ## convert J/m2 to w/m2
      unit<-"w/m2"
    }else if (var=="tp"){
      var_daily_dat<-var_daily_dat*2*1000  ## convert m/12hour to mm/day
      unit<-"mm/day"
    }else if (var=="t2m"){
      var_daily_dat<-var_daily_dat-273.16  ## k to deg C
      unit<-"deg C"
    }
    exportnc(var_daily_dat,var,unit,outfile=paste0("./PROJECT/EOS_SM/Climate/ERA/",
                                                   var,"/",y,".daily_",var,".nc"),year=y)
  }
#}

#get_climate_date(era_f1[uid],era_f2[uid],vars[uid])


