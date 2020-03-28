### autumn SM or sm before the multi year average EOS
library(ncdf4)
library(raster)
library(abind)
library(parallel)
#get the mutltiyear 
setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
load('./CSIF/mean_pheno.RData')
#dim(mean_eos)<-c(720,360)

cl<-makeCluster(getOption("cl.cores",8))
##get_sm_data
sm_f<-list.files("./SMAP/fourday/",full.names = T)
sm<-array(NA,dim=c(360*720,92*3))
for (i in 1:3){
  ncin<-nc_open(sm_f[i])
  dat<-ncvar_get(ncin,varid="SM_AM")
  for (j in 1:92){
    temp<-t(apply(dat[,,j],2,rev))
    sm[,(i-1)*92+j]<-temp
  }
}

dim(sm)<-c(720,360,92,3)
mean_sm<-parApply(cl,sm,c(1,2,3),mean,na.rm=T)
annual_mean_sm<-parApply(cl,mean_sm,c(1,2),mean,na.rm=T)

get_eos_mean<-function(vec){
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
  if(is.na(vec[185])){
    return(c(NA,NA,NA))
  }
  end<-min(184,vec[185]/4+92)
  sm15<-get_mean(vec,end-15/4,end)
  sm30<-get_mean(vec,end-30/4,end)
  sm60<-get_mean(vec,end-60/4,end)
  return(c(sm15,sm30,sm60))
}

sm_eos<-abind(mean_sm,mean_sm,mean_eos,along=3)

mean_pre_eos_sm<-parApply(cl,sm_eos,c(1,2),get_eos_mean)

##### save the mean_pre_eos_sm
save(mean_pre_eos_sm,annual_mean_sm,file="./SMAP/analysis/mean_pre_eos_sm.RData")

