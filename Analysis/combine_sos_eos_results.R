### combine sos_eos
sosdate<-array(NA,dim=c(17,720,360))
eosdate<-array(NA,dim=c(17,720,360))
thresh<-array(NA,dim=c(720,360))

files<-list.files("/rigel/glab/users/zy2309/PROJECT/EOS_SM/CSIF/pheno/",full.names = T)

for (i in 1:36){
  load(files[i])
  sosdate[,,(i*10-9):(i*10)]<-sos_all
  eosdate[,,(i*10-9):(i*10)]<-eos_all
  thresh[,(i*10-9):(i*10)]<-threshold[1,,]
}




#image(sosdate[5,,])

# generate the land and crop mask -----------------------------------------
library(ncdf4)
# library(raster)
# export_nc<-function(dat,outfile,varname){
#   latmin<- -90
#   latmax<- 90
#   latd<- 0.5
#   lonmin<- -180
#   lonmax<- 180
#   lond<- 0.5
#
#   lat<- seq(latmin+latd/2,latmax-latd/2,latd)
#   long<-seq(lonmin+lond/2,lonmax-lond/2,lond)
#
#   dimlat<-ncdim_def('latitude','deg',lat)
#   dimlong<-ncdim_def('longitude','deg',long)
#   ncvar<-ncvar_def(varname,'NA',list(dimlong,dimlat),-9999,longname=varname,prec='float',compression=9)
#
#   if (file.exists(outfile)){
#     file.remove(outfile)
#   }
#   ncout<-nc_create(outfile,list(ncvar))
#   ncvar_put(ncout,varid=ncvar,dat)
#   nc_close(ncout)
# }
#
# ncin<-nc_open('/Users/yzhang/Project/SIF_phenology/data/mcd12c1_landcover1_majority.nc')
# lccmg<-ncvar_get(ncin,varid='landcover1')
# ras_lc <- raster(lccmg)
#
# ## aggregate by a factor of 3, with "modal"
# m <- aggregate(ras_lc, fact = 10, fun = modal, na.rm = TRUE)
# m_mat<-getValues(m)
# dim(m_mat)<-c(360,720)
# #m_mat<-apply(t(m_mat),1,rev)
# m_mat<-t(m_mat)
# outfile<-"/Users/yzhang/Project/EOS_SM/Global_mcd12c1_landcover1_majority.nc"
# export_nc(m_mat,outfile,'lc1')

# ncin<-nc_open("/rigel/glab/users/zy2309/PROJECT/EOS_SM/other_data/Global_mcd12c1_landcover1_majority.nc")
# ncin<-nc_open("/Users/yzhang/Project/EOS_SM/Global_mcd12c1_landcover1_majority.nc")
# lc1<-ncvar_get(ncin,varid="lc1")
# lc_type<-c("ENF","EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","CNV")
# crop<-!(lc1==0|lc1==15|lc1==12)
# crop[crop==0]<-NA
# save(crop,file="/Users/yzhang/Project/EOS_SM/land_mask.RData")

load("/rigel/glab/users/zy2309/PROJECT/EOS_SM/other_data/land_mask.RData")
for (i in 1:17){
  sosdate[i,,]<-sosdate[i,,]*crop
  eosdate[i,,]<-eosdate[i,,]*crop
}
save(sosdate,eosdate,thresh,file="/rigel/glab/users/zy2309/PROJECT/EOS_SM/CSIF/pheno_results.RData")

#image(eosdate[4,,],zlim=c(-100,300))


# get the multi year avearge sos and EOS ----------------------------------

mean_sos<-apply(sosdate,c(2,3),median,na.rm=T)
mean_eos<-apply(eosdate,c(2,3),median,na.rm=T)

sd_sos<-apply(sosdate,c(2,3),sd,na.rm=T)
sd_eos<-apply(eosdate,c(2,3),sd,na.rm=T)
save(mean_sos,mean_eos,sd_sos,sd_eos,file="/rigel/glab/users/zy2309/PROJECT/EOS_SM/CSIF/mean_pheno.RData")




