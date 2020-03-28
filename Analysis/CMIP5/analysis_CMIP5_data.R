args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ### id for different rcp and CMIP5 model combinations

library(raster)
library(abind)

setwd("/rigel/glab/users/zy2309/")
cmip5_models<-c("CESM1-CAM5","CCSM4","GFDL-ESM2M","NorESM1-ME","NorESM1-M",
                "MIROC-ESM-CHEM","MIROC-ESM","IPSL-CM5B-LR","IPSL-CM5A-MR","IPSL-CM5A-LR",
                "inmcm4","HadGEM2-ES","HadGEM2-CC","GISS-E2-R-CC","GISS-E2-R","GISS-E2-H-CC",
                "GISS-E2-H","GFDL-ESM2G","CanESM2","CESM1-BGC","BNU-ESM","bcc-csm1-1","bcc-csm1-1-m")

rcps<-c("rcp45","rcp85")
vars<-c("tas","mrsos")

projlat<-CRS("+proj=longlat +datum=WGS84 +no_defs")

temp_thr<-10.76810+273.15  ## to deg K
sm_thr<-0.1056108*100  ### convert to kg/m2
## also load the area percentage for 0.5 degree pixels.
landarea_no_arctic<-raster("./DATA/land_data/0.5deg.area.no_arctic.tif")
load("./PROJECT/EOS_SM/CSIF/mean_pheno.RData")
veg_area<-raster(t(mean_eos!=-999)[360:1,])
extent(veg_area)<-c(-180,180,-90,90)


### function for reprojection
reproject2HD<-function(mat,nlat,ext_corr){
  var_year<-raster(t(mat[,nlat:1])) 
  extent(var_year)<-ext_corr
  projection(var_year)<-projlat
  var_year_HD<-projectRaster(from=var_year,to=landarea_no_arctic)
  return(var_year_HD)
}
## global land area
global_land<-133197693
veg_land<-66258632
###### get two calculations (1) get the area change for each year
get_area_for_each_year<-function(model_data){
  ### the 
  prec_lim_area<-rep(NA,95)
  for (i in 1:95){
    ### precipitation limited regions. temp>temp_thr | sm<sm_thr
    prec_lim<-model_data[[1]][[i]]>temp_thr|model_data[[2]][[i]]<sm_thr
    #prec_list[[i]]<-prec_lim
    ## calculate the preci limited regions
    prec_lim_area[i]<-cellStats(prec_lim*landarea_no_arctic/1e6*veg_area,sum,na.rm=T)
  }
  #prec_stack<-stack(prec_list)
  return(prec_lim_area)
}


### get the average of three models
### need to convert the resolution into 0.5 degree.

load_HD_data<-function(model){
  load(paste0("./PROJECT/EOS_SM/CMIP5/",rcp,".",model,".MAT_EOS_SM.RData"))
  MAT_hd_list<-list()
  pre_SM_hd_list<-list()
  nlat<-dim(MAT)[2]
  for (i in 1:95){
    MAT_year_HD<-reproject2HD(MAT[,,i],nlat,ext_corr)
    pre_SM_year_HD<-reproject2HD(pre_SM[,,i],nlat,ext_corr)
    MAT_hd_list[[i]]<-MAT_year_HD
    pre_SM_hd_list[[i]]<-pre_SM_year_HD
  }
  MAT_hd_stack<-stack(MAT_hd_list)
  pre_SM_hd_stack<-stack(pre_SM_hd_list)
  return(list(MAT_hd_stack,pre_SM_hd_stack))
}


rcp<-rcps[uid]
model_data<-list()

### read the data using a loop
for (i in 1:length(cmip5_models)){
  model_data[[i]]<-load_HD_data(cmip5_models[i])
}

ensemble_temp<-(model_data[[1]][[1]]+model_data[[2]][[1]]+model_data[[3]][[1]]+
                  model_data[[4]][[1]]+model_data[[5]][[1]]+model_data[[6]][[1]]+
                  model_data[[7]][[1]]+model_data[[8]][[1]]+model_data[[9]][[1]]+
                  model_data[[10]][[1]]+model_data[[11]][[1]]+model_data[[12]][[1]]+
                  model_data[[13]][[1]]+model_data[[14]][[1]]+model_data[[15]][[1]]+
                  model_data[[16]][[1]]+model_data[[17]][[1]]+model_data[[18]][[1]]+
                  model_data[[19]][[1]]+model_data[[20]][[1]]+model_data[[21]][[1]]+
                  model_data[[22]][[1]]+model_data[[23]][[1]])/23
ensemble_sm<-(model_data[[1]][[2]]+model_data[[2]][[2]]+model_data[[3]][[2]]+
                model_data[[4]][[2]]+model_data[[5]][[2]]+model_data[[6]][[2]]+
                model_data[[7]][[2]]+model_data[[8]][[2]]+model_data[[9]][[2]]+
                model_data[[10]][[2]]+model_data[[11]][[2]]+model_data[[12]][[2]]+
                model_data[[13]][[2]]+model_data[[14]][[2]]+model_data[[15]][[2]]+
                model_data[[16]][[2]]+model_data[[17]][[2]]+model_data[[18]][[2]]+
                model_data[[19]][[2]]+model_data[[20]][[2]]+model_data[[21]][[2]]+
                model_data[[22]][[2]]+model_data[[23]][[2]])/23
ensemble<-list(ensemble_temp,ensemble_sm)

area_all<-array(NA,dim=c(95,24))

area_all[,1]<-get_area_for_each_year(ensemble)
for (i in 1:23){
  area_all[,i+1]<-get_area_for_each_year(model_data[[i]])
}

#save(area_all,file=paste0("./PROJECT/EOS_SM/analysis/CMIP5/",rcp,".area_stat.RData"))
### also calculate the regions for individual models.


### get the (2) maps for current 2006-2015    2031-2040    2061-2070  2091-2100
mat_ave<-list()
sm_ave<-list()
precip_lim<-list()

mat_ave[[1]]<-calc(ensemble_temp[[1:10]],mean,na.rm=T)
mat_ave[[2]]<-calc(ensemble_temp[[26:35]],mean,na.rm=T)
mat_ave[[3]]<-calc(ensemble_temp[[56:65]],mean,na.rm=T)
mat_ave[[4]]<-calc(ensemble_temp[[86:95]],mean,na.rm=T)

sm_ave[[1]]<-calc(ensemble_sm[[1:10]],mean,na.rm=T)
sm_ave[[2]]<-calc(ensemble_sm[[26:35]],mean,na.rm=T)
sm_ave[[3]]<-calc(ensemble_sm[[56:65]],mean,na.rm=T)
sm_ave[[4]]<-calc(ensemble_sm[[86:95]],mean,na.rm=T)

for (i in 1:4){
  precip_lim[[i]]<-(mat_ave[[i]]>temp_thr|sm_ave[[i]]<sm_thr)
}

save(precip_lim,area_all,mat_ave,sm_ave,
     file=paste0("./PROJECT/EOS_SM/analysis/CMIP5/",rcp,".area_stat.RData"))

