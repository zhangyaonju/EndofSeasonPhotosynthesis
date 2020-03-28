args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ### id for different rcp and CMIP5 model combinations

# analyze the cmip5 dataset -----------------------------------------------

library(ncdf4)
library(raster)
library(abind)
library(lubridate)

setwd("/rigel/glab/users/zy2309/")
cmip5_models<-c("CESM1-CAM5","CCSM4","GFDL-ESM2M","NorESM1-ME","NorESM1-M",
                "MIROC-ESM-CHEM","MIROC-ESM","IPSL-CM5B-LR","IPSL-CM5A-MR","IPSL-CM5A-LR",
                "inmcm4","HadGEM2-ES","HadGEM2-CC","GISS-E2-R-CC","GISS-E2-R","GISS-E2-H-CC",
                "GISS-E2-H","GFDL-ESM2G","CanESM2","CESM1-BGC","BNU-ESM","bcc-csm1-1","bcc-csm1-1-m")

#### "MPI-ESM-MR","MPI-ESM-LR", does not have mrsos
rcps<-c("rcp45","rcp85")
vars<-c("tas","mrsos")

projlat<-CRS("+proj=longlat +datum=WGS84 +no_defs")

# # check the files for each model ------------------------------------------
# for (i in 1:length(cmip5_models)){
#   stat<-c()
#   for (r in 1:2){
#     for (v in 1:2){
#       stat[r*2-2+v]<-length(list.files(paste0("./DATA/CMIP5/",cmip5_models[i]),
#                                        pattern=paste0(vars[v],"_.*",cmip5_models[i],"_",rcps[r],".*")))
#     }
#   }
#   print(paste(cmip5_models[i],paste(stat,collapse = "  ")))
# }


# process functions --------------------------------------------------------

### get the year from the file name
st_end_from_name<-function(file){
  st_end_string<-strsplit(basename(file),"_")[[1]][6]
  st_ym<-interval(ymd(paste0(substr(st_end_string,1,6),"01")),ymd("20060101"))
  st_diff<-st_ym %/% months(1)
  end_ym<-interval(ymd("21001201"),ymd(paste0(substr(st_end_string,8,13),"01")))
  end_diff<-end_ym %/% months(1)
  num_mon<-interval(ymd(paste0(substr(st_end_string,1,6),"01")),
                    ymd(paste0(substr(st_end_string,8,13),"01")))
  st_end<-num_mon %/% months(1) +1
  ######### get the n range
  return((1+max(st_diff,0)):(st_end-max(end_diff,0)))
}



read_data<-function(files,var){
  model_name<-strsplit(basename(files[1]),"_")[[1]][3]
  ncin<-nc_open(files[1])
  st_end<-st_end_from_name(files[1])
  var_dat<-ncvar_get(ncin,varid=ncin$var[[var]])[,,st_end]
  latdat<-ncvar_get(ncin,varid=ncin$var[["lat_bnds"]])
  londat<-ncvar_get(ncin,varid=ncin$var[["lon_bnds"]])
  nlat<-dim(latdat)[2]
  nlon<-dim(londat)[2]
  
  if (length(files)>1){
    for (i in 2:length(files)){
      ncin<-nc_open(files[i])
      st_end<-st_end_from_name(files[i])
      temp<-ncvar_get(ncin,varid=ncin$var[[var]])[,,st_end]
      var_dat<-abind(var_dat,temp,along=3)
    }
  }
  # check if the data has 95 years of months
  if(dim(var_dat)[3]==95*12){
    print(paste(model_name,"good range"))
  }
  
  #### rearrange the data and extent 
  if (latdat[2,1]<0){
    ext<-c(londat[2,nlon/2]-360,londat[2,nlon/2],latdat[2,1],latdat[1,nlat])
    latdat<-latdat[,2:(nlat-1)]
    var_dat<-abind(var_dat[(nlon/2+1):nlon,2:(nlat-1),],var_dat[1:(nlon/2),2:(nlat-1),],along=1)
  }else{
    ext<-c(londat[2,nlon/2]-360,londat[2,nlon/2],latdat[1,nlat],latdat[2,1])
    latdat<-latdat[,(nlat-1):2]
    var_dat<-abind(var_dat[(nlon/2+1):nlon,(nlat-1):2,],var_dat[1:(nlon/2),(nlat-1):2,],along=1)
  }
  return(list(var_dat,ext))
}
## get preseason soil moisture  ### input is two years of sm data, and the average EOS date
load("./PROJECT/EOS_SM/CSIF/mean_pheno.RData")
### resample it to the resolution of the CMIP5 model

get_mat_EOS_sm<-function(model,rcp){
  ### two functions for calculating the EOS sm
  preseason_sm<-function(vec){
    ### first 24 is the monthly mrsos, the 25th is the eos date in yday
    if(is.na(vec[25])){
      return(NA)
    }
    end<-min(24,vec[25]/30.42+12)
    #t0<-get_mean(vec,end-15/30.42,end)
    t1<-get_mean(vec,end-1,end)
    #t2<-get_mean(vec,end-60/8,end)
    #t3<-get_mean(vec,end-90/8,end)
    return(t1)
  }
  
  get_mean<-function(dat,s1,e1){
    s1<-max(1,s1)
    sint<-ceiling(s1-0.5)
    eint<-floor(e1-0.5)
    sum_dat<-((0.5+sint-s1)*dat[sint]+(1.5-sint+s1)*dat[sint+1])*
      (0.5+sint-s1)/2+
      ((e1-0.5-eint)*dat[eint+2]+(2.5+eint-e1)*dat[eint+1])*
      (e1-0.5-eint)/2
    return(sum_dat/(e1-s1))
  }
  
  tas_f<-list.files(paste0("./DATA/CMIP5/",model),pattern=paste0("tas_.*",rcp,".*"),full.names = T)
  mrsos_f<-list.files(paste0("./DATA/CMIP5/",model),pattern=paste0("mrsos_.*",rcp,".*"),full.names = T)
  #get the data using a loop
  tas<-read_data(tas_f,"tas")
  tas_dat<-tas[[1]]
  ext<-dim(tas_dat)[1:2]
  dim(tas_dat)<-c(ext,12,95)
  MAT<-apply(tas_dat,c(1,2,4),mean,na.rm=T)
  
  sm<-read_data(mrsos_f,"mrsos")
  sm_dat<-sm[[1]]
  nlat<-dim(sm_dat)[2]
  # process the mean EOS date to the same resolution
  fm_mask<-raster(t(sm_dat[,,1]))
  
  ext_corr<-sm[[2]]
  extent(fm_mask)<-ext_corr
  
  projection(fm_mask)<-projlat
  Rmean_eos<-raster(t(mean_eos)[360:1,])
  extent(Rmean_eos)<-c(-180,180,-90,90)
  projection(Rmean_eos)<-projlat
  meaneos_proj<-projectRaster(from = Rmean_eos,to=fm_mask)
  ### convert back to matrix for calculation
  mat_eos<-t(as.matrix(meaneos_proj))[,nlat:1]
  
  pre_SM<-array(NA,dim=c(ext,95))
  for (i in 2006:2100){
    ## two years of data
    if(i==2006){
      cdat<-abind(sm_dat[,,((i-2006)*12+1):((i-2005)*12)],sm_dat[,,((i-2006)*12+1):((i-2005)*12)],
                  mat_eos,along=3)
    }else{
      cdat<-abind(sm_dat[,,((i-2007)*12+1):((i-2005)*12)],mat_eos,along=3)
    }
    pre_SM[,,i-2005]<-apply(cdat,c(1,2),preseason_sm)
  }
  outfile<-paste0("./PROJECT/EOS_SM/CMIP5/",rcp,".",model,".MAT_EOS_SM.RData")
  save(MAT,pre_SM,ext_corr,file=outfile)
}

rcp<-rcps[uid%%2+1]
model<-cmip5_models[ceiling(uid/2)]

get_mat_EOS_sm(model,rcp)


