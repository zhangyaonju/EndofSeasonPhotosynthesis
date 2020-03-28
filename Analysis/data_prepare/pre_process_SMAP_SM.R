args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)
library(raster)
library(lubridate)
library(ncdf4)

#library("hdf5r")
# Soil moisture is retrieved over land targets on the descending (AM) SMAP half-orbits 
# when the SMAP spacecraft is travelling from North to South, while the SMAP instruments 
# are operating in the nominal mode.  The L3_SM_P product represents soil moisture 
# retrieved over the entre UTC day. Retrievals are performed but flagged as questionable 
# over urban areas, mountainous areas with high elevation variability, and areas with high 
# ( &gt; 5 kg/m**2) vegetation water content; for retrievals using the high-resolution radar, 
# cells in the nadir region are also flagged. Retrievals are inhibited for permanent snow/ice,
# frozen ground, and excessive static or transient open water in the cell, and for excessive 
# RFI in the sensor data.";
proj<-CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
lat_proj<-CRS("+proj=longlat +datum=WGS84 +no_defs")

convert2raster<-function(dat){
  ras_dat<-raster(dat)
  extent(ras_dat)<-c(-17367530.44516138,17367530.44516138,-7314540.08900637,7314540.08900637)
  projection(ras_dat)<-proj
  # extent -180,180,-85.0445, 85.0445
  dat_lonlat<-projectRaster(ras_dat,res=c(0.5,0.5),crs=lat_proj,method="bilinear")
  return(dat_lonlat)
}

### get the monthly smap surface soil moisture.

# read_hdf<-function(file){
#   h5in<-H5File$new(file,mode="r")
#   am_sm<-h5in[["Soil_Moisture_Retrieval_Data_AM/soil_moisture"]]
#   am_dat<-am_sm$read()
#   am_dat[am_dat<=-9999]<-NA
#   #image(am_dat)
#   
#   pm_sm<-h5in[["Soil_Moisture_Retrieval_Data_PM/soil_moisture_pm"]]
#   pm_dat<-pm_sm$read()
#   pm_dat[pm_dat<=-9999]<-NA
#   #image(pm_dat)
#   dat<-abind(am_dat,pm_dat,along=3)
#   return(dat)
# }

exportnc<-function(dat,y,outfile){
  dimx<-ncdim_def(name = "longitude", unit="",vals=seq(-179.75,179.75,0.5))
  dimy<-ncdim_def(name = "latitude", unit="",vals=seq(89.75,-89.75,-0.5))
  if (dim(dat)[3]==92){
    dimt<-ncdim_def(name = "day",units = "",vals = 1:92*4-3+(y-2000)*1000)
  }else if (dim(dat)[3]==122){
    dimt<-ncdim_def(name = "day",units = "",vals = 1:122*3-2+(y-2000)*1000)
  }else if(dim(dat)[3]==12){
    dimt<-ncdim_def(name = "month",units = "",vals = 1:12+(y-2000)*100)
  }
  var_tre<-ncvar_def(name = "SM_AM",units = "cm3/cm3",dim = list(dimy,dimx,dimt),
                     missval = -9999,longname = "SMAP soil moisutre at 6:00am",compression=9)
  dat[is.na(dat)]<- -9999
  ncout<-nc_create(filename = outfile,vars = list(var_tre))
  ncvar_put(nc = ncout,varid = var_tre,vals = dat)
  nc_close(ncout)
}


##### read in the tiff data and aggreate to monthly

setwd("/rigel/glab/users/")
amf<-list.files("./Datasets/SMAP_L3_daily/",pattern="*.SM_AM.tif$",full.names = T)
#pmf<-list.files("./Datasets/SMAP_L3_daily/",pattern="*.SM_PM.tif$",full.names = T)

get_dat_to_doy<-function(fi){
  fdate<-substr(basename(fi),14,21)
  ftime<-strptime(fdate,format = "%Y%m%d")
  ydoy<-year(ftime)*1000+yday(ftime)
}

all_m<-unique(substr(basename(amf),14,19))

ras_frame<-raster(array(1,dim=c(360,720)))
extent(ras_frame)<-c(-180,180,-90,90)
projection(ras_frame)<-lat_proj

#### get three day resolution data for 6:00
aggreate_reproj<-function(y,interval){
  year_f<-amf[substr(basename(amf),14,17)==y]
  doy<-get_dat_to_doy(year_f)%%1000
  obs<-ceiling(365/interval)
  day_dat<-array(NA,dim=c(360,720,obs))
  for (d in 1:obs){
    trf<-year_f[doy<=d*interval&doy>(d-1)*interval]
    if (length(trf)>0){
      if (length(trf)==1){
        mean_moisture<-raster(trf)
        mean_moisture[mean_moisture<= -9999]<-NA
      }else{
        thr_stack<-stack(trf)
        thr_stack[thr_stack<= -9999]<-NA
        mean_moisture<-calc(thr_stack,mean,na.rm=T)
      }
      extent(mean_moisture)<-c(-17367530.44516138,17367530.44516138,-7314540.08900637,7314540.08900637)
      projection(mean_moisture)<-proj
      # extent -180,180,-85.0445, 85.0445
      dat_lonlat<-projectRaster(mean_moisture,ras_frame,method="bilinear")
      day_dat[,,d]<-as.matrix(dat_lonlat)
    }
  }
  # export the data to nc
  if (interval==3){
    outf<-paste0("./zy2309/PROJECT/EOS_SM/SMAP/threeday/",y,".AM_SM_3day.nc")
  }else if(interval==4){
    outf<-paste0("./zy2309/PROJECT/EOS_SM/SMAP/fourday/",y,".AM_SM_4day.nc")
  }
  exportnc(day_dat,y,outfile=outf)
}

mon_s<-c(1,11,21,31,41,51,61,71,82,92,102,112,123)
three_f<-list.files("./zy2309/PROJECT/EOS_SM/SMAP/threeday/",pattern=".nc",full.names = T)
aggreatemon<-function(y){
  year_f<-three_f[substr(basename(three_f),1,4)==y]
  ncin<-nc_open(year_f)
  dat<-ncvar_get(ncin,varid="SM_AM")
  sm_mon<-array(NA,dim=c(360,720,12))
  for (i in 1:12){
    sm_mon[,,i]<-apply(dat[,,mon_s[i]:(mon_s[i+1]-1)],c(1,2),mean,na.rm=T)
  }
  outf<-paste0("./zy2309/PROJECT/EOS_SM/SMAP/mon/",y,".AM_SM_mon.nc")
  exportnc(sm_mon,y,outfile=outf)
}

aggreate_reproj(uid+2014,3)
aggreate_reproj(uid+2014,4)
aggreatemon(uid+2014)

