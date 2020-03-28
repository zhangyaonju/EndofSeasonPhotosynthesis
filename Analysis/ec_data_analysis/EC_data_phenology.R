### process the EC data as the figure 3
#
# this code is used to process the flux data to analyze the climate controls on the End of photosynthesis.
# this code need to load the phenology retrieval functions first.
# several key processes are included. 1) check the data quality and plot the data quality for each site year
#                                     2) retrieve phenology for sites with more than 5 year
#                                     3) get the preseason climate based on the multi year average EOP
#                                     4) calculate the climate and EOP correlation.
#                                     5) get the maximum correlation with T and P
#                                     6) get the SM from SMAP for all sites as a backup
#
# @ YAO ZHANG 2019

library(tidyverse)
library(lubridate)

# check site data quality -------------------------------------------------
## do a site selection based on the rules below.
## 1. not cropland, not managed grassland
## 2. at least 5 years of data
## 3. global extent
### site list for tier 1 subset
setwd("/Users/yzhang/")
all_sites<-read.csv("./Data/FLUXNET_subset_tier1_DD/Site_information_tier1.csv",stringsAsFactors = F)
## remove sites with less than 5 years
all_sites_5y<-all_sites[all_sites$END_YR-all_sites$START_YR>=5,]
## remove cropland
all_sites_5y_no_crop<-all_sites_5y[all_sites_5y$IGBP!="CRO",]
## remove managed grassland
mana_grass<-c("AT-Neu","CH-Cha","CH-Fru","CH-Oe1","De-Gri","DK-ZaH","US-Cop","US-Var")
all_sites_5y_no_crop_no_managed_grass<-all_sites_5y_no_crop[!all_sites_5y_no_crop$SITE_ID %in% mana_grass,]

all_site_list<-data.frame(site_id=all_sites_5y_no_crop_no_managed_grass$SITE_ID,
                          eos_obs_years=NA,
                          variability=NA)
files<-list.files("./Data/FLUXNET_subset_tier1_DD/Data/",full.names = T)
### check the data quality
for (i in 1:length(all_site_list$site_id)){
  dat<-read.csv(files[substr(basename(files),5,10)==all_site_list$site_id[i]],stringsAsFactors = F)
  dat<-dat %>% mutate(.,year=floor(TIMESTAMP/10000),yearmonth=floor(TIMESTAMP/100))
  dat$NEE_VUT_REF_QC[dat$NEE_VUT_REF_QC< -9000|dat$GPP_NT_VUT_REF< -9000]<-0
  # get the monthly data quality
  yearmonth_data_quality<-dat %>% 
    group_by(yearmonth) %>% 
    dplyr::summarize(avg=mean(NEE_VUT_REF_QC,na.rm=T),month_gpp=mean(GPP_NT_VUT_REF))
  
  ### the the multi year average growing season based on GPP
  monthly_gpp<-yearmonth_data_quality %>%
    group_by(month=yearmonth%%100) %>%
    dplyr::summarize(msc_gpp=mean(month_gpp[avg>0.5]))
  
  ## get the month that eos ends
  decrease<-diff(c(monthly_gpp$msc_gpp[12],monthly_gpp$msc_gpp,monthly_gpp$msc_gpp))<0
  
  small_month<-which(monthly_gpp$msc_gpp<max(min(monthly_gpp$msc_gpp,na.rm = T),0)*0.7+
                       max(monthly_gpp$msc_gpp,na.rm = T)*0.3&decrease)
  eos_month<-(small_month[small_month>which.max(monthly_gpp$msc_gpp)][1]-1)%%12+1
  
  ### the mean seasonal average >1 g is considered as growing season, get the number of months without valid
  ##  observations during the growing season for each year
  year_gs_summary<-yearmonth_data_quality %>%
    dplyr::filter(., yearmonth%%100 %in% which(monthly_gpp$msc_gpp>1)) %>%
    group_by(year=floor(yearmonth/100)) %>%
    dplyr::summarize(avg=mean(avg>0.5))
  
  year_eos_summary<-yearmonth_data_quality %>%
    dplyr::filter(., yearmonth%%100 %in% ((eos_month+(-2:0)-1)%%12+1)) %>%
    group_by(year=floor(yearmonth/100)) %>%
    dplyr::summarize(avg=mean(avg>0.5))
  
  month_qa<-yearmonth_data_quality %>%
    mutate(.,month=yearmonth%%100) %>%
    group_by(month) %>%
    dplyr::summarise(m_year_avg = mean(avg[avg!=0]))
  
  yearmonth_data_quality$myear_avg<-month_qa$m_year_avg
  miss_month<-yearmonth_data_quality %>% 
    group_by(year=floor(yearmonth/100)) %>% 
    dplyr::summarise(missingmonths = sum(myear_avg-avg>0.5|avg==0))
  
  #yearly_data_quality<-dat %>% group_by(year) %>% summarise(avg=mean(NEE_VUT_REF_QC))
  print(paste(all_site_list$site_id[i],"   ",sum(year_eos_summary$avg==1),"    ",
              sum(miss_month$missingmonths>2)))
  good_year<-year_eos_summary$year[year_eos_summary$avg>0.5&year_gs_summary$avg>0.5]
  all_site_list$eos_obs_years[i]<-length(good_year)
  all_site_list$variability[i]<-min(monthly_gpp$msc_gpp,na.rm=T)>0.3*max(monthly_gpp$msc_gpp,na.rm=T)&
    max(monthly_gpp$msc_gpp,na.rm=T)>5
  if ("SWC_F_MDS_1" %in% names(dat)){
    good_obs<-dat %>% 
      dplyr::select(.,TIMESTAMP,TA_F,SW_IN_POT,SW_IN_F,VPD_F,P_F,PA_F,WS_F,LE_F_MDS,H_F_MDS,
                    NEE_VUT_REF,NEE_VUT_REF_QC,RECO_NT_VUT_REF,RECO_DT_VUT_REF,GPP_NT_VUT_REF,GPP_DT_VUT_REF,
                    SWC_F_MDS_1,
                    SWC_F_MDS_1_QC)
  }else{
    good_obs<-dat %>% 
      dplyr::select(.,TIMESTAMP,TA_F,SW_IN_POT,SW_IN_F,VPD_F,P_F,PA_F,WS_F,LE_F_MDS,H_F_MDS,
                    NEE_VUT_REF,NEE_VUT_REF_QC,RECO_NT_VUT_REF,RECO_DT_VUT_REF,GPP_NT_VUT_REF,GPP_DT_VUT_REF)
  }
  
  #if (all_site_list$eos_obs_years[i]>=5){
  save(good_obs,good_year,year_eos_summary,eos_month,monthly_gpp,
       file=paste0("./Project/EOS_SM/FLUX_analysis/selected_sites/",
                            all_site_list$site_id[i],".RData"))
  #}
}
save(all_site_list,file="./Project/EOS_SM/FLUX_analysis/all_site_list_stat.RData")

# graph the sites data with data quality. ---------------------------------
### get corrected month

setwd("/Users/yzhang/")
plot_graph_with_QA<-function(site_id){
  load(paste0("./Project/EOS_SM/FLUX_analysis/selected_sites/",site_id,".RData"))
  ### get the number of years to determine the figure size
  good_obs<-good_obs %>% 
    mutate(.,year=floor(TIMESTAMP/10000),
           yearmonth=floor(TIMESTAMP/100),
           date=ymd(TIMESTAMP))
  good_obs$NEE_VUT_REF_QC[good_obs$NEE_VUT_REF_QC< -9000|good_obs$GPP_NT_VUT_REF< -9000]<-0
  good_obs[good_obs< -9990]<-NA
  years<-unique(good_obs$year)
  nrows<-ceiling(length(years)/5)
  ymax<-max(monthly_gpp$msc_gpp,na.rm=T)*1.5
  pdf(paste0("./Project/EOS_SM/FLUX_analysis/graph_sites_QA/",site_id,".pdf"),
      width=8,height=3*nrows+1)
  par(mfrow=c(nrows,1),mar=c(2,3,1,2.5),oma=c(1,1,1,1),mgp=c(3,0.5,0))
  for (i in 1:nrows){
    ### get the range of y
    row_dat<-good_obs %>%
      dplyr::filter(., year %in% years[i*5-4:0])
    bad_years<-setdiff(years[i*5-4:0], good_year)

    plot(NA,xlim=c(min(row_dat$date),ymd(paste0(years[i*5-4]+5,"0101"))),
         ylim=c(0, ymax),axes=F,xlab="",ylab="")
    box()
    axis.Date(1,x = row_dat$date)
    axis(2,las=2)
    
    mtext(side=2,line=2,expression(paste("GPP_NT (g C"^-2,"d"^-1,")",sep="")))
    ### plot the bad years, and the eos months
    if (length(bad_years)>0){
      for (y in bad_years){
        polygon(c(ymd(paste0(y,"0101")),ymd(paste0(y,"1231")),ymd(paste0(y,"1231")),ymd(paste0(y,"0101"))),
                c(0,0,ymax,ymax),border=NA,col=adjustcolor("grey40",0.2))
      }
    }
    ## and the eos months
    eos_months<-(eos_month+(-2:0)-1)%%12+1
    for (y in 1:5){
      for (j in 1:3){
        st_date<-ymd(paste0(years[i*5-4:0][y],formatC(eos_months[j],width=2,flag = "0"),"01"))
        polygon(c(st_date,st_date+30,st_date+30,st_date),
                c(0,0,ymax,ymax),border=NA,col=adjustcolor("yellow",0.4))
      }
    }
    ### plot the gpp NT data
    points(row_dat$date,row_dat$GPP_NT_VUT_REF,col="blue",cex=0.3,lwd=0.3)
    abline(h=max(min(monthly_gpp$msc_gpp,na.rm=T),0)*0.7+max(monthly_gpp$msc_gpp,na.rm=T)*0.3)
    ### add the NEE QA
    
    lines(row_dat$date,(row_dat$NEE_VUT_REF_QC+3)*ymax/4,col="red")
    abline(h=4*ymax/4,lty=3)
  }
  dev.off()
}

for (k in 1:length(all_site_list$site_id)){
  plot_graph_with_QA(all_site_list$site_id[k])
}
              
#### check AU-Dry DE-Akm
#### altogether 50 sites

# get phenology for each site ---------------------------------------------

load("~/Documents/Project/EOS_SM/FLUX_analysis/all_site_list_stat.RData")
sites_used<-all_site_list$site_id[all_site_list$eos_obs_years>4&!all_site_list$variability]

## get the phenology for each site.
source("./Documents/GitHub/EOS_soil_moisture/data_analysis/ec_data_analysis/pheno_retrieval.R")

map_pheno<-function(dat, t, y, w, thr, method ){
  t<- dat %>% pull(t)
  y<- dat %>% pull(y)
  w<- dat %>% pull(w)
  pheno<-method(y,t,w,thr)
  #pred<-doubleLog.zhang(pheno,t=2011+(0:364)/365)
  
  return(data.frame(eos=pheno[1]))
}

retrieve_phen_prepare<-function(rdata){
  load(rdata)
  sid<-substr(basename(rdata),1,6)
  good_obs<-good_obs %>% 
    mutate(., GPP = GPP_NT_VUT_REF,doy = yday(ymd(TIMESTAMP))) 
  
  # if the peak growing season is in the mid year 4-9, year starts from doy1 to doy 365
  if (which.max(monthly_gpp$msc_gpp)<=9&which.max(monthly_gpp$msc_gpp)>=4){
    offset=0
  }else if (eos_month<8){### for south hemisphere or dryland sites, the phenological year is shifted by 120 days.
    offset=120
  }else{     ### for south hemisphere or dryland sites, the phenological year is shifted by -240 days.
    offset=-240
  }
  data_pheno<-data.frame(t=floor(good_obs$TIMESTAMP/10000)+good_obs$doy/365,
                         y=good_obs$GPP,
                         w=good_obs$NEE_VUT_REF_QC*0.95+0.05,
                         year=year(ymd(good_obs$TIMESTAMP)+offset),
                         yearmonth=floor(good_obs$TIMESTAMP/100),
                         date=ymd(good_obs$TIMESTAMP))

  ### only keep the valid years
  data_pheno<-data_pheno %>%
    filter(year %in% good_year)
  data_pheno$y[data_pheno$y< -999]<-NA
  yearmonth_data_quality<-data_pheno %>% 
    group_by(yearmonth) %>% 
    dplyr::summarize(qa_avg=mean(w,na.rm=T))
  years<-unique(data_pheno$year)
  
  #### if months's dataquality is 0 during the non growing season, the GPP for this month is set to zero
  gs_months<-which(monthly_gpp$msc_gpp>max(0,min(monthly_gpp$msc_gpp,na.rm=T))*0.7+
                     min(monthly_gpp$msc_gpp,na.rm=T)*0.3)
  
  data_pheno$y[data_pheno$yearmonth %in% yearmonth_data_quality$yearmonth[yearmonth_data_quality$qa_avg==0.05]&
                 !data_pheno$yearmonth%%100 %in% gs_months]<-0
  data_pheno$y[data_pheno$y< -3|data_pheno$y>30]<-0
  
  smoothed_val<-data_pheno %>% as.tibble %>%
    nest(-year) %>% 
    mutate(fitted_dat = map(data, ~ splinefit(t=.$t,y=.$y,w=.$w))) %>% 
    .$fitted_dat
  
  data_pheno$fitted<-unlist(smoothed_val)
  thr<-get_MR(smoothed_val)
  
  results<-data_pheno %>% 
    split(.$year) %>% 
    map_dfr(map_pheno,t="t",y='y',w='w',thr,method=pheno_wSpline,.id='year') %>%
    mutate(.,site=sid)
  
  #### start the graph
  all_years<-unique(floor(good_obs$TIMESTAMP/10000))
  nrows<-ceiling(length(all_years)/5)
  ymax<-max(monthly_gpp$msc_gpp,na.rm=T)*1.5
  pdf(paste0("/Users/yzhang/Project/EOS_SM/FLUX_analysis/graph_sites_pheno/",sid,'.pdf'),
      width=8,height=3*nrows+1)

  par(mfrow=c(nrows,1),mar=c(2,3,1,2.5),oma=c(1,1,1,1),mgp=c(3,0.5,0))
  for (i in 1:nrows){
    ### get the range of y
    row_dat<-data_pheno %>%
      dplyr::filter(., year %in% all_years[i*5-4:0])
    
    plot(NA,xlim=c(ymd(paste0(all_years[i*5-4],"0101")),ymd(paste0(all_years[i*5-4]+5,"0101"))),
         ylim=c(0, ymax),axes=F,xlab="",ylab="")
    box()
    mtext(side=2,line=2,expression(paste("GPP_NT (g C"^-2,"d"^-1,")",sep="")))
    axis.Date(1,x = row_dat$date)
    axis(2,las=2)
    points(row_dat$date,row_dat$y,col="blue",cex=0.3,lwd=0.3)
    lines(row_dat$date,row_dat$fitted,col="darkgreen",lwd=2)
    abline(h=thr[3],lty=2)
    #### plot the eos date
    row_eos<-results$eos[results$year %in% all_years[i*5-4:0]]
    row_eos_date<-as_date(date_decimal(row_eos))
    abline(v=row_eos_date,col="red",lwd=0.7)
  }
  dev.off()
  
  return(results)
}

spline_res<-list()
site_list<-list.files('./Project/EOS_SM/FLUX_analysis/selected_sites/',full.names = T)
for (s in 1:length(sites_used)){
  print(site_list[substr(basename(site_list),1,6)==sites_used[s]])
  retrive_results<-retrieve_phen_prepare(site_list[substr(basename(site_list),1,6)==sites_used[s]])
  spline_res[[s]]<-retrive_results
}

save(spline_res,file="./Project/EOS_SM/FLUX_analysis/eos_date_all_50_sites.RData")




# get preseason climate ---------------------------------------------------

# get 15 30 60 90 pre date var average
get_var_pre<-function(var,tin,eos){
  fourdate<-c(15,30,60,90)
  var_mean<-c()
  eos_id<-which(tin==eos)
  for (i in 1:4){
    var_mean[i]<-mean(var[max(eos_id-fourdate[i],1):(eos_id)],na.rm=T)
  }
  return(var_mean)
}

####### get the preseason temerpature and precipitation, as well as for sm if exist.
get_preseason_T_P_SM<-function(rdata,eos_date){
  load(rdata)
  sid<-substr(basename(rdata),1,6)
  good_obs[good_obs< -9990]<-NA
  climate<-good_obs %>% as.tibble(.) %>%
    mutate(date=ymd(good_obs$TIMESTAMP),year=floor(good_obs$TIMESTAMP/10000))
  ## get the mean annual temeprature
  mat<-climate %>%
    group_by(year) %>%
    dplyr::summarise(annual_T = mean(TA_F)) %>%
    .$annual_T %>% mean(na.rm=T)
  ## get the mean eos date
  mean_eos<-round(median(yday(date_decimal(eos_date$eos)),na.rm=T))
  ### get var average for a specfic period and year
  climate_pre_EOS<-as.data.frame(array(NA,dim=c(dim(eos_date)[1],12)))
  names(climate_pre_EOS)<-c(paste0("Ta",0:3),paste0("Prec",0:3),paste0("SM",0:3))
  for (ty in 1:length(eos_date$year)){
    year_eos<-as_date(date_decimal(mean_eos/365+as.numeric(eos_date$year[ty])))
    climate_pre_EOS[ty,1:4]<-get_var_pre(climate$TA_F,climate$date,year_eos)
    climate_pre_EOS[ty,5:8]<-get_var_pre(climate$P_F,climate$date,year_eos)*c(15,30,60,90)
    if("SWC_F_MDS_1" %in% names(climate)){
      climate_pre_EOS[ty,9:12]<-get_var_pre(climate$SWC_F_MDS_1,climate$date,year_eos)
    }
  }
  climate_with_EOS<-cbind(eos_date,climate_pre_EOS,mat)
  return(climate_with_EOS)
}

#### load the eos date
load("./Project/EOS_SM/FLUX_analysis/eos_date_all_50_sites.RData")
preseason_var<-list()
site_list<-list.files('./Project/EOS_SM/FLUX_analysis/selected_sites/',full.names = T)
for (s in 1:length(sites_used)){
  print(site_list[substr(basename(site_list),1,6)==sites_used[s]])
  site_preseason<-get_preseason_T_P_SM(site_list[substr(basename(site_list),1,6)==sites_used[s]],spline_res[[s]])
  preseason_var[[s]]<-site_preseason
}
save(preseason_var,file="./Project/EOS_SM/FLUX_analysis/eos_climate_all_50_sites.RData")


# calculate the correlation for each site ---------------------------------

all_site_summary<-as.data.frame(array(NA,dim=c(50,19)))
names(all_site_summary)<-c("site",paste0("r_T_EOS_",0:3),paste0("r_P_EOS_",0:3),
                           paste0("p_T_EOS_",0:3),paste0("p_P_EOS_",0:3),"MAT","Pre_SM")
all_site_summary$site<-sites_used
for (i in 1:length(sites_used)){
  site_dat<-preseason_var[[i]]
  for (v in 1:8){
    cor_spearman<-cor.test(site_dat$eos-as.numeric(site_dat$year),site_dat[,3+v])
    all_site_summary[i,v+1]<-cor_spearman$estimate
    all_site_summary[i,v+9]<-cor_spearman$p.value
  }
  all_site_summary$MAT[i]<-site_dat$mat[1]
  all_site_summary$Pre_SM[i]<-mean(site_dat$SM1,na.rm=T)
}

maximum_corrlation<-as.data.frame(array(NA,dim=c(50,9)))
names(maximum_corrlation)<-c("site","r_max_T_EOS","p_max_T_EOS","length_max_T_EOS",
                             "r_max_P_EOS","p_max_P_EOS","length_max_P_EOS","MAT","Pre_SM")
maximum_corrlation$site<-sites_used
for (i in 1:length(sites_used)){
  for (v in 1:2){
    max_id<-which.max(abs(all_site_summary[i,-2:1+v*4]))
    maximum_corrlation[i,-1:1+v*3]<-c(all_site_summary[i,c(v*4-3+max_id,v*4+5+max_id)],max_id)
  }
  maximum_corrlation$MAT[i]<-all_site_summary$MAT[i]
  maximum_corrlation$Pre_SM[i]<-all_site_summary$Pre_SM[i]
}

save(all_site_summary,maximum_corrlation,
     file="./Project/EOS_SM/FLUX_analysis/eos_climate_max_corr_50_sites.RData")



load("./Project/EOS_SM/FLUX_analysis/eos_climate_max_corr_50_sites.RData")
# add SMAP soil moisture to the analysis ----------------------------------
library(raster)
sites_used<-all_site_summary$site
all_sites<-read.csv("./Data/FLUXNET_subset_tier1_DD/Site_information_tier1.csv",stringsAsFactors = F)
used_sites_loc<-all_sites %>% as_tibble %>%
  filter (.$SITE_ID %in% sites_used)
#### load the EOP SMAP data
load("./Project/EOS_SM/analysis/mean_pre_eos_sm.RData")
mean_eos_sm<-raster(t(mean_pre_eos_sm[2,,360:1]))
extent(mean_eos_sm)<-c(-180,180,-90,90)
### extract pre EOP SM using CSIF EOS 
loc<-as.matrix(cbind(used_sites_loc$LOCATION_LONG,used_sites_loc$LOCATION_LAT))
maximum_corrlation$SM_SMAP_CSIF_EOS<-extract(x = mean_eos_sm,loc)
# save(loc,file="./Project/EOS_SM/FLUX_analysis/site_50_loc.RData")

#### load the ts of SM data and get the preseason SM
load("./Project/EOS_SM/FLUX_analysis/site_50_4day_4year_sm.RData")
load("./Project/EOS_SM/FLUX_analysis/eos_date_all_50_sites.RData")

dim(ts_smap_sm)<-c(92,4,50)
mean_seasonal_sm<-apply(ts_smap_sm,c(1,3),mean,na.rm=T)

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

for (s in 1:50){
  eos<-round(median(spline_res[[s]]$eos-as.numeric(spline_res[[s]]$year),na.rm=T)*365)
  maximum_corrlation$SM_SMAP_site_EOS[s]<-get_mean(c(mean_seasonal_sm[,s],mean_seasonal_sm[,s]),
                                                   (eos/4-30/4+92),eos/4+92)
  maximum_corrlation$years_used[s]<-sum(!is.na(spline_res[[s]]$eos))
  maximum_corrlation$igbp[s]<-used_sites_loc$IGBP[s]
}
save(all_site_summary,maximum_corrlation,
     file="./Project/EOS_SM/FLUX_analysis/eos_climate_max_corr_50_sites.RData")


### load the time series of SM for each site !!!!!!
### this should be run on the server
if (F){
  library(raster)
  library(abind)
  library(ncdf4)
  setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/SMAP/")
  fourday_f<-list.files("./fourday/",full.names = T)
  load("./site_50_loc.RData")
  ts_smap_sm<-array(NA,dim=c(92*4,50))

  # read the sm data and generate a stack
  ncin<-nc_open(fourday_f[1])
  sm<-ncvar_get(ncin,varid='SM_AM')
  for (i in 2:length(fourday_f)){
    ncin<-nc_open(fourday_f[i])
    temp<-ncvar_get(ncin,varid='SM_AM')
    sm<-abind(sm,temp,along=3)
  }
  sm<-sm[360:1,,]
  ### convert the loc into id in coordinate
  coor_id<-cbind(round((loc[,2]+89.75)*2+1),
                 round(loc[,1]+179.75)*2+1)
  for (i in 1:50){
    ts_smap_sm[,i]<-sm[coor_id[i,1],coor_id[i,2],]
  }
  save(ts_smap_sm,file="./analysis/site_50_4day_4year_sm.RData")
}

# sites_SM_Stat<-as.data.frame(sites_used)
# sites_SM_Stat$sm_obs<-NA
# ### years of SM obs
# for (i in 1:50){
#   load(paste0("./Project/EOS_SM/FLUX_analysis/selected_sites/",sites_used[i],".RData"))
#   sites_SM_Stat$sm_obs[i]<-sum(good_obs$SWC_F_MDS_1_QC==1,na.rm=T)/365
# }
# 
