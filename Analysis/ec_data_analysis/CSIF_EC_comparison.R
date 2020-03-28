
library(tidyverse)
library(lubridate)
setwd("/Users/yzhang/")
load("/Users/yzhang/Project/PAR_limit/FLUX_analysis/site_list.RData")
all_sites<-read.csv("./Data/FLUXNET_subset_tier1_DD/Site_information_tier1.csv",stringsAsFactors = F)
load("./Project/EOS_SM/FLUX_analysis/all_site_list_stat.RData")
sites_used<-all_site_list$site_id[all_site_list$eos_obs_years>4&!all_site_list$variability]

south_sites<-c("AU-DaP", "AU-Stp")#, "ZA-Kru")
### sites less than 5 years.
l5_sites<-c("US-Blo","CA-Man")
north_sites<-sites_used[sites_used %in% m_site_list$site_id& !(sites_used %in% l5_sites)]

csif_dat<-read.csv("/Users/yzhang/Project/SIF_phenology/retrieve_site/sites166_with_CSIF_daily_clear_BISE.csv",
                   header = T,stringsAsFactors = F)
# csif_dat<-read.csv("/Users/yzhang/Project/SIF_phenology/retrieve_site/sites166_with_CSIF_clear.csv",
#                    header = T,stringsAsFactors = F)[-c(1:92),]
csif_dat[csif_dat>10000]<-0
##### this is the sites used for comparison with EC sites
sites_csif_comp<-c(south_sites,as.character(north_sites))

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

retrieve_phen_prepare<-function(sid,method=pheno_wSpline){
  load(paste0("./Project/EOS_SM/FLUX_analysis/selected_sites/",
              sid,".RData"))
  good_obs<-good_obs %>% 
    mutate(., GPP = GPP_NT_VUT_REF,doy = yday(ymd(TIMESTAMP))) 
  #### retrieve CSIF_SITE for comparison
  csif_site<-cbind(rep(2001:2017,each=92)+rep(((1:92)*4-2)/365,17),
                   as_date(strptime(rep(2001:2017,each=92)*1000+rep(c((0:90)*4+2,365),17),format="%Y%j")),
                   subset(csif_dat,select=gsub("-",'.',sid)))
  
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
  csif_pheno<-data.frame(t=csif_site[,1],
                         y=csif_site[,3],
                         w=rep(1,dim(csif_site)[1]),
                         year=year(csif_site[,2]+offset),
                         date=csif_site[,2])
  
  ### only keep the valid years
  good_year <-good_year[good_year>2001]
  ###valid years
  num_years<-length(good_year)
  
  data_pheno<-data_pheno %>%
    filter(year %in% good_year)
  data_pheno$y[data_pheno$y< -999]<-NA
  yearmonth_data_quality<-data_pheno %>% 
    group_by(yearmonth) %>% 
    dplyr::summarize(qa_avg=mean(w,na.rm=T))
  years<-unique(data_pheno$year)
  years<-years[years>2000]
  ##     for csif
  csif_pheno<-csif_pheno %>%
    filter(year %in% good_year)
  
  #### if months's dataquality is 0 during the non growing season, the GPP for this month is set to zero
  gs_months<-which(monthly_gpp$msc_gpp>max(0,min(monthly_gpp$msc_gpp,na.rm=T))*0.7+
                     min(monthly_gpp$msc_gpp,na.rm=T)*0.3)
  
  data_pheno$y[data_pheno$yearmonth %in% yearmonth_data_quality$yearmonth[yearmonth_data_quality$qa_avg==0.05]&
                 !data_pheno$yearmonth%%100 %in% gs_months]<-0
  data_pheno$y[data_pheno$y< -3|data_pheno$y>30]<-0
  
  #### first smoothing and get thr
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
  
  ### for csif
  smoothed_csif<-csif_pheno %>% as.tibble %>%
    nest(-year) %>% 
    mutate(fitted_dat = map(data, ~ splinefit(t=.$t,y=.$y,w=.$w))) %>% 
    .$fitted_dat
  csif_pheno$fitted<-unlist(smoothed_csif)
  thr_csif<-get_MR(smoothed_csif)
  
  results_csif<-csif_pheno %>% 
    split(.$year) %>% 
    map_dfr(map_pheno,t="t",y='y',w='w',thr_csif,method=pheno_wSpline,.id='year') %>%
    mutate(.,site=sid)
  
  scalar<-thr[2]/thr_csif[2]
  
  all_years<-unique(floor(good_obs$TIMESTAMP/10000))
  all_years<-all_years[all_years>2000]
  nrows<-ceiling(length(all_years)/5)
  ymax<-max(monthly_gpp$msc_gpp,na.rm=T)*1.5
  pdf(paste0("./Project/EOS_SM/FLUX_analysis/graph_ts_csif/",sid,".pdf"),
      width=8,height=3*nrows+1)
  par(mfrow=c(nrows,1),mar=c(2,3,1,2.5),oma=c(1,1,1,1),mgp=c(3,0.5,0))
  for (i in 1:nrows){
    ### get the range of y
    row_dat<-data_pheno %>%
      dplyr::filter(., year %in% all_years[i*5-4:0])
    row_csif<-csif_pheno %>%
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
    abline(v=row_eos_date,col="green",lwd=0.7)
    
    ##### plot csif
    points(row_csif$date,row_csif$y*scalar,col="red",cex=0.3,lwd=0.3)
    lines(row_csif$date,row_csif$fitted*scalar,col="brown4",lwd=2)
    row_eos_csif<-results_csif$eos[results_csif$year %in% all_years[i*5-4:0]]
    row_eos_csif_date<-as_date(date_decimal(row_eos_csif))
    abline(v=row_eos_csif_date,col="coral",lwd=0.7)
  }
  dev.off()
  results<-results[!is.na(results$eos),]
  results[results[,1]==results[,2],]<-NA
  all_comp<-merge(results,results_csif,by.x="year",by.y="year")
  
  return(all_comp)
}

spline_res<-list()
site_list<-list.files('./Project/EOS_SM/FLUX_analysis/selected_sites/',full.names = T)
for (s in 1:length(sites_csif_comp)){
  print(sites_csif_comp[s])
  retrive_results<-retrieve_phen_prepare(sites_csif_comp[s])
  spline_res[[s]]<-retrive_results
}

save(spline_res,file="./Project/EOS_SM/FLUX_analysis/eos_date_comparison_CSIF.RData")


all<-bind_rows(spline_res, .id = "column_label")
par(mfrow=c(1,2))
plot(all$eos.x %% 1, all$eos.y %% 1,xlim=c(0,1),ylim=c(0,1))
cor.test(all$eos.x %% 1, all$eos.y %% 1)
abline(0,1)
ano<-list()
for (i in 1:length(spline_res)){
  dat<-spline_res[[i]]
  gpp_eos_mean<-mean(dat$eos.x %% 1,na.rm=T)
  gpp_eos_ano<-(dat$eos.x %% 1)-gpp_eos_mean
  csif_eos_mean<-mean(dat$eos.y %% 1,na.rm=T)
  csif_eos_ano<-(dat$eos.y %% 1)-csif_eos_mean
  site_dat<-data.frame(y=dat$year,
                       gpp_eos=dat$eos.x,
                       gpp_eos_ano,
                       csif_eos=dat$eos.y,
                       csif_eos_ano,
                       sid = dat$site.x)
  ano[[i]]<-site_dat
}
all<-bind_rows(ano, .id = "column_label")
plot(all$gpp_eos_ano,all$csif_eos_ano,xlim=c(-0.3,0.3),ylim=c(-0.3,0.3))
cor.test(all$gpp_eos_ano,all$csif_eos_ano)
abline(0,1)
save(spline_res,all,file="./Project/EOS_SM/FLUX_analysis/eos_date_comparison_CSIF.RData")
###      bise    clear 
### all  0.9495 0.9497
### ano  0.7256 0.6789


