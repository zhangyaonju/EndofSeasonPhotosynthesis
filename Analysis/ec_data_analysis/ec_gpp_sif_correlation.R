library(tidyverse)
library(lubridate)
setwd("~/")
load("~/Documents/Project/PAR_limit/FLUX_analysis/site_list.RData")
all_sites<-read.csv("./Data/FLUXNET_subset_tier1_DD/Site_information_tier1.csv",stringsAsFactors = F)
load("~/Documents/Project/EOS_SM/FLUX_analysis/all_site_list_stat.RData")
sites_used<-as.character(all_site_list$site_id[all_site_list$eos_obs_years>4&!all_site_list$variability])

south_sites<-c("AU-DaP", "AU-Stp")#, "ZA-Kru")
### sites less than 5 years.
l5_sites<-c("US-Blo","CA-Man")
north_sites<-sites_used[sites_used %in% m_site_list$site_id& !(sites_used %in% l5_sites)]

select_sites<-c(north_sites,south_sites)

csif_dat<-read.csv("~/Documents/Project/SIF_phenology/retrieve_site/sites166_with_CSIF_daily_clear_BISE.csv",
                   header = T,stringsAsFactors = F)

# get the correlation between
site_list<-list.files('~/Documents/Project/EOS_SM/FLUX_analysis/selected_sites/',full.names = T)
site_stat<-as.data.frame(array(NA,dim=c(length(select_sites),3)))
names(site_stat)<-c("site_id","cor.coef","years_used")
  
site_data<-list()
csif_years<-2001:2017
for (i in 1:length(select_sites)){
  site_file<-site_list[which(substr(basename(site_list),1,6)==select_sites[i])]
  load(site_file)
  
  good_obs<-good_obs %>% 
    mutate(., GPP = GPP_NT_VUT_REF,doy = yday(ymd(TIMESTAMP))) 
  data_pheno<-data.frame(t=floor(good_obs$TIMESTAMP/10000)+good_obs$doy/365,
                         y=good_obs$GPP,
                         w=good_obs$NEE_VUT_REF_QC*0.95+0.05,
                         year=year(ymd(good_obs$TIMESTAMP)),
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
  # intersect years
  ava_years<-years[years>2000]
  ### convert daily to 4day
  day4data<-as.data.frame(array(NA,dim=c(92*length(ava_years),2)))
  csif_sitedata<-csif_dat[,c(1,which(names(csif_dat)==gsub("-",".",select_sites[i])))]
  for (k in 1:length(ava_years)){
    ec_year_dat<-data_pheno$y[floor(data_pheno$t)==ava_years[k]]
    first91<-ec_year_dat[1:364]
    dim(first91)<-c(4,91)
    day4data[((k-1)*92+1):(k*92),1]<-ava_years[k]*1000+1:92*4-3
    day4data[((k-1)*92+1):(k*92-1),2]<-apply(first91,2,mean,na.rm=T)
    day4data[k*92,2]<-mean(ec_year_dat[365:length(ec_year_dat)],na.rm=T)
  }
  ### combine the data for comparison
  merged_csif_gpp<-merge(day4data,csif_sitedata,by.x="V1",by.y='DATE')
  cort<-cor.test(merged_csif_gpp[,2],merged_csif_gpp[,3])
  ###
  site_data[[i]]<-merged_csif_gpp
  site_stat$site_id[i]<-select_sites[i]
  site_stat$cor.coef[i]<-cort$estimate
  site_stat$years_used[i]<-paste0(ava_years,collapse=',')
}
save(site_data,site_stat,file="~/Documents/Project/EOS_SM/analysis/revision/CSIF_GPP_comparison.RData")





