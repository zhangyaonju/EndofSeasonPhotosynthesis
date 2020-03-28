#### objectives
##  1. get the trend of EOP for the past two decades
##  2. trend and fluctuation along the othegonal
##  3. contribution based on precipitation change and temperature change

load("/rigel/glab/users/zy2309/PROJECT/EOS_SM/CSIF/pheno_results.RData")
## get the trend and the significance of the eosdate
library(zyp)

mk_trend<-function(vec){
  if (sum(is.na(vec))>10){
    return(c(NA,NA,NA,NA))
  }
  x<-1:length(vec)
  sen<-zyp.sen(vec~x)
  mk<-MannKendall(vec)
  return(c(sen$coefficients[2],mk$sl,sd(vec,na.rm=T),mean(vec,na.rm=T)))
}

mk_zscore<-function(vec){
  if (sum(is.na(vec))>10){
    return(c(NA,NA,NA,NA))
  }
  x<-1:length(vec)
  zscore<-(vec-mean(vec,na.rm=T))/sd(vec,na.rm=T)
  sen<-zyp.sen(vec~x)
  mk<-MannKendall(vec)
  return(c(sen$coefficients[2],mk$sl,sd(vec,na.rm=T),mean(zscore,na.rm=T)))
}

eop_trend<-apply(eosdate,c(2,3),mk_trend)
save(eop_trend,file="/rigel/glab/users/zy2309/PROJECT/EOS_SM/analysis/trend/trend_and_SD_EOP.RData")

#### grab the climate variables with strongest correlations.
era_f<-list.files("/rigel/glab/users/zy2309/PROJECT/EOS_SM/Climate/eos_ERA/",full.names = T)
pre_s_length<-'/rigel/glab/users/zy2309/PROJECT/EOS_SM/analysis/max_corr_climate.RData'

pre_eos_temp<-array(NA,dim=c(17,720,360))
pre_eos_prec<-array(NA,dim=c(17,720,360))
select_ind<-function(vec){
  return(vec[vec[5]])
}
library(abind)
load(pre_s_length)
ind_p<-max_mat_p[2,,]
ind_t<-max_mat_t[2,,]
for (i in 1:17){
  load(era_f[i])
  tp_dat<-abind(tp_eos,ind_p,along=1)
  pre_eos_prec[i,,]<-apply(tp_dat,c(2,3),select_ind)
  t2m_dat<-abind(t2m_eos,ind_p,along=1)
  pre_eos_temp[i,,]<-apply(t2m_dat,c(2,3),select_ind)
}
## calcualte the trend of these climate variables.
tp_trend<-apply(pre_eos_prec,c(2,3),mk_trend)
t2m_trend<-apply(pre_eos_temp,c(2,3),mk_trend)
t2m_zscore<-apply(pre_eos_temp,c(2,3),mk_zscore)
save(tp_trend,t2m_trend,t2m_zscore,file="/rigel/glab/users/zy2309/PROJECT/EOS_SM/analysis/trend/trend_and_SD_climate.RData")

## combine with P and T trend to get the contributions from these two factors.









