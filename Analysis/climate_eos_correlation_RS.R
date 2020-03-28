#### get correlation between the EOS date and the climate variable RS
####   calcualte the 4day Ta Rn
args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ## uid for the four different lags

library(abind)
library(ppcor)
setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
cli_f<-list.files("./Climate/eos_RS2/",full.names = T)

get_corr<-function(vec){
  if(sum(!is.na(vec[1:14]))<8|sum(!is.na(vec[15:28]))<8){
    return(c(NA,NA))
  }
  dim(vec)<-c(14,2)
  co<-cor.test(vec[,1],vec[,2])
  return(c(co$estimate,co$p.value))
}

get_pcorr4<-function(vec){
  dim(vec)<-c(14,5)
  vec_nona<-vec[complete.cases(vec),]
  if (dim(vec_nona)[1]<8|length(unique(vec_nona[,3]))==1){
    return(rep(NA,8))
  }
  p_cor1<-pcor.test(vec_nona[,1],vec_nona[,2],vec_nona[,3:5])      # lst_day
  p_cor2<-pcor.test(vec_nona[,1],vec_nona[,3],vec_nona[,c(2,4,5)]) # lst_night
  p_cor3<-pcor.test(vec_nona[,1],vec_nona[,4],vec_nona[,c(2,3,5)]) # prec
  p_cor4<-pcor.test(vec_nona[,1],vec_nona[,5],vec_nona[,c(2,3,4)]) # PAR
  
  return(c(p_cor1$estimate,p_cor1$p.value,p_cor2$estimate,p_cor2$p.value,
           p_cor3$estimate,p_cor3$p.value,p_cor4$estimate,p_cor4$p.value))
}


get_pcorr3<-function(vec){
  dim(vec)<-c(14,4)
  vec_nona<-vec[complete.cases(vec),]
  if (dim(vec_nona)[1]<8|length(unique(vec_nona[,3]))==1){
    return(rep(NA,6))
  }
  p_cor1<-pcor.test(vec_nona[,1],vec_nona[,2],vec_nona[,3:4])    # lst_mean
  p_cor2<-pcor.test(vec_nona[,1],vec_nona[,3],vec_nona[,c(2,4)]) # prec
  p_cor3<-pcor.test(vec_nona[,1],vec_nona[,4],vec_nona[,c(2,3)]) # par

  return(c(p_cor1$estimate,p_cor1$p.value,p_cor2$estimate,p_cor2$p.value,
           p_cor3$estimate,p_cor3$p.value))
}

get_correlation<-function(lag){
  load("./CSIF/pheno_results.RData")
  ### lag represent preseason length 15 day 1-, 2-, 3- months
  lst_day_dat<-array(NA,dim=c(14,720,360))
  lst_night_dat<-array(NA,dim=c(14,720,360))
  par_dat<-array(NA,dim=c(14,720,360))
  tp_dat<-array(NA,dim=c(14,720,360))
  for (i in 1:14){
    load(cli_f[i])
    lst_day_dat[i,,]<-lst_day_eos[lag,,]
    lst_night_dat[i,,]<-lst_night_eos[lag,,]
    tp_dat[i,,]<-tp_eos[lag,,]
    par_dat[i,,]<-par_eos[lag,,]
  }
  lst_mean<-(lst_day_dat+lst_night_dat)/2
  ## combine with the eos date, and calculate the correltion and partial correlation
  lst_mean_all<-abind(lst_mean,eosdate[3:16,,],along=1)
  corr_lst_mean<-apply(lst_mean_all,c(2,3),get_corr)
  lst_day_all<-abind(lst_day_dat,eosdate[3:16,,],along=1)
  corr_lst_day<-apply(lst_day_all,c(2,3),get_corr)
  lst_night_all<-abind(lst_night_dat,eosdate[3:16,,],along=1)
  corr_lst_night<-apply(lst_night_all,c(2,3),get_corr)
  tp_all<-abind(tp_dat,eosdate[3:16,,],along=1)
  corr_tp<-apply(tp_all,c(2,3),get_corr)
  par_all<-abind(par_dat,eosdate[3:16,,],along=1)
  corr_par<-apply(par_all,c(2,3),get_corr)

  save(corr_lst_day,corr_lst_night,corr_lst_mean,corr_tp,corr_par,
       file=paste0("./analysis/RS_corr2/",lag-1,".mon_corr.RData"))
  
  ### combine all and calculate the partial correlation
  climate_dat<-abind(eosdate[3:16,,],lst_day_dat,lst_night_dat,tp_dat,par_dat,along=1)
  pcorr_all<-apply(climate_dat,c(2,3),get_pcorr4)
  pcorr_lst_day4<-pcorr_all[1:2,,]
  pcorr_lst_night4<-pcorr_all[3:4,,]
  pcorr_tp4<-pcorr_all[5:6,,]
  pcorr_par4<-pcorr_all[7:8,,]
  climate_dat<-abind(eosdate[3:16,,],lst_mean,tp_dat,par_dat,along=1)
  pcorr_all<-apply(climate_dat,c(2,3),get_pcorr3)
  pcorr_lst_mean3<-pcorr_all[1:2,,]
  pcorr_tp3<-pcorr_all[3:4,,]
  pcorr_par3<-pcorr_all[5:6,,]
  
  save(pcorr_lst_day4,pcorr_lst_night4,pcorr_tp4,pcorr_par4,
       pcorr_lst_mean3,pcorr_tp3,pcorr_par3,file=paste0("./analysis/RS_corr2/",lag-1,".mon_pcorr.RData"))
}

get_correlation(uid)

