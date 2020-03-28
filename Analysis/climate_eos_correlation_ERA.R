#### get correlation between the EOS date and the climate variable ERA
####   calcualte the 4day Ta Rn
args = commandArgs(trailingOnly=F)
print(args)
myargs <-sub('-','',args[length(args)])
print(myargs)

uid <- as.numeric(myargs)  ## uid for the four different lags

library(abind)
library(ppcor)
setwd("/rigel/glab/users/zy2309/PROJECT/EOS_SM/")
cli_f<-list.files("./Climate/eos_ERA/",full.names = T)

get_corr<-function(vec){
  if(sum(!is.na(vec[1:17]))<8|sum(!is.na(vec[18:34]))<8){
    return(c(NA,NA))
  }
  dim(vec)<-c(17,2)
  co<-cor.test(vec[,1],vec[,2])
  return(c(co$estimate,co$p.value))
}

get_pcorr<-function(vec){
  dim(vec)<-c(17,4)
  vec_nona<-vec[complete.cases(vec),]
  if (dim(vec_nona)[1]<8|length(unique(vec_nona[,3]))==1){
    return(rep(NA,6))
  }
  tryCatch({
    p_cor1<-pcor.test(vec_nona[,1],vec_nona[,2],vec_nona[,3:4])     # t2m
    p_cor2<-pcor.test(vec_nona[,1],vec_nona[,3],vec_nona[,c(2,4)]) # tp
    p_cor3<-pcor.test(vec_nona[,1],vec_nona[,4],vec_nona[,c(2,3)]) # par
    
    return(c(p_cor1$estimate,p_cor1$p.value,p_cor2$estimate,p_cor2$p.value,p_cor3$estimate,p_cor3$p.value))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(rep(NA,6))
}

get_correlation<-function(lag){
  load("./CSIF/pheno_results.RData")
  ### lag represent preseason length 15 day 1-, 2-, 3- months
  t2m_dat<-array(NA,dim=c(17,720,360))
  par_dat<-array(NA,dim=c(17,720,360))
  tp_dat<-array(NA,dim=c(17,720,360))
  for (i in 1:17){
    load(cli_f[i])
    t2m_dat[i,,]<-t2m_eos[lag,,]
    tp_dat[i,,]<-tp_eos[lag,,]
    par_dat[i,,]<-par_eos[lag,,]
  }
  ## combine with the eos date, and calculate the correltion and partial correlation
  t2m_all<-abind(t2m_dat,eosdate,along=1)
  corr_t2m<-apply(t2m_all,c(2,3),get_corr)
  tp_all<-abind(tp_dat,eosdate,along=1)
  corr_tp<-apply(tp_all,c(2,3),get_corr)
  par_all<-abind(par_dat,eosdate,along=1)
  corr_par<-apply(par_all,c(2,3),get_corr)
  
  save(corr_t2m,corr_tp,corr_par,file=paste0("./analysis/ERA_corr/",lag-1,".mon_corr.RData"))
  
  ### combine all and calculate the partial correlation
  climate_dat<-abind(eosdate,t2m_dat,tp_dat,par_dat,along=1)
  pcorr_all<-apply(climate_dat,c(2,3),get_pcorr)
  pcorr_t2m<-pcorr_all[1:2,,]
  pcorr_tp<-pcorr_all[3:4,,]
  pcorr_par<-pcorr_all[5:6,,]
  
  save(pcorr_t2m,pcorr_tp,pcorr_par,file=paste0("./analysis/ERA_corr/",lag-1,".mon_pcorr.RData"))
}

get_correlation(uid)


