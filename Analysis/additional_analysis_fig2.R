##### analyze the threshold for determining the water and temeprature limitation/
library(raster)
library(ncdf4)
#library(rpart)
library(boot)
library(e1071)
setwd("/Users/yzhang/")
load("./Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")

load("./Project/EOS_SM/analysis/max_corr_climate.RData")
#load("./Project/EOS_SM/analysis/pre_season_sm.RData")
load("./Project/EOS_SM/analysis/mean_annual_t.RData")
load("./Project/EOS_SM/analysis/mean_pheno.RData")
load("./Project/EOS_SM/analysis/mean_pre_eos_sm.RData")

mean_eos[,1:180]<-mean_eos[,1:180]+180
#sm_pre<-t(as.matrix(sm_pre30))[,360:1]
sm_pre<-mean_pre_eos_sm[2,,]
sm_pre[sm_pre>1]<-NA
#
ncin<-nc_open("./Project/EOS_SM/Global_mcd12c1_landcover1_majority.nc")
lc_data<-ncvar_get(ncin,"lc1")
nc_close(ncin)

lc_data1<- lc_data
lc_data1[lc_data>=1&lc_data<=5]<-1  #forest
lc_data1[lc_data>=6&lc_data<=8]<-2  # woody land
lc_data1[lc_data>=9&lc_data<=11|lc_data==14]<-3 # grass wet sav
lc_data1[lc_data==12]<-NA #cropland
lc_data1[lc_data1>=4|lc_data<1]<-NA
#lc_ras<-raster(lc_data)

max_t<-t(as.matrix(max_val_t))[,360:1]
max_p<-t(as.matrix(max_val_p))[,360:1]

all<-data.frame(as.vector(sm_pre),as.vector(mean_annual_T),as.vector(max_t),as.vector(max_p),
           as.factor(as.vector(max_t)>0),as.factor(as.vector(max_p)<0),
           as.factor((as.vector(max_t)>0)*2+(as.vector(max_p)>0)*1-1))

### get the threshold using the decision tree algorithm.
names(all)<-c("sm","mean_T","t_corr","p_corr","t_lim","p_lim","tp_lim")

# for the temperature 
all_temp<-all[complete.cases(all)&abs(all[,3])>0.41,c(1,2,5)]
all_prec<-all[complete.cases(all)&abs(all[,4])>0.41,c(1,2,6)]
all_tp<-all[complete.cases(all)&(all[,4]>0.41|all[,3]>0.41)&all[,7]!=-1&all[,7]!=2,c(1,2,7)]
all_tp[,3]<-as.factor(all_tp[,3]==1)
if(F){
# bootstrapping to calculate the mean, sd of threshold --------------------

temp_thresh<-array(NA,dim=c(2000,2))
prec_thresh<-array(NA,dim=c(2000,2))
tp_thresh<-array(NA,dim=c(2000,2))
### cost function
cost_fun<-function(xy,all){
  costv<-sum(((all[,1]>xy[1]&all[,2]<xy[2])*1-(as.numeric(all[,3])-1))^2,na.rm=T)
  return(costv)
}

for (i in 1:2000){
  print(i)
  boot.temp.sample = sample(dim(all_temp)[1],replace=T)
  optM<-optim(par=c(0.1,10),cost_fun,gr="L-BGGS-B",all=all_temp[boot.temp.sample,])
  temp_thresh[i,]<-optM$par
  boot.prec.sample = sample(dim(all_prec)[1],replace=T)
  optM<-optim(par=c(0.1,10),cost_fun,gr="L-BGGS-B",all=all_prec[boot.prec.sample,])
  prec_thresh[i,]<-optM$par
  boot.tp.sample = sample(dim(all_tp)[1],replace=T)
  optM<-optim(par=c(0.1,10),cost_fun,gr="L-BGGS-B",all=all_tp[boot.tp.sample,])
  tp_thresh[i,]<-optM$par
}

save(temp_thresh,prec_thresh,tp_thresh,file="./Project/EOS_SM/analysis/bootstrapping_threshold_sm.RData")
####### new bootstrapping analysis

sd_by_bootstrap<-function(dat){
  # quat<-quantile(dat,c(0.05,0.95))
  # dat<-dat[dat<quat[2]]
  # dat<-dat[dat>quat[1]]
  median.fun <- function(data, idx){
    med<-median(data[idx])
    return(med)
  }
  tryCatch({
    results <- boot(data=dat, statistic=median.fun, R=10000)
    ci_95<-boot.ci(results,conf=0.95, type="bca")
    return(ci_95$bca[4:5])
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  return(rep(NA,2))
}



#apply(temp_thresh,2,mean)
quantile(temp_thresh[,1],c(0.05,0.5,0.95),na.rm=T)
quantile(temp_thresh[,2],c(0.05,0.5,0.95),na.rm=T)

#apply(prec_thresh,2,mean,na.rm=T)
quantile(prec_thresh[,1],c(0.05,0.5,0.95),na.rm=T)
quantile(prec_thresh[,2],c(0.05,0.5,0.95),na.rm=T)

#apply(tp_thresh,2,mean,na.rm=T)
quantile(tp_thresh[,1],c(0.05,0.5,0.95),na.rm=T)
quantile(tp_thresh[,2],c(0.05,0.5,0.95),na.rm=T)
# tp_quantile
# 5%        50%        95% 
# -0.1837891  1.5424194  2.4693634 ????
# 5%      50%      95% 
# 10.11730 10.66717 11.21735 

}


# get the accruacy of forest grassland ------------------------------------
all_lc<-data.frame(as.vector(lc_data1),as.vector(max_t),as.vector(max_p),
                as.factor((as.vector(max_t)>0)*2+(as.vector(max_p)>0)*1-1))
tp_lc<-all_lc[complete.cases(all_lc)&(all_lc[,2]>0.41|all_lc[,3]>0.41)&all_lc[,4]!=-1&all_lc[,4]!=2,c(1,4)]
lc_pred<-tp_lc[,1]==3
caret::confusionMatrix(as.factor(lc_pred),as.factor(tp_lc[,2]==0))

### accu 0.7471  kappa: 0.4953

# use the SVM to get the decision tree ------------------------------------
coef_SVM<-as.data.frame(array(NA,dim=c(3,3)))
names(coef_SVM)<-c("beta0",'beta1','beta2')

svmfit_t = svm(t_lim~.,dat=all_temp,kernel="linear",cost=1,gamma=0.5,scale=FALSE)
beta = t(svmfit_t$coefs) %*% svmfit_t$SV
beta0 = svmfit_t$rho
coef_SVM[1,]<-c(beta0,beta)

plot(all_temp$sm,all_temp$mean_T,col=as.numeric(all_temp$t_lim)+2,cex=0.1)
abline(beta0/beta[2],-beta[1]/beta[2])
pred <- predict(svmfit_t,all_temp[,c(1,2)])
caret::confusionMatrix(pred,all_temp$t_lim)
#### acc 0.9064  Kappa : 0.7948 
#p
svmfit_p = svm(p_lim~.,dat=all_prec,kernel="linear",cost=10,scale=FALSE)
beta = t(svmfit_p$coefs) %*% svmfit_p$SV
beta0 = svmfit_p$rho
coef_SVM[2,]<-c(beta0,beta)
pred <- predict(svmfit_p,all_prec[,c(1,2)])
caret::confusionMatrix(pred,all_prec$p_lim)

#### Accuracy : 0.8538 Kappa : 0.6936  
# tp
svmfit_tp = svm(tp_lim~.,dat=all_tp,kernel="linear",cost=10,scale=FALSE)
beta = t(svmfit_tp$coefs) %*% svmfit_tp$SV
beta0 = svmfit_tp$rho
coef_SVM[3,]<-c(beta0,beta)

pred <- predict(svmfit_tp,all_tp[,c(1,2)])
caret::confusionMatrix(pred,all_tp$tp_lim)
#### Accuracy : 0.9125 Kappa : 0.8249 

save(coef_SVM,file="./Project/EOS_SM/analysis/SVM_decision.RData")

get_accustat<-function(dat,thresh){
  pred<-as.factor(dat[,1]>thresh[1]&dat[,2]<thresh[2])
  return(caret::confusionMatrix(pred,dat[,3]))
}

get_accustat(all_temp,c(0.10204089,11.70457)) ### Accuracy : 0.8824  kappa : 0.7425
get_accustat(all_prec,c(0.1122440,9.202018))  ##  Accuracy : 0.8013  Kappa : 0.599  
get_accustat(all_tp,c(0.1063477,10.66717))    ##  Accuracy : 0.8836  Kappa : 0.7668
# using the RPART decision tree -------------------------------------------

# for (i in 1:2000){
#   # boot.temp.sample = sample(dim(all_temp)[1],replace=T)
#   # fit_temp<-rpart(t_lim~mean_T+sm,method="class",data=all_temp[boot.temp.sample,],
#   #                 control=rpart.control(maxdepth = 2))
#   # temp_thresh[i,]<-fit_temp$splits[c(1,4),4]
#   print(i)
#   boot.prec.sample = sample(dim(all_prec)[1],replace=T)
#   fit_prec<-rpart(p_lim~mean_T+sm,method="class",data=all_prec,#[boot.prec.sample,],
#                   control=rpart.control(maxdepth = 2))
#   if (dim(fit_prec$splits)[1]!=5){
#     next
#   }
#   prec_thresh[i,]<-fit_prec$splits[c(1,4),4]
#   
#   boot.tp.sample = sample(dim(all_tp)[1],replace=T)
#   tp_fit<-rpart(tp_lim~mean_T+sm,method="class",data=all_tp[boot.tp.sample,],
#                 control=rpart.control(maxdepth = 2))
#   if (dim(tp_fit$splits)[1]!=5){
#     next
#   }
#   tp_thresh[i,]<-tp_fit$splits[c(1,4),4]
# }
