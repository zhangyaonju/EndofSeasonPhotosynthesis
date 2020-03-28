### get the separate mat and pssm
cli_file<-list.files("~/Documents/Project/EOS_SM/mean_annual_T/",full.names = T)
dat<-array(NA,dim=c(720,360,9))
MAT<<-array(NA,dim=c(720,360,3))
for (y in 2001:2009){
  cli_f<-cli_file[substr(basename(cli_file),1,4)==y]
  load(cli_f)
  dat[,,y-2000]<-mean_t
}
MAT[,,1]<-apply(dat,c(1,2),mean,na.rm=T)

dat<-array(NA,dim=c(720,360,9))
for (y in 2009:2017){
  cli_f<-cli_file[substr(basename(cli_file),1,4)==y]
  load(cli_f)
  dat[,,y-2008]<-mean_t
}
MAT[,,2]<-apply(dat,c(1,2),mean,na.rm=T)

dat<-array(NA,dim=c(720,360,17))
for (y in 2001:2017){
  cli_f<-cli_file[substr(basename(cli_file),1,4)==y]
  load(cli_f)
  dat[,,y-2000]<-mean_t
}
MAT[,,3]<-apply(dat,c(1,2),mean,na.rm=T)


#### this is the functions to get the maximum preseason length, only need to be run onces.


## this figure only shows the correlation with t and p, the pre-season length and the partial correlation
## is shown in the supporting information
## the remote sensing dataset based analysis is also shown in supporting information. 

get_the_max<-function(vec){
  if(sum(is.na(vec))>1){
    return(c(NA,NA))
  }
  ind<-which.max(abs(vec))
  val<-vec[ind]
  return(c(val,ind))
}

to_raster<-function(dat){
  max_val<-raster(apply(dat[1,,],1,rev))
  extent(max_val)<-c(-180,180,-90,90)
  return(max_val)
}

## ERA dataset
temp_corr<-array(NA,dim=c(720,360,4))
prec_corr<-array(NA,dim=c(720,360,4))

era_files<-list.files("~/Documents/Project/EOS_SM/analysis/ERA_corr_sep/",pattern="_corr.1.",full.names = T)
for (i in 1:4){
  load(era_files[i])
  temp_corr[,,i]<-corr_t2m[1,,]
  prec_corr[,,i]<-corr_tp[1,,]
}

max_mat_t1<-apply(temp_corr,c(1,2),get_the_max)
max_mat_p1<-apply(prec_corr,c(1,2),get_the_max)

max_val_t1<-to_raster(max_mat_t1)
max_val_p1<-to_raster(max_mat_p1)


era_files<-list.files("~/Documents/Project/EOS_SM/analysis/ERA_corr_sep/",pattern="_corr.2.",full.names = T)
for (i in 1:4){
  load(era_files[i])
  temp_corr[,,i]<-corr_t2m[1,,]
  prec_corr[,,i]<-corr_tp[1,,]
}

max_mat_t2<-apply(temp_corr,c(1,2),get_the_max)
max_mat_p2<-apply(prec_corr,c(1,2),get_the_max)

max_val_t2<-to_raster(max_mat_t2)
max_val_p2<-to_raster(max_mat_p2)

save(max_val_p1,max_val_t1,max_val_p2,max_val_t2,
     max_mat_t1,max_mat_p1,max_mat_t1,max_mat_p2,MAT,
     file="~/Documents/Project/EOS_SM/analysis/max_corr_climate_sep.RData")

### load the PSSM data
library(raster)
library("e1071")
load("~/Documents/Project/EOS_SM/analysis/max_corr_climate_sep.RData")
load("~/Documents/Project/EOS_SM/analysis/mean_pre_eos_sm.RData")
image(max_val_p1,zlim=c(-1,1))
image(max_val_p2,zlim=c(-1,1))
image(max_val_t1,zlim=c(-1,1))
image(max_val_t2,zlim=c(-1,1))

sm_pre<-mean_pre_eos_sm[2,,]


max_t1<-t(as.matrix(max_val_t1))[,360:1]
max_p1<-t(as.matrix(max_val_p1))[,360:1]

max_t2<-t(as.matrix(max_val_t2))[,360:1]
max_p2<-t(as.matrix(max_val_p2))[,360:1]

all<-data.frame(as.vector(sm_pre),as.vector(MAT[,,1]),as.vector(MAT[,,2]),
                as.vector(max_t1),as.vector(max_p1),
                as.vector(max_t2),as.vector(max_p2),
                as.factor((as.vector(max_t1)>0)*2+(as.vector(max_p1)>0)*1-1),
                as.factor((as.vector(max_t2)>0)*2+(as.vector(max_p2)>0)*1-1))
names(all)<-c("sm","mean_T1",'mean_T2',"t_corr1","p_corr1","t_corr2","p_corr2","tp_lim1","tp_lim2")

all_tp<-all[complete.cases(all)&(all[,4]>0.55|all[,5]>0.55|all[,6]>0.55|all[,7]>0.55)&all[,8]!=-1&all[,8]!=2&all[,9]!=-1&all[,9]!=2,c(1,2,3,8,9)]
all_tp[,4]<-as.factor(all_tp[,4]!=1)
all_tp[,5]<-as.factor(all_tp[,5]!=1)
#plot(all_tp[,1],all_tp[,2],col=all_tp[,3])

load("~/Documents/Project/EOS_SM/analysis/SVM_decision.RData")
pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,2]-coef_SVM[3,1]/coef_SVM[3,3]>0)
plot(all_tp[,1],all_tp[,2],col=pred)
caret::confusionMatrix(pred,all_tp$tp_lim1)   ### accuracy 0.8196, kappa  0.638    0.05 0.8579; 0.7158
# 8671  7976    7299                                     ### 0.1 0.9012  0.8016       0.01 8805       0.7596

#coef_SVM<-as.data.frame(array(NA,dim=c(3,3)))
#names(coef_SVM)<-c("beta0",'beta1','beta2')

svmfit_t1 = svm(tp_lim1~.,dat=all_tp[,c(1,2,4)],kernel="linear",cost=1,gamma=0.5,scale=FALSE)
beta = t(svmfit_t$coefs) %*% svmfit_t$SV
beta0 = svmfit_t$rho
coef_SVM[1,]<-c(beta0,beta)


pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,3]-coef_SVM[3,1]/coef_SVM[3,3]>0)
caret::confusionMatrix(pred,all_tp$tp_lim2)   ### accuracy 0.8624, kappa  0.722   0.05 0.8944 0.7883
 #8774   8179    781                                   ### 0.1 0.915 0.8298       0.01 0.9113 0.8222

svmfit_t2 = svm(tp_lim2~.,dat=all_tp[,c(1,3,5)],kernel="linear",cost=1,gamma=0.5,scale=FALSE)
beta = t(svmfit_t2$coefs) %*% svmfit_t2$SV
beta0 = svmfit_t2$rho
coef_SVM[2,]<-c(beta0,beta)

# cross validation scenerio
pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,2]-coef_SVM[3,1]/coef_SVM[3,3]>0)
caret::confusionMatrix(pred,all_tp$tp_lim2)    ## 0.8633  kappa 0.724     0.05 8953; 0.7902
# 8776  8191    7832                                         ### 0.1 0.9157 0.8311      0.01 0.9122 0.8241

pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,3]-coef_SVM[3,1]/coef_SVM[3,3]>0)
caret::confusionMatrix(pred,all_tp$tp_lim1)    ## 0.8186  kappa 0.636     0.05 8565    0.7131
 #8671  7964     728                                      ### 0.1 0.9004 0.8001      0.01 0.879   0.7566
##### start the graph.

save(coef_SVM,file="~/Documents/Project/EOS_SM/analysis/SVM_separate.RData")

load("~/Documents/Project/EOS_SM/analysis/SVM_separate.RData")
load("~/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")
pdf("~/Dropbox/YAOZHANG/paper/=2019_EOS_SM/PNAS_r1/figures/FigS_separate_peiod_SVM.pdf",height=8.9,width=6)

par(oma=c(1,1,1,1),mgp=c(3,0.5,0))
par(fig=c(0,0.8,0.5,1),mar=c(3,4,1,0))
all<-data.frame(as.vector(sm_pre),as.vector(MAT[,,1]),as.vector(MAT[,,2]),
                as.vector(max_t1),as.vector(max_p1),
                as.vector(max_t2),as.vector(max_p2),
                (as.vector(max_t1)>0)*2+(as.vector(max_p1)>0)*1-1,
                (as.vector(max_t2)>0)*2+(as.vector(max_p2)>0)*1-1)
names(all)<-c("sm","mean_T1",'mean_T2',"t_corr1","p_corr1","t_corr2","p_corr2","tp_lim1","tp_lim2")

all_tp<-all[complete.cases(all)&(all[,4]>0.55|all[,5]>0.55|all[,6]>0.55|all[,7]>0.55)&all[,8]!=-1&all[,8]!=2&all[,9]!=-1&all[,9]!=2,]

col_ano<-adjustcolor(GMT_panoply[(all_tp[,8])*8+4],alpha.f = 0.6)
plot(all_tp[,1],all_tp[,2],col=col_ano,cex=0.3,pch=15,ylim=c(-20,35),xlim=c(0,0.6),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
### add the threshold
abline(coef_SVM[1,1]/coef_SVM[1,3], -coef_SVM[1,2]/coef_SVM[1,3],col="seagreen3",lwd=2)
abline((coef_SVM[1,1] - 1)/coef_SVM[1,3], -coef_SVM[1,2]/coef_SVM[1,3],col="seagreen3", lty = 2)
abline((coef_SVM[1,1] + 1)/coef_SVM[1,3], -coef_SVM[1,2]/coef_SVM[1,3],col="seagreen3", lty = 2)

abline(coef_SVM[3,1]/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black",lwd=2)
abline((coef_SVM[3,1] - 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)
abline((coef_SVM[3,1] + 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)

mtext(side=1,line=1.8,expression(paste("Pre-EOP soil moisture (m"^3," m"^-3,")",sep="")))
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=3,LETTERS[1],cex=2,font=2,padj=-6,las=2)
#text(0.888,7.5, labels = expression(italic(R)[paste("EOP,Ta",sep="")]), xpd = NA, srt = -90,cex=1.1)   
axis(1,at=0:6/10,tck = -0.02)
axis(2,at=(-2:3)*10,tck = -0.02,las=2)
box()

par(fig=c(0,0.98,0.5,1),mar=c(3,4,1,3),new=T)
plot(max_val_p1, legend.only=TRUE, col=GMT_panoply[c(4,12)],horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=c(-.5,0.5), label=c("Prec. limited",'Ta limated'),
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))



par(fig=c(0,0.8,0,0.5),mar=c(3,4,1,0),new=T)

col_ano<-adjustcolor(GMT_panoply[(all_tp[,9])*8+4],alpha.f = 0.6)
plot(all_tp[,1],all_tp[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-20,35),xlim=c(0,0.6),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
### add the threshold
abline(coef_SVM[2,1]/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3],col="seagreen3",lwd=2)
abline((coef_SVM[2,1] - 1)/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3],col="seagreen3", lty = 2)
abline((coef_SVM[2,1] + 1)/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3],col="seagreen3", lty = 2)

abline(coef_SVM[3,1]/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black",lwd=2)
abline((coef_SVM[3,1] - 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)
abline((coef_SVM[3,1] + 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)

mtext(side=1,line=1.8,expression(paste("Pre-EOP soil moisture (m"^3," m"^-3,")",sep="")))
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=3,LETTERS[2],cex=2,font=2,padj=-6,las=2)
#text(0.888,7.5, labels = expression(italic(R)[paste("EOP,Ta",sep="")]), xpd = NA, srt = -90,cex=1.1)   
axis(1,at=0:6/10,tck = -0.02)
axis(2,at=(-2:3)*10,tck = -0.02,las=2)
box()

par(fig=c(0,0.98,0,0.5),mar=c(3,4,1,3),new=T)
plot(max_val_p1, legend.only=TRUE, col=GMT_panoply[c(4,12)],horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=c(-.5,0.5), label=c("Prec. limited",'Ta limated'),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=0.8))
dev.off()

pred <- as.factor(all_tp[,1]*coef_SVM[1,2]/coef_SVM[1,3]+all_tp[,2]-coef_SVM[1,1]/coef_SVM[1,3]>0)
caret::confusionMatrix(as.factor(all_tp$tp_lim1==0),pred)  # 0.728 0.8635
pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,2]-coef_SVM[3,1]/coef_SVM[3,3]>0)
caret::confusionMatrix(as.factor(all_tp$tp_lim1==0),pred)  # 0.7299 0.8655


pred <- as.factor(all_tp[,1]*coef_SVM[2,2]/coef_SVM[2,3]+all_tp[,3]-coef_SVM[2,1]/coef_SVM[2,3]>0)
caret::confusionMatrix(pred,as.factor(all_tp$tp_lim2==0)) # 0.7758 0.8882
pred <- as.factor(all_tp[,1]*coef_SVM[3,2]/coef_SVM[3,3]+all_tp[,3]-coef_SVM[3,1]/coef_SVM[3,3]>0)
caret::confusionMatrix(as.factor(all_tp$tp_lim2==0),pred) # 0.781 0.8908
