

library(raster)
library(ncdf4)
#setwd("C://Users/shazh/Dropbox/")
setwd("/Users/yzhang/Documents/")
load("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")

load("./Project/EOS_SM/analysis/max_corr_climate.RData")
#load("./Project/EOS_SM/analysis/pre_season_sm.RData") 
load("./Project/EOS_SM/analysis/mean_annual_t.RData")
load("./Project/EOS_SM/analysis/mean_pheno.RData")
load("./Project/EOS_SM/analysis/mean_pre_eos_sm.RData")
load("./Project/EOS_SM/analysis/SVM_decision.RData")
##EOS climate and SM together
### adjust EOS and shift south hemisphere 180 day
mean_eos[,1:180]<-mean_eos[,1:180]+180
sm_pre<-t(as.matrix(sm_pre30))[,360:1]

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

sm_pre<-mean_pre_eos_sm[2,,]

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/=2019_EOS_SM/figures_PNAS/Fig2_double_limit_SVM.pdf",height=13.5,width=6)

par(oma=c(1,1,0,1),mgp=c(3,0.5,0))
par(fig=c(0,0.7,0.66,0.93),mar=c(3,4,0,0))
all<-cbind(as.vector(lc_data1),as.vector(sm_pre),as.vector(mean_annual_T),as.vector(max_t),
           as.vector(mean_eos))
all_na<-all[complete.cases(all)&abs(all[,4])>0.41,]
col_ano<-adjustcolor(GMT_panoply[ceiling((all_na[,4]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-20,35),xlim=c(0,0.6),
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
text(0.888,7.5, labels = expression(italic(R)[paste("EOP,Ta",sep="")]), xpd = NA, srt = -90,cex=1.1)   
axis(1,at=0:6/10,tck = -0.02)
axis(2,at=(-2:3)*10,tck = -0.02,las=2)
box()

ymax=1.5

par(fig=c(0,0.7,0.93,1),mar=c(0,4,1,0),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,2]> 0.6|negativeall[,2]< 0,2]=NA
positiveall[positiveall[,2]> 0.6|positiveall[,2]< 0,2]=NA


plot(NA,xlim=c(0,0.6),ylim=c(0,ymax),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")

histpos<-hist(positiveall[,2],breaks=(0:60)*0.01,plot=F)
rect(xleft=histpos$mids-1.5/500, ybottom=0, xright=histpos$mids+1.5/500, ytop=histpos$counts/1000,
     col=adjustcolor(GMT_panoply[15],alpha.f=0.7),border="white")
histneg<-hist(negativeall[,2],breaks=(0:60)*0.01,plot=F)
rect(xleft=histneg$mids-1.5/500, ybottom=0, xright=histneg$mids+1.5/500, ytop=histneg$counts/1000,
     col=adjustcolor(GMT_panoply[2],alpha.f=0.7),border="white")


arrows(x0 = 0.11,x1 = 0.05, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
arrows(x0 = 0.13,x1 = 0.19, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
text(0.02,ymax*0.8,"dry",xpd = NA)
text(0.22,ymax*0.8,"wet",xpd = NA)
text(0.35,ymax*0.62,expression("RED:"),col="red",pos=2,xpd = NA,cex=0.78)
text(0.32,ymax*0.59,expression("higher T"[air] %=>%"delayed EOP"),pos=4,xpd = NA,cex=0.78)
text(0.35,ymax*0.85,expression("BLUE:"),col="blue",pos=2,xpd = NA,cex=0.78)
text(0.32,ymax*0.82,expression("higher T"[air] %=>%"advanced EOP"),pos=4,xpd = NA,cex=0.78)
####
par(fig=c(0.7,0.9,0.66,0.93),mar=c(3,0,1,2),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,3]>35|negativeall[,3]< -20,3]=NA
positiveall[positiveall[,3]>35|positiveall[,3]< -20,3]=NA

plot(NA,xlim=c(0,ymax),ylim=c(-20,35),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
histpos<-hist(positiveall[,3],breaks=-20:35,plot=F)
rect(xleft=0, ybottom=histpos$mids-0.3, xright=histpos$counts/1000, ytop=histpos$mids+0.3,
     col=adjustcolor(GMT_panoply[15],alpha.f=0.7),border="white")
histneg<-hist(negativeall[,3],breaks=-20:35,plot=F)
rect(xleft=0, ybottom=histneg$mids-0.3, xright=histneg$counts/1000, ytop=histneg$mids+0.3,
     col=adjustcolor(GMT_panoply[2],alpha.f=0.7),border="white")

arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 8,y1=-2,length = 0.1)
arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 12,y1=22,length = 0.1)
text(ymax*0.8,-9.33,"cold",xpd = NA,srt = -90)
text(ymax*0.8,29.33,"warm",xpd = NA,srt = -90)

par(fig=c(0.2,1,0.66,0.93),mar=c(3,4,1,3),new=T)
plot(max_val_p, legend.only=TRUE, col=GMT_panoply,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))




##### prec
par(oma=c(0,1,1,1))
par(fig=c(0,0.7,0.33,0.6),mar=c(3,4,0,0),new=T)
all<-cbind(as.vector(lc_data1),as.vector(sm_pre),as.vector(mean_annual_T),as.vector(max_p),
           as.vector(mean_eos))
all_na<-all[complete.cases(all)&abs(all[,4])>0.412,]
col_ano<-adjustcolor(GMT_panoply[ceiling((all_na[,4]+1)*8)],alpha.f = 0.6)
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-20,35),xlim=c(0,0.6),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
### add the threshold

abline(coef_SVM[2,1]/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3],col="seagreen3",lwd=2)
abline((coef_SVM[2,1] - 1)/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3], col="seagreen3",lty = 2)
abline((coef_SVM[2,1] + 1)/coef_SVM[2,3], -coef_SVM[2,2]/coef_SVM[2,3], col="seagreen3",lty = 2)

abline(coef_SVM[3,1]/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black",lwd=2)
abline((coef_SVM[3,1] - 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)
abline((coef_SVM[3,1] + 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)


mtext(side=1,line=1.8,expression(paste("Pre-EOP soil moisture (m"^3," m"^-3,")",sep="")))
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
text(0.888,7.5, labels = expression(italic(R)[paste("EOP,Prec.",sep="")]), xpd = NA, srt = -90,cex=1.1)   
mtext(side=2,line=3,LETTERS[2],cex=2,font=2,padj=-6,las=2)
axis(1,at=0:6/10,tck = -0.02)
axis(2,at=-2:3*10,tck = -0.02,las=2)
box()
par(fig=c(0,0.7,0.6,0.66),mar=c(0,4,1,0),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,2]> 0.6|negativeall[,2]< 0,2]=NA
positiveall[positiveall[,2]> 0.6|positiveall[,2]< 0,2]=NA

plot(NA,xlim=c(0,0.6),ylim=c(0,ymax),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
histpos<-hist(positiveall[,2],breaks=(0:60)*0.01,plot=F)
rect(xleft=histpos$mids-1.5/500, ybottom=0, xright=histpos$mids+1.5/500, ytop=histpos$counts/1000,
     col=adjustcolor(GMT_panoply[15],alpha.f=0.7),border="white")
histneg<-hist(negativeall[,2],breaks=(0:60)*0.01,plot=F)
rect(xleft=histneg$mids-1.5/500, ybottom=0, xright=histneg$mids+1.5/500, ytop=histneg$counts/1000,
     col=adjustcolor(GMT_panoply[2],alpha.f=0.7),border="white")

arrows(x0 = 0.11,x1 = 0.05, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
arrows(x0 = 0.13,x1 = 0.19, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
text(0.02,ymax*0.8,"dry",xpd = NA)
text(0.22,ymax*0.8,"wet",xpd = NA)
text(0.35,ymax*0.62,expression("RED:"),col="red",pos=2,xpd = NA,cex=0.78)
text(0.32,ymax*0.59,expression("higher Prec." %=>%"delayed EOP"),pos=4,xpd = NA,cex=0.78)
text(0.35,ymax*0.85,expression("BLUE:"),col="blue",pos=2,xpd = NA,cex=0.78)
text(0.32,ymax*0.82,expression("higher Prec." %=>%"advanced EOP"),pos=4,xpd = NA,cex=0.78)
####
par(fig=c(0.7,0.9,0.33,0.6),mar=c(3,0,1,2),new=T)
negativeall<-all_na[all_na[,4]<0,]
positiveall<-all_na[all_na[,4]>0,]
negativeall[negativeall[,3]>35|negativeall[,3]< -20,3]=NA
positiveall[positiveall[,3]>35|positiveall[,3]< -20,3]=NA

plot(NA,xlim=c(0,ymax),ylim=c(-20,35),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
histpos<-hist(positiveall[,3],breaks=-20:35,plot=F)
rect(xleft=0, ybottom=histpos$mids-0.3, xright=histpos$counts/1000, ytop=histpos$mids+0.3,
     col=adjustcolor(GMT_panoply[15],alpha.f=0.7),border="white")
histneg<-hist(negativeall[,3],breaks=-20:35,plot=F)
rect(xleft=0, ybottom=histneg$mids-0.3, xright=histneg$counts/1000, ytop=histneg$mids+0.3,
     col=adjustcolor(GMT_panoply[2],alpha.f=0.7),border="white")

arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 8,y1=-2,length = 0.1)
arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 12,y1=22,length = 0.1)
text(ymax*0.8,-9.33,"cold",xpd = NA,srt = -90)
text(ymax*0.8,29.33,"warm",xpd = NA,srt = -90)

par(fig=c(0.2,1,0.33,0.6),mar=c(3,4,1,3),new=T)
plot(max_val_p, legend.only=TRUE, col=GMT_panoply,horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(
       mgp=c(3,0.2,0),tck=0.3,
       cex.axis=0.8))
#dev.off()
#######

biome_col<-adjustcolor(c("darkgreen","blue","red"),alpha.f = 0.4)

par(fig=c(0,0.7,0,0.27),mar=c(3,4,0,0),new=T)
all<-cbind(as.vector(lc_data1),as.vector(sm_pre),as.vector(mean_annual_T),as.vector(max_p),
           as.vector(mean_eos))
all_na<-all[complete.cases(all),]#&abs(all[,4])>0.412,]
col_ano<-biome_col[all_na[,1]]
plot(all_na[,2],all_na[,3],col=col_ano,cex=0.3,pch=15,ylim=c(-20,35),xlim=c(0,0.6),
     xaxs="i",yaxs="i",ylab="",xlab="",axes=F)
mtext(side=1,line=1.8,expression(paste("Pre-EOP soil moisture (m"^3," m"^-3,")",sep="")))
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))
mtext(side=2,line=3,LETTERS[3],cex=2,font=2,padj=-6,las=2)
axis(1,at=0:6/10,tck = -0.02)
axis(2,at=(-2:3)*10,tck = -0.02,las=2)
box()

par(fig=c(0,0.7,0.27,0.33),mar=c(0,4,1,0),new=T)
forest<-all_na[all_na[,1]==1,]
woodland<-all_na[all_na[,1]==2,]
grassland<-all_na[all_na[,1]==3,]

forest[forest[,2]> 0.6|forest[,2]< 0,2]=NA
woodland[woodland[,2]> 0.6|woodland[,2]< 0,2]=NA
grassland[grassland[,2]> 0.6|grassland[,2]< 0,2]=NA


plot(NA,xlim=c(0,0.6),ylim=c(0,ymax),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")

histfor<-hist(forest[,2],breaks=(0:60)*0.01,plot=F)
histwoo<-hist(woodland[,2],breaks=(0:60)*0.01,plot=F)
histgra<-hist(grassland[,2],breaks=(0:60)*0.01,plot=F)

rect(xleft=histgra$mids-1.5/500, ybottom=0, xright=histgra$mids+1.5/500, ytop=histgra$counts/1000,
     col=adjustcolor("red",alpha.f = 0.7),border="white")
rect(xleft=histwoo$mids-1.5/500, ybottom=0, xright=histwoo$mids+1.5/500, ytop=histwoo$counts/1000,
     col=adjustcolor("blue",alpha.f = 0.7),border="white")
rect(xleft=histfor$mids-1.5/500, ybottom=0, xright=histfor$mids+1.5/500, ytop=histfor$counts/1000,
     col=adjustcolor("darkgreen",alpha.f = 0.7),border="white")

arrows(x0 = 0.11,x1 = 0.05, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
arrows(x0 = 0.13,x1 = 0.19, y0 = ymax*0.8,y1=ymax*0.8,length = 0.1)
text(0.02,ymax*0.8,"dry",xpd = NA)
text(0.22,ymax*0.8,"wet",xpd = NA)
####
par(fig=c(0.7,0.9,0,0.27),mar=c(3,0,1,2),new=T)

plot(NA,xlim=c(0,ymax),ylim=c(-20,35),axes=F,xaxs='i',yaxs="i",xlab="",ylab="")
forest<-all_na[all_na[,1]==1,]
woodland<-all_na[all_na[,1]==2,]
grassland<-all_na[all_na[,1]==3,]

forest[forest[,3]> 35|forest[,3]< -20,3]=NA
woodland[woodland[,3]> 35|woodland[,3]< -20,3]=NA
grassland[grassland[,3]> 35|grassland[,3]< -20,3]=NA

histfor<-hist(forest[,3],breaks=-20:35,plot=F)
histwoo<-hist(woodland[,3],breaks=-20:35,plot=F)
histgra<-hist(grassland[,3],breaks=-20:35,plot=F)
rect(xleft=0, ybottom=histgra$mids-0.3, xright=histgra$counts/1000, ytop=histgra$mids+0.3,
     col=adjustcolor("red",alpha.f = 0.7),border="white")
rect(xleft=0, ybottom=histwoo$mids-0.3, xright=histwoo$counts/1000, ytop=histwoo$mids+0.3,
     col=adjustcolor("blue",alpha.f = 0.7),border="white")
rect(xleft=0, ybottom=histfor$mids-0.3, xright=histfor$counts/1000, ytop=histfor$mids+0.3,
     col=adjustcolor("darkgreen",alpha.f = 0.7),border="white")



arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 8,y1=-2,length = 0.1)
arrows(x0 = ymax*0.8,x1 = ymax*0.8, y0 = 12,y1=22,length = 0.1)
text(ymax*0.8,-9.33,"cold",xpd = NA,srt = -90)
text(ymax*0.8,29.33,"warm",xpd = NA,srt = -90)


par(fig=c(0.2,1,0,0.27),mar=c(3,4,1,3),new=T)
plot(max_val_p, legend.only=TRUE, col=rep(c("darkgreen","blue","red"),each=20),horizontal=F,zlim=c(0.5,3.5),
     legend.width=1, legend.shrink=0.75,
     axis.args=list(at=1:3,label=c("Forest","Woodland","Grassland"),
                    mgp=c(3,0.2,0),tck=0.3,
                    cex.axis=0.8))

dev.off()





