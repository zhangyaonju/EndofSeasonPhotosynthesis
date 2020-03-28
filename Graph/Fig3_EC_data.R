
### retrieve SM from the SMAP dataset.
## replace no obs with SMAP soil moisture.

setwd("/Users/yzhang/Documents/")
library(raster)
#library(e1071)
# graph the FLUXNET data --------------------------------------------------
load("./Project/EOS_SM/FLUX_analysis/eos_climate_all_50_sites.RData")
load("./Project/EOS_SM/FLUX_analysis/eos_climate_max_corr_50_sites.RData")
load("./Project/EOS_SM/analysis/SVM_decision.RData")
load("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")

interpolate<-function(c1,c2,c3){
  mid_col<-colorRampPalette(c(c1,c2,c3),interpolate = "linear",space="Lab")(10)
  return(mid_col)
}
interpolate2<-function(c1,c2){
  ramp1<-colorRampPalette(c(c1,"white"),interpolate = "linear",space="Lab")(7)[1:5]
  ramp2<-colorRampPalette(c("white",c2),interpolate = "linear",space="Lab")(7)[3:7]
  return(c(ramp1,ramp2))
}
yellow_red<-interpolate("#fff314","#fff4db","#ff024a")
blue_green<-interpolate("#300bff","#99faff","#5eff0d")

colorsquare<-array(NA,dim=c(10,10))
for (i in 1:10){
  colorsquare[i,]<-interpolate2(blue_green[i],yellow_red[i])
}
rgbsquare<-col2rgb(colorsquare)
dim(rgbsquare)<-c(3,10,10)
r<-raster(rgbsquare[1,,])
g<-raster(rgbsquare[2,,])
b<-raster(rgbsquare[3,,])
rgbstack<-stack(r,g,b)

maximum_corrlation$EOP_SM<-maximum_corrlation$Pre_SM/100
maximum_corrlation$EOP_SM[is.na(maximum_corrlation$EOP_SM)]<-
  maximum_corrlation$SM_SMAP_site_EOS[is.na(maximum_corrlation$EOP_SM)]


# create the site list ----------------------------------------------------

load("./Project/EOS_SM/FLUX_analysis/eos_date_comparison_CSIF.RData")
all_sites<-read.csv("./Data/FLUXNET_subset_tier1_DD/Site_information_tier1.csv",stringsAsFactors = F)
si_list<-all_sites[match(all_site_summary$site,all_sites$SITE_ID),]
si_list$years<-NA
si_list$validate<-NA
for (i in 1:50){
  sid<-preseason_var[[i]]$site[1]
  years_used<-preseason_var
  si_list$years[which(si_list$SITE_ID==sid)]<-paste(preseason_var[[i]]$year,collapse=", ")
  if (sid %in% all$sid){
    si_list$validate[which(si_list$SITE_ID==sid)]<-T
  }else{
    si_list$validate[which(si_list$SITE_ID==sid)]<-F
  }
}

write.csv(si_list,"./Dropbox/YAOZHANG/paper/2019_EOS_SM/figures/Table_S1.csv",row.names = F)
# # seprate data using SVM --------------------------------------------------
# 
# 
# sepa<-maximum_corrlation[,c(11,8)]
# sepa$tp_lim<-(maximum_corrlation$r_max_T_EOS>0)*2+(maximum_corrlation$r_max_P_EOS>0)*1
# sep_used<-sepa[sepa$tp_lim>=1&sepa$tp_lim<=2,]
# sep_used$tp_lim<-as.factor(sep_used$tp_lim)
# svmfit_tp = svm(tp_lim~.,dat=sep_used,kernel="linear",cost=10,scale=FALSE)
# beta = t(svmfit_tp$coefs) %*% svmfit_tp$SV
# beta0 = svmfit_tp$rho

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/=2019_EOS_SM/figures_PNAS/Fig3_EC_analysis.pdf",width=6,height=6)
par(mar=c(4,4,1,1),mgp=c(3,0.5,0),fig=c(0,1,0,1))
plot(NA,xlim=c(0,0.7),ylim=c(-5,30),axes=F,xlab="",ylab="",xaxs="i",yaxs="i")
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
mtext(side=1,line=2,expression(paste("Pre-EOP soil moisture (m"^3," m"^-3,")",sep="")))
mtext(side=2,line=2.2,expression(paste("Mean annual Temperature (",degree,"C)",sep="")))

maximum_corrlation$col_ano<-colorsquare[cbind(11-ceiling((maximum_corrlation$r_max_T_EOS+1)*5),
                                              ceiling((maximum_corrlation$r_max_P_EOS+1)*5))]
forest<-maximum_corrlation[maximum_corrlation$igbp %in% c("MF", "ENF", "DBF", "EBF"),]
forest<-forest[order(forest$p_max_P_EOS,decreasing = T),]
points(forest$SM_SMAP_site_EOS,forest$MAT,
       col=forest$col_ano,pch=15,cex=forest$years_used/10+0.5)

woodland<-maximum_corrlation[maximum_corrlation$igbp %in% c("OSH", "WSA"),]
points(woodland$SM_SMAP_site_EOS,woodland$MAT,
       col=woodland$col_ano,pch=16,cex=woodland$years_used/10+0.5)

grassland<-maximum_corrlation[maximum_corrlation$igbp %in% c("GRA", "SAV"),]
points(grassland$SM_SMAP_site_EOS,grassland$MAT,
       col=grassland$col_ano,pch=17,cex=grassland$years_used/10+0.5)

# xci<-c(0.1033739, 0.1056108, 0.1075379)
# yci<-c(10.30316, 10.76810, 11.47712)
# lines(c(xci[2],0.7),c(yci[2],yci[2]),lty=1,lwd=1.5,col='darkmagenta')
# polygon(c(xci[3],0.7,0.7,xci[3]),c(yci[1],yci[1],yci[3],yci[3]),
#         col=adjustcolor("darkmagenta",alpha.f = 0.5),border=NA)
# lines(c(xci[2],xci[2]),c(-20,yci[2]),lty=1,lwd=1.5,col='darkmagenta')
# polygon(c(xci[1],xci[1],xci[3],xci[3]),c(-20,yci[3],yci[3],-20),
#         col=adjustcolor("darkmagenta",alpha.f = 0.5),border=NA)

abline(coef_SVM[3,1]/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black",lwd=2)
abline((coef_SVM[3,1] - 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)
abline((coef_SVM[3,1] + 1)/coef_SVM[3,3], -coef_SVM[3,2]/coef_SVM[3,3],col="black", lty = 2)

# abline(beta0/beta[2], -beta[1]/beta[2],col="black",lwd=2)
# abline((beta0 - 1)/beta[2], -beta[1]/beta[2],col="black", lty = 2)
# abline((beta0 + 1)/beta[2], -beta[1]/beta[2],col="black", lty = 2)

polygon(c(0.2,0.7,3,3),c(30,12,50,50),col="white",border=NA)

legend(0.17,29,c("Grassland","Woodland","Forest"),col=rep("black",3),pch=c(2,1,0),bty="n")
legend(0.34,29,c("5","10","15"),pch=c(2,2,2),col=rep("black",3),pt.cex=c(1:3)*5/10+0.5,bty="n")

par(fig=c(0.75,0.95,0.75,0.95),mgp=c(3,0.2,0),new=T)
plotRGB(rgbstack,xaxs="i",yaxs="i")

arrows(x0 = 0.55,x1 = 0.95, y0 = 0.45, y1=0.05,length = 0.1)
arrows(x0 = 0.45,x1 = 0.05, y0 = 0.55, y1=0.95,length = 0.1)
text(0.8,0.3,"P limited",cex=0.7,srt=315)
text(0.3,0.8,"T limited",cex=0.7,srt=315)

axis(1,tck=-0.01,at=0:4/4,label=-2:2/2,cex.axis=0.7)
axis(2,las=2,tck=-0.01,at=0:4/4,label=-2:2/2,cex.axis=0.7)
mtext(side=1,line=1.5,expression(paste(italic(R),""["EOP,"],""[Prec.],sep="")))
mtext(side=2,line=1.5,expression(paste(italic(R),""["EOP,"],""[Ta],sep="")))
box()
dev.off()
