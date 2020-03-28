### plot the future with precipitation limitation.


### for both RCP45 and RCP85

library(raster)
## load the color for the figures
#setwd("C:/Users/shazh/Dropbox/")
setwd("/Users/yzhang/Documents/")
load("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")
coast<-shapefile("./Dropbox/YAOZHANG/code/R-code/tools/ne110coast/ne_110m_coastline.shp")

load("./Project/EOS_SM/analysis/mean_pheno.RData")
veg_area<-raster(t(mean_eos!=-999)[360:1,])
extent(veg_area)<-c(-180,180,-90,90)

landarea_no_arctic<-raster("./Data/land/0.5deg.area.no_arctic.tif")
veg_land<-66258632
## for RCP45
load("./Project/EOS_SM/analysis/CMIP5_svm/rcp45.area_stat_Median.RData")
load("./Project/EOS_SM/analysis/CMIP5_svm/rcp45.change_maps_Median.RData")
pre_plot_rcp45<-calc(stack(ens_precip_lim),sum,na.rm=T)*veg_area
area_all_rcp45<-area_all
mean_rcp45<-area_all[,1]*100/veg_land
sd_rcp45<-apply(area_all[,2:24]/veg_land,1,sd)*100
load("./Project/EOS_SM/analysis/CMIP5_svm/rcp85.area_stat_Median.RData")
load("./Project/EOS_SM/analysis/CMIP5_svm/rcp85.change_maps_Median.RData")
pre_plot_rcp85<-calc(stack(ens_precip_lim),sum,na.rm=T)*veg_area
area_all_rcp85<-area_all
mean_rcp85<-area_all[,1]*100/veg_land
sd_rcp85<-apply(area_all[,2:24]/veg_land,1,sd)*100

lat_label<-expression(paste("50",degree,"S",sep=""),
                      paste("25",degree,"S",sep=""),
                      paste("0",degree,sep=""),
                      paste("25",degree,"N",sep=""),
                      paste("50",degree,"N",sep=""),
                      paste("75",degree,"N",sep=""))
lon_label<-expression(paste("180",degree,"W",sep=""),
                      paste("90",degree,"W",sep=""),
                      paste("0",degree,sep=""),
                      paste("90",degree,"E",sep=""),
                      paste("180",degree,"E",sep=""))



pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/=2019_EOS_SM/figures_PNAS/Fig4_svm.pdf",height=8,width=7)
par(oma=c(1,1,1,1),mgp=c(3,0.5,0))
par(fig=c(0,1,0.7,1),mar=c(0,4,1,1))
image(pre_plot_rcp45,zlim=c(0,4),col=c("grey",blue_orange[12:9]),axes=F,
      xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
#axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=3,LETTERS[1],cex=2,font=2,padj=-3.5,las=2)
box()
#par(fig=c(0,0.2,0.7,0.9),mar=c(1,4,0,1),new=T)
plot(pre_plot_rcp45, legend.only=TRUE, col=rep(rev(c(blue_orange[12:9],"grey")),each=10),
     horizontal=F,zlim=c(0,5),
     legend.width=3, legend.shrink=0.7,smallplot=c(0.13,0.15,0.1,0.6),
     axis.args=list(at=0:4+0.5,label=c("no limitation","2010s","2030s","2060s","2090s"),
                    mgp=c(3,0.2,0),tck=0,
                    cex.axis=1))
par(fig=c(0,1,0.4,0.7),mar=c(1,4,0,1),new=T)
image(pre_plot_rcp85,zlim=c(0,4),col=c("grey",blue_orange[12:9]),axes=F,
      xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=3,LETTERS[2],cex=2,font=2,padj=-3.5,las=2)
box()
par(fig=c(0,1,0,0.4),mar=c(1,4,2,1),new=T)
plot(NA,xlim=c(2006,2100),ylim=c(40,100),yaxs="i",xlab="",ylab="",axes=F)
axis(1,tck=-0.01)
axis(2,las=2,tck=-0.01)
mtext(side=2,line=2.1,"Precip. limited region (%)")
box()
lines(2006:2100,mean_rcp45,lwd=2,col="blue")
# polygon(c(2006:2100,2100:2006),c(mean_rcp45-sd_rcp45,rev(mean_rcp45+sd_rcp45)),
#         border = NA,col = adjustcolor("blue",alpha.f = 0.4))
for (i in 1:23){
  lines(2006:2100,area_all_rcp45[,1+i]/veg_land*100,lwd=0.5,col=adjustcolor("blue",alpha.f = 0.3))
}

lines(2006:2100,mean_rcp85,lwd=2,col="red")
# polygon(c(2006:2100,2100:2006),c(mean_rcp85-sd_rcp85,rev(mean_rcp85+sd_rcp85)),
#         border = NA,col = adjustcolor("red",alpha.f = 0.4))
for (i in 1:23){
  lines(2006:2100,area_all_rcp85[,1+i]/veg_land*100,lwd=0.5,col=adjustcolor("red",alpha.f = 0.3))
}
mtext(side=2,line=3,LETTERS[3],cex=2,font=2,padj=-4.5,las=2)

legend("topleft",c("RCP4.5","RCP8.5"),col=c("blue","red"),lwd=c(2,2),lty=c(1,1))
dev.off()

#plot(NA,xlim=c(2006,2100),ylim=c(20,80),yaxs="i",xlab="",ylab="",axes=F)


