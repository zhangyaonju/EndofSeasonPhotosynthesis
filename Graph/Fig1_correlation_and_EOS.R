
# graph the correlation of climate and EOS --------------------------------

library(raster)
## load the color for the figures
#setwd("C:/Users/shazh/Dropbox/")
setwd("/Users/yzhang/Documents/")
load("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData")
coast<-shapefile("/Users/yzhang/Dropbox/YAOZHANG/code/R-code/tools/ne110coast/ne_110m_coastline.shp")

if(F){ 
  #### this is the functions to get the maximum preseason length, only need to be run onces.
era_files<-list.files("./Project/EOS_SM/analysis/ERA_corr/",pattern="_corr",full.names = T)
rs_files<-list.files("./Project/EOS_SM/analysis/RS_corr/",pattern="_corr",full.names = T)

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

## ERA dataset
temp_corr<-array(NA,dim=c(720,360,4))
prec_corr<-array(NA,dim=c(720,360,4))

##remote sensing dataset
lst_corr<-array(NA,dim=c(720,360,4))
pre_corr<-array(NA,dim=c(720,360,4))
for (i in 1:4){
  load(era_files[i])
  temp_corr[,,i]<-corr_t2m[1,,]
  prec_corr[,,i]<-corr_tp[1,,]
  load(rs_files[i])
  lst_corr[,,i]<-corr_lst_mean[1,,]
  pre_corr[,,i]<-corr_tp[1,,]
}
to_raster<-function(dat){
  max_val<-raster(apply(dat[1,,],1,rev))
  extent(max_val)<-c(-180,180,-90,90)
  return(max_val)
}

max_mat_t<-apply(temp_corr,c(1,2),get_the_max)
max_mat_p<-apply(prec_corr,c(1,2),get_the_max)
max_mat_RSt<-apply(lst_corr,c(1,2),get_the_max)
max_mat_RSp<-apply(pre_corr,c(1,2),get_the_max)
max_val_t<-to_raster(max_mat_t)
max_val_p<-to_raster(max_mat_p)
max_val_RSt<-to_raster(max_mat_RSt)
max_val_RSp<-to_raster(max_mat_RSp)

save(max_val_p,max_val_t,max_val_RSp,max_val_RSt,
     max_mat_t,max_mat_p,max_mat_RSt,max_mat_RSp,
     file="./Project/EOS_SM/analysis/max_corr_climate.RData")
}


# The graph functions start here!!! ---------------------------------------
load("./Project/EOS_SM/analysis/mean_pre_eos_sm.RData")
load("./Project/EOS_SM/analysis/max_corr_climate.RData")
#load("./Project/EOS_SM/analysis/pre_season_sm.RData")
load("./Project/EOS_SM/analysis/mean_annual_t.RData")
mean_raster_T<-raster(t(mean_annual_T[,360:1]))
extent(mean_raster_T)<-c(-180,180,-90,90)

mean_eos_sm<-raster(t(mean_pre_eos_sm[2,,360:1]))
extent(mean_eos_sm)<-c(-180,180,-90,90)
mean_eos_sm[mean_eos_sm>0.6]<- 0.6
masked_T<-mask(mean_raster_T,mask=max_val_p)
#### need to convert the correlation into significance. at 0.1 0.05 0.01
# for n= 17, df =15, r = 0.412 0.482 0.606
ano_color1<-rep.int(GMT_panoply[c(1,3,5,8,9,12,14,16)],c(39,13,7,41,41,7,13,39))
ano_color2<-rep.int(GMT_panoply[c(1,3,5,8,9,12,14,16)],c(34,13,7,46,46,7,13,34))
##  for n=14, df =12, r= 0.458 0.532	0.661
crit1<-c(-0.606,-0.482,-0.412,0,0.412,0.482,0.606)
crit2<-c(-0.661,-0.532,-0.458,0,0.458,0.532,0.661)

get_corr_code<-function(dat,crit){
  dat[dat< crit[1]]<- crit[1]-0.02
  dat[dat> crit[1]& dat< crit[2]]<- crit[2]-0.02
  dat[dat> crit[2]& dat< crit[3]]<- crit[3]-0.02
  dat[dat> crit[3]& dat< crit[4]]<- crit[4]-0.02
  dat[dat> crit[4]& dat< crit[5]]<- crit[5]-0.02
  dat[dat> crit[5]& dat< crit[6]]<- crit[6]-0.02
  dat[dat> crit[6]& dat< crit[7]]<- crit[7]-0.02
  dat[dat> crit[7]]<- crit[7]+0.02
  return(dat)
}
max_val_t_fig<-get_corr_code(max_val_t,crit1)
max_val_p_fig<-get_corr_code(max_val_p,crit1)
max_val_RSt_fig<-get_corr_code(max_val_RSt,crit2)
max_val_RSp_fig<-get_corr_code(max_val_RSp,crit2)

##### stat the corr and preseason
get_pre_season<-function(cdat,crit){ ### crit should have 5 values
  stat_matrix<-array(NA,dim=c(4,4))  ### the row for the preseason, col for the significance
  dim(cdat)<-c(2,259200)
  cdat<-t(cdat)
  cdat<-cdat[complete.cases(cdat),]
  for (i in 1:4){
    temp<-cdat[cdat[,1]< crit[1+i]&cdat[,1]>=crit[i],2]
    stat_matrix[,i]<-as.numeric(table(temp))
  }
  colnames(stat_matrix)=c("A","B","C","D")
  rownames(stat_matrix)=c("15","30","60","90")
  return(stat_matrix)
}

stat_t<-get_pre_season(max_mat_t,crit=c(-1,-0.412,0,0.412,1))
stat_p<-get_pre_season(max_mat_p,c(-1,-0.412,0,0.412,1))
stat_RSt<-get_pre_season(max_mat_RSt,c(-1,-0.458,0,0.458,1))
stat_RSp<-get_pre_season(max_mat_RSp,c(-1,-0.458,0,0.458,1))

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
                      

pdf("/Users/yzhang/Dropbox/YAOZHANG/paper/=2019_EOS_SM/figures_PNAS/Fig1_correlation_corrected_Pre_season.pdf",height=9,width=12)
### a
par(mar=c(4,4,1.5,1),oma=c(1,1,0.5,1),mgp=c(3,0.3,0))
par(fig=c(0,0.5,0.66,1),mar=c(4,4,1,1))
image(max_val_t_fig,zlim=c(-1,1),col=ano_color1,axes=F,xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[1],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,expression(italic(R)["EOP,Ta"]~"(ERA-Interim)"))
box()
par(fig=c(0.07,0.17,0.745,0.86),mar=c(1,1,0,0),new=T)
barplot(stat_t/1000,col=c("grey80","grey60","grey50","grey40"),border="white",space=0.1,
        xlab="",axes=F,xaxt="n",ylim=c(0,12))
axis(2,at=0:3*5,las=2,tck=-0.02,cex.axis=0.7)
text( c(0.5,1.6,2.7,3.8)+0.1, par("usr")[1]-0.5, labels = c('-*',"-","+","+*"), cex.axis=0.7, xpd = TRUE)
#mtext(side=2,line=1,expression("counts ("%*%"1000)"),cex=0.8)
#axis(1,tck=-0.01)

### b
par(fig=c(0.5,1,0.66,1),mar=c(4,4,1,1),new=T)
image(max_val_p_fig,zlim=c(-1,1),col=ano_color1,axes=F,xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[2],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,expression(italic(R)["EOP,Prec."]~"(ERA-Interim)"))
box()
text(30,-92, expression(italic(R)),xpd=T)
par(fig=c(0.57,0.67,0.745,0.86),mar=c(1,1,0,0),new=T)
barplot(stat_p/1000,col=c("grey80","grey60","grey50","grey40"),border="white",space=0.1,
        xlab="",axes=F,xaxt="n",ylim=c(0,12))
axis(2,at=0:3*5,las=2,tck=-0.02,cex.axis=0.7)
text( c(0.5,1.6,2.7,3.8)+0.1, par("usr")[1]-0.5, labels = c('-*',"-","+","+*"), cex.axis=0.7, xpd = TRUE)


### c for LST mean and EOS
par(fig=c(0,0.5,0.33,0.66),mar=c(4,4,1,1),new=T)
image(max_val_RSt_fig,zlim=c(-1,1),col=ano_color2,axes=F,xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[3],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,expression(italic(R)["EOP,LST"]~"(MODIS)"))
box()
par(fig=c(0.07,0.17,0.415,0.53),mar=c(1,1,0,0),new=T)
barplot(stat_RSt/1000,col=c("grey80","grey60","grey50","grey40"),border="white",space=0.1,
        xlab="",axes=F,xaxt="n",ylim=c(0,12))
axis(2,at=0:3*5,las=2,tck=-0.02,cex.axis=0.7)
text( c(0.5,1.6,2.7,3.8)+0.1, par("usr")[1]-0.5, labels = c('-*',"-","+","+*"), cex.axis=0.7, xpd = TRUE)

### d for mswep prec and EOS
par(fig=c(0.5,1,0.33,0.66),mar=c(4,4,1,1),new=T)
image(max_val_RSp_fig,zlim=c(-1,1),col=ano_color2,axes=F,xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[4],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,expression(italic(R)["EOP,Prec."]~"(GPCP)"))
text(30,-92, expression(italic(R)),xpd=T)
box()
par(fig=c(0.57,0.67,0.415,0.53),mar=c(1,1,0,0),new=T)
barplot(stat_RSp/1000,col=c("grey80","grey60","grey50","grey40"),border="white",space=0.1,
        xlab="",axes=F,xaxt="n",ylim=c(0,12))
axis(2,at=0:3*5,las=2,tck=-0.02,cex.axis=0.7)
#axis(1,cex.axis=0.7)
text( c(0.5,1.6,2.7,3.8)+0.1, par("usr")[1]-0.5, labels = c('-*',"-","+","+*"), cex.axis=0.7, xpd = TRUE)
#text(1,col=NA,at=1:4-0.5,cex.axis=0.7,lab=c('-*',"-","+","+*"))

### e mean annual temperature
par(fig=c(0,0.5,0,0.33),mar=c(4,4,1,1),new=T)
image(masked_T,zlim=c(-20,30),col=rep(LST_color[1:10*25-10],each=10),axes=F,
      xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[5],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,"MAT")
text(130,-92, expression(paste(degree,"C",sep="")),xpd=T,pos=4)
box()

### f pre season soil water
par(fig=c(0.5,1,0,0.33),mar=c(4,4,1,1),new=T)
image(mean_eos_sm*100,zlim=c(0,60),col=rep(rev(blue_orange)[c(3,6,7:12)],each=10),axes=F,
      xlim=c(-180,180),ylim=c(-60,85),ylab="",xlab="")
lines(coast,col='grey30',lwd=0.5)
axis(1,at=c(-2:2*90),label=lon_label,tck=-0.01)
axis(2,at=c(-2:3*25),label=lat_label,las=2,tck=-0.01)
mtext(side=2,line=2.3,LETTERS[6],cex=2,font=2,padj=-3.5,las=2)
mtext(side=3,line=0,"Pre-EOP surface soil moisture")
text(130,-92, expression(paste("cm"^3," cm"^-3,sep="")),xpd=T,pos=4)
box()


par(fig=c(0,1,0.71,1),new=T,oma=c(0,2,2,2),mar=c(0,0,0,0))
plot(max_val_p_fig, legend.only=TRUE, col=ano_color1,horizontal=T,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.6,
     axis.args=list(at=round(crit1,2),
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))

par(fig=c(0,1,0.38,0.66),new=T,oma=c(0,2,2,2),mar=c(0,0,0,0))
plot(max_val_p_fig, legend.only=TRUE, col=ano_color2,horizontal=T,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.6,
     axis.args=list(at=round(crit2,2),
       mgp=c(3,0.2,0),tck=1,
       cex.axis=1))

par(fig=c(0.0,0.5,0.02,0.33),new=T,mar=c(4,4,1,1),oma=c(1,1,0.5,1))
plot(max_val_p_fig, legend.only=TRUE, col=rep(LST_color[1:10*25-10],each=10),horizontal=T,zlim=c(-20,30),
     legend.width=1, legend.shrink=0.7,
     axis.args=list(at=-2:3*10,
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))


par(fig=c(0.5,1,0.02,0.33),new=T,mar=c(4,4,1,1),oma=c(1,1,0.5,1))
plot(max_val_p_fig, legend.only=TRUE, col=rep(rev(blue_orange)[c(3,6,7:12)],each=10),
     horizontal=T,zlim=c(0,0.6),
     legend.width=1, legend.shrink=0.7,
     axis.args=list(at=0:4*0.15,
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))

dev.off()





