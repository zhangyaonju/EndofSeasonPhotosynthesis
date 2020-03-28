### comparison between GRACE and SMAP 
## get the grace data
library(raster)
library(ncdf4)
library(abind)

grace_f<-list.files("~/Documents/Data/GRACE/",pattern=".nc",full.names = T)

grace_dat<-array(NA,dim=c(360,180,12*4))

for (i in 1:length(grace_f)){
  # grace time
  obs_t<-substr(basename(grace_f[i]),15,21)
  obs_n<-(as.numeric(substr(obs_t,1,4))-2014)*12+round(as.numeric(substr(obs_t,5,7))/30.5)
  print(paste(obs_t,obs_n))
  ncin<-nc_open(grace_f[i])
  grace_tws<-ncvar_get(ncin,varid="lwe_thickness")
  grace_dat[,,obs_n]<-grace_tws
}
#shift the longitude
grace_lwe<-abind(grace_dat[181:360,,],grace_dat[1:180,,],along=1)

save(grace_lwe,file='~/Documents/Data/GRACE/GRACE_2014_2017.RData')
load('~/Documents/Data/GRACE/GRACE_2014_2017.RData')

#### get the smap data for analysis
smap_f<-list.files('~/Documents/Data/SMAP/',full.names = T)
smap_dat<-array(NA,dim=c(360,180,12*4))
for (i in 1:3){
  # smap time
  obs_t<-substr(basename(smap_f[i]),1,4)
  ncin<-nc_open(smap_f[i])
  smap_sm<-ncvar_get(ncin,varid="SM_AM")
  for (j in 1:12){
    smap_r<-raster(smap_sm[,,j])
    smap_low<-aggregate(smap_r,2)
    smap_mat<-t(apply(as.matrix(smap_low),2,rev))
    smap_dat[,,i*12+j]<-smap_mat
  }
}
#shift the longitude


#comparison ## 1. regression in time
get_correlation<-function(vec){
  n_obs<-length(vec)/2
  if(sum(!is.na(vec[1:n_obs]))<5|sum(!is.na(vec[(n_obs+1):(2*n_obs)]))<5){
    return(c(NA,NA))
  }
  cort<-cor.test(vec[1:n_obs],vec[(n_obs+1):(n_obs*2)])
  return(c(cort$estimate,cort$p.value))
}

dat<-abind(grace_lwe,smap_dat,along=3)
cor_ana<-apply(dat,c(1,2),get_correlation)

dim(cor_ana)

image(cor_ana[1,,],zlim=c(-1,1))
### plot the graphs.
load('~/Dropbox/YAOZHANG/code/R-code/tools/RCOLOR.RData')
cor_ef<-raster(t(cor_ana[1,,180:1]))
extent(cor_ef)<-c(-180,180,-90,90)
pval<-raster(t(cor_ana[2,,180:1]<0.1))
extent(pval)<-c(-180,180,-90,90)


shp<-shapefile("~/Documents/Data/GIS_data/global/land_no_greenland_antarctica/ne_110m_land.shp")
pdf("~/Dropbox/YAOZHANG/paper/=2019_EOS_SM/PNAS_r1/figures/SMAP_GRACE_correlation.pdf",width=8,height=4)
par(oma=c(1,1,1,1),mar=c(2,2,1,5))
par(fig=c(0,1,0,1))
image(cor_ef,zlim=c(-1,1),col=blue_orange,ylim=c(-55,85),axes=F)
axis(1,tck=-0.02)
axis(2,las=2,tck=-0.02)
box()
lines(shp,col='grey70',lw=0.5)

plot(cor_ef, legend.only=TRUE, col=blue_orange,
     horizontal=F,zlim=c(-1,1),
     legend.width=1, legend.shrink=0.7,
     axis.args=list(at=c(-2:2)/2,
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))
dev.off()

### get the anomaly for the four seasons.
grace_lwe<-grace_lwe[,,13:48]
dim(grace_lwe)<-c(360,180,12,3)
grace_mean<-apply(grace_lwe,c(1,2),mean,na.rm=T)
mean_spring_grace<-apply(grace_lwe[,,3:5,],c(1,2),mean,na.rm=T)
mean_summer_grace<-apply(grace_lwe[,,6:8,],c(1,2),mean,na.rm=T)
mean_autumn_grace<-apply(grace_lwe[,,9:11,],c(1,2),mean,na.rm=T)
mean_winter_grace<-apply(grace_lwe[,,c(1,2,12),],c(1,2),mean,na.rm=T)

smap<-smap_dat[,,1:36]
smap_mean<-apply(smap,c(1,2),mean,na.rm=T)
dim(smap)<-c(360,180,12,3)
mean_spring_smap<-apply(smap[,,3:5,],c(1,2),mean,na.rm=T)
mean_summer_smap<-apply(smap[,,6:8,],c(1,2),mean,na.rm=T)
mean_autumn_smap<-apply(smap[,,9:11,],c(1,2),mean,na.rm=T)
mean_winter_smap<-apply(smap[,,c(1,2,12),],c(1,2),mean,na.rm=T)

### graph the comparison
graph_matrix<-function(dat,rang,colramp,msk){
  dat[dat<rang[1]]<-rang[1]
  dat[dat>rang[2]]<-rang[2]
  rdat<-raster(t(dat[,180:1]))
  extent(rdat)<-c(-180,180,-90,90)
  
  rdat<-mask(rdat,msk)
  
  image(rdat,zlim=rang,col=colramp,ylim=c(-55,85),axes=F)
  axis(1,tck=-0.02)
  axis(2,las=2,tck=-0.02)
  box()
  lines(shp,col='grey70',lw=0.5)
}

grace_msk<-rasterize(shp,cor_ef)



pdf("~/Dropbox/YAOZHANG/paper/=2019_EOS_SM/PNAS_r1/figures/SMAP_GRACE_four_season.pdf",width=8,height=9)
par(mar=c(2,2,1,1),oma=c(3,1,1,1),mgp=c(3,0.3,0))
par(fig=c(0,0.5,0.75,1))
graph_matrix(mean_spring_grace-grace_mean, c(-0.25,0.25),ncc_red_cyan,msk=grace_msk)
mtext(side=3,line=0.5,"GRACE LWET (m)")
text(-180,-45,LETTERS[1],pos=4,font=2,cex=1.5)
par(fig=c(0,0.5,0.5,0.75),new=T)
graph_matrix(mean_summer_grace-grace_mean, c(-0.25,0.25),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[2],pos=4,font=2,cex=1.5)
par(fig=c(0,0.5,0.25,0.5),new=T)
graph_matrix(mean_autumn_grace-grace_mean, c(-0.25,0.25),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[3],pos=4,font=2,cex=1.5)
par(fig=c(0,0.5,0,0.25),new=T)
graph_matrix(mean_winter_grace-grace_mean, c(-0.25,0.25),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[4],pos=4,font=2,cex=1.5)

par(fig=c(0.5,1,0.75,1),new=T)
graph_matrix(mean_spring_smap-smap_mean, c(-0.1,0.1),ncc_red_cyan,msk=grace_msk)
mtext(side=3,line=0.5,"SMAP SSM (m3/m3)")
text(-180,-45,LETTERS[5],pos=4,font=2,cex=1.5)
par(fig=c(0.5,1,0.5,0.75),new=T)
graph_matrix(mean_summer_smap-smap_mean, c(-0.1,0.1),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[6],pos=4,font=2,cex=1.5)
par(fig=c(0.5,1,0.25,0.5),new=T)
graph_matrix(mean_autumn_smap-smap_mean, c(-0.1,0.1),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[7],pos=4,font=2,cex=1.5)
par(fig=c(0.5,1,0,0.25),new=T)
graph_matrix(mean_winter_smap-smap_mean, c(-0.1,0.1),ncc_red_cyan,msk=grace_msk)
text(-180,-45,LETTERS[8],pos=4,font=2,cex=1.5)

par(oma=c(1,1,1,1))
par(fig=c(0,0.5,0,0.25),new=T)
plot(cor_ef, legend.only=TRUE, col=blue_orange,
     horizontal=T,zlim=c(-0.25,0.25),
     legend.width=1, legend.shrink=0.9,
     axis.args=list(at=c(-3:3)/20,
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))
par(fig=c(0.5,1,0,0.25),new=T)
plot(cor_ef, legend.only=TRUE, col=blue_orange,
     horizontal=T,zlim=c(-0.1,0.1),
     legend.width=1, legend.shrink=0.9,
     axis.args=list(at=c(-2:2)/20,
                    mgp=c(3,0.2,0),tck=1,
                    cex.axis=1))
dev.off()
