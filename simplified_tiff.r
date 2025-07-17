#function for volcano plot
plot_volcanoS<-function(dat_non_sig,dat_sig,x_lim,y_lim,xlab_seq,ylab_seq,type){
  plot(dat_non_sig$log2FoldChange, -log10(dat_non_sig$padj),xlim=x_lim,ylim=y_lim,
       main=NA, xlab=NA,xaxt="n",bty="n", yaxt="n",
       ylab=NA,col=mcol[1],
       pch=20, cex=0.2)
  axis(side = 1,at = xlab_seq,labels = NA,tck=-0.02,lwd=0.5)
  axis(side = 2,labels = NA,at = ylab_seq,las=1,tck=-0.02,lwd=0.5)
  abline(v=0, col="black", lty=3, lwd=0.5)
  abline(v=-2, col="black", lty=4, lwd=0.5)
  abline(v=2, col="black", lty=4, lwd=0.5)
  abline(h=-log10(0.001), col="black", lty=3, lwd=0.5)
  if (type=="pro"){
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=20, col=mcol[2], cex=0.2)
  } else if (type=="snc") {
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=21, 
           bg=scol[dat_sig$Type2], cex=0.2,lwd=0.2)
  }
}
mcol<-c("black","red")
scol<-brewer.pal(9, "Set1")
pcol<-brewer.pal(4,"Paired")

tiff("IBC_volcano_1.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#frozen  nonBC vs Healthy
res<-data.frame(results(FFPE_dds,contrast = c("Disease","FBC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-20,20),c(0,20),seq(-20,20,5),seq(0,20,5),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","FBC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,25),c(0,30),seq(-10,25,5),seq(0,30,10),"snc")
dev.off()

tiff("IBC_volcano_2.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE IBC vs BC
res<-data.frame(results(FFPE_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,15),seq(-10,10,5),seq(0,15,3),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-25,10),c(0,9),seq(-25,10,5),seq(0,9,3),"snc")
dev.off()

tiff("IBC_volcano_3.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#PBMC IBC vs Healthy
res<-data.frame(results(PBMC_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,60),seq(-10,10,5),seq(0,60,10),"pro")

res<-data.frame(results(PBMC_snc_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,60),seq(-10,10,5),seq(0,60,10),"snc")
dev.off()

tiff("IBC_volcano_4.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#PBMC IBC vs nonIBC
res<-data.frame(results(PBMC_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,60),seq(-10,10,5),seq(0,60,10),"pro")

res<-data.frame(results(PBMC_snc_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-15,15),c(0,60),seq(-15,15,5),seq(0,60,10),"snc")
dev.off()

tiff("IBC_volcano_5.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#PBMC nonIBC vs Healthy
res<-data.frame(results(PBMC_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-5,5),c(0,9),seq(-5,5,5),seq(0,9,3),"pro")

res<-data.frame(results(PBMC_snc_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-5,5),c(0,30),seq(-5,5,5),seq(0,30,10),"snc")
dev.off()

tiff("IBC_volcano_6.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#Plasma IBC vs Healthy
res<-data.frame(results(Plasma_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,50),seq(-30,30,10),seq(0,50,10),"pro")

res<-data.frame(results(Plasma_snc_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-20,40),c(0,30),seq(-20,40,10),seq(0,30,10),"snc")
dev.off()

tiff("IBC_volcano_7.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#Plasma IBC vs nonIBC
res<-data.frame(results(Plasma_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,40),seq(-30,30,10),seq(0,40,10),"pro")

res<-data.frame(results(Plasma_snc_dds,contrast = c("Disease","IBC","BC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,9),seq(-30,30,10),seq(0,9,3),"snc")
dev.off()

tiff("IBC_volcano_8.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#Plasma BC vs Healthy
res<-data.frame(results(Plasma_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,30),seq(-30,30,10),seq(0,30,10),"pro")

res<-data.frame(results(Plasma_snc_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,9),seq(-30,30,10),seq(0,9,3),"snc")
dev.off()

tiff("HR_nonIBC_volcano_1.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE HR+ vs HR-
res<-data.frame(results(FFPE_BC_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-25,25),c(0,40),seq(-25,25,5),seq(0,40,10),"pro")

res<-data.frame(results(FFPE_BC_snc_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,9),seq(-10,10,5),seq(0,9,3),"snc")
dev.off()

tiff("HR_IBC_volcano_1.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE HR+ vs HR-
res<-data.frame(results(FFPE_IBC_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,25),c(0,12),seq(-10,25,5),seq(0,12,3),"pro")

res<-data.frame(results(FFPE_IBC_snc_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"snc")
dev.off()

tiff("HR_nonIBC_volcano_2.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#PBMC HR+ vs HR-
res<-data.frame(results(PBMC_BC_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"pro")

res<-data.frame(results(PBMC_BC_snc_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"snc")
dev.off()

tiff("HR_IBC_volcano_2.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#PBMC HR+ vs HR-
res<-data.frame(results(PBMC_IBC_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"pro")

res<-data.frame(results(PBMC_IBC_snc_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"snc")
dev.off()

tiff("HR_nonIBC_volcano_3.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#Plasma HR+ vs HR-
res<-data.frame(results(Plasma_BC_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,25),c(0,9),seq(-10,25,5),seq(0,9,3),"pro")

res<-data.frame(results(Plasma_BC_snc_HR,contrast = c("HR","BCP","BCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,6),seq(-10,10,5),seq(0,6,3),"snc")
dev.off()

tiff("HR_IBC_volcano_3.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#Plasma HR+ vs HR-
res<-data.frame(results(Plasma_IBC_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-15,30),c(0,20),seq(-15,30,5),seq(0,20,5),"pro")

res<-data.frame(results(Plasma_IBC_snc_HR,contrast = c("HR","IBCP","IBCN"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,25),c(0,15),seq(-10,25,5),seq(0,15,3),"snc")
dev.off()

postscript("legend.eps")
plot_volcano(non_sig,sig,c(-10,25),c(0,15),seq(-10,25,5),seq(0,15,3),
             "sncRNA (Plasma IBC: HR+ vs HR-)","snc",-10,15)
dev.off()

tiff("tissue_volcano_1.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE  nonBC vs Forzen nonBC
res<-data.frame(results(FFPE_dds,contrast = c("Disease","BC","FBC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,30),c(0,9),seq(-30,30,10),seq(0,9,3),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","BC","FBC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,10),c(0,40),seq(-10,10,5),seq(0,40,10),"snc")
dev.off()

tiff("tissue_volcano_2.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE  nonBC vs Forzen nonBC
res<-data.frame(results(FFPE_dds,contrast = c("Disease","IBC","FBC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,20),c(0,20),seq(-10,20,10),seq(0,20,5),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","IBC","FBC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-30,10),c(0,75),seq(-30,10,10),seq(0,75,15),"snc")
dev.off()

tiff("tissue_volcano_3.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE  nonBC vs Forzen nonBC
res<-data.frame(results(FFPE_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,30),c(0,12),seq(-10,30,10),seq(0,12,3),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","IBC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-10,30),c(0,30),seq(-10,30,10),seq(0,30,10),"snc")
dev.off()

tiff("tissue_volcano_4.tiff",units = "px",res = 300,height=500,width=800)
par(mfrow=c(1,2),mar = c(1, 1, 0.1, 0.1))
#FFPE  nonBC vs Forzen nonBC
res<-data.frame(results(FFPE_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-20,30),c(0,15),seq(-20,30,10),seq(0,15,5),"pro")

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","BC","Healthy"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcanoS(non_sig,sig,c(-20,30),c(0,20),seq(-20,30,10),seq(0,20,10),"snc")
dev.off()


#pheatmap

res<-data.frame(results(FFPE_dds,contrast = c("Disease","IBC","nonIBC"),alpha = 0.05))
res<-cbind(res,pro_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)

res<-data.frame(results(FFPE_snc_dds,contrast = c("Disease","IBC","nonIBC"),alpha = 0.05))
res<-cbind(res,snc_dat[,c(1,89)])
res<-res[complete.cases(res),]
res$log10padj<-log10(1/res$padj)
sig<-rbind(sig,subset(res, padj<=0.001 & abs(log2FoldChange)>=2))

mat<-rbind(assay(FFPE_vst),assay(FFPE_snc_vst))
mat<-merge(mat,sig,by=0)
row.names(mat)<-mat$Row.names
mat[,2:23]<-mat[,2:23]-rowMeans(mat[,2:23])
anno<-data.frame(row.names=colnames(mat)[10:23],
                 "Disease"=c(rep("nonIBC",4),rep("IBC",10)),
                 "HR"=c("N","P","N","P","P","N","N",rep("P",6),"N"))
heat1<-pheatmap(as.matrix(mat[,c(10,12,11,13,15,16,23,14,17:22)]), annotation_col = anno,cluster_rows = F,
         cluster_cols = F,labels_row = mat$Name, show_colnames = F)



figure <- multi_panel_figure(width = 90, height = 60, columns = 3,rows=1) %>%
  fill_panel("FFPE_IBC_nonIBC.tiff",column=1:2) %>%
  fill_panel(heat1,column=3)


