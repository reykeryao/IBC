rm(list=ls())
gc()
library (DESeq2)
library(genefilter)
library(gplots)
library(RColorBrewer)
library(shape)
library(ggfortify)
library(pheatmap)
library("BiocParallel")
#library("ComplexHeatmap")
library(circlize)
library(scales)
library(ggpubr)
library(VennDiagram)
register(MulticoreParam(4))
library(tidyr)
library(multipanelfigure)

#function for volcano plot
plot_volcano<-function(dat_non_sig,dat_sig,x_lim,y_lim,xlab_seq,
                       ylab_seq,title,type,cex=1){
  plot(dat_non_sig$log2FoldChange, -log10(dat_non_sig$padj),xlim=x_lim,ylim=y_lim,
       main=title, xlab=bquote(log[2](FoldChange)),xaxt="n",bty="n", yaxt="n",
       ylab=bquote("-"*log[10]*"(padj)"),,cex=cex,cex.lab=0.7,cex.main=0.7,
       pch=21, bg="black")
  axis(side = 1,at = xlab_seq,labels = xlab_seq,cex.axis=0.7)
  axis(side = 2,labels = ylab_seq,at = ylab_seq,las=1,cex.axis=0.7)
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-2, col="black", lty=4, lwd=2.0)
  abline(v=2, col="black", lty=4, lwd=2.0)
  abline(h=-log10(0.001), col="black", lty=3, lwd=2.0)
  if (type=="pro"){
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=21, bg="red", cex=cex)
  } else if (type=="snc") {
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=21, col="black",
           bg=scol[dat_sig$Type2], cex=cex,lwd=0.5)
    if(dim(dat_sig)[1]>0){
      legend("topleft",legend = levels(dat_sig$Type2),col=scol,bty="n",pch=20,cex=0.5)
    }
  }
}
mcol<-c("black","red")
scol<-brewer.pal(9, "Set1")
pcol<-brewer.pal(4,"Paired")
#simplified volcano for TIFF output
plot_volcanoS<-function(dat_non_sig,dat_sig,x_lim,y_lim,xlab_seq,ylab_seq,type,cex=0.5){
  plot(dat_non_sig$log2FoldChange, -log10(dat_non_sig$padj),xlim=x_lim,ylim=y_lim,
       main=NA, xlab=NA,xaxt="n",bty="n", yaxt="n",
       ylab=NA,col=mcol[1],
       pch=16, cex=cex)
  axis(side = 1,at = xlab_seq,labels = NA,tck=-0.02,lwd=0.5)
  axis(side = 2,labels = NA,at = ylab_seq,las=1,tck=-0.02,lwd=0.5)
  abline(v=0, col="black", lty=3, lwd=0.5)
  abline(v=-2, col="black", lty=4, lwd=0.5)
  abline(v=2, col="black", lty=4, lwd=0.5)
  abline(h=-log10(0.001), col="black", lty=3, lwd=0.5)
  if (type=="pro"){
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=16, col=mcol[2], cex=cex)
  } else if (type=="snc") {
    points(dat_sig$log2FoldChange,-log10(dat_sig$padj), pch=16, 
           bg=scol[dat_sig$Type2], cex=cex,lwd=0.2)
    if(dim(dat_sig)[1]>0){
      legend("topleft",legend = levels(dat_sig$Type2),col=scol,cex=0.5)
    }
  }
}

dat<-read.delim("../../../combined_repeats.counts")
dat$Type1[dat$Type1=="RC"]<-"Other"
dat$Type1[dat$Type1=="Retroposon"]<-"Other"
dat$Type1[dat$Type1=="Unknown"]<-"Other"
dat<-dat[c(-5060,-5066,-5084,-5290,-5292),]
rownames(dat)<-paste(dat$ID,dat$Type1,sep="@")
dat<-dat[c(1:3,8:11,4:7,22,24,23,25,13,14,21,12,15:20,
           42:54,36,39,40,37,38,41,27,28,35,26,29:34,
           71:83,65,68,69,66,67,70,56,57,64,55,58:63)]
dat$Type1<-as.factor(dat$Type1)
dat$Type2<-as.factor(dat$Type2)
coldata <- data.frame(c(rep("Healthy",4),rep("F_nonIBC",4),rep("nonIBC",4),rep("IBC",10),
                        rep(c(rep("Healthy",13),rep("nonIBC",6),rep("IBC",10)),2)),
                      c(rep("H",4),rep("NA",4),
                        rep(c("N","P"),each=2),
                        c(rep("N",3),rep("P",7)),
                        rep(c(rep("H",13),
                              c("N","N","N","P","P","P"),
                              c("N","N","N",rep("P",7))),2)),
                      c(rep("Forzen Tissue",8),rep("FFPE",14),rep("PBMC",29),rep("Plasma",29)),
                      colnames(dat)[4:83],
                      row.names=colnames(dat)[4:83])
colnames(coldata) <- c("Disease","HR","Tissue","Name")
coldata$Disease<-factor(coldata$Disease,levels=c("Healthy","F_nonIBC","nonIBC","IBC"))
coldata$HR<-factor(coldata$HR,levels=c("H","NA","N","P"))
#rep DE HR
# creat deseq objects for HR analysis using the same size factors, 
FFPE_HR <- DESeqDataSetFromMatrix(countData = dat[,12:25],
                                  colData = coldata[9:22,],design = ~Disease+HR+Disease:HR)
PBMC_HR <- DESeqDataSetFromMatrix(countData = dat[,39:54],
                                  colData = coldata[36:51,],design = ~Disease+HR+Disease:HR)
Plasma_HR <- DESeqDataSetFromMatrix(countData = dat[,68:83],
                                    colData = coldata[65:80,],design = ~Disease+HR+Disease:HR)
sizeFactors(FFPE_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/mrna_ffpe_HR_deseq_object.rds"))
sizeFactors(PBMC_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/mrna_pbmc_HR_deseq_object.rds"))
sizeFactors(Plasma_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/mrna_plasma_HR_deseq_object.rds"))
FFPE_HR <- DESeq(FFPE_HR,parallel = T)
PBMC_HR <- DESeq(PBMC_HR,parallel = T)
Plasma_HR <- DESeq(Plasma_HR,parallel = T)
saveRDS(FFPE_HR,"../../../BigWig/DESeq2_objexts/repeats_ffpe_HR_deseq_object.rds")
saveRDS(PBMC_HR,"../../../BigWig/DESeq2_objexts/repeats_pbmc_HR_deseq_object.rds")
saveRDS(Plasma_HR,"../../../BigWig/DESeq2_objexts/repeats_plasma_HR_deseq_object.rds")

#mature miRNA
dat<-read.delim("../../../mature_miRNA.counts")
rownames(dat)<-dat$ID
dat<-dat[c(1,6:9,2:5,20,22,21,23,11,12,19,10,13:18,
           40:52,34,37,38,35,36,39,25,26,33,24,27:32,
           69:81,63,66,67,64,65,68,54,55,62,53,56:61)]
# creat deseq objects for HR analysis using the same size factors, 
FFPE_HR <- DESeqDataSetFromMatrix(countData = dat[,10:23],
                                  colData = coldata[9:22,],design = ~Disease+HR+Disease:HR)
PBMC_HR <- DESeqDataSetFromMatrix(countData = dat[,37:52],
                                  colData = coldata[36:51,],design = ~Disease+HR+Disease:HR)
Plasma_HR <- DESeqDataSetFromMatrix(countData = dat[,66:81],
                                    colData = coldata[65:80,],design = ~Disease+HR+Disease:HR)
sizeFactors(FFPE_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_ffpe_HR_deseq_object.rds"))
sizeFactors(PBMC_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_pbmc_HR_deseq_object.rds"))
sizeFactors(Plasma_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_plasma_HR_deseq_object.rds"))
FFPE_HR <- DESeq(FFPE_HR,parallel = T)
PBMC_HR <- DESeq(PBMC_HR,parallel = T)
Plasma_HR <- DESeq(Plasma_HR,parallel = T)
saveRDS(FFPE_HR,"../../../BigWig/DESeq2_objexts/matureMiRNA_ffpe_HR_deseq_object.rds")
saveRDS(PBMC_HR,"../../../BigWig/DESeq2_objexts/matureMiRNA_pbmc_HR_deseq_object.rds")
saveRDS(Plasma_HR,"../../../BigWig/DESeq2_objexts/matureMiRNA_plasma_HR_deseq_object.rds")

##non-mature miRNA
#mature miRNA
dat<-read.delim("../../../combined_non_mature_miRNA.counts")
dat<-dat[c(1,6:9,2:5,20,22,21,23,11,12,19,10,13:18,
           40:52,34,37,38,35,36,39,25,26,33,24,27:32,
           69:81,63,66,67,64,65,68,54,55,62,53,56:61)]
# creat deseq objects for HR analysis using the same size factors, 
FFPE_HR <- DESeqDataSetFromMatrix(countData = dat[,10:23],
                                  colData = coldata[9:22,],design = ~Disease+HR+Disease:HR)
PBMC_HR <- DESeqDataSetFromMatrix(countData = dat[,37:52],
                                  colData = coldata[36:51,],design = ~Disease+HR+Disease:HR)
Plasma_HR <- DESeqDataSetFromMatrix(countData = dat[,66:81],
                                    colData = coldata[65:80,],design = ~Disease+HR+Disease:HR)
sizeFactors(FFPE_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_ffpe_HR_deseq_object.rds"))
sizeFactors(PBMC_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_pbmc_HR_deseq_object.rds"))
sizeFactors(Plasma_HR) <- sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/snc_plasma_HR_deseq_object.rds"))
FFPE_HR <- DESeq(FFPE_HR,parallel = T)
PBMC_HR <- DESeq(PBMC_HR,parallel = T)
Plasma_HR <- DESeq(Plasma_HR,parallel = T)
saveRDS(FFPE_HR,"../../../BigWig/DESeq2_objexts/non_matureMiRNA_ffpe_HR_deseq_object.rds")
saveRDS(PBMC_HR,"../../../BigWig/DESeq2_objexts/non_matureMiRNA_pbmc_HR_deseq_object.rds")
saveRDS(Plasma_HR,"../../../BigWig/DESeq2_objexts/non_matureMiRNA_plasma_HR_deseq_object.rds")


#tRF
dat<-read.delim("../../../tRNA_frag.counts")
rownames(dat)<-paste(dat$ID,dat$Frag,sep="@")
dat<-dat[,-2]
dat<-dat[c(1,6:9,2:5,20,22,21,23,11,12,19,10,13:18,
           40:52,34,37,38,35,36,39,25,26,33,24,27:32,
           69:81,63,66,67,64,65,68,54,55,62,53,56:61)]
# creat deseq objects for HR analysis using the same size factors, 
#sizefactor using tRF object
FFPE_HR <- DESeqDataSetFromMatrix(countData = dat[,10:23],
                                  colData = coldata[9:22,],design = ~Disease+HR+Disease:HR)
PBMC_HR <- DESeqDataSetFromMatrix(countData = dat[,37:52],
                                  colData = coldata[36:51,],design = ~Disease+HR+Disease:HR)
Plasma_HR <- DESeqDataSetFromMatrix(countData = dat[,66:81],
                                    colData = coldata[65:80,],design = ~Disease+HR+Disease:HR)
SF<-sizeFactors(readRDS("../../../BigWig/DESeq2_objexts/DESeq_tRF.rds"))
sizeFactors(FFPE_HR) <- SF[c(19,21,20,22,10,11,18,9,12:17)]
sizeFactors(PBMC_HR) <- SF[c(33,36,37,34,35,38,24,25,32,23,26:31)]
sizeFactors(Plasma_HR) <- SF[c(62,65,66,63,64,67,53,54,61,52,55:60)]
FFPE_HR <- DESeq(FFPE_HR,parallel = T)
PBMC_HR <- DESeq(PBMC_HR,parallel = T)
Plasma_HR <- DESeq(Plasma_HR,parallel = T)
saveRDS(FFPE_HR,"../../../BigWig/DESeq2_objexts/tRF_ffpe_HR_deseq_object.rds")
saveRDS(PBMC_HR,"../../../BigWig/DESeq2_objexts/tRF_pbmc_HR_deseq_object.rds")
saveRDS(Plasma_HR,"../../../BigWig/DESeq2_objexts/tRF_plasma_HR_deseq_object.rds")

#read in miRNA (mature) and repeat DESeq2 obj
miRNA<-readRDS("../../../BigWig/DESeq2_objexts/DESeq_mature_miRNA.rds")
rep<-readRDS("../../../BigWig/DESeq2_objexts/DESeq_Rep_by_mRNA.rds")
#FigS3 sup
pdf("S3_sup.pdf",height=7,width=6)
par(mfrow=c(2,2))
res<-results(miRNA,contrast = c("Type","FFPE IBC","FFPE non-IBC"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-5,25),c(0,30),seq(-5,25,5),seq(0,30,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)
tmp<-readRDS("../../../BigWig/DESeq2_objexts/matureMiRNA_ffpe_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-results(tmp,contrast = comp)
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-8,6),c(0,3),seq(-8,6,2),seq(0,3,1),"Mature miRNA","pro",.7)

res<-data.frame(data.frame(results(rep,contrast = c("Type","FFPE IBC","FFPE non-IBC"))))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-25,25),c(0,20),seq(-25,25,5),seq(0,20,10),"Repeat elements","snc",.7)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into="Name",sep = "@",extra = "drop")
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

tmp<-readRDS("../../../BigWig/DESeq2_objexts/repeats_ffpe_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-results(tmp,contrast = comp)
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-8,8),c(0,3),seq(-8,8,2),seq(0,3,1),"Repeat elements","snc",.7)
dev.off()

#FigS8 sup
pdf("S8_sup.pdf",height=7,width=12)
par(mfrow=c(2,4))
res<-results(miRNA,contrast = c("Type","PBMC IBC","PBMC Healthy"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-5,25),c(0,30),seq(-5,25,5),seq(0,30,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)

res<-results(miRNA,contrast = c("Type","PBMC IBC","PBMC non-IBC"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-5,25),c(0,30),seq(-5,25,5),seq(0,30,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)
tmp<-readRDS("../../../BigWig/DESeq2_objexts/matureMiRNA_pbmc_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-results(tmp,contrast = comp)
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-6,6),c(0,3),seq(-6,6,2),seq(0,3,1),"Mature miRNA","pro",.7)
res<-results(miRNA,contrast = c("Type","PBMC non-IBC","PBMC Healthy"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-20,5),c(0,20),seq(-20,5,5),seq(0,20,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)

res<-data.frame(results(rep,contrast = c("Type","PBMC IBC","PBMC Healthy")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-3,15),c(0,9),seq(-3,15,3),seq(0,9,3),"Repeat elements","snc",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

res<-data.frame(results(rep,contrast = c("Type","PBMC IBC","PBMC non-IBC")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-30,30),c(0,20),seq(-30,30,10),seq(0,20,5),"Repeat elements","snc",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

tmp<-readRDS("../../../BigWig/DESeq2_objexts/repeats_pbmc_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-data.frame(results(tmp,contrast = comp))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-6,6),c(0,3),seq(-6,6,2),seq(0,3,1),"Mature miRNA","pro",.7)

res<-data.frame(results(rep,contrast = c("Type","PBMC non-IBC","PBMC Healthy")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-20,20),c(0,15),seq(-20,20,5),seq(0,15,3),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)
dev.off()

#FigS12 sup
pdf("S12_sup.pdf",height=7,width=12)
par(mfrow=c(2,4))
res<-results(miRNA,contrast = c("Type","Plasma IBC","Plasma Healthy"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-10,40),c(0,40),seq(-10,40,10),seq(0,40,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)

res<-results(miRNA,contrast = c("Type","Plasma IBC","Plasma non-IBC"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-10,30),c(0,15),seq(-10,30,10),seq(0,15,3),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)
tmp<-readRDS("../../../BigWig/DESeq2_objexts/matureMiRNA_plasma_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-results(tmp,contrast = comp)
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-10,30),c(0,6),seq(-10,30,10),seq(0,6,3),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)

res<-results(miRNA,contrast = c("Type","Plasma non-IBC","Plasma Healthy"))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
plot_volcano(non_sig,sig,c(-25,25),c(0,20),seq(-25,25,5),seq(0,20,10),"Mature miRNA","pro",.7)
text(sig$log2FoldChange,-log10(sig$padj),labels = rownames(sig),pos=2,cex=0,3)


res<-data.frame(results(rep,contrast = c("Type","Plasma IBC","Plasma Healthy")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-40,40),c(0,150),seq(-40,40,10),seq(0,150,30),"Repeat elements","snc",.7)
#text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

res<-data.frame(results(rep,contrast = c("Type","Plasma IBC","Plasma non-IBC")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-40,40),c(0,100),seq(-40,40,10),seq(0,100,25),"Repeat elements","snc",.7)
#text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

tmp<-readRDS("../../../BigWig/DESeq2_objexts/repeats_plasma_HR_deseq_object.rds")
comp <- list(c("HR_P_vs_N","DiseaseIBC.HRP"))
res<-data.frame(results(tmp,contrast = comp))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-30,30),c(0,15),seq(-30,30,10),seq(0,15,3),"Repeat elements","snc",.7)
#text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)

res<-data.frame(results(rep,contrast = c("Type","Plasma non-IBC","Plasma Healthy")))
res<-res[complete.cases(res$padj),]
non_sig<-subset(res, padj>0.001 | abs(log2FoldChange)<2)
sig<-subset(res, padj<=0.001 & abs(log2FoldChange)>=2)
sig$Name<-rownames(sig)
sig<-separate(sig,Name,into=c("Name","Type2"),sep = "@",extra = "drop")
sig$Type2<-factor(sig$Type2)
plot_volcano(non_sig,sig,c(-30,30),c(0,80),seq(-30,30,10),seq(0,80,20),"Repeat elements","snc",.7)
#text(sig$log2FoldChange,-log10(sig$padj),labels = sig$Name,pos=2,cex=0,3)
dev.off()