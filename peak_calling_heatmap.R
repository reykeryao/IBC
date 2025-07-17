rm(list=ls())
gc()
library (DESeq2)
#library(genefilter)
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
library(readxl)
setwd("/Users/junyao/Documents/NGS/IBC/nc_RNA_paper/")
dat<-read.delim("../IBC_final/IBC_git/all_IBC_cleanedup.counts")
dat<-dat[c(1:33,37:65,69:89)]
colnames(dat)[34:36]<-paste0("PBMC_H1",1:3)
colnames(dat)[63:65]<-paste0("Plasma_H1",1:3)
dat$Type<-as.factor(dat$Type)
dat$Type2<-as.factor(dat$Type2)

peakxl<-data.frame(read_excel("../plasma_peakCall/macs2_high_conf_peak_new.xlsx",sheet = 1))
rownames(peakxl)<-peakxl$HeatID
peakxl<-peakxl[,c(29:57,14)]
colnames(peakxl)<-c(colnames(dat)[53:81],"Type")
sum_reads<-colSums(dat[,53:81])/1e6
peakxl[,1:29]<-t(t(peakxl[,1:29])/sum_reads)
peakxl[peakxl==0]<-2^ -5
peakxl[,1:29]<-log10(peakxl[,1:29])

#plasma heatmap
dat1<-read.delim("../IBC_final/IBC_git/Plasma.DESeq.norm.counts")
dat1<-dat1[,c(1:14,17:35)]
colnames(dat1)[1]<-"ID"
tRF<-read.delim("tRF/tRF.DEseq.norm.counts")
tRF<-tRF[,c(85,86,54:63,66:84)]
dat1[dat1==0]<-2^-5
tRF[tRF==0]<-2^-5
dat1[,5:33]<-log10(dat1[,5:33])
tRF[,3:31]<-log10(tRF[,3:31])
tRF$Type2=paste(tRF$ID,tRF$Type,sep="-")
  
pro_list<-c("CDC42EP1","TRBJ1_1","TRBJ1_6","MSH6",
            "TREML1","ARMCX6","GP9","MRFAP1",
            "AL358473.2", "AL031005.1", "C15orf54", 
            "AC007036.3", "AP001324.1")
set1<-dat1[match(pro_list,dat1$Name),]
rownames(set1)<-set1$Name
set1<-set1[,5:33]
set1$Type<-c(rep("Protein coding",8),rep("lncRNA",5))
tRF_list<-c("LeuTAG-tRH5","GlnCTG-tRH5","AlaAGC-tRF3", 
            "ArgTCG-tRF3", "AsnGTT-tRF3")
set2<-tRF[match(tRF_list,tRF$Type2),]
rownames(set2)<-set2$ID
set2<-set2[,c(3:31,2)]
colnames(set2)<-colnames(set1)
snc_list<-c("SCARNA4","SNORA2A","SNORA8",  "SNORA21",
            "SNORA44","SNORA64","SNORA72", "SNORD14C", 
            "SNORD32A", "SNORD95","SNORD118","MIRLET7G", 
            "MIR26B", "MIR30C1","MIR101", "MIR144", 
            "MIR150", "MIR199A2", "MIR339", "MIR422A",
            "MIR1273D", "MIR1302", "MIR4454","MIR4512"
            )


set3<-dat1[match(snc_list,dat1$Name),]
rownames(set3)<-set3$Name
set3<-set3[,5:33]
set3$Type<-c(rep("snoRNA",11),rep("miRNA",13))

names(set1)<-names(set2)<-names(set3)<-names(peakxl)
set1<-rbind(set1,set2,set3,peakxl)
set1<-set1[,c(1:13,14,17:18,15:16,19,21:22,29,20,23:28,30)]
set1_type<-c(rep("Healthy",13),rep(c("HR-","HR+"),each=3),rep("HR-",3),
             rep("HR+",7))

pdf("plasma_heat.pdf",height=15,width=8)
scol<-colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(50)
heatmap.2(as.matrix(set1[,1:29]),labRow = rownames(set1),dendrogram = "none",
          scale="none",margins = c(10, 20),Rowv = F,
          density.info = "none",trace="none",symm=F,symbreaks=F,keysize=1,Colv = F,symkey=F,
          col =scol)
dev.off()

peakxl<-peakxl[,c(1:13,14,17:18,15:16,19,21:22,29,20,23:28,30)]
pdf("plasma_heatB.pdf",height=10,width=8)
scol<-colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(50)
heatmap.2(as.matrix(peakxl[,1:29]),labRow = peakxl$Type,dendrogram = "none",
          scale="none",margins = c(10, 20),Rowv = F,
          density.info = "none",trace="none",symm=F,symbreaks=F,keysize=1,Colv = F,symkey=F,
          col =scol)
dev.off()

## get new DE list after removing H11-H13, and new mapping
setwd("/stor/work/Lambowitz/yaojun/Work/IBC_test/BigWig/")
readExcel = function(filename, tibble=FALSE) {
  require(readxl)
  sheets = readxl::excel_sheets(filename)
  x = lapply(sheets, function(.) {
    readxl::read_excel(filename, sheet=., na=c("", "NA", "NaN"))
  })
  if (!tibble) {x = lapply(x, as.data.frame)}
  names(x) = sheets
  return(x)
}
deseqResults = readExcel("deseq_comparison_results.xlsx")
Plasma1<-deseqResults$Plasma_IBC_H
Plasma2<-deseqResults$Plasma_IBC_nonIBC
pro1<-Plasma1[Plasma1$padj<0.001,]
pro1<-pro1[abs(pro1$log2FoldChange)>=2,]
pro2<-Plasma2[Plasma2$padj<0.001,]
pro2<-pro2[abs(pro2$log2FoldChange)>=2,]
pro<-merge(pro1,pro2,by=1:3,all=T)

'''
# output gene id list for IGV
write.table(pro$id[pro$Type2=="Protein coding"],"IGV/Plasma/protein.id",quote=F,row.names=F,col.names=F)
write.table(pro$id[pro$Type2=="Antisense"|pro$Type2=="lincRNA"|pro$Type2=="Other lncRNA"],"IGV/Plasma/lncRNA.id",quote=F,row.names=F,col.names=F)
write.table(pro$id[pro$Type2=="Antisense"|pro$Type2=="lincRNA"|pro$Type2=="Other lncRNA"|pro$Type2=="Pseudogene"],"IGV/Plasma/lncRNA.id",quote=F,row.names=F,col.names=F)
write.table(pro$id[pro$Type2=="scaRNA"|pro$Type2=="snoRNA"|pro$Type2=="snRNA"],"IGV/Plasma/sncRNA.id",quote=F,row.names=F,col.names=F)
'''

#relaxed miRNA standard
miRNA1<-Plasma1[abs(Plasma1$log2FoldChange)>=2 & 
                  Plasma1$Type2=="miRNA",]
miRNA2<-Plasma2[abs(Plasma2$log2FoldChange)>=2 & 
                  Plasma2$Type2=="miRNA",]
miRNA<-merge(miRNA1,miRNA2,by=1:3)
write.table(paste0(sub("_","-",miRNA$id),":"),
            "IGV/Plasma/miRNA.id",quote=F,row.names=F,col.names=F)
write.table(paste0(sub("_","-",pro$id[pro$Type2=="miRNA"]),":"),
            "IGV/Plasma/miRNA_sig.id",quote=F,row.names=F,col.names=F)

tmp<-deseqResults$FFPE_IBC_nonIBC
