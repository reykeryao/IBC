rm(list=ls())
gc()
library(tidyverse)
library(tidyr)
library(readxl)
readExcel<-function(filename, tibble=FALSE) {
  require(readxl)
  sheets = readxl::excel_sheets(filename)
  x = lapply(sheets, function(.) {
    readxl::read_excel(filename, sheet=., na=c("", "NA", "NaN"))
  })
  if (!tibble) {x = lapply(x, as.data.frame)}
  names(x) = sheets
  return(x)
}

setwd("/stor/work/Lambowitz/yaojun/Work/IBC_test/BigWig/DESeq2_objexts/")
Sheet<-readExcel("deseq_sheet.xlsx")
## ----------------FFPE---------------
## FFPE IBC vs nonIBC

dat<-read.delim("../../../cleaned.combined_with_Rep_miRNA_processed.counts")
dat<-dat[dat$Type%in%c("sncRNA","miRNA"),]
dat$PlasmaH_rep<-apply(dat[,69:81],1,function(x){sum(x>0)>5})
dat$PlasmaB_rep<-apply(dat[,63:68],1,function(x){sum(x>0)>2})
dat$PlasmaI_rep<-apply(dat[,53:62],1,function(x){sum(x>0)>4})
dat<-dat[,c(1,82:86)]

dds<-readRDS("../../../BigWig/DESeq2_objexts/sncrna_plasma_deseq_object.rds")
res1<-data.frame(results(dds,contrast=c("Disease","IBC","Healthy")))
res2<-data.frame(results(dds,contrast=c("Disease","IBC","nonIBC")))
norm_counts<-data.frame(counts(dds,normalized=T))
dat<-cbind(dat,norm_counts)
res1<-cbind(res1,dat)
res2<-cbind(res2,dat)
res1<-res1[complete.cases(res1),]
res2<-res2[complete.cases(res2),]
### IBC vs Healthy complete case and significant
res3<-res1[abs(res1$log2FoldChange)>=2 & res1$Type=="miRNA",]
### IBC vs nonIBC complete case and significant
res4<-res2[abs(res2$log2FoldChange)>=2& res2$Type=="miRNA",]
### IBC vs Healthy and nonIBC
res5<-merge(res3,res4[,1:7],by="ID")
res5<-res5[res5$log2FoldChange.x*res5$log2FoldChange.y>0,]
### up in IBC
res6<-res5[res5$log2FoldChange.x>0 & res5$PlasmaI_rep,]
res7<-res5[res5$log2FoldChange.x<0 & res5$PlasmaH_rep & res5$PlasmaB_rep,]
