rm(list=ls())
gc()
library(tidyverse)
library(tidyr)
library(openxlsx)
library(DESeq2)
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

dds<-readRDS("FFPE/sncrna_ffpe_deseq_object.rds")
Counts<-data.frame(counts(dds,normalized=T))
dds<-readRDS("PBMC/sncrna_pbmc_deseq_object.rds")
Counts<-cbind(Counts,data.frame(counts(dds,normalized=T)))
dds<-readRDS("Plasma/sncrna_plasma_deseq_object.rds")
Counts<-cbind(Counts,data.frame(counts(dds,normalized=T)))
SNAR<-Counts["SNAR",]
SNAR<-data.frame(t(SNAR))
SNAR$Type<-c(rep("FFPE:non-IBC",4),rep("FFPE:IBC",10),
               rep("Frozen tissue:Healthy",4),rep("Frozen tissue:non-IBC",4),
               rep("PBMC:Healthy",13),rep("PBMC:non-IBC",6),
               rep("PBMC:IBC",10),rep("Plasma:Healthy",13),
               rep("Plasma:non-IBC",6),rep("Plasma:IBC",10))
SNAR$Type<-factor(SNAR$Type,levels=c("Frozen tissue:Healthy","Frozen tissue:non-IBC",
                                     "FFPE:non-IBC","FFPE:IBC","PBMC:Healthy",
                                     "PBMC:non-IBC","PBMC:IBC","Plasma:Healthy",
                                     "Plasma:non-IBC","Plasma:IBC"))
## in "/stor/work/Lambowitz/yaojun/Work/IBC_test/BigWig/DESeq2_objexts"
pdf("SNAR_IBC.pdf")
par(mar=c(10,6,2,1))
boxplot(SNAR$SNAR~SNAR$Type , col="gray80",ylim=c(0,800), 
        main="SNAR",
        xaxt = "n",axes=F,ylab=NA , xlab=NA,outline=FALSE)
axis(1, at=1:10, labels = FALSE)
axis(2, at=seq(0,800,200), labels = seq(0,800,200),las=2)
axis(2, at=400,labels = "DESeq2 normalized counts",adj=0.5,line = 3,tick=FALSE)
stripchart(SNAR ~ Type,data = SNAR,
           method = "jitter",pch = 19,
           vertical = TRUE,add = TRUE)
axis(1, at=1:10, labels = levels(SNAR$Type),las=2)
dev.off()