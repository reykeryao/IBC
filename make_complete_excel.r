rm(list=ls())
gc()
library(tidyverse)
library(readxl)
library(DESeq2)
library(openxlsx)
library(tidyr)
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

setwd("/stor/work/Lambowitz/yaojun/Work/IBC_test/Document/IBC_final/IBC_git/")
DE_path="/stor/work/Lambowitz/yaojun/Work/IBC_test/BigWig/DESeq2_objexts/"
Ct_path="/stor/work/Lambowitz/yaojun/Work/IBC_test/"
#process dat, set FLAGS for dat
dat<-read.delim(paste0(Ct_path,"cleaned.combined_with_Rep_miRNA_processed.counts"))
#remove tRNA from dat
dat<-dat[grep("tRNA",dat$Type,invert = T),]
tRF<-read.delim(paste0(Ct_path,"tRNA_frag.counts"))
rownames(tRF)<-paste(tRF$ID,tRF$Frag,sep=":")
tRF$Type<-"tRNA"
tRF$Type2<-paste0("tRNA:",tRF$Frag)
tRF<-tRF[,-2]
dat<-rbind(dat,tRF)
dat<-dat[,c(1,20,22,21,23,11,12,19,10,13:18,6:9,2:5,40:52,34,37,38,35,36,39,25,26,33,24,
          27:32,69:81,63,66,67,64,65,68,54,55,62,53,56:61,82,83)]



dat$Rep_FFPEB<-apply(dat[,2:5],1,function(x){sum(x>0)>=2})
dat$Rep_FFPEB_HRN<-apply(dat[,2:3],1,function(x){sum(x>0)>=1})
dat$Rep_FFPEB_HRP<-apply(dat[,4:5],1,function(x){sum(x>0)>=1})
dat$Rep_FFPEI<-apply(dat[,6:15],1,function(x){sum(x>0)>=5})
dat$Rep_FFPEI_HRN<-apply(dat[,6:8],1,function(x){sum(x>0)>=2})
dat$Rep_FFPEI_HRP<-apply(dat[,9:15],1,function(x){sum(x>0)>=4})
dat$Rep_BCH<-apply(dat[,16:19],1,function(x){sum(x>0)>=2})
dat$Rep_BC<-apply(dat[,20:23],1,function(x){sum(x>0)>=2})
dat$Rep_PBMCH<-apply(dat[,24:36],1,function(x){sum(x>0)>=7})
dat$Rep_PBMCB<-apply(dat[,37:42],1,function(x){sum(x>0)>=3})
dat$Rep_PBMCB_HRN<-apply(dat[,37:39],1,function(x){sum(x>0)>=2})
dat$Rep_PBMCB_HRP<-apply(dat[,40:42],1,function(x){sum(x>0)>=2})
dat$Rep_PBMCI<-apply(dat[,43:52],1,function(x){sum(x>0)>=5})
dat$Rep_PBMCI_HRN<-apply(dat[,43:45],1,function(x){sum(x>0)>=2})
dat$Rep_PBMCI_HRP<-apply(dat[,46:52],1,function(x){sum(x>0)>=4})
dat$Rep_PlasmaH<-apply(dat[,53:65],1,function(x){sum(x>0)>=7})
dat$Rep_PlasmaB<-apply(dat[,66:71],1,function(x){sum(x>0)>=3})
dat$Rep_PlasmaB_HRN<-apply(dat[,66:68],1,function(x){sum(x>0)>=2})
dat$Rep_PlasmaB_HRP<-apply(dat[,69:71],1,function(x){sum(x>0)>=2})
dat$Rep_PlasmaI<-apply(dat[,72:81],1,function(x){sum(x>0)>=5})
dat$Rep_PlasmaI_HRN<-apply(dat[,72:74],1,function(x){sum(x>0)>=2})
dat$Rep_PlasmaI_HRP<-apply(dat[,75:81],1,function(x){sum(x>0)>=4})
dat<-dat[,c(1,82:105)]

# A Breast tissue. Create combined normalized counts
mRNAObj<-readRDS(paste0(DE_path,"FFPE/mrna_ffpe_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"FFPE/lnc_ffpe_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"FFPE/sncrna_ffpe_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"FFPE/repeats_ffpe_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"FFPE/trna_ffpe_deseq_object.rds"))
  
Counts<-data.frame(counts(mRNAObj,normalized=T))
Counts_tmp<-data.frame(counts(lncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(sncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(repObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(tRFObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
#combined with FLAGs
Counts<-cbind(Counts,dat[,1:11])
Counts<-Counts[,c(23:25,1:22,26:33)]
#create results
contrast_nameA<-c("F_nonIBC",rep("IBC",3),rep("nonIBC",2))
contrast_nameB<-c("Healthy","nonIBC",rep(c("Healthy","F_nonIBC"),2))
repro_col<-c(33,32,29,26,29,32,29,33,26,32,26,33)
for (i in 1:6){
  # get results
  contra<-c("Disease",contrast_nameA[i],contrast_nameB[i])
  res<-data.frame(results(mRNAObj,contrast=contra))
  res_tmp<-data.frame(results(lncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:25),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:33)]
  if (i==1){
    Sheet<-list(res)
  } else {
    Sheet<-c(Sheet,list(res))
  }
}
# FFPE HR
mRNAObj<-readRDS(paste0(DE_path,"FFPE/mrna_ffpe_HR_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"FFPE/lnc_ffpe_HR_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"FFPE/sncrna_ffpe_HR_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"FFPE/repeats_ffpe_HR_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"FFPE/trna_ffpe_HR_deseq_object.rds"))

comp <- list(c("HR_P_vs_N"),c("HR_P_vs_N","DiseaseIBC.HRP"),
             c("Disease_IBC_vs_nonIBC"),c("Disease_IBC_vs_nonIBC","DiseaseIBC.HRP"))
repro_col=c(28,27,31,30,30,27,31,28)
for (i in 1:4){
  # get results
  res<-data.frame(results(mRNAObj,contrast=comp[i]))
  res_tmp<-data.frame(results(lncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:25),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:33)]
  Sheet<-c(Sheet,list(res))
}

# B PBMC. Create combined normalized counts
mRNAObj<-readRDS(paste0(DE_path,"PBMC/mrna_pbmc_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"PBMC/lnc_pbmc_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"PBMC/sncrna_pbmc_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"PBMC/repeats_pbmc_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"PBMC/trna_pbmc_deseq_object.rds"))

Counts<-data.frame(counts(mRNAObj,normalized=T))
Counts_tmp<-data.frame(counts(lncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(sncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(repObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(tRFObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
#combined with FLAGs
Counts<-cbind(Counts,dat[,c(1:3,12:18)])
Counts<-Counts[,c(30:32,1:29,33:39)]
#create results
contrast_nameA<-c(rep("IBC",2),"nonIBC")
contrast_nameB<-c("Healthy","nonIBC","Healthy")
repro_col<-c(37,33,37,34,34,33)

for (i in 1:3){
  # get results
  contra<-c("Disease",contrast_nameA[i],contrast_nameB[i])
  res<-data.frame(results(mRNAObj,contrast=contra))
  res_tmp<-data.frame(results(lncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:32),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:40)]
  Sheet<-c(Sheet,list(res))
}
# PBMC HR
#add mature RNA, tRF, and repeat obj
mRNAObj<-readRDS(paste0(DE_path,"PBMC/mrna_pbmc_HR_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"PBMC/lnc_pbmc_HR_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"PBMC/sncrna_pbmc_HR_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"PBMC/repeats_pbmc_HR_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"PBMC/trna_pbmc_HR_deseq_object.rds"))
comp <- list(c("HR_P_vs_N"),c("HR_P_vs_N","DiseaseIBC.HRP"),
             c("Disease_IBC_vs_nonIBC"),c("Disease_IBC_vs_nonIBC","DiseaseIBC.HRP"))
repro_col=c(36,35,39,38,38,35,39,36)
for (i in 1:4){
  res<-data.frame(results(mRNAObj,contrast=comp[i]))
  res_tmp<-data.frame(results(lncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:32),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:40)]
  Sheet<-c(Sheet,list(res))
}

# C Plasma. Create combined normalized counts
mRNAObj<-readRDS(paste0(DE_path,"Plasma/mrna_plasma_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"Plasma/lnc_plasma_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"Plasma/sncrna_plasma_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"Plasma/repeats_plasma_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"Plasma/trna_plasma_deseq_object.rds"))

Counts<-data.frame(counts(mRNAObj,normalized=T))
Counts_tmp<-data.frame(counts(lncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(sncObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(repObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
Counts_tmp<-data.frame(counts(tRFObj,normalized=T))
Counts<-rbind(Counts,Counts_tmp)
#combined with FLAGs
Counts<-cbind(Counts,dat[,c(1:3,19:25)])
Counts<-Counts[,c(30:32,1:29,33:39)]
#create results
contrast_nameA<-c(rep("IBC",2),"nonIBC")
contrast_nameB<-c("Healthy","nonIBC","Healthy")
repro_col<-c(37,33,37,34,34,33)

for (i in 1:3){
  # get results
  contra<-c("Disease",contrast_nameA[i],contrast_nameB[i])
  res<-data.frame(results(mRNAObj,contrast=contra))
  res_tmp<-data.frame(results(lncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=contra))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:32),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:40)]
  Sheet<-c(Sheet,list(res))
}

# Plasma HR
mRNAObj<-readRDS(paste0(DE_path,"Plasma/mrna_plasma_HR_deseq_object.rds"))
lncObj<-readRDS(paste0(DE_path,"Plasma/lnc_plasma_HR_deseq_object.rds"))
sncObj<-readRDS(paste0(DE_path,"Plasma/sncrna_plasma_HR_deseq_object.rds"))
repObj<-readRDS(paste0(DE_path,"Plasma/repeats_plasma_HR_deseq_object.rds"))
tRFObj<-readRDS(paste0(DE_path,"Plasma/trna_plasma_HR_deseq_object.rds"))
comp <- list(c("HR_P_vs_N"),c("HR_P_vs_N","DiseaseIBC.HRP"),
             c("Disease_IBC_vs_nonIBC"),c("Disease_IBC_vs_nonIBC","DiseaseIBC.HRP"))
repro_col=c(36,35,39,38,38,35,39,36)
for (i in 1:4){
  # get results
  res<-data.frame(results(mRNAObj,contrast=comp[i]))
  res_tmp<-data.frame(results(lncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(sncObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(repObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res_tmp<-data.frame(results(tRFObj,contrast=comp[i]))
  res<-rbind(res,res_tmp)
  res<-cbind(res,Counts[,c(c(1:32),repro_col[2*i-1],repro_col[2*i])])
  res<-res[,c(7:9,1:6,10:40)]
  Sheet<-c(Sheet,list(res))
}
names(Sheet)<- c("Frozen_nonIBCvsFrozen_Healthy","FFPE_IBCvsFFPE_nonIBC",
                 "FFPE_IBCvsFrozen_Healthy","FFPE_IBCvsFrozen_nonIBC",
                 "FFPE_nonIBCvsFrozen_H","FFPE_nonIBCvsFrozen_nonIBC",
                 "FFPE_nonIBC|HRPvsHRN","FFPE_IBC|HRPvsHRN",
                 "HRN|FFPE_IBCvsFFPE_nonIBC","HRP(FFPE_IBCvsFFPE_nonIBC",
                 "PBMC_IBCvsPBMC_Healthy","PBMC_IBCvsPBMC_nonIBC",
                 "PBMC_nonIBCvsPBMC_Healthy",
                 "PBMC_nonIBC|HRPvsHRN","PBMC_IBC|HRPvsHRN",
                 "HRN|PBMC_IBCvsPBMC_nonIBC","HRP|PBMC_IBCvsPBMC_nonIBC",
                 "Plasma_IBCvsPlasma_Healthy","Plasma_IBCvsPlasma_nonIBC",
                 "Plasma_nonIBCvsPlasma_H",
                 "Plasma_nonIBC|HRPvsHRN","Plasma_IBC|HRPvsHRN",
                 "HRN|Plasma_IBCvsPLasma_nonIBC","HRP|Plasma_IBCvsPlasma_nonIBC")  

write.xlsx(Sheet, file = paste0(DE_path,"deseq_sheet.xlsx"))

