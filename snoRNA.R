# map name of comparason to res or HR
# 01. Frozen tissue; nonIBC vs Healthy (4:4)
# 02 FFPE; IBC vs nonIBC (10:4)
# 03 FFPE IBC vs Frozen Healthy (10:4)
# 04 FFPE IBC vs Frozen nonIBC (10:4)
# 05 FFPE nonIBC vs Frozen Healthy (4:4)
# 06 FFPE nonIBC vs Frozen nonIBC (4:4)
# 07 PBMC IBC vs Healthy (10:16)
# 08 PBMC IBC vs nonIBC (10:6)
# 09 PBMC nonIBC vs Healthy (6:16)
# 10 Plasma IBC vs Healthy (10:16)
# 11 Plasma IBC vs nonIBC (10:6)
# 12 Plasma nonIBC vs Healthy (6:16)
# HR results HR+ vs HR-
# HR_01, FFPE nonIBC
# HR_02, FFPE IBC
# HR_03, PBMC nonIBC
# HR_04, PBMC IBC
# HR_05, Plasma nonIBC
# HR_06, Plasma IBC

snc<-subset(comp_res[[10]],Type2=="snoRNA" | Type2=="scaRNA")
snc$Type2<-droplevels(snc$Type2)
sig<-subset(snc, padj<=0.001 & abs(log2FoldChange)>=2)

snc<-subset(comp_res[[11]],Type2=="snoRNA" | Type2=="scaRNA")
snc$Type2<-droplevels(snc$Type2)
sig1<-subset(snc, padj<=0.001 & abs(log2FoldChange)>=2)

snc<-subset(comp_res[[12]],Type2=="snoRNA" | Type2=="scaRNA")
snc$Type2<-droplevels(snc$Type2)
sig2<-subset(snc, padj<=0.001 & abs(log2FoldChange)>=2)

setdiff(intersect(sig$Name[sig$log2FoldChange>0],sig1$Name[sig1$log2FoldChange>0]),sig2$Name[sig2$log2FoldChange>0])
setdiff(intersect(sig$Name[sig$log2FoldChange<0],sig1$Name[sig1$log2FoldChange<0]),sig2$Name[sig2$log2FoldChange<0])

