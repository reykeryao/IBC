#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(ggrepel)
library(MASS)
library(matrixStats)
library(scales)
library(tidyr)
library(WriteXLS)

source("load_grch38_all.R")
source("genenamize.R")
source("coding_lengths.R")

genes$gene_name = gsub("CCN6", "WISP3", genes$gene_name)
genes$gene_name = gsub("ARMH4", "C14orf37", genes$gene_name)
names(fracExon) = gsub("CCN6", "WISP3", names(fracExon))
names(fracExon) = gsub("ARMH4", "C14orf37", names(fracExon))

allCounts = read.table("cleaned.combined_with_Rep_miRNA_processed.counts.gz",
                       sep="\t", header=TRUE, row.names=1, check.names=FALSE)
fixSampleNames = function(s) {
    s = gsub("BC", "FT", s)
    s = gsub("^(FT|FFPE|PBMC|Plasma)", "\\1-", s)
    s = gsub("^FT-(\\d)", "FT-B\\1", s)
    s = gsub("^(FFPE|PBMC|Plasma)-(\\d+)$", "\\1-I\\2", s)
    s = gsub("-([BHI])([1-9])$", "-\\10\\2", s)
    s = gsub("-", "_", s)
    return(s)
}
colnames(allCounts) = fixSampleNames(colnames(allCounts))

exCounts = read.table("transcriptome.ENSG.counts.gz",
                      sep="\t", header=TRUE, row.names=1, check.names=FALSE)
colnames(exCounts) = gsub("\\.", "_", colnames(exCounts))
exCounts = exCounts[ , colnames(exCounts) %in% colnames(allCounts)]
exCountsT = read.table("transcriptome.ENST.counts.gz",
                       sep="\t", header=TRUE, row.names=1, check.names=FALSE)

interGenes = intersect(rownames(exCounts), rownames(allCounts))

exRatios = exCounts[interGenes, ] / allCounts[interGenes, colnames(exCounts)]
exRatios[exRatios > 1] = 1
exRatios = geneNamize(exRatios, allCounts, "ID")

inRatios = 1 - exRatios

## =============================================================================
## g -> (ag) / ((a-1)g + 1) = t
## -----------------------------------------------------------------------------
intronicDepthRatio = function(g, t, max.=10, min.=0.01) {
    if (is.vector(g) && length(g) > 1 && length(g) == length(t)) {
        return(mapply(intronicDepthRatio, g, t))
    }
    if (is.na(g) || is.na(t)) {return(NA)}
    if (g == 0) {return(min.)}
    if (t == 1) {return(max.)}
    out = try(
        uniroot(function(a) {t - a*g/((a-1)*g+1)}, interval=c(0, 1000))$root,
        silent = TRUE
    )
    if (!is(out, "try-error")) {
        return(max(min., min(max., out)))
    } else {
        return(NA)
    }
}

## -----------------------------------------------------------------------------
ffpeIdrs = intronicDepthRatio(
    1 - fracExon[rownames(exRatios)],
    rowMeans(inRatios[ , grep("FFPE", colnames(inRatios))], na.rm=TRUE)
)
pbmcIdrs = intronicDepthRatio(
    1 - fracExon[rownames(exRatios)],
    rowMeans(inRatios[ , grep("PBMC", colnames(inRatios))], na.rm=TRUE)    
)

## -----------------------------------------------------------------------------
patternToName = c(FT_H = "Frozen Healthy",
                  FT_B = "Frozen non-IBC",
                  FFPE_B = "FFPE non-IBC",
                  FFPE_I = "FFPE IBC",
                  PBMC_H = "PBMC Healthy",
                  PBMC_B = "PBMC non-IBC",
                  PBMC_I = "PBMC IBC",
                  Plasma_H = "Plasma Healthy",
                  Plasma_B = "Plasma non-IBC",
                  Plasma_I = "Plasma IBC")

groupwiseIntronicFracs = list()
for (samplePattern in names(patternToName)) {
    gif = data.frame("genomic intron fraction" = 1 - fracExon[rownames(exRatios)],
                     "transcriptome intronic fraction" = rowMeans(inRatios[ ,
                         grepl(paste0("^", samplePattern), colnames(inRatios))
                     ], na.rm=TRUE),
                     length = log10(geneLens[rownames(inRatios)]),
                     check.names = FALSE)
    groupwiseIntronicFracs[[samplePattern]] = gif
}
idrs = sapply(groupwiseIntronicFracs, function(.) {
    intronicDepthRatio(.$`genomic intron fraction`,
                       .$`transcriptome intronic fraction`)
})
rownames(idrs) = rownames(groupwiseIntronicFracs[[1]])

# -----------------------------------------------------------------------------
write.table(data.frame(gene=rownames(idrs), idrs),
            "intronic_depth_ratio.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
R.utils::gzip("intronic_depth_ratio.tsv", overwrite=TRUE)
