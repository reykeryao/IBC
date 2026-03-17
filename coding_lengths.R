## -----------------------------------------------------------------------------
tltn = gtf %>%
    filter(feature == "CDS") %>%
    mutate(trx_id = gsub('^.*=transcript:(.*?)(;|$).*', '\\1', attribute)) %>%
    mutate(strsign = c("-"=-1, "+"=+1)[strand]) %>%
    mutate(cdsstart = pmin(start * strsign, end * strsign)) %>%
    group_by(trx_id) %>%
    summarize(cdsstart = min(cdsstart)) %>%
    mutate(gene_id = trx[trx_id, "gene_id"]) %>%
    group_by(gene_id) %>%
    summarize(cdsstart = abs(min(cdsstart))) %>%
    as.data.frame
rownames(tltn) = tltn$gene_id
tltn$distance = abs(tltn$cdsstart - genes[rownames(tltn), "tss"])
tltn$gene_name = genes[rownames(tltn), "gene_name"]
## -----------------------------------------------------------------------------

rnalen = gtf %>%
    filter(feature == "exon") %>%
    mutate(trx_id = gsub('^.*=transcript:(.*?)(;|$).*', '\\1', attribute)) %>%
    mutate(gene_id = trx[trx_id, "gene_id"]) %>%
    mutate(len = end-start) %>%
    group_by(trx_id) %>%
    summarize(len=sum(len), gene_id=paste(unique(gene_id), collapse=",")) %>%
    group_by(gene_id) %>%
    summarize(len=max(len)) %>%
    as.data.frame
rownames(rnalen) = rnalen$gene_id
rnalen$gene_name = genes[rownames(rnalen), "gene_name"]

## -----------------------------------------------------------------------------
geneLens = geneNamize(data.frame(row.names = genes$gene_id,
                                 length = genes$end-genes$start), genes)
geneLens = structure(geneLens$length, names=rownames(geneLens))

rnaLens = structure(rnalen$len, names=rownames(geneNamize(rnalen, rnalen)))

fracIntron = (geneLens-rnaLens[names(geneLens)]) / geneLens
fracExon = 1 - fracIntron
