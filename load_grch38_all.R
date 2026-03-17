gtfFile = 'Homo_sapiens.GRCh38.109.gff3.gz'
if (!file.exists(gtfFile)) {
    download.file(paste0(
        "https://download.rgd.mcw.edu/pub/data_release/GFF3/Ensembl/",
        gtfFile
    ), destfile=gtfFile)
}
gtf = read.table(gtfFile,
                 sep='\t', row.names=NULL, header=FALSE,
                 comment.char='#', quote='',
                 check.names=FALSE, stringsAsFactors=FALSE)
colnames(gtf) = c(
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
)
gtf$gene_id = gsub('^.*=gene:(.*?)(;|$).*', '\\1', gtf$attribute)
gtf$gene_name = gsub('^.*Name=(.*?)(;|$).*', '\\1', gtf$attribute)
gtf$biotype = gsub('^.*biotype=(.*?)(;|$).*', '\\1', gtf$attribute)

## gtf = gtf[grepl('\\d+$|^[X|Y]$', gtf$seqname), ]

genes = gtf[grepl('^ID=gene:', gtf$attribute), ]
genes$tss = genes$start
genes[genes$strand == '-', 'tss'] = genes[genes$strand == '-', 'end']
rownames(genes) = genes$gene_id

trx = gtf[grepl('^ID=transcript:', gtf$attribute), ]
rownames(trx) = gsub('^.*=transcript:(.*?)(;|$).*', '\\1', trx$attribute)
trx$tss = trx$start
trx[trx$strand == '-', 'tss'] = trx[trx$strand == '-', 'end']
