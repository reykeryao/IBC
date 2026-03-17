geneNamize = function(df, gtf, name='gene_name') {
    geneNames = gtf[rownames(df), name]
    names(geneNames) = rownames(df)
    geneNames[is.na(geneNames)] = names(geneNames)[is.na(geneNames)]
    geneNameCounts = table(geneNames)
    uniqueGeneNames = geneNames[geneNames %in%
                                names(which(geneNameCounts == 1))]
    geneNameDegenerate = !(geneNames %in% uniqueGeneNames)
    geneNames[geneNameDegenerate] = names(geneNames)[geneNameDegenerate]
    rownames(df) = geneNames
    return(df)
}
