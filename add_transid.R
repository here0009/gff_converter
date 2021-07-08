#!/usr/bin/env Rscript

# Generate a tab file of ensembl id and refseq id
# Usage: Rscript add_transid.R gene_id.tsv gene_id2.tsv
# cols of input file :
# NAME  NCBI    ENSEML
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)
# supported types ENSEMBLTRANS, REFSEQ
input_file = args[1]
output_file = args[2]
gene_id <- read.table(file = input_file, sep = "\t", col.names = c("Name", "NCBI", "ENSEMBL"))

ensembl_ids <- unlist(gene_id[3])
trans_id = tryCatch(
    mapIds(org.Hs.eg.db, keys=ensembl_ids, keytype="ENSEMBLTRANS", column="REFSEQ"), 
    error = function(e) NULL)
# mapIds may raise error, try to ignore it.
counts = 0
print(gene_id[1,1])
for (i in 1:nrow(gene_id)){
    name = gene_id[i, 1]
    ncbi = gene_id[i, 2]
    ensembl = gene_id[i, 3]
    # print(c(name, ncbi, ensembl))
    if (ncbi == "None" && ensembl != "None" && !is.na(trans_id[ensembl])){
        print(c(name, ncbi, ensembl, trans_id[ensembl]))
        gene_id[i, 2] = trans_id[ensembl]
        counts = counts + 1
    }
    else if (ncbi != "None" && ensembl == "None" && !is.na(trans_id[ncbi])){
        print(c(name, ncbi, ensembl, trans_id[ncbi]))
        gene_id[i, 3] = trans_id[ncbi]
        counts = counts + 1
    }
}

sprintf("Add %s lines", counts)
write.table(gene_id, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)

