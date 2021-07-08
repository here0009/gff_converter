#!/usr/bin/env Rscript
# Use R libarry org.Hs.eg.db to get the refseq_id from ensemb_id if refseq_id is missing, or vice vera.
# Supported types ENSEMBLTRANS, REFSEQ
# Input : A tsv file including name, refseq_id and ensembl_id of transcripts, missing values will be represented as 'None'
# cols of input file :
# Name  NCBI    ENSEMBL
# Usage: 
# Rscript add_transid.R <input_file> <output_file> <trans_type>
# Example1: 
# Rscript add_transid.R test/ensembl_gene_ids.tsv test/ensembl_gene_ids2.tsv e2n
# Example2: 
# Rscript add_transid.R test/ncbi_gene_ids.tsv test/ncbi_gene_ids2.tsv n2e
# Output: 
# tsv file with added refseq_ids and ensembl_ids

library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)

input_file = args[1]  # input_file = "test/ncbi_gene_ids.tsv"
output_file = args[2] # output_file = "test/ncbi_gene_ids2.tsv"
trans_type = args[3]  # e2n: ENSEMBL to NCBI; n2e: NCBI to ENSEMBL; both
gene_ids_table <- read.table(file = input_file, sep = "\t", header = TRUE)
counts = 0

if (trans_type == 'e2n' || trans_type == 'both'){
    ensembl_ids <- unlist(gene_ids_table[gene_ids_table[,2] == "None", ][,3])
    quries = length(ensembl_ids)
    trans_ids = tryCatch(mapIds(org.Hs.eg.db, keys=ensembl_ids, keytype="ENSEMBLTRANS", column="REFSEQ"), error = function(e) NULL)
    # mapIds may raise error, try to ignore it.
    if (!is.null(trans_ids)){
        for (i in 1:nrow(gene_ids_table)){
            name = gene_ids_table[i, 1]
            ncbi = gene_ids_table[i, 2]
            ensembl = gene_ids_table[i, 3]
            if (ncbi == "None" && ensembl != "None" && !is.na(trans_ids[ensembl])){
                # print(c(name, ncbi, ensembl, trans_ids[ensembl]))
                gene_ids_table[i, 2] = trans_ids[ensembl]
                counts = counts + 1
            }
        }
    }
}

if (trans_type == 'n2e' || trans_type == 'both'){
    gene_ids_table$NCBI2 <- sapply(strsplit(as.character(gene_ids_table$NCBI), "\\."), "[", 1) 
    # get the trimmed ncbi_id so we can use org.Hs.eg.db to do the query
    # e.g."NM_001256128.1" will be trimmed to "NM_001256128"
    ncbi_ids <- unlist(gene_ids_table[gene_ids_table[,3] == "None", ][,4])
    quries = length(ncbi_ids)
    trans_ids = tryCatch(mapIds(org.Hs.eg.db, keys=ncbi_ids, keytype="REFSEQ", column="ENSEMBLTRANS"), error = function(e) NULL)
    if (!is.null(trans_ids)){
        for (i in 1:nrow(gene_ids_table)){
            name = gene_ids_table[i, 1]
            ncbi = gene_ids_table[i, 2]
            ensembl = gene_ids_table[i, 3]
            ncbi2 = gene_ids_table[i, 4]
            if (ncbi != "None" && ensembl == "None" && !is.na(trans_ids[ncbi2])){
                # print(c(name, ncbi, ensembl, trans_ids[ncbi2]))
                gene_ids_table[i, 3] = trans_ids[ncbi2]
                counts = counts + 1
            }
        }
    }
}

sprintf("There are %s quries, add %s trans_ids", quries, counts)
write.table(gene_ids_table[, c("Name","NCBI","ENSEMBL")], file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
