#!/usr/bin/env Rscript
# Use R libarry org.Hs.eg.db to translate NCBI refseq_id to ENSEMBL transcript_id if refseq_id is missing, or vice vera.
# Supported types ENSEMBLTRANS, REFSEQ
# Input : A tsv file including name, refseq_id and ensembl_id of transcripts, missing values will be represented as 'None'
# cols of input file :
# Name  NCBI    ENSEMBL
# Usage: 
# Rscript org_Hs_eg_db_translate.R input_file output_file trans_type
# Example:
# Test Data: 
# Rscript org_Hs_eg_db_translate.R test/ensembl_name_id.tsv test/R_ensembl_ids.tsv e2n
# Rscript org_Hs_eg_db_translate.R test/ncbi_name_id.tsv test/R_ncbi_name_id.tsv n2e
# Real Data:
# Rscript org_Hs_eg_db_translate.R gff_out/ncbi_name_id.tsv gff_out/R_trans_ncbi_name_id.tsv n2e
# Rscript org_Hs_eg_db_translate.R gff_out/ensembl_name_id.tsv gff_out/R_trans_ensembl_name_id.tsv e2n
# Output: 
# tsv file with tranlated refseq_ids and ensembl_ids


library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)

input_file = args[1]  
output_file = args[2] 
trans_type = args[3]  # e2n: ENSEMBL to NCBI; n2e: NCBI to ENSEMBL; both
gene_ids_table <- read.table(file = input_file, sep = "\t", header = TRUE)
total_lines = length(gene_ids_table[, 1])
trans_counts = 0

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
                trans_counts = trans_counts + 1
            }
        }
    }
}

if (trans_type == 'n2e' || trans_type == 'both'){
    # gene_ids_table$NCBI2 <- sapply(strsplit(as.character(gene_ids_table$NCBI), "\\."), "[", 1) 
    # get the trimmed ncbi_id so we can use org.Hs.eg.db to do the query
    # e.g."NM_001256128.1" will be trimmed to "NM_001256128" , already trimmed at the table
    # ncbi_ids <- unlist(gene_ids_table[gene_ids_table[,3] == "None", ][,4])
    ncbi_ids <- unlist(gene_ids_table[gene_ids_table[,3] == "None", ][,2])
    quries = length(ncbi_ids)
    trans_ids = tryCatch(mapIds(org.Hs.eg.db, keys=ncbi_ids, keytype="REFSEQ", column="ENSEMBLTRANS"), error = function(e) NULL)
    if (!is.null(trans_ids)){
        for (i in 1:nrow(gene_ids_table)){
            name = gene_ids_table[i, 1]
            ncbi = gene_ids_table[i, 2]
            ensembl = gene_ids_table[i, 3]
            # ncbi2 = gene_ids_table[i, 4]
            if (ncbi != "None" && ensembl == "None" && !is.na(trans_ids[ncbi])){
                # print(c(name, ncbi, ensembl, trans_ids[ncbi2]))
                gene_ids_table[i, 3] = trans_ids[ncbi]
                trans_counts = trans_counts + 1
            }
        }
    }
}

sprintf("There are %s total lines, %s quries, and %s got trans_ids.", total_lines, quries, trans_counts)
write.table(gene_ids_table[, c("Name","NCBI","ENSEMBL")], file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
