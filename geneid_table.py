#!/usr/bin/env python3

"""
generate a table contain NCBI, ENSEMBL and NAME of transcripts and mRNA, cols of the file:
#NAME  NCBI    ENSEML
usage:
python3 geneid_table.py data/gene2ensembl ensembl_test_output.gff gene_id.txt ENSEMBL
"""
from GFF import GffRecord
import sys


geneid_table = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]
gff_type = sys.argv[4]
assert gff_type in {'NCBI', 'ENSEMBL'}
if gff_type == 'NCBI':
    source_col = 3
    dest_col = 4
elif gff_type == 'ENSEMBL':
    source_col = 4
    dest_col = 3
# cols of the gene2ENSEMBL file
#tax_id GeneID  Ensembl_gene_identifier RNA_nucleotide_accession.version        Ensembl_rna_identifier  protein_accession.version       Ensembl_protein_identifier

trans_dict = {'None':'None'}
with open(geneid_table) as geneid_table_fhand:
    for line in geneid_table_fhand:
        # print(line)
        if line.startswith('#'):
            continue
        lst = line.split('\t')
        # print(lst)
        ncbi_id = lst[3]
        ensembl_id = lst[4].split('.')[0]
        if gff_type == 'NCBI':
            trans_dict[ncbi_id] = ensembl_id
        elif gff_type == 'ENSEMBL':
            trans_dict[ensembl_id] = ncbi_id

# total records and records got trans id
line_counts = 0
trans_line_counts = 0
output_file_fhand = open(output_file, 'w')
output_file_fhand.write('\t'.join(['#Name', 'NCBI', 'ENSEMBL']) + '\n')
with open(input_file) as input_file_fhand:
    for line in input_file_fhand:
        lst = line.split('\t')
        if len(lst) < 9:
            continue
        gff_record = GffRecord(line, gff_type)
        if gff_record.type in {'mRNA', 'transcript'}:
            line_counts += 1
            name = gff_record.name if gff_record.name else 'None'
            source_id = gff_record.id if gff_record.id else 'None'
            trans_id = trans_dict.get(gff_record.id, 'None')
            trans_line_counts += trans_id != 'None'
            if gff_type == 'NCBI':
                w_line = '\t'.join([name, source_id, trans_id])
            else:
                w_line = '\t'.join([name, trans_id, source_id])
            output_file_fhand.write(w_line + '\n')

output_file_fhand.write('\n')
output_file_fhand.close()
print(f'There are {line_counts} records {trans_line_counts} got trans id')
            



