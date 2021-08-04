import pandas as pd
import os.path

bed_file = 'get_coord/gene_list.bed'
fusion_file = 'data/CosmicFusionExport_GRCh37.tsv'
# output_file = 'get_coord/fusion_gene_cosmic.tsv'
output_dir = 'cosmic_fusion_report_grch37'
bed_fhand = open(bed_file)
fusion_fhand = open(fusion_file)
# output_fhand = open(output_file, 'w')
PRIMARY_SITE = {'haematopoietic_and_lymphoid_tissue'}
FUSION_TYPE = {'Observed mRNA'}

def get_query_names(bed_fhand):
    """
    Get the query gene names from the bed file
    """
    query_names = set()
    for line in bed_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        query_names.add(line.split('\t')[4])
    return query_names

def get_fusion_table(fusion_fhand):
    """
    Get the fusion table from the tsv file
    """
    fusion_table = pd.read_csv(fusion_fhand, sep='\t')
    if PRIMARY_SITE:
        fusion_table = fusion_table[fusion_table['PRIMARY_SITE'].isin(PRIMARY_SITE)]
    if FUSION_TYPE:
        fusion_table = fusion_table[fusion_table['FUSION_TYPE'].isin(FUSION_TYPE)]
    return fusion_table

query_names = get_query_names(bed_fhand)
print(query_names)
fusion_table = get_fusion_table(fusion_fhand)
print(fusion_table.head())
end5_fusion_table = fusion_table[fusion_table["5'_GENE_NAME"].isin(query_names)]
end3_fusion_table = fusion_table[fusion_table["3'_GENE_NAME"].isin(query_names)]
one_end_fusion_table = fusion_table[fusion_table["3'_GENE_NAME"].isin(query_names) | fusion_table["5'_GENE_NAME"].isin(query_names)]
two_ends_fusion_table = fusion_table[fusion_table["3'_GENE_NAME"].isin(query_names) & fusion_table["5'_GENE_NAME"].isin(query_names)]
print(end5_fusion_table.head())
print(end3_fusion_table.head())
two_ends_fusion_table.to_csv(os.path.join(output_dir, 'two_ends.tsv'), sep = '\t', index = False)
one_end_fusion_table.to_csv(os.path.join(output_dir, 'one_end.tsv'), sep = '\t', index = False)
group_cols = ["FUSION_ID","TRANSLOCATION_NAME","5'_CHROMOSOME","5'_STRAND","5'_GENE_ID","5'_GENE_NAME","5'_LAST_OBSERVED_EXON","5'_GENOME_START_FROM","5'_GENOME_START_TO","5'_GENOME_STOP_FROM","5'_GENOME_STOP_TO","3'_CHROMOSOME","3'_STRAND","3'_GENE_ID","3'_GENE_NAME","3'_FIRST_OBSERVED_EXON","3'_GENOME_START_FROM","3'_GENOME_START_TO","3'_GENOME_STOP_FROM","3'_GENOME_STOP_TO"]
# group_cols2 = ["FUSION_ID","TRANSLOCATION_NAME","5'_GENE_NAME","3'_GENE_NAME"]
one_end_counts = one_end_fusion_table.groupby(by = group_cols).size().reset_index(name='counts')
one_end_counts.to_csv(os.path.join(output_dir, 'one_end_count.tsv'), sep = '\t', index = False)
two_end_counts = two_ends_fusion_table.groupby(by = group_cols).size().reset_index(name='counts')
two_end_counts.to_csv(os.path.join(output_dir, 'two_ends_count.tsv'), sep = '\t', index = False)
