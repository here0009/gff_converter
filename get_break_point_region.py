#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get the breakpoint region from the input cosmic file
"""

import pandas as pd
import re
from get_bed_file import get_gff, get_bed, get_bed_record
gff_style = 'ENSEMBL'
extra_bp = 10
CHOOSEN_TYPES = {'mRNA', 'transcript'}
CHOOSEN_SUB_TYPES = {'CDS', '3_UTR', '5_UTR', 'intron'}
OUTPUT_BED_SUB_TYPES = {'CDS'}
FUSION_BREAK_POINT_SUBTYPES = {'CDS', '3_UTR', '5_UTR', 'intron'}


input_fusion_file = 'cosmic_fusion_report_grch37/one_end_count.tsv'
extra_paterner_file = 'cosmic_fusion_report_grch37/extra_parterner.tsv'
ensemble_id_file = 'cosmic_fusion_report_grch37/ensemble_id.tsv'
gene_name_id = 'get_coord/trans_gene_id.tsv'  # gene_name, ncbi_id, ensembl_id in the input gene list
input_gff = 'gff_out/ENSEMBL_output.gff'
gff_fhand = open(input_gff)
break_point_output = 'cosmic_fusion_report_grch37/fusion_parterner_break_point.bed'
break_point2_output = 'cosmic_fusion_report_grch37/fusion_parterner_break_point2.bed'
fusion_parterner_bed = 'cosmic_fusion_report_grch37/fusion_parterner.bed'
fusion_parterner_fhand = open(fusion_parterner_bed, 'w')
break_point_fhand = open(break_point_output, 'w')
break_point2_fhand = open(break_point2_output, 'w')


def get_ensembl_ids(input_fusion_file, extra_paterner_file, ensemble_id_file, counts_threhold):
    """
    get the ensemble ids from fusion_counts_table of column table, where the counts are greater than the threshold.
    FUSION_ID	TRANSLOCATION_NAME	5'_CHROMOSOME	5'_STRAND	5'_GENE_ID	5'_GENE_NAME	5'_LAST_OBSERVED_EXON	5'_GENOME_START_FROM	5'_GENOME_START_TO	5'_GENOME_STOP_FROM	5'_GENOME_STOP_TO	3'_CHROMOSOME	3'_STRAND	3'_GENE_ID	3'_GENE_NAME	3'_FIRST_OBSERVED_EXON	3'_GENOME_START_FROM	3'_GENOME_START_TO	3'_GENOME_STOP_FROM	3'_GENOME_STOP_TO	counts
    Translocation_Name: 
    ENST00000254108.7(FUS):r.1_904_ENST00000442448.1(ERG):r.1141_5034
    """
    fusion_counts_table = pd.read_csv(input_fusion_file, sep='\t')
    fusion_counts_table = fusion_counts_table[fusion_counts_table['counts'] >= counts_threhold]
    fusion_parterner_regex = re.compile('(ENST\d+)\.?\d?\(([A-Za-z0-9]+)\)')
    fusion_counts_table['fusion_parterners'] = fusion_counts_table['TRANSLOCATION_NAME'].apply(lambda x: fusion_parterner_regex.findall(x))
    fusion_counts_table[fusion_counts_table['fusion_parterners'].str.len() != 2].to_csv(extra_paterner_file, sep='\t', index = False)
    fusion_counts_table = fusion_counts_table[fusion_counts_table['fusion_parterners'].str.len() == 2]
    fusion_counts_table["5'_GENE_NAME"] = fusion_counts_table["fusion_parterners"].apply(lambda x: x[0][1])
    fusion_counts_table["5'_ENSEMBL_ID"] = fusion_counts_table["fusion_parterners"].apply(lambda x: x[0][0])
    fusion_counts_table["3'_GENE_NAME"] = fusion_counts_table["fusion_parterners"].apply(lambda x: x[1][1])
    fusion_counts_table["3'_ENSEMBL_ID"] = fusion_counts_table["fusion_parterners"].apply(lambda x: x[1][0])
    fusion_counts_table.to_csv(ensemble_id_file, sep='\t', index = False)
    return fusion_counts_table


ensembl_id_table = get_ensembl_ids(input_fusion_file, extra_paterner_file, ensemble_id_file, 5)
query_id_dict = pd.Series(ensembl_id_table["5'_GENE_NAME"].values, index=ensembl_id_table["5'_ENSEMBL_ID"]).to_dict()
query_id_dict.update(pd.Series(ensembl_id_table["3'_GENE_NAME"].values, index=ensembl_id_table["3'_ENSEMBL_ID"]).to_dict())

name_id_table = pd.read_csv(gene_name_id, sep='\t')
genelist_ensembl_ids = set(name_id_table['ENSEMBL'].tolist())

id_gff_dict = get_gff(query_id_dict, gff_fhand, gff_style, CHOOSEN_TYPES, CHOOSEN_SUB_TYPES) # 1st parameter is query_id_dict or set
len_query_ensembl_ids, len_found_ids = len(query_id_dict), len(id_gff_dict)
not_found_fusion_ids = list(set(query_id_dict.keys()) - set(id_gff_dict.keys()))
print(f'There are {len_query_ensembl_ids} quries, found {len_found_ids}, the ids not found are:')
print(not_found_fusion_ids)
query_id2 = set(query_id_dict.keys()) - genelist_ensembl_ids  # the ids in fusion but not in the input gene list
query_id2_dict = {k: query_id_dict[k] for k in query_id2}
id_gff2_dict = {k: id_gff_dict[k] for k in query_id2} # the ids in fusion but not in the input gene list
print('The ensembl_id that not in list are:')
print([(query_id2_dict[id], id) for id in query_id2_dict])
get_bed(query_id2_dict, id_gff_dict, fusion_parterner_fhand, extra_bp, OUTPUT_BED_SUB_TYPES)


def get_break_point_bed(ensembl_id, start_from, stop_to, id_gff_dict, output_subtypes, output_gff_set, output_gff_list):
    """
    get the break point gff record list
    """
    
    if ensembl_id in id_gff_dict:
        gff_record = id_gff_dict[ensembl_id]
        children_list = []
        for _type in output_subtypes:
            children_list.extend(gff_record.get_children(_type))
        children_list.sort(key=lambda x: int(x.start))
        # print(children_list)
        for child_record in children_list:
            if int(child_record.start) <= start_from <= int(child_record.end) or int(child_record.start) <= stop_to <= int(child_record.end):
                if child_record not in output_gff_set:
                    output_gff_set.add(child_record)
                    output_gff_list.append(child_record)


def get_break_point_region(fusion_table, id_gff_dict, break_point_fhand, extra_bp, output_subtypes):
    """
    get the break point region of start_from and end_to
    """
    bed_header = '\t'.join(['#chrom', 'start', 'end', 'strand', 'name', 'rna_id', 'sub_type'])
    break_point_fhand.write(bed_header + '\n')
    output_gff_list = []
    output_gff_set = set()
    for _, row in fusion_table.iterrows():
        get_break_point_bed(row["5'_ENSEMBL_ID"], row["5'_GENOME_START_FROM"], row["5'_GENOME_STOP_TO"], id_gff_dict, output_subtypes, output_gff_set, output_gff_list)
        get_break_point_bed(row["3'_ENSEMBL_ID"], row["3'_GENOME_START_FROM"], row["3'_GENOME_STOP_TO"], id_gff_dict, output_subtypes, output_gff_set, output_gff_list)
    # print(output_gff_list)
    for record in output_gff_list:
        name = id_gff_dict[record.parent].name
        bed_record = get_bed_record(record, name, extra_bp)
        break_point_fhand.write(bed_record + '\t' + record.type + '\n')


get_break_point_region(ensembl_id_table, id_gff_dict, break_point_fhand, extra_bp, FUSION_BREAK_POINT_SUBTYPES)
get_break_point_region(ensembl_id_table, id_gff2_dict, break_point2_fhand, extra_bp, FUSION_BREAK_POINT_SUBTYPES)
