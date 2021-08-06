#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage: find_fusion_cosmic_table.py [-h] -i input_bed -f fusion_table -o output_dir

Find the Gene Fusion Record in Cosmic Table of the given BED file, Group the Same Record Togther

optional arguments:
  -h, --help            show this help message and exit
  -i input_bed, --input input_bed
                        Input BED File (default: None)
  -f fusion_table, --fusion_table fusion_table
                        Fusion Table (default: None)
  -o output_dir, --output_dir output_dir
                        Output Directory (default: None)

Example:
    python3 find_fusion_cosmic_table.py -i get_coord/red_gene/red_gene_list.bed -f data/CosmicFusionExport_GRCh37.tsv -o cosmic_fusion_report_grch37/red_gene

    python3 find_fusion_cosmic_table.py -i get_coord/gene_list.bed -f data/CosmicFusionExport_GRCh37.tsv -o cosmic_fusion_report_grch37/all_gene

"""


import pandas as pd
import os.path
import argparse
from typing import NamedTuple, TextIO
import os

class Arags(NamedTuple):
    """
    Command Line arguments
    """
    bed_file:TextIO
    fusion_table:TextIO
    output_dir:str



def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Find the Gene Fusion Record in Cosmic Table of the given BED file, Group the Same Record Togther", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input',
                        metavar='input_bed',
                        type=argparse.FileType('r'),
                        help='Input BED File',
                        required=True
                        )

    parser.add_argument('-f', '--fusion_table',
                        metavar='fusion_table',
                        type=argparse.FileType('r'),
                        help='Fusion Table',
                        required=True
                        )
    parser.add_argument('-o', '--output_dir',
                        metavar='output_dir',
                        help='Output Directory',
                        required=True
                        )
    args = parser.parse_args()
    return Arags(bed_file=args.input, 
                fusion_table=args.fusion_table,
                output_dir=args.output_dir
                )


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

def get_fusion_table(fusion_fhand, PRIMARY_SITE, FUSION_TYPE):
    """
    Get the fusion table from the tsv file
    """
    fusion_table = pd.read_csv(fusion_fhand, sep='\t')
    if PRIMARY_SITE:
        fusion_table = fusion_table[fusion_table['PRIMARY_SITE'].isin(PRIMARY_SITE)]
    if FUSION_TYPE:
        fusion_table = fusion_table[fusion_table['FUSION_TYPE'].isin(FUSION_TYPE)]
    # save the data type as int or string, otherwise it will be converted to float
    fusion_table = fusion_table.astype({"5'_LAST_OBSERVED_EXON":int,"5'_GENOME_START_FROM":int,"5'_GENOME_START_TO":int,"5'_GENOME_STOP_FROM":int,"5'_GENOME_STOP_TO":int,"3'_FIRST_OBSERVED_EXON":int,"3'_GENOME_START_FROM":int,"3'_GENOME_START_TO":int,"3'_GENOME_STOP_FROM":int,"3'_GENOME_STOP_TO":int}) 
    return fusion_table

def get_fusion_report(fusion_table, query_names, output_dir):
    """
    Get the fusion record and the counts of fusion record in tsv file
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"The output dir is {output_dir}")
    one_end_fusion_table = fusion_table[fusion_table["3'_GENE_NAME"].isin(query_names) | fusion_table["5'_GENE_NAME"].isin(query_names)]
    two_ends_fusion_table = fusion_table[fusion_table["3'_GENE_NAME"].isin(query_names) & fusion_table["5'_GENE_NAME"].isin(query_names)]
    # aviod the trailing zeros by '{0:g}'.format()
    two_ends_fusion_table.to_csv(os.path.join(output_dir, 'two_ends.tsv'), sep = '\t', index = False, float_format = '{0:g}'.format)
    one_end_fusion_table.to_csv(os.path.join(output_dir, 'one_end.tsv'), sep = '\t', index = False, float_format = '{0:g}'.format)
    group_cols = ["FUSION_ID","TRANSLOCATION_NAME","5'_CHROMOSOME","5'_STRAND","5'_GENE_ID","5'_GENE_NAME","5'_LAST_OBSERVED_EXON","5'_GENOME_START_FROM","5'_GENOME_START_TO","5'_GENOME_STOP_FROM","5'_GENOME_STOP_TO","3'_CHROMOSOME","3'_STRAND","3'_GENE_ID","3'_GENE_NAME","3'_FIRST_OBSERVED_EXON","3'_GENOME_START_FROM","3'_GENOME_START_TO","3'_GENOME_STOP_FROM","3'_GENOME_STOP_TO"]
    one_end_counts = one_end_fusion_table.groupby(by = group_cols).size().reset_index(name='counts')
    one_end_counts.to_csv(os.path.join(output_dir, 'one_end_count.tsv'), sep = '\t', index = False, float_format = '{0:g}'.format)
    print(f"There are {one_end_counts.shape[0]} one end fusions")
    two_ends_counts = two_ends_fusion_table.groupby(by = group_cols).size().reset_index(name='counts')
    two_ends_counts.to_csv(os.path.join(output_dir, 'two_ends_count.tsv'), sep = '\t', index = False, float_format = '{0:g}'.format)
    print(f"There are {two_ends_counts.shape[0]} two ends fusions")

def main():
    PRIMARY_SITE = {'haematopoietic_and_lymphoid_tissue'}
    FUSION_TYPE = {'Observed mRNA'}
    args = get_args()
    query_names = get_query_names(args.bed_file)
    print("The names of the query genes are:")
    print(query_names)
    fusion_table = get_fusion_table(args.fusion_table, PRIMARY_SITE, FUSION_TYPE)
    print("The head of the fusion table is:")
    print(fusion_table.head())
    get_fusion_report(fusion_table, query_names, args.output_dir)


if __name__ == '__main__':
    main()
