#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage: get_break_point_region.py [-h] -f fusion_file -g gff_file -t trans_table -o output_dir

Find the break point region based on the input fusion File and gff file

optional arguments:
  -h, --help            show this help message and exit
  -f fusion_file, --fusion_file fusion_file
                        Input Fusion File, Cosmic Record (default: None)
  -g gff_file, --gff_file gff_file
                        Input GFF File (default: None)
  -t trans_table, --trans_table trans_table
                        TransTable of Gene Name, NCBI ID, ENSEMBL ID (default: None)
  -o output_dir, --output_dir output_dir
                        Output Directory (default: None)
''  # gene_name, ncbi_id, ensembl_id in the input gene list
    # input_gff = ''
Example:
    python3 get_break_point_region.py -f cosmic_fusion_report_grch37/red_gene/one_end_count.tsv -g gff_out/ENSEMBL_output.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37/red_gene
    python3 get_break_point_region.py -f cosmic_fusion_report_grch37/all_gene/one_end_count.tsv -g gff_out/ENSEMBL_output.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37

    python3 get_break_point_region.py -f cosmic_fusion_report_grch37/all_gene/one_end_count.tsv -g gff_out/ENSEMBL_output.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37/test
"""


from GFF import GffRecord
import pandas as pd
from typing import NamedTuple, TextIO
import re
from get_bed_file import get_gff, get_bed, get_bed_record
import argparse
import os.path
import os

class Arags(NamedTuple):
    """
    Command Line arguments
    """
    fusion_file:TextIO
    gff_file:TextIO
    trans_table:TextIO
    output_dir:str


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Find the break point region based on the input fusion File and gff file", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-f', '--fusion_file',
                        metavar='fusion_file',
                        type=argparse.FileType('r'),
                        help='Input Fusion File, Cosmic Record',
                        required=True
                        )

    parser.add_argument('-g', '--gff_file',
                        metavar='gff_file',
                        type=argparse.FileType('r'),
                        help='Input GFF File',
                        required=True
                        )
    parser.add_argument('-t', '--trans_table',
                        metavar='trans_table',
                        type=argparse.FileType('r'),
                        help='TransTable of Gene Name, NCBI ID, ENSEMBL ID',
                        required=True
                        )
    parser.add_argument('-o', '--output_dir',
                        metavar='output_dir',
                        help='Output Directory',
                        required=True
                        )
    args = parser.parse_args()
    return Arags(
                fusion_file=args.fusion_file,
                gff_file=args.gff_file,
                trans_table=args.trans_table,
                output_dir=args.output_dir
                )


def get_ensembl_ids(input_fusion_file, output_dir, counts_threhold):
    """
    get the ensemble ids from fusion_counts_table of column table, where the counts are greater than the threshold.
    FUSION_ID	TRANSLOCATION_NAME	5'_CHROMOSOME	5'_STRAND	5'_GENE_ID	5'_GENE_NAME	5'_LAST_OBSERVED_EXON	5'_GENOME_START_FROM	5'_GENOME_START_TO	5'_GENOME_STOP_FROM	5'_GENOME_STOP_TO	3'_CHROMOSOME	3'_STRAND	3'_GENE_ID	3'_GENE_NAME	3'_FIRST_OBSERVED_EXON	3'_GENOME_START_FROM	3'_GENOME_START_TO	3'_GENOME_STOP_FROM	3'_GENOME_STOP_TO	counts
    Translocation_Name: 
    ENST00000254108.7(FUS):r.1_904_ENST00000442448.1(ERG):r.1141_5034
    """
    extra_paterner_file = os.path.join(output_dir, 'extra_parterner.tsv') # store the fusion got parterners not equal to 2.
    # ensemble_id_file = os.path.join(output_dir, 'ensemble_ids.txt')
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
    return fusion_counts_table


def get_break_point_bed(fusion_id, ensembl_id, chrom, break_points, id_gff_dict, output_subtypes, output_gff_set, output_gff_list):
    """
    get the break point gff record list
    """
    if len(break_points) == 0:
        print(f"{fusion_id},{ensembl_id} have no break points!")
        return
    if ensembl_id not in id_gff_dict:
        print(f"{fusion_id}, {ensembl_id} not in gff record!")
        return
    gff_record = id_gff_dict[ensembl_id]
    # if gff_record.seqid != chrom:
    #     print(f"{fusion_id},{ensembl_id} and gff_record are not on the same chromosome!")
    #     return
    children_list = []
    for _type in output_subtypes:
        children_list.extend(gff_record.get_children(_type))
    if len(children_list) == 0:
        print(f"{gff_record.id},{gff_record.name} have no children!")
        return
    children_list.sort(key=lambda x: int(x.start))
    # print("The Fusions are:")
    # print(fusion_id, break_points)
    # print(children_list)
    idx_child = 0
    idx_breaks = 0
    # break point is before the 1st child, add the 500bp seq around break point to the output list
    while idx_breaks < len(break_points) and break_points[idx_breaks] < int(children_list[0].start): 
        _s, _e = str(break_points[idx_breaks] - 250), str(break_points[idx_breaks] + 250)
        idx_breaks += 1
        _gff = GffRecord('\t'.join([gff_record.seqid, gff_record.source, "BreakPoint", _s, _e, gff_record.score, gff_record.strand, gff_record.phase, gff_record.atrb.string]), gff_record.source)
        if _gff not in output_gff_set:
            output_gff_list.append(_gff)
            output_gff_set.add(_gff)
        idx_breaks += 1
    # break point is in the middle of some children, add the gffrecord of chilren to the output list
    while idx_breaks < len(break_points):
        while idx_child < len(children_list) and break_points[idx_breaks] > int(children_list[idx_child].end):
            idx_child += 1
        if idx_child >= len(children_list):
            break
        if int(children_list[idx_child].start) <= break_points[idx_breaks] <= int(children_list[idx_child].end):
            if children_list[idx_child] not in output_gff_set:
                output_gff_list.append(children_list[idx_child])
                output_gff_set.add(children_list[idx_child])
        idx_breaks += 1
    # break point is after the last child, add the 500bp seq around break point to the output list
    while idx_breaks < len(break_points) and break_points[idx_breaks] > int(children_list[-1].end):
        _s, _e = str(break_points[idx_breaks] - 250), str(break_points[idx_breaks] + 250)
        idx_breaks += 1
        _gff = GffRecord('\t'.join([gff_record.seqid, gff_record.source, "BreakPoint", _s, _e, gff_record.score, gff_record.strand, gff_record.phase, gff_record.atrb.string]), gff_record.source)
        if _gff not in output_gff_set:
            output_gff_list.append(_gff)
            output_gff_set.add(_gff)
        idx_breaks += 1



def get_break_point_region(fusion_table, id_gff_dict, break_point_fhand, extra_bp, output_subtypes):
    """
    get the break point region of start_from and end_to
    """
    bed_header = '\t'.join(['#chrom', 'start', 'end', 'strand', 'name', 'rna_id', 'sub_type'])
    break_point_fhand.write(bed_header + '\n')
    output_gff_list = []
    output_gff_set = set()
    for _, row in fusion_table.iterrows():
        fusion_id = row['FUSION_ID']
        break_points = sorted(list(set([row["5'_GENOME_START_FROM"], row["5'_GENOME_START_TO"], row["5'_GENOME_STOP_FROM"], row["5'_GENOME_STOP_TO"]])))
        get_break_point_bed(fusion_id, row["5'_ENSEMBL_ID"], row["5'_CHROMOSOME"], break_points, id_gff_dict, output_subtypes, output_gff_set, output_gff_list)
        break_points = sorted(list(set([row["3'_GENOME_START_FROM"], row["3'_GENOME_START_TO"], row["3'_GENOME_STOP_FROM"], row["3'_GENOME_STOP_TO"]])))
        get_break_point_bed(fusion_id, row["3'_ENSEMBL_ID"], row["3'_CHROMOSOME"], break_points, id_gff_dict, output_subtypes, output_gff_set, output_gff_list)
    # print(output_gff_list)
    for record in output_gff_list:
        name = id_gff_dict[record.parent].name
        bed_record = get_bed_record(record, name, extra_bp)
        break_point_fhand.write(bed_record + '\t' + record.type + '\n')


def main():
    """
    Main
    """
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)
    gff_style = 'ENSEMBL'
    extra_bp = 10
    fusion_record_threshold = 5
    CHOOSEN_TYPES = {'mRNA', 'transcript'}
    CHOOSEN_SUB_TYPES = {'CDS', '3_UTR', '5_UTR', 'intron'}
    OUTPUT_BED_SUB_TYPES = {'CDS'}
    FUSION_BREAK_POINT_SUBTYPES = {'CDS', '3_UTR', '5_UTR', 'intron'}
    ensembl_id_table = get_ensembl_ids(args.fusion_file, args.output_dir, fusion_record_threshold)
    # get the query ensemble ids in ensembl_id_table
    query_id_dict = pd.Series(ensembl_id_table["5'_GENE_NAME"].values, index=ensembl_id_table["5'_ENSEMBL_ID"]).to_dict()
    query_id_dict.update(pd.Series(ensembl_id_table["3'_GENE_NAME"].values, index=ensembl_id_table["3'_ENSEMBL_ID"]).to_dict())
    # get the ensembl ids in input gene list
    name_id_table = pd.read_csv(args.trans_table, sep='\t')
    genelist_ensembl_ids = set(name_id_table['ENSEMBL'].tolist())
    # get the gff_record of query ids
    id_gff_dict = get_gff(query_id_dict, args.gff_file, gff_style, CHOOSEN_TYPES, CHOOSEN_SUB_TYPES)
    len_query_ensembl_ids, len_found_ids = len(query_id_dict), len(id_gff_dict)
    not_found_fusion_ids = list(set(query_id_dict.keys()) - set(id_gff_dict.keys()))
    print(f'There are {len_query_ensembl_ids} quries, found {len_found_ids}, the ids not found are:')
    print(not_found_fusion_ids)
    query_id2 = set(query_id_dict.keys()) - genelist_ensembl_ids  # the ids in fusion file but not in the input gene list
    query_id2_dict = {k: query_id_dict[k] for k in query_id2}
    print('The ensembl_id that not in the input gene list are:')
    print([(query_id2_dict[id], id) for id in query_id2_dict])
    # export the fusion parterner information to a bed file
    fusion_parterner_bed = os.path.join(args.output_dir, 'fusion_parterner.bed')
    fusion_parterner_fhand = open(fusion_parterner_bed, 'w')
    get_bed(query_id2_dict, id_gff_dict, fusion_parterner_fhand, extra_bp, OUTPUT_BED_SUB_TYPES)
    # export the region contain break point to a bed file
    break_point_all = os.path.join(args.output_dir, 'break_point_all.bed') # all of the break points
    break_point_all_fhand = open(break_point_all, 'w')
    get_break_point_region(ensembl_id_table, id_gff_dict, break_point_all_fhand, extra_bp, FUSION_BREAK_POINT_SUBTYPES)
    # id_gff2_dict = {k: id_gff_dict[k] for k in query_id2} # the ids in fusion but not in the input gene list
    # break_point_partial = os.path.join(args.output_dir, 'break_point_partial.bed') # break points the not contain the already selected transcripts.
    # break_point_partial_fhand = open(break_point_partial, 'w')
    # get_break_point_region(ensembl_id_table, id_gff2_dict, break_point_partial_fhand, extra_bp, FUSION_BREAK_POINT_SUBTYPES)


if __name__ == '__main__':
    main()
