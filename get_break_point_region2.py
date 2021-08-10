#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage: get_break_point_region.py [-h] -f fusion_file -g gff_file -t trans_table -o output_dir

Find the break point region based on the input fusion File and gff file, get the break point region -+100bp

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
    python3 get_break_point_region2.py -f cosmic_fusion_report_grch37/red_gene/one_end_count.tsv -g gff_out/ENSEMBL_output_exon.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37/red_gene -i 10
    
    python3 get_break_point_region2.py -f cosmic_fusion_report_grch37/all_gene/one_end_count.tsv -g gff_out/ENSEMBL_output_exon.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37/all_gene -i 10

    python3 get_break_point_region2.py -f cosmic_fusion_report_grch37/all_gene/one_end_count.tsv -g gff_out/ENSEMBL_output_exon.gff -t get_coord/trans_gene_id.tsv -o cosmic_fusion_report_grch37/test -i 1
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
    fusion_threshold:int


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
    parser.add_argument('-i', '--fusion_threshold',
                        metavar='fusion_threshold',
                        help='Fusion Threshold',
                        type=int,
                        required=True
                        )
    args = parser.parse_args()
    return Arags(
                fusion_file=args.fusion_file,
                gff_file=args.gff_file,
                trans_table=args.trans_table,
                output_dir=args.output_dir,
                fusion_threshold=args.fusion_threshold
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


def add_record(record, record_list, ei_num, fusion_id, bp_start, bp_end, partner):
    """
    add the record to the record_set, ei_num is the number of exon/intron in the record.
    """
    record_list.append((record, ei_num, fusion_id, bp_start, bp_end, partner))


def get_break_point_bed(fusion_id, ensembl_id, chrom, bp_start, bp_end, partner, id_gff_dict, output_subtypes, output_gff_list):
    """
    get the break point gff record list
    """
    def check_break_point(gff_record):
        """
        check if the break point is on the border of the gff_record
        """
        if (partner == 'five_partner' and gff_record.strand == '+') or (partner == 'three_partner' and gff_record.strand == '-'):
            if (gff_record.type == 'exon' and int(gff_record.end) != bp_start) or (gff_record.type == 'intron' and int(gff_record.start) - 1 != bp_start):
                print(f"The break point of fusion_id: {fusion_id}, {ensembl_id} is not on the border, check it again")
                print(f'Break point is {bp_start}')
                print(gff_record)
        elif (partner == 'five_partner' and gff_record.strand == '-') or (partner == 'three_partner' and gff_record.strand == '+'):
            if (gff_record.type == 'exon' and int(gff_record.start) != bp_end) or (gff_record.type == 'intron' and int(gff_record.end) + 1 != bp_end):
                print(f"The break point of {fusion_id}, {ensembl_id} is not on the border, check it again")
                print(f'Break point is {bp_end}')
                print(gff_record)

    chrom_convert = {23:'X', 24:'Y'} # convert the chrom in cosmic record to the chrom in gff file
    if ensembl_id not in id_gff_dict:
        print(f"{fusion_id}, {ensembl_id} not in gff record!")
        return
    gff_record = id_gff_dict[ensembl_id]
    seqid = gff_record.seqid[3:] if gff_record.seqid.startswith('chr') else gff_record.seqid
    if str(seqid).lower() != str(chrom).lower() and seqid != chrom_convert.get(chrom, None):
        print(f"{fusion_id},{ensembl_id} and gff_record are not on the same chromosome!")
        print(f"gff_record: seqid: {seqid}, fusion record chrom: {chrom}")
        return
    children_list = []
    for _type in output_subtypes:
        children_list.extend(gff_record.get_children(_type))
    if len(children_list) == 0:
        print(f"{gff_record.id},{gff_record.name} have no children!")
        return
    children_list.sort(key=lambda x: int(x.start))
    # break point is before the 1st child, add the 500bp seq around break point to the output list
    if int(children_list[0].start) > bp_end:
        _s, _e = str(bp_start - 250), str(bp_end + 250)
        _atrb = f'Parent={gff_record.id};'
        _gff = GffRecord('\t'.join([gff_record.seqid, gff_record.source, "BreakPoint", _s, _e, gff_record.score, gff_record.strand, gff_record.phase, _atrb]), gff_record.source)
        add_record(_gff, output_gff_list, 0, fusion_id, bp_start, bp_end, partner)
     # break point is in the middle of some children, add the gffrecord of chilren to the output list
    exons, introns = 0, 0
    for gff_record in children_list:
        if gff_record.type == 'intron':
            introns += 1
        elif gff_record.type == 'exon':
            exons += 1
        if bp_start <= int(gff_record.start) <= bp_end or bp_start <= int(gff_record.end) <= bp_end:
            if gff_record.type == 'intron':
                _num = introns
            else:
                _num = exons
            check_break_point(gff_record)
            add_record(gff_record, output_gff_list, _num, fusion_id, bp_start, bp_end, partner)
    # break point is after the last child, add the 500bp seq around break point to the output list
    if int(children_list[-1].end) < bp_start:
        _s, _e = str(bp_start - 250), str(bp_end + 250)
        _atrb = f'Parent={gff_record.id};'
        _gff = GffRecord('\t'.join([gff_record.seqid, gff_record.source, "BreakPoint", _s, _e, gff_record.score, gff_record.strand, gff_record.phase, _atrb]), gff_record.source)
        add_record(_gff, output_gff_list, -1, fusion_id, bp_start, bp_end, partner)

def get_break_point_region(fusion_table, id_gff_dict, break_point_fhand, extra_bp, output_subtypes):
    """
    get the break point region of start_from and end_to
    """
    adj_length = 10
    bed_header = '\t'.join(['#chrom', 'start', 'end', 'strand', 'name', 'rna_id', 'sub_type',  'ei_num', 'fusion_id', 'break_point_start', 'break_point_end', 'partner'])
    break_point_fhand.write(bed_header + '\n')
    output_gff_list = [] #(gff_reocrd, ei_num, fusion_id, bp_start, bp_end, partner)
    # output_gff_set = set()
    for _, row in fusion_table.iterrows():
        fusion_id = row['FUSION_ID']
        if row["5'_STRAND"] == '+':
            bp_start, bp_end = row["5'_GENOME_STOP_FROM"], row["5'_GENOME_STOP_TO"] + adj_length
        else:
            bp_start, bp_end = row["5'_GENOME_START_FROM"] - adj_length, row["5'_GENOME_START_TO"]
        get_break_point_bed(fusion_id, row["5'_ENSEMBL_ID"], row["5'_CHROMOSOME"], bp_start, bp_end, 'five_partner', id_gff_dict, output_subtypes, output_gff_list)
        if row["3'_STRAND"] == '+':
            bp_start, bp_end = row["3'_GENOME_START_FROM"] - adj_length, row["3'_GENOME_START_TO"]
        else:
            bp_start, bp_end = row["3'_GENOME_STOP_FROM"], row["3'_GENOME_STOP_TO"] + adj_length
        get_break_point_bed(fusion_id, row["3'_ENSEMBL_ID"], row["3'_CHROMOSOME"], bp_start, bp_end, 'three_partner', id_gff_dict, output_subtypes, output_gff_list)
    # print(output_gff_list)
    for record, ei_num, fusion_id, bp_start, bp_end, part in output_gff_list:
        if not record.parent:
            print(record, ei_num, fusion_id, bp_start, bp_end, part)
            continue
        name = id_gff_dict[record.parent].name
        bed_record = get_bed_record(record, name, extra_bp)
        _record = '\t'.join([bed_record, record.type, str(ei_num), str(fusion_id), str(bp_start), str(bp_end), part])
        break_point_fhand.write(_record + '\n')


def main():
    """
    Main
    """
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)
    gff_style = 'ENSEMBL'
    extra_bp = 10
    CHOOSEN_TYPES = {'mRNA', 'transcript'}
    CHOOSEN_SUB_TYPES = {'exon', 'intron'}
    OUTPUT_BED_SUB_TYPES = {'CDS'}
    FUSION_BREAK_POINT_SUBTYPES = {'exon', 'intron'}
    ensembl_id_table = get_ensembl_ids(args.fusion_file, args.output_dir, args.fusion_threshold)
    query_id_dict = pd.Series(ensembl_id_table["5'_GENE_NAME"].values, index=ensembl_id_table["5'_ENSEMBL_ID"]).to_dict()
    query_id_dict.update(pd.Series(ensembl_id_table["3'_GENE_NAME"].values, index=ensembl_id_table["3'_ENSEMBL_ID"]).to_dict())
    print('The ensembl_id and gene_name of queries are:')
    print([(query_id_dict[id], id) for id in query_id_dict])
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
    break_point_all = os.path.join(args.output_dir, 'break_point_all.tsv') # all of the break points
    break_point_all_fhand = open(break_point_all, 'w')
    get_break_point_region(ensembl_id_table, id_gff_dict, break_point_all_fhand, extra_bp, FUSION_BREAK_POINT_SUBTYPES)


if __name__ == '__main__':
    main()
