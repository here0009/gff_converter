#!/usr/bin/env python3
"""
usage: get_coord.py [-h] -i input_file -o output_bed -g gff_file [-s gff_style] [-e extra_bp]

Get coordinates of the input gene names based on GFF file, and get the output bed file

optional arguments:
  -h, --help            show this help message and exit
  -i input_file, --input_file input_file
                        Input Gene Name File (default: None)
  -o output_bed, --output_file output_bed
                        Output BED File (default: None)
  -g gff_file, --gff_file gff_file
                        Input Gff File (default: None)
  -s gff_style, --style gff_style
                        Gff style, "ENSEMBL" or "NCBI" or "UCSC" (default: "UCSC")
  -e extra_bp, --extra_bp extra_bp

Example:
Test:
    python3 get_coord.py -i get_coord/gene_list.txt -o get_coord/output_hg19_trimmed.bed -g get_coord/hg19.gff  -e 10
    python3 get_coord.py -i get_coord/gene_list.txt -o get_coord/output_hg19_UCSC.bed -g hg19_gff/hg19.ncbiRefSeq.gtf  -e 10
    python3 get_coord.py -i get_coord/gene_list.txt -o get_coord/output_hg19_refGene.bed -g hg19_gff/hg19.refGene.gtf  -e 10

Note:
    The gff file used for get coordinates is downloaded from UCSC.
    http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/
"""


import argparse
from os import truncate
from typing import Counter, NamedTuple, OrderedDict, TextIO
from GFF import GffRecord, UCSC_TRANSCRIPT


class Arags(NamedTuple):
    """
    Command Line arguments
    """
    input_file:TextIO
    gff_file:TextIO
    output_file:TextIO
    gff_style:str
    extra_bp:int


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Get coordinates of the input gene names based on GFF file, and get the output bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input_file',
                        metavar='input_file',
                        type=argparse.FileType('r'),
                        help='Input Gene Name File',
                        required=True
                        )

    parser.add_argument('-o', '--output_file',
                        metavar='output_bed',
                        type=argparse.FileType('w'),
                        help='Output BED File',
                        required=True
                        )
    parser.add_argument('-g', '--gff_file', 
                        metavar='gff_file', 
                        type=argparse.FileType('r'), 
                        help='Input Gff File', 
                        required=True
                        )
    parser.add_argument('-s', '--style',
                        metavar='gff_style',
                        choices=['ENSEMBL', 'NCBI', 'UCSC'],
                        help='Gff style, "ENSEMBL" or "NCBI"' or 'UCSC',
                        default='UCSC'
                        )
    parser.add_argument('-e', '--extra_bp',
                        metavar='extra_bp',
                        type=int)
    args = parser.parse_args()
    return Arags(input_file=args.input_file, 
                output_file=args.output_file, 
                gff_file=args.gff_file, 
                gff_style=args.style, 
                extra_bp=args.extra_bp
                )


def get_query_name(input_fhand):
    """
    Get the gene in an orderdict
    """
    query_counts = 0
    query_name_dict = OrderedDict()
    for line in input_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        query_counts += 1
        gene_name = line.upper()  # input gene name may be lower case
        query_name_dict[gene_name] = []
    print(query_name_dict)
    return query_name_dict

def get_gff(query_name_dict, gff_fhand, gff_style):
    """
    Get the gff file
    """
    for line in gff_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        gff_record = UCSC_TRANSCRIPT(line, gff_style)
        if gff_record.type == 'transcript' and gff_record.name in query_name_dict:
            query_name_dict[gff_record.name].append(gff_record)
    query_counts = Counter()
    for query_name in query_name_dict:
        query_counts[query_name] = len(query_name_dict[query_name])
    print('The transcript ids for each query are:')
    query_names = list(query_counts.keys())
    query_names.sort(key = lambda x: query_counts[x] , reverse=True)
    for key in query_names:
        print('\t'.join([key, str(query_counts[key])]))

def get_bed(query_name_dict, output_fhand, extra_bp):
    """
    Get the bed file
    """
    bed_header = '\t'.join(['#chrom', 'start', 'end', 'strand', 'name', 'rna_id'])
    output_fhand.write(bed_header + '\n')
    not_found_names = []
    for key, gff_lst in query_name_dict.items():
        if len(gff_lst) == 0:
            not_found_names.append(key)
        for gff_record in gff_lst:
            bed_record = '\t'.join([gff_record.seqid, 
                                    str(int(gff_record.start) - extra_bp), 
                                    str(int(gff_record.end) + extra_bp + 1), 
                                    gff_record.strand, 
                                    gff_record.name,
                                    gff_record.id])
            output_fhand.write(bed_record + '\n')
    not_found_counts = len(not_found_names)
    total_queries = len(query_name_dict)
    print(f"There are {total_queries} total_queries")
    print(f'{not_found_counts} not found, the gene names not found are:')
    print(not_found_names)
        

def main():
    args = get_args()
    # name_gff_dict = get_dict(args.gff_file, args.gff_style)
    # get_bed_file(name_gff_dict, args.input_file, args.output_file, args.extra_bp)
    query_name_dict = get_query_name(args.input_file)
    get_gff(query_name_dict, args.gff_file, args.gff_style)
    get_bed(query_name_dict, args.output_file, args.extra_bp)


if __name__ == '__main__':
    main()


