#!/usr/bin/env python3
"""
usage: get_bed_file.py [-h] -i input_file -o output_bed -g gff_file [-s gff_style] [-e extra_bp]

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
                        Gff style, "ENSEMBL" or "NCBI" (default: NCBI)
  -e extra_bp, --extra_bp extra_bp

Example:
    python3 get_bed_file.py -i get_coord/gene_list.bed -o get_coord/refGene3.bed -g hg19_gff/UCSC_refGene_output.gff -e 10
    python3 get_bed_file.py -i get_coord/gene_list2.bed -o get_coord/ncbiRefSeq.bed -g hg19_gff/UCSC_ncbiRefSeq_output.gff -e 10

Note:
    The gff file used for get coordinates is downloaded from UCSC.
    http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/
"""


import argparse
from typing import NamedTuple, OrderedDict, TextIO
from GFF import RNA

# CHOOSEN_SUB_TYPES = {'3_UTR', '5_UTR'}


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
                        choices=['ENSEMBL', 'NCBI'],
                        help='Gff style, "ENSEMBL" or "NCBI"',
                        default='NCBI'
                        )
    parser.add_argument('-e', '--extra_bp',
                        metavar='extra_bp',
                        type=int,
                        default=0)
    args = parser.parse_args()
    return Arags(input_file=args.input_file, 
                output_file=args.output_file, 
                gff_file=args.gff_file, 
                gff_style=args.style,
                extra_bp=args.extra_bp
                )


def trimm_id(ref_id):
    """
    trimm the version number after ncbi refseq id or ensembl id
    """
    return ref_id.split('.')[0]

def get_bed_record(gff_record, name, extra_bp):
    """
    Get the bed record based on gff_record
    """
    bed_record = '\t'.join([gff_record.seqid, 
                                    str(int(gff_record.start) - extra_bp), 
                                    str(int(gff_record.end) + extra_bp + 1), 
                                    gff_record.strand, 
                                    name,
                                    trimm_id(gff_record.parent)])
    return bed_record

def get_query_id(input_fhand):
    """
    Get the gene in an orderdict
    Input file format
    #chrom	start	end	strand	name	rna_id
    chr10	27035515	27150027	-	ABI1	NM_005470
    """
    query_counts = 0
    query_id_dict = OrderedDict()
    for line in input_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        query_counts += 1
        lst = line.split('\t')
        gene_name, gene_id = lst[4], lst[5]  # input gene name may be lower case
        query_id_dict[gene_id] = gene_name
    return query_id_dict

def get_gff(query_id_dict, gff_fhand, gff_style, CHOOSEN_TYPES, CHOOSEN_SUB_TYPES):
    """
    Get the gff file
    """
    id_gff_dict = dict() # id:gff_record
    children_list = []
    for line in gff_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        gff_record = RNA(line, gff_style)
        if gff_record.type in CHOOSEN_TYPES and gff_record.id in query_id_dict:
            id_gff_dict[gff_record.id] = gff_record
        # elif gff_record.type == 'CDS' and gff_record.parent in query_id_dict:
        elif gff_record.type in CHOOSEN_SUB_TYPES and gff_record.parent in query_id_dict:
            children_list.append(gff_record)
    # print(id_gff_dict)
    # print(children_list)
    for gff_cds in children_list:
        id_gff_dict[gff_cds.parent].add_children(gff_cds)
    return id_gff_dict


def get_bed(query_id_dict, id_gff_dict, output_fhand, extra_bp, OUTPUT_BED_SUB_TYPES):
    """
    Get the bed file
    """
    bed_header = '\t'.join(['#chrom', 'start', 'end', 'strand', 'name', 'rna_id'])
    output_fhand.write(bed_header + '\n')
    not_found_names = []
    no_subtype_names = []
    for id in query_id_dict:
        if id not in id_gff_dict:
            not_found_names.append((id, query_id_dict[id]))
            continue
        gff_record = id_gff_dict[id]
        children_list = []
        for _type in OUTPUT_BED_SUB_TYPES:
            children_list.extend(gff_record.get_children(_type))
        if len(children_list) == 0:
            no_subtype_names.append((gff_record.id, gff_record.name))
        children_list.sort(key=lambda x: int(x.start))
        for child_record in children_list:
            bed_record = get_bed_record(child_record, gff_record.name, extra_bp)
            output_fhand.write(bed_record + '\n')

    not_found_counts = len(not_found_names)
    total_queries = len(query_id_dict)
    no_subtype_counts = len(no_subtype_names)
    print(f"There are {total_queries} total_queries")
    print(f'{not_found_counts} not found, the gene names not found are:')
    print(not_found_names)
    print(f'{no_subtype_counts} no subtypes, the gene names without subtypes are:')
    print(no_subtype_names)


def main():
    CHOOSEN_TYPES = {'transcript', 'mRNA'}
    CHOOSEN_SUB_TYPES = {'CDS'}
    OUTPUT_BED_SUB_TYPES = {'CDS'}
    args = get_args()
    query_id_dict = get_query_id(args.input_file)
    # print(query_id_dict)
    id_gff_dict = get_gff(query_id_dict, args.gff_file, args.gff_style, CHOOSEN_TYPES, CHOOSEN_SUB_TYPES)
    # print(id_gff_dict)
    get_bed(query_id_dict, id_gff_dict, args.output_file, args.extra_bp, OUTPUT_BED_SUB_TYPES)


if __name__ == '__main__':
    main()




