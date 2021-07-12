#!/usr/bin/env python3
"""
usage: gff_converter.py [-h] -i input_gff -o output_gff [-a assembly_report] -t output_name_id_table -s gff_style [--add_intron] [--add_utr]

Convert NCBI/ENSEMBL style Gff file to a concise UCSC style

optional arguments:
  -h, --help            show this help message and exit
  -i input_gff, --input input_gff
                        Input Gff File (default: None)
  -o output_gff, --output output_gff
                        Output Gff File (default: None)
  -a assembly_report, --assembly assembly_report
                        Assembly Report File (default: data/GCF_000001405.25_GRCh37.p13_assembly_report.txt)
  -t output_name_id_table, --table output_name_id_table
                        Output Name Id table File (default: None)
  -s gff_style, --style gff_style
                        Gff style, "ENSEMBL" or "NCBI" (default: None)
  --add_intron          Add intron to gff file (default: False)
  --add_utr             Add UTR to gff file (default: False)

Example:

Test Data:
    python3 gff_converter.py -i test/ncbi_test.gff -o test/ncbi_test_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t test/ncbi_name_id.tsv -s NCBI  --add_intron --add_utr

    python3 gff_converter.py -i test/ensembl_test.gff3 -o test/ensembl_test_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t test/ensembl_name_id.tsv -s ENSEMBL  --add_intron --add_utr

Real Data:
    python3 gff_converter.py -i gff_data/GCF_000001405.25_GRCh37.p13_genomic.gff -o gff_out/NCBI_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t  gff_out/ncbi_name_id.tsv -s NCBI --add_intron --add_utr

    python3 gff_converter.py -i gff_data/Homo_sapiens.GRCh37.87.Ensembl.gff3  -o gff_out/ENSEMBL_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t  gff_out/ensembl_name_id.tsv -s ENSEMBL --add_intron --add_utr
"""


import argparse
from typing import NamedTuple, TextIO
from collections import OrderedDict
from GFF import RNA, GffRecord, TransTables, SELECTED_TYPES, SELECTED_SUB_TYPES, WRITE_SUB_TYPES, TYPE_TRANS_DICT, ATRB_TRANS_DICT

class Arags(NamedTuple):
    """
    Command Line arguments
    """
    gff_file:TextIO
    out_file:TextIO
    assembly_report_file:TextIO
    name_id_table:TextIO
    gff_style:str
    add_intron:bool
    add_utr:bool


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Convert NCBI/ENSEMBL style Gff file to a concise UCSC style", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input',
                        metavar='input_gff',
                        type=argparse.FileType('r'),
                        help='Input Gff File',
                        required=True
                        )

    parser.add_argument('-o', '--output',
                        metavar='output_gff',
                        type=argparse.FileType('w'),
                        help='Output Gff File',
                        required=True
                        )
    parser.add_argument('-a', '--assembly', 
                        metavar='assembly_report', 
                        type=argparse.FileType('r'), 
                        help='Assembly Report File', 
                        required=True
                        )
    parser.add_argument('-t', '--table',
                        metavar='output_name_id_table',
                        type=argparse.FileType('w'),
                        help='Output Name Id table File',
                        required=True
                        )
    parser.add_argument('-s', '--style',
                        metavar='gff_style',
                        choices=['ENSEMBL', 'NCBI'],
                        help='Gff style, "ENSEMBL" or "NCBI"',
                        required=True
                        )
    parser.add_argument('--add_intron',
                        help='Add intron to gff file',
                        action='store_true'
                        )
    parser.add_argument('--add_utr',
                        help='Add UTR to gff file',
                        action='store_true'
                        )
    args = parser.parse_args()
    return Arags(gff_file=args.input, 
                out_file=args.output, 
                assembly_report_file=args.assembly, 
                name_id_table=args.table,
                gff_style=args.style, 
                add_intron=args.add_intron, 
                add_utr=args.add_utr
                )


def get_seqid_trans_table(assembly_report_file):
    """
    return seqid_conversion dict based on assembly report file
    """
    chroms = set([str(i) for i in range(1, 24)]) | set(['X', 'Y', 'MT'])
    ncbi_to_ucsc = dict()
    ensembl_to_ucsc = dict()
    for line in assembly_report_file:
        if line.startswith('#'):
            continue
        lst = line.strip().split('\t')
        seqid, refseq, genbank, ucsc = lst[0], lst[4], lst[6], lst[9]
        if seqid in chroms:
            ensembl_to_ucsc[seqid] = ucsc
        else:
            ensembl_to_ucsc[refseq] = ucsc
        ncbi_to_ucsc[genbank] = ucsc
    return ncbi_to_ucsc, ensembl_to_ucsc


def get_trans_tables(seqid_trans_dict, type_trans_dict, atrb_trans_dict):
    """
    Get the tables for name transision
    """
    trans_tables = TransTables()
    trans_tables.add_table('seqid', seqid_trans_dict)
    trans_tables.add_table('type', type_trans_dict)
    trans_tables.add_table('atrb', atrb_trans_dict)
    return trans_tables

def read_gff(gff_file, gff_style, trans_tables):
    """
    1. For each gff record, if it belongs to selected types or subtypes, do the transformation and trimming
    2. read selected_types into an OrderDict, selected_sub_types into a list
    3. add selected_sub_type if its parent is in the OrderDict
    4. return the OrderDict
    """
    id_gff_dict = OrderedDict([])  # gff_record.id : rna
    subtype_gff_list = []

    for line in gff_file:  # add selected type to id_list and id_gff_dict
        if not line or line.startswith('#'):
            continue
        curr_record = GffRecord(line, gff_style)
        if curr_record.type in SELECTED_TYPES:
            curr_record = RNA(curr_record.convert(
                trans_tables).string, gff_style)
            if curr_record.id is not None:
                id_gff_dict[curr_record.id] = curr_record
        elif curr_record.type in SELECTED_SUB_TYPES and curr_record.parent is not None:
            curr_record = curr_record.convert(trans_tables)
            subtype_gff_list.append(curr_record)

    for gff in subtype_gff_list:  # add subtype as children if the parent is in id_list
        if gff.parent in id_gff_dict:
            id_gff_dict[gff.parent].add_children(gff)
    
    return id_gff_dict


def write_table(gff_style, name_id_table, gff):
    """
    Write name id table file, the head is:
    Name   NCBI   ENSEMBL
    """
    name = gff.name if gff.name else 'None'
    source_id = gff.id if gff.id else 'None'
    if gff_style == "NCBI":
        name_id_table.write('\t'.join([name, source_id, 'None']) + '\n')
    else:
        name_id_table.write('\t'.join([name, 'None', source_id]) + '\n')

def write_gff(out_file, add_intron, add_utr, gff_style, name_id_table, id_gff_dict):
    """
    write gff record to the output file
    """
    line_counts = 0
    selected_type_counts = 0
    for rna in id_gff_dict.values():
        out_file.write(rna.string + ';\n')
        write_table(gff_style, name_id_table, rna)
        selected_type_counts += 1
        sub_gff_list = []
        if 'intron' not in rna.children and add_intron:
            rna.add_intron()
        if add_utr and '5_UTR' not in rna.children and '3_UTR' not in rna.children:
            rna.add_utrs()
        for sub_type in WRITE_SUB_TYPES:
            sub_gff_list.extend(rna.get_children(sub_type))
        sub_gff_list.sort(key=lambda x: int(x.start))
        for _gff in sub_gff_list:
            out_file.write(_gff.string + ';\n')
            line_counts += 1
        out_file.write('\n')
    print(f"Written {selected_type_counts} seqs and {line_counts} lines")

def main():
    """
    Main
    """
    args = get_args()
    seqid_dicts = get_seqid_trans_table(args.assembly_report_file)
    if args.gff_style == 'NCBI':
        seqid_trans_dict = seqid_dicts[0]
    elif args.gff_style == 'ENSEMBL':
        seqid_trans_dict = seqid_dicts[1]
    trans_tables = get_trans_tables(
        seqid_trans_dict, TYPE_TRANS_DICT, ATRB_TRANS_DICT)
    id_rna_dict = read_gff(args.gff_file, args.gff_style, trans_tables.tables)
    args.name_id_table.write('\t'.join(['Name', 'NCBI', 'ENSEMBL']) + '\n') # header for name_id_table
    write_gff(args.out_file, args.add_intron, args.add_utr,
              args.gff_style, args.name_id_table, id_rna_dict)
    args.assembly_report_file.close()
    args.gff_file.close()
    args.out_file.close()
    args.name_id_table.close()

if __name__ == "__main__":
    main()
