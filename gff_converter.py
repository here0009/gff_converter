#!/usr/bin/env python3
"""
Input:
Output:
Usage:
    gff_converter.py [-h] -i File -o File [-a File] -s str [--add_intron]
Example:
    python3 gff_converter.py -i ensembl_test.gff3 -o ensembl_test_out2.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -s ENSEMBL --add_intron
    python3 gff_converter.py -i ensembl_test.gff3 -o ensembl_test_out2.gff -s ENSEMBL --add_intron
    python3 gff_converter.py ncbi_test.gff ncbi_test_out.gff data/GCF_000001405.25_GRCh37.p13_assembly_report.txt NCBI
"""


import argparse
from typing import NamedTuple, TextIO
from GFF import GffRecord, GffAttributes, SELECTED_TYPES, SELECTED_SUB_TYPES, SELECTED_ATTRIBUTES


class Arags(NamedTuple):
    """
    Command Line arguments
    """
    gff_file:TextIO
    out_file:TextIO
    assembly_report_file:TextIO
    gff_style:str
    add_intron:bool


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(description="Convert NCBI/ENSEMBL style GFF file to a concise UCSC style", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='File', type=argparse.FileType('r'), help='Input GFF File', required=True)
    parser.add_argument('-o', '--output', metavar='File', type=argparse.FileType('w'), help='Output File', required=True)
    parser.add_argument('-a', '--assembly', metavar='File', type=argparse.FileType('r'), help='Assembly Report File', default='data/GCF_000001405.25_GRCh37.p13_assembly_report.txt')
    parser.add_argument('-s', '--style', metavar='str', choices=['ENSEMBL', 'NCBI'], help='GFF style', required=True)
    parser.add_argument('--add_intron', help='Add intron to gff file', action='store_true')
    args = parser.parse_args()
    return Arags(gff_file=args.input, out_file=args.output, assembly_report_file=args.assembly, gff_style=args.style, add_intron=args.add_intron)


def seqid_conversion(assembly_report_file, gff_style):
    """
    return seqid_conversion dict based on assembly report file
    gff_style includs ESEMBLE AND NCBI
    """
    chroms = set([str(i) for i in range(1, 24)]) | set(['X', 'Y', 'MT'])
    refseq_to_ucsc = dict()
    genbank_to_ucsc = dict()
    for line in assembly_report_file:
        if line.startswith('#'):
            continue
        lst = line.strip().split('\t')
        seqid, refseq, genbank, ucsc = lst[0], lst[4], lst[6], lst[9]
        if seqid in chroms:
            genbank_to_ucsc[seqid] = ucsc
        else:
            genbank_to_ucsc[refseq] = ucsc
        refseq_to_ucsc[genbank] = ucsc
    if gff_style == 'ENSEMBL':
        return genbank_to_ucsc
    elif gff_style == 'NCBI':
        return refseq_to_ucsc

args = get_args()
seqid_trans_dict = seqid_conversion(args.assembly_report_file, args.gff_style)
line_counts = 0
seq_counts = 0
gff_records = []
pre_id = None
pre_exon_record = None

for line in args.gff_file:
    if not line or line.startswith('#'):
        continue
    curr_record = GffRecord(line, args.gff_style)
    if curr_record.id and (curr_record.type in SELECTED_TYPES):
        pre_id = curr_record.id
        curr_record = curr_record.convert(seqid_trans_dict)
        if gff_records:
            for _record in gff_records:
                gff_records.sort(key = lambda x : int(x.start))
                args.out_file.write(_record.string + '\n')
                line_counts += 1
        gff_records = []
        if seq_counts != 0:
            args.out_file.write('\n')
        args.out_file.write(curr_record.string + '\n') # seperate line
        seq_counts += 1
        line_counts += 1
    elif curr_record.parent == pre_id:
        if curr_record.type in SELECTED_SUB_TYPES:
            curr_record = curr_record.convert(seqid_trans_dict)
            gff_records.append(curr_record)

        elif args.add_intron and curr_record.type == 'exon':
            curr_record = curr_record.convert(seqid_trans_dict)
            if pre_exon_record is not None and  pre_exon_record.type == 'exon' and curr_record.parent != None and curr_record.parent == pre_exon_record.parent:  # add intron record
                _start = str(int(pre_exon_record.end) + 1)
                _end = str(int(curr_record.start) - 1)
                _atrb = f'Parent={curr_record.parent}'
                if int(_end) >= int(_start):  # do not forget use int for comparision
                    gff_records.append(GffRecord('\t'.join([curr_record.seqid, curr_record.source, 'intron', _start, _end, '.', curr_record.strand, '.', _atrb]), args.gff_style))
            pre_exon_record = curr_record

if gff_records:
    gff_records.sort(key = lambda x : int(x.start))
    for _record in gff_records:
        args.out_file.write(_record.string + '\n')
        line_counts += 1

args.out_file.write('\n')
args.assembly_report_file.close()
args.gff_file.close()
args.out_file.close()
print (f"Written total {seq_counts} seqs and {line_counts} lines")
