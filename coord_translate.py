#!/usr/bin/env python3

"""
usage: table_translate.py [-h] -i input.tsv -o output.tsv [-r File] [-t str]

Translate NCBI Refseq ID to ENSEMBL Transcript ID, or vice versa, based on a transformation table provided by NCBI    

optional arguments:
  -h, --help            show this help message and exit
  -i input.tsv, --input input.tsv
                        Input tsv file, 3 cols of NAME/NCBI/ENSEML (default: None)
  -o output.tsv, --output output.tsv
                        Output tsv file, 3 cols of NAME/NCBI/ENSEML (default: None)
  -r File, --ref File   Translation Table (default: data/homo_sapiens_gene2ensembl)
  -t str, --trans_type str
                        Transformation Style, can be: e2n:ENSEMBL to NCBI; n2e: NCBI to ENSEMBL; both: e2n and n2e    
                        (default: both)
Example: 
   python3 coord_translate.py -n test/ncbi_test.gff -e test/ensembl_test.gff3 -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -o test/trans_coord_id.tsv
   python3 coord_translate.py -n gff_data/GCF_000001405.25_GRCh37.p13_genomic.gff -e gff_data/Homo_sapiens.GRCh37.87.Ensembl.gff3 -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -o trans_coord_id.tsv
"""


import argparse
from collections import defaultdict
from typing import NamedTuple, OrderedDict, TextIO
from GFF import GffRecord, RNA, SELECTED_TYPES, TYPE_TRANS_DICT, ATRB_TRANS_DICT
from gff_converter import get_seqid_trans_table, get_trans_tables, read_gff

SELECTED_SUB_TYPES = {'exon'}  # only select exons

class Arags(NamedTuple):
    """
    Command Line arguments
    """
    ensembl_file: TextIO
    ncbi_file: TextIO
    assembly_report_file: TextIO
    output_file: TextIO


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Match NCBI Refseq ID and ENSEMBL Transcript ID based on the coordinates. Transcripts/mRNAs that got the same seqid, start, end, strand and exones  are considered to be identical.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-n', '--ncbi_gff',
                        metavar='ncbi.gff',
                        type=argparse.FileType('r'),
                        help='NCBI Gff File',
                        required=True
                        )
    parser.add_argument('-e', '--ensembl_gff',
                        metavar='ensembl.gff',
                        type=argparse.FileType('r'),
                        help='ENSEMBL Gff File',
                        required=True
                        )
    parser.add_argument('-a', '--assembly_report',
                        metavar='assembly_report',
                        type=argparse.FileType('r'),
                        help='assembly report file',
                        default='data/GCF_000001405.25_GRCh37.p13_assembly_report.txt',
                        # required=True
                        )
    parser.add_argument('-o', '--output',
                        metavar='output.tsv',
                        type=argparse.FileType('w'),
                        help="output tsv File, cols: seqid;start;end;strand;ncbi_id;ncbi_name;ensembl_id;ensembl_name",
                        required=True
                        )
    args = parser.parse_args()
    return Arags(ensembl_file=args.ensembl_gff,
                 ncbi_file=args.ncbi_gff,
                 assembly_report_file=args.assembly_report,
                 output_file=args.output
                 )


def test_get_coord_dict(id_gff_dict, test_dict):
    """
    test how many mRNA/transcript got the same start, end and exon information
    """
    for rna in id_gff_dict.values():
        exons = sorted(rna.get_coordinates('exon'))
        relative_exon = []
        start = int(rna.start)
        for s, e in exons:
            relative_exon.extend([s - start, e -start])
        primary_key = '\t'.join([rna.seqid, rna.start, rna.end, rna.strand])
        exon_feature = ' '.join(str(pos) for pos in relative_exon)
        test_dict[(primary_key, exon_feature)].append(rna)
    for k, v in test_dict:
        if len(v) > 1:
            print(k, v)


def get_coord_dict(id_gff_dict, gff_style, coord_dict):
    """
    generate relative exon position as a feature of transcript/mRNA, set seqid, start, end, strand and exon as keys of transcript/mRNA.
    [ncbi_id, ncbi_name, ensembl_id, ensembl_name] as vals
    """
    for rna in id_gff_dict.values():
        exons = sorted(rna.get_coordinates('exon'))
        relative_exon = []
        start = int(rna.start)
        for s, e in exons:
            relative_exon.extend([s - start, e - start])
        primary_key = '\t'.join([rna.seqid, rna.start, rna.end, rna.strand])
        exon_coordinates = ' '.join(str(pos) for pos in relative_exon)
        key = tuple([primary_key, exon_coordinates])
        if key not in coord_dict:
            coord_dict[key] = ['None'] * 4
        if gff_style == 'NCBI':
            coord_dict[key][0] = rna.id if rna.id is not None else 'None'
            coord_dict[key][1] = rna.name if rna.name is not None else 'None'
        elif gff_style == 'ENSEMBL':
            coord_dict[key][2] = rna.id if rna.id is not None else 'None'
            coord_dict[key][3] = rna.name if rna.name is not None else 'None'


def write_coord_dict(coord_dict, output_file):
    """
    output table cols
    seqid start end strand ncbi_id ncbi_name ensembl_id ensembl_name
    """
    ncbi_id_counts = 0
    ensembl_id_counts = 0
    common_id_counts = 0
    line_counts = len(coord_dict)
    for key, val in coord_dict.items():
        ncbi_id_counts += val[0] != 'None'
        ensembl_id_counts += val[2] != 'None'
        common_id_counts += val[0] != 'None' and val[2] != 'None'
        output_file.write(key[0] + '\t' + '\t'.join(val) + '\n')
    print(f'There are {line_counts} lines, {ncbi_id_counts} ncbi_ids, {ensembl_id_counts} ensembl_ids, {common_id_counts} got common ids.')

def main():
    """
    Main
    """
    args = get_args()

    ncbi_seqid_trans_dict, ensembl_seqid_trans_dict = get_seqid_trans_table(
        args.assembly_report_file)
    ncbi_trans_tables = get_trans_tables(
        ncbi_seqid_trans_dict, TYPE_TRANS_DICT, ATRB_TRANS_DICT)
    ensembl_trans_tables = get_trans_tables(
        ensembl_seqid_trans_dict, TYPE_TRANS_DICT, ATRB_TRANS_DICT)
    
    # read mRNA and transcript into dictionary of {id:gff_record}
    # add its exon to the mRNA or transcript to generate relative exon position as a feature to differ mRNA/transcript that got the same start and end positions.
    ncbi_gff_dict = read_gff(args.ncbi_file, 'NCBI', ncbi_trans_tables.tables) 
    ensembl_gff_dict = read_gff(args.ensembl_file, 'ENSEMBL', ensembl_trans_tables.tables)
    
    # test how many rna/transcript got the same seqid, strand, start, end, and exons
    # test_dict = defaultdict(list)
    # test_get_coord_dict(ncbi_gff_dict, test_dict)
    # test_get_coord_dict(ensembl_gff_dict, test_dict)
    
    coord_dict = OrderedDict()
    get_coord_dict(ncbi_gff_dict, 'NCBI', coord_dict)
    get_coord_dict(ensembl_gff_dict, 'ENSEMBL', coord_dict)

    write_coord_dict(coord_dict, args.output_file)

    args.ncbi_file.close()
    args.ensembl_file.close()
    args.assembly_report_file.close()
    args.output_file.close()


if __name__ == "__main__":
    main()
