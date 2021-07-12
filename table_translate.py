#!/usr/bin/env python3

"""
table_translate.py [-h] -i input.tsv -o output.tsv -r translation_table [-t trans_type]

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

Test Data:
    python3 table_translate.py -i test/ncbi_name_id.tsv -o test/trans_ncbi_name_id.tsv  -r data/homo_sapiens_gene2ensembl -t n2e 
    python3 table_translate.py -i test/ensembl_name_id.tsv -o test/trans_ensembl_name_id.tsv -r data/homo_sapiens_gene2ensembl -t e2n  

Real Data:
    python3 table_translate.py -i gff_out/ncbi_name_id.tsv -o gff_out/trans_ncbi_name_id.tsv  -r data/homo_sapiens_gene2ensembl -t n2e 
    python3 table_translate.py -i gff_out/ensembl_name_id.tsv -o gff_out/trans_ensembl_name_id.tsv -r data/homo_sapiens_gene2ensembl -t e2n  
"""


import argparse
from typing import NamedTuple, TextIO

class Arags(NamedTuple):
    """
    Command Line arguments
    """
    input_file: TextIO
    output_file: TextIO
    ref_table: TextIO
    trans_type: str


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Translate NCBI Refseq ID to ENSEMBL Transcript ID, or vice versa, based on a transformation table provided by NCBI", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input',
                        metavar='input.tsv',
                        type=argparse.FileType('r'),
                        help='Input tsv file, 3 cols of NAME/NCBI/ENSEML',
                        required=True
                        )
    parser.add_argument('-o', '--output',
                        metavar='output.tsv',
                        type=argparse.FileType('w'),
                        help='Output tsv file, 3 cols of NAME/NCBI/ENSEML',
                        required=True
                        )
    parser.add_argument('-r', '--ref',
                        metavar='translation_table',
                        type=argparse.FileType('r'),
                        help='Translation Table',
                        required=True
                        )
    parser.add_argument('-t', '--trans_type',
                        metavar='trans_type',
                        choices=['e2n', 'n2e', 'both'],
                        help='Transformation Style, can be: e2n:ENSEMBL to NCBI; n2e: NCBI to ENSEMBL; both: e2n and n2e',
                        default='both'
                        )
    args = parser.parse_args()
    return Arags(input_file=args.input,
                 output_file=args.output,
                 ref_table=args.ref,
                 trans_type=args.trans_type
                 )


args = get_args()


def trimm_id(ref_id):
    """
    trimm the version number after ncbi refseq id or ensembl id
    """
    return ref_id.split('.')[0]


def get_dicts(ref_table):
    """
    get the dictonary from ref table
    cols of the gene2ENSEMBL file:
    #tax_id GeneID  Ensembl_gene_identifier RNA_nucleotide_accession.version        Ensembl_rna_identifier  protein_accession.version       Ensembl_protein_identifier
    """
    n2e_dict = {'None':'None'}  # ncbi to ensembl dictionary
    e2n_dict = {'None':'None'} # ensembl to ncbi dictionary
    for line in ref_table:
        if line.startswith('#'):
            continue
        lst = line.strip().split('\t')
        ncbi_id = trimm_id(lst[3])
        ensembl_id = trimm_id(lst[4])
        n2e_dict[ncbi_id] = ensembl_id
        e2n_dict[ensembl_id] = ncbi_id
    return n2e_dict, e2n_dict

def translate_ids(input_file, output_file, trans_type, n2e_dict, e2n_dict):
    """
    translate the ids in the input file if needed, and write it to the output file
    """
    # total line total quries and the quries got translated
    total_lines = 0
    qurey_lines = 0
    trans_lines = 0
    output_file.write('\t'.join(['Name', 'NCBI', 'ENSEMBL']) + '\n')  # header

    for line in input_file:
        total_lines += 1
        if not line or total_lines == 1:
            continue
        name, ncbi_id, ensembl_id = line.strip().split('\t')
        if trans_type in {'n2e', 'both'} and ncbi_id != 'None' and ensembl_id == 'None':
            ensembl_id = n2e_dict.get(ncbi_id, 'None')
            qurey_lines += 1
            trans_lines += ensembl_id != 'None'
        if trans_type in {'e2n', 'both'} and ensembl_id != 'None' and ncbi_id == 'None':
            ncbi_id = e2n_dict.get(ensembl_id, 'None')
            qurey_lines += 1
            trans_lines += ncbi_id != 'None'
        output_file.write('\t'.join([name, ncbi_id, ensembl_id]) + '\n')

    print(f'There are {total_lines - 1} lines in the file,  {qurey_lines} quries and {trans_lines} got trans_ids.')

def main():
    """
    Main
    """
    n2e_dict, e2n_dict = get_dicts(args.ref_table)
    translate_ids(args.input_file, args.output_file, args.trans_type, n2e_dict, e2n_dict)
    args.output_file.write('\n')
    args.input_file.close()
    args.ref_table.close()
    args.output_file.close()

if __name__ == "__main__":
    main()
