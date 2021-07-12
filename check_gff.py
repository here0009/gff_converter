#!/usr/bin/env python3
"""
To check if a gff file is valid as the format:
chr1	NCBI	mRNA	21877763	21904903	.	+	.	ID=NM_001369803; name=ALPL;
chr1	NCBI	5_UTR	21877763	21877921	.	+	.	Parent=NM_001369803;
chr1	NCBI	intron	21877922	21880470	.	+	.	Parent=NM_001369803;
chr1	NCBI	5_UTR	21880471	21880574	.	+	.	Parent=NM_001369803;
chr1	NCBI	CDS	21880575	21880635	.	+	0	Parent=NM_001369803;

1. every seq got its own id and name
2. every subseq got its parent, and is directy placed under its parent
3. if we join the coordinates of subseq end by end, its equal to the coordinates of the seq


usage: check_gff.py [-h] -i input.gff [-s gff_style]

Test if the gff file is of correct format

optional arguments:
  -h, --help            show this help message and exit
  -i input.gff, --input input.gff
                        Input Gff File (default: None)
  -s gff_style, --gff_style gff_style
                        The style of input gff file (default: None)

Example:

Test Data:
    python3 check_gff.py -i test/ncbi_test_output.gff -s NCBI
    python3 check_gff.py -i test/ensembl_test_output.gff -s ENSEMBL

Real Data:
    python3 check_gff.py -i gff_out/NCBI_output.gff -s NCBI
    python3 check_gff.py -i gff_out/ENSEMBL_output.gff -s ENSEMBL
"""

from GFF import GffRecord, SELECTED_SUB_TYPES, SELECTED_TYPES, RNA
import argparse
from typing import NamedTuple, TextIO
SELECTED_SUB_TYPES = {'CDS', '5_UTR', '3_UTR', 'intron'}


class Arags(NamedTuple):
    """
    Command Line arguments
    """
    input_file: TextIO
    gff_style: str


def get_args():
    """ 
    Get command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Test if the gff file is of correct format", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input',
                        metavar='input.gff',
                        type=argparse.FileType('r'),
                        help='Input Gff File',
                        required=True
                        )
    parser.add_argument('-s', '--gff_style',
                        metavar='gff_style',
                        choices=['UCSC', 'NCBI', 'ENSEMBL'],
                        help='The style of input gff file'
                        )
    args = parser.parse_args()
    return Arags(input_file=args.input,
                gff_style=args.gff_style
                )




def check_gff(input_fhand, gff_style):
    sub_coords = []
    rna = None
    line_counts = 0
    seq_counts = 0
    invalid_counts = 0
    for line in input_fhand:
        if not line or line.startswith('#'):
            continue
        if line == '\n': # blank line 
            # continue
            if rna and sub_coords:
                start, end = int(rna.start), int(rna.end)
                if rna.id == 'None' or rna.name == 'None':
                    invalid_counts += 1
                    print("The following seq do not got name or id:")
                    print(rna)
                if sub_coords[0][0] != start or sub_coords[-1][-1] != end:
                    invalid_counts += 1
                    print("The following rna and its subsequneces got inconsistent start or ending")
                    print(rna)
                    print(sub_coords)
                for i in range(len(sub_coords) - 1):
                    if sub_coords[i][1] + 1 != sub_coords[i + 1][0]:
                        invalid_counts += 1
                        print("The following rna got wrong subseq")
                        print(rna)
                        print(sub_coords)
                        print(sub_coords[i][1], sub_coords[i + 1][0])
            continue

        gff_record = GffRecord(line, gff_style)
        if gff_record.type in SELECTED_TYPES:
            rna = RNA(line, gff_style)
            sub_coords = []
            seq_counts += 1
        elif gff_record.type in SELECTED_SUB_TYPES:
            line_counts += 1
            gff = GffRecord(line, gff_style)
            if gff.parent != rna.id:
                print("The following subseq got wrong parent")
                print(gff)
                invalid_counts += 1
            else:
                sub_coords.append((int(gff.start), int(gff.end)))
    print(f"Checked {seq_counts} seqs and {line_counts} lines, {invalid_counts} of them is invalid")
            
def main():
    """
    Main
    """
    args = get_args()
    input_hand = args.input_file
    gff_style = args.gff_style
    check_gff(input_hand, gff_style)


if __name__ == "__main__":
    main()


