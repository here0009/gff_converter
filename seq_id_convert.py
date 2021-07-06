#!/usr/bin/env python3
import sys
assembly_report_file = sys.argv[1]
input_fhand = open(assembly_report_file)
chroms = set([str(i) for i in range(1, 24)]) | set(['X', 'Y', 'MT'])
refseq_to_ucsc = dict()
genbank_to_ucsc = dict()
for line in input_fhand:
    if line.startswith('#'):
        continue
    lst = line.split('\t')
    seqid, refseq, genbank, ucsc = lst[0], lst[4], lst[6], lst[8]
    if seqid in chroms:
        refseq_to_ucsc[seqid] = ucsc
    else:
        refseq_to_ucsc[refseq] = ucsc
    genbank_to_ucsc[genbank] = ucsc
for k, v in refseq_to_ucsc.items():
    print(k, v)
for k, v in genbank_to_ucsc.items():
    print(k, v)
