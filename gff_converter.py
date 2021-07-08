#!/usr/bin/env python3
"""
Input:
Output:
Usage:
    python3 gff_converter.py <input_file> <output_file> <assembly report> <gff_file_type>
Example: 
    python3 gff_converter.py ncbi_test.gff ncbi_test_out.gff data/GCF_000001405.25_GRCh37.p13_assembly_report.txt NCBI
"""
import sys
from GFF import GffRecord, GffAttributes

gff_file = sys.argv[1]
out_file = sys.argv[2]
assembly_report_file = sys.argv[3]
gff_type = sys.argv[4].upper()
assert gff_type in {'ENSEMBL', 'NCBI'}
NO_INTRON = True  # did the gff file include intron
SELECTED_TYPES = {'mRNA':'mRNA', 'transcript':'transcript'}  # col3 : col8 ID_type
SELECTED_SUB_TYPES = {'CDS':'CDS', 'intron':'intron', 'five_prime_UTR':'5_UTR', 'three_prime_UTR':'3_UTR'}
SELECTED_ATTRIBUTES = {'mRNA':['ID','Name', 'Parent'], 'transcript':['ID','Name', 'Parent'], 'CDS':['Parent'], 'intron':['Parent'], 'five_prime_UTR':['Parent'], 'three_prime_UTR':['Parent'], 'exon':['Parent']}
if gff_type == 'ENSEMBL':
    ATRB_ID_SEP = ':'  # e.g. ID=transcript:ENST00000423372;Parent=gene:ENSG00000237683
elif gff_type == 'NCBI':
    ATRB_ID_SEP = '-' # e.g. ID=exon-NM_001005277.1-1;Parent=rna-NM_001005277.1


class GffAttributes:
    """
    attributes for gff record
    """
    def __init__(self,  attributes, gff_type):
        self.string = attributes
        self.gff_type = gff_type
        self.atrb_dict = {}
        if self.string != '':
            self._lst = attributes.split(';')
            for _atrb in self._lst:
                if '=' not in _atrb:
                    continue
                _key, _val = _atrb.split('=')
                self.atrb_dict[_key] = _val
    
    def __str__(self):
        return ';'.join([f'{k}={v}'for k, v in self.atrb_dict.items()])
    
    def trimm_atrb(self, feature):
        if feature not in SELECTED_ATTRIBUTES:
            return GffAttributes('', gff_type)
        res = []
        if self.gff_type == 'NCBI' and feature in SELECTED_TYPES: # use the gene name if got one
            if  'Parent' != 'None':
                self.atrb_dict['Name'] = self.atrb_dict['Parent']
                del self.atrb_dict['Parent']

        for k in SELECTED_ATTRIBUTES[feature]:
            if k in self.atrb_dict:
                v = self.atrb_dict[k]
                idx = v.find(ATRB_ID_SEP)
                if idx != -1:
                    v = v[idx + 1:]
                res.append(f'{k}={v}')
        return GffAttributes(';'.join(res), gff_type)


class GffRecord:
    """
    gff record got 8 cols:
    seqid source type start end score strand phase attributes
    """
    def __init__(self, records, gff_type):
        self.string = records
        self.records = records.strip().split('\t')
        self.seqid = self.records[0]
        self.source = gff_type
        self.type = self.records[2]
        self.start = self.records[3]
        self.end = self.records[4]
        self.score = self.records[5]
        self.strand = self.records[6]
        self.phase = self.records[7] 
        self.attributes = GffAttributes(self.records[8], gff_type)
        self.id = self.attributes.atrb_dict.get('ID', None)
        self.parent = self.attributes.atrb_dict.get('Parent', None)
        self.name = self.attributes.atrb_dict.get('Name', None)

    def __str__(self):
        return self.string
    
    def __repr__(self):
        return self.string
        
    def convert(self):
        convert_seqid = seqid_trans_dict[self.seqid]
        convert_type = self.type
        convert_atrb = self.attributes.trimm_atrb(self.type).string
        if self.type in SELECTED_TYPES:
            convert_type = SELECTED_TYPES[self.type]
        elif self.type in SELECTED_SUB_TYPES:
            convert_type = SELECTED_SUB_TYPES[self.type]
        return  GffRecord('\t'.join([convert_seqid, self.source, convert_type, self.start, self.end, self.score, self.strand, self.phase, convert_atrb]), gff_type)


def seqid_conversion(assembly_report_file, gff_type):
    """
    return seqid_conversion dict based on assembly report file
    gff_type includs ESEMBLE AND NCBI
    """
    input_fhand = open(assembly_report_file)
    chroms = set([str(i) for i in range(1, 24)]) | set(['X', 'Y', 'MT'])
    refseq_to_ucsc = dict()
    genbank_to_ucsc = dict()
    for line in input_fhand:
        if line.startswith('#'):
            continue
        lst = line.strip().split('\t')
        seqid, refseq, genbank, ucsc = lst[0], lst[4], lst[6], lst[9]
        if seqid in chroms:
            genbank_to_ucsc[seqid] = ucsc
        else:
            genbank_to_ucsc[refseq] = ucsc
        refseq_to_ucsc[genbank] = ucsc
    if gff_type == 'ENSEMBL':
        return genbank_to_ucsc
    elif gff_type == 'NCBI':
        return refseq_to_ucsc


seqid_trans_dict = seqid_conversion(assembly_report_file, gff_type)
line_counts = 0
seq_counts = 0
gff_records = []
pre_id = None
input_fhand = open(gff_file)
output_fhand = open(out_file, 'w')
pre_exon_record = None

for line in input_fhand:
    if not line or line.startswith('#'):
        continue
    curr_record = GffRecord(line, gff_type)
    if curr_record.id and (curr_record.type in SELECTED_TYPES):
        pre_id = curr_record.id
        curr_record = curr_record.convert()
        if gff_records:
            for _record in gff_records:
                gff_records.sort(key = lambda x : int(x.start))
                output_fhand.write(_record.string + '\n')
                line_counts += 1
        gff_records = []
        if seq_counts != 0:
            output_fhand.write('\n')
        output_fhand.write(curr_record.string + '\n') # seperate line
        seq_counts += 1
        line_counts += 1
    elif curr_record.parent == pre_id:
        if curr_record.type in SELECTED_SUB_TYPES:
            curr_record = curr_record.convert()
            gff_records.append(curr_record)

        elif NO_INTRON and curr_record.type == 'exon':
            curr_record = curr_record.convert()
            if pre_exon_record is not None and  pre_exon_record.type == 'exon' and curr_record.parent != None and curr_record.parent == pre_exon_record.parent:  # add intron record
                _start = str(int(pre_exon_record.end) + 1)
                _end = str(int(curr_record.start) - 1)
                _atrb = f'Parent={curr_record.parent}'
                if int(_end) >= int(_start):  # do not forget use int for comparision
                    gff_records.append(GffRecord('\t'.join([curr_record.seqid, curr_record.source, 'intron', _start, _end, '.', curr_record.strand, '.', _atrb]), gff_type))
            pre_exon_record = curr_record

if gff_records:
    gff_records.sort(key = lambda x : int(x.start))
    for _record in gff_records:
        output_fhand.write(_record.string + '\n')
        line_counts += 1

output_fhand.write('\n')
output_fhand.close()
input_fhand.close()
print (f"Written total {seq_counts} seqs and {line_counts} lines")


