#! /usr/env/bin
import sys
import os


gff_file = sys.argv[1]
out_file = sys.argv[2]
gff_type = sys.argv[3].upper()
if gff_type == 'ENSEMBLE':
    SELECTED_TYPES = {'mRNA'}  # col3 : col8 ID_type
    SELECTED_SUB_TYPES = {'CDS', 'intron', 'five_prime_UTR', 'three_prime_UTR'}
    SELECTED_ATTRIBUTES = {'mRNA':['ID','Name'], 'CDS':['Parent'], 'intron':['Parent'], 'five_prime_UTR':['Parent'], 'three_prime_UTR':['Parent'], 'exon':['Parent']}
    ATRB_ID_SEP = ':'  # e.g. ID=transcript:ENST00000423372;Parent=gene:ENSG00000237683
elif gff_type == 'NCBI':
    SELECTED_TYPES = {'mRNA'}
    SELECTED_SUB_TYPES = {'CDS', 'intron', 'five_prime_UTR', 'three_prime_UTR'}
    SELECTED_ATTRIBUTES = ['ID', 'Parent']
    ATRB_ID_SEP = '-' # e.g. ID=exon-NM_001005277.1-1;Parent=rna-NM_001005277.1
# assembly_report_file = sys.argv[4]


# class GffTrans:
#     """
#     Tran
#     """

class GffAttributes:
    """
    attributes for gff record
    """
    def __init__(self,  attributes):
        self.string = attributes
        self.atrb_dict = {}
        if self.string != '':
            self._lst = attributes.split(';')
            for _atrb in self._lst:
                _key, _val = _atrb.split('=')
                self.atrb_dict[_key] = _val
    
    def __str__(self):
        return ';'.join([f'{k}={v}'for k, v in self.atrb_dict.items()])
    
    def trimm_atrb(self, feature):
        if feature not in SELECTED_ATTRIBUTES:
            return GffAttributes('')
        res = []
        for k in SELECTED_ATTRIBUTES[feature]:
            if k in self.atrb_dict:
                v = self.atrb_dict[k]
                idx = v.find(ATRB_ID_SEP)
                if idx != -1:
                    v = v[idx + 1:]
                res.append(f'{k}={v}')
        return GffAttributes(';'.join(res))


class GffRecord:
    """
    gff record got 8 cols:
    seqid source type start end score strand phase attributes
    """
    def __init__(self, records):
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
        self.attributes = GffAttributes(self.records[8])
        self.id = self.attributes.atrb_dict.get('ID', None)
        self.parent = self.attributes.atrb_dict.get('Parent', None)
        self.name = self.attributes.atrb_dict.get('Name', None)
        
    
    def __str__(self):
        return self.string
    
    def __repr__(self):
        return self.string
    
    def convert_id(self):
        return 'chr' + self.seqid
    
    def convert(self):
        return  GffRecord('\t'.join([self.convert_id(), self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes.trimm_atrb(self.type).string]))


limits = 10000
line_counts = 0
seq_counts = 0
gff_records = []
pre_id = None
output_fhand = open(out_file, 'w')
pre_exon_record = None
with open(gff_file) as fhand:
    
    for line in fhand:
        if not line or line.startswith('#'):
            continue
        curr_record = GffRecord(line)
        # print(pre_id, curr_record)
        if curr_record.id and (curr_record.type in SELECTED_TYPES):
            pre_id = curr_record.id
            curr_record = curr_record.convert()
            # print(curr_record, gff_records)
            if gff_records:
                for _record in gff_records:
                    gff_records.sort(key = lambda x : int(x.start))
                    output_fhand.write(_record.string + '\n')
            gff_records = []
            if seq_counts != 0:
                output_fhand.write('\n')
            output_fhand.write(curr_record.string + '\n') # seperate line
            seq_counts += 1
            line_counts += 1
        elif curr_record.parent == pre_id:
            if curr_record.type in SELECTED_SUB_TYPES:
                curr_record = curr_record.convert()
                # print('selected sub feature', curr_record)
                gff_records.append(curr_record)
                line_counts += 1
            elif curr_record.type == 'exon':
                curr_record = curr_record.convert()
                # print('pre:', pre_exon_record)
                # print('curr:', curr_record)
                if pre_exon_record is not None and  pre_exon_record.type == 'exon' and curr_record.parent != None and curr_record.parent == pre_exon_record.parent:  # add intron record
                    _start = str(int(pre_exon_record.end) + 1)
                    _end = str(int(curr_record.start) - 1)
                    _atrb = f'Parent={curr_record.parent}'
                    # print('start', 'end', 'atrb', _start, _end, _atrb, _end >= _start)
                    if int(_end) >= int(_start):  # do not forget use int for comparision
                        # output_fhand.write('\t'.join([curr_record.seqid, curr_record.source, 'intron', _start, _end, '.', curr_record.strand, '.', _atrb]) + '\n')
                        gff_records.append(GffRecord('\t'.join([curr_record.seqid, curr_record.source, 'intron', _start, _end, '.', curr_record.strand, '.', _atrb])))
                pre_exon_record = curr_record
                line_counts += 1
    # if line_counts > limits:
    #     break

if gff_records:
    gff_records.sort(key = lambda x : int(x.start))
    for _record in gff_records:
        # print(_record.string)
        output_fhand.write(_record.string + '\n')
output_fhand.write('\n')
print (f"Written total {seq_counts} seqs {line_counts} records")
output_fhand.close()
