#!/usr/bin/env python3
"""
Classes and Parameters for GFF conversion
"""

SELECTED_TYPES = {'mRNA':'mRNA', 'transcript':'transcript'}  # col3 : col8 ID_type
SELECTED_SUB_TYPES = {'CDS':'CDS', 'intron':'intron', 'five_prime_UTR':'5_UTR', 'three_prime_UTR':'3_UTR'}
SELECTED_ATTRIBUTES = {'mRNA':['ID','Name', 'Parent'], 'transcript':['ID','Name', 'Parent'], 'CDS':['Parent'], 'intron':['Parent'], 'five_prime_UTR':['Parent'], 'three_prime_UTR':['Parent'], 'exon':['Parent']}

class GffAttributes:
    """
    attributes for gff record
    """
    def __init__(self,  attributes, gff_style):
        self.string = attributes
        self.gff_style = gff_style
        self.atrb_dict = {}
        if self.string != '':
            self._lst = attributes.split(';')
            for _atrb in self._lst:
                if '=' not in _atrb:
                    continue
                _key, _val = _atrb.split('=')
                self.atrb_dict[_key] = _val
        if self.gff_style == 'ENSEMBL':
            self.atrb_id_sep = ':'  # e.g. ID=transcript:ENST00000423372;Parent=gene:ENSG00000237683
        elif self.gff_style == 'NCBI':
            self.atrb_id_sep = '-' # e.g. ID=exon-NM_001005277.1-1;Parent=rna-NM_001005277.1
    
    def __str__(self):
        return ';'.join([f'{k}={v}'for k, v in self.atrb_dict.items()])
    
    def trimm_atrb(self, feature):
        if feature not in SELECTED_ATTRIBUTES:
            return GffAttributes('', self.gff_style)
        res = []
        if self.gff_style == 'NCBI' and feature in SELECTED_TYPES: # use the gene name if got one
            if  'Parent' != 'None':
                self.atrb_dict['Name'] = self.atrb_dict['Parent']
                del self.atrb_dict['Parent']
        for k in SELECTED_ATTRIBUTES[feature]:
            if k in self.atrb_dict:
                v = self.atrb_dict[k]
                idx = v.find(self.atrb_id_sep)
                if idx != -1:
                    v = v[idx + 1:]
                res.append(f'{k}={v}')
        return GffAttributes(';'.join(res), self.gff_style)


class GffRecord:
    """
    gff record got 8 cols:
    seqid source type start end score strand phase attributes
    """
    def __init__(self, records, gff_style):
        self.string = records
        self.records = records.strip().split('\t')
        self.seqid = self.records[0]
        self.source = gff_style
        self.type = self.records[2]
        self.start = self.records[3]
        self.end = self.records[4]
        self.score = self.records[5]
        self.strand = self.records[6]
        self.phase = self.records[7] 
        self.attributes = GffAttributes(self.records[8], gff_style)
        self.id = self.attributes.atrb_dict.get('ID', None)
        self.parent = self.attributes.atrb_dict.get('Parent', None)
        self.name = self.attributes.atrb_dict.get('Name', None)
        
    def __str__(self):
        return self.string
    
    def __repr__(self):
        return self.string
        
    def convert(self, seqid_trans_dict):
        convert_seqid = seqid_trans_dict[self.seqid]
        convert_type = self.type
        convert_atrb = self.attributes.trimm_atrb(self.type).string
        if self.type in SELECTED_TYPES:
            convert_type = SELECTED_TYPES[self.type]
        elif self.type in SELECTED_SUB_TYPES:
            convert_type = SELECTED_SUB_TYPES[self.type]
        return  GffRecord('\t'.join([convert_seqid, self.source, convert_type, self.start, self.end, self.score, self.strand, self.phase, convert_atrb]), self.source)
