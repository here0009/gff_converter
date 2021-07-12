#!/usr/bin/env python3
"""
Classes and Parameters for GFF conversion
"""

SELECTED_TYPES = {'mRNA', 'transcript'}  # col3 : col8 ID_type
SELECTED_SUB_TYPES = {'CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'intron'}
WRITE_SUB_TYPES = {'CDS',  '5_UTR', '3_UTR', 'intron'}
SELECTED_ATTRIBUTES =  {'mRNA':['ID','Name', 'Parent'], 'transcript':['ID','Name', 'Parent'], 'CDS':['Parent'], 'intron':['Parent'], 'five_prime_UTR':['Parent'], 'three_prime_UTR':['Parent'], 'exon':['Parent']}
WRITE_ATTRIBUTES = {'mRNA':['ID','Name'], 'transcript':['ID','Name'], 'CDS':['Parent'], 'intron':['Parent'], 'five_prime_UTR':['Parent'], 'three_prime_UTR':['Parent'], 'exon':['Parent']}
TYPE_TRANS_DICT = {'CDS': 'CDS',  'five_prime_UTR': '5_UTR', 'three_prime_UTR': '3_UTR',
                   'intron': 'intron', 'mRNA': 'mRNA', 'transcript': 'transcript', 'exon': 'exon'}
ATRB_TRANS_DICT = {'ID': 'ID', 'Name': 'name', 'Parent': 'Parent'}


class TransTables():
    """
    Store the dictionaries for name transformation
    """
    def __init__(self):
        self.tables = dict()

    def add_table(self, name, dictionary):
        self.tables[name] = dictionary

    def get_table(self, name):
        return self.tables.get(name, None)

class GffAttributes:
    """
    attributes for gff record
    """
    def __init__(self,  attributes, gff_style):
        self.string = attributes
        self.gff_style = gff_style
        self.atrb_dict = {}
        if self.string != '':
            self.lst = attributes.split(';')
            for atrb in self.lst:
                atrb = atrb.strip()  # remove extra spaces if any
                if '=' not in atrb:
                    continue
                key, val = atrb.split('=')
                self.atrb_dict[key] = val
        if self.gff_style == 'ENSEMBL':
            self.atrb_id_sep = ':'  # e.g. ID=transcript:ENST00000423372;Parent=gene:ENSG00000237683
        elif self.gff_style == 'NCBI':
            self.atrb_id_sep = '-' # e.g. ID=exon-NM_001005277.1-1;Parent=rna-NM_001005277.1
    
    def __str__(self):
        return self.string
    
    def __repr__(self):
        return self.string
    
    def trimm_atrb(self, type_name, atrb_trans_dict):
        # type_name is name of type or subtype
        if type_name not in SELECTED_ATTRIBUTES:
            return GffAttributes('', self.gff_style)
        res = []
        if self.gff_style == 'NCBI' and type_name in SELECTED_TYPES: # NCBI use "NM_*"  or "NR_* "as name, replace it with the gene name if got one
            if 'Parent' in self.atrb_dict:
                self.atrb_dict['Name'] = self.atrb_dict['Parent']
                del self.atrb_dict['Parent']
        for k in WRITE_ATTRIBUTES[type_name]:
            if k in self.atrb_dict:
                v = self.atrb_dict[k]
                idx = v.find(self.atrb_id_sep)
                if idx != -1:
                    v = v[idx + 1:]
                if k == 'ID' or k == 'Parent':  # remove the version after dot for id and parent
                    v = v.split('.')[0]
                res.append(f'{atrb_trans_dict[k]}={v}')  # transfer Name to name
        return GffAttributes('; '.join(res), self.gff_style)


class GffRecord:
    """
    gff record got 9 cols:
    seqid source type start end score strand phase attributes
    """
    def __init__(self, record, gff_style):
        self.string = record
        lst = record.strip().split('\t')
        assert len(lst) == 9
        self.seqid = lst[0]
        self.source = gff_style
        self.type = lst[2]
        self.start = lst[3]
        self.end = lst[4]
        self.score = lst[5]
        self.strand = lst[6]
        self.phase = lst[7] 
        self.attributes = GffAttributes(lst[8], gff_style)
        self.id = self.attributes.atrb_dict.get('ID', None)
        self.parent = self.attributes.atrb_dict.get('Parent', None)
        self.name = self.attributes.atrb_dict.get(
            'Name', None) or self.attributes.atrb_dict.get('name', None) # Attributes 'Name' will be converted to 'name'
        self.children = dict()
        
    def __str__(self):
        return self.string
    
    def __repr__(self):
        return self.string
    
    def add_children(self, gff_record):
        if gff_record.type not in self.children:
            self.children[gff_record.type] = []
        self.children[gff_record.type].append(gff_record)
    
    def get_children(self, sub_type):
        return self.children.get(sub_type, [])
    
    def get_coordinates(self, sub_type):
        return [(int(gff.start), int(gff.end)) for gff in self.get_children(sub_type)]
    
    def add_coordinates(self, sub_type, coordinates):
        """
        add coordinates as gff record to self.children[sub_type]
        """
        for start, end in coordinates:
            gff = GffRecord('\t'.join([self.seqid, self.source, sub_type, str(start), str(end), '.', self.strand, '.', f'Parent={self.id}']), self.source)
            self.add_children(gff)
        
    def convert(self, tables):
        """
        convert to seqid, type,  and trimm atrb according to tables
        """
        if 'seqid' in tables and self.seqid in tables['seqid']:
            seqid = tables['seqid'][self.seqid]
        else:
            seqid = self.seqid
        if 'type' in tables and self.type in tables['type']:
            convert_type = tables['type'][self.type]
        else:
            convert_type = self.type
        if 'atrb' in tables:
            atrb = self.attributes.trimm_atrb(self.type, tables['atrb']).string
        else:
            atrb = self.atrb.string
        return  GffRecord('\t'.join([seqid, self.source, convert_type, self.start, self.end, self.score, self.strand, self.phase, atrb]), self.source)


class RNA(GffRecord):
    """
    For mRNA, transcripts, we can add introns or utrs to it
    """ 
    def get_introns(self):
        introns = []
        exons = sorted(self.get_coordinates('exon'))
        for i in range(len(exons) - 1):
            if exons[i + 1][0] - 1 >= exons[i][1] + 1:
                introns.append((exons[i][1] + 1, exons[i + 1][0] - 1))
        return introns
                
    def get_utrs(self):
        left_utrs = []
        right_utrs = []
        cds = sorted(self.get_coordinates('CDS'))
        exons = sorted(self.get_coordinates('exon'))
        assert self.strand in '-+'
        if cds and exons:
            cds_start = cds[0][0]
            for s, e in exons:
                if e < cds_start:
                    left_utrs.append((s, e))
                elif s < cds_start:
                    left_utrs.append((s, cds_start - 1))
                else:
                    break
            cds_end = cds[-1][-1]
            for s, e in exons[::-1]:
                if s > cds_end:
                    right_utrs.append((s, e))
                elif e > cds_end:
                    right_utrs.append((cds_end + 1, e))
                else:
                    break
        # return five_prime_utrs and three_prime_utrs
        if self.strand == '+':
            return left_utrs, right_utrs
        else: 
            return right_utrs, left_utrs
    
    def add_intron(self):
        intron_coords = self.get_introns()
        self.add_coordinates('intron', intron_coords)

    def add_utrs(self):
        five_prime_utrs, three_prime_utrs = self.get_utrs()
        self.add_coordinates('5_UTR', five_prime_utrs)
        self.add_coordinates('3_UTR', three_prime_utrs)
