from collections import defaultdict
from typing import OrderedDict


bed_file = 'get_coord/final_output.bed'
fusion_file = 'get_coord/cosmic_fusion.tsv'
output_file = 'get_coord/fusion_gene_list.tsv'
bed_fhand = open(bed_file)
fusion_fhand = open(fusion_file)
output_fhand = open(output_file, 'w')


def get_query_names(bed_fhand):
    """
    Get the query gene names from the bed file
    """
    query_names = OrderedDict()
    for line in bed_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        query_names[line.split('\t')[4]] = []
    return query_names

def get_fusion_names(fusion_fhand):
    """
    Get the fusion gene names from the fusion file
    """
    fusion_genes = defaultdict(list)
    for line in fusion_fhand:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fusion_name = line.split('\t')[0]
        lst = fusion_name.split('-')  # partern of fusion gene
        for part in lst:
            if '_' in part:
                name, id = part.split('_')
            else:
                name = part
            fusion_genes[name].append(fusion_name)
    return fusion_genes

query_names = get_query_names(bed_fhand)
fusion_genes = get_fusion_names(fusion_fhand)
output_fhand.write('\t'.join(['#Name', 'Fusion']) + '\n')
for query in query_names:
    if query in fusion_genes:
        query_names[query] = fusion_genes[query]
        fusions = ';'.join(query_names[query])
        output_fhand.write('\t'.join([query, fusions]) + '\n')
print(query_names)



    
