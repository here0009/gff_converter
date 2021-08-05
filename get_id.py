import pandas as pd
ref_table = pd.read_csv('data/Gene_transcript_new.list', sep = '\t')
print(ref_table.head())
query_names = pd.read_csv('test/gene_list_605.tsv', sep = '\t')['#Name'].to_list()
# print(query_ids)
print(len(query_names))
ref_dict = pd.Series(ref_table['ID'].values, index=ref_table['#Name']).to_dict()
print(ref_dict)
out_put_fhand = open('test/gene_id_list_605.tsv', 'w')
out_put_fhand.write('#Name\tID\n')
not_found = []
for query_name in query_names:
    if query_name in ref_dict:
        out_put_fhand.write('\t'.join([query_name, ref_dict[query_name]])+ '\n')
    else:
        not_found.append(query_name)
out_put_fhand.close()
print("The following names are not found:")
print(not_found)
