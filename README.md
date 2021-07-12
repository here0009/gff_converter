# GFF Converter

## Introduction

GFF Converter is used to convert NCBI or ENSEMBL style gff3 files into a more concise UCSC style gff file with the script `gff_converter.py`

NCBI style input data:
```
NC_000001.10	BestRefSeq	gene	65419	71585	.	+	.	ID=gene-OR4F5;Dbxref=GeneID:79501,HGNC:HGNC:14825;Name=OR4F5;description=olfactory receptor family 4 subfamily F member 5;gbkey=Gene;gene=OR4F5;gene_biotype=protein_coding
NC_000001.10	BestRefSeq	mRNA	65419	71585	.	+	.	ID=rna-NM_001005484.2;Parent=gene-OR4F5;Dbxref=GeneID:79501,Genbank:NM_001005484.2,HGNC:HGNC:14825;Name=NM_001005484.2;gbkey=mRNA;gene=OR4F5;product=olfactory receptor family 4 subfamily F member 5;tag=RefSeq Select;transcript_id=NM_001005484.2
NC_000001.10	BestRefSeq	exon	65419	65433	.	+	.	ID=exon-NM_001005484.2-1;Parent=rna-NM_001005484.2;Dbxref=GeneID:79501,Genbank:NM_001005484.2,HGNC:HGNC:14825;gbkey=mRNA;gene=OR4F5;product=olfactory receptor family 4 subfamily F member 5;tag=RefSeq Select;transcript_id=NM_001005484.2
```

ENSEMBL style input data:
```
1	ensembl_havana	gene	69091	70008	.	+	.	ID=gene:ENSG00000186092;Name=OR4F5;biotype=protein_coding;description=olfactory receptor%2C family 4%2C subfamily F%2C member 5 [Source:HGNC Symbol%3BAcc:14825];gene_id=ENSG00000186092;logic_name=ensembl_havana_gene;version=4
1	ensembl_havana	mRNA	69091	70008	.	+	.	ID=transcript:ENST00000335137;Parent=gene:ENSG00000186092;Name=OR4F5-001;biotype=protein_coding;ccdsid=CCDS30547.1;havana_transcript=OTTHUMT00000003223;havana_version=1;tag=basic;transcript_id=ENST00000335137;version=3
1	ensembl_havana	exon	69091	70008	.	+	.	Parent=transcript:ENST00000335137;Name=ENSE00002319515;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;exon_id=ENSE00002319515;rank=1;version=1
```

UCSC style output data；
```
chr1	NCBI	mRNA	65419	71585	.	+	.	ID=NM_001005484; name=OR4F5;
chr1	NCBI	5_UTR	65419	65433	.	+	.	Parent=NM_001005484;
chr1	NCBI	intron	65434	65519	.	+	.	Parent=NM_001005484;
chr1	NCBI	5_UTR	65520	65564	.	+	.	Parent=NM_001005484;
chr1	NCBI	CDS	65565	65573	.	+	0	Parent=NM_001005484;
chr1	NCBI	intron	65574	69036	.	+	.	Parent=NM_001005484;
chr1	NCBI	CDS	69037	70008	.	+	0	Parent=NM_001005484;
chr1	NCBI	3_UTR	70009	71585	.	+	.	Parent=NM_001005484;
```

It will also generate a table of (name, ncbi_id, ensembl_id) for each mRNA/transcript in  the input file. 
The missing value will be represented as 'None', and can be queried using `table_translate.py` or `org_Hs_eg_db_translate.R`.

The format of the output table is 
```
Name	NCBI	ENSEMBL
DDX11L1	NR_046018	ENST00000456328
WASH7P	NR_024540	None
OR4F5	NM_001005484	ENST00000641515
OR4F29	NM_001005221	ENST00000426406
```

`check_gff.py` can check if the generated gff file is valid:

1. Is every seq got its own id and name
2. Is every subseq got its parent, and is directy placed under its parent
3. If we join the coordinates of subseq end by end, its equal to the coordinates of the seq

## Data Preprocess and Preparation

* The input file should met the [gff3 file specifications](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/)
* Assembly report file shoule be provided for the seqid converion, it can be downloaded from [NCBI FTP site](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/)
* You can translate ensembl transcript ids to ncbi refseq ids, or vice versa based on the table [gene2ensembl](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) provided by NCBI
* A translation table for Homo sapiens can be generated by the following command:
`awk 'BEGIN { FS = "\t" } ; {if(NR == 1 || $1=="9606") print $0}' gene2ensembl > homo_sapiens_gene2ensembl`
* You can also make the translation using the R libray [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

## Install

You should install `python3` to use `gff_converter.py`, `coord_translate.py` and `table_translate.py`.
You should install `R` and the R library `org.Hs.eg.db` to use `org_Hs_eg_db_translate.R`.
A `yaml` file is provided to build a conda environment named `gfftools` by `conda/mamba` using the following command:

```bash
mamba env create -f envs/environment.yml
```

## Usage

The scripts can used with the following command:

```bash
gff_converter.py [-h] -i input_gff -o output_gff [-a assembly_report] -t output_name_id_table -s gff_style [--add_intron] [--add_utr]

table_translate.py [-h] -i input.tsv -o output.tsv -r translation_table [-t trans_type]

Rscript org_Hs_eg_db_translate.R input_file output_file trans_type

coord_translate.py [-h] -n ncbi.gff -e ensembl.gff -a assembly_report -o output.tsv

check_gff.py [-h] -i input.gff [-s gff_style]
```


## Test Data

A small test data is provided in the folder `test` for testing. You can use the following command to test the scripts:

```bash
python3 gff_converter.py -i test/ncbi_test.gff -o test/ncbi_test_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t test/ncbi_name_id.tsv -s NCBI  --add_intron --add_utr

python3 table_translate.py -i test/ncbi_name_id.tsv -o test/trans_ncbi_name_id.tsv  -r data/homo_sapiens_gene2ensembl -t n2e

Rscript org_Hs_eg_db_translate.R test/ncbi_name_id.tsv test/R_ncbi_name_id.tsv n2e

python3 coord_translate.py -n test/ncbi_test.gff -e test/ensembl_test.gff3 -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -o test/trans_coord_id.tsv

python3 check_gff.py -i test/ensembl_test_output.gff -s ENSEMBL
```


## Report

The scripts have been tested using the follwing files:

* NCBI style gff file : `GCF_000001405.25_GRCh37.p13_genomic.gff`
* ENSEMBL style gff file: `Homo_sapiens.GRCh37.87.Ensembl.gff3`

Script  for the testing:

```bash
python3 gff_converter.py -i gff_data/GCF_000001405.25_GRCh37.p13_genomic.gff -o gff_out/NCBI_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t  gff_out/ncbi_name_id.tsv -s NCBI --add_intron --add_utr
python3 gff_converter.py -i gff_data/Homo_sapiens.GRCh37.87.Ensembl.gff3  -o gff_out/ENSEMBL_output.gff -a data/GCF_000001405.25_GRCh37.p13_assembly_report.txt -t  gff_out/ensembl_name_id.tsv -s ENSEMBL --add_intron --add_utr

python3 table_translate.py -i gff_out/ncbi_name_id.tsv -o gff_out/trans_ncbi_name_id.tsv  -r data/homo_sapiens_gene2ensembl -t n2e 
python3 table_translate.py -i gff_out/ensembl_name_id.tsv -o gff_out/trans_ensembl_name_id.tsv -r data/homo_sapiens_gene2ensembl -t e2n

Rscript org_Hs_eg_db_translate.R gff_out/ncbi_name_id.tsv gff_out/R_trans_ncbi_name_id.tsv n2e
Rscript org_Hs_eg_db_translate.R gff_out/ensembl_name_id.tsv gff_out/R_trans_ensembl_name_id.tsv e2n

python3 check_gff.py -i gff_out/NCBI_output.gff -s NCBI
python3 check_gff.py -i gff_out/ENSEMBL_output.gff -s ENSEMBL
```

The following table is the Number Quries and TransId generated by `table_tranlate.py` and `org_Hs_eg_db_translate.R`
Method | TransType | Quries | TransId 
-------|-------|-------|-------
table_tranlate.py | n2e | 68879 | 43056
table_tranlate.py | e2n | 95160 | 40702
org_Hs_eg_db_translate.R | n2e | 68879 | 12638
org_Hs_eg_db_translate.R | e2n | 95160 | 9438
