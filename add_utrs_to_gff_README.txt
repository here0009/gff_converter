add_utrs_to_gff.py 
******************

Application to add explicit five_prime_UTR and three_prime_UTR features (as 
inferred from the gene, mRNA, exon and CDS features) to GFF3 data. The input 
file should be a GFF3 or gzipped GFF3. Note the app was developed based on 
the GFF3 formatting conventions used by NCBI, and may not work with all 
GFF3 files.

Example command:
add_utrs_to_gff.py input.gff3 > output_with_utrs.gff3

Requires python 2.7 or newer.

A sample input and output file is provided:
  /test/input/test.gff3
  /test/baseline/test.gff3 -- modified output

Any questions and comments should be addressed to:
info@ncbi.nlm.nih.gov
