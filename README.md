# Split long senqences by barcodes

Long sequences might contain barcodes. This tool looks for barcodes occurences 
using fuzzy search, and split each sequence by the right end of each detected barcode. 

## Example of arguments:
python main.py -i PAQ39826_pass_barcode07_8074fd48_d90929d6_0.fastq 
-bar "AAGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACCT, GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAAT"

