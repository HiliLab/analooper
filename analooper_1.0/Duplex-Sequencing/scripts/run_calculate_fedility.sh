#!/bin/bash

cd $1

rm $1.res.txt

#python ~/tools/Duplex-Sequencing/extract_tag_and_motifs_from_single_read.py --infile1 $1.assembled.fastq  --outfile1 $1.processed.fasta  --barcode_length 14 --spacer_length 5 > tag_extraction.txt
#python ~/tools/Duplex-Sequencing/id_duplex_pairs.py --infile $1.processed.fasta --outfile $1.processed.txt  --barcode_length 14 >  tb_duplex_pairing.log 2>&1
python  ~/tools/Duplex-Sequencing/calculate_fidelity.py  --infile $1.processed.txt --outfile $1.res.txt  --motif_length 5

