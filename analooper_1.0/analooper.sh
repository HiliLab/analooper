#!/bin/bash
# analooper bash script
# Hili Lab
# version 1.0
# 2018
#
# Step 1: Setup variables for run:
clear
# set defaults for the options
CODON=5
LENGTH=114
BARCODE=14
SPACER=5
# parse the options
while getopts 'b:c:l:s:' opt ; do
case $opt in
b) BARCODE=$OPTARG ;;
c) CODON=$OPTARG ;;
l) LENGTH=$OPTARG ;;
s) SPACER=$OPTARG ;;
esac
done
# skip processed options
shift $((OPTIND-1))
# check for mandatory positional parameters
if [ $# -lt 3 ]; then
echo "Usage: $0 [options] outputFolder Read1 Read2"
echo "Options: -b barcode_length (def: $BARCODE) | -c codon_length (def: $CODON) | -l duplex_length (def: $LENGTH)| -s spacer_length (def: $SPACER)"
exit 1
fi
DIR="$1"
LEFT="$2"
RIGHT="$3"
# assembly of terms
echo "Assembling with barcode_length=$BARCODE, codon_length=$CODON, duplex_length=$LENGTH, and spacer_length=$SPACER"
# Step 2: Use Trimmomatic
cd "$(dirname "$0")"
java -jar ./Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $2 $3  $1/filtered_R1.paired.fq  $1/filtered_R1.unpaired.fq  $1/filtered_R2.paired.fq  $1/filtered_R2.unpaired.fq CROP:$LENGTH MINLEN:$LENGTH
# Step 3: Use FLASH
./FLASH-1.2.11/flash -m $LENGTH -M $LENGTH -d $1 $1/filtered_R1.paired.fq $1/filtered_R2.paired.fq --interleaved-output
# Step 4: Use Extract Tag script from Duplex Sequencing Software Package
python ./Duplex-Sequencing/extract_tag_and_motifs_from_single_read.py --infile1 $1/out.extendedFrags.fastq --outfile1 $1/processed.fastq --barcode_length $BARCODE --spacer_length $SPACER > $1/tag_extraction.txt
# Step 5: Use Duplex Pair Identification script from Duplex Sequencing Package
python ./Duplex-Sequencing/id_duplex_pairs.py --infile $1/processed.fastq --outfile $1/processed.txt --barcode_length $BARCODE > $1/tb_duplex_pairing.log 2>&1
# Step 6: Determine Fidelity of LOOPER
python ./LOOPER-fidelity/calculate_fidelity.py --infile $1/processed.txt --outfile $1/fidelity.txt --motif_length $CODON
echo "Fidelity results are in: $DIR"
