#!/bin/bash
cd $1
java -jar /usr/local/trimmomatic/latest/trimmomatic-0.32.jar PE -threads 4 -phred33 $1\_R1.fastq.gz $1\_R2.fastq.gz  $1\.filtered_R1.paired.fq  $1\.filtered_R1.unpaired.fq  $1\.filtered_R2.paired.fq  $1\.filtered_R2.unpaired.fq CROP:114 MINLEN:114
