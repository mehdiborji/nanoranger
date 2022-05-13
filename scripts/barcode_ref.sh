#!/bin/bash

ref_fasta=$1
genome_dir=$2
echo $ref_fasta
echo $genome_dir

#module load gcc/6.2.0
#module load star/2.7.9a

STAR \
--runMode genomeGenerate \
--runThreadN 4 \
--genomeDir $2 \
--genomeFastaFiles $1 \
--genomeSAindexNbases 6 \
--genomeChrBinNbits 7 \
--limitGenomeGenerateRAM 92000000000
