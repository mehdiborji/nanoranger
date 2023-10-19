#!/bin/bash

input_fastq=$1
genome_dir=$2
out_name=$3

echo $input_fastq
echo $genome_dir
echo $out_name

#module load gcc/6.2.0
#module load star/2.7.9a

STAR \
--runThreadN $4 \
--readFilesIn $input_fastq \
--genomeDir $genome_dir \
--alignIntronMax 1 \
--outFileNamePrefix $out_name \
--outSAMmode NoQS \
--outSAMattributes AS nM MD \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 0 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--scoreGenomicLengthLog2scale 0 \
--scoreDelBase -1 \
--scoreDelOpen 0 \
--scoreInsOpen 0 \
--scoreInsBase -1 \
--seedSearchStartLmax 4 \
--seedSearchStartLmaxOverLread 0.9 \
--alignEndsType EndToEnd \
--readNameSeparator space \
--readFilesCommand zcat


#--outSAMtype BAM Unsorted \

#module load gcc/9.2.0
#module load samtools/1.14

mv "$out_name"Aligned.out.sam "$out_name".sam
#samtools sort -@$4 -o "$out_name".bam "$out_name"Aligned.out.bam
#samtools index -@$4 "$out_name".bam

rm "$out_name"*out*