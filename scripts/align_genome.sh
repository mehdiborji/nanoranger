#!/bin/bash

echo 'cores =' $1
echo 'ref =' $2
echo 'infile=' $3
echo 'outdir =' $4
echo 'sample =' $5

minimap2 -aY --eqx -x splice -t $1  --secondary=no --sam-hit-only $2 $3 > $4/$5_genome.sam

samtools sort -@$1 -o $4/$5_genome.bam $4/$5_genome.sam
rm $4/$5_genome.sam
samtools index -@$1 $4/$5_genome.bam
echo 'number of genome aligned reads = ' $(samtools view -@$1 -c $4/$5_genome.bam)

#samtools view -h -b -F 2308 -@16 $3/$5.bam > $3/$5_pri.bam
#samtools view -@16 -c $3/$5_pri.bam
#samtools index -@16 $3/$5_pri.bam