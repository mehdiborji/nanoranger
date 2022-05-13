#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:00:00
#SBATCH -p priority
#SBATCH -o poreranger_job_%A.out

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3
echo 'cores =' $4

python ~/poreranger/pipeline.py --s \
--i $1 --o $2 --e $3 --c $4 \
--g ~/refs/GRCh38.primary_assembly.genome.fa.gz
