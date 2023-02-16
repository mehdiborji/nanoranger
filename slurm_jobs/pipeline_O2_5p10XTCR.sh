#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:30:00
#SBATCH -p priority
#SBATCH -o poreranger_job_%A.out

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3

python ~/poreranger/pipeline.py --c 16 --i $1 --o $2 --e $3 --s --m 5p10XTCR --human
