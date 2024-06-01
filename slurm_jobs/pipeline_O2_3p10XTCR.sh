#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH -t 1:30:00
#SBATCH -p short
#SBATCH -o poreranger_job_%A.out
#SBATCH --account=chen_fec176

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3

python ~/nanoranger/pipeline.py --c 16 --i $1 --o $2 --e $3 --m 3p10XTCR
