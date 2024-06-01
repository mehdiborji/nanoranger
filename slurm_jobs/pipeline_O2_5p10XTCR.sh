#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 0:60:00
#SBATCH -p priority
#SBATCH -o poreranger_5p10XTCR_job_%A.out
#SBATCH --account=chen_fec176

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3
echo 'trans_ref =' $4
echo 'mixcr_species =' $5

python ~/nanoranger/pipeline.py --c 16 --i $1 --o $2 --e $3 --m 5p10XTCR --s --t $4 --x $5
