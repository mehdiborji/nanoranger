#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 0:20:00
#SBATCH -p short
#SBATCH -o poreranger_3pXCR_slideseq_job_%A.out
#SBATCH --account=chen_fec176

echo 'inputfq =' $1
echo 'outdir =' $2
echo 'sample =' $3
echo 'barcodes =' $4
echo 'trans_ref =' $5
echo 'mixcr_species =' $6

python ~/nanoranger/pipeline.py --c 16 --i $1 --o $2 --e $3 --m 3pXCR_slideseq --b $4 --t $5 --x $6 --s

# old $6 >> --h for human otherwise mouse references are selected